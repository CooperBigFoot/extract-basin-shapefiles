import logging
import shutil
import subprocess
from pathlib import Path

import geopandas as gpd
import rasterio
import whitebox
from pyproj import CRS

# Configure logging
logger = logging.getLogger(__name__)


def detect_optimal_utm_zone(bounds: tuple[float, float, float, float]) -> int:
    """Detect optimal UTM zone for given geographic bounds.

    Args:
        bounds: Geographic bounds as (min_lon, min_lat, max_lon, max_lat)

    Returns:
        EPSG code for optimal UTM zone

    Raises:
        ValueError: If bounds are invalid or in polar regions
    """
    min_lon, min_lat, max_lon, max_lat = bounds

    # Validate bounds
    if not (-180 <= min_lon <= 180 and -180 <= max_lon <= 180):
        raise ValueError(f"Invalid longitude bounds: {min_lon}, {max_lon}")
    if not (-90 <= min_lat <= 90 and -90 <= max_lat <= 90):
        raise ValueError(f"Invalid latitude bounds: {min_lat}, {max_lat}")
    if min_lon > max_lon or min_lat > max_lat:
        raise ValueError("Invalid bounds: min values must be less than max values")

    # Check for polar regions (where UTM is not appropriate)
    if max_lat > 84 or min_lat < -80:
        raise ValueError(
            "Data extends into polar regions where UTM is not suitable. "
            "Consider using a polar stereographic projection."
        )

    # Calculate center point for zone determination
    center_lon = (min_lon + max_lon) / 2
    center_lat = (min_lat + max_lat) / 2

    # Calculate UTM zone number (1-60)
    utm_zone = int((center_lon + 180) / 6) + 1

    # Handle edge cases
    if utm_zone > 60:
        utm_zone = 60
    elif utm_zone < 1:
        utm_zone = 1

    # Determine hemisphere
    if center_lat >= 0:
        # Northern hemisphere: EPSG:326XX
        epsg_code = 32600 + utm_zone
        hemisphere = "North"
    else:
        # Southern hemisphere: EPSG:327XX
        epsg_code = 32700 + utm_zone
        hemisphere = "South"

    # Check if data spans multiple zones (warn if > 6 degrees longitude)
    lon_span = max_lon - min_lon
    if lon_span > 6:
        logger.warning(
            f"Data spans {lon_span:.1f} degrees longitude (multiple UTM zones). "
            f"Using zone {utm_zone}{hemisphere[0]} for center point."
        )

    logger.info(f"Detected optimal UTM zone: {utm_zone}{hemisphere[0]} (EPSG:{epsg_code})")

    return epsg_code


def fix_whitebox_geotiff(input_path: Path, output_path: Path, target_crs: CRS = None, big_tiff: bool = False) -> Path:
    """Fix corrupted GeoTIFF metadata from WhiteboxTools output.

    Args:
        input_path: Path to input raster with corrupted metadata
        output_path: Path for output raster with fixed metadata
        target_crs: Target CRS to apply. If None, attempts to preserve existing CRS
        big_tiff: Whether to create BigTIFF format for large files
    """
    try:
        # Determine CRS to use
        if target_crs is not None:
            crs_string = target_crs.to_string()
        else:
            # Try to read existing CRS from input file
            try:
                with rasterio.open(input_path) as src:
                    if src.crs is not None:
                        crs_string = src.crs.to_string()
                    else:
                        logger.warning(f"No CRS found in {input_path}, using EPSG:4326 as fallback")
                        crs_string = "EPSG:4326"
            except Exception:
                logger.warning(f"Could not read CRS from {input_path}, using EPSG:4326 as fallback")
                crs_string = "EPSG:4326"

        # Use GDAL to rebuild clean metadata
        cmd = [
            "gdal_translate",
            "-co",
            "COMPRESS=LZW",
            "-co",
            "TILED=YES",
        ]

        # Add BigTIFF option if requested
        if big_tiff:
            cmd.extend(["-co", "BIGTIFF=YES"])

        cmd.extend(
            [
                "-a_srs",
                crs_string,
                str(input_path),
                str(output_path),
            ]
        )

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug(f"Fixed GeoTIFF metadata with CRS {crs_string}: {output_path}")
        return output_path

    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to fix GeoTIFF metadata: {e.stderr}")
        raise RuntimeError(f"GDAL translate failed: {e.stderr}") from e


def extract_basin_shapefiles(
    dem_path: str | Path,
    gauges_shapefile_path: str | Path,
    output_directory: str | Path,
    breach_dist: int = 5,
    stream_extract_threshold: int = 100,
    snap_dist: float = 1000,
    use_dinf: bool = True,
    target_utm_zone: int = None,
    big_tiff: bool = False,
) -> str:
    """Extract watershed basins for gauge points.

    Args:
        dem_path: Path to input DEM raster file.
        gauges_shapefile_path: Path to input gauges point shapefile.
        output_directory: Path to output directory for final results.
        breach_dist: Maximum distance for depression breaching. Defaults to 5.
        stream_extract_threshold: Flow accumulation threshold for stream
            extraction. Defaults to 100.
        snap_dist: Maximum distance for snapping pour points to streams.
            Defaults to 1000.
        use_dinf: Whether to use D-infinity algorithm (more accurate) or D8.
            Defaults to True.
        target_utm_zone: Optional UTM zone EPSG code (e.g., 32642). If None,
            will auto-detect optimal zone for geographic data. Defaults to None.
        big_tiff: Whether to create BigTIFF format for large files (>4GB).
            Defaults to False.

    Returns:
        Path to final watersheds shapefile.
    """
    # Convert paths to Path objects
    dem_path = Path(dem_path)
    gauges_path = Path(gauges_shapefile_path)
    output_dir = Path(output_directory)

    # Validate inputs
    _validate_inputs(dem_path, gauges_path, breach_dist, snap_dist)

    logger.info("Starting watershed delineation process")

    # Setup working environment
    temp_dir = _setup_working_environment(output_dir)

    try:
        # Process watershed delineation with CRS optimization
        temp_dem, temp_gauges = _validate_and_harmonize_crs(dem_path, gauges_path, temp_dir, target_utm_zone)

        # Pre-process DEM
        logger.info("Pre-processing DEM (this may take a while for large DEMs)")
        filled_dem = _condition_dem_optimized(temp_dir, temp_dem, breach_dist, big_tiff)
        flow_pointer, flow_accum, d8_pointer = _calculate_flow_analysis_optimized(
            temp_dir, filled_dem, use_dinf, big_tiff
        )

        # Extract and enhance streams
        streams_raster, stream_links = _extract_and_enhance_streams(
            temp_dir, flow_pointer, flow_accum, stream_extract_threshold, big_tiff
        )

        # Snap all points with fallback methods
        snapped_points = _snap_points_with_fallback(temp_dir, temp_gauges, streams_raster, flow_accum, snap_dist)

        # Process watersheds
        gauges_gdf = gpd.read_file(snapped_points)
        total_gauges = len(gauges_gdf)
        logger.info(f"Processing {total_gauges} gauges")

        # Process all watersheds using watershed function with D8 pointer
        all_watersheds = _process_watersheds_with_watershed_function(gauges_gdf, d8_pointer, temp_dir, big_tiff)

        # Process results
        if all_watersheds is not None and not all_watersheds.empty:
            # Add gauge attributes
            final_watersheds = _add_gauge_attributes(all_watersheds, gauges_gdf)

            # Calculate areas
            final_watersheds["AREA_KM2"] = final_watersheds.geometry.area / 1_000_000

            # Save final outputs
            final_output = _save_final_outputs(final_watersheds, snapped_points, output_dir)

            logger.info(f"Completed! Total watersheds: {len(final_watersheds)}")
            return str(final_output)
        else:
            raise RuntimeError("No watersheds were successfully processed")

    except Exception as e:
        logger.error(f"Watershed delineation failed: {e}")
        raise RuntimeError(f"Watershed delineation failed: {e}") from e

    finally:
        _cleanup_temp_files(temp_dir)


def _condition_dem_optimized(temp_dir: Path, temp_dem: Path, breach_dist: int, big_tiff: bool = False) -> Path:
    """Optimized DEM conditioning with improved depression handling."""
    logger.info("Conditioning DEM with optimized settings")

    # Get CRS from input DEM for metadata fixing
    with rasterio.open(temp_dem) as src:
        dem_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True
    wbt.max_procs = -1  # Use all available processors

    try:
        # Primary method: Use fill_depressions_wang_and_liu (fastest for large areas)
        filled_dem_raw = temp_dir / "dem_filled_raw.tif"
        wbt.fill_depressions_wang_and_liu(dem=str(temp_dem), output=str(filled_dem_raw))

        # Fix the metadata corruption from WhiteboxTools
        filled_dem = temp_dir / "dem_filled.tif"
        fix_whitebox_geotiff(filled_dem_raw, filled_dem, dem_crs, big_tiff)

        # Clean up raw file
        filled_dem_raw.unlink()

        return filled_dem

    except Exception as e:
        logger.error(f"Wang & Liu depression filling failed: {e}")
        # Fallback to Planchon & Darboux method
        logger.warning("Using Planchon & Darboux depression filling as fallback")
        return _condition_dem_planchon_darboux(wbt, temp_dem, temp_dir, breach_dist, dem_crs, big_tiff)


def _condition_dem_planchon_darboux(
    wbt: whitebox.WhiteboxTools, temp_dem: Path, temp_dir: Path, breach_dist: int, dem_crs: CRS, big_tiff: bool = False
) -> Path:
    """Improved fallback DEM conditioning using Planchon & Darboux method."""
    try:
        # Use more robust Planchon & Darboux depression filling
        filled_dem_raw = temp_dir / "dem_filled_pd_raw.tif"
        wbt.fill_depressions(dem=str(temp_dem), output=str(filled_dem_raw))

        # Fix metadata after filling
        filled_dem = temp_dir / "dem_filled_pd.tif"
        fix_whitebox_geotiff(filled_dem_raw, filled_dem, dem_crs, big_tiff)

        # Clean up raw file
        filled_dem_raw.unlink()

        return filled_dem

    except Exception as e:
        logger.error(f"Planchon & Darboux filling failed: {e}")
        # Final fallback to original breach method
        logger.warning("Using original breach method as final fallback")
        return _condition_dem_original_fixed(wbt, temp_dem, breach_dist, temp_dir, dem_crs, big_tiff)


def _condition_dem_original_fixed(
    wbt: whitebox.WhiteboxTools, temp_dem: Path, breach_dist: int, temp_dir: Path, dem_crs: CRS, big_tiff: bool = False
) -> Path:
    """Original DEM conditioning method with metadata fixing."""
    try:
        filled_dem_raw = temp_dir / "dem_filled_orig_raw.tif"
        wbt.fill_single_cell_pits(dem=str(temp_dem), output=str(filled_dem_raw))

        # Fix metadata after filling
        filled_dem = temp_dir / "dem_filled_orig.tif"
        fix_whitebox_geotiff(filled_dem_raw, filled_dem, dem_crs, big_tiff)

        breached_dem_raw = temp_dir / "dem_breached_raw.tif"
        wbt.breach_depressions_least_cost(
            dem=str(filled_dem),
            output=str(breached_dem_raw),
            dist=breach_dist,
            fill=True,
        )

        # Fix metadata after breaching
        breached_dem = temp_dir / "dem_breached.tif"
        fix_whitebox_geotiff(breached_dem_raw, breached_dem, dem_crs, big_tiff)

        # Clean up raw files
        filled_dem_raw.unlink()
        breached_dem_raw.unlink()

        return breached_dem

    except Exception as e:
        raise RuntimeError(f"DEM conditioning failed: {e}") from e


def _calculate_flow_analysis_optimized(
    temp_dir: Path, filled_dem: Path, use_dinf: bool = True, big_tiff: bool = False
) -> tuple[Path, Path, Path]:
    """Enhanced flow analysis with D-infinity option and metadata fixing."""
    algorithm = "D-infinity" if use_dinf else "D8"
    logger.info(f"Calculating flow direction and accumulation using {algorithm}")

    # Get CRS from input DEM for metadata fixing
    with rasterio.open(filled_dem) as src:
        dem_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True
    wbt.max_procs = -1

    try:
        if use_dinf:
            # Use D-infinity for more accurate flow routing
            flow_pointer_raw = temp_dir / "dinf_pointer_raw.tif"
            wbt.d_inf_pointer(dem=str(filled_dem), output=str(flow_pointer_raw))

            # Fix metadata for flow pointer
            flow_pointer = temp_dir / "dinf_pointer.tif"
            fix_whitebox_geotiff(flow_pointer_raw, flow_pointer, dem_crs, big_tiff)

            # D-infinity flow accumulation
            flow_accum_raw = temp_dir / "dinf_flow_accum_raw.tif"
            wbt.d_inf_flow_accumulation(i=str(flow_pointer), output=str(flow_accum_raw))

            # Fix metadata for flow accumulation
            flow_accum = temp_dir / "dinf_flow_accum.tif"
            fix_whitebox_geotiff(flow_accum_raw, flow_accum, dem_crs, big_tiff)

            # ALSO create D8 flow pointer for watershed function (required)
            d8_pointer_raw = temp_dir / "d8_pointer_raw.tif"
            wbt.d8_pointer(dem=str(filled_dem), output=str(d8_pointer_raw))

            d8_pointer = temp_dir / "d8_pointer.tif"
            fix_whitebox_geotiff(d8_pointer_raw, d8_pointer, dem_crs, big_tiff)

        else:
            # Fallback to D8 method
            flow_pointer_raw = temp_dir / "d8_pointer_raw.tif"
            wbt.d8_pointer(dem=str(filled_dem), output=str(flow_pointer_raw))

            # Fix metadata for flow pointer
            flow_pointer = temp_dir / "d8_pointer.tif"
            fix_whitebox_geotiff(flow_pointer_raw, flow_pointer, dem_crs, big_tiff)

            # D8 flow accumulation
            flow_accum_raw = temp_dir / "d8_flow_accum_raw.tif"
            wbt.d8_flow_accumulation(i=str(flow_pointer), output=str(flow_accum_raw), pntr=True)

            # Fix metadata for flow accumulation
            flow_accum = temp_dir / "d8_flow_accum.tif"
            fix_whitebox_geotiff(flow_accum_raw, flow_accum, dem_crs, big_tiff)

            d8_pointer = flow_pointer  # Same as flow_pointer for D8

        # Clean up raw files
        flow_pointer_raw.unlink()
        flow_accum_raw.unlink()
        if use_dinf:
            d8_pointer_raw.unlink()

        return flow_pointer, flow_accum, d8_pointer

    except Exception as e:
        if use_dinf:
            logger.error(f"D-infinity flow analysis failed: {e}")
            logger.warning("Falling back to D8 flow analysis")
            return _calculate_flow_analysis_optimized(temp_dir, filled_dem, use_dinf=False, big_tiff=big_tiff)
        else:
            raise RuntimeError(f"Flow analysis failed: {e}") from e


def _extract_and_enhance_streams(
    temp_dir: Path, flow_pointer: Path, flow_accum: Path, threshold: int, big_tiff: bool = False
) -> tuple[Path, Path]:
    """Extract and enhance stream network with stream links."""
    logger.info("Extracting and enhancing stream network")

    # Get CRS from flow accumulation for metadata fixing
    with rasterio.open(flow_accum) as src:
        flow_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    try:
        # Extract basic streams
        streams_raster_raw = temp_dir / "streams_raw.tif"
        wbt.extract_streams(
            flow_accum=str(flow_accum),
            output=str(streams_raster_raw),
            threshold=threshold,
        )

        # Fix metadata for streams
        streams_raster = temp_dir / "streams.tif"
        fix_whitebox_geotiff(streams_raster_raw, streams_raster, flow_crs, big_tiff)

        # Create stream links for better network topology
        stream_links_raw = temp_dir / "stream_links_raw.tif"
        try:
            wbt.stream_link_identifier(
                d8_pntr=str(flow_pointer), streams=str(streams_raster), output=str(stream_links_raw)
            )

            # Fix metadata for stream links
            stream_links = temp_dir / "stream_links.tif"
            fix_whitebox_geotiff(stream_links_raw, stream_links, flow_crs, big_tiff)

            # Clean up raw files
            streams_raster_raw.unlink()
            stream_links_raw.unlink()

            logger.info("Enhanced stream network with stream links")

        except Exception as e:
            logger.warning(f"Stream link enhancement failed: {e}")
            # Clean up raw file and use basic streams
            streams_raster_raw.unlink()
            stream_links = streams_raster  # Use basic streams as fallback

        return streams_raster, stream_links

    except Exception as e:
        raise RuntimeError(f"Stream extraction failed: {e}") from e


def _snap_points_with_fallback(temp_dir: Path, gauges: Path, streams: Path, flow_accum: Path, snap_dist: float) -> Path:
    """Snap points to streams with multiple fallback methods."""
    logger.info("Snapping gauge points to streams with fallback methods")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    # Try primary method: Jenson snap
    try:
        snapped_points = temp_dir / "gauges_snapped_jenson.shp"
        wbt.jenson_snap_pour_points(
            pour_pts=str(gauges),
            streams=str(streams),
            output=str(snapped_points),
            snap_dist=snap_dist,
        )

        # Verify the output exists and has features
        if snapped_points.exists():
            test_gdf = gpd.read_file(snapped_points)
            if not test_gdf.empty:
                logger.info("Successfully snapped points using Jenson method")
                return snapped_points

        logger.warning("Jenson snap produced no results, trying alternative method")

    except Exception as e:
        logger.warning(f"Jenson snap failed: {e}")

    # Fallback method: Snap to flow accumulation peaks
    try:
        snapped_points_alt = temp_dir / "gauges_snapped_flowaccum.shp"
        wbt.snap_pour_points(
            pour_pts=str(gauges), flow_accum=str(flow_accum), output=str(snapped_points_alt), snap_dist=snap_dist
        )

        # Verify the output
        if snapped_points_alt.exists():
            test_gdf = gpd.read_file(snapped_points_alt)
            if not test_gdf.empty:
                logger.info("Successfully snapped points using flow accumulation method")
                return snapped_points_alt

        logger.warning("Flow accumulation snap also failed")

    except Exception as e:
        logger.warning(f"Flow accumulation snap failed: {e}")

    # Final fallback: Return original points with warning
    logger.warning("All snapping methods failed, using original gauge locations")
    logger.warning("This may result in less accurate watershed delineation")

    # Copy original points to expected output location
    fallback_points = temp_dir / "gauges_original.shp"
    original_gdf = gpd.read_file(gauges)
    original_gdf.to_file(fallback_points)

    return fallback_points


def _process_watersheds_with_watershed_function(
    gauges_gdf: gpd.GeoDataFrame, d8_pointer: Path, temp_dir: Path, big_tiff: bool = False
) -> gpd.GeoDataFrame:
    """Process all watersheds using the watershed function with robust error handling."""
    logger.info("Processing watersheds with watershed function")

    # Get CRS from D8 pointer for metadata fixing
    with rasterio.open(d8_pointer) as src:
        pointer_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    try:
        # Save gauges to temporary file
        gauges_path = temp_dir / "all_gauges.shp"
        gauges_gdf.to_file(gauges_path)

        # Use watershed function with D8 pointer
        watersheds_raster_raw = temp_dir / "watersheds_raw.tif"

        # Check if the D8 pointer file exists and is valid
        if not d8_pointer.exists():
            raise FileNotFoundError(f"D8 pointer file not found: {d8_pointer}")

        logger.info(f"Using D8 pointer: {d8_pointer}")
        wbt.watershed(
            d8_pntr=str(d8_pointer),
            pour_pts=str(gauges_path),
            output=str(watersheds_raster_raw),
        )

        # Check if watershed output was created
        if not watersheds_raster_raw.exists():
            raise RuntimeError("Watershed function failed to create output file")

        # Fix metadata for watershed raster
        watersheds_raster = temp_dir / "watersheds.tif"
        fix_whitebox_geotiff(watersheds_raster_raw, watersheds_raster, pointer_crs, big_tiff)

        # Convert to vector
        watersheds_vector = temp_dir / "watersheds.shp"
        wbt.raster_to_vector_polygons(i=str(watersheds_raster), output=str(watersheds_vector))

        if watersheds_vector.exists():
            watershed_gdf = gpd.read_file(watersheds_vector)
            if not watershed_gdf.empty:
                logger.info(f"Successfully processed {len(watershed_gdf)} watersheds")
                return watershed_gdf

        logger.warning("Watershed processing produced no results, trying individual processing")
        return _process_watersheds_individually(gauges_gdf, d8_pointer, temp_dir, big_tiff)

    except Exception as e:
        logger.error(f"Watershed processing failed: {e}")
        logger.warning("Falling back to individual watershed processing")
        return _process_watersheds_individually(gauges_gdf, d8_pointer, temp_dir, big_tiff)


def _process_watersheds_individually(
    gauges_gdf: gpd.GeoDataFrame, d8_pointer: Path, temp_dir: Path, big_tiff: bool = False
) -> gpd.GeoDataFrame:
    """Process watersheds individually as fallback method."""
    logger.info("Processing watersheds individually")

    # Get CRS from D8 pointer for metadata fixing
    with rasterio.open(d8_pointer) as src:
        pointer_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    all_watersheds = []

    for idx, gauge in gauges_gdf.iterrows():
        try:
            logger.info(f"Processing individual watershed {idx + 1}/{len(gauges_gdf)}")

            # Save individual gauge
            individual_gauge = temp_dir / f"gauge_{idx}.shp"
            single_gauge_gdf = gpd.GeoDataFrame([gauge], crs=gauges_gdf.crs)
            single_gauge_gdf.to_file(individual_gauge)

            # Process individual watershed
            individual_watershed_raw = temp_dir / f"watershed_{idx}_raw.tif"
            wbt.watershed(
                d8_pntr=str(d8_pointer),
                pour_pts=str(individual_gauge),
                output=str(individual_watershed_raw),
            )

            if individual_watershed_raw.exists():
                # Fix metadata
                individual_watershed = temp_dir / f"watershed_{idx}.tif"
                fix_whitebox_geotiff(individual_watershed_raw, individual_watershed, pointer_crs, big_tiff)

                # Convert to vector
                individual_vector = temp_dir / f"watershed_{idx}.shp"
                wbt.raster_to_vector_polygons(i=str(individual_watershed), output=str(individual_vector))

                if individual_vector.exists():
                    watershed_gdf = gpd.read_file(individual_vector)
                    if not watershed_gdf.empty:
                        # Add FID to match with gauge
                        watershed_gdf["FID"] = gauge.get("FID", idx + 1)
                        all_watersheds.append(watershed_gdf)
                        logger.info(f"Successfully processed watershed {idx + 1}")

                # Clean up individual files
                individual_watershed_raw.unlink(missing_ok=True)
                individual_watershed.unlink(missing_ok=True)
                individual_gauge.unlink(missing_ok=True)
                individual_vector.unlink(missing_ok=True)

        except Exception as e:
            logger.warning(f"Failed to process individual watershed {idx + 1}: {e}")
            continue

    if all_watersheds:
        combined_watersheds = gpd.pd.concat(all_watersheds, ignore_index=True)
        logger.info(f"Successfully processed {len(combined_watersheds)} watersheds individually")
        return combined_watersheds
    else:
        logger.error("No watersheds could be processed individually")
        return gpd.GeoDataFrame()


def _add_gauge_attributes(watersheds_gdf: gpd.GeoDataFrame, gauges_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Add gauge attributes to watersheds."""
    logger.info("Adding gauge attributes to watersheds")

    try:
        # Ensure both have FID columns for joining
        if "FID" not in gauges_gdf.columns:
            gauges_gdf = gauges_gdf.reset_index()
            gauges_gdf["FID"] = gauges_gdf.index + 1

        # Join with gauge attributes (drop geometry from gauges to avoid conflicts)
        gauge_attrs = gauges_gdf.drop(columns=["geometry"])
        final_watersheds = watersheds_gdf.merge(gauge_attrs, on="FID", how="left")

        logger.info(f"Successfully added attributes to {len(final_watersheds)} watersheds")
        return final_watersheds

    except Exception as e:
        logger.warning(f"Failed to add gauge attributes: {e}")
        return watersheds_gdf


def _validate_inputs(dem_path: Path, gauges_path: Path, breach_dist: int, snap_dist: float) -> None:
    """Validate input parameters and file existence."""
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM file not found: {dem_path}")
    if not gauges_path.exists():
        raise FileNotFoundError(f"Gauges shapefile not found: {gauges_path}")
    if breach_dist <= 0:
        raise ValueError("breach_dist must be positive")
    if snap_dist <= 0:
        raise ValueError("snap_dist must be positive")


def _setup_working_environment(output_dir: Path) -> Path:
    """Setup working environment and create necessary directories."""
    logger.info("Setting up working environment")
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir / "temp_processing"
    temp_dir.mkdir(exist_ok=True)
    return temp_dir


def _validate_and_harmonize_crs(
    dem_path: Path, gauges_path: Path, temp_dir: Path, target_utm_zone: int = None
) -> tuple[Path, Path]:
    """Validate input data and harmonize coordinate reference systems with UTM optimization."""
    logger.info("Validating input data and optimizing CRS for hydrological analysis")

    try:
        # Read DEM metadata
        with rasterio.open(dem_path) as dem_src:
            dem_crs = dem_src.crs
            dem_bounds = dem_src.bounds

        # Read gauges
        gauges_gdf = gpd.read_file(gauges_path)
        logger.info(f"Loaded {len(gauges_gdf)} gauge records")

        # Determine target CRS
        target_crs = _determine_optimal_crs(dem_crs, dem_bounds, gauges_gdf, target_utm_zone)

        # Reproject DEM if needed
        temp_dem_path = temp_dir / "input_dem.tif"
        if dem_crs != target_crs:
            logger.info(f"Reprojecting DEM from {dem_crs} to {target_crs}")
            _reproject_raster(dem_path, temp_dem_path, target_crs)
        else:
            logger.info("DEM already in target CRS, copying to temp directory")
            shutil.copy2(dem_path, temp_dem_path)

        # Reproject gauges if needed
        temp_gauges_path = temp_dir / "gauges_reprojected.shp"
        if gauges_gdf.crs != target_crs:
            logger.info(f"Reprojecting gauges from {gauges_gdf.crs} to {target_crs}")
            gauges_gdf = gauges_gdf.to_crs(target_crs)
        else:
            logger.info("Gauges already in target CRS")

        gauges_gdf.to_file(temp_gauges_path)

        return temp_dem_path, temp_gauges_path

    except Exception as e:
        raise RuntimeError(f"CRS harmonization failed: {e}") from e


def _determine_optimal_crs(
    dem_crs: CRS, dem_bounds: rasterio.coords.BoundingBox, gauges_gdf: gpd.GeoDataFrame, target_utm_zone: int = None
) -> CRS:
    """Determine optimal CRS for hydrological analysis."""

    # If target UTM zone is manually specified
    if target_utm_zone is not None:
        logger.info(f"Using manually specified UTM zone: EPSG:{target_utm_zone}")
        target_crs = CRS.from_epsg(target_utm_zone)

        # Validate the manual zone makes sense
        _validate_manual_utm_zone(dem_bounds, gauges_gdf, target_utm_zone)

        return target_crs

    # Check if current CRS is suitable for hydrological analysis
    if _is_suitable_projected_crs(dem_crs):
        logger.info(f"Current CRS {dem_crs} is suitable for hydrological analysis")
        return dem_crs

    # Auto-detect optimal UTM zone for geographic data
    if dem_crs.is_geographic:
        logger.info("Data in geographic coordinates, detecting optimal UTM zone")

        # Get combined bounds in geographic coordinates
        combined_bounds = _get_combined_geographic_bounds(dem_bounds, gauges_gdf, dem_crs)

        # Detect optimal UTM zone
        optimal_epsg = detect_optimal_utm_zone(combined_bounds)
        target_crs = CRS.from_epsg(optimal_epsg)

        logger.info(f"Auto-selected UTM zone for optimal hydrological analysis: {target_crs}")
        return target_crs

    else:
        # Data is in some other projected system - assess if suitable
        logger.warning(
            f"Data in projected CRS {dem_crs}, but suitability unclear. "
            "Consider specifying target_utm_zone for optimal results."
        )
        return dem_crs


def _is_suitable_projected_crs(crs: CRS) -> bool:
    """Check if CRS is suitable for hydrological analysis (projected, metric units)."""
    if crs.is_geographic:
        return False

    # Check if units are metric
    try:
        axis_unit = crs.axis_info[0].unit_name.lower()
        return axis_unit in ["metre", "meter", "m"]
    except Exception:
        return False


def _get_combined_geographic_bounds(
    dem_bounds: rasterio.coords.BoundingBox, gauges_gdf: gpd.GeoDataFrame, dem_crs: CRS
) -> tuple[float, float, float, float]:
    """Get combined bounds of DEM and gauges in geographic coordinates."""

    # Convert DEM bounds to geographic if needed
    if dem_crs.is_geographic:
        dem_geo_bounds = (dem_bounds.left, dem_bounds.bottom, dem_bounds.right, dem_bounds.top)
    else:
        # Reproject DEM bounds to geographic
        from pyproj import Transformer

        transformer = Transformer.from_crs(dem_crs, CRS.from_epsg(4326), always_xy=True)
        min_x, min_y = transformer.transform(dem_bounds.left, dem_bounds.bottom)
        max_x, max_y = transformer.transform(dem_bounds.right, dem_bounds.top)
        dem_geo_bounds = (min_x, min_y, max_x, max_y)

    # Get gauges bounds in geographic coordinates
    if gauges_gdf.crs.is_geographic:
        gauges_bounds = gauges_gdf.total_bounds
    else:
        gauges_geo = gauges_gdf.to_crs(4326)
        gauges_bounds = gauges_geo.total_bounds

    # Combine bounds
    combined_bounds = (
        min(dem_geo_bounds[0], gauges_bounds[0]),  # min_lon
        min(dem_geo_bounds[1], gauges_bounds[1]),  # min_lat
        max(dem_geo_bounds[2], gauges_bounds[2]),  # max_lon
        max(dem_geo_bounds[3], gauges_bounds[3]),  # max_lat
    )

    return combined_bounds


def _validate_manual_utm_zone(
    dem_bounds: rasterio.coords.BoundingBox, gauges_gdf: gpd.GeoDataFrame, target_utm_zone: int
) -> None:
    """Validate that manually specified UTM zone is reasonable for the data."""

    try:
        # Extract zone number and hemisphere from EPSG code
        if 32601 <= target_utm_zone <= 32660:
            zone_num = target_utm_zone - 32600
            hemisphere = "North"
        elif 32701 <= target_utm_zone <= 32760:
            zone_num = target_utm_zone - 32700
            hemisphere = "South"
        else:
            logger.warning(f"UTM zone EPSG:{target_utm_zone} is not standard. Proceeding anyway.")
            return

        # Get approximate longitude range for this UTM zone
        zone_center_lon = (zone_num - 1) * 6 - 180 + 3
        zone_min_lon = zone_center_lon - 3
        zone_max_lon = zone_center_lon + 3

        # Check if data is reasonably close to this UTM zone
        # (allow some flexibility since UTM zones can handle data outside their strict bounds)
        combined_bounds = _get_combined_geographic_bounds(dem_bounds, gauges_gdf, CRS.from_epsg(4326))
        data_center_lon = (combined_bounds[0] + combined_bounds[2]) / 2

        lon_distance = abs(data_center_lon - zone_center_lon)
        if lon_distance > 12:  # More than 2 zones away
            logger.warning(
                f"Manual UTM zone {zone_num}{hemisphere[0]} may not be optimal. "
                f"Data center longitude {data_center_lon:.1f}° is {lon_distance:.1f}° "
                f"from zone center {zone_center_lon:.1f}°"
            )

        # Check hemisphere
        data_center_lat = (combined_bounds[1] + combined_bounds[3]) / 2
        if hemisphere == "North" and data_center_lat < -10:
            logger.warning(f"Using Northern hemisphere UTM zone for data centered at {data_center_lat:.1f}° latitude")
        elif hemisphere == "South" and data_center_lat > 10:
            logger.warning(f"Using Southern hemisphere UTM zone for data centered at {data_center_lat:.1f}° latitude")

    except Exception as e:
        logger.warning(f"Could not validate manual UTM zone: {e}")


def _reproject_raster(input_path: Path, output_path: Path, target_crs: CRS) -> None:
    """Reproject raster to target CRS using rasterio with cubic resampling for DEMs."""
    from rasterio.warp import Resampling, calculate_default_transform, reproject

    with rasterio.open(input_path) as src:
        transform, width, height = calculate_default_transform(src.crs, target_crs, src.width, src.height, *src.bounds)

        kwargs = src.meta.copy()
        kwargs.update({"crs": target_crs, "transform": transform, "width": width, "height": height})

        with rasterio.open(output_path, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.cubic,  # Better for elevation data
                )


def _save_final_outputs(final_watersheds: gpd.GeoDataFrame, snapped_points: Path, output_dir: Path) -> Path:
    """Save final watershed and gauge shapefiles with CRS information."""
    logger.info("Saving final outputs")

    try:
        final_output = output_dir / "watersheds.shp"
        final_watersheds.to_file(final_output)

        gauge_output = output_dir / "gauges_snapped.shp"
        snapped_gauges = gpd.read_file(snapped_points)
        snapped_gauges.to_file(gauge_output)

        # Save CRS information to a readme file
        crs_info_file = output_dir / "crs_info.txt"
        with open(crs_info_file, "w") as f:
            f.write("Coordinate Reference System Information\n")
            f.write("=====================================\n\n")
            f.write(f"All output files use CRS: {final_watersheds.crs}\n")
            f.write(f"CRS Name: {final_watersheds.crs.name}\n")
            f.write(f"EPSG Code: {final_watersheds.crs.to_epsg()}\n\n")
            f.write("This CRS was selected for optimal hydrological analysis accuracy.\n")

        logger.info(f"Final watersheds saved to: {final_output}")
        logger.info(f"Snapped gauges saved to: {gauge_output}")
        logger.info(f"CRS information saved to: {crs_info_file}")

        return final_output

    except Exception as e:
        raise RuntimeError(f"Failed to save outputs: {e}") from e


def _cleanup_temp_files(temp_dir: Path) -> None:
    """Clean up temporary processing files."""
    logger.info("Cleaning up temporary files")
    try:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
            logger.info("Temporary files cleaned up successfully")
    except Exception as e:
        logger.warning(f"Failed to clean up temporary files: {e}")
