import logging
import shutil
import subprocess
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio
import rasterio.features
import whitebox
from pyproj import CRS, Transformer

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

    if not (-180 <= min_lon <= 180 and -180 <= max_lon <= 180):
        raise ValueError(f"Invalid longitude bounds: {min_lon}, {max_lon}")
    if not (-90 <= min_lat <= 90 and -90 <= max_lat <= 90):
        raise ValueError(f"Invalid latitude bounds: {min_lat}, {max_lat}")
    if min_lon > max_lon or min_lat > max_lat:
        raise ValueError("Invalid bounds: min values must be less than max values")

    if max_lat > 84 or min_lat < -80:
        raise ValueError(
            "Data extends into polar regions where UTM is not suitable. "
            "Consider using a polar stereographic projection."
        )

    center_lon = (min_lon + max_lon) / 2
    center_lat = (min_lat + max_lat) / 2

    utm_zone = int((center_lon + 180) / 6) + 1

    if utm_zone > 60:
        utm_zone = 60
    elif utm_zone < 1:
        utm_zone = 1

    if center_lat >= 0:
        epsg_code = 32600 + utm_zone
        hemisphere = "North"
    else:
        epsg_code = 32700 + utm_zone
        hemisphere = "South"

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
        if target_crs is not None:
            crs_string = target_crs.to_string()
        else:
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

        cmd = [
            "gdal_translate",
            "-co",
            "COMPRESS=LZW",
            "-co",
            "TILED=YES",
        ]

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
    target_utm_zone: int = None,
    big_tiff: bool = False,
    keep_temp_dir_on_fail: bool = False,
) -> str:
    """Extract watershed basins for gauge points using D8 flow analysis.

    Args:
        dem_path: Path to input DEM raster file.
        gauges_shapefile_path: Path to input gauges point shapefile.
        output_directory: Path to output directory for final results.
        breach_dist: Maximum distance for depression breaching. Defaults to 5.
        stream_extract_threshold: Flow accumulation threshold for stream
            extraction. Defaults to 100.
        snap_dist: Maximum distance for snapping pour points to streams.
            Defaults to 1000.
        target_utm_zone: Optional UTM zone EPSG code (e.g., 32642). If None,
            will auto-detect optimal zone for geographic data. Defaults to None.
        big_tiff: Whether to create BigTIFF format for large files (>4GB).
            Defaults to False.
        keep_temp_dir_on_fail: Whether to preserve temporary directory on failure
            for debugging purposes. Defaults to False.

    Returns:
        Path to final watersheds shapefile.
    """
    dem_path = Path(dem_path)
    gauges_path = Path(gauges_shapefile_path)
    output_dir = Path(output_directory)

    _validate_inputs(dem_path, gauges_path, breach_dist, snap_dist)

    logger.info("Starting watershed delineation process using D8 flow analysis")

    temp_dir = _setup_working_environment(output_dir)

    try:
        temp_dem, temp_gauges = _validate_and_harmonize_crs(dem_path, gauges_path, temp_dir, target_utm_zone)

        logger.info("Pre-processing DEM")
        filled_dem = _condition_dem(temp_dir, temp_dem, breach_dist, big_tiff)

        flow_pointer, flow_accum = _calculate_flow_analysis(temp_dir, filled_dem, big_tiff)

        streams_raster = _extract_streams(temp_dir, flow_accum, stream_extract_threshold, big_tiff)

        snapped_points = _snap_points(temp_dir, temp_gauges, streams_raster, snap_dist)

        gauges_gdf = gpd.read_file(snapped_points)
        total_gauges = len(gauges_gdf)
        logger.info(f"Processing {total_gauges} gauges")

        all_watersheds, gauge_mapping = _process_watersheds(gauges_gdf, flow_pointer, temp_dir, big_tiff)

        if all_watersheds is not None and not all_watersheds.empty:
            final_watersheds = _map_gauge_attributes_by_value(all_watersheds, gauge_mapping)
            final_watersheds["AREA_KM2"] = final_watersheds.geometry.area / 1_000_000

            final_output = _save_final_outputs(final_watersheds, snapped_points, output_dir)

            logger.info(f"Completed! Total watersheds: {len(final_watersheds)}")
            return str(final_output)
        else:
            raise RuntimeError("No watersheds were successfully processed")

    except Exception as e:
        logger.error(f"Watershed delineation failed: {e}")
        if keep_temp_dir_on_fail:
            logger.info(f"Preserving temporary directory for debugging: {temp_dir}")
        raise RuntimeError(f"Watershed delineation failed: {e}") from e

    finally:
        _cleanup_temp_files(temp_dir, keep_on_fail=keep_temp_dir_on_fail)


def _condition_dem(temp_dir: Path, temp_dem: Path, breach_dist: int, big_tiff: bool = False) -> Path:
    """DEM conditioning using least-cost breaching."""
    logger.info("Conditioning DEM using least-cost breaching")

    with rasterio.open(temp_dem) as src:
        dem_crs = src.crs
        # Get pixel size (assuming square pixels)
        pixel_size = abs(src.transform[0])  # meters per pixel

    # Convert breach distance from meters to cells
    search_dist_cells = max(1, int(breach_dist / pixel_size))

    logger.info(
        f"Using breach search distance of {search_dist_cells} cells ({breach_dist}m at {pixel_size}m resolution)"
    )

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True
    wbt.max_procs = -1

    breached_dem_raw = temp_dir / "dem_breached_raw.tif"
    wbt.breach_depressions_least_cost(
        dem=str(temp_dem),
        output=str(breached_dem_raw),
        dist=search_dist_cells,
        fill=True,
    )

    breached_dem = temp_dir / "dem_filled.tif"
    fix_whitebox_geotiff(breached_dem_raw, breached_dem, dem_crs, big_tiff)

    breached_dem_raw.unlink()
    return breached_dem


def _calculate_flow_analysis(temp_dir: Path, filled_dem: Path, big_tiff: bool = False) -> tuple[Path, Path]:
    """Calculate D8 flow direction and accumulation."""
    logger.info("Calculating D8 flow direction and accumulation")

    with rasterio.open(filled_dem) as src:
        dem_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True
    wbt.max_procs = -1

    flow_pointer_raw = temp_dir / "d8_pointer_raw.tif"
    wbt.d8_pointer(dem=str(filled_dem), output=str(flow_pointer_raw))

    flow_pointer = temp_dir / "d8_pointer.tif"
    fix_whitebox_geotiff(flow_pointer_raw, flow_pointer, dem_crs, big_tiff)

    flow_accum_raw = temp_dir / "d8_flow_accumulation_raw.tif"
    wbt.d8_flow_accumulation(i=str(flow_pointer), output=str(flow_accum_raw), pntr=True)

    flow_accum = temp_dir / "d8_flow_accum.tif"
    fix_whitebox_geotiff(flow_accum_raw, flow_accum, dem_crs, big_tiff)

    flow_pointer_raw.unlink()
    flow_accum_raw.unlink()

    return flow_pointer, flow_accum


def _extract_streams(temp_dir: Path, flow_accum: Path, threshold: int, big_tiff: bool = False) -> Path:
    """Extract stream network."""
    logger.info("Extracting stream network")

    with rasterio.open(flow_accum) as src:
        flow_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    streams_raster_raw = temp_dir / "streams_raw.tif"
    wbt.extract_streams(
        flow_accum=str(flow_accum),
        output=str(streams_raster_raw),
        threshold=threshold,
    )

    streams_raster = temp_dir / "streams.tif"
    fix_whitebox_geotiff(streams_raster_raw, streams_raster, flow_crs, big_tiff)

    streams_raster_raw.unlink()
    return streams_raster


def _snap_points(temp_dir: Path, gauges: Path, streams: Path, snap_dist: float) -> Path:
    """Snap points to streams."""
    logger.info("Snapping gauge points to streams")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    snapped_points = temp_dir / "gauges_snapped.shp"
    wbt.jenson_snap_pour_points(
        pour_pts=str(gauges),
        streams=str(streams),
        output=str(snapped_points),
        snap_dist=snap_dist,
    )

    if snapped_points.exists():
        test_gdf = gpd.read_file(snapped_points)
        if not test_gdf.empty:
            logger.info("Successfully snapped points")
            return snapped_points

    logger.warning("Snapping failed, using original gauge locations")
    fallback_points = temp_dir / "gauges_original.shp"
    original_gdf = gpd.read_file(gauges)
    original_gdf.to_file(fallback_points)
    return fallback_points


def _process_watersheds(
    gauges_gdf: gpd.GeoDataFrame, d8_pointer: Path, temp_dir: Path, big_tiff: bool = False
) -> tuple[gpd.GeoDataFrame, dict]:
    """Process watersheds using UnnestBasins tool.

    Returns:
        Tuple of (watersheds GeoDataFrame, gauge mapping dictionary)
    """
    logger.info("Processing watersheds using UnnestBasins tool")

    with rasterio.open(d8_pointer) as src:
        pointer_crs = src.crs

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = True

    # Create gauge mapping dictionary before processing
    # UnnestBasins assigns VALUE = record_index + 1 (1-based indexing)
    gauge_mapping = {}
    for idx, (_, gauge_row) in enumerate(gauges_gdf.iterrows()):
        value_id = idx + 1  # UnnestBasins uses 1-based indexing
        gauge_mapping[value_id] = gauge_row.to_dict()

    logger.info(f"Created gauge mapping for {len(gauge_mapping)} gauges with VALUE IDs 1-{len(gauge_mapping)}")

    # Save gauges without FID column - order is critical for UnnestBasins
    gauges_path = temp_dir / "all_gauges.shp"
    gauges_gdf.to_file(gauges_path)

    watersheds_unnested_base = temp_dir / "watersheds_unnested.tif"
    wbt.unnest_basins(d8_pntr=str(d8_pointer), pour_pts=str(gauges_path), output=str(watersheds_unnested_base))

    # Find all unnested basin files
    file_pattern = f"{watersheds_unnested_base.stem}_*.tif"
    file_list = list(temp_dir.glob(file_pattern))

    if not file_list:
        raise RuntimeError("UnnestBasins function failed to create output files")

    logger.info(f"Found {len(file_list)} unnested basin files to process")

    # Process each unnested basin file
    all_basin_parts = []
    for idx, ws_file in enumerate(sorted(file_list)):
        logger.info(f"Processing unnested basin file {idx + 1}/{len(file_list)}: {ws_file.name}")

        ws_raster_fixed = temp_dir / f"watersheds_unnested_{idx + 1}_fixed.tif"
        fix_whitebox_geotiff(ws_file, ws_raster_fixed, pointer_crs, big_tiff)

        ws_vector = temp_dir / f"watersheds_unnested_{idx + 1}.shp"

        try:
            wbt.raster_to_vector_polygons(i=str(ws_raster_fixed), output=str(ws_vector))
            if ws_vector.exists():
                ws_shp_tmp = gpd.read_file(ws_vector)
                if not ws_shp_tmp.empty:
                    ws_shp_tmp["geometry"] = ws_shp_tmp["geometry"].buffer(0)
                    all_basin_parts.append(ws_shp_tmp)
        except Exception as e:
            logger.warning(f"Failed to convert basin {idx + 1}: {e}")

        ws_raster_fixed.unlink(missing_ok=True)
        ws_file.unlink(missing_ok=True)

    if all_basin_parts:
        ws_shp = pd.concat(all_basin_parts, ignore_index=True)
        logger.info(f"Successfully processed {len(ws_shp)} total watershed polygons")
        return ws_shp, gauge_mapping
    else:
        raise RuntimeError("No valid watersheds were processed")


def _map_gauge_attributes_by_value(watersheds_gdf: gpd.GeoDataFrame, gauge_mapping: dict) -> gpd.GeoDataFrame:
    """Map gauge attributes to watersheds using VALUE field from UnnestBasins.

    Args:
        watersheds_gdf: Watershed polygons with VALUE field from UnnestBasins
        gauge_mapping: Dictionary mapping VALUE IDs to gauge attributes

    Returns:
        Watersheds GeoDataFrame with gauge attributes added
    """
    logger.info("Mapping gauge attributes to watersheds using VALUE field")

    # Validate VALUE field exists
    if "VALUE" not in watersheds_gdf.columns:
        raise RuntimeError("VALUE field not found in watershed data - UnnestBasins may have failed")

    # Get unique VALUES and validate mapping
    watershed_values = set(watersheds_gdf["VALUE"].dropna().astype(int))
    available_values = set(gauge_mapping.keys())

    logger.info(f"Found {len(watershed_values)} unique watershed VALUE IDs: {sorted(watershed_values)}")
    logger.info(f"Available gauge mappings for VALUE IDs: {sorted(available_values)}")

    # Check for missing mappings
    missing_values = watershed_values - available_values
    if missing_values:
        logger.warning(
            f"WARNING: {len(missing_values)} watersheds have VALUE IDs without gauge mappings: {sorted(missing_values)}"
        )

    # Check for unused mappings
    unused_values = available_values - watershed_values
    if unused_values:
        logger.info(
            f"INFO: {len(unused_values)} gauges did not generate watersheds (VALUE IDs: {sorted(unused_values)})"
        )

    # Create result DataFrame with gauge attributes
    result_rows = []
    successful_mappings = 0
    failed_mappings = 0

    for idx, row in watersheds_gdf.iterrows():
        try:
            value_id = int(row["VALUE"])

            if value_id in gauge_mapping:
                # Get gauge attributes
                gauge_attrs = gauge_mapping[value_id].copy()

                # Combine watershed geometry with gauge attributes
                result_row = gauge_attrs.copy()
                result_row["geometry"] = row["geometry"]
                result_row["VALUE"] = value_id

                result_rows.append(result_row)
                successful_mappings += 1
            else:
                logger.warning(f"No gauge mapping found for watershed VALUE {value_id}")
                failed_mappings += 1

        except (ValueError, TypeError) as e:
            logger.warning(f"Invalid VALUE field in watershed {idx}: {row.get('VALUE', 'N/A')} - {e}")
            failed_mappings += 1

    if not result_rows:
        raise RuntimeError("No valid watershed-gauge mappings were created")

    # Create final GeoDataFrame
    result_gdf = gpd.GeoDataFrame(result_rows, crs=watersheds_gdf.crs)

    logger.info(f"Successfully mapped {successful_mappings} watersheds to gauges")
    if failed_mappings > 0:
        logger.warning(f"Failed to map {failed_mappings} watersheds")

    # Validate 1:1 mapping
    gauge_counts = result_gdf.groupby("VALUE").size()
    duplicate_values = gauge_counts[gauge_counts > 1]
    if not duplicate_values.empty:
        logger.warning(f"WARNING: Found duplicate VALUE assignments: {duplicate_values.to_dict()}")
    else:
        logger.info("✓ Confirmed 1:1 watershed-to-gauge mapping (no duplicates)")

    return result_gdf


def _add_gauge_attributes(watersheds_gdf: gpd.GeoDataFrame, gauges_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """DEPRECATED: This function has been replaced by _map_gauge_attributes_by_value().

    The spatial join approach caused duplicate gauge assignments.
    Use VALUE-based mapping instead for reliable 1:1 relationships.
    """
    logger.warning(
        "DEPRECATED: _add_gauge_attributes() should not be called. Use _map_gauge_attributes_by_value() instead."
    )
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
    shutil.rmtree(temp_dir, ignore_errors=True)  # Clean previous runs
    temp_dir.mkdir(exist_ok=True)
    return temp_dir


def _validate_and_harmonize_crs(
    dem_path: Path, gauges_path: Path, temp_dir: Path, target_utm_zone: int = None
) -> tuple[Path, Path]:
    """Validate input data and harmonize coordinate reference systems with UTM optimization."""
    logger.info("Validating input data and optimizing CRS for hydrological analysis")

    try:
        with rasterio.open(dem_path) as dem_src:
            dem_crs = dem_src.crs
            dem_bounds = dem_src.bounds

        gauges_gdf = gpd.read_file(gauges_path)
        logger.info(f"Loaded {len(gauges_gdf)} gauge records")

        target_crs = _determine_optimal_crs(dem_crs, dem_bounds, gauges_gdf, target_utm_zone)

        temp_dem_path = temp_dir / "input_dem.tif"
        if dem_crs != target_crs:
            logger.info(f"Reprojecting DEM from {dem_crs} to {target_crs}")
            _reproject_raster(dem_path, temp_dem_path, target_crs)
        else:
            logger.info("DEM already in target CRS, copying to temp directory")
            shutil.copy2(dem_path, temp_dem_path)

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

    if target_utm_zone is not None:
        logger.info(f"Using manually specified UTM zone: EPSG:{target_utm_zone}")
        target_crs = CRS.from_epsg(target_utm_zone)

        _validate_manual_utm_zone(dem_bounds, gauges_gdf, target_utm_zone)

        return target_crs

    if _is_suitable_projected_crs(dem_crs):
        logger.info(f"Current CRS {dem_crs} is suitable for hydrological analysis")
        return dem_crs

    if dem_crs.is_geographic:
        logger.info("Data in geographic coordinates, detecting optimal UTM zone")

        combined_bounds = _get_combined_geographic_bounds(dem_bounds, gauges_gdf, dem_crs)

        optimal_epsg = detect_optimal_utm_zone(combined_bounds)
        target_crs = CRS.from_epsg(optimal_epsg)

        logger.info(f"Auto-selected UTM zone for optimal hydrological analysis: {target_crs}")
        return target_crs

    else:
        logger.warning(
            f"Data in projected CRS {dem_crs}, but suitability unclear. "
            "Consider specifying target_utm_zone for optimal results."
        )
        return dem_crs


def _is_suitable_projected_crs(crs: CRS) -> bool:
    """Check if CRS is suitable for hydrological analysis (projected, metric units)."""
    if crs.is_geographic:
        return False

    try:
        axis_unit = crs.axis_info[0].unit_name.lower()
        return axis_unit in ["metre", "meter", "m"]
    except Exception:
        return False


def _get_combined_geographic_bounds(
    dem_bounds: rasterio.coords.BoundingBox, gauges_gdf: gpd.GeoDataFrame, dem_crs: CRS
) -> tuple[float, float, float, float]:
    """Get combined bounds of DEM and gauges in geographic coordinates."""

    if dem_crs.is_geographic:
        dem_geo_bounds = (dem_bounds.left, dem_bounds.bottom, dem_bounds.right, dem_bounds.top)
    else:
        transformer = Transformer.from_crs(dem_crs, CRS.from_epsg(4326), always_xy=True)
        min_x, min_y = transformer.transform(dem_bounds.left, dem_bounds.bottom)
        max_x, max_y = transformer.transform(dem_bounds.right, dem_bounds.top)
        dem_geo_bounds = (min_x, min_y, max_x, max_y)

    if gauges_gdf.crs.is_geographic:
        gauges_bounds = gauges_gdf.total_bounds
    else:
        gauges_geo = gauges_gdf.to_crs(4326)
        gauges_bounds = gauges_geo.total_bounds

    combined_bounds = (
        min(dem_geo_bounds[0], gauges_bounds[0]),
        min(dem_geo_bounds[1], gauges_bounds[1]),
        max(dem_geo_bounds[2], gauges_bounds[2]),
        max(dem_geo_bounds[3], gauges_bounds[3]),
    )

    return combined_bounds


def _validate_manual_utm_zone(
    dem_bounds: rasterio.coords.BoundingBox, gauges_gdf: gpd.GeoDataFrame, target_utm_zone: int
) -> None:
    """Validate that manually specified UTM zone is reasonable for the data."""

    try:
        if 32601 <= target_utm_zone <= 32660:
            zone_num = target_utm_zone - 32600
            hemisphere = "North"
        elif 32701 <= target_utm_zone <= 32760:
            zone_num = target_utm_zone - 32700
            hemisphere = "South"
        else:
            logger.warning(f"UTM zone EPSG:{target_utm_zone} is not standard. Proceeding anyway.")
            return

        zone_center_lon = (zone_num - 1) * 6 - 180 + 3

        combined_bounds = _get_combined_geographic_bounds(dem_bounds, gauges_gdf, CRS.from_epsg(4326))
        data_center_lon = (combined_bounds[0] + combined_bounds[2]) / 2

        lon_distance = abs(data_center_lon - zone_center_lon)
        if lon_distance > 12:
            logger.warning(
                f"Manual UTM zone {zone_num}{hemisphere[0]} may not be optimal. "
                f"Data center longitude {data_center_lon:.1f}° is {lon_distance:.1f}° "
                f"from zone center {zone_center_lon:.1f}°"
            )

        data_center_lat = (combined_bounds[1] + combined_bounds[3]) / 2
        if hemisphere == "North" and data_center_lat < -10:
            logger.warning(f"Using Northern hemisphere UTM zone for data centered at {data_center_lat:.1f}° latitude")
        elif hemisphere == "South" and data_center_lat > 10:
            logger.warning(f"Using Southern hemisphere UTM zone for data centered at {data_center_lat:.1f}° latitude")

    except Exception as e:
        logger.warning(f"Could not validate manual UTM zone: {e}")


logger = logging.getLogger(__name__)


def _reproject_raster(input_path: Path, output_path: Path, target_crs: CRS) -> None:
    """Efficiently reproject raster using GDAL with VRT for memory efficiency."""

    logger.info("Reprojecting raster using GDAL VRT method")

    # Step 1: Create VRT with target projection (no data copying)
    vrt_path = output_path.with_suffix(".vrt")

    vrt_cmd = [
        "gdalwarp",
        "-of",
        "VRT",
        "-t_srs",
        target_crs.to_string(),
        "-r",
        "bilinear",
        "-multi",
        "-wo",
        "NUM_THREADS=ALL_CPUS",
        "-wm",
        "2GB",
        str(input_path),
        str(vrt_path),
    ]

    try:
        subprocess.run(vrt_cmd, capture_output=True, text=True, check=True)

        # Step 2: Convert VRT to compressed GeoTIFF
        translate_cmd = [
            "gdal_translate",
            "-of",
            "GTiff",
            "-co",
            "COMPRESS=LZW",
            "-co",
            "TILED=YES",
            "-co",
            "BLOCKXSIZE=512",
            "-co",
            "BLOCKYSIZE=512",
            "-co",
            "BIGTIFF=IF_SAFER",
            str(vrt_path),
            str(output_path),
        ]

        subprocess.run(translate_cmd, capture_output=True, text=True, check=True)

        # Clean up VRT
        vrt_path.unlink(missing_ok=True)

        logger.info(f"Reprojection completed: {output_path}")

    except subprocess.CalledProcessError as e:
        # Clean up on failure
        vrt_path.unlink(missing_ok=True)
        raise RuntimeError(f"GDAL reprojection failed: {e.stderr}") from e


def _save_final_outputs(final_watersheds: gpd.GeoDataFrame, snapped_points: Path, output_dir: Path) -> Path:
    """Save final watershed and gauge shapefiles with CRS information."""
    logger.info("Saving final outputs")

    try:
        final_output = output_dir / "watersheds.shp"
        final_watersheds.to_file(final_output)

        gauge_output = output_dir / "gauges_snapped.shp"
        shutil.copy(snapped_points, gauge_output)  # Copy the snapped points shapefile and its sidecars

        crs_info_file = output_dir / "crs_info.txt"
        with open(crs_info_file, "w") as f:
            f.write("Coordinate Reference System Information\n")
            f.write("=====================================\n\n")
            f.write(f"All output files use CRS: {final_watersheds.crs.to_wkt(pretty=True)}\n")
            if final_watersheds.crs and final_watersheds.crs.to_epsg():
                f.write(f"CRS Name: {final_watersheds.crs.name}\n")
                f.write(f"EPSG Code: {final_watersheds.crs.to_epsg()}\n\n")
            f.write("This CRS was selected for optimal hydrological analysis accuracy.\n")

        logger.info(f"Final watersheds saved to: {final_output}")
        logger.info(f"Snapped gauges saved to: {gauge_output}")
        logger.info(f"CRS information saved to: {crs_info_file}")

        return final_output

    except Exception as e:
        raise RuntimeError(f"Failed to save outputs: {e}") from e


def _cleanup_temp_files(temp_dir: Path, keep_on_fail: bool = False) -> None:
    """Clean up temporary processing files.

    Args:
        temp_dir: Path to temporary directory
        keep_on_fail: Whether to preserve temp files (used when called after failure)
    """
    if keep_on_fail:
        logger.info(f"Temporary files preserved at: {temp_dir}")
        return

    logger.info("Cleaning up temporary files")
    try:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
            logger.info("Temporary files cleaned up successfully")
    except Exception as e:
        logger.warning(f"Failed to clean up temporary files: {e}")
