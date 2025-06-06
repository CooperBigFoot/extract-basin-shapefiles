import logging
import shutil
from pathlib import Path

import geopandas as gpd
import rasterio
import whitebox

# Configure logging
logger = logging.getLogger(__name__)


def extract_basin_shapefiles(
    dem_path: str | Path,
    gauges_shapefile_path: str | Path,
    output_directory: str | Path,
    breach_dist: int = 5,
    stream_extract_threshold: int = 100,
    snap_dist: float = 1000,
) -> str:
    """Extract watershed basin shapefiles from DEM and gauge locations.

    Args:
        dem_path: Path to input DEM raster file.
        gauges_shapefile_path: Path to input gauges point shapefile.
        output_directory: Path to output directory for final results.
        breach_dist: Maximum distance for depression breaching. Defaults to 5.
        stream_extract_threshold: Flow accumulation threshold for stream
            extraction. Defaults to 100.
        snap_dist: Maximum distance for snapping pour points to streams.
            Defaults to 1000.

    Returns:
        Path to final watersheds shapefile.

    Raises:
        FileNotFoundError: If input files do not exist.
        RuntimeError: If watershed delineation fails.
        ValueError: If input parameters are invalid.
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
    wbt = _initialize_whitebox_tools(temp_dir)

    try:
        # Process watershed delineation
        temp_dem, temp_gauges = _validate_and_harmonize_crs(dem_path, gauges_path, temp_dir)

        breached_dem = _condition_dem(wbt, temp_dem, breach_dist)
        flow_pointer, flow_accum = _calculate_flow_analysis(wbt, breached_dem)
        snapped_points = _extract_streams_and_snap_points(
            wbt, flow_accum, temp_gauges, stream_extract_threshold, snap_dist
        )

        watershed_files = _delineate_watersheds(wbt, flow_pointer, snapped_points)
        watershed_gdfs = _convert_to_vector_polygons(wbt, watershed_files)
        final_watersheds = _merge_and_attribute_watersheds(watershed_gdfs, snapped_points)

        # Save final outputs
        final_output = _save_final_outputs(final_watersheds, snapped_points, output_dir)
        _create_hillshade(wbt, temp_dem, output_dir)

        logger.info(f"Watershed delineation completed. Total watersheds: {len(final_watersheds)}")
        return str(final_output)

    except Exception as e:
        logger.error(f"Watershed delineation failed: {e}")
        raise RuntimeError(f"Watershed delineation failed: {e}") from e

    finally:
        _cleanup_temp_files(temp_dir)


def _validate_inputs(
    dem_path: Path,
    gauges_path: Path,
    breach_dist: int,
    snap_dist: float,
) -> None:
    """Validate input parameters and file existence.

    Args:
        dem_path: Path to DEM file.
        gauges_path: Path to gauges shapefile.
        breach_dist: Depression breaching distance.
        snap_dist: Snapping distance.

    Raises:
        FileNotFoundError: If input files do not exist.
        ValueError: If parameters are invalid.
    """
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM file not found: {dem_path}")

    if not gauges_path.exists():
        raise FileNotFoundError(f"Gauges shapefile not found: {gauges_path}")

    if breach_dist <= 0:
        raise ValueError("breach_dist must be positive")

    if snap_dist <= 0:
        raise ValueError("snap_dist must be positive")


def _setup_working_environment(output_dir: Path) -> Path:
    """Setup working environment and create necessary directories.

    Args:
        output_dir: Output directory path.

    Returns:
        Path to temporary processing directory.
    """
    logger.info("Setting up working environment")
    output_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = output_dir / "temp_processing"
    temp_dir.mkdir(exist_ok=True)

    return temp_dir


def _initialize_whitebox_tools(temp_dir: Path) -> whitebox.WhiteboxTools:
    """Initialize WhiteboxTools with proper configuration.

    Args:
        temp_dir: Temporary directory for processing.

    Returns:
        Configured WhiteboxTools instance.
    """
    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = False  # Use logging instead
    return wbt


def _validate_and_harmonize_crs(dem_path: Path, gauges_path: Path, temp_dir: Path) -> tuple[Path, Path]:
    """Validate input data and harmonize coordinate reference systems.

    Args:
        dem_path: Path to DEM file.
        gauges_path: Path to gauges shapefile.
        temp_dir: Temporary directory for processing.

    Returns:
        Tuple of (temp_dem_path, temp_gauges_path).

    Raises:
        RuntimeError: If CRS harmonization fails.
    """
    logger.info("Validating input data and harmonizing CRS")

    try:
        # Read DEM CRS
        with rasterio.open(dem_path) as dem_src:
            dem_crs = dem_src.crs

        # Read and reproject gauges if necessary
        gauges_gdf = gpd.read_file(gauges_path)

        if gauges_gdf.crs != dem_crs:
            logger.info(f"Reprojecting gauges from {gauges_gdf.crs} to {dem_crs}")
            gauges_gdf = gauges_gdf.to_crs(dem_crs)

        # Save files to temporary directory
        temp_gauges_path = temp_dir / "gauges_reprojected.shp"
        gauges_gdf.to_file(temp_gauges_path)

        temp_dem_path = temp_dir / "input_dem.tif"
        shutil.copy2(dem_path, temp_dem_path)

        return temp_dem_path, temp_gauges_path

    except Exception as e:
        raise RuntimeError(f"CRS harmonization failed: {e}") from e


def _condition_dem(wbt: whitebox.WhiteboxTools, temp_dem: Path, breach_dist: int) -> Path:
    """Hydrologically condition the DEM by filling pits and breaching."""
    logger.info("Conditioning DEM")

    try:
        # Fill single-cell pits
        filled_dem = temp_dem.parent / "dem_filled.tif"
        wbt.fill_single_cell_pits(dem=str(temp_dem), output=str(filled_dem))

        # Breach depressions using least-cost method
        breached_dem = temp_dem.parent / "dem_breached.tif"
        wbt.breach_depressions_least_cost(
            dem=str(filled_dem),
            output=str(breached_dem),
            dist=breach_dist,
            fill=True,
            fill_pits=True,
        )

        return breached_dem

    except Exception as e:
        raise RuntimeError(f"DEM conditioning failed: {e}") from e


def _calculate_flow_analysis(wbt: whitebox.WhiteboxTools, breached_dem: Path) -> tuple[Path, Path]:
    """Calculate flow direction and accumulation grids.

    Args:
        wbt: WhiteboxTools instance.
        breached_dem: Path to conditioned DEM file.

    Returns:
        Tuple of (flow_pointer_path, flow_accum_path).

    Raises:
        RuntimeError: If flow analysis fails.
    """
    logger.info("Calculating flow direction and accumulation")

    try:
        # Calculate D8 flow pointer
        flow_pointer = breached_dem.parent / "d8_pointer.tif"
        wbt.d8_pointer(dem=str(breached_dem), output=str(flow_pointer))

        # Calculate D8 flow accumulation
        flow_accum = breached_dem.parent / "d8_flow_accum.tif"
        wbt.d8_flow_accumulation(input=str(flow_pointer), output=str(flow_accum), pntr=True)

        return flow_pointer, flow_accum

    except Exception as e:
        raise RuntimeError(f"Flow analysis failed: {e}") from e


def _extract_streams_and_snap_points(
    wbt: whitebox.WhiteboxTools,
    flow_accum: Path,
    temp_gauges: Path,
    stream_extract_threshold: int,
    snap_dist: float,
) -> Path:
    """Extract streams and snap pour points to stream network.

    Args:
        wbt: WhiteboxTools instance.
        flow_accum: Path to flow accumulation raster.
        temp_gauges: Path to temporary gauges shapefile.
        stream_extract_threshold: Threshold for stream extraction.
        snap_dist: Maximum distance for snapping pour points.

    Returns:
        Path to snapped points shapefile.

    Raises:
        RuntimeError: If stream extraction or point snapping fails.
    """
    logger.info("Extracting streams and snapping pour points")

    try:
        # Extract streams from flow accumulation
        streams_raster = flow_accum.parent / "streams.tif"
        wbt.extract_streams(
            flow_accum=str(flow_accum),
            output=str(streams_raster),
            threshold=stream_extract_threshold,
        )

        # Snap pour points to streams
        snapped_points = flow_accum.parent / "gauges_snapped.shp"
        wbt.jenson_snap_pour_points(
            pour_pts=str(temp_gauges),
            streams=str(streams_raster),
            output=str(snapped_points),
            snap_dist=snap_dist,
        )

        return snapped_points

    except Exception as e:
        raise RuntimeError(f"Stream extraction or point snapping failed: {e}") from e


def _delineate_watersheds(wbt: whitebox.WhiteboxTools, flow_pointer: Path, snapped_points: Path) -> list[Path]:
    """Delineate watersheds using unnested basins approach.

    Args:
        wbt: WhiteboxTools instance.
        flow_pointer: Path to flow pointer raster.
        snapped_points: Path to snapped points shapefile.

    Returns:
        List of paths to watershed raster files.

    Raises:
        RuntimeError: If watershed delineation fails.
    """
    logger.info("Delineating watersheds")

    try:
        # Use unnest_basins to handle nested watersheds properly
        base_output = flow_pointer.parent / "watersheds_unnested.tif"
        wbt.unnest_basins(
            d8_pntr=str(flow_pointer),
            pour_pts=str(snapped_points),
            output=str(base_output),
        )

        # Find all unnested watershed raster files
        watershed_files = list(flow_pointer.parent.glob("watersheds_unnested_*.tif"))

        # If no individual files, check for single output
        if not watershed_files and base_output.exists():
            watershed_files = [base_output]

        if not watershed_files:
            raise RuntimeError("No watershed files were generated")

        return watershed_files

    except Exception as e:
        raise RuntimeError(f"Watershed delineation failed: {e}") from e


def _convert_to_vector_polygons(wbt: whitebox.WhiteboxTools, watershed_files: list[Path]) -> list[gpd.GeoDataFrame]:
    """Convert watershed rasters to vector polygons.

    Args:
        wbt: WhiteboxTools instance.
        watershed_files: List of watershed raster file paths.

    Returns:
        List of GeoDataFrames containing watershed polygons.

    Raises:
        RuntimeError: If raster to vector conversion fails.
    """
    logger.info("Converting watersheds to vector polygons")

    try:
        watershed_gdfs = []
        for i, raster_file in enumerate(watershed_files):
            vector_file = raster_file.parent / f"watershed_{i}.shp"

            # Convert raster to vector polygons
            wbt.raster_to_vector_polygons(i=str(raster_file), output=str(vector_file))

            # Read the resulting shapefile
            if vector_file.exists():
                gdf = gpd.read_file(vector_file)
                if not gdf.empty:
                    gdf["FID"] = i + 1  # Assign watershed ID
                    watershed_gdfs.append(gdf)

        if not watershed_gdfs:
            raise RuntimeError("No watersheds were successfully converted to vectors")

        return watershed_gdfs

    except Exception as e:
        raise RuntimeError(f"Vector conversion failed: {e}") from e


def _merge_and_attribute_watersheds(watershed_gdfs: list[gpd.GeoDataFrame], snapped_points: Path) -> gpd.GeoDataFrame:
    """Merge watershed polygons and add gauge attributes.

    Args:
        watershed_gdfs: List of watershed GeoDataFrames.
        snapped_points: Path to snapped points shapefile.

    Returns:
        Final GeoDataFrame with attributed watersheds.

    Raises:
        RuntimeError: If merging or attribution fails.
    """
    logger.info("Merging watersheds and adding attributes")

    try:
        # Merge all watershed polygons
        merged_watersheds = gpd.pd.concat(watershed_gdfs, ignore_index=True)

        # Read snapped gauges to get attributes
        snapped_gauges = gpd.read_file(snapped_points)

        # Add FID to snapped gauges if not present (WhiteboxTools starts from 0)
        if "FID" not in snapped_gauges.columns:
            snapped_gauges["FID"] = range(1, len(snapped_gauges) + 1)
        else:
            snapped_gauges["FID"] = snapped_gauges["FID"] + 1

        # Join watershed polygons with gauge attributes
        final_watersheds = merged_watersheds.merge(snapped_gauges.drop(columns=["geometry"]), on="FID", how="left")

        # Calculate watershed areas in km²
        final_watersheds["AREA_KM2"] = final_watersheds.geometry.area / 1_000_000

        return final_watersheds

    except Exception as e:
        raise RuntimeError(f"Watershed merging and attribution failed: {e}") from e


def _save_final_outputs(
    final_watersheds: gpd.GeoDataFrame,
    snapped_points: Path,
    output_dir: Path,
) -> Path:
    """Save final watershed and gauge shapefiles.

    Args:
        final_watersheds: Final watersheds GeoDataFrame.
        snapped_points: Path to snapped points shapefile.
        output_dir: Output directory.

    Returns:
        Path to final watersheds shapefile.

    Raises:
        RuntimeError: If saving outputs fails.
    """
    logger.info("Saving final outputs")

    try:
        # Save final watersheds shapefile
        final_output = output_dir / "watersheds.shp"
        final_watersheds.to_file(final_output)

        # Save final gauges shapefile
        gauge_output = output_dir / "gauges_snapped.shp"
        snapped_gauges = gpd.read_file(snapped_points)
        snapped_gauges.to_file(gauge_output)

        logger.info(f"Final watersheds saved to: {final_output}")
        logger.info(f"Snapped gauges saved to: {gauge_output}")

        return final_output

    except Exception as e:
        raise RuntimeError(f"Failed to save outputs: {e}") from e


def _create_hillshade(wbt: whitebox.WhiteboxTools, temp_dem: Path, output_dir: Path) -> None:
    """Create hillshade raster for visualization.

    Args:
        wbt: WhiteboxTools instance.
        temp_dem: Path to temporary DEM file.
        output_dir: Output directory.

    Raises:
        RuntimeError: If hillshade creation fails.
    """
    logger.info("Creating hillshade for visualization")

    try:
        hillshade_output = output_dir / "hillshade.tif"
        wbt.hillshade(
            dem=str(temp_dem),
            output=str(hillshade_output),
            azimuth=315.0,
            altitude=30.0,
            zfactor=1.0,
        )
        logger.info(f"Hillshade saved to: {hillshade_output}")

    except Exception as e:
        logger.warning(f"Hillshade creation failed: {e}")
        # Don't raise - hillshade is optional for visualization


def _cleanup_temp_files(temp_dir: Path) -> None:
    """Clean up temporary processing files.

    Args:
        temp_dir: Temporary directory to remove.
    """
    logger.info("Cleaning up temporary files")

    try:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
    except Exception as e:
        logger.warning(f"Failed to clean up temporary files: {e}")

