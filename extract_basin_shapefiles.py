import contextlib
import logging
import multiprocessing as mp
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    batch_size: int = 50,
    n_workers: int = None,
    use_watershed_function: bool = True,
) -> str:
    """Optimized watershed extraction with parallel processing and batching.

    Args:
        dem_path: Path to input DEM raster file.
        gauges_shapefile_path: Path to input gauges point shapefile.
        output_directory: Path to output directory for final results.
        breach_dist: Maximum distance for depression breaching. Defaults to 5.
        stream_extract_threshold: Flow accumulation threshold for stream
            extraction. Defaults to 100.
        snap_dist: Maximum distance for snapping pour points to streams.
            Defaults to 1000.
        batch_size: Number of watersheds to process in each batch. Defaults to 50.
        n_workers: Number of parallel workers. Defaults to CPU count - 1.
        use_watershed_function: Use watershed function instead of unnest_basins.

    Returns:
        Path to final watersheds shapefile.
    """
    # Convert paths to Path objects
    dem_path = Path(dem_path)
    gauges_path = Path(gauges_shapefile_path)
    output_dir = Path(output_directory)

    # Validate inputs
    _validate_inputs(dem_path, gauges_path, breach_dist, snap_dist)

    logger.info("Starting optimized watershed delineation process")

    # Setup working environment
    temp_dir = _setup_working_environment(output_dir)

    # Set number of workers
    if n_workers is None:
        n_workers = max(1, mp.cpu_count() - 1)

    logger.info(f"Using {n_workers} parallel workers with batch size {batch_size}")

    try:
        # Process watershed delineation with optimizations
        temp_dem, temp_gauges = _validate_and_harmonize_crs(dem_path, gauges_path, temp_dir)

        # Pre-process DEM once (most expensive operations)
        logger.info("Pre-processing DEM (this may take a while for large DEMs)")
        breached_dem = _condition_dem_optimized(temp_dir, temp_dem, breach_dist)
        flow_pointer, flow_accum = _calculate_flow_analysis_optimized(temp_dir, breached_dem)

        # Extract streams once
        streams_raster = _extract_streams_optimized(temp_dir, flow_accum, stream_extract_threshold)

        # Snap all points at once
        snapped_points = _snap_points_optimized(temp_dir, temp_gauges, streams_raster, snap_dist)

        # Split gauges into batches and process in parallel
        gauges_gdf = gpd.read_file(snapped_points)
        total_gauges = len(gauges_gdf)
        logger.info(f"Processing {total_gauges} gauges in batches of {batch_size}")

        # Create batches
        batches = []
        for i in range(0, total_gauges, batch_size):
            batch_gauges = gauges_gdf.iloc[i : i + batch_size].copy()
            batch_gauges["batch_id"] = i // batch_size
            batches.append((batch_gauges, i // batch_size))

        logger.info(f"Created {len(batches)} batches for processing")

        # Process batches in parallel
        if use_watershed_function:
            all_watersheds = _process_batches_with_watershed_function(
                batches, flow_pointer, flow_accum, temp_dir, n_workers
            )
        else:
            all_watersheds = _process_batches_parallel(batches, flow_pointer, temp_dir, n_workers)

        # Merge all watershed results
        if all_watersheds:
            final_watersheds = gpd.pd.concat(all_watersheds, ignore_index=True)

            # Add gauge attributes
            final_watersheds = _add_gauge_attributes(final_watersheds, gauges_gdf)

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


def _condition_dem_optimized(temp_dir: Path, temp_dem: Path, breach_dist: int) -> Path:
    """Optimized DEM conditioning with better parameters."""
    logger.info("Conditioning DEM with optimized settings")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = False
    wbt.max_procs = -1  # Use all available processors

    try:
        # Use fill_depressions_wang_and_liu (faster than breach for large areas)
        filled_dem = temp_dir / "dem_filled.tif"
        wbt.fill_depressions_wang_and_liu(dem=str(temp_dem), output=str(filled_dem))

        return filled_dem

    except Exception:
        # Fallback to original method if optimized version fails
        logger.warning("Optimized DEM conditioning failed, using original method")
        return _condition_dem_original(wbt, temp_dem, breach_dist)


def _condition_dem_original(wbt: whitebox.WhiteboxTools, temp_dem: Path, breach_dist: int) -> Path:
    """Original DEM conditioning method as fallback."""
    try:
        filled_dem = temp_dem.parent / "dem_filled.tif"
        wbt.fill_single_cell_pits(dem=str(temp_dem), output=str(filled_dem))

        breached_dem = temp_dem.parent / "dem_breached.tif"
        wbt.breach_depressions_least_cost(
            dem=str(filled_dem),
            output=str(breached_dem),
            dist=breach_dist,
            fill=True,
        )
        return breached_dem
    except Exception as e:
        raise RuntimeError(f"DEM conditioning failed: {e}") from e


def _calculate_flow_analysis_optimized(temp_dir: Path, filled_dem: Path) -> tuple[Path, Path]:
    """Optimized flow analysis with parallel processing."""
    logger.info("Calculating flow direction and accumulation (optimized)")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = False
    wbt.max_procs = -1

    try:
        # Use D8 flow pointer (fastest option)
        flow_pointer = temp_dir / "d8_pointer.tif"
        wbt.d8_pointer(dem=str(filled_dem), output=str(flow_pointer))

        # Flow accumulation with optimized settings
        flow_accum = temp_dir / "d8_flow_accum.tif"
        wbt.d8_flow_accumulation(i=str(flow_pointer), output=str(flow_accum), pntr=True)

        return flow_pointer, flow_accum

    except Exception as e:
        raise RuntimeError(f"Flow analysis failed: {e}") from e


def _extract_streams_optimized(temp_dir: Path, flow_accum: Path, threshold: int) -> Path:
    """Extract streams with optimized settings."""
    logger.info("Extracting stream network")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = False

    try:
        streams_raster = temp_dir / "streams.tif"
        wbt.extract_streams(
            flow_accum=str(flow_accum),
            output=str(streams_raster),
            threshold=threshold,
        )
        return streams_raster
    except Exception as e:
        raise RuntimeError(f"Stream extraction failed: {e}") from e


def _snap_points_optimized(temp_dir: Path, gauges: Path, streams: Path, snap_dist: float) -> Path:
    """Snap all points to streams at once."""
    logger.info("Snapping all gauge points to streams")

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(temp_dir)
    wbt.verbose = False

    try:
        snapped_points = temp_dir / "gauges_snapped.shp"
        wbt.jenson_snap_pour_points(
            pour_pts=str(gauges),
            streams=str(streams),
            output=str(snapped_points),
            snap_dist=snap_dist,
        )
        return snapped_points
    except Exception as e:
        raise RuntimeError(f"Point snapping failed: {e}") from e


def _process_batch_watersheds(args) -> gpd.GeoDataFrame:
    """Process a batch of watersheds - designed for multiprocessing."""
    batch_gauges, batch_id, flow_pointer_path, temp_dir_path, use_watershed_func = args

    # Create batch-specific working directory
    batch_dir = Path(temp_dir_path) / f"batch_{batch_id}"
    batch_dir.mkdir(exist_ok=True)

    wbt = whitebox.WhiteboxTools()
    wbt.work_dir = str(batch_dir)
    wbt.verbose = False

    try:
        # Save batch gauges to temporary file
        batch_gauges_path = batch_dir / f"batch_gauges_{batch_id}.shp"
        batch_gauges.to_file(batch_gauges_path)

        if use_watershed_func:
            # Use watershed function (faster for multiple points)
            watersheds_raster = batch_dir / f"watersheds_batch_{batch_id}.tif"
            wbt.watershed(
                d8_pntr=str(flow_pointer_path),
                pour_pts=str(batch_gauges_path),
                output=str(watersheds_raster),
            )

            # Convert to vector
            watersheds_vector = batch_dir / f"watersheds_batch_{batch_id}.shp"
            wbt.raster_to_vector_polygons(i=str(watersheds_raster), output=str(watersheds_vector))

            if watersheds_vector.exists():
                watershed_gdf = gpd.read_file(watersheds_vector)
                if not watershed_gdf.empty:
                    watershed_gdf["batch_id"] = batch_id
                    return watershed_gdf
        else:
            # Use unnest_basins (original method)
            watersheds_raster = batch_dir / f"watersheds_unnested_batch_{batch_id}.tif"
            wbt.unnest_basins(
                d8_pntr=str(flow_pointer_path),
                pour_pts=str(batch_gauges_path),
                output=str(watersheds_raster),
            )

            # Find all generated watershed files
            watershed_files = list(batch_dir.glob(f"watersheds_unnested_batch_{batch_id}_*.tif"))
            if not watershed_files and watersheds_raster.exists():
                watershed_files = [watersheds_raster]

            # Convert each to vector and combine
            batch_gdfs = []
            for i, raster_file in enumerate(watershed_files):
                vector_file = batch_dir / f"watershed_{batch_id}_{i}.shp"
                wbt.raster_to_vector_polygons(i=str(raster_file), output=str(vector_file))

                if vector_file.exists():
                    gdf = gpd.read_file(vector_file)
                    if not gdf.empty:
                        gdf["FID"] = i + (batch_id * 1000)  # Unique FID across batches
                        gdf["batch_id"] = batch_id
                        batch_gdfs.append(gdf)

            if batch_gdfs:
                return gpd.pd.concat(batch_gdfs, ignore_index=True)

        return gpd.GeoDataFrame()  # Return empty if nothing processed

    except Exception as e:
        logger.error(f"Batch {batch_id} processing failed: {e}")
        return gpd.GeoDataFrame()

    finally:
        # Clean up batch directory
        with contextlib.suppress(Exception):
            shutil.rmtree(batch_dir)


def _process_batches_with_watershed_function(
    batches: list, flow_pointer: Path, flow_accum: Path, temp_dir: Path, n_workers: int
) -> list[gpd.GeoDataFrame]:
    """Process batches using the watershed function (recommended)."""
    logger.info("Processing batches with watershed function")

    # Prepare arguments for multiprocessing
    args_list = [(batch_gauges, batch_id, flow_pointer, temp_dir, True) for batch_gauges, batch_id in batches]

    all_watersheds = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all batches
        future_to_batch = {executor.submit(_process_batch_watersheds, args): args[1] for args in args_list}

        # Collect results as they complete
        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                result = future.result()
                if not result.empty:
                    all_watersheds.append(result)
                    logger.info(f"Completed batch {batch_id} ({len(result)} watersheds)")
                else:
                    logger.warning(f"Batch {batch_id} produced no watersheds")
            except Exception as e:
                logger.error(f"Batch {batch_id} failed: {e}")

    return all_watersheds


def _process_batches_parallel(
    batches: list, flow_pointer: Path, temp_dir: Path, n_workers: int
) -> list[gpd.GeoDataFrame]:
    """Process batches using unnest_basins (fallback method)."""
    logger.info("Processing batches with unnest_basins")

    # Prepare arguments for multiprocessing
    args_list = [(batch_gauges, batch_id, flow_pointer, temp_dir, False) for batch_gauges, batch_id in batches]

    all_watersheds = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all batches
        future_to_batch = {executor.submit(_process_batch_watersheds, args): args[1] for args in args_list}

        # Collect results as they complete
        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                result = future.result()
                if not result.empty:
                    all_watersheds.append(result)
                    logger.info(f"Completed batch {batch_id} ({len(result)} watersheds)")
                else:
                    logger.warning(f"Batch {batch_id} produced no watersheds")
            except Exception as e:
                logger.error(f"Batch {batch_id} failed: {e}")

    return all_watersheds


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

        return final_watersheds

    except Exception as e:
        logger.warning(f"Failed to add gauge attributes: {e}")
        return watersheds_gdf


# Include original helper functions for validation, setup, etc.
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


def _validate_and_harmonize_crs(dem_path: Path, gauges_path: Path, temp_dir: Path) -> tuple[Path, Path]:
    """Validate input data and harmonize coordinate reference systems."""
    logger.info("Validating input data and harmonizing CRS")

    try:
        with rasterio.open(dem_path) as dem_src:
            dem_crs = dem_src.crs

        gauges_gdf = gpd.read_file(gauges_path)

        if gauges_gdf.crs != dem_crs:
            logger.info(f"Reprojecting gauges from {gauges_gdf.crs} to {dem_crs}")
            gauges_gdf = gauges_gdf.to_crs(dem_crs)

        temp_gauges_path = temp_dir / "gauges_reprojected.shp"
        gauges_gdf.to_file(temp_gauges_path)

        temp_dem_path = temp_dir / "input_dem.tif"
        shutil.copy2(dem_path, temp_dem_path)

        return temp_dem_path, temp_gauges_path

    except Exception as e:
        raise RuntimeError(f"CRS harmonization failed: {e}") from e


def _save_final_outputs(final_watersheds: gpd.GeoDataFrame, snapped_points: Path, output_dir: Path) -> Path:
    """Save final watershed and gauge shapefiles."""
    logger.info("Saving final outputs")

    try:
        final_output = output_dir / "watersheds.shp"
        final_watersheds.to_file(final_output)

        gauge_output = output_dir / "gauges_snapped.shp"
        snapped_gauges = gpd.read_file(snapped_points)
        snapped_gauges.to_file(gauge_output)

        logger.info(f"Final watersheds saved to: {final_output}")
        logger.info(f"Snapped gauges saved to: {gauge_output}")

        return final_output

    except Exception as e:
        raise RuntimeError(f"Failed to save outputs: {e}") from e


def _cleanup_temp_files(temp_dir: Path) -> None:
    """Clean up temporary processing files."""
    logger.info("Cleaning up temporary files")
    try:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
    except Exception as e:
        logger.warning(f"Failed to clean up temporary files: {e}")
