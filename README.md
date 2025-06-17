# Watershed Basin Shapefile Extraction

This script provides an automated solution for delineating watershed basins from a Digital Elevation Model (DEM) and a set of gauge locations. It automates the hydrological conditioning of the DEM, stream network extraction, and watershed delineation using the `WhiteboxTools` library with intelligent CRS optimization.

---

## 📖 Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Function Parameters](#function-parameters)
- [Processing Steps](#processing-steps)
- [Outputs](#outputs)
- [Features](#features)

---

## Requirements

The script requires the following Python libraries:

- `geopandas`
- `rasterio`
- `whitebox`
- `pandas`
- `pyproj`
- `pathlib` (built-in)
- `logging` (built-in)

---

## ⚙️ Installation

1. **Clone the repository:**

    ```bash
    git clone https://github.com/CooperBigFoot/extract-basin-shapefiles.git
    cd extract-basin-shapefiles
    ```

2. **Create the virtual environment and install dependencies using uv:**

    ```bash
    # Create the virtual environment
    uv venv

    # Install the dependencies from pyproject.toml
    uv sync
    ```

---

## Usage

The script provides two main functions:

1. `extract_basin_shapefiles` - Complete workflow starting from raw DEM
2. `extract_basin_shapefiles_from_filled_dem` - Workflow starting from pre-filled/conditioned DEM

### Basic Example

```python
from extract_basin_shapefiles import extract_basin_shapefiles, extract_basin_shapefiles_from_filled_dem

# Basic usage with automatic UTM zone detection
result = extract_basin_shapefiles(
    dem_path="path/to/dem.tif",
    gauges_shapefile_path="path/to/gauges.shp",
    output_directory="path/to/output"
)

# Advanced usage with custom parameters
result = extract_basin_shapefiles(
    dem_path="path/to/dem.tif",
    gauges_shapefile_path="path/to/gauges.shp", 
    output_directory="path/to/output",
    breach_dist=10,
    stream_extract_threshold=50,
    snap_dist=500,
    target_utm_zone=32642,  # Manual UTM zone specification
    big_tiff=True
)

# Usage with pre-filled DEM (skips conditioning step)
result = extract_basin_shapefiles_from_filled_dem(
    filled_dem_path="path/to/filled_dem.tif",
    gauges_shapefile_path="path/to/gauges.shp",
    output_directory="path/to/output"
)
```

---

## Function Parameters

### `extract_basin_shapefiles`

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `dem_path` | `str` or `Path` | Path to the input DEM raster file. | *Required* |
| `gauges_shapefile_path` | `str` or `Path` | Path to the input point shapefile containing gauge locations. | *Required* |
| `output_directory` | `str` or `Path` | Path to the directory where the final shapefiles will be saved. | *Required* |
| `breach_dist` | `int` | The maximum distance (in cells) for breaching depressions in the DEM. | `5` |
| `stream_extract_threshold` | `int` | The flow accumulation threshold for extracting the stream network. | `100` |
| `snap_dist` | `float` | The maximum distance (in map units) to snap gauge points to streams. | `1000` |
| `target_utm_zone` | `int` | Optional UTM zone EPSG code. If None, auto-detects optimal zone. | `None` |
| `big_tiff` | `bool` | Whether to create BigTIFF format for large files (>4GB). | `False` |
| `keep_temp_dir_on_fail` | `bool` | Whether to preserve temporary files on failure for debugging. | `False` |

### `extract_basin_shapefiles_from_filled_dem`

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `filled_dem_path` | `str` or `Path` | Path to the pre-filled/conditioned DEM raster file. | *Required* |
| `gauges_shapefile_path` | `str` or `Path` | Path to the input point shapefile containing gauge locations. | *Required* |
| `output_directory` | `str` or `Path` | Path to the directory where the final shapefiles will be saved. | *Required* |
| `stream_extract_threshold` | `int` | The flow accumulation threshold for extracting the stream network. | `100` |
| `snap_dist` | `float` | The maximum distance (in map units) to snap gauge points to streams. | `1000` |
| `target_utm_zone` | `int` | Optional UTM zone EPSG code. If None, auto-detects optimal zone. | `None` |
| `big_tiff` | `bool` | Whether to create BigTIFF format for large files (>4GB). | `False` |
| `keep_temp_dir_on_fail` | `bool` | Whether to preserve temporary files on failure for debugging. | `False` |

---

## 🔄 Processing Steps

The watershed extraction process follows these main steps:

### 1. **Input Validation & Setup**

- Validates input files and parameters
- Creates temporary working directories
- Sets up logging for process tracking

### 2. **CRS Optimization & Harmonization**

- **Automatic UTM zone detection**: Analyzes data bounds to select optimal UTM zone
- **Manual UTM specification**: Option to specify target UTM zone
- **CRS validation**: Ensures projected coordinate system with metric units
- Reprojects DEM and gauge points to harmonized CRS for accurate analysis

### 3. **DEM Conditioning** (skipped for filled DEM workflow)

- Uses least-cost breaching algorithm to remove depressions
- Configurable breach distance for different terrain types
- Creates hydrologically correct DEM for flow analysis

### 4. **Flow Analysis**

- Calculates D8 flow direction (pointer) using all available CPU cores
- Computes flow accumulation from flow direction
- Uses WhiteboxTools with maximum processor utilization

### 5. **Stream Network Extraction**

- Extracts stream network based on flow accumulation threshold
- Creates binary raster of stream pixels
- Higher thresholds = simpler stream networks

### 6. **Point Snapping**

- Snaps gauge points to nearest stream pixels using Jenson algorithm
- Fallback to original locations if snapping fails
- Ensures watershed outlets are properly positioned

### 7. **Watershed Delineation**

- Uses WhiteboxTools `UnnestBasins` function for multiple watershed processing
- Processes all gauge points simultaneously for efficiency
- Maintains 1:1 relationship between gauges and watersheds

### 8. **Results Integration**

- Maps gauge attributes to watersheds using VALUE-based linking
- Calculates watershed areas in km²
- Validates 1:1 gauge-to-watershed relationships

### 9. **Output Generation**

- Saves final watershed polygons shapefile
- Saves snapped gauge points shapefile
- Creates CRS information file for reference

### 10. **Cleanup**

- Removes temporary processing files (unless debugging)
- Preserves only final output files

---

## Outputs

The function generates the following output files in the specified output directory:

| File | Description |
| :--- | :--- |
| `watersheds.shp` | Final watershed polygons with gauge attributes and calculated areas |
| `gauges_snapped.shp` | Gauge points snapped to the stream network |
| `crs_info.txt` | Information about the coordinate reference system used |

### Watershed Attributes

The output watershed shapefile includes:

- **Original gauge attributes**: All attributes from the input gauge shapefile
- **VALUE**: Unique identifier linking watershed to original gauge
- **AREA_KM2**: Calculated watershed area in square kilometers

---

## 🚀 Features

### **Intelligent CRS Management**

- Automatic UTM zone detection for optimal hydrological analysis
- Geographic bounds analysis to select appropriate projection
- Manual UTM zone override capability
- Cross-hemisphere data handling

### **Robust Processing**

- Comprehensive error handling and logging
- Fallback methods for algorithm failures
- Input validation and bounds checking
- Temporary file management with cleanup

### **Flexible Workflow Options**

- Complete workflow from raw DEM
- Pre-filled DEM workflow for efficiency
- Configurable processing parameters
- BigTIFF support for large datasets

### **Quality Assurance**

- 1:1 gauge-to-watershed relationship validation
- Duplicate detection and reporting
- Area calculations and validation
- Comprehensive logging and status reporting

---

## 📊 Usage Tips

- **For large DEMs**: Enable `big_tiff=True` and ensure sufficient disk space
- **For detailed watersheds**: Use lower `stream_extract_threshold` values (e.g., 50-100)
- **For coarse analysis**: Use higher `stream_extract_threshold` values (e.g., 500-1000)
- **For debugging**: Set `keep_temp_dir_on_fail=True` to examine intermediate files
- **For efficiency**: Use pre-filled DEM workflow if you have multiple gauge sets for the same area
- **For accuracy**: Let the system auto-detect UTM zones unless you have specific projection requirements
