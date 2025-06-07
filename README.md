# Watershed Basin Shapefile Extraction

This script provides an optimized, parallel processing solution for delineating watershed basins from a Digital Elevation Model (DEM) and a set of gauge locations. It automates the hydrological conditioning of the DEM, stream network extraction, and watershed delineation using the `WhiteboxTools` library with advanced batching and multiprocessing capabilities.

---

## 📖 Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Function Parameters](#function-parameters)
- [Processing Steps](#processing-steps)
- [Outputs](#outputs)
- [Performance Features](#performance-features)

---

## Requirements

The script requires the following Python libraries:

- `geopandas`
- `rasterio`
- `whitebox`
- `multiprocessing` (built-in)
- `concurrent.futures` (built-in)
- `pathlib` (built-in)

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

The main function to use is `extract_basin_shapefiles`. This function takes paths to a DEM raster and a gauge location shapefile and performs the watershed analysis with optimized parallel processing.

To use it, follow the instructions in the `extract_basin_shapefiles` notebook.

### Basic Example

```python
from extract_basin_shapefiles import extract_basin_shapefiles

# Basic usage with default parameters
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
    batch_size=100,
    n_workers=8,
    stream_extract_threshold=50,
    snap_dist=500,
    use_watershed_function=True
)
```

---

## Function Parameters

The `extract_basin_shapefiles` function accepts the following arguments:

| Parameter | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `dem_path` | `str` or `Path` | Path to the input DEM raster file. | *Required* |
| `gauges_shapefile_path` | `str` or `Path` | Path to the input point shapefile containing gauge locations. | *Required* |
| `output_directory` | `str` or `Path` | Path to the directory where the final shapefiles will be saved. | *Required* |
| `breach_dist` | `int` | The maximum distance (in pixels) for breaching depressions in the DEM. | `5` |
| `stream_extract_threshold` | `int` | The flow accumulation threshold for extracting the stream network. Lower values = more detailed stream network. | `100` |
| `snap_dist` | `float` | The maximum distance (in map units) to snap gauge points to the extracted stream network. | `1000` |
| `batch_size` | `int` | Number of watersheds to process in each batch for parallel processing. | `50` |
| `n_workers` | `int` | Number of parallel workers to use. If `None`, uses CPU count - 1. | `None` |
| `use_watershed_function` | `bool` | Use optimized watershed function instead of unnest_basins method. Recommended for better performance. | `True` |

---

## 🔄 Processing Steps

The watershed extraction process follows these main steps:

### 1. **Input Validation & Setup**

- Validates input files and parameters
- Creates temporary working directories
- Sets up logging and multiprocessing environment

### 2. **CRS Harmonization**

- Reads DEM and gauge shapefiles
- Reprojects gauge points to match DEM coordinate system
- Copies inputs to temporary working directory

### 3. **DEM Conditioning**

- **Optimized approach**: Uses `fill_depressions_wang_and_liu` for faster processing
- **Fallback approach**: Uses traditional fill single cell pits + breach depressions
- Removes spurious sinks and conditions the DEM for flow analysis

### 4. **Flow Analysis**

- Calculates D8 flow direction (pointer) using all available CPU cores
- Computes flow accumulation from flow direction
- Uses optimized WhiteboxTools settings for maximum performance

### 5. **Stream Network Extraction**

- Extracts stream network based on flow accumulation threshold
- Creates binary raster of stream pixels
- Higher thresholds = simpler stream networks

### 6. **Point Snapping**

- Snaps all gauge points to the nearest stream pixels
- Uses Jenson snap pour points algorithm
- Ensures watershed outlets are located on the stream network

### 7. **Watershed Delineation (Parallel Processing)**

- Splits gauge points into batches for parallel processing
- **Method A** (default): Uses `watershed` function for multiple points simultaneously
- **Method B** (fallback): Uses `unnest_basins` for individual watershed processing
- Processes batches across multiple CPU cores

### 8. **Results Integration**

- Merges watershed results from all batches
- Adds original gauge attributes to watershed polygons
- Calculates watershed areas in km²
- Ensures unique identifiers across all watersheds

### 9. **Output Generation**

- Saves final watershed polygons shapefile
- Saves snapped gauge points shapefile
- Provides summary statistics and completion status

### 10. **Cleanup**

- Removes temporary processing files
- Cleans up batch-specific directories
- Preserves only final output files

---

## Outputs

The function generates the following output files in the specified output directory:

| File | Description |
| :--- | :--- |
| `watersheds.shp` | Final watershed polygons with gauge attributes and calculated areas |
| `gauges_snapped.shp` | Gauge points snapped to the stream network |

### Watershed Attributes

The output watershed shapefile includes:

- **Original gauge attributes**: All attributes from the input gauge shapefile
- **FID**: Unique identifier for each watershed
- **AREA_KM2**: Calculated watershed area in square kilometers
- **batch_id**: Batch processing identifier (for debugging)

---

## 🚀 Performance Features

### **Parallel Processing**

- Utilizes multiple CPU cores for watershed delineation
- Automatic detection of available processors
- Configurable batch sizes to optimize memory usage

### **Optimized Algorithms**

- Uses fastest available WhiteboxTools functions
- Wang & Liu depression filling algorithm for large DEMs
- D8 flow analysis with maximum processor utilization

### **Memory Management**

- Batch processing prevents memory overflow with large datasets
- Temporary file cleanup to manage disk space
- Efficient CRS handling and data harmonization

### **Error Handling**

- Comprehensive error logging and recovery
- Fallback methods for algorithm failures
- Individual batch failure isolation

### **Scalability**

- Handles datasets from dozens to thousands of gauge points
- Configurable parameters for different system capabilities
- Progress tracking and status reporting

---

## 📊 Performance Tips

- **For large DEMs**: Increase `batch_size` and ensure sufficient RAM
- **For many gauge points**: Reduce `batch_size` to prevent memory issues  
- **For detailed watersheds**: Lower `stream_extract_threshold` values
- **For faster processing**: Keep `use_watershed_function=True` (default)
- **For system optimization**: Set `n_workers` to match your CPU capabilities
