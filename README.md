MODIS-VIIRS-fires 🔥🛰️

This repository provides data and code for processing and clustering fire detections from MODIS and VIIRS satellite instruments, with applications in identifying persistent combustion sources such as flares, burn pits, and wildfires near military and civilian locations.

Reference paper is: Exposures to combustion sources near military operations in Iraq and Afghanistan using satellite observations

🔍 Overview

Fire detections are analyzed from:

- MODIS (Moderate Resolution Imaging Spectroradiometer)
- VIIRS (Visible Infrared Imaging Radiometer Suite)

Available for download from NASA FIRMS https://firms.modaps.eosdis.nasa.gov/download/

💻 Code:
- Reads and processes fire data (e.g., MODIS MCD14ML, VIIRS VNP14IMG).
- Converts data to spatial format and transforms to local UTM projections.
- Clusters detections using hierarchical density based scanning (HDBSCAN) to identify persistent sources.
- Supports comparison across years, regions, and satellite products.

📦 Features
- Annual spatial clustering of fire detections.
- Adjustable minPts values for tuning.
- Parallel processing with foreach and doParallel.

📈 Data Inputs
- MODIS .csv (e.g.,`fire_archive_M-C61_*.csv`)
- VIIRS .csv (e.g., `fire_archive_SV-C2_*.csv`)

⚙️ Requirements
- R packages: sf, dbscan, tidyverse, doParallel, foreach
