## MODIS and VIIRS Fires Analysis
## Exposures to combustion sources near military operations in Iraq and Afghanistan using satellite observations
## Journal of Exposure Science and Environmental Epidemiology
## Author: Meredith Franklin
## Date: October 2, 2025
## -----
## Steps:
## MODIS C6.1 and VIIRS (data downloaded from https://firms.modaps.eosdis.nasa.gov/download/)
## Iraq: project latitude and longitude to UTM 38N
## Afghanistan: project latitude and longitude to UTM 42N
## Djibouti: project latitude and longitude to UTM 38N
## Split data into annual subsets
## Perform annual clustering with HDBSCAN 
## Set minPts based on tuning

library(dbscan)
library(sf)
library(doParallel)
library(foreach)

## Read data
## ----- Iraq
iraq_modis<-read.csv("fire_archive_M-C61_576385.csv")
## ----- Afghanistan
af_modis <- read.csv('fire_archive_M-C61_576384.csv')
## ----- Djibouti
dji_modis<-read.csv("DL_FIRE_M-C61_587727/fire_archive_M-C61_587727.csv")
dji_viirs_npp<-read.csv("fire_archive_SV-C2_587731.csv")

## ----- HDBSCAN function
hdbscan_by_year <- function(data, 
                            acq_date_col = "acq_date",
                            coords = c("longitude", "latitude"),
                            crs = 4326,
                            utm_crs = NULL,
                            minPts_fun = NULL,
                            years = NULL,
                            cores = 4) {
    
# --- Convert acq_date to year/month/day for annual splits
    data <- data %>%
      mutate(acq_date = as.character(!!sym(acq_date_col))) %>%
      separate(col = acq_date, into = c("year", "month", "day"), sep = "-", remove = FALSE)
    
# --- Convert to sf for CRS 
    if (!inherits(data, "sf")) {
      data <- st_as_sf(data, coords = coords, crs = crs, remove = FALSE)
    }
    
# --- Optionally reproject to UTM 
    if (!is.null(utm_crs)) {
      data <- st_transform(data, crs = utm_crs)
    }
    
# --- Get years 
    if (is.null(years)) {
      years <- sort(unique(data$year))
    }
    
# --- Set up parallel
    num_cores <- min(cores, parallel::detectCores() - 1)
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    clust_list <- foreach(i = years, .packages = c("dbscan", "sf", "dplyr"), .combine = 'c') %dopar% {
      cat("\nProcessing year:", i, "\n")
      data_i <- data[data$year == i, ]
      
      if (nrow(data_i) == 0) {
        cat(i, "No data available, skipping...\n")
        return(NULL)
      }
      
# --- Dynamic or fixed minPts and apply hdbscan
      minPts_value <- if (!is.null(minPts_fun)) minPts_fun(i) else 3
      
      clust_prob <- tryCatch({
        clust <- hdbscan(st_coordinates(data_i), minPts = minPts_value)
        
        result <- data.frame(
          data[data$year == i, ],
          cluster = clust$cluster,
          member_prob = clust$membership_prob
        )
        
        cat(i, "\nNoise fires:", sum(result$cluster == 0),
            "\nClustered fires:", sum(result$cluster != 0),
            "\nClusters:", max(result$cluster, na.rm = TRUE), "\n")
        
        list(i = result)
      }, error = function(e) {
        cat("Error processing year:", i, "-", conditionMessage(e), "\n")
        NULL
      })
      
      gc()
      return(clust_prob)
    }
    
    stopCluster(cl)
    
    names(clust_list) <- years
    all_clust <- do.call("rbind", clust_list)
    
    return(all_clust)
  }

## ----- Iraq
iraq_modis_clusters <- hdbscan_by_year(
  data = iraq_modis,
  coords = c("longitude", "latitude"),
  crs = 4326,
  utm_crs = 32638,
  acq_date_col = "acq_date",
  minPts_fun = function(y) {
    if (y == 2003) return(6) else return(5)
  },
  years = 2002:2012,
  cores = 4
)

##-----Afghanistan

afg_modis_clusters <- hdbscan_by_year(
  data = af_modis,
  coords = c("longitude", "latitude"),
  crs = 4326,
  utm_crs = 32642,  # UTM Zone 42N for Afghanistan
  acq_date_col = "acq_date",
  minPts_fun = function(y) {3},  
  years = 2002:2012,
  cores = 4
)

##-----Djibouti

dji_modis_clusters <- hdbscan_by_year(
  data = dji_modis,
  coords = c("longitude", "latitude"),
  crs = 4326,
  utm_crs = 32638,  # UTM Zone 38N for Djibouti
  acq_date_col = "acq_date",
  minPts_fun = function(y) {5},  
  years = 2016:2023,
  cores = 4
)

dji_viirs_clusters <- hdbscan_by_year(
  data = dji_viirs_npp,
  coords = c("longitude", "latitude"),
  crs = 4326,
  utm_crs = 32638,  # UTM Zone 38N for Djibouti
  acq_date_col = "acq_date",
  minPts_fun = function(y) {5},  
  years = 2016:2023,
  cores = 4
)

# Optional: filter out noise points (cluster == 0) 
iraq_modis_clusters_no_noise <- iraq_modis_clusters %>% filter(cluster !=0)
afg_modis_clusters_no_noise <- afg_modis_clusters %>% filter(cluster !=0)
dji_modis_clusters_no_noise <- dji_modis_clusters %>% filter(cluster != 0)
dji_viirs_clusters_no_noise <- dji_viirs_clusters %>% filter(cluster != 0)


