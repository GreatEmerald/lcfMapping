# Markov paper
# Process WUR change data

# Access libraries
library(sf)
library(pbapply)
library(ranger)

source("utils/filterBands.R")
source("utils/loadData.R")
source("utils/extractDates.R")
source("utils/dataManagement.R")
source("utils/harmonicsFunctions.R")
source("RFfunction.R")

# Link to data folder (for yourself)
linkData <- "data/"


# Link to landsat gpkg of wur change dataset 
linkRawValidation = paste0(linkData, "/raw/WURChange20152019_Landsat8_TS.gpkg")
nameBands <- st_layers(linkRawValidation)

# Read in b2 to filter
b2 = st_read(linkRawValidation, nameBands$name[2])
st_geometry(b2)=NULL

# change colnames dates
dates = extractDates()
NewColDates = paste0("X", gsub("-", ".", dates))
colnames(b2)[4:194] = NewColDates

# Filter blue band
b2Filtered = filterBands(b2, smoothLoessPlot, dates)

# Apply on other bands
b1 <- st_read(linkRawValidation, nameBands$name[1])
b3 <- st_read(linkRawValidation, nameBands$name[3])
b4 <- st_read(linkRawValidation, nameBands$name[4])
b5 <- st_read(linkRawValidation, nameBands$name[5])
b6 <- st_read(linkRawValidation, nameBands$name[6])
b7 <- st_read(linkRawValidation, nameBands$name[7])
st_geometry(b1) = NULL
st_geometry(b3) = NULL
st_geometry(b4) = NULL
st_geometry(b5) = NULL
st_geometry(b6) = NULL
st_geometry(b7) = NULL
colnames(b1)[4:194] = NewColDates
colnames(b3)[4:194] = NewColDates
colnames(b4)[4:194] = NewColDates
colnames(b5)[4:194] = NewColDates
colnames(b6)[4:194] = NewColDates
colnames(b7)[4:194] = NewColDates

mean(is.na(b2))
mean(is.na(b2Filtered))
 
# nothing saved from this, need to rerun the filter on b2 next time

# Apply b2 filter to other bands
b2Matrix = as.matrix(b2Filtered[,NewColDates])

ApplyFilter = function(b) 
{
    bFiltered = b
    temp = as.matrix(bFiltered)[,NewColDates]
    temp[is.na(b2Matrix)] = NA
    bFiltered[,NewColDates] = temp
    return(bFiltered)
}
b1Filtered = ApplyFilter(b1)
b3Filtered = ApplyFilter(b3)
b4Filtered = ApplyFilter(b4)
b5Filtered = ApplyFilter(b5)
b6Filtered = ApplyFilter(b6)
b7Filtered = ApplyFilter(b7)

mean(is.na(b7[,NewColDates]))
mean(is.na(b7Filtered[,NewColDates]))
mean(is.na(b2Filtered))

b1FilteredSF <- DFtoSF(b1Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b2Filtered = cbind("location_id" = b2$location_id,
                   "sample_x" = b2$sample_x,
                   "sample_y" = b2$sample_y, b2Filtered)
b2FilteredSF <- DFtoSF(b2Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b3FilteredSF <- DFtoSF(b3Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b4FilteredSF <- DFtoSF(b4Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b5FilteredSF <- DFtoSF(b5Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b6FilteredSF <- DFtoSF(b6Filtered, coords = c("sample_x","sample_y"), validation = TRUE)
b7FilteredSF <- DFtoSF(b7Filtered, coords = c("sample_x","sample_y"), validation = TRUE)

# Save as one gpkg with mulitple layers
st_write(b1FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b1")
st_write(b2FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b2")
st_write(b3FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b3")
st_write(b4FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b4")
st_write(b5FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b5")
st_write(b6FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b6")
st_write(b7FilteredSF, paste0(linkData,"processed/WURchangeFiltered.gpkg"), "b7")

## Calc VIs
# Convert to numeric
b4Filtered = as.data.frame(sapply(b4Filtered[,NewColDates], as.numeric))
b5Filtered = as.data.frame(sapply(b5Filtered[,NewColDates], as.numeric))

# Calculate NDVI
ndvi = (b5Filtered - b4Filtered) / (b5Filtered + b4Filtered)
plot(as.numeric(apply(ndvi, 2, function(x){mean(x, na.rm = TRUE)})), ylab="mean ndvi")

# Save ndvi as gpkg
temp = cbind("location_id" = b1Filtered$location_id,
             "sample_x" = b1Filtered$sample_x,
             "sample_y" = b1Filtered$sample_y,
             ndvi)
ndviSF <- DFtoSF(temp, coords = c("sample_x","sample_y"), validation = TRUE) # first source
st_write(ndviSF, paste0(linkData, "/processed/WURchangeVIs.gpkg"), "NDVI")


##
# Get temporal harmonic metrics

# Apply function to get the harmonics of NDVI
test = t(pbapply(as.matrix(ndvi), 1, getHarmonics))
HarmMetrics = cbind(b1Filtered[, c("location_id","sample_x", "sample_y")], test)

# Change colnames
names(HarmMetrics)[4:(length(HarmMetrics))]= c("min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save harmonics 
HarmMetricsSF <- DFtoSF(HarmMetrics, coords = c("sample_x","sample_y"), validation = TRUE)
st_write(HarmMetricsSF, paste0(linkData, "/processed/WURchangeHarmonics.gpkg"), "NDVI")


HarmMetrics = st_read(paste0(linkData, "/processed/WURchangeHarmonics.gpkg"), "NDVI")
wurChangeCSV = read.csv(paste0(linkData, "/raw/reference_global_100m_orig&change_year2015-2019_20210407.csv"))
sum(wurChangeCSV$location_id %in% HarmMetrics$location_id)
sum(HarmMetrics$location_id %in% wurChangeCSV$location_id)
change2015 = subset(wurChangeCSV, dataYear=="2015")
change2016 = subset(wurChangeCSV, dataYear=="2016")
change2017 = subset(wurChangeCSV, dataYear=="2017")
change2018= subset(wurChangeCSV, dataYear=="2018")

dataVali = loadChangeValidationData()
#loadChangeValidationData()


### 
# Now split time series into yearly datasets
b4 = st_read(paste0(linkData, "processed/WURchangeFiltered.gpkg"), "b4")
b5 = st_read(paste0(linkData, "processed/WURchangeFiltered.gpkg"), "b5")
st_geometry(b4)=NULL
st_geometry(b5)=NULL

# 2015 #
# 68 observations / dates
b42015 = b4[,colnames(b4)[grepl("2014|2015|2016", colnames(b4))]]
b52015 = b5[,colnames(b5)[grepl("2014|2015|2016", colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b42015, as.numeric))
b5Temp = as.data.frame(sapply(b52015, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl("2014|2015|2016",dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics, cl=parallel::detectCores()))
HarmMetrics = cbind(b4[, c("location_id","sample_x", "sample_y")], HarmMetrics)

# change colnames
names(HarmMetrics)[4:(length(HarmMetrics))]= c("min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2015 features as GPKG
temp = DFtoSF(HarmMetrics, coords = c("sample_x","sample_y"), validation = TRUE)
st_write(temp, paste0(linkData, "/processed/2015/WURchangeHarmonics.gpkg"), "NDVI")

# 2016 #
# 69 observations / dates
b42016 = b4[,colnames(b4)[grepl("2015|2016|2017", colnames(b4))]]
b52016 = b5[,colnames(b5)[grepl("2015|2016|2017", colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b42016, as.numeric))
b5Temp = as.data.frame(sapply(b52016, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl("2015|2016|2017",dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))
HarmMetrics = cbind(b4[, c("location_id","sample_x", "sample_y")], HarmMetrics)

# change colnames
names(HarmMetrics)[4:(length(HarmMetrics))]= c("min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2016 features as GPKG
temp = DFtoSF(HarmMetrics, coords = c("sample_x","sample_y"), validation = TRUE)
st_write(temp, paste0(linkData, "/processed/2016/WURchangeHarmonics.gpkg"), "NDVI")

# 2017 #
# 69 observations / dates
b42017 = b4[,colnames(b4)[grepl("2016|2017|2018", colnames(b4))]]
b52017 = b5[,colnames(b5)[grepl("2016|2017|2018", colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b42017, as.numeric))
b5Temp = as.data.frame(sapply(b52017, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl("2016|2017|2018",dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))
HarmMetrics = cbind(b4[, c("location_id","sample_x", "sample_y")], HarmMetrics)

# change colnames
names(HarmMetrics)[4:(length(HarmMetrics))]= c("min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2017 features as GPKG
temp = DFtoSF(HarmMetrics, coords = c("sample_x","sample_y"), validation = TRUE)
st_write(temp, paste0(linkData, "/processed/2017/WURchangeHarmonics.gpkg"), "NDVI")

# 2018 #
# 68 observations / dates
b42018 = b4[,colnames(b4)[grepl("2017|2018|2019", colnames(b4))]]
b52018 = b5[,colnames(b5)[grepl("2017|2018|2019", colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b42018, as.numeric))
b5Temp = as.data.frame(sapply(b52018, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl("2017|2018|2019",dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))
HarmMetrics = cbind(b4[, c("location_id","sample_x", "sample_y")], HarmMetrics)

# change colnames
names(HarmMetrics)[4:(length(HarmMetrics))]= c("min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2018 features as GPKG
temp = DFtoSF(HarmMetrics, coords = c("sample_x","sample_y"), validation = TRUE)
st_write(temp, paste0(linkData, "/processed/2018/WURchangeHarmonics.gpkg"), "NDVI")
