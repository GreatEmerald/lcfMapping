# Markov paper
# Split time series into yearly forecasts

# Access libraries and functions
library(sf)
library(pbapply)
library(ranger)
source("utils/extractDates.R")
source("utils/dataManagement.R")
source("utils/loadData.R")
source("utils/harmonicsFunctions.R")
source("RFfunction.R")

# Link to data folder (adapt for yourself)
linkData <- "data/"


## Apply on Validation set ##

# Yearly data consists of TS from year adjacent years + year itself
# Example: 2015 consists of data from 2014, 2015 and 2016
# In this script implemented with grepl function

b4 = st_read(paste0(linkData, "processed/WURvalidationFiltered.gpkg"), "b4")
b5 = st_read(paste0(linkData, "processed/WURvalidationFiltered.gpkg"), "b5")
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
dates = dates[grepl("2014|2015|2016",dates)] # important to run before getHarmonics below
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics, cl=parallel::detectCores()))

# adjust temporal harmonics
coordsData = read.csv(paste0(linkData, "processed/WURvalidationIDcoords.csv"))
HarmMetrics = cbind(coordsData, HarmMetrics)

# change colnames
names(HarmMetrics)[2:(length(HarmMetrics))]= c("x", "y", "min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2015 features as GPKG
# should create 2015 folder in data/processed folder
for (year in 2015:2018)
    if (!dir.exists(paste0("data/processed/", year)))
        dir.create(paste0("data/processed/", year))
temp = DFtoSF(HarmMetrics)
st_write(temp, "data/processed/2015/WURvalidationHarmonics.gpkg", "NDVI")
harmonics2015 = st_read("data/processed/2015/WURvalidationHarmonics.gpkg", "NDVI")
st_geometry(harmonics2015)=NULL


# load data with function (test)
dataTrain = loadTrainingData()
dataVali2015 = loadValidationData("2015")


# 2016 #
# 69 observations / dates
string2016 = "2015|2016|2017"
b4Temp = b4[,colnames(b4)[grepl(string2016, colnames(b4))]]
b5Temp = b5[,colnames(b5)[grepl(string2016, colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b4Temp, as.numeric))
b5Temp = as.data.frame(sapply(b5Temp, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl(string2016,dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))

# adjust temporal harmonics
coordsData = read.csv(paste0(linkData, "processed/WURvalidationIDcoords.csv"))
HarmMetrics = cbind(coordsData, HarmMetrics)

# change colnames
names(HarmMetrics)[2:(length(HarmMetrics))]= c("x", "y", "min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2016 features as GPKG
temp = DFtoSF(HarmMetrics)
st_write(temp, "data/processed/2016/WURvalidationHarmonics.gpkg", "NDVI")
harmonics2016 = st_read("data/processed/2016/WURvalidationHarmonics.gpkg", "NDVI")
st_geometry(harmonics2016)=NULL

# load data with function
dataTrain = loadTrainingData()
dataVali = loadValidationData("2016")

# 2017 #
# 69 observations / dates
string2017 = "2016|2017|2018"
b4Temp = b4[,colnames(b4)[grepl(string2017, colnames(b4))]]
b5Temp = b5[,colnames(b5)[grepl(string2017, colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b4Temp, as.numeric))
b5Temp = as.data.frame(sapply(b5Temp, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl(string2017,dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))

# adjust temporal harmonics
coordsData = read.csv(paste0(linkData, "processed/WURvalidationIDcoords.csv"))
HarmMetrics = cbind(coordsData, HarmMetrics)

# change colnames
names(HarmMetrics)[2:(length(HarmMetrics))]= c("x", "y", "min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2017 features as GPKG
temp = DFtoSF(HarmMetrics)
st_write(temp, "data/processed/2017/WURvalidationHarmonics.gpkg", "NDVI")
harmonics2017 = st_read("data/processed/2017/WURvalidationHarmonics.gpkg", "NDVI")
st_geometry(harmonics2017)=NULL

# load data with function
dataTrain = loadTrainingData()
dataVali = loadValidationData("2017")

# 2018 #
# 68 observations / dates
string2018 = "2017|2018|2019"
b4Temp = b4[,colnames(b4)[grepl(string2018, colnames(b4))]]
b5Temp = b5[,colnames(b5)[grepl(string2018, colnames(b5))]]

# calc ndvi
b4Temp = as.data.frame(sapply(b4Temp, as.numeric))
b5Temp = as.data.frame(sapply(b5Temp, as.numeric))
ndvi = (b5Temp - b4Temp) / (b5Temp + b4Temp)
rm(b4Temp)
rm(b5Temp)

# calc temporal harmonics
dates = extractDates()
dates = dates[grepl(string2018,dates)]
HarmMetrics = t(pbapply(as.matrix(ndvi), 1, getHarmonics))

# adjust temporal harmonics
coordsData = read.csv(paste0(linkData, "processed/WURvalidationIDcoords.csv"))
HarmMetrics = cbind(coordsData, HarmMetrics)

# change colnames
names(HarmMetrics)[2:(length(HarmMetrics))]= c("x", "y", "min", "max", "intercept", "co", 
                                               "si", "co2", "si2", "trend", "phase1", 
                                               "amplitude1", "phase2", "amplitude2")
names(HarmMetrics)

# Save (ndvi) 2018 features as GPKG
temp = DFtoSF(HarmMetrics)
st_write(temp, "data/processed/2018/WURvalidationHarmonics.gpkg", "NDVI")
#harmonics2018 = st_read("data/processed/2018/WURvalidationHarmonics.gpkg", "NDVI")
#st_geometry(harmonics2018)=NULL
