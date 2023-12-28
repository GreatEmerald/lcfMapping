# Markov paper
# Run the baseline RF model

source("utils/loadData.R")
source("RFfunctionNew.R")

linkData = "data/"
outDir = file.path(linkData, "output", "wurChange")

if (!dir.exists(outDir))
    dir.create(outDir, recursive=TRUE)

dataTrain = loadTrainingData()
# train2015 = loadChangeTrainingData("2015")
# train2016 = loadChangeTrainingData("2016")
# train2017 = loadChangeTrainingData("2017")
# train2018 = loadChangeTrainingData("2018")
# temp = subset(train2015, sample_id %in% train2016$sample_id)
# temp = subset(temp, sample_id %in% train2017$sample_id)
# temp = subset(temp, sample_id %in% train2018$sample_id)
# train2015 = subset(train2015, sample_id %in% temp$sample_id)
# train2016 = subset(train2016, sample_id %in% temp$sample_id)
# train2017 = subset(train2017, sample_id %in% temp$sample_id)
# train2018 = subset(train2018, sample_id %in% temp$sample_id)
# listChange = list("2015" = train2015,
#                   "2016" = train2016,
#                   "2017" = train2017,
#                   "2018" = train2018)
change2015 = loadChangeValidationData(year = "2015")
change2016 = loadChangeValidationData(year = "2016")
change2017 = loadChangeValidationData(year = "2017")
change2018 = loadChangeValidationData(year = "2018")
temp = subset(change2015, location_id %in% change2016$location_id)
temp = subset(temp, location_id %in% change2017$location_id)
temp = subset(temp, location_id %in% change2018$location_id)
change2015 = subset(change2015, sample_id %in% temp$sample_id)
change2016 = subset(change2016, sample_id %in% temp$sample_id)
change2017 = subset(change2017, sample_id %in% temp$sample_id)
change2018 = subset(change2018, sample_id %in% temp$sample_id)
years = list("2015" = change2015,
             "2016" = change2016,
             "2017" = change2017,
             "2018" = change2018)


# Run baseline RF model
listDFs = runRandomForest(train=dataTrain, years=years,
                         features=loadFeaturesNames(), PredictType="quantiles")
# Save results
lapply(1:length(listDFs), function(i) {
    result = listDFs[[i]]
    #names(result) = loadClassNames()
    write.csv(result, paste0(outDir, "/predictions-", names(years)[[i]], "-median.csv"), row.names = F)
})

# Check statistics
for (i in 1:length(listDFs))
{
    Prediction = listDFs[[i]]
    names(Prediction) = loadClassNames()
    print(getStats(Prediction, years[[i]][loadClassNames()]))
}

## Run Recurrent code

listDFs = runRecurrentRF(train=dataTrain, years=years,
                         features=loadFeaturesNames(), PredictType="quantiles")
# Save results
lapply(1:length(listDFs), function(i) {
    result = listDFs[[i]]
    write.csv(result, paste0(outDir, "/predictions-", names(years)[[i]], "-median-recurrent.csv"), row.names = F)
})

rec2015 = listDFs[[1]]
rec2016 = listDFs[[2]]
rec2017 = listDFs[[3]]
rec2018 = listDFs[[4]]

write.csv(rec2015, paste0(linkData, "/output/wurChange/predictions-2015-median-recurrent.csv"), row.names = F)
write.csv(rec2016, paste0(linkData, "/output/wurChange/predictions-2016-median-recurrent.csv"), row.names = F)
write.csv(rec2017, paste0(linkData, "/output/wurChange/predictions-2017-median-recurrent.csv"), row.names = F)
write.csv(rec2018, paste0(linkData, "/output/wurChange/predictions-2018-median-recurrent.csv"), row.names = F)

