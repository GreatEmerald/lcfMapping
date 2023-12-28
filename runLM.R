# Markov paper
# Postprocess basic RF using a linear regression
library(pbapply)

source("utils/loadData.R")

linkData = "data/"
inDir = file.path(linkData, "output", "wurChange")
outDir = file.path(linkData, "output", "linear")

if (!dir.exists(outDir))
    dir.create(outDir)

# Load all RF files into one big data.frame
inFiles = paste("predictions", 2015:2018, "median.csv", sep="-")
RFlist = lapply(file.path(inDir, inFiles), read.csv)
# Add a year and ID column
RFlist = lapply(1:length(RFlist), function(i){
    result = RFlist[[i]]
    result$year = (2015:2018)[i]
    result$ID = 1:nrow(result)
    return(result)
})
# Collapse into a single DF
RFdata = do.call(rbind, RFlist)

# Test
#plot(grassland~year, RFdata[RFdata$ID == 1,])
#fitted(lm(grassland~year, RFdata[RFdata$ID == 1,]))
#RFdata[RFdata$ID == 1,]

RunLM = function(RFslice)
{
    outDF = RFslice[c("year", "ID")]
    for (class in loadClassNames())
    {
        # Predict linear trend in the four years
        result = fitted(lm(as.formula(paste0(class, "~year")), RFslice))
        # Clamp to valid range
        result[result > 100] = 100
        result[result < 0] = 0
        outDF[[class]] = result
    }
    return(outDF)
}

# Run on all points - will take a while
LMlist = pbtapply(RFdata, list(RFdata$ID), RunLM, cl=parallel::detectCores())
LMout = do.call(rbind, LMlist)
# Rescale results to 100%
LMout[loadClassNames()] = LMout[loadClassNames()] / rowSums(LMout[loadClassNames()], na.rm=T) * 100

# Unpack back into the same structure as input
LMoutlist = split(LMout, LMout$year)
# Save to disk
lapply(LMoutlist, function(LMout) {
    year = LMout$year[1]
    LMout$ID = NULL
    LMout$year = NULL
    write.csv(LMout, paste0(outDir, "/smooth", year, ".csv"), row.names = FALSE)
})
