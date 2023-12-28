# Markov paper
# 29/11/2021
# Read in data

writeDates = function(b1Landsat){
  
  string = colnames(b1Landsat)[4:194]
  stringSub = substr(string, 2, 11)
  stringFinal = gsub("\\.", "-", stringSub)
  DFdates = data.frame(date=as.Date(stringFinal))
  write.csv(DFdates, paste0(linkData, "processed/dates.csv"), row.names=FALSE)
  
  return(TRUE)
}

extractDates = function(){
  
  DateCSV = read.csv(file="data/processed/dates.csv")
  dates = as.Date(DateCSV$date)
  
  return(dates)
}