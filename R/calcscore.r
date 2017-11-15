# Calculates biotic index scores for individual samples.
# @description Calculates chosen biotic index on individual samples,
# based on vectors of abundances and taxon names passed as arguments.
#
# @param abundances An integer vector representing the abundances of taxa
# in a sample. Absences should be represented as \code{NA} rather than
# zero values but it will process data with zeros.
# @param taxonlist A string vector containing taxon names (family level)
# present in all the samples being analysed.
# @param index A string representing the choice of index value. Defaults to
# BMWP.
# @return A list containing the different components of each score. List
# dimensions depend on the index calculated.

calcscore<-function(abundances, taxonlist){

  # recombine taxon list and abundances
  sampledata<-as.data.frame(cbind.data.frame(taxonlist, abundances))

  # find rows where taxa are present (>0 abundance)
  taxapresent<-sampledata[sampledata[,2] > 0,]

  # check that there are any taxa present in the sample, then extract scores
  if (nrow(taxapresent)!=0){
    samplescores<-extractrows(taxapresent, indextable)
    scorelist<-c(sum(samplescores$BMWP_CR, na.rm=TRUE))

  return(scorelist)
  # this closes the nrows if
  }
  # this is the function closure
}
