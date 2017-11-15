#' Calculate the Costa Rican BMWP biotic index for freshwater invertebrate samples
#' @description Calculates the Costa Rican version of the BMWP freshwater invertebrate biotic index.
#'
#' @param df A dataframe containing list of taxon names and their abundances
#' in samples, along with sample identifiers.  Default format is for taxon
#' names to be in the first column and sample abundances in subsequent
#' columns with identifers as column headers. See built-in \code{\link{almond}}
#'  dataset for an example. If data are in the transposed format i.e taxa as
#'  columns and samples as rows, the \code{\link{transposedata}} function can
#'  be used prior to calculation.
#'
#' @param type Indicates type of data being processed. Options are "num" for
#' numeric data, "log" for integer log abundance categories (1-5) or "alpha"
#' for alphabetic abundance categories (A-E). Default value is "num".
#'
#' @return A data frame consisting of columns of index values with samples in
#' rows.
#' @export calcBMWP_CR
#' @examples
#'
#' # calculate the Costa Rican version of the BMWP index for
#' # the built in 'almond' dataset
#' # 'type' does not have to specified as default is used
#' # ("num")
#'
#' calcBMWP_CR(almond)
#'
#' # example of processing data in alphabetic log abundance categories
#' # using the 'type' argument
#'
#' # 'braidburn' dataset contains alphabetic log category data
#' # see ?braidburn for details
#'
#' # calculate the Costa Rican BMWP index for this dataset
#'
#' calcBMWP_CR(braidburn, "alpha")
#'
#' # example of processing data in numeric log abundance categories
#' # using the 'type' argument
#'
#' # 'greenburn' dataset contains numeric log category data
#' # see ?greenburn for details
#'
#' # calculate the Costa Rican BMWP index for this dataset
#'
#' calcBMWP_CR(greenburn, "log")
#'
calcBMWP_CR<-function(df, type="num"){

  # explicitly specify index here
#  index<-"BMWP_CR"

  # check that a correct type has been specified
  TYPES<-c("num", "log", "alpha")
  datatype<-pmatch(type, TYPES)
  if (is.na(datatype))
    stop("Invalid data type specified")


  # if abundances are recorded as alphabetic categories, convert them
  if (type=="alpha"){
    df<-convertalpha(df)
  }
  # if abundances are recorded as log categories, convert them
  if (type=="log"){
    # check maximum value is 5
    maxabund<-max(df[,2:ncol(df)], na.rm=TRUE)
    if (maxabund>5){
        stop("Maximum value is > 5; not log categories")
    }
    df<-convertlog(df)
  }

  # check for and combine oligochaete families to class level

  # set up vector of oligochaete taxa
  families<-c("Lumbricidae", "Lumbriculidae", "Enchytraeidae", "Haplotaxidae", "Naididae", "Tubificidae", "Oligochaeta")

  # create logical vector of rows with oligochaetes
  present<- df[,1] %in% families

  # extract non oligochaete rows
  rest<-df[!present,]

  # subset rows with worms present
  worms<-df[present,]

  if (nrow(worms)!=0){
    # convert taxon to character for replacement
    worms[,1]<-as.character(worms[,1])

    # if there is more than one row of worms
    if(nrow(worms)>1){
      # if there is more than one sample
      if (ncol(worms)>2){
        # sum abundance across all oligochaetes and add to first row
        worms[1,-1]<-colSums(worms[,-1],na.rm=TRUE)
      } else {
        worms[1,-1]<-sum(worms[,2], na.rm=FALSE)
      }
    }
    # add taxon string back in
    worms[1,1]<-"Oligochaeta"

    # convert back to factor
    worms[,1]<-as.factor(worms[,1])

    # just take first row (sum)
    worms<-worms[1,]

    # recombine with rest of taxa
    df<-rbind(rest, worms)
  }

  # separate out sample taxon list
  taxonlist<-df[,1]

  # and samples
  samples<-df[,2:ncol(df)]

  # if only one sample present, need to process differently
  if (ncol(df)==2){

    # calculate scores
    output<-calcscore(samples, taxonlist=taxonlist)

    # transpose output
    output<-t(output)

    # add on sample name to row
    row.names(output)<-names(df[2])
  }
  # if there is more than one sample, apply across columns
  else{
    output<-apply(samples, 2, calcscore, taxonlist=taxonlist)
  }

  # only need to bind rows for multiple samples
  if (ncol(df)>2){
    output<-rbind(output)
    output<-t(output)
  }

  # add column names depending on index
  colnames(output)<-c("BMWP_CR")

  # add on sample identifier column and set row names to null
  output<-as.data.frame(cbind.data.frame(row.names(output), output))
  row.names(output)<-NULL
  colnames(output)[1]<-"Sample"

  return(output)
}
