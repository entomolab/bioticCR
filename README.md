bioticCR
========

### Introduction

The bioticCR package calculates the Costa Rican version of the BMWP
freshwater invertebrate biotic index.

### Installation

In order to install the development version, you need to first install
the 'devtools' package.

    install.packages("devtools")

Then use the 'install\_github' function within this package to load the
development version from GitHub.

    devtools::install_github('robbriers/bioticCR')

### Input format

The data are assumed to be stored in a dataframe with the first column
containing a list of taxon names and subsequent columns containing
sample headers and individual taxon abundances. Absences should be
recorded as 'NA' rather than zero (so NA or blank cells in a csv file
imported via read.csv), although the calculations will work with zero
abundances as well, but this is not recommended. In addition to actual
abundances, log categories (1-5) or alphabetic log categories (A-E) can
also be processed, through specification of the data type (see below).

If you have data that are stored in a 'tidy' format [sensu
Wickham](https://www.jstatsoft.org/article/view/v059i10), then these can
be converted to the correct format for this package using the
[reshape2](https://cran.r-project.org/package=reshape2) package, or the
[tidyr](https://cran.r-project.org/package=tidyr) package prior to
calculation.

An example of the default layout of data can be seen by accessing the
built-in 'almond' dataset.

    # load library
    library(bioticCR)

    # show the format of the built-in almond dataset
    head(almond)

### Functions

There is one main function 'calcBMWP\_CR'. This calculates the value of
the Costa Rican version of BMWP for the sample data. The first argument
to the function is the object containing the data. A second (optional)
argument relates to the format of the dataset.The options are "num" for
numeric abundances, "log" for integer log categories (1-5) or "alpha"
for alphabetic log abundance categories (A-E). If the data are in the
default format (actual integer abundances) then it can be omitted. The
use of these is best illustrated through examples using the built-in
datasets ('almond', 'braidburn' and 'greenburn').

    # calculate the Costa Rican BMWP index for the River Almond dataset
    # 'type' does not have to specified as the default is used
    # ("num")

    calcBMWP_CR(almond)

To process data in either integer or alphabetic log categories, the
'type' argument should specify either "log" or "alpha". An example is
shown below.

    # example of processing data in alphabetic log abundance categories
    # using the 'type' argument

    # 'braidburn' dataset contains alphabetic log category data
    # see ?braidburn for details

    calcBMWP_CR(braidburn, "alpha")
