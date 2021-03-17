# LOAD LIBRARIES

``` r
library(reshape2) # melt(); citation("reshape2")
library(dplyr)    # %>%;    citation("dplyr")
library(plyr)     # ddply;  citation("plyr") # must be loaded after dplyr
library(rlang)

library(pander)   # pandoc.table
```

# DEFINE FUNCTIONS

## Standard error

The `se()` function is used to calculate the `standard error` of the
mean. `x` represents a numeric vector containing replicates of a given
measurement. The `na.rm = TRUE` argument indicates that by default this
function ignores missing values.

``` r
# x:  numeric

se <- function(x, ...) {
  sd(x, ...) / sqrt(length(x))
}
```

``` r
se(c(0.647603, 0.547048, 0.529873, NA, 0.908040, 0.835195))
```

    ## [1] NA

``` r
se(c(0.647603, 0.547048, 0.529873, NA, 0.908040, 0.835195), na.rm = TRUE)
```

    ## [1] 0.06965193

## Gaussian error propagation

The following function is based on article published on January 22, 2015
by Lee Pang in the R bloggers website
(<https://www.r-bloggers.com/easy-error-propagation-in-r/>). It allows
convenient integration with libraries of the tidyverse and the dplyr
grammar. It builds on the mutate() function to apply the chain rule on
any given mathematical formula.

To illustrate its working let’s assume one needs to propagate the error
on the following *Z* calculation:

The Gaussian propagated error *d**Z* can be calculated by applying the
chain-rule to *Z*, as follows:

The `mutate_with_error()` receives any given formula, construct the
respective formula for error propagation and returns the results of both
calculations. For this two arguments must be input to the
`mutate_with_error()` function. First, the `.data` argument receives a
`data.frame` containing the mean values to be used in the calculation
(i.e. *X* and *Y* in the example above) plus the standard error of these
means (i.e. *d**X* and *d**Y*) as individual columns. This nomenclature
is important: all columns containing standard errors must be named with
*d* appended to its respective mean values column. Then `f` receives a
`formula` object indicating the calculation to be done. *Vide* below for
more details on the structure of the datasets.

Inside the `mutate_with_error()` function we have the `exprs` object
which is a list of the two calculations to be done - the value of
interest and the propagated error. First, the `deparse(f[[3]])` command
transform the right-hand side of the formula `f` into a `character`
string. Then, a new `character` string is constructed containing the
full right-hand side of the formula that will be used to calculated the
propagated error (*vide* below). Finally, the left-hand side of the new
error propagation formula is created by appending the character *d* to
the original formula left-hand side (i.e. *Z* becomes *d**Z*). The
`mutate_with_error()` run these commands and return the results of the
calculation defined by the formula and its associated propagated error
appended to `.data`.

``` r
# .data: data.frame
#     f: formula
mutate_with_error = function(.data, f) {
  require(dplyr)
  
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    
    # expression to compute new variable errors
    sapply(all.vars(f[[3]]), function(v) {
      dfdp = deparse(D(f[[3]], v))
      sprintf('(d%s*(%s))^2', v, dfdp)
    }) %>%
      paste(collapse='+') %>%
      sprintf('sqrt(%s)', .)
  )
  
  names(exprs) = c(
    deparse(f[[2]]),                 # unchanged names
    sprintf('d%s', deparse(f[[2]]))  # 'd' appended to the names
  )
  
  # MATHEMAGICS!
  .data[names(exprs)] <- lapply(exprs , function(x) { eval(parse(text= x), envir = .data) })
  #.data %>% mutate_(.dots=exprs) #### -> mutate_() is deprecated, figure out how to fix.
  
  .data
}
```

``` r
data.frame(X = c(0.647, 0.547, 0.529, 0.908, 0.835), Y = c(1.072, 0.905, 0.877, 1.503, 1.383)) %>%
  summarise(dX = se(X), dY = se(Y), X = mean(X), Y = mean(Y)) %>%
  mutate_with_error(Z ~ (X-Y)/(X+Y)^2) %>%
  select(X, dX, Y, dY, Z, dZ) %>%
  pandoc.table()
```

|   X    |   dX    |   Y   |   dY   |    Z    |   dZ    |
|:------:|:-------:|:-----:|:------:|:-------:|:-------:|
| 0.6932 | 0.07639 | 1.148 | 0.1264 | -0.1342 | 0.03859 |

# LOAD DATA

The sections below demonstrates the procedure for calculating
*K*<sub>*s*</sub>, *K*<sub>*d*</sub> and *R**G**R*<sup>*S**T**R*</sup>
using data from Ishihara et al., 2017.

Datasets from the original publication are available online at \[LINK\].
Datasets are expected to contain replicate measurements as rows,
metadata information and measurement values as columns. Data from
Ishihara et al., 2017 is composed by replicate measurements from three
independent experiments (`Exp`) done in various genotypes (`Genotype`)
at different time points (`time`). It contains measurements of 13C
labelling alanine and serine enrichments in proteins (`Prot_Ala` and
`Prot_Ser`) and in free amino acids (`Free_Ala` and `Free_Ser`) plus
glucose enrichment in cell walls (`Glc`). \[NEED TO CHECK THE TEXT
ITSELF, PHRASE BETTER\]

-   All missing values should be set as `NA` in your original table
-   `CSV` files are a convenient way of loading data to R.

``` r
ishihara2017_data <- read.csv2("data/KDKSRGR_13C.enrichment[3].csv")

# check data structure
str(ishihara2017_data)
```

    ## 'data.frame':    117 obs. of  8 variables:
    ##  $ Genotype: chr  "Ang0" "Bsch2" "Bu2" "Col0" ...
    ##  $ Exp     : chr  "Exp1" "Exp1" "Exp1" "Exp1" ...
    ##  $ time    : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Prot_Ala: num  0.00827 0.00999 0 0.01054 0 ...
    ##  $ Prot_Ser: num  0.01137 0.01236 0.00932 0.01314 0 ...
    ##  $ Free_Ala: num  0.0126 0.0129 0.0119 0.0114 0.0164 ...
    ##  $ Free_Ser: num  0.00757 0.00637 0.00664 0.00782 0.01192 ...
    ##  $ Glc     : num  0.0361 0.0422 0.0343 0.0363 0.0344 ...

# PREPARE DATA

There are many different ways in which R can prepare datasets.
Nevertheless, the important bits are: the data must be summarised -
i.e. the mean and standard error of the mean should be used below which
must be placed in individual columns with names as indicated above.

``` r
vars_prot <- c("Prot_Ala", "Prot_Ser")
vars_free <- c("Free_Ala", "Free_Ser")
vars_glc  <- c("Glc")
vars_id   <- c("Genotype", "time")

ishihara2017_summary <- ishihara2017_data %>%
  #ddply(.(Genotype, time), function(each.data) {
  ddply(as.quoted(vars_id), function(each.data) {
    ldply(c(vars_prot, vars_free, vars_glc), function(each_var) {
      rbind(
        data.frame(variable = each_var,              value = mean(each.data[,each_var], na.rm = T)),
        data.frame(variable = paste0("d", each_var), value =   se(each.data[,each_var], na.rm = T))
      )
    })
  }) %>%
  #dcast(Genotype + time ~ variable)
  dcast(paste(paste(vars_id, collapse = " + "), "variable", sep = " ~ "))
```

``` r
ishihara2017_summary
```

# PULSE EXPERIMENTS

## Estimation of KS

We calculate the average enrichment of free labelled alanine and serine
at ZT24 and then the *K*<sub>*S*</sub> - using Gaussian error
propagation (REF). Data table is cast back into wide format prior to
calculations.

``` r
ishihara_KS <- ishihara2017_summary %>%
  subset(time == 24) %>%
  mutate_with_error( SAt1t2 ~ (Free_Ala + Free_Ser) / 2 ) %>%  # calculate the correction factor
  mutate_with_error( KS_Ala ~ (Prot_Ala / SAt1t2) ) %>%        # calculate KS for Alanine
  mutate_with_error( KS_Ser ~ (Prot_Ser / SAt1t2) ) %>%        # calculate KS for Serine
  select(
    -all_of(c( # remove unwanted columns
      "time",
      c(vars_prot, vars_free, vars_glc, "SAt1t2"),
      paste0("d", c(vars_prot, vars_free, vars_glc, "SAt1t2"))
    ))
  )

ishihara_KS
```

## Estimation of RGRSTR

Cases in which calculations involve comparisons between time points
demands reorganization of the data set we are using. Below we define a
logic to create a individual column for values of each time point. This
will be achieved by melting the original data set, renaming the variable
names (i.e. append the time to their names) and finally re-casting the
tidy data set into a wide data set; as follows:

``` r
ishihara_pulse_data <- ishihara2017_summary %>%
  subset(time %in% c(0, 24)) %>%
  melt(id.vars = vars_id) %>%
  mutate(variable = paste(variable, time, sep = "_")) %>%
  select(-time) %>%           # remove unwanted columns
  dcast(Genotype ~ variable)

ishihara_pulse_data
```

``` r
ishihara_RGRp <- ishihara_pulse_data %>%
  mutate_with_error( RGRp ~ Glc_24 - Glc_0 ) %>%
  select(c("Genotype", "RGRp", "dRGRp"))

ishihara_RGRp
```

## Estimation of KD

Calculation of *K*<sub>*d*<sub>*p*</sub></sub> can then be easily
achieved by joining the previously calculated datasets containing
*K*<sub>*s*</sub> and *R**G**R*<sub>*p*</sub><sup>*S**T**R*</sup>
results and apply the *K*<sub>*d*<sub>*p*</sub></sub> formula.

``` r
ishihara_KDp <- join(ishihara_KS, ishihara_RGRp, by = "Genotype") %>%
  mutate_with_error( KDp_Ala ~ KS_Ala - RGRp ) %>%
  mutate_with_error( KDp_Ser ~ KS_Ser - RGRp ) %>%
  select(c("Genotype", "KDp_Ala", "dKDp_Ala", "KDp_Ser", "dKDp_Ser"))

ishihara_KDp
```

# CHASE EXPERIMENTS

``` r
ishihara_chase_data <- ishihara2017_summary %>%
  subset(time %in% c(24, 120)) %>%
  melt(id.vars = vars_id) %>%
  mutate(variable = paste(variable, time, sep = "_")) %>%
  select(-time) %>%           # remove unwanted columns
  dcast(Genotype ~ variable)

ishihara_chase_data
```

## Estimation of KS

Note that results are in \#/day
:(120*h* − 24*h*) ÷ 24*h*/*d**a**y* = 4*d**a**y**s*; This value has to
be input as ‘numeric’ and not as an object given the behaviour of
`mutate_with_error()`. (This can be solved dynamically, to implement in
the future.)

``` r
# ishihara_KSloss <- ishihara_chase_data %>%
#   mutate_with_error( KSloss_Ala ~ (log(Prot_Ala_120) - log(Prot_Ala_24)) / (4) ) %>%
#   mutate_with_error( KSloss_Ser ~ (log(Prot_Ser_120) - log(Prot_Ser_24)) / (4) ) %>%
#   # remove unwanted columns
#   select(c("Genotype", "KSloss_Ala", "dKSloss_Ala", "KSloss_Ser", "dKSloss_Ser"))

ishihara_KSloss <- ishihara_chase_data %>%
  mutate_with_error( KSloss_Ala ~ -((log(Prot_Ala_120) - log(Prot_Ala_24))) / (4) ) %>%
  mutate_with_error( KSloss_Ser ~ -((log(Prot_Ser_120) - log(Prot_Ser_24))) / (4) ) %>%
  # remove unwanted columns
  select(c("Genotype", "KSloss_Ala", "dKSloss_Ala", "KSloss_Ser", "dKSloss_Ser"))

ishihara_KSloss
```

## Estimation of RGRSTR

Similarly to what was done for the pulse data, the
*R**G**R*<sub>*c*</sub><sup>*S**T**R*</sup> calculation also demands
re-organization of the original dataset.

Note that results are in \#/day
:(120*h* − 24*h*) ÷ 24*h*/*d**a**y* = 4*d**a**y**s*; This value has to
be input as ‘numeric’ and not as an object given the behaviour of
`mutate_with_error()`.

``` r
ishihara_RGRc <- ishihara_chase_data %>%
  mutate_with_error( RGRc ~ 1 - (Glc_120 / Glc_24)^(1/4) ) %>% 
  select(c("Genotype", "RGRc", "dRGRc"))

ishihara_RGRc
```

## Estimation of KD

``` r
# ishihara_KDc <- join(ishihara_KSloss, ishihara_RGRc, by = "Genotype") %>%
#   mutate_with_error( KDc_Ala ~ -KSloss_Ala - RGRc ) %>%
#   mutate_with_error( KDc_Ser ~ -KSloss_Ser - RGRc ) %>%
#   select(c("Genotype", "KDc_Ala", "dKDc_Ala", "KDc_Ser", "dKDc_Ser"))

ishihara_KDc <- join(ishihara_KSloss, ishihara_RGRc, by = "Genotype") %>%
  mutate_with_error( KDc_Ala ~ KSloss_Ala - RGRc ) %>%
  mutate_with_error( KDc_Ser ~ KSloss_Ser - RGRc ) %>%
  select(c("Genotype", "KDc_Ala", "dKDc_Ala", "KDc_Ser", "dKDc_Ser"))

ishihara_KDc
```
