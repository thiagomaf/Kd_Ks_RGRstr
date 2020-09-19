# LOAD LIBRARIES
[Very briefly indicate why each library is used and include explicity bibliographic citations]
```{r}
library(reshape2) # melt(); citation("reshape2")
library(dplyr)    # %>%;    citation("dplyr")
library(plyr)     # ddply;  citation("plyr") # must be loaded after dplyr

library(pander)   # pandoc.table
```

# DEFINE FUNCTIONS
## Standard error
The `se()` function is used to calculate the `standard error` of the mean. `x` represents a numeric vector containing replicates of a given measurement. The `na.rm = TRUE` argument indicates that by default this function ignores missing values.
```{r echo=TRUE}
# x:  numeric

se <- function(x, na.rm = TRUE) {
  sd(x, na.rm) / sqrt(length(x))
}
```

```{r echo=TRUE}
se(c(0.647603, 0.547048, 0.529873, NA, 0.908040, 0.835195))
```

## Gaussian error propagation
The following function is based on article published on January 22, 2015 by Lee Pang in the R bloggers website (https://www.r-bloggers.com/easy-error-propagation-in-r/). It allows convenient integration with libraries of the tidyverse and the dplyr grammar. It builds on the mutate() function to apply the chain rule on any given mathematical formula.

To illustrate its working let's assume one needs to propagate the error on the following calculation:
$$
Z = (X-Y)/(X+Y)^2
$$
![Z](images/eq_Z_12px.png?raw=true "Z")

For this case, two arguments must be input to the `mutate_with_error()` function. First, the `.data` argument receives a `data.frame` containing the mean `x` and `y` values as individual columns plus the standard error of these means as `dx` and `dy` columns. This nomenclature is important: all columns containing standard errors must be named with `d` appended to its respective mean values column. Then `f` receives a `function` object indicating the calculation to be done. _Vide_ below for more details on the structure of the datasets.

Inside the `mutate_with_error()` function we have the `exprs` object which is a list of commands to be executed. First, the `deparse(f[[3]])` command transform the right-hand side of the formula `f` into a `character` string. Then, a new `character` string is constructed containing the full right-hand side of the formula that will be used to calculated the propagated error (_vide_ below). Finally, the left-hand side of the new error propagation formula is created by appending the character `d` to the original formula left-hand side (i.e. `Z` becomes `dZ`). The `mutate_with_error()` apply these commands and return the results of the calculation defined by the formula and its associated propagated error appended to `.data` and returns the results.

$$
dZ = sqrt((dX*(1/(X + Y)^2 - (X - Y) * (2 * (X + Y))/((X + Y)^2)^2))^2+(dY*(-(1/(X + Y)^2 + (X - Y) * (2 * (X + Y))/((X + Y)^2)^2)))^2)
$$
![dZ](images/eq_dZ_12px.png?raw=true "dZ")

```{r}
# .data: data.frame
#     f: function
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
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  
  #.data %>% mutate_(.dots=exprs) #### -> mutate_() is deprecated, figure out how to fix.
  .data[names(exprs)] <- lapply(exprs , function(x) { eval(parse(text= x), envir = .data) })
  
  .data
}
```

```{r, results='asis'}
data.frame(X = c(0.647, 0.547, 0.529, 0.908, 0.835), Y = c(1.072, 0.905, 0.877, 1.503, 1.383)) %>%
  summarise(dX = se(X), dY = se(Y), X = mean(X), Y = mean(Y)) %>%
  mutate_with_error(Z ~ (X-Y)/(X+Y)^2) %>%
  select(X, dX, Y, dY, Z, dZ) %>%
  pandoc.table()
```


# LOAD DATA

* All missing values should be set as NA in your original table
* CSV files are a convenient way of loading data to R.

```{r}
#################################################################################
vars_prot <- c("Prot_Ala", "Prot_Ser")
vars_free <- c("Free_Ala", "Free_Ser")
vars_id   <- c("Genotype", "time")
#################################################################################
```

```{r}
# LOAD ALL NEEDED DATA
ishihara2017_data <- read.csv2(paste0(main.folder, "older/KDKSRGR_13C.enrichment[3].csv"))

# check data structure
str(ishihara2017_data)
```

* "Genotype", "Exp", "time" - are table keys, i.e. identify unique sample replicates (id.vars)
* "Prot_Ala", "Prot_Ser"    - measurements of 13C labelled alananine and serine enrichments in proteins
* "Free_Ala", "Free_Ser"    - measurements of 13C labelled alananine and serine enrichments in free amino acids
* "Glc"                     - measurements of 13C labelled glucose enrichment in cell walls

```{r}
# print head of data
head(ishihara2017_data)
```


# PREPARE DATA
We first must transform data tables in the _wide_ format to the _tidy_ format. Check https://garrettgman.github.io/tidying/ for details. In addition to that, missing values are removed.
We then proceed to summarize the biological replicates by calculating the mean and the standard error of the mean for each genotype in each time point.
```{r}
ishihara2017_melt <- melt(ishihara2017_data, id.vars = c("Genotype", "Exp", "time")) %>%
  ddply(.(Genotype, variable, time), function(.each) {
    data.frame(value = mean(.each$value, na.rm = TRUE), dvalue = se(.each$value, na.rm = TRUE))
    }) %>%
  filter(complete.cases(.)) # remove any soapy missing values; NA

head(ishihara2017_melt)
#################################################################################
```

# ESTIMATION OF PROTEIN SYNTHESIS RATES (KS)
We calculate the average enrichment of free labelled alanine and serine at ZT24 and then the KS - using Gaussian error propagation (REF).
Data table is cast back into wide format prior to calculations.

```{r}
# NOT SURE MELTING IS NEEDED
ishihara_KS <- ishihara2017_melt %>%
  subset(time == 24) %>%                                                      # select data only from ZT24
  subset(variable %in% c(vars_prot, vars_free)) %>%                           # select only amino acids data
  ddply(as.quoted(vars_id[which(vars_id != "time")]), function(.each) {
        # Transform to wide format
        .means  <- .each %>% 
          dcast(Genotype + time ~ variable, value.var =  "value", na.rm = TRUE) # only mean
        .errors <- .each %>%
          dcast(Genotype + time ~ variable, value.var = "dvalue", na.rm = TRUE) # only se
        
        # This crazy stuff renames all columns which are not keys ('vars') and 
        # append "d" in front of it. This is important for the error propagation
        # calculation ahead.
        colnames(.errors)[which(!colnames(.errors) %in% vars_id)] <- 
          paste0("d", colnames(.errors)[which(!colnames(.errors) %in% vars_id)])
        
        # join both tables
        join(.means, .errors, by = vars_id[which(vars_id != "time")])
    }) %>%
  mutate_with_error( SAt1t2 ~ (Free_Ala + Free_Ser) / 2 ) %>%                 # calculate the correction factor
  mutate_with_error( KS_Ala ~ (Prot_Ala / SAt1t2) ) %>%                       # calculate KS for Alanine
  mutate_with_error( KS_Ser ~ (Prot_Ser / SAt1t2) ) %>%                       # calculate KS for Serine
  select(-c(c(vars_prot, vars_free, "SAt1t2", "time"), paste0("d", c(vars_prot, vars_free, "SAt1t2"))))   # remove unwanted columns

ishihara_KS
```

# ESTIMATION OF PROTEIN DEGRADATION RATES (KD)
## Pulse
### Estimation of pulse RGRSTR
```{r}
ishihara_RGRp <- ishihara2017_melt %>%
  subset(time %in% c(0, 24)) %>%                                            # select data only from ZT24
  subset(variable %in% c(vars_prot, vars_free, "Glc")) %>%                           # select only amino acids data
  ddply(as.quoted(vars_id[which(vars_id != "time")]), function(.each) {
      # Transform to wide format
      .means  <- .each %>% 
        dcast(Genotype ~ variable + time, value.var =  "value", na.rm = TRUE) # only mean
      .errors <- .each %>%
        dcast(Genotype ~ variable + time, value.var = "dvalue", na.rm = TRUE) # only se
      
      # This crazy stuff renames all columns which are not keys ('vars') and 
      # append "d" in front of it. This is important for the error propagation
      # calculation ahead.
      colnames(.errors)[which(!colnames(.errors) %in% vars_id)] <- 
        paste0("d", colnames(.errors)[which(!colnames(.errors) %in% vars_id)])
      
      # join both tables
      join(.means, .errors, by = vars_id[which(vars_id != "time")])
  }) %>%
  mutate_with_error( RGRp ~ Glc_24 - Glc_0 ) %>%   # calculate KSloss for Serine
  select(c("Genotype", "RGRp", "dRGRp"))

ishihara_RGRp
```

### Estimation of pulse KD
```{r}
ishihara_KDp <- join(ishihara_KS, ishihara_RGRp, by = "Genotype") %>%
  mutate_with_error( KDp_Ala ~ KS_Ala - RGRp ) %>%
  mutate_with_error( KDp_Ser ~ KS_Ser - RGRp ) %>%
  select(c("Genotype", "KDp_Ala", "dKDp_Ala", "KDp_Ser", "dKDp_Ser"))

ishihara_KDp
```


## Chase
### Estimation of chase RGRSTR
```{r}
var_glc <- "Glc"
vars_id <- c("Genotype", "time")

ishihara_RGRc <- ishihara2017_melt %>%
  subset(time %in% c(24, 120)) %>%
  subset(variable == var_glc) %>%
  ddply(as.quoted(vars_id[which(vars_id != "time")]), function(.each) {
      # Transform to wide format
      .means  <- .each %>% 
        dcast(Genotype ~ variable + time, value.var =  "value", na.rm = TRUE) # only mean
      .errors <- .each %>%
        dcast(Genotype ~ variable + time, value.var = "dvalue", na.rm = TRUE) # only se
      
      # This crazy stuff renames all columns which are not keys ('vars') and 
      # append "d" in front of it. This is important for the error propagation
      # calculation ahead.
      colnames(.errors)[which(!colnames(.errors) %in% vars_id)] <- 
        paste0("d", colnames(.errors)[which(!colnames(.errors) %in% vars_id)])
      
      # join both tables
      join(.means, .errors, by = vars_id[which(vars_id != "time")])
  }) %>%
  #mutate_with_error( RGRc ~ (Glc_120 / Glc_24)^(1/4) - 1 ) %>%      # calculate RGR[STR]; t = 4 days
  mutate_with_error( RGRc ~ -((Glc_120 / Glc_24)^(1/4) - 1) ) %>%      # calculate RGR[STR]; t = 4 days
  # mutate_with_error( RGRc ~ 1 - (Glc_120 / Glc_24)^(1/4) ) %>%      # calculate RGR[STR]; t = 4 days
  select(c("Genotype", "RGRc", "dRGRc"))
  
ishihara_RGRc
```

### Estimation of chase KD
```{r}
# NOT SURE MELTING IS NEEDED
ishihara_KSloss <- ishihara2017_melt %>%
  subset(time %in% c(120, 24)) %>%                                            # select data only from ZT24
  subset(variable %in% c(vars_prot, vars_free)) %>%                           # select only amino acids data
  ddply(as.quoted(vars_id[which(vars_id != "time")]), function(.each) {
      # Transform to wide format
      .means  <- .each %>% 
        dcast(Genotype ~ variable + time, value.var =  "value", na.rm = TRUE) # only mean
      .errors <- .each %>%
        dcast(Genotype ~ variable + time, value.var = "dvalue", na.rm = TRUE) # only se
      
      # This crazy stuff renames all columns which are not keys ('vars') and 
      # append "d" in front of it. This is important for the error propagation
      # calculation ahead.
      colnames(.errors)[which(!colnames(.errors) %in% vars_id)] <- 
        paste0("d", colnames(.errors)[which(!colnames(.errors) %in% vars_id)])
      
      # join both tables
      join(.means, .errors, by = vars_id[which(vars_id != "time")])
  }) %>%
  # mutate_with_error( KSloss_Ala ~ (log(Prot_Ala_120/Free_Ala_120) - log(Prot_Ala_24/Free_Ala_24)) / (4) ) %>%    # calculate KSloss for Alanine
  # mutate_with_error( KSloss_Ser ~ (log(Prot_Ser_120/Free_Ser_120) - log(Prot_Ser_24/Free_Ser_24)) / (4) ) %>%    # calculate KSloss for Serine
  mutate_with_error( KSloss_Ala ~ (log(Prot_Ala_120) - log(Prot_Ala_24)) / (4) ) %>%                             # calculate KSloss for Alanine
  mutate_with_error( KSloss_Ser ~ (log(Prot_Ser_120) - log(Prot_Ser_24)) / (4) ) %>%                             # calculate KSloss for Serine
  select(c("Genotype", "KSloss_Ala", "dKSloss_Ala", "KSloss_Ser", "dKSloss_Ser"))                                # remove unwanted columns

ishihara_KSloss
```

```{r}
join(ishihara_KSloss, ishihara_RGRc, by = "Genotype") %>%
  mutate_with_error( KDc_Ala ~ -KSloss_Ala - RGRc ) %>%
  mutate_with_error( KDc_Ser ~ -KSloss_Ser - RGRc ) %>%
  #subset(Genotype == "Col0") %>%
  select(c("Genotype", "KDc_Ala", "dKDc_Ala", "KDc_Ser", "dKDc_Ser"))
```
