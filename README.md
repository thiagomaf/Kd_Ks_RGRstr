# SETUP
```{r}
#################################################################################
# DEFINE WORKING DIRECTORY
#################################################################################
where.am.I <- "casa"
#where.am.I <- "MPI"
# working directory
if (where.am.I == "casa") {
  main.folder   <- "D:/MPI-MP/Dados/Hiro/"
  path          <- "Paper2020 - methods/R"
} else if (where.am.I == "MPI") {
  main.folder   <- "//mpimp-golm/USER/HOMES/Moraes/Dados/Hiro/"
  path          <- "Paper2020 - methods/R"
}
# set working directory
setwd(paste0(main.folder, path))
library(reshape2) # melt(); citation("reshape2")
library(dplyr)    # %>%;    citation("dplyr")
library(plyr)     # ddply;  citation("plyr") # must be loaded after dplyr
se <- function(x, na.rm = FALSE) {
  if(na.rm == TRUE) {  x <- as.vector(na.exclude(x)) }
  
  sd(x) / sqrt(length(x))
}
mutate_with_error = function(.data, f) {
  # this is based on work done by someone else and published in:
  # https://www.r-bloggers.com/easy-error-propagation-in-r/
  
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
  
  .data %>%
    # the standard evaluation alternative of mutate()
    mutate_(.dots=exprs)
}
vars_prot <- c("Prot_Ala", "Prot_Ser")
vars_free <- c("Free_Ala", "Free_Ser")
vars_id   <- c("Genotype", "time")
#################################################################################
```

# LOAD DATA

* All missing values should be set as NA in your original table
* CSV files are a convenient way of loading data to R.

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
We then proceed to summarise the biological replicates by calculating the mean and the standard error of the mean for each genotype in each time point.
```{r}
ishihara2017_melt <- melt(ishihara2017_data, id.vars = c("Genotype", "Exp", "time")) %>%
  #filter(value != 0) %>%       #include or not?
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
