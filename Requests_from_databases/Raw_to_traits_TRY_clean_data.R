# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# This code cleans the raw output from TRY into a trait x species matrix usable in the exploratories. #
# It is based on a script provided by J. Kattge & S. Tautenhahn, adapted and extended by M. Neyret 

# It goes through the following steps:
# 1. Check and format species names
# 2. Identify (and exclude) duplicates
# 3. Identify (and exclude) non-mature plants
# 4. Identify (and exclude) non-healthy plants
# 5. Identify (and exclude) non-standard explosition (i.e. shade or experiments)
# 6. Filter non-standard trait measurements
# 7. ? Identify (and exclude) measurements with non-standard units
# 8. Identify (and exclude) measurements from different climate zones
# 9. Check for outliers, incl. removing data sources with too many outliers
# 10. Then, for each trait, displays trait distribution with and without excluded data
# 11. Finally, get average value for each trait/species

# Note: For all steps excluding data, we make a list of all observations to exclude
# They will actually be excluded only at steps 10-11

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 0. Load libraries and data ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(readr)
library(data.table)
library(stringr)

setwd("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/")
TRYdata_original <- fread("Input_Data/11741.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)
setDT(TRYdata_original)

TRYdata <- TRYdata_original

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 1. Check and format species names ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# For now, all species should already be matched with exploratories species, so this part is not needed.
# For future updates it is also possible to rematch species to strict genus + species denomination by
# uncommenting the following lines.

# Species names should be only two words (genus + species)
TRYdata[, temp.word_length := str_count(AccSpeciesName, "\\w+")] # So we split by blank space, etc. to check

test <- TRYdata[temp.word_length > 2, unique(AccSpeciesName)] # Check for the different cases
TRYdata[, new.AccSpeciesName := AccSpeciesName] # Create new columns
TRYdata[temp.word_length > 2, new.AccSpeciesName := word(AccSpeciesName, 1, 2, sep = " ")] # Keep only 2 first words of the name

TRYdata <- TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 2. Identify (and exclude) duplicates ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# First we check that no ObsDataID is duplicated
length(TRYdata$ObsDataID) == length(unique(TRYdata$ObsDataID)) # -> OK

# Then there is apparently entries that were entered multiple times
obs_exclude_duplicate <- unique(TRYdata[DataName == "Duplicate data", ObservationID])

# So we remove all the entries with corresponding ObservationID
TRYdata = TRYdata[!(ObservationID %in% obs_exclude_duplicate),]



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 3. Identify  non-mature plants ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

TRYdata[DataID %in% c(413), unique(OrigValueStr)] # Check different levels and complete if needed below ("dummy")

TRYdata[, temp.juvenile := OrigValueStr == "dummy" |
  OrigValueStr == "immature" |
  OrigValueStr == "junvenil" |
  OrigValueStr == "juvenil (3 years)" |
  OrigValueStr == "juvenile" |
  OrigValueStr == "Juvenile" |
  OrigValueStr == "juvenile, 11-14 weeks" |
  OrigValueStr == "juvenile, 6 weeks" |
  OrigValueStr == "juveniles" |
  OrigValueStr == "sapling" |
  OrigValueStr == "Sapling (1 - 10 y)" |
  OrigValueStr == "saplings" |
  OrigValueStr == "seedling" |
  OrigValueStr == "Seedling" |
  OrigValueStr == "Seedling (0 - 1 y)" |
  OrigValueStr == "seedlings" |
  OrigValueStr == "seedlings, < 1/2 year" |
  OriglName == "Developmental stage" & OrigValueStr == "S" |
  OriglName == "Seedlings (True/False)" & OrigValueStr == "T" |
  OriglName == "Juvenile" & OrigValueStr == "Y"]

obs_exclude_juvenile <- unique(TRYdata[temp.juvenile == TRUE, ObservationID])

TRYdata <- TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns


# %%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 4. Identify non-healthy plants and organs ####
# %%%%%%%%%%%%%%%#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

TRYdata[DataID %in% c(3057, 1961), unique(OrigValueStr)] # Check different levels and complete below in 'dummy' if needed
TRYdata[DataID %in% c(3057, 1961), temp.unhealthy := OrigValueStr %in% c("Dead", "dummy")]

obs_exclude_dead <- unique(TRYdata[temp.unhealthy == TRUE, ObservationID])

TRYdata <- TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 5. Identify non-standard exposition ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# check how much data it would take out

TRYdata[, temp.exposition := OrigValueStr %in% c(
  "growth-chamber", "growth-chamber", "Climate chamber",
  "Climate Chamber", "Climate chamber, non-limiting conditions, (cf. dataset reference)",
  "climate chambers", "controlled environment room", "experimental", "experimental treatment", "GH", "Glass house", "Glasshouse", "Glasshouse experiment",
  "Greehouse", "Green house", "greenhouse", "Greenhouse", "Greenhouse plants", "Controlled climate chamber", "climate chamber", "open-top chamber",
  "open-sided growth chamber",
  "Greenhouse, grrowth container", "Greenhouse, Indiana University", "Greenhouse: highlight_highpH_competition", "Greenhouse: highlight_highpH_nocompetition",
  "Greenhouse: highlight_lowpH_competition", "Greenhouse: highlight_lowpH_nocompetition", "Greenhouse: lowleight_lowpH_competition", "Greenhouse: lowlight_highpH_competition", "Greenhouse: lowlight_highpH_nocompetition",
  "Greenhouse: lowlight_lowpH_nocompetition", "groth chamber", "growth chamber", "Growth chamber", "Growth Chamber", "Growth chamber, -N", "Growth chamber, +N", "growth chambers", "growth_chamber", "mesocosm", "Shade - Natural environment"
) |
  OriglName == "Artificial growth conditions (G=greenhouse, C=growth chamber)" & OrigValueStr == "G" |
  OriglName == "growingCondition" & OrigValueStr == "GH" |
  OriglName == "Natural / Greenhouse" & OrigValueStr == "G"]

TRYdata[DataID %in% c(327, 308, 210) & !temp.exposition, unique(OrigValueStr)] # Check different levels and complete above in 'dummy' if needed

obs_exclude_exposition <- unique(TRYdata[temp.exposition == TRUE, ObservationID])

TRYdata <- TRYdata[, .SD, .SDcols = colnames(TRYdata)[!grepl("temp.", colnames(TRYdata))]] # Removes temporary columns

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 6. Filter non-standard trait measurements ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Check the different modalities of DataName for each trait
tapply(
  TRYdata[TraitName != "" & DataID != "", ]$DataName,
  TRYdata[TraitName != "" & DataID != "", ]$TraitName,
  unique
)

# We check trait by trait and remove the Data which are not standard
remove_DataName <- c()

# Trait: `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included`
# DataName: "SLA: petiole  included"    "SLA: petiole included (2)" "SLA: petiole included (1)"
## -> all should be ok

# Trait:  `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded`
# DataName: "SLA: undefined if petiole in- or excluded"  "SLA: undefined if petiole in- or excluded (1)" "Leaf mass per area (LMA)"
## -> all should be ok

# Trait: `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`
# DataName: "Leaf dry matter content per leaf water-saturated mass (LDMC)"
#           "Leaf water content per leaf water-saturated mass (LWC, 1-LDMC)"  ** !!!! ???? Remove **
#           "LDMC min"       ** Remove **
#           "LDMC max"       ** Remove **
#           "Leaf dry fraction"
#           "Leaf dry matter concentration  (NIRS based)"
#           "Leaf dry matter content (min)"
#           "Leaf dry matter content (max)"
#           "Leaf lamina dry matter content (LDMC)"

remove_DataName <- c(remove_DataName, "Leaf water content per leaf water-saturated mass (LWC, 1-LDMC)", "LDMC min", "LDMC max")

# Trait: $`Leaf nitrogen (N) content per leaf dry mass`
# DataName: "Leaf nitrogen content per dry mass (Nmass)"
#           "N amount at 15 N measurement"
#           "Leaf senescent nitrogen (N) content per dry mass" ** Remove **
#           "Leaf nitrogen concentration  (NIRS based)"

remove_DataName <- c(remove_DataName, "Leaf senescent nitrogen (N) content per dry mass")

# Trait: `Leaf phosphorus (P) content per leaf dry mass`
# DataName: "Leaf phosphorus content per dry mass (Pmass)"
#           "Leaf senescent phosphorus (P) content per dry mass" ** Remove **

remove_DataName <- c(remove_DataName, "Leaf senescent phosphorus (P) content per dry mass")

# Trait: `Mycorrhiza type`
# DataName: "Mycorrhizal type"
#           "Original term for mycorrhizal type given by Selivanov"
#           "Mycorrhiza: Arbuscular mycorrhizal (AM) fungi"
#           "Mycorrhiza: Ectomycorrhizal (ECM) fungi"
#           "Mycorrhiza: No Mycorrhizal (NM) fungi"
#           "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi"
#           "Mycorrhiza: Orchid mycorrhizal (ORM) fungi"
#           "Mycorrhiza: Arbutoid mycorrhizal (ARB) fungi"
#           "Mycorrhiza: Presence of any non-AM mycorrhizal fungus (binary)"
#           "Mycorrhiza: Presence of any  alternative to AM mycorrhizal fungus (binary)"
#           "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)" ** Remove **
#           "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)"         ** Remove **
#           "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)"      ** Remove **
#           "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model"                                                   ** Remove **
#           "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model"                                               ** Remove **
#           "Mycorrhiza: Inferred binary AM state under our best HRM-model (i.e. \"Yes\" if AM_retained_likelihood > 0.5)"                           ** Remove **
#           "Mycorrhiza type according to Maherali"

remove_DataName <- c(
  remove_DataName, "Mycorrhiza: Likelihood the species was in the 'Stable AM Loss' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species was in the 'Labile' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species was in the 'Stable AM' state inferred under our best HRM-model (SI Figure 1, Extended Methods)",
  "Mycorrhiza: Likelihood the species had lost AM interactions under our best HRM-model",
  "Mycorrhiza: Likelihood the species had retained AM interactions under our best HRM-model",
  "Mycorrhiza: Inferred binary AM state under our best HRM-model (i.e. \"Yes\" if AM_retained_likelihood > 0.5)"
)

# Trait: Plant height vegetative`
#           "Plant height vegetative"                                           "Vegetative plant height, not elongated"
#           "Flowering plant height, heighest leaf elongated" ** Remove **       "Flowering plant height, heighest leaf not elongated"
#           "Average tree height"                                                "Height at 20 Years"
#           "Maximum plant height" ** Remove **                                  "Plant height observed"
#           "Plant height (unspecified if vegetative or reproductive)"           "Plant height of aquatic plants"
#           "Height of seedlings"                                                "Maximum height"
#           "Maximum height min"  ** Remove **                                   "Maximum height max"   ** Remove **
#           "Maximum height extreme"    ** Remove **                              "Plant height from flora"
#           "Plant height vegetative min."                             "          Plant height vegetative max."
#           "Vegetative plant height elongated"   ** Remove **

remove_DataName <- c(
  remove_DataName,
  "Flowering plant height, heighest leaf elongated",
  "Maximum plant height",
  "Maximum height min",
  "Maximum height extreme",
  "Maximum height max",
  "Vegetative plant height elongated"
)

# Trait: Seed dry mass`
#           "Seed dry mass"                               "Seed mass"
#           "Seed mass original value: mean"              "Seed mass min"  ** Remove **
#           "Seed mass max"    ** Remove **               "Seed mass original value: min" ** Remove **
#           "Seed mass original value: max"  ** Remove **

remove_DataName <- c(
  remove_DataName,
  "Seed mass min",
  "Seed mass max",
  "Seed mass original value: min",
  "Seed mass original value: max"
)
TRYdata[TraitName == "Seed dry mass", table(DataName)]


# Trait: Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)`
#           "Wood density; stem specific density; wood specific gravity"
# --> OK

# Trait: Mycorrhizal infection intensity
#           "Intensity of mycorrhizal infection"  
#           "Root length colonisation by vesicular-arbucular mycorrhizal fungi"

# --> OK

TRYdata[TraitName == 'Mycorrhizal infection intensity', table(DataName)]
remove_DataName <- c(
  remove_DataName,
  "Root length colonisation by vesicular-arbucular mycorrhizal fungi"
)


ObsDataID_exclude_nonstandard <- unique(TRYdata[DataName %in% remove_DataName, ObsDataID])



### I'm also renaming mycorrhiza into homogeneous names
TRYdata[
  TraitName == "Mycorrhiza type" & !(ObsDataID %in% ObsDataID_exclude_nonstandard),
  StdValue_myco := ifelse(OrigValueStr %in% c("ECTO", "ecto", "E.ch.ect.", "Ecto", "EC", "ectomycorrhizal", "Mycorrhiza: Ectomycorrhizal (ECM) fungi", "EM", "ec?", "Ectomycorrhiza"), "-Ectomycorrhiza",
    ifelse(OrigValueStr %in% c("Arbutoid", "AbtM", "Ecto arbut.", "E.t.ect.arb.", "Mycorrhiza: Arbutoid mycorrhizal (ARB) fungi"), "-Arbutoid mycorrhiza",
      ifelse(OrigValueStr %in% c("VAM", "AMNM", "NM/AM", "?va", "VA", "Mycorrhiza: Arbuscular mycorrhizal (AM) fungi", "AM", "arbuscular", "Gram-AM", "Forb-AM", "Ph.th.end."), "-Arbuscular mycorrhiza",
        # ifelse(OrigValueStr %in% c("AMNM", "NM/AM"), "-Weak arbuscular mycorrhiza",
        ifelse(OrigValueStr %in% c("AM + EM", "EC/AM"), "-Ecto and/or arbuscular mycorrhiza",
          ifelse(OrigValueStr %in% c("ericoid", "Mycorrhiza: Ericoid Mycorrhiza (ER) fungi", "Ericoid", "E.t.ect.er.", "ErM", "EM", "Ecto er.", "Endo er.", "ER", "ERICOID", "pyroloid"), "-Ericoid mycorrhiza",
            ifelse(OrigValueStr %in% c("Orchid", "OrM", "orchid", "E.t.end.", "Mycorrhiza: Orchid mycorrhizal (ORM) fungi"), "-Orchid mycorrhiza",
              ifelse(OrigValueStr %in% c("DS", "Ps.end.", "Ectendo", "Endo", "ecto, ectendo", "Endo-unidentified", "Ecto, endo, VA"), "-Other or undefined endomycorrhiza",
                ifelse(OrigValueStr %in% c("NM", 0, "Non", "no", "absent", "Mycorrhiza: No Mycorrhizal (NM) fungi", "N"), "-Non mycorrhizal",
                  NA
                )
              )
            )
          )
        )
        # )
      )
    )
  )
]



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 7. ? Identify measurements with non-standard units ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# tapply(
#  TRYdata[TraitName != "" & DataID != "", ]$UnitName,
#  TRYdata[TraitName != "" & DataID != "", ]$TraitName,
#  unique
# )

# -> We can use the StdValue column, which already standardize values, for all traits (!!! to check if new traits are added)
# This solves this issue


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 8. Identify measurements from different climate zones ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

## We currently decided to keep all datasets in the analysis, whatever their geographic origin is.
## You can change this by uncommenting the following lines

histogram(TRYdata[TRYdata$DataID == 59, StdValue, by = Dataset ]$StdValue, breaks = 100) # Latitude: bw 45 et 55 
histogram(TRYdata[TRYdata$DataID == 60, StdValue, by = Dataset ]$StdValue, breaks = 100) # Longitude: -5-30 (Europe) or 5-25 (Central Europe)
histogram(TRYdata[TRYdata$DataID == 61, StdValue, by = Dataset ]$StdValue, breaks = 100) # Altitude: < 1500

data_exclude_lat = TRYdata[TRYdata$DataID == 59 & (StdValue < 45 | StdValue> 55),]$ObservationID
data_exclude_lon = TRYdata[TRYdata$DataID == 60 &  (StdValue < 5 | StdValue > 25),]$ObservationID 
data_exclude_alt = TRYdata[TRYdata$DataID == 61 & (StdValue > 1500) ,]$ObservationID

obs_exclude_geography = unique(c(data_exclude_lat, data_exclude_lon, data_exclude_alt))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####  9. Average by contributor, then exclude outliers ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# For this section, we need to remove all previous identified sources of bias to effectively identify "outliers".
#TRYdata_all_data = copy(TRYdata)
TRYdata = copy(TRYdata_all_data)

TRYdata <- TRYdata[!(ObservationID %in%
  c(
   # obs_exclude_geography,
   #obs_exclude_exposition,
    obs_exclude_dead,
    obs_exclude_juvenile,
    obs_exclude_duplicate
  ) |
  ObsDataID %in% ObsDataID_exclude_nonstandard), ]

### We then average data for each contributor, to avoid too high weights from one dataset to bias the data

## 1. For quantitative traits -> take the mean of all estimates, removing min/max values when calculating the mean
# For most traits, the dominant value type is trait mean or single measurements, and we have enough data to kepp only those and remove the maximum/inimum values.
# However, this is not the case for Myco infection intensity, for which maximum values are provided in most cases, so we keep these values
TRYdata[TraitName == 'Mycorrhizal infection intensity', ValueKindName := tolower(ValueKindName)]

TRYdata_contributors <- TRYdata[TraitName != "Mycorrhiza type" & !(ValueKindName %in% c('Minimum', 'Maximum','Low', 'High')), 
                                list(StdValue = mean(StdValue, na.rm = T)), by = c(
  "new.AccSpeciesName", "LastName", "FirstName",
  "DatasetID", "Dataset", "SpeciesName", "AccSpeciesID", "AccSpeciesName",
  "TraitID", "TraitName", "UnitName", "Reference"
)]

# To calculate the range we do not exclude min/max
TRYdata_contributors_min_max <- TRYdata[TraitName != "Mycorrhiza type",  
                                list(Min = min(StdValue, na.rm = T),
                                     Max = max(StdValue, na.rm = T)),
                                by = c(
                                  "new.AccSpeciesName", "LastName", "FirstName",
                                  "DatasetID", "Dataset", "SpeciesName", "AccSpeciesID", "AccSpeciesName",
                                  "TraitID", "TraitName", "UnitName", "Reference"
                                )]
TRYdata_contributors_min_max[, Range := ifelse(Min != Max, paste(round(as.numeric(Min), 2), round(as.numeric(Max), 2), sep = '___'), round(as.numeric(Min), 2))][, Range := gsub("Inf___-Inf", NA, Range)]

TRYdata_contributors = merge(TRYdata_contributors, TRYdata_contributors_min_max)

## For myco traits --> choose the most common, if it's twice as common as the 2nd most common, otherwise describe as 'Mixed'
# Need to correct one binary database
TRYdata[Dataset == "Mycorrhizal Association Database" & TraitName == "Mycorrhiza type", OrigValueStr := ifelse(OrigValueStr == "Yes", DataName, NA)]
choose_myco <- function(data) {
  data = data[!is.na(data) & data != 'NA']
  x <- table(data)
  if (length(x) == 0) {
    value = NA
    prop = NA
  }
  else {
    # if ("-Weak arbuscular mycorrhiza" %in% names(x) & "-Arbuscular mycorrhiza" %in% names(x)) {
    #    x["-Arbuscular mycorrhiza"] <- x["-Arbuscular mycorrhiza"] + x["-Weak arbuscular mycorrhiza"]
    #    x <- x[names(x) != "-Weak arbuscular mycorrhiza"]
    #  }
    if (length(x) == 1) {
      value =names(x)
      prop = '100 %'
    }
    else {
      X <- sort(x, decreasing = TRUE)
      if (X[1] >= 2 * X[2]) {
        value = names(X[1])
        prop = paste(round(X[1]/sum(X), 2)*100, '%')
      }
      else {
        value = 'Mixed'
        prop = 'NA'
      }
    }
  }
  return(list(value = value, proportion = prop))
}

TRYdata_contributors_myco = TRYdata[TraitName == "Mycorrhiza type", list(StdValue = as.character(choose_myco(StdValue_myco)$value),
                                                                         Range = paste(unique(StdValue_myco), collapse = '___'),
                                                                         Min = NA,
                                                                         Max = NA), by = c(
  "new.AccSpeciesName", "LastName", "FirstName",
  "DatasetID", "Dataset", "SpeciesName", "AccSpeciesID", "AccSpeciesName",
  "TraitID", "TraitName", "UnitName", "Reference"
)]

TRYdata_contributors = rbind(TRYdata_contributors, TRYdata_contributors_myco)

## We also want to have an idea of how many measurements were made.
TRYdata[, Replicates2 := Replicates]
# If Replicates = '' but StdValue or StdValue_myco non NA, we consider that there is 1 replicate
TRYdata[!(is.na(StdValue) & is.na(StdValue_myco)) & Replicates == '', Replicates2 := 1]
# If ValueKindName is single value, we consider that it's one replicate.
TRYdata[ValueKindName == 'Single', Replicates2 := 1]
# If Replicate == Originat str value, it probably means that it's an error -> replace by 1.
TRYdata[Replicates == OrigValueStr, Replicates2 := 1]
# If StdValue = NA, we can't count it -> replace by 0
TRYdata[is.na(StdValue) & is.na(StdValue_myco), Replicates2 := 0]
# There is also a category "Number of replicates" in DataName column, which we need to re-merge into the main dataset to use it
Replicates_per_Observation = TRYdata[DataName == 'Number of replicates', list(Rep = unique(OrigValueStr)), by = c('ObservationID')]

TRYdata[TraitName != '' , Replicates2 := sapply(ObservationID, function(x){if (length(Replicates_per_Observation[ObservationID == x,]$Rep)>0){Replicates_per_Observation[ObservationID == x,]$Rep}
  else (NA)})]
TRYdata[, Replicates2 := as.numeric(Replicates2)]

TRYdata_contributors_n_measurements <- TRYdata[!(ValueKindName %in% c('Minimum', 'Maximum','Low', 'High')),  
                                               list(n_replicates = sum(Replicates2)),
                                               by = c(
                                                 "new.AccSpeciesName", "LastName", "FirstName",
                                                 "DatasetID", "Dataset", "SpeciesName", "AccSpeciesID", "AccSpeciesName",
                                                 "TraitID", "TraitName", "UnitName", "Reference"
                                               )]


### Merge into 1 dataset

TRYdata_contributors_all <- merge(TRYdata_contributors,
                              TRYdata_contributors_n_measurements
)

# identify outliers
is.outlier <- function(x) {
  print(x)
  outliers <- boxplot.stats(x)$out
  is.outlier <- x %in% outliers
  return(is.outlier)
}

TRYdata_contributors_all[TraitName != "Mycorrhiza type", is.outlier := is.outlier(as.numeric(StdValue)), by = c("TraitName", "new.AccSpeciesName")]
TRYdata_contributors_all[TraitName == "Mycorrhiza type", is.outlier := FALSE, by = c("TraitName", "new.AccSpeciesName")]

TRYdata_contributors_all = TRYdata_contributors_all[TraitName != '',]
TRYdata_contributors_all = TRYdata_contributors_all[!(TraitName != "Mycorrhiza type" & UnitName == ''),]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### 11. Finally, get average value for each trait/species ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

### Number of datasets per value
N_contrib = TRYdata_contributors_all[is.outlier == FALSE, 
                                     list(N_contrib = length(unique(Dataset)),
                                          References = paste(unique(Reference), collapse = '___')),
                                     by = c("TraitName", "new.AccSpeciesName",'AccSpeciesID')]

###  For quantitative traits ###
Final_data <- TRYdata_contributors_all[TraitName != "Mycorrhiza type" & !is.outlier, 
                                   list(StdValue = as.character(mean(as.numeric(StdValue), na.rm = T)),
                                        Std_std =  as.character(sd(as.numeric(StdValue), na.rm = T)),
                                        Reps = sum(n_replicates, na.rm = T),
                                        Min = min(Min, na.rm = T),
                                        Max = max(Max, na.rm = T)),
                                   by = c("TraitName", "new.AccSpeciesName", "UnitName", 'AccSpeciesID')]

Final_data[, Range := ifelse(Min != Max, paste(round(as.numeric(Min), 2), round(as.numeric(Max), 2), sep = '___'), round(as.numeric(Min), 2))][, Range := gsub("Inf___-Inf", NA, Range)]

### We, again, need to aggregate for mycorrhiza
Final_data_myco = TRYdata_contributors_all[TraitName == "Mycorrhiza type", 
                                       list(StdValue = as.character(choose_myco(StdValue)$value),
                                            Std_std = as.character(choose_myco(StdValue)$prop),
                                            Reps = sum(n_replicates),
                                            Min = NA,
                                            Max = NA,
                                            Range = sapply(Range, function(x){paste(unique(unlist(strsplit(x, '___'))), collapse = "___")})),
                                  by = c("TraitName", "new.AccSpeciesName", "UnitName",  'AccSpeciesID')
]

Final_data <- rbind(
  Final_data,
  Final_data_myco
)

Final_data = merge(Final_data, N_contrib)


Final_data[Final_data$Std_std == 0, Std_std := NA]


Final_data[TraitName == 'Mycorrhiza type' & StdValue == "Unknown", StdValue := NA]
Final_data = Final_data[TraitName !='',]

Final_data$Range = as.character(Final_data$Range)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####       12. Export dataset, raw and by species         ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# Fill in the final dataset with all species
#write.csv(Final_data, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits.csv") #6845
#write.csv(Final_data, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits_onlyCEurope.csv") #5314
#write.csv(Final_data, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Try_request/Final_data/TRY_plant_traits_with_shade.csv") #6892
nrow(Final_data)

