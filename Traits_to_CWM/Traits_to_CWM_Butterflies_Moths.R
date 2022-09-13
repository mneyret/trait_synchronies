# This script takes as input the abundances and species-level traits of moth and butterflies species found 
# in the Exploratories grasslands and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.

library(data.table)
library(reshape2)
library(vegan)
library(taxize)
library(FD)
library(readr)
library(textclean)
library(readxl)
library(dplyr)
library(traitdataform)
library(betapart)

setwd("~/Data")

figures_path = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/'
cwm_path = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data'
matched_traits_path = '~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Matched_trait_datasets/'

# ############### #
#### Load data ####
# ############### #

### Abundances datasets from the Exploratories
## Raw diversity
allsp <- fread("Abundances/210112_EP_species_diversity_GRL_BEXIS.txt") # https://www.bexis.uni-jena.de/ddm/data/Showdata/27706
allsp$Species = gsub('_$', '', allsp$Species ) # Remove _ if last character

## Species information
fgs <- fread("Abundances/210112_EP_species_info_GRL_BEXIS.txt") # https://www.bexis.uni-jena.de/ddm/data/Showdata/27707
fgs$Species = gsub(' $', '', fgs$Species ) # Correct species names
fgs$Species = gsub(' ', '_', fgs$Species ) 

Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE) # Merge abundance and information datasets
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))] # Correct plot ID format

Abundances_lepi = Abundance_all[Group_broad == 'Lepidoptera',][, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species] # Select only lepidoptera and homogenise taxonomy

### Traits from multipe databases
# Non-exploratories databases
butterflies_UKtraits = fread("Traits/Arthropods/5b5a13b6-2304-47e3-9c9d-35237d1232c6/data/ecological_traits_format.csv") #https://catalogue.ceh.ac.uk/documents/5b5a13b6-2304-47e3-9c9d-35237d1232c6
butterflies_EUMaghreb = data.table(read_excel("Traits/Arthropods/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", skip = 1)) #https://butterflytraits.github.io/European-Butterfly-Traits/index.html
butterflies_EUMaghreb_traits = data.table(read_excel("Traits/Arthropods/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", sheet = 2)) #https://butterflytraits.github.io/European-Butterfly-Traits/index.html

# From Exploratories
butterflies_Boerschig = data.table(read_excel("Traits/Arthropods/Butterflies_Boerschig2013.xlsx")) # Extracted from https://ars.els-cdn.com/content/image/1-s2.0-S1439179113001199-mmc1.pdf
colnames(butterflies_Boerschig)[2:11] = paste('B', colnames(butterflies_Boerschig)[2:11], sep = '_')

moths_traits_morpho_std = fread("Traits/Arthropods/23926_2_Dataset/23926_2_data.csv") # https://www.bexis.uni-jena.de/ddm/data/Showdata/23926
moths_traits_morpho = dcast.data.table(moths_traits_morpho_std, scientificName ~traitName, value.var = 'traitValue', fun.aggregate = mean)[, species := scientificName]
moths_traits_lifeH = fread("Traits/Arthropods/21228_2_Dataset/21228_2_data.csv") #https://www.bexis.uni-jena.de/ddm/data/Showdata/21228
moths_traits = merge.data.table(moths_traits_morpho, moths_traits_lifeH, by = 'species', all = T)

# Homogenise taxonomy
moths_traits[, scientific_name_gbif := get_gbif_taxonomy(species)$scientificName, by = species]
butterflies_UKtraits[, scientific_name_gbif := get_gbif_taxonomy(scientific_name)$scientificName, by = scientific_name][is.na(scientific_name_gbif) , scientific_name_gbif:= scientific_name]
butterflies_EUMaghreb[, scientific_name_gbif := get_gbif_taxonomy(`Taxa name`)$scientificName, by = `Taxa name`][is.na(scientific_name_gbif) , scientific_name_gbif:= `Taxa name`]
butterflies_EUMaghreb_traits[, scientific_name_gbif := get_gbif_taxonomy(Taxon)$scientificName, by = Taxon][is.na(scientific_name_gbif) , scientific_name_gbif:= Taxon]
butterflies_Boerschig[, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species][is.na(scientific_name_gbif) , scientific_name_gbif:= Species]

# ############################ #
#### Trait-by-trait merge ####
# ############################ #
# The different databases code each trait differently. We recode the traits in a semi-continuous scale common to all databases.

### Voltinism

# butterflies_UKtraits: 1, more or partial
# butterflies_EUMaghreb: number of generations per year, min and max
# moths_traits: semivoltine (i.e., all individuals must undergo two periods of hiber- nation to complete their development), (2) strictly univoltine, and (3) multivoltine

# --> code as 0.5 if semivoltine, 1 if strictly uni, 1.5 if 1 or more, 2 if more
Test_voltinism = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('obligate_univoltine', 'partial_generation', 'obligate_multivoltine', 'scientific_name_gbif'),],
                       butterflies_EUMaghreb[,.SD, .SDcols = c('Vol_biennial_0.5', 'Vol_univoltine_1', 'Vol_uni+partial2_1.5', 'Vol_bivoltine_2', 'Vol_multivoltine_3', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T),
                       moths_traits[,.SD, .SDcols = c('voltinism',  'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)
types = colnames(Test_voltinism)[2:10]
Test_voltinism[,  cases := paste(types[!is.na(.SD)], collapse = '-'), .SDcols = types, by = scientific_name_gbif] # Create all combinations of traits in the different databases to see if they match

# Semivoltine
Test_voltinism[voltinism == 1, voltinism_use := 0.5]
Test_voltinism[Vol_biennial_0.5 == 1, voltinism_use := 0.5]

# Strict univoltine
Test_voltinism[, temp := (voltinism == 2 | is.na(voltinism)) & (obligate_univoltine == 1 | is.na(obligate_univoltine)) & (Vol_univoltine_1 == 1 | is.na(Vol_univoltine_1)) &
                 all(is.na(c(Vol_biennial_0.5, partial_generation ,   obligate_multivoltine, `Vol_uni+partial2_1.5` , Vol_bivoltine_2  ,     Vol_multivoltine_3))) &
                 !all(is.na(c(voltinism, obligate_univoltine, Vol_univoltine_1))),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, voltinism_use := 1]

# Strict multivoltine
Test_voltinism[, temp := (voltinism == 3 | is.na(voltinism)) & (obligate_multivoltine == 1 | is.na(obligate_multivoltine)) &  (Vol_bivoltine_2 == 1 | Vol_multivoltine_3 == 1 | is.na(Vol_bivoltine_2) | is.na(Vol_multivoltine_3))&
                 all(is.na(c(obligate_univoltine,partial_generation, Vol_biennial_0.5,Vol_univoltine_1, `Vol_uni+partial2_1.5`))) &
                   !all(is.na(c(voltinism, obligate_multivoltine, Vol_bivoltine_2, Vol_multivoltine_3))),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, voltinism_use := 2]

# If the databases say it can be both univoltine and multivoltine, then we take 1.5
Test_voltinism[, temp := 
                         (obligate_univoltine == 1 & obligate_multivoltine == 1 | (is.na(obligate_univoltine) & is.na(obligate_multivoltine))) &
                            (Vol_univoltine_1 == 1 & (Vol_bivoltine_2 == 1 | Vol_multivoltine_3 == 1) | (is.na(Vol_univoltine_1) & is.na(Vol_bivoltine_2)))&
                 is.na(voltinism_use) & cases !='',
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, voltinism_use := 1.5]

# If it's mostly uni + half generation, then we give the score 1.25
Test_voltinism[, temp := cases %in% c(
  'obligate_univoltine-partial_generation','obligate_univoltine-Vol_univoltine_1-Vol_uni+partial2_1.5', 
  'Vol_univoltine_1-Vol_uni+partial2_1.5', 'obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_uni+partial2_1.5',
  "obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_bivoltine_2", 'bligate_univoltine-partial_generation-Vol_univoltine_1','obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2') |
    (cases == 'obligate_univoltine-partial_generation-voltinism' & voltinism == 2),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, voltinism_use := 1.25]


# For the rest, I keep data of Mangels et al in priority, then UE database, then UK database
Test_voltinism[is.na(voltinism_use) & voltinism == 3, voltinism_use := 2]

Test_voltinism[, temp := (cases %in% c(
  "obligate_multivoltine-Vol_univoltine_1-Vol_bivoltine_2-Vol_multivoltine_3",'obligate_multivoltine-Vol_univoltine_1-Vol_multivoltine_3',
  'bligate_multivoltine-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3','obligate_univoltine-obligate_multivoltine-Vol_bivoltine_2-Vol_multivoltine_3',
  'obligate_univoltine-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3','Vol_bivoltine_2-Vol_multivoltine_3-voltinism-NA-NA-NA-NA-NA-NA-NA',
  'obligate_univoltine-partial_generation-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_bivoltine_2-Vol_multivoltine_3','obligate_univoltine-Vol_multivoltine_3', 'partial_generation-obligate_multivoltine',
  'partial_generation-obligate_multivoltine-Vol_univoltine_1-Vol_bivoltine_2-Vol_multivoltine_3','Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'obligate_multivoltine-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'obligate_multivoltine-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'obligate_univoltine-obligate_multivoltine-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'partial_generation-obligate_multivoltine-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3',
  'bligate_multivoltine-Vol_biennial_0.5-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2-Vol_multivoltine_3')) | (cases %in% c('obligate_univoltine-partial_generation-voltinism', 'obligate_univoltine-Vol_multivoltine_3') & voltinism == 2),
  by = scientific_name_gbif]
Test_voltinism[temp == TRUE & is.na(voltinism_use), voltinism_use := 2]

Test_voltinism[, temp := cases %in% c(
  "obligate_univoltine-Vol_univoltine_1-Vol_bivoltine_2", 'Vol_uni+partial2_1.5-Vol_bivoltine_2'),
  by = scientific_name_gbif]
Test_voltinism[temp == TRUE & is.na(voltinism_use), voltinism_use := 1.5]


### Specificity 
# butterflies_UKtraits:  specificity (mono/oligo genus/ oligo family /poly)
# butterflies_EUMaghreb: specificity HPS_monophage	HPS_oligophag (within one genus)	HPS_oligophag (within one family )	HPS_polyphag
# moths_traits: (1) narrow specialists (host plants within one plant genus), (2) moderate specialists (host plants within one plant family), (3) moderate generalists (host plants recorded from two to four families), and (4) wide generalists (host plants in five or more families)

Test_generalism = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('specificity', 'scientific_name_gbif'),],
                             butterflies_EUMaghreb[,.SD, .SDcols = c('HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T),
                        moths_traits[,.SD, .SDcols = c('feeding_niche',  'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)

Test_generalism[, c('HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag') :=
                  lapply(.SD, as.numeric), .SDcols = c('HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag')]
Test_generalism[specificity == '', specificity := NA]

# One species
Test_generalism[, temp := (specificity == 'Monophagous' | is.na(specificity)) & (is.na(HPS_monophage) | HPS_monophage == 1) & (feeding_niche == 1 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, HPS_monophage)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, generalism_use := 1]

# One genus
Test_generalism[, temp := (specificity  == 'Oligophagous (Genus)' | is.na(specificity)) & (is.na(`HPS_oligophag (within one genus)`) | `HPS_oligophag (within one genus)` == 1) & (feeding_niche == 1 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, `HPS_oligophag (within one genus)`)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, generalism_use := 2]

# One family
Test_generalism[, temp := (specificity  == 'Oligophagous (Family)' | is.na(specificity)) & (is.na(`HPS_oligophag (within one family )`) | `HPS_oligophag (within one family )` == 1) & (feeding_niche == 2 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, `HPS_oligophag (within one family )`)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, generalism_use := 3]

# More than one family
Test_generalism[, temp := (specificity  == 'Polyphagous' | is.na(specificity)) & (is.na(HPS_polyphag) | HPS_polyphag == 1) & (feeding_niche %in% 3:4 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, HPS_polyphag)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, generalism_use := 4]


# Other cases
Test_generalism[ is.na(generalism_use) & !is.na(feeding_niche), generalism_use := ifelse(feeding_niche == 1, 2,
                                                                       ifelse(feeding_niche == 2, 3,
                                                                              4))]


Test_generalism[ is.na(generalism_use) & `HPS_oligophag (within one genus)` == 1,    generalism_use := 2]
Test_generalism[ is.na(generalism_use) & `HPS_oligophag (within one family )` == 1,    generalism_use := 3]
Test_generalism[ is.na(generalism_use) & `Polyphagous` == 1,    generalism_use := 4]


### Wing size / wingspan / weight 
# Mangels: size = wingspan, weight
# UK: estimated_dry_mass, forewing_minimum"  forewing_maximum
# UE: average wingspan, max forewing

Test_wing_weight = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('forewing_maximum', 'scientific_name_gbif', 'estimated_dry_mass'),],
                               butterflies_EUMaghreb[, list(FoL = mean(as.numeric(c(FoL_var_male_average, FoL_HR_average)), na.rm = T)), by = 'scientific_name_gbif'], by = 'scientific_name_gbif', all = T),
                               moths_traits[,.SD, .SDcols = c('size',  'wing_length', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)

# We impute wing size data as we can assume all dimensions to be very well correlated
miced_size = data.table(mice::complete(mice(Test_wing_weight)))

# Check that imputed size distribution
check_distribution = rbind(Test_wing_weight[, list(forewing_maximum, scientific_name_gbif, type = 'raw')],
                           miced_size[, list(forewing_maximum, scientific_name_gbif, type = 'imputed')])
ggplot(check_distribution, aes(forewing_maximum, fill = type)) + geom_histogram() +
  facet_wrap(~type, ncol = 1)

### Overwintering stage 
# Mangels: hibernation: 1 egg 2 larvae 3 pupa 4 adult
# UK: overwintering_egg" "overwintering_larva" "overwintering_pupa" "overwintering_adult
# UE: "OvS_egg_E" "OvS_larva_L" "OvS_pupa_P" "OvS_adult_A"

Test_overwintering = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c("overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult", 'scientific_name_gbif'),],
                               butterflies_EUMaghreb[, .SD, .SDcols = c("OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A", 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T),
                         moths_traits[,.SD, .SDcols = c('hibernation', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)

Test_overwintering[, c("overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult","OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A","hibernation") :=
                     lapply(.SD, function(x){
                       x = gsub('[c?]', '', x)
                       return(ifelse(is.na(x) | x == "" | x == 0 | x == 'NA', NA, x))}), .SDcols = c("overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult","OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A","hibernation"),
                   by = 'scientific_name_gbif']


# Egg
Test_overwintering[, temp := (overwintering_egg == 1 | is.na(overwintering_egg)) & 
                            (OvS_egg_E == 1  | is.na(OvS_egg_E)) & 
                           (hibernation == 1 | is.na(hibernation)) &
                    all(is.na(c(overwintering_larva, overwintering_pupa, overwintering_adult, OvS_larva_L, OvS_pupa_P, OvS_adult_A))|
                          c(overwintering_larva, overwintering_pupa, overwintering_adult, OvS_larva_L, OvS_pupa_P, OvS_adult_A) == 2) &
                  !(all(is.na(c(overwintering_egg, OvS_egg_E, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, wintering_stage := 1]

# Larvae
Test_overwintering[, temp := (overwintering_larva == 1 | is.na(overwintering_larva)) & 
                     (OvS_larva_L == 1  | is.na(OvS_larva_L)) & 
                     (hibernation == 2 | is.na(hibernation)) &
                     all(is.na(c(overwintering_pupa, overwintering_egg, overwintering_adult, OvS_pupa_P, OvS_egg_E, OvS_adult_A))|
                           c(overwintering_pupa, overwintering_egg, overwintering_adult, OvS_pupa_P, OvS_egg_E, OvS_adult_A) == 2)&
                     !(all(is.na(c(overwintering_larva, OvS_larva_L, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, wintering_stage := 2]

# Puppa
Test_overwintering[, temp := (overwintering_pupa == 1 | is.na(overwintering_pupa)) & 
                     (OvS_pupa_P == 1  | is.na(OvS_pupa_P)) & 
                     (hibernation == 3 | is.na(hibernation)) &
                     all(is.na(c(overwintering_larva, overwintering_egg, overwintering_adult, OvS_larva_L, OvS_egg_E, OvS_adult_A))|
                           c(overwintering_larva, overwintering_egg, overwintering_adult, OvS_larva_L, OvS_egg_E, OvS_adult_A) ==2) &
                     !(all(is.na(c(overwintering_pupa, OvS_pupa_P, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, wintering_stage := 3]

# Adult
Test_overwintering[, temp := (overwintering_adult == 1 | is.na(overwintering_adult)) & 
                     (OvS_adult_A == 1  | is.na(OvS_adult_A)) & 
                     (hibernation == 4 | is.na(hibernation)) &
                     all(is.na(c(overwintering_larva, overwintering_egg, overwintering_pupa, OvS_larva_L, OvS_egg_E, OvS_pupa_P))|
                           c(overwintering_larva, overwintering_egg, overwintering_pupa, OvS_larva_L, OvS_egg_E, OvS_pupa_P) == 2) &
                     !(all(is.na(c(overwintering_adult, OvS_adult_A, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, wintering_stage := 4]


# When datasets disagree, we use data from Mangels et al. which was tailored to Exploratories species
Test_overwintering[is.na(wintering_stage), wintering_stage := hibernation]

# for the rest take maximum stage
Test_overwintering[is.na(wintering_stage) & (overwintering_egg   == 1   | OvS_egg_E   == 1), wintering_stage := 1]
Test_overwintering[is.na(wintering_stage) & (overwintering_larva == 1 | OvS_larva_L == 1), wintering_stage := 2]
Test_overwintering[is.na(wintering_stage) & (overwintering_pupa  == 1  | OvS_pupa_P  == 1), wintering_stage := 3]
Test_overwintering[is.na(wintering_stage) & (overwintering_adult == 1 | OvS_adult_A == 1), wintering_stage := 4]


### Flight period 
Test_flight = merge(butterflies_UKtraits[,list(Flight_max = sum(.SD, na.rm = T)), .SDcols = c("ad_jan", "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec'), by = scientific_name_gbif],
                         butterflies_EUMaghreb_traits[, .SD, .SDcols = c("FMo_Max", 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)
Test_flight[Flight_max == 0, Flight_max := NA]
mice_flight = complete(mice(Test_flight))

ggplot(Test_flight , aes(Flight_max, FMo_Max)) + geom_point()

### Put all traits together together 
Merged_traits = Reduce(function(...) merge(..., by = 'scientific_name_gbif', all = TRUE), 
                       list(
                         data.table(mice_flight)[,.SD, .SDcols = c('scientific_name_gbif', 'Flight_max')],
                         Test_overwintering[,.SD, .SDcols = c('scientific_name_gbif','wintering_stage')],
                         data.table(miced_size)[,.SD, .SDcols = c('scientific_name_gbif','forewing_maximum')],
                         Test_generalism[,.SD, .SDcols = c('scientific_name_gbif','generalism_use')],
                         Test_voltinism[,.SD, .SDcols = c('scientific_name_gbif','voltinism_use')],
                         butterflies_Boerschig[, .SD, .SDcols = c('scientific_name_gbif', 'B_Egg_number')]
                       ))

### A few species are very abundant in a few plots but currently don't have data, leading to very low trait coverage.
# We fill them with additional data for some abundant species
# data from http://lepiforum.org
# and http://www.pyrgus.de
Lythria_purpuraria = data.table('scientific_name_gbif'= 'Lythria purpuraria', 'Flight_max' = 6, 'wintering_stage' = 3, 'forewing_maximum' = 15, 'generalism_use' = 1, 'voltinism_use' = 2, 'B_Egg_number' = NA)
Zygaena_carniolica = data.table('scientific_name_gbif'= 'Zygaena carniolica', 'Flight_max' = 2, 'wintering_stage' = 2, 'forewing_maximum' = 13, 'generalism_use' = 3, 'voltinism_use' = 1, 'B_Egg_number' = NA)

Merged_traits = rbind(Merged_traits, rbind(Lythria_purpuraria, Zygaena_carniolica))
Merged_traits[scientific_name_gbif == 'Aphantopus hyperantus', generalism_use := 3] #http://www.pyrgus.de/Aphantopus_hyperantus_en.html

Merged_traits[, Size.log := log(forewing_maximum)]
Merged_traits[, Eggs.log := log(B_Egg_number)]
Merged_traits[, Flight.log := log(Flight_max)]

# ############# #
#### Output ####
# ############# #

### Species-level PCA

# Check the existence of a species-level slow-fast axis
pca_but_sp = dudi.pca(Merged_traits[gsub(' ', '_', scientific_name_gbif) %in% Abundance_all$Species & complete.cases(Merged_traits), 
                                    list('Flight_period.log' = Flight.log, 
                                         'Wintering_stage' = wintering_stage, 
                                         Size.log, 
                                         'Generalism' = generalism_use, 
                                         'Volitnism' = voltinism_use)], scannf = FALSE, nf = 2)

gg_but = fviz_pca(pca_but_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                     col.ind = "steelblue",
                     fill.ind = "white",
                     col.var = "black")
ggsave(gg_but, file = paste(figures_path, 'species_pca_but.pdf', sep = ''), width = 6, height = 5)


### Trait data coverage
Merged_traits[!is.na(scientific_name_gbif) & scientific_name_gbif %in% Abundances_lepi$scientific_name_gbif,-1][, lapply(.SD, function(x){length(x[!is.na(x) & x != 'NA'])})]
butterfly_CC = check_coverage(Merged_traits[!is.na(scientific_name_gbif),], Abundances_lepi[!is.na(value) & value != 0,], c('wintering_stage', 'Flight.log','forewing_maximum', 'generalism_use', 'voltinism_use', 'B_Egg_number'),
                              'scientific_name_gbif', 'scientific_name_gbif')
butterfly_CC_B = check_coverage(butterflies_Boerschig[!is.na(scientific_name_gbif),], Abundances_butterflies[!is.na(value) & value != 0,], colnames(boerschig_traits)[2:11],
                                'Species', 'Species')
butterfly_CC[, -c(1,2)][, lapply(.SD, range)]
butterfly_CC_B[, -c(1,2)][, lapply(.SD, median)]

### Community-weighted means
# Weighted
butterfly_cwm = my_cwm(Merged_traits[!is.na(scientific_name_gbif),], Abundances_lepi[!is.na(value) & value != 0,], c('Flight.log', 'wintering_stage', 'Size.log', 'generalism_use', 'voltinism_use', 'Eggs.log'),
                       'scientific_name_gbif', 'scientific_name_gbif')

pca_butterflies = dudi.pca(butterfly_cwm[,1:5],scannf = FALSE, nf = 2)
fviz_pca_biplot(pca_butterflies)
fviz_pca_biplot(pca_butterflies, axes = c(1,3))

cor.test(pca_butterflies$l1$RS1, env_data_lui[Plot %in% butterfly_cwm$Plot,]$LUI)

write.csv(butterfly_cwm[ , list(
  "but_Generalism" = generalism_use,
  "but_Size" = Size.log,
  "but_GenYear" = voltinism_use,
  "but_hibernation" = wintering_stage,
  "but_flight" = Flight.log,
  'but_Eggs' = Eggs.log,
  "Plot" = Plot   ,
  "Year" = Year), ], paste(cwm_path, "CWM_Butterflies.csv", sep = ''))


# Unweighted: changes with species presence-absence, not abundance
Abundances_lepi_presence_absence = Abundances_lepi[, list(value = sum(value, na.rm = T)), by = list(Plot,Year, scientific_name_gbif)]
Abundances_lepi_presence_absence[value>1, value := 1]
butterfly_cwm_noweight = my_cwm(Merged_traits[!is.na(scientific_name_gbif),], Abundances_lepi_presence_absence[!is.na(value) & value != 0,], c('Flight.log', 'wintering_stage', 'Size.log', 'generalism_use', 'voltinism_use', 'Eggs.log'),
                       'scientific_name_gbif', 'scientific_name_gbif')

write.csv(butterfly_cwm_noweight[ , list(
  "but_Generalism" = generalism_use,
  "but_Size" = Size.log,
  "but_GenYear" = voltinism_use,
  "but_hibernation" = wintering_stage,
  "but_flight" = Flight.log,
  'but_Eggs' = Eggs.log,
  "Plot" = Plot   ,
  "Year" = Year), ], paste(cwm_path, "CWM_Butterflies_noweight.csv", sep = ''))



#######################################
# Check turnover accross LUI gradient #
#######################################

data_lui <- fread("Environment/LUI_input_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = dcast(Abundances_lepi[!is.na(scientific_name_gbif) & value>0,],  Plot~scientific_name_gbif, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)



#### Save matched trait dataset ####

# Short version, only custom codes
write.csv(Merged_traits, paste(matched_traits_path, 'Matched_lepidoptera.csv', sep = ''))

#### Complete version, with verbatim codes and names
#
#Merged_traits_melt = melt.data.table(Merged_traits, id.vars = 'scientific_name_gbif', variable.name = "traitName", value.name = 'traitValue')
#
## UK database
#merged_UK = merge.data.table(Merged_traits_melt, butterflies_UKtraits[, c('scientific_name','obligate_univoltine', 'partial_generation', 'obligate_multivoltine', 
#                                                                          'specificity',
#                                                                          'forewing_maximum', 'estimated_dry_mass', 
#                                                                          "overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult", 
#                                                                          'ad_jan', "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec','scientific_name_gbif')] , by = 'scientific_name_gbif')
#
#merged_UK[traitName == 'voltinism_use', 'Generation_time' :=    paste(c('obligate_univoltine', 'partial_generation', 'obligate_multivoltine')[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c('obligate_univoltine', 'partial_generation', 'obligate_multivoltine'), by = scientific_name_gbif] 
#merged_UK[traitName == 'Flight_max', 'Flight_month_adult' :=       paste(c('ad_jan', "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec')[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c('ad_jan', "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec'), by = scientific_name_gbif] 
#merged_UK[traitName == 'wintering_stage', 'Wintering_stage' :=  paste(c("overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult")[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c("overwintering_egg", "overwintering_larva" ,"overwintering_pupa" ,"overwintering_adult"), by = scientific_name_gbif] 
#merged_UK[traitName == 'forewing_maximum', 'Forewing_maximum' := forewing_maximum] 
#merged_UK[traitName == 'forewing_maximum', 'Estimated_dry_mass' := estimated_dry_mass] 
#merged_UK[traitName == 'generalism_use', 'Specificity' :=  specificity, by = scientific_name_gbif] 
#
#merged_UK_melt = melt.data.table(merged_UK[, list(traitName, Generation_time, Flight_month_adult, Wintering_stage, Specificity, Forewing_maximum, Estimated_dry_mass, verbatimScientificName = scientific_name, scientific_name_gbif)], id.vars = c('traitName', 'scientific_name_gbif','verbatimScientificName'), variable.name = 'verbatimTraitname', value.name = 'verbatimTraitValue')
#merged_UK_melt = unique(merged_UK_melt[!is.na(verbatimTraitValue),])
#merged_UK_melt[, Source := 'Cook et al. 2021, UK butterflies trait database']
#
## EU database
#merged_EU = merge.data.table(Merged_traits_melt, butterflies_EUMaghreb[, c('Taxa name','Vol_biennial_0.5', 'Vol_univoltine_1', 'Vol_uni+partial2_1.5', 'Vol_bivoltine_2', 'Vol_multivoltine_3',
#                                                                           'HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag',
#                                                                           'FoL_var_male_average', 'FoL_HR_average',
#                                                                           "OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A", 
#                                                                           'scientific_name_gbif')] , by = 'scientific_name_gbif')
#merged_EU = merge.data.table(merged_EU, butterflies_EUMaghreb_traits[, list(`Taxa name` = Taxon, FMo_Max, scientific_name_gbif)], by = c('Taxa name', 'scientific_name_gbif'))
#
#merged_EU[traitName == 'voltinism_use', 'Generation_time' :=    paste(c('Vol_biennial_0.5', 'Vol_univoltine_1', 'Vol_uni+partial2_1.5', 'Vol_bivoltine_2', 'Vol_multivoltine_3')[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c('Vol_biennial_0.5', 'Vol_univoltine_1', 'Vol_uni+partial2_1.5', 'Vol_bivoltine_2', 'Vol_multivoltine_3'), by = scientific_name_gbif] 
#merged_EU[traitName == 'Flight_max', 'fMo_Max' := FMo_Max] 
#merged_EU[traitName == 'wintering_stage', 'Wintering_stage' :=  paste(c("OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A")[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c("OvS_egg_E", "OvS_larva_L" ,"OvS_pupa_P" ,"OvS_adult_A"), by = scientific_name_gbif] 
#merged_EU[traitName == 'forewing_maximum', 'foL_var_male_average' := FoL_var_male_average] 
#merged_EU[traitName == 'forewing_maximum', 'foL_HR_average' := FoL_HR_average] 
#merged_EU[traitName == 'generalism_use', 'Specificity':=    paste(c('HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag')[!is.na(.SD) & .SD != 'NA' & .SD != ''], collapse = '; '), .SDcols = c( 'HPS_monophage', 'HPS_oligophag (within one genus)', 'HPS_oligophag (within one family )', 'HPS_polyphag'), by = scientific_name_gbif] 
#
#merged_EU_melt = melt.data.table(merged_EU[, list(traitName,Generation_time, fMo_Max, Wintering_stage, foL_var_male_average, foL_HR_average, Specificity, verbatimScientificName = `Taxa name`, scientific_name_gbif)], id.vars = c('traitName', 'scientific_name_gbif','verbatimScientificName'), variable.name = 'verbatimTraitname', value.name = 'verbatimTraitValue')
#merged_EU_melt = unique(merged_EU_melt[!is.na(verbatimTraitValue),])
#merged_EU_melt[, Source := 'Middleton-Welling et al. 2020, Europe and maghreb butterflies trait database']
#
#
## Moths traits from Mangels et al. (Exploratories)
#merged_moth_traits = merge.data.table(Merged_traits_melt[, list(scientific_name_gbif,  traitName, traitValue)], 
#                                      moths_traits[, c('species',
#                                                                           'voltinism',
#                                                                           'feeding_niche',
#                                                                           'size',  'wing_length',
#                                                                           "hibernation", 
#                                                                           'scientific_name_gbif')] , by = 'scientific_name_gbif')
#
#merged_moth_traits[traitName == 'voltinism_use',    'Voltinism' :=   voltinism] 
#merged_moth_traits[traitName == 'wintering_stage',  'Hibernation' :=  hibernation] 
#merged_moth_traits[traitName == 'forewing_maximum', 'Size' := size] 
#merged_moth_traits[traitName == 'forewing_maximum', 'Wing_length' := wing_length] 
#merged_moth_traits[traitName == 'generalism_use',   'Feeding_niche':= feeding_niche] 
#
#merged_moth_traits_melt = melt.data.table(merged_moth_traits[, list(traitName, Voltinism, Hibernation, Size, Wing_length, Feeding_niche, verbatimScientificName = `species`, scientific_name_gbif)], id.vars = c('traitName', 'scientific_name_gbif','verbatimScientificName'), variable.name = 'verbatimTraitname', value.name = 'verbatimTraitValue')
#merged_moth_traits_melt = unique(merged_moth_traits_melt[!is.na(verbatimTraitValue),])
#merged_moth_traits_melt[, Source := 'Mangels et al.']
#
#
## Put al datasources together
#Merged_traits_all = rbindlist(list(
#   merge.data.table(Merged_traits_melt, merged_UK_melt, id.vars = c('scientific_name_gbif', 'traitName')),
#   merge.data.table(Merged_traits_melt, merged_moth_traits_melt, id.vars = c('scientific_name_gbif', 'traitName')),
#   merge.data.table(Merged_traits_melt, merged_EU_melt, id.vars = c('scientific_name_gbif', 'traitName'), all.x = T)))
#
#
#Merged_traits_all = Merged_traits_all[!is.na(Merged_traits_all$traitValue)]
#Merged_traits_all[traitName == 'B_Egg_number', c('verbatimTraitname', 'Source') := list('Egg number', 'Boerschig et al. 2013')]
#Merged_traits_all = Merged_traits_all[!(traitName %in% c('Size.log', 'Flight.log', 'Eggs.log'))]
#Merged_traits_all[scientific_name_gbif == 'Lythria purpuraria', Source =  c('http://lepiforum.org, http://www.pyrgus.de')]
#Merged_traits_all[scientific_name_gbif == 'Zygaena carniolica', Source =  c('http://lepiforum.org, http://www.pyrgus.de')]
#Merged_traits_all[scientific_name_gbif == 'Aphantopus hyperantus', Source =  c('http://lepiforum.org, http://www.pyrgus.de')]

