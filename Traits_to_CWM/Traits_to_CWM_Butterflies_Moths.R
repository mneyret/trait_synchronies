# This script takes as input the abundances and species-level traits of moth and butterflies species
# and outputs a CWM matrix averaged for all considered years.
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

setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/")

# Abundances
## Raw diversity
allsp <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_diversity_GRL_BEXIS.txt")
allsp$Species = gsub('_$', '', allsp$Species ) # Remove _ if last character
## Species information
fgs <- fread("//Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_info_GRL_BEXIS.txt")
fgs$Species = gsub(' $', '', fgs$Species )
fgs$Species = gsub(' ', '_', fgs$Species )
Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE)
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))]

Abundances_lepi = Abundance_all[Group_broad == 'Lepidoptera',][, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species]

lepi_sp = unique(c(Abundance_all[Group_broad == 'Lepidoptera', get_gbif_taxonomy(Species)$scientificName, by = Species]$V1,
            moth_abundance_night[, get_gbif_taxonomy(Species)$scientificName, by = Species]$V1))

### Sorting out the multiple trait databases
# non-BE
butterflies_UKtraits = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/5b5a13b6-2304-47e3-9c9d-35237d1232c6/data/ecological_traits_format.csv")
butterflies_EUMaghreb = data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", skip = 1))
butterflies_EUMaghreb_traits = data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Arthropods/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", sheet = 2))

# From exploratories
butterflies_Boerschig = data.table(read_excel("Traits/Arthropods/Butterflies_Boerschig2013.xlsx"))
colnames(butterflies_Boerschig)[2:11] = paste('B', colnames(butterflies_Boerschig)[2:11], sep = '_')

moths_traits_morpho_std = fread("Traits/Arthropods/23926_2_Dataset/23926_2_data.csv")
moths_traits_morpho = dcast.data.table(moths_traits_morpho_std, scientificName ~traitName, value.var = 'traitValue', fun.aggregate = mean)[, species := scientificName]
moths_traits_lifeH = fread("Traits/Arthropods/21228_2_Dataset/21228_2_data.csv")
moths_traits = merge.data.table(moths_traits_morpho, moths_traits_lifeH, by = 'species', all = T)

# Homogenise taxonomy
moths_traits[, scientific_name_gbif := get_gbif_taxonomy(species)$scientificName, by = species]
butterflies_UKtraits[, scientific_name_gbif := get_gbif_taxonomy(scientific_name)$scientificName, by = scientific_name][is.na(scientific_name_gbif) , scientific_name_gbif:= scientific_name]
butterflies_EUMaghreb[, scientific_name_gbif := get_gbif_taxonomy(`Taxa name`)$scientificName, by = `Taxa name`][is.na(scientific_name_gbif) , scientific_name_gbif:= `Taxa name`]
butterflies_EUMaghreb_traits[, scientific_name_gbif := get_gbif_taxonomy(Taxon)$scientificName, by = Taxon][is.na(scientific_name_gbif) , scientific_name_gbif:= Taxon]
butterflies_Boerschig[, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species][is.na(scientific_name_gbif) , scientific_name_gbif:= Species]

##### Trait-by-trait merge ####
### Voltinism ####
# butterflies_UKtraits: 1, more or partial
# butterflies_EUMaghreb: number of generations per year, min and max
# moths_traits: semivoltine (i.e., all individuals must undergo two periods of hiber- nation to complete their development), (2) strictly univoltine, and (3) multivoltine

# --> code as 0.5 if semivoltine, 1 if strictly uni, 1.5 if 1 or more, 2 if more
Test_voltinism = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('obligate_univoltine', 'partial_generation', 'obligate_multivoltine', 'scientific_name_gbif'),],
                       butterflies_EUMaghreb[,.SD, .SDcols = c('Vol_biennial_0.5', 'Vol_univoltine_1', 'Vol_uni+partial2_1.5', 'Vol_bivoltine_2', 'Vol_multivoltine_3', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T),
                       moths_traits[,.SD, .SDcols = c('voltinism',  'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)
types = colnames(Test_voltinism)[2:10]
Test_voltinism[,  cases := paste(types[!is.na(.SD)], collapse = '-'), .SDcols = types, by = scientific_name_gbif]

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


# For the rest, FOR now, I keep data of Mangels et al > UE > UK
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

# Special cases
View(Test_voltinism[is.na(voltinism_use) & cases != '',])
Test_voltinism[is.na(voltinism_use), table(cases)]


### Specificity ####
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


### Wing size / wingspan / weight ####
# Mangels: size = wingspan, weight
# UK: estimated_dry_mass, forewing_minimum"  forewing_maximum
# UE: average wingspan, max forewing

Test_wing_weight = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('forewing_maximum', 'scientific_name_gbif', 'estimated_dry_mass'),],
                               butterflies_EUMaghreb[, list(FoL = mean(as.numeric(c(FoL_var_male_average, FoL_HR_average)), na.rm = T)), by = 'scientific_name_gbif'], by = 'scientific_name_gbif', all = T),
                               moths_traits[,.SD, .SDcols = c('size',  'wing_length', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)


miced_size = data.table(mice::complete(mice(Test_wing_weight)))

check_distribution = rbind(Test_wing_weight[, list(forewing_maximum, scientific_name_gbif, type = 'raw')],
                           miced_size[, list(forewing_maximum, scientific_name_gbif, type = 'imputed')])
ggplot(check_distribution, aes(forewing_maximum, fill = type)) + geom_histogram() +
  facet_wrap(~type, ncol = 1)

### Overwintering stage ####
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


# Consider that the BE dataset is right
Test_overwintering[is.na(wintering_stage), wintering_stage := hibernation]

# for the rest take maximum stage
Test_overwintering[is.na(wintering_stage) & (overwintering_egg   == 1   | OvS_egg_E   == 1), wintering_stage := 1]
Test_overwintering[is.na(wintering_stage) & (overwintering_larva == 1 | OvS_larva_L == 1), wintering_stage := 2]
Test_overwintering[is.na(wintering_stage) & (overwintering_pupa  == 1  | OvS_pupa_P  == 1), wintering_stage := 3]
Test_overwintering[is.na(wintering_stage) & (overwintering_adult == 1 | OvS_adult_A == 1), wintering_stage := 4]


### Flight period ####
Test_flight = merge(butterflies_UKtraits[,list(Flight_max = sum(.SD, na.rm = T)), .SDcols = c("ad_jan", "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec'), by = scientific_name_gbif],
                         butterflies_EUMaghreb_traits[, .SD, .SDcols = c("FMo_Max", 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)
Test_flight[Flight_max == 0, Flight_max := NA]
mice_flight = complete(mice(Test_flight))

ggplot(Test_flight , aes(Flight_max, FMo_Max)) + geom_point()

### Put everything together ####

Merged_traits = Reduce(function(...) merge(..., by = 'scientific_name_gbif', all = TRUE), 
                       list(
                         data.table(mice_flight)[,.SD, .SDcols = c('scientific_name_gbif', 'Flight_max')],
                         Test_overwintering[,.SD, .SDcols = c('scientific_name_gbif','wintering_stage')],
                         data.table(miced_size)[,.SD, .SDcols = c('scientific_name_gbif','forewing_maximum')],
                         Test_generalism[,.SD, .SDcols = c('scientific_name_gbif','generalism_use')],
                         Test_voltinism[,.SD, .SDcols = c('scientific_name_gbif','voltinism_use')],
                         butterflies_Boerschig[, .SD, .SDcols = c('scientific_name_gbif', 'B_Egg_number')]
                       ))

# Fill with additional data for some abundant species
# data from http://lepiforum.org
# http://www.pyrgus.de
Lythria_purpuraria = data.table('scientific_name_gbif'= 'Lythria purpuraria', 'Flight_max' = 6, 'wintering_stage' = 3, 'forewing_maximum' = 15, 'generalism_use' = 1, 'voltinism_use' = 2, 'B_Egg_number' = NA)
Zygaena_carniolica = data.table('scientific_name_gbif'= 'Zygaena carniolica', 'Flight_max' = 2, 'wintering_stage' = 2, 'forewing_maximum' = 13, 'generalism_use' = 3, 'voltinism_use' = 1, 'B_Egg_number' = NA)

Merged_traits = rbind(Merged_traits, rbind(Lythria_purpuraria, Zygaena_carniolica))
Merged_traits[scientific_name_gbif == 'Aphantopus hyperantus', generalism_use := 3] #http://www.pyrgus.de/Aphantopus_hyperantus_en.html

Merged_traits[scientific_name_gbif %in% Abundances_lepi[ Plot == 'AEG24' & value > 0, unique(scientific_name_gbif)],]
Merged_traits[, Size.log := log(forewing_maximum)]
Merged_traits[, Eggs.log := log(B_Egg_number)]
Merged_traits[, Flight.log := log(Flight_max)]

# Species-level PCA
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
ggsave(gg_but, file = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/species_pca_but.pdf', width = 6, height = 5)


# Coverage
Merged_traits[!is.na(scientific_name_gbif) & scientific_name_gbif %in% Abundances_lepi$scientific_name_gbif,-1][, lapply(.SD, function(x){length(x[!is.na(x) & x != 'NA'])})]
butterfly_CC = check_coverage(Merged_traits[!is.na(scientific_name_gbif),], Abundances_lepi[!is.na(value) & value != 0,], c('wintering_stage', 'Flight.log','forewing_maximum', 'generalism_use', 'voltinism_use', 'B_Egg_number'),
                              'scientific_name_gbif', 'scientific_name_gbif')
butterfly_CC_B = check_coverage(butterflies_Boerschig[!is.na(scientific_name_gbif),], Abundances_butterflies[!is.na(value) & value != 0,], colnames(boerschig_traits)[2:11],
                                'Species', 'Species')
butterfly_CC[, -c(1,2)][, lapply(.SD, range)]
butterfly_CC_B[, -c(1,2)][, lapply(.SD, median)]

# CWM
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
  "Year" = Year), ], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Butterflies.csv")


# CWM without weight
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
  "Year" = Year), ], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Butterflies_noweight.csv")


### Check turnover
data_lui <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Environment/LUI_input_data/LUI_standardized_global.txt")
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

library(betapart)
comm.test = dcast(Abundances_lepi[!is.na(scientific_name_gbif) & value>0,],  Plot~scientific_name_gbif, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)


