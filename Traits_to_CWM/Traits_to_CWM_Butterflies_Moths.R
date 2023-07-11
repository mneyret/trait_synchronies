# This script takes as input the abundances and species-level traits of moth and butterflies species found 
# in the Exploratories grasslands and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.


#!!! Gbif version: Pre-Nov 2022 (taxonomy matching has been affected by the following update)

# ******************* #
#### 1. Load data ####
# ****************** #
set.seed(101)
Abundances_lepi = Abundance_all[Group_broad == 'Lepidoptera',]
Abundances_lepi[, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species] # Select only lepidoptera and homogenise taxonomy
Abundances_lepi[is.na(scientific_name_gbif), scientific_name_gbif := gsub('_', ' ', Species)]

### Traits from multipe databases
# Non-exploratories databases
butterflies_UKtraits = fread("Data/Trait_data/UK_butterflies/data/ecological_traits_format.csv") # https://catalogue.ceh.ac.uk/documents/5b5a13b6-2304-47e3-9c9d-35237d1232c6
butterflies_EUMaghreb = data.table(read_excel("Data/Trait_data/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", skip = 1)) #https://butterflytraits.github.io/European-Butterfly-Traits/index.html
butterflies_EUMaghreb_traits = data.table(read_excel("Data/Trait_data/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx",  sheet = 2)) #https://butterflytraits.github.io/European-Butterfly-Traits/index.html
#= data.table(read_excel("Data/Trait_data/European_Maghreb_Butterfly_Trait_data_v1.2/European_Maghreb_Butterfly_Trait_data_v1.2.xlsx", sheet = 2)) #https://butterflytraits.github.io/European-Butterfly-Traits/index.html

# From Exploratories
butterflies_Boerschig = data.table(read_excel("Data/Trait_data/Butterflies_Boerschig2013.xlsx")) # Extracted from https://ars.els-cdn.com/content/image/1-s2.0-S1439179113001199-mmc1.pdf
colnames(butterflies_Boerschig)[2:11] = paste('B', colnames(butterflies_Boerschig)[2:11], sep = '_')

moths_traits_morpho_std = fread("Data/Trait_data/23926_2_Dataset/23926_2_data.csv") # https://www.bexis.uni-jena.de/ddm/data/Showdata/23926
moths_traits_morpho = dcast.data.table(moths_traits_morpho_std, scientificName ~traitName, value.var = 'traitValue', fun.aggregate = mean)[, species := scientificName]
moths_traits_lifeH = fread("Data/Trait_data/21228_2_Dataset/21228_2_data.csv") #https://www.bexis.uni-jena.de/ddm/data/Showdata/21228
moths_traits = merge.data.table(moths_traits_morpho, moths_traits_lifeH, by = 'species', all = T)

#Homogenise taxonomy
moths_traits[, scientific_name_gbif := get_gbif_taxonomy(species)$scientificName, by = species]
moths_traits[is.na(scientific_name_gbif),scientific_name_gbif := gsub('_', ' ', species)]


butterflies_UKtraits[, scientific_name_gbif := get_gbif_taxonomy(scientific_name)$scientificName, by = scientific_name][is.na(scientific_name_gbif) , scientific_name_gbif:= scientific_name]
butterflies_EUMaghreb[, scientific_name_gbif := get_gbif_taxonomy(`Taxa name`)$scientificName, by = `Taxa name`][is.na(scientific_name_gbif) , scientific_name_gbif:= `Taxa name`]
butterflies_EUMaghreb_traits[, scientific_name_gbif := get_gbif_taxonomy(Taxon)$scientificName, by = Taxon][is.na(scientific_name_gbif) , scientific_name_gbif:= Taxon]
butterflies_Boerschig[, scientific_name_gbif := get_gbif_taxonomy(Species)$scientificName, by = Species][is.na(scientific_name_gbif) , scientific_name_gbif:= Species]

# ************************************** #
#### 2. Merge traits across databases ####
# ************************************** #

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
Test_voltinism[voltinism == 1, Voltinism_use := 0.5]
Test_voltinism[Vol_biennial_0.5 == 1, Voltinism_use := 0.5]

# Strict univoltine
Test_voltinism[, temp := (voltinism == 2 | is.na(voltinism)) & (obligate_univoltine == 1 | is.na(obligate_univoltine)) & (Vol_univoltine_1 == 1 | is.na(Vol_univoltine_1)) &
                 all(is.na(c(Vol_biennial_0.5, partial_generation ,   obligate_multivoltine, `Vol_uni+partial2_1.5` , Vol_bivoltine_2  ,     Vol_multivoltine_3))) &
                 !all(is.na(c(voltinism, obligate_univoltine, Vol_univoltine_1))),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, Voltinism_use := 1]

# Strict multivoltine
Test_voltinism[, temp := (voltinism == 3 | is.na(voltinism)) & (obligate_multivoltine == 1 | is.na(obligate_multivoltine)) &  (Vol_bivoltine_2 == 1 | Vol_multivoltine_3 == 1 | is.na(Vol_bivoltine_2) | is.na(Vol_multivoltine_3))&
                 all(is.na(c(obligate_univoltine,partial_generation, Vol_biennial_0.5,Vol_univoltine_1, `Vol_uni+partial2_1.5`))) &
                   !all(is.na(c(voltinism, obligate_multivoltine, Vol_bivoltine_2, Vol_multivoltine_3))),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, Voltinism_use := 2]

# If the databases say it can be both univoltine and multivoltine, then we take 1.5
Test_voltinism[, temp := 
                         (obligate_univoltine == 1 & obligate_multivoltine == 1 | (is.na(obligate_univoltine) & is.na(obligate_multivoltine))) &
                            (Vol_univoltine_1 == 1 & (Vol_bivoltine_2 == 1 | Vol_multivoltine_3 == 1) | (is.na(Vol_univoltine_1) & is.na(Vol_bivoltine_2)))&
                 is.na(Voltinism_use) & cases !='',
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, Voltinism_use := 1.5]

# If it's mostly uni + half generation, then we give the score 1.25
Test_voltinism[, temp := cases %in% c(
  'obligate_univoltine-partial_generation','obligate_univoltine-Vol_univoltine_1-Vol_uni+partial2_1.5', 
  'Vol_univoltine_1-Vol_uni+partial2_1.5', 'obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_uni+partial2_1.5',
  "obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_bivoltine_2", 'bligate_univoltine-partial_generation-Vol_univoltine_1','obligate_univoltine-partial_generation-Vol_univoltine_1-Vol_uni+partial2_1.5-Vol_bivoltine_2') |
    (cases == 'obligate_univoltine-partial_generation-voltinism' & voltinism == 2),
               by = scientific_name_gbif]
Test_voltinism[temp == TRUE, Voltinism_use := 1.25]


# For the rest, I keep data of Mangels et al in priority, then UE database, then UK database
Test_voltinism[is.na(Voltinism_use) & voltinism == 3, Voltinism_use := 2]

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
Test_voltinism[temp == TRUE & is.na(Voltinism_use), Voltinism_use := 2]

Test_voltinism[, temp := cases %in% c(
  "obligate_univoltine-Vol_univoltine_1-Vol_bivoltine_2", 'Vol_uni+partial2_1.5-Vol_bivoltine_2'),
  by = scientific_name_gbif]
Test_voltinism[temp == TRUE & is.na(Voltinism_use), Voltinism_use := 1.5]


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
Test_generalism[temp == TRUE, Generalism_use := 1]

# One genus
Test_generalism[, temp := (specificity  == 'Oligophagous (Genus)' | is.na(specificity)) & (is.na(`HPS_oligophag (within one genus)`) | `HPS_oligophag (within one genus)` == 1) & (feeding_niche == 1 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, `HPS_oligophag (within one genus)`)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, Generalism_use := 2]

# One family
Test_generalism[, temp := (specificity  == 'Oligophagous (Family)' | is.na(specificity)) & (is.na(`HPS_oligophag (within one family )`) | `HPS_oligophag (within one family )` == 1) & (feeding_niche == 2 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, `HPS_oligophag (within one family )`)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, Generalism_use := 3]

# More than one family
Test_generalism[, temp := (specificity  == 'Polyphagous' | is.na(specificity)) & (is.na(HPS_polyphag) | HPS_polyphag == 1) & (feeding_niche %in% 3:4 | is.na(feeding_niche)) &
                  !(all(is.na(c(specificity, HPS_polyphag)))), by = scientific_name_gbif]
Test_generalism[temp == TRUE, Generalism_use := 4]


# Other cases
Test_generalism[ is.na(Generalism_use) & !is.na(feeding_niche), Generalism_use := ifelse(feeding_niche == 1, 2,
                                                                       ifelse(feeding_niche == 2, 3,
                                                                              4))]

Test_generalism[ is.na(Generalism_use) & `HPS_oligophag (within one genus)` == 1,    Generalism_use := 2]
Test_generalism[ is.na(Generalism_use) & `HPS_oligophag (within one family )` == 1,    Generalism_use := 3]
Test_generalism[ is.na(Generalism_use) & specificity == "Polyphagous",    Generalism_use := 4]


### Wing size / wingspan / weight 
# Mangels: size = wingspan, weight
# UK: estimated_dry_mass, forewing_minimum"  forewing_maximum
# UE: average wingspan, max forewing

Test_wing_weight = merge(merge(butterflies_UKtraits[,.SD, .SDcols = c('forewing_maximum', 'scientific_name_gbif', 'estimated_dry_mass'),],
                               butterflies_EUMaghreb[, list(FoL = mean(as.numeric(c(FoL_var_male_average, FoL_HR_average)), na.rm = T)), by = 'scientific_name_gbif'], by = 'scientific_name_gbif', all = T),
                               moths_traits[,.SD, .SDcols = c('size',  'wing_length', 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)


Test_wing_weight[scientific_name_gbif %in% wrong_species,]
all_data[scientific_name_gbif %in% wrong_species & traitName == 'forewing_maximum',]
plot(Test_wing_weight$forewing_maximum, Test_wing_weight$wing_length )

# We impute wing size data as we can assume all dimensions to be very well correlated
#miced_size = Test_wing_weight#data.table(mice::complete(mice(Test_wing_weight)))
miced_size = data.table(mice::complete(mice(Test_wing_weight)))

# Check that imputed size distribution
check_distribution = merge.data.table(Test_wing_weight[, list(forewing_maximum, scientific_name_gbif, type = 'raw')],
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
Test_overwintering[temp == TRUE, Wintering_stage := 1]

# Larvae
Test_overwintering[, temp := (overwintering_larva == 1 | is.na(overwintering_larva)) & 
                     (OvS_larva_L == 1  | is.na(OvS_larva_L)) & 
                     (hibernation == 2 | is.na(hibernation)) &
                     all(is.na(c(overwintering_pupa, overwintering_egg, overwintering_adult, OvS_pupa_P, OvS_egg_E, OvS_adult_A))|
                           c(overwintering_pupa, overwintering_egg, overwintering_adult, OvS_pupa_P, OvS_egg_E, OvS_adult_A) == 2)&
                     !(all(is.na(c(overwintering_larva, OvS_larva_L, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, Wintering_stage := 2]

# Puppa
Test_overwintering[, temp := (overwintering_pupa == 1 | is.na(overwintering_pupa)) & 
                     (OvS_pupa_P == 1  | is.na(OvS_pupa_P)) & 
                     (hibernation == 3 | is.na(hibernation)) &
                     all(is.na(c(overwintering_larva, overwintering_egg, overwintering_adult, OvS_larva_L, OvS_egg_E, OvS_adult_A))|
                           c(overwintering_larva, overwintering_egg, overwintering_adult, OvS_larva_L, OvS_egg_E, OvS_adult_A) ==2) &
                     !(all(is.na(c(overwintering_pupa, OvS_pupa_P, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, Wintering_stage := 3]

# Adult
Test_overwintering[, temp := (overwintering_adult == 1 | is.na(overwintering_adult)) & 
                     (OvS_adult_A == 1  | is.na(OvS_adult_A)) & 
                     (hibernation == 4 | is.na(hibernation)) &
                     all(is.na(c(overwintering_larva, overwintering_egg, overwintering_pupa, OvS_larva_L, OvS_egg_E, OvS_pupa_P))|
                           c(overwintering_larva, overwintering_egg, overwintering_pupa, OvS_larva_L, OvS_egg_E, OvS_pupa_P) == 2) &
                     !(all(is.na(c(overwintering_adult, OvS_adult_A, hibernation)))), by = scientific_name_gbif]
Test_overwintering[temp == TRUE, Wintering_stage := 4]


# When datasets disagree, we use data from Mangels et al. which was tailored to Exploratories species
Test_overwintering[is.na(Wintering_stage), Wintering_stage := hibernation]

# for the rest take maximum stage
Test_overwintering[is.na(Wintering_stage) & (overwintering_egg   == 1   | OvS_egg_E   == 1), Wintering_stage := 1]
Test_overwintering[is.na(Wintering_stage) & (overwintering_larva == 1 | OvS_larva_L == 1), Wintering_stage := 2]
Test_overwintering[is.na(Wintering_stage) & (overwintering_pupa  == 1  | OvS_pupa_P  == 1), Wintering_stage := 3]
Test_overwintering[is.na(Wintering_stage) & (overwintering_adult == 1 | OvS_adult_A == 1), Wintering_stage := 4]


### Flight period 
Test_flight = merge(butterflies_UKtraits[,list(Flight_max = sum(.SD, na.rm = T)), .SDcols = c("ad_jan", "ad_feb", "ad_mar", "ad_apr","ad_may", 'ad_jun', 'ad_jul', 'ad_agu', 'ad_sep', 'ad_oct', 'ad_nov', 'ad_dec'), by = scientific_name_gbif],
                    butterflies_EUMaghreb_traits[, .SD, .SDcols = c("FMo_Max", 'scientific_name_gbif')], by = 'scientific_name_gbif', all = T)
Test_flight[Flight_max == 0, Flight_max := NA]
mice_flight = data.table(complete(mice(Test_flight)))

ggplot(Test_flight , aes(Flight_max, FMo_Max)) + geom_point()

### Put all traits together together 
Merged_traits = Reduce(function(...) merge(..., by = 'scientific_name_gbif', all = TRUE), 
                       list(
                         mice_flight[!is.na(scientific_name_gbif),.SD, .SDcols = c('scientific_name_gbif', 'Flight_max')],
                         Test_overwintering[!is.na(scientific_name_gbif),.SD, .SDcols = c('scientific_name_gbif','Wintering_stage')],
                                miced_size[!is.na(scientific_name_gbif),.SD, .SDcols = c('scientific_name_gbif','forewing_maximum')],
                                       Test_generalism[!is.na(scientific_name_gbif),.SD, .SDcols = c('scientific_name_gbif','Generalism_use')],
                                              Test_voltinism[!is.na(scientific_name_gbif),.SD, .SDcols = c('scientific_name_gbif','Voltinism_use')],
                                                    butterflies_Boerschig[!is.na(scientific_name_gbif), .SD, .SDcols = c('scientific_name_gbif', 'B_Egg_number')]
                       ))


# Check main species for which we don't have data
Abundances_lepi_agg = Abundances_lepi[Plot %in% plots & value >0 &!(scientific_name_gbif %in% Merged_traits$scientific_name_gbif), list(value = sum(value, na.rm = T)), by = list(Species,scientific_name_gbif)]

get_gbif_taxonomy('Polyommatus icarus')
Abundances_lepi_agg[!(scientific_name_gbif %in% Merged_traits$scientific_name_gbif) & value>0, ]

Merged_traits[grepl('icarus', scientific_name_gbif),]

### A few species are very abundant in a few plots but currently don't have data, leading to very low trait coverage.
# We fill them with additional data for some abundant species
# data from http://lepiforum.org
# and http://www.pyrgus.de
Lythria_purpuraria = data.table('scientific_name_gbif'= 'Lythria purpuraria', 'Flight_max' = 6, 'Wintering_stage' = 3, 'forewing_maximum' = 15, 'Generalism_use' = 1, 'Voltinism_use' = 2, 'B_Egg_number' = NA)
Zygaena_carniolica = data.table('scientific_name_gbif'= 'Zygaena carniolica', 'Flight_max' = 2, 'Wintering_stage' = 2, 'forewing_maximum' = 13, 'Generalism_use' = 3, 'Voltinism_use' = 1, 'B_Egg_number' = NA)

Merged_traits = rbind(Merged_traits, rbind(Lythria_purpuraria, Zygaena_carniolica))
Merged_traits[scientific_name_gbif == 'Aphantopus hyperantus', Generalism_use := 3] #http://www.pyrgus.de/Aphantopus_hyperantus_en.html

Merged_traits[, logSize := log(forewing_maximum)]
Merged_traits[, Eggs.log := log(B_Egg_number)]
Merged_traits[, logFlight := log(Flight_max)]

Merged_traits = Merged_traits[!is.na(scientific_name_gbif), lapply(.SD, mean, na.rm = T), by = scientific_name_gbif]

#Merged_traits = fwrite('Code/Data/Temporary_data/Lepidoptera_traits_matched.csv')

Merged_traits = fread('Code/Data/Temporary_data/Lepidoptera_traits_matched.csv')

Merged_traits[, logSize := Size.log]
Merged_traits[, Wintering_stage := wintering_stage]
Merged_traits[, logFlight := log(Flight_max)]
Merged_traits[, Generalism_use := generalism_use]
Merged_traits[, Voltinism_use := voltinism_use]

Merged_traits = Merged_traits[BE_species != '', list(logSize,Wintering_stage,logFlight,Generalism_use,Voltinism_use,BE_species)]


## Save species-matched trait data

Melted_traits = melt.data.table(Merged_traits, id.vars = 'BE_species', variable.name = 'traitName', value.name = 'traitValue')
# Add info

#logFlight      Wintering_stage Generalism_use  Voltinism_use  
traitUnits = c("month (log-transformed)","unitless","unitless","unitless",'mm (log-transformed)')
traitDescription = c("Length of the adult flying period",
                     "Wintering stage coded as egg = 0, larva = 1, pupa = 2, adult = 3",
                     "Feeding generalism coded as monophagous (feeding on 1 species) = 1, oligophagous (feeding on one genus) = 2, oligophagous (feeding on one family) = 3, polyphagous (more than one family) = 4)",
                     "Number of generations per year, coded as: semivoltine = 0.5, strictly univoltine = 1, univoltine + sometimes half a generation = 1.25, uni or multivoltine = 1.5, strictly multivoltine = 2",
                     'Body size: maximum length of the forewing')

traitDataID = c("NA","Bexis ID 21228","Bexis ID 21228","Bexis ID 21228", "Bexis ID 21228")

traitRef = c("doi:10.1007/s10531-017-1411-z, 10.1016/j.baae.2013.09.002,https://doi.org/10.5285/5b5a13b6-2304-47e3-9c9d-35237d1232c6,https://doi.org/10.1038/s41597-020-00697-7",
             "doi:10.1007/s10531-017-1411-z, 10.1016/j.baae.2013.09.002,https://doi.org/10.5285/5b5a13b6-2304-47e3-9c9d-35237d1232c6,https://doi.org/10.1038/s41597-020-00697-7",
             "doi:10.1007/s10531-017-1411-z, 10.1016/j.baae.2013.09.002,https://doi.org/10.5285/5b5a13b6-2304-47e3-9c9d-35237d1232c6,https://doi.org/10.1038/s41597-020-00697-7",
             "doi:10.1007/s10531-017-1411-z, 10.1016/j.baae.2013.09.002,https://doi.org/10.5285/5b5a13b6-2304-47e3-9c9d-35237d1232c6,https://doi.org/10.1038/s41597-020-00697-7",
             'https://doi.org/10.5285/5b5a13b6-2304-47e3-9c9d-35237d1232c6, Waring, P., Townsend, M., Lewington, R., 2017, Field guide to the Moths of Great Britain and Ireland, Third Edition, Bloomsbury Wildlife, London')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('logFlight',    'Wintering_stage', 'Generalism_use',  'Voltinism_use', 'logSize'  )

Melted_traits_all = add_info(Melted_traits[traitName %in% c('logFlight', 'Wintering_stage', 'logSize', 'Generalism_use', 'Voltinism_use', 'Eggs.log')], traitRefs = traitRef, traitDataIDs = traitDataID, traitDescriptions = traitDescription, traitUnits = traitUnits, traitsOnly = TRUE)

fwrite(Melted_traits_all, "Data/Temporary_data/Lepidoptera_traits_with_info.csv")
#fwrite(Abundances_lepi, "Data/Temporary_data/Lepidoptera_abundance_gbif.csv")



# ************************** #
#### 2. Species-level PCA ####
# ************************** #

# Check the existence of a species-level slow-fast axis
pca_but_sp = dudi.pca(Merged_traits[#gsub(' ', '_', scientific_name_gbif) %in% Abundance_all$Species & 
  complete.cases(Merged_traits), 
                                    list('Flight_period.log' = logFlight, 
                                         'Wintering_stage' = Wintering_stage, 
                                         logSize, 
                                         'Generalism' = Generalism_use, 
                                         'Volitnism' = Voltinism_use)], scannf = FALSE, nf = 2)

gg_but = fviz_pca(pca_but_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                     col.ind = "steelblue",
                     fill.ind = "white",
                     col.var = "black")
ggsave(gg_but, file = 'Results/species_pca_but.pdf', width = 6, height = 5)


# ********************************** #
#### 3. Community-weighted traits ####
# ********************************** #

Merged_traits[#!is.na(scientific_name_gbif) & scientific_name_gbif %in% Abundances_lepi$scientific_name_gbif,-1][
  , lapply(.SD, function(x){length(x[!is.na(x) & x != 'NA'])})]
butterfly_CC = check_coverage(Merged_traits[!is.na(BE_species),], Abundances_lepi[!is.na(value) & value != 0,], c('Wintering_stage', 'logSize', 'logFlight','Generalism_use', 'Voltinism_use'),
                              'BE_species', 'Species')

butterfly_CC[, -c(1,2)][, lapply(.SD, median)]
butterfly_cwm = my_cwm(Merged_traits[!is.na(BE_species),], Abundances_lepi[!is.na(value) & value != 0,], c('logFlight', 'Wintering_stage', 'logSize', 'Generalism_use', 'Voltinism_use'),
                       'BE_species', 'Species')


# Melt and merge
CWM_CC_butterflies = merge.data.table(melt.data.table(butterfly_cwm[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                      melt.data.table(butterfly_CC[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
#logFlight      Wintering_stage Generalism_use  Voltinism_use  


CWM_CC_butterflies = add_info(CWM_CC_butterflies, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_butterflies, "Data/CWM_data/CWM_butterflies.csv")

#write.csv(butterfly_cwm[ , list(
#  "but_Generalism" = Generalism_use,
#  "but_Size" = logSize,
#  "but_GenYear" = Voltinism_use,
#  "but_hibernation" = Wintering_stage,
#  "but_flight" = logFlight,
#  'but_Eggs' = Eggs.log,
#  "Plot" = Plot   ,
#  "Year" = Year), ], paste(cwm_path, "CWM_Butterflies.csv", sep = ''))
#

# Unweighted: changes with species presence-absence, not abundance
Abundances_lepi_presence_absence = Abundances_lepi[, list(value = sum(value, na.rm = T)), by = list(Plot,Year, Species)]
Abundances_lepi_presence_absence[value>1, value := 1]
butterfly_CC_noweight = check_coverage(Merged_traits[!is.na(BE_species),], Abundances_lepi_presence_absence[!is.na(value) & value != 0,], c('Wintering_stage', 'logFlight','logSize',  'Generalism_use', 'Voltinism_use'),
                              'BE_species', 'Species')

butterfly_cwm_noweight = my_cwm(Merged_traits[!is.na(BE_species),], Abundances_lepi_presence_absence[!is.na(value) & value != 0,], c('logFlight', 'Wintering_stage', 'logSize', 'Generalism_use', 'Voltinism_use'),
                       'BE_species', 'Species')


# Melt and merge
CWM_CC_butterflies_noweight = merge.data.table(melt.data.table(butterfly_cwm_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                      melt.data.table(butterfly_CC_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info
CWM_CC_butterflies_noweight = add_info(CWM_CC_butterflies_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_butterflies_noweight, "Data/CWM_data/CWM_butterflies_noweight.csv")

#######################################
# Check turnover accross LUI gradient #
#######################################

data_lui <- fread("Environment/LUI_input_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = dcast(Abundances_lepi[!is.na(Species) & value>0,],  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)

