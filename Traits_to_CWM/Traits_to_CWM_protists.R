# This script takes as input the abundances and species-level traits of PROTIST species
# and outputs a CWM matrix averaged for all considered years.
library(car)
library(readr)
library(data.table)
library(reshape2)
library(vegan)
library(Hmisc)

Cercozoa_traits <- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Protists/FunctionalTraitsTable_Size_Dumack.xlsx"))
Cercozoa_abundance <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Protists/Cercozoa\ copie2.txt")

colnames(Cercozoa_traits) = capitalize(tolower(colnames(Cercozoa_traits)))
colnames(Cercozoa_traits) = gsub(' ', '_', colnames(Cercozoa_traits))
colnames(Cercozoa_traits) = gsub('[()\\./]', '', colnames(Cercozoa_traits))

### Reformat abundance data
Cercozoa_abundance_format = melt.data.table(Cercozoa_abundance, id.vars = colnames(Cercozoa_abundance)[c(1:3,595:605)], 
                     value.name = 'value',
                     variable.name = 'Plot_Year')
Cercozoa_abundance_format[, c('Plot', 'Year') := list(paste(substr(Plot_Year, 1, 3), substr(Plot_Year, 5, 6), sep = ''),
                     substr(Plot_Year, 8, 9))]

Cercozoa_traits[, Size := as.numeric(Recode(Sizes,
                                    " '≤10µm' = 1;
                                     '≥51µm' = 4;
                                     '11µm-30µm' = 2;
                                     '31µm-50µm'  = 3"
                                     ))]
# Recoding the main traits
Cercozoa_traits[, c('locomotion_code',
                     'locomotion_num',
                     'morpho_code',
                     'nutrition_code'):=  list(
                       dplyr::recode(Locomotion, !!!c(locomotion_unknown = NA,
                                                      freely_swimming = 'freely_swimming',
                                                      'non_motile' = 'non_motile',
                                                      'gliding' = 'gliding')),
                       dplyr::recode(Locomotion, !!!c(locomotion_unknown = NA,
                                                      freely_swimming = 2,
                                                      'non_motile' = 0,
                                                      'gliding' = 1)),
                       dplyr::recode(gsub('[/()]', '', Morphology), !!!c(
                         'NA' = NA,
                         endoparasite = 'endo',
                         testate_cell_silica = 'testate',
                         testate_cell_organic = 'testate',
                         naked_amoeboflagellate = 'naked',
                         naked_amoeba = 'naked',
                         naked_flagellate = 'naked')),
                       dplyr::recode(Nutrition, 
                                     !!!c('NA' = NA,
                                          bacterivore = 'bacterial_cons',
                                          animal_parasite = 'secondary_cons',
                                          plant_parasite = 'primary_cons',
                                          nutrition_unknown = NA,
                                          eukaryvore = 'secondary_cons',
                                          omnivore = 'secondary_cons',
                                          not_plant_parasite = 'secondary_cons',
                                          autotroph = 'primary_prod')))]

Cercozoa_traits[, Naked_amoeba := as.numeric(Morphology == 'naked_amoeba')]
Cercozoa_traits[, Naked := as.numeric(morpho_code == 'naked')]
trait_names = c("morpho_code",
                "nutrition_code","Size",  'Naked_amoeba')

#### Aggregate abundance to Genus_2
Cercozoa_abundance_format = Cercozoa_abundance_format[, list(value = sum(as.numeric(value))), by = c("Order", "Family","Genus", 
                                                                                                "Plot", "Year") ]



# Species-level strategies
# Bacterial consumers
data_bact = Cercozoa_traits[nutrition_code == 'bacterial_cons' & Genus %in% Cercozoa_abundance_format[value>0 & grepl('F', Plot), unique(Genus)] & Sizes != 'NA',]
plot_bact = ggplot(data_bact, aes(Sizes, morpho_code)) + ylab('Morphology') + xlab('Size') +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) + theme_bw()
cor.test(data_bact$Naked, data_bact$Size)

# Secondary consumers
data_sec = Cercozoa_traits[nutrition_code == 'secondary_cons' & Genus %in% Cercozoa_abundance_format[value>0 & grepl('F', Plot), unique(Genus)] & Sizes != 'NA',]
plot_sec = ggplot(data_sec, aes(Sizes, morpho_code)) + ylab('Morphology') + xlab('Size') +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) + theme_bw()
cor.test(data_sec$Naked, data_sec$Size)

cor.test(data_sec$Size, data_sec$Naked_amoeba)


# Cover
CC_Protists = check_coverage(unique(Cercozoa_traits[,.SD, .SDcols = !c('Family')]), 
                             Cercozoa_abundance_format, trait_names, 'Genus', 'Genus')
CC_Protists_prim_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'primary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'primary_cons',], trait_names[trait_names != 'nutrition_code']
                                , 'Genus', 'Genus')
CC_Protists_bacterial_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'bacterial_cons',.SD, .SDcols = !c('Family')]), 
                                        Cercozoa_abundance_format[nutrition_code == 'bacterial_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')
CC_Protists_sec_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'secondary_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')


# CWM
CWM_Protists <- my_cwm(unique(Cercozoa_traits[ ,.SD, .SDcols = c('Genus',trait_names)]), 
                                Cercozoa_abundance_format, 
                                trait_names, 'Genus', 'Genus')

test = merge(env_data_lui, CWM_Protists[, list(P_patho = mean(nutrition_code_primary_cons/(nutrition_code_primary_prod+nutrition_code_primary_cons +nutrition_code_bacterial_cons+nutrition_code_secondary_cons)),
                                               Naked_amoeba_1), by = Plot], by = 'Plot')

cor.test(test$Naked_amoeba_1, test$LUI)

CWM_Protists_bact<- my_cwm(unique(Cercozoa_traits[nutrition_code %in% c('bacterial_cons') ,.SD, .SDcols = c('Genus',trait_names)]), 
                                Cercozoa_abundance_format, 
                                trait_names, 'Genus', 'Genus')
CWM_Protists_sec_cons <- my_cwm(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = c('Genus',trait_names)]), 
                            Cercozoa_abundance_format,
                            trait_names, 'Genus', 'Genus')


cor.test(CWM_Protists_bact$morpho_code_naked, CWM_Protists_bact$Size)
cor.test(CWM_Protists_sec_cons$morpho_code_naked, CWM_Protists_sec_cons$Size)
cor.test(CWM_Protists_sec_cons$morpho_code_naked, CWM_Protists_sec_cons$Size)

test = merge(env_data_lui, CWM_Protists_sec_cons[, mean(Naked_amoeba_1), by = Plot], by = 'Plot')
cor.test(test$LUI, test$V1)

write.csv(CWM_Protists[, list("P_patho" = nutrition_code_primary_cons/(nutrition_code_primary_prod+nutrition_code_primary_cons +nutrition_code_bacterial_cons+nutrition_code_secondary_cons),
                              "P_naked" = morpho_code_naked / (morpho_code_naked + morpho_code_testate),
                              "P_Size" = Size,
                                        "Plot" = Plot,
                                        "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists.csv')
write.csv(CWM_Protists_bact[, list("Pb_Naked" = morpho_code_naked / (morpho_code_testate+morpho_code_naked),
                                        "Pb_Size" = Size,
                                        "Plot" = Plot,
                                        "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists_bact.csv')

write.csv(CWM_Protists_sec_cons[, list("Ps_Naked" = morpho_code_naked / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Size" = Size,
                                       "Plot" = Plot,
                                       "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists_sec_cons.csv')
