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

Cercozoa_traits[, Size_max := as.numeric(Recode(Sizes,
                                    " '≤10µm' = 10;
                                     '≥51µm' = 70;
                                     '11µm-30µm' = 30;
                                     '31µm-50µm'  = 50"
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


trait_names = c("morpho_code",
                "nutrition_code","Size_max" )

#### Aggregate abundance to Genus_2
Cercozoa_abundance_format = Cercozoa_abundance_format[, list(value = sum(as.numeric(value))), by = c("Order", "Family","Genus", 
                                                                                                "Plot", "Year") ]

# Add log size
#Cerco_traits_estimates$log_size = log(Cerco_traits_estimates$size_estimates)
#trait_names = c(trait_names, 'log_size' )


# Genus-level PCA
data_bact= Cercozoa_traits[nutrition_code == 'bacterial_cons',.SD, .SDcols = c('locomotion_code', 'morpho_code', 'nutrition_code', 'Size_max_inclplasmodia')]
data_bact$locomotion_code = factor(data_bact$locomotion_code)
data_bact$morpho_code = factor(data_bact$morpho_code)
data_bact$nutrition_code = factor(data_bact$nutrition_code)
pca_bact = dudi.mix(data_bact[complete.cases(data_bact),-3], , scannf = FALSE, nf = 4)
draw_dudi_mix(Data = pca_bact, c(1, 3), save = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Ind_level/Prot_bact13bis.pdf', 
              new_var_names = c('Freely_swimming', 'Gliding', 'Morpho:naked', 'Morpho:testate', 'Size'))
data_secondary= Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = c('locomotion_code', 'morpho_code', 'nutrition_code', 'Size_max_inclplasmodia')]
data_secondary$locomotion_code = factor(data_secondary$locomotion_code)
data_secondary$morpho_code = factor(data_secondary$morpho_code)
data_secondary$nutrition_code = factor(data_secondary$nutrition_code)
pca_secondary = dudi.mix(data_secondary[complete.cases(data_secondary),-3], , scannf = FALSE, nf = 4)
draw_dudi_mix(Data = pca_secondary, c(1, 2), save = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Ind_level/Prot_secondary12.pdf', 
              new_var_names = c('Freely_swimming', 'Gliding', 'Non-motile','Endo', 'Morpho:naked', 'Morpho:testate', 'Size'))

# Cover
CC_Protists = check_coverage(unique(Cercozoa_traits[,.SD, .SDcols = !c('Family')]), 
                             Cercozoa_abundance_format, trait_names, 'Genus', 'Genus')
CC_Protists_prim_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'primary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'primary_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')
CC_Protists_bacterial_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'bacterial_cons',.SD, .SDcols = !c('Family')]), 
                                        Cercozoa_abundance_format[nutrition_code == 'bacterial_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')
CC_Protists_sec_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'secondary_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')

CWM_Protists <- my_cwm(unique(Cercozoa_traits[, .SD, .SDcols = c('Genus',trait_names)]),
                               Cercozoa_abundance_format, trait_names, 'Genus', 'Genus')
CWM_Protists_prim_bact<- my_cwm(unique(Cercozoa_traits[nutrition_code %in% c('primary_cons','bacterial_cons') ,.SD, .SDcols = c('Genus',trait_names)]), 
                                Cercozoa_abundance_format, 
                                trait_names, 'Genus', 'Genus')
CWM_Protists_bacterial_cons <- my_cwm(unique(Cercozoa_traits[nutrition_code == 'bacterial_cons',.SD, .SDcols = c('Genus',trait_names)]), 
                                             Cercozoa_abundance_format, 
                                      trait_names, 'Genus', 'Genus')
CWM_Protists_sec_cons <- my_cwm(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = c('Genus',trait_names)]), 
                            Cercozoa_abundance_format,
                            trait_names, 'Genus', 'Genus')

write.csv(CWM_Protists[, list("P_Naked" = morpho_code_naked / (morpho_code_testate+morpho_code_naked+morpho_code_NA+morpho_code_endo),
                              "P_Teste" = morpho_code_testate / (morpho_code_testate+morpho_code_naked+morpho_code_NA+morpho_code_endo),
                              "P_Endo" = morpho_code_endo / (morpho_code_testate+morpho_code_naked+morpho_code_NA+morpho_code_endo),
                              "P_Bact" = nutrition_code_bacterial_cons / (nutrition_code_NA +  nutrition_code_bacterial_cons +nutrition_code_primary_cons +nutrition_code_secondary_cons),
                              #"P_Parasite" = nutrition_code_not_plant_parasite/ (nutrition_code_NA +  nutrition_code_bacterial_cons +nutrition_code_primary_cons +nutrition_code_secondary_cons),
                              "P_Primary" = nutrition_code_primary_cons/ (nutrition_code_NA +  nutrition_code_bacterial_cons +nutrition_code_primary_cons +nutrition_code_secondary_cons),
                              "P_Secondary" = nutrition_code_secondary_cons/ (nutrition_code_NA +  nutrition_code_bacterial_cons +nutrition_code_primary_cons +nutrition_code_secondary_cons),
                              "P_Size" = Size_max,
                              "Plot" = Plot,
                              "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists.csv')
write.csv(CWM_Protists_prim_bact[, list("Ppb_Naked" = morpho_code_naked / (morpho_code_testate+morpho_code_naked+morpho_code_endo),
                                        "Ppb_Teste" = morpho_code_testate / (morpho_code_testate+morpho_code_naked+morpho_code_endo),
                                        "Ppb_Endo" = morpho_code_endo / (morpho_code_testate+morpho_code_naked+morpho_code_endo),
                                        "Ppb_Primary" = nutrition_code_primary_cons/ (nutrition_code_bacterial_cons +nutrition_code_primary_cons),
                                        "Ppb_Size" = Size_max,
                                        "Plot" = Plot,
                                        "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists_prim_bact.csv')
write.csv(CWM_Protists_bacterial_cons[, list("Pb_Teste" = morpho_code_testate / (morpho_code_testate+morpho_code_naked),
                                             "Pb_Size" = Size_max,
                                             "Plot" = Plot,
                                             "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists_bacterial_cons.csv')
write.csv(CWM_Protists_sec_cons[, list("Ps_Teste" = morpho_code_testate / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Endo" = morpho_code_endo / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Naked" = morpho_code_naked / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Size" = Size_max,
                                       "Plot" = Plot,
                                       "Year" = Year)], '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Protists_sec_cons.csv')