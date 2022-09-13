# This script takes as input the abundances and species-level traits of protist taxa
# and outputs a CWM matrix averaged for all considered years.
library(car)
library(readr)
library(data.table)
library(reshape2)
library(vegan)
library(Hmisc)

setwd("~/Data")

figures_path = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/'
cwm_path = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data'
matched_traits_path = '~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Matched_trait_datasets/'


# ############### #
#### Load data ####
# ############### #

# Size data provided by Kenneth Dumack, 2021
Sizes = list('Filoreta'='≥51µm', 'Gromia'='≥51µm', 'Plasmodiophora'='≤10µm', 'Woronina'='≤10µm', 'Polymyxa-lineage_X'='≤10µm', 'Ligniera'='≤10µm',
             'Sorosphaerula'='≤10µm', 'Sorosphaera'='≤10µm', 'Polymyxa'='≤10µm', 'Spongospora_lineage_X'='≤10µm', 'Spongospora'='≤10µm', 'Phagomyxida_XX'='≤10µm',
             'Phagomyxa'='≤10µm', 'Hillenburgia'='≤10µm', 'Maullinia'='≤10µm', 'Arachnula'='≥51µm', 'Leptophryidae_X'='≥51µm', 'Leptophrys'='≥51µm',
             'Theratromyxa'='≥51µm', 'Platyreta'='≥51µm', 'Vernalophrys'='≥51µm', 'Planctomyxa'='31µm-50µm', 'Arachnomyxa'='31µm-50µm',
             'Thalassomyxa-lineage_X'='≥51µm', 'Vampyrellidae_X'='31µm-50µm', 'Vampyrella'='31µm-50µm', 'Thalassomyxa'='≥51µm', 'Placopodidae_X'='11µm-30µm',
             'Placopus'='11µm-30µm', 'Novel-clade-10_XX'='≤10µm', 'Aquavolonida_XX'='≤10µm', 'Aquavolon'='≤10µm', 'Tremulida_XX'='≤10µm', 'Tremula'='≤10µm',
             'Massisteria'='≤10µm', 'Minimassisteria'='≤10µm', 'Mesofila'='11µm-30µm', 'Nanofila'='≤10µm', 'Limnofila'='≤10µm', 'Clathrulina'='31µm-50µm',
             'Hedriocystis'='11µm-30µm', 'Metromonadea_XXX'='11µm-30µm', 'Metopion'='≤10µm', 'Metromonas'='≤10µm', 'Micrometopion'='≤10µm', 'Minorisa'='≤10µm',
             'Bigelowiella'='≤10µm', 'Norrisiella'='≤10µm', 'Chlorarachnion'='11µm-30µm', 'Partenskyella'='≤10µm', 'Gymnochlora'='≤10µm', 'Amorphochlora'='≤10µm',
             'Lotharella'='≤10µm', 'Leptogromia'='≥51µm', 'Imbricatea_Novel-clade-3_X'='≥51µm', 'Assulina'='≥51µm', 'Euglypha'='≥51µm', 'Paulinella'='11µm-30µm',
             'Ovulinata'='11µm-30µm', 'Micropyxidiella'='≤10µm', 'Trachelocorythion'='11µm-30µm', 'Sphenoderia'='≥51µm', 'Trinematidae_X'='≥51µm',
             'Trinema'='31µm-50µm', 'Corythion'='31µm-50µm', 'Cyphoderiidae_X'='≥51µm', 'Cyphoderia'='≥51µm', 'Tracheleuglypha'='≥51µm', 'Placocista'='≥51µm',
             'Discomonadida_XX'='11µm-30µm', 'Discomonadidae_X'='11µm-30µm', 'Discomonas'='11µm-30µm', 'Krakenida_XX'='≤10µm', 'Krakenidae_X'='≤10µm',
             'Kraken'='≤10µm', 'Variglissida_XX'='≤10µm', 'Clautriavia'='11µm-30µm', 'Quadricilia'='11µm-30µm', 'Nudifilidae_X'='≤10µm', 'Nudifila'='≤10µm',
             'Auranticordis'='≥51µm', 'Pseudopirsonia'='≤10µm', 'Abollifer'='11µm-30µm', 'Cyranomonas'='≤10µm', 'Spongomonadida_XX'='≤10µm',
             'Spongomonadidae_X'='≤10µm', 'Spongomonas'='≤10µm', 'Peregrinia'='11µm-30µm', 'Esquamula'='≤10µm', 'Penardeugenia'='11µm-30µm',
             'Thaumatomonas'='11µm-30µm', 'Thaumatomastix'='11µm-30µm', 'Thaumatospina'='11µm-30µm', 'Ovaloplaca'='11µm-30µm', 'Scutellomonas'='11µm-30µm',
             'Allas'='11µm-30µm', 'Reckertia'='11µm-30µm', 'Cavernomonas'='≤10µm', 'Cercomonas'='≤10µm', 'Eocercomonas'='≤10µm', 'Filomonas'='≤10µm',
             'Neocercomonas'='11µm-30µm', 'Paracercomonadidae_X'='11µm-30µm', 'Brevimastigomonas'='≤10µm', 'Phytocercomonas'='≤10µm', 'Metabolomonas'='≤10µm',
             'Nucleocercomonas'='≤10µm', 'Paracercomonas'='≤10µm', 'Allapsidae_X'='≤10µm', 'Allapsidae_Group_Te'='≤10µm', 'Allantion'='≤10µm', 
             'Allapsa'='≤10µm', 'Teretomonas'='≤10µm', 'Viridiraptoridae_X'='11µm-30µm', 'Viridiraptor'='11µm-30µm', 'Orciraptor'='11µm-30µm', 
             'Bodomorphidae_X'='≤10µm', 'Bodomorpha'='≤10µm', 'Dujardinidae_X'='≤10µm', 'Dujardina'='≤10µm', 'Proleptomonadidae_X'='≤10µm',
             'Proleptomonas'='≤10µm', 'Sandonidae_X'='≤10µm', 'Sandonidae_Clade_N_F_A'='≤10µm', 'Flectomonas'='≤10µm', 'Mollimonas'='≤10µm', 'Neoheteromita'='≤10µm',
             'Sandona'='≤10µm', 'Pansomonadida_XX'='≤10µm', 'Agitatidae_X'='≤10µm', 'Aurigamonas'='≤10µm', 'Agitata'='≤10µm', 'Sainouridea_XX'='≤10µm',
             'Sainouridae_X'='≤10µm', 'Cholamonas'='≤10µm', 'Sainouron'='≤10µm', 'Acantholus'='≤10µm', 'Homocognatis'='≤10µm', 'Helkesimastigidae_X'='≤10µm',
             'Helkesimastix'='≤10µm', 'Guttulinopsidae_X'='≤10µm', 'Olivorum'='≤10µm', 'Guttulinopsis'='≤10µm', 'Rosculus'='≤10µm', 'Puppisaman'='≤10µm',
             'Phaeodarea_XX'='≥51µm', 'Metusettidae_X'='≥51µm', 'Gazeletta'='≥51µm', 'Challengeriidae_X'='≥51µm', 'Chalengeron'='≥51µm', 'Protcystis'='≥51µm',
             'Entocannula'='≥51µm', 'Porospathidae_X'='≥51µm', 'Porospathis'='≥51µm', 'Tuscaroridae_X'='≥51µm', 'Tuscaretta'='≥51µm', 'Aulosphaerida_X'='≥51µm',
             'Auloscena'='≥51µm', 'Aulosphaera'='≥51µm', 'Sagosphaeridae_X'='≥51µm', 'Sagoscena'='≥51µm', 'Medusettidae_X'='≥51µm', 'Medusetta'='≥51µm',
             'Coelodendridae_X'='≥51µm', 'Coelodendrum'='≥51µm', 'Conchariidae_X'='≥51µm', 'Conchellium'='≥51µm', 'Aulacanthidae_X'='≥51µm', 'Aulacantha'='≥51µm',
             'Auloceros'='≥51µm', 'Aulographis'='≥51µm', 'Phaeodinidae_X'='≥51µm', 'Phaeodina'='≥51µm', 'Cryothecomonas-lineage_X'='≤10µm', 'Cryothecomonadidae_X'='≤10µm',
             'Cryothecomonas'='≤10µm', 'Protaspa-lineage_X'='31µm-50µm', 'Protaspididae_X'='31µm-50µm', 'Protaspa'='31µm-50µm', 'Rhogostoma-lineage_X'='≤10µm', 'Rhogostomidae_X'='≤10µm',
             'Capsellina'='≤10µm', 'Rhogostoma'='≤10µm', 'Sacciforma'='11µm-30µm', 'Ventricleftida_CCW10-lineage_X'='31µm-50µm', 'Ventricleftida_XX'='31µm-50µm', 'Ventricleftida_CCW10-lineage_X'='31µm-50µm',
             'Ventrifissuridae_X'='31µm-50µm', 'Verrucomonas'='31µm-50µm', 'Ventrifissura'='31µm-50µm', 'Ebria'='11µm-30µm', 'Botuliforma'='31µm-50µm', 'Thecofilosea_Novel-clade-4_X'='31µm-50µm',
             'Pseudodifflugia'='31µm-50µm', 'Rhizaspididae_X'='31µm-50µm', 'Rhizaspis'='31µm-50µm', 'Fisculla'='11µm-30µm', 'Lecythium'='31µm-50µm', 'Diaphoropodon'='31µm-50µm', 'Trachyrhizium'='11µm-30µm')



Cercozoa_abundance <- data.table(read_excel("Abundances/Protists/Data_Sheet_3_Contrasting_Responses_of_Protistan_Plant_Parasites_and_Phagotrophs_to_Ecosystems_Land_Management_and_Soil_Properties/TableS5.xlsx"))
# Accessible from https://www.frontiersin.org/articles/10.3389/fmicb.2020.01823/full, Table S5
# Also includes trait data based on https://github.com/Kenneth-Dumack/Functional-Traits-Cercozoa-Endomyxa

colnames(Cercozoa_abundance) = as.character(Cercozoa_abundance[1,])
Cercozoa_abundance = Cercozoa_abundance[-1,]

### Reformat abundance data
Cercozoa_abundance_format = melt.data.table(Cercozoa_abundance, id.vars = colnames(Cercozoa_abundance)[c(1:3,595:605)], 
                     value.name = 'value',
                     variable.name = 'Plot_Year')
Cercozoa_abundance_format[, c('Plot', 'Year') := list(paste(substr(Plot_Year, 1, 3), substr(Plot_Year, 5, 6), sep = ''),
                     substr(Plot_Year, 8, 9))]


Cercozoa_traits = unique(Cercozoa_abundance_format[, list(Genus, nutrition, morphology, locomotion, Plot_Year, Plot, Year)])
Cercozoa_traits[, Size := Sizes[Genus]]

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
                       dplyr::recode(locomotion, !!!c(locomotion_unknown = NA,
                                                      '0' = NA,
                                                      freely_swimming = 'freely_swimming',
                                                      'non-motile_endoparasite' = 'non_motile',
                                                      'creeping/gliding' = 'gliding')),
                       dplyr::recode(locomotion, !!!c(locomotion_unknown = NA,
                                                      freely_swimming = 2,
                                                      '0' = NA,
                                                      'non-motile_endoparasite' = 0,
                                                      'creeping/gliding' = 1)),
                       dplyr::recode(gsub('[/()]', '', morphology), !!!c(
                         'NA' = NA,
                         morphology_unknown = NA,
                         '0' = NA,
                         endoparasite = 'endo',
                         "testate_amoeboflagellate/amoeba_(silica)"  = 'testate',
                         'testate_amoeboflagellate/amoeba_(organic/agglutinated)' = 'testate',
                         "flagellate/intracellular" = 'naked',
                         testate_cell_organic = 'testate',
                         naked_amoeboflagellate = 'naked',
                         naked_amoeba = 'naked',
                         naked_flagellate = 'naked')),
                       dplyr::recode(nutrition, 
                                     !!!c('NA' = NA,
                                          "0" = NA,
                                          bacterivore = 'bacterial_cons',
                                          animal_parasite = 'secondary_cons',
                                          plant_parasite = 'primary_cons',
                                          nutrition_unknown = NA,
                                          eukaryvore = 'secondary_cons',
                                          omnivore = 'secondary_cons',
                                          parasite_not_plant = 'secondary_cons',
                                          autotroph = 'primary_prod')))]

Cercozoa_traits[, Naked_amoeba := as.numeric(morphology == 'naked_amoeba')]
Cercozoa_traits[, Naked := as.numeric(morpho_code == 'naked')]
trait_names = c("morpho_code",
                "nutrition_code","Size",  'Naked_amoeba')

#### Aggregate abundance to Genus
Cercozoa_abundance_format[, nutrition_code:=  dplyr::recode(nutrition, 
                                                            !!!c('NA' = NA,
                                                                 bacterivore = 'bacterial_cons',
                                                                 animal_parasite = 'secondary_cons',
                                                                 plant_parasite = 'primary_cons',
                                                                 nutrition_unknown = NA,
                                                                 eukaryvore = 'secondary_cons',
                                                                 omnivore = 'secondary_cons',
                                                                 not_plant_parasite = 'secondary_cons',
                                                                 autotroph = 'primary_prod'))]


# ############## #
#### Output ####
# ############## #


# Species-level strategies
# Bacterial consumers
data_bact = Cercozoa_traits[nutrition_code == 'bacterial_cons' & Genus %in% Cercozoa_abundance_format[value>0 & grepl('G', Plot), unique(Genus)] & Sizes != 'NA',]
plot_bact = ggplot(data_bact, aes(Sizes, morpho_code)) + ylab('Morphology') + xlab('Size') +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) + theme_bw()
cor.test(data_bact$Naked, data_bact$Size)

# Secondary consumers
data_sec = Cercozoa_traits[nutrition_code == 'secondary_cons' & Genus %in% Cercozoa_abundance_format[value>0 & grepl('G', Plot), unique(Genus)] & Sizes != 'NA',]
plot_sec = ggplot(data_sec, aes(Sizes, morpho_code)) + ylab('Morphology') + xlab('Size') +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) + theme_bw()
cor.test(data_sec$Naked, data_sec$Size)

cor.test(data_sec$Size, data_sec$Naked_amoeba)


### Coverage
CC_Protists = check_coverage(unique(Cercozoa_traits[,.SD, .SDcols = !c('Family')]), 
                             Cercozoa_abundance_format, trait_names, 'Genus', 'Genus')
CC_Protists_prim_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'primary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'primary_cons',], trait_names[trait_names != 'nutrition_code']
                                , 'Genus', 'Genus')
CC_Protists_bacterial_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'bacterial_cons',.SD, .SDcols = !c('Family')]), 
                                        Cercozoa_abundance_format[nutrition_code == 'bacterial_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')


Cercozoa_traits[nutrition_code == 'secondary_cons', length(Size[!is.na(Size)])]
Cercozoa_abundance_format[nutrition_code == 'secondary_cons',]
CC_Protists_sec_cons <- check_coverage(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = !c('Family')]), 
                                Cercozoa_abundance_format[nutrition_code == 'secondary_cons',], trait_names[trait_names != 'nutrition_code'], 'Genus', 'Genus')


### CWM
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
                                        "Year" = Year)], paste(cwm_path, 'CWM_Protists.csv',sep = ''))
write.csv(CWM_Protists_bact[, list("Pb_Naked" = morpho_code_naked / (morpho_code_testate+morpho_code_naked),
                                        "Pb_Size" = Size,
                                        "Plot" = Plot,
                                        "Year" = Year)],  paste(cwm_path, 'CWM_Protists_bact.csv',sep = ''))

write.csv(CWM_Protists_sec_cons[, list("Ps_Naked" = morpho_code_naked / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Size" = Size,
                                       "Plot" = Plot,
                                       "Year" = Year)],  paste(cwm_path, 'CWM_Protists_sec_cons.csv',sep = ''))


# Non-weighted community traits

Cercozoa_abundance_format_presence_absence = Cercozoa_abundance_format[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Genus)]
Cercozoa_abundance_format_presence_absence[value>1, value := 1]

CWM_Protists_noweight <- my_cwm(unique(Cercozoa_traits[ ,.SD, .SDcols = c('Genus',trait_names)]), 
                                Cercozoa_abundance_format_presence_absence, 
                       trait_names, 'Genus', 'Genus')
CWM_Protists_bact_noweight <- my_cwm(unique(Cercozoa_traits[nutrition_code %in% c('bacterial_cons') ,.SD, .SDcols = c('Genus',trait_names)]), 
                                     Cercozoa_abundance_format_presence_absence, 
                           trait_names, 'Genus', 'Genus')
CWM_Protists_sec_cons_noweight <- my_cwm(unique(Cercozoa_traits[nutrition_code == 'secondary_cons',.SD, .SDcols = c('Genus',trait_names)]), 
                                         Cercozoa_abundance_format_presence_absence,
                                trait_names, 'Genus', 'Genus')

write.csv(CWM_Protists_noweight[, list("P_patho" = nutrition_code_primary_cons/(nutrition_code_primary_prod+nutrition_code_primary_cons +nutrition_code_bacterial_cons+nutrition_code_secondary_cons),
                              "P_naked" = morpho_code_naked / (morpho_code_naked + morpho_code_testate),
                              "P_Size" = Size,
                              "Plot" = Plot,
                              "Year" = Year)], paste(cwm_path, "CWM_Protists_noweight.csv",sep = ''))
write.csv(CWM_Protists_bact_noweight[, list("Pb_Naked" = morpho_code_naked / (morpho_code_testate+morpho_code_naked),
                                   "Pb_Size" = Size,
                                   "Plot" = Plot,
                                   "Year" = Year)], paste(cwm_path, "CWM_Protists_bact_noweight.csv",sep = ''))

write.csv(CWM_Protists_sec_cons_noweight[, list("Ps_Naked" = morpho_code_naked / (morpho_code_testate+ morpho_code_naked + morpho_code_endo),
                                       "Ps_Size" = Size,
                                       "Plot" = Plot,
                                       "Year" = Year)], paste(cwm_path, "CWM_Protists_sec_cons_noweight.csv",sep = ''))


# ############## #
#### Turnover ####
# ############## #
data_lui <- fread("Environment/LUI_input_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years

data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

library(betapart)
comm.test_patho = dcast(Cercozoa_abundance_format[nutrition_code == 'primary_cons', list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Genus)],  Plot~Genus, value.var = 'value', fill = 0)
comm.test_bact = dcast(Cercozoa_abundance_format[nutrition_code == 'bacterial_cons', list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Genus)],  Plot~Genus, value.var = 'value', fill = 0)
comm.test_sec = dcast(Cercozoa_abundance_format[nutrition_code == 'secondary_cons', list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Genus)],  Plot~Genus, value.var = 'value', fill = 0)

rownames(comm.test_patho)= comm.test_patho$Plot
rownames(comm.test_bact)= comm.test_bact$Plot
rownames(comm.test_sec)= comm.test_sec$Plot
comm.test_patho = comm.test_patho[,-1]
comm.test_bact = comm.test_bact[,-1]
comm.test_sec = comm.test_sec[,-1]

beta.multi.abund(comm.test_patho)
beta.multi.abund(comm.test_bact)
beta.multi.abund(comm.test_sec)


comm_patho_min_max = matrix(c(colSums(comm.test_patho[min_lui_plots,]),colSums(comm.test_patho[max_lui_plots,])), nrow = 2)
comm_bact_min_max = matrix(c(colSums(comm.test_bact[min_lui_plots,]),colSums(comm.test_bact[max_lui_plots,])), nrow = 2)
comm_sec_min_max = matrix(c(colSums(comm.test_sec[min_lui_plots,]),colSums(comm.test_sec[max_lui_plots,])), nrow = 2)


beta.multi.abund(comm_patho_min_max)
beta.multi.abund(comm_bact_min_max)
beta.multi.abund(comm_sec_min_max)


