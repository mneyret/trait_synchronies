# This script takes as input the abundances and species-level traits of bat species found 
# in the Exploratories grasslands and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.

library(splitstackshape)
library(plyr)
library(stringr)
library(data.table)
library(readxl)
library(traitdataform)
library(ade4)
library(factoextra)
library(FD)
library(mice)
library(betapart)

setwd("~/Data")

figures_path = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Figures/'
cwm_path = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data'
matched_traits_path = '/Users/Margot/Desktop/Research/Senckenberg/Documents/Papers/Traits/Matched_traits_datasets/'

# ############### #
#### Load data ####
# ############### #

### Abundances 
## Raw diversity
allsp <- fread("Abundances/210112_EP_species_diversity_GRL_BEXIS.txt") # https://www.bexis.uni-jena.de/ddm/data/Showdata/27706
allsp$Species = gsub('_$', '', allsp$Species ) # Clean up species names
## Species information
fgs <- fread("Abundances/210112_EP_species_info_GRL_BEXIS.txt") # https://www.bexis.uni-jena.de/ddm/data/Showdata/27707
fgs$Species = gsub(' $', '', fgs$Species )
fgs$Species = gsub(' ', '_', fgs$Species ) # Clean up species names

Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE) # Merge abundance and info datasets
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))] # Homogenise up plot ID


# Subset the data by bats actually found in Germany, based on: https://www.eurobats.org/sites/default/files/documents/pdf/National_Reports/Inf.MoP7_.20-National%20Implementation%20Report%20of%20Germany.pdf
german_bats = c('Barbastella_barbastellus' ,
                'Eptesicus_nilssonii' ,
                'Eptesicus_serotinus' ,
                'Myotis_alcathoe',
                'Myotis_bechsteinii' ,
                'Myotis_brandtii',
                'Myotis_dasycneme' ,
                'Myotis_daubentonii',
                'Myotis_emarginatus' ,
                'Myotis_myotis' ,
                'Myotis_mystacinus' ,
                'Myotis_nattereri',
                'Nyctalus_leisleri' ,
                'Nyctalus_noctula',
                'Pipistrellus_kuhlii' ,
                'Pipistrellus_nathusii',
                'Pipistrellus_pipistrellus' ,
                'Pipistrellus_pygmaeus',
                'Plecotus_auritus' ,
                'Plecotus_austriacus',
                'Rhinolophus_ferrumequinum' ,
                'Rhinolophus_hipposideros',
                'Vespertilio_murinus'
)


# Using data from Conena et al. 2021 for morphological traits
morpho_traits <- fread("Traits/Bats/doi_10.5061_dryad.tmpg4f4xm__v4/conenna_et_al_2021_trait_data_csv.csv", header=TRUE) # https://datadryad.org/stash/dataset/doi:10.5061/dryad.tmpg4f4xm
morpho_traits[, Species := gsub(' ', '_', Binomial2019)]


# Life cycle traits from Wilkinson et al. 2002
lifecycle_traits = data.table(read_excel('Traits/Bats/Life_cycle_Wilkinson2002.xlsx')) #https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1046%2Fj.1474-9728.2002.00020.x&file=ACEL_020_sm_table.doc
gbif = get_gbif_taxonomy(lifecycle_traits$Species)
lifecycle_traits$Lifespan = as.numeric(lifecycle_traits$Lifespan)
lifecycle_traits$Offspring = as.numeric(lifecycle_traits$Offspring)
lifecycle_traits$Female_mass = as.numeric(lifecycle_traits$Female_mass)

# Homogenizing species names
lifecycle_traits[Species == 'Myotis_daubentoni', Species := 'Myotis_daubentonii']
lifecycle_traits[Species == 'Myotis_bechsteini', Species := 'Myotis_bechsteinii']
lifecycle_traits[Species == 'Myotis_brandti', Species := 'Myotis_brandtii']

# Using data from Michela et al. 2014 for generation time
gen_time <- data.table(read_excel("Traits/Bats/doi_10.5061_dryad.gd0m3__v1/Generation\ Lenght\ for\ Mammals.xlsx")) # https://natureconservation.pensoft.net/articles.php?id=1343&element_type=5&display_type=element&element_id=31
gen_time[, Species := gsub(' ', '_', Scientific_name)]
gen_time = gen_time[Species %in% german_bats,]



# ################ #
#### Data merge ####
# ################ #

all_traits = merge.data.table(morpho_traits, lifecycle_traits, by = 'Species', all = T)
all_traits = merge.data.table(all_traits, gen_time[, list(GenLen_m14 = as.numeric(GenerationLength_d),
                                                          Long_m14 =   as.numeric(Max_longevity_d),
                                                          Mass_m14 =   as.numeric(AdultBodyMass_g),
                                                          Species)], by = 'Species', all = T)

all_traits[, Mass.log := log(body.mass)]

# Species present in the exploratories
explo_species = c("Barbastella_barbastellus"  ,"Myotis_myotis"            , "Myotis_nattereri"    ,     "Myotis_sp"                ,
                  "Nyctaloid"                , "Nyctalus_leisleri"       ,  "Nyctalus_noctula"   ,      "Pipistrellus_nathusii"    ,
                  "Pipistrellus_pipistrellus", "Pipistrellus_pygmaeus"   ,  "Plecotus_sp") 

nyctaloid_species = c('Nyctalus_noctula', 'Vespertilio_murinus', 'Eptesicus_serotinus', 'Nyctalus_leisleri', 'Eptesicus_nilssonii') # These are the german bat species that are classified as nyctaloid
myotis_species = gsub(' ', '_', german_bats[grepl('Myotis', german_bats)]) # These are the german Myotis
plecotus_species = gsub(' ', '_',german_bats[grepl('Plecotus', german_bats)]) # These are the german Plecotus

bat_traits = c('Mass.log','GenLen_m14', 'body.mass','Female_mass', 'forearm.length', 'aspect.ratio', 'wing.loading', 'peak.f', 'duration', 'Lifespan', 'Offspring')
Bat_traits = all_traits 
Bat_traits[, RWL := wing.loading/(body.mass^1/3)] # ref: https://www.nature.com/articles/s41598-019-41125-0

# Trait subset to use
bat_traits2 = c('Mass.log', 'Lifespan', 'Offspring')

# For all aggregate species (nyctaloid, Myotis, Plecotus) we average trat data cross corresponding German species
Nyctaloid_traits = Bat_traits[Species %in% nyctaloid_species, lapply(.SD, mean), .SDcols = bat_traits2][, Species := 'Nyctaloid']
Myotis_traits = Bat_traits[Species %in% myotis_species, lapply(.SD, mean, na.rm = T), .SDcols = bat_traits2][, Species := 'Myotis_sp']
Plecotus_traits = Bat_traits[Species %in% plecotus_species, lapply(.SD, mean, na.rm = T), .SDcols = bat_traits2][, Species := 'Plecotus_sp']

# merge aggregate species back with all trait data
Bat_traits_full = rbindlist(list(Bat_traits[, .SD, .SDcols = c(bat_traits2, 'Species')], Nyctaloid_traits, Myotis_traits, Plecotus_traits), use.names=TRUE)


# ############# #
#### Output ####
# ############# #

### Save matched trait dataset
write.csv(x = Bat_traits_full, file = paste(matched_traits_path, 'Matched_bats.csv', sep = ''))

### Check trait coverage
Bats_CC = check_coverage(Bat_traits_full, Abundance_all[ Group_broad == 'Bats',], bat_traits2, 'Species', 'Species')

### Species-level PCA
pca_bats_sp = dudi.pca(mice::complete(mice(Bat_traits_full[Species %in% explo_species, ..bat_traits2])), scannf = FALSE, nf = 2)
gg_bats_sp = fviz_pca(pca_bats_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                       col.ind = "steelblue",
                       fill.ind = "white",
                       col.var = "black")
ggsave(gg_bats_sp, file = 'species_pca_bats.pdf', width = 6, height = 5)


### Calculate CWM 
# weighted by abundance
Abundance_bats = Abundance_all[Species %in% Bat_traits_full$Species,]

CWM_bats = my_cwm(Bat_traits_full, Abundance_bats, bat_traits2, 'Species', 'Species')

fwrite(CWM_bats[, list(
  "bat_mass"= Mass.log,
 # "bat_aspect"= aspect.ratio,
  "bat_lifespan"= Lifespan,
  "bat_offspring"= Offspring   ,
  Plot,
  Year
)], paste(cwm_path,"CWM_Bats.csv", sep = ''))


# Non-weighted community traits, i.e. depends only on species presence/absence

Abundance_bats_presence_absence = Abundance_bats[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)]
Abundance_bats_presence_absence[value>1, value := 1]

CWM_bats_noweight = my_cwm(Bat_traits_full, Abundance_bats_presence_absence, bat_traits2, 'Species', 'Species')

fwrite(CWM_bats_noweight[, list(
  "bat_mass"= Mass.log,
  "bat_lifespan"= Lifespan,
  "bat_offspring"= Offspring   ,
  Plot,
  Year
)], paste(cwm_path, "CWM_bats_noweight.csv", sep = ''))


#######################################
# Check turnover accross LUI gradient #
#######################################

data_lui <- fread("Environment/LUI_input_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = dcast(Abundance_bats[!is.na(Species), list(value = sum(value, na.rm = T)), by = list(Plot, Species)],  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)

