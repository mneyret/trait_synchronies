#### Bats
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


setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/")
Abundances_all = fread("Abundances/Dataset_clean.txt")


# Bats found in Germany, from: https://www.eurobats.org/sites/default/files/documents/pdf/National_Reports/Inf.MoP7_.20-National%20Implementation%20Report%20of%20Germany.pdf
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
setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bats/doi_10.5061_dryad.tmpg4f4xm__v4/")
morpho_traits <- fread("conenna_et_al_2021_trait_data_csv.csv", header=TRUE) # Load file with only species included in analyses
morpho_traits[, Species := gsub(' ', '_', Binomial2019)]

#german_bats %in% morpho_traits$Species --> All species present in dataset

# Life cycle traits from https://onlinelibrary.wiley.com/doi/10.1046/j.1474-9728.2002.00020.x
lifecycle_traits = data.table(read_excel('/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bats/Life_cycle_Wilkinson2002.xlsx'))
gbif = get_gbif_taxonomy(lifecycle_traits$Species)
lifecycle_traits$Lifespan = as.numeric(lifecycle_traits$Lifespan)
lifecycle_traits$Offspring = as.numeric(lifecycle_traits$Offspring)
lifecycle_traits$Female_mass = as.numeric(lifecycle_traits$Female_mass)
ggplot(lifecycle_traits, aes(Offspring, log(Female_mass), color = Species %in% german_bats)) + geom_point() + geom_smooth(method = 'lm')
ggplot(lifecycle_traits, aes(Offspring, Lifespan, color = Species %in% german_bats)) + geom_point() + geom_smooth(method = 'lm')

lifecycle_traits[Species == 'Myotis_daubentoni', Species := 'Myotis_daubentonii']
lifecycle_traits[Species == 'Myotis_bechsteini', Species := 'Myotis_bechsteinii']
lifecycle_traits[Species == 'Myotis_brandti', Species := 'Myotis_brandtii']


# Using data from Michela et al. 2014 for generation time
gen_time <- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bats/doi_10.5061_dryad.gd0m3__v1/Generation\ Lenght\ for\ Mammals.xlsx"))
gen_time[, Species := gsub(' ', '_', Scientific_name)]
gen_time = gen_time[Species %in% german_bats,]

all_traits = merge.data.table(morpho_traits, lifecycle_traits, by = 'Species', all = T)
all_traits = merge.data.table(all_traits, gen_time[, list(GenLen_m14 = as.numeric(GenerationLength_d),
                                                          Long_m14 =   as.numeric(Max_longevity_d),
                                                          Mass_m14 =   as.numeric(AdultBodyMass_g),
                                                          Species)], by = 'Species', all = T)

all_traits[ , Mass.log := log(body.mass)]
explo_species = c("Barbastella_barbastellus"  ,"Myotis_myotis"            , "Myotis_nattereri"    ,     "Myotis_sp"                ,
                  "Nyctaloid"                , "Nyctalus_leisleri"       ,  "Nyctalus_noctula"   ,      "Pipistrellus_nathusii"    ,
                  "Pipistrellus_pipistrellus", "Pipistrellus_pygmaeus"   ,  "Plecotus_sp") 
nyctaloid_species = c('Nyctalus_noctula', 'Vespertilio_murinus', 'Eptesicus_serotinus', 'Nyctalus_leisleri', 'Eptesicus_nilssonii')
myotis_species = gsub(' ', '_', german_bats[grepl('Myotis',german_bats)])
plecotus_species = gsub(' ', '_',german_bats[grepl('Plecotus',german_bats)])

bat_traits = c('Mass.log','GenLen_m14', 'body.mass','Female_mass', 'forearm.length', 'aspect.ratio', 'wing.loading', 'peak.f', 'duration', 'Lifespan', 'Offspring')
Bat_traits = all_traits 
Bat_traits[, RWL := wing.loading/(body.mass^1/3)] # check ref in https://www.nature.com/articles/s41598-019-41125-0
bat_traits2 = c('Mass.log', 'Lifespan', 'Offspring')

Nyctaloid_traits = Bat_traits[Species %in% nyctaloid_species, lapply(.SD, mean), .SDcols = bat_traits2][, Species := 'Nyctaloid']
Myotis_traits = Bat_traits[Species %in% myotis_species, lapply(.SD, mean, na.rm = T), .SDcols = bat_traits2][, Species := 'Myotis_sp']
Plecotus_traits = Bat_traits[Species %in% plecotus_species, lapply(.SD, mean, na.rm = T), .SDcols = bat_traits2][, Species := 'Plecotus_sp']

Bat_traits_full = rbindlist(list(Bat_traits[, .SD, .SDcols = c(bat_traits2, 'Species')], Nyctaloid_traits, Myotis_traits, Plecotus_traits), use.names=TRUE)
Bats_CC = check_coverage(Bat_traits_full, Abundances_all[ Group_broad == 'Bats',], bat_traits2, 'Species', 'Species')

Bat_traits_full[Species %in% explo_species,]
Bats_CC[, -c(1,2)][,lapply(.SD, min)]
Bats_CC[, -c(1,2)][,lapply(.SD, median)]
Bats_CC[, -c(1,2)][,lapply(.SD, max)]

# Species-level
pca_bats = dudi.pca(Bat_traits_full[Species %in% german_bats, ..bat_traits2],  scannf = FALSE, nf = 2)
fviz_pca(pca_bats)


CWM_bats = my_cwm(Bat_traits_full, Abundances_all, bat_traits2, 'Species', 'Species')

pca = dudi.pca(complete(mice(CWM_bats[, lapply(.SD, mean), .SDcols = bat_traits2, by = Plot][, -1])), , scannf = FALSE, nf = 2)
cwm = CWM_bats
tot_pca = fviz_pca(pca)
tot_pca
env_data_lui[, RD := scale(Fertil.) + scale(Mowing)]
quanti.coord <- supcol(pca, data.frame(scale(env_data_lui[Plot %in% cwm$Plot, c( 'LUI')])))$cosup * pca$eig[1:2]^2
tot_pca12 <-fviz_add(tot_pca, quanti.coord, axes = c(1, 2), "arrow", color = "blue", linetype = "solid", repel = T,
                     addlabel = T
)
tot_pca12
cor.test(pca$l1$RS1, env_data_lui[Plot %in% CWM_bats$Plot,]$LUI)

fwrite(CWM_bats[, list(
  "bat_mass"= Mass.log,
 # "bat_aspect"= aspect.ratio,
  "bat_lifespan"= Lifespan,
  "bat_offspring"= Offspring   ,
  Plot,
  Year
)], "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Bats.csv")


CWM_bats_use = merge.data.table(CWM_bats_use, env_data_lui, by = 'Plot')
CWM_bats_use[, (bat_traits) := lapply(.SD, function(x){
  mod = lm(x ~ Clay   +       pH + Soil.depth + Temperature  + Precipitation      +  TWI  +  Grassland + Bulk.density)
  return(residuals(mod))
}), .SDcols = bat_traits]

pca_bats = dudi.pca(CWM_bats_use[, .SD, .SDcols = bat_traits2[-c(4,6)]] , scannf = FALSE, nf = 3)

quanti.coord <- supcol(pca_bats, data.frame(scale(CWM_bats_use[, c('Fertil.', 'LUI', 'Mowing','Grazing')])))$cosup * pca_bats$eig[1:3]^2
pca_plot = fviz_pca(pca_bats, habillage = substr(CWM_bats_use$Plot, 1, 1), axes = c(1,2))
fviz_add(pca_plot, quanti.coord, axes = c(1, 2), "arrow", color = "black", linetype = "solid",# repel = T,
         addlabel = T)

# check : The effect of local land use and loss of forests on bats and nocturnal insects
