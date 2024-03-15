# This script demonstrate the analysis conducted for the manuscript Neyret, M., Le Provost, G., Boesing, A.L. et al. A slow-fast trait continuum at the whole community level in relation to land-use intensification. Nat Commun 15, 1251 (2024).

# This script takes as input the abundances and species-level traits of bat species found 
# in the Exploratories grasslands and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.

# ****************** #
#### 1. Load data ####
# ****************** #

Abundance_bats = Abundance_all[Group_broad == 'Bats',]

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
morpho_traits <- fread("Data/Trait_data/conenna_et_al_2021_trait_data_csv.csv", header=TRUE) # https://datadryad.org/stash/dataset/doi:10.5061/dryad.tmpg4f4xm
morpho_traits[, Species := gsub(' ', '_', Binomial2019)]


# Life cycle traits from Wilkinson et al. 2002
lifecycle_traits = data.table(read_excel('Data/Trait_data/Life_cycle_Wilkinson2002.xlsx')) #https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1046%2Fj.1474-9728.2002.00020.x&file=ACEL_020_sm_table.doc
gbif = get_gbif_taxonomy(lifecycle_traits$Species)
lifecycle_traits$lifespan = as.numeric(lifecycle_traits$lifespan)
lifecycle_traits$Number_offspring = as.numeric(lifecycle_traits$Number_Offspring)
lifecycle_traits$Female_mass = as.numeric(lifecycle_traits$Female_mass)

# Homogenizing species names
lifecycle_traits[Species == 'Myotis_daubentoni', Species := 'Myotis_daubentonii']
lifecycle_traits[Species == 'Myotis_bechsteini', Species := 'Myotis_bechsteinii']
lifecycle_traits[Species == 'Myotis_brandti', Species := 'Myotis_brandtii']

# Using data from Michela et al. 2014 for generation time
gen_time <- data.table(read_excel("Data/Trait_data/Generation\ Lenght\ for\ Mammals.xlsx")) # https://natureconservation.pensoft.net/articles.php?id=1343&element_type=5&display_type=element&element_id=31
gen_time[, Species := gsub(' ', '_', Scientific_name)]
gen_time = gen_time[Species %in% german_bats,]



# ********************* #
#### 2. Trait merge ####
# ******************** #

Bat_traits = merge.data.table(morpho_traits, lifecycle_traits, by = 'Species', all = T)
Bat_traits = merge.data.table(Bat_traits, gen_time[, list(Gen_length = as.numeric(GenerationLength_d),
                                                          Species)], by = 'Species', all = T)

Bat_traits[, logBody_mass := log(body.mass)]
Bat_traits[, c('Peak_freq', 'Lifespan','Number_offspring', 'Forearm_length', 'Aspect_ratio', 'Wing_loading', 'Duration' )  := list(peak.f, Lifespan,Offspring, forearm.length, aspect.ratio, wing.loading, duration)]
Bat_traits[, Number_offspring := as.numeric(Number_offspring)]
Bat_traits[, Lifespan := as.numeric(Lifespan)]

# Species present in the exploratories
explo_species = c("Barbastella_barbastellus"  ,"Myotis_myotis"            , "Myotis_nattereri"    ,     "Myotis_sp"                ,
                  "Nyctaloid"                , "Nyctalus_leisleri"       ,  "Nyctalus_noctula"   ,      "Pipistrellus_nathusii"    ,
                  "Pipistrellus_pipistrellus", "Pipistrellus_pygmaeus"   ,  "Plecotus_sp") 

nyctaloid_species = c('Nyctalus_noctula', 'Vespertilio_murinus', 'Eptesicus_serotinus', 'Nyctalus_leisleri', 'Eptesicus_nilssonii') # These are the german bat species that are classified as nyctaloid
myotis_species = gsub(' ', '_', german_bats[grepl('Myotis', german_bats)]) # These are the german Myotis
plecotus_species = gsub(' ', '_',german_bats[grepl('Plecotus', german_bats)]) # These are the german Plecotus

# Trait subset to use
all_traits = c("Forearm_length", "Aspect_ratio", "Wing_loading", "Peak_freq", "Duration", "Lifespan", "Number_offspring", 'Gen_length','logBody_mass')
bat_traits2 = c('logBody_mass', 'Lifespan', 'Number_offspring')
# For all aggregate species (nyctaloid, Myotis, Plecotus) we average trat data cross corresponding German species
Nyctaloid_traits = Bat_traits[Species %in% nyctaloid_species, lapply(.SD, mean), .SDcols = all_traits][, Species := 'Nyctaloid']
Myotis_traits = Bat_traits[Species %in% myotis_species, lapply(.SD, mean, na.rm = T), .SDcols = all_traits][, Species := 'Myotis_sp']
Plecotus_traits = Bat_traits[Species %in% plecotus_species, lapply(.SD, mean, na.rm = T), .SDcols = all_traits][, Species := 'Plecotus_sp']

# merge aggregate species back with all trait data
Bat_traits_full = rbindlist(list(Bat_traits[, .SD, .SDcols = c(all_traits, 'Species')], Nyctaloid_traits, Myotis_traits, Plecotus_traits), use.names=TRUE)


### Save matched trait dataset

traitUnits = c("mm","unitless","N/m3","kHz",'ms',"year","year-1",'day','g')
traitDescription = c("Distance from elbow to wrist, used as a proxy for body size",
                     "Aspect ratio is defined as in Norberg & Rayner (1987) as B^2/S where B is the wingspan (measured from tip to tip of the fully opened wings) and S is the wing area (as the combined surface of wings, patagium (membrame included between the legs and the tail) and body, excluding the surface of the head protruding from the line of attachment of the wings to the body)",
                     "Wing loading is defined as in Norberg & Rayner (1987) as Mg/S where S is the wing area (as the combined surface of wings, patagium (membrame included between the legs and the tail) and body, excluding the surface of the head protruding from the line of attachment of the wings to the body), and Mg is the weight (mass times gravitational acceleration g)",
                     "Frequency of maximum energy,  obtained for the harmonic with maximum energy based on Monadjem et al., 2010",
                     "Duration of the echolocation call",
                     "Maximum observed lifespan",
                     'Number of offspring per year',
                     "Generation year: average time between two consecutive generations",
                     "Body mass (log)")

traitDataID = c("NA","NA","NA","NA","NA","NA","NA","NA","NA")

traitRef = c("https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1890/08-1494.1",
             "https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1086/368289, https://doi.org/10.1098/rspb.2018.1222",
             "https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1086/368289, https://doi.org/10.1098/rspb.2018.1222",
             "https://doi.org/10.1111/geb.13278, original data is from Collen, A. (2012). The evolution of echolocation in bats: a comparative approach. Doctoral dissertation, University College London",
             "https://doi.org/10.1111/geb.13278, original data is from Collen, A. (2012). The evolution of echolocation in bats: a comparative approach. Doctoral dissertation, University College London",
             "https://doi.org/10.1046/j.1474-9728.2002.00020.x, original data is from Schober & Grimmberger (1997), Neverly (1987), Stebbings & Griffith (1986), Haensel (1994)",
             "https://doi.org/10.1046/j.1474-9728.2002.00020.x, original data is from Schober & Grimmberger (1997), Neverly (1987), Stebbings & Griffith (1986), Haensel (1994)",
             "https://doi.org/10.3897/natureconservation.5.5734, original data is from https://doi.org/10.1890/08-1494.1 and",
             "https://doi.org/10.1111/geb.13278")

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) =  c("Forearm_length", "Aspect_ratio", "Wing_loading", "Peak_freq", "Duration", "Lifespan", "Number_offspring", 'Gen_length','logBody_mass')

Bat_traits_melt = melt.data.table(Bat_traits_full[Species %in% Abundance_bats$Species, .SD, .SDcols = c("Forearm_length", "Aspect_ratio", "Wing_loading", "Peak_freq", "Duration", "Lifespan", "Number_offspring", 'Gen_length','logBody_mass', 'Species')], variable.name = 'traitName', value.name = 'traitValue')
Bat_traits_info = add_info(Bat_traits_melt, traitRef, traitDataID, traitDescription, traitUnits, c('19849, 19850 synthesised in 27707'))

fwrite(Bat_traits_info, "Data/Temporary_data/Bat_traits.csv")

# ************************** #
#### 3. Species-level PCA ####
# ************************** #

### Species-level PCA
pca_bats_sp = dudi.pca(mice::complete(mice(Bat_traits_full[Species %in% explo_species, ..all_traits])), scannf = FALSE, nf = 2)
gg_bats_sp = fviz_pca(pca_bats_sp, title = '', repel = T, geom = 'point', alpha = 0.3,
                       col.ind = "steelblue",
                       fill.ind = "white",
                       col.var = "black")
ggsave(gg_bats_sp, file = 'species_pca_bats.pdf', width = 6, height = 5)


# ********************************** #
#### 4. Community-weighted traits ####
# ********************************** #

### Check trait coverage
Bats_CC = check_coverage(Bat_traits_full, Abundance_bats, all_traits, 'Species', 'Species')

### Calculate CWM 
Bats_CWM = my_cwm(Bat_traits_full, Abundance_bats, all_traits, 'Species', 'Species')

# Melt and merge
CWM_CC_bats = merge.data.table(melt.data.table(Bats_CWM[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                       melt.data.table(Bats_CC[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


# Add info

# "forearm.length"          "aspect.ratio"            "wing.loading"            "peak.freq"                  "duration"               
#  "lifespan"                "Number_Offspring"               "Gen_length" "logBody_mass"    
traitUnits = c("mm","unitless","N/m3","kHz",'ms',"year","year-1",'day','g (log-transformed)')
traitDescription = c("Distance from elbow to wrist, used as a proxy for body size",
                     "Aspect ratio is defined as in Norberg & Rayner (1987) as B^2/S where B is the wingspan (measured from tip to tip of the fully opened wings) and S is the wing area (as the combined surface of wings, patagium (membrame included between the legs and the tail) and body, excluding the surface of the head protruding from the line of attachment of the wings to the body)",
                     "Wing loading is defined as in Norberg & Rayner (1987) as Mg/S where S is the wing area (as the combined surface of wings, patagium (membrame included between the legs and the tail) and body, excluding the surface of the head protruding from the line of attachment of the wings to the body), and Mg is the weight (mass times gravitational acceleration g)",
                     "Frequency of maximum energy,  obtained for the harmonic with maximum energy based on Monadjem et al., 2010",
                     "Duration of the echolocation call",
                     "Maximum observed lifespan",
                     'Number of offspring per year',
                     "Generation year: average time between two consecutive generations",
                     "Body mass")

traitDataID = c("NA","NA","NA","NA","NA","NA","NA","NA","NA")

traitRef = c("https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1890/08-1494.1",
             "https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1086/368289, https://doi.org/10.1098/rspb.2018.1222",
             "https://doi.org/10.1111/geb.13278, original data is from https://doi.org/10.1086/368289, https://doi.org/10.1098/rspb.2018.1222",
             "https://doi.org/10.1111/geb.13278, original data is from Collen, A. (2012). The evolution of echolocation in bats: a comparative approach. Doctoral dissertation, University College London",
             "https://doi.org/10.1111/geb.13278, original data is from Collen, A. (2012). The evolution of echolocation in bats: a comparative approach. Doctoral dissertation, University College London",
             "https://doi.org/10.1046/j.1474-9728.2002.00020.x, original data is from Schober & Grimmberger (1997), Neverly (1987), Stebbings & Griffith (1986), Haensel (1994)",
             "https://doi.org/10.1046/j.1474-9728.2002.00020.x, original data is from Schober & Grimmberger (1997), Neverly (1987), Stebbings & Griffith (1986), Haensel (1994)",
             "https://doi.org/10.3897/natureconservation.5.5734, original data is from https://doi.org/10.1890/08-1494.1 and",
             "https://doi.org/10.1111/geb.13278")

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c("Forearm_length",   "Aspect_ratio"  ,   "Wing_loading"   ,  "Peak_freq" , "Duration"        , "Lifespan"   ,      "Number_offspring", "Gen_length" ,      "logBody_mass")

CWM_CC_bats = add_info(CWM_CC_bats, traitRef, traitDataID, traitDescription, traitUnits, c('19849, 19850 synthesised in 27707'))


fwrite(CWM_CC_bats, "Data/CWM_data/CWM_bats.csv")

# ************************************** #
#### 5. Non-weighted community traits ####
# ************************************** #

Abundance_bats_presence_absence = Abundance_bats[, list(value = sum(value, na.rm = T), Year = 'NA'), by = list(Plot, Species)]
Abundance_bats_presence_absence[value>1, value := 1]

### Check trait coverage
Bats_CC_noweight = check_coverage(Bat_traits_full, Abundance_bats_presence_absence, all_traits, 'Species', 'Species')

### Calculate CWM 
Bats_CWM_noweight = my_cwm(Bat_traits_full, Abundance_bats_presence_absence, all_traits, 'Species', 'Species')
Bats_CWM_noweight$Year = as.character(Bats_CWM_noweight$Year)
Bats_CWM_noweight$Year = 'NA'

# Melt and merge
CWM_CC_bats_noweight = merge.data.table(melt.data.table(Bats_CWM_noweight, id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                        melt.data.table(Bats_CC_noweight, id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))


CWM_CC_bats_noweight = add_info(CWM_CC_bats_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('19849, 19850 synthesised in 27707'))
fwrite(CWM_CC_bats_noweight, "Data/CWM_data/CWM_bats_noweight.csv")


#######################
# Turnover (Table S6) #
#######################

data_lui <- fread("Data/Environment_function_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = dcast(Abundance_bats[!is.na(Species), list(value = sum(value, na.rm = T)), by = list(Plot, Species)],  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]

beta.multi.abund(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)

