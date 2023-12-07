# This script demonstrate the analysis conducted for the manuscript "A fast-slow trait continuum at the level of entire communities" by Neyret et al. 
# Author: Margot Neyret - Please get in touch if you have questions.

# This script takes as input the abundances and species-level traits of plant taxa
# and outputs a CWM matrix averaged for all considered years.

# The code start with trait data already matched from TRY.
# Please check https://www.bexis.uni-jena.de/ddm/data/Showdata/27610 for the procedure and dataset

# ************************** #
#### 3. Species-level PCA ####
# ************************** #

Species_trait_data_clean = fread("Data/Trait_data/Plants/MAIN_Plant_traits_final_all.csv", sep = ";")
# Available here: https://www.bexis.uni-jena.de/ddm/data/Showdata/27610

# Individual PCAs
Data_cast = dcast.data.table(unique(Species_trait_data_clean),
                             scientificName ~ traitName,
                             value.var = 'traitValue')
Data_cast[, c(
  "Fine_roots_diameter",
  "Height",
  "LDMC" ,
  "LeafN",
  "LeafP",
  "Myco_intensity",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "SLA_all"    ,
  "SLA_with_petiole"   ,
  "SSD"     ,
  "Seed_mass"     ,
  "Specific_root_length"
) := lapply(.SD, as.numeric), .SDcols = c(
  "Fine_roots_diameter",
  "Plant_height_vegetative",
  "Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC" ,
  "Leaf_nitrogen_N_content_per_leaf_dry_mass",
  "Leaf_phosphorus_P_content_per_leaf_dry_mass",
  "Mycorrhizal_infection_intensity",
  "Root_tissue_density" ,
  "Root_weight_ratio"  ,
  "Rooting_depth"   ,
  "Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded"    ,
  "Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included"   ,
  "Stem_specific_density_SSD_or_wood_density_stem_dry_mass_per_stem_fresh_volume"     ,
  "Seed_dry_mass"     ,
  "Specific_root_length"
)]


# We need to exclude tree saplings from grassland CWM measures, so here is a list of some tree species found in grasslands
list_trees = c(
  "Acer_campestre",
  "Acer_platanoides" ,
  "Acer_pseudoplatanus",
  "Acer_sp",
  "Alnus_glutinosa",
  "Betula_pendula",
  "Carpinus_betulus",
  "Crataegus_monogyna",
  "Crataegus_sp" ,
  "Fagus_sylvatica",
  "Fraxinus_excelsior"  ,
  "Malus_sylvestris" ,
  "Picea_abies" ,
  "Pinus_sylvestris",
  "Populus_sp" ,
  "Populus_tremula" ,
  "Prunus_avium"   ,
  "Prunus_sp" ,
  "Quercus_petraea"   ,
  "Quercus_robur" ,
  "Quercus_rubra"  ,
  "Quercus_sp"    ,
  "Sorbus_aria" ,
  "Sorbus_aucuparia" ,
  "Sorbus_torminalis",
  "Tilia_platyphyllos"  ,
  "Tilia_sp"       ,
  "Ulmus_minor"      ,
  "Crataegus_laevigata"    ,
  "Abies_alba"     ,
  "Aesculus_hippocastanum",
  "Alnus_incana" ,
  "Larix_decidua" ,
  "Malus_sp"   ,
  "Prunus_padus" ,
  "Prunus_serotina" ,
  "Pseudotsuga_menziesii" ,
  "Pyrus_communis",
  "Pyrus_pyraster" ,
  "Robinia_pseudoacacia",
  "Salix_caprea" ,
  "Salix_cinerea",
  "Salix_sp"  ,
  "Taxus_baccata",
  "Tilia_cordata"  ,
  "Ulmus_glabra" ,
  "Ulmus_laevis"  ,
  "Alnus_sp"    ,
  "Betula_sp"    ,
  "Larix_sp"  ,
  "Populus_nigra"  ,
  "Ulmus_sp" ,
  "Carya_ovata"  ,
  "Pinus_mugo"   ,
  "Pinus_cf_mugo" ,
  'Tree_seedling'
)
list_bushes = c("Hedera_helix" ,
                "Juniperus_communis",
                'Clematis_vitalba',
                'Sambucus_nigra',
                'Corylus_avellana')

pca_AG = dudi.pca(mice::complete(mice(Data_cast[!(scientificName %in% c(list_trees, list_bushes)),][, list(#"Height",
  Root_tissue_density,
  LDMC ,
  LeafN,
  LeafP,
  SLA_with_petiole ,
  SSD ,
  Seed_mass)])), scannf = FALSE, nf = 2)

gg_plants = fviz_pca(pca_AG, title = '', repel = T, geom = 'point', alpha = 0.3,
                     col.ind = "steelblue",
                     fill.ind = "white",
                     col.var = "black")
ggsave(gg_plants, file = 'Results/species_pca_plants.pdf', width = 6, height = 5)

# ************************************ #
#### 3. Calculate coverage and CWM  ####
# ************************************ #

plant_traits =   c("Fine_roots_diameter","Height","LDMC" ,"LeafN","LeafP","Myco_intensity","Root_tissue_density" ,"Root_weight_ratio"  ,"Rooting_depth"   ,"SLA_all"    ,"SLA_with_petiole"   ,"SSD"     ,"Seed_mass"     ,"Specific_root_length")
### Check trait coverage
Plants_CC = check_coverage(Data_cast, Abundances_plants, plant_traits, 'scientificName', 'Species')

### Calculate CWM 
Plants_CWM = my_cwm(Data_cast, Abundances_plants, plant_traits, 'scientificName', 'Species')

# Melt and merge
CWM_CC_plants = merge.data.table(melt.data.table(Plants_CWM[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                 melt.data.table(Plants_CC[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

# Add info

traitUnits = c('mm2/mg','mm2/mg','g/g' ,'mg/g' ,'mg/g' ,'%/100','m','mg','cm3', 'mm','mg/cm3',
  '', 'cm', 'm/g')
traitDescription = c('Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_petiole_included',
                     'Leaf_area_per_leaf_dry_mass_specific_leaf_area_SLA_or_1LMA:_undefined_if_petiole_is_in-_or_excluded',
                     'Leaf_dry_mass_per_leaf_fresh_mass_leaf_dry_matter_content_LDMC' ,'Leaf_nitrogen_N_content_per_leaf_dry_mass' ,
                     'Leaf_phosphorus_P_content_per_leaf_dry_mass' ,'Mycorrhizal_infection_intensity',
                     'Plant vegetative height','Seed dry mass','SSD', 'Diameter of fine roots','Root tissue density',
                     'Ratio of above:below mass', 'Rooting depth', 'Specific root length')

traitDataID = c('Bexis ID 27586','Bexis ID 27586','Bexis ID 27586' ,'Bexis ID 27586' ,'Bexis ID 27586' ,
                'Bexis ID 27586','Bexis ID 27586','Bexis ID 27586','Bexis ID 27586', 
                'Bexis ID 26587','Bexis ID 26587',
                'Bexis ID 26587', 'Bexis ID 26587', 'Bexis ID 26587')

traitRef = c('TRY, doi: 10.1111/gcb.14904','TRY, doi: 10.1111/gcb.14904','TRY, doi: 10.1111/gcb.14904' ,'TRY, doi: 10.1111/gcb.14904' ,'TRY, doi: 10.1111/gcb.14904' ,'TRY, doi: 10.1111/gcb.14904','TRY, doi: 10.1111/gcb.14904','TRY, doi: 10.1111/gcb.14904','TRY, doi: 10.1111/gcb.14904', 'https://doi.org/10.1111/oik.07874,  https://doi.org/10.1111/1365-2745.13862','https://doi.org/10.1111/oik.07874,  https://doi.org/10.1111/1365-2745.13862',
             'https://doi.org/10.1111/oik.07874,  https://doi.org/10.1111/1365-2745.13862', 'https://doi.org/10.1111/oik.07874,  https://doi.org/10.1111/1365-2745.13862', 'https://doi.org/10.1111/oik.07874,  https://doi.org/10.1111/1365-2745.13862')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('SLA_with_petiole','SLA_all','LDMC' ,'LeafN' ,'LeafP' ,'Mycorrhizal_inf_int','Height','Seed_mass','SSD', 'Fine_roots_diameter','Root_tissue_density',
                                                                                       'Root_weight_ratio', 'Rooting_depth', 'Specific_root_length')


CWM_CC_plants = add_info(CWM_CC_plants, traitRef, traitDataID, traitDescription, traitUnits, c('27386 synthesised in 27707'))

fwrite(CWM_CC_plants, "Data/CWM_data/CWM_plants.csv")



# *************************************** #
#### 3. Non-weighted community traits  ####
# *************************************** #
Abundances_plants_presenceabsence = copy(Abundances_plants)
Abundances_plants_presenceabsence[value>1, value :=1]
### Check trait coverage
Plants_CC_noweight = check_coverage(Data_cast, Abundances_plants_presenceabsence, plant_traits, 'scientificName', 'Species')

### Calculate CWM 
Plants_CWM_noweight = my_cwm(Data_cast, Abundances_plants_presenceabsence, plant_traits, 'scientificName', 'Species')

# Melt and merge
CWM_CC_plants_noweight = merge.data.table(melt.data.table(Plants_CWM_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                 melt.data.table(Plants_CC_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))

# Add info
CWM_CC_plants_noweight = add_info(CWM_CC_plants_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('27386 synthesised in 27707'))

fwrite(CWM_CC_plants_noweight, "Data/CWM_data/CWM_plants_noweight.csv")

# ***************** #
#### 5. Turnover ####
# ***************** #

data_lui <- fread("Data/Environment_function_data/LUI_standardized_global.txt") # from https://www.bexis.uni-jena.de/lui/LUICalculation/index; new components, standardised, global, all regions, all years
data_lui = data_lui[Year > 2007 & Year <= 2018, list(LUI = mean(LUI)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
min_lui_plots = data_lui[rank(LUI) <= 10,Plot]
max_lui_plots = data_lui[rank(LUI) > 140,Plot]

comm.test = dcast(Abundances_grasslands[Year < 2019, list(value = sum(value, na.rm = T), Year = 'NA'), by = c('Plot', 'Species')],  Plot~Species, value.var = 'value', fill = 0)
rownames(comm.test)= comm.test$Plot
comm.test = comm.test[,-1]
comm.test[comm.test>0]=1

beta.multi.abund(comm.test)
beta.multi(comm.test)

comm_min_max = matrix(c(colSums(comm.test[min_lui_plots,]),colSums(comm.test[max_lui_plots,])), nrow = 2)
beta.multi.abund(comm_min_max)
