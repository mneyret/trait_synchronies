# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### Get trait species from Exploratories and from gROOT ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

library(traitdataform)

#### First, merge data from the Exploratories ####
AG_BG_traits <- read_delim("~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Below-ground_From_Joana/26587_Above- and below-ground trait data for Grasslands EPs plants at the species level measured in controlled conditions 2017-2018_1.1.7/26587.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
Fine_traits <- read_delim("~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Below-ground_From_Joana/26546_Fine root and mycorrhizal traits of 82 grassland species measured in a greenhouse experiment on sand, 2018_1.1.2/26546.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
t_ag_bg = tpl.get(gsub('_', ' ', unique(AG_BG_traits$scientificName)))
t_fine = tpl.get(gsub('_', ' ', unique(Fine_traits$species)))

roottraitlist <- as.thesaurus(
  biomass_allocation = as.trait("biomass_allocation", expectedUnit = "", valueType = "numeric",
                    identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Body_length"), 
  SRL = as.trait("SRL", expectedUnit = "m/g", valueType = "numeric",
                            identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Antenna_length"),
  AD = as.trait("AD", expectedUnit = "mm", valueType = "numeric",
                              identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  D_first_order = as.trait("D_first_order", expectedUnit = "microm", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  hairlength = as.trait("hairlength", expectedUnit = "microm", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  CF = as.trait("CF", expectedUnit = "%/100", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  col = as.trait("col", expectedUnit = "%/100", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  hairincidence = as.trait("hairincidence", expectedUnit = "%/100", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  SRSA = as.trait("SRSA", expectedUnit = "m2/g", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  LDMC = as.trait("LDMC", expectedUnit = "%/100", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  SLA = as.trait("SLA", expectedUnit = "cm2/g", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length"),
  RTD = as.trait("RTD", expectedUnit = "g/cm3", valueType = "numeric",
                identifier = "http://t-sita.cesab.org/BETSI_vizInfo.jsp?trait=Femur_length")
)

astrait = as.traitdata(Fine_traits, traits = colnames(Fine_traits)[2:13], taxa = 'species')
std_roots = standardize(astrait, thesaurus = roottraitlist)

std_roots    = data.frame(std_roots, stringsAsFactors=FALSE); std_roots$measurementID = as.numeric( std_roots$measurementID) ; std_roots$dataset = 'Bexis_26546'
AG_BG_traits = data.frame(AG_BG_traits, stringsAsFactors=FALSE); AG_BG_traits$measurementID = as.numeric( AG_BG_traits$measurementID); AG_BG_traits$dataset = 'Bexis_26587'


all_roots_B = rbind(std_roots, AG_BG_traits[,-5])
unique(all_roots_B$traitName)

tpl_all_roots_B = tpl.get(gsub('_', ' ', unique(all_roots_B$scientificName)))
all_roots_B$traitName =  as.character(revalue(all_roots_B$traitName, 
                                 c(col = 'Root_mycorrhizal colonization')))
all_roots_B[is.na(all_roots_B$scientificNameStd),]$scientificNameStd = as.character(all_roots_B[is.na(all_roots_B$scientificNameStd),]$scientificName)

unique(all_roots_B$traitName)

#### Try to compare w gRoot ####

# !!!!! Will have to subset for localisation !!!!! #
GRooTAggregateSpeciesVersion <- read_csv("~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Below-ground_from-gROOT/DataFiles/GRooTAggregateSpeciesVersion.csv")
#GRooTGRooTFullVersion <- read_csv("~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/Below-ground_from-gROOT/DataFiles/GRooTFullVersion.csv")

all_species_names$name.synonym

GRooTAggregateSpeciesVersion$scientificName = paste(GRooTAggregateSpeciesVersion$genusTNRS, GRooTAggregateSpeciesVersion$speciesTNRS)
groot_species = unique(GRooTAggregateSpeciesVersion$scientificName)

groot_species = groot_species[!(groot_species %in% c("Anona squamosa"))]
sp_names_std = get_gbif_taxonomy(groot_species)
std_names = sp_names_std$scientificNameStd
names(std_names) = sp_names_std$scientificName

GRooTAggregateSpeciesVersion$scientificNameStd = as.character(std_names[GRooTAggregateSpeciesVersion$scientificName])
GRooTAggregateSpeciesVersion[is.na(GRooTAggregateSpeciesVersion$scientificNameStd),]$scientificNameStd = GRooTAggregateSpeciesVersion[is.na(GRooTAggregateSpeciesVersion$scientificNameStd),]$scientificName

GRooTAggregateSpeciesVersion$dataset = "groot"
GRooTAggregateSpeciesVersion$traitValue = GRooTAggregateSpeciesVersion$medianSpecies
# We will use the following traits: Root_mycorrhizal colonization, Rooting_depth, Root_tissue_density,  Specific_root_length
# Will add : Fine_roots_diameter, hairlength

colnames(GRooTAggregateSpeciesVersion)
colnames(all_roots_B)


all_roots_BE_groot = merge(all_roots_B[ all_roots_B$traitName %in% c('Root_mycorrhizal colonization', 'Rooting_depth', 'Root_tissue_density',  'Specific_root_length'), c('scientificName', 'scientificNameStd', "traitName", "traitValue") ],
                           GRooTAggregateSpeciesVersion[GRooTAggregateSpeciesVersion$traitName %in% c('Root_mycorrhizal colonization', 'Rooting_depth', 'Root_tissue_density',  'Specific_root_length'), c('scientificName', 'scientificNameStd', "traitName", "traitValue") ],
                           by = c('scientificNameStd', "traitName")#,
                           #all = TRUE
                           )

setDT(all_roots_BE_groot)

all_roots_BE_groot[, c('traitValue.x2',  'traitValue.y2' ) := list(as.vector(scale(traitValue.x)), as.vector(scale(traitValue.y))),
                   by = c('traitName')]

all_roots_BE_groot[, cor(traitValue.x2, traitValue.y2, use = 'pairwise.complete.obs' ),
                   by = c('traitName')]
ggplot(all_roots_BE_groot, aes(traitValue.x2, traitValue.y2)) +
  geom_point() + facet_wrap(~traitName) + geom_smooth(method = "lm")

# Not so good !




#### Check how it fits with the abundance dataset ####

Abundances_plants = Abundances[Abundances$Group_broad == 'Plant',]
Abundances_plants$Species = gsub("_aggr.", "", Abundances_plants$Species)
Abundances_plants$Species = gsub("_agg\\.", "", Abundances_plants$Species)
Abundances_plants$Species = gsub("_agg", "", Abundances_plants$Species)
#      - Species with uncertain ID, but resembling a species (cf.) can be attributed this species' traits
Abundances_plants$Species = gsub("_cf_", " ", Abundances_plants$Species)
#      - Homogenise "sp" and "spec"
Abundances_plants$Species = gsub("_spec", "_sp", Abundances_plants$Species)
Abundances_plants$Species = gsub("_sp\\.", "_sp", Abundances_plants$Species)
Abundances_plants$Species = gsub("\\.", "-", Abundances_plants$Species)
Abundances_plants$Species = gsub('_', ' ', Abundances_plants$Species)
Abundances_plants$Species = gsub('Bromus hordeaceus.incl B commutatus.', 'Bromus hordeaceus', Abundances_plants$Species)
Abundances_plants$Species = gsub('Neottia nidus avis', 'Neottia nidus-avis', Abundances_plants$Species)
Abundances_plants$Species = gsub('Impatiens noli tangere', 'Impatiens noli-tangere', Abundances_plants$Species)
Abundances_plants$Species = gsub('Buphtalmum salicifolium', 'Buphthalmum salicifolium', Abundances_plants$Species)
Abundances_plants$Species = gsub('Silene flos.cuculi', 'Silene flos-cuculi', Abundances_plants$Species)
Abundances_plants$Species = gsub('Athyrium filix femina', 'Athyrium filix-femina', Abundances_plants$Species)
Abundances_plants$Species = gsub('Capsella bursa.pastoris', 'Capsella bursa-pastoris', Abundances_plants$Species)
Abundances_plants$Species = gsub('Dryopteris filix mas', 'Dryopteris filix-mas', Abundances_plants$Species)
Abundances_plants$Species = gsub('Mycelis muralis', 'Lactuca muralis', Abundances_plants$Species)
Abundances_plants$Species = gsub('Persicaria lapathifolium', 'Persicaria lapathifolia', Abundances_plants$Species)
Abundances_plants$Species = gsub('Primula elatior-veris', 'Primula elatior', Abundances_plants$Species) # We will average over the two species after getting the traits
#Abundances_plants$Species = c(Abundances_plants$Species, 'Primula veris')
Abundances_plants$Species = gsub('Rubus fruticosus corylifolius', 'Rubus fruticosus', Abundances_plants$Species) # We will average over the two species after getting the traits
#Abundances_plants$Species = c(Abundances_plants$Species, 'Rubus corylifolius')
Abundances_plants$Species = gsub("Viola reichenbachiana riviniana", 'Viola reichenbachiana', Abundances_plants$Species) # We will average over the two species after getting the traits
#Abundances_plants$Species = c(Abundances_plants$Species, "Viola riviniana")
Abundances_plants$Species = gsub("Ribes uva crispa", 'Ribes uva-crispa', Abundances_plants$Species)
Abundances_plants$Species = gsub("Valerianella officinalis", 'Valeriana officinalis', Abundances_plants$Species) # I guess that was probably a typo 

plant_BE = c(unique(Abundances_plants$Species), "Viola riviniana", 'Rubus corylifolius', 'Primula veris')
tpl_BE = tpl.get(plant_BE)


unique(tpl_all_roots_B$name[(tpl_all_roots_B$name %in% tpl_BE$name)])
tpl_all_roots_B$name[is.na(tpl_all_roots_B$name)]
unique(tpl_BE$name[tpl_BE$name %in% tpl_all_roots_B$name])
tpl_BE[is.na(tpl_BE$name),]



#### Check with the new species list ###
New_sp_Ralph <- read_excel("~/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Plants/EPs_Relevee_2008_2016_2019_SpName_Ralph_Bolliger.xlsx")
New_sp_Ralph$NewVersion_Dataset_26106_TPL[!(New_sp_Ralph$NewVersion_Dataset_26106_TPL %in% gsub(' ', '_', tpl_BE$name))]





Abundances_plants$species_check = Abundances_plants$scientificNameStd 
all_roots_B$species_check = all_roots_B$scientificNameStd 
all_roots_B_cast = dcast(all_roots_B, species_check ~ traitName, value.var = 'traitValue', fun.aggregate = mean, na.rm = T)

prop = my_cwm(data.table(all_roots_B_cast), data.table(Abundances_plants), colnames(all_roots_B_cast)[2:19], return_only_cover = TRUE)

test$name
t2$original.search[grepl('sp', t2$original.search)]
t_fine$name
t_ag_bg$name
sort(t_ag_bg[!(t_ag_bg$name %in% t2$name) ,]$original.search)
t_fine[!(t_fine$name %in% t_ag_bg$name) ,]$name

length(test$name)
length(unique(all_roots_B$scientificName))
length(t2$name)
length(unique(Abundances_plants$Species))
