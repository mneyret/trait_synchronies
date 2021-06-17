#### Bats
library(splitstackshape)
library(plyr)
library(stringr)


# Using code from Conena et al. 2021
setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bats/doi_10.5061_dryad.tmpg4f4xm__v4/")

trait.dat <- fread("conenna_et_al_2021_trait_data_csv.csv", header=TRUE) # Load file with only species included in analyses
trait.dat$Binomialtest<- trait.dat$Binomial2019
trait.dat$Binomialtest<-gsub(" ", "_", trait.dat$Binomialtest)
trait.dat$Binomial.tree <- revalue(trait.dat$Binomialtest, c( #change old names to new ones (IUCN standards)
  "Austronomus_australis" = "Tadarida_australis",
  "Austronomus_kuboriensis"= "Tadarida_kuboriensis",
  "Baeodon_alleni" = "Rhogeessa_alleni",
  "Baeodon_gracilis"= "Rhogeessa_gracilis",
  "Chaerephon_jobimena"= "Tadarida_jobimena",
  "Dermanura_azteca" = "Dermanura_aztecus",
  "Dermanura_cinerea" = "Dermanura_cinereus",
  "Dermanura_glauca" = "Dermanura_glaucus",
  "Dermanura_gnoma" = "Dermanura_gnomus",
  "Dermanura_rosenbergi"= "Dermanura_rosenbergii",
  "Dermanura_tolteca"= "Dermanura_toltecus",
  "Diclidurus_isabella"= "Diclidurus_isabellus",
  "Gardnerycteris_crenulatum"= "Mimon_crenulatum",
  "Hypsugo_anthonyi"= "Pipistrellus_anthonyi",
  "Hypsugo_joffrei"= "Pipistrellus_joffrei",
  "Hypsugo_kitcheneri"= "Pipistrellus_kitcheneri",
  "Hypsugo_lophurus"= "Pipistrellus_lophurus",
  "Hypsugo_macrotis"= "Pipistrellus_macrotis",
  "Hypsugo_savii"= "Pipistrellus_savii",
  "Hypsugo_vordermanni"= "Pipistrellus_vordermanni",
  "Lophostoma_occidentalis"= "Lophostoma_aequatorialis",
  "Mormopterus_kalinowskii"= "Nyctinomops kalinowskii",
  #"Mormopterus_ridei"=  "Mormopterus_planiceps",
  "Ozimops_loriae"= "Mormopterus_loriae",
  "Paremballonura_atrata"= "Emballonura_atrata",
  "Paremballonura_tiavato"= "Emballonura_tiavato",
  "Perimyotis_subflavus"= "Pipistrellus_subflavus",
  "Pipistrellus_anchietae"= "Hypsugo_anchietae",
  #"Rhinolophus_microglobosus"= "Rhinolophus_stheno",
  "Vampyriscus_bidens"= "Vampyressa_bidens",
  "Vampyriscus_brocki"= "Vampyressa_brocki",
  "Vampyriscus_nymphaea"= "Vampyressa_nymphaea"))
trait.dat[, Species := Binomialtest]

#######################
##### IMPUTATIONS #####
#######################

# Randomly select 100 of 10000 trees from Phylogeny. In this case, here are provided the IDs for the ones used in the analyses.
id <- c(2656,3721,5728,9080,2017,8980,9442,6604,6287,618,2058,1764,6862,3837,7688,4970,7165,9903,3794,7760,9329,
        2117,6503,1253,2666,3852,134,3814,8673,3394,4807,5978,4920,1857,8246,6662,7914,1076,7210,4097,8177,6445,
        7797,5507,5274,7859,233,4750,7288,6894,4753,8569,4359,2435,703,990,3146,5157,6582,4045,9074,2919,4563,
        3304,6468,2564,4754,7612,837,8693,3367,8335,3442,3314,4729,8856,8578,3870,7713,9531,4312,7068,3968,3227,
        7508,2010,7051,1207,2434,1421,2375,584,6364,8682,7716,7898,4510,4062,8030,5990)
phylo.trees<-read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_nexus.trees") # Reference for the tree: Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree:
#Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS biology, 17(12), e3000494.


ImputedDatasets_list<-list(0)

for (i in 1:10) {
  print(i)
  t <- id[[i]]
  tree.n<-phylo.trees[[t]]
  tip.lab <- tree.n[["tip.label"]]
  tip.lab <- data.frame(tip.lab)
  split.lab <- concat.split(tip.lab, "tip.lab", sep = "_", structure = "compact")
  split.lab.bat <- subset(split.lab, tip.lab_4 == "CHIROPTERA")
  split.lab.bat$Binomial.tree <- paste(split.lab.bat$tip.lab_1, split.lab.bat$tip.lab_2, sep = "_")
  split.lab.bat <- merge(split.lab.bat, trait.dat, by = "Binomial.tree")
  tree.n[["tip.label"]] <- str_trim(gsub("[A-Z]{2,}","",tree.n[["tip.label"]]))
  tree.n[["tip.label"]] <- substr(tree.n[["tip.label"]],1,nchar(tree.n[["tip.label"]])-2)
  tree.trimmed<-ape::drop.tip(tree.n, tree.n$tip[!tree.n$tip %in% split.lab.bat$Binomial.tree])
  tree.matrix <- ape::cophenetic.phylo(tree.trimmed)
  pca <- prcomp(tree.matrix, scale = TRUE)
  eigenvect.pca <- unclass(pca$rotation[, 1:5])
  eigenvect.pca <- as.data.frame(eigenvect.pca)
  eigenvect.pca$Binomial <- row.names(eigenvect.pca)
  eigenvect.pca <- eigenvect.pca[order(eigenvect.pca$Binomial),]
  eigenvect.pca1 <- eigenvect.pca %>% dplyr::select(PC1, PC2, PC3, PC4, PC5, Binomial)
  data2 <- merge(trait.dat, eigenvect.pca1, by.x="Binomial.tree", by.y="Binomial", all.x=TRUE)
  data3 <- data2[order(data2$Binomial2019),]
  data <- data3 %>% dplyr::select(Binomial2019,  PC1, PC2, PC3, PC4, PC5,
                                  body.mass, forearm.length, aspect.ratio, wing.loading, peak.f, duration, IUCN.ID.2019)
  data$IUCN.ID.2019 <- as.factor(data$IUCN.ID.2019)
  data$Binomial2019 <- NULL
  data$IUCN.ID.2019 <- NULL
  
  # Imputing
  imp.data <- mice(data, m = 25,  maxit=50, method = c( "pmm", "pmm", "pmm",                                  #"lda","lda","lda",
                                                        "pmm", "pmm", "pmm", "pmm", "pmm", "pmm", "pmm",
                                                        "pmm"))
  complete.imp <- mice::complete(imp.data, action = "long", include = FALSE)
  complete.imp$.id <- as.factor(complete.imp$.id)
  complete.imp$.imp <- as.factor(complete.imp$.imp)
  data <- data3 %>% dplyr::select(Binomial2019,  PC1, PC2, PC3, PC4, PC5,
                                  body.mass, forearm.length, aspect.ratio, wing.loading, peak.f, duration, IUCN.ID.2019)
  data$IUCN.ID.2019 <- as.factor(data$IUCN.ID.2019)
  data$rn <- seq(1,915,1)
  data$rn <- as.factor(data$rn)
  complete.imp <- merge(complete.imp, data[, c("rn", "Binomial2019", "IUCN.ID.2019")], by.x = ".id", by.y = "rn",)
  ImputedDatasets_list[[i]] <- data.frame(complete.imp)
}


whole_data = data.table(Binomial2019 = character(0),
                        body.mass  = numeric(0),
                        forearm.length  = numeric(0),
                        aspect.ratio = numeric(0),
                        wing.loading = numeric(0),
                        peak.f = numeric(0),
                        duration = numeric(0))
for (i in 1: length(ImputedDatasets_list)){
  data = data.table(copy(ImputedDatasets_list[[i]]))
  whole_data = rbind(whole_data,
                     data[, lapply(.SD, mean), .SDcols = c('body.mass' ,'forearm.length' ,'aspect.ratio' ,'wing.loading' ,'peak.f', 'duration'), by = Binomial2019])
  
}

all_trait_data = whole_data[, lapply(.SD, mean), .SDcols = c('body.mass' ,'forearm.length' ,'aspect.ratio' ,'wing.loading' ,'peak.f', 'duration'), by = Binomial2019]


all_trait_data = trait.dat
all_trait_data[, Species := gsub(' ', '_', Binomial2019)]

# Keep only bats found in Germany, from: https://www.eurobats.org/sites/default/files/documents/pdf/National_Reports/Inf.MoP7_.20-National%20Implementation%20Report%20of%20Germany.pdf
german_bats = c('Barbastella barbastellus' ,
                'Eptesicus nilssonii' ,
                'Eptesicus serotinus' ,
                'Myotis alcathoe',
                'Myotis bechsteinii' ,
                'Myotis brandtii',
                'Myotis dasycneme' ,
                'Myotis daubentonii',
                'Myotis emarginatus' ,
                'Myotis myotis' ,
                'Myotis mystacinus' ,
                'Myotis nattereri',
                'Nyctalus leisleri' ,
                'Nyctalus noctula',
                'Pipistrellus kuhlii' ,
                'Pipistrellus nathusii',
                'Pipistrellus pipistrellus' ,
                'Pipistrellus pygmaeus',
                'Plecotus auritus' ,
                'Plecotus austriacus',
                'Rhinolophus ferrumequinum' ,
                'Rhinolophus hipposideros',
                'Vespertilio murinus'
                )

explo_species = c("Barbastella_barbastellus"  ,"Myotis_myotis"            , "Myotis_nattereri"    ,     "Myotis_sp"                ,
                  "Nyctaloid"                , "Nyctalus_leisleri"       ,  "Nyctalus_noctula"   ,      "Pipistrellus_nathusii"    ,
                  "Pipistrellus_pipistrellus", "Pipistrellus_pygmaeus"   ,  "Plecotus_sp") 
nyctaloid_species = c('Nyctalus_noctula', 'Vespertilio_murinus', 'Eptesicus_serotinus', 'Nyctalus_leisleri', 'Eptesicus_nilssonii')
myotis_species = gsub(' ', '_', german_bats[grepl('Myotis',german_bats)])
plecotus_species = gsub(' ', '_',german_bats[grepl('Plecotus',german_bats)])

bat_traits = c('body.mass', 'forearm.length', 'aspect.ratio', 'wing.loading', 'peak.f', 'duration')
Bat_traits = data.table(mice::complete(mice(all_trait_data[, .SD, .SDcols = c('Species', bat_traits)])))
Bat_traits[, RWL := wing.loading/(body.mass^1/3)] # check ref in https://www.nature.com/articles/s41598-019-41125-0
bat_traits2 = c('body.mass', 'forearm.length', 'aspect.ratio', 'wing.loading',
               'peak.f', 'duration', 'RWL')

Nyctaloid_traits = Bat_traits[Species %in% nyctaloid_species, lapply(.SD, mean), .SDcols = bat_traits2][, Species := 'Nyctaloid']
Myotis_traits = Bat_traits[Species %in% myotis_species, lapply(.SD, mean), .SDcols = bat_traits2][, Species := 'Myotis_sp']
Plecotus_traits = Bat_traits[Species %in% plecotus_species, lapply(.SD, mean), .SDcols = bat_traits2][, Species := 'Plecotus_sp']

Bat_traits_full = rbindlist(list(Bat_traits,Nyctaloid_traits,Myotis_traits,Plecotus_traits), use.names=TRUE)

CWM_bats = my_cwm(Bat_traits_full, Abundances_all, bat_traits2, 'Species', 'Species')

fwrite(CWM_bats[, list(
  "bat_mass"= body.mass,
  "bat_forearm"= forearm.length,
  "bat_aspect"= aspect.ratio,
  "bat_WL"= wing.loading,
  "bat_peakf"= peak.f   ,
  "bat_duration"= duration,
  "bat_RWL"= RWL,
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
