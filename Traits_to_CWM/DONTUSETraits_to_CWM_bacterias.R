# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### In this script, we transform all the microbial data into one CWM matrix. #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
library(factoextra)
library(ggvegan)
library(readr)
library(data.table)
library(readxl)
library(ade4)
library(vegan)
library(factoextra)
library(FactoMineR)
library(fungaltraits)
library(mice)
library(ggcorrplot)
install.packages('fastDummies')
library(fastDummies)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####  New code with bacterial trait database  ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

Bact_data = read_csv('/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Bacteria/Bact_data_24866_25066.csv')
Trait_bacterias <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Bacteria/Bacteria_traits_order.csv")

setDT(Bact_data)
setDT(Trait_bacterias)
Abundance_order <- Bact_data[, list(value = sum(Read_count)), by = list(Plot = Plot_ID, Year = Year, Std_Order = Std_Order)]

## Transform data frames
traits_all <- c(#"metabolism.value", 
  "anaerobic_score.mean", "is.ammoni_pthw.value", "is.nitrif.pthw.value",
  "is.motility.value", "d2_up.mean","d1_up.mean", "volume_est", "doubling_h.mean", "sporulation.value")

Trait_bacterias[, volume_est := d2_up.mean*d1_up.mean]
                
## Check distributions
hist(log(Trait_bacterias$d2_up.mean))
hist(log(Trait_bacterias$d1_up.mean))
hist(log(Trait_bacterias$doubling_h.mean))
# Much better with log !
Trait_bacterias[, c('log_d2_up.mean', 'log_d1_up.mean', 'log_doubling_h.mean', "log_volume") := lapply(.SD, log), .SDcols =  c('d2_up.mean', 'd1_up.mean', 'doubling_h.mean', 'volume_est')]

corr <- round(cor(Trait_bacterias[, c('d2_up.mean', 'doubling_h.mean', 'sporulation.prop', 'd1_up.mean')], use = "pairwise.complete.obs"),1)
p.mat <- cor_pmat(Trait_bacterias[, c('d2_up.mean', 'doubling_h.mean', 'sporulation.prop', 'd1_up.mean')])
ggcorrplot(corr, method = "circle", type = "lower", p.mat = p.mat)


#traits_all = c(traits_all, 'log_d2_up.mean', 'log_d1_up.mean', 'log_doubling_h.mean', "log_volume")
CC_order <- check_coverage(Trait_bacterias, Abundance_order, traits_all, "Std_Order", "Std_Order")
CWM_order <- my_cwm(Trait_bacterias, Abundance_order, traits_all, "Std_Order", "Std_Order")


CWM_aggr = CWM_order[, lapply(.SD, mean), by = Plot, .SDcols = colnames(CWM_order)[3:13]]
colnames(CWM_aggr)[c(3:10, 12, 14, 16)]-1
dudi1 = dudi.pca(CWM_aggr[, c(4:9, 12, 14, 16)-1], col.w = c(1, 1, 1/4, 1/4, 1/4, 1/4, 1, 1, 1), scannf = FALSE, nf = 3)
dudi = dudi.pca(CWM_aggr[, c(3, 4, 5, 16, 12, 14)-1],  scannf = FALSE, nf = 3)

fviz_pca_biplot(dudi)
fviz_pca(dudi1)
fviz_pca_var(dudi)

fviz_pca_biplot(dudi, c(1,3), habillage = as.factor(gsub('[0-9]*', '', CWM_aggr$Plot)))


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# # Stuff to think about:
# #   - no sure which traits to include
# #   - funfun database: very little data on white rot/brown rot; only species name (we have mostly OTUs)
# #   - % sapro, decomposers, etc: % of the total or % of not-unknown ?
# #
# #
# #
#
#
# setwd("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis")
# Abundances = fread("Data/Raw_data/Abundances/Dataset_clean.txt")
#
#
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# #### Import already existing soil microbial indices ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# #microb_soil_prop_2014 <- read_delim("Data/Raw_data/Microbes/20251.txt",
# #                       "\t", escape_double = FALSE, trim_ws = TRUE)
# #microb_soil_prop_2011 <- read_delim("Data/Raw_data/Microbes/20251.txt",
# #                     "\t", escape_double = FALSE, trim_ws = TRUE)
#
#
# #microb_soil_prop = setDT(merge(microb_soil_prop_2011[, c('Year', 'EP_Plot_ID', 'fungi_bacteria', 'Ratio_Cmic_Nmic',"gram_positive",  "gram_negative" )], microb_soil_prop_2014[, c('Year', 'EP_Plot_ID', 'fungi_bacteria', "Ratio_Cmic_Nmic", "gram_positive",  "gram_negative")], by = c('EP_Plot_ID')))
# #soil_microbial_cmw = microb_soil_prop[, list(fungi_bacteria = (fungi_bacteria.x + fungi_bacteria.y)/2,
# #                                             Cmic_Nmic = (Ratio_Cmic_Nmic.x + Ratio_Cmic_Nmic.y)/2,
# #                                             gram_pos_neg = (gram_positive.x/gram_negative.x + gram_positive.y/gram_negative.y)/2
# #                                             ), by = list(Plot = EP_Plot_ID)]
#
#
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ####      Classify fungi based on trophic class      ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
all_fungi = Abundances[Abundances$Group_broad == "soilfungi",]
#
# fungi_groups = all_fungi[, list(prop = sum(value)), by = list(Fun_group_broad, Plot = Plot_bexis)]
# #fungi_groups_cast = dcast(fungi_groups, Plot ~ Fun_group_broad)
#
# #!!!!!!! Maybe replace line below w/o "other":
# #fungi_groups_cast[, c('soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont')] =
# #  fungi_groups_cast[, c('soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont')] / (rowSums(fungi_groups_cast[, c('soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont', 'soilfungi.other')] ))
#
# #soil_microbial_cmw = merge(soil_microbial_cmw, fungi_groups_cast[, c('Plot', 'soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont')], by = 'Plot')
#
# # Calculate the % of white-rot fungi
# # fungal_data = fungal_traits()[, c('trait_fg','species','speciesMatched', 'Genus')]
# # fungal_data[!is.na(fungal_data$trait_fg) & fungal_data$trait_fg == "White Rot",]
#
#
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ####     Add copio/oligo and nitrifying bacterias    ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# #all_bacteria = setDT(Abundances[Abundances$Group_broad == "bacteria.RNA",])
# #all_bacteria$Plot = all_bacteria$Plot_bexis
# #all_bacteria_bygroup = all_bacteria[Year == 2011, sum(value), by = list(Plot =Plot_bexis, Group_fine = Group_fine)]
#
# #bacteria_details_2008 <- setDT(read_delim("Data/Raw_data/Bacteria/19526.txt",
#                     "\t", escape_double = FALSE, trim_ws = TRUE))
#
# #bacteria_details_2008_byphylym = bacteria_details_2008[!grep('W', plotid), sum(Data), by = list(Plot = plotid, Group_fine = Phylum)]
#
# # Let's attribute copio/oligo characteristics
# # Proteobacteria and Bacteroidetes =copiotrophic bacteria.
# # Acidobacteria, Actinobacteria, Gemmatimonadetes, Planctomycetes, and Verrucomicrobia, = oligotrophic bacteria,
# #oli_cop =  c(Bacteroidetes = 'copio', #**
#             Proteobacteria =  'copio',
#             Deltaproteobacteria = 'oligo',
#             Betaproteobacteria  = 'copio',
#             Alphaproteobacteria = 'copio',  #*
#             Gammaproteobacteria   = 'copio', #***
#             Acidobacteria = 'oligo',
#             Gemmatimonadetes = 'oligo',
#             Planctomycetes = 'oligo',
#             Verrucomicrobia = 'oligo')
#
# bacteria_details_2008 = bacteria_details_2008[!grep('W', plotid), functional_type :=
#                         ifelse(Phylum %in% names(oli_cop), oli_cop[Phylum],  oli_cop[Class])]
#
# ratio_oli_cop = bacteria_details_2008[!is.na(functional_type) , list(oligo_copio = sum(Data[functional_type == "oligo"], na.rm = T)/sum(Data[functional_type == "copio"], na.rm =T)),
#                                       , by = list(Plot =plotid)
#                                         ]
#
# #all_bacteria_bygroup = all_bacteria_bygroup[, functional_type := ifelse(Group_fine %in%names(oli_cop), oli_cop[Group_fine], NA)]
# #ratio_oli_cop = all_bacteria_bygroup[,list(oligo_copio = sum(V1[functional_type == "oligo"], na.rm = T)/sum(V1[functional_type == "copio"], na.rm =T))
#                                                                   , by = list(Plot =Plot)]
#
# soil_microbial_cmw = merge(soil_microbial_cmw, ratio_oli_cop, all = T )
#
# #### Calculate % of nitro bacterias
#
# bacteria_details_2008_nitri = bacteria_details_2008[!grep('W', plotid),][
#                                                            Genus %in% c('Nitrosomonas', 'Nitrosococcus', 'Nitrobacter' , 'Nitrospira'),
#                                                            list(nitro = sum(Data)), by = list(Plot = plotid)]
# bacteria_details_2008_nitri = merge(bacteria_details_2008_nitri, bacteria_details_2008[!grep('W', plotid),][ ,
#                                                            list(tot = sum(Data)), by = list(Plot = plotid)])
# bacteria_details_2008_nitri[, prop_nitri := nitro/tot]
#
# soil_microbial_cmw = merge(soil_microbial_cmw, bacteria_details_2008_nitri[, .SD, .SDcols = c('Plot', 'prop_nitri')], all = T )
#
#
# # Check distributions
# soil_microbial_cmw[, sapply(.SD, hist), .SDcols = colnames(soil_microbial_cmw)[2:8]]
#
# soil_microbial_cmw = soil_microbial_cmw[ ,Plot := ifelse(nchar(Plot)== 5, Plot, paste(substr(Plot, 1, 3),substr(Plot, 4, 4), sep = '0'))][order(Plot),]
#


######### CHECK ########@
soil_microbial_cmw[, 2:9] <- mice::complete(mice(soil_microbial_cmw[, 2:9]))


write.csv(soil_microbial_cmw, "Data/CWM_data/CWM_microbes.csv")
