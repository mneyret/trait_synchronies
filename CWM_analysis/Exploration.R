
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
####          Start looking at some results          ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
library(readxl)
library(data.table)
library(ade4)
library(factoextra)
library(ggfortify)
library(autoplot)
library(ggcorrplot)
library(lavaan)
library(semPlot)
library(ghibli)
install.packages("wesanderson")
library(wesanderson)
my_palette <- ghibli_palette("LaputaMedium", direction = -1)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
####  *** Prepare data  ***  ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

#### Parameters ####
setwd("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/")

env_corr <- FALSE
region_corr <- FALSE

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
#### *** Correction for region and/or environment *** ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

corrections = function(CWM_data, Env_data = NA, env_corr = FALSE, region_corr = FALSE, variables = NA, env_variables = NA){
  CWM_data_copy = copy(CWM_data)
  if (env_corr == TRUE & is.na(env_data)){
    stop('Please provide environmental data')
  }
 
   if (is.na(variables)){ 
    variables = colnames(CWM_data_copy)[!(colnames(CWM_data_copy) %in% c('Plot', 'Year'))]
   }
  
  if (is.na(env_variables)){ 
    env_variables = colnames(Env_data)[!(colnames(Env_data) %in% c('Plot', 'Year'))]
  }
  
  print(region_corr)
  if (region_corr == TRUE){
  CWM_data_copy[, (variables) := lapply(.SD, function(x){
    if (region_corr == TRUE){
      mod = lm(x~substr(Plot, 1, 1))
      return(residuals(mod))
    }
  }),
  .SDcols = variables]
    }
  
  if (env_corr == TRUE){
    CWM_data_copy = merge(CWM_data_copy, Env_data[, .SD, .SDcols = c('Plot',env_variables)], by = 'Plot')
    
    CWM_data_copy[, (variables) := lapply(.SD, function(x){
      if (region_corr == TRUE){
        mod = lm(x~., CWM_data_copy[, ..env_variables])
        return(residuals(mod))
      }
    }),
    .SDcols = variables]
  }
  return(CWM_data_copy)
}


# %%%%%%%%%%%%%%%%%%%%%%%%% #
#### *** Trait data *** ####
# %%%%%%%%%%%%%%%%%%%%%%%%% #

# Plants
CWM_plants_raw <- rbind(fread("Data/CWM_data/27608.csv"), fread("Data/CWM_data/27606.csv"))
CWM_plants <- dcast.data.table(CWM_plants_raw, Plot + Year ~ traitName, value.var = "CWM")
CWM_plants_AG <- CWM_plants[grepl("G", Plot), .SD, .SDcols = c("Plot", "Year", "Height", "LDMC", "LeafN", "LeafP", "Seed_mass", "SLA", "SSD")]
CWM_plants_BG <- CWM_plants[grepl("G", Plot), .SD, .SDcols = c("Plot", "Year", "Specific_root_length", "Fine_roots_diameter", "Mycorrhizal_inf_int", "Root_tissue_density", "Root_weight_ratio", "Rooting_depth")]

traits_plantAG <- c("Height", "LDMC", "LeafN", "Seed_mass", "SLA")
traits_plantBG <- c("Specific_root_length", "Fine_roots_diameter", "Mycorrhizal_inf_int", "Root_tissue_density", "Root_weight_ratio", "Rooting_depth")
CWM_plants_AG_agg <- CWM_plants_AG[Year < 2015, lapply(.SD, mean), by = c("Plot"), .SDcols = traits_plantAG]
CWM_plants_BG_agg <- CWM_plants_BG[Year < 2015, lapply(.SD, mean), by = c("Plot"), .SDcols = traits_plantBG]


CWM_plants_AG_agg = corrections(CWM_plants_AG_agg, region_corr = region_corr, env_corr = env_corr)
CWM_plants_BG_agg = corrections(CWM_plants_BG_agg, region_corr = region_corr, env_corr = env_corr)



# Bacterias
CWM_bacterias_GO_complete <- fread("Data/CWM_data/CWM_bacterias_GO.csv")

CWM_bacterias <- CWM_bacterias_GO_complete[, list(
  Plot = Plot,
  Year = Year,
  bact_AOAB = as.numeric(AOAB),
  bact_anaer = log(anaerobic_score.mean),
  bact_nitr = Npathways.value_Nitrifier / Npathways.value_Denitrifier,
  bact_motil = is.motility.value_yes / (is.motility.value_no + is.motility.value_yes),
  bact_d1 = d1_up.mean,
  # d2 = log_d2,
  bact_Size = log_volume,
  bact_Time = log(doubling_h.mean),
  bact_spor = sporulation.value_sporulation / (sporulation.value_sporulation + sporulation.value_not_sporulation),
  bact_gensize = genome_size.mean,
  FB = fungi_bacteria,
  CN = Ratio_Cmic_Nmic,
  fun_spore_size = spore_size,
  AMF = AMF,
  bact_long = as.numeric(est_elongation) # ,
  #pathotroph = soilfungi.pathotroph#,
  # symbionts =  soilfungi.symbiont,
  # decomposer = soilfungi.decomposer#,
  # oligo = oli_copio.value_oligo/(oli_copio.value_oligo +oli_copio.value_copio)#,
  # O.C.ratio = oli_copio.value_oligo/oli_copio.value_copio
)]
traits_bacteria <- colnames(CWM_bacterias)[!(colnames(CWM_bacterias) %in% c("Plot", "Year"))]
CWM_bacterias_agg <- CWM_bacterias[, lapply(.SD, mean, na.rm = T), by = c("Plot"), .SDcols = traits_bacteria]

CWM_bacterias_agg = corrections(CWM_bacterias_agg, region_corr = region_corr, env_corr = env_corr)

# Above-ground Arthropods
CWM_Arthropods_above_herb <- fread("Data/CWM_data/CWM_Arthropods_above_herb.csv")
CWM_Arthropods_above_carni <- fread("Data/CWM_data/CWM_Arthropods_above_carni.csv")

CWM_Arthropods_above_herb <- CWM_Arthropods_above_herb[, list(
  Plot = Plot,
  Year = Year,
  aH_Dispersal = Dispersal_ability,
  aH_Size = logBody_Size,
  aH_Generalism = Feeding_generalism,
  aH_Chewers = Feeding_mode_c / (Feeding_mode_c + Feeding_mode_s),
  aH_Stratum_herb = Stratum_use_simple_h / (Stratum_use_simple_s + Stratum_use_simple_h + Stratum_use_simple_t)
)]
CWM_Arthropods_above_carni <- CWM_Arthropods_above_carni[, list(
  Plot = Plot,
  Year = Year,
  aC_Dispersal = Dispersal_ability,
  aC_Size = logBody_Size,
  #      aC_Generalism =  Feeding_generalism,
  aC_Extraint = Feeding_mode_e / (Feeding_mode_c + Feeding_mode_s + Feeding_mode_e),
  aC_Stratum_herb = Stratum_use_simple_h / (Stratum_use_simple_s + Stratum_use_simple_h + Stratum_use_simple_t)
)]

CWM_Arthropods_above_herb = corrections(CWM_Arthropods_above_herb, region_corr = region_corr, env_corr = env_corr)
CWM_Arthropods_above_carni = corrections(CWM_Arthropods_above_carni, region_corr = region_corr, env_corr = env_corr)

traits_arthropod_herb <- colnames(CWM_Arthropods_above_herb)[!(colnames(CWM_Arthropods_above_herb) %in% c("Year", "Plot"))]
traits_arthropod_carni <- colnames(CWM_Arthropods_above_carni)[!(colnames(CWM_Arthropods_above_carni) %in% c("Year", "Plot"))]

# Mites
CWM_Mites <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_mites_all.csv")
CWM_Mites_soil <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_mites_soil.csv")
CWM_Mites_litter_surface <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_mites_litter_surface.csv")

CWM_Mites <- CWM_Mites[, list(
  Plot = Plot,
  Year = Year,
  mit_Surface = Surface_1,
  mit_Litter = Litter_1,
  mit_Soil_1 = Soil_1,
  mit_Sex = Reproduction_sex,
  mit_Mass = Mass,
  mit_Size = Size,
  mit_Cuticule = CuticuleCalcium_1,
  mit_DaystoAdult = DaystoAdult
)]
CWM_Mites_soil <- CWM_Mites_soil[, list(
  Plot = Plot,
  Year = Year,
  mitS_Sex = Reproduction_sex,
  mitS_Mass = Mass,
  mitS_DaystoAdult = DaystoAdult
)]
CWM_Mites_litter_surface <- CWM_Mites_litter_surface[, list(
  Plot = Plot,
  Year = Year,
  mitL_Sex = Reproduction_sex,
  mitL_Mass = Mass,
  mitL_Size = Size,
  mitL_Cuticule = CuticuleCalcium_1,
  mitL_DaystoAdult = DaystoAdult
)]

CWM_Mites = corrections(CWM_Mites, region_corr = region_corr, env_corr = env_corr)
CWM_Mites_soil = corrections(CWM_Mites_soil, region_corr = region_corr, env_corr = env_corr)
CWM_Mites_litter_surface = corrections(CWM_Mites_litter_surface, region_corr = region_corr, env_corr = env_corr)

traits_mites <- colnames(CWM_Mites)[!(colnames(CWM_Mites) %in% c("Year", "Plot"))]
traits_mites_soil <- colnames(CWM_Mites_soil)[!(colnames(CWM_Mites_soil) %in% c("Year", "Plot"))]
traits_mites_litter_surface <- colnames(CWM_Mites_litter_surface)[!(colnames(CWM_Mites_litter_surface) %in% c("Year", "Plot"))]

# Collembola
CWM_Collembola <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_coll_cwm_all.csv")
CWM_Collembola_epi <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_coll_cwm_epi.csv")
CWM_Collembola_eue <- fread("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_coll_cwm_eue.csv")


CWM_Collembola <- CWM_Collembola[, list(
  Plot = Plot,
  Year = Year,
  col_Size = Size,
  col_Ocelli = Ocelli,
  col_Pigment = Pigment_1,
  col_Furca = Furca_1,
  col_PAO = PAO_1,
  col_PSO = PSO_1,
  col_Asp = Asp_1,
  col_Scales = Scales_1
)]
CWM_Collembola_epi <- CWM_Collembola_epi[, list(
  Plot = Plot,
  Year = Year,
  colL_Size = Size,
  colL_Ocelli = Ocelli,
  colL_Pigment = Pigment_1,
  colL_PAO = PAO_1,
  colL_Asp = Asp_1,
  colL_Scales = Scales_1
)]
CWM_Collembola_eue <- CWM_Collembola_eue[, list(
  Plot = Plot,
  Year = Year,
  colS_Size = Size,
  colS_Furca = Furca_1,
  colS_PAO = PAO_1,
  colS_PSO = PSO_1,
  colS_Asp = Asp_1,
  colS_Scales = Scales_1
)]

CWM_Collembola = corrections(CWM_Collembola, region_corr = region_corr, env_corr = env_corr)
CWM_Collembola_epi = corrections(CWM_Collembola_epi, region_corr = region_corr, env_corr = env_corr)
CWM_Collembola_eue = corrections(CWM_Collembola_eue, region_corr = region_corr, env_corr = env_corr)


traits_collembola <- colnames(CWM_Collembola)[!(colnames(CWM_Collembola) %in% c("Year", "Plot"))]
traits_collembola_epi <- colnames(CWM_Collembola_epi)[!(colnames(CWM_Collembola_epi) %in% c("Year", "Plot"))]
traits_collembola_eue <- colnames(CWM_Collembola_eue)[!(colnames(CWM_Collembola_eue) %in% c("Year", "Plot"))]


# Protists
CWM_protists_cons2 <- fread("Data/CWM_data/CWM_Protists_sec_cons.csv")
CWM_protists_prim_bact <- fread("Data/CWM_data/CWM_Protists_prim_bact.csv")

CWM_protists_cons2 <- CWM_protists_cons2[, list(
  p2_Testate = morpho_code_testate,
  p2_Size = size_estimates,
  Plot = Plot,
  Year = Year
)]
CWM_protists_prim_bact <- CWM_protists_prim_bact[, list(
  p1_Testate = morpho_code_testate / (morpho_code_NA + morpho_code_naked + morpho_code_testate),
  p1_Parasite = morpho_code_endo + (morpho_code_endo + morpho_code_NA + morpho_code_naked + morpho_code_testate),
  p1_Size = size_estimates,
  Plot = Plot,
  Year = Year
)]

traits_protists_cons2 <- colnames(CWM_protists_cons2)[!(colnames(CWM_protists_cons2) %in% c("Year", "Plot"))]
traits_protists_bact <- colnames(CWM_protists_prim_bact)[!(colnames(CWM_protists_prim_bact) %in% c("Year", "Plot"))]

CWM_protists_cons2_agg <- CWM_protists_cons2[grepl("G", Plot), lapply(.SD, mean), .SDcols = traits_protists_cons2, by = Plot]
CWM_protists_prim_bact_agg <- CWM_protists_prim_bact[grepl("G", Plot), lapply(.SD, mean), .SDcols = traits_protists_bact, by = Plot]

CWM_protists_cons2_agg = corrections(CWM_protists_cons2_agg, region_corr = region_corr, env_corr = env_corr)
CWM_protists_prim_bact_agg = corrections(CWM_protists_prim_bact_agg, region_corr = region_corr, env_corr = env_corr)

# %%%%%%%%%%%%%%%%%%%%%%%%% #
#### *** Environment *** ####
# %%%%%%%%%%%%%%%%%%%%%%%%% #
# Environmental data
env_data <- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Environment/env_data.xlsx"))
lui_data <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Environment/LUI_input_data/LUI_standardized_global.txt")
functions_data <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Functions/may2019_grassland_functions.csv")
functions_data[, Bulk.density := mean(c(Bulk.density2011, Bulk.density2014), na.rm = T), by = "Plot"]
functions_data[, Plot := Plotn]


# %%%%%%%%%%%%%%%%%% #
#### *** PCAs *** ####
# %%%%%%%%%%%%%%%%%% #

axes_direction_pca <- function(pca, cwm, traitnames, conditions) {
  # This functions reorient PCA axes so that the direction is reproducible and interpretable
  if (length(conditions) != length(traitnames)) {
    print("Number of axes, traits, and conditions should be equal")
  }
  else {
    cond_positive <- ifelse(gsub(" ", "", conditions) == ">0", TRUE, FALSE)
    I <- length(traitnames)
    for (i in 1:I) {
      if ((pca$li[cwm[get(traitnames[i]) == max(get(traitnames[i])), Plot], i] > 0) != cond_positive[i]) {
        print(paste("Changed direction of axis", i))
        pca$li[, i] <- -pca$li[, i]
        pca$co[, i] <- -pca$co[, i]
      }
    }
  }
  return(pca)
}


#### Environmental PCA ####
env_data <- data.table(mice::complete(mice::mice(env_data)))
lui_data <- lui_data[, c("Grazing", "Mowing", "Fertil.", "LUI") := lapply(.SD, as.numeric), .SDcols = c("G_STD", "M_STD", "F_STD", "LUI")]
lui_data[, Plot := ifelse(nchar(PLOTID) == 5, PLOTID, paste(substr(PLOTID, 1, 3), "0", substr(PLOTID, 4, 4), sep = ""))]
lui_data <- lui_data[Year <= 2015 & Year >= 2008, lapply(.SD, mean), .SDcols = c("Grazing", "Mowing", "Fertil.", "LUI"), by = Plot]
env_data_lui <- merge(env_data[, c("Clay", "pH", "Soil.depth", "Tmean.Annual", "Precip.Annual", "TWI", "Plot")], lui_data, by = "Plot")
env_data_lui[, c("Clay", "pH", "Soil.depth", "Tmean.Annual", "Precip.Annual", "TWI", "Grazing", "Mowing", "Fertil.") :=
  lapply(.SD, scale), .SDcols = c("Clay", "pH", "Soil.depth", "Tmean.Annual", "Precip.Annual", "TWI", "Grazing", "Mowing", "Fertil.")]

colnames(env_data_lui) <- c(
  "Plot", "Clay", "pH", "Soil.depth", "Temperature",
  "Precipitation", "TWI", "Grazing", "Mowing", "Fertil.", "LUI"
)

env_data_lui <- merge(env_data_lui, functions_data[, c("Plot", "Bulk.density")])
cor_env = cor(env_data_lui[,-1])
pmat = cor_pmat(env_data_lui[,-1])
corrpot_env_no_corr = ggcorrplot(cor_env, p.mat = pmat, method = 'circle', 'lower', hc.order	 = T)

env_data_lui = corrections(env_data_lui, env_corr = F, region_corr = region_corr)
cor_env = cor(env_data_lui[,-1])
pmat = cor_pmat(env_data_lui[,-1])
corrpot_env = ggcorrplot(cor_env, p.mat = pmat, method = 'circle', 'lower', hc.order	 = T)


pca_env <- my_dudi_pca(env_data_lui, c("Clay", "pH", "Bulk.density", "Soil.depth", "Temperature", "Precipitation", "TWI"))

pca_env <- axes_direction_pca(pca_env, env_data_lui, c("Temperature", "Fertil."), c("<0", ">0"))

plot_env <- fviz_pca(pca_env,
  habillage = substr(env_data_lui$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Plants - AG"
)
fviz_pca(pca_env,
  habillage = substr(env_data$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Plants - AG", axes = c(2, 3)
)

ggsave(plot = plot_env, paste("Results/plot_env__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)

#-##########################-#
#### BIG PCA - on traits ####
#-##########################-#

# here we just try to put all traits together and see what happens
all_cwm <- merge(CWM_plants_AG_agg[grepl("G", Plot), .SD, .SDcols = c("Plot", traits_plantAG)],
  CWM_plants_BG_agg[grepl("G", Plot), .SD, .SDcols = c("Plot", traits_plantBG)],
  by = "Plot", all = T
)
all_cwm <- merge(all_cwm, CWM_bacterias_agg[, .SD, .SDcols = c("Plot", traits_bacteria)], all = T, by = "Plot")
all_cwm <- merge(all_cwm, CWM_protists_prim_bact_agg[grepl("G", Plot), .SD, .SDcols = c("Plot", traits_protists_bact)], all = T)
all_cwm <- merge(all_cwm, CWM_protists_cons2_agg[grepl("G", Plot), .SD, .SDcols = c("Plot", traits_protists_cons2)], all = T)
all_cwm <- merge(all_cwm, CWM_Arthropods_above_herb[, .SD, .SDcols = c("Plot", traits_arthropod_herb)], all = T)
all_cwm <- merge(all_cwm, CWM_Arthropods_above_carni[, .SD, .SDcols = c("Plot", traits_arthropod_carni)], all = T)
all_cwm <- merge(all_cwm, CWM_Mites_soil[, .SD, .SDcols = c("Plot", traits_mites_soil)], all = T)
all_cwm <- merge(all_cwm, CWM_Mites_litter_surface[, .SD, .SDcols = c("Plot", traits_mites_litter_surface)], all = T)
all_cwm <- merge(all_cwm, CWM_Collembola_epi[, .SD, .SDcols = c("Plot", traits_collembola_epi)], all = T)
all_cwm <- merge(all_cwm, CWM_Collembola_eue[, .SD, .SDcols = c("Plot", traits_collembola_eue)], all = T)

all_cwm <- merge(all_cwm, env_data_lui, all = T)

# Fill missing values
library(mice)
all_cwm <- data.table(mice::complete(mice(all_cwm, m = 20, method = "cart")))

# exclude_traits = c('bact_Width','bact_Length', 'bact_logWidth', 'bact_logLength','bact_logSize', 'bact_logDoubling',
#                   'prot1_log_size', 'prot1_morpho_code_testate', 'prot2_log_size')
exclude_traits <- NA
use_traits <- c(
  traits_plantAG,
  traits_plantBG,
  traits_bacteria,
  traits_protists_bact,
  traits_protists_cons2,
  traits_arthropod_herb,
  traits_arthropod_carni,
  traits_collembola_epi,
  traits_collembola_eue,
  traits_mites_soil,
  traits_mites_litter_surface
) # [-(exclude_traits)]
groups <- c(
  "Plants - AG", "Plants - BG", "Bacteria - fungi", "Protists, bacterivores", "Protists, predators", "Arthropods, herbivores", "Arthropods, carnivores",
  "Collembola epi", "Collembola eue", "Mites soil", "Mites litter + surface"
)
n_traits <- c(
  length(traits_plantAG),
  length(traits_plantBG),
  length(traits_bacteria),
  length(traits_protists_bact),
  length(traits_protists_cons2),
  length(traits_arthropod_herb),
  length(traits_arthropod_carni),
  length(traits_collembola_epi),
  length(traits_collembola_eue),
  length(traits_mites_soil),
  length(traits_mites_litter_surface)
)

arrow_colors <- rep(
  c(groups),
  c(n_traits)
) # [  -exclude_traits]
arrow_colors <- factor(arrow_colors, levels = c(
  "Plants - AG", "Plants - BG", "Bacteria - fungi", "Protists, bacterivores", "Protists, predators", "Arthropods, herbivores", "Arthropods, carnivores",
  "Collembola epi", "Collembola eue", "Mites soil", "Mites litter + surface"
))

com <- scale(data.frame(all_cwm[complete.cases(all_cwm), .SD, .SDcols = c(use_traits, colnames(env_data_lui)[-1])]))
rownames(com) <- all_cwm[complete.cases(all_cwm), ]$Plot

# if (env_corr == 'env_corr'){
#  com2 = data.table(com)
#  com2 =  com2[, lapply(.SD, function(x) {
#    mod <- lm(x ~ Clay + pH + Soil.depth + Precip.Annual + Tmean.Annual + TWI, com2)
#    print(length(residuals(mod)))
#    return(residuals(mod))
#  }), .SDcols = colnames(com2)[1:41]]
#  com2 = data.frame(com2)
#  rownames(com2) = rownames(com)
# }

# cccorr = cor(com2)
# pmat = cor_pmat(com2)
# ggcorrplot(cc, p.mat = pmat)
# a = apply(pmattabme, 2, function(x){length(x[x<0.05])})
# keep_traits = a>15

tot_pca <- dudi.pca(com[, use_traits],
  col.w = rep(1 / n_traits, times = n_traits),
  scannf = FALSE, nf = 3
)

library(RColorBrewer)
my_cols <- c(
  "black",
  "#71475F", "#c09bb1",
  "#536179",
  "#CC3E40", "#812224",
  "black",
  "#F8C9A0", "#F3933F",
  "#8cb369", "#386641",
  "#72ACFD", "#2E4BB2",
  "black", "black"
)

# Without env
tot_pca12 <- fviz_pca(tot_pca,
  col.var = arrow_colors, palette = my_cols, geom = c("point"), habillage = substr(rownames(tot_pca$li), 1, 1),
  repel = T
)
tot_pca13 <- fviz_pca(tot_pca,
  col.var = arrow_colors, palette = my_cols, geom = c("point"), habillage = substr(rownames(tot_pca$li), 1, 1),
  repel = T, axes = c(1, 3)
)

ggsave(plot = tot_pca12, paste("Results/plot_all_traits12__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
ggsave(plot = tot_pca13, paste("Results/plot_all_traits13__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
#ggsave(plot = tot_pca14, paste("Results/plot_all_traits14__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
#ggsave(plot = tot_pca34, paste("Results/plot_all_traits34__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)

# With env as sup var
tot_pca12 <- fviz_pca_var(tot_pca, col.var = arrow_colors, palette = my_cols[-c(1, 7, 14)], geom = c("arrow", "text"), repel = T, alpha = 0.2, max.overlaps = 10)
tot_pca12 <- fviz_add(tot_pca12, quanti.coord, color = "black", geom = "arrow", repel = T, lwd = 2, labelsize = 5, linetype = "longdash")
tot_pca13 <- fviz_pca_var(tot_pca, col.var = arrow_colors, palette = my_cols[-c(1, 7, 14)], geom = c("arrow", "text"), repel = T, axes = c(1, 3))
tot_pca13 <- fviz_add(tot_pca13, axes = c(1, 3), quanti.coord, color = "black", geom = "arrow", repel = T, lwd = 2, labelsize = 5, linetype = "longdash")
tot_pca23 <- fviz_pca_var(tot_pca, col.var = arrow_colors, palette = my_cols[-c(1, 7, 14)], geom = c("arrow", "text"), repel = T, axes = c(2, 3))
tot_pca23 <- fviz_add(tot_pca23, axes = c(2, 3), quanti.coord, color = "black", geom = "arrow", repel = T, lwd = 2, labelsize = 5, linetype = "longdash")
#tot_pca14 <- fviz_pca_var(tot_pca, col.var = arrow_colors, palette = my_cols[-c(1, 7, 14)], geom = c("arrow", "text"), repel = T, axes = c(1, 4))
#tot_pca14 <- fviz_add(tot_pca14, axes = c(1, 4), quanti.coord, color = "black", geom = "arrow", repel = T, lwd = 2, labelsize = 5, linetype = "longdash")
#tot_pca34 <- fviz_pca_var(tot_pca, col.var = arrow_colors, palette = my_cols[-c(1, 7, 14)], geom = c("arrow", "text"), repel = T, axes = c(3, 4))
#tot_pca34 <- fviz_add(tot_pca34, axes = c(3, 4), quanti.coord, color = "black", geom = "arrow", repel = T, lwd = 2, labelsize = 5, linetype = "longdash")

ggsave(plot = tot_pca12, paste("Results/plot_all_env_traits12__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
ggsave(plot = tot_pca13, paste("Results/plot_all_env_traits13__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
ggsave(plot = tot_pca14, paste("Results/plot_all_env_traits14__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)
ggsave(plot = tot_pca34, paste("Results/plot_all_env_traits34__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 10)

correlation_matrix <- cbind(data.table(com), data.table(-tot_pca$li))

corr_p <- cor_pmat(correlation_matrix, sig.level = 0.01)
corr <- round(cor(correlation_matrix), 1)

corr2 <- corr
corr2[corr_p > 0.01] <- 0

corrplot <- ggcorrplot(corr,
  hc.order = TRUE, p.mat = corr_p,
  type = "lower", insig = "blank"
)
ggsave(plot = corrplot, paste("Results/corrplot__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 13, height = 15)

D <- dist(correlation_matrix)
names(D) <- rownames(com)
hc <- hclust(dist(correlation_matrix))
plot(hc)
#-################-#
#### Group PCAs ####
#-################-#

#### __Plants AG ####
pca_plants_AG <- my_dudi_pca(CWM_plants_AG_agg, traits_plantAG)

# Turn the ACP to have always fast/disturbed on the right
pca_plants_AG_stand <- axes_direction_pca(pca_plants_AG, CWM_plants_AG_agg, c("SLA", "Height", "Seed_mass"), c(">0", "<0", ">0"))

# Calculating sup var loadings and rotating them if needed too
quanti.coord <- supcol(pca_plants_AG, scale(env_data_lui[, -1]))$cosup * pca_plants_AG$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_plants_AG$li[, i] == -pca_plants_AG_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}

plot_plants_AG <- fviz_pca(pca_plants_AG_stand,
  habillage = substr(CWM_plants_AG_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Plants - AG", repel = T, geom.var = c("arrow", "text")
)
plot_plants_AG_12 <- fviz_add(plot_plants_AG, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T, addlabel = T)

plot_plants_AG_13 <- fviz_pca(pca_plants_AG_stand,
  habillage = substr(CWM_plants_AG_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Plants - AG", axes = c(1, 3), geom.var = c("arrow", "text")
)
plot_plants_AG_13 <- fviz_add(plot_plants_AG_13, quanti.coord, axes = c(1, 3), "arrow", color = "seagreen", linetype = "solid", repel = T, addlabel = T)

ggsave(plot = plot_plants_AG_12, paste("Results/plot_plants_AG_12__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
ggsave(plot = plot_plants_AG_13, paste("Results/plot_plants_AG_13__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
axis_plant_fast <- pca_plants_AG_stand$li$Axis1
axis_plant_dist <- pca_plants_AG_stand$li$Axis2

cor(pca_plants_AG_stand$li$Axis1, CWM_plants_AG_agg$SLA)# env_data_lui$LUI)

## Same w or w/o env_corr
# Plants AG axis 1: aquisitive/fast -> conservative/slow (= less resources)
# Plants AG axis 2: small -> big (= less disturbance?)
# Plants AG axis 3: increased seed mass

#### __Plants BG ####
pca_plants_BG <- my_dudi_pca(CWM_plants_BG_agg, plantBG_traits)
pca_plants_BG_stand <- axes_direction_pca(pca_plants_BG, CWM_plants_BG_agg, c("Specific_root_length", "Rooting_depth", "Rooting_depth"), c(">0", ">0", ">0"))

quanti.coord <- supcol(pca_plants_BG, scale(env_data_lui[, -1]))$cosup * pca_plants_BG$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_plants_BG$li[, i] == -pca_plants_BG_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}

plot_plants_BG <- fviz_pca(pca_plants_BG_stand,
  habillage = substr(CWM_plants_BG_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Plants - BG", repel = T
)
plot_plants_BG_12 <- fviz_add(plot_plants_BG, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

plot_plants_BG_13 <- fviz_pca(pca_plants_BG_stand,
  habillage = substr(CWM_plants_BG_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4],
  geom = c("point"), xlim = c(-5, 5), title = "Plants - BG", axes = c(1, 3)
)
plot_plants_BG_13 <- fviz_add(plot_plants_BG_13, quanti.coord, axes = c(1, 3), "arrow", color = "seagreen", linetype = "solid", repel = T)

ggsave(plot = plot_plants_BG_12, paste("Results/plot_plants_BG_12__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
ggsave(plot = plot_plants_BG_13, paste("Results/plot_plants_BG_13__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)


axis_plantBG_fast <- pca_plants_BG_stand$li$Axis1
axis_plantBG_dist <- pca_plants_BG_stand$li$Axis2

# Plants BG axis 1: aquisitive/fast -> conservative/slow (= less resources)
### (Plants BG axis 1: collaboration? see Bergmann 2020 -> does not really work any more)
# Plants BG axis 3: collaboration


#### __Bacterias ####

pca_bacterias <- my_dudi_pca(CWM_bacterias_agg, traits_bacteria)

pca_bacterias_stand <- axes_direction_pca(pca_bacterias, CWM_bacterias_agg, c("FB", "bact_spor", "bact_Size"), c("<0", "<0", ">0"))

quanti.coord <- supcol(pca_bacterias, scale(env_data_lui[, -1]))$cosup * pca_bacterias$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_bacterias$li[, i] == -pca_bacterias_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}

plot_bacterias <- fviz_pca_biplot(pca_bacterias_stand,
  habillage = substr(CWM_bacterias_agg$Plot, 1, 1), , palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Soil", repel = T
)
plot_bacterias_12 <- fviz_add(plot_bacterias, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T, addlabel = T)

plot_bacterias_13 <- fviz_pca(pca_bacterias_stand,
  habillage = substr(CWM_bacterias_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Soil", axes = c(1, 3), repel = T
)
plot_bacterias_13 <- fviz_add(plot_bacterias_13, quanti.coord, axes = c(1, 3), "arrow", color = "seagreen", linetype = "solid", repel = T, addlabel = T)

ggsave(plot = plot_bacterias_12, paste("Results/plot_bacterias_12__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
ggsave(plot = plot_bacterias_13, paste("Results/plot_bacterias_13__region",region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)


# Bacterias axis 1: aquisitive/fast -> conservative/slow (= less resources)
# Bacterias axis 2: resistance to disturbance OR size if using LOG)
# Bacterias axis 3: size (only w/o LOG)

axis_bactfungi_fast_wet <- pca_bacterias_stand$li$Axis1

#### __Protists ####
# merging 2 leves bc not enough variables to run PCA
CWM_protists <- merge(CWM_protists_cons2_agg, CWM_protists_prim_bact_agg, by = "Plot")

pca_protists2 <- my_dudi_pca(CWM_protists, c(traits_protists_bact, traits_protists_cons2), NF = 2)
pca_protists2_stand <- axes_direction_pca(pca_protists2, CWM_protists, c("p1_Parasite"), c(">0"))

quanti.coord <- supcol(pca_protists2, scale(env_data_lui[, -1]))$cosup * pca_protists2$eig^2
for (i in 1:ncol(quanti.coord)) {
  if (sum(pca_protists2$li[, i] == -pca_protists2_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_protists2 <- fviz_pca(pca_protists2_stand,
  habillage = substr(CWM_protists_cons2_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Protists 2 consumers", repel = T, c(1,2)
)
plot_protists2_12 <- fviz_add(plot_protists2, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)
ggsave(plot = plot_protists2_12, paste("Results/plot_protists_cons2__region", region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)


axis_protists1_fast <- pca_protists2_stand$li$Axis1
# axis_protists1_fast = pca_protists2_stand$li$Axis1

### axis 1: r->K

# cons bact : bact
# pca_protists_bact = my_dudi_pca(CWM_protists_cons_bact_agg, 2:4)
# pca_protists_bact = axes_direction_pca(pca_protists_bact, CWM_protists_cons_bact_agg, c('Size'), c('>0'))

# quanti.coord <- supcol(pca_protists_bact, scale(env_data_lui[, -1]))$cosup * pca_protists_bact$eig[1:3]^2

# plot_protists_bact = fviz_pca(pca_protists_bact, habillage = substr(CWM_protists_cons_bact_agg$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
#                              geom = c("point"),  xlim = c(-5, 5), title = 'Protists bacterivores', repel =T)
# plot_protists_bact_12 = fviz_add(plot_protists_bact, quanti.coord, axes = c(1, 2), 'arrow', color = 'seagreen', linetype = 'solid', repel = T)

# ggsave(plot = plot_protists_bact_12, 'Results/plot_protists_bact.pdf',width = 5, height = 4)
### axis 1: r->K
### axis 2: mobility

# summary(step(lm(pca_protists_bact$li$Axis1 ~ TotalFertilization + TotalGrazing  + Tmean.Annual + Precip.Annual + TWI, env_data_lui)))
# summary(step(lm(pca_protists_bact$li$Axis2 ~ TotalFertilization + TotalGrazing  + Tmean.Annual + Precip.Annual + TWI, env_data_lui)))


# cons1 : not doable
# pca_protists1 = my_dudi_pca(CWM_protists_cons1_agg, 2:4)
# plot_protists1 = fviz_pca(pca_protists1,# habillage = substr(CWM_protists_cons1_agg$Plot, 1, 1),
#                          geom = c("point"),  xlim = c(-5, 5), title = 'Soil')
# ggsave(plot = plot_protists1, 'Results/plot_protists_cons1.pdf', width = 5, height = 4)


#### __Arthropods ####
### ______AG, herbivores ####
pca_Arthropods_herb <- my_dudi_pca(CWM_Arthropods_above_herb, traits_arthropod_herb)
pca_Arthropods_herb_stand <- axes_direction_pca(pca_Arthropods_herb, CWM_Arthropods_above_herb, c("aH_Generalism", "aH_Size"), c(">0", "<0"))

quanti.coord <- supcol(pca_Arthropods_herb, scale(env_data_lui[Plot %in% CWM_Arthropods_above_herb$Plot, -1]))$cosup * pca_Arthropods_herb$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_Arthropods_herb$li[, i] == -pca_Arthropods_herb_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_Arthropods_herb <- fviz_pca(pca_Arthropods_herb_stand,
  habillage = substr(CWM_Arthropods_above_herb$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Arthropods_herb", repel = T
)
plot_Arthropods_herb_12 <- fviz_add(plot_Arthropods_herb, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

plot_Arthropods_herb_13 <- fviz_pca(pca_Arthropods_herb_stand,
  habillage = substr(CWM_Arthropods_above_herb$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Arthropods_herb", axes = c(1, 3), repel = T
)
plot_Arthropods_herb_13 <- fviz_add(plot_Arthropods_herb_13, quanti.coord, axes = c(1, 3), "arrow", color = "seagreen", linetype = "solid", repel = T)

ggsave(plot = plot_Arthropods_herb_12, paste("Results/plot_Arthropods_herb_12__region", region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
ggsave(plot = plot_Arthropods_herb_13, paste("Results/plot_Arthropods_herb_13__region", region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)

axis_arH_fast <- pca_Arthropods_herb_stand$li$Axis1
axis_arH_disturbed <- pca_Arthropods_herb_stand$li$Axis2

# ______AG, secondary consumers ####
pca_Arthropods_carni <- my_dudi_pca(CWM_Arthropods_above_carni, traits_arthropod_carni, NF = 2)
pca_Arthropods_carni_stand <- axes_direction_pca(pca_Arthropods_carni, CWM_Arthropods_above_carni, c("aC_Size", "aC_Dispersal"), c("<0", ">0"))

quanti.coord <- supcol(pca_Arthropods_carni, scale(env_data_lui[Plot %in% CWM_Arthropods_above_herb$Plot, -1]))$cosup * pca_Arthropods_carni$eig[1:2]^2
for (i in 1:ncol(quanti.coord)) {
  if (sum(pca_Arthropods_carni$li[, i] == -pca_Arthropods_carni_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}

plot_Arthropods_carni <- fviz_pca(pca_Arthropods_carni_stand,
  habillage = substr(CWM_Arthropods_above_carni$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Arthropods_carni", repel = T
)
plot_Arthropods_carni_12 <- fviz_add(plot_Arthropods_carni, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

plot_Arthropods_carni_13 <- fviz_pca(pca_Arthropods_carni_stand,
  habillage = substr(CWM_Arthropods_above_carni$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Arthropods_carni", repel = T
)
plot_Arthropods_carni_13 <- fviz_add(plot_Arthropods_carni_13, quanti.coord, axes = c(1, 3), "arrow", color = "seagreen", linetype = "solid", repel = T)

ggsave(plot = plot_Arthropods_carni_12, paste("Results/plot_Arthropods_carni_12__region", region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)
ggsave(plot = plot_Arthropods_carni_13, paste("Results/plot_Arthropods_carni_13__region", region_corr, '_env', env_corr ,".pdf", sep = ''), width = 5, height = 4)


summary(step(lm(pca_Arthropods_carni$li$Axis1 ~ TotalFertilization + TotalGrazing + Tmean.Annual + Precip.Annual + TWI, env_data_lui[Plot %in% CWM_Arthropods_above_herb$Plot, ])))
summary(step(lm(pca_Arthropods_carni$li$Axis2 ~ TotalFertilization + TotalGrazing + TotalMowing + Tmean.Annual + Precip.Annual + TWI, env_data_lui[Plot %in% CWM_Arthropods_above_herb$Plot, ])))

axis_arC_fast <- pca_Arthropods_carni_stand$li$Axis2
axis_arC_disturbed <- pca_Arthropods_carni_stand$li$Axis1


# ______Mites ####
# Soil mites
pca_Mites_soil <- my_dudi_pca(CWM_Mites_soil, traits_mites_soil)
pca_Mites_soil_stand <- axes_direction_pca(pca_Mites_soil, CWM_Mites_soil, c("mitS_DaystoAdult", "mitS_Mass"), c("<0", ">0"))

quanti.coord <- supcol(pca_Mites_soil, scale(env_data_lui[Plot %in% rownames(pca_Mites_soil$li), -1]))$cosup * pca_Mites_soil$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_Mites_soil$li[, i] == -pca_Mites_soil_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_Mites_soil <- fviz_pca(pca_Mites_soil_stand,
  habillage = substr(CWM_Mites_soil$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Mites", repel = T
)
plot_Mites_soil_12 <- fviz_add(plot_Mites_soil, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

axis_mitS_fast <- pca_Mites_soil_stand$li$Axis1

# Litter and surface mites
pca_Mites_litter <- my_dudi_pca(CWM_Mites_litter_surface, traits_mites_litter_surface)
pca_Mites_litter_stand <- axes_direction_pca(pca_Mites_litter, CWM_Mites_litter_surface, c("mitL_DaystoAdult", "mitL_Mass"), c("<0", ">0"))

quanti.coord <- supcol(pca_Mites_litter, scale(env_data_lui[Plot %in% rownames(pca_Mites_litter$li), -1]))$cosup * pca_Mites_litter$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_Mites_litter$li[, i] == -pca_Mites_litter_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_Mites_litter <- fviz_pca(pca_Mites_litter_stand,
  habillage = substr(CWM_Mites_litter_surface$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Mites", repel = T
)
plot_Mites_litter_12 <- fviz_add(plot_Mites_litter, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

axis_mitL_slow <- pca_Mites_litter_stand$li$Axis1
axis_mitL_repro <- pca_Mites_litter_stand$li$Axis2



# ______Collembola ####
# Surface and litter
pca_coll_epi <- my_dudi_pca(CWM_Collembola_epi, traits_collembola_epi)
pca_coll_epi_stand <- axes_direction_pca(pca_coll_epi, CWM_Collembola_epi, c("colL_Size", "colL_Pigment"), c("<0", ">0"))

quanti.coord <- supcol(pca_coll_epi, scale(env_data_lui[Plot %in% rownames(pca_coll_epi$li), -1]))$cosup * pca_coll_epi$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_coll_epi$li[, i] == -pca_coll_epi_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_coll_epi <- fviz_pca(pca_coll_epi_stand,
  habillage = substr(CWM_Collembola_epi$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Mites", repel = T
)
plot_coll_epi_12 <- fviz_add(plot_coll_epi, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)


# Soil
pca_coll_eue <- my_dudi_pca(CWM_Collembola_eue, traits_collembola_eue)
pca_coll_eue_stand <- axes_direction_pca(pca_coll_eue, CWM_Collembola_eue, c("colS_Furca"), c(">0"))

quanti.coord <- supcol(pca_coll_eue, scale(env_data_lui[Plot %in% rownames(pca_coll_eue$li), -1]))$cosup * pca_coll_eue$eig[1:3]^2
for (i in 1:3) {
  if (sum(pca_coll_eue$li[, i] == -pca_coll_eue_stand$li[, i]) != 0) {
    quanti.coord[, i] <- -quanti.coord[, i]
  }
}
plot_coll_eue <- fviz_pca(pca_coll_eue_stand,
  habillage = substr(CWM_Collembola_eue$Plot, 1, 1), palette = my_palette, col.var = my_palette[4], alpha.ind = 1,
  geom = c("point"), xlim = c(-5, 5), title = "Mites", repel = T
)
plot_coll_eue_12 <- fviz_add(plot_coll_eue, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

axis_coll_eue_dist <- pca_coll_eue_stand$li$Axis1


#-##############################-# 
#### Big PCA - on predefined axes ####
#-#############################
All_axes <- data.table(
  arC_disturbed = axis_arC_disturbed,
  mitS_fast = axis_mitS_fast,
  protists1_fast = axis_protists1_fast,
  arC_fast = axis_arC_fast,
  arH_disturbed = axis_arH_disturbed,
  arH_fast = axis_arH_fast,
  bactfungi_fast_wet = axis_bactfungi_fast_wet,
  coll_eue_disp = axis_coll_eue_disp,
  mitL_repro = axis_mitL_repro,
  mitL_size = axis_mitL_size,
  plant_dist = axis_plant_dist,
  plant_fast = axis_plant_fast,
  plantBG_dist = axis_plantBG_dist,
  plantBG_fast = axis_plantBG_fast
)


all_pca <- dudi.pca(All_axes, scannf = FALSE, nf = 4)
fviz_pca(all_pca)
fviz_pca(all_pca, c(2,3))

quanti.coord <-supcol(all_pca, scale(env_data_lui[,-1]))$cosup * all_pca$eig[1:4]^2

plot_all_pca <- fviz_pca(all_pca)
plot_all_pca <- fviz_add(plot_all_pca, quanti.coord, axes = c(1, 2), "arrow", color = "seagreen", linetype = "solid", repel = T)

plot_all_pca_13 <- fviz_pca(all_pca, c(2,3))
plot_all_pca_13 <- fviz_add(plot_all_pca_13, quanti.coord, axes = c(2, 3), "arrow", color = "seagreen", linetype = "solid", repel = T, c(1,3))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
##### *** Prepare SEM analyses *** ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

library(lavaan)

pcas <- list(
  "pca_plants_AG" = pca_plants_AG,
  "pca_plants_BG" = pca_plants_BG,
  "pca_bacterias" = pca_bacterias,
  "pca_protists_bact" = pca_protists_bact,
  "pca_protists2" = pca_protists2,
  "pca_Arthropods_herb" = pca_Arthropods_herb,
  "pca_Arthropods_carni" = pca_Arthropods_carni
)
#  'pca_Arthropods_ground' = pca_Arthropods_ground,
#   'pca_Arthropods_ground_herb' = pca_Arthropods_ground_herb,
#   'pca_Arthropods_ground_carni' = pca_Arthropods_ground_carni,
#   'pca_Arthropods_above_herb' = pca_Arthropods_above_herb,
#   'pca_Arthropods_above_carni' = pca_Arthropods_above_carni)

dd <- Reduce(
  function(x, y) {
    merge(x, y, all = T)
  },
  lapply(
    names(pcas),
    function(X) {
      x <- get(X)
      x <- x$li[, 1:3]
      colnames(x) <- paste(gsub("pca_", "", X), c(1, 2, 3), sep = "_")
      x[, "Plot"] <- as.character(rownames(x))
      return(x)
    }
  )
)
setDT(dd)

# If using only plants, bacterias
# colnames(dd) = c("Plot", "AG_fast_slow", 'AG_size', 'AG3', 'BG_collab', 'BG_fast_slow', 'AG_BG', 'Bact_fast_slow',
#                 'Bact_disturbance', 'Bact_size', "LUI")
colnames(dd) <- c(
  "Plot",
  "AG_fastslow", "AG_size", "AG_seed",
  "BG_fastslow", "BG_fastslow2", "BG_collab",
  "Bact_fastslow", "Bact_disturbance", "Bact_size",
  "Prot_bact_size", "Prot_bact_2", "Prot_bact_3",
  "Prot_2_size", "Prot_2_2", "Prot_2_3",
  "Arthro_herb_feeding", "Arthro_herb_size", "Arthro_herb3",
  "Arthro_carni_feeding", "Arthro_carni_size", "Arthro_carni3"
)

dd <- merge.data.table(dd, env_data_lui, by = "Plot")
dd[, Region := substr(Plot, 1, 1)]
histogram(dd$TotalMowing)
dd[, Disturbance := scale(TotalMowing) + scale(TotalGrazing)]

corr <- round(cor(dd[, -1], use = "pairwise.complete.obs"), 1)
p.mat <- cor_pmat(dd[, -1])
ggcorr <- ggcorrplot(corr, method = "circle", type = "lower", p.mat = p.mat, colors = c("#CD4F38FF", "white", "#3D4F7DFF"))

ggsave(plot = ggcorr, file = paste("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/axes_correlations__region", region_corr, '_env', env_corr ,".pdf", sep = ''))

# %%%%%%%%%%%%%%%%%%% #
##### *** SEM *** ####
# %%%%%%%%%%%%%%%%%%% #

# https://benwhalley.github.io/just-enough-r/cfa-sem.html

#### Below-ground ####
# Everything together
# Old models
{
  c(
    "Plot", "Env_severity", "env_management", "env3",
    "AG_fastslow", "AG_size", "AG_seed",
    "BG_fastslow", "BG_fastslow2", "BG_collab",
    "Bact_fastslow", "Bact_disturbance", "Bact_size",
    "Prot_bact_size", "Prot_bact_2", "Prot_bact_3",
    "Prot_2_size", "Prot_2_2", "Prot_2_3",
    "Arthro_herb_feeding", "Arthro_herb_size", "Arthro_herb3",
    "Arthro_carni_feeding", "Arthro_carni_size", "Arthro_carni3",
    "LUI"
  )

  "AG_fastslow ~ LUI
            AG_size ~ LUI
            
            BG_fastslow ~ AG_fastslow + LUI
            BG_fastslow2 ~  AG_fast_slow  + LUI
            BG_collab ~ AG_size + LUI
            
            Bact_fast_slow ~  AG_fast_slow + BG_fastslow + LUI
            Bact_disturbance ~ AG_size + BG_collab+ LUI
            Bact_size ~ AG_size + BG_collab+ LUI
            
            Prot_2_size ~ AG_fastslow + AG_size + BG_collab+ Bact_disturbance +Bact_size+ LUI
         
            AG_fast_slow ~~ 0*AG_size
            BG_collab ~~ 0*BG_fastslow2
            BG_fastslow ~~ 0*BG_fastslow2
            Bact_fast_slow ~~ 0*Bact_disturbance+ 0*Bact_size
            Bact_disturbance ~~ 0*Bact_size
"

  test_model <- sem(
    "AG_fastslow ~  LUI + Env_severity
             AG_size ~ Env_severity 
             
             BG_fastslow ~ LUI + AG_size
             BG_fastslow2 ~ AG_size + Env_severity

             Bact_fastslow ~  AG_fastslow + BG_fastslow + Env_severity  
             Bact_disturbance ~ AG_fastslow +BG_fastslow+ Env_severity + LUI + BG_fastslow2

             Prot_bact_size ~ BG_fastslow + Bact_disturbance +  LUI + Env_severity
             Prot_2_size ~   Bact_disturbance +Bact_fastslow

             AG_fastslow ~~ 0*AG_size
             BG_fastslow ~~ 0*BG_fastslow2
             Bact_fastslow ~ 0* Bact_disturbance
             Env_severity ~~ 0*LUI
    

",
    data = dd[-1]
  )
  test_model
  summary(test_model)
  data.table(modificationindices(test_model))[mi > 5, ][order(mi, decreasing = T), ]
  fitmeasures(test_model, c("cfi", "rmsea", "rmsea.ci.upper", "bic"))
  # RMSEA < 0.05, CFI > 0.9, BIC


  summary(test_model)

  #' Fast/slow', 'Size', 'fast/slow', BG_fastslow2, 'fast/slow',  'disturbance', 'Size', 'Size', 'LUI', 'Environment'
  #  x = c(0,        1,       0,          1,              2,                3,            2.5,       2.5,     2,      4)
  #  y = c(3.5,        3.5,       3,         3 ,              2.5,                2.5,             2,    1.5,      4,     4 )/2
  x <- c(0, 0, 1, 1, 2, 2, 3, 3, 0, 1)
  y <- c(1, 0.5, 1, 0.5, 1, 0.5, 0.5, 0.5, 2, 2)
  ly <- matrix(c(x, y), ncol = 2)

  sp <- semPaths(test_model,
    nodeLabels = c("Fast/slow", "Size", "fast/slow", "?", "fast/slow", "disturbance", "Size", "Size", "LUI", "Environment"),
    # nodeLabels = 1:8,
    # style = "lisrel",
    # layout="spring",
    layout = ly,
    what = "std",
    # label.cex=2,
    node.width = 2,
    edge.label.cex = 1 # ,
    # asize = 1.4#,
    # filetype = 'pdf',
    #  filename = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/semplot_BG',
  )

  sp$graphAttributes$Edges$curve

  # By replacing this vector, we add curve to our plot
  sp$graphAttributes$Edges$curve <- c(0, 0, 0.5, 0, 0, 0, 1, -1.5, 0, 1, 0, 0, -0.1, 0.2, 0, 0, 1, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  plot(sp)

  # Then we can plot manipulated p with plot()-function and see the curvature
}
#### New models including latent variables ####

### Below-ground, Fast-slow
model_with_latent <- sem("
             Fertility_int =~  Tmean.Annual+ TWI   +Precip.Annual
             Fertility_added =~   TotalFertilization 

             AG_fastslow ~ Fertility_int + Fertility_added 

             BG_fastslow ~  Fertility_int + Fertility_added 
             Bact_fastslow ~ AG_fastslow + BG_fastslow + Fertility_int 
             
             Prot_bact_size ~ BG_fastslow + BG_fastslow + AG_fastslow +  Fertility_added        + Tmean.Annual  + Precip.Annual
             Prot_2_size ~    BG_fastslow + BG_fastslow  + Prot_bact_size + Fertility_int        +Tmean.Annual

            Tmean.Annual ~~   Precip.Annual
            Tmean.Annual ~~ TWI
                         ",
  data = dd
)

summary(model_with_latent)
data.table(modificationindices(model_with_latent))[mi > 5, ][order(mi, decreasing = T), ]
ly <- matrix(c(
  c(1, 2, 3, 4, 2, 1, 3, 2.5, 1.5, 2, 3),
  c(3, 3, 3, 3, 1, 0, 1.5, 1.5, 0.5, 2, 2)
),
ncol = 2
)
semPaths(model_with_latent,
  what = "std",
  layout = ly,
  filetype = "pdf",
  filename = "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Semplot_latent_BG_fast-slow",
)


### Below-ground, Disturbance
model_with_latent_BG_all <- sem("
             Fertility_int =~  Tmean.Annual + TWI  
             Fertility_added =~   TotalFertilization 

             AG_fastslow ~ Fertility_int + Fertility_added 

             BG_fastslow ~  Fertility_int + Fertility_added + Tmean.Annual
             Bact_fastslow ~ AG_fastslow + BG_fastslow + Fertility_int 
             
             Prot_bact_size ~ BG_fastslow + BG_fastslow + AG_fastslow +  Fertility_added        + Tmean.Annual 
             Prot_2_size ~    BG_fastslow + BG_fastslow  + Prot_bact_size + Fertility_int        +Tmean.Annual

             Tmean.Annual ~~ TWI
            
             TotalFertilization ~~ TotalMowing
             TotalFertilization ~~ TotalGrazing
             TotalMowing ~~ TotalGrazing

             AG_size ~ TotalGrazing  
             BG_fastslow2  ~  TotalGrazing + TotalMowing 
             Bact_disturbance ~ AG_size  + TotalMowing  + BG_fastslow2
             Prot_bact_size ~ Bact_disturbance +  TotalGrazing + TotalMowing  + BG_fastslow2
             Prot_2_size ~   Bact_disturbance   
                         ",
  data = dd
)
summary(model_with_latent_BG_all)
data.table(modificationindices(model_with_latent_BG_all))[mi > 5, ][order(mi, decreasing = T), ]

ly <- matrix(c(
  c(1, 2, 2.5, 3.5, 4.5, 2, 4),
  c(3, 3, 2, 1, 0, 4, 4)
),
ncol = 2
)
semPaths(model_with_latent_BG_all, what = "std") # ,
#      layout = ly,
#       filetype = 'pdf',
#       filename = '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Semplot_latent_BG_size',
#      )



### Below-ground, both pathways
model_with_latent2 <- sem("
             AG_size ~ TotalGrazing  
             BG_fastslow2  ~  TotalGrazing + TotalMowing 
             Bact_disturbance ~ AG_size  + TotalMowing  + BG_fastslow2
             Prot_bact_size ~ Bact_disturbance +  TotalGrazing + TotalMowing  + BG_fastslow2
             Prot_2_size ~   Bact_disturbance   
                         ",
  data = dd
)
summary(model_with_latent2)
data.table(modificationindices(model_with_latent2))[mi > 5, ][order(mi, decreasing = T), ]
ly <- matrix(c(
  c(1, 2, 2.5, 3.5, 4.5, 2, 4),
  c(3, 3, 2, 1, 0, 4, 4)
),
ncol = 2
)
semPaths(model_with_latent2,
  what = "std",
  layout = ly,
  filetype = "pdf",
  filename = "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Semplot_latent_BG_size",
)






#### Above-ground ####
# Full model

#### New models including latent variables ####

### Above-ground, Fast-slow
model_with_latent_AG_fastslow <- sem("
            Fertility_int =~    Tmean.Annual+ TWI  
            Fertility_added =~   TotalFertilization 

            AG_fastslow ~ Fertility_int + Fertility_added 
            #Arthro_herb_size  ~  AG_fastslow  + Fertility_added 
            Arthro_herb_feeding ~  AG_fastslow  + Fertility_added 
            Arthro_carni_feeding ~  Arthro_herb_feeding  + Fertility_int + Fertility_added 
                         ",
  data = dd
)

summary(model_with_latent_AG_fastslow)
data.table(modificationindices(model_with_latent_AG_fastslow))[mi > 5, ][order(mi, decreasing = T), ]
ly <- matrix(c(
  c(1, 2, 3, 1, 3, 2, 1.5, 2.5),
  c(3, 3, 3, 1, 0, 0.5, 2, 2)
),
ncol = 2
)
semPaths(model_with_latent_AG_fastslow,
  what = "std",
  layout = ly,
  filetype = "pdf",
  filename = "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Semplot_latent_AG_fast-slow",
)



### Above-ground, Size
model_with_latent_AG_size <- sem("
            AG_size ~ TotalGrazing + TotalMowing
            Arthro_herb_size  ~   AG_size + TotalGrazing + TotalMowing
            Arthro_carni_size  ~  Arthro_herb_size + AG_size + TotalGrazing + TotalMowing
                         ",
  data = dd
)

summary(model_with_latent_AG_size)
data.table(modificationindices(model_with_latent_AG_size))[mi > 5, ][order(mi, decreasing = T), ]
ly <- matrix(c(
  c(1, 2, 3, 1, 3),
  c(3, 2, 1, 4, 4)
),
ncol = 2
)
semPaths(model_with_latent_AG_size,
  what = "std",
  layout = ly,
  filetype = "pdf",
  filename = "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Results/Semplot_latent_AG_size",
)



### SEM FD ###

FD_Arthropods <- fread(file = "Data/CWM_data/FD_Arthropods.csv")
FD_Arthropods_ground <- fread(file = "Data/CWM_data/FD_Arthropods_ground.csv")
FD_Arthropods_above <- fread(file = "Data/CWM_data/FD_Arthropods_above.csv")
FD_plants_AG <- fread("Data/CWM_data/FD_plants_AG.csv")
FD_plants_BG <- fread("Data/CWM_data/FD_plants_BG.csv")
