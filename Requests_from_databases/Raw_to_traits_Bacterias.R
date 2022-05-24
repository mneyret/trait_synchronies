# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Trying to match BE bacterial data with online database #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
library(readr)
library(data.table)
library(mice)
library(traitdataform)
library(ggcorrplot)
library(dplyr)
library(ade4)
library(factoextra)

setwd('~/Desktop/Research/Senckenberg/')
# ############### #
#### Load data ####
# ############### #

# Traits
condensed_species_NCBI <- setDT(read_delim("~/Data/Traits/Bacteria/From_trait_database/bacteria-archaea-traits-master/output/condensed_species_NCBI.csv",
  ";",
  escape_double = FALSE, col_types = cols(
    carbon_substrates = col_character(), sporulation = col_character(),
    rRNA16S_genes = col_character(),
    cell_shape = col_character(), d1_up = col_number(),
    d2_up = col_number(), doubling_h = col_number(),
    gram_stain = col_character(), motility = col_character(),
    pathways = col_character(), pathways.prop = col_character(), 
    range_tmp = col_character(), range_salinity= col_character(), genome_size =col_number()
  ),
  trim_ws = TRUE
))
# Old dataset (because traitdatafrom does not work anymore)
#Bact_data_old = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Bact_traits_all_old.csv")

#condensed_species_NCBI = merge.data.table(condensed_species_NCBI[, .SD, .SDcol = colnames(condensed_species_NCBI)[!colnames(condensed_species_NCBI) %in% c( 'species','genus'   , 'family','order',  'class',  'phylum', 'superkingdom')]],
#                               Bact_data_old[, .SD, .SDcol =colnames(Bact_data_old)[colnames(Bact_data_old) %in% c( 'species_tax_id','species','genus', 'family','order',  'class',  'phylum', 'superkingdom', "Std_Order", "Std_Genus")]],
#                               by = 'species_tax_id')


corr <- round(cor(condensed_species_NCBI[, c('d1_up', 'd2_up', 'doubling_h')], use = "pairwise.complete.obs"),1)
p.mat <- cor_pmat(condensed_species_NCBI[, c('d1_up', 'd2_up', 'doubling_h')])
ggcorrplot(corr, method = "circle", type = "lower", p.mat = p.mat)

ggplot(condensed_species_NCBI, aes(log(d2_lo), x = motility)) + geom_boxplot()

# Abundance
Bact_data_2011 <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/From_Exploratories/24866.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Bact_data_2014 <- read_delim("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/From_Exploratories/25066.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

Bact_data_2011$Year <- 2011
Bact_data_2014$Year <- 2014
Bact_data <- rbind(Bact_data_2011, Bact_data_2014)
setDT(Bact_data)
Bact_data[, c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Sequence_variant") :=
  tstrsplit(Taxonomy, ", ", type.convert = TRUE, fixed = TRUE)]

# ########################## #
#### Standardize taxonomy ####
# ########################## #

# Function to automatically try to standardize taxonomy (DO ONLY ONCE)
try_standardize_taxo <- function(taxo) {
  taxo_unique <- unique(taxo)
 # print("a")
  gbif_taxo <- get_gbif_taxonomy(taxo_unique[taxo_unique != "" & !is.na(taxo_unique)])
  get_gbif_taxonomy(taxo_unique[taxo_unique != "" & !is.na(taxo_unique)])
  gbif_ScientifiNameSTD <- gbif_taxo$scientificNameStd
 # print("b")
  names(gbif_ScientifiNameSTD) <- gbif_taxo$scientificName
 # print("c")
  taxo_standardized <- gbif_ScientifiNameSTD[taxo]
  taxo_standardized[is.na(taxo_standardized)] <- taxo[is.na(taxo_standardized)]
  return(taxo_standardized)
}

## Standardize Abundance data -> takes some time, do only once!
# Bact_data[!(Order %in% c("bacterium enrichment culture clone Anammox_55",
#                         "bacterium enrichment culture clone SRAO_37",
#                         "bacterium enrichment culture clone Anammox_3")), Std_Order := lapply(.SD, try_standardize_taxo), .SDcols = c("Order")]
# Bact_data[(Order %in% c("bacterium enrichment culture clone Anammox_55",
#                          "bacterium enrichment culture clone SRAO_37",
#                          "bacterium enrichment culture clone Anammox_3")), Std_Order := Order]
# Bact_data[!(Genus %in% c("" ,"bacterium enrichment culture clone auto67_4W",
#                         "bacterium enrichment culture clone auto10_4W",    "Estrella",
#                         "bacterium enrichment culture clone auto112_4W",   "Diplosphaera",
#                         "bacterium enrichment culture clone Anammox_55",   "bacterium enrichment culture clone B30(2011)",
#                         "bacterium enrichment culture clone heteroC52_4W", "bacterium enrichment culture clone SRAO_37",
#                         "Dechlorosoma",                                    "bacterium enrichment culture clone auto79_4W",
#                         "bacterium enrichment culture clone auto73_4W",    "bacterium enrichment culture clone auto8_4W",
#                         "bacterium enrichment culture clone Anammox_3",    "bacterium enrichment culture clone BBMC-4",
#                         "bacterium enrichment culture clone heteroA49_4W", "bacterium enrichment culture clone JCA3",
#                         "bacterium enrichment culture clone SRAO_22" )), Std_Genus := lapply(.SD, try_standardize_taxo), .SDcols = c("Genus")]
# Bact_data[(Genus %in% c("" ,"bacterium enrichment culture clone auto67_4W",
#                          "bacterium enrichment culture clone auto10_4W",    "Estrella",
#                          "bacterium enrichment culture clone auto112_4W",   "Diplosphaera",
#                          "bacterium enrichment culture clone Anammox_55",   "bacterium enrichment culture clone B30(2011)",
#                          "bacterium enrichment culture clone heteroC52_4W", "bacterium enrichment culture clone SRAO_37",
#                          "Dechlorosoma",                                    "bacterium enrichment culture clone auto79_4W",
#                          "bacterium enrichment culture clone auto73_4W",    "bacterium enrichment culture clone auto8_4W",
#                          "bacterium enrichment culture clone Anammox_3",    "bacterium enrichment culture clone BBMC-4",
#                          "bacterium enrichment culture clone heteroA49_4W", "bacterium enrichment culture clone JCA3",
#                          "bacterium enrichment culture clone SRAO_22" )), Std_Genus := Genus]
# 
# write_csv(Bact_data, "/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/Raw_data/Bacteria/Bact_data_24866_25066.csv")
Bact_data = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Bact_data_24866_25066.csv")

## Standardize Abundance data -> takes some time, do only once!
#condensed_species_NCBI[, Std_Order := lapply(.SD, try_standardize_taxo), .SDcols = c("order")]
#condensed_species_NCBI[!(genus %in% c('Agrobacterium', 'Nitrosovibrio','Pimelobacter', 'Methanosaeta')), Std_Genus:= lapply(.SD, try_standardize_taxo), .SDcols = c("genus")]
#condensed_species_NCBI[(genus %in% c('Agrobacterium', 'Nitrosovibrio','Pimelobacter', 'Methanosaeta')), Std_Genus := genus]

#write_csv(condensed_species_NCBI, "/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Bact_traits_all.csv")

Bact_traits_all = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Bact_traits_all.csv")


# How much of order & genus present in the trait dataset?
Bact_data[Genus %in% Bact_traits_all[!is.na(motility),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
Bact_data[Genus %in% Bact_traits_all[!is.na(sporulation),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
Bact_data[Genus %in% Bact_traits_all[!is.na(d1_up),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
Bact_data[Genus %in% Bact_traits_all[!is.na(d1_lo),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]

unique(Bact_data$Genus[!Bact_data$Genus %in% Bact_traits_all$Std_Genus])
unique(Bact_data$Std_Order[Bact_data$Std_Order %in% Bact_traits_all$Std_Order])
unique(Bact_data$Std_Order[!Bact_data$Std_Order %in% Bact_traits_all$Std_Order])


# ########################## #
#### Trait manipulation ####
# ########################## #

## Copio/oligo
# Add % of copio and oligotrophs
copio_groups = c('Actinobacteria', 'Betaproteobacteria', 'Gammaproteobacteria', 'Proteobacteria', 'Bacteroidetes')
oligo_groups = c('Acidobacteria', 'Verrucomicrobia', 'Planctomycetes')

Bact_traits_all[, oli_copio := ifelse(phylum %in% copio_groups | class  %in% copio_groups, 'copio',
                                      ifelse(phylum %in% oligo_groups | class  %in% oligo_groups, 'oligo',
                                             NA))]

traits <- c('oli_copio',"genome_size",  "d2_up", "d1_up","doubling_h")


#### Condensing traits to different taxonomic levels
# The idea here is to use genus data if genus data is available.
# If the data at the genus level is not available, and if the data at the order level is reliable enough, 
# THEN we use order data for the corresponding genera.

find_main_trait <- function(Traitvalue, Traitname, threshold, Std_Genus) {
  Res <- list()
  for (i in 1:ncol(Traitvalue)) {
    traitvalue <- Traitvalue[[i]]
    traitname <- Traitname[i]

    # For qualitative traits
    if (is.character(traitvalue)) {
      print(traitname)
      if (length(unique(c(traitvalue, NA))) > 1) {
        # If there is more than one possibility: we keep the maximum, but only if the proportion is higher than the threshold (e.g. 60%)
        Tab <- table(traitvalue) / sum(table(traitvalue))
        Trait <- names(Tab[Tab == max(Tab)])[1]
        names(Trait) <- paste(traitname, ".value", sep = "")[1]
        Prop <- max(Tab)
        names(Prop) <- paste(traitname, ".prop", sep = "")
        if (Prop < threshold) {
          Trait <- "NA"
        }
      }
      else {
        Prop <- 0
        Trait <- "NA"
      }
      res <- list(Prop, Trait)
      names(res) <- c(paste(traitname, c(".prop", ".value"), sep = ""))
    }


    # For quantitative traits
    if (is.numeric(traitvalue)) {
      Mean <- mean(traitvalue, na.rm = T)
      names(Mean) <- paste(traitname, ".mean", sep = "")
      Sd <- sd(traitvalue, na.rm = T)
      names(Sd) <- paste(traitname, ".sd", sep = "")
      Cv <- Sd / Mean

      names(Cv) <- paste(traitname, ".cv", sep = "")
      if (!is.na(Cv) & Cv > 1) {
        Mean <- as.numeric(NA)
      }
      res <- list(Mean, Sd, Cv)
     # print(unlist(res))
      names(res) <- c(paste(traitname, c(".mean", ".sd", ".cv"), sep = ""))
    }

    Res <- append(Res, res)
  }
  return(Res)
}

Bact_traits_all[, lapply(.SD, function(x){length(x[!is.na(x)])/length(x)})]
### What should be the trophic level for trait aggregation?
#library(cati)


Trait_order <- Bact_traits_all[, find_main_trait(.SD, traits, 0.6),  .SDcols = traits, by = Std_Order]
Trait_genus <- Bact_traits_all[, find_main_trait(.SD, traits, 0.6),  .SDcols = traits, by = Std_Genus]

traits_all <- c(  "d2_up.mean","d1_up.mean", "doubling_h.mean", 
                'oli_copio.value', 'genome_size.mean')

#### Match with abundance data ####
# There are some ?errors? in the Abundance dataset, with genera being associated to two orders. We need to correct this first, avoiding 
# the cases where it is only due to "uncultured" or "unidentified" genera.
genus_order = tapply(Bact_data$Std_Order, Bact_data$Std_Genus, unique)
genus_with_two_orders = genus_order[sapply(genus_order, length)==2 & !(grepl('uncultured', names(genus_order)))]
names_genus_with_two_orders = names(genus_with_two_orders)
corr_genus_with_two_orders = sapply(names_genus_with_two_orders, function(x){Bact_traits_all[Std_Genus == x, unique(Std_Order)]})

Bact_data_corrected = copy(Bact_data)
Bact_data_corrected[Std_Genus %in% names(corr_genus_with_two_orders),]$Std_Order = sapply(Bact_data_corrected[Std_Genus %in% names(corr_genus_with_two_orders),Std_Genus],
      function(x){
        print(corr_genus_with_two_orders[[x]])
      })

# To fit with abundance data, I need to have a name for each genus, incl. only ID to order level.
# For genus names that are not in the trait datasets, I replace the Genus by the Order.
Bact_data_corrected[!(Std_Genus %in% Bact_traits_all$Std_Genus) 
                      , Std_Genus := paste(Std_Order, 'genus')]
Bact_data_corrected[is.na(Std_Genus) & !is.na(Std_Order)
                    , Std_Genus := paste(Std_Order, 'genus')]

# Now I create a new trait table with all genera and corresponding orders
New_traits_genus_order = data.table(unique(Bact_data_corrected[!is.na(Std_Genus), c('Std_Genus', 'Std_Order')]))
New_traits_genus_order[, (traits_all) := NA]
New_traits_genus_order = melt.data.table(New_traits_genus_order, variable.name = "trait", id.var = c('Std_Genus', 'Std_Order'))
New_traits_genus_order = New_traits_genus_order[!grepl('uncultured', Std_Order),]
New_traits_genus_order$value = unlist(as.numeric(New_traits_genus_order$value))

for (i in 1:nrow(New_traits_genus_order)){
 # print('_______')
  genus = New_traits_genus_order$Std_Genus[i]
  order = New_traits_genus_order$Std_Order[i]
  trait = New_traits_genus_order$trait[i]
  if (genus %in% Trait_genus$Std_Genus){
    if (!is.na(Trait_genus[Std_Genus == genus, ..trait])){
      New_traits_genus_order$value[i] = unlist(Trait_genus[Std_Genus == genus, ..trait])
    }
  }
  else{
  if (!is.na(order) & order %in% Trait_order$Std_Order) {
    New_traits_genus_order$value[i] = unlist(Trait_order[Std_Order == order, ..trait])
  }
 }
}


Trait_genus_order = dcast.data.table(New_traits_genus_order, Std_Order + Std_Genus ~ trait, value.var = 'value')
Trait_genus_order = Trait_genus_order[!grepl('uncultivated', Std_Order) & !grepl('unidentified', Std_Order) & !grepl('metagenome', Std_Order) , ]

Trait_genus_order[, c( 'd2_up.mean', 'd1_up.mean','doubling_h.mean', 'genome_size.mean' ) :=
                    lapply(.SD, as.numeric), 
            .SDcols = c( 'd2_up.mean', 'd1_up.mean','doubling_h.mean', 'genome_size.mean'  )]

# Aggregate Abundance data to Genus
Abundance_genus = Bact_data_corrected[, list(value = sum(Read_count, na.rm = T)), by = c( Plot = 'Plot_ID', Std_Genus = 'Std_Genus', Year = 'Year')]
Abundance_genus[, Plot := Plot_ID]


### Transform data frames
# Check distribution
Trait_genus[, est_volume := d1_up.mean*3.14*d1_up.mean*d2_up.mean/4]
Trait_genus_order[, est_volume := d1_up.mean*3.14*d1_up.mean*d2_up.mean/4]
Trait_genus[, est_elongation := d1_up.mean/d2_up.mean]
Trait_genus_order[, est_elongation := d1_up.mean/d2_up.mean]

Trait_genus[, c('log_d1', 'log_d2', 'log_volume', 'log_doubling_H') := lapply(.SD, log), .SDcols = c('d1_up.mean', 'd2_up.mean', 'est_volume',  'doubling_h.mean')]
Trait_genus_order[, c('log_d1', 'log_d2', 'log_volume', 'log_doubling_H') := lapply(.SD, log), .SDcols = c('d1_up.mean', 'd2_up.mean', 'est_volume',  'doubling_h.mean')]

traits_all = unique(c(traits_all, 'log_d1', 'log_d2', 'log_volume', 'log_doubling_H', 'est_elongation'))


# Check species-level
traits_test = Trait_genus_order[Std_Genus %in% unique(Abundance_genus$Std_Genus), list(log_volume, log_doubling_H, est_elongation, genome_size.mean)]
pca_mic = dudi.pca(complete(mice(traits_test)), , scannf = FALSE, nf = 2)
fviz_pca(pca_mic)

# Solution 1: use only genus data
#CC_genus  <- check_coverage(Trait_genus, Abundance_genus, traits_all, "Std_Genus", "Std_Genus")

# Solution 2: use order values when needed
CC_genus_order  <- check_coverage(Trait_genus_order, Abundance_genus, traits_all, "Std_Genus", "Std_Genus")
Trait_genus_order[, lapply(.SD, function(x){length(x[!is.na(x)])})]


write_csv(Trait_genus_order, "~/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Bacteria_traits_genus_order.csv")
write_csv(Trait_genus, "~/Desktop/Research/Senckenberg//Data/Traits/Bacteria/Bacteria_traits_genus.csv")


# Calculate CWM
CWM_bacterias_GO <- my_cwm(Trait_genus_order, Abundance_genus, traits_all, "Std_Genus", "Std_Genus")
CWM_bacterias_GO[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1,3), 0, substr(Plot, 4,4), sep = ''))]


# Calculate CWM without weight

Abundance_genus_presence_absence = Abundance_genus[, list(value = sum(value), Year = 'NA'), by = c( "Plot_ID", "Std_Genus", "Plot", "taxo", "plot_year")]
Abundance_genus_presence_absence[value >1, value := 1]
CWM_bacterias_GO_noweight <- my_cwm(Trait_genus_order, Abundance_genus_presence_absence, traits_all, "Std_Genus", "Std_Genus")
CWM_bacterias_GO_noweight[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1,3), 0, substr(Plot, 4,4), sep = ''))]


# Check correlations
#melt_GO = melt.data.table(CWM_bacterias_GO, id.var = c('Plot', 'Year'))

#melt_al = merge.data.table(melt_GO, melt_G, by = c('Plot', 'Year', 'variable'))
#melt_al[, c('value.x','value.y') := lapply(.SD, function(x){as.numeric(scale(x))}),
#        by = variable, .SDcols = c('value.x','value.y')] 
#ggplot(melt_al, aes(value.x , value.y)) + geom_point() + geom_smooth(method = 'lm') +
#  facet_wrap(~variable)
#
#melt_al[, cor(value.x, value.y), by = variable]


# We also add FB ratio from corresponding datasets
microb_soil_prop_2011 <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Others/20250.txt")
microb_soil_prop_2014 <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Bacteria/Others/20251.txt")

microb_soil_prop = rbindlist(list(microb_soil_prop_2011[, c('Year', 'EP_Plot_ID', 'fungi_bacteria', 'Ratio_Cmic_Nmic',"gram_positive",  "gram_negative" )], 
                                  microb_soil_prop_2014[, c('Year', 'EP_Plot_ID', 'fungi_bacteria', "Ratio_Cmic_Nmic", "gram_positive",  "gram_negative")]))
microb_soil_prop[, Plot := ifelse(nchar(EP_Plot_ID) == 5, EP_Plot_ID, paste(substr(EP_Plot_ID, 1, 3), '0', substr(EP_Plot_ID, 4, 4), sep = ''))]


# And the proportion of each type of fungi
## Raw diversity
allsp <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_diversity_GRL_BEXIS.txt")
allsp$Species = gsub('_$', '', allsp$Species ) # Remove _ if last character
## Species information
fgs <- fread("//Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/210112_EP_species_info_GRL_BEXIS.txt")
fgs$Species = gsub(' $', '', fgs$Species )
fgs$Species = gsub(' ', '_', fgs$Species )
Abundance_all <- merge.data.table(allsp, fgs, by ="Species", all.x=TRUE)
Abundance_all[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), '0', substr(Plot, 4, 4), sep = ''))]
Abundances_fungi = Abundance_all[Group_broad == "soilfungi",]
table(Abundances_fungi[, c('Fun_group_broad', 'Fun_group_fine')])
Abundances_fungi[, sum(value), by = Group_fine]

Prop_fun_group = Abundances_fungi[, list(value = sum(value)), by = c('Plot', 'Fun_group_broad', 'Year')]
Prop_fun_group = dcast.data.table(Prop_fun_group, Year+ Plot~Fun_group_broad, value.var = 'value')
Prop_fun_group_final = Prop_fun_group[, c('total', 'Year') := list(rowSums(.SD), as.numeric(Year)), .SD = c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont')]
Prop_fun_group_final[, c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont') := lapply(.SD, function(x){x/total}), .SD = c('soilfungi.decomposer', 'soilfungi.other', 'soilfungi.pathotroph', 'soilfungi.symbiont')]

microb_soil_prop2 = merge(microb_soil_prop, Prop_fun_group_final, by = c('Plot', 'Year'), all = T)

CWM_bacterias_GO_complete = merge(CWM_bacterias_GO, microb_soil_prop2[, c('Plot','Year','fungi_bacteria', 'Ratio_Cmic_Nmic', 'soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont')], by = c('Plot', 'Year'), all = T)
CWM_bacterias_GO_complete_noweight = merge(CWM_bacterias_GO_noweight, microb_soil_prop2[, c('Plot','Year','fungi_bacteria', 'Ratio_Cmic_Nmic', 'soilfungi.decomposer', 'soilfungi.pathotroph', 'soilfungi.symbiont')], by = c('Plot', 'Year'), all = T)


# Select traits
CWM_bacterias_selection = CWM_bacterias_GO_complete[grepl('G', Plot), list(Plot = Plot,
                                                           Year = Year,
                                                           mic_Bvolume = log_volume,
                                                          mic_Bgenome_size = genome_size.mean/1000,
                                                           mic_FB = fungi_bacteria,
                                                           mic_Fpathotroph = soilfungi.pathotroph,
                                                           mic_Fsymbionts =  soilfungi.symbiont,
                                                           mic_Fdecomposer = soilfungi.decomposer,
                                                          mic_O.C.ratio = oli_copio.value_oligo/oli_copio.value_copio
                                                           )]
CWM_bacterias_selection_noweight = CWM_bacterias_GO_complete_noweight[grepl('G', Plot), list(Plot = Plot,
                                                                           Year = Year,
                                                                           mic_Bvolume = log_volume,
                                                                           mic_Bgenome_size = genome_size.mean/1000,
                                                                           mic_FB = fungi_bacteria,
                                                                           mic_Fpathotroph = soilfungi.pathotroph,
                                                                           mic_Fsymbionts =  soilfungi.symbiont,
                                                                           mic_Fdecomposer = soilfungi.decomposer,
                                                                           mic_O.C.ratio = oli_copio.value_oligo/oli_copio.value_copio
)]
#CWM_bacterias_selection = data.table(mice::complete(mice(CWM_bacterias_selection)))
#CWM_bacterias_selection_noweight = data.table(mice::complete(mice(CWM_bacterias_selection_noweight)))

write_csv(CWM_bacterias_selection, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Microbes.csv')
write_csv(CWM_bacterias_selection_noweight, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_Microbes_noweight.csv')

CWM_bacterias_agg = CWM_bacterias_selection[, lapply(.SD, mean, na.rm = T), .SDcols = colnames(CWM_bacterias_selection)[!colnames(CWM_bacterias_selection) %in% c('Plot', 'Year')], by = 'Plot']
bact_pca = dudi.pca(CWM_bacterias_agg[,-1], scannf = FALSE, nf = 3)
fviz_pca(bact_pca, habillage= substr(CWM_bacterias_agg$Plot, 1, 1))
fviz_pca(bact_pca, axes = c(2, 3), habillage= substr(CWM_bacterias_agg$Plot, 1, 1))

for (r in c('A', 'H', 'S')){
  CWM_bacterias_agg = CWM_bacterias_selection[grepl(r, Plot), lapply(.SD, mean, na.rm = T), .SDcols = colnames(CWM_bacterias_selection)[!colnames(CWM_bacterias_selection) %in% c('Plot', 'Year')], by = 'Plot']
  bact_pca = dudi.pca(CWM_bacterias_agg[,-1], scannf = FALSE, nf = 3)
  plot(fviz_pca(bact_pca, title = r))
}


hist(CWM_bacterias_G_complete$is.motility.value_yes)

#write_csv(CWM_bacterias_G_complete, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_bacterias_G.csv')
write_csv(CWM_bacterias_GO_complete, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_bacterias_GO.csv')
write_csv(CWM_bacterias_GO_complete_noweight, '/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/Data/CWM_data/CWM_bacterias_GO_noweight.csv')

