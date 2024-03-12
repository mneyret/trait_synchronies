# This script demonstrate the analysis conducted for the manuscript "A fast-slow trait continuum at the level of entire communities" by Neyret et al. 
# Author: Margot Neyret - Please get in touch if you have questions.

# This script takes as input the OTU abundances, bacteiral traits and community-level fungal trait information
# and outputs a matched trait dataset, a CWM matrix for all considered years and a species-level PCA.

# Gbif version: Pre-Nov 2022 (taxonomy matching might change afterwards)

# ***************************** #
#### 1. Load and merge data ####
# **************************** #

### Traits from Madin et al.
condensed_species_NCBI <- setDT(read_delim("Data/Trait_data/condensed_species_NCBI.csv",
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

# Abundance
Bact_data_2011 <- read_delim("Data/Abundance_data/24866_1_Dataset/24866_1_data.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Bact_data_2014 <- read_delim("Data/Abundance_data/25066_1_Dataset/25066_1_data.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

Bact_data_2011$Year <- 2011
Bact_data_2014$Year <- 2014
Bact_data <- rbind(Bact_data_2011, Bact_data_2014)
setDT(Bact_data)
Bact_data[, c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Sequence_variant") :=
  tstrsplit(Taxonomy, ", ", type.convert = TRUE, fixed = TRUE)]


## Standardize taxonomy ##

## Standardize Abundance data -> takes some time, do only once! 

Bact_data[Order != '', Std_Order := get_gbifid_(Order)[[1]]$order[1], by = Order ]
Bact_data[Genus != '', Std_Genus := get_gbifid_(Genus)[[1]]$genus[1], by = Genus ]

#write_csv(Bact_data, "Data/Abundance_data/Bact_data_24866__25066.csv")
#Bact_data = fread("Data/Temporary_data/Bact_data_24866__25066.csv")

## Standardize trait data -> takes some time, do only once!

condensed_species_NCBI[order != '', Std_Order := get_gbifid_(order)[[1]]$order[1], by = order ]
condensed_species_NCBI[genus != '', Std_Genus := get_gbifid_(genus)[[1]]$genus[1], by = genus ]

#write_csv(condensed_species_NCBI, "Traits/Bacteria/Temp_data/Bact_traits_all.csv")
#Bact_traits_all = fread("Data/Temporary_data/Bact_traits_all.csv")

# Check how many  orders & genera present in the trait dataset?
# Bact_data[Genus %in% Bact_traits_all[!is.na(motility),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
# Bact_data[Genus %in% Bact_traits_all[!is.na(sporulation),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
# Bact_data[Genus %in% Bact_traits_all[!is.na(d1_up),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
# Bact_data[Genus %in% Bact_traits_all[!is.na(d1_lo),]$Std_Genus, sum(Read_count)]/Bact_data[, sum(Read_count)]
#
# unique(Bact_data$Genus[!Bact_data$Genus %in% Bact_traits_all$Std_Genus])
# unique(Bact_data$Std_Order[Bact_data$Std_Order %in% Bact_traits_all$Std_Order])
# unique(Bact_data$Std_Order[!Bact_data$Std_Order %in% Bact_traits_all$Std_Order])


# *************************** #
#### 2. Trait manipulation ####
# *************************** #

## Copiotrophs/oligotrophs
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

## Match with abundance data ##

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
Trait_genus[, est_volume := d1_up.mean*d1_up.mean*3.14*d2_up.mean/4]
Trait_genus_order[, est_volume := d1_up.mean*3.14*d1_up.mean*d2_up.mean/4]
Trait_genus[, Elongation := d1_up.mean/d2_up.mean]
Trait_genus_order[, Elongation := d1_up.mean/d2_up.mean]
Trait_genus_order[, Genome_size := genome_size.mean/1000] # Convert to M bd

Trait_genus[, c('logLength', 'logDiameter', 'logVolume', 'logDoubling_time') := lapply(.SD, log), .SDcols = c('d1_up.mean', 'd2_up.mean', 'est_volume',  'doubling_h.mean')]
Trait_genus_order[, c('logLength', 'logDiameter', 'logVolume', 'logDoubling_time') := lapply(.SD, log), .SDcols = c('d1_up.mean', 'd2_up.mean', 'est_volume',  'doubling_h.mean')]

traits_all = c("d2_up.mean", "d1_up.mean" , "doubling_h.mean", "oli_copio.value", "Genome_size", 'logLength', 'logDiameter', 'logVolume', 'logDoubling_time', 'Elongation')


# Save species-matched trait data

traitUnits = c("microm (log-transformed)","microm (log-transformed)",'microm^3 (log-transformed)', 'Mbp', "hour (log-transformed)",'unitless', 'unitless')
traitDescription = c("Largest length",
                     "Largest diameter",
                     "Volume of the organism calculated as d1*3.14*d1*d2/4",
                     'Genome size',
                     "Minimum doubling time in hours",
                     "Cell elongation calculated as d1/d2",
                     "Copiotrophic or oligotrophic bacteria")
traitDataID = c('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')

traitRef = c("https://doi.org/10.1038/s41597-020-0497-4",
             "https://doi.org/10.1038/s41597-020-0497-4",
             "Calculated manually from https://doi.org/10.1038/s41597-020-0497-4",
             'https://doi.org/10.1038/s41597-020-0497-4',
             'https://doi.org/10.1038/s41597-020-0497-4',
             'Calculated manually from https://doi.org/10.1038/s41597-020-0497-4',
             '(Actinobacteria, Betaproteobacteria, Gammaproteobacteria, Proteobacteria, Bacteroidetes were considered as copiotrophic; Acidobacteria, Verrucomicrobia, Planctomycetes as oligotrophic)')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'oli_copio.value')


Bacteria_trait_melt = melt.data.table(unique(Trait_genus_order[, .SD, .SDcols = c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'oli_copio.value', 'Std_Genus')]), id.vars = c( 'Std_Genus'), , variable.name = 'traitName', value.name = 'traitValue')
Bacteria_trait_info = add_info(Bacteria_trait_melt[traitName %in%  c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'oli_copio.value')], traitRefs = traitRef, traitDataIDs = traitDataID, traitDescriptions = traitDescription, traitUnits = traitUnits, traitsOnly = TRUE)
Bacteria_trait_info[traitName == 'oli_copio.value', traitName := 'Oligo_copio' ]
fwrite(Bacteria_trait_info, "Data/Temporary_data/Bacteria_traits.csv")


# ************************** #
#### 2. Species-level PCA ####
# ************************** #

traits_test = Trait_genus_order[Std_Genus %in% unique(Abundance_genus$Std_Genus), list(logVolume, logDoubling_time, Elongation, Genome_size)]
pca_mic = dudi.pca(complete(mice(traits_test)), , scannf = FALSE, nf = 2)
fviz_pca(pca_mic)

# ********************************** #
#### 3. Community-weighted traits ####
# ********************************** #

# Use order values when needed
CC_genus_order  <- check_coverage(Trait_genus_order, Abundance_genus, traits_all, "Std_Genus", "Std_Genus")

# Add to match CWM categorical variables
CC_genus_order[, OC_ratio := oli_copio.value]
Trait_genus_order[, lapply(.SD, function(x){length(x[!is.na(x)])})]

# Calculate CWM
CWM_bacterias_GO <- my_cwm(Trait_genus_order, Abundance_genus, traits_all, "Std_Genus", "Std_Genus")
CWM_bacterias_GO[, OC_ratio := oli_copio.value_oligo/oli_copio.value_copio]

# Melt and merge
CWM_CC_bacterias = merge.data.table(melt.data.table(CWM_bacterias_GO[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                       melt.data.table(CC_genus_order[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))
CWM_CC_bacterias[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1,3), 0, substr(Plot, 4,4), sep = ''))]

CWM_CC_bacterias = CWM_CC_bacterias[traitName %in% c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'OC_ratio')]
# Add info

traitUnits = c("microm (log-transformed)","microm (log-transformed)",'microm^3 (log-transformed)', 'Mbp', "hour (log-transformed)",'unitless', 'unitless')
traitDescription = c("Largest length",
                     "Largest diameter",
                     "Volume of the organism calculated as d1*3.14*d1*d2/4",
                     'Genome size',
                     "Minimum doubling time in hours",
                     "Cell elongation calculated as d1/d2",
                     "Ratio of copiotrophic:oligotrophic bacteria (Actinobacteria, Betaproteobacteria, Gammaproteobacteria, Proteobacteria, Bacteroidetes were considered as copiotrophic; Acidobacteria, Verrucomicrobia, Planctomycetes as oligotrophic)")
traitDataID = c('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')

traitRef = c("https://doi.org/10.1038/s41597-020-0497-4",
             "https://doi.org/10.1038/s41597-020-0497-4",
             "Calculated manually from https://doi.org/10.1038/s41597-020-0497-4",
             'https://doi.org/10.1038/s41597-020-0497-4',
             'https://doi.org/10.1038/s41597-020-0497-4',
             'Calculated manually from https://doi.org/10.1038/s41597-020-0497-4',
             'Calculated manually')

names(traitRef) = names(traitDataID) = names(traitDescription) = names(traitUnits) = c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'OC_ratio')


CWM_CC_bacterias = add_info(CWM_CC_bacterias, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_bacterias, "Data/CWM_data/CWM_microbes.csv")


# ************************************** #
#### 4. Non-weighted community traits ####
# ************************************** #

# To have non-weighted traits we need to rarefy first
Abundance_genus_agg = Abundance_genus[, list(value= sum(value)), by = list(Plot, taxo)]

bact.cast<-dcast.data.table(Abundance_genus[, list(value = value, 'taxo_Year' = paste(taxo, Year, sep = '__'), Plot)],taxo_Year~Plot,value.var="value",fill=0)
bact.cast_matrix<-as.matrix(bact.cast[,-1])
rownames(bact.cast_matrix)<-bact.cast$taxo_Year
bact_spec<-otu_table(bact.cast_matrix, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
bact_phylo<-phyloseq(bact_spec) # conversion step two for getting a phyloseq object
Bacteria_rarefied <- rarefy_even_depth(bact_phylo,
                                       sample.size = min(sample_sums(bact_phylo)),
                                       rngseed = 1, replace = FALSE, trimOTUs = TRUE,
                                       verbose = TRUE) #sample size should be the smallest number of sequences per #sample
## Transform to usable dataframe/datatable format in the long format
Bacteria_rarefied = data.frame(Bacteria_rarefied)
Abundance_genus_presence_absence = melt.data.table(data.table(Bacteria_rarefied)[, c('taxo', 'Year') := tstrsplit(rownames(Bacteria_rarefied), '__')],
                                         id.var = c('taxo', 'Year'), value.name = 'value', variable.name = 'Plot')

Abundance_genus_presence_absence[value >1, value := 1]
Abundance_genus_presence_absence[, sum(value), by = taxo]


# Use order values when needed
CC_genus_order_noweight  <- check_coverage(Trait_genus_order, Abundance_genus_presence_absence, traits_all, "Std_Genus", "taxo")
# Add to match CWM categorical variables
CC_genus_order_noweight[, OC_ratio := oli_copio.value]

# Calculate CWM
CWM_bacterias_GO_noweight <- my_cwm(Trait_genus_order, Abundance_genus_presence_absence, traits_all, "Std_Genus", "taxo")
CWM_bacterias_GO_noweight[, OC_ratio := oli_copio.value_oligo/oli_copio.value_copio]

# Melt and merge
CWM_CC_bacterias_noweight = merge.data.table(melt.data.table(CWM_bacterias_GO_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitValue'),
                                    melt.data.table(CC_genus_order_noweight[, Year := as.character(Year)], id.vars = c('Plot', 'Year'), variable.name = "traitName", value.name = 'traitCoverage'))
CWM_CC_bacterias_noweight[, Plot := ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1,3), 0, substr(Plot, 4,4), sep = ''))]

CWM_CC_bacterias_noweight = CWM_CC_bacterias_noweight[traitName %in% c('logLength', 'logDiameter', 'logVolume', 'Genome_size', 'logDoubling_time', 'Elongation', 'OC_ratio')]
# Add info
CWM_CC_bacterias_noweight = add_info(CWM_CC_bacterias_noweight, traitRef, traitDataID, traitDescription, traitUnits, c('21446, 21447, 21448, 21449, 24690, 25306 synthesised in 27707'))

fwrite(CWM_CC_bacterias_noweight, "Data/CWM_data/CWM_microbes_noweight.csv")

