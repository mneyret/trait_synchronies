# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Trying to match BE fungal data with Funfun database #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")

install.packages("traitdataform")
library(traitdataform)
library(fungaltraits)
library(data.table)
library(readxl)
library(ggplot2)
setwd("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Fungi")

### Load data
# Trait data
#Funfun_data = setDT(fungal_traits())
#setDT(Funfun_data)
for (i in unique(Funfun_data$speciesMatched)){
  print(i)
  taxo = try(get_gbif_taxonomy(i))
  if (class(taxo) == "try-error") next
  sp = try(taxo$scientificName)
  o = try(taxo$order)
  f = try(taxo$family)
  if (class(sp) == "try-error" | is.null(sp)) sp = i 
  if (class(f) == "try-error" | is.null(f)) f = NA
  if (class(o) == "try-error" | is.null(o)) o = NA
  Funfun_data[speciesMatched == i, c('scientificName', 'Family', 'Order') := list(sp, f, o)]
}

#fwrite(Funfun_data, 'Funfun_data_taxo.csv')
Funfun_data =fread( 'Funfun_data_taxo.csv')
Baessler <- data.table(read_excel("~/Desktop/Research/Senckenberg/Data/Traits/Fungi/Baessler_JApplEcol_2014.xlsx"))
Baessler[, scientificName := lapply(Species, function(x){
          y = get_gbif_taxonomy(x)
          return(y$scientificName)})]


### Abundance data
#Total
Abundances_all = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Dataset_clean.txt")
lookup_table = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/26473_2_Dataset/26473_2_data.csv')
Abundances_fungi = Abundances_all[Group_broad == 'soilfungi',]
Abundances_fungi[, ASV := gsub('soilf_OTU', 'ASV', Species)]
Abundances_fungi = merge(Abundances_fungi[,-c(1,2)], lookup_table[,-2], by = 'ASV', all.x = T)
#data_AMF <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/19786.txt")
#data_fungi <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/21048.txt")
#data_fungi = data_fungi[nchar(Plotid) < 6,]
#data_fungi[, Plot := ifelse(nchar(Plotid) == 5, Plotid, paste(substr(Plotid, 1, 3), 0, substr(Plotid, 4, 4), sep = ''))]
#data_AMF[, Plot := ifelse(nchar(Plotid) == 5, Plotid, paste(substr(Plotid, 1, 3), 0, substr(Plotid, 4, 4), sep = ''))]
#all_species = data.table(rbind(data_AMF[, c("Plot","Abundance", "OTU",  "kingdom", "pyhlum", "class", "order", "family","genus", "species")],
#                               data_fungi[, c("Plot","Abundance", "OTU",  "kingdom", "pyhlum", "class", "order", "family","genus", "species")])
#)
#all_species[, scientificName:= '']
Abundances_fungi[, scientificName:= get_gbif_taxonomy(Species)$scientificName, by = Species]

fungal_traits  = c(
  'extension_rate','fruiting_body_size',
  'guild_fg', 'tissue_c', 'spore_size','spore_length', 
  'fruiting_body_size', 'total_genes', 'tissue_cn', 'tissue_cp')

#Size_data_Luang_2020 <- data.table(read_excel("/Users/Margot/Desktop/Research/Senckenberg/Data/Traits/Fungi/Size_data_Luang_2020.xlsx"))
#Size_data_fungi = Size_data_Luang_2020[Domain == 'Fungi' & Order != '',]
#Size_data_fungi[, c('Phylum', 'Class', 'Order','Family', 'Genus', 'Species') :=
#                  list(gsub('p__','', Phylum),
#                       gsub('c__','', Class),
#                       gsub('o__','', Order),
#                       gsub('f__','', Family),
#                       gsub('g__','', Genus),
#                       gsub('s__','', Species))]
#Size_data_fungi[, Size := .SD, .SDcols = c('Propagule size (μm)')]
#ggplot(Size_data_fungi, aes(Size, Order)) + geom_boxplot()

#order_in_data = sort(unique(c(all_species$family, all_species$order)))
#order_in_data = order_in_data[!grepl('ceae', order_in_data)]
#order_in_size = sort(unique(Size_data_fungi$Order))

#all_species[, Order := ifelse(order %in% order_in_size, order,
#                              ifelse(family %in% order_in_size, family,
#                                     NA))]
#all_species = merge(all_species, Size_data_fungi[, list(Size = mean(Size, na.rm = T)), by = Order], all.x = T)

#size_cwm = all_species[, list(spore_size=sum(Abundance*Size, na.rm = T)/sum(Abundance, na.rm = T)), by = Plot]

## Try analysis at genus level
 # Check that genus is the accepted one
Abundances_fungi[, genus := tstrsplit(scientificName, ' ')[1]]
Funfun_data[, genus := tstrsplit(scientificName, ' ')[1]]

Funfun_genus = merge.data.table(
                     x = Funfun_data[, lapply(.SD, function(x){mean(x, na.rm = T)}), .SDcols = fungal_traits[c(1,3:9)], by = c('genus')],
                     y = Funfun_data[, lapply(.SD, function(x){paste(unique(x[!is.na(x) & x != 'NULL']), collapse = '_')}), .SDcols = fungal_traits[2], by = c('genus')])


# Check coverage at the genus level
CC = check_coverage(Funfun_genus, Abundances_fungi[!(genus %in% c('', 'unclassified')),], unique(fungal_traits), 'genus',  'genus')

round(CC[, lapply(.SD, mean),.SDcols = fungal_traits]*100)
CC[, lapply(.SD, range),.SDcols = fungal_traits]
#==> Not very promising

## Try analysis at Family level
# Check that genus is the accepted one
Funfun_family = merge.data.table(
  x = Funfun_data[, lapply(.SD, function(x){mean(x, na.rm = T)}), .SDcols = fungal_traits[c(1,3:9)], by = c('Family')],
  y = Funfun_data[, lapply(.SD, function(x){paste(unique(x[!is.na(x) & x != 'NULL']), collapse = '_')}), .SDcols = fungal_traits[2], by = c('Family')])

# Check coverage at the genus level
CC = check_coverage(Funfun_family, Abundances_fungi[!(genus %in% c('', 'unclassified')),], unique(fungal_traits), 'Family',  'Family')

round(CC[, lapply(.SD, mean),.SDcols = fungal_traits]*100)
CC[, lapply(.SD, range),.SDcols = fungal_traits]
#==> Not very promising