# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# Trying to match BE fungal data with Funfun database #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

#install.packages("devtools")
#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("traitecoevo/fungaltraits")

#install.packages("traitdataform")
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
Funfun_data = fread( 'Funfun_data_taxo.csv')
Fungaltraits = data.table(read_excel('FungalTraits1.2_ver_16Dec_2020.xlsx'))
#Baessler <- data.table(read_excel("~/Desktop/Research/Senckenberg/Data/Traits/Fungi/Baessler_JApplEcol_2014.xlsx"))
#Baessler[, scientificName := lapply(Species, function(x){
#          y = get_gbif_taxonomy(x)
#          return(y$scientificName)})]


### Abundance data
#Total
Abundances_all = fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Dataset_clean.txt")
lookup_table = fread('/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/26473_2_Dataset/26473_2_data.csv')
Abundances_fungi = Abundances_all[Group_broad == 'soilfungi',]
Abundances_fungi[, ASV := gsub('soilf_OTU', 'ASV', Species)]
Abundances_fungi = merge(Abundances_fungi[,-c(1,2)], lookup_table[,-2], by = 'ASV', all.x = T)
data_AMF <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/19786.txt")
#data_fungi <- fread("/Users/Margot/Desktop/Research/Senckenberg/Data/Abundances/Fungi/21048.txt")
#data_fungi = data_fungi[nchar(Plotid) < 6,]
#data_fungi[, Plot := ifelse(nchar(Plotid) == 5, Plotid, paste(substr(Plotid, 1, 3), 0, substr(Plotid, 4, 4), sep = ''))]
data_AMF[, Plot := ifelse(nchar(Plotid) == 5, Plotid, paste(substr(Plotid, 1, 3), 0, substr(Plotid, 4, 4), sep = ''))]
all_species = data.table(rbind(data_AMF[, c("Plot","Abundance", "OTU",  "kingdom", "pyhlum", "class", "order", "family","genus", "species")],
                               data_fungi[, c("Plot","Abundance", "OTU",  "kingdom", "pyhlum", "class", "order", "family","genus", "species")])
)
#all_species[, scientificName:= '']
Abundances_fungi[, scientificName := get_gbif_taxonomy(Species)$scientificName, by = Species]

fungal_traits  = c(
  'extension_rate',#'growth_form_fg',
  #'guild_fg',
   'tissue_c', 'spore_size','spore_length', 'spore_width', 'spore_shape',
  'fruiting_body_size', 'total_genes', 'tissue_cn')

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

## What families/orders are mostly represented for each trait
Funfun_data[, genus := tstrsplit(scientificName, ' ')[1]]
Funfun_data[, spore_shape := spore_length/spore_width]

coltraits = c(fungal_traits, "Family", "Order", "genus")
Funfun_data01 = copy(Funfun_data[, .SD, .SDcols = coltraits])
Funfun_data01[, fungal_traits] = Funfun_data01[,lapply(.SD,function(x) ifelse(is.na(x),0,1)),.SDcols=fungal_traits]
Funfun_data01$guild = Funfun_data$guild_fg


## Try analysis at genus level
 # Check that genus is the accepted one
Abundances_fungi[, genus := tstrsplit(scientificName, ' ')[1]]

ggplot(Funfun_data[!is.na(fruiting_body_size),], aes(fruiting_body_size, x = Genus, fill = Order)) + geom_boxplot()


Funfun_genus = #merge.data.table(
                      Funfun_data[, lapply(.SD, function(x){mean(x, na.rm = T)}), .SDcols = fungal_traits, by = c('genus')]
                   #  y = Funfun_data[, lapply(.SD, function(x){paste(unique(x[!is.na(x) & x != 'NULL']), collapse = '_')}), .SDcols = fungal_traits[2], by = c('genus')]
                  #   )

Funfun_genus[genus %in%  Abundances_fungi[!(genus %in% c('', 'unclassified')),genus],unique(genus)]

# Check coverage at the genus level
CC = check_coverage(Funfun_genus, Abundances_fungi[Order %in% c('Agaricales', 'Russulales', 'Boletales')
  ,], unique(fungal_traits), 'genus',  'Genus')

round(CC[, lapply(.SD, mean),.SDcols = fungal_traits]*100)
CC[, lapply(.SD, range),.SDcols = fungal_traits]
#==> Not very promising


Funfun_genus[, spore_shape := spore_length/spore_width]
CWM_fungi = my_cwm(Funfun_genus, Abundances_fungi[ Order %in% c('Agaricales', 'Russulales', 'Boletales') & Genus != 'unclassified',], unique(fungal_traits), 'genus',  'Genus')

CWM_fungi1 = data.table(complete(mice(CWM_fungi)))
CWM_fungi = CWM_fungi1[, lapply(.SD, mean), by = 'Plot', .SDcols = c('tissue_cn', 'spore_size', 'spore_shape', 'fruiting_body_size')]

pca = dudi.pca(CWM_fungi[, c('tissue_cn', 'spore_size', 'spore_shape', 'fruiting_body_size')])

quanti.coord <- supcol(pca, env_data_lui[, c('Fertil.', 'LUI', 'Mowing','Grazing')])$cosup * pca$eig[1:2]^2
pcaplot = fviz_pca(pca)
pcaplot = fviz_add(pcaplot, quanti.coord, axes = c(1, 2), "arrow", color = "black", linetype = "solid"#,# repel = T,
                      #addlabel = T
)

cor(pca$li[,1], env_data_lui$LUI)

# Other test with the FungalTraits dataset
test_sp = Fungaltraits[!is.na(Ectomycorrhiza_exploration_type_template), GENUS]
test_sp = Fungaltraits[!is.na(Decay_type_template), GENUS]
Fungaltraits[, Explo := ifelse(Ectomycorrhiza_exploration_type_template == 'contact',0,
                               ifelse(grepl('short-distance', Ectomycorrhiza_exploration_type_template), 1,
                                      ifelse(grepl('medium-distance', Ectomycorrhiza_exploration_type_template), 2,
                                             ifelse(Ectomycorrhiza_exploration_type_template == 'long-distance',3,
                                                    NA))))]
#Fungaltraits[, Genus := get_gbif_taxonomy(GENUS)$genus, by = GENUS]

Abundances_fungi[Fun_group_fine == "Saprotroph" & Genus %in% test_sp, sum(value)] /Abundances_fungi[Fun_group_fine == "Saprotroph", sum(value)]

CWM_explo = my_cwm(Fungaltraits, Abundances_fungi[Fun_group_fine == "EMF" ,], c('Explo', 'primary_photobiont'), 'GENUS',  'Genus')
Abundances_fungi[Genus %in% Fungaltraits[!is.na(Explo), GENUS], sum(value)] /Abundances_fungi[, sum(value)]
