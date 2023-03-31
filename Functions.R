my_cwm <- function(Traits0, Abundances0, trait_names, trait_taxo, abundance_taxo) {
  # Traits must be a LARGE data table with one column called 'taxo' (= species or whatever) and other colnames including trait_names
  # Abundance must be a LONG data table with one column called 'taxo' (= species or whatever), one "value" (= abundance), one "Plot" and one "Year"
  Traits <- Traits0
  Abundances <- Abundances0
  Traits[, taxo := .SD, .SDcols = trait_taxo]
  Abundances[, taxo := .SD, .SDcols = abundance_taxo]
  common_spe <- intersect(Traits$taxo, Abundances$taxo)
  common_spe = common_spe[!is.na(common_spe)]
  
  # n species in traits not in ab:
  print(paste(length(unique(Traits$taxo[!(Traits$taxo %in% common_spe)])), "species in traits but not found in abundance", collapse = " "))
  print(paste(length(unique(Abundances$taxo[!(Abundances$taxo %in% common_spe)])), "species in abundance but not found in traits", collapse = " "))
  
  Traits_mat <-Traits[taxo %in% common_spe, ..trait_names]
  char_col = colnames(Traits_mat)[sapply(Traits_mat, class) %in% c('character', 'factor')]
  
  Traits_mat[, (char_col) := lapply(.SD, function(x){
    x = as.character(x)
    x[is.na(x)] = 'NA'
    x[is.nan(x)] = 'NA'
    return(as.factor(x))}), .SDcols = char_col]
  Traits_mat = data.frame(Traits_mat)
  rownames(Traits_mat) = Traits[taxo %in% common_spe,]$taxo
  Traits_mat = Traits_mat[order(rownames(Traits_mat)),]
  
  Abundances[, plot_year := paste(Plot, Year, sep = '_')]
  Abundances_cast = dcast(Abundances[taxo %in% common_spe,], plot_year ~ taxo, value.var = 'value', fun.aggregate = sum)
  Abundances_mat <- as.matrix(Abundances_cast[, -1])
  Abundances_mat = Abundances_mat[rowSums(Abundances_mat, na.rm = T) != 0,]
  Abundances_mat <- Abundances_mat[, order(colnames(Abundances_mat))]
  rownames(Abundances_mat) = Abundances_cast[rowSums(Abundances_cast[,-1], na.rm = T) != 0,]$plot_year
  
  cwm = functcomp(Traits_mat, Abundances_mat, CWM.type = 'all')
  cwm$Plot = sapply(rownames(cwm), function(x){unlist(strsplit(x, '_'))[[1]]})
  cwm$Year = as.numeric(sapply(rownames(cwm), function(x){unlist(strsplit(x, '_'))[[2]]}))
  
  return(data.table(cwm))
  
#  Traits2 <- Traits[, .SD, .SDcols = c("taxo", trait_names)]
 # Abundances2 <- Abundances[, .SD, .SDcols = c("Year", "taxo", "value", "Plot")]
  
  #traits_ab <- merge.data.table(Traits2, Abundances2, by = "taxo", all = TRUE)
  
  #traits_ab = traits_ab[value > 0, ]
  #traits_ab$ID = 1:nrow(traits_ab)
  #traits_ab[, value2 := ifelse(is.na(taxo), 0, value)]
  #traits_ab_cover = traits_ab
  
  #MAX = traits_ab_cover[, list(MAX =  sum(value2, na.rm = T)/sum(value, na.rm = T)), by = c('Plot',"Year")]
  
  #trait_type = unlist(traits_ab_cover[, lapply(.SD, is.numeric), .SDcols = trait_names])
  #trait_quanti = names(trait_type[trait_type])
  #print(trait_quanti)
  #print("_______")
  #trait_quali = names(trait_type[!trait_type])
  #print(trait_type)
  #print(trait_quali)
  #old_col = copy(colnames(traits_ab))
  
  # For qualitative traits: create dummy variables
  
 # if (length(trait_quali)>0){
#  traits_ab_cover[, (trait_quali) := lapply(.SD, function(x){ifelse(x == 'NA', NA, x)}), .SDcols = trait_quali]
#  dummy_columns(traits_ab_cover, select_columns = trait_quali, ignore_na = TRUE)
#  dummycols = colnames(traits_ab_cover)[!(colnames(traits_ab_cover) %in% old_col) ]
#  traits = c(trait_quanti, dummycols)}
#  if (length(trait_quali)==0){traits = trait_quanti}
#  #print(traits)
#  cwm <- traits_ab_cover[,
#                   lapply(.SD, function(x) {
#                     sum(x[!is.na(x)] * value2[!is.na(x)], na.rm = T) / sum(value2[!is.na(x)], na.rm = T)
#                   }),
#                   by = c("Plot", "Year"),
#                   .SDcols = traits
#  ][order(Plot), ]
  
  
}

check_coverage <- function(Traits0, Abundances0, trait_names, trait_taxo, abundance_taxo) {
  # Traits must be a LARGE data table with one column called 'taxo' (= species or whatever) and other colnames including trait_names
  # Abundance must be a LONG data table with one column called 'taxo' (= species or whatever), one "value" (= abundance), one "Plot" and one "Year"
  Traits <- Traits0
  Abundances <- Abundances0
  Traits[, taxo := .SD, .SDcols = trait_taxo]
  Abundances[, taxo := .SD, .SDcols = abundance_taxo]
  common_spe <- intersect(Traits$taxo, Abundances$taxo)

  # n species in traits not in ab:
  print(paste(length(unique(Traits$taxo[!(Traits$taxo %in% common_spe)])), "species in traits but not found in abundance", collapse = " "))
  print(paste(length(unique(Abundances$taxo[!(Abundances$taxo %in% common_spe)])), "species in abundance but not found in traits", collapse = " "))

  Traits2 <- Traits[, .SD, .SDcols = c("taxo", trait_names)]
  Abundances2 <- Abundances[, .SD, .SDcols = c("Year", "taxo", "value", "Plot")]

  traits_ab <- merge.data.table(Traits2, Abundances2, by = "taxo", all = TRUE)

  traits_ab_cover <- traits_ab
  traits_ab_cover <- traits_ab_cover[value > 0, ]
  traits_ab_cover$ID <- 1:nrow(traits_ab_cover)
  traits_ab_cover[, value2 := ifelse(is.na(taxo) | taxo == ' genus', 0, value)]
  #MAX = traits_ab_cover[, list(MAX =  sum(value2, na.rm = T)/sum(value, na.rm = T)), by = c('Plot',"Year")]
  
  # For quantitative traits
  traits_ab_cover2 <- cbind(
    traits_ab_cover[, lapply(.SD, function(x) {
      x = unlist(x)
       x = ifelse(is.null(x) | (is.na(x) |  (x == "NA")), 0, 1)
   #   x[is.null(x)] <- 0
  #    x[(is.na(x) |  (x == "NA"))] <- 0
  #    x[x != 0] <- 1
   #   x <- as.numeric(x)
      return(x)
    }),
    by = ID,
    .SDcols = trait_names
    ],
    traits_ab_cover[, .SD, .SDcols = c("taxo", "Plot", "Year", "value", 'value2')]
  )
  
  

  final <- traits_ab_cover2[, lapply(.SD, function(x) {
    sum(x * value2, na.rm = T)/sum(value, na.rm = T)
  }),
  .SDcols = trait_names, by = c('Plot', 'Year')][order(Plot), ]
  return(final)
}

add_info = function(data, traitRefs, traitDataIDs, traitDescriptions, traitUnits, AbundanceID, traitsOnly = FALSE){
  data[, traitUnit := dplyr::recode(traitName, !!!traitUnits)]
  data[, traitDescription := dplyr::recode(traitName, !!!traitDescriptions)]
  data[, traitDataRef := dplyr::recode(traitName, !!!traitRefs)]
  data[, TraitDataID := dplyr::recode(traitName, !!!traitDataIDs)]
  if (traitsOnly == FALSE){
  data[, AbundanceDataID := AbundanceID]}
  return(data)
}

get_gbif = function(x){
  if(is.na(x)){
    return(NA)
  } else {
    if (grepl('sp[ec]*[\\.]*$',x)){
      x_genus = gsub('sp[ec]*[\\.]*$','', x)
      x_genus = gsub(' ', '', x_genus)
      genus = try(get_gbif_taxonomy(x_genus)$scientificName)
      if (class(genus) == 'try-error') {return(paste(x_genus, 'sp.'))}
      return(paste(genus, 'sp.'))
    }
    else{
      names = get_gbif_taxonomy(unique(x))
      if(names$warnings =="No matching species concept! "){
        return(names$verbatimScientificName)
      } else{
        return(names$scientificName)
      }}}
}

calculate_coeff = function(DATA){
  mod = 'value ~ Fertil 
         value ~ Mowing
         value ~ Grazing
      
      Mowing ~~ Fertil 
      Grazing ~~ Fertil 
      Grazing ~~ Mowing
'
  sem_components = sem(mod, data = DATA)
  res = as.list(partable(sem_components)$est[1:3])
  names(res) = partable(sem_components)$rhs[1:3]
  
  lui_mod = lm(value ~ LUI, data = DATA)
  coeff = coefficients(summary(lui_mod))[2,]
  res$estimate = coeff['Estimate']
  res$p = coeff['Pr(>|t|)']
  
  res$ci = paste(round(confint(lui_mod)[2,], 2), collapse = ' - ')
  
  lui_mod_uncorr = lm(value_uncorr ~ LUI, data = DATA)
  coeff = coefficients(summary(lui_mod_uncorr))[2,]
  res$estimate_uncorr = coeff['Estimate']
  res$p_uncorr = coeff['Pr(>|t|)']
  res$ci_uncorr = paste(round(confint(lui_mod_uncorr)[2,], 2), collapse = ' - ')
  
  return(res)
}



#### Functions ####

