
# CWM calculation
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
}

load_env_data = function(){
  data_lui_raw <- fread("Data/Environment_function_data/LUI_standardized_global.txt") 
  
  data_lui = data_lui_raw[Year > 2007 & Year <= 2018, list(LUI = mean(LUI),
                                                           Mowing = mean(M_STD),
                                                           Fertil = mean(F_STD),
                                                           Grazing = mean(G_STD)), by = list(Plot = ifelse(nchar(PLOTID) == 5,PLOTID, paste(substr(PLOTID, 1, 3), '0', substr(PLOTID, 4, 4), sep = '')))]
  data_lui[, c('Fertilisation', 'Disturbance') := list( Fertil, sqrt(Mowing/mean(Mowing) + Grazing/mean(Grazing)))]
  
  env_data <- data.table(read_excel("Data/Environment_function_data/31018_5_Dataset/31018_5_data.xslx", sheet = 1))
  env_data[, Plot := EP_PlotID]
  env_data_lui = merge.data.table(data_lui, env_data[, c('Plot', 'LII',   "Soil.pH", "Soil.depth", "Soil.clay.content" , "Soil.sand.content", "TWI" , 'Grassland.1000')], by = 'Plot')
  
  Temp = fread("Data/Environment_function_data/Ta_200_2008_2018_3a7aa69b37636254/plots.csv")
  Precip = fread("Data/Environment_function_data/precipitation_radolan_95a5c5f55a798133/plots.csv")
  Climate_data <- merge.data.table(Temp[, plotID := gsub('f', '', plotID)], Precip, by = c('plotID', 'datetime'))
  Climate_data[, Plot := plotID]
  Climate <- Climate_data[datetime>=2008 & datetime <= 2018, list("Mean_Temp" = mean(Ta_200, na.rm = T), "Mean_precip" = mean(precipitation_radolan, na.rm = T)), by = list("Plot" = plotID)]
  
  # Merge all environmental data
  env_data_lui = merge.data.table(env_data_lui, Climate, by = 'Plot')
  env_data_lui$Region = substr(env_data_lui$Plot, 1, 1)
  
  # Selection of environmental variables: removal of highly correlated variables (Soil.depth, Soil.sand.content, Mean_precip)
  env_var = c("Soil.pH", "Soil.clay.content", "TWI", 'Mean_Temp')
  env_data_lui[, (c(env_var, 'LUI', 'Mowing', 'Fertil', 'Grazing')) := lapply(.SD, function(x){as.numeric(scale(x))}), .SDcols = c(env_var, 'LUI', 'Mowing', 'Fertil', 'Grazing')]
  return(list(env_data, env_data_lui, env_var, data_lui_raw))}


load_eco_functions = function(){
  # Load datasets
  additional_info <- fread('Data/Environment_function_data/27087_21_Dataset/synthesis_grassland_function_metadata_ID27087.csv')
  additional_info <- additional_info[, .(ColumnName, AggregatedColumnName, codedYear)]
  setnames(additional_info, old = c("ColumnName", "codedYear"), new = c("variable", "Year"))
  original_synth_func <- fread('Data/Environment_function_data/27087_21_Dataset/27087_21_data.csv')
  
  # get function-year combinations as columns
  # make format even longer to obtain a column "Year" and a column "function"
  synth_func <- melt.data.table(original_synth_func, id.vars = c("Plot", "Plotn", "Explo", "Year"))
  sum(is.na(synth_func$value))
  # delete all missing function-year combinations (without excluding NA values)
  synth_func <- synth_func[!value %in% "NM"]
  sum(is.na(synth_func$value))
  
  synth_func <- merge(synth_func, additional_info, by = c("variable", "Year"))
  synth_func <- synth_func[, .(AggregatedColumnName, Plot, Plotn, Explo, value)]
  synth_func <-dcast.data.table(synth_func, Plot + Plotn + Explo ~ AggregatedColumnName, value.var = "value")
  synth_func[, colnames(synth_func)[4:87] := lapply(.SD, as.numeric), .SDcols = colnames(synth_func)[4:87]]
  setkey(synth_func, 'Plot')
  
  ## Also the respiration dataset
  Resp = fread('Data/Environment_function_data/26908_4_Dataset/26908_4_data.csv')
  Resp[, Plot := ifelse(nchar(EP_Plotid) == 5, EP_Plotid, paste(substr(EP_Plotid, 1, 3), '0', substr(EP_Plotid, 4, 4), sep = ''))]
  
  EF_clean = synth_func[, list(Plot = ifelse(nchar(Plot) == 5, Plot, paste(substr(Plot, 1, 3), 0, substr(Plot, 4, 4),sep = '')),
                               Dung.decomposition = as.numeric(scale(dung.removal)),
                               Litter.decomposition =  as.numeric(scale(Litter.decomposition)),
                               Root.decomposition=  as.numeric(scale(Root.decomposition)),
                               Biomass.production = as.numeric(scale(Biomass)),
                               'Urease'= as.numeric(scale(sqrt(Urease))),
                               'DEA'= as.numeric(scale(sqrt(DEA))),
                               'Potential.nitrification'= as.numeric(scale(sqrt(Potential.nitrification2011)))+
                                 as.numeric(scale(Potential.nitrification2014)),
                               'nifH'= as.numeric(scale(sqrt(nifH))),
                               'amoA_AOB'= as.numeric(scale(sqrt(amoA_AOB.2011 + amoA_AOB.2016))),
                               'amoA_AOA'= as.numeric(scale(sqrt(amoA_AOA.2011 + amoA_AOA.2016))),
                               'nxrA_NS'=  as.numeric(scale(sqrt(nxrA_NS))),
                               'n16S_NB'=  as.numeric(scale(sqrt(`16S_NB`))),
                               'beta_Glucosidase' = as.numeric(scale(sqrt(beta_Glucosidase))),
                               'N_Acetyl_beta_Glucosaminidase' = as.numeric(scale(sqrt(N_Acetyl_beta_Glucosaminidase))),
                               'Xylosidase' = as.numeric(scale(sqrt(Xylosidase))),
                               'Phosphatase' = as.numeric(scale(log(Phosphatase)))
  )
  ]
  EF_clean = merge.data.table(EF_clean, Resp[!grepl('W', Plot),  list(Respiration = mean(Rs)), by = Plot])
  return(EF_clean)}
# Function to check how much % of the individuals (or cover) has trait data available

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


# Function to match trait data with metadata
add_info = function(data, traitRefs, traitDataIDs, traitDescriptions, traitUnits, AbundanceID, traitsOnly = FALSE){
  data[, traitUnit := dplyr::recode(traitName, !!!traitUnits)]
  data[, traitDescription := dplyr::recode(traitName, !!!traitDescriptions)]
  data[, traitDataRef := dplyr::recode(traitName, !!!traitRefs)]
  data[, TraitDataID := dplyr::recode(traitName, !!!traitDataIDs)]
  if (traitsOnly == FALSE){
  data[, AbundanceDataID := AbundanceID]}
  return(data)
}

# Function to get gbif taxonomy (wrapper around get_gif_taxonomy function)
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

# Function to calculate the coefficients of the SEM linking each trait to the fertilisation, mowing and grazing components of Land-use intensity
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


## Function to correct CWM and function data for region and/or environmental covariates
corrections <- function(CWM_data, Env_data = NULL, env_corr = FALSE,  variables = NULL, env_variables = NULL) {
  CWM_data_copy <- CWM_data[complete.cases(CWM_data),]
  if (env_corr == TRUE & is.null(Env_data)) {
    stop("Please provide environmental data")
  }
  if (is.null(variables)) {
    variables <- colnames(CWM_data_copy)[!(colnames(CWM_data_copy) %in% c("Plot", "Year"))]
  }
  
  if (is.null(env_variables)) {
    env_variables <- colnames(Env_data)[!(colnames(Env_data) %in% c("Plot", "Year"))]
  }
  
  if (env_corr == TRUE) {
    CWM_data_copy <- merge(CWM_data_copy, Env_data[, .SD, .SDcols = c("Plot", env_variables)], by = "Plot")
    
    CWM_data_copy[, (variables) := lapply(.SD, function(x) {
      mod <- lm(x ~ ., CWM_data_copy[, .SD, .SDcols = env_variables])
      return(residuals(mod))
    }),
    .SDcols = variables
    ]
    CWM_data_copy = CWM_data_copy[, .SD, .SDcols = c(variables, 'Plot')]
  }
  return(CWM_data_copy)
}

## Function to correct CWM and function data for region and/or environmental covariates, long format
corrections_long <- function(value, env_data, env_corr = FALSE) {
  if (env_corr == TRUE) {
    mod <- lm(value ~ ., env_data)
    return(residuals(mod))
  } else {
    return(value)
  }
}

## Function to calculate individual responses to fertilisation, mowing and grazing
calculate_coeff_new = function(DATA){
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
  parameters = model_parameters(lui_mod)
  res$estimate = parameters$Coefficient[2]
  res$p = parameters$p[2]
  res$ci_low =  round(parameters$CI_low[2], 2)
  res$ci_high = round(parameters$CI_high[2], 2)
  res$ci =  paste(res$ci_low,res$ci_high, sep = ' - ')
  res$n =  nobs(lui_mod)
  res$r2 =  rsq(lui_mod)
  return(res)
}

## Function to run PCA, renaming rows as needed
my_dudi_pca <- function(CWM_matrix, col_to_use, NF = 3) {
  matrix <- data.frame(CWM_matrix[complete.cases(CWM_matrix[, ..col_to_use])])
  rownames(matrix) <- matrix$Plot
  pca <- dudi.pca(matrix[, col_to_use], scannf = FALSE, nf = NF)
  return(pca)
}

## Function to scale between 0 and 1 (used for colors mostly)
scale01 = function(x){
  return((x-min(x))/(max(x)-min(x)))
}

## Function to change the direction of PCA axis (for reproductibility of axes direction)
axes_direction_pca <- function(pca, cwm, traitnames, conditions) {
  if (length(conditions) != length(traitnames)) {
    print("Number of axes, traits, and conditions should be equal")
  }
  else {
    cond_positive <- ifelse(gsub(" ", "", conditions) == ">0", TRUE, FALSE)
    I <- length(traitnames)
    for (i in 1:I) {
      if ((pca$li[cwm[get(traitnames[i]) == max(get(traitnames[i])), Plot], i] > 0) != cond_positive[i]) {
        pca$li[, i] <- -pca$li[, i]
        pca$co[, i] <- -pca$co[, i]
      }}}
  return(pca)
}

## Function to print the large table with individual trait responses
print_table = function(Hypo_table){
  Table_use = Hypo_table[,.SD, .SDcols = c("Trait_short", "Expected_direction", 'Est','CI', 'r2', 'Padj', 'Fit_exp')]
  H_kable <- kable_paper(kbl(Table_use, escape = FALSE))
  H_kable <- column_spec(H_kable, which(colnames(Table_use) == "Est"),
                         background = rgb(pal.fnc(scale01(Hypo_table$Est)), maxColorValue=255),
                         color = 'white',
                         bold  = Hypo_table$Signif == "YES",
                         italic  = !Hypo_table$Signif == "no"
  )
  H_kable <- column_spec(H_kable, which(colnames(Table_use) == "Fit_exp"),
                         background = ifelse(is.na(Hypo_table$Fit_exp) | Hypo_table$Fit_exp == 'Inconclusive', 'grey', 'white')
  )
  
  H_kable = gsub( '+', '&#43;', H_kable, fixed = T)
  H_kable = gsub( '@@@', '&#', H_kable)
  H_kable
}

# Function to auto run PCAs
run_group_pca = function(cwm, traits, trait_direction, direction = '>0', plot = TRUE, env_corr, Labels = NA, annot = NA, naxes = 2){
  
  data = data.frame(cwm[, .SD, .SDcols = c('LUI', traits)])
  rownames(data) = cwm$Plot
  data = data[complete.cases(data),]
  pca_stand = PCA(data, graph=FALSE, quanti.sup = 1)
  
  # Store the direction of the expected fast_slow axis (1 is left to right, -1 is right to left)
  direction_axis = numeric()
  pca_to_plot = pca_stand
  
  for (i in 1:length(direction)){
    if (pca_stand$var$coord[trait_direction[i], i] > 0 & direction[i]  == '>0' |
        pca_stand$var$coord[trait_direction[i], i] < 0 & direction[i]  == '<0'){
      direction_axis = c(direction_axis, 1)
    } else {
      direction_axis = c(direction_axis,-1)
    }
    pca_to_plot$var$coord[,i] = pca_to_plot$var$coord[,i]*direction_axis[i]
    pca_to_plot$ind$coord[,i] = pca_to_plot$ind$coord[,i]*direction_axis[i]
    pca_to_plot$quanti.sup$coord[,i] = pca_to_plot$quanti.sup$coord[,i]*direction_axis[i]
  }
  
  naxes = length(traits)
  
  corr_unsignificant = 0
  for (trait in traits){
    cor = cor.test(unlist(data[, trait]), pca_stand$ind$coord[,1])
    if (cor$p.value > 0.05){
      corr_unsignificant = corr_unsignificant+1
    }
  }
  
  keep_criteria = data.table('More_than_random' = pca_stand$eig[1] / (1/length(traits)),
                             'Prop.traits.40' =   length(pca_stand$var$cor[,1][abs(pca_stand$var$cor[,1])>0.4])/length(pca_stand$var$cor[,1]),
                             'Prop.traits.25' =   length(pca_stand$var$cor[,1][abs(pca_stand$var$cor[,1])>0.25])/length(pca_stand$var$cor[,1]),
                             'Cor_unsignificant' =  corr_unsignificant)
  keep_criteria[, enough_support := More_than_random >2 & Prop.traits.40 >= 0.60]
  keep_criteria[, partial := More_than_random >2 & Prop.traits.25 >= 0.60 & Cor_unsignificant <= 1]
  keep_criteria[, LUI.Axis1 := round(cor.test(data$LUI, pca_stand$ind$coord[,1])$p.value, 3)]
  keep_criteria[, LUI.Axis2 := round(cor.test(data$LUI, pca_stand$ind$coord[,1])$p.value, 3)]
  names(Labels)  = traits
  
  if (is.na(annot)){
    Criteria_label = ifelse(keep_criteria$enough_support, 'Criteria met', 'Criteria not met!')
    if (keep_criteria$partial == TRUE & keep_criteria$enough_support != TRUE){
      Criteria_label = paste('Criteria partially met', sep = ' ')
    }
    
    Criteria_label = paste(Criteria_label, '\n Axis1 x LUI: p =', keep_criteria$LUI.Axis1,
                           '\n Axis2 x LUI: p =', keep_criteria$LUI.Axis2, sep = '')
  }else{
    Criteria_label =  annot
  }
  
  rownames(pca_to_plot$var$coord) = rownames(pca_to_plot$var$cor) = rownames(pca_to_plot$var$cos2) = Labels
  
  if (env_corr == TRUE){
    plot_pca_12 = fviz_pca(pca_to_plot, title = '', geom = 'point', geom.var = c("arrow", 'text'), col.quanti.sup = 'brown3', col.var = 'grey30', 
                           labelsize = 4, repel = T)
  }
  if (env_corr == FALSE){
    plot_pca_12 = fviz_pca(pca_to_plot, title = '', geom = 'point', geom.var = c("arrow", 'text'), col.quanti.sup = 'brown3', col.var = 'grey30', habillage = factor(substr(rownames(pca_to_plot$ind$coord), 1, 1)), 
                           labelsize = 5, repel = T) + 
      scale_shape_discrete( breaks = c('A', 'H', 'S'), labels = c('South', 'Central', 'North'), name = 'Region')+
      scale_color_viridis(discrete = TRUE, begin = 0.3, breaks = c('A', 'H', 'S'), labels = c('South', 'Central', 'North'), name = 'Region')+
      theme(legend.position = 'none') 
    
  }
  res = list(PCA = pca_stand, plot12 = plot_pca_12, criteria = keep_criteria, direction = direction_axis)
  return(res)
}

# Function to format thetable testing all the hypotheses
return_hypo_result_color = function(Signif, Direction){
  res = NA
  if (Signif == 'no'){
    res = 'Inconclusive' } 
  if (Signif == 'yes'){
    if (Direction %in% c('a', 'b')) {res = Direction}
    if (Direction %in% c('yes', 'no')) {res = ifelse(Direction == "yes", "@@@10004;", "X")}  }
  if (Signif == 'YES'){
    if (Direction %in% c('a', 'b')) {res = toupper(Direction)}
    if (Direction %in% c('yes', 'no')) {res = ifelse(Direction == "yes", "@@@9989;", "@@@10060;")} }
  
  return(res)}

## Function to convert a SEM into a text variable to facilitate SEM handling
autowrite_SEM_full = function(model_raw, trophic_levels = Trophic_levels, LUI = T, Data = dd){
  # This function takes a model formula without parameters
  # and creates all the parameter names for direct, indirect and total effects
  model_init <- sem(model_raw,
                    data = Data)
  coefficients = data.table(parameterEstimates(model_init, standardized = T))
  coefficients = coefficients[lhs != rhs,]
  coefficients = coefficients[op == '~',]
  
  coefficients_lm = coefficients[op == '~',]
  exo_var = unique(coefficients_lm$rhs[!(coefficients_lm$rhs %in% coefficients_lm$lhs)])
  
  effect_list = list()
  model_lines = list()
  
  for (group in unique(coefficients$lhs)){
    # Write down main model line
    RHS = coefficients[lhs == group, ]
    right_hand = paste(paste(RHS$rhs, RHS$lhs, sep = '__'), '*', RHS$rhs, collapse = ' + ', sep = '')
    lm_line = paste(group, '~', right_hand)
    model_lines = append(model_lines, lm_line)
    
    for (exo in exo_var){
      # Direct effect
      direct = paste('Direct__', exo, '__', group, '__', trophic_levels[group], sep = '')
      
      direct_rhs = paste(' := ', exo, '__', group, sep = '')
      direct_line = paste(direct, direct_rhs , sep = '')
      model_lines = append(model_lines, direct_line)
      effect_list = append(effect_list, direct)
      
      # Indirect effect
      RHS_indirect = RHS[!(rhs %in% exo_var),]
      if (nrow(RHS_indirect) > 0){
        indirect = paste('Indirect__', exo, '__', group, '__', trophic_levels[group], sep = '')
        indirect_rhs = paste( 
          paste('Total__', exo, '__', RHS_indirect[, rhs], '__', trophic_levels[RHS_indirect[, rhs]], sep = ''), 
          '*',
          paste(RHS_indirect[, rhs], '__', group, sep = ''), sep = '', collapse = ' + ')
        indirect_line = paste( indirect , indirect_rhs, sep = ' := ')
        model_lines = append(model_lines, indirect_line)
        effect_list = append(effect_list, indirect)
        
        # Total effect
        total = paste('Total__', exo, '__', group, '__', trophic_levels[group], sep = '')
        
        total_rhs = paste(direct, indirect, sep = ' + ')
        total_line = paste(total, total_rhs, sep = ' := ')
        
        model_lines = append(model_lines, total_line)
        effect_list = append(effect_list, total)
      }
      else {
        total = paste('Total__', exo, '__', group, '__', trophic_levels[group], sep = '')
        total_line = paste(total, direct, sep = ' := ')
        model_lines = append(model_lines, total_line)
        effect_list = append(effect_list, total)
      }} }
  
  effect_list = unlist(effect_list)
  overall_effects = list()
  for (trophic_lvl in unique(trophic_levels[coefficients$lhs])){
    for (type in c('Direct', 'Indirect', 'Total')){
      #  print(type)
      if (!(type == 'Indirect' & trophic_lvl == 0)){
        for (exo in exo_var){
          effect = paste(type, exo, trophic_lvl, sep = '__')
          relevant_effects = effect_list[grepl(type, effect_list) & grepl(trophic_lvl, effect_list) & grepl(exo, effect_list)]
          rhs = paste('(', 
                      paste(relevant_effects, collapse = ' + '),
                      ')/', length(relevant_effects))
          overall_effects = append(overall_effects, paste(effect, rhs, sep = ':='))
        }
      }}
  }
  whole_model = paste(c(unlist(model_lines), unlist(overall_effects)), collapse = ' \n ')
  return(whole_model)
}



