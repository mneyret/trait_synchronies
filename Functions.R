my_cwm <- function(Traits0, Abundances0, trait_names, trait_taxo, abundance_taxo) {
  #library(fastDummies)
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


#my_functcomp =function (x, a, CWM.type = c("dom", "all"), bin.num = NULL) 
#{
#  if (!is.matrix(x) & !is.data.frame(x)) {
#    stop("'x' must be a matrix or a data frame.", "\n")
#  }  else {x <- data.frame(x)}
#  if (!is.matrix(a))   {stop("'a' must be a matrix.", "\n")}
#  if (is.null(row.names(x)))  { stop("'x' must have row names.", "\n")
#    } else {x.n <- row.names(x)}
#  if (is.null(colnames(a))) {stop("'a' must have column names.", "\n")
#    } else {a.n <- colnames(a)}
#  s.x <- dim(x)[1]
#  s.a <- dim(a)[2]
#  if (s.x != s.a) 
#    stop("Different number of species in 'x' and 'a'.", "\n")
#  if (any(x.n != a.n)) 
#    stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
#         "\n")
#  com <- dim(a)[1]
#  t <- dim(x)[2]
#  com.names <- row.names(a)
#  sp.names <- row.names(x)
#  tr.names <- names(x)
#  a[which(is.na(a))] <- 0
#  CWM.type <- match.arg(CWM.type)
#  is.bin <- function(k) all(k[!is.na(k)] %in% c(0, 1))
#  bin.var <- rep(NA, t)
#  names(bin.var) <- tr.names
#  for (i in 1:t) bin.var[i] <- is.bin(x[, i])
#  if (!all(bin.var[bin.num])) 
#    stop("'bin.num' points to non-binary variables.\n")
#  bin.var[bin.num] <- FALSE
#  type <- sapply(x, data.class)
#  type[type %in% c("numeric", "integer")] <- "C"
#  type[type == "ordered"] <- "O"
#  type[type == "factor"] <- "N"
#  type[bin.var] <- "B"
#  sum.a <- apply(a, 1, sum)
#  a <- a/sum.a
#  a <- t(a)
#  a <- data.frame(a)
#  temp <- list()
#  for (i in 1:t) {
#    if (type[i] == "C") {
#      vec <- numeric(com)
#      for (j in 1:com) vec[j] <- weighted.mean(x[, i], 
#                                               a[, j], na.rm = T)
#      temp[[i]] <- matrix(vec, com, 1, dimnames = list(com.names, 
#                                                       tr.names[i]))
#    }
#    else {
#      x[, i] <- as.factor(x[, i])
#      fac <- data.frame()
#      which.dom <- rep(NA, com)
#     # for (k in 1:com) {
#    #    temp2 <- tapply(a[, k], x[, i], sum)
#    #    fac <- rbind(fac, temp2)
#    #    which.dom[k] <- sample(levels(x[, i])[which(fac[k, 
#    #    ] == max(fac[k, ]))], size = 1)
#    #  }
#      colnames(fac) <- paste(tr.names[i], "_", levels(x[, 
#                                                        i]), sep = "")
#      rownames(fac) <- com.names
#    #  which.dom <- data.frame(which.dom)
#    #  colnames(which.dom) <- tr.names[i]
#    #  rownames(which.dom) <- com.names
#     # if (CWM.type == "dom") 
#    #    temp[[i]] <- which.dom
#    #  if (CWM.type == "all") 
#        temp[[i]] <- fac
#    }
#  }
#  temp <- data.frame(temp)
#  return(temp)
#}
#
#check_coverage <- function(Traits0, Abundances0, trait_names, trait_taxo, abundance_taxo) {
#  # Traits must be a LARGE data table with one column called 'taxo' (= species or whatever) and other colnames including trait_names
#  # Abundance must be a LONG data table with one column called 'taxo' (= species or whatever), one "value" (= abundance), one "Plot" and one "Year"
#  Traits <- Traits0
#  Abundances <- Abundances0
#  Traits[, taxo := .SD, .SDcols = trait_taxo]
#  Abundances[, taxo := .SD, .SDcols = abundance_taxo]
#  common_spe <- intersect(Traits$taxo, Abundances$taxo)
#
#  # n species in traits not in ab:
#  print(paste(length(unique(Traits$taxo[!(Traits$taxo %in% common_spe)])), "species in traits but not found in abundance", collapse = " "))
#  print(paste(length(unique(Abundances$taxo[!(Abundances$taxo %in% common_spe)])), "species in abundance but not found in traits", collapse = " "))
#
#  Traits2 <- Traits[, .SD, .SDcols = c("taxo", trait_names)]
#  Abundances2 <- Abundances[, .SD, .SDcols = c("Year", "taxo", "value", "Plot")]
#
#  traits_ab <- merge.data.table(Traits2, Abundances2, by = "taxo", all = TRUE)
#
#  traits_ab_cover <- traits_ab
#  traits_ab_cover <- traits_ab_cover[value > 0, ]
#  traits_ab_cover$ID <- 1:nrow(traits_ab_cover)
#  traits_ab_cover[, value2 := ifelse(is.na(taxo) | taxo == ' genus', 0, value)]
#  #MAX = traits_ab_cover[, list(MAX =  sum(value2, na.rm = T)/sum(value, na.rm = T)), by = c('Plot',"Year")]
#  
#  # For quantitative traits
#  traits_ab_cover2 <- cbind(
#    traits_ab_cover[, lapply(.SD, function(x) {
#      x = unlist(x)
#       x = ifelse(is.null(x) | (is.na(x) |  (x == "NA")), 0, 1)
#   #   x[is.null(x)] <- 0
#  #    x[(is.na(x) |  (x == "NA"))] <- 0
#  #    x[x != 0] <- 1
#   #   x <- as.numeric(x)
#      return(x)
#    }),
#    by = ID,
#    .SDcols = trait_names
#    ],
#    traits_ab_cover[, .SD, .SDcols = c("taxo", "Plot", "Year", "value", 'value2')]
#  )
#  
#  
#
#  final <- traits_ab_cover2[, lapply(.SD, function(x) {
#    sum(x * value2, na.rm = T)/sum(value, na.rm = T)
#  }),
#  .SDcols = trait_names, by = c('Plot', 'Year')][order(Plot), ]
#  return(final)
#}
#
#my_FD <- function(Traits, Abundances, trait_names) {
#  common_spe <- intersect(Traits$species_check, Abundances[Abundances$value > 0, ]$species_check)
#  Ab <- Abundances[species_check %in% common_spe, ][, list(value = sum(value)), by = c("Plot", "species_check")][order(species_check), ]
#  print(head(Ab))
#  Ab <- Ab[value > 0, ]
#  Ab <- dcast.data.table(Ab, Plot ~ species_check)
#  Traits <- Traits[Traits$species_check %in% colnames(Ab[, -1]), ][, (trait_names) := lapply(.SD, as.numeric), .SDcols = trait_names]
#  tt <- data.frame(Traits[order(species_check), ..trait_names])
#  tt <- mice::complete(mice(tt))
#  rownames(tt) <- sort(Traits$species_check)
#  print(nrow(tt))
#  print(ncol(Ab))
#  a <- dbFD(tt, Ab[, -1], corr = "lingoes", calc.CWM = F, calc.FDiv = T, m = 3)
#  return(a)
#}
#
#corr_env_cwm <- function(env_data, CWM, col_to_use) {
#  env_data <- env_data[env_data$Plot %in% CWM$Plot, ]
#  CWM <- CWM[CWM$Plot %in% env_data$Plot, ..col_to_use]
#  print(nrow(CWM))
#  print(nrow(env_data))
#
#  CWM_corr =  CWM[, lapply(.SD, function(x) {
#    mod <- lm(x ~ ., data = env_data[, c("Clay", "pH",'Soil.depth', "Precip.Annual", "Tmean.Annual", "TWI")])
#    print(length(residuals(mod)))
#    return(residuals(mod))
#  }), .SDcols = colnames(CWM)]
#  
#  return(CWM_corr)
#}
#
#my_dudi_pca <- function(CWM_matrix, col_to_use, NF = 3) {
#  matrix <- data.frame(CWM_matrix[complete.cases(CWM_matrix[, ..col_to_use])])
#  rownames(matrix) <- matrix$Plot
#  pca <- dudi.pca(matrix[, col_to_use], scannf = FALSE, nf = NF)
#  return(pca)
#}
#
#draw_dudi_mix = function(Data, axes, save = NA, new_var_names = NA){
#  library(ggrepel)
#eig = round(Data$eig/sum(Data$eig), 2)[axes]
#LI = Data$l1[axes]
#colnames(LI) = c('X', 'Y')
#CO = Data$c1[axes]*(eig*sum(Data$eig))^2
#colnames(CO) = c('X', 'Y')
#if (!length(new_var_names)>1){
#  print(rownames(CO))
#  new_var_names = rownames(CO)
#}
#gg = ggplot(CO, aes(xend = X, yend = Y, x = 0, y = 0)) + theme_minimal() +
#  geom_point(data =LI, aes(x = X, y = Y), color = 'grey50', alpha = 0.5) +
#  geom_hline(yintercept = 0, lwd = 0.2) + geom_vline(xintercept = 0, lwd = 0.2) +
#  geom_segment(arrow = arrow(length = unit(0.10, "inches"))) +
#  xlab(paste('FAMD', axes[1],  '=', eig[1]*100,'%')) +  ylab(paste('FAMD', axes[2],  '=', eig[2]*100,'%')) +
#  geom_text_repel(data = CO,  aes(x = 1.1*X, y = 1.1*Y, label = new_var_names))
#if (!is.na(save)){
#  ggsave(gg, file = save, width = 5, height = 5)
#}
#return(gg)
#}
#

