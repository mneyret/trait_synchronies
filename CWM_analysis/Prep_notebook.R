# Library loading; creation of functions; parameters
library(readxl)
library(data.table)
library(ade4)
library(factoextra)
library(ggfortify)
library(ggcorrplot)
library(lavaan)
library(semPlot)
library(ghibli)
library(mice)
library(kableExtra)

setwd("/Users/Margot/Desktop/Research/Senckenberg/Project_Ecosystem_strat/Analysis/")

#### Functions ####
# Function to correct for region and/or environment
corrections <- function(CWM_data, Env_data = NA, env_corr = FALSE, region_corr = FALSE, variables = NA, env_variables = NA) {
  CWM_data_copy <- copy(CWM_data)
  if (env_corr == TRUE & is.na(env_data)) {
    stop("Please provide environmental data")
  }
  
  if (is.na(variables)) {
    variables <- colnames(CWM_data_copy)[!(colnames(CWM_data_copy) %in% c("Plot", "Year"))]
  }
  
  if (is.na(env_variables)) {
    env_variables <- colnames(Env_data)[!(colnames(Env_data) %in% c("Plot", "Year"))]
  }
  
  print(region_corr)
  if (region_corr == TRUE) {
    CWM_data_copy[, (variables) := lapply(.SD, function(x) {
      if (region_corr == TRUE) {
        mod <- lm(x ~ substr(Plot, 1, 1))
        return(residuals(mod))
      }
    }),
    .SDcols = variables
    ]
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

# Function to change the direction of PCA axis (for reproductibility of axes direction)
axes_direction_pca <- function(pca, cwm, traitnames, conditions) {
  if (length(conditions) != length(traitnames)) {
    print("Number of axes, traits, and conditions should be equal")
  }
  else {
    cond_positive <- ifelse(gsub(" ", "", conditions) == ">0", TRUE, FALSE)
    I <- length(traitnames)
    for (i in 1:I) {
      if ((pca$li[cwm[get(traitnames[i]) == max(get(traitnames[i])), Plot], i] > 0) != cond_positive[i]) {
        #   print(paste("Changed direction of axis", i))
        pca$li[, i] <- -pca$li[, i]
        pca$co[, i] <- -pca$co[, i]
      }
    }
  }
  return(pca)
}

# Function to draw pca
my_dudi_pca <- function(CWM_matrix, col_to_use, NF = 3) {
  matrix <- data.frame(CWM_matrix[complete.cases(CWM_matrix[, ..col_to_use])])
  rownames(matrix) <- matrix$Plot
  pca <- dudi.pca(matrix[, col_to_use], scannf = FALSE, nf = NF)
  return(pca)
}

