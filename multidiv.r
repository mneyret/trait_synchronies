# From: https://github.com/eric-allan/multidiversity/blob/master/README.md


#### a function to calculate multifunctionality
#### threshold can be a proportion of the maximum e.g. 0.5 or can be the median or mean
#### it can also be a vector of the same length as the number of diversities to
#### allow different thresholds for different diversities
#### threshold = FALSE calculates the mean of the scaled diversities or functions
#### scaling can be by any function specified in "sc", i.e. max, mean, sd etc., max is the default
#### if the maximum should be a mean of the top n values specify sc = "maxx", mx = n
#### centering by the mean is possible with cent = TRUE, to take z-scores of the diversities/processes, use sc="sd" and cent = TRUE 
#### "by" specifies a vector of the same length as nrow(x) to be used to split the data and use different
#### thresholds for the groups in "by"
#### "weights" allows different weightings for the different diversities/functions and should be a vector of the same length 
#### as the number of diversities/functions (shorter vectors will be recycled)
#### to weight a diversity as 0 and drop it from the calculation, code the weight as "NA"
#### the default is weights = 1: all diversities weighted the same

### a function to calculate max based on several values (e.g. a mean of top 5) the number of values is specified by mx
maxx <- function(x, n, ...){     
return(mean(rev(sort(x))[1:n], ...))}

### the multidiversity function
multidiv <- function(x, threshold=FALSE, sc = "max", mx = 1, cent = FALSE, by =FALSE, weights = 1){

result <- matrix(nrow=nrow(x), ncol=2)
weights <- rep_len(weights, length.out=ncol(x))   ### expand weights

### split the data on "by"
if(any(by!=FALSE)){

xs <- split(x, by)
xst <- list()

for(i in 1:length(unique(by))){
if(mx > 1){
xst[[i]] <- scale(xs[[i]], scale = apply(xs[[i]], 2, maxx, mx, na.rm=T), center = cent)}
else{
xst[[i]] <- scale(xs[[i]], scale = apply(xs[[i]], 2, match.fun(sc), na.rm=T), center = cent) 
}}
x.stand <- do.call("rbind", xst)
}

### otherwise standardise globally
else{
if(mx > 1){
x.stand <- scale(x, scale = apply(x, 2, maxx, mx, na.rm=T), center = cent)}
else{
x.stand <- scale(x, scale = apply(x, 2, match.fun(sc), na.rm=T), center = cent)
}}

#### sum diversities measured on each plot
#### NA weighted diversities are dropped from the calculation
x2 <- sweep(x, 2, weights, "*")  ## remove diversities with NA threshold from calc. of how many measured
gm <- apply(x2, 1, function(x)(sum(complete.cases(x))))


### no threshold: average standardised values
if(FALSE %in% threshold){  ### prevent error message if threshold is vector
x.stand2 <- sweep(x.stand, 2, weights, "*")
m <- apply(x.stand2, 1, mean, na.rm=T)
}

else{

if (any(c("median","mean") %in% threshold)){  ### mean or median threshold

tf <- match.fun(threshold)

x.thresh <- apply(x.stand, 2, function(x)(1*(x > tf(x, na.rm=T))))
weights2 <- matrix(rep(weights, nrow(x.thresh)), nrow=nrow(x.thresh), byrow=TRUE)
x.thresh <- x.thresh*weights2   ### multiply by weights
gg <- apply(x.thresh, 1, sum, na.rm=T)
m <- gg/gm

}

else{ 

x.thresh <- 1*sweep(x.stand, 2, threshold, ">")  ### does each variable pass threshold?

weights2 <- matrix(rep(weights, nrow(x.thresh)), nrow=nrow(x.thresh), byrow=TRUE)
x.thresh <- x.thresh*weights2   ### multiply by weights

gg <- apply(x.thresh, 1, sum,na.rm=T)
m <- gg/gm
}}

result[, 1] <- m
result[, 2] <- gm

colnames(result)<-c("m", "groups measured")

return(result)
}