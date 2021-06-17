#####Calculate microbes richness#####

require(data.table)
require(phyloseq)

source("R/SCRIPTS UTILES/BE_plots_zero.R")

set.seed(1)

##############fungi dataset############
f11<-fread("Exploratories/Data/GRASSLANDS/FungiJuly2020/26470.txt")
length(unique(f11$Plotid))*length(unique(f11$OTU)) #miss zeros too, but too large to add
f14<-fread("Exploratories/Data/GRASSLANDS/FungiJuly2020/26471.txt")
f17<-fread("Exploratories/Data/GRASSLANDS/FungiJuly2020/26472.txt")
finfo<-fread("Exploratories/Data/GRASSLANDS/FungiJuly2020/26473.txt")

f11$Year<-2011; f11$DataID<-26470
f14$Year<-2014; f14$DataID<-26471
f17$Year<-2017; f17$DataID<-26472

soilf<-rbindlist(list(f11,f14,f17))
rm(f11,f14,f17)

soilf<-data.table(BEplotZeros(soilf,"Plotid",plotnam = "Plot"))
setnames(soilf,c("Plotid","OTU","Abundance"),c("Plot_bexis","Species","value"))
soilf$type<-"OTU_number"

#prepare OTU information data
finfo$Seq<-finfo$Kingdom<-NULL
summary(as.factor(finfo$Guild))
finfo<-finfo[OTU %in% unique(soilf$Species)] #remove forest species
finfo$Trophic_level<-finfo$Guild
finfo[Trophic_level %in% c("AMF","EMF","Lichen","EricoidM","OrchidM"),Trophic_level:="soilfungi.symbiont"]
finfo[Trophic_level %in% c("Saprotroph"),Trophic_level:="soilfungi.decomposer"]
finfo[Trophic_level %in% c("Pathogen","Parasite","Epiphyte"),Trophic_level:="soilfungi.pathotroph"]
finfo[Trophic_level %in% c("unknown","Endophyte"),Trophic_level:="soilfungi.other"]
summary(as.factor(finfo$Trophic_level))

finfo$Fun_group_broad<-finfo$Trophic_level
setnames(finfo,"Guild","Fun_group_fine")
setnames(finfo,"Phylum","Group_fine")
finfo$Group_broad<-"soilfungi"
finfo$Class<-finfo$Order<-finfo$Family<-finfo$Genus<-finfo$Species<-NULL

#merge
soilf<-merge(soilf,finfo,by.x="Species",by.y="OTU")


#####calculate richness (no rarefaction)
# rigroup<-soilf
# rigroup<-unique(rigroup[,.(value,Plot,Trophic_level,Species,Year)])
# rigroup$Species<-NULL
# set(rigroup,j="value",i=which(rigroup$value>0),value=1)
# rigroup<-dcast(rigroup,Plot~Trophic_level+Year,fun.aggregate=sum,value.var="value",na.rm=F)

#####calculate richness with rarefaction 2011######
fun11<-dcast.data.table(soilf[Year=="2011"],Species~Plot,value.var="value",fill=0)
sp11<-fun11$Species
fun11<-as.matrix(fun11[,-1])
rownames(fun11)<-sp11

spec11<-otu_table(fun11, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
phylo11<-phyloseq(spec11) # conversion step two for getting a phyloseq object

set.seed(1)
rarefied11 <- rarefy_even_depth(phylo11) #sanple size should be the smallest number of sequences per sample, alternatively, you can use 0.9*smallest number to also shuffle this sample a bit.
rarefied11 <- as.data.frame(rarefied11)


# rarefied11 <- melt.data.table(rarefied11,variable.name = "Plot")
# richrar11 <- rarefied11[,.(rich11=sum(value>1)),by=Plot]
# plot(richrar11$rich11~
#        (rigroup$soilfungi.pathotroph_2011+
#           rigroup$soilfungi.other_2011+
#           rigroup$soilfungi.pathotroph_2011+
#           rigroup$soilfungi.decomposer_2011)) #good correlation!

#melt and add trophic group
rarefied11$Species<-row.names(rarefied11)
rarefied11<-data.table(melt(rarefied11,variable.name = "Plot"))
rarefied11<-merge(rarefied11,finfo,by.x="Species",by.y="OTU")

#calculate richness
richrar11 <- rarefied11[,.(rich11=sum(value>1)),by=Plot]
rarefied11<-unique(rarefied11[,.(value,Plot,Trophic_level,Species)])
rarefied11$Species<-NULL
set(rarefied11,j="value",i=which(rarefied11$value>0),value=1)
rarefied11<-dcast(rarefied11,Plot~Trophic_level,fun.aggregate=sum,value.var="value",na.rm=F)
setnames(rarefied11,names(rarefied11)[2:ncol(rarefied11)],
         paste(colnames(rarefied11)[2:ncol(rarefied11)],"2011",sep="_"))

#####calculate richness with rarefaction 2014######
fun14<-dcast.data.table(soilf[Year=="2014"],Species~Plot,value.var="value",fill=0)
sp14<-fun14$Species
fun14<-as.matrix(fun14[,-1])
rownames(fun14)<-sp14

spec14<-otu_table(fun14, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
phylo14<-phyloseq(spec14) # conversion step two for getting a phyloseq object

set.seed(1)
rarefied14 <- rarefy_even_depth(phylo14) #sanple size should be the smallest number of sequences per sample, alternatively, you can use 0.9*smallest number to also shuffle this sample a bit.
rarefied14 <- as.data.frame(rarefied14)

#melt and add trophic group
rarefied14$Species<-row.names(rarefied14)
rarefied14<-data.table(melt(rarefied14,variable.name = "Plot"))
rarefied14<-merge(rarefied14,finfo,by.x="Species",by.y="OTU")

#calculate richness
richrar14 <- rarefied14[,.(rich14=sum(value>1)),by=Plot]
rarefied14<-unique(rarefied14[,.(value,Plot,Trophic_level,Species)])
rarefied14$Species<-NULL
set(rarefied14,j="value",i=which(rarefied14$value>0),value=1)
rarefied14<-dcast(rarefied14,Plot~Trophic_level,fun.aggregate=sum,value.var="value",na.rm=F)
setnames(rarefied14,names(rarefied14)[2:ncol(rarefied14)],
         paste(colnames(rarefied14)[2:ncol(rarefied14)],"2014",sep="_"))


#####calculate richness with rarefaction 2017######
fun17<-dcast.data.table(soilf[Year=="2017"],Species~Plot,value.var="value",fill=0)
sp17<-fun17$Species
fun17<-as.matrix(fun17[,-1])
rownames(fun17)<-sp17

spec17<-otu_table(fun17, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
phylo17<-phyloseq(spec17) # conversion step two for getting a phyloseq object

set.seed(1)
rarefied17 <- rarefy_even_depth(phylo17) #sanple size should be the smallest number of sequences per sample, alternatively, you can use 0.9*smallest number to also shuffle this sample a bit.
rarefied17 <- as.data.frame(rarefied17)

#melt and add trophic group
rarefied17$Species<-row.names(rarefied17)
rarefied17<-data.table(melt(rarefied17,variable.name = "Plot"))
rarefied17<-merge(rarefied17,finfo,by.x="Species",by.y="OTU")

#calculate richness
richrar17 <- rarefied17[,.(rich17=sum(value>1)),by=Plot]
rarefied17<-unique(rarefied17[,.(value,Plot,Trophic_level,Species)])
rarefied17$Species<-NULL
set(rarefied17,j="value",i=which(rarefied17$value>0),value=1)
rarefied17<-dcast(rarefied17,Plot~Trophic_level,fun.aggregate=sum,value.var="value",na.rm=F)
setnames(rarefied17,names(rarefied17)[2:ncol(rarefied17)],
         paste(colnames(rarefied17)[2:ncol(rarefied17)],"2017",sep="_"))

#########merge all
allrarfun<-merge(rarefied11,rarefied14,by="Plot")
allrarfun<-merge(allrarfun,rarefied17,by="Plot")

rm(list=ls()[! ls() %in% c("allrarfun")])

######################

source("R/SCRIPTS UTILES/BE_plots_zero.R")

######Protists#######
prot<-fread("Exploratories/Data/GRASSLANDS/200220_EP_species_diversity_GRL.txt")
spinfo<-fread("Exploratories/Data/GRASSLANDS/200220_EP_species_info_GRL.txt")
prot<-merge(prot,spinfo,by="Species")
rm(spinfo)
prot<-prot[Group_broad %in% "Protists"]

prot<-unique(prot[,.(value,Plot,Trophic_level,Species,Year)])
prot$Species<-NULL
set(prot,j="value",i=which(prot$value>0),value=1)
prot<-dcast(prot,Plot~Trophic_level+Year,fun.aggregate=sum,value.var="value",na.rm=F)

######################

######Bacteria (2011)#######
bac11<-fread("Exploratories/Data/GRASSLANDS/TEXTfiles/190319_Update/24866.txt")
bac11$Taxonomy<-NULL
bac11<-dcast.data.table(bac11,Sequence_variant~Plot_ID,value.var="Read_count",fill=0)
sp11<-bac11$Sequence_variant
bac11<-as.matrix(bac11[,Sequence_variant:=NULL])
rownames(bac11)<-sp11

spec11<-otu_table(bac11, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
phylo11<-phyloseq(spec11) # conversion step two for getting a phyloseq object

set.seed(1)
rarefied11 <- rarefy_even_depth(phylo11) #sanple size should be the smallest number of sequences per sample, alternatively, you can use 0.9*smallest number to also shuffle this sample a bit.
rarefied11 <- data.table(rarefied11)

rarefied11 <- melt.data.table(rarefied11,variable.name = "Plot")
richrar11 <- rarefied11[,.(rich11=sum(value>1)),by=Plot]


#compare with sikorski rarefaction
# load("Exploratories/Data/GRASSLANDS/bacterialRichness2011grassland_BExIS_ID_24866_new.Rdata")
# Richness2011$Plot0<-row.names(Richness2011)
# Richness2011$Plot0<-gsub("SG","SEG",Richness2011$Plot0)
# Richness2011$Plot0<-gsub("AG","AEG",Richness2011$Plot0)
# Richness2011$Plot0<-gsub("HG","HEG",Richness2011$Plot0)
# Richness2011<-data.table(BEplotNonZeros(Richness2011,"Plot0",plotnam = "Plot"))
# Richness2011<-merge(Richness2011,richrar11,by="Plot")
# 
# plot(Richness2011$rich11~Richness2011$Chao1)
# plot(Richness2011$rich11~Richness2011$Observed)
# plot(Richness2011$rich11~Richness2011$ACE) #not identical but well correlated, not due to seed

bac14<-fread("Exploratories/Data/GRASSLANDS/TEXTfiles/190319_Update/25066.txt")
bac14$Taxonomy<-NULL
bac14<-dcast.data.table(bac14,Sequence_variant~Plot_ID,value.var="Read_count",fill=0)
sp14<-bac14$Sequence_variant
bac14<-as.matrix(bac14[,Sequence_variant:=NULL])
rownames(bac14)<-sp14

spec14<-otu_table(bac14, taxa_are_rows = TRUE) # conversion step one for getting a phyloseq object
phylo14<-phyloseq(spec14) # conversion step two for getting a phyloseq object

set.seed(1)
rarefied14 <- rarefy_even_depth(phylo14) #sanple size should be the smallest number of sequences per sample, alternatively, you can use 0.9*smallest number to also shuffle this sample a bit.
rarefied14 <- data.table(rarefied14)

rarefied14 <- melt.data.table(rarefied14,variable.name = "Plot")
richrar14 <- rarefied14[,.(rich14=sum(value>1)),by=Plot]

richrar11<-data.table(BEplotZeros(richrar11,"Plot",plotnam = "Plot"))
richrar14<-data.table(BEplotZeros(richrar14,"Plot",plotnam = "Plot"))

bac<-merge(richrar11,richrar14,by="Plot",all=T)
setnames(bac,names(bac),c("Plot","bacteria_2011","bacteria_2017"))

rm(list=ls()[! ls() %in% c("allrarfun","prot","bac")])
######################

microbes<-merge(allrarfun,prot,by="Plot")
microbes<-merge(microbes,bac,by="Plot")

write.table(microbes,"Exploratories/Data/GRASSLANDS/200704_microbesRichness.txt",row.names=F)
