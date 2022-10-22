# code consists of 6 main parts
#1)loading/ installing packages
#2) analysis of the wing length relation with latitude
#3) analysis of phylogenetic patterns
#4) analysis of the wing length  relation with mean annual temperature
#5) analysis of the relationship between temperature/ latitude and the wing length in individual genera
#6) final plotting of figures+addendum for the additional supplementary wing

# to run the code, please put all the files "main.txt", "df_temp.txt", "cranston_names.csv"
# as well as files c_cadi_dna, c_cadiv_dna, c_coi_dna, c_18s_dna and c_s28_dna into your working directory
#additionaly put the files c_cadi_dna, c_cadiv_dna, c_coi_dna, c_18s_dna and c_s28_dna into the
# apex package folder (after installing the package) to run the analysis
# warning messages while running the tidying on the mblm results arises while 
#mblm package is not fully supported by the newer  broom package version, but we can not find any output difference

#1)LOADING/ INSTALLING PACKAGES
#Packages installation
#install.packages(usdm)
#install.packages(adegenet)
#install.packages(phytools)
#install.packages(biomaRt)
#install.packages(rentrez)
#install.packages(biomartr)
#install.packages(msa)
#install.packages(devtools)
#install.packages(tidyverse)
#install.packages(stats)
#install.packages(fitdistrplus)
#install.packages(devtools)
#install.packages(tidyverse)
#install.packages(phangorn)
#install.packages(seqinr)
#install.packages(geiger)
#install.packages(apex)
#install.packages(ncdf4)
#install.packages(RColorBrewer)
#install.packages(ggplot2)
#install.packages(raster)
#install.packages(sp)
#install.packages(tidyverse)
#install.packages(ggmap)
#install.packages(dplyr)
#install.packages(ggeffects)
#install.packages(raster)
#install.packages(sp)
#install.packages(rgdal)
#install.packages("diversitree")
#install.packages("caper")
#install.packages(plyr)
#install.packages(dplyr)
#install.packages(reshape2)
#install.packages(nlme)
#install.packages(visreg)
#install.packages(car)
#install.packages(ggplot2)
#install.packages(gridExtra)
#install.packages(FD)
#install.packages(mblm)
#install.packages(trend)
#install.packages(broom)
#install.packages(ggfortify)
#install.packages(DataCombine)
#install.packages(graphics)
#install.packages(Kendall)
#install.packages(boot)
#install.packages(modifiedmk)
#install.packages(vegan)
#install.packages(gtools)
#install.packages(vegan)
#install.packages(codyn)
#install.packages("lmeInfo")
#install.packages(effectsize)
#install.packages(sjPlot)
#install.packages(sjmisc)
#install.packages(sjlabelled)
#install.packages(msa)
#install.packages(nlme)
#install.packages(stats)
#install.packages(mblm)
#install.packages(trend)
#install.packages(ggpubr)
#install.packages("rgeos")
#install.packages("rnaturalearth")
#install.packages("sf")
#install.packages("rnaturalearthdata")
#install.packages(mapview)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.14")
#install.packages("adegenet", dep=TRUE)
#install.packages ("biomartr")
#install.packages ("installr")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#("BiocManager")
#BiocManager::install("biomaRt")
#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("msa")

# load packages
require(usdm)
require(adegenet)
require(phytools)
library(biomaRt)
library (rentrez)
library (biomartr)
library(msa)
library (tidyverse)
library(stats)
library(fitdistrplus)
library(devtools)
library(phangorn)
library(seqinr)
library(geiger)
library(apex)
library(ncdf4)
library(RColorBrewer)
require(ggplot2)
library(raster)
library(sp)
library(ggmap)
require(dplyr)
require(ggeffects)
library(raster)
library(sp)
library(rgdal)
library("diversitree")
library("caper")
require(plyr)
require(dplyr)
require(reshape2)
require(nlme)
require(visreg)
require(car)
require(gridExtra)
require(FD)
require(mblm)
require(trend)
require(broom)
require(ggfortify)
library(DataCombine)
require(graphics)
require(Kendall)
require(boot)
require(modifiedmk)
require(vegan)
library(gtools)
require(codyn)
require("lmeInfo")
require(effectsize)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
require(nlme)
library(ggpubr)
library("rgeos")
library("rnaturalearth")
library("rnaturalearthdata")
library(mapview)
require(sf)

#########################################

# 2) ANALYSIS OF THE WING LENGTH RELATION WITH LATITUDE
#load main dataset
size= read.delim("main.txt",sep="\t",header=TRUE,fill = TRUE)
binnedSamples <- cut(size$Lat_un, breaks = c(-50,-40,-30,-20,-10,0, 10, 20, 30, 40, 50, 60, 70,80) )#intervals of coordinates
size$interval=binnedSamples
binnedSamples1 <- cut(size$Lat1_cor, breaks = c(0, 10, 20, 30, 40, 50, 60, 70,80) )#intervals of coordinates

size$interval_simplified=binnedSamples1#interval for latitude with - removed to collapse N and S hemisphere, they are coded as separate variable "hemisphere"

#split the data set for the non-simplified latitude and get average wing lengths per interval
dt2=aggregate(size[, 7], list(size$interval),median)#median wing length per lat interval
vector1<-c(-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70)
dt2$lat2=vector1 #replace factorial interval with numeric interval

#remove genera with less then 60 measurments (specimens measured) from the analysis 
#[otherwise too many levels of factor for gls models to handle]
b=table(size$Group.4)#count number of measurments per genus
b=as.data.frame(b)
bx=b[order(-b$Freq), ]#order genera by the number of specimens
b$Group.4=b$Var1#rename column in dataset to match the name of the genus column in the main datset
size1=merge(size,b,by = c("Group.4")) #merge table with number of specimes per genus with main dataset
size1=size1[order(-size1$Freq), ] #order the table by number of specimens in the genus
size2=size1[c(1:2643),]#remove genera with less than 60 observations

#plot fig 3
pex1=ggplot(size2,aes(x=Lat1_cor,y=wl,colour=as.factor(hemisphere)))+geom_smooth(method="lm",se=.99,lwd=1.5)#temperature plot
pex1=pex1+xlab("Geographic latitude [Â°]")+ylab("Wing length [Âµm]")
pex1=pex1+theme_classic()
pex1=pex1+theme(legend.position = "none")
pex1=pex1+scale_color_manual(values = c("#FF9933", "#000000"))
pex1=pex1+theme(axis.line = element_line(size = 2, linetype=1,colour = "grey"))
pex1=pex1 + theme(axis.ticks = element_blank())
pex1_1=pex1+ theme(text = element_text(size = 30))

# plot same with all data points, supplementary figure 2A
pex3=ggplot(size2,aes(x=Lat1_cor,y=wl,colour=as.factor(hemisphere)))+geom_smooth(method="lm",se=.99,lwd=1.5)+geom_point()#temperature plot
pex3=pex3+xlab("Geographic latitude [Â°]")+ylab("Wing length [Âµm]")
pex3=pex3+theme_classic()
pex3=pex3+theme(legend.position = "none")
pex3=pex3+scale_color_manual(values = c("#FF9933", "#000000"))
pex3=pex3+theme(axis.line = element_line(size = 2, linetype=1,colour = "grey"))
pex3=pex3 + theme(axis.ticks = element_blank())
pex3=pex3+ theme(text = element_text(size = 30))

# main gls model
md1<- gls(wl~Lat1_cor+hemisphere+Group.4+Lat1_cor*Group.4+Lat1_cor*hemisphere,data=size2, correlation=corAR1(), method="ML",na.action = na.omit)
mod.gls.1.1 <- update(md1, correlation=corARMA(p=1))#best fit
#we also have tested this models, you can run them, but mod.gls.3 can take several hours to update
#mod.gls.3 <- update(md1, correlation=corARMA(p=3))
#mod.gls.0 <- update(md1, correlation=NULL)
#anova(md1, mod.gls.3, mod.gls.1.1,mod.gls.0) # AR(2) vs AR(1)
an2=Anova(mod.gls.1.1)#anova on the best fit model with latitude
plot(density(resid(md1, type='pearson')))
lines(density(resid(mod.gls.1.1, type='pearson')), col='red')

#testing git commits

#########################################
# 3) ANALYSIS OF PHYLOGENETIC PATTERNS

#Phylogenetic analysis based on the sequences published in:
#Cranston, P. S., Hardy, N. B., & Morse, G. E. (2012). A dated molecular phylogeny for the Chironomidae (Diptera). Systematic Entomology, 37(1), 172-188.


#put the files  c_cadi_dna, c_cadiv_dna, c_coi_dna, c_18s_dna and c_s28_dna into the
# apex package folder to run the analysis
files <- dir(system.file(package="apex"),patter="c_", full=FALSE)# this will change on your computer base on directory
files 
#reading the 5 abovemntioned files as multiFASTA object
x <- read.multiFASTA(files, add.gaps=TRUE)
names(x@dna) # dispplay names of the genes
#concatenate 5 disparate genes into a single file for the phylogenetic analysis
y <- concatenate(x)
#gett the tree from 5 genes
trees <- getTree(x)
#get the tree from concatenated data
al_dna_phyDat <- phyDat(y,type = "DNA", levels = NULL)
#calculate distance matrix
dna_dist <- dist.ml(al_dna_phyDat, model="JC69")
#build neighbor joining tree
dna_NJ  <- NJ(dna_dist)
# find optimal tree; might run for a while, hours to a day or two depending on a PC
fit <- pml(dna_NJ, al_dna_phyDat,maxit=50)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic",maxit=50)
#bootstrap the fitJC; might run for a while, hours to a day or two depending on a PC
bs <- bootstrap.pml(fitJC, bs=100, optNni=FALSE, multicore=FALSE, control = pml.control(trace=0))
#consensus tree from bootstrap
cnet <- consensusNet(bs, .5)
#root first tree
x=root(cnet, 158, r = TRUE)
x=drop.tip(x, "Conochironomus  sp.", root.edge = 0) #remove conochironomus, due to insufficient molecular data
#upload size data for the taxa in the tree

data=read.csv("cranston_names.csv",sep=",")
size=data[,c(2,6)]
#break size into the intervals for the phylogenetic analysis
binnedSamples <- cut(data$wl, breaks = c(0,1000,2000,3000,4000,5000,6000) )#intervals of coordinates
data$interval.w=binnedSamples
size=data[,c(2,14)]

size$int=as.numeric(size$interval.w)

size$int[is.na(size$int)] = 0


#list with interval brakes
test <- setNames(size$int, size$name)#generate named vector required for contmp
body.size<-as.data.frame(data)[,6] #choosing the 3rd column of dataset X2#


#create a reconstruction of ancestral character state for size
obj<-contMap(x,test,plot=FALSE)
#plot reconstruction on a tree
tr=plot(obj,type="phylogram",leg.txt="wing length (mm)",lwd=2,
        mar=c(4,2,4,2))
obj1=obj
drop.tip.contMap(obj1, tip)
tr1=plot(obj,type="fan",leg.txt="wing length (mm)",lwd=2,
        mar=c(4,2,4,2),tiplabels=FALSE)
#write.tree(x, file="100bs_tr19042021.tre") #save the tree
#adjust the tree plot
title(main="Chironomidae phylogenetic tree")
axis(1)
title(xlab="Time from the root")
#calculate correlation between the trait distribution (size) and phylogenetic structure
phylosig(x,body.size,method="lambda",test=TRUE)
phenogram(x, test, spread.labels = TRUE)
fancyTree(x, type = "phenogram95", x = type)
#########################################
#4 ANALYSIS OF THE WING LENGTH RELATION WITH MEAN ANNUAL TEMPERATURE
#load dataset with bioclimatic variables including mean annual temperature. Bioclimatic vars. from https://worldclim.org/
size_temp= read.delim("df_temp.txt",sep="\t",header=TRUE,fill = TRUE)
#create a coordinate object
my.sf.point <- st_as_sf(x = size_temp,
                        coords = c("Longitude","Lat_un"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#subset only extant taxa, removing some fossils contained in dataset
size_temp=subset(size_temp,Geological.age=="rec")
#remove genera with less then 60 measurements
b1=table(size_temp$Group.4)
b1=as.data.frame(b1)
by=b1[order(-b1$Freq), ]
b1$Group.4=b1$Var1
size_temp1=merge(size_temp,b1,by = c("Group.4"))
size_temp1=size_temp1[order(-size_temp1$Freq), ]
size_temp1=size_temp1[c(1:2643),]#with more than 60 obs
#as mean annual temperature is given in worldclim data as Â°C*10 it has to be adjusted for the calculations
size_temp1$temp=size_temp1$Temp/10 

#plot wing length vs mean annual temperature
pex=ggplot(size_temp1,aes(x=temp,y=wl,colour=as.factor(hemisphere)))+geom_smooth(method="lm",se=.99,lwd=1.5)#temperature plot
pex=pex+xlab("Annual Mean Temperature [CÂ°]")+ylab("Wing length [Âµm]")
pex=pex+theme_classic()
pex=pex+theme(legend.position = "none")
pex=pex+scale_color_manual(values = c("#FF9933", "#000000"))
pex=pex+theme(axis.line = element_line(size = 2, linetype=1,colour = "grey"))
pex=pex + theme(axis.ticks = element_blank())
pex=pex+ theme(text = element_text(size = 30))

#Supplementary figure 2B
pex4=ggplot(size_temp1,aes(x=temp,y=wl,colour=as.factor(hemisphere)))+geom_smooth(method="lm",se=.99,lwd=1.5)+geom_point()#temperature plot
pex4=pex4+xlab("Annual Mean Temperature [CÂ°]")+ylab("Wing length [Âµm]")
pex4=pex4+theme_classic()
pex4=pex4+theme(legend.position = "none")
pex4=pex4+scale_color_manual(values = c("#FF9933", "#000000"))
pex4=pex4+theme(axis.line = element_line(size = 2, linetype=1,colour = "grey"))
pex4=pex4 + theme(axis.ticks = element_blank())
pex4=pex4+ theme(text = element_text(size = 30))

#intercation between temp and genera removed due to the singular fit of the high level random effect
md1_1<- gls(wl~temp+hemisphere+Group.4+temp*hemisphere,data=size_temp1, correlation=corAR1(), method="ML",na.action = na.omit)
mod.gls.1.1x <- update(md1_1, correlation=corARMA(p=1))#best fit

an3=Anova(mod.gls.1.1x)#anova on the best fit model with latitude

#plot the distribution of the specimens from dataset on the map of the world
base_world <- map_data("world")
#Figure 2
sizep=ggplot() +
  geom_polygon(data=base_world, aes(x=long, y=lat, group=group)) +
  geom_point(data=size_temp, aes(x=Longitude, y=Lat_un, colour=wl), size=3, alpha=I(0.7))+xlab("longitude,Â°")+ylab("latitude,Â°")

sizep=sizep+scale_color_gradientn(colours = rainbow(5))+theme_classic()
sizep=sizep+ theme(text = element_text(size = 30))

#########################################
#5 ANALYSIS OF THE RELATIONSHIP BETWEEN TEMPERATURE/ LATITUDE AND THE WING LENGTH IN INDIVIDUAL GENERA
#We extract the data for the separate genera, with over 60 records, and conductin
#a analysis of wing length vs temperature/ latitude for each of them
#A) wing vs latitude
theme_set(theme_pubr())#theme for the images
#rename the dataset with temperature and latitude for the analysis below
size3=size_temp1
#median based linear model for the wing size vs latitude for all the genera toogether
mb1=mblm(wl~Lat1_cor,data=size3)
summary (mb1)
#extract individual genera as separate datasets

mk<- with(size3,
            by(size3, Group.4,
               function(x) mk.test(size3$wl)))

tidy<- with(size3,
          by(size3, Group.4,
             function(x) mblm(wl~Lat1_cor,data=x)))

sapply(tidy,coef)


Cricotopus=subset(size3,Group.4=="Cricotopus")#1
Cricotopus=Cricotopus[order(-Cricotopus$Lat1_cor), ]#order by latitude (normalized) to allow mk test on the wl~lat
Dicrotendipes=subset(size3,Group.4=="Dicrotendipes")#2
Dicrotendipes=Dicrotendipes[order(-Dicrotendipes$Lat1_cor), ]
Micropsectra=subset(size3,Group.4=="Micropsectra")#3
Micropsectra=Micropsectra[order(-Micropsectra$Lat1_cor), ]
Tanytarsus=subset(size3,Group.4=="Tanytarsus")#4
Tanytarsus=Tanytarsus[order(-Tanytarsus$Lat1_cor), ]
Chironomus=subset(size3,Group.4=="Chironomus")#5
Chironomus=Chironomus[order(-Chironomus$Lat1_cor), ]
Stenochironomus=subset(size3,Group.4=="Stenochironomus")#6
Stenochironomus=Stenochironomus[order(-Stenochironomus$Lat1_cor), ]
Limnophyes=subset(size3,Group.4=="Limnophyes")#7
Limnophyes=Limnophyes[order(-Limnophyes$Lat1_cor), ]
Paratanytarsus=subset(size3,Group.4=="Paratanytarsus")#8
Paratanytarsus=Paratanytarsus[order(-Paratanytarsus$Lat1_cor), ]
Parachironomus=subset(size3,Group.4=="Parachironomus")#9
Parachironomus=Parachironomus[order(-Parachironomus$Lat1_cor), ]
Bryophaenocladius=subset(size3,Group.4=="Bryophaenocladius")#10
Bryophaenocladius=Bryophaenocladius[order(-Bryophaenocladius$Lat1_cor), ]
Rheotanytarsus=subset(size3,Group.4=="Rheotanytarsus")#11
Rheotanytarsus=Rheotanytarsus[order(-Rheotanytarsus$Lat1_cor), ]
Cladotanytarsus=subset(size3,Group.4=="Cladotanytarsus")#12
Cladotanytarsus=Cladotanytarsus[order(-Cladotanytarsus$Lat1_cor), ]
Conchapelopia=subset(size3,Group.4=="Conchapelopia")#13
Conchapelopia=Conchapelopia[order(-Conchapelopia$Lat1_cor), ]
Chaetocladius=subset(size3,Group.4=="Chaetocladius")#14
Chaetocladius=Chaetocladius[order(-Chaetocladius$Lat1_cor), ]
Stictochironomus=subset(size3,Group.4=="Stictochironomus")#15
Stictochironomus=Stictochironomus[order(-Stictochironomus$Lat1_cor), ]
Pseudosmittia=subset(size3,Group.4=="Pseudosmittia")#16
Pseudosmittia=Pseudosmittia[order(-Pseudosmittia$Lat1_cor), ]
Tanypus=subset(size3,Group.4=="Tanypus")#17
Tanypus=Tanypus[order(-Tanypus$Lat1_cor), ]
Polypedilum=subset(size3,Group.4=="Polypedilum")#18
Polypedilum=Polypedilum[order(-Polypedilum$Lat1_cor), ]
Zavrelimyia=subset(size3,Group.4=="Zavrelimyia")#19
Zavrelimyia=Zavrelimyia[order(-Zavrelimyia$Lat1_cor), ]
Ablabesmyia=subset(size3,Group.4=="Ablabesmyia")#20
Ablabesmyia=Ablabesmyia[order(-Ablabesmyia$Lat1_cor), ]
Corynoneura=subset(size3,Group.4=="Corynoneura")#21
Corynoneura=Corynoneura[order(-Corynoneura$Lat1_cor), ]
Orthocladius=subset(size3,Group.4=="Orthocladius")#22
Orthocladius=Orthocladius[order(-Orthocladius$Lat1_cor), ]
Eukiefferiella=subset(size3,Group.4=="Eukiefferiella")#23
Eukiefferiella=Eukiefferiella[order(-Eukiefferiella$Lat1_cor), ]
Cladopelma=subset(size3,Group.4=="Cladopelma")#24
Cladopelma=Cladopelma[order(-Cladopelma$Lat1_cor), ]

#Calculate Theil-Sen's slopes of the wing vs latitude for each genus

mk_Chironomus=tidy(mk.test(Chironomus$wl))
tidy_Chironomus=tidy(mblm(wl~Lat1_cor,data=Chironomus))

mk_Cricotopus=tidy(mk.test(Cricotopus$wl))
tidy_Cricotopus=tidy(mblm(wl~Lat1_cor,data=Cricotopus))

mk_Dicrotendipes=tidy(mk.test(Dicrotendipes$wl))
tidy_Dicrotendipes=tidy(mblm(wl~Lat1_cor,data=Dicrotendipes))

mk_Micropsectra=tidy(mk.test(Micropsectra$wl))
tidy_Micropsectra=tidy(mblm(wl~Lat1_cor,data=Micropsectra))

mk_Tanytarsus=tidy(mk.test(Tanytarsus$wl))
tidy_Tanytarsus=tidy(mblm(wl~Lat1_cor,data=Tanytarsus))

mk_Stenochironomus=tidy(mk.test(Stenochironomus$wl))
tidy_Stenochironomus=tidy(mblm(wl~Lat1_cor,data=Stenochironomus))

mk_Paratanytarsus=tidy(mk.test(Paratanytarsus$wl))
tidy_Paratanytarsus=tidy(mblm(wl~Lat1_cor,data=Paratanytarsus))

mk_Conchapelopia=tidy(mk.test(Conchapelopia$wl))
tidy_Conchapelopia=tidy(mblm(wl~Lat1_cor,data=Conchapelopia))

mk_Bryophaenocladius=tidy(mk.test(Bryophaenocladius$wl))
tidy_Bryophaenocladius=tidy(mblm(wl~Lat1_cor,data=Bryophaenocladius))

mk_Parachironomus=tidy(mk.test(Parachironomus$wl))
tidy_Parachironomus=tidy(mblm(wl~Lat1_cor,data=Parachironomus))


mk_Cladotanytarsus=tidy(mk.test(Cladotanytarsus$wl))
tidy_Cladotanytarsus=tidy(mblm(wl~Lat1_cor,data=Cladotanytarsus))

mk_Chaetocladius=tidy(mk.test(Chaetocladius$wl))
tidy_Chaetocladius=tidy(mblm(wl~Lat1_cor,data=Chaetocladius))


mk_Limnophyes=tidy(mk.test(Limnophyes$wl))
tidy_Limnophyes=tidy(mblm(wl~Lat1_cor,data=Limnophyes))

mk_Rheotanytarsus=tidy(mk.test(Rheotanytarsus$wl))
tidy_Rheotanytarsus=tidy(mblm(wl~Lat1_cor,data=Rheotanytarsus))

mk_Ablabesmyia=tidy(mk.test(Ablabesmyia$wl))
tidy_Ablabesmyia=tidy(mblm(wl~Lat1_cor,data=Ablabesmyia))


mk_Tanypus=tidy(mk.test(Tanypus$wl))
tidy_Tanypus=tidy(mblm(wl~Lat1_cor,data=Tanypus))

mk_Stictochironomus=tidy(mk.test(Stictochironomus$wl))
tidy_Stictochironomus=tidy(mblm(wl~Lat1_cor,data=Stictochironomus))


mk_Cladopelma=tidy(mk.test(Cladopelma$wl))
tidy_Cladopelma=tidy(mblm(wl~Lat1_cor,data=Cladopelma))

mk_Corynoneura=tidy(mk.test(Corynoneura$wl))
tidy_Corynoneura=tidy(mblm(wl~Lat1_cor,data=Corynoneura))

mk_Zavrelimyia=tidy(mk.test(Zavrelimyia$wl))
tidy_Zavrelimyia=tidy(mblm(wl~Lat1_cor,data=Zavrelimyia))


mk_Eukiefferiella=tidy(mk.test(Eukiefferiella$wl))
tidy_Eukiefferiella=tidy(mblm(wl~Lat1_cor,data=Eukiefferiella))


mk_Orthocladius=tidy(mk.test(Orthocladius$wl))
tidy_Orthocladius=tidy(mblm(wl~Lat1_cor,data=Orthocladius))

mk_Pseudosmittia=tidy(mk.test(Pseudosmittia$wl))
tidy_Pseudosmittia=tidy(mblm(wl~Lat1_cor,data=Pseudosmittia))
mk_Polypedilum=tidy(mk.test(Polypedilum$wl))
tidy_Polypedilum=tidy(mblm(wl~Lat1_cor,data=Polypedilum))

#bind output of Theil-Sen's slope analysis into a single dataset

tidyy=bind_rows(
  tidy_Cricotopus	,
  tidy_Dicrotendipes	,
  tidy_Micropsectra	,
  tidy_Tanytarsus	,
  tidy_Chironomus	,
  tidy_Stenochironomus	,
  tidy_Limnophyes	,
  tidy_Paratanytarsus	,
  tidy_Parachironomus	,
  tidy_Bryophaenocladius,
  tidy_Rheotanytarsus	,
  tidy_Cladotanytarsus	,
  tidy_Conchapelopia	,
  tidy_Chaetocladius	,
  tidy_Stictochironomus	,
  tidy_Pseudosmittia,
  tidy_Tanypus	,
  tidy_Polypedilum,
  tidy_Zavrelimyia	,
  tidy_Ablabesmyia	,
  tidy_Corynoneura	,
  tidy_Orthocladius,
  tidy_Eukiefferiella	,
  tidy_Cladopelma)


name=c( "Cricotopus",       "Dicrotendipes" ,    "Micropsectra" ,     "Tanytarsus"  ,
        "Chironomus" ,       "Stenochironomus" ,  "Limnophyes"   ,     "Paratanytarsus" ,
        "Parachironomus" ,   "Bryophaenocladius", "Rheotanytarsus",    "Cladotanytarsus",
        "Conchapelopia",     "Chaetocladius",     "Stictochironomus",  "Pseudosmittia",
        "Tanypus" ,          "Polypedilum"  ,     "Zavrelimyia",       "Ablabesmyia",
        "Corynoneura" ,      "Orthocladius" ,     "Eukiefferiella" ,   "Cladopelma"  )
#add names of the genera to the dataset with the Theil-Sen's slope output
#join the names to the outputs
tidyy=subset(tidyy,term=="Lat1_cor" )
df1=cbind(name,tidyy)
#bind Mann-Kendall analysis output into a table test 
tidmk=bind_rows(
  mk_Cricotopus	,
  mk_Dicrotendipes	,
  mk_Micropsectra	,
  mk_Tanytarsus	,
  mk_Chironomus	,
  mk_Stenochironomus	,
  mk_Limnophyes	,
  mk_Paratanytarsus	,
  mk_Parachironomus	,
  mk_Bryophaenocladius,
  mk_Rheotanytarsus	,
  mk_Cladotanytarsus	,
  mk_Conchapelopia	,
  mk_Chaetocladius	,
  mk_Stictochironomus	,
  mk_Pseudosmittia,
  mk_Tanypus	,
  mk_Polypedilum,
  mk_Zavrelimyia	,
  mk_Ablabesmyia	,
  mk_Corynoneura	,
  mk_Orthocladius,
  mk_Eukiefferiella	,
  mk_Cladopelma)
df2=cbind(name,tidmk)

#sort significant relationship wing length vs temp, so we can plot significant and 
#non-significant relationships

a <- df1 %>%
  filter( p.value>0.05) %>%
  mutate(final = ifelse(p.value >0.05, "insignificant"))
b <- df1 %>%
  filter( p.value<=0.05) %>%
  mutate(final = ifelse(p.value<=0.05, "significant"))
df3=rbind(a,b)

#Plot the Thei-Sen's slopes for each genus analyzed
g3=ggplot(df3, aes(as.factor(name),estimate))+ geom_point(size=5,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("[Âµm per Â° of latitude]")
g3=g3+ theme(axis.text.y = element_text(face = "italic"))
g3=g3+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

#calculate rate of wing size change based on the mean size of the genus representatives
size4=aggregate(size3[, 19], list(size3$Group.4),mean)
size4$name=size4$Group.1
size_change<-merge(df3,size4,by = c("name"), all.x = TRUE, all.y = TRUE)
size_change$mean_wl=size_change$x
summary(lm(size_change$estimate~size_change$x))
abline(lm(size_change$estimate~size_change$x))
size4$Group.4=size4$Group.1
size5=size4
#rate of wing size change vs mean wing size
size_slope<-merge(size5,df3,by = c("name"), all.x = TRUE, all.y = TRUE)
#model of the rate of wing size change based on the mean size of the genus representatives
mdx<- gls(estimate~x,data=size_slope, correlation=corAR1(), method="ML",na.action = na.omit)
mod.x1 <- update(mdx, correlation=corARMA(p=1))
an1=Anova(mod.x1)
#plot rate of wing size change based on the mean size of the genus representatives
pex_2=ggplot(size_slope,aes(x=x,y=estimate))+geom_point(aes(size = 2))+geom_smooth(method="lm",se=FALSE,lwd=1.5, colour = "#FF9966")+xlab("wing length, Ã‚Âµm")+ylab("Sen's slope of the wing lenght shift")#temperature plot
pex_2=pex_2+xlab("Wing length[Âµm]")+ylab("Theil-Sen's slope, wing lenght")#temperature plot
pex_2=pex_2+theme_classic()
pex_2=pex_2+theme(legend.position = "none")
pex_2=pex_2+theme(axis.line = element_line(size = 2, linetype=1,colour = "grey"))
pex_2=pex_2 + theme(axis.ticks = element_blank())
pex_2=pex_2+ theme(text = element_text(size = 30))

#B Wing vs mean annual temperature (code same as above,
#just uses temperature instead of latitude as a predictor)
Cricotopus=subset(size3,Group.4=="Cricotopus")#1
Cricotopus=Cricotopus[order(-Cricotopus$temp), ]#order by latitude (normalized) to allow mk test on the wl~lat
Dicrotendipes=subset(size3,Group.4=="Dicrotendipes")#2
Dicrotendipes=Dicrotendipes[order(-Dicrotendipes$temp), ]
Micropsectra=subset(size3,Group.4=="Micropsectra")#3
Micropsectra=Micropsectra[order(-Micropsectra$temp), ]
Tanytarsus=subset(size3,Group.4=="Tanytarsus")#4
Tanytarsus=Tanytarsus[order(-Tanytarsus$temp), ]
Chironomus=subset(size3,Group.4=="Chironomus")#5
Chironomus=Chironomus[order(-Chironomus$temp), ]
Stenochironomus=subset(size3,Group.4=="Stenochironomus")#6
Stenochironomus=Stenochironomus[order(-Stenochironomus$temp), ]
Limnophyes=subset(size3,Group.4=="Limnophyes")#7
Limnophyes=Limnophyes[order(-Limnophyes$temp), ]
Paratanytarsus=subset(size3,Group.4=="Paratanytarsus")#8
Paratanytarsus=Paratanytarsus[order(-Paratanytarsus$temp), ]
Parachironomus=subset(size3,Group.4=="Parachironomus")#9
Parachironomus=Parachironomus[order(-Parachironomus$temp), ]
Bryophaenocladius=subset(size3,Group.4=="Bryophaenocladius")#10
Bryophaenocladius=Bryophaenocladius[order(-Bryophaenocladius$temp), ]
Rheotanytarsus=subset(size3,Group.4=="Rheotanytarsus")#11
Rheotanytarsus=Rheotanytarsus[order(-Rheotanytarsus$temp), ]
Cladotanytarsus=subset(size3,Group.4=="Cladotanytarsus")#12
Cladotanytarsus=Cladotanytarsus[order(-Cladotanytarsus$temp), ]
Conchapelopia=subset(size3,Group.4=="Conchapelopia")#13
Conchapelopia=Conchapelopia[order(-Conchapelopia$temp), ]
Chaetocladius=subset(size3,Group.4=="Chaetocladius")#14
Chaetocladius=Chaetocladius[order(-Chaetocladius$temp), ]
Stictochironomus=subset(size3,Group.4=="Stictochironomus")#15
Stictochironomus=Stictochironomus[order(-Stictochironomus$temp), ]
Pseudosmittia=subset(size3,Group.4=="Pseudosmittia")#16
Pseudosmittia=Pseudosmittia[order(-Pseudosmittia$temp), ]
Tanypus=subset(size3,Group.4=="Tanypus")#17
Tanypus=Tanypus[order(-Tanypus$temp), ]
Polypedilum=subset(size3,Group.4=="Polypedilum")#18
Polypedilum=Polypedilum[order(-Polypedilum$temp), ]
Zavrelimyia=subset(size3,Group.4=="Zavrelimyia")#19
Zavrelimyia=Zavrelimyia[order(-Zavrelimyia$temp), ]
Ablabesmyia=subset(size3,Group.4=="Ablabesmyia")#20
Ablabesmyia=Ablabesmyia[order(-Ablabesmyia$temp), ]
Corynoneura=subset(size3,Group.4=="Corynoneura")#21
Corynoneura=Corynoneura[order(-Corynoneura$temp), ]
Orthocladius=subset(size3,Group.4=="Orthocladius")#22
Orthocladius=Orthocladius[order(-Orthocladius$temp), ]
Eukiefferiella=subset(size3,Group.4=="Eukiefferiella")#23
Eukiefferiella=Eukiefferiella[order(-Eukiefferiella$temp), ]
Cladopelma=subset(size3,Group.4=="Cladopelma")#24
Cladopelma=Cladopelma[order(-Cladopelma$temp), ]


data=size3[,c(19,33)]
data1=na.omit(data)
#median based linear model, all genera vs temperature
mb2=mblm(wl~temp,data=data1)

mk_Chironomus=tidy(mk.test(Chironomus$wl))
tidy_Chironomus=tidy(mblm(wl~temp,data=Chironomus))

mk_Cricotopus=tidy(mk.test(Cricotopus$wl))
tidy_Cricotopus=tidy(mblm(wl~temp,data=Cricotopus))

mk_Dicrotendipes=tidy(mk.test(Dicrotendipes$wl))
tidy_Dicrotendipes=tidy(mblm(wl~temp,data=Dicrotendipes))

mk_Micropsectra=tidy(mk.test(Micropsectra$wl))
tidy_Micropsectra=tidy(mblm(wl~temp,data=Micropsectra))

mk_Tanytarsus=tidy(mk.test(Tanytarsus$wl))
tidy_Tanytarsus=tidy(mblm(wl~temp,data=Tanytarsus))

mk_Stenochironomus=tidy(mk.test(Stenochironomus$wl))
tidy_Stenochironomus=tidy(mblm(wl~temp,data=Stenochironomus))

mk_Paratanytarsus=tidy(mk.test(Paratanytarsus$wl))
tidy_Paratanytarsus=tidy(mblm(wl~temp,data=Paratanytarsus))

mk_Conchapelopia=tidy(mk.test(Conchapelopia$wl))
tidy_Conchapelopia=tidy(mblm(wl~temp,data=Conchapelopia))

mk_Bryophaenocladius=tidy(mk.test(Bryophaenocladius$wl))
tidy_Bryophaenocladius=tidy(mblm(wl~temp,data=Bryophaenocladius))

mk_Parachironomus=tidy(mk.test(Parachironomus$wl))
tidy_Parachironomus=tidy(mblm(wl~temp,data=Parachironomus))


mk_Cladotanytarsus=tidy(mk.test(Cladotanytarsus$wl))
tidy_Cladotanytarsus=tidy(mblm(wl~temp,data=Cladotanytarsus))

mk_Chaetocladius=tidy(mk.test(Chaetocladius$wl))
tidy_Chaetocladius=tidy(mblm(wl~temp,data=Chaetocladius))


mk_Limnophyes=tidy(mk.test(Limnophyes$wl))
tidy_Limnophyes=tidy(mblm(wl~temp,data=Limnophyes))

mk_Rheotanytarsus=tidy(mk.test(Rheotanytarsus$wl))
tidy_Rheotanytarsus=tidy(mblm(wl~temp,data=Rheotanytarsus))

mk_Ablabesmyia=tidy(mk.test(Ablabesmyia$wl))
tidy_Ablabesmyia=tidy(mblm(wl~temp,data=Ablabesmyia))


mk_Tanypus=tidy(mk.test(Tanypus$wl))
tidy_Tanypus=tidy(mblm(wl~temp,data=Tanypus))

mk_Stictochironomus=tidy(mk.test(Stictochironomus$wl))
tidy_Stictochironomus=tidy(mblm(wl~temp,data=Stictochironomus))


mk_Cladopelma=tidy(mk.test(Cladopelma$wl))
tidy_Cladopelma=tidy(mblm(wl~temp,data=Cladopelma))

mk_Corynoneura=tidy(mk.test(Corynoneura$wl))
tidy_Corynoneura=tidy(mblm(wl~temp,data=Corynoneura))

mk_Zavrelimyia=tidy(mk.test(Zavrelimyia$wl))
tidy_Zavrelimyia=tidy(mblm(wl~temp,data=Zavrelimyia))


mk_Eukiefferiella=tidy(mk.test(Eukiefferiella$wl))
tidy_Eukiefferiella=tidy(mblm(wl~temp,data=Eukiefferiella))


mk_Orthocladius=tidy(mk.test(Orthocladius$wl))
tidy_Orthocladius=tidy(mblm(wl~temp,data=Orthocladius))

mk_Pseudosmittia=tidy(mk.test(Pseudosmittia$wl))
tidy_Pseudosmittia=tidy(mblm(wl~temp,data=Pseudosmittia))
mk_Polypedilum=tidy(mk.test(Polypedilum$wl))
tidy_Polypedilum=tidy(mblm(wl~temp,data=Polypedilum))

#bind output of mblm models

tidyy=bind_rows(
  tidy_Cricotopus	,
  tidy_Dicrotendipes	,
  tidy_Micropsectra	,
  tidy_Tanytarsus	,
  tidy_Chironomus	,
  tidy_Stenochironomus	,
  tidy_Limnophyes	,
  tidy_Paratanytarsus	,
  tidy_Parachironomus	,
  tidy_Bryophaenocladius,
  tidy_Rheotanytarsus	,
  tidy_Cladotanytarsus	,
  tidy_Conchapelopia	,
  tidy_Chaetocladius	,
  tidy_Stictochironomus	,
  tidy_Pseudosmittia,
  tidy_Tanypus	,
  tidy_Polypedilum,
  tidy_Zavrelimyia	,
  tidy_Ablabesmyia	,
  tidy_Corynoneura	,
  tidy_Orthocladius,
  tidy_Eukiefferiella	,
  tidy_Cladopelma)





name=c( "Cricotopus",       "Dicrotendipes" ,    "Micropsectra" ,     "Tanytarsus"  ,
        "Chironomus" ,       "Stenochironomus" ,  "Limnophyes"   ,     "Paratanytarsus" ,
        "Parachironomus" ,   "Bryophaenocladius", "Rheotanytarsus",    "Cladotanytarsus",
        "Conchapelopia",     "Chaetocladius",     "Stictochironomus",  "Pseudosmittia",
        "Tanypus" ,          "Polypedilum"  ,     "Zavrelimyia",       "Ablabesmyia",
        "Corynoneura" ,      "Orthocladius" ,     "Eukiefferiella" ,   "Cladopelma"  )

tidyy=subset(tidyy,term=="temp" )
df1=cbind(name,tidyy)


#bind mk test 
tidmk=bind_rows(
  mk_Cricotopus	,
  mk_Dicrotendipes	,
  mk_Micropsectra	,
  mk_Tanytarsus	,
  mk_Chironomus	,
  mk_Stenochironomus	,
  mk_Limnophyes	,
  mk_Paratanytarsus	,
  mk_Parachironomus	,
  mk_Bryophaenocladius,
  mk_Rheotanytarsus	,
  mk_Cladotanytarsus	,
  mk_Conchapelopia	,
  mk_Chaetocladius	,
  mk_Stictochironomus	,
  mk_Pseudosmittia,
  mk_Tanypus	,
  mk_Polypedilum,
  mk_Zavrelimyia	,
  mk_Ablabesmyia	,
  mk_Corynoneura	,
  mk_Orthocladius,
  mk_Eukiefferiella	,
  mk_Cladopelma)
df2=cbind(name,tidmk)

#sort significant relationship wl vs temp

a <- df1 %>%
  filter( p.value>0.05) %>%
  mutate(final = ifelse(p.value >0.05, "insignificant"))
b <- df1 %>%
  filter( p.value<=0.05) %>%
  mutate(final = ifelse(p.value<=0.05, "significant"))
df3=rbind(a,b)

#plot relationships between wing length and temperature
g4=ggplot(df3, aes(as.factor(name),estimate))+ geom_point(size=5,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("[Âµm per Â°C]")
g4=g4+ theme(axis.text.y = element_text(face = "italic"))
g4=g4+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))

#########################################
#6) FINAL PLOTTING OF FIGURES
#Figure 1 - photo of the Chironomidae wing, not produced by code
# Figure 2, map with body size
sizep
#Figure 3
pex1_1
#Figure 4
grid.arrange(pex, pex_2,nrow=1) 
#Figure 5
grid.arrange( g3,g4,nrow=1)
#supplementary fig 1
tr=plot(obj,type="phylogram",leg.txt="wing length (mm)",lwd=2,
        mar=c(4,2,4,2))
#supplementary figure 2
grid.arrange( pex3,pex4,nrow=1) 



#ADDENDUM
################################
#code produces supplementary image showing lack of the
#signif. relationship between size and rate of wing length shift in log-normalaized dataset

library(ggpubr)
theme_set(theme_pubr())
size3=size2
#mblm function for all the genera above 60 records

##########################################
#######################################
##########################################
size3$wl=log10(size3$wl)+1
Cricotopus=subset(size3,Group.4=="Cricotopus")#1
Cricotopus=Cricotopus[order(-Cricotopus$Lat1_cor), ]#order by latitude (normalized) to allow mk test on the wl~lat
Dicrotendipes=subset(size3,Group.4=="Dicrotendipes")#2
Dicrotendipes=Dicrotendipes[order(-Dicrotendipes$Lat1_cor), ]
Micropsectra=subset(size3,Group.4=="Micropsectra")#3
Micropsectra=Micropsectra[order(-Micropsectra$Lat1_cor), ]
Tanytarsus=subset(size3,Group.4=="Tanytarsus")#4
Tanytarsus=Tanytarsus[order(-Tanytarsus$Lat1_cor), ]
Chironomus=subset(size3,Group.4=="Chironomus")#5
Chironomus=Chironomus[order(-Chironomus$Lat1_cor), ]
Stenochironomus=subset(size3,Group.4=="Stenochironomus")#6
Stenochironomus=Stenochironomus[order(-Stenochironomus$Lat1_cor), ]
Limnophyes=subset(size3,Group.4=="Limnophyes")#7
Limnophyes=Limnophyes[order(-Limnophyes$Lat1_cor), ]
Paratanytarsus=subset(size3,Group.4=="Paratanytarsus")#8
Paratanytarsus=Paratanytarsus[order(-Paratanytarsus$Lat1_cor), ]
Parachironomus=subset(size3,Group.4=="Parachironomus")#9
Parachironomus=Parachironomus[order(-Parachironomus$Lat1_cor), ]
Bryophaenocladius=subset(size3,Group.4=="Bryophaenocladius")#10
Bryophaenocladius=Bryophaenocladius[order(-Bryophaenocladius$Lat1_cor), ]
Rheotanytarsus=subset(size3,Group.4=="Rheotanytarsus")#11
Rheotanytarsus=Rheotanytarsus[order(-Rheotanytarsus$Lat1_cor), ]
Cladotanytarsus=subset(size3,Group.4=="Cladotanytarsus")#12
Cladotanytarsus=Cladotanytarsus[order(-Cladotanytarsus$Lat1_cor), ]
Conchapelopia=subset(size3,Group.4=="Conchapelopia")#13
Conchapelopia=Conchapelopia[order(-Conchapelopia$Lat1_cor), ]
Chaetocladius=subset(size3,Group.4=="Chaetocladius")#14
Chaetocladius=Chaetocladius[order(-Chaetocladius$Lat1_cor), ]
Stictochironomus=subset(size3,Group.4=="Stictochironomus")#15
Stictochironomus=Stictochironomus[order(-Stictochironomus$Lat1_cor), ]
Pseudosmittia=subset(size3,Group.4=="Pseudosmittia")#16
Pseudosmittia=Pseudosmittia[order(-Pseudosmittia$Lat1_cor), ]
Tanypus=subset(size3,Group.4=="Tanypus")#17
Tanypus=Tanypus[order(-Tanypus$Lat1_cor), ]
Polypedilum=subset(size3,Group.4=="Polypedilum")#18
Polypedilum=Polypedilum[order(-Polypedilum$Lat1_cor), ]
Zavrelimyia=subset(size3,Group.4=="Zavrelimyia")#19
Zavrelimyia=Zavrelimyia[order(-Zavrelimyia$Lat1_cor), ]
Ablabesmyia=subset(size3,Group.4=="Ablabesmyia")#20
Ablabesmyia=Ablabesmyia[order(-Ablabesmyia$Lat1_cor), ]
Corynoneura=subset(size3,Group.4=="Corynoneura")#21
Corynoneura=Corynoneura[order(-Corynoneura$Lat1_cor), ]
Orthocladius=subset(size3,Group.4=="Orthocladius")#22
Orthocladius=Orthocladius[order(-Orthocladius$Lat1_cor), ]
Eukiefferiella=subset(size3,Group.4=="Eukiefferiella")#23
Eukiefferiella=Eukiefferiella[order(-Eukiefferiella$Lat1_cor), ]
Cladopelma=subset(size3,Group.4=="Cladopelma")#24
Cladopelma=Cladopelma[order(-Cladopelma$Lat1_cor), ]

require(broom)
require(mblm)
require(trend)
#################################
###########################

mb1=mblm(wl~Lat1_cor,data=size3)
#Chironomus

mk_Chironomus=tidy(mk.test(Chironomus$wl))
tidy_Chironomus=tidy(mblm(wl~Lat1_cor,data=Chironomus))

mk_Cricotopus=tidy(mk.test(Cricotopus$wl))
tidy_Cricotopus=tidy(mblm(wl~Lat1_cor,data=Cricotopus))

mk_Dicrotendipes=tidy(mk.test(Dicrotendipes$wl))
tidy_Dicrotendipes=tidy(mblm(wl~Lat1_cor,data=Dicrotendipes))

mk_Micropsectra=tidy(mk.test(Micropsectra$wl))
tidy_Micropsectra=tidy(mblm(wl~Lat1_cor,data=Micropsectra))

mk_Tanytarsus=tidy(mk.test(Tanytarsus$wl))
tidy_Tanytarsus=tidy(mblm(wl~Lat1_cor,data=Tanytarsus))

mk_Stenochironomus=tidy(mk.test(Stenochironomus$wl))
tidy_Stenochironomus=tidy(mblm(wl~Lat1_cor,data=Stenochironomus))

mk_Paratanytarsus=tidy(mk.test(Paratanytarsus$wl))
tidy_Paratanytarsus=tidy(mblm(wl~Lat1_cor,data=Paratanytarsus))

mk_Conchapelopia=tidy(mk.test(Conchapelopia$wl))
tidy_Conchapelopia=tidy(mblm(wl~Lat1_cor,data=Conchapelopia))

mk_Bryophaenocladius=tidy(mk.test(Bryophaenocladius$wl))
tidy_Bryophaenocladius=tidy(mblm(wl~Lat1_cor,data=Bryophaenocladius))

mk_Parachironomus=tidy(mk.test(Parachironomus$wl))
tidy_Parachironomus=tidy(mblm(wl~Lat1_cor,data=Parachironomus))


mk_Cladotanytarsus=tidy(mk.test(Cladotanytarsus$wl))
tidy_Cladotanytarsus=tidy(mblm(wl~Lat1_cor,data=Cladotanytarsus))

mk_Chaetocladius=tidy(mk.test(Chaetocladius$wl))
tidy_Chaetocladius=tidy(mblm(wl~Lat1_cor,data=Chaetocladius))


mk_Limnophyes=tidy(mk.test(Limnophyes$wl))
tidy_Limnophyes=tidy(mblm(wl~Lat1_cor,data=Limnophyes))

mk_Rheotanytarsus=tidy(mk.test(Rheotanytarsus$wl))
tidy_Rheotanytarsus=tidy(mblm(wl~Lat1_cor,data=Rheotanytarsus))

mk_Ablabesmyia=tidy(mk.test(Ablabesmyia$wl))
tidy_Ablabesmyia=tidy(mblm(wl~Lat1_cor,data=Ablabesmyia))


mk_Tanypus=tidy(mk.test(Tanypus$wl))
tidy_Tanypus=tidy(mblm(wl~Lat1_cor,data=Tanypus))

mk_Stictochironomus=tidy(mk.test(Stictochironomus$wl))
tidy_Stictochironomus=tidy(mblm(wl~Lat1_cor,data=Stictochironomus))


mk_Cladopelma=tidy(mk.test(Cladopelma$wl))
tidy_Cladopelma=tidy(mblm(wl~Lat1_cor,data=Cladopelma))

mk_Corynoneura=tidy(mk.test(Corynoneura$wl))
tidy_Corynoneura=tidy(mblm(wl~Lat1_cor,data=Corynoneura))

mk_Zavrelimyia=tidy(mk.test(Zavrelimyia$wl))
tidy_Zavrelimyia=tidy(mblm(wl~Lat1_cor,data=Zavrelimyia))


mk_Eukiefferiella=tidy(mk.test(Eukiefferiella$wl))
tidy_Eukiefferiella=tidy(mblm(wl~Lat1_cor,data=Eukiefferiella))


mk_Orthocladius=tidy(mk.test(Orthocladius$wl))
tidy_Orthocladius=tidy(mblm(wl~Lat1_cor,data=Orthocladius))

mk_Pseudosmittia=tidy(mk.test(Pseudosmittia$wl))
tidy_Pseudosmittia=tidy(mblm(wl~Lat1_cor,data=Pseudosmittia))
mk_Polypedilum=tidy(mk.test(Polypedilum$wl))
tidy_Polypedilum=tidy(mblm(wl~Lat1_cor,data=Polypedilum))

#bind output of mblm models

tidyy=bind_rows(
  tidy_Cricotopus	,
  tidy_Dicrotendipes	,
  tidy_Micropsectra	,
  tidy_Tanytarsus	,
  tidy_Chironomus	,
  tidy_Stenochironomus	,
  tidy_Limnophyes	,
  tidy_Paratanytarsus	,
  tidy_Parachironomus	,
  tidy_Bryophaenocladius,
  tidy_Rheotanytarsus	,
  tidy_Cladotanytarsus	,
  tidy_Conchapelopia	,
  tidy_Chaetocladius	,
  tidy_Stictochironomus	,
  tidy_Pseudosmittia,
  tidy_Tanypus	,
  tidy_Polypedilum,
  tidy_Zavrelimyia	,
  tidy_Ablabesmyia	,
  tidy_Corynoneura	,
  tidy_Orthocladius,
  tidy_Eukiefferiella	,
  tidy_Cladopelma)





name=c( "Cricotopus",       "Dicrotendipes" ,    "Micropsectra" ,     "Tanytarsus"  ,
        "Chironomus" ,       "Stenochironomus" ,  "Limnophyes"   ,     "Paratanytarsus" ,
        "Parachironomus" ,   "Bryophaenocladius", "Rheotanytarsus",    "Cladotanytarsus",
        "Conchapelopia",     "Chaetocladius",     "Stictochironomus",  "Pseudosmittia",
        "Tanypus" ,          "Polypedilum"  ,     "Zavrelimyia",       "Ablabesmyia",
        "Corynoneura" ,      "Orthocladius" ,     "Eukiefferiella" ,   "Cladopelma"  )

tidyy=subset(tidyy,term=="Lat1_cor" )
df1=cbind(name,tidyy)


#bind mk test 
tidmk=bind_rows(
  mk_Cricotopus	,
  mk_Dicrotendipes	,
  mk_Micropsectra	,
  mk_Tanytarsus	,
  mk_Chironomus	,
  mk_Stenochironomus	,
  mk_Limnophyes	,
  mk_Paratanytarsus	,
  mk_Parachironomus	,
  mk_Bryophaenocladius,
  mk_Rheotanytarsus	,
  mk_Cladotanytarsus	,
  mk_Conchapelopia	,
  mk_Chaetocladius	,
  mk_Stictochironomus	,
  mk_Pseudosmittia,
  mk_Tanypus	,
  mk_Polypedilum,
  mk_Zavrelimyia	,
  mk_Ablabesmyia	,
  mk_Corynoneura	,
  mk_Orthocladius,
  mk_Eukiefferiella	,
  mk_Cladopelma)
df2=cbind(name,tidmk)

#sort significant relationship wl vs temp

a <- df1 %>%
  filter( p.value>0.05) %>%
  mutate(final = ifelse(p.value >0.05, "insignificant"))
b <- df1 %>%
  filter( p.value<=0.05) %>%
  mutate(final = ifelse(p.value<=0.05, "significant"))
df3=rbind(a,b)




size4=aggregate(size3[, 19], list(size3$Group.4),mean)#mean wing length per fossil amber deposit
#aggregate reduced dataset, average wing lenght per genus
size4$name=size4$Group.1
size_change<-merge(df3,size4,by = c("name"), all.x = TRUE, all.y = TRUE)
size_change$mean_wl=size_change$x
summary(lm(size_change$estimate~size_change$x))

#supplementary figure 3
plot(size_change$estimate~size_change$x,pch=19,xlab="mean wing size, log",ylab="size shift (log)")
abline(lm(size_change$estimate~size_change$x))

size4$Group.4=size4$Group.1
size5=size4




#rate of wing size change vs mean wing size
size_slope<-merge(size5,df3,by = c("name"), all.x = TRUE, all.y = TRUE)

mdy<- gls(estimate~x,data=size_slope, correlation=corAR1(), method="ML",na.action = na.omit)

mod.y1 <- update(mdy, correlation=corARMA(p=1))#best fit (still best fit)
anx=Anova(mod.y1)

size3 %>%
  group_by(spp) %>%tally()


