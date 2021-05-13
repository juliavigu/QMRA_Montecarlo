#############ANALYSIS OF MICROBIAL SAMPLES##########
#############JULIA VILA GUILERA- #################
###########UCL PHD PROJECT- 2021##########

#### PREPARING DATASET TO WORK WITH ####

#Import data
setwd("~/Desktop/OBJ2. full QMRA analysis")
dataset<-read.csv("codedsamplesr.csv")

#install needed packages
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("MASS")
library(MASS)
install.packages("viridis")
library(viridis)
install.packages("hrbrthemes")
library(hrbrthemes)

#Create Geomean function
geomean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#### QUICK LOOK AT RAW DATA ####
#Count number of samples collected, by type
#Village level samples
sum(with(dataset, HHID=="village" & dataset$Csamplesource==9)) #school and AWC soil
sum(with(dataset, HHID=="village" & (dataset$Csamplesource==1 ))) #school and AWC handpumps
sum(with(dataset,  dataset$Csamplesource==10 | dataset$Csamplesource==12)) #nalas and ponds

#Household-level samples
sum(with(dataset, HHID!="village" & dataset$Csamplesource==9 | dataset$Csamplesource==11)) #house soil or floor swab
sum(with(dataset, HHID!="village" & dataset$Csamplesource==1)) #house HP
sum(with(dataset,  dataset$Csamplesource==5)) #house wells
sum(with(dataset,  HHID!="village" & dataset$Csamplesource==6)) #house borewells
sum(with(dataset,  dataset$Csamplesource==2)) #house stored HP
sum(with(dataset,  dataset$Csamplesource==3)) #house stored well
sum(with(dataset,  dataset$Csamplesource==4)) #house stored borewell
sum(with(dataset,  dataset$Csamplesource==7)) #child's hands
sum(with(dataset,  dataset$Csamplesource==8)) #mothers hands


#### DATA CLEANING ####
#To perform analysis, drop errors, transform TMTC into max detection limit (300CFU/) transform non-detects into lowest detection limit (1CFU/100mL)
cleands<- subset(dataset, (Methoderror!=4 & Methoderror!=3)) #Method error 4= damaged sample, Method error 3=tmtc and too undiluted, we exclude them from analysis
cleands$Village[cleands$Village=="DevdaSath"]<-"Devdasath" #fixing some small issues with names
cleands$Block[cleands$Block=="Ghatol "]<-"Ghatol" #fixing some small issues with names

#### DESCRIPTIVE SUMMARY ANALYSIS OF MICROBIAL LOADS ####
## Faecal contamination loads and comparisons ##

### Number of samples by type (cleaned data, usable for analysis)
#Village level samples
sum(with(cleands, HHID=="village" & Csamplesource==9)) #school soil
sum(with(cleands, HHID=="village" & Csamplesource==1)) #school HP
sum(with(cleands,  Csamplesource==10)) #nalas

#Household-level samples
sum(with(cleands, HHID!="village" & Csamplesource==9)) #house soil
sum(with(cleands, HHID!="village" & Csamplesource==1)) #house HP
sum(with(cleands, Csamplesource==5)) #house wells
sum(with(cleands,  Csamplesource==6)) #house borewells
sum(with(cleands,  Csamplesource==2)) #house stored HP
sum(with(cleands,  Csamplesource==3)) #house stored well
sum(with(cleands,  Csamplesource==4)) #house stored borewell
sum(with(cleands,  Csamplesource==7)) #child's hands
sum(with(cleands,  Csamplesource==8)) #mothers hands

#Percentage of positive samples
cleands %>% group_by(Csamplesource) %>% summarise(sum(CFU_100mL==0))

#Now imputate censored data to be able to work with distributions
cleands$CFU_100mL[cleands$CFU_100mL==0]<-1 #substituting non-detects for minimum detection limit of 1cfu
cleands$Log10CFU<-log10(cleands$CFU_100mL) #recalculate logarithmic scale

#Logarithmic E.coli categories-frequency of samples
#dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU==0))
cleands %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU<=1))
cleands %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>1 & Log10CFU<=2))
cleands %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>2 & Log10CFU<=3))
cleands %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>3 & Log10CFU<=4))
cleands %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>4))

##Comparing HP contamination across villages
cleands %>%
  group_by(Village)%>%
  filter(Csamplesource==1)%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
         lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)

cleands %>%
  group_by(Block)%>%
  filter(Csamplesource==1)%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
       lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
       upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)

#HP significant differences between blocks?
tempdataset<-subset(cleands, (Block=="Ghatol" | Block=="Kushalgarh")& Csamplesource==1, na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)

##Comparing surface water contamination across villages and blocks
cleands %>%
  group_by(Village)%>%
  filter(Csamplesource==10)%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
         lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)

cleands %>%
  group_by(Block)%>%
  filter(Csamplesource==10)%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
         lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)

#Surface water significant differences between blocks?
tempdataset<-subset(cleands, (Block=="Ghatol" | Block=="Kushalgarh") & (Csamplesource==10), na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)

##Comparing school soil contamination across villages
dataset %>%
  group_by(Village)%>%
  filter(Csamplesource==9 & HHID=="village")%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
         lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)
cleands %>%
  group_by(Block)%>%
  filter(Csamplesource==9 & HHID=="village")%>%
  summarise(n=n(), 
            gm=geomean(Log10CFU), 
            sd=sd(Log10CFU))%>%
  mutate(se = sd / sqrt(n),
         lower.ci = gm - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = gm + qt(1 - (0.05 / 2), n - 1) * se)

#School soil significant differences between blocks?
tempdataset<-subset(cleands, (Block=="Ghatol" | Block=="Kushalgarh")& Csamplesource==9, na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)

##Comparing Sample types (t-test comparisons)
aggregate(cleands[, 15], list(cleands$Csamplesource), geomean, na.rm=TRUE)
aggregate(cleands[, 15], list(cleands$Csamplesource), sd, na.rm=TRUE)

#Between water source comparisons
#Handpump samples (N=68) vs Borewell samples (N=7)
tempdataset<-subset(cleands, Csamplesource==1 | Csamplesource==6, na.rm=TRUE)
table(tempdataset$Csamplesource)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump samples (N=68) vs open well samples (N=7)
tempdataset<-subset(cleands, Csamplesource==1 | Csamplesource==5, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell samples (N=7) vs well samples (N=7)
tempdataset<-subset(cleands, Csamplesource==5 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs well samples (N=7)
tempdataset<-subset(cleands, Csamplesource==5 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs handpump samples (N=68)
tempdataset<-subset(cleands, Csamplesource==1 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs borewell samples (N=7)
tempdataset<-subset(cleands, Csamplesource==6 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==3 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs stored borewelll samples (N=3)
tempdataset<-subset(cleands, Csamplesource==4 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=24) vs stored HP samples (N=37)
tempdataset<-subset(cleands, Csamplesource==2 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell sample (N=7) vs stored HP samples (N=37)
tempdataset<-subset(cleands, Csamplesource==2 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell sample (N=7) vs stored borewell samples (N=3)
tempdataset<-subset(cleands, Csamplesource==4 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell sample (N=7) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==3 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump sample (N=68) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==3 | Csamplesource==1, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump sample (N=68) vs stored borwell samples (N=3)
tempdataset<-subset(cleands, Csamplesource==4 | Csamplesource==1, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump sample (N=68) vs stored hanpdump samples (N=37)
tempdataset<-subset(cleands, Csamplesource==2 | Csamplesource==1, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Well sample (N=7) vs stored hanpdump samples (N=37)
tempdataset<-subset(cleands, Csamplesource==5 | Csamplesource==2, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Well sample (N=7) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==5 | Csamplesource==3, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Well sample (N=7) vs stored borewell samples (N=3)
tempdataset<-subset(cleands, Csamplesource==5 | Csamplesource==4, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Stored borewell (N=3) vs stored handmpum samples (N=37)
tempdataset<-subset(cleands, Csamplesource==2 | Csamplesource==4, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Stored borewell sample (N=3) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==3 | Csamplesource==4, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Stored handmpump sample (N=37) vs stored well samples (N=3)
tempdataset<-subset(cleands, Csamplesource==2 | Csamplesource==3, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)

###Comparing different sample types
#Create temp dataset where sample types are grouped into domestic water sources (1) vs stored water (2)
# vs hanswabs (7) vs streams (10) vs soil samples (9)
tempdataset<-cleands
tempdataset<-tempdataset%>% mutate(Csamplesource=recode(Csamplesource, '5'=1L,'6'=1L,'3'=2L,'4'=2L,'8'=7L))

#Geomean, SD and N for each sample type?
table(tempdataset$Csamplesource)
aggregate(tempdataset[, 15], list(tempdataset$Csamplesource), geomean, na.rm=TRUE)
aggregate(tempdataset[, 15], list(tempdataset$Csamplesource), sd, na.rm=TRUE)

#Significant differences between sample types?
#Source vs stored
compareds<-subset(tempdataset, Csamplesource==1 | Csamplesource==2, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Source vs handswab
compareds<-subset(tempdataset, Csamplesource==1 | Csamplesource==7, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Source vs stream
compareds<-subset(tempdataset, Csamplesource==1 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Source vs soil
compareds<-subset(tempdataset, Csamplesource==1 | Csamplesource==9, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Stored vs handswab
compareds<-subset(tempdataset, Csamplesource==2 | Csamplesource==7, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Stored vs stream
compareds<-subset(tempdataset, Csamplesource==2 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Stored vs soil
compareds<-subset(tempdataset, Csamplesource==2 | Csamplesource==9, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Handswab vs stream
compareds<-subset(tempdataset, Csamplesource==7 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Handswab vs soil
compareds<-subset(tempdataset, Csamplesource==7 | Csamplesource==9, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)
#Stream vs soil
compareds<-subset(tempdataset, Csamplesource==9 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=compareds, conf.int=TRUE)


#ECOLI CONTAMINATION DISTRIBUTIONS
#BOXPLOT WITH DIFERENT ENV. SAMPLES
#Create temp dataset where sample types are grouped into domestic water sources (1) vs stored water (2)
# vs hanswabs (7) vs streams (10) vs soil samples (9)
tempdataset<-cleands
tempdataset<-tempdataset %>% 
  mutate(Csamplesource=factor(Csamplesource))%>%
  mutate(Csamplesource=recode(Csamplesource, '5'="1",'6'="1",'3'="2",'4'="2",'8'="7"))%>%
  filter(Csamplesource=="1" | Csamplesource=="2" | Csamplesource=="7" | Csamplesource=="9" | Csamplesource=="10")

tempdataset %>%
  ggplot(aes (x=reorder(Csamplesource, Log10CFU, FUN=median), y=Log10CFU, color=Csamplesource)) +
  geom_boxplot() +
  scale_color_viridis(discrete = TRUE, begin=0.4, end=0.8, alpha=1, option="inferno") +
  geom_jitter(size=0.4, alpha=0.9) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("E. coli Log10 CFU/exposure unit")+
  xlab("")+
  scale_x_discrete(labels=c("Source water ","Stored  water","Hands","Streams", "Soil")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

#BOXPLOT WITH WATER SAMPLES
#Create temp dataset where sample types are grouped into 
# hp (1), shp(2), well (5), swell (3), borewell (6), sborewell (4), nala(10)
tempdataset<-cleands
tempdataset<-tempdataset %>% 
  mutate(Csamplesource=factor(Csamplesource))%>%
  filter(Csamplesource=="1" | Csamplesource=="2" | Csamplesource=="3" | 
           Csamplesource=="4" | Csamplesource=="5"| Csamplesource=="6" |
           Csamplesource=="10") %>%
  mutate(Csamplesource=recode(Csamplesource, '1'="A",'2'="B",'6'="C",'4'="D",'5'="E", '3'="F",'10'="G"))
#re-order in the desired order
tempdataset$Csamplesource<-factor(tempdataset$Csamplesource, levels=c("A", "B", "C", "D", "E", "F", "G"))

tempdataset %>%
  ggplot(aes (x=Csamplesource, y=Log10CFU, color=Csamplesource)) +
  geom_boxplot() +
  scale_color_viridis(discrete = TRUE, begin=0.4, end=0.8, alpha=1, option="inferno") +
  geom_jitter(size=0.4, alpha=0.9) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("E. coli Log10 CFU/100 mL")+
  xlab("")+
  scale_x_discrete(labels=c("Hand-pump","Stored hand-pump","Borewell","Stored borewell", "Well", "Stored well", "Streams")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))


###################################
#### RISK ASSESSMENT ANALYSIS ####
##################################

#Install some extra packages
install.packages("fitdistrplus")
library(fitdistrplus)
install.packages("mc2d")
library(mc2d)

###########   MONTECARLO MODEL ###########
#Prepare the clean data for the different pathways
hp<-subset(cleands, Csamplesource==1)
shp<-subset(cleands, Csamplesource==2)
borewell<-subset(cleands, Csamplesource==6)
sborewell<-subset(cleands, Csamplesource==4)
well<-subset(cleands, Csamplesource==5)
swell<-subset(cleands, Csamplesource==3)
nala<-subset(cleands, Csamplesource==10)
soil<-subset(cleands, Csamplesource==9)
childhand<-subset(cleands, Csamplesource==7)
caregiverhand<-subset(cleands, Csamplesource==8)
sourcewater<-subset(cleands, Csamplesource==1 | Csamplesource==5 | Csamplesource==6)
storedwater<-subset(cleands, Csamplesource==2 | Csamplesource==3 | Csamplesource==4)

####MODEL FOR ALL PATHWAYS

#E. COLI CONCENTRATIONS IN THE DIFFERENT ENVIRONMENTAL MEDIA
sourcewaterdistr<-fitdistr(sourcewater$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
storedwaterdistr<-fitdistr(storedwater$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
naladistr<-fitdistr(nala$CFU_100mL/100, "lognormal")# ln N(logmean, logsd) CFU/mL 
childhanddistr<-fitdistr(childhand$CFU_100mL, "lognormal") #ln N(logmean, logsd) CFU/hand
caregiverhanddistr<-fitdistr(caregiverhand$CFU_100mL, "lognormal") #ln N(logmean, logsd) CFU/hand
soildistr<-fitdistr(soil$CFU_100mL*1000, "lognormal") #ln N(logmean, logsd) CFU/mg

#SIMULATIONS ARE ASSIGNED RANDOM E.COLI CONCENTRATION AS PER THE DISTRIBUTIONS FITTED
ndvar(10000)
sourcewaterconc<-mcstoc(rlnorm,type="VU",meanlog=sourcewaterdistr$estimate['meanlog'],sdlog=sourcewaterdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
storedwaterconc<-mcstoc(rlnorm,type="VU",meanlog=storedwaterdistr$estimate['meanlog'],sdlog=storedwaterdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
nalaconc<-mcstoc(rlnorm,type="VU",meanlog=naladistr$estimate['meanlog'],sdlog=naladistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
childhandconc<-mcstoc(rlnorm,type="VU",meanlog=childhanddistr$estimate['meanlog'],sdlog=childhanddistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
caregiverhandconc<-mcstoc(rlnorm,type="VU",meanlog=caregiverhanddistr$estimate['meanlog'],sdlog=caregiverhanddistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
soilconc<-mcstoc(rlnorm,type="VU",meanlog=soildistr$estimate['meanlog'],sdlog=soildistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability

#DISTRIBUTION PARAMETERS FOR THE INGESTION OF DIFFERENT EXPOSURE PATHWAYS, PER AGE GROUP
#1-ingestion of drinking water (ml/day)
#0-5 months
meandrinkingwateringestion05<-62
sddrinkingwateringestion05<-161
#6-11 months
meandrinkingwateringestion611<-269
sddrinkingwateringestion611<-437
#12-23 months
meandrinkingwateringestion1223<-146
sddrinkingwateringestion1223<-255

#2-ingestion of bathing water (ml/day)
#0-5 months
meanbathingwateringestion05<-0
sdbathingwateringestion05<-0
#6-11 months
meanbathingwateringestion611<-0
sdbathingwateringestion611<-0
#12-23 months
meanbathingwateringestion1223<-38
sdbathingwateringestion1223<-35

#3-ingestion of soil (mg/day)
#0-5 months
meansoilingestion05<-19
sdsoilingestion05<-2
#6-11 months
meansoilingestion611<-90
sdsoilingestion611<-2
#12-23 months
meansoilingestion1223<-89
sdsoilingestion1223<-2

#4i5-frequency of hand-to-mouth contact (handfraction/day)
#0-5 months
shape05<-1.58
scale05<-58.22
#6-11 months
shape611<-1.58
scale611<-58.22
#12-23 months
shape1223<-2.66
scale1223<-56.62


#SIMULATIONS ARE ASSIGNED RANDOM INGESTION VOLUMES AS PER THE DISTRIBUTIONS FOUND
#1-Drinking water (stored water)
#0-5 months
storedwatering05<-mcstoc(rlnorm,type="VU",meanlog=log(meandrinkingwateringestion05),sdlog=log(sddrinkingwateringestion05), seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#6-11 months
storedwatering611<-mcstoc(rlnorm,type="VU",meanlog=log(meandrinkingwateringestion611),sdlog=log(sddrinkingwateringestion611), seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#12-23 months
storedwatering1223<-mcstoc(rlnorm,type="VU",meanlog=log(meandrinkingwateringestion1223),sdlog=log(sddrinkingwateringestion1223), seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#2-Bathing water 
#0-5 months
#bathingwatering05<-no ingestion
#6-11 months
#bathingwatering611<-noingestion
#12-23 months
bathingwatering1223<-mcstoc(rlnorm,type="VU",meanlog=log(meanbathingwateringestion1223),sdlog=log(sddrinkingwateringestion1223), seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#3-ingestion of soil
#0-5 months
soiling05<-mcstoc(rlnorm,type="VU",meanlog=log(meansoilingestion05),sdlog=log(sdsoilingestion05), seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)
#6-11 months
soiling611<-mcstoc(rlnorm,type="VU",meanlog=log(meansoilingestion611),sdlog=log(sdsoilingestion611), seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)
#12-23 months
soiling1223<-mcstoc(rlnorm,type="VU",meanlog=log(meansoilingestion1223),sdlog=log(sdsoilingestion1223), seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)

#4 i 5-frequency of hand contacts
#0-5 months
handcontacts05<-mcstoc(rweibull, type="VU", shape=shape05, scale=scale05, seed=41)
#6-11 months
handcontacts611<-mcstoc(rweibull, type="VU", shape=shape611, scale=scale611, seed=41)
#12-23 months
handcontacts1223<-mcstoc(rweibull, type="VU", shape=shape1223, scale=scale1223, seed=41)

#CALCULATE DAILY DOSES FOR EACH PATHWAY, PER AGE GROUPS
#Dose-response parameters that depend on the pathogens
#********input parameters***********#
ratio<-0.00005 #e.coli to pathogen ratio. E.coli0157=0.08, campylobacter=0.66, rotavirus=0.00005
alpha<-0.2531 #bpoisson pathogen parameter. E.coli0157=0.2099, campylobacter=0.145, rotavirus=0.2531
n50<-6.16 #bpoisson pathogen parameter. E.coli0157=1120, campylobacter=895.52, rotavirus=6.16

#Doses
#1-Drinking water (stored water)
#0-5 months
storedwaterdose05<-storedwaterconc*ratio*storedwatering05
#6-11 months
storedwaterdose611<-storedwaterconc*ratio*storedwatering611
#12-23 months
storedwaterdose1223<-storedwaterconc*ratio*storedwatering1223

#2-Bathing water 
#0-5 months
bathingwaterdose05<-0 #noingestion
#6-11 months
bathingwaterdose611<-0 #noingestion
#12-23 months
bathingwaterdose1223<-nalaconc*ratio*bathingwatering1223

#3-Soil
#0-5 months
soildose05<-soilconc*ratio*soiling05
#6-11 months
soildose611<-soilconc*ratio*soiling611
#12-23 months
soildose1223<-soilconc*ratio*soiling1223

#4-Child's hands
fractionHT<-0.067 #fraction of hand residue transferred per hand contact
fractionCH05<- (1-0.04) #fraction of contacts that are with child's hands
fractionCH611<- (1-0.46) #fraction of contacts that are with child's hands
fractionCH1223<- (1-0.23) #fraction of contacts that are with child's hands
meanhoursawake05<-10.3 #mean hours awake  for 0-5
meanhoursawake611<-11.3 #mean hours awake for 6-11
meanhoursawake1223<-11.4 #mean hours awake for 12-23
#0-5 months
childhanddose05<-childhandconc*ratio*handcontacts05*meanhoursawake05*fractionHT*fractionCH05
#6-11 months
childhanddose611<-childhandconc*ratio*handcontacts611*meanhoursawake611*fractionHT*fractionCH611
#12-23 months
childhanddose1223<-childhandconc*ratio*handcontacts1223*meanhoursawake1223*fractionHT*fractionCH1223

#5-Caregiver's hands
fractionHT<-0.067 #fraction of hand residue transferred per hand contact
fractionMH05<- (0.04) #fraction of contacts that are with mothers's hands
fractionMH611<- (0.46) #fraction of contacts that are with mothers's hands
fractionMH1223<- (0.23) #fraction of contacts that are with mothers's hands
meanhoursawake05<-10.3 #mean hours awake  for 0-5
meanhoursawake611<-11.3 #mean hours awake for 6-11
meanhoursawake1223<-11.4 #mean hours awake for 12-23
#0-5 months
caregiverhanddose05<-caregiverhandconc*ratio*handcontacts05*meanhoursawake05*fractionHT*fractionMH05
#6-11 months
caregiverhanddose611<-caregiverhandconc*ratio*handcontacts611*meanhoursawake611*fractionHT*fractionMH611
#12-23 months
caregiverhanddose1223<-caregiverhandconc*ratio*handcontacts1223*meanhoursawake1223*fractionHT*fractionMH1223

#CALCULATE DAILY INFECTION RISKS FOR EACH PATHWAY, PER AGE GROUP
#1-Drinking water (stored water)
#0-5 months
storedwaterrisk05<- 1-(1+(storedwaterdose05/n50)*(2^(1/alpha)-1))^-alpha
storedwaterinfrisk05<-mc(storedwaterrisk05)
#6-11 months
storedwaterrisk611<- 1-(1+(storedwaterdose611/n50)*(2^(1/alpha)-1))^-alpha
storedwaterinfrisk611<-mc(storedwaterrisk611)
#12-23 months
storedwaterrisk1223<- 1-(1+(storedwaterdose1223/n50)*(2^(1/alpha)-1))^-alpha
storedwaterinfrisk1223<-mc(storedwaterrisk1223)

#2- Bathing water 
#0-5 months
bathingwaterinfrisk05<-0
#6-11 months
bathingwaterinfrisk611<-0
#12-23 months
bathingwaterrisk1223<- 1-(1+(bathingwaterdose1223/n50)*(2^(1/alpha)-1))^-alpha
bathingwaterinfrisk1223<-mc(bathingwaterrisk1223)

#3- Soil
#0-5 months
soilrisk05<-1-(1+(soildose05/n50)*(2^(1/alpha)-1))^-alpha
soilinfrisk05<-mc(soilrisk05)
#6-11 months
soilrisk611<-1-(1+(soildose611/n50)*(2^(1/alpha)-1))^-alpha
soilinfrisk611<-mc(soilrisk611)
#12-23 months
soilrisk1223<- 1-(1+(soildose1223/n50)*(2^(1/alpha)-1))^-alpha
soilinfrisk1223<-mc(soilrisk1223)

#4- Child's hands
#0-5 months
childhandrisk05<-1-(1+(childhanddose05/n50)*(2^(1/alpha)-1))^-alpha
childhandinfrisk05<-mc(childhandrisk05)
#6-11 months
childhandrisk611<-1-(1+(childhanddose611/n50)*(2^(1/alpha)-1))^-alpha
childhandinfrisk611<-mc(childhandrisk611)
#12-23 months
childhandrisk1223<-1-(1+(childhanddose1223/n50)*(2^(1/alpha)-1))^-alpha
childhandinfrisk1223<-mc(childhandrisk1223)

#5- Caregiver's hands
#0-5 months
caregiverhandrisk05<-1-(1+(caregiverhanddose05/n50)*(2^(1/alpha)-1))^-alpha
caregiverhandinfrisk05<-mc(caregiverhandrisk05)
#6-11 months
caregiverhandrisk611<-1-(1+(caregiverhanddose611/n50)*(2^(1/alpha)-1))^-alpha
caregiverhandinfrisk611<-mc(caregiverhandrisk611)
#12-23 months
caregiverhandrisk1223<-1-(1+(caregiverhanddose1223/n50)*(2^(1/alpha)-1))^-alpha
caregiverhandinfrisk1223<-mc(caregiverhandrisk1223)

#OUTPUT RESULTS OF DAILY INFECTION RISKS
#1-Drinking water (stored water)
#0-5 months
print(storedwaterinfrisk05)
quantile(storedwaterinfrisk05, c(0.1,0.25,0.5,0.75, 0.9))
#6-11 months
print(storedwaterinfrisk611)
quantile(storedwaterinfrisk611, c(0.1,0.25,0.5,0.75, 0.9))
#12-23 months
print(storedwaterinfrisk1223)
quantile(storedwaterinfrisk1223, c(0.1,0.25,0.5,0.75, 0.9))


#2-Bathing water
#0-5 months
#no ingestion
#6-11 months
#no ingestion
#12-23 months
print(bathingwaterinfrisk1223)
quantile(bathingwaterinfrisk1223, c(0.1,0.25,0.5,0.75, 0.9))

#3-Soil
#0-5 months
print(soilinfrisk05)
quantile(soilinfrisk05, c(0.1,0.25,0.5,0.75, 0.9))
#6-11 months
print(soilinfrisk611)
quantile(soilinfrisk611, c(0.1,0.25,0.5,0.75, 0.9))
#12-23 months
print(soilinfrisk1223)
quantile(soilinfrisk1223, c(0.1,0.25,0.5,0.75, 0.9))

#4-Child's hands
#0-5 months
print(childhandinfrisk05)
quantile(childhandinfrisk05, c(0.1,0.25,0.5,0.75, 0.9))
#6-11 months
print(childhandinfrisk611)
quantile(childhandinfrisk611, c(0.1,0.25,0.5,0.75, 0.9))
#12-23 months
print(childhandinfrisk1223)
quantile(childhandinfrisk1223, c(0.1,0.25,0.5,0.75, 0.9))

#5-Caregiver's hands
#0-5 months
print(caregiverhandinfrisk05)
quantile(caregiverhandinfrisk05, c(0.1,0.25,0.5,0.75, 0.9))
#6-11 months
print(caregiverhandinfrisk611)
quantile(caregiverhandinfrisk611, c(0.1,0.25,0.5,0.75, 0.9))
#12-23 months
print(caregiverhandinfrisk1223)
quantile(caregiverhandinfrisk1223, c(0.1,0.25,0.5,0.75, 0.9))


#Visualise results with boxplot
boxplot(storedwaterinfrisk05, at=1, xlim=c(0,13), ylim=c(0,1), ylab="Daily infection risk (%)", 
       col="orangered", outline=FALSE)
boxplot(storedwaterinfrisk611, at=2, add=TRUE,  col="sienna2")
boxplot(storedwaterinfrisk1223, at=3, add=TRUE, outcex=0.1, col="tan1")
boxplot(bathingwaterinfrisk1223, at=4, add=TRUE,  col="tan1")
boxplot(soilinfrisk05, at=5, add=TRUE, outline=FALSE,  col="orangered")
boxplot(soilinfrisk611, at=6, add=TRUE, outline=FALSE, col="sienna2")
boxplot(soilinfrisk1223, at=7, add=TRUE, outline=FALSE, col="tan1")
boxplot(childhandinfrisk05, at=8, add=TRUE, outline=FALSE, col="orangered")
boxplot(childhandinfrisk611, at=9, add=TRUE, outline=FALSE, col="sienna2")
boxplot(childhandinfrisk1223, at=10, add=TRUE, outline=FALSE, col="tan1")
boxplot(caregiverhandinfrisk05, at=11, add=TRUE, outline=FALSE, col="orangered")
boxplot(caregiverhandinfrisk611, at=12, add=TRUE, outline=FALSE, col="sienna2")
boxplot(caregiverhandinfrisk1223, at=13, add=TRUE, outline=FALSE, col="tan1")
axis(side = 1, at = c(2,4,6,9,12), cex.axis=0.5, las=2, labels = c("Drinking water","Bathing water", "Soil ingestion", "Mouthing own hands", "Mouthing caregiver's hands"))
#legend("topright", fill = c("sienna1", "sienna2", "sienna3", legend = c("0-5 months","6-11 months","12-23 months"), horiz = F)



#MODEL FOR A SINGLE PATHWAY AND EXPOSURE
#Read in the microbial ***INPUT DATA***
inputdata<-sourcewater$CFU_100mL
unittransform<-0.01 #to change from cfu/100ml to cfu/ml. Only for water samples
ratio<-0.66 #e.coli to pathogen ratio. E.coli0157=0.08, campylobacter=0.66, rotavirus=0.00005
alpha<-0.145 #bpoisson pathogen parameter. E.coli0157=0.2099, campylobacter=0.145, rotavirus=0.2531
n50<-1120 #bpoisson pathogen parameter. E.coli0157=1120, campylobacter=895.52, rotavirus=6.16
n<-365 #number of yearly exposures
#Log-normally distributed exposures
meaningestion<-146 #mean mL, grams or hand fraction ingested
sdingestion<-216 #sd
#Weibull
shape<-1.58
scale<-58
#handsfreq<-rweibull(n=10000, shape=shape, scale=scale)
#childhands<-handsfreq

#***INPUT DATA***

#Convert e.coli to pathogen concentrations 
conc<-inputdata

#Fit distributions to microbial concentration data
fit<-fitdistr(conc, "lognormal")

#Generate mcnodes with pathogen concentration and ingestion parameters
ndvar(10000)

c<-mcstoc(rlnorm,type="VU",meanlog=fit$estimate['meanlog'],sdlog=fit$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
i<-mcstoc(rlnorm,type="VU",meanlog=log(meaningestion),sdlog=log(sdingestion), seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#Calculate daily dose of exposure
dose<-(c*ratio*unittransform)*(i)

#Calculate daily infection risk following the beta-poisson function
risk<- 1-(1+(dose/n50)*(2^(1/alpha)-1))^-alpha
infrisk<-mc(risk)
sourcewaterinfrisk<-infrisk
#Calculate annual infection risk
riskannual<-(1-(1-risk)^n)
infriskannual<-mc(riskannual)

#Visualise results
print(infrisk)
print(infriskannual)
plot(density(infrisk$risk))
plot(density(infriskannual$riskannual))

boxplot(storedwaterinfrisk, at=1, xlim=c(0, 3), ylab="Daily Campylobacter infection risk", 
        main="Daily infection risk", col=(c("darkorange")))
boxplot(sourcewaterinfrisk, at=2, add=TRUE, col=(c("darkorange4")))
axis(side = 1, at = c(1,2), labels = c("storedwater","sourcewater"))
legend("topright", fill = rainbow(3, s = 0.4), legend = c("0-5 months","6-11 months","12-23 months"), horiz = F)

plot.mc <- function(infrisk, prec=0.001, stat = c("median","mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm=TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint=TRUE, xlim=NULL,ylim=NULL)
  

##############################
# Define Beta poisson model
bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}
expon<-function(k,d){
  return((1-(exp(1)^(-k*d))))
}

#Pathogen dose response parameters
pathogen<-boxplot(
  org=c("ecoli0157","campylobacter","rotavirus"),
  alpha=c(0.2099, 0.1450, 0.2531),
  n50=c(1120.00,895.52, 6.16),
  pdi=c(NA,0.3, 0.5), #only used for calculation of disease risk
  k=c(0.001715176,NA, NA), #only used in the exponential model, not for betapoisson
  ratio=c(0.08, 0.66, 0.00005),
  equation=c("bpois", "bpois", "bpois")
)
#Simulate infection risk for each pathogen
#Run 10000 montecarlo simulations for each pathogen to estimate infection risk given the 
#parameter distributions for c i t and d

simulator<-function(patho){ #for each pathogen
  ratio<-as.numeric(patho["ratio"]) #store the ratio info for each pathogen
  conc<-convertdata(ratio) #pathogen concentrations in samples are stored in "conc"
  #mcstoc is a function that evaluates 3 ratios for each mcnode (i.e. each mc object). These ratios are
  #the variability ratio the uncertainty ratio and the overall uncertinty ratio.
  #it is assumed that the concentration of pathogens has variability
  c<-mcstoc(rlnorm,type="V",meanlog=mean(log(conc)),sdlog=sd(log(conc))) #distribution of pathogen concentration with "V"ariability
  #it is assumed that exposure has variability and uncertainty
  #children ingest between 200 and 900 ml of water every day
  i<-mcstoc(rlnorm,type="VU",meanlog=log(3),sdlog=log(3)) #lognorm distr of ingestion of pathogens, with a min and max, with Variability and Uncertainty (VU)
  #time and duration of exposure 1 day
  t<-1
  
  #Dose (ingestion in ml/hour)  C, I , T represent distributions and a montecarlo simulator 
  #will be used to sample from those distributions to output a distribution for the dose 
  d<-c*i*t
  
  #MC
  ndvar(10000)
  
  riski<-if(patho["equation"]=="expon"){
  k<-as.numeric(patho["k"])
  expon(k,d)
  } else if(pathogen["equation"]=="bpois"){
    n50<-as.numeric(patho["n50"])
    alpha<-as.numeric(patho["alpha"])
    bpois(n50,alpha,d)
    }
  
  #riskd <- riski *as.numeric(path["pdi"])
  
  #risk infection per event
  #if (path["org"] != "ecoli"){
  infrisk<-mc(c,d,riski)
  print(patho["org"])
  #print(infrisk)
  summary(infrisk)
  #plot(infrisk)
  
  #plot kernel density of risk
  plot(density(infrisk$riski))
  #hist(infrisk)
  #}
  
  n<-365
  infriskperyear<-(1-(1-riski)^n)
  annualinfrisk<-mc(c,d,infriskperyear)
  summary(annualinfrisk)
 
  #risk disease per event 
  #disrisk<-mc(c,d,riskd)
  #print(path["org"])
  #print(diskrisk)
  #summary(diskrisk)
  #plot(diskrisk)
}

apply(pathogen,1,simulator)

## Exploratory data analysis of distributions
#Distribution of pathogen loads
hp<-subset(cleands, cleands$Csamplesource==9)
hist(hp$Log10CFU, breaks=30)

hist(cleands$CFU_100mL)
cleands$Log10CFU<-log10(cleands$CFU_100mL)
hist(cleands$Log10CFU, breaks=30)


#Subset for stored handpump samples only (n=37)
shp<-dataset[dataset$Csamplesource == 2,c("CFU_100mL","Log10CFU")]
#only for positive samples...
shp<-shp[!(apply(shp,1,function(y) any(y==0))),]
#fem transformacio ln de la concentracio de ecolis
shp$LnCFU<-log(shp$CFU_100mL)
#investiguem la distribucio de les shp mostres en forma ln
head(shp$LnCFU)
x <- seq(0,max(shp$LnCFU),length=700) #creem vector q determina quan hist breaks
hst<-hist(shp$LnCFU, breaks=10)
#fitegem la lognorm distribution a les shp i ens guardem els valors del fit 
fit_params <- fitdistr(shp$LnCFU,"lognormal")
fit <- dlnorm(x, fit_params$estimate['meanlog'], fit_params$estimate['sdlog'])

#comparem el lognorm fit aamb les dades
plot(x, fit, type="l", ylab="Density",
     xlab="X", ylim=c(0,max(hst$density)), xlim=c(0,10))
title(main = "Density histogram with lognormal fit")
lines(hst$mid, hst$density, type="l", col="red")
legend(8,0.15,legend=c("Fit","Data"),lty=c(1,1),col=c("black","red"))

ggplot(data=shp) + 
  geom_histogram(binwidth = 0.8, aes(LnCFU)) + 
  stat_function(fun = dlnorm, args = list(meanlog = 1.619, sdlog = 0.0448), 
                colour = "red")

wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)

##############DATASET WITH (TMTC & TFTC)##############
#Number of samples
for (i in 1:12) {print(sum(dataset$Csamplesource==i))}

#GeoMean, max, min, and Nsamples of all sample sources
aggregate(dataset[, 15], list(dataset$Csamplesource), summary, na.rm=TRUE)
aggregate(dataset[, 15], list(dataset$Csamplesource), geomean, na.rm=TRUE)
aggregate(dataset[, 15], list(dataset$Csamplesource), mean, na.rm=TRUE)
aggregate(dataset[, 15], list(dataset$Csamplesource), sd, na.rm=TRUE)

#Comparison source and stored water
source<-subset(dataset, (dataset$Csamplesource==1 |dataset$Csamplesource==5| dataset$Csamplesource==6))
stored<-subset(dataset, (dataset$Csamplesource==2 |dataset$Csamplesource==3| dataset$Csamplesource==4))
mean(source$Log10CFU)
mean(stored$Log10CFU)
geomean(stored$Log10CFU)
geomean(source$Log10CFU)
geomean(stored$Log10CFU)/geomean(source$Log10CFU)
mean(stored$Log10CFU)/mean(source$Log10CFU) #difference between stored and source water contamination level
sd(source$Log10CFU)
sd(stored$Log10CFU)
max(stored$Log10CFU)
min(stored$Log10CFU)

#Log categories
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU==0))
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>0 & Log10CFU<=1))
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>1 & Log10CFU<=2))
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>2 & Log10CFU<=3))
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>3 & Log10CFU<=4))
dataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>4))

#Log categories for overall stored water
stored %>% summarise(sum(Log10CFU==0))
stored %>% summarise(sum(Log10CFU>0 & Log10CFU<=1))
stored %>% summarise(sum(Log10CFU>1 & Log10CFU<=2))
stored %>% summarise(sum(Log10CFU>2 & Log10CFU<=3))
stored %>% summarise(sum(Log10CFU>3 & Log10CFU<=4))
stored %>% summarise(sum(Log10CFU>4))

#############DATASET WITHout ERRORS (TMTC & TFTC)- FOR SOIL##############
#Datset without errors
Gooddataset<- subset(dataset, (dataset$Methoderror!=4 & dataset$Methoderror!=1 & dataset$Methoderror!=3))
housegooddat<-subset(Gooddataset, Gooddataset$HHID!="village")
housedat<-subset(dataset, dataset$HHID!="village")

#Number of samples without errors
for (i in 1:12) {print(sum(Gooddataset$Csamplesource==i))}
#Number of household samples without errors
for (i in 1:12) {print(sum(housegooddat$Csamplesource==i))}

#Mean, max, min, and Nsamples of all sample sources without errors
aggregate(Gooddataset[, 15], list(Gooddataset$Csamplesource), geomean, na.rm=TRUE)
aggregate(Gooddataset[, 15], list(Gooddataset$Csamplesource), sd, na.rm=TRUE)
aggregate(Gooddataset[, 15], list(Gooddataset$Csamplesource), summary, na.rm=TRUE)

#Mean, max, min, and Nsamples of household sample sources without errors
aggregate(housegooddat[, 15], list(housegooddat$Csamplesource), mean, na.rm=TRUE)
aggregate(housegooddat[, 15], list(housegooddat$Csamplesource), geomean, na.rm=TRUE)
aggregate(housegooddat[, 15], list(housegooddat$Csamplesource), sd, na.rm=TRUE)
aggregate(housegooddat[, 15], list(housegooddat$Csamplesource), summary, na.rm=TRUE)

#Log categories for household samples without errors
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU==0))
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>0 & Log10CFU<=1))
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>1 & Log10CFU<=2))
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>2 & Log10CFU<=3))
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>3 & Log10CFU<=4))
housegooddat %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>4))

##########DATASET WITHOUT FAULTY HAND SAMPLES ONLY##############
#Datset without errors
nofaultdataset<- subset(dataset, (dataset$Methoderror!=4))

#Number of samples without faulty samples
for (i in 7:8) {print(sum(nofaultdataset$Csamplesource==i))}

#Mean, max, min, stdev and Nsamples of all sample sources
aggregate(nofaultdataset[, 15], list(nofaultdataset$Csamplesource), mean, na.rm=TRUE)
aggregate(nofaultdataset[, 15], list(nofaultdataset$Csamplesource), geomean, na.rm=TRUE)
aggregate(nofaultdataset[, 15], list(nofaultdataset$Csamplesource), sd, na.rm=TRUE)
aggregate(nofaultdataset[, 15], list(nofaultdataset$Csamplesource), summary, na.rm=TRUE)

#Log categories
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU==0))
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>0 & Log10CFU<=1))
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>1 & Log10CFU<=2))
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>2 & Log10CFU<=3))
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>3 & Log10CFU<=4))
nofaultdataset %>% group_by(Csamplesource) %>% summarise(sum(Log10CFU>4))



###*******COMPARING LOCATIONS*****

tempdataset<-subset(dataset, Csamplesource==3 | dataset$Csamplesource==4, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
tempdataset<-subset(dataset, Csamplesource==1, na.rm=TRUE)
wilcox.test(tempdataset$Log10CFU, data=tempdataset, conf.int=TRUE)

#Comparing Ghatol and Kushalgarh (t-test comparisons) 

#Handpump samples (N=68)
tempdataset<-subset(dataset, (Block=="Ghatol" | Block=="Kushalgarh")& Csamplesource==1, na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)
#Nala samples (N=25)
tempdataset<-subset(dataset, (Block=="Ghatol" | Block=="Kushalgarh")& Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)
#Stored water from handpump samples (N=37)
tempdataset<-subset(dataset, (Block=="Ghatol" | Block=="Kushalgarh")& Csamplesource==2, na.rm=TRUE)
wilcox.test(Log10CFU~Block, data=tempdataset, conf.int=TRUE)

#Comparing Sample types (t-test comparisons) 
#between Source comparisons
#Handpump samples (N=68) vs Borewell samples (N=7)
tempdataset<-subset(dataset, Csamplesource==1 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump samples (N=68) vs open well samples (N=7)
tempdataset<-subset(dataset, Csamplesource==1 | Csamplesource==5, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell samples (N=7) vs well samples (N=7)
tempdataset<-subset(dataset, Csamplesource==5 | Csamplesource==6, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=25) vs well samples (N=7)
tempdataset<-subset(dataset, Csamplesource==5 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Nala samples (N=25) vs handpump samples (N=68)
tempdataset<-subset(dataset, Csamplesource==1 | Csamplesource==10, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)

#between stored compariosns
#Handpump stored samples (N=37) vs Borewell stored samples (N=3)
tempdataset<-subset(dataset, Csamplesource==2 | Csamplesource==4, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Handpump stored samples (N=37) vs open well stored samples (N=3)
tempdataset<-subset(dataset, Csamplesource==2 | Csamplesource==3, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Borewell stored samples (N=3) vs well stored samples (N=3)
tempdataset<-subset(dataset, Csamplesource==4 | Csamplesource==3, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)

#Stored vs source comparison
#handpump samples (N=68) vs stored handpump samples (N=37)
tempdataset<-subset(dataset, Csamplesource==1 | Csamplesource==2, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#Well samples (N=7) vs stored well samples (N=3)
tempdataset<-subset(dataset, Csamplesource==3 | Csamplesource==5, na.rm=TRUE)
wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)
#boreWell samples (N=7) vs stored borewell samples (N=3)
tempdataset<-subset(dataset, Csamplesource==4 | Csamplesource==6, na.rm=TRUE)


####### Probabilistic model Montecarlo
install.packages("ggplot2")
install.packages("MASS")
library(ggplot2)
library(MASS)
#Subset for stored handpump samples only (n=37)
shp<-dataset[dataset$Csamplesource == 2,c("CFU_100mL","Log10CFU")]
#only for positive samples...
shp<-shp[!(apply(shp,1,function(y) any(y==0))),]
#fem transformacio ln de la concentracio de ecolis
shp$LnCFU<-log(shp$CFU_100mL)
#investiguem la distribucio de les shp mostres en forma ln
head(shp$LnCFU)
x <- seq(0,max(shp$LnCFU),length=700) #creem vector q determina quan hist breaks
hst<-hist(shp$LnCFU, breaks=10)
#fitegem la lognorm distribution a les shp i ens guardem els valors del fit 
fit_params <- fitdistr(shp$LnCFU,"lognormal")
fit <- dlnorm(x, fit_params$estimate['meanlog'], fit_params$estimate['sdlog'])

#comparem el lognorm fit aamb les dades
plot(x, fit, type="l", ylab="Density",
     xlab="X", ylim=c(0,max(hst$density)), xlim=c(0,10))
title(main = "Density histogram with lognormal fit")
lines(hst$mid, hst$density, type="l", col="red")
legend(8,0.15,legend=c("Fit","Data"),lty=c(1,1),col=c("black","red"))

ggplot(data=shp) + 
  geom_histogram(binwidth = 0.8, aes(LnCFU)) + 
  stat_function(fun = dlnorm, args = list(meanlog = 1.619, sdlog = 0.0448), 
                colour = "red")

wilcox.test(Log10CFU~Csamplesource, data=tempdataset, conf.int=TRUE)