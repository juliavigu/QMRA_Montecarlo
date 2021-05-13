#Install packages
install.packages("fitdistrplus")
library(fitdistrplus)
install.packages("mc2d")
library(mc2d)
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

#Import data
setwd("~/Desktop/OBJ2. full QMRA analysis")
dataset<-read.csv("codedsamplesr.csv")

#Clean data
cleands<- subset(dataset, (Methoderror!=4 & Methoderror!=3)) #Method error 4= damaged sample, Method error 3=tmtc and too undiluted, we exclude them from analysis
cleands$Village[cleands$Village=="DevdaSath"]<-"Devdasath" #fixing some small issues with names
cleands$Block[cleands$Block=="Ghatol "]<-"Ghatol" #fixing some small issues with names
#Now imputate censored data
cleands$CFU_100mL[cleands$CFU_100mL==0]<-1 #substituting non-detects for minimum detection limit of 1cfu
cleands$Log10CFU<-log10(cleands$CFU_100mL) #recalculate logarithmic scale


###########   MONTECARLO MODEL ###########
#Prepare the clean data for the different pathways
hp<-subset(cleands, Csamplesource==1)
shp<-subset(cleands, Csamplesource==2)
bw<-subset(cleands, Csamplesource==6)
sbw<-subset(cleands, Csamplesource==4)
w<-subset(cleands, Csamplesource==5)
sw<-subset(cleands, Csamplesource==3)
nala<-subset(cleands, Csamplesource==10)

####MODEL FOR ALL PATHWAYS
ndvar(10000) #number of samples for the variability dimension
nsu(100) #number of samples for the uncertainty dimension

#E. COLI CONCENTRATIONS IN THE DIFFERENT ENVIRONMENTAL MEDIA
#Define e. coli concentration distributions
hpdistr<-fitdistr(hp$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
shpdistr<-fitdistr(shp$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
bwdistr<-fitdistr(bw$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
sbwdistr<-fitdistr(sbw$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
wdistr<-fitdistr(w$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
swdistr<-fitdistr(sw$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
naladistr<-fitdistr(nala$CFU_100mL/100, "lognormal")# ln N(logmean, logsd) CFU/mL 

#SIMULATIONS ARE ASSIGNED RANDOM E.COLI CONCENTRATION AS PER THE DISTRIBUTIONS FITTED
#create a matrix of 10000*100 simulations, to which random values from the e.coli distributions are assigned
hpconc<-mcstoc(rlnorm,type="VU",meanlog=hpdistr$estimate['meanlog'],sdlog=hpdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
shpconc<-mcstoc(rlnorm,type="VU",meanlog=shpdistr$estimate['meanlog'],sdlog=shpdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
bwconc<-mcstoc(rlnorm,type="VU",meanlog=bwdistr$estimate['meanlog'],sdlog=bwdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
sbwconc<-mcstoc(rlnorm,type="VU",meanlog=sbwdistr$estimate['meanlog'],sdlog=sbwdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
wconc<-mcstoc(rlnorm,type="VU",meanlog=wdistr$estimate['meanlog'],sdlog=wdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
swconc<-mcstoc(rlnorm,type="VU",meanlog=swdistr$estimate['meanlog'],sdlog=swdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
nalaconc<-mcstoc(rlnorm,type="VU",meanlog=naladistr$estimate['meanlog'],sdlog=naladistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability

#DISTRIBUTION PARAMETERS FOR THE INGESTION OF DIFFERENT EXPOSURE PATHWAYS, PER AGE GROUP
#1-ingestion of drinking water (ml/day)
#0-5 months
logmeandrinkingwateringestion05<-4.12
logsddrinkingwateringestion05<-1.23
#6-11 months
logmeandrinkingwateringestion611<-5.59
logsddrinkingwateringestion611<-0.85
#12-23 months
logmeandrinkingwateringestion1223<-4.98
logsddrinkingwateringestion1223<-0.89

#2-ingestion of bathing water (ml/day)
#0-5 months
#logmeanbathingwateringestion05<-0
#logsdbathingwateringestion05<-0
#6-11 months
#logmeanbathingwateringestion611<-0
#logsdbathingwateringestion611<-0
#12-23 months
logmeanbathingwateringestion1223<-3.64
logsdbathingwateringestion1223<-0.72


#SIMULATIONS ARE ASSIGNED RANDOM INGESTION VOLUMES AS PER THE DISTRIBUTIONS FOUND
#1-Drinking water 
#0-5 months
drinkingwatering05<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion05,sdlog=logsddrinkingwateringestion05, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#6-11 months
drinkingwatering611<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion611,sdlog=logsddrinkingwateringestion611, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#12-23 months
drinkingwatering1223<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion1223,sdlog=logsddrinkingwateringestion1223, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#2-Bathing water 
#0-5 months
#bathingwatering05<-no ingestion
#6-11 months
#bathingwatering611<-noingestion
#12-23 months
bathingwatering1223<-mcstoc(rlnorm,type="VU",meanlog=logmeanbathingwateringestion1223,sdlog=logsddrinkingwateringestion1223, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#CALCULATE DAILY DOSES FOR EACH PATHWAY, PER AGE GROUPS
#Dose-response parameters that depend on the pathogens
#********input parameters***********#
#************************************
#************************************
ratio<-0.08 #e.coli to pathogen ratio. E.coli0157=0.08, campylobacter=0.66, rotavirus=0.00005
alpha<-0.145 #bpoisson pathogen parameter. E.coli0157=0.2099, campylobacter=0.145, rotavirus=0.2531
n50<-1120 #bpoisson pathogen parameter. E.coli0157=1120, campylobacter=895.52, rotavirus=6.16
#*****************************
#*****************************
#*****************************
#Doses
#1-Drinking water (hp)
#0-5 months
hpdose05<-hpconc*ratio*drinkingwatering05
#6-11 months
hpdose611<-hpconc*ratio*drinkingwatering611
#12-23 months
hpdose1223<-hpconc*ratio*drinkingwatering1223

#2-Drinking water (shp)
#0-5 months
shpdose05<-shpconc*ratio*drinkingwatering05
#6-11 months
shpdose611<-shpconc*ratio*drinkingwatering611
#12-23 months
shpdose1223<-shpconc*ratio*drinkingwatering1223

#3-Drinking water (bw)
#0-5 months
bwdose05<-bwconc*ratio*drinkingwatering05
#6-11 months
bwdose611<-bwconc*ratio*drinkingwatering611
#12-23 months
bwdose1223<-bwconc*ratio*drinkingwatering1223

#4-Drinking water (sbw)
#0-5 months
sbwdose05<-sbwconc*ratio*drinkingwatering05
#6-11 months
sbwdose611<-sbwconc*ratio*drinkingwatering611
#12-23 months
sbwdose1223<-sbwconc*ratio*drinkingwatering1223

#5-Drinking water (w)
#0-5 months
wdose05<-wconc*ratio*drinkingwatering05
#6-11 months
wdose611<-wconc*ratio*drinkingwatering611
#12-23 months
wdose1223<-wconc*ratio*drinkingwatering1223

#6-Drinking water (sw)
#0-5 months
swdose05<-swconc*ratio*drinkingwatering05
#6-11 months
swdose611<-swconc*ratio*drinkingwatering611
#12-23 months
swdose1223<-swconc*ratio*drinkingwatering1223

#7-Bathing water (nala)
#0-5 months
naladose05<-0 #noingestion
#6-11 months
naladose611<-0 #noingestion
#12-23 months
naladose1223<-nalaconc*ratio*bathingwatering1223


#CALCULATE DAILY INFECTION RISKS FOR EACH PATHWAY, PER AGE GROUP
#1-Drinking water (hp)
#0-5 months
hprisk05<- 1-(1+(hpdose05/n50)*(2^(1/alpha)-1))^-alpha
hpinfrisk05<-mc(hprisk05)
#6-11 months
hprisk611<- 1-(1+(hpdose611/n50)*(2^(1/alpha)-1))^-alpha
hpinfrisk611<-mc(hprisk611)
#12-23 months
hprisk1223<- 1-(1+(hpdose1223/n50)*(2^(1/alpha)-1))^-alpha
hpinfrisk1223<-mc(hprisk1223)

#2-Drinking water (shp)
#0-5 months
shprisk05<- 1-(1+(shpdose05/n50)*(2^(1/alpha)-1))^-alpha
shpinfrisk05<-mc(shprisk05)
#6-11 months
shprisk611<- 1-(1+(shpdose611/n50)*(2^(1/alpha)-1))^-alpha
shpinfrisk611<-mc(shprisk611)
#12-23 months
shprisk1223<- 1-(1+(shpdose1223/n50)*(2^(1/alpha)-1))^-alpha
shpinfrisk1223<-mc(shprisk1223)

#3-Drinking water (bw)
#0-5 months
bwrisk05<- 1-(1+(bwdose05/n50)*(2^(1/alpha)-1))^-alpha
bwinfrisk05<-mc(bwrisk05)
#6-11 months
bwrisk611<- 1-(1+(bwdose611/n50)*(2^(1/alpha)-1))^-alpha
bwinfrisk611<-mc(bwrisk611)
#12-23 months
bwrisk1223<- 1-(1+(bwdose1223/n50)*(2^(1/alpha)-1))^-alpha
bwinfrisk1223<-mc(bwrisk1223)

#4-Drinking water (sbw)
#0-5 months
sbwrisk05<- 1-(1+(sbwdose05/n50)*(2^(1/alpha)-1))^-alpha
sbwinfrisk05<-mc(sbwrisk05)
#6-11 months
sbwrisk611<- 1-(1+(sbwdose611/n50)*(2^(1/alpha)-1))^-alpha
sbwinfrisk611<-mc(sbwrisk611)
#12-23 months
sbwrisk1223<- 1-(1+(sbwdose1223/n50)*(2^(1/alpha)-1))^-alpha
sbwinfrisk1223<-mc(sbwrisk1223)

#5-Drinking water (w)
#0-5 months
wrisk05<- 1-(1+(wdose05/n50)*(2^(1/alpha)-1))^-alpha
winfrisk05<-mc(wrisk05)
#6-11 months
wrisk611<- 1-(1+(wdose611/n50)*(2^(1/alpha)-1))^-alpha
winfrisk611<-mc(wrisk611)
#12-23 months
wrisk1223<- 1-(1+(wdose1223/n50)*(2^(1/alpha)-1))^-alpha
winfrisk1223<-mc(wrisk1223)

#6-Drinking water (sw)
#0-5 months
swrisk05<- 1-(1+(swdose05/n50)*(2^(1/alpha)-1))^-alpha
swinfrisk05<-mc(swrisk05)
#6-11 months
swrisk611<- 1-(1+(swdose611/n50)*(2^(1/alpha)-1))^-alpha
swinfrisk611<-mc(swrisk611)
#12-23 months
swrisk1223<- 1-(1+(swdose1223/n50)*(2^(1/alpha)-1))^-alpha
swinfrisk1223<-mc(swrisk1223)

#7- Bathing water 
#0-5 months
nalainfrisk05<-0
#6-11 months
nalainfrisk611<-0
#12-23 months
nalarisk1223<- 1-(1+(naladose1223/n50)*(2^(1/alpha)-1))^-alpha
nalainfrisk1223<-mc(nalarisk1223)


#CALCULATE 2-YEAR INFECTION RISKS FOR EACH PATHWAY, PER AGE GROUP
#SCENARIO 1: 608 DAYS DRINKING FROM BOREWELL, NO NALA
#Drinking water (bw)
n<-608
sc1bw2yearrisk<-1-(1-bwrisk611)^n
sc1bw2yearinfrisk<-mc(sc1bw2yearrisk)

#Drinking water (sbw)
n<-608
sc1sbw2yearrisk<-1-(1-sbwrisk611)^n
sc1sbw2yearinfrisk<-mc(sc1sbw2yearrisk)

#SCENARIO 2: 608 DAYS DRINKING FROM HP, 52 DAYS OF NALA
#Drinking water (hp)
n<-608
sc2hp2yearrisk<-1-(1-hprisk611)^n
sc2hp2yearinfrisk<-mc(sc2hp2yearrisk)

#Drinking water (shp)
n<-608
sc2shp2yearrisk<-1-(1-shprisk611)^n
sc2shp2yearinfrisk<-mc(sc2shp2yearrisk)

#Bathing water
n<-52
sc2nala2yearrisk<-1-(1-nalarisk1223)^n
sc2nala2yearinfrisk<-mc(sc2nala2yearrisk)

#SCENARIO 3: 430 DAYS DRINKING FROM HP, 178 DAYS DRINKING FROM WELL, 16 DAYS NALA
#Drinking water (hp+well)
nhp<-430
nwell<-178
sc32hpwellyearrisk<-1-(((1-hprisk611)^nhp)*((1-wrisk611)^nwell))
sc32hpwellyearinfrisk<-mc(sc32hpwellyearrisk)

#Drinking water (shp+swell)
nhp<-430
nwell<-178
sc3shpswell2yearrisk<-1-(((1-hprisk611)^nhp)*((1-wrisk611)^nwell))
sc3shpswell2yearinfrisk<-mc(sc3shpswell2yearrisk)

#Bathing water
n<-17
sc3nala2yearrisk<-1-(1-nalarisk1223)^n
sc3nala2yearinfrisk<-mc(sc3nala2yearrisk)


#OUTPUT RESULTS OF DAILY INFECTION RISKS
#1-Drinking water (hp)
#0-5 months
quantile(hpinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(hpinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(hpinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#2-Drinking water (shp)
#0-5 months
quantile(shpinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(shpinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(shpinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#3-Drinking water (bw)
#0-5 months
quantile(bwinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(bwinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(bwinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#4-Drinking water (sbw)
#0-5 months
quantile(sbwinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(sbwinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(sbwinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#4-Drinking water (w)
#0-5 months
quantile(winfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(winfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(winfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#5-Drinking water (sw)
#0-5 months
quantile(swinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
quantile(swinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
quantile(swinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#7-Bathing water
#0-5 months
#no ingestion
#6-11 months
#no ingestion
#12-23 months
quantile(nalainfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#OUTPUT RESULTS OF 2-YEAR INFECTION RISKS
#Scenario1
quantile(sc1bw2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
quantile(sc1sbw2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#Scenario2
quantile(sc2hp2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
quantile(sc2shp2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
quantile(sc2nala2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#Scenario3
quantile(sc32hpwellyearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
quantile(sc3shpswell2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
quantile(sc3nala2yearinfrisk, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))





