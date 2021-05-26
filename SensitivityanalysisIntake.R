#Sensitivity- Ingestion
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
cleands<- subset(dataset, (Methoderror!=4 & Methoderror!=3)) #Method error 4= damaged sample, could not count, Method error 3=tmtc and too undiluted-->we exclude them from analysis
cleands$Village[cleands$Village=="DevdaSath"]<-"Devdasath" #fixing some small issues with varnames
cleands$Block[cleands$Block=="Ghatol "]<-"Ghatol" #fixing some small issues with varnames
#Now imputate censored data
cleands$CFU_100mL[cleands$CFU_100mL==0]<-1 #substituting non-detects for minimum detection limit of 1cfu
cleands$Log10CFU<-log10(cleands$CFU_100mL) #recalculate logarithmic scale to include imputed non-detects


###########   MONTECARLO MODEL ###########
#E. coli concentrations found in each exposure pathway
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

####MODEL
ndvar(10000) #number of samples for the variability dimension. 10000 iterations along the prob density functions
nsu(100) #number of samples for the uncertainty dimension. 100 iterations along the uncertatainty distribution

#E. COLI CONCENTRATIONS DISTRIBUTIONS
#Define E.coli concentration distributions in each exposure pathway
sourcewaterdistr<-fitdistr(sourcewater$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
storedwaterdistr<-fitdistr(storedwater$CFU_100mL/100, "lognormal") # ln N(logmean, logsd) CFU/mL 
naladistr<-fitdistr(nala$CFU_100mL/100, "lognormal")# ln N(logmean, logsd) CFU/mL 
childhanddistr<-fitdistr(childhand$CFU_100mL, "lognormal") #ln N(logmean, logsd) CFU/hand
caregiverhanddistr<-fitdistr(caregiverhand$CFU_100mL, "lognormal") #ln N(logmean, logsd) CFU/hand
soildistr<-fitdistr(soil$CFU_100mL*1000, "lognormal") #ln N(logmean, logsd) CFU/mg

#SIMULATIONS ARE ASSIGNED RANDOM E.COLI CONCENTRATION AS PER THE DISTRIBUTIONS FITTED
#mcstoc creates a matrix of 10000*100 simulations, to which random values from the e.coli distributions are assigned
sourcewaterconc<-mcstoc(rlnorm,type="VU",meanlog=sourcewaterdistr$estimate['meanlog'],sdlog=sourcewaterdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
storedwaterconc<-mcstoc(rlnorm,type="VU",meanlog=storedwaterdistr$estimate['meanlog'],sdlog=storedwaterdistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
nalaconc<-mcstoc(rlnorm,type="VU",meanlog=naladistr$estimate['meanlog'],sdlog=naladistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
childhandconc<-mcstoc(rlnorm,type="VU",meanlog=childhanddistr$estimate['meanlog'],sdlog=childhanddistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
caregiverhandconc<-mcstoc(rlnorm,type="VU",meanlog=caregiverhanddistr$estimate['meanlog'],sdlog=caregiverhanddistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability
soilconc<-mcstoc(rlnorm,type="VU",meanlog=soildistr$estimate['meanlog'],sdlog=soildistr$estimate['sdlog'], seed=41) #distribution of pathogen concentration with "V"ariability

#INGESTION DISTRIBUTION PARAMETERS, PER AGE GROUP

#******************ALTER SENSITIVITY PARAMETERS***************
sensitivitylogchange<- +1
sensitivityweibullchange<-+10
#******************ALTERING SENSITIVITY PARAMETERS***************
#*
#1-ingestion of drinking water (ml/day)
#0-5 months
logmeandrinkingwateringestion05<-3.37+sensitivitylogchange
logsddrinkingwateringestion05<-1.23
#6-11 months
logmeandrinkingwateringestion611<-5.24+sensitivitylogchange
logsddrinkingwateringestion611<-0.85
#12-23 months
logmeandrinkingwateringestion1223<-4.58+sensitivitylogchange
logsddrinkingwateringestion1223<-0.89

#2-ingestion of bathing water (ml/day)
#0-5 months
#logmeanbathingwateringestion05<-0
#logsdbathingwateringestion05<-0
#6-11 months
#logmeanbathingwateringestion611<-0
#logsdbathingwateringestion611<-0
#12-23 months
logmeanbathingwateringestion1223<-3.49+sensitivitylogchange
logsdbathingwateringestion1223<-0.55

#3-ingestion of soil (mg/day)
#0-5 months
logmeansoilingestion05<-2.73+sensitivitylogchange
logsdsoilingestion05<-0.69
#6-11 months
logmeansoilingestion611<-4.26+sensitivitylogchange
logsdsoilingestion611<-0.69
#12-23 months
logmeansoilingestion1223<-4.25+sensitivitylogchange
logsdsoilingestion1223<-0.69

#4i5-frequency of hand-to-mouth contact (handfraction/day)
#0-5 months
shape05<-1.58
scale05<-58.22+sensitivityweibullchange
#6-11 months
shape611<-1.58
scale611<-58.22+sensitivityweibullchange
#12-23 months
shape1223<-2.66
scale1223<-56.62+sensitivityweibullchange


#SIMULATIONS ARE ASSIGNED RANDOM INGESTION VOLUMES AS PER THE DISTRIBUTIONS FOUND
#1-Drinking water (stored water)
#0-5 months
storedwatering05<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion05,sdlog=logsddrinkingwateringestion05, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#6-11 months
storedwatering611<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion611,sdlog=logsddrinkingwateringestion611, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)
#12-23 months
storedwatering1223<-mcstoc(rlnorm,type="VU",meanlog=logmeandrinkingwateringestion1223,sdlog=logsddrinkingwateringestion1223, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#2-Bathing water 
#0-5 months
#bathingwatering05<-no ingestion
#6-11 months
#bathingwatering611<-noingestion
#12-23 months
bathingwatering1223<-mcstoc(rlnorm,type="VU",meanlog=logmeanbathingwateringestion1223,sdlog=logsddrinkingwateringestion1223, seed=41) #lognorm distr of daily ingestion of pathogens (ml/day), with a min and max, with Variability and Uncertainty (VU)

#3-ingestion of soil
#0-5 months
soiling05<-mcstoc(rlnorm,type="VU",meanlog=logmeansoilingestion05,sdlog=logsdsoilingestion05, seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)
#6-11 months
soiling611<-mcstoc(rlnorm,type="VU",meanlog=logmeansoilingestion611,sdlog=logsdsoilingestion611, seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)
#12-23 months
soiling1223<-mcstoc(rlnorm,type="VU",meanlog=logmeansoilingestion1223,sdlog=logsdsoilingestion1223, seed=41) #lognorm distr of daily ingestion of pathogens (mg/day), with Variability and Uncertainty (VU)

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
#************************************
#************************************
ratio<-0.08 #e.coli to pathogen ratio. E.coli0157=0.08, campylobacter=0.66, rotavirus=0.00005
alpha<-0.2099 #bpoisson pathogen parameter. E.coli0157=0.2099, campylobacter=0.145, rotavirus=0.2531
n50<-1120 #bpoisson pathogen parameter. E.coli0157=1120, campylobacter=895.52, rotavirus=6.16
#*****************************
#*****************************
#*****************************

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

#CALCULATE 2-YEAR INFECTION RISKS FOR EACH PATHWAY, PER AGE GROUP
#1-Drinking water (stored water)
n<-608
storedwater2yearrisk<-1-(1-storedwaterrisk611)^n
storedwater2yearinfrisk<-mc(storedwater2yearrisk)

#2-Bathing water
n<-52
bathingwater2yearrisk<-1-(1-bathingwaterrisk1223)^n
bathingwater2yearinfrisk<-mc(bathingwater2yearrisk)

#3-Soil
n<-517
soil2yearrisk<-1-(1-soilrisk611)^n
soil2yearinfrisk<-mc(soil2yearrisk)

#Mouthing childs hands
n<-730
childhand2yearrisk<-1-(1-childhandrisk611)^n
childhand2yearinfrisk<-mc(childhand2yearrisk)

#Mouthing caregiver's hands
n<-730
caregiverhand2yearrisk<-1-(1-caregiverhandrisk611)^n
caregiverhand2yearinfrisk<-mc(caregiverhand2yearrisk)

#OUTPUT RESULTS OF DAILY INFECTION RISKS
#1-Drinking water (stored water)
#0-5 months
#print(storedwaterinfrisk05)
#quantile(storedwaterinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
#print(storedwaterinfrisk611)
#quantile(storedwaterinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
print(storedwaterinfrisk1223)
quantile(storedwaterinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#2-Bathing water
#0-5 months
#no ingestion
#6-11 months
#no ingestion
#12-23 months
print(bathingwaterinfrisk1223)
quantile(bathingwaterinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#3-Soil
#0-5 months
#print(soilinfrisk05)
#quantile(soilinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
#print(soilinfrisk611)
#quantile(soilinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
print(soilinfrisk1223)
quantile(soilinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#4-Child's hands
#0-5 months
#print(childhandinfrisk05)
#quantile(childhandinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
#print(childhandinfrisk611)
#quantile(childhandinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
print(childhandinfrisk1223)
quantile(childhandinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#5-Caregiver's hands
#0-5 months
#print(caregiverhandinfrisk05)
#quantile(caregiverhandinfrisk05, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#6-11 months
#print(caregiverhandinfrisk611)
#quantile(caregiverhandinfrisk611, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))
#12-23 months
print(caregiverhandinfrisk1223)
quantile(caregiverhandinfrisk1223, c(0.025,0.1,0.25,0.5,0.75, 0.9, 0.975))

#OUTPUT RESULTS OF 2-YEAR INFECTION RISKS
#summary(storedwater2yearinfrisk)
#summary(bathingwater2yearinfrisk)
#summary(soil2yearinfrisk)
#summary(childhand2yearinfrisk)
#summary(caregiverhand2yearinfrisk)



