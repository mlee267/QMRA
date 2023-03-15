library(mc2d)
library(truncnorm)
ndvar(7e3) 
ndunc(7e3) 

###TABLE
results <- data.frame(matrix(nrow = 8, ncol = 12))
colnames(results) <- c("Age_Group", "Treatment", "Sampling_date", "Dose CFU","L95", "U95","Sub Clincal","L95", "U95","Clinical", "L95", "U95")
results$Age_Group<- rep(c(rep("50-61",4), rep("16-21",4)))
results$Treatment<- c(("Pre-Heat Shock"),("Post-Heat Shock"), ("Pre-Chlorine Dioxide Shock"), ("Post-Chlorine Dioxide Shock"),("Pre-Contionus Treatment"),("Post-Continous Treatment1"),("Post-Continous Treatment2"),("Post-Continous Treatment3"))
results$Sampling_date<-c(("May 17,2021"),("May 19,2021"), ("June 16,2021"), ("August 9,2021"),("February 2,2022"),("March 2,2022"),("April 1,2022"),("May 4,2022"))

### Pre-Heat Shock 
# breathing rate, m^3 air/min for (>50 to 61 population)
B <- mcstoc(runif, "U", min = 0.013, max = 0.017) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Pre-Heat Shock Sampling May 17, 2021)
# convert to L from 100 mL
original<-c((10*c(11629,7350,64690,12224,19226,1062,5549,99,373,436,12,53,231,940,2273,267,4770,6676,689,872,36,231,760,507)))
HS3Pre<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  HS3Pre[ ,i] <- rep_len(resample, ndvar())}
HS3Pre <- mcdata(HS3Pre, type="VU")
summary(HS3Pre)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
HS3Pre_dose<- ((wat_conv * HS3Pre)/Conver)
summary(HS3Pre_dose)
##Sub-clinical risk
rHS3PreSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PHS3Pre<- 1 -exp(-rHS3PreSubClin * HS3Pre_dose)
summary(PHS3Pre)
##Clinical Risk CSI
rHS3PreCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSIHS3Pre<- 1-exp(-rHS3PreCSI * HS3Pre_dose)
summary(CSIHS3Pre)
#Filling in Table
results[1, 4] <- median(apply(unmc(HS3Pre_dose), 2, mean))
results[1, 5] <- quantile(apply(unmc(HS3Pre_dose), 2, mean), probs =0.025)
results[1, 6] <- quantile(apply(unmc(HS3Pre_dose), 2, mean), probs = 0.975)
results[1, 7] <- median(apply(unmc(PHS3Pre), 2, mean))
results[1, 8] <- quantile(apply(unmc(PHS3Pre), 2, mean), probs = 0.025)
results[1, 9] <- quantile(apply(unmc(PHS3Pre), 2, mean), probs = 0.975)
results[1, 10] <- median(apply(unmc(CSIHS3Pre), 2, mean))
results[1, 11] <- quantile(apply(unmc(CSIHS3Pre), 2, mean), probs = 0.025)
results[1, 12] <- quantile(apply(unmc(CSIHS3Pre), 2, mean), probs = 0.975)


###Post-Heat Shock 
# breathing rate, m^3 air/min for (>50 to 61 population)
B <- mcstoc(runif, "U", min = 0.013, max = 0.017) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Post-Heat Shock Sampling May 19, 2021)
# convert to L from 100 mL
original<-c(10*c(79,4,1163,53,13,99,231,10,7,9,373,10),rep(0,12))
HS3Post<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  HS3Post[ ,i] <- rep_len(resample, ndvar())}
HS3Post <- mcdata(HS3Post, type="VU")
summary(HS3Post)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
HS3Post_dose<- ((wat_conv * HS3Post)/Conver)
summary(HS3Post_dose)
##Sub-clinical risk
rHS3PostSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PHS3Post<- 1-exp(-rHS3PostSubClin * HS3Post_dose)
summary(PHS3Post)
##Clinical Risk CSI
rHS3PostCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSIHS3Post<- 1-exp(-rHS3PostCSI * HS3Post_dose)
summary(CSIHS3Post)

results[2, 4] <- median(apply(unmc(HS3Post_dose), 2, mean))
results[2, 5] <- quantile(apply(unmc(HS3Post_dose), 2, mean), probs =0.025)
results[2, 6] <- quantile(apply(unmc(HS3Post_dose), 2, mean), probs = 0.975)
results[2, 7] <- median(apply(unmc(PHS3Post), 2, mean))
results[2, 8] <- quantile(apply(unmc(PHS3Post), 2, mean), probs = 0.025)
results[2, 9] <- quantile(apply(unmc(PHS3Post), 2, mean), probs = 0.975)
results[2, 10] <- median(apply(unmc(CSIHS3Post), 2, mean))
results[2, 11] <- quantile(apply(unmc(CSIHS3Post), 2, mean), probs = 0.025)
results[2, 12] <- quantile(apply(unmc(CSIHS3Post), 2, mean), probs = 0.975)


##Pre-chlorine dioxide 
# breathing rate, m^3 air/min for (>50 to 61 population)
B <- mcstoc(runif, "U", min = 0.013, max = 0.017) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Pre-Chlorine Dioxide Sampling June 16, 2021)
# convert to L from 100 mL
original<-c(10*c(22726,13677,22726,14598,15725,12898,6081,8723,19226,2273,12898,3849,13677,3501,17178,22726,19226,15,17178,14598,2580,17178),rep(0,2))
CDS1Pre<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CDS1Pre[ ,i] <- rep_len(resample, ndvar())}
CDS1Pre <- mcdata(CDS1Pre, type="VU")
summary(CDS1Pre)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CDS1Pre_dose<- ((wat_conv * CDS1Pre)/Conver)
summary(CDS1Pre_dose)
##sub-clinical risk
rCDS1PreSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCDS1Pre<- 1-exp(-rCDS1PreSubClin * CDS1Pre_dose)
summary(PCDS1Pre)
##Clinical Risk CSI
rCDS1PreCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICDS1Pre<- 1-exp(-rCDS1PreCSI * CDS1Pre_dose)
summary(CSICDS1Pre)

results[3, 4] <- median(apply(unmc(CDS1Pre_dose), 2, mean))
results[3, 5] <- quantile(apply(unmc(CDS1Pre_dose), 2, mean), probs =0.025)
results[3, 6] <- quantile(apply(unmc(CDS1Pre_dose), 2, mean), probs = 0.975)
results[3, 7] <- median(apply(unmc(PCDS1Pre), 2, mean))
results[3, 8] <- quantile(apply(unmc(PCDS1Pre), 2, mean), probs = 0.025)
results[3, 9] <- quantile(apply(unmc(PCDS1Pre), 2, mean), probs = 0.975)
results[3, 10] <- median(apply(unmc(CSICDS1Pre), 2, mean))
results[3, 11] <- quantile(apply(unmc(CSICDS1Pre), 2, mean), probs = 0.025)
results[3, 12] <- quantile(apply(unmc(CSICDS1Pre), 2, mean), probs = 0.975)


##Post-Chlorine Dioxide 
# breathing rate, m^3 air/min for (>50 to 61 population)
B <- mcstoc(runif, "U", min = 0.013, max = 0.017) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Post-Chlorine Dioxide Sampling August 8,2021)
# convert to L from 100 mL
original<-c(10*c(22,13,7,53,16,258,16,2,36,9,4,8,12,12,4,31,3,99), rep(0,6))
CD1Post<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CD1Post[ ,i] <- rep_len(resample, ndvar())}
CD1Post <- mcdata(CD1Post, type="VU")
summary(CD1Post)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CD1Post_dose<- ((wat_conv * CD1Post)/Conver)
summary(CD1Post_dose)
##sub-clinical risk
rCD1PostSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCD1Post<- 1- exp(-rCD1PostSubClin * CD1Post_dose)
summary(PCD1Post)
##Clinical Risk CSI
rCD1PostCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICD1Post<- 1- exp(-rCD1PostCSI * CD1Post_dose)
summary(CSICD1Post)

results[4, 4] <- median(apply(unmc(CD1Post_dose), 2, mean))
results[4, 5] <- quantile(apply(unmc(CD1Post_dose), 2, mean), probs =0.025)
results[4, 6] <- quantile(apply(unmc(CD1Post_dose), 2, mean), probs = 0.975)
results[4, 7] <- median(apply(unmc(PCD1Post), 2, mean))
results[4, 8] <- quantile(apply(unmc(PCD1Post), 2, mean), probs = 0.025)
results[4, 9] <- quantile(apply(unmc(PCD1Post), 2, mean), probs = 0.975)
results[4, 10] <- median(apply(unmc(CSICD1Post), 2, mean))
results[4, 11] <- quantile(apply(unmc(CSICD1Post), 2, mean), probs = 0.025)
results[4, 12] <- quantile(apply(unmc(CSICD1Post), 2, mean), probs = 0.975)



####PRE chlorine dioxide continuous
# breathing rate, m^3 air/min for (16 to <21 population)
B <- mcstoc(runif, "U", min = 0.012, max = 0.016) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Pre-Continuous Chlorine Dioxide Sampling February 2,2022)
# convert to L from 100 mL
original<-c(10*c(3501.00 ,	328.10 ,	759.60 ,	258.00,	317.50 ,	5.20 ,	8.40 ,	85.40 ,	26.40 ,	12.30 ,	1.10 ,	2.00 ,	41.60  ,	15.50 ,	47.40 ,	6.60 ,	239.60 ,	1.10  ,	3.50 ,	7.40 ), rep (0,4))
CTPre<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CTPre[ ,i] <- rep_len(resample, ndvar())}
CTPre <- mcdata(CTPre, type="VU")
summary(CTPre)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CTPre_dose<- ((wat_conv * CTPre)/Conver)
summary(CTPre_dose)
##subclinical risk
rCTPreSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCDS1Pre<- 1-exp(-rCTPreSubClin * CTPre_dose)
summary(PCDS1Pre)
##Clinical Risk CSI
rCTPreCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICTS1Pre<- 1-exp(-rCTPreCSI* CTPre_dose)
summary(CSICTS1Pre)

results[5, 4] <- median(apply(unmc(CTPre_dose), 2, mean))
results[5, 5] <- quantile(apply(unmc(CTPre_dose), 2, mean), probs =0.025)
results[5, 6] <- quantile(apply(unmc(CTPre_dose), 2, mean), probs = 0.975)
results[5, 7] <- median(apply(unmc(PCDS1Pre), 2, mean))
results[5, 8] <- quantile(apply(unmc(PCDS1Pre), 2, mean), probs = 0.025)
results[5, 9] <- quantile(apply(unmc(PCDS1Pre), 2, mean), probs = 0.975)
results[5, 10] <- median(apply(unmc(CSICTS1Pre), 2, mean))
results[5, 11] <- quantile(apply(unmc(CSICTS1Pre), 2, mean), probs = 0.025)
results[5, 12] <- quantile(apply(unmc(CSICTS1Pre), 2, mean), probs = 0.975)


##Post-Chlorine dioxide continuous  
# breathing rate, m^3 air/min for (16 to <21 population)
B <- mcstoc(runif, "U", min = 0.012, max = 0.016) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Post-Continuous Chlorine Dioxide 1 Sampling March 2, 2022)
# convert to L from 100 mL
original<-c(10*c(10.40 ,	1.10 ,	448.90 ,	85.40 ,	221.90 ,	5.20 ,	1.00 ,	1.10 ,	98.90 ,	1.10 ,	10.60  ,	14.10 ,	8.40 ,	8.40),rep(0,10))
CT1Post<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CT1Post[ ,i] <- rep_len(resample, ndvar())}
CT1Post <- mcdata(CT1Post, type="VU")
summary(CT1Post)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CT1Post_dose<- ((wat_conv * CT1Post)/Conver)
summary(CT1Post_dose)
##subclinical risk
rCT1PostSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCT1Post<- 1- exp(-rCT1PostSubClin * CT1Post_dose)
summary(PCT1Post)
##Clinical Risk CSI
rCT1PostCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICT1Post<- 1- exp(-rCT1PostCSI * CT1Post_dose)
summary(CSICT1Post)

results[6, 4] <- median(apply(unmc(CT1Post_dose), 2, mean))
results[6, 5] <- quantile(apply(unmc(CT1Post_dose), 2, mean), probs =0.025)
results[6, 6] <- quantile(apply(unmc(CT1Post_dose), 2, mean), probs = 0.975)
results[6, 7] <- median(apply(unmc(PCT1Post), 2, mean))
results[6, 8] <- quantile(apply(unmc(PCT1Post), 2, mean), probs = 0.025)
results[6, 9] <- quantile(apply(unmc(PCT1Post), 2, mean), probs = 0.975)
results[6, 10] <- median(apply(unmc(CSICT1Post), 2, mean))
results[6, 11] <- quantile(apply(unmc(CSICT1Post), 2, mean), probs = 0.025)
results[6, 12] <- quantile(apply(unmc(CSICT1Post), 2, mean), probs = 0.975)

##Post-chlorine dioxide continuous post 2 
# breathing rate, m^3 air/min for (16 to <21 population)
B <- mcstoc(runif, "U", min = 0.012, max = 0.016) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Post-Continuous Chlorine Dioxide 2 Sampling April 1, 2022)
# convert to L from 100 mL
original<-c(10*c(6081.00 ,	3614.00 ,	2580.00 ,	3175.00 ,	9049.00 ,	5223.00 ,	41.60 ,	4.70  ,	16.90 ,	26.40 ,	1.00 ,	4.70  ,	105.70 ), rep(0,11))
CT2Post<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CT2Post[ ,i] <- rep_len(resample, ndvar())}
CT2Post <- mcdata(CT2Post, type="VU")
summary(CT2Post)

# convert MPN to CFU (1 MPN =1.1 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CT2Post_dose<- ((wat_conv * CT2Post)/Conver)
summary(CT2Post_dose)
##subclinical risk
rCT2PostSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCT2Post<- 1- exp(-rCT2PostSubClin * CT2Post_dose)
summary(PCT2Post)
##Clinical Risk CSI
rCT2PostCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICT2Post<- 1- exp(-rCT2PostCSI * CT2Post_dose)
summary(CSICT2Post)

results[7, 4] <- median(apply(unmc(CT2Post_dose), 2, mean))
results[7, 5] <- quantile(apply(unmc(CT2Post_dose), 2, mean), probs =0.025)
results[7, 6] <- quantile(apply(unmc(CT2Post_dose), 2, mean), probs = 0.975)
results[7, 7] <- median(apply(unmc(PCT2Post), 2, mean))
results[7, 8] <- quantile(apply(unmc(PCT2Post), 2, mean), probs = 0.025)
results[7, 9] <- quantile(apply(unmc(PCT2Post), 2, mean), probs = 0.975)
results[7, 10] <- median(apply(unmc(CSICT2Post), 2, mean))
results[7, 11] <- quantile(apply(unmc(CSICT2Post), 2, mean), probs = 0.025)
results[7, 12] <- quantile(apply(unmc(CSICT2Post), 2, mean), probs = 0.975)

##POST chlorine dioxide continuous post 3 Population 
# breathing rate, m^3 air/min for (16 to <21 population)
B <- mcstoc(runif, "U", min = 0.012, max = 0.016) 
#showertime, min
tsh <- mcstoc(rtruncnorm, "U", a = 0, mean = 7.8, sd = 0.02) 
# Aerosol concentrations, no./m^3 air (conventional shower fixtures)
Caer1 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.30)
Caer2 <- mcstoc(rlnorm, "U", meanlog = 17.5, sdlog = 0.17)
Caer3 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer4 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer5 <- mcstoc(rlnorm, "U", meanlog = 19.4, sdlog = 0.35)
Caer6 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer7 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer8 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer9 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
Caer10 <- mcstoc(rlnorm, "U", meanlog = 20.0, sdlog = 0.31)
# Aerosol volumes, m^3 water 
Vaer1 <- (4 / 3) * pi * (((1 / 2) * (10 ^ -6)) ^ 3)
Vaer2 <- (4 / 3) * pi * (((2 / 2) * (10 ^ -6)) ^ 3)
Vaer3 <- (4 / 3) * pi * (((3 / 2) * (10 ^ -6)) ^ 3)
Vaer4 <- (4 / 3) * pi * (((4 / 2) * (10 ^ -6)) ^ 3)
Vaer5 <- (4 / 3) * pi * (((5 / 2) * (10 ^ -6)) ^ 3)
Vaer6 <- (4 / 3) * pi * (((6 / 2) * (10 ^ -6)) ^ 3)
Vaer7 <- (4 / 3) * pi * (((7 / 2) * (10 ^ -6)) ^ 3)
Vaer8 <- (4 / 3) * pi * (((8 / 2) * (10 ^ -6)) ^ 3)
Vaer9 <- (4 / 3) * pi * (((9 / 2) * (10 ^ -6)) ^ 3)
Vaer10 <- (4 / 3) * pi * (((10 / 2) * (10 ^ -6)) ^ 3)
# Percentage of total Legionella per size fraction
F1 <- 17.5 / 100
F2 <- 16.39 / 100
F3 <- 15.56 / 100
F4 <- 6.67 / 100
F5 <- 3.89 / 100
F6 <- 2.5 / 100
F7 <- 2.78 / 100
F8 <- 5.00 / 100
F9 <- 5.28 / 100
F10 <- 3.89 / 100
# Deposition efficiency of aerosols by size category
D1 <- mcstoc(runif, "U", min = 0.23, max = 0.25)
D2 <- mcstoc(runif, "U", min = 0.40, max = 0.53)
D3 <- mcstoc(runif, "U", min = 0.36, max = 0.62)
D4 <- mcstoc(runif, "U", min = 0.29, max = 0.61)
D5 <- mcstoc(runif, "U", min = 0.19, max = 0.52)
D6 <- mcstoc(runif, "U", min = 0.10, max = 0.40)
D7 <- mcstoc(runif, "U", min = 0.06, max = 0.29)
D8 <- mcstoc(runif, "U", min = 0.03, max = 0.19)
D9 <- mcstoc(runif, "U", min = 0.01, max = 0.12)
D10 <- mcstoc(runif, "U", min = 0.01, max = 0.06)
# Total aerosol concentration, m^3 water/m^3 air (conventional)
CVaer_conv <- (Caer1 * Vaer1) + (Caer2 * Vaer2) + (Caer3 * Vaer3) +(Caer4 * Vaer4) + (Caer5 * Vaer5) + (Caer6 * Vaer6) + (Caer7 * Vaer7) +(Caer8 * Vaer8) + (Caer9 * Vaer9) + (Caer10 * Vaer10)
# F: Inhalable proportion of total Legionella 
FD <- (F1 * D1) + (F2 * D2) + (F3 * D3) + (F4 * D4) + (F5 * D5) + (F6 * D6) +(F7 * D7) + (F8 * D8) + (F9 * D9) + (F10 * D10)
# Equivalent inhaled aerosol volume, L (conventional)
wat_conv <- B * tsh * (CVaer_conv * 1e3) * FD

#Culturable L.pneumophila concentrations values from (Post-Continuous Chlorine Dioxide 3 Sampling May 4, 2022)
# convert to L from 100 mL
original<-c(10*c(	155.00  ,	27.50,	104.00),rep(0,21))
CT3Post<-matrix(nrow = ndvar(), ncol = ndunc())
for (i in 1:ndunc()) {
  resample<- sample(original, length(original), replace= TRUE) 
  CT3Post[ ,i] <- rep_len(resample, ndvar())}
CT3Post <- mcdata(CT3Post, type="VU")
summary(CT3Post)

# convert MPN to CFU (1 MPN =1.2 CFU)
Conver<- mcstoc(rtruncnorm, "U", a = 0, mean = 1.1, sd = 0.16)
CT3Post_dose<- ((wat_conv * CT3Post)/Conver)
summary(CT3Post_dose)
##subclinical risk
rCT3PostSubClin<-mcstoc(rlnorm, type="U", meanlog= -2.93, sdlog =0.49)
PCT3Post<- 1- exp(-rCT3PostSubClin * CT3Post_dose)
summary(PCT3Post)
##Clinical Risk CSI
rCT3PostCSI<-mcstoc(rlnorm, type="U", meanlog = -9.69, sdlog= 0.30)
CSICT3Post<- 1- exp(-rCT3PostCSI * CT3Post_dose)
summary(CSICT3Post)

results[8, 4] <- median(apply(unmc(CT3Post_dose), 2, mean))
results[8, 5] <- quantile(apply(unmc(CT3Post_dose), 2, mean), probs =0.025)
results[8, 6] <- quantile(apply(unmc(CT3Post_dose), 2, mean), probs = 0.975)
results[8, 7] <- median(apply(unmc(PCT3Post), 2, mean))
results[8, 8] <- quantile(apply(unmc(PCT3Post), 2, mean), probs = 0.025)
results[8, 9] <- quantile(apply(unmc(PCT3Post), 2, mean), probs = 0.975)
results[8, 10] <- median(apply(unmc(CSICT3Post), 2, mean))
results[8, 11] <- quantile(apply(unmc(CSICT3Post), 2, mean), probs = 0.025)
results[8, 12] <- quantile(apply(unmc(CSICT3Post), 2, mean), probs = 0.975)

library(data.table)
write_csv(results,"C:\\Users\\monica\\Desktop\\QMRA Results\\R Code\\Results_Culture.csv")


