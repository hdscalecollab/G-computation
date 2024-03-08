#G Computation Simulation Code for Effect Modification by Sex
#Manuscript "Simulating the Impact of Green Space Exposure on Cardiometabolic Biomarkers in a
#Diverse Population Living in San Diego, California: A G-Computation Application"

#Adapted from code provided by Hernán MA, Robins JM (2020). Causal Inference: What If. Boca Raton: Chapman & Hall/CRC.
#R code by Joy Shi and Sean McGrath
#Original code can be found here: https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/

#Version 3.6.2024

# Contact: Anais Teyton, ateyton@ucsd.edu

#########################################
# Calculating Risk Differences by Age
# Example provided is for Blood Glucose Level Outcome
# Can be substituted for any of the other 6 biomarkers
#########################################

#Libraries
library(readxl)
library(tidyr)
library(dplyr)
library(boot)
library("ggplot2") 
library("ggthemes") 

#########################################
# 65 Years and over for Simulated Exposure
# <65 Years old exposed to true NDVI value
#########################################

#Set working directory and upload clean sample data
setwd("ENTER WORKING DIRECTORY HERE")
pq <- read_excel("ENTER DATASET HERE")
pq = subset(pq, select = c(pq_study_id_unique, ndvi_mean_400m, cholesterol_result, ldl_calculated_result, triglyceride_result, hemoglobin_a1c_result, glucose_mgdl, sys, dia, sex, incomecat, racecat, age, educcat, lat, bmi, hdl_result, m1_ws_cm, mets) )

#DATA DESCRIPTION PROVIDED IN CODE #1
pq$sex <- as.factor(pq$sex)
pq$incomecat <- as.factor(pq$incomecat)
pq$racecat <- as.factor(pq$racecat)
pq$educcat <-as.factor(pq$educcat)
pq$lat <- as.factor(pq$lat)

pq=rename(pq, intervention_b = ndvi_mean_400m)

#Creating standardization function
standardization <- function(data, indices) {
  
  # create a dataset with 11 copies of each subject
  d <- data[indices,]
  d$interv <- -1 # 1st copy: equal to original one
  
  d1 <- d # 2nd copy (Minimum): intervention set to 1, outcome to missing
  d1$interv <- 1
  d1$intervention_b [d1$age >= 65] <- -.0815979
  d1$glucose_mgdl <- NA
  
  d2 <- d # 3rd copy (10th percentile): intervention set to 2, outcome to missing
  d2$interv <- 2
  d2$intervention_b [d2$age >= 65]  <- .086611807346344
  d2$glucose_mgdl <- NA
  
  d3 <- d # 4th copy (20th percentile): intervention set to 3, outcome to missing
  d3$interv <- 3
  d3$intervention_b [d3$age >= 65]  <- .1140746735036373
  d3$glucose_mgdl <- NA
  
  d4 <- d # 5th copy (30th percentile): intervention set to 4, outcome to missing
  d4$interv <- 4
  d4$intervention_b [d4$age >= 65]  <- .133780911564827
  d4$glucose_mgdl <- NA
  
  d5 <- d # 6th copy (40th percentile): intervention set to 5, outcome to missing
  d5$interv <- 5
  d5$intervention_b [d5$age >= 65]  <- .146772064268589
  d5$glucose_mgdl <- NA
  
  d6 <- d # 7th copy (50th percentile): intervention set to 6, outcome to missing
  d6$interv <- 6
  d6$intervention_b [d6$age >= 65]  <- .1665423661470413
  d6$glucose_mgdl <- NA
  
  d7 <- d # 8th copy (60th percentile): intervention set to 7, outcome to missing
  d7$interv <- 7
  d7$intervention_b  [d7$age >= 65]  <- .1879632100462914
  d7$glucose_mgdl <- NA
  
  d8 <- d # 9th copy (70th percentile): intervention set to 8, outcome to missing
  d8$interv <- 8
  d8$intervention_b  [d8$age >= 65]  <- .2076895833015442
  d8$glucose_mgdl <- NA
  
  d9 <- d # 10th copy (80th percentile): intervention set to 9, outcome to missing
  d9$interv <- 9
  d9$intervention_b [d9$age >= 65]   <- .2252189218997955
  d9$glucose_mgdl <- NA
  
  d10 <- d # 11th copy (90th percentile): intervention set to 10, outcome to missing
  d10$interv <- 10
  d10$intervention_b [d10$age >= 65]  <- .259517103433609
  d10$glucose_mgdl <- NA
  
  d11 <- d # 12th copy (Maximum): intervention set to 11, outcome to missing
  d11$interv <- 11
  d11$intervention_b [d11$age >= 65]  <- .3624528
  d11$glucose_mgdl <- NA
  
  d.onesample <- rbind(d, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11) # combining datasets
  
  
  # linear regression model to estimate mean outcome conditional on treatment and confounders
  model <- glm(glucose_mgdl ~ intervention_b + age + intervention_b*age + lat + racecat + incomecat + educcat + sex, data=d.onesample)
  
  summary(model)
  
  d.onesample$predicted_meanY <- predict(model, d.onesample)
  
  # estimate mean outcome in each of the groups
  # estimate difference in mean outcome in each of the 'exposure' groups (deciles + max NDVI) vs. 'unexposed' group (min NDVI)
  
  return(c(
    mean(d.onesample$predicted_meanY[d.onesample$interv == -1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 2]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 3]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 4]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 5]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 6]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 7]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 8]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 9]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 10]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 11]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 2])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 3])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 4])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 5])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 6])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 7])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 8])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 9])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 10])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 11])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1])
  ))
  
}

# bootstrap
results <- boot(data = pq, statistic = standardization, R = 1000)

# generating confidence intervals
se <- c(sd(results$t[, 1]),
        sd(results$t[, 2]),
        sd(results$t[, 3]),
        sd(results$t[, 4]),
        sd(results$t[, 5]),
        sd(results$t[, 6]),
        sd(results$t[, 7]),
        sd(results$t[, 8]),
        sd(results$t[, 9]),
        sd(results$t[, 10]),
        sd(results$t[, 11]),
        sd(results$t[, 12]),
        sd(results$t[, 13]),
        sd(results$t[, 14]),
        sd(results$t[, 15]),
        sd(results$t[, 16]),
        sd(results$t[, 17]),
        sd(results$t[, 18]),
        sd(results$t[, 19]),
        sd(results$t[, 20]),
        sd(results$t[, 21]),
        sd(results$t[, 22]))

mean <- results$t0
ll <- mean - qnorm(0.975) * se
ul <- mean + qnorm(0.975) * se

bootstrap_gluc_ageover65 <-
  data.frame(cbind(
    c("Observed",
      "Minimum",
      "10th Percentile",
      "20th Percentile",
      "30th Percentile",
      "40th Percentile",
      "50th Percentile",
      "60th Percentile",
      "70th Percentile",
      "80th Percentile",
      "90th Percentile",
      "Maximum",
      "RD 10th - Min",
      "RD 20th - Min",
      "RD 30th - Min",
      "RD 40th - Min",
      "RD 50th - Min",
      "RD 60th - Min",
      "RD 70th - Min",
      "RD 80th - Min",
      "RD 90th - Min",
      "RD Max - Min"),
    mean,
    se,
    ll,
    ul
  ))
bootstrap_gluc_ageover65

#########################################
# <65 Years for Simulated Exposure
# 65 Years and over exposed to true NDVI value
#########################################

#Set working directory and upload clean sample data
setwd("ENTER WORKING DIRECTORY HERE")
pq <- read_excel("ENTER DATASET HERE")
pq = subset(pq, select = c(pq_study_id_unique, ndvi_mean_400m, cholesterol_result, ldl_calculated_result, triglyceride_result, hemoglobin_a1c_result, glucose_mgdl, sys, dia, sex, incomecat, racecat, age, educcat, lat, bmi, hdl_result, m1_ws_cm, mets) )

#DATA DESCRIPTION PROVIDED IN CODE #1
pq$sex <- as.factor(pq$sex)
pq$incomecat <- as.factor(pq$incomecat)
pq$racecat <- as.factor(pq$racecat)
pq$educcat <-as.factor(pq$educcat)
pq$lat <- as.factor(pq$lat)

pq=rename(pq, intervention_b = ndvi_mean_400m)

#Creating standardization function
standardization <- function(data, indices) {
  
  # create a dataset with 11 copies of each subject
  d <- data[indices,]
  d$interv <- -1 # 1st copy: equal to original one
  
  d1 <- d # 2nd copy (Minimum): intervention set to 1, outcome to missing
  d1$interv <- 1
  d1$intervention_b [d1$age < 65] <- -0815979
  d1$glucose_mgdl <- NA
  
  d2 <- d # 3rd copy (10th percentile): intervention set to 2, outcome to missing
  d2$interv <- 2
  d2$intervention_b [d2$age < 65]  <- .086611807346344
  d2$glucose_mgdl <- NA
  
  d3 <- d # 4th copy (20th percentile): intervention set to 3, outcome to missing
  d3$interv <- 3
  d3$intervention_b [d3$age < 65]  <- .1140746735036373
  d3$glucose_mgdl <- NA
  
  d4 <- d # 5th copy (30th percentile): intervention set to 4, outcome to missing
  d4$interv <- 4
  d4$intervention_b [d4$age < 65]  <- .133780911564827
  d4$glucose_mgdl <- NA
  
  d5 <- d # 6th copy (40th percentile): intervention set to 5, outcome to missing
  d5$interv <- 5
  d5$intervention_b [d5$age < 65]  <- .146772064268589
  d5$glucose_mgdl <- NA
  
  d6 <- d # 7th copy (50th percentile): intervention set to 6, outcome to missing
  d6$interv <- 6
  d6$intervention_b [d6$age < 65]  <- .1665423661470413
  d6$glucose_mgdl <- NA
  
  d7 <- d # 8th copy (60th percentile): intervention set to 7, outcome to missing
  d7$interv <- 7
  d7$intervention_b  [d7$age < 65]  <- .1879632100462914
  d7$glucose_mgdl <- NA
  
  d8 <- d # 9th copy (70th percentile): intervention set to 8, outcome to missing
  d8$interv <- 8
  d8$intervention_b  [d8$age < 65]  <- .2076895833015442
  d8$glucose_mgdl <- NA
  
  d9 <- d # 10th copy (80th percentile): intervention set to 9, outcome to missing
  d9$interv <- 9
  d9$intervention_b [d9$age < 65]   <- .2252189218997955
  d9$glucose_mgdl <- NA
  
  d10 <- d # 11th copy (90th percentile): intervention set to 10, outcome to missing
  d10$interv <- 10
  d10$intervention_b [d10$age < 65]  <- .259517103433609
  d10$glucose_mgdl <- NA
  
  d11 <- d # 12th copy (Maximum): intervention set to 11, outcome to missing
  d11$interv <- 11
  d11$intervention_b [d11$age < 65]  <- .3624528
  d11$glucose_mgdl <- NA
  
  d.onesample <- rbind(d, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11) # combining datasets
  
  
  # linear regression model to estimate mean outcome conditional on treatment and confounders
  model <- glm(glucose_mgdl ~ intervention_b + age + intervention_b*age+  lat + racecat + incomecat + educcat + sex, data=d.onesample)
  
  summary(model)
  
  d.onesample$predicted_meanY <- predict(model, d.onesample)
  
  # estimate mean outcome in each of the groups
  # estimate difference in mean outcome in each of the 'exposure' groups (deciles + max NDVI) vs. 'unexposed' group (min NDVI)
  
  return(c(
    mean(d.onesample$predicted_meanY[d.onesample$interv == -1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 2]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 3]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 4]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 5]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 6]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 7]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 8]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 9]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 10]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 11]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 2])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 3])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 4])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 5])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 6])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 7])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 8])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 9])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 10])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1]),
    mean(d.onesample$predicted_meanY[d.onesample$interv == 11])-
      mean(d.onesample$predicted_meanY[d.onesample$interv == 1])
  ))
  
}

# bootstrap
results <- boot(data = pq, statistic = standardization, R = 1000)

# generating confidence intervals
se <- c(sd(results$t[, 1]),
        sd(results$t[, 2]),
        sd(results$t[, 3]),
        sd(results$t[, 4]),
        sd(results$t[, 5]),
        sd(results$t[, 6]),
        sd(results$t[, 7]),
        sd(results$t[, 8]),
        sd(results$t[, 9]),
        sd(results$t[, 10]),
        sd(results$t[, 11]),
        sd(results$t[, 12]),
        sd(results$t[, 13]),
        sd(results$t[, 14]),
        sd(results$t[, 15]),
        sd(results$t[, 16]),
        sd(results$t[, 17]),
        sd(results$t[, 18]),
        sd(results$t[, 19]),
        sd(results$t[, 20]),
        sd(results$t[, 21]),
        sd(results$t[, 22]))

mean <- results$t0
ll <- mean - qnorm(0.975) * se
ul <- mean + qnorm(0.975) * se

bootstrap_gluc_ageunder65 <-
  data.frame(cbind(
    c("Observed",
      "Minimum",
      "10th Percentile",
      "20th Percentile",
      "30th Percentile",
      "40th Percentile",
      "50th Percentile",
      "60th Percentile",
      "70th Percentile",
      "80th Percentile",
      "90th Percentile",
      "Maximum",
      "RD 10th - Min",
      "RD 20th - Min",
      "RD 30th - Min",
      "RD 40th - Min",
      "RD 50th - Min",
      "RD 60th - Min",
      "RD 70th - Min",
      "RD 80th - Min",
      "RD 90th - Min",
      "RD Max - Min"),
    mean,
    se,
    ll,
    ul
  ))
bootstrap_gluc_ageunder65