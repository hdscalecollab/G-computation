#G Computation Simulation Code for Effect Modification by Sex
#Manuscript "Simulating the Impact of Green Space Exposure on Cardiometabolic Biomarkers in a
#Diverse Population Living in San Diego, California: A G-Computation Application"

#Appending dataframes together into final dataset
# Example provided is for Blood Glucose Level Outcome
# Can be substituted for any of the other 6 biomarkers

#Version 3.6.2024

# Contact: Anais Teyton, ateyton@ucsd.edu

#Libraries
library(openxlsx)

# Create a list of all data frames
df_list <- list(bootstrap_gluc, bootstrap_gluc_male, bootstrap_gluc_female, bootstrap_gluc_hisp, bootstrap_gluc_nothisp, bootstrap_gluc_inc1, bootstrap_gluc_inc2, bootstrap_gluc_inc3, bootstrap_gluc_ageover65, bootstrap_gluc_ageunder65)

# Add a "name" column to each dataframe
bootstrap_gluc <- bootstrap_gluc %>% mutate(name = "bootstrap_gluc")
bootstrap_gluc_male <- bootstrap_gluc_male %>% mutate(name = "bootstrap_gluc_male")
bootstrap_gluc_female <- bootstrap_gluc_female %>% mutate(name = "bootstrap_gluc_female")
bootstrap_gluc_hisp <- bootstrap_gluc_hisp %>% mutate(name = "bootstrap_gluc_hisp")
bootstrap_gluc_nothisp <- bootstrap_gluc_nothisp %>% mutate(name = "bootstrap_gluc_nothisp")
bootstrap_gluc_inc1 <- bootstrap_gluc_inc1 %>% mutate(name = "bootstrap_gluc_inc1")
bootstrap_gluc_inc2 <- bootstrap_gluc_inc2 %>% mutate(name = "bootstrap_gluc_inc2")
bootstrap_gluc_inc3 <- bootstrap_gluc_inc3 %>% mutate(name = "bootstrap_gluc_inc3")
bootstrap_gluc_ageover65 <- bootstrap_gluc_ageover65 %>% mutate(name = "bootstrap_gluc_ageover65")
bootstrap_gluc_ageunder65 <- bootstrap_gluc_ageunder65 %>% mutate(name = "bootstrap_gluc_ageunder65")

gluc_FinalResults <- bind_rows(
  bootstrap_gluc,
  bootstrap_gluc_male,
  bootstrap_gluc_female,
  bootstrap_gluc_hisp,
  bootstrap_gluc_nothisp,
  bootstrap_gluc_inc1,
  bootstrap_gluc_inc2,
  bootstrap_gluc_inc3,
  bootstrap_gluc_ageover65,
  bootstrap_gluc_ageunder65
)

gluc_FinalResults$mean <- round(as.numeric(gluc_FinalResults$mean), 2)
gluc_FinalResults$ll <- round(as.numeric(gluc_FinalResults$ll), 2)
gluc_FinalResults$ul <- round(as.numeric(gluc_FinalResults$ul), 2)

write.xlsx(gluc_FinalResults, "SAVE FINAL RESULTS HERE")

