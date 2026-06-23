library(lubridate)
library(tidyverse)
library(hms)
# data cleaning
df <- read_csv("MMF-Final-3.csv")

# Finding missed observations and fixing errors
df <- df %>% 
  add_row(
    Infection_number = 34,
    Date = "18/06/2026",
    `Person Infected` = "Abalo", 
    Time = as_hms("18:10:00"), 
    `Infected by` = "Thenuja"
  )
df <- df %>%
  mutate(`Infected by` = case_when(
    `Person Infected` == "Kimberley" ~ "Started",  # Change Time to 13:55 if Person is Disebo
    TRUE               ~ `Infected by`     # Keep everything else exactly the same
  ))
df <- df %>%
  mutate(`Infected by` = case_when(
    `Person Infected` == "Morgan" ~ "Started",  # Change Time to 13:55 if Person is Disebo
    TRUE               ~ `Infected by`     # Keep everything else exactly the same
  ))
df <- df %>%
  mutate(`Infected by` = case_when(
    `Person Infected` == "Mandie" ~ "Started",  # Change Time to 13:55 if Person is Disebo
    TRUE               ~ `Infected by`     # Keep everything else exactly the same
  ))



df <- df %>%
  mutate(Date = case_when(
    `Person Infected` == "Gebrekiros" ~ "18/06/2026",  # Change Time to 13:55 if Person is Disebo
    TRUE               ~ Date     # Keep everything else exactly the same
  ))
df <- df %>%
  mutate(Time = case_when(
    `Person Infected` == "Disebo" ~ as_hms("13:58:00"),  # Change Time to 13:55 if Person is Disebo
    TRUE               ~ Time     # Keep everything else exactly the same
  ))
# Adding Location, 0 for AIMS and 1 for Empire
df<- df %>% mutate("Location" = c(1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0))

# Adding doctor visits to a specific column, 0 for asymptomatic, 1 for symptomatic
df<- df %>% mutate("Doctor" = c(1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

# Adding duplicate Infections for the network

df <- df %>% 
  add_row(
    Infection_number = 35,
    Date = "16/06/2026",
    `Person Infected` = "Kimberley", 
    Time = as_hms("18:05:00"), 
    `Infected by` = "Cebile", 
    Location = 1, 
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 36,
    Date = "17/06/2026",
    `Person Infected` = "Pendo", 
    Time = as_hms("10:15:00"), 
    `Infected by` = "Cebile", 
    Location = 1, 
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 36,
    Date = "18/06/2026",
    `Person Infected` = "Ethel", 
    Time = as_hms("18:22:00"), 
    `Infected by` = "Mario", 
    Location = 0, 
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 37,
    Date = "18/06/2026",
    `Person Infected` = "Cebile", 
    Time = as_hms("13:52:00"), 
    `Infected by` = "Boikanyo", 
    Location = 1, 
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 38,
    Date = "19/06/2026",
    `Person Infected` = "Gebrekiros", 
    Time = as_hms("15:10:00"), 
    `Infected by` = "Thulisile", 
    Location = 1,
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 39,
    Date = "18/06/2026",
    `Person Infected` = "Ethel", 
    Time = as_hms("00:47:00"), 
    `Infected by` = "Ennie", 
    Location = 0,
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 40,
    Date = "19/06/2026",
    `Person Infected` = "Kisaka", 
    Time = as_hms("12:22:00"), 
    `Infected by` = "Mija", 
    Location = 1,
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 41,
    Date = "19/06/2026",
    `Person Infected` = "Haron", 
    Time = as_hms("15:04:00"), 
    `Infected by` = "Mija", 
    Location = 1,
    Doctor = 0
  )
df <- df %>% 
  add_row(
    Infection_number = 42,
    Date = "18/06/2026",
    `Person Infected` = "Vix", 
    Time = as_hms("15:00:00"), 
    `Infected by` = "Thenuja", 
    Location = 0, 
    Doctor = 0  )


# 2. Save the objects as CSV files, run only the corrective code you want for the specific data set you would like
write_csv(df, "MMF-Final+Duplicates2.csv")

mean(table(MMF_Final_Locations_Doctors$`Infected by`))
