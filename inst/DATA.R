#===============================================================================
#
#  PROGRAM: DATA.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: NHANES data preparation for our motivating and supplementary
#           analyses on self-reported race, socioeconomic status, lead
#           lead exposure, and telomere length.
#
#  INPUT:   NHANES 1999-2000 (XXXX_A) and 2001-2002 (XXXX_B) data files:
#
#             - TELO_A, TELO_B - Telomere Length
#             - DEMO_A, DEMO_B - Demographics
#             - OCQ_A,  OCQ_B  - Occupation Questionnaire
#             - FSQ_A,  FSQ_B  - Food Security Questionnaire
#             - HOQ_A,  HOQ_B  - Housing Questionnaire
#             - HIQ_A,  HIQ_B  - Health Insurance Questionnaire
#             - LAB_A,  LAB_B  - CBC with 5-Part Differential
#             - LEAD_A, LEAD_B - Cadmium, Lead, Mercury, Cotinine
#
#  OUTPUT:  NHANES.rda
#
#           Final analytic dataset for Non-Hispanic White and Non-Hispanic
#           Black participants with measured telomere length in the 1999-2002
#           waves of the Health and Nutrition Examination Survey (NHANES).
#
#  UPDATED: 2024-05-25
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(nhanesA, here, tidyverse, update = FALSE)

#--- READ-IN DATA --------------------------------------------------------------

TELO_A <- nhanes("TELO_A")   #-- 1999-2000 - Telomere Length

TELO_B <- nhanes("TELO_B")   #-- 2001-2002 - Telomere Length

DEMO_A <- nhanes("DEMO")     #-- 1999-2000 - Demographics

DEMO_B <- nhanes("DEMO_B")   #-- 2001-2002 - Demographics

OCQ_A  <- nhanes("OCQ")      #-- 1999-2000 - Occupation Categories

OCQ_B  <- nhanes("OCQ_B")    #-- 2001-2002 - Occupation Categories

FSQ_A  <- nhanes("FSQ")      #-- 1999-2000 - Food Security

FSQ_B  <- nhanes("FSQ_B")    #-- 2001-2002 - Food Security

HOQ_A  <- nhanes("HOQ")      #-- 1999-2000 - Housing Characteristics

HOQ_B  <- nhanes("HOQ_B")    #-- 2001-2002 - Housing Characteristics

HIQ_A  <- nhanes("HIQ")      #-- 1999-2000 - Health Insurance Characteristics

HIQ_B  <- nhanes("HIQ_B")    #-- 2001-2002 - Health Insurance Characteristics

LAB_A  <- nhanes("LAB25")    #-- 1999-2000 - CBC with 5-Part Differential

LAB_B  <- nhanes("L25_B")    #-- 2001-2002 - CBC with 5-Part Differential

LEAD_A <- nhanes("LAB06")    #-- 1999-2000 - Cadmium, Lead, Mercury, Cotinine

LEAD_B <- nhanes("L06_B")    #-- 2001-2002 - Cadmium, Lead, Mercury, Cotinine

#=== DATA PRE-PROCESSING =======================================================

#--- TELOMERE LENGTH -----------------------------------------------------------

TELO <- bind_rows(TELO_A, TELO_B) |>

  mutate(lTELOMEAN = log(TELOMEAN))

#--- DEMOGRAPHICS --------------------------------------------------------------

DEMO <- bind_rows(DEMO_A, DEMO_B) |>

  mutate(

    #- Probability of Selection

    iWTMEC4YR = 1 / WTMEC4YR,

    #- Race

    RACE_2CAT = case_when(

      RIDRETH1 == "Non-Hispanic White" ~ 0,

      RIDRETH1 == "Non-Hispanic Black" ~ 1,

      TRUE ~ NA_real_),

    #- Age

    AGE = RIDAGEYR,

    #- Sex

    SEX = factor(RIAGENDR),

    #- Education

    EDUC_3CAT = case_when(

      DMDEDUC2 %in% c("Less Than 9th Grade",

        "9-11th Grade (Includes 12th grade with no diploma)",

        "High School Grad/GED or Equivalent") ~ 0,

      DMDEDUC2 == "Some College or AA degree" ~ 1,

      DMDEDUC2 == "College Graduate or above" ~ 2,

      TRUE ~ NA_real_) |>

      factor(0:2, c("High School or GED", "Some College", "College Graduate")),

    #- Marital Status

    MARTL_3CAT = case_when(

      DMDMARTL == "Never married" ~ 0,

      DMDMARTL %in% c("Widowed", "Divorced", "Separated") ~ 1,

      DMDMARTL %in% c("Married", "Living with partner") ~ 2,

      TRUE ~ NA_real_) |>

      factor(0:2, c("Never Married", "Widowed/Divorced/Separated",

        "Married/Living with Partner")),

    #- Household Size

    HHSIZE_5CAT = case_when(

      DMDHHSIZ == 1 ~ 0,

      DMDHHSIZ == 2 ~ 1,

      DMDHHSIZ == 3 ~ 2,

      DMDHHSIZ == 4 ~ 3,

      DMDHHSIZ %in% 5:7 ~ 4,

      TRUE ~ NA_real_) |>

      factor(0:4, c("1 Person", paste(c(2:4, "5+"), "People"))),

    #- Household Income

    HHINC_5CAT = case_when(

      INDHHINC %in% c("$     0 to $ 4,999", "$ 5,000 to $ 9,999",

        "$10,000 to $14,999", "$15,000 to $19,999", 13) ~ 0,

      INDHHINC %in% c("$20,000 to $24,999", "$25,000 to $34,999") ~ 1,

      INDHHINC %in% c("$35,000 to $44,999", "$45,000 to $54,999") ~ 2,

      INDHHINC %in% c("$55,000 to $64,999", "$65,000 to $74,999") ~ 3,

      INDHHINC == "$75,000 and Over" ~ 4,

      TRUE ~ NA_real_) |>

      factor(0:4, c("$0 - $20,000", "$20,000 - $35,000", "$35,000 - $55,000",

        "$55,000 - $75,000", "$75,000+")),

    #- Poverty Income Ratio

    PIR_3CAT = cut(INDFMPIR, right = FALSE, breaks = c(-Inf, 1.3, 3.5, Inf),

      labels = c("< 1.3", "1.3 <= PIR < 3.5", ">= 3.5")) |>

      factor())

#--- OCCUPATION CATEGORIES -----------------------------------------------------

OCQ <- bind_rows(OCQ_A, OCQ_B) |>

  mutate(

    #- Align 'Longest' and 'Recent' Work

    OCD390 = if_else(OCQ390G == "Same as current occupation", OCD240, OCD390),

    #- Employment Status

    EMPSTAT_4CAT = case_when(

      (OCQ150 == "Working at a job or business," &

        OCQ180 >= 35 & OCQ180 < 77777) |

          (OCQ150 == "With a job or business but not at work," &

             OCQ180 >= 35 & OCQ180 < 77777) ~ 0,

      (OCQ150 == "Working at a job or business," &  OCQ180 > 0 & OCQ180 < 35) |

        (OCQ150 == "With a job or business but not at work," &

           OCQ180 > 0 & OCQ180 < 35) ~ 1,

      (OCQ150 %in% c("Looking for work, or",

        "Not working at a job or business?") & OCQ380 == "Retired") ~ 2,

      (OCQ150 %in% c("Looking for work, or",

        "Not working at a job or business?")) ~ 3,

      TRUE ~ NA_real_) |>

      factor(0:3, c("Full-Time", "Part-Time", "Retired", "Not Working")),

    #- Occupation Category

    OCC_5CAT = case_when(OCQ390G == "Never worked" ~ 0,

      OCD390 %in% c(11, 17:21, 23:24, 26:27, 32:40) ~ 1,

      OCD390 %in% c(28:31, 41) ~ 2,

      OCD390 %in% c(8, 10, 12:16, 22) ~ 3,

      OCD390 %in% c(1:7, 9, 25) ~ 4,

      OCD390 %in% c(98) ~ NA_real_,

      TRUE ~ NA_real_) |>

      factor(0:4, c("No Work", "Low Blue Collar", "High Blue Collar",

        "Low White Collar", "High White Collar")))

#--- FOOD SECURITY -------------------------------------------------------------

FSQ <- bind_rows(FSQ_A, FSQ_B) |>

  mutate(

    #- WIC Utilization

    WIC_2CAT = case_when(

      is.na(FSD160) ~ 0,

      FSD160 == "No" ~ 0,

      FSD160 == "Yes" ~ 1,

      TRUE ~ NA_real_) |>

      factor(0:1, c("No WIC", "Received WIC")),

    #- Food Security

    FDSEC_3CAT = case_when(HHFDSEC == 1 ~ 0, HHFDSEC == 2 ~ 1,

      HHFDSEC %in% 3:4 ~ 2, TRUE ~ NA_real_) |>

      factor(0:2, c("Food Secure", "Marginally Food Secure", "Food Insecure")))

#--- HOUSING CHARACTERISTICS ---------------------------------------------------

HOQ <- bind_rows(HOQ_A, HOQ_B) |>

  mutate(

    #- Home Ownership

    OWNHOME_2CAT = case_when(

      HOQ065 %in% c("Rented", "Other arrangement") ~ 0,

      HOQ065 == "Owned or being bought" ~ 1,

      TRUE ~ NA_real_) |>

      factor(0:1, c("Does Not Own Home", "Owns Home")),

    #- Housing Type

    HOD_4CAT = case_when(

      HOD010 == "A one family house detached from any other house," ~ 0,

      HOD010 == "A one family house attached to one or more houses," ~ 1,

      HOD010 == "An apartment," ~ 2,

      HOD010 %in% c("A mobile home or trailer,", "Something else,",

        "Dormitory? ") ~ 3,

      TRUE ~ NA_real_) |>

      factor(0:3, c("Family Home Detached", "Family Home Attached",

        "Apartment", "Other")))

#--- HEALTH INSURANCE CHARACTERISTICS ------------------------------------------

HIQ <- bind_rows(HIQ_A, HIQ_B) |>

  mutate(

    #- Health Insurance

    HIQ_2CAT = case_when(

      HID010 == "No" ~ 0,

      HID010 == "Yes" ~ 1,

      TRUE ~ NA_real_) |>

      factor(0:1, c("Not Insured", "Insured")))

#--- CBC WITH 5-PART DIFFERENTIAL ----------------------------------------------

LAB <- bind_rows(LAB_A, LAB_B)

#--- CADMIUM, LEAD, MERCURY, COTININE ------------------------------------------

LEAD <- bind_rows(LEAD_A, LEAD_B) |>

  mutate(

    LBXBPB_LOD = if_else(LBXBPB < 0.3, 0.3 / sqrt(2), LBXBPB)) #-- LOD / sqrt(2)

#=== FINAL ANALYTIC DATASET ====================================================

NHANES <- TELO |>

  #-- MERGE DATASETS

  left_join(DEMO, "SEQN") |>

  left_join(OCQ,  "SEQN") |>

  left_join(FSQ,  "SEQN") |>

  left_join(HOQ,  "SEQN") |>

  left_join(HIQ,  "SEQN") |>

  left_join(LAB,  "SEQN") |>

  left_join(LEAD, "SEQN") |>

  #-- FILTER ANALYTIC SAMPLE

  filter(!is.na(RACE_2CAT), !is.na(TELOMEAN)) |>

  #-- MAKE REFUSED/UNKNOWN/MISSING AN EXPLICIT FACTOR LEVEL

  mutate(across(where(is.factor),

    function(x) fct_na_value_to_level(x, "Refused/Unknown"))) |>

  #-- KEEP NECESSARY VARIABLES

  select(SEQN, iWTMEC4YR, WTMEC4YR, SDMVPSU, SDMVSTRA, TELOMEAN, lTELOMEAN,

    RACE_2CAT, AGE, SEX, EDUC_3CAT, MARTL_3CAT, PIR_3CAT, EMPSTAT_4CAT,

    OCC_5CAT, HHSIZE_5CAT, HHINC_5CAT, HIQ_2CAT, WIC_2CAT, FDSEC_3CAT, HOD_4CAT,

    OWNHOME_2CAT, LBXWBCSI, LBXLYPCT, LBXMOPCT, LBXNEPCT, LBXEOPCT, LBXBAPCT,

    LBXBPB_LOD)

save(NHANES, file = here("data", "NHANES.rda"))

#=== END =======================================================================
