#' Race, SES, and Telomere Length Data
#'
#' National Health and Nutrition Examination Survey (NHANES) data on race,
#' socioeconomic status, and leukocyte telomere length from the 1999-2000
#' and 2001-2002 survey waves.
#'
#' @name NHANES
#'
#' @docType data
#'
#' @usage data(NHANES)
#'
#' @details
#'
#' Our initial sample consisted of 7,839 participants in the 1999-2002 NHANES
#' waves with laboratory measures recorded, including telomere length,
#' \code{lTELOMEAN}, which was assayed via quantitative polymerase chain
#' reaction (PCR; see Cawthorn, 2002). Our primary endpoint is the
#' log-transformed mean ratio of an individual's telomere length to a standard
#' reference DNA sample across all leukocyte cell types (mean T/S),
#' \code{TELOMEAN}. We focus on the 1999-2002 NHANES waves, as they featured
#' 4-year adjusted survey weights, \code{WTMEC4YR}, designed for aggregating
#' data across cohorts. Among the initial 7,839 participants, 5,308 (67.7\%)
#' self-identified as either non-Hispanic White or non-Hispanic Black.
#' Excluding those participants without our outcome of interest, our final
#' analytic sample contained 5,298 Non-Hispanic White or Non-Hispanic Black
#' identifying participants with measured telomere length. Race,
#' \code{RACE_2CAT}, is our variable of interest. We further included study
#' participant age, sex, and blood cell composition to account for known
#' differences in these factors, as well as twelve indicators of SES. Ten of
#' these, namely marital status, education level, household income, insurance
#' status, Special Supplemental Nutrition Program for Women, Infants, and
#' Children (WIC) usage, household size, home ownership, home type, food
#' security status, and an individual’s poverty income ratio (PIR), were
#' extracted directly from the NHANES demographic and occupation questionnaires.
#' Occupation category was constructed by mapping occupation group codes in the
#' NHANES 1999-2002 occupation questionnaire to the national statistics
#' socioeconomic job classifications, and employment status was derived from
#' three occupational measures: type of work done last week, hours worked last
#' week at all jobs, and main reason for not working last week (see Rehkopf et
#' al., 2008, Rose et al., 2005).
#'
#' @format
#' A dataset with 5,298 observations (rows) of 29 variables (columns):
#' \describe{
#'    \item{SEQN}{Numeric: Respondent Sequence Number}
#'    \item{iWTMEC4YR}{Numeric: 1/WTMEC4YR (Full Sample 4 Year Probability of
#'                     Selection)}
#'    \item{WTMEC4YR}{Numeric: Full Sample 4 Year Interview Weight}
#'    \item{SDMVPSU}{Numeric: Masked Variance Pseudo-PSU}
#'    \item{SDMVSTRA}{Numeric: Masked Variance Pseudo-Stratum}
#'    \item{TELOMEAN}{Numeric: Mean T/S Ratio (See Details)}
#'    \item{lTELOMEAN}{Numeric: log(TELOMEAN)}
#'    \item{RACE_2CAT}{Numeric: 0 = Non-Hispanic White, 1 = Non-Hispanic Black
#'                     (0/1 Coded for Current Functionality)}
#'    \item{AGE}{Numeric: Age at Screening (Years)}
#'    \item{SEX}{Factor w/ 2 Levels: Self-Reported Sex - Male, Female}
#'    \item{EDUC_3CAT}{Factor w/ 3 Levels: Education - High School or GED,
#'                     Some College, College Graduate}
#'    \item{MARTL_3CAT}{Factor w/ 3 Levels: Marital Status - Never Married,
#'                      Widowed/Divorced/Separated, Married/Living with Partner}
#'    \item{HHSIZE_3CAT}{Factor w/ 5 Levels: Household Size - 1 Person,
#'                       2 People, 3 People, 4 People, 5+ People}
#'    \item{HHINC_5CAT}{Factor w/ 5 Levels: Annual Household Income -
#'                      $0 - $20,000, $20,000 - $35,000, $35,000 - $55,000,
#'                      $55,000 - $75,000, $75,000+}
#'    \item{PIR}{Factor w/ 3 Levels: Family Poverty-Income Ratio Category -
#'               < 1.3, 1.3 <= PIR < 3.5, >= 3.5}
#'    \item{EMPSTAT_4CAT}{Factor w/ 4 Levels: Employment Status - Full-Time,
#'                        Part-Time, Retired, Not Working}
#'    \item{OCC_5CAT}{Factor w/ 5 Levels: Occupation Category - No Work,
#'                    Low Blue Collar, High Blue Collar, Low White Collar,
#'                    High White Collar}
#'    \item{WIC_2CAT}{Factor w/ 2 Levels: WIC Utilization - No WIC,
#'                    Received WIC}
#'    \item{FDSEC_3CAT}{Factor w/ 3 Levels: Food Security Status - Food Secure,
#'                      Marginally Food Secure, Food Insecure}
#'    \item{HOD_4CAT}{Factor w/ 4 Levels: Home Type - Family Home Detached,
#'                    Family Home Attached, Apartment, Other}
#'    \item{OWNHOME_2CAT}{Factor w/ 2 Levels: Home Ownership -
#'                        Does Not Own Home, Owns Home}
#'    \item{HIQ_2CAT}{Factor w/ 2 Levels: Insurance Status - Not Insured,
#'                    Insured}
#'    \item{LBXWBCSI}{Numeric: White Blood Cell Count (SI)}
#'    \item{LBXLYPCT}{Numeric: Lymphocyte Percent (\%)}
#'    \item{LBXMOPCT}{Numeric: Monocyte Percent (\%)}
#'    \item{LBXNEPCT}{Numeric: Segmented Neutrophils Percent (\%)}
#'    \item{LBXEOPCT}{Numeric: Eosinophils Percent (\%)}
#'    \item{LBXBAPCT}{Numeric: Basophils Percent (\%)}
#'    \item{LBXBPB_LOD}{Numeric: Blood Lead Concentration (ug/dL; LOD = 0.3
#'                      ug/dL; Imputed by LOD / sqrt(2))}
#' }
#'
#' @keywords datasets
#'
#' @references
#'
#' Richard M Cawthon. Telomere measurement by quantitative pcr. Nucleic acids
#' research, 30(10):e47–e47, 2002.
#'
#' David H Rehkopf, Lisa F Berkman, Brent Coull, and Nancy Krieger. The
#' non-linear risk of mortality by income level in a healthy population: Us
#' national health and nutrition examination survey mortality follow-up cohort,
#' 1988–2001. BMC Public Health, 8(1):1–11, 2008.
#'
#' David Rose, David J Pevalin, and Karen O’Reilly. The National Statistics
#' Socio-economic Classification: origins, development and use. Palgrave
#' Macmillan, 2005.
#'
#' @source <https://www.cdc.gov/nchs/nhanes/index.htm>
#'
#' @examples
#' data(NHANES)
NULL
