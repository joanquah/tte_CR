# Clean data for analysis
library(tidyr)
library(knitr)
library(dplyr)
library(ggplot2)
library(car)  # For deltaMethod for CI
library(skimr)
library(lubridate)
library(data.table)
library(survival)
library(survminer)
library(gtsummary) #For baseline table
library(nnet) # For multinomial log regression
library(tableone)
library(stringr) # For string manipulation
library(cobalt) # SMD
library(mice)

# From clean data
infection_types_index <- inf_add
baseline_outcomes_index <- df_baseline_all
ast_all <- all_ast_merged 
ast_all_index <- all_ast_merged_index
anti_treat
anti_treat_index <- anti_treat_index_new
all_vap_bsi
vap_bsi_index

length(unique(ast_all$recordid)) #9369
length(unique(baseline_outcomes_index$recordid)) #9678
length(unique(f02$recordid)) #9171

f07a <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07aNew_DATA_2025-04-07_1345.csv", header = TRUE, sep = ",")
f07e <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07eNew_DATA_2025-04-07_1346.csv", header = TRUE, sep = ",")
f07m <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07mNew_DATA_2025-04-07_1347.csv", header = TRUE, sep = ",")
f02 <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F02New_DATA_2025-04-07_1344.csv", header = TRUE, sep = ",")

# Run correction code

checkast <-ast_all%>% select(recordid,spec_date, org_names_all, pathogen_group,ris_Carbapenems, Meropenem)
length(unique(ast_all$recordid)) #9369

CR <- checkast %>%
  mutate(recordid = ifelse(grepl("_BSI|_VAP", recordid), sub("_BSI|_VAP", "", recordid), recordid)) %>%
  mutate(susceptibility = case_when(
    ris_Carbapenems %in% c(2, 3) ~ "Susceptible",
    ris_Carbapenems == 1 ~ "Resistant",
    ris_Carbapenems %in% c(4, 5) ~ "Unknown"
  )) %>%
  mutate(organism = case_when(
    # Assign "CRAB" only if org_names_all matches AND susceptibility is NOT "Susceptible"
    org_names_all %in% c("Acinetobacter", "acinetobacter spp", "acinetobacter baumini") & 
      susceptibility != "Susceptible" ~ "CRAB",
    
    # Assign "CRE" only if org_names_all matches AND susceptibility is NOT "Susceptible"
    org_names_all %in% c("K. pneumoniae", "E. coli", 
                         "Serratia", "Proteus", "Providencia", "Morganella", 
                         "Enterobacter", "Enterobacter aerogenes", "Enterobacter cloacae", 
                         "Klebsiella", "Klebsiella acidogenes", "Klebsiella oxytocin",
                         "Citrobacter") & 
      susceptibility != "Susceptible" ~ "CRE",
    
    # Assign "CRE" only if org_names_all="Pseudomonas" AND susceptibility is NOT "Susceptible"
    org_names_all %in% c("Pseudomonas", "pseudomonas spp", "pseudomonas aeruginosa") & 
      susceptibility != "Susceptible" ~ "CRPAE",
    
    # Default to NA if no match or susceptibility is "Susceptible"
    TRUE ~ NA_character_
  )) %>%
  filter(organism == "CRE"|organism =="CRAB"|organism =="CRPAE") %>%
  #filter(!(susceptibility == "Susceptible"))%>%
  filter(susceptibility =="Resistant")

length(unique(CR$recordid)) # 3610  # 3780  #9369  #3859 CNS  # 3078 CR

outcome <- baseline_outcomes_index %>% 
  select(recordid, date_enrolment,age_new, age_group, country_income, country, country_ab, sex, hpd_adm_date, hpd_hosp_date, inf_onset,
         comorbidities_Chalson, sofa_score,adm_ward_types_new, infection_types, adm_ward_types_new, UTI_source,
         ho_discharge_date, ho_dischargestatus, d28_date,d28_status,d28_death_date,first28_death, 
         mortality,mortality_date, score_respiration,score_coagulation,score_liver,score_CVS,score_CNV, score_renal,
  ) %>%
  mutate(infection_types = ifelse(infection_types %in% c("Hospital-acquired BSI", "Healthcare-associated BSI"), "BSI", infection_types))%>%
  filter(age_group==1)%>% 
  mutate(country_income = ifelse(country_ab == "CN", "Upper middle income", country_income))%>%
  mutate(country_income = ifelse(country_ab == "IN", "Lower middle income", country_income)) %>%
  mutate(adm_date = pmin(hpd_adm_date, hpd_hosp_date, na.rm = TRUE)) 

length(unique(outcome$recordid)) #7882 adult with baseline and outcome data

CR2 <- CR %>% inner_join(outcome, by = "recordid") %>% distinct() %>%
  mutate(
    inf_onset = as.Date(inf_onset),
    spec_date = as.Date(spec_date),
    #ast_date = as.Date(ast_date),
    mortality_date = as.Date(mortality_date)) %>%  #2451 adult inpatients who has CR
  filter(!(mortality_date == inf_onset & !is.na(mortality_date))) %>% 
  filter(!(mortality_date == (inf_onset + 1))| is.na(mortality_date))%>%   #2283 adult pts who had CR and did not die within 1 day of infection onset
  select(recordid, inf_onset, spec_date, ho_discharge_date, ho_dischargestatus,d28_date,d28_status,mortality, mortality_date, 
         everything())%>%  
  filter(!(is.na(ho_discharge_date) & is.na(ho_dischargestatus) & is.na(mortality) & is.na(mortality_date) & is.na(d28_date) & is.na(d28_status)))

length(unique(CR2$recordid)) #2667 #2843  #7110  #2260 adult inpatient who did not die within 1 day and has baseline and outcme data

# Data correction to change 
# remove PK-002-A0-0357 from the dataset
CR2 <- CR2[CR2$recordid != "PK-002-A0-0357",]
length(unique(CR2$recordid)) #2259 adult inpatient who did not die within 1 day and has baseline and outcme data

abx_all <- anti_treat %>% 
  filter(!is.na(anti_start)) %>%
  filter(!is.na(anti_names)) %>%
  mutate(anti_new = case_when(
    anti_names %in% c("Imipenem", "Meropenem", "Doripenem") ~ "Carbapenem",
    anti_names %in% c("Colistin", "Polymyxin B") ~ "Polymyxin",
    anti_names %in% c("Ampicillin/sulbactam", "Cefoperazone/sulbactam") ~ "Sulbactam",
    anti_names == "Ceftazidime/avibactam" ~ "Ceftazidime/avibactam",
    anti_names == "Tigecycline"  ~ "Tigecycline",
    TRUE ~ anti_names
  )) 

abx <- abx_all %>%
  group_by(recordid) %>%
  filter(any(anti_new %in% c("Polymyxin", "Sulbactam", "Ceftazidime/avibactam"))) %>%
  ungroup()

CR_abx2 <- abx %>% inner_join(CR2, by = "recordid") %>% distinct()
length(unique(CR_abx2$recordid)) #1439   #1538  #2251  # 1429 received studied drug (polymyxin or sulbactam or ceftazidime avibactam)

excluded_names <- c(
  "Amphotericin B", "Other_Fluconazole", "Caspofungin","Other_Voriconazole", "Anidulafungin", 'Voriconazole',"Fluconazole","Posaconazole",
  "Clindamycin","Teicoplanin", "Chloramphenicol", "Linezolid", "Daptomycin","Metronidazole", "Nitrofurantoin", "Fusidic acid", "Vancomycin", "Rifampicin", 
  "Azithromycin", "Clarithromycin", "Erythromycin", NA, 
  "Amoxicillin-clavulanate", "Cloxacillin", "Flucloxacillin","Piperacillin-tazobactam","Piperacillin/tazobactam", 
  "Ampicillin","Amoxicillin","Benzylpenicillin","Oxacillin", "Penicillin", 
  "Ceftazidime", "Ceftriaxone", "Cefuroxime",  "Cefepime","Cefoxitin", "Cefotaxime","Cefazolin", 
  "Cefixime", "Cefpirome", "Ceftazidime/clavulanic acid","Cephalexin",
  "Co-trimoxazole", "Trimethoprim","Trimethoprim/sulfamethoxazole",
  "Ciprofloxacin", "Levofloxacin" , "Ofloxacin", "Sitafloxacin","Moxifloxacin",
  "Doxycycline" , "Tetracycline", "Minocycline", "Gentamicin", "Amikacin" , "Streptomycin",
  "Ceftolozane/tazobactum", "Ceftolozane/tazobactam","Ertapenem")


CR_abx2 <- CR_abx2[!is.na(CR_abx2$anti_names) & !CR_abx2$anti_names %in% excluded_names, ] 

length(unique(CR_abx2$recordid)) #1439   #1538 #2251  #1556 had qualifying antibiotics # 1429 received studied drug (polymyxin or sulbactam or ceftazidime avibactam)

 # # Define the antibiotics of interest
# abx_to_exclude <- c("Ceftazidime/avibactam", "Colistin", "Polymyxin B")
# #,"Cefoperazone/sulbactam")
# 
# # Identify recordid to exclude
# recordid_to_exclude <- CR_abx2 %>%
#   mutate(delay = as.numeric(difftime(anti_start, inf_onset, units = "days"))) %>%
#   filter(anti_names %in% abx_to_exclude & delay <= -5) %>%
#   distinct(recordid)

CR_names <- CR_abx2 %>% 
  mutate(spec_day = as.integer(difftime(spec_date, inf_onset, units = "days")))%>%
  relocate(spec_day, .after = spec_date) %>% 
  #filter(spec_day<6)
  filter(spec_day<6) %>%    #1270
  mutate(onset4 = inf_onset + 4)%>% relocate(onset4, .after = spec_date) %>%
  mutate(onset5 = inf_onset + 5) %>% relocate(onset5, .after = onset4) %>%
  select(recordid, organism, adm_date, inf_onset, spec_date, onset4, onset5, anti_names, anti_new, anti_start, anti_end, susceptibility)%>%
  mutate (delay = as.numeric(difftime(anti_start, inf_onset, units = "days")))%>%
  filter(delay >= -5 | anti_new == "Carbapenem") %>% #1241
  # filter(!recordid %in% recordid_to_exclude$recordid) %>%
  filter(delay <= 5) %>%  #1195
  mutate(dur = as.integer(difftime(anti_end, anti_start, units = "days"))) %>%
  mutate (dur = dur+1) %>%
  # filter(dur>2) %>%
  filter(onset4 >= anti_start & onset4 < anti_end) %>%
  # filter(onset5 >= anti_start & onset5 < anti_end) %>%
  mutate(anti_names = case_when(
    recordid == "BN-001-A0-0028" & anti_start == ymd("2024-01-27") & anti_end == ymd("2024-02-09") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "MG-001-A0-0001" & anti_start == ymd("2023-03-31") & anti_end == ymd("2023-04-10") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "MG-001-A0-0039" & anti_start == ymd("2024-08-08") & anti_end == ymd("2024-08-22") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "MY-004-A0-0001" & anti_start == ymd("2023-03-29") & anti_end == ymd("2023-04-03") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "NP-003-A0-0038" & anti_start == ymd("2024-06-25") & anti_end == ymd("2024-07-02") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "NP-003-A0-0051" & anti_start == ymd("2024-10-24") & anti_end == ymd("2024-11-03") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "PH-001-A0-0051" & anti_start == ymd("2023-09-25") & anti_end == ymd("2023-10-01") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "TH-001-A0-0110" & anti_start == ymd("2023-01-10") & anti_end == ymd("2023-01-11") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "TH-001-A0-0153" & anti_start == ymd("2023-02-10") & anti_end == ymd("2023-02-18") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "TH-001-A0-0391" & anti_start == ymd("2023-07-27") & anti_end == ymd("2023-07-31") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "TH-005-A0-0001" & anti_start == ymd("2024-03-05") & anti_end == ymd("2024-03-22") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0080" & anti_start == ymd("2024-03-18") & anti_end == ymd("2024-03-26") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0080" & anti_start == ymd("2024-04-11") & anti_end == ymd("2024-04-16") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0099" & anti_start == ymd("2024-05-17") & anti_end == ymd("2024-05-25") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0120" & anti_start == ymd("2024-07-26") & anti_end == ymd("2024-08-02") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0123" & anti_start == ymd("2024-08-07") & anti_end == ymd("2024-08-28") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0128" & anti_start == ymd("2024-08-22") & anti_end == ymd("2024-08-31") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0149" & anti_start == ymd("2024-12-02") & anti_end == ymd("2024-12-04") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0159" & anti_start == ymd("2024-12-15") & anti_end == ymd("2024-12-24") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "NP-003-A0-0038" & anti_start == ymd("2024-06-25") & anti_end == ymd("2024-07-02") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "NP-003-A0-0045" & anti_start == ymd("2024-09-16") & anti_end == ymd("2024-09-25") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0120" & anti_start == ymd("2024-07-26") & anti_end == ymd("2024-08-02") & anti_names == "Colistin" ~ "Neb Col",
    recordid == "VN-005-A0-0155" & anti_start == ymd("2024-12-04") & anti_end == ymd("2024-12-12") & anti_names == "Colistin" ~ "Neb Col",
    TRUE ~ anti_names
  ))%>%
  group_by(recordid) %>%
  mutate(
    anti_onset4   = paste(unique(anti_names[onset4 >= anti_start & onset4 <= anti_end]), collapse = " + "),
    anti_onset5   = paste(unique(anti_names[onset5 >= anti_start & onset5 <= anti_end]), collapse = " + ")
  ) %>%
  mutate(
    anti_onset4_s   = paste(unique(anti_new[onset4 >= anti_start & onset4 <= anti_end]), collapse = " + "),
  ) %>%
  ungroup()%>% 
  distinct()  #929
   #903 qualified ITT  #849 initiated studied drug within 5 days from infection onset

length(unique(CR_names$recordid)) #1290. For onset 4 988 for ITT, 976 for PP; For onset5 94 for ITT, 941 for PP   # 849  #957 after fixing 1 line
#CR_names <- CR_names %>% select (recordid, organism, susceptibility, anti_ast, anti_onset4,anti_onset5) %>% distinct

CR_names_4 <- CR_names %>%
  mutate(arm = case_when(
    grepl("Ceftazidime/avibactam", anti_onset4) ~ "Ceftazidime/avibactam based",
    anti_onset4 %in% c("Polymyxin B", "Colistin") ~ "Polymyxin monotherapy",
    anti_onset4 %in% c("Cefoperazone/sulbactam", "Ampicillin/sulbactam") ~ "Sulbactam based",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset4) ~ "Polymyxin Sulbactam",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset4) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset4) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset4) ~ "Sulbactam based",
    TRUE ~ NA_character_
  )) %>%
  mutate(arm2 = case_when(
    anti_onset4 %in% c("Polymyxin B", "Colistin") ~ "Polymyxin monotherapy",
    anti_onset4 %in% c("Cefoperazone/sulbactam", "Ampicillin/sulbactam") ~ "Sulbactam monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset4) ~ "Polymyxin Sulbactam",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset4) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset4) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset4) ~ "Sulbactam combination",
    TRUE ~ NA_character_
  )) %>%
  filter(grepl("Colistin|Polymyxin B|Ampicillin/sulbactam|Cefoperazone/sulbactam|Ceftazidime/avibactam", anti_onset4))%>%
  mutate(arm = ifelse(anti_onset4 == "Polymyxin B + Colistin", "Polymyxin monotherapy", arm)) %>%
  mutate(arm = ifelse(anti_onset4 == "Cefoperazone/sulbactam + Ampicillin/sulbactam", "Sulbactam based", arm)) %>%
  mutate(arm2 = ifelse(anti_onset4 == "Cefoperazone/sulbactam + Ampicillin/sulbactam", "Sulbactam monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset4 == "Cefoperazone/sulbactam + Neb Col", "Sulbactam based", arm)) %>%
  mutate(arm2 = ifelse(anti_onset4 == "Cefoperazone/sulbactam + Neb Col", "Sulbactam monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset4 == "Colistin + Neb Col", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset4 == "Colistin + Neb Col", "Polymyxin monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset4 == "Polymyxin B + Neb Col", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset4 == "Polymyxin B + Neb Col", "Polymyxin monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset4 == "Neb Col + Polymyxin B", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset4 == "Neb Col + Polymyxin B", "Polymyxin monotherapy", arm2)) %>%
  
  mutate(arm3 = case_when(
    grepl("Ceftazidime/avibactam", anti_onset4) ~ "Ceftazidime/avibactam based",
    anti_onset4 %in% c("Polymyxin B", "Colistin","Polymyxin B + Colistin", 
                       "Polymyxin B + Ampicillin/sulbactam", "Ampicillin/sulbactam + Polymyxin B",
                       "Ampicillin/sulbactam + Colistin","Colistin + Ampicillin/sulbactam",
                       "Colistin + Cefoperazone/sulbactam","Cefoperazone/sulbactam + Colistin",
                       "Colistin + Ampicillin/sulbactam + Cefoperazone/sulbactam") ~ "Polymyxin monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Meropenem|Imipenem|Tigecycline|Fosfomycin", anti_onset4) ~ "Polymyxin combination",
    TRUE ~ NA_character_
  )) %>%
  mutate(anti_onset4_s = str_split(anti_onset4_s, " \\+ ") %>% 
           lapply(sort) %>% 
           sapply(paste, collapse = " + ")) %>%
  mutate(number = case_when(
    str_count(anti_onset4_s, "\\+") == 0 ~ "Monotherapy (1 drug)",  # No '+' means single drug
    str_count(anti_onset4_s, "\\+") == 1 ~ "Dual therapy (2 drugs)", # One '+' means two drugs
    str_count(anti_onset4_s, "\\+") >= 2 ~ "Triple or more (≥3 drugs)" # Two or more '+' means 3+ drugs
  )) %>%
  distinct() 

table(CR_names_4$arm) # 90 ceftazidime/avibactam,  125 polymyxin monotherapy,  333 polymyxin combination
length(unique(CR_names_4$recordid)) #316 #336  #310  #982  #753  #796 had studied drug on index date


#####
# ICU and MV
icu_mv <- baseline_outcomes_index %>% 
  select (recordid,"icu_hd_ap","icu_hd_ap_1","icu_hd_ap_1_2","icu_hd_ap_1_1",            
          "icu_hd_ap_1_3","icu_hd_ap_2_2","icu_hd_ap_2_1","icu_hd_ap_2_3",           
          "icu_hd_ap_3_2","icu_hd_ap_3_1","icu_hd_ap_3_3","icu_hd_ap_4_2",           
          "icu_hd_ap_4_1","icu_hd_ap_4_3","icu_hd_ap_5_2","icu_hd_ap_5_1",            
          "icu_hd_ap_5_3","mv_ap","mv_ap_1","mv_ap_1_2",               
          "mv_ap_1_1","mv_ap_1_3","mv_ap_2_2","mv_ap_2_1",                
          "mv_ap_2_3","mv_ap_3_2", "mv_ap_3_1","mv_ap_3_3",                
          "mv_ap_4_2","mv_ap_4_1", "mv_ap_4_3","mv_ap_5_2",               
          "mv_ap_5_1","mv_ap_5_3")


CR_combined <- left_join(CR_names_4, icu_mv, by = "recordid") %>% distinct() %>%
  mutate(
    # ICU status at onset4
    icu_at_onset4 = case_when(
      is.na(onset4) ~ NA_integer_,
      (!is.na(icu_hd_ap_1_2) & onset4 >= icu_hd_ap_1_2 & 
         (onset4 <= icu_hd_ap_1_3 | (icu_hd_ap_1_1 == 1 & is.na(icu_hd_ap_1_3)))) ~ 1,
      (!is.na(icu_hd_ap_2_2) & onset4 >= icu_hd_ap_2_2 & 
         (onset4 <= icu_hd_ap_2_3 | (icu_hd_ap_2_1 == 1 & is.na(icu_hd_ap_2_3)))) ~ 1,
      (!is.na(icu_hd_ap_3_2) & onset4 >= icu_hd_ap_3_2 & 
         (onset4 <= icu_hd_ap_3_3 | (icu_hd_ap_3_1 == 1 & is.na(icu_hd_ap_3_3)))) ~ 1,
      (!is.na(icu_hd_ap_4_2) & onset4 >= icu_hd_ap_4_2 & 
         (onset4 <= icu_hd_ap_4_3 | (icu_hd_ap_4_1 == 1 & is.na(icu_hd_ap_4_3)))) ~ 1,
      (!is.na(icu_hd_ap_5_2) & onset4 >= icu_hd_ap_5_2 & 
         (onset4 <= icu_hd_ap_5_3 | (icu_hd_ap_5_1 == 1 & is.na(icu_hd_ap_5_3)))) ~ 1,
      TRUE ~ 0
    ),
    # Ventilator status at onset4
    vent_at_onset4 = case_when(
      is.na(onset4) ~ NA_integer_,
      (!is.na(mv_ap_1_2) & onset4 >= mv_ap_1_2 & 
         (onset4 <= mv_ap_1_3 | (mv_ap_1_1 == 1 & is.na(mv_ap_1_3)))) ~ 1,
      (!is.na(mv_ap_2_2) & onset4 >= mv_ap_2_2 & 
         (onset4 <= mv_ap_2_3 | (mv_ap_2_1 == 1 & is.na(mv_ap_2_3)))) ~ 1,
      (!is.na(mv_ap_3_2) & onset4 >= mv_ap_3_2 & 
         (onset4 <= mv_ap_3_3 | (mv_ap_3_1 == 1 & is.na(mv_ap_3_3)))) ~ 1,
      (!is.na(mv_ap_4_2) & onset4 >= mv_ap_4_2 & 
         (onset4 <= mv_ap_4_3 | (mv_ap_4_1 == 1 & is.na(mv_ap_4_3)))) ~ 1,
      (!is.na(mv_ap_5_2) & onset4 >= mv_ap_5_2 & 
         (onset4 <= mv_ap_5_3 | (mv_ap_5_1 == 1 & is.na(mv_ap_5_3)))) ~ 1,
      TRUE ~ 0
    ),
    # Length of stay before onset4
    los_onset4 = as.integer(difftime(onset4, adm_date, units = "days")),
    # Sum up number of days in ICU before onset4
    iculos_onset4 = rowSums(
      cbind(
        pmax(pmin(onset4, icu_hd_ap_1_3, na.rm = TRUE) - icu_hd_ap_1_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, icu_hd_ap_2_3, na.rm = TRUE) - icu_hd_ap_2_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, icu_hd_ap_3_3, na.rm = TRUE) - icu_hd_ap_3_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, icu_hd_ap_4_3, na.rm = TRUE) - icu_hd_ap_4_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, icu_hd_ap_5_3, na.rm = TRUE) - icu_hd_ap_5_2, 0, na.rm = TRUE)
      ), na.rm = TRUE
    ),
    # Mechanical ventilation duration before onset4
    mvdur_onset4 = rowSums(
      cbind(
        pmax(pmin(onset4, mv_ap_1_3, na.rm = TRUE) - mv_ap_1_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, mv_ap_2_3, na.rm = TRUE) - mv_ap_2_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, mv_ap_3_3, na.rm = TRUE) - mv_ap_3_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, mv_ap_4_3, na.rm = TRUE) - mv_ap_4_2, 0, na.rm = TRUE),
        pmax(pmin(onset4, mv_ap_5_3, na.rm = TRUE) - mv_ap_5_2, 0, na.rm = TRUE)
      ), na.rm = TRUE
    )
  )%>%
  select(-c(icu_hd_ap, icu_hd_ap_1, icu_hd_ap_1_2, icu_hd_ap_1_1, icu_hd_ap_1_3, icu_hd_ap_2_2, icu_hd_ap_2_1, icu_hd_ap_2_3, icu_hd_ap_3_2, icu_hd_ap_3_1, icu_hd_ap_3_3, 
            icu_hd_ap_4_2, icu_hd_ap_4_1, icu_hd_ap_4_3, icu_hd_ap_5_2, icu_hd_ap_5_1, icu_hd_ap_5_3, 
            mv_ap, mv_ap_1, mv_ap_1_2, mv_ap_1_1, mv_ap_1_3, mv_ap_2_2, mv_ap_2_1, mv_ap_2_3, mv_ap_3_2, mv_ap_3_1, mv_ap_3_3, 
            mv_ap_4_2, mv_ap_4_1, mv_ap_4_3, mv_ap_5_2, mv_ap_5_1, mv_ap_5_3))


CR_targeted <- CR_combined %>% 
  select(recordid, inf_onset, spec_date, onset4, delay, arm, arm2, arm3, number, 
         anti_onset4, anti_onset4_s, anti_new, anti_start, organism, icu_at_onset4, vent_at_onset4, 
         los_onset4, iculos_onset4, mvdur_onset4, susceptibility, organism) %>%
  filter(anti_new != "Carbapenem") %>%  
  group_by(recordid) %>%
  filter(delay == min(delay)) %>%  
  ungroup()%>%
  select (recordid, inf_onset, onset4, delay, arm,  arm2, arm3,
          number, anti_onset4, anti_onset4_s, organism,  icu_at_onset4, vent_at_onset4, los_onset4, iculos_onset4, mvdur_onset4) %>%
  distinct()

# join CR_targeted and outcome dataframe, follow CR_targeted recordid
CR_final <- left_join(CR_targeted, outcome, by = "recordid")%>%
  mutate(
    # Follow-up days from dayzero
    fup_day_onset4 = as.integer(pmin(d28_date, d28_death_date, mortality_date, na.rm = TRUE) - onset4),
    # 21-day mortality from dayzero
    mortality_date_final = pmin(d28_death_date, mortality_date, na.rm = TRUE),
    mort_21d_onset4 = ifelse(!is.na(mortality_date_final) & (mortality_date_final - onset4) <= 21, 1, 0),
    mort_14d_onset4 = ifelse(!is.na(mortality_date_final) & (mortality_date_final - onset4) <= 14, 1, 0),
    # Mortality days from dayzero, assign 0 if no mortality
    mortday_onset4 = ifelse(is.na(mortality_date_final), 0, as.integer(mortality_date_final - onset4))) %>% distinct()

CR_mono <- ast_all_index %>% select (recordid, org_combined) %>% filter (recordid %in% CR_final$recordid)%>%
  mutate(mono_poly = ifelse(grepl(",", org_combined), "polymicrobial", "monomicrobial"))%>% distinct()

# 
# length(unique(CR_final$recordid)) #335  #271
# length(unique(CR_mono$recordid)) #331
# # Which recordid is missing in CR_mono
# CR_final[!CR_final$recordid %in% CR_mono$recordid, "recordid"]

#merge CR_final and CR_mono
CR_final <- merge(CR_final, CR_mono, by = "recordid", all.x = TRUE) %>% distinct()
# For recordid CN-001-A0-0032, CN-002-A0-0005, IN-004-A0-0009, assign mono_poly as monomicrobial

#CR_final[CR_final$recordid %in% c("CN-001-A0-0032", "CN-002-A0-0005","IN-004-A0-0007", "IN-004-A0-0009"), "mono_poly"] <- "monomicrobial"
#CR_final[CR_final$recordid %in% c("CN-001-A0-0032", "CN-002-A0-0005","IN-004-A0-0007", "IN-004-A0-0009"), "org_combined"] <- "K. pneumoniae"

length(unique(CR_final$recordid)) #335  #982  #723

CR_final%>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms")

CR_combined_new <- checkast %>%
  left_join(CR_final %>% select(recordid, onset4), by = "recordid")%>%
  filter(spec_date <= onset4) %>%
  filter(!grepl("Coagulase-negative staphylococci|Chryseobacterium|Corynebacterium|Ochrobactrum|Abiotrophia", org_names_all)) %>%
  filter(!grepl("Coagulase-negative staphylococci (CoNS)", org_names_all, fixed = TRUE))%>%
  distinct()

unique(CR_combined_new$org_names_all)

CR_combined_new <- CR_combined_new %>%
  group_by(recordid) %>%
  summarise(
    org_combined_new = paste(unique(org_names_all), collapse = ", "),
    .groups = "drop"
  )%>%
  mutate(monopoly = ifelse(grepl(",", org_combined_new), "polymicrobial", "monomicrobial"))%>% distinct() 

CR_final <- CR_final %>%
  left_join(CR_combined_new, by= "recordid")%>%
  mutate(aci = ifelse(grepl("Acinetobacter|acinetobacter baumini", org_combined_new), 1, 0)) %>%
  mutate(pae= ifelse(grepl("Pseudomonas", org_combined_new), 1, 0)) %>%
  mutate(ent= ifelse(grepl("Enterobacter|Serratia|K. pneumoniae|Klebsiella|Aeromonas|Enterobacter|E. coli|Proteus|Providencia|Morganella|Proteus|Serratia|Citrobacter|Aeromonas", org_combined_new), 1, 0))%>%
  group_by(recordid) %>%
  mutate(
    cr_aci = ifelse(any(organism == "CRAB", na.rm = TRUE), 1, 0),
    cr_ent = ifelse(any(organism == "CRE", na.rm = TRUE), 1, 0),
    cr_pae = ifelse(any(organism == "CRPAE", na.rm = TRUE), 1, 0)
  ) %>%
  ungroup()

CR_complete_outcome <- CR_final %>% filter(!is.na(d28_date) | !is.na(d28_status) | !is.na(d28_death_date) | !is.na(mortality_date))
CR_incomplete_outcome <- CR_final %>% filter(is.na(d28_date) & is.na(d28_status) & is.na(d28_death_date) & is.na(mortality_date))
unique(CR_incomplete_outcome$recordid) # "CN-002-A0-0001" "CN-002-A0-0027" "KH-002-A0-0179"
# "CN-001-A0-0020" "TH-002-A0-0344"

length(unique(CR_final$recordid)) #501  #511 #532  #723
length(unique(CR_final$recordid)) #668 -> 654 (26 apr)  #723 with tigecycline and fosfomycin
nrow(CR_final) #779

CR_final_drop <- CR_final %>% filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(CR_final_drop$recordid)) #638 wihtout tigecycline and fosfomycin  #638
CR_final <- CR_final_drop

#######################################################################################################################################
################# STOP here ####################
#######################################################################################################################################
colnames(CR_final)

# CRAB <- CR_final %>% filter(organism == "CRAB")
# length(unique(CRAB$recordid)) #487  #471
# nrow(CRAB) #487
# CRE <- CR_final %>% filter(organism == "CRE")
# length(unique(CRE$recordid)) #221  #231 
# nrow(CRE) #221
# CRPAE <- CR_final %>% filter(organism == "CRPAE")
# length(unique(CRPAE$recordid)) # 89  #77 
# nrow(CRPAE) #89

CRAB_drop <- CR_final_drop %>% filter(organism == "CRAB")
length(unique(CRAB_drop$recordid)) #429  
nrow(CRAB_drop) 
CRE_drop <- CR_final_drop %>% filter(organism == "CRE")
length(unique(CRE_drop$recordid)) #187 
nrow(CRE_drop) 
CRPAE_drop <- CR_final_drop %>% filter(organism == "CRPAE")
length(unique(CRPAE_drop$recordid)) #65
nrow(CRPAE_drop) 

CRAB <- CR_final %>% filter(cr_aci==1)
#%>% filter(monopoly == "monomicrobial")
length(unique(CRAB$recordid)) #487 (370 mono)   #445 (no tige)
nrow(CRAB) #487   #445
#429 ITT 425 PP
#430 5 days
#466

CRE_CRPAE <- CR_final %>% filter(cr_pae==1 | cr_ent==1)
# %>% filter(monopoly == "monomicrobial")
length(unique(CRE_CRPAE$recordid)) #316 (221 mono)  #264 (no tige)
nrow(CRE_CRPAE) #316
# 246
# 274

CR_final_mono <- CR_final %>% filter(monopoly == "monomicrobial") 
length(unique(CR_final_mono$recordid)) #570
nrow(CR_final_mono) 

#skim(CRAB)
#skim(CRE_CRPAE)

###########################################################################################################################################
check <- CRE_CRPAE %>% select(recordid, onset4, org_combined_new, anti_onset4, arm, arm2, arm3, aci, ent, pae, cr_aci, cr_pae, cr_ent)%>%
  filter(!is.na(arm3))%>%
  filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(check$recordid)) #261 (with tige), #211 (no tige)

crab <- CRAB %>% filter(!is.na(arm2))
#%>% filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(crab$recordid)) #443 (no tige)  vs 485 (with tige)
#427
#428 day 5
#462

table(crab$arm, crab$mort_21d_onset4)

crab %>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRAB") 

crab_mono <- CRAB %>% filter(!is.na(arm2))%>% filter(monopoly =="monomicrobial")
length(unique(crab_mono$recordid)) #345 #358
 
crab_bsi <-crab %>% filter(infection_types == "BSI")
length(unique(crab_bsi$recordid)) #144
crab_vap <- crab %>% filter(infection_types == "VAP")
length(unique(crab_vap$recordid)) # 318

table(crab_bsi$arm)
table(crab_vap$arm)
###====================================#####
cre_crpae <- CRE_CRPAE %>% filter(!is.na(arm3)) 
# %>%
# filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(cre_crpae$recordid)) #261 (with tige), #211 (no tige)  #246
# 196
#216

cre_crpae_mono <- CRE_CRPAE %>% filter(!is.na(arm3))%>% filter(monopoly =="monomicrobial")
length(unique(cre_crpae_mono$recordid)) #162

cre_crpae %>% count(arm3) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRE and CRPAE")

table(cre_crpae$arm3, cre_crpae$mort_21d_onset4)

cre_crpae_bsi <-cre_crpae %>% filter(infection_types == "BSI")
length(unique(cre_crpae_bsi$recordid))  #127
cre_crpae_vap <- cre_crpae %>% filter(infection_types == "VAP")
length(unique(cre_crpae_vap$recordid)) #89

table(cre_crpae_bsi$arm3)
table(cre_crpae_vap$arm3)

