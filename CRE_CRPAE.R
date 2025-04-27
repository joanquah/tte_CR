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

test <- f07e %>% filter (recordid =="MY-005-A0-0173") %>% select (recordid, route_antib_39)

# From clean data
infection_types_index <- inf_add
baseline_outcomes_index <- df_baseline_all
ast_all <- all_ast_merged 
ast_all_index <- all_ast_merged_index
anti_treat
anti_treat_index <- anti_treat_index_new
all_vap_bsi
vap_bsi_index

f07a <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07aNew_DATA_2025-03-05_1010.csv", header = TRUE, sep = ",")
f07e <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07eNew_DATA_2025-03-05_1011.csv", header = TRUE, sep = ",")
f07m <- read.csv("C:/Users/Joan/Downloads/data/raw data/ACORNHAI-F07mNew_DATA_2025-03-05_1011.csv", header = TRUE, sep = ",")

# Run correction code

# Source (to think of how to overcome this)
# source <- f07a %>% select (recordid, bsi_pri_sec, bsi_pri_sec_1___1, "bsi_pri_sec_1___2", "bsi_pri_sec_1___3","bsi_pri_sec_1___4",
#                            "bsi_pri_sec_1___5", "bsi_pri_sec_1___6", "bsi_pri_sec_1___7" , "bsi_pri_sec_1_other")%>%
#   rename(
#     SSTI= bsi_pri_sec_1___1, Pulmonary = bsi_pri_sec_1___2, GI = bsi_pri_sec_1___3,
#     UTI = bsi_pri_sec_1___4, SSI = bsi_pri_sec_1___5, Unknown = bsi_pri_sec_1___6,
#     Other = bsi_pri_sec_1___7, Other_freetext = bsi_pri_sec_1_other) %>%
#   filter(recordid %in% CR_final$recordid)%>%
#   distinct()
# Use slice_min instead of summarise (onset=min(onset))

CR <- ast_all %>% 
  mutate(pathogen_group = ifelse(org_names_all %in% c("Klebsiella", "Citrobacter"), "GNR", pathogen_group)) %>%
  mutate(recordid = ifelse(grepl("_BSI|_VAP", recordid), sub("_BSI|_VAP", "", recordid), recordid)) %>%
  select(recordid, spec_date,ast_date,  org_names_all,
         #inf_onset,, infection_types, org_combined, 
         Meropenem, Doripenem, Ertapenem, Imipenem, ris_Carbapenems,
         org, org_oth, pathogen_group) %>%
  mutate(susceptibility = case_when(
    ris_Carbapenems %in% c(2, 3) ~ "Susceptible",
    ris_Carbapenems == 1 ~ "Resistant",
    ris_Carbapenems %in% c(4, 5) ~ "Unknown"
  )) %>%
  filter(pathogen_group == "GNB") %>%
  filter(!(org_names_all == "Stenotrophomonas")) %>%
  filter(!(susceptibility == "Susceptible"))%>%
  # filter(susceptibility =="Resistant") %>%
  mutate(organism = case_when(
    org_names_all == "Acinetobacter" ~ "CRAB",
    org_names_all %in% c("Pseudomonas", "K. pneumoniae", "E. coli", 
                         "Serratia", "Proteus", "Enterobacter") ~ "CRE_CRPAE",
    TRUE ~ NA_character_  
  ))

length(unique(CR$recordid)) # 3610

outcome <- baseline_outcomes_index %>% 
  select(recordid, date_enrolment,age_new, age_group, country_income, country, country_ab, sex, hpd_adm_date, hpd_hosp_date, inf_onset,
         comorbidities_Chalson, sofa_score,adm_ward_types_new, infection_types, adm_ward_types_new, UTI_source,
         ho_discharge_date, ho_dischargestatus, d28_date,d28_status,d28_death_date,first28_death, 
         mortality,mortality_date, score_respiration,score_coagulation,score_liver,score_CVS,score_CNV, score_renal, fbis_score,
         ) %>%
  mutate(infection_types = ifelse(infection_types %in% c("Hospital-acquired BSI", "Healthcare-associated BSI"), "BSI", infection_types))%>%
  filter(age_group==1)%>% 
  mutate(country_income = ifelse(country_ab == "CN", "Upper middle income", country_income))%>%
  mutate(country_income = ifelse(country_ab == "IN", "Lower middle income", country_income)) %>%
  mutate(adm_date = pmin(hpd_adm_date, hpd_hosp_date, na.rm = TRUE)) 

length(unique(outcome$recordid)) #7498

CR2 <- CR %>% inner_join(outcome, by = "recordid") %>% distinct() %>%
  mutate(
    inf_onset = as.Date(inf_onset),
    spec_date = as.Date(spec_date),
    ast_date = as.Date(ast_date),
    mortality_date = as.Date(mortality_date)) %>%  
  filter(!(mortality_date == inf_onset & !is.na(mortality_date))) %>% 
  filter(!(mortality_date == (inf_onset + 1))| is.na(mortality_date))%>%   #1639 adult pts who had CR and did not die within 1 day of infection onset
  select(recordid, inf_onset, spec_date, ast_date, ho_discharge_date, ho_dischargestatus,d28_date,d28_status,mortality, mortality_date, 
         everything())%>%
filter(!(is.na(ho_discharge_date) & is.na(ho_dischargestatus) & is.na(mortality) & is.na(mortality_date) & is.na(d28_date) & is.na(d28_status)))

length(unique(CR2$recordid)) #2667

# Check for incomplete data-> all cleared, no missing outcome 
# Hashtag line above first
#CR_incomplete <- CR2 %>% filter(is.na(ho_discharge_date) & is.na(ho_dischargestatus) & is.na(mortality) & is.na(mortality_date) & is.na(d28_date) & is.na(d28_status))
#baseline_outcome_incomplete <- outcome %>% filter(recordid %in% CR_incomplete$recordid) %>% filter (is.na(sofa_score))
#CR_incomplete <- baseline_outcomes_index %>% filter(recordid %in% baseline_outcome_incomplete$recordid)
#CR_incomplete2 <- baseline_outcomes_index %>% filter(recordid %in% CR_incomplete$recordid) %>% select(recordid, score_respiration, score_coagulation,score_CVS, score_liver,score_CVS,score_CNV, score_renal, sofa_score)

abx_all <- anti_treat %>% 
  select (-vap_epi, -bsi_epi, -anti_group_ent, -anti_group_pa_ab) %>%
  filter(!is.na(anti_start)) %>%
  filter(!is.na(anti_names)) %>%
  mutate(anti_new = case_when(
    anti_names %in% c("Ertapenem", "Imipenem", "Meropenem", "Doripenem") ~ "Carbapenem",
    anti_names %in% c("Colistin", "Polymyxin B") ~ "Polymyxin",
    anti_names %in% c("Ampicillin/sulbactam", "Cefoperazone/sulbactam") ~ "Sulbactam",
    anti_names == "Ceftazidime/avibactam" ~ "Ceftazidime/avibactam",
    anti_names %in% c("Amikacin", "Gentamicin") ~ "Aminoglycoside",
    anti_names %in% c("Tigecycline", "Minocycline") ~ "Tetracycline",
    TRUE ~ anti_names
  )) 

abx <- abx_all %>%
  group_by(recordid) %>%
  filter(any(anti_new %in% c("Polymyxin", "Sulbactam", "Ceftazidime/avibactam"))) %>%
  ungroup()


#Subjects with missing antibiotics data 
# CR_missing_abx <- CR2 %>% anti_join(abx, by = "recordid") %>% select(recordid) %>% distinct()  # need site to fill
#check_f07e_complete <- f07e %>% filter(recordid %in% CR$recordid) %>% select(recordid, f07e_antibiotic_complete)
#length(unique(CR_missing_abx$recordid))  #62

CR_abx2 <- abx %>% inner_join(CR2, by = "recordid") %>% distinct()
length(unique(CR_abx2$recordid)) #1439

excluded_names <- c(
  "Amphotericin B", "Other_Fluconazole", "Caspofungin","Other_Voriconazole", "Anidulafungin", 'Voriconazole',"Fluconazole","Posaconazole",
  "Clindamycin","Teicoplanin", "Chloramphenicol", "Linezolid", "Daptomycin","Metronidazole", "Nitrofurantoin", "Fusidic acid", "Vancomycin", "Rifampicin", 
  "Azithromycin", "Clarithromycin", "Erythromycin", NA, 
  "Amoxicillin-clavulanate", "Cloxacillin", "Flucloxacillin","Piperacillin-tazobactam","Piperacillin/tazobactam", 
  "Ampicillin","Amoxicillin","Benzylpenicillin","Oxacillin", "Penicillin", 
  "Ceftazidime", "Ceftriaxone", "Cefuroxime",  "Cefepime","Cefoxitin", "Cefotaxime","Cefazolin", 
  "Cefixime", "Cefpirome", "Ceftazidime/clavulanic acid","Cephalexin")
#"Doxycycline" , "Tetracycline", "Minocycline",
  #"Ciprofloxacin", "Levofloxacin" , "Ofloxacin", "Sitafloxacin","Moxifloxacin",
  #"Co-trimoxazole", "Trimethoprim","Trimethoprim/sulfamethoxazole")
  #"Aztreonam", "Gentamicin", "Amikacin" , "Streptomycin", 
  #"Ceftolozane/tazobactum", "Ceftolozane/tazobactam"

CR_abx2 <- CR_abx2[!is.na(CR_abx2$anti_names) & !CR_abx2$anti_names %in% excluded_names, ] 

length(unique(CR_abx2$recordid)) #1439

CR_names <- CR_abx2 %>% 
  mutate(onset4 = inf_onset + 4) %>% relocate(onset4, .after = ast_date) %>%
  mutate(onset5 = inf_onset + 5) %>% relocate(onset5, .after = onset4) %>%
  select(recordid, organism, adm_date, inf_onset, spec_date, ast_date, onset4, onset5, anti_names, anti_new, anti_start, anti_end, susceptibility)%>%
  distinct() %>%
  mutate (delay = as.numeric(difftime(anti_start, inf_onset, units = "days")))%>%
  filter(delay >= -5) %>%
  filter(delay <= 5) %>%
  mutate(dur = as.integer(difftime(anti_end, anti_start, units = "days"))) %>%
  mutate (dur = dur+1) %>%
  #filter(dur>2) %>%
  filter(onset4 >= anti_start & onset4 < anti_end) %>%
  #filter(onset5 >= anti_start & onset5 < anti_end) %>%
  group_by(recordid) %>%
  mutate(
    anti_ast = paste(unique(anti_names[ast_date >= anti_start & ast_date <= anti_end]), collapse = " + "),
    anti_onset4   = paste(unique(anti_names[onset4 >= anti_start & onset4 <= anti_end]), collapse = " + "),
    anti_onset5   = paste(unique(anti_names[onset5 >= anti_start & onset5 <= anti_end]), collapse = " + ")
  ) %>%
  ungroup()%>% 
  mutate(spec_day = as.integer(difftime(spec_date, inf_onset, units = "days")))%>%
  relocate(spec_day, .after = spec_date) %>%
  distinct() %>%
  filter(spec_day<6)

length(unique(CR_names$recordid)) #For onset 4 988 for ITT, 976 for PP; For onset5 94 for ITT, 941 for PP 
#CR_names <- CR_names %>% select (recordid, organism, susceptibility, anti_ast, anti_onset4,anti_onset5) %>% distinct

CR_names_4 <- CR_names %>%
  filter(organism == "CRE_CRPAE") %>%
  mutate(arm = case_when(
    anti_onset4 %in% c("Polymyxin B", "Colistin","Polymyxin B + Colistin") ~ "Polymyxin monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset4) & grepl("Meropenem|Imipenem|Tigecycline|Fosfomycin|Ciprofloxacin|Levofloxacin|Moxifloxacin|Aztreonam|Trimethoprim/sulfamethoxazole|Ampicllin/sulbactam|Minocycline|Amikacin", anti_onset4) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam", anti_onset4) ~ "Cefoperazone/sulbactam based",
    grepl("Ceftazidime/avibactam", anti_onset4) ~ "Ceftazidime/avibactam based",
    TRUE ~ NA_character_
  )) %>%
  filter(grepl("Colistin|Polymyxin B|Ceftazidime/avibactam|Cefoperazone/sulbactam", anti_onset4))%>%
  mutate(arm = ifelse(anti_onset4 == "Polymyxin B + Ampicillin/sulbactam", "Polymyxin combination", arm)) %>%
  mutate(arm = ifelse(anti_onset4 == "Ampicillin/sulbactam + Polymyxin B", "Polymyxin combination", arm)) %>%
  mutate(arm = ifelse(anti_onset4 == "Colistin + Ampicillin/sulbactam", "Polymyxin combination", arm)) %>%
  filter(!(is.na(arm)))

length(unique(CR_names_4$recordid)) #316

CR_names_5 <- CR_names %>%
  filter(organism == "CRE_CRPAE") %>%
  mutate(arm = case_when(
    anti_onset5 %in% c("Polymyxin B", "Colistin","Polymyxin B + Colistin") ~ "Polymyxin monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Meropenem|Imipenem|Tigecycline|Fosfomycin|Ciprofloxacin|Levofloxacin|Moxifloxacin|Aztreonam|Trimethoprim/sulfamethoxazole|Ampicllin/sulbactam|Minocycline|Amikacin", anti_onset5) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam", anti_onset5) ~ "Cefoperazone/sulbactam based",
    grepl("Ceftazidime/avibactam", anti_onset5) ~ "Ceftazidime/avibactam based",
    TRUE ~ NA_character_
  )) %>%
  filter(grepl("Colistin|Polymyxin B|Ceftazidime/avibactam|Cefoperazone/sulbactam", anti_onset5))%>%
  mutate(arm = ifelse(anti_onset5 == "Polymyxin B + Ampicillin/sulbactam", "Polymyxin combination", arm)) %>%
  mutate(arm = ifelse(anti_onset5 == "Ampicillin/sulbactam + Polymyxin B", "Polymyxin combination", arm)) %>%
  mutate(arm = ifelse(anti_onset5 == "Colistin + Ampicillin/sulbactam", "Polymyxin combination", arm)) %>%
  filter(!(is.na(arm)))

#Distribution of study arms
CR_names_4 %>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms")
#377/100/90/60 (5) VS 405/106/85/53 (4)


length(unique(CR_names_4$recordid)) #316
length(unique(CR_names_5$recordid)) #316

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

CR_combined <- left_join(CR_names_4, icu_mv, by = "recordid") %>% distinct()
CR_combined <- left_join(CR_names_5, icu_mv, by = "recordid") %>% distinct()

# remove columns "Meropenem", "Doripenem", "Ertapenem","Imipenem","ris_Carbapenems","Colistin","Polymyxin B","ris_Polymyxins", "Fosfomycin","Tigecycline","ris_Glycylcyclines","Ceftazidime/avibactam","Piperacillin/tazobactam","Ceftazidime","Moxifloxacin","Minocycline","Ceftolozane/tazobactam","Ampicillin","Ceftriaxone"
#select(-c("Meropenem", "Doripenem", "Ertapenem","Imipenem","ris_Carbapenems","Colistin","Polymyxin B","ris_Polymyxins", "Fosfomycin",
#          "Tigecycline","ris_Glycylcyclines","Ceftazidime/avibactam","Piperacillin/tazobactam","Ceftazidime","Moxifloxacin","Minocycline",
#         "Ceftolozane/tazobactam","Ampicillin","Ceftriaxone", "ris_Phosphonics", "Cefoperazone/sulbactam", "Ampicillin/sulbactam"))%>%

CR_combined <- CR_combined %>%
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

CR_combined <- CR_combined %>%
  mutate(
    # ICU status at onset5
    icu_at_onset5 = case_when(
      is.na(onset5) ~ NA_integer_,
      (!is.na(icu_hd_ap_1_2) & onset5 >= icu_hd_ap_1_2 & 
         (onset5 <= icu_hd_ap_1_3 | (icu_hd_ap_1_1 == 1 & is.na(icu_hd_ap_1_3)))) ~ 1,
      (!is.na(icu_hd_ap_2_2) & onset5 >= icu_hd_ap_2_2 & 
         (onset5 <= icu_hd_ap_2_3 | (icu_hd_ap_2_1 == 1 & is.na(icu_hd_ap_2_3)))) ~ 1,
      (!is.na(icu_hd_ap_3_2) & onset5 >= icu_hd_ap_3_2 & 
         (onset5 <= icu_hd_ap_3_3 | (icu_hd_ap_3_1 == 1 & is.na(icu_hd_ap_3_3)))) ~ 1,
      (!is.na(icu_hd_ap_4_2) & onset5 >= icu_hd_ap_4_2 & 
         (onset5 <= icu_hd_ap_4_3 | (icu_hd_ap_4_1 == 1 & is.na(icu_hd_ap_4_3)))) ~ 1,
      (!is.na(icu_hd_ap_5_2) & onset5 >= icu_hd_ap_5_2 & 
         (onset5 <= icu_hd_ap_5_3 | (icu_hd_ap_5_1 == 1 & is.na(icu_hd_ap_5_3)))) ~ 1,
      TRUE ~ 0
    ),
    # Ventilator status at onset4
    vent_at_onset5 = case_when(
      is.na(onset5) ~ NA_integer_,
      (!is.na(mv_ap_1_2) & onset5 >= mv_ap_1_2 & 
         (onset5 <= mv_ap_1_3 | (mv_ap_1_1 == 1 & is.na(mv_ap_1_3)))) ~ 1,
      (!is.na(mv_ap_2_2) & onset5 >= mv_ap_2_2 & 
         (onset5 <= mv_ap_2_3 | (mv_ap_2_1 == 1 & is.na(mv_ap_2_3)))) ~ 1,
      (!is.na(mv_ap_3_2) & onset5 >= mv_ap_3_2 & 
         (onset5 <= mv_ap_3_3 | (mv_ap_3_1 == 1 & is.na(mv_ap_3_3)))) ~ 1,
      (!is.na(mv_ap_4_2) & onset5 >= mv_ap_4_2 & 
         (onset5 <= mv_ap_4_3 | (mv_ap_4_1 == 1 & is.na(mv_ap_4_3)))) ~ 1,
      (!is.na(mv_ap_5_2) & onset5 >= mv_ap_5_2 & 
         (onset5 <= mv_ap_5_3 | (mv_ap_5_1 == 1 & is.na(mv_ap_5_3)))) ~ 1,
      TRUE ~ 0
    ),
    # Length of stay before onset4
    los_onset5 = as.integer(difftime(onset5, adm_date, units = "days")),
    # Sum up number of days in ICU before onset4
    iculos_onset5 = rowSums(
      cbind(
        pmax(pmin(onset5, icu_hd_ap_1_3, na.rm = TRUE) - icu_hd_ap_1_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_2_3, na.rm = TRUE) - icu_hd_ap_2_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_3_3, na.rm = TRUE) - icu_hd_ap_3_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_4_3, na.rm = TRUE) - icu_hd_ap_4_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_5_3, na.rm = TRUE) - icu_hd_ap_5_2, 0, na.rm = TRUE)
      ), na.rm = TRUE
    ),
    # Mechanical ventilation duration before onset4
    mvdur_onset5 = rowSums(
      cbind(
        pmax(pmin(onset5, mv_ap_1_3, na.rm = TRUE) - mv_ap_1_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, mv_ap_2_3, na.rm = TRUE) - mv_ap_2_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, mv_ap_3_3, na.rm = TRUE) - mv_ap_3_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, mv_ap_4_3, na.rm = TRUE) - mv_ap_4_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, mv_ap_5_3, na.rm = TRUE) - mv_ap_5_2, 0, na.rm = TRUE)
      ), na.rm = TRUE
    )
  )%>%
  select(-c(icu_hd_ap, icu_hd_ap_1, icu_hd_ap_1_2, icu_hd_ap_1_1, icu_hd_ap_1_3, icu_hd_ap_2_2, icu_hd_ap_2_1, icu_hd_ap_2_3, icu_hd_ap_3_2, icu_hd_ap_3_1, icu_hd_ap_3_3, 
            icu_hd_ap_4_2, icu_hd_ap_4_1, icu_hd_ap_4_3, icu_hd_ap_5_2, icu_hd_ap_5_1, icu_hd_ap_5_3, 
            mv_ap, mv_ap_1, mv_ap_1_2, mv_ap_1_1, mv_ap_1_3, mv_ap_2_2, mv_ap_2_1, mv_ap_2_3, mv_ap_3_2, mv_ap_3_1, mv_ap_3_3, 
            mv_ap_4_2, mv_ap_4_1, mv_ap_4_3, mv_ap_5_2, mv_ap_5_1, mv_ap_5_3))

CR_targeted <- CR_combined %>% 
  select(recordid, inf_onset, onset4, delay, arm, anti_new, anti_start, organism, icu_at_onset4, vent_at_onset4, los_onset4, iculos_onset4, mvdur_onset4, susceptibility, organism) %>%
  filter(anti_new != "Carbapenem") %>%  
  group_by(recordid) %>%
  filter(delay == min(delay)) %>%  
  ungroup()%>%
  select (recordid, inf_onset, onset4, delay, arm,  organism,  icu_at_onset4, vent_at_onset4, los_onset4, iculos_onset4, mvdur_onset4) %>%
  distinct()
  # filter(n() > 1 & n_distinct(anti_start) > 1) 
  # mortality_date_final, infection_types, mono_poly)

CR_targeted <- CR_combined %>% 
  select(recordid, inf_onset, onset5, delay, arm, anti_new, anti_start, organism, icu_at_onset5, vent_at_onset5, los_onset5, iculos_onset5, mvdur_onset5, susceptibility, organism) %>%
  filter(anti_new != "Carbapenem") %>%  
  group_by(recordid) %>%
  filter(delay == min(delay)) %>%  
  ungroup()%>%
  select (recordid, inf_onset, onset5, delay, arm,  organism,  icu_at_onset5, vent_at_onset5, los_onset5, iculos_onset5, mvdur_onset5) %>%
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

CR_final <- left_join(CR_targeted, outcome, by = "recordid")%>%
  mutate(
    # Follow-up days from dayzero
    fup_day_onset5 = as.integer(pmin(d28_date, d28_death_date, mortality_date, na.rm = TRUE) - onset5),
    # 21-day mortality from dayzero
    mortality_date_final = pmin(d28_death_date, mortality_date, na.rm = TRUE),
    mort_21d_onset5 = ifelse(!is.na(mortality_date_final) & (mortality_date_final - onset5) <= 21, 1, 0),
    mort_14d_onset5 = ifelse(!is.na(mortality_date_final) & (mortality_date_final - onset5) <= 14, 1, 0),
    # Mortality days from dayzero, assign 0 if no mortality
    mortday_onset5 = ifelse(is.na(mortality_date_final), 0, as.integer(mortality_date_final - onset5))) %>% distinct()

CR_mono <- ast_all_index %>% select (recordid, org_combined) %>% filter (recordid %in% CR_final$recordid)%>%
  mutate(mono_poly = ifelse(grepl(",", org_combined), "polymicrobial", "monomicrobial"))%>% distinct()


length(unique(CR_final$recordid)) #316 (4), 316(5)
length(unique(CR_mono$recordid)) #318 (4), 313(5)
# Which recordid is missing in CR_mono
CR_final[!CR_final$recordid %in% CR_mono$recordid, "recordid"]

#merge CR_final and CR_mono
CR_final <- merge(CR_final, CR_mono, by = "recordid", all.x = TRUE) %>% distinct()
# For recordid CN-001-A0-0032, CN-002-A0-0005, IN-004-A0-0009, assign mono_poly as monomicrobial
CR_final[CR_final$recordid %in% c("CN-001-A0-0032", "CN-002-A0-0005", "IN-004-A0-0009"), "mono_poly"] <- "monomicrobial"


CR_final <- CR_final %>%
  left_join(CR_sofa_sum, by = "recordid") %>%
  left_join(CR_sofa_median, by = "recordid") %>%
  left_join(qsofa, by = "recordid") %>%
  left_join(cmb, by = "recordid")%>%
  left_join(cre, by = "recordid") %>%
  #for country_income2, change to 2 groups upper middle income and high income, low income and low middle income
  mutate(country_income2 = ifelse(country_income == "Low income" | country_income == "Lower middle income", "Low income and Lower middle income", "Upper middle income and High income")) %>%
  relocate(country_income2, .after = country_income) 

# Convert columns to appropriate classes
CR_final <- CR_final %>%
  mutate(delay = ifelse(delay < 0, 0, delay))%>%
  mutate(delay_group = ifelse(delay == 0, 0, 1)) %>%
  mutate(
    # Continuous variables as numeric
    across(c(age_new, sofa_score, comorbidities_Chalson,  
             los_onset4, iculos_onset4, mvdur_onset4, mortday_onset4, delay, 
             sofa_score_sum, sofa_score_median, fbis_score,fup_day_onset4, qsofa,crea), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset4, vent_at_onset4, 
             mort_21d_onset4, mono_poly, UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, 
             delay_group, country_income2), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined, arm, organism,ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

# Convert columns to appropriate classes
CR_final <- CR_final %>%
  mutate(delay = ifelse(delay < 0, 0, delay))%>%
  mutate(delay_group = ifelse(delay == 0, 0, 1)) %>%
  mutate(
    # Continuous variables as numeric
    across(c(age_new, sofa_score, comorbidities_Chalson,  
             los_onset5, iculos_onset5, mvdur_onset5, mortday_onset5, delay, 
             sofa_score_sum, sofa_score_median, fbis_score,fup_day_onset5, qsofa,crea), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset5, vent_at_onset5, 
             mort_21d_onset5, mono_poly, UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, 
             delay_group, country_income2), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined, arm, organism,ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

CR_final <- CR_final %>% 
  # remove these participants TEMPRORAILY due to missing admission details and outcome data (To prioritize)
  filter(recordid != "CN-002-A0-0005")%>%
  filter(recordid != "CN-002-A0-0044")%>%
  filter(recordid != "CN-002-A0-0048")%>%
  filter(recordid != "CN-002-A0-0056")%>%
  filter(recordid != "CN-002-A0-0058")%>%
  filter(recordid != "CN-002-A0-0065")
# "CN-002-A0-0005,CN-002-A0-0044,CN-002-A0-0048, CN-002-A0-0058	, CN-002-A0-0065, CN-002-A0-0056
length(unique(CR_final$recordid)) #310 (4), 311(5) 

CR_complete_outcome <- CR_final %>%
  filter(
    !is.na(d28_date) | 
      !is.na(d28_status) | 
      !is.na(d28_death_date) | 
      !is.na(mortality_date)
  )

CR_incomplete_outcome <- CR_final %>%
  filter(
    is.na(d28_date) & 
      is.na(d28_status) & 
      is.na(d28_death_date) & 
      is.na(mortality_date)
  )

####=================================================###

# CR_new_4 <- CR_new %>% select (recordid, organism, susceptibility, anti_onset4) %>% distinct %>%
#   filter(grepl("Sulbactam|Ceftazidime/avibactam|Polymyxin", anti_onset4))%>%
#   filter(organism == "CRE_CRPAE")
# CR_new_5 <- CR_new %>% select (recordid, organism, susceptibility, anti_onset5) %>% distinct %>%
#   filter(grepl("Sulbactam|Ceftazidime/avibactam|Polymyxin", anti_onset5))
# CR_ast <- CR_new %>% select (recordid, organism, susceptibility, anti_ast) %>% distinct %>%
  # filter(grepl("Sulbactam|Ceftazidime/avibactam|Polymyxin", anti_ast))

# length(unique(CR_ast$recordid)) #853
# length(unique(CR_new_4$recordid)) #372
# length(unique(CR_new_5$recordid)) #859
# length(unique(CR_names_4$recordid)) #324

# Check ceftazidime avibactam ST
# CR_ceftazidime_avibactam <- CR %>% filter(anti_new == "Ceftazidime/avibactam")
# ceftazidime_avibactam <- f07m %>% select(recordid, org, spec_date, ast_date, susceptibility_7, res_7,std_7, mic_7,mic_ug_7, mic_mg_7)%>%
#   filter(!is.na(susceptibility_7))%>%
#   filter(recordid %in% CR_ceftazidime_avibactam$recordid)%>%
#   distinct()
# 
# CR_polymyxin  <- CR %>% filter(anti_new == "Polymyxin")
# coli_poly <- f07m %>% select(recordid, org, spec_date, ast_date, susceptibility_26, res_26,std_26, mic_26,mic_ug_26, mic_mg_26,
#                              susceptibility_45, res_45,std_45, mic_45,mic_ug_45, mic_mg_45)%>%
#   filter(recordid %in% CR_polymyxin$recordid)%>%
#   # remove rows when susceptibility_26 is NA or suscetibility_45 is NA
#   filter(!is.na(susceptibility_26) | !is.na(susceptibility_45))%>%
#   distinct()
# 
# # List all recordid in CR_polymyxin tha is not avalable in coli_poly dataframe
# CR_missing_polymyxinST <- CR_polymyxin %>% anti_join(coli_poly, by = "recordid") %>% distinct()  # need site to fill
# 
# CR_sulbactam  <- CR %>% filter(anti_new == "Sulbactam") %>%
#   mutate(spec_date = as.Date(spec_date))  # Convert spec_date to Date type
# sulbactam <- f07m %>%
#   mutate(spec_date = as.Date(spec_date)) %>%
#   inner_join(CR_abx_simp, by = c("recordid", "spec_date", "org")) %>%  # Match rows exactly 
#   select(recordid, org, spec_date, susceptibility_4, res_4,std_4, mic_4,mic_ug_4, mic_mg_4, susceptibility_12, res_12,std_12, mic_12,mic_ug_12, mic_mg_12)%>%
#   filter(!is.na(susceptibility_4) | !is.na(susceptibility_12))%>%
#   distinct()
# 
# CR_targeted_sulbactam <- CR_targeted %>%
#   filter(targetedabx=="Sulbactam")
# 
# sulbactam <- f07m %>% select(recordid, org, spec_date, ast_date, susceptibility_4, res_4,std_4, mic_4,mic_ug_4, mic_mg_4,
#                              susceptibility_12, res_12,std_12, mic_12,mic_ug_12, mic_mg_12)%>%
#   filter(recordid %in% CR_sulbactam$recordid & 
#            spec_date %in% CR_sulbactam$spec_date & 
#            org %in% CR_sulbactam$org) %>%
#   # remove rows when susceptibility_26 is NA or suscetibility_45 is NA
#   filter(!is.na(susceptibility_4) | !is.na(susceptibility_12))%>%
#   distinct() %>%
#   filter(recordid %in% CR_targeted_sulbactam$recordid)
# 
# # Name all recordid with anti_names== Colistin and ris_Polymyxins==5 in CR dataframe
# CR_ceftaz_avi_missing <- CR %>% filter(anti_names == "Ceftazidime/avibactam" & (`Ceftazidime/avibactam` == 5 | `Ceftazidime/avibactam` == 1)) # Use back ticks
# 
# # Check sofa date
# sofa_date <- sofa_all %>% select (recordid, sofa_res_date, sofa_co_date, sofa_liv_date, sofa_ren)%>%
#   filter(recordid %in% CR$recordid)
# 
# # Which recordid is available in CR but not available in sofa_date
# # CR_missing_sofa <- CR %>% anti_join(sofa_date, by = "recordid") %>% distinct()  # need site to fill
# 
# length(unique(CR$recordid))
# length(unique(sofa_date$recordid))


#Distribution of study arms
CR_final %>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms")

# Export to excel
# library(writexl)
# write_xlsx(CR_combined, "CR_combined.xlsx")


############################ 
## TARGET TRIAL EMULATION ##
############################
# cci
cmb <- baseline_outcomes_index %>% select (recordid, hpd_admreason, cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd, cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab, cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___hivwa, cmb_comorbidities___hivna, cmb_comorbidities___mlr, cmb_comorbidities___mal, cmb_comorbidities___mst, cmb_comorbidities___mld, cmb_comorbidities___liv, cmb_comorbidities___pep, cmb_comorbidities___renal, cmb_comorbidities___tub, cmb_comorbidities___oth)%>%
  filter(recordid %in% CR_final$recordid)

# BSI only-immunosuppresion status and Pitts
comorb <- f07a %>% select (recordid, inf_onset, hiv, end, ins, am, cc, pt, child_c, neu,hae, sot,
                           temp, temp_1, res_rate, res_rate_1, heart_rate, heart_rate_1, sys_bp, sys_bp_1, men_sta, 
                           acute, iv_agents, mv,ca)%>%
  filter(recordid %in% CR_final$recordid) 

# qsofa
qsofa <- severity_new %>% filter(recordid %in% CR_final$recordid) %>% mutate(qsofa = severity_score)%>% select (-severity_score, -severity_score_scale, -adm_ward_types_ori, -adm_ward_types_new)
qsofa_detailed <- severity %>% filter(recordid %in% CR_final$recordid) %>% 
  select (recordid,redcap_event_name, ser_gcs_under15, ser_rr_22up, ser_sbp_under100, ser_abnormal_temp_adult, hai_icu48days, hai_surg,mic_bloodcollect,mic_rec_antibiotic, 
          ser_gcs_under15_new,ser_gcs_under15_new,ser_rr_22up_new, ser_sbp_under100_new)

# SOFA
colnames(sofa_all) <- make.names(colnames(sofa_all), unique = TRUE)

CR_sofa_all <- sofa_all %>% filter(recordid %in% CR_final$recordid)
CR_sofa_new <- sofa_new %>% filter(recordid %in% CR_final$recordid)
median_sofa <- median(CR_sofa_new$sofa_score, na.rm = TRUE)

# Ignore NA and Sum Available Scores (the total SOFA score for that patient may be underestimated, but it might still be valid depending on your study goal, This method is acceptable if the missingness is random and the missing values are not due to a systematic issue.)
CR_sofa_new$sofa_score_sum <- rowSums(CR_sofa_new[, c("score_respiration", "score_coagulation", 
                                                      "score_liver", "score_CVS", 
                                                      "score_CNV", "score_renal")], na.rm = TRUE)
CR_sofa_sum <- CR_sofa_new %>% select(recordid, sofa_score_sum)
# Impute missing values in each SOFA sub-score (e.g., score_respiration, score_coagulation, etc.) using mean or median:
CR_sofa_median_imputed <- CR_sofa_new %>% select(-sofa_score, -sofa_score_sum) %>%
  mutate(across(starts_with("score_"), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  rowwise() %>%
  mutate(sofa_score_median = sum(c_across(starts_with("score_"))))
CR_sofa_median <- CR_sofa_median_imputed %>% select(recordid, sofa_score_median)

CR_final <- CR_final %>%
  left_join(CR_sofa_sum, by = "recordid") %>%
  left_join(CR_sofa_median, by = "recordid") %>%
  left_join(qsofa, by = "recordid") %>%
  left_join(cmb, by = "recordid") 

# # Multiple imputation
# library(mice)
# md.pattern(CR_sofa_new)
# imputed_data <- mice(CR_sofa_new, m = 5, method = "pmm", seed = 123) # m=5 creates 5 imputed dataset, method pmm uses predictive mean matching, seed 123 ensures reproducibility
# summary(imputed_data)
# CR_sofa_imputed <- complete(imputed_data, action = "long") 
# CR_sofa_imputed <- CR_sofa_imputed %>%
#   mutate(sofa_score_mice = rowSums(select(., starts_with("score_")), na.rm = TRUE))
# pooled_model <- with(imputed_data, lm(sofa_score ~ age_new + sex + country_income))
# summary(pool(pooled_model))
# To be continued using mice

# Convert columns to appropriate classes
CR_final <- CR_final %>%
  # if delay is negative value, assign 0
  mutate(delay = ifelse(delay < 0, 0, delay))%>%
  mutate(delay_group = ifelse(delay == 0, 0, 1)) %>%
  mutate(
    # Continuous variables as numeric
    across(c(age_new, sofa_score, comorbidities_Chalson,  
             los_onset4, iculos_onset4, mvdur_onset4, mortday_onset4, delay, 
             sofa_score_sum, sofa_score_median, fbis_score,fup_day_onset4, qsofa), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset4, vent_at_onset4, 
             mort_21d_onset4, mono_poly, UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, 
             delay_group), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined, arm, organism,ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

# Check structure
str(CR_final)

library(gtsummary)
# Stratify age into 18-40,41-60,51-60, 61-80, 81-100
CR_final$age_group <- cut(CR_final$age_new, breaks = c(17, 40, 60, 80, 100), labels = c("18-40", "41-60", "61-80", "81-100"))
CR_final$cci_group <- cut(CR_final$comorbidities_Chalson, 
                          breaks = c(-1, 0, 2, 4, Inf), 
                          labels = c("0", "1-2", "3-4", "5 and above"))
skim (CR_final)


# Baseline demographic table with p-value, for abx_outcome, intervention is Arm, outcome is mortality
table1 <- CR_final %>% 
  select(arm, age_new, age_group,sex, country_income, country_income2,
         comorbidities_Chalson, cci_group,sofa_score_sum, sofa_score_median, qsofa, adm_ward_types_new,infection_types, 
         icu_at_onset4, vent_at_onset4, los_onset4, iculos_onset4, mvdur_onset4, mono_poly,
         cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
         cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
         cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
         cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
         cmb_comorbidities___mal, cmb_comorbidities___mlr,
         delay,delay_group, crea)%>% 
  #sofa_score_sum_group,
  tbl_summary(
    by = arm,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})", # Show median (IQR) instead of mean (SD)
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany" # Show missing data only if present
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",  # Kruskal-Wallis test for continuous variables
    all_categorical() ~ "chisq.test" # Chi-squared test for categorical variables
  ))
  
print(table1)




# Crude analysis
# Fit logistic model and estimate 21 day mortality risks
model1 <- glm(mort_21d_onset4~ arm, data = CR_final, family = "binomial")
summary(model1)

# Generate new dataset to store estimated risks
new_data <- data.frame(arm = levels(CR_final$arm))

# list unique CR_targeted_abx_name
levels(CR_final$arm)

# Obtained estimated risks based on model
new_data$risk <- predict(model1, newdata = new_data, type = "response")
new_data

# Estimate risk difference and risk ratio from logistic model
#new_data$risk[2] - new_data$risk[1]
#new_data$risk[2] / new_data$risk[1]

# Get predictions with standard errors
pred <- predict(model1, newdata = new_data, type = "link", se.fit = TRUE)

# Convert log-odds to probability scale
new_data$risk <- plogis(pred$fit)

# Compute confidence intervals (logit scale)
new_data$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
new_data$upper <- plogis(pred$fit + 1.96 * pred$se.fit)

# View the final dataset with estimated risks and confidence intervals
print(new_data)

# Horizontal Error Bar Plot (Whisker and Dot Plot)
ggplot(new_data, aes(y = arm, x = risk)) +
  geom_point(size = 4) +  # Dot for point estimate
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +  # Horizontal error bar
  labs(title = "21-Day Mortality Risk by Treatment",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal()

# Expand dataset based on survival time and fits pooled logistic regression models to estimate risks over time
# Creating survtime variable that represents the duration (in days) of follow-up
CR_final$survtime <-ifelse(CR_final$mortday_onset4 == 0, 21, CR_final$mortday_onset4)

# Expanding the dataset using the survtime variable
CR_final.surv <- uncount(CR_final, survtime,.remove = F)

# Creating variables for time
CR_final.surv <- CR_final.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()

# Creating variable for timesq
CR_final.surv$timesq <- CR_final.surv$time^2

# Creating event variable
CR_final.surv$event <- ifelse(
  CR_final.surv$time == CR_final.surv$survtime-1 & CR_final.surv$mort_21d_onset4 == 1, 1, 0)

# Fitting a pooled logistic regression model with time, treatment and product term between treatment and time
#fit.pool1 <-  glm(event ~ CR_targeted_abx_name + time + CR_targeted_abx_name*time, 
#                  data = CR_final.surv, 
#                  family = "binomial")
#summary(fit.pool1)
#exp(fit.pool1$coefficients)

# Fitting a pooled logistic regression model with time (linear and quadratic terms), treatment and product terms between treatment and time
#fit.pool2 <-  glm(event ~ CR_targeted_abx_name + time + timesq + CR_targeted_abx_name*time + CR_targeted_abx_name*timesq, 
#                  data = CR_final.surv, 
#                 family = "binomial")
#summary(fit.pool2)
#exp(fit.pool2$coefficients)

# pooled logistic regression 
# Fitting pooled logistic regression model with treatment and time (linear and quadratic terms)
fit.pool3 <- glm(event ~ CR_targeted_abx_name + time + timesq,  
                 data = CR_final.surv, 
                 family = "binomial")  # Without product term
summary(fit.pool3)
exp(fit.pool3$coefficients)

# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0 <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                       time =seq(0,21), timesq=seq(0,21)^2)
results1 <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                       time =seq(0,21), timesq=seq(0,21)^2)
results2 <- data.frame(CR_targeted_abx_name = "Ceftazidime/avibactam",
                       time =seq(0,21), timesq=seq(0,21)^2)

# Obtain predicted hazards from pooled logistic regression model
results0$hazard0 <- predict(fit.pool3, results0, type="response")
results1$hazard1 <- predict(fit.pool3, results1, type="response")
results2$hazard2 <- predict(fit.pool3, results2, type="response")

# Estimate survival probabilities from hazards
# S(t) = cumulative product of (1 - h(t))
results0$surv0 <- cumprod(1-results0$hazard0)
results1$surv1 <- cumprod(1-results1$hazard1)
results2$surv2 <- cumprod(1-results2$hazard2)

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
results0$risk0 <- 1 - results0$surv0
results1$risk1 <- 1 - results1$surv1
results2$risk2 <- 1 - results2$surv2

#View(results0)
#View(results1)
#View(results2)

# Constructing risk curves
# Combine results for each treatment group into a single dataset
results_combined <-merge(results0, results1, by=c("time", "timesq"))
results_combined <-merge(results_combined, results2, by=c("time", "timesq"))

# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results_combined$time_updated <- results_combined$time + 1
results_combined <- results_combined %>% add_row(time_updated=0, risk0=0, risk1=0, risk2=0)%>%
  arrange(time_updated)

# Creating plot
ggplot(results_combined, aes(x=time_updated)) + 
  geom_step(aes(y=risk0, color="Polymyxin monotherapy")) +
  geom_step(aes(y=risk1, color="Polymyxin + Carbapenem")) +
  geom_step(aes(y=risk2, color="Ceftazidime/avibactam")) +
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks=seq(0,21,1)) +
  ylab("Cumulative Incidence of mortality") +
  labs(colour = "Treatment") +
  theme_bw() + 
  theme(legend.position="bottom")

#  Risk at end of follow-up
risk0 <-results0[results0$time==21,]$risk0
risk1 <-results1[results1$time==21,]$risk1
risk2 <-results2[results2$time==21,]$risk2
risk0
risk1
risk2

# Risk ratio and risk difference
risk1 - risk0
risk1 / risk0
head(results_combined)
# Cox proportional hazards model for comparison

# Fit Cox model
fit.cox <- coxph(Surv(survtime, event) ~ arm, data = CR_final.surv)
summary(fit.cox)
exp(fit.cox$coefficients)

# This is to show that pooled logistic regression can approximate cox.

# To generate confidence intervals for risk estimates at time = 20
# Obtain predicted hazards and standard errors
pred0 <- predict(fit.pool3, results0, type = "link", se.fit = TRUE)
pred1 <- predict(fit.pool3, results1, type = "link", se.fit = TRUE)
pred2 <- predict(fit.pool3, results2, type = "link", se.fit = TRUE)

# Convert logit scale hazards to probabilities
results0$hazard0 <- plogis(pred0$fit)
results1$hazard1 <- plogis(pred1$fit)
results2$hazard2 <- plogis(pred2$fit)

# Compute standard errors of hazards on probability scale using delta method
results0$se_hazard0 <- pred0$se.fit * results0$hazard0 * (1 - results0$hazard0)
results1$se_hazard1 <- pred1$se.fit * results1$hazard1 * (1 - results1$hazard1)
results2$se_hazard2 <- pred2$se.fit * results2$hazard2 * (1 - results2$hazard2)

# Compute cumulative survival probabilities
results0$surv0 <- cumprod(1 - results0$hazard0)
results1$surv1 <- cumprod(1 - results1$hazard1)
results2$surv2 <- cumprod(1 - results2$hazard2)

# Compute cumulative standard errors using Greenwood’s formula
results0$se_surv0 <- results0$surv0 * sqrt(cumsum((results0$se_hazard0 / (1 - results0$hazard0))^2))
results1$se_surv1 <- results1$surv1 * sqrt(cumsum((results1$se_hazard1 / (1 - results1$hazard1))^2))
results2$se_surv2 <- results2$surv2 * sqrt(cumsum((results2$se_hazard2 / (1 - results2$hazard2))^2))

# Compute risk estimates
results0$risk0 <- 1 - results0$surv0
results1$risk1 <- 1 - results1$surv1
results2$risk2 <- 1 - results2$surv2

# Compute confidence intervals for risks
results0$lower_CI_risk0 <- 1 - (results0$surv0 * exp(1.96 * results0$se_surv0 / results0$surv0))
results0$upper_CI_risk0 <- 1 - (results0$surv0 * exp(-1.96 * results0$se_surv0 / results0$surv0))

results1$lower_CI_risk1 <- 1 - (results1$surv1 * exp(1.96 * results1$se_surv1 / results1$surv1))
results1$upper_CI_risk1 <- 1 - (results1$surv1 * exp(-1.96 * results1$se_surv1 / results1$surv1))

results2$lower_CI_risk2 <- 1 - (results2$surv2 * exp(1.96 * results2$se_surv2 / results2$surv2))
results2$upper_CI_risk2 <- 1 - (results2$surv2 * exp(-1.96 * results2$se_surv2 / results2$surv2))

risk0 <- results0[results0$time == 21, c("risk0", "lower_CI_risk0", "upper_CI_risk0")]
risk1 <- results1[results1$time == 21, c("risk1", "lower_CI_risk1", "upper_CI_risk1")]
risk2 <- results2[results2$time == 21, c("risk2", "lower_CI_risk2", "upper_CI_risk2")]

print(risk0)
print(risk1)
print(risk2)

#1. extract hazard estimates and their standard errors from predict().
#2. Convert hazards to cumulative survival probabilities using the product rule.
#3. Use Greenwood’s formula to compute standard errors for survival.
#4. Convert survival estimates to risk estimates and compute confidence intervals.

##########
library(nnet)
# Fit multinomial logistic regression model
model2 <- multinom (arm ~ age_new + I(age_new*age_new) + sex + 
                 as.factor(country_income)+ comorbidities_Chalson + I(comorbidities_Chalson* comorbidities_Chalson ) + 
                 sofa_score_sum  + I(sofa_score_sum * sofa_score_sum) + 
                 infection_types + delay + I(delay*delay) + icu_at_onset4 + vent_at_onset4 + adm_ward_types_new +
                 los_onset4 + I(los_onset4* los_onset4) + mono_poly, 
               data=CR_final)
summary(model2) # 

model2a <- multinom (CR_targeted_abx_name ~ age_new + sex + 
                  as.factor(country_income)+ comorbidities_Chalson  + 
                  sofa_score_sum + infection_types,
                data=CR_final)
summary(model2a) # 

model2b <- multinom (CR_targeted_abx_name ~ age_new  + sex + 
                  as.factor(country_income)+ as.factor(country_ab) + comorbidities_Chalson  + 
                  sofa_score_sum + adm_ward_types_new + icu_at_ast + vent_at_ast + 
                  infection_types + delay  + mvdur_ast + 
                  los_ast  + iculos_ast + as.factor(org_names_all),
                data=CR_final)
summary(model2b) #

model2c <- multinom (CR_targeted_abx_name ~ age_new + as.factor(country_income) + 
                  sofa_score_sum + infection_types + as.factor(org_names_all),
                data=CR_final)
summary(model2c) #

#Use the predicted probabilities from this model to estimate nonstabilized IP weights
# Generate predicted probabilities for each treatment arm
CR_final$prob_matrix <- predict(model2, newdata = CR_final, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.
head(CR_final$prob_matrix) 

# Try without this line
# Convert treatment variable to character (to match column names in prob_matrix)
# CR_final$CR_targeted_abx_name <- as.character(CR_final$CR_targeted_abx_name)

# For each patient, select the probability corresponding to the treatment they actually received.
# Extract the probability corresponding to the observed treatment for each row
CR_final$prob <- mapply(function(i, treatment) {
  CR_final$prob_matrix[i, treatment]
}, i = seq_len(nrow(CR_final)), treatment = CR_final$arm)  #This assigns prob[i] as the probability of the observed treatment for individual i.

# Compute weight (nonstabilized IPW)
CR_final$w2 <- 1 / CR_final$prob

# Before truncation
hist(CR_final$w2)

# Truncating weights at 99th percentile
CR_final$w2[CR_final$w2 > quantile(CR_final$w2, 0.99)] <- quantile(CR_final$w2, 0.99)

# Check the distribution of the nonstabilized weights
summary(CR_final$w2)  # Check distribution of weights
sd(CR_final$w2)

# After truncation
hist(CR_final$w2[CR_final$w2 <= quantile(CR_final$w2, 0.99)])
###
# use MSM with nonstabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter error due to non-integer number of successes as a result of weights
# msm.w <- geeglm(mort_21d_ast ~ CR_targeted_abx_name, data=CR_final, weights=w2, id=recordid, family=binomial())
msm.w <- glm(mort_21d_onset4 ~ arm, data=CR_final, family=binomial(), weights=w2)
summary(msm.w)
exp(coef(msm.w))  # Convert log-odds to odds ratios
exp(confint(msm.w))  # 95% confidence interval for odds ratios

# Outputting risks
new_data <- data.frame(arm = levels(CR_final$arm))
new_data$risk <- predict(msm.w, newdata = new_data, type = "response")
new_data

# Obtain predicted risks and standard errors
predictions <- predict(msm.w, newdata = new_data, type = "response", se.fit = TRUE)

# Calculate confidence intervals (95% CI)
conf_int_lower <- predictions$fit - 1.96 * predictions$se.fit
conf_int_upper <- predictions$fit + 1.96 * predictions$se.fit

# Add these confidence intervals to the new_data dataframe
new_data$conf_int_lower <- conf_int_lower
new_data$conf_int_upper <- conf_int_upper

# Output the new data with risks and confidence intervals
new_data

# Horizontal Error Bar Plot (Whisker and Dot Plot)
ggplot(new_data, aes(y = arm, x = risk)) +
  geom_point(size = 3, color = "blue") +  # Dots for the predicted risks
  geom_errorbarh(aes(xmin = conf_int_lower, xmax = conf_int_upper), height = 0.2) +  # Horizontal error bars
  labs(title = "21-Day Mortality Risk by Treatment",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 14))

### Stop here

# Estimation of stabilized ip weights 
# Fit a multinomial logistic regression model for treatment assignment
stab <- multinom(arm~ 1, data = CR_final)
summary(stab)

# Obtain predicted probabilities (p.num) for each treatment arm
CR_final$p.num <- predict(stab, newdata = CR_final, type = "probs")

# Create stabilized weights for each observation based on the treatment assigned
# For each arm, calculate the stabilized weights based on the formula:
# sw = p(num) / p(predicted) for the arm that was actually received
CR_final$sw <- NA
CR_final$sw[CR_final$CR_targeted_abx_name == "Polymyxin monotherapy"] <- CR_final$p.num[, "Polymyxin monotherapy"] / CR_final$prob[, "Polymyxin monotherapy"]
CR_final$sw[CR_final$CR_targeted_abx_name == "Polymyxin + Carbapenem"] <- CR_final$p.num[, "Polymyxin + Carbapenem"] / CR_final$prob[, "Polymyxin + Carbapenem"]
CR_final$sw[CR_final$CR_targeted_abx_name == "Ceftazidime/avibactam"] <- CR_final$p.num[, "Ceftazidime/avibactam"] / CR_final$prob[, "Ceftazidime/avibactam"]

# Check the distribution of stabilized weights
summary(CR_final$sw)
sd(CR_final$sw)
# As we have seen, stabilizing the inverse probability weights decreases the range of the weight distribution, which in turn results in more efficient estimates from the outcome model in the parametric setting.

# Truncating weights at 99th percentile
CR_final$sw[CR_final$sw > quantile(CR_final$sw, 0.99)] <- quantile(CR_final$sw, 0.99)


# MSM with stabilized weights
#msm.sw <- geeglm(mort_21d_ast ~ CR_targeted_abx_name, data=CR_final, weights=sw,id=recordid,family=binomial())
msm.sw <- glm(mort_21d_ast ~ CR_targeted_abx_name, data=CR_final, family=binomial(), weights=sw)
summary(msm.sw)

# Outputting risks
new_data3 <- data.frame(CR_targeted_abx_name = levels(CR_final$CR_targeted_abx_name))
new_data3$risk <- predict(msm.sw, newdata=new_data3, type="response")
new_data3

#results compare when using stabilized versus nonstabilized weights are nearly the same but not identical (when unrounded), which is expected of a parametrically estimated outcome model.

# Check data
#sum(is.na(CR_final$mort_21d_ast))  # Count NAs
#sum(is.infinite(CR_final$mort_21d_ast))  # Count Inf values
#sum(is.nan(CR_final$mort_21d_ast))  # Count NaN values
#sum(is.na(CR_final$w2))
#sum(is.infinite(CR_final$w2))
#sum(is.nan(CR_final$w2))
#table(CR_final$CR_targeted_abx_name)
#sum(is.na(CR_final$recordid))
#CR_final <- CR_final %>% filter(!is.na(recordid))
#summary(CR_final$w2)
#sum(CR_final$w2 <= 0)  # Count weights that are zero or negative
#table(CR_final$mort_21d_ast)
unique(CR_final$CR_targeted_abx_name)
#quantile(CR_final$w2, probs = c(0, 0.25, 0.5, 0.75, 1))
#length(unique(CR_final$recordid)) == nrow(CR_final)

### Estimating the causal effect and constructing risk curves for a survival outcome via an IP weighted pooled logistic model
# Converting wide data to long form
CR_final$survtime <-ifelse(CR_final$mortday_ast == 0, 21, CR_final$mortday_ast)
# Expanding the dataset using the survtime variable
CR_final.surv <- uncount(CR_final, survtime,.remove = F)
# Creating variables for time
CR_final.surv <- CR_final.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()
# Creating variable for timesq
CR_final.surv$timesq <- CR_final.surv$time^2
# Creating event variable
CR_final.surv$event <- ifelse(
  CR_final.surv$time == CR_final.surv$survtime-1 & CR_final.surv$mort_21d_ast == 1, 1, 0)

# Fitting a pooled logistic regression model with time, treatment and product term between treatment and time, with nonstabilized weight
fit.pool.w2 <-  glm(event ~ 
                      CR_targeted_abx_name + time + timesq,
                      weights = w2,
                    data = CR_final.surv, 
                    family = "binomial")
summary(fit.pool.w2)
exp(fit.pool.w2$coefficients)
# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0.w2 <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                          time =seq(0,20), timesq=seq(0,20)^2)
results1.w2 <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                          time =seq(0,20), timesq=seq(0,20)^2)
results2.w2 <- data.frame(CR_targeted_abx_name = "Ceftazidime/avibactam",
                          time =seq(0,20), timesq=seq(0,20)^2)

# Obtain predicted hazards from pooled logistic regression model
results0.w2$hazard0 <- predict(fit.pool.w2, results0.w2, type="response")
results1.w2$hazard1 <- predict(fit.pool.w2, results1.w2, type="response")
results2.w2$hazard2 <- predict(fit.pool.w2, results2.w2, type="response")

# Estimate survival probabilities from hazards
results0.w2$surv0 <- cumprod(1-results0.w2$hazard0)
results1.w2$surv1 <- cumprod(1-results1.w2$hazard1)
results2.w2$surv2 <- cumprod(1-results2.w2$hazard2)

# Estimate risks from survival probabilities
results0.w2$risk0 <- 1 - results0.w2$surv0
results1.w2$risk1 <- 1 - results1.w2$surv1
results2.w2$risk2 <- 1 - results2.w2$surv2

# Risks at end of follow-up
risk0.w2 <-results0.w2[results0.w2$time==20,]$risk0
risk1.w2 <-results1.w2[results1.w2$time==20,]$risk1
risk2.w2 <-results2.w2[results2.w2$time==20,]$risk2
risk0.w2
risk1.w2
risk2.w2

# Risk ratio and risk difference
risk1.w2 - risk0.w2
risk1.w2 / risk0.w2

### Construct marginal cumulative incidence (risk) curves for all-cause mortality, by treatment group
# Combine results for each treatment group into a single dataset
results_combined <-merge(results0.w2, results1.w2, by=c("time", "timesq"))
results_combined <-merge(results_combined, results2.w2, by=c("time", "timesq"))
# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results_combined$time_updated <- results_combined$time + 1
results_combined <- results_combined %>% add_row(time_updated=0, risk0=0, risk1=0, risk2=0)%>%
  arrange(time_updated)

# Creating plot 
ggplot(results_combined, aes(x=time_updated)) + 
  geom_line(aes(y=risk0, color="Polymyxin monotherapy")) +
  geom_line(aes(y=risk1, color="Polymyxin + Carbapenem")) +
  geom_line(aes(y=risk2, color="Ceftazidime/avibactam")) +
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks=seq(0,21,1)) +
  ylab("Cumulative Incidence of mortality") +
  labs(colour = "Treatment") +
  theme_bw() + 
  theme(legend.position="bottom")

### STOP HERE
# Fit IP-weighted pooled logistic regression, with stabilized weights
fit.pool.sw <-  glm(event ~ 
                      CR_targeted_abx_name + time + timesq + CR_targeted_abx_name*time + CR_targeted_abx_name*timesq, 
                    weights = sw,
                    data = CR_final.surv, 
                    family = "binomial")
summary(fit.pool.sw)
exp(fit.pool.sw$coefficients)


fit.pool.sw.simp <-  glm(event ~ 
                           CR_targeted_abx_name + time + timesq, # without product term
                         weights = sw,
                         data = CR_final.surv, 
                         family = "binomial")
summary(fit.pool.sw.simp)
exp(fit.pool.sw.simp$coefficients)

# Cox proportional hazards model for comparison
fit.cox <- coxph(Surv(survtime, event) ~ CR_targeted_abx_name, data = CR_final.surv,
                 weights = w2)
summary(fit.cox)
exp(fit.cox$coefficients)

# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0.sw <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                          time =seq(0,20), timesq=seq(0,20)^2)
results1.sw <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                          time =seq(0,20), timesq=seq(0,20)^2)

# Obtain predicted hazards from pooled logistic regression model
results0.sw$hazard0 <- predict(fit.pool.sw, results0.sw, type="response")
results1.sw$hazard1 <- predict(fit.pool.sw, results1.sw, type="response")

# Estimate survival probabilities from hazards
results0.sw$surv0 <- cumprod(1-results0.sw$hazard0)
results1.sw$surv1 <- cumprod(1-results1.sw$hazard1)

# Estimate risks from survival probabilities
results0.sw$risk0 <- 1 - results0.sw$surv0
results1.sw$risk1 <- 1 - results1.sw$surv1

# Risks at end of follow-up
risk0.sw <-results0.sw[results0.sw$time==20,]$risk0
risk1.sw <-results1.sw[results1.sw$time==20,]$risk1
risk0.sw
risk1.sw

# Risk ratio and risk difference
risk1.sw - risk0.sw
risk1.sw / risk0.sw

### Construct marginal cumulative incidence (risk) curves for all-cause mortality, by treatment group
# Combine results for each treatment group into a single dataset
results_combined <-merge(results0.sw, results1.sw, by=c("time", "timesq"))

# Create a new "time" variable to reflect the fact that risks start at 0
# and are estimated at the end of each time interval
results_combined$time_updated <- results_combined$time + 1
results_combined <- results_combined %>% add_row(time_updated=0, risk0=0, risk1=0)%>%
  arrange(time_updated)

# Creating plot 
ggplot(results_combined, aes(x=time_updated)) + 
  geom_line(aes(y=risk0, color="Polymyxin monotherapy")) +
  geom_line(aes(y=risk1, color="Polymyxin + Carbapenem")) +
  geom_line(aes(y=risk2, color="Ceftazidime/avibactam")) +
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks=seq(0,21,1)) +
  ylab("Cumulative Incidence of mortality") +
  labs(colour = "Treatment") +
  theme_bw() + 
  theme(legend.position="bottom")


################################################
### Obtain percentile-based bootstrapped 95% CIs for each quantity ###

# Create input list of eds (eligible persons)
CR_final_ids <-data.frame(id = unique(CR_final.surv$recordid)) #resampling individual, select all person time of that recordid, resampled into bootstrap samples


# Create a function to obtain risks, RD, and RR from each bootstrap sample
risk.boot <- function(data, indices) {
  # Select individuals into each boostrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the boostrapped sample
  d <- left_join(boot.ids, CR_final.surv, by = c("id" = "recordid"))
  
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      CR_targeted_abx_name + time + timesq + CR_targeted_abx_name*time + CR_targeted_abx_name*timesq, 
                    weights = d$sw,
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  #Include all time points under each treatment level
  results0 <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("CR_targeted_abx_name", "time", "timesq")
  colnames(results1) <- c("CR_targeted_abx_name", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  #Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  
  # Estimate survival probabilityes S(t) from hazads, h(t)
  # S (t) = cumulative product of (1 - h(t))
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph$rd <- graph$risk1 - graph$risk0
  graph$rr <- graph$risk1 / graph$risk0
  return(c(graph$risk0[which(graph$time==20)], 
           graph$risk1[which(graph$time==20)], 
           graph$rd[which(graph$time==20)], 
           graph$rr[which(graph$time==20)]))
}

# Run 2 boostrap samples 
library(boot)
set.seed(123)
risk.results <- boot(data = CR_final_ids, statistic = risk.boot, R = 2)

# Print point estimates from the original data
head(risk.results$t0)


# 95% CI for risk in polmyxin monotherapy group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =1)

# 95% CI for risk in polmyxin + carbapenem group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =2)

# 95% CI for risk difference
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =3)

# 95% CI for risk ratio
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =4)


### Create a parametric cumulative incidnce (risk) plot that includes 95% CIs

# Estimate CIs for plot using bootstrapping

### Polymyxin monotherapy group

# Create a function to obtain risk in polymyxin monotherapy group at each time t from each bootstrap sample
risk.boot.0 <- function(data, indices) {
  #Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, CR_final.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      CR_targeted_abx_name + time + timesq + CR_targeted_abx_name*time + CR_targeted_abx_name*timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("CR_targeted_abx_name", "time", "timesq")
  colnames(results1) <- c("CR_targeted_abx_name", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk0)
  
}

# Run 2 boostrap samples
set.seed(123)
risk.results0 <- boot(data = CR_final_ids, statistic = risk.boot.0, R = 50)

# Combine relevant boostrapped results into a dataframe
risk.boot.results.0 <- data.frame(cbind(risk0 = risk.results0$t0,
                                        t(risk.results0$t)))
# Format boostrapped results for plotting
risk.boot.graph.0 <- data.frame (cbind(time = seq(0,20),
                                       mean.0 = risk.boot.results.0$risk0),
                                 ll.0 = (apply((risk.boot.results.0)[,-1], 1, quantile, probs = 0.025)),
                                 ul.0 = (apply((risk.boot.results.0)[,-1], 1, quantile, probs = 0.975)))                                      


### Polymyxin combination group
# Create a function to obtain risk in polymyxin combination group at each time t from each bootstrap sample
risk.boot.1 <- function(data, indices) {
  #Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, CR_final.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      CR_targeted_abx_name + time + timesq + CR_targeted_abx_name*time + CR_targeted_abx_name*timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(CR_targeted_abx_name = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(CR_targeted_abx_name = "Polymyxin + Carbapenem",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("CR_targeted_abx_name", "time", "timesq")
  colnames(results1) <- c("CR_targeted_abx_name", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk1)
  
}

# Run 2 boostrap samples
set.seed(123)
risk.results1 <- boot(data = CR_final_ids, statistic = risk.boot.1, R = 500)

# Combine relevant boostrapped results into a dataframe
risk.boot.results.1 <- data.frame(cbind(risk1 = risk.results1$t0,
                                        t(risk.results1$t)))
# Format boostrapped results for plotting
risk.boot.graph.1 <- data.frame (cbind(time = seq(0,20),
                                       mean.1 = risk.boot.results.1$risk1),
                                 ll.1 = (apply((risk.boot.results.1)[,-1], 1, quantile, probs = 0.025)),
                                 ul.1 = (apply((risk.boot.results.1)[,-1], 1, quantile, probs = 0.975)))

# Prepare data
risk.boot.graph.pred <- merge(risk.boot.graph.0, risk.boot.graph.1, 
                              by = "time")

# Edit data frame to reflect that risks are estimated at the END of each interval
risk.boot.graph.pred$time_0 <- risk.boot.graph.pred$time + 1
zero <-data.frame(cbind(0,0,0,0,0,0,0,0))
zero <-setNames(zero, names(risk.boot.graph.pred))
risk.boot.graph <- rbind(zero, risk.boot.graph.pred)

# Create plot
plot.plr.ci <- ggplot(risk.boot.graph,
                      aes(x=time_0)) +
  geom_line(aes(y =mean.1,
                color = "Polymyxin + Carbapenem"),
            size =1.5) +
  geom_ribbon(aes(ymin = ll.1, ymax = ul.1,
                  fill = "Polymyxin + Carbapenem"), alpha = 0.4) +
  geom_line(aes(y =mean.0,
                color = "Polymyxin monotherapy"),
            size =1.5) +
  geom_ribbon(aes(ymin = ll.0, ymax = ul.0,
                  fill = "Polymyxin monotherapy"), alpha = 0.4) + 
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks=seq(0,21,1)) +
  ylab("Cumulative Incidence of mortality") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  font("xlab",size=14)+
  font("ylab",size=14)+
  font("legend.text",size=10)+
  scale_color_manual(values = c("Polymyxin monotherapy" = "blue",
                                "Polymyxin + Carbapenem" = "red")) +
  scale_fill_manual(values = c("Polymyxin monotherapy" = "blue",
                               "Polymyxin + Carbapenem" = "red"))
# Plot
plot.plr.ci


test <- risk.boot.graph %>%
  mutate(valid_CI = ll.1 <= mean.1 & mean.1 <= ul.1) %>%
  filter(!valid_CI)

###########################################################################################
### Covariate balance plot
# matching
if (!require("MatchIt")) install.packages("MatchIt")
library(MatchIt)
if (!require("cobalt")) install.packages("cobalt")
library(cobalt)
install.packages("Hmisc")  # If not already installed
library(Hmisc)  # Load the package to access wtd.var()

# Use the MatchIt package to perform matching
match_data <- matchit(formula = CR_targeted_abx_name ~ survtime +
                        age_new + sex + country_income + comorbidities_Chalson +
                        sofa_score_sum + infection_types + org_names_all + delay + los_ast,
                      data = CR_final,
                      method = "nearest")
# The above specifications ensure that matchit() performs exact matching in a 1:1 ratio
match_data

# Generate matched dataset (this will only contain one row per person, corresponding to time 0).
CR_final_matched <- match.data(match_data)

# Check number of treated and untreated prior to matching
table(CR_final$CR_targeted_abx_name)

# Check number of treated and untreated in matched sample
table(CR_final_matched$CR_targeted_abx_name)

# Restrict long-formatted dataset to subset of individuals who were matched
ids <- CR_final_matched$recordid
final_matched <- CR_final[CR_final$recordid %in% ids,]

### Construct a covariate balance plot ###
# Create subsets of data, according to CR_targeted_abx_name
CR_final_matched_0 <- CR_final_matched[CR_final_matched$CR_targeted_abx_name == "Polymyxin monotherapy",]
CR_final_matched_1 <- CR_final_matched[CR_final_matched$CR_targeted_abx_name == "Polymyxin + Carbapenem",]

# List variables to include in the plot
varlist <- c("age_new", "sex", "country_income", "comorbidities_Chalson", "sofa_score_sum", "infection_types", "org_names_all", "mono_poly", "delay", "los_ast")

# Create function to calculate mean difference, or standardized mean difference for age
meanfctn <- function(x){
  if(x == "age_new"){
    t0 <- CR_final_matched_0[[x]]
    t1 <- CR_final_matched_1[[x]]
    md <- (mean(t1) - mean(t0)) / sd(t1)
  } else {
    t0 <- CR_final_matched_0[[x]]
    t1 <- CR_final_matched_1[[x]]
    md <- mean(t1) - mean(t0)}
  return(c(var =x, md=md))
}

cov_plot_vars <- lapply(varlist, meanfctn) %>% do.call(rbind,.) %>% as.data.frame() 
cov_plot_vars$md <- as.numeric(cov_plot_vars$md)

# Plot mean differences
cov_plot <- ggplot(data = cov_plot_vars) + 
  geom_point(aes(x = md, y = var), color = "steelblue") + 
  geom_vline(xintercept = 0) +
  xlim(-1, 1) +
  labs( y = "Covariates", x = "Mean difference", title = "Covariate balance plot in matched population")  

cov_plot

# Generate covariate balance plot with *standardized* mean differences 
love.plot(match_data, 
          binary = "std", 
          title = "Covariate balance plot",
          sample.names = c("Unmatched", "Matched"),
          addl = ~ age_new + sex + country_income + comorbidities_Chalson + sofa_score_sum + infection_types + org_names_all + mono_poly + delay + los_ast,
          drop.distance = TRUE,
          grid = TRUE)


### Construct a covariate balance plot for the weighted population ###
# Create subsets of data, according to CR_targeted_abx_name
CR_final_weighted_0 <- CR_final[CR_final$CR_targeted_abx_name == "Polymyxin monotherapy",]
CR_final_weighted_1 <- CR_final[CR_final$CR_targeted_abx_name == "Polymyxin + Carbapenem",]
CR_final_weighted_2 <- CR_final[CR_final$CR_targeted_abx_name == "Ceftazidime/avibactam",]

# List variables to include in the plot
varlist <- c("age_new", "sex", "country_income", "comorbidities_Chalson", "sofa_score_sum", "infection_types", "org_names_all", "mono_poly", "delay", "los_ast")

# Create function to calculate mean difference, or standardized mean difference for age
meanfctn <- function(x){
  if(x == "age_new"){
    t0 <-CR_final_weighted_0[[x]]
    t1 <- CR_final_weighted_1[[x]]
    t2 <- CR_final_weighted_2[[x]]
    md <- (mean(t1) - mean(t0)) / sd(t1)
  } else {
    t0 <- CR_final_weighted_0[[x]]
    t1 <- CR_final_weighted_1[[x]]
    t2 <- CR_final_weighted_2[[x]]
    md <- mean(t1) - mean(t0)}
  return(c(var =x, md=md))
}

# Calculate mean differences for covariates (SMD for age)
wmean_fctn <- function(x){
  if(x == "age_new"){
    md <- (weighted.mean(CR_final_weighted_1[[x]], CR_final_weighted_1$sw) - 
             weighted.mean(CR_final_weighted_0[[x]], CR_final_weighted_0$sw)) / 
      sqrt(wtd.var(CR_final_weighted_1[[x]], CR_final_weighted_1$sw))  # Corrected weights
  } else {
    t0 <- weighted.mean(CR_final_weighted_0[[x]], CR_final_weighted_0$sw)
    t1 <- weighted.mean(CR_final_weighted_1[[x]], CR_final_weighted_1$sw)
    md <- t1 - t0
  }
  return(c(var = x, md = md))
}

# Create the covariate plot
covplot_w <- lapply(varlist, wmean_fctn) %>% do.call(rbind,.) %>% as.data.frame()
covplot_w$md <- as.numeric(covplot_w$md)

# Plot it
covplot_weighted <- ggplot(data = covplot_w) +
  geom_point(aes(x = md, y = var), color = "steelblue") + scale_x_continuous(limits = c(-0.1, 0.1)) +
  geom_vline(xintercept = 0) +
  labs(y = "Covariates", x = "Mean Difference", title = "Covariate Balance Plot")

covplot_weighted


#######################################################################
fit.pool.sw.simp <-  glm(event ~ 
                           CR_targeted_abx_name + time + timesq, # without product term
                         weights = sw,
                         data = CR_final.surv, 
                         family = "binomial")
summary(fit.pool.sw.simp)
exp(fit.pool.sw.simp$coefficients)

# Cox proportional hazards model for comparison
fit.cox <- coxph(Surv(survtime, event) ~ CR_targeted_abx_name, data = CR_final.surv,
                 weights = sw)
summary(fit.cox)
exp(fit.cox$coefficients)


coef_summary <- summary(fit.pool.sw.simp)
beta <- coef_summary$coefficients[, "Estimate"]
se <- coef_summary$coefficients[, "Std. Error"]

# Extract coefficients and covariance matrix
beta <- coef(fit.pool.sw.simp)
vcov_mat <- vcov(fit.pool.sw.simp)

# Linear predictors
eta_arm1 <- beta[1]
eta_arm2 <- beta[1] + beta[2]
eta_arm3 <- beta[1] + beta[3]

# Standard errors of linear predictors
se_eta_arm1 <- sqrt(vcov_mat[1, 1])
se_eta_arm2 <- sqrt(vcov_mat[1, 1] + vcov_mat[2, 2] + 2 * vcov_mat[1, 2])
se_eta_arm3 <- sqrt(vcov_mat[1, 1] + vcov_mat[3, 3] + 2 * vcov_mat[1, 3])

# Confidence intervals on the log-odds scale
ci_eta_arm1 <- c(
  eta_arm1 - 1.96 * se_eta_arm1,
  eta_arm1 + 1.96 * se_eta_arm1
)
ci_eta_arm2 <- c(
  eta_arm2 - 1.96 * se_eta_arm2,
  eta_arm2 + 1.96 * se_eta_arm2
)
ci_eta_arm3 <- c(
  eta_arm3 - 1.96 * se_eta_arm3,
  eta_arm3 + 1.96 * se_eta_arm3
)

# Transform to probability scale
risk_arm1 <- 1 / (1 + exp(-eta_arm1))
risk_arm2 <- 1 / (1 + exp(-eta_arm2))
risk_arm3 <- 1 / (1 + exp(-eta_arm3))

ci_risk_arm1 <- 1 / (1 + exp(-ci_eta_arm1))
ci_risk_arm2 <- 1 / (1 + exp(-ci_eta_arm2))
ci_risk_arm3 <- 1 / (1 + exp(-ci_eta_arm3))

# Display risks and confidence intervals
list(
  risk_arm1 = c(risk_arm1, ci_risk_arm1),
  risk_arm2 = c(risk_arm2, ci_risk_arm2),
  risk_arm3 = c(risk_arm3, ci_risk_arm3)
)

# Create a data frame for plotting
plot_data <- data.frame(
  Arm = c("Carbapenem + Polymyxin", "Polymxin monotherapy", "Ceftazidime avibactam monotherapy"),
  Risk = c(0.4053872, 0.4319106, 0.3505323),
  CI_Low = c(0.3308195, 0.3185274, 0.2439864),
  CI_High = c(0.4845912, 0.5529080, 0.4744089)
)

ggplot(plot_data, aes(x = Risk, y = Arm)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = CI_Low, xmax = CI_High), height = 0.3) +
  xlab("21-day risk of mortality") +
  ylab("Treatment Arm") +
  ggtitle("Carbapenem-resistant enterobacterales") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.1)) +
  theme_minimal(base_size = 14)

####################################################################################
skim(CR_sofa_all)
# Check data
sofa_date <- sofa_all %>% select (recordid, sofa_res_date, sofa_co_date, sofa_liv_date, sofa_ren)
checkdate <- merge (CR, sofa_date, by = "recordid", all.x = TRUE) %>%
  select(recordid, 
         #infection_types, inf_onset,  spec_day,
         spec_date, ast_date, sofa_res_date, sofa_co_date, sofa_liv_date, sofa_ren) %>%
  filter (recordid %in% CR_final$recordid) %>%
  distinct()

checkast <-ast_all %>% select(recordid,spec_date,org, ast_date, org_names, org_names_all, ris_Carbapenems)
checkastindex <-ast_all_index %>% select(recordid,infection_types,spec_date,org_names_all, org_combined) 



##########
############################ 
## TARGET TRIAL EMULATION ##
############################
# cci
cmb <- baseline_outcomes_index %>% select (recordid, hpd_admreason, cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd, cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab, cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___hivwa, cmb_comorbidities___hivna, cmb_comorbidities___mlr, cmb_comorbidities___mal, cmb_comorbidities___mst, cmb_comorbidities___mld, cmb_comorbidities___liv, cmb_comorbidities___pep, cmb_comorbidities___renal, cmb_comorbidities___tub, cmb_comorbidities___oth)%>%
  filter(recordid %in% CR_final$recordid)

# BSI only-immunosuppresion status and Pitts
comorb <- f07a %>% select (recordid, inf_onset, hiv, end, ins, am, cc, pt, child_c, neu,hae, sot,
                           temp, temp_1, res_rate, res_rate_1, heart_rate, heart_rate_1, sys_bp, sys_bp_1, men_sta, 
                           acute, iv_agents, mv,ca)%>%
  filter(recordid %in% CR_final$recordid) 

# qsofa
qsofa <- severity_new %>% filter(recordid %in% CR_final$recordid) %>% mutate(qsofa = severity_score)%>% select (-severity_score, -severity_score_scale, -adm_ward_types_ori, -adm_ward_types_new)
qsofa_detailed <- severity %>% filter(recordid %in% CR_final$recordid) %>% 
  select (recordid,redcap_event_name, ser_gcs_under15, ser_rr_22up, ser_sbp_under100, ser_abnormal_temp_adult, hai_icu48days, hai_surg,mic_bloodcollect,mic_rec_antibiotic, 
          ser_gcs_under15_new,ser_gcs_under15_new,ser_rr_22up_new, ser_sbp_under100_new)

# SOFA
colnames(sofa_all) <- make.names(colnames(sofa_all), unique = TRUE)

CR_sofa_all <- sofa_all %>% filter(recordid %in% CR_final$recordid)
CR_sofa_new <- sofa_new %>% filter(recordid %in% CR_final$recordid)
median_sofa <- median(CR_sofa_new$sofa_score, na.rm = TRUE)

# Ignore NA and Sum Available Scores (the total SOFA score for that patient may be underestimated, but it might still be valid depending on your study goal, This method is acceptable if the missingness is random and the missing values are not due to a systematic issue.)
CR_sofa_new$sofa_score_sum <- rowSums(CR_sofa_new[, c("score_respiration", "score_coagulation", 
                                                      "score_liver", "score_CVS", 
                                                      "score_CNV", "score_renal")], na.rm = TRUE)
CR_sofa_sum <- CR_sofa_new %>% select(recordid, sofa_score_sum)
# Impute missing values in each SOFA sub-score (e.g., score_respiration, score_coagulation, etc.) using mean or median:
CR_sofa_median_imputed <- CR_sofa_new %>% select(-sofa_score, -sofa_score_sum) %>%
  mutate(across(starts_with("score_"), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  rowwise() %>%
  mutate(sofa_score_median = sum(c_across(starts_with("score_"))))
CR_sofa_median <- CR_sofa_median_imputed %>% select(recordid, sofa_score_median)

pathogen <- CR_final %>% select(recordid, spec_date, org_combined, mono_poly) %>% distinct() 
pathogen2 <- checkast %>% select(recordid, spec_date, org_names) %>% distinct()
pathogen <- pathogen %>% left_join(pathogen2, by = c("recordid", "spec_date")) %>% distinct()
# For recordid VN-004-A0-0083, assign org_names= org_combined
pathogen[pathogen$recordid == "VN-004-A0-0083", "org_names"] <- "K. pneumoniae"
pathogen[pathogen$recordid == "TH-001-A0-0130", "org_names"] <- "Pseudomonas"

cre <- sofa_all %>% select (recordid, sofa_ren,sofa_ren_2, sofa_cre_1, sofa_cre_2)%>%
  filter (recordid %in% CR_final$recordid)  #only 23 missing out of 310

# Compute the median of non-NA values in sofa_cre_2 and sofa_cre_1 * 88.4
median_crea <- cre %>%
  mutate(temp_crea = ifelse(is.na(sofa_cre_2), sofa_cre_1 * 88.4, sofa_cre_2)) %>%
  pull(temp_crea) %>%
  median(na.rm = TRUE)

# Mutate new column 'crea' based on conditions
cre <- cre %>%
  mutate(crea = case_when(
    !is.na(sofa_cre_2) ~ sofa_cre_2,                   # If sofa_cre_2 is not NA, use it
    !is.na(sofa_cre_1) ~ sofa_cre_1 * 88.4,            # Else, use sofa_cre_1 * 88.4
    TRUE ~ median_crea                                 # Else, use median imputation
  ))

# For recordid PH-002-A0-0036, assugn crea as median_crea
cre[cre$recordid == "PH-002-A0-0036", "crea"] <- median_crea
cre[cre$recordid == "CN-002-A0-0026", "crea"] <- 119

cre <- cre %>%  select(recordid, crea)

CR_final <- CR_final %>%
  left_join(CR_sofa_sum, by = "recordid") %>%
  left_join(CR_sofa_median, by = "recordid") %>%
  left_join(qsofa, by = "recordid") %>%
  left_join(cmb, by = "recordid")%>%
  left_join(cre, by = "recordid") %>%
  #for country_income2, change to 2 groups upper middle income and high income, low income and low middle income
  mutate(country_income2 = ifelse(country_income == "Low income" | country_income == "Lower middle income", "Low income and Lower middle income", "Upper middle income and High income")) %>%
  relocate(country_income2, .after = country_income) 

# Convert columns to appropriate classes
CR_final <- CR_final %>%
  mutate(delay = ifelse(delay < 0, 0, delay))%>%
  mutate(delay_group = ifelse(delay == 0, 0, 1)) %>%
  mutate(
    # Continuous variables as numeric
    across(c(age_new, sofa_score, comorbidities_Chalson,  
             los_onset4, iculos_onset4, mvdur_onset4, mortday_onset4, delay, 
             sofa_score_sum, sofa_score_median, fbis_score,fup_day_onset4, qsofa,crea), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset4, vent_at_onset4, 
             mort_21d_onset4, mono_poly, UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, 
             delay_group, country_income2), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined, arm, organism,ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

# Sensitivity analysis
# Monomicrobial only
CR_final <- CR_final %>% filter(mono_poly == "monomicrobial")
CR_final <- CR_final %>% filter(infection_types == "VAP")
CR_final <- CR_final %>% filter(infection_types == "VAP")

# Fit multinomial logistic regression model
model2 <- multinom (arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      sofa_score_sum  + infection_types + 
                       icu_at_onset4 + vent_at_onset4+
                      # icu_at_onset5 + vent_at_onset5+
                     + mono_poly,
                    data=CR_final)
summary(model2) 

# Only adjust for age
model2h <- multinom (arm ~ age_new, 
                    data=CR_final)


# Drop charlson
model2g <- multinom (arm ~ age_new + sex + country_income2 + 
                      sofa_score_sum  + infection_types + icu_at_onset4 + vent_at_onset4 + 
                      mono_poly, 
                    data=CR_final)

# using iculos_onset4, mvdur_onset4"  
model2f <- multinom (arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      sofa_score_sum  + infection_types + icu_at_onset4 + vent_at_onset4 + 
                      mono_poly + iculos_onset4 + mvdur_onset4,
                    data=CR_final)

# Using selected comorbidities
model2e <- multinom (arm ~ age_new + sex + country_income2 + 
                       cmb_comorbidities___cpd + 
                       cmb_comorbidities___cog +
                       cmb_comorbidities___diab + 
                       cmb_comorbidities___diad +
                       cmb_comorbidities___onc +
                       cmb_comorbidities___mst +  
                       cmb_comorbidities___liv +
                       cmb_comorbidities___renal + 
                      sofa_score_sum  + infection_types + icu_at_onset4 + vent_at_onset4 + 
                      mono_poly, 
                    data=CR_final)

# Using qsofa
model2d <- multinom (arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      qsofa  + infection_types + icu_at_onset4 + vent_at_onset4 + 
                      mono_poly, 
                    data=CR_final)


# add variables
model2c <- multinom (arm ~ age_new + sex +  country_income+ comorbidities_Chalson + 
                      sofa_score_sum  + infection_types + icu_at_onset4 + vent_at_onset4 + 
                       + delay_group  + los_onset4 + mono_poly + crea + adm_ward_types_new, 
                    data=CR_final)

# using sofa median
model2b <- multinom (arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      sofa_score_median + infection_types + icu_at_onset4 + vent_at_onset4 + 
                      mono_poly, 
                    data=CR_final)
# using quadratic
model2a <- multinom (arm ~ age_new + I(age_new*age_new) + sex + country_income2 + 
                       comorbidities_Chalson + I(comorbidities_Chalson*comorbidities_Chalson) +
                      sofa_score_sum + I(sofa_score_sum*sofa_score_sum) + 
                        infection_types + icu_at_onset4 + vent_at_onset4 + mono_poly, 
                    data=CR_final)
# test model
model2 <- multinom (arm ~ age_new + I(age_new*age_new) + sex + 
                      as.factor(country_income)+ comorbidities_Chalson + I(comorbidities_Chalson* comorbidities_Chalson) + 
                      sofa_score_sum  + I(sofa_score_sum * sofa_score_sum) + 
                      infection_types + delay + I(delay*delay) + icu_at_onset4 + vent_at_onset4 + adm_ward_types_new +
                      los_onset4 + I(los_onset4* los_onset4) + mono_poly + cre + I(crea*crea), 
                    data=CR_final)

#Use the predicted probabilities from this model to estimate nonstabilized IP weights
# Generate predicted probabilities for each treatment arm
CR_final$prob_matrix <- predict(model2, newdata = CR_final, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.

# Extract the probability corresponding to the observed treatment 
CR_final$prob <- mapply(function(i, treatment) {
  CR_final$prob_matrix[i, treatment]
}, i = seq_len(nrow(CR_final)), treatment = CR_final$arm)  # assign prob[i] as the probability of the observed treatment for individual i.

# Compute weight (nonstabilized IPW)
CR_final$w2 <- 1 / CR_final$prob
# Truncating weights at 99th percentile
CR_final$w2[CR_final$w2 > quantile(CR_final$w2, 0.99)] <- quantile(CR_final$w2, 0.99)
# Check the distribution of the nonstabilized weights
summary(CR_final$w2) 
sd(CR_final$w2)

###
# use MSM with nonstabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter error due to non-integer number of successes as a result of weights
# msm.w <- geeglm(mort_21d_onset5 ~ arm, data=CR_final, weights=w2, id=recordid, family=binomial())
# msm.w <- glm(mort_21d_onset5 ~ arm, data=CR_final, family=binomial(), weights=w2)
msm.w <- glm(mort_21d_onset4 ~ arm, data=CR_final, family=binomial(), weights=w2)
summary(msm.w)
#exp(coef(msm.w))  # Convert log-odds to odds ratios
#exp(confint(msm.w))  # 95% confidence interval for odds ratios

# Outputting risks
new_data <- data.frame(arm = levels(CR_final$arm))
new_data$risk <- predict(msm.w, newdata = new_data, type = "response")
new_data

# Obtain predicted risks and standard errors
predictions <- predict(msm.w, newdata = new_data, type = "response", se.fit = TRUE)

# Calculate confidence intervals (95% CI)
conf_int_lower <- predictions$fit - 1.96 * predictions$se.fit
conf_int_upper <- predictions$fit + 1.96 * predictions$se.fit

# Add these confidence intervals to the new_data dataframe
new_data$conf_int_lower <- conf_int_lower
new_data$conf_int_upper <- conf_int_upper

# Output the new data with risks and confidence intervals
new_data

# Horizontal Error Bar Plot (Whisker and Dot Plot)
ggplot(new_data, aes(y = arm, x = risk)) +
  geom_point(size = 3, color = "blue") +  # Dots for the predicted risks
  geom_errorbarh(aes(xmin = conf_int_lower, xmax = conf_int_upper), height = 0.2) +  # Horizontal error bars
  labs(title = "21-Day Mortality Risk by Treatment",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 14))

CR_abx <- CR_names_4 %>% filter (recordid %in% CR_final$recordid) %>% select(recordid, arm, anti_onset4) %>% distinct()
