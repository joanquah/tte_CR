#### Sensitivity analysis - day 5####
checkast <-ast_all%>% select(recordid,spec_date, org_names_all, pathogen_group,ris_Carbapenems, Meropenem)

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
         #ast_date, 
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
  mutate(onset4 = inf_onset + 4)%>% relocate(onset4, .after = spec_date) %>%
  mutate(onset5 = inf_onset + 5) %>% relocate(onset5, .after = onset4) %>%
  select(recordid, organism, adm_date, inf_onset, spec_date, onset4, onset5, anti_names, anti_new, anti_start, anti_end, susceptibility)%>%
  mutate (delay = as.numeric(difftime(anti_start, inf_onset, units = "days")))%>%
  filter(delay >= -5 | anti_new == "Carbapenem") %>%
  # filter(!recordid %in% recordid_to_exclude$recordid) %>%
  filter(delay <= 5) %>%
  mutate(dur = as.integer(difftime(anti_end, anti_start, units = "days"))) %>%
  mutate (dur = dur+1) %>%
  filter(dur>2) %>%
  # filter(onset4 >= anti_start & onset4 < anti_end) %>%
  filter(onset5 >= anti_start & onset5 < anti_end) %>%
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
  mutate(spec_day = as.integer(difftime(spec_date, inf_onset, units = "days")))%>%
  relocate(spec_day, .after = spec_date) %>%
  distinct() %>%
  filter(spec_day<6)  #903 qualified ITT  #849 initiated studied drug within 5 days from infection onset

length(unique(CR_names$recordid)) #1290. For onset 4 988 for ITT, 976 for PP; For onset5 94 for ITT, 941 for PP   # 850
#CR_names <- CR_names %>% select (recordid, organism, susceptibility, anti_ast, anti_onset4,anti_onset5) %>% distinct

CR_names_5 <- CR_names %>%
  mutate(arm = case_when(
    grepl("Ceftazidime/avibactam", anti_onset5) ~ "Ceftazidime/avibactam based",
    anti_onset5 %in% c("Polymyxin B", "Colistin") ~ "Polymyxin monotherapy",
    anti_onset5 %in% c("Cefoperazone/sulbactam", "Ampicillin/sulbactam") ~ "Sulbactam based",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset5) ~ "Polymyxin Sulbactam",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset5) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset5) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset5) ~ "Sulbactam based",
    TRUE ~ NA_character_
  )) %>%
  mutate(arm2 = case_when(
    anti_onset5 %in% c("Polymyxin B", "Colistin") ~ "Polymyxin monotherapy",
    anti_onset5 %in% c("Cefoperazone/sulbactam", "Ampicillin/sulbactam") ~ "Sulbactam monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset5) ~ "Polymyxin Sulbactam",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset5) ~ "Polymyxin combination",
    grepl("Cefoperazone/sulbactam|Ampicillin/sulbactam", anti_onset5) & grepl("Meropenem|Imipenem|Doripenem|Tigecycline|Fosfomycin", anti_onset5) ~ "Sulbactam combination",
    TRUE ~ NA_character_
  )) %>%
  filter(grepl("Colistin|Polymyxin B|Ampicillin/sulbactam|Cefoperazone/sulbactam|Ceftazidime/avibactam", anti_onset5))%>%
  mutate(arm = ifelse(anti_onset5 == "Polymyxin B + Colistin", "Polymyxin monotherapy", arm)) %>%
  mutate(arm = ifelse(anti_onset5 == "Cefoperazone/sulbactam + Ampicillin/sulbactam", "Sulbactam based", arm)) %>%
  mutate(arm2 = ifelse(anti_onset5 == "Cefoperazone/sulbactam + Ampicillin/sulbactam", "Sulbactam monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset5 == "Cefoperazone/sulbactam + Neb Col", "Sulbactam based", arm)) %>%
  mutate(arm2 = ifelse(anti_onset5 == "Cefoperazone/sulbactam + Neb Col", "Sulbactam monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset5 == "Colistin + Neb Col", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset5 == "Colistin + Neb Col", "Polymyxin monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset5 == "Polymyxin B + Neb Col", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset5 == "Polymyxin B + Neb Col", "Polymyxin monotherapy", arm2)) %>%
  mutate(arm = ifelse(anti_onset5 == "Neb Col + Polymyxin B", "Polymyxin monotherapy", arm)) %>%
  mutate(arm2 = ifelse(anti_onset5 == "Neb Col + Polymyxin B", "Polymyxin monotherapy", arm2)) %>%
  
  mutate(arm3 = case_when(
    grepl("Ceftazidime/avibactam", anti_onset5) ~ "Ceftazidime/avibactam based",
    anti_onset5 %in% c("Polymyxin B", "Colistin","Polymyxin B + Colistin", 
                       "Polymyxin B + Ampicillin/sulbactam", "Ampicillin/sulbactam + Polymyxin B",
                       "Ampicillin/sulbactam + Colistin","Colistin + Ampicillin/sulbactam",
                       "Colistin + Cefoperazone/sulbactam","Cefoperazone/sulbactam + Colistin",
                       "Colistin + Ampicillin/sulbactam + Cefoperazone/sulbactam") ~ "Polymyxin monotherapy",
    grepl("Polymyxin B|Colistin", anti_onset5) & grepl("Meropenem|Imipenem|Tigecycline|Fosfomycin", anti_onset5) ~ "Polymyxin combination",
    TRUE ~ NA_character_
  )) %>%
  distinct() 

table(CR_names_4$arm) # 90 ceftazidime/avibactam,  125 polymyxin monotherapy,  333 polymyxin combination
length(unique(CR_names_4$recordid)) #316 #336  #310  #982  #753  #723 had studied drug on index date

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

CR_combined <- left_join(CR_names_5, icu_mv, by = "recordid") %>% distinct() %>%
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
    # Ventilator status at onset5
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
    # Length of stay before onset5
    los_onset5 = as.integer(difftime(onset5, adm_date, units = "days")),
    # Sum up number of days in ICU before onset5
    iculos_onset5 = rowSums(
      cbind(
        pmax(pmin(onset5, icu_hd_ap_1_3, na.rm = TRUE) - icu_hd_ap_1_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_2_3, na.rm = TRUE) - icu_hd_ap_2_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_3_3, na.rm = TRUE) - icu_hd_ap_3_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_4_3, na.rm = TRUE) - icu_hd_ap_4_2, 0, na.rm = TRUE),
        pmax(pmin(onset5, icu_hd_ap_5_3, na.rm = TRUE) - icu_hd_ap_5_2, 0, na.rm = TRUE)
      ), na.rm = TRUE
    ),
    # Mechanical ventilation duration before onset5
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
  select(recordid, inf_onset, spec_date, onset5, delay, arm, arm2, arm3, 
         anti_onset5, anti_new, anti_start, organism, icu_at_onset5, vent_at_onset5, 
         los_onset5, iculos_onset5, mvdur_onset5, susceptibility, organism) %>%
  filter(anti_new != "Carbapenem") %>%  
  group_by(recordid) %>%
  filter(delay == min(delay)) %>%  
  ungroup()%>%
  select (recordid, inf_onset, onset5, delay, arm,  arm2, arm3,
          anti_onset5, organism,  icu_at_onset5, vent_at_onset5, los_onset5, iculos_onset5, mvdur_onset5) %>%
  distinct()

# join CR_targeted and outcome dataframe, follow CR_targeted recordid
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
  left_join(CR_final %>% select(recordid, onset5), by = "recordid")%>%
  filter(spec_date <= onset5) %>%
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

CR_final_drop <- CR_final %>% filter(!grepl("Tigecycline|Fosfomycin", anti_onset5))
length(unique(CR_final_drop$recordid)) #638 wihtout tigecycline and fosfomycin  #638
CR_final <- CR_final_drop

#######################################################################################################################################
################# STOP here ####################
#######################################################################################################################################

colnames(CR_final)

CRAB_drop <- CR_final_drop %>% filter(organism == "CRAB")
length(unique(CRAB_drop$recordid)) #429  
nrow(CRAB_drop) 
CRE_drop <- CR_final_drop %>% filter(organism == "CRE")
length(unique(CRE_drop$recordid)) #187 
nrow(CRE_drop) 
CRPAE_drop <- CR_final_drop %>% filter(organism == "CRPAE")
length(unique(CRPAE_drop$recordid)) #65
nrow(CRPAE_drop) 

CRAB <- CR_final %>% filter(cr_aci==1)%>%
   #remove arm="Ceftazidime/avibactam based"
  filter(!(arm=="Ceftazidime/avibactam based"))
#%>% filter(monopoly == "monomicrobial")
length(unique(CRAB$recordid)) #487 (370 mono)   #445 (no tige)
nrow(CRAB) #487   #445
#429 ITT 425 PP

CRE_CRPAE <- CR_final %>% filter(cr_pae==1 | cr_ent==1)
# %>% filter(monopoly == "monomicrobial")
length(unique(CRE_CRPAE$recordid)) #316 (221 mono)  #264 (no tige)
nrow(CRE_CRPAE) #316
# 246

CR_final_mono <- CR_final %>% filter(monopoly == "monomicrobial") 
length(unique(CR_final_mono$recordid)) #570
nrow(CR_final_mono) 

#skim(CRAB)
#skim(CRE_CRPAE)

###########################################################################################################################################
check <- CRE_CRPAE %>% select(recordid, onset5, org_combined_new, anti_onset5, arm, arm2, arm3, aci, ent, pae, cr_aci, cr_pae, cr_ent)%>%
  filter(!is.na(arm3))%>%
  filter(!grepl("Tigecycline|Fosfomycin", anti_onset5))
length(unique(check$recordid)) #261 (with tige), #211 (no tige)

crab <- CRAB %>% filter(!is.na(arm2))
#%>% filter(!grepl("Tigecycline|Fosfomycin", anti_onset5))
length(unique(crab$recordid)) #443 (no tige)  vs 485 (with tige)
#427

crab %>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRAB") 

crab_mono <- CRAB %>% filter(!is.na(arm2))%>% filter(monopoly =="monomicrobial")
length(unique(crab_mono$recordid)) #345

cre_crpae <- CRE_CRPAE %>% filter(!is.na(arm3)) 
# %>%
# filter(!grepl("Tigecycline|Fosfomycin", anti_onset5))
length(unique(cre_crpae$recordid)) #261 (with tige), #211 (no tige)  #246
# 196

cre_crpae %>% count(arm3) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRE and CRPAE")

####
  
  
  # cci
  cmb <- baseline_outcomes_index %>% select (recordid, hpd_admreason, cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd, cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab, cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___hivwa, cmb_comorbidities___hivna, cmb_comorbidities___mlr, cmb_comorbidities___mal, cmb_comorbidities___mst, cmb_comorbidities___mld, cmb_comorbidities___liv, cmb_comorbidities___pep, cmb_comorbidities___renal, cmb_comorbidities___tub, cmb_comorbidities___oth,
  )%>%
  filter(recordid %in% CR_final$recordid)

# BSI only-immunosuppresion status and Pitts
comorb <- f07a %>% select (recordid, inf_onset, hiv, end, ins, am, cc, pt, child_c, neu,hae, sot,
                           temp, temp_1, res_rate, res_rate_1, heart_rate, heart_rate_1, sys_bp, sys_bp_1, men_sta, 
                           acute, iv_agents, mv,ca)%>%
  filter(recordid %in% CR_final$recordid) 

# qsofa
qsofa <- severity_new %>% filter(recordid %in% CR_final$recordid) %>% mutate(qsofa = severity_score)%>% select (-severity_score, -severity_score_scale, -adm_ward_types_ori, -adm_ward_types_new)
qsofa_detailed <- severity %>% filter(recordid %in% CR_final$recordid) %>% 
  select (recordid,redcap_event_name, ser_gcs_under15, ser_rr_22up, ser_sbp_under100, ser_abnormal_temp_adult, hai_icu48days, hai_surg, , hai_have_med_device___vent,mic_bloodcollect,mic_rec_antibiotic, 
          ser_gcs_under15_new,ser_gcs_under15_new,ser_rr_22up_new, ser_sbp_under100_new)

CR_icu_ven <- qsofa_detailed %>% select(recordid, hai_icu48days,hai_surg, hai_have_med_device___vent)%>% filter(recordid %in% CR_final$recordid)
CR_icu <- f02 %>% select (recordid, f02_deleted, f02_infected_episode_complete, hai_icu48days, hai_have_med_device___vent)%>%
  #filter(recordid %in% CR_final$recordid)%>%
  filter(f02_infected_episode_complete == 2)%>%
  # filter rows when f02_deleted is not Y
  select(-f02_infected_episode_complete, -f02_deleted) %>%
  group_by(recordid) %>%
  summarise(
    hai_icu48days = if_else(any(hai_icu48days == "Y"), "Y", "N"),
    hai_have_med_device___vent = if_else(any(hai_have_med_device___vent == 1), 1L, 0L)
  )%>%
  distinct()%>%
  mutate(hai_icu48days = ifelse(hai_icu48days == "Y", 1, 0)) 

cre <- sofa_all %>% select (recordid, sofa_ren,sofa_ren_2, sofa_cre_1, sofa_cre_2)%>%
  filter (recordid %in% CR_final$recordid)  #only 23 missing out of 310
# Combine crea1 and crea2
cre <- cre %>%
  mutate(crea = case_when(
    !is.na(sofa_cre_2) ~ sofa_cre_2,                   # If sofa_cre_2 is not NA, use it
    !is.na(sofa_cre_1) ~ sofa_cre_1 * 88.4,            # Else, use sofa_cre_1 * 88.4
    TRUE ~ NA_real_                           # Else, use NA
  ))%>%
  select(recordid, crea)
cre[cre$recordid == "CN-002-A0-0026", "crea"] <- 119
cre[cre$recordid == "PH-002-A0-0036", "crea"] <- NA_real_

# SOFA
colnames(sofa_all) <- make.names(colnames(sofa_all), unique = TRUE)
CR_sofa_all <- sofa_all %>% filter(recordid %in% CR_final$recordid)
CR_sofa_new <- sofa_new %>% filter(recordid %in% CR_final$recordid)
# median_sofa <- median(CR_sofa_new$sofa_score, na.rm = TRUE)

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
CR_final[CR_final$recordid == "PK-003-A0-0093", "score_CVS"] <- 0
CR_final[CR_final$recordid == "PK-003-A0-0100", "score_CVS"] <- 0
CR_final[CR_final$recordid == "PK-003-A0-0105", "score_CVS"] <- 0
CR_final[CR_final$recordid == "MY-003-A0-0085", "score_CVS"] <- 0
CR_final[CR_final$recordid == "PK-003-A0-0319", "score_CVS"] <- 0
CR_final[CR_final$recordid == "MY-003-A0-0060", "score_CVS"] <- 0
CR_final[CR_final$recordid == "SL-001-A0-0050", "score_CVS"] <- 0


length(unique(CR_final$recordid)) #982

CR_final <- CR_final %>%
  #left_join(CR_sofa_sum, by = "recordid") %>%
  #left_join(CR_sofa_median, by = "recordid") %>%
  left_join(qsofa, by = "recordid") %>%
  left_join(cmb, by = "recordid")%>%
  left_join(cre, by = "recordid") %>%
  left_join(CR_icu, by = "recordid") %>%
  #for country_income2, change to 2 groups upper middle income and high income, low income and low middle income
  mutate(country_income2 = ifelse(country_income == "Low income" | country_income == "Lower middle income", "Low income and Lower middle income", "Upper middle income and High income")) %>%
  relocate(country_income2, .after = country_income) 
# %>%  select(-spec_date) %>% distinct()


# pathogen <- CR_final %>% select(recordid, org_combined, mono_poly) %>% distinct() 
# pathogen2 <- checkast %>% select(recordid, org_names) %>% distinct()
# pathogen <- pathogen %>% left_join(pathogen2, by = c("recordid", "spec_date")) %>% distinct()
# # For recordid VN-004-A0-0083, assign org_names= org_combined
# pathogen[pathogen$recordid == "VN-004-A0-0083", "org_names"] <- "K. pneumoniae"
# pathogen[pathogen$recordid == "TH-001-A0-0130", "org_names"] <- "Pseudomonas"
# pathogen[pathogen$recordid == "CN-001-A0-0032", "org_names"] <- "K. pneumoniae"
# pathogen[pathogen$recordid == "CN-002-A0-0005", "org_names"] <- "K. pneumoniae"
# pathogen[pathogen$recordid == "IN-004-A0-0007", "org_names"] <- "K. pneumoniae"
# pathogen[pathogen$recordid == "IN-004-A0-0009", "org_names"] <- "K. pneumoniae"

# Compute the median of non-NA values in sofa_cre_2 and sofa_cre_1 * 88.4
# median_crea <- cre %>%
#   mutate(temp_crea = ifelse(is.na(sofa_cre_2), sofa_cre_1 * 88.4, sofa_cre_2)) %>%
#   pull(temp_crea) %>%
#   median(na.rm = TRUE)
# 
# # Mutate new column 'crea' based on conditions
# cre <- cre %>%
#   mutate(crea = case_when(
#     !is.na(sofa_cre_2) ~ sofa_cre_2,                   # If sofa_cre_2 is not NA, use it
#     !is.na(sofa_cre_1) ~ sofa_cre_1 * 88.4,            # Else, use sofa_cre_1 * 88.4
#     TRUE ~ median_crea                                 # Else, use median imputation
#   ))
# 
# # For recordid PH-002-A0-0036, assugn crea as median_crea
# cre[cre$recordid == "PH-002-A0-0036", "crea"] <- median_crea
# cre[cre$recordid == "CN-002-A0-0026", "crea"] <- 119
# 
# cre <- cre %>%  select(recordid, crea)


length(unique(CR_final$recordid)) #859 (26 apr)  #740 if select resistant only in CR

CR_final_1row <-CR_final %>% select(-organism) %>% distinct() 
length(unique(CR_final_1row$recordid)) #753
nrow(CR_final_1row) #753

# filter rows with recordid more than 1 row
check <- CR_final_1row %>%
  group_by(recordid) %>%
  filter(n() > 1) %>%
  ungroup()

### Multiple imputation ###
# Step 1: Define the variables with missing data to impute and predictor variables
vars_to_impute <- c("score_respiration", "score_liver", "crea")
predictors <- c(
  "icu_at_onset5", "vent_at_onset5", "los_onset5", "iculos_onset5", "mvdur_onset5", "age_new", 
  "score_respiration", "score_coagulation", "score_liver", "score_CVS", "score_CNV","score_renal", 
  "qsofa", "hpd_admreason","comorbidities_Chalson",
  "cmb_comorbidities___none", "cmb_comorbidities___aids", "cmb_comorbidities___onc",
  "cmb_comorbidities___cpd", "cmb_comorbidities___cog", "cmb_comorbidities___rheu", "cmb_comorbidities___dem",
  "cmb_comorbidities___diab", "cmb_comorbidities___diad", "cmb_comorbidities___hop", "cmb_comorbidities___hivwa",
  "cmb_comorbidities___hivna", "cmb_comorbidities___mlr", "cmb_comorbidities___mal", "cmb_comorbidities___mst",
  "cmb_comorbidities___mld", "cmb_comorbidities___liv", "cmb_comorbidities___pep", "cmb_comorbidities___renal",
  "cmb_comorbidities___tub", "cmb_comorbidities___oth", 
  "crea", "hai_icu48days", "hai_have_med_device___vent"
)

# Step 2: Initialize imputation setup
ini <- mice(CR_final_1row, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix

# Step 3: Set up methods and predictor matrix
meth[] <- "" # Do not impute unless specified
meth[vars_to_impute] <- "pmm"
pred[,] <- 0 #This zeroes out the entire predictor matrix (i.e., no variables predict anything by default)
pred[vars_to_impute, predictors] <- 1  # activates only chosen predictors to predict your chosen

# Step 4: Perform imputation (only impute 4 variables, use only the listed predictors, run 20 imputations ie create 20 complete datasets with diff imputed values)
imp <- mice(CR_final_1row, method = meth, predictorMatrix = pred, m = 20, seed = 123)

# Step 5: Extract pooled values using Rubin’s rule and average imputations
# Convert to long format and compute mean across imputations
imputed_values <- complete(imp, "long") %>%
  group_by(.id) %>%
  summarise(across(all_of(vars_to_impute), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(across(c(score_respiration, score_liver, crea), round)) %>%
  rename_with(~ paste0(., "_imp"), all_of(vars_to_impute))

# Step 6: Bind imputed values back to original data
CR_final_1row <- CR_final_1row %>%
  bind_cols(imputed_values %>% select(-.id))

CR_sofa_imputed <- CR_final_1row %>%
  select(recordid, score_respiration_imp, score_respiration, score_CVS, score_liver_imp,score_liver, 
         score_coagulation,score_CNV, score_renal, sofa_score, crea, crea_imp) %>%
  # new column is the sum of score_respiration_imp, score_CVS_imp, score_liver_imp, score_coagulation, score_CNV, score_renal
  mutate(sofa_imputed = rowSums(select(., score_respiration_imp, score_CVS, score_liver_imp, score_coagulation, score_CNV, score_renal), na.rm = TRUE))
skim(CR_sofa_imputed)


# Convert columns to appropriate classes
CR_final <- CR_final_1row %>%
  mutate(
    malignancy = if_else(cmb_comorbidities___mst == 1 | cmb_comorbidities___onc == 1, 1, 0),
    diabetes   = if_else(cmb_comorbidities___diad == 1 | cmb_comorbidities___diab == 1, 1, 0),
    liver      = if_else(cmb_comorbidities___liv == 1 | cmb_comorbidities___mld == 1, 1, 0),
    renal      = if_else(cmb_comorbidities___renal == 1, 1, 0)) %>%
  mutate(sofa_imp = rowSums(select(., score_respiration_imp, score_CVS, score_liver_imp, score_coagulation, score_CNV, score_renal), na.rm = TRUE)) %>%
  mutate(delay = ifelse(delay < 0, 0, delay))%>%
  mutate(delay_group = ifelse(delay == 0, 0, 1)) %>%
  mutate(
    # Continuous variables as numeric
    across(c(age_new, sofa_score, comorbidities_Chalson,  
             los_onset5, iculos_onset5, mvdur_onset5, mortday_onset5, delay, 
             #sofa_score_sum, sofa_score_median, 
             sofa_imp,
             score_respiration, score_respiration_imp, score_coagulation, score_liver, score_liver_imp, score_CVS, score_CNV, score_renal,
             fup_day_onset5, qsofa,crea, crea_imp), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset5, vent_at_onset5, 
             mort_21d_onset5, 
             #mono_poly, 
             UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, malignancy, diabetes, liver, renal,
             delay_group, country_income2, hai_icu48days, hai_have_med_device___vent, monopoly,
             aci, pae, ent, cr_aci, cr_ent, cr_pae), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined_new, arm, arm2,arm3,
             #organism,
             ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

skim(CR_final)


#### For CRAB ####
# Function for analysis
estimate_mortality <- function(data, treatment_var, weight_var) {
  library(ggplot2)
  library(dplyr)
  library(nnet)
  library(patchwork)
  
  # Create treatment column dynamically
  data <- data %>%
    mutate(treatment = !!sym(treatment_var))
  
  # ---- Crude model ----
  crude_model <- glm(mort_21d_onset5 ~ treatment, data = data, family = "binomial")
  crude_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_crude <- predict(crude_model, newdata = crude_newdata, type = "link", se.fit = TRUE)
  crude_newdata$risk <- plogis(pred_crude$fit)
  crude_newdata$lower <- plogis(pred_crude$fit - 1.96 * pred_crude$se.fit)
  crude_newdata$upper <- plogis(pred_crude$fit + 1.96 * pred_crude$se.fit)
  crude_newdata$analysis <- "Crude"
  
  # ---- Multinomial logistic regression for inverse probability weights ----
  multinom_model <- multinom(as.formula(paste(treatment_var, "~ age_new  + sex + country_income2 + 
                                                comorbidities_Chalson  +
                                             sofa_imp + infection_types + 
                                              los_onset5 + hai_icu48days + hai_have_med_device___vent + monopoly +
                                              delay + ent + pae")), 
                             data = data)
  
  prob_matrix <- predict(multinom_model, newdata = data, type = "probs")
  prob <- mapply(function(i, trt) prob_matrix[i, trt], seq_len(nrow(data)), data$treatment)
  
  # Unstabilized weights
  data$weight_unstabilized <- 1 / prob
  w2_99 <- quantile(data$weight_unstabilized, 0.99)
  data$weight_unstabilized <- pmin(data$weight_unstabilized, w2_99)
  
  # Stabilized weights
  marginal_probs <- prop.table(table(data$treatment))
  data$weight_stabilized <- mapply(function(trt, p) marginal_probs[trt] / p, data$treatment, prob)
  sw_99 <- quantile(data$weight_stabilized, 0.99)
  data$weight_stabilized <- pmin(data$weight_stabilized, sw_99)
  
  # Choose which weight to use
  data$final_weight <- if (weight_var == "w2") data$weight_unstabilized else data$weight_stabilized
  
  # ---- MSM model using final weights ----
  msm_model <- glm(mort_21d_onset5 ~ treatment, data = data, family = "binomial", weights = final_weight)
  adjusted_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_adj <- predict(msm_model, newdata = adjusted_newdata, type = "link", se.fit = TRUE)
  adjusted_newdata$risk <- plogis(pred_adj$fit)
  adjusted_newdata$lower <- plogis(pred_adj$fit - 1.96 * pred_adj$se.fit)
  adjusted_newdata$upper <- plogis(pred_adj$fit + 1.96 * pred_adj$se.fit)
  adjusted_newdata$analysis <- "Adjusted"
  
  # ---- Combine results and convert to percent ----
  result_table <- bind_rows(crude_newdata, adjusted_newdata) %>%
    mutate(across(c(risk, lower, upper), ~ .x * 100))
  
  # Order treatment and analysis for plotting
  result_table <- result_table %>%
    mutate(treatment = factor(treatment, levels = rev(unique(data$treatment))),
           analysis = factor(analysis, levels = c("Crude", "Adjusted")))
  
  # ---- Summary of weights ----
  weight_summary <- data %>%
    summarise(
      mean   = mean(final_weight, na.rm = TRUE),
      median = median(final_weight, na.rm = TRUE),
      min    = min(final_weight, na.rm = TRUE),
      max    = max(final_weight, na.rm = TRUE),
      sd     = sd(final_weight, na.rm = TRUE)
    )
  print(weight_summary)
  
  # ---- Plot ----
  plot <- ggplot(result_table, aes(x = risk, y = treatment)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    facet_wrap(~ analysis, ncol = 1, scales = "free_y") +
    scale_x_continuous(name = "Estimated 21-Day Mortality Risk (%)", limits = c(0, 100)) +
    labs(y = "Treatment", title = paste("21-Day Mortality Risk by", treatment_var, "with", weight_var)) +
    theme_minimal(base_size = 14)
  
  print(plot)
  return(result_table)
}

#estimate_mortality(CR_aida, "arm2", "w2")  # Uses truncated unstabilized weights
#estimate_mortality(CR_aida, "arm2", "sw")  # Uses truncated stabilized weights

# # Using CR_final (514 CRAB, include polymicrobial)
# table(CR_final$arm) # 4 arms. Polymyxin mono 119, Polymyxin combi 108, Polymyxin sulbactam 95, sulbactam based 192
# table(CR_final$arm2) # Drop. 5 arms. Polymyxin mono 119, Polymyxin combi 108, Polymyxin sulbactam 95, Sulbactam combi 23, Sulbactam mono 169
# table(CR_final$arm3)

#estimate_mortality(CR_final, "arm", "w2")  # residual confounding +++ as polymyxin monotherapy superior to combination
#estimate_mortality(CR_final, "arm", "sw")  

#estimate_mortality(CR_final, "number", "w2")
#estimate_mortality(CR_final, "number", "sw")

crab$arm <- droplevels(crab$arm)
table(crab$arm)  
estimate_mortality(crab, "arm", "w2")


###########
#### For CRE and CRPAE ####
# Function for analysis
estimate_mortality <- function(data, treatment_var, weight_var) {
  library(ggplot2)
  library(dplyr)
  library(nnet)
  library(patchwork)
  
  # Create treatment column dynamically
  data <- data %>%
    mutate(treatment = !!sym(treatment_var))
  
  # ---- Crude model ----
  crude_model <- glm(mort_21d_onset5 ~ treatment, data = data, family = "binomial")
  crude_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_crude <- predict(crude_model, newdata = crude_newdata, type = "link", se.fit = TRUE)
  crude_newdata$risk <- plogis(pred_crude$fit)
  crude_newdata$lower <- plogis(pred_crude$fit - 1.96 * pred_crude$se.fit)
  crude_newdata$upper <- plogis(pred_crude$fit + 1.96 * pred_crude$se.fit)
  crude_newdata$analysis <- "Crude"
  
  # ---- Multinomial logistic regression for inverse probability weights ----
  multinom_model <- multinom(as.formula(paste(treatment_var, "~ age_new  + sex + country_income2 + 
                                                comorbidities_Chalson  +
                                              sofa_imp + infection_types + 
                                              los_onset5 + hai_icu48days + hai_have_med_device___vent + monopoly +
                                              delay + aci")), 
                             data = data)
  
  prob_matrix <- predict(multinom_model, newdata = data, type = "probs")
  prob <- mapply(function(i, trt) prob_matrix[i, trt], seq_len(nrow(data)), data$treatment)
  
  # Unstabilized weights
  data$weight_unstabilized <- 1 / prob
  w2_99 <- quantile(data$weight_unstabilized, 0.99)
  data$weight_unstabilized <- pmin(data$weight_unstabilized, w2_99)
  
  # Stabilized weights
  marginal_probs <- prop.table(table(data$treatment))
  data$weight_stabilized <- mapply(function(trt, p) marginal_probs[trt] / p, data$treatment, prob)
  sw_99 <- quantile(data$weight_stabilized, 0.99)
  data$weight_stabilized <- pmin(data$weight_stabilized, sw_99)
  
  # Choose which weight to use
  data$final_weight <- if (weight_var == "w2") data$weight_unstabilized else data$weight_stabilized
  
  # ---- MSM model using final weights ----
  msm_model <- glm(mort_21d_onset5 ~ treatment, data = data, family = "binomial", weights = final_weight)
  adjusted_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_adj <- predict(msm_model, newdata = adjusted_newdata, type = "link", se.fit = TRUE)
  adjusted_newdata$risk <- plogis(pred_adj$fit)
  adjusted_newdata$lower <- plogis(pred_adj$fit - 1.96 * pred_adj$se.fit)
  adjusted_newdata$upper <- plogis(pred_adj$fit + 1.96 * pred_adj$se.fit)
  adjusted_newdata$analysis <- "Adjusted"
  
  # ---- Combine results and convert to percent ----
  result_table <- bind_rows(crude_newdata, adjusted_newdata) %>%
    mutate(across(c(risk, lower, upper), ~ .x * 100))
  
  # Order treatment and analysis for plotting
  result_table <- result_table %>%
    mutate(treatment = factor(treatment, levels = rev(unique(data$treatment))),
           analysis = factor(analysis, levels = c("Crude", "Adjusted")))
  
  # ---- Summary of weights ----
  weight_summary <- data %>%
    summarise(
      mean   = mean(final_weight, na.rm = TRUE),
      median = median(final_weight, na.rm = TRUE),
      min    = min(final_weight, na.rm = TRUE),
      max    = max(final_weight, na.rm = TRUE),
      sd     = sd(final_weight, na.rm = TRUE)
    )
  print(weight_summary)
  
  # ---- Plot ----
  plot <- ggplot(result_table, aes(x = risk, y = treatment)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    facet_wrap(~ analysis, ncol = 1, scales = "free_y") +
    scale_x_continuous(name = "Estimated 21-Day Mortality Risk (%)", limits = c(0, 100)) +
    labs(y = "Treatment", title = paste("21-Day Mortality Risk by", treatment_var, "with", weight_var)) +
    theme_minimal(base_size = 14)
  
  print(plot)
  return(result_table)
}

cre_crpae$arm3 <- droplevels(cre_crpae$arm3)
table(cre_crpae$arm3)
estimate_mortality(cre_crpae, "arm3", "w2")
