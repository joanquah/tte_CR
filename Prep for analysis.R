# Preparing for analysis
### Preparing data for analaysis
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
CR_sofa_new[CR_sofa_new$recordid == "PK-003-A0-0093", "score_CVS"] <- 0
CR_sofa_new[CR_sofa_new$recordid == "PK-003-A0-0100", "score_CVS"] <- 0
CR_sofa_new[CR_sofa_new$recordid == "PK-003-A0-0105", "score_CVS"] <- 0


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
length(unique(CR_final$recordid)) #335  #982

CR_final <- CR_final %>%
  left_join(CR_sofa_sum, by = "recordid") %>%
  left_join(CR_sofa_median, by = "recordid") %>%
  left_join(qsofa, by = "recordid") %>%
  left_join(cmb, by = "recordid")%>%
  left_join(cre, by = "recordid") %>%
  left_join(CR_icu, by = "recordid") %>%
  left_join(CR_combined_new, by= "recordid") %>%
  #for country_income2, change to 2 groups upper middle income and high income, low income and low middle income
  mutate(country_income2 = ifelse(country_income == "Low income" | country_income == "Lower middle income", "Low income and Lower middle income", "Upper middle income and High income")) %>%
  relocate(country_income2, .after = country_income) 
# %>%  select(-spec_date) %>% distinct()

CR_final<- CR_final %>%
  # if hai_icu48days is "Y" , assign 1, else 0
  mutate(hai_icu48days = ifelse(hai_icu48days == "Y", 1, 0)) 

### Multiple imputation ###
# Step 1: Define the variables with missing data to impute and predictor variables
vars_to_impute <- c("score_respiration", "score_liver", "crea")
predictors <- c(
  "icu_at_onset4", "vent_at_onset4", "los_onset4", "iculos_onset4", "mvdur_onset4", "age_new", 
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
ini <- mice(CR_final, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix

# Step 3: Set up methods and predictor matrix
meth[] <- "" # Do not impute unless specified
meth[vars_to_impute] <- "pmm"
pred[,] <- 0 #This zeroes out the entire predictor matrix (i.e., no variables predict anything by default)
pred[vars_to_impute, predictors] <- 1  # activates only chosen predictors to predict your chosen

# Step 4: Perform imputation (only impute 4 variables, use only the listed predictors, run 20 imputations ie create 20 complete datasets with diff imputed values)
imp <- mice(CR_final, method = meth, predictorMatrix = pred, m = 20, seed = 123)

# Step 5: Extract pooled values using Rubin’s rule and average imputations
# Convert to long format and compute mean across imputations
imputed_values <- complete(imp, "long") %>%
  group_by(.id) %>%
  summarise(across(all_of(vars_to_impute), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(across(c(score_respiration, score_liver, crea), round)) %>%
  rename_with(~ paste0(., "_imp"), all_of(vars_to_impute))

# Step 6: Bind imputed values back to original data
CR_final <- CR_final %>%
  bind_cols(imputed_values %>% select(-.id))

CR_sofa_imputed <- CR_final %>%
  select(recordid, score_respiration_imp, score_respiration, score_CVS, score_liver_imp,score_liver, 
         score_coagulation,score_CNV, score_renal, sofa_score, crea, crea_imp) %>%
  # new column is the sum of score_respiration_imp, score_CVS_imp, score_liver_imp, score_coagulation, score_CNV, score_renal
  mutate(sofa_imputed = rowSums(select(., score_respiration_imp, score_CVS, score_liver_imp, score_coagulation, score_CNV, score_renal), na.rm = TRUE))
skim(CR_sofa_imputed)



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


# Convert columns to appropriate classes
CR_final <- CR_final %>%
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
             los_onset4, iculos_onset4, mvdur_onset4, mortday_onset4, delay, 
             sofa_score_sum, sofa_score_median, sofa_imp,
             score_respiration, score_respiration_imp, score_coagulation, score_liver, score_liver_imp, score_CVS, score_CNV, score_renal,
             fbis_score,fup_day_onset4, qsofa,crea, crea_imp), as.numeric),
    
    # Binary variables as factor
    across(c(sex, adm_ward_types_new, infection_types, icu_at_onset4, vent_at_onset4, 
             mort_21d_onset4, 
             #mono_poly, 
             UTI_source,first28_death, mortality,
             cmb_comorbidities___none, cmb_comorbidities___aids, cmb_comorbidities___onc, cmb_comorbidities___cpd,
             cmb_comorbidities___cog, cmb_comorbidities___rheu, cmb_comorbidities___dem, cmb_comorbidities___diab,
             cmb_comorbidities___diad, cmb_comorbidities___hop, cmb_comorbidities___oth, cmb_comorbidities___tub, cmb_comorbidities___renal,
             cmb_comorbidities___pep, cmb_comorbidities___liv, cmb_comorbidities___mld, cmb_comorbidities___mst,cmb_comorbidities___hivwa, cmb_comorbidities___hivna,
             cmb_comorbidities___mal, cmb_comorbidities___mlr, malignancy, diabetes, liver, renal,
             delay_group, country_income2, hai_icu48days, hai_have_med_device___vent, monopoly), as.factor),
    
    # Categorical variables as factor
    across(c(country_income, country, country_ab, org_combined_new, arm, arm2,arm3,number,
             organism,ho_dischargestatus, hpd_admreason,d28_status), as.factor),
    
    # Date variables
    across(c(adm_date, date_enrolment,hpd_adm_date, hpd_hosp_date), ymd)
  ) 

CR_final <- CR_final %>%
  #relocate(spec_date_new, .after = spec_date) %>%
  relocate(monopoly, .after = mono_poly) %>%
  relocate(org_combined_new, .after = org_combined) %>%
  # select(-spec_date) %>%
  distinct() %>%
  mutate(aci = ifelse(grepl("Acinetobacter", org_combined_new), 1, 0))

length(unique(CR_final$recordid))  #982

CRAB <- CR_final %>% filter(organism == "CRAB") 
length(unique(CRAB$recordid)) #487
nrow(CRAB) #487
CRE_CRPAE <- CR_final %>% filter(organism == "CRE_CRPAE")
length(unique(CRE_CRPAE$recordid)) #315
nrow(CRE_CRPAE) #315

# filter rows with recordid more than 1 row
check <- CRE_CRPAE %>%
  group_by(recordid) %>%
  filter(n() > 1) %>%
  ungroup()

# filter rows with organism is not NA
CR_final_mono <- CR_final %>% filter(!is.na(organism)) #973
length(unique(check$recordid)) #888
nrow(check) #973

check <- CR_final %>% select(recordid, org_combined_new, organism, monopoly, infection_types) %>%
  filter(!is.na(organism))%>%
  filter(monopoly == "monomicrobial")

check2 <- CR %>% filter (recordid %in% check$recordid) 
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

# Levofloxacin
#levo_st <- f07m %>% select(recordid, org, spec_date, ast_date, susceptibility_36, res_36,std_36, mic_36,mic_ug_36, mic_mg_36)
#mino_st <- f07m %>% select(recordid, org, spec_date, ast_date, susceptibility_39, res_39,std_39, mic_39,mic_ug_39, mic_mg_39)
  #   filter(!is.na(susceptibility_7))%>%
  #   filter(recordid %in% CR_ceftazidime_avibactam$recordid)%>%
  #   distinct()


#Distribution of study arms
# CR_final %>%
#   count(arm) %>%
#   arrange(desc(n)) %>%
#   kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms")

# Export to excel
# library(writexl)
# write_xlsx(CR_final, "CR_final_CRAB.xlsx")
# 
# # Export as Rdata
# save(CR_final, file = "CR_final_CRAB.RData")
