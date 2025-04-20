CRAB <- CR_final %>% filter(cr_aci==1)
#%>% filter(monopoly == "monomicrobial")
length(unique(CRAB$recordid)) #487 (370 mono)   #445 (no tige)
nrow(CRAB) #487   #445

CRE_CRPAE <- CR_final %>% filter(cr_pae==1 | cr_ent==1)
# %>% filter(monopoly == "monomicrobial")
length(unique(CRE_CRPAE$recordid)) #316 (221 mono)  #264 (no tige)
nrow(CRE_CRPAE) #316

CR_final_mono <- CR_final %>% filter(monopoly == "monomicrobial") 
length(unique(CR_final_mono$recordid)) #591   #532 (no tige)
nrow(CR_final_mono) #591

#skim(CRAB)
#skim(CRE_CRPAE)

check <- CRE_CRPAE %>% select(recordid, onset4, org_combined_new, anti_onset4, arm, arm2, arm3, aci, ent, pae, cr_aci, cr_pae, cr_ent)%>%
  filter(!is.na(arm3))%>%
  filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(check$recordid)) #261 (with tige), #211 (no tige)

crab <- CRAB %>%
  filter(!is.na(arm2))
#%>%
#  filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(crab$recordid)) #443 (no tige)  vs 485 (with tige)

crab %>%
  count(arm) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRAB") 

crab_mono <- CRAB %>%
  filter(!is.na(arm2))%>%
  filter(monopoly =="monomicrobial")
length(unique(crab_mono$recordid)) #345

cre_crpae <- CRE_CRPAE %>%
  filter(!is.na(arm3)) 
# %>%
# filter(!grepl("Tigecycline|Fosfomycin", anti_onset4))
length(unique(cre_crpae$recordid)) #261 (with tige), #211 (no tige)

cre_crpae %>% count(arm3) %>%
  arrange(desc(n)) %>%
  kable(col.names = c("Antibiotic", "Count"), caption = "Distribution of study arms for CRE and CRPAE")


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
  crude_model <- glm(mort_21d_onset4 ~ treatment, data = data, family = "binomial")
  crude_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_crude <- predict(crude_model, newdata = crude_newdata, type = "link", se.fit = TRUE)
  crude_newdata$risk <- plogis(pred_crude$fit)
  crude_newdata$lower <- plogis(pred_crude$fit - 1.96 * pred_crude$se.fit)
  crude_newdata$upper <- plogis(pred_crude$fit + 1.96 * pred_crude$se.fit)
  crude_newdata$analysis <- "Crude"
  
  # ---- Multinomial logistic regression for inverse probability weights ----
  multinom_model <- multinom(as.formula(paste(treatment_var, "~ age_new + sex + country_income2 + comorbidities_Chalson + 
                                              diabetes + malignancy + renal + liver +
                                              sofa_imp + infection_types + 
                                              los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
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
  msm_model <- glm(mort_21d_onset4 ~ treatment, data = data, family = "binomial", weights = final_weight)
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

crab_mono$arm <- droplevels(crab_mono$arm)
estimate_mortality(crab_mono, "arm", "w2")

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
  crude_model <- glm(mort_21d_onset4 ~ treatment, data = data, family = "binomial")
  crude_newdata <- data.frame(treatment = levels(factor(data$treatment)))
  pred_crude <- predict(crude_model, newdata = crude_newdata, type = "link", se.fit = TRUE)
  crude_newdata$risk <- plogis(pred_crude$fit)
  crude_newdata$lower <- plogis(pred_crude$fit - 1.96 * pred_crude$se.fit)
  crude_newdata$upper <- plogis(pred_crude$fit + 1.96 * pred_crude$se.fit)
  crude_newdata$analysis <- "Crude"
  
  # ---- Multinomial logistic regression for inverse probability weights ----
  multinom_model <- multinom(as.formula(paste(treatment_var, "~ age_new + sex + country_income2 + comorbidities_Chalson + 
                                              diabetes + malignancy + renal + liver +
                                              sofa_imp + infection_types + 
                                              los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
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
  msm_model <- glm(mort_21d_onset4 ~ treatment, data = data, family = "binomial", weights = final_weight)
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

#cre_crpae$arm3 <- droplevels(cre_crpae$arm3)
#table(cre_crpae$arm3)
estimate_mortality(cre_crpae, "arm3", "w2")







##############################################################################################################
############################################### STOP HERE ####################################################
##############################################################################################################
# Using CR_final but only 4 cleaner arms
CR_final_4arm <- CR_final %>% 
  filter(anti_onset4 %in% c(
    "Colistin", "Colistin + Meropenem", "Colistin + Doripenem", "Colistin + Imipenem",
    "Colistin + Neb Col", "Polymyxin B", "Neb Col + Polymyxin B", "Meropenem + Polymyxin B",
    "Ampicillin/sulbactam", "Cefoperazone/sulbactam", "Ampicillin/sulbactam + Cefoperazone/sulbactam",
    "Cefoperazone/sulbactam + Neb Col", "Ampicillin/sulbactam + Polymyxin B",
    "Ampicillin/sulbactam + Colistin", "Cefoperazone/sulbactam + Polymyxin B",
    "Cefoperazone/sulbactam + Colistin", "Colistin + Ampicillin/sulbactam",
    "Colistin + Cefoperazone/sulbactam", "Polymyxin B + Ampicillin/sulbactam",
    "Polymyxin B + Cefoperazone/sulbactam"
  ))
CR_final_4arm$arm <- droplevels(CR_final_4arm$arm)
CR_final_4arm$arm2 <- droplevels(CR_final_4arm$arm2)  # poly combi 77, poly mono 119, poly sul 64, sul 169
table(CR_final_4arm$arm)  # 77 Polymyxin carbapenem, 119 polymyxin mono, 64 polymyxin sulbactam, 169 sulbactam monotherapy
length(unique(CR_final_4arm$recordid))  # 429
estimate_mortality(CR_final_4arm, "arm2", "w2")
estimate_mortality(CR_final_4arm, "arm2", "sw")

# Using CR_final_mono
CR_final_mono <- CR_final %>% filter(monopoly == "monomicrobial")  #396 instead of 514
table(CR_final_mono$arm2) # Poly combi 71, poly mono 92, poly sul 77, sul 156 (14 combi and 142 mono)
length(unique(CR_final_mono$recordid))  #396
estimate_mortality(CR_final_mono, "arm", "w2")  
estimate_mortality(CR_final_mono, "arm", "sw") 

estimate_mortality(CR_final_mono, "arm2", "w2")
estimate_mortality(CR_final_mono, "arm2", "sw")

# Using CR_final strictly Sulbactam mono VS Polymyxin mono VS Polymyxin + Carbapenem only
CR_final_strict <- CR_final %>% 
  filter(anti_onset4 %in% c("Colistin", "Colistin + Meropenem", "Colistin + Doripenem", "Colistin + Imipenem",
                          "Colistin + Neb Col","Polymyxin B", "Neb Col + Polymyxin B",  "Meropenem + Polymyxin B", "Neb Col + Polymyxin B",
                          "Ampicillin/sulbactam","Cefoperazone/sulbactam", "Ampicillin/sulbactam + Cefoperazone/sulbactam",
                          "Cefoperazone/sulbactam + Neb Col"))%>%
  filter(arm2=="Polymyxin combination" | arm2 == "Polymyxin monotherapy" | arm2 == "Sulbactam monotherapy")
CR_final_strict$arm <- droplevels(CR_final_strict$arm)
CR_final_strict$arm2 <- droplevels(CR_final_strict$arm2)
table(CR_final_strict$arm)  # poly mero 77, poly mono 119, sulbactam mono 169
table(CR_final_strict$arm2)
length(unique(CR_final_strict$recordid)) # 365
estimate_mortality(CR_final_strict, "arm2", "w2")
estimate_mortality(CR_final_strict, "arm2", "sw")

# Using CR_final_mono strictly Sulbactam mono (141) VS Polymyxin mono (92) VS Polymyxin + Carbapenem (54) only
CR_final_mono_strict <- CR_final_mono %>% 
  filter(anti_onset4 %in% c("Colistin", "Colistin + Meropenem", "Colistin + Doripenem", "Colistin + Imipenem",
                            "Colistin + Neb Col","Polymyxin B", "Neb Col + Polymyxin B",  "Meropenem + Polymyxin B", "Neb Col + Polymyxin B",
                            "Ampicillin/sulbactam","Cefoperazone/sulbactam", "Ampicillin/sulbactam + Cefoperazone/sulbactam",
                            "Cefoperazone/sulbactam + Neb Col"))%>%
  filter(arm2=="Polymyxin combination" | arm2 == "Polymyxin monotherapy" | arm2 == "Sulbactam monotherapy")
CR_final_mono_strict$arm <- droplevels(CR_final_mono_strict$arm)
CR_final_mono_strict$arm2 <- droplevels(CR_final_mono_strict$arm2)
table(CR_final_mono_strict$arm)
table(CR_final_mono_strict$arm2)
estimate_mortality(CR_final_mono_strict, "arm2", "w2")
estimate_mortality(CR_final_mono_strict, "arm2", "sw")

CRAB_strict <- CRAB %>% 
  filter(anti_onset4 %in% c("Colistin", "Colistin + Meropenem", "Colistin + Doripenem", "Colistin + Imipenem",
                                                  "Colistin + Neb Col","Polymyxin B", "Neb Col + Polymyxin B",  "Meropenem + Polymyxin B", "Neb Col + Polymyxin B",
                                                  "Ampicillin/sulbactam","Cefoperazone/sulbactam", "Ampicillin/sulbactam + Cefoperazone/sulbactam",
                                                  "Cefoperazone/sulbactam + Neb Col"))%>%
  filter(arm=="Polymyxin combination" | arm == "Polymyxin monotherapy" | arm == "Sulbactam based")
CRAB_strict$arm <- droplevels(CRAB_strict$arm)
table(CRAB_strict$arm)
estimate_mortality(CRAB_strict, "arm", "w2")

# Using CR_final_4arm but BSI only (sample size too small for subgroup analysis)
CR_final_4arm_bsi <- CR_final_4arm %>%
  filter (infection_types == "BSI")
table(CR_final_4arm_bsi$arm)
table(CR_final_4arm_bsi$arm2) 
estimate_mortality(CR_final_4arm_bsi, "arm2", "w2")
estimate_mortality(CR_final_4arm_bsi, "arm2", "sw")

# Using CR_final_4arm but monomicrobial only (sample size too small for subgroup analysis)
CR_final_4arm_mono <- CR_final_4arm %>%
  filter (monopoly == "monomicrobial")
table(CR_final_4arm_bsi$arm)
table(CR_final_4arm_bsi$arm2)  # 25 poly combi, 43 poly mono, 27 poly sul, 49 sul mono
estimate_mortality(CR_final_4arm_bsi, "arm2", "w2")
estimate_mortality(CR_final_4arm_bsi, "arm2", "sw")

table(CR_final$anti_onset4)

# Emulate Treat-GNB
CR_treatgnb <- CR_final %>%
  mutate(
    arm_treatgnb = case_when(
      # Rule 1: Tigecycline + (Colistin or Polymyxin B)
      grepl("Tigecycline", anti_onset4, ignore.case = TRUE) & 
        (grepl("Colistin", anti_onset4, ignore.case = TRUE) | 
           grepl("Polymyxin B", anti_onset4, ignore.case = TRUE)) ~ "polymyxin tigecycline",
      
      # Rule 2: Sulbactam + (Colistin or Polymyxin B)
      (grepl("Ampicillin/sulbactam", anti_onset4, ignore.case = TRUE) | 
         grepl("Cefoperazone/sulbactam", anti_onset4, ignore.case = TRUE)) & 
        (grepl("Colistin", anti_onset4, ignore.case = TRUE) | 
           grepl("Polymyxin B", anti_onset4, ignore.case = TRUE)) ~ "polymyxin sulbactam",
      
      # Rule 3: Monotherapy from `arm`
      arm == "Polymyxin monotherapy" ~ "polymyxin monotherapy",
      
      # Rule 4: Carbapenem + (Colistin or Polymyxin B), and no other drugs
      (grepl("Meropenem|Imipenem|Doripenem", anti_onset4, ignore.case = TRUE)) & 
        (grepl("Colistin|Polymyxin B", anti_onset4, ignore.case = TRUE)) &
        !grepl("Ampicillin/sulbactam|Cefoperazone/sulbactam|Tigecycline|Neb Col", anti_onset4, ignore.case = TRUE) ~ "polymyxin carbapenem",
      
      # Default: NA (or you can assign "other" if preferred)
      TRUE ~ NA_character_
    )
  ) %>%
  relocate(
    arm_treatgnb,
    .after = arm
  )

# Drop arm_treatgnb if it is NA
CR_treatgnb <- CR_treatgnb %>%
  filter(!is.na(arm_treatgnb))
table(CR_treatgnb$arm_treatgnb)  # 77 poly penem, 119 poly mono, 92 poly sul, 20 poly tige
estimate_mortality(CR_treatgnb, "arm_treatgnb", "w2")
estimate_mortality(CR_treatgnb, "arm_treatgnb", "sw")

