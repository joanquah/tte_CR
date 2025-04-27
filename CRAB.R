CR_fin <- crab %>% select (recordid, onset4,arm, arm2, anti_onset4, number, infection_types,mort_21d_onset4,
                               monopoly,org_combined_new, aci)
#levo <- CR_fin %>% filter(grepl("Levofloxacin", anti_onset4)) 
#mino <- CR_fin %>% filter(grepl("Minocycline", anti_onset4)) 
fosfo <- CR_fin %>% filter(grepl("Fosfomycin", anti_onset4))  #16
tige <- CR_fin %>% filter(grepl("Tigecycline", anti_onset4)) #23
nebcol <- CR_fin %>% filter(grepl("Neb Col", anti_onset4)) # 8
sulbacta <- CR_fin %>% filter(grepl("sulbactam", anti_onset4)) # 287

CR_fin <- CR_fin %>%
  mutate(pae = ifelse(grepl("Pseudomonas", org_combined_new), 1, 0))%>%
  mutate(enterobacterales = ifelse(grepl("K. pneumoniae|E. coli|Klebsiella|Enterobacter", org_combined_new), 1, 0))
############################ 
## TARGET TRIAL EMULATION ##
############################

# Run preparing data for analysis
#########################################################################################
# Crude analysis
# Fit logistic model and estimate 21 day mortality risks
model1 <- glm(mort_21d_onset4~ arm, data = crab, family = "binomial")
summary(model1)

# Generate new dataset to store estimated risks
new_data <- data.frame(arm = levels(crab$arm))

# list unique CR_targeted_abx_name
levels(crab$arm)

# Obtained estimated risks based on model
new_data$risk <- predict(model1, newdata = new_data, type = "response")
new_data

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

############################################################################################
# Fit multinomial logistic regression model
model2 <- multinom (arm ~ age_new  + sex + country_income2 + 
                      comorbidities_Chalson  +
                      diabetes + malignancy + renal + liver +
                      qsofa + infection_types + 
                      los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
                      delay + ent + pae,
                    data=crab)
summary(model2) # 


#Use the predicted probabilities from this model to estimate nonstabilized IP weights
# Generate predicted probabilities for each treatment arm
crab$prob_matrix <- predict(model2, newdata = crab, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.
head(crab$prob_matrix) 

# For each patient, select the probability corresponding to the treatment they actually received.
# Extract the probability corresponding to the observed treatment for each row
crab$prob <- mapply(function(i, treatment) {
  crab$prob_matrix[i, treatment]
}, i = seq_len(nrow(crab)), treatment = crab$arm)  #This assigns prob[i] as the probability of the observed treatment for individual i.

# Compute weight (nonstabilized IPW)
crab$w2 <- 1 / crab$prob

# Before truncation
# hist(crab$w2)

# Truncating weights at 99th percentile
crab$w2[crab$w2 > quantile(crab$w2, 0.99)] <- quantile(crab$w2, 0.99)

# Check the distribution of the nonstabilized weights
summary(crab$w2)  # Check distribution of weights
sd(crab$w2)

# After truncation
# hist(crab$w2[crab$w2 <= quantile(crab$w2, 0.99)])
###
# use MSM with nonstabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter error due to non-integer number of successes as a result of weights
# msm.w <- geeglm(mort_21d_ast ~ CR_targeted_abx_name, data=crab, weights=w2, id=recordid, family=binomial())
msm.w <- glm(mort_21d_onset4 ~ arm, data=crab, family=binomial(), weights=w2)
summary(msm.w)
exp(coef(msm.w))  # Convert log-odds to odds ratios
exp(confint(msm.w))  # 95% confidence interval for odds ratios

# Outputting risks
new_data <- data.frame(arm = levels(crab$arm))
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

#####################################################################################################

# checkast <-ast_all_index %>% select(recordid,spec_date, org_names_all, pathogen_group_combined,
#                                     `Amikacin`, `Gentamicin`, `Fosfomycin`, `Minocycline`, `Tigecycline`, `Levofloxacin`,
#                                     `Ampicillin/sulbactam`, `Cefoperazone/sulbactam`) %>%
#   filter(recordid %in% crab$recordid)
# 
# amikacin <- CR_names_4 %>% filter(anti_new== "Amikacin") %>% left_join(checkast %>% select(recordid, spec_date, `Amikacin`), by = c("recordid", "spec_date")) %>% distinct()
# gentamicin <- CR_names_4 %>% filter(anti_new== "Gentamicin") %>% left_join(checkast %>% select(recordid, spec_date, `Gentamicin`), by = c("recordid", "spec_date")) %>% distinct()
# tigecycline <- CR_names_4 %>% filter(anti_new== "Tigecycline") %>% left_join(checkast %>% select(recordid, spec_date, `Tigecycline`), by = c("recordid", "spec_date")) %>% distinct()
# levofloxacin <- CR_names_4 %>% filter(anti_new== "Levofloxacin") %>% left_join(checkast %>% select(recordid, spec_date, `Levofloxacin`), by = c("recordid", "spec_date")) %>% distinct()
# minocyline <- CR_names_4 %>% filter(anti_new== "Minocycline") %>% left_join(checkast %>% select(recordid, spec_date, `Minocycline`), by = c("recordid", "spec_date")) %>% distinct()
# fosfomycin <- CR_names_4 %>% filter(anti_new== "Fosfomycin") %>% left_join(checkast %>% select(recordid, spec_date, `Fosfomycin`), by = c("recordid", "spec_date")) %>% distinct()
# 
# ####################################################################################################
# # Show data
# colnames(crab)
# crab_column <- crab %>% 
#   select (recordid, arm,arm_group, number, age_new, sex, country_income2,
#           comorbidities_Chalson, sofa_score, sofa_score_sum, sofa_score_median, qsofa,
#           infection_types, delay, delay_group, crea,adm_ward_types_new,
#           icu_at_onset4, vent_at_onset4,iculos_onset4,mvdur_onset4, mono_poly,
#           country, mort_21d_onset4,prob, w2)
# skim(crab_column)          
# 
#                                 