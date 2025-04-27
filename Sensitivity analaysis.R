# Sensitivity analysis
# Monomicrobial only
CR_final <- CR_final %>% filter(mono_poly == "monomicrobial")
CR_final <- CR_final %>% filter(infection_types == "VAP")
CR_final <- CR_final %>% filter(infection_types == "BSI")
#Add UTI source to BSI

model2i <- multinom (arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                       sofa_imp  + infection_types + 
                       hai_icu48days + hai_have_med_device___vent+
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
model2c <- multinom (arm ~ age_new + sex +  country_income + comorbidities_Chalson + 
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
CR_final$prob_matrix <- predict(model2i, newdata = CR_final, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.

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
