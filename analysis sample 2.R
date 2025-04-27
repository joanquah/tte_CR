### Arm2 ####
# Crude analysis
# Fit logistic model and estimate 21 day mortality risks
model1 <- glm(mort_21d_onset4~ arm2, data = CR_final, family = "binomial")
new_data <- data.frame(arm2 = levels(CR_final$arm2))
# levels(CR_final$arm2)
new_data$risk <- predict(model1, newdata = new_data, type = "response")
pred <- predict(model1, newdata = new_data, type = "link", se.fit = TRUE)
new_data$risk <- plogis(pred$fit)
new_data$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
new_data$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
print(new_data)
ggplot(new_data, aes(y = arm2, x = risk)) +
  geom_point(size = 4) +  # Dot for point estimate
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +  # Horizontal error bar
  labs(title = "21-Day Mortality Risk by Treatment (Unadjusted analysis)",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal()

# Adjusted analysis
# Fit multinomial logistic regression model
model2 <- multinom (arm2 ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      # cmb_comorbidities___diad + cmb_comorbidities___diab + 
                      # cmb_comorbidities___mst + cmb_comorbidities___onc + 
                      sofa_imp  + infection_types + los_onset4 +
                      hai_icu48days + hai_have_med_device___vent + delay,
                      # icu_at_onset4 + vent_at_onset4
                      # icu_at_onset5 + vent_at_onset5+
                      #+ monopoly,
                    data=CR_final)
# Generate predicted probabilities for each treatment arm
CR_final$prob_matrix <- predict(model2, newdata = CR_final, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.
CR_final$prob <- mapply(function(i, treatment) {
  CR_final$prob_matrix[i, treatment]
}, i = seq_len(nrow(CR_final)), treatment = CR_final$arm2)  # assign prob[i] as the probability of the observed treatment for individual i.
CR_final$w2 <- 1 / CR_final$prob
CR_final$w2[CR_final$w2 > quantile(CR_final$w2, 0.99)] <- quantile(CR_final$w2, 0.99)
summary(CR_final$w2) 
sd(CR_final$w2)
hist(CR_final$w2, breaks = 30, main = "Distribution of Weights", xlab = "Weight")



###
# use MSM with nonstabilized weights
msm.w <- glm(mort_21d_onset4 ~ arm2, data=CR_final, family=binomial(), weights=w2)
# use MSM with stabilized weights
msm.w <- glm(mort_21d_onset4 ~ arm2, data=CR_final, family=binomial(), weights=sw)
#exp(coef(msm.w))  # Convert log-odds to odds ratios
#exp(confint(msm.w))  # 95% confidence interval for odds ratios
# Outputting risks
new_data1 <- data.frame(arm2 = levels(CR_final$arm2))
new_data1$risk <- predict(msm.w, newdata = new_data, type = "response")
new_data1

# Obtain predicted risks and standard errors
predictions <- predict(msm.w, newdata = new_data1, type = "response", se.fit = TRUE)

# Calculate confidence intervals (95% CI)
conf_int_lower <- predictions$fit - 1.96 * predictions$se.fit
conf_int_upper <- predictions$fit + 1.96 * predictions$se.fit

# Add these confidence intervals to the new_data dataframe
new_data1$lower <- conf_int_lower
new_data1$upper <- conf_int_upper

# Output the new data with risks and confidence intervals
new_data1

# Horizontal Error Bar Plot (Whisker and Dot Plot)
ggplot(new_data1, aes(y = arm2, x = risk)) +
  geom_point(size = 3, color = "blue") +  # Dots for the predicted risks
  geom_errorbarh(aes(xmin = conf_int_lower, xmax = conf_int_upper), height = 0.2) +  # Horizontal error bars
  labs(title = "21-Day Mortality Risk by Treatment (Adjusted analysis)",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 14))

library(patchwork)

# Create individual plots
p1 <- ggplot(new_data, aes(y = arm2, x = risk)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  labs(title = "21-Day Mortality Risk by Treatment (Crude analysis)",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +
  theme_minimal()

p2 <- ggplot(new_data1, aes(y = arm2, x = risk)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  labs(title = "21-Day Mortality Risk by Treatment (Adjusted analysis)",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +
  theme_minimal()

# Combine side by side or stacked
p1 / p2  # stacked vertically
# p1 | p2  # side by side

# ### Number ### Doesnt really make sense because choice of antibiotic is more important than number of antibiotics
# # Crude analysis
# # Fit logistic model and estimate 21 day mortality risks
# model1 <- glm(mort_21d_onset4~ number, data = CR_final_mono, family = "binomial")
# summary(model1)
# 
# # Generate new dataset to store estimated risks
# new_data <- data.frame(number = levels(CR_final_mono$number))
# 
# # list unique CR_targeted_abx_name
# levels(CR_final_mono$number)
# 
# # Obtained estimated risks based on model
# new_data$risk <- predict(model1, newdata = new_data, type = "response")
# 
# # Get predictions with standard errors
# pred <- predict(model1, newdata = new_data, type = "link", se.fit = TRUE)
# 
# # Convert log-odds to probability scale
# new_data$risk <- plogis(pred$fit)
# 
# # Compute confidence intervals (logit scale)
# new_data$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
# new_data$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
# 
# # View the final dataset with estimated risks and confidence intervals
# print(new_data)
# 
# # Horizontal Error Bar Plot (Whisker and Dot Plot)
# ggplot(new_data, aes(y = number, x = risk)) +
#   geom_point(size = 4) +  # Dot for point estimate
#   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +  # Horizontal error bar
#   labs(title = "21-Day Mortality Risk by Treatment",
#        x = "Estimated Risk (95% CI)", y = "Treatment") +
#   xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
#   theme_minimal()
# 
# 
# 
# # Fit multinomial logistic regression model
# model2 <- multinom (number ~ age_new + sex + country_income + 
#                       comorbidities_Chalson + 
#                       sofa_imp + infection_types +
#                       # hai_icu48days + hai_have_med_device___vent + 
#                       # crea_imp + delay +
#                       icu_at_onset4 + vent_at_onset4, 
#                       # los_onset4, 
#                       # icu_at_onset5 + vent_at_onset5+
#                       # monopoly,
#                     data=CR_final_mono)
# summary(model2) 
# 
# #Use the predicted probabilities from this model to estimate nonstabilized IP weights
# # Generate predicted probabilities for each treatment arm
# CR_final_mono$prob_matrix <- predict(model2, newdata = CR_final_mono, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.
# 
# # Extract the probability corresponding to the observed treatment 
# CR_final_mono$prob <- mapply(function(i, treatment) {
#   CR_final_mono$prob_matrix[i, treatment]
# }, i = seq_len(nrow(CR_final_mono)), treatment = CR_final_mono$number)  # assign prob[i] as the probability of the observed treatment for individual i.
# 
# # Compute weight (nonstabilized IPW)
# CR_final_mono$w2 <- 1 / CR_final_mono$prob
# # Truncating weights at 99th percentile
# CR_final_mono$w2[CR_final_mono$w2 > quantile(CR_final_mono$w2, 0.99)] <- quantile(CR_final_mono$w2, 0.99)
# # Check the distribution of the nonstabilized weights
# summary(CR_final_mono$w2) 
# sd(CR_final_mono$w2)
# 
# ###
# # use MSM with nonstabilized weights
# options(warn=-1) # Need to suppress warning or else geeglm will encounter error due to non-integer number of successes as a result of weights
# # msm.w <- geeglm(mort_21d_onset5 ~ arm, data=CR_final, weights=w2, id=recordid, family=binomial())
# # msm.w <- glm(mort_21d_onset5 ~ arm, data=CR_final, family=binomial(), weights=w2)
# msm.w <- glm(mort_21d_onset4 ~ number, data=CR_final_mono, family=binomial(), weights=w2)
# summary(msm.w)
# #exp(coef(msm.w))  # Convert log-odds to odds ratios
# #exp(confint(msm.w))  # 95% confidence interval for odds ratios
# 
# # Outputting risks
# new_data <- data.frame(number = levels(CR_final_mono$number))
# new_data$risk <- predict(msm.w, newdata = new_data, type = "response")
# new_data
# 
# # Obtain predicted risks and standard errors
# predictions <- predict(msm.w, newdata = new_data, type = "response", se.fit = TRUE)
# 
# # Calculate confidence intervals (95% CI)
# conf_int_lower <- predictions$fit - 1.96 * predictions$se.fit
# conf_int_upper <- predictions$fit + 1.96 * predictions$se.fit
# 
# # Add these confidence intervals to the new_data dataframe
# new_data$conf_int_lower <- conf_int_lower
# new_data$conf_int_upper <- conf_int_upper
# 
# # Output the new data with risks and confidence intervals
# new_data
# 
# # Horizontal Error Bar Plot (Whisker and Dot Plot)
# ggplot(new_data, aes(y = number, x = risk)) +
#   geom_point(size = 3, color = "blue") +  # Dots for the predicted risks
#   geom_errorbarh(aes(xmin = conf_int_lower, xmax = conf_int_upper), height = 0.2) +  # Horizontal error bars
#   labs(title = "21-Day Mortality Risk by Treatment",
#        x = "Estimated Risk (95% CI)", y = "Treatment") +
#   xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
#   theme_minimal() +  # Minimal theme
#   theme(axis.text.y = element_text(size = 12), 
#         axis.text.x = element_text(size = 12), 
#         axis.title = element_text(size = 14))
# 
