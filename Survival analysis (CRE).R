# Check analysis
# Plot survival curve

# cre_crpae <- cre_crpae %>% select (recordid, arm3,age_new,sex, country_income, country_income2, 
#                                    comorbidities_Chalson, sofa_imp,qsofa,crea_imp,
#                                    malignancy, diabetes, liver, renal, 
#                                    infection_types,
#                                    fup_day_onset4, mort_21d_onset4,  mort_14d_onset4, mortday_onset4,
#                                    monopoly,aci, pae, ent, delay_group, delay,  
#                                    hai_icu48days, hai_have_med_device___vent, icu_at_onset4,vent_at_onset4, los_onset4
# )

# Crude analysis
# Fit logistic model and estimate 21 day mortality risks
model1 <- glm(mort_21d_onset4~ arm3, data = cre_crpae, family = "binomial")
summary(model1)

# Generate new dataset to store estimated risks
new_data <- data.frame(arm3 = levels(cre_crpae$arm3))

# list unique arm3
levels(cre_crpae$arm3)

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
ggplot(new_data, aes(y = arm3, x = risk)) +
  geom_point(size = 4) +  # Dot for point estimate
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +  # Horizontal error bar
  labs(title = "21-Day Mortality Risk by Treatment",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal()

# Data correction to change 
#For cre_crpae dataframe, assign d28_death_date as 2024-11-25 and fup_day_onset4 as 13 for recordid CN-002-A0-0033
cre_crpae$fup_day_onset4[cre_crpae$recordid == "CN-002-A0-0033"] <- 13
cre_crpae$mortday_onset4[cre_crpae$recordid == "CN-002-A0-0033"] <- 13
cre_crpae$d28_death_date[cre_crpae$recordid == "CN-002-A0-0033"] <- as.Date("2024-11-25")
cre_crpae$fup_day_onset4[cre_crpae$recordid == "TH-005-A0-0123"] <- 8
cre_crpae$mort_date_final[cre_crpae$recordid == "TH-005-A0-0123"]  <- as.Date("2024-06-23")
cre_crpae$d28_date[cre_crpae$recordid == "TH-005-A0-0123"]  <- as.Date("2024-06-23")
cre_crpae$d28_death_date[cre_crpae$recordid == "TH-005-A0-0123"]  <- as.Date("2024-07-09")
cre_crpae$mortday_onset4[cre_crpae$recordid == "TH-005-A0-0123"] <- 12
# remove PK-002-A0-0357 from the dataset
cre_crpae <- cre_crpae[cre_crpae$recordid != "PK-002-A0-0357",]

# Expand dataset based on survival time and fits pooled logistic regression models to estimate risks over time
# Creating survtime variable that represents the duration (in days) of follow-up
cre_crpae$survtime <-ifelse(cre_crpae$mortday_onset4 == 0, 21, cre_crpae$mortday_onset4)

# Expanding the dataset using the survtime variable
cre_crpae.surv <- uncount(cre_crpae, survtime,.remove = F)

# Creating variables for time
cre_crpae.surv <- cre_crpae.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()

# Creating variable for timesq
cre_crpae.surv$timesq <- cre_crpae.surv$time^2

# Creating event variable
cre_crpae.surv$event <- ifelse(
  cre_crpae.surv$time == cre_crpae.surv$survtime-1 & cre_crpae.surv$mort_21d_onset4 == 1, 1, 0)

# Fitting a pooled logistic regression model with time, treatment and product term between treatment and time
#fit.pool1 <-  glm(event ~ arm3 + time + arm3*time, 
#                  data = cre_crpae.surv, 
#                  family = "binomial")
#summary(fit.pool1)
#exp(fit.pool1$coefficients)

# Fitting a pooled logistic regression model with time (linear and quadratic terms), treatment and product terms between treatment and time
#fit.pool2 <-  glm(event ~ arm3 + time + timesq + arm3*time + arm3*timesq, 
#                  data = cre_crpae.surv, 
#                 family = "binomial")
#summary(fit.pool2)
#exp(fit.pool2$coefficients)

# pooled logistic regression 
# Fitting pooled logistic regression model with treatment and time (linear and quadratic terms)
fit.pool3 <- glm(event ~ arm3 + time + timesq,  
                 data = cre_crpae.surv, 
                 family = "binomial")  # Without product term
summary(fit.pool3)
exp(fit.pool3$coefficients)

# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0 <- data.frame(arm3 = "Polymyxin monotherapy", 
                       time =seq(0,21), timesq=seq(0,21)^2)
results1 <- data.frame(arm3 = "Polymyxin combination",
                       time =seq(0,21), timesq=seq(0,21)^2)
results2 <- data.frame(arm3 = "Ceftazidime/avibactam based",
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
  geom_step(aes(y=risk1, color="Polymyxin combination")) +
  geom_step(aes(y=risk2, color="Ceftazidime/avibactam based")) +
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

# Merge confidence intervals from results0, results1, results2
results_combined <- results0 %>%
  select(time, lower_CI_risk0, upper_CI_risk0) %>%
  merge(results1 %>% select(time, lower_CI_risk1, upper_CI_risk1), by = "time") %>%
  merge(results2 %>% select(time, lower_CI_risk2, upper_CI_risk2), by = "time") %>%
  merge(results_combined, by = "time") %>%
  mutate(time_updated = time + 1) %>%
  arrange(time_updated)

ggplot(results_combined, aes(x = time_updated)) + 
  # Ceftazidime/avibactam based (risk0)
  geom_ribbon(aes(ymin = lower_CI_risk0, ymax = upper_CI_risk0, 
                  fill = "#00468B"), alpha = 0.2) +
  geom_step(aes(y = risk0, color = "Ceftazidime/avibactam based"), linewidth = 1) +
  
  # Polymyxin combination (risk1)
  geom_ribbon(aes(ymin = lower_CI_risk1, ymax = upper_CI_risk1, 
                  fill = "#ED0000"), alpha = 0.2) +
  geom_step(aes(y = risk1, color = "Polymyxin combination"), linewidth = 1) +
  
  # Polymyxin monotherapy (risk2)
  geom_ribbon(aes(ymin = lower_CI_risk2, ymax = upper_CI_risk2, 
                  fill = "#42B540"), alpha = 0.2) +
  geom_step(aes(y = risk2, color = "Polymyxin monotherapy"), linewidth = 1) +
  
  # Lancet color scheme
  scale_color_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  scale_fill_manual(
    values = c(
      "#00468B", "#ED0000", "#42B540"
    ),
    guide = "none"  # Hide fill legend (redundant with color)
  ) +
  
  # Axis/labels
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 1)) +
  ylab("Cumulative Incidence of Mortality") +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# First ensure time_updated alignment
# Add time_updated=0 row with NA CIs if missing
results_combined <- results_combined %>% 
  add_row(time_updated = 0, risk0 = 0, risk1 = 0, risk2 = 0,
          lower_CI_risk0 = NA, upper_CI_risk0 = NA,
          lower_CI_risk1 = NA, upper_CI_risk1 = NA,
          lower_CI_risk2 = NA, upper_CI_risk2 = NA) %>% 
  arrange(time_updated)

# Create the final plot
ggplot(results_combined, aes(x = time_updated)) +
  # Ceftazidime/avibactam based
  geom_ribbon(aes(ymin = lower_CI_risk0, ymax = upper_CI_risk0, 
                  fill = "Ceftazidime/avibactam based"), alpha = 0.2) +
  geom_step(aes(y = risk0, color = "Ceftazidime/avibactam based"), linewidth = 1) +
  
  # Polymyxin combination
  geom_ribbon(aes(ymin = lower_CI_risk1, ymax = upper_CI_risk1,
                  fill = "Polymyxin combination"), alpha = 0.2) +
  geom_step(aes(y = risk1, color = "Polymyxin combination"), linewidth = 1) +
  
  # Polymyxin monotherapy
  geom_ribbon(aes(ymin = lower_CI_risk2, ymax = upper_CI_risk2,
                  fill = "Polymyxin monotherapy"), alpha = 0.2) +
  geom_step(aes(y = risk2, color = "Polymyxin monotherapy"), linewidth = 1) +
  
  # Color and fill scales
  scale_color_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  
  # Axis and labels
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 1)) +
  ylab("Cumulative Incidence of Mortality") +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Ensure proper time alignment
  coord_cartesian(xlim = c(0, 21))

###############################################################
# Cox proportional hazards model for comparison

# Fit Cox model
fit.cox <- coxph(Surv(survtime, event) ~ arm3, data = cre_crpae.surv)
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

############################################################################
library(nnet)
# Fit multinomial logistic regression model

model2 <- multinom (arm3 ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                      sofa_imp + infection_types + 
                      los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
                      delay + aci,
                    data=cre_crpae)
summary(model2) 

model2 <- multinom (arm3 ~ age_new + I(age_new*age_new) + sex + 
                      as.factor(country_income)+ comorbidities_Chalson + I(comorbidities_Chalson* comorbidities_Chalson ) + 
                      sofa_score_sum  + I(sofa_score_sum * sofa_score_sum) + 
                      infection_types + delay + I(delay*delay) + icu_at_onset4 + vent_at_onset4 + 
                      los_onset4 + I(los_onset4* los_onset4) + mono_poly, 
                    data=cre_crpae)
summary(model2) # 

model2a <- multinom (arm3 ~ age_new + sex + 
                       as.factor(country_income)+ comorbidities_Chalson  + 
                       sofa_score_sum + infection_types,
                     data=cre_crpae)
summary(model2a) # 

model2b <- multinom (arm3 ~ age_new  + sex + 
                       as.factor(country_income)+ as.factor(country_ab) + comorbidities_Chalson  + 
                       sofa_score_sum + adm_ward_types_new + icu_at_ast + vent_at_ast + 
                       infection_types + delay  + mvdur_ast + 
                       los_ast  + iculos_ast + as.factor(org_names_all),
                     data=cre_crpae)
summary(model2b) #

model2c <- multinom (arm3 ~ age_new + as.factor(country_income) + 
                       sofa_score_sum + infection_types + as.factor(org_names_all),
                     data=cre_crpae)
summary(model2c) #

#Use the predicted probabilities from this model to estimate nonstabilized IP weights
# Generate predicted probabilities for each treatment arm
cre_crpae$prob_matrix <- predict(model2, newdata = cre_crpae, type = "probs") # returns a matrix where each row represents an individual, and each column represents the predicted probability of receiving each treatment.
head(cre_crpae$prob_matrix) 

# Try without this line
# Convert treatment variable to character (to match column names in prob_matrix)
# cre_crpae$arm3 <- as.character(cre_crpae$arm3)

# For each patient, select the probability corresponding to the treatment they actually received.
# Extract the probability corresponding to the observed treatment for each row
cre_crpae$prob <- mapply(function(i, treatment) {
  cre_crpae$prob_matrix[i, treatment]
}, i = seq_len(nrow(cre_crpae)), treatment = cre_crpae$arm3)  #This assigns prob[i] as the probability of the observed treatment for individual i.

# Compute weight (nonstabilized IPW)
cre_crpae$w2 <- 1 / cre_crpae$prob

# Before truncation
hist(cre_crpae$w2)

# Truncating weights at 99th percentile
cre_crpae$w2[cre_crpae$w2 > quantile(cre_crpae$w2, 0.99)] <- quantile(cre_crpae$w2, 0.99)

# Check the distribution of the nonstabilized weights
summary(cre_crpae$w2)  # Check distribution of weights
sd(cre_crpae$w2)

# After truncation
hist(cre_crpae$w2[cre_crpae$w2 <= quantile(cre_crpae$w2, 0.99)])
###
# use MSM with nonstabilized weights
options(warn=-1) # Need to suppress warning or else geeglm will encounter error due to non-integer number of successes as a result of weights
# msm.w <- geeglm(mort_21d_ast ~ arm3, data=cre_crpae, weights=w2, id=recordid, family=binomial())
msm.w <- glm(mort_21d_onset4 ~ arm3, data=cre_crpae, family=binomial(), weights=w2)
summary(msm.w)
exp(coef(msm.w))  # Convert log-odds to odds ratios
exp(confint(msm.w))  # 95% confidence interval for odds ratios

# Outputting risks
new_data <- data.frame(arm3 = levels(cre_crpae$arm3))
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
ggplot(new_data, aes(y = arm3, x = risk)) +
  geom_point(size = 3, color = "blue") +  # Dots for the predicted risks
  geom_errorbarh(aes(xmin = conf_int_lower, xmax = conf_int_upper), height = 0.2) +  # Horizontal error bars
  labs(title = "21-Day Mortality Risk by Treatment",
       x = "Estimated Risk (95% CI)", y = "Treatment") +
  xlim(0, 1) +  # Ensure risk estimates range from 0 to 1
  theme_minimal() +  # Minimal theme
  theme(axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 14))

#################################### Stop here#################################

# Estimation of stabilized ip weights 
# Fit a multinomial logistic regression model for treatment assignment
stab <- multinom(arm3~ 1, data = cre_crpae)
summary(stab)

# Obtain predicted probabilities (p.num) for each treatment arm
cre_crpae$p.num <- predict(stab, newdata = cre_crpae, type = "probs")

# Create stabilized weights for each observation based on the treatment assigned
# For each arm, calculate the stabilized weights based on the formula:
# sw = p(num) / p(predicted) for the arm that was actually received
cre_crpae$sw <- NA
cre_crpae$sw[cre_crpae$arm3 == "Polymyxin monotherapy"] <- cre_crpae$p.num[, "Polymyxin monotherapy"] / cre_crpae$prob[, "Polymyxin monotherapy"]
cre_crpae$sw[cre_crpae$arm3 == "Polymyxin combination"] <- cre_crpae$p.num[, "Polymyxin combination"] / cre_crpae$prob[, "Polymyxin combination"]
cre_crpae$sw[cre_crpae$arm3 == "Ceftazidime/avibactam based"] <- cre_crpae$p.num[, "Ceftazidime/avibactam based"] / cre_crpae$prob[, "Ceftazidime/avibactam based"]

# Check the distribution of stabilized weights
summary(cre_crpae$sw)
sd(cre_crpae$sw)
# As we have seen, stabilizing the inverse probability weights decreases the range of the weight distribution, which in turn results in more efficient estimates from the outcome model in the parametric setting.

# Truncating weights at 99th percentile
cre_crpae$sw[cre_crpae$sw > quantile(cre_crpae$sw, 0.99)] <- quantile(cre_crpae$sw, 0.99)


# MSM with stabilized weights
#msm.sw <- geeglm(mort_21d_ast ~ arm3, data=cre_crpae, weights=sw,id=recordid,family=binomial())
msm.sw <- glm(mort_21d_ast ~ arm3, data=cre_crpae, family=binomial(), weights=sw)
summary(msm.sw)

# Outputting risks
new_data3 <- data.frame(arm3 = levels(cre_crpae$arm3))
new_data3$risk <- predict(msm.sw, newdata=new_data3, type="response")
new_data3

#results compare when using stabilized versus nonstabilized weights are nearly the same but not identical (when unrounded), which is expected of a parametrically estimated outcome model.

################################################################################################################################################
### Estimating the causal effect and constructing risk curves for a survival outcome via an IP weighted pooled logistic model
# Converting wide data to long form
cre_crpae$survtime <-ifelse(cre_crpae$mortday_onset4 == 0, 21, cre_crpae$mortday_onset4)
# Expanding the dataset using the survtime variable
cre_crpae.surv <- uncount(cre_crpae, survtime,.remove = F)
# Creating variables for time
cre_crpae.surv <- cre_crpae.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()
# Creating variable for timesq
cre_crpae.surv$timesq <- cre_crpae.surv$time^2
# Creating event variable
cre_crpae.surv$event <- ifelse(
  cre_crpae.surv$time == cre_crpae.surv$survtime-1 & cre_crpae.surv$mort_21d_onset4 == 1, 1, 0)

# Fitting a pooled logistic regression model with time, treatment and product term between treatment and time, with nonstabilized weight
fit.pool.w2 <-  glm(event ~ 
                      arm3 + time + timesq,
                    weights = w2,
                    data = cre_crpae.surv, 
                    family = "binomial")
summary(fit.pool.w2)
exp(fit.pool.w2$coefficients)
# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0.w2 <- data.frame(arm3 = "Polymyxin monotherapy", 
                          time =seq(0,21), timesq=seq(0,21)^2)
results1.w2 <- data.frame(arm3 = "Polymyxin combination",
                          time =seq(0,21), timesq=seq(0,21)^2)
results2.w2 <- data.frame(arm3 = "Ceftazidime/avibactam based",
                          time =seq(0,21), timesq=seq(0,21)^2)

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
risk0.w2 <-results0.w2[results0.w2$time==21,]$risk0
risk1.w2 <-results1.w2[results1.w2$time==21,]$risk1
risk2.w2 <-results2.w2[results2.w2$time==21,]$risk2
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
  geom_step(aes(y=risk0, color="Polymyxin monotherapy")) +
  geom_step(aes(y=risk1, color="Polymyxin combination")) +
  geom_step(aes(y=risk2, color="Ceftazidime/avibactam based")) +
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks=seq(0,21,1)) +
  ylab("Cumulative Incidence of mortality") +
  labs(colour = "Treatment") +
  theme_bw() + 
  theme(legend.position="bottom")

# To generate confidence intervals for risk estimates at time = 20
# Obtain predicted hazards and standard errors
pred0.w2 <- predict(fit.pool.w2, results0.w2, type = "link", se.fit = TRUE)
pred1.w2 <- predict(fit.pool.w2, results1.w2, type = "link", se.fit = TRUE)
pred2.w2 <- predict(fit.pool.w2, results2.w2, type = "link", se.fit = TRUE)

# Convert logit scale hazards to probabilities
results0.w2$hazard0 <- plogis(pred0.w2$fit)
results1.w2$hazard1 <- plogis(pred1.w2$fit)
results2.w2$hazard2 <- plogis(pred2.w2$fit)

# Compute standard errors of hazards on probability scale using delta method
results0.w2$se_hazard0 <- pred0.w2$se.fit * results0.w2$hazard0 * (1 - results0.w2$hazard0)
results1.w2$se_hazard1 <- pred1.w2$se.fit * results1.w2$hazard1 * (1 - results1.w2$hazard1)
results2.w2$se_hazard2 <- pred2.w2$se.fit * results2.w2$hazard2 * (1 - results2.w2$hazard2)

# Compute cumulative survival probabilities
results0.w2$surv0 <- cumprod(1 - results0.w2$hazard0)
results1.w2$surv1 <- cumprod(1 - results1.w2$hazard1)
results2.w2$surv2 <- cumprod(1 - results2.w2$hazard2)

# Compute cumulative standard errors using Greenwood’s formula
results0.w2$se_surv0 <- results0.w2$surv0 * sqrt(cumsum((results0.w2$se_hazard0 / (1 - results0.w2$hazard0))^2))
results1.w2$se_surv1 <- results1.w2$surv1 * sqrt(cumsum((results1.w2$se_hazard1 / (1 - results1.w2$hazard1))^2))
results2.w2$se_surv2 <- results2.w2$surv2 * sqrt(cumsum((results2.w2$se_hazard2 / (1 - results2.w2$hazard2))^2))

# Compute risk estimates
results0.w2$risk0 <- 1 - results0.w2$surv0
results1.w2$risk1 <- 1 - results1.w2$surv1
results2.w2$risk2 <- 1 - results2.w2$surv2

# Compute confidence intervals for risks
results0.w2$lower_CI_risk0 <- 1 - (results0.w2$surv0 * exp(1.96 * results0.w2$se_surv0 / results0.w2$surv0))
results0.w2$upper_CI_risk0 <- 1 - (results0.w2$surv0 * exp(-1.96 * results0.w2$se_surv0 / results0.w2$surv0))
results1.w2$lower_CI_risk1 <- 1 - (results1.w2$surv1 * exp(1.96 * results1.w2$se_surv1 / results1.w2$surv1))
results1.w2$upper_CI_risk1 <- 1 - (results1.w2$surv1 * exp(-1.96 * results1.w2$se_surv1 / results1.w2$surv1))
results2.w2$lower_CI_risk2 <- 1 - (results2.w2$surv2 * exp(1.96 * results2.w2$se_surv2 / results2.w2$surv2))
results2.w2$upper_CI_risk2 <- 1 - (results2.w2$surv2 * exp(-1.96 * results2.w2$se_surv2 / results2.w2$surv2))

# Merge confidence intervals from results0, results1, results2
results_combined <- results0.w2 %>%
  select(time, lower_CI_risk0, upper_CI_risk0) %>%
  merge(results1.w2 %>% select(time, lower_CI_risk1, upper_CI_risk1), by = "time") %>%
  merge(results2.w2 %>% select(time, lower_CI_risk2, upper_CI_risk2), by = "time") %>%
  merge(results_combined, by = "time") %>%
  mutate(time_updated = time + 1) %>%
  arrange(time_updated)

# Create the final plot
ggplot(results_combined, aes(x = time_updated)) +
  # Ceftazidime/avibactam based
  geom_ribbon(aes(ymin = lower_CI_risk0, ymax = upper_CI_risk0, 
                  fill = "Polymyxin monotherapy"), alpha = 0.2) +
  geom_step(aes(y = risk0, color = "Polymyxin monotherapy"), linewidth = 1) +
  
  # Polymyxin combination
  geom_ribbon(aes(ymin = lower_CI_risk1, ymax = upper_CI_risk1,
                  fill = "Polymyxin combination"), alpha = 0.2) +
  geom_step(aes(y = risk1, color = "Polymyxin combination"), linewidth = 1) +
  
  # Polymyxin monotherapy
  geom_ribbon(aes(ymin = lower_CI_risk2, ymax = upper_CI_risk2,
                  fill = "Ceftazidime/avibactam based"), alpha = 0.2) +
  geom_step(aes(y = risk2, color = "Ceftazidime/avibactam based"), linewidth = 1) +
  
  # Color and fill scales
  scale_color_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  
  # Axis and labels
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 1)) +
  ylab("Cumulative Incidence of Mortality") +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# First ensure time_updated alignment
# Add time_updated=0 row with NA CIs if missing
results_combined <- results_combined %>% 
  add_row(time_updated = 0, risk0 = 0, risk1 = 0, risk2 = 0,
          lower_CI_risk0 = NA, upper_CI_risk0 = NA,
          lower_CI_risk1 = NA, upper_CI_risk1 = NA,
          lower_CI_risk2 = NA, upper_CI_risk2 = NA) %>% 
  arrange(time_updated)

# Create the final plot
ggplot(results_combined, aes(x = time_updated)) +
  # Ceftazidime/avibactam based
  geom_ribbon(aes(ymin = lower_CI_risk0, ymax = upper_CI_risk0, 
                  fill = "Polymyxin monotherapy"), alpha = 0.2) +
  geom_step(aes(y = risk0, color = "Polymyxin monotherapy"), linewidth = 1) +
  
  # Polymyxin combination
  geom_ribbon(aes(ymin = lower_CI_risk1, ymax = upper_CI_risk1,
                  fill = "Polymyxin combination"), alpha = 0.2) +
  geom_step(aes(y = risk1, color = "Polymyxin combination"), linewidth = 1) +
  
  # Polymyxin monotherapy
  geom_ribbon(aes(ymin = lower_CI_risk2, ymax = upper_CI_risk2,
                  fill = "Ceftazidime/avibactam based"), alpha = 0.2) +
  geom_step(aes(y = risk2, color = "Ceftazidime/avibactam based"), linewidth = 1) +
  
  # Color and fill scales
  scale_color_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c(
      "Ceftazidime/avibactam based" = "#00468B",
      "Polymyxin combination" = "#ED0000",
      "Polymyxin monotherapy" = "#42B540"
    )
  ) +
  
  # Axis and labels
  xlab("Time (days)") +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 1)) +
  ylab("Cumulative Incidence of Mortality") +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.8, "cm")
  ) +
  
  # Ensure proper time alignment
  coord_cartesian(xlim = c(0, 21))


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


