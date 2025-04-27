############################ 
## CRAB ##
############################
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
CR_final$survtime <-ifelse(CR_final$mort_21d_onset4 == 0, 21, CR_final$mortday_onset4)

# Expanding the dataset using the survtime variable
CR_final.surv <- uncount(CR_final, survtime,.remove = F)

# Creating variables for time
CR_final.surv <- CR_final.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()

# Creating variable for timesq
CR_final.surv$timesq <- CR_final.surv$time^2

# Creating event variable
CR_final.surv$event <- ifelse(
  CR_final.surv$time == CR_final.surv$survtime-1 & CR_final.surv$mort_21d_onset4 == 1, 1, 0)

# Pooled logistic regression 
# Fitting pooled logistic regression model with treatment and time (linear and quadratic terms)
fit.pool3 <- glm(event ~ arm + time + timesq,  
                 data = CR_final.surv, 
                 family = "binomial")  # Without product term
summary(fit.pool3)
exp(fit.pool3$coefficients)

# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0 <- data.frame(arm = "Ceftazidime/avibactam based", 
                       time =seq(0,21), timesq=seq(0,21)^2)
results1 <- data.frame(arm = "Polymyxin combination",
                       time =seq(0,21), timesq=seq(0,21)^2)
results2 <- data.frame(arm = "Polymyxin monotherapy",
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
  geom_step(aes(y=risk0, color="Ceftazidime/avibactam based")) +
  geom_step(aes(y=risk1, color="Polymyxin combination")) +
  geom_step(aes(y=risk2, color="Polymyxin monotherapy")) +
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
    
    

# Cox proportional hazards model for comparison

# Fit Cox model
# fit.cox <- coxph(Surv(survtime, event) ~ arm, data = CR_final.surv)
# summary(fit.cox)
# exp(fit.cox$coefficients)

# This is to show that pooled logistic regression can approximate cox.



##########
library(nnet)
# Fit multinomial logistic regression model
# model2 <- multinom (arm ~ age_new + I(age_new*age_new) + sex + 
#                  as.factor(country_income)+ comorbidities_Chalson + I(comorbidities_Chalson* comorbidities_Chalson ) + 
#                  sofa_score_sum  + I(sofa_score_sum * sofa_score_sum) + 
#                  infection_types + delay + I(delay*delay) + icu_at_onset4 + vent_at_onset4 + adm_ward_types_new +
#                  los_onset4 + I(los_onset4* los_onset4) + mono_poly, 
#                data=CR_final)
# summary(model2) # 

# Fit multinomial logistic regression model
model2 <- multinom (arm ~ age_new + sex + country_income + comorbidities_Chalson + 
                      sofa_imp  + infection_types + 
                      hai_icu48days + hai_have_med_device___vent +
                      #icu_at_onset4 + vent_at_onset4+
                      # icu_at_onset5 + vent_at_onset5+
                      + monopoly,
                    data=CR_final)
summary(model2) 

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

### Estimating the causal effect and constructing risk curves for a survival outcome via an IP weighted pooled logistic model
# Converting wide data to long form
CR_final$survtime <-ifelse(CR_final$mort_21d_onset4 == 0, 21, CR_final$mortday_onset4)
# Expanding the dataset using the survtime variable
CR_final.surv <- uncount(CR_final, survtime,.remove = F)
# Creating variables for time
CR_final.surv <- CR_final.surv %>% group_by(recordid) %>% mutate(time=row_number()-1) %>% ungroup()
# Creating variable for timesq
CR_final.surv$timesq <- CR_final.surv$time^2
# Creating event variable
CR_final.surv$event <- ifelse(
  CR_final.surv$time == CR_final.surv$survtime-1 & CR_final.surv$mort_21d_onset4 == 1, 1, 0)

# Fitting a pooled logistic regression model with time, treatment and product term between treatment and time, with nonstabilized weight
fit.pool.w2 <-  glm(event ~ 
                      arm + time + timesq,
                    weights = w2,
                    data = CR_final.surv, 
                    family = "binomial")
summary(fit.pool.w2)
exp(fit.pool.w2$coefficients)
# Computing risks at each time point of follow-up
# Create datasets to store results
# Include all time points under each treatment level
results0.w2 <- data.frame(arm = "Ceftazidime/avibactam based", 
                          time =seq(0,21), timesq=seq(0,21)^2)
results1.w2 <- data.frame(arm = "Polymyxin combination",
                          time =seq(0,21), timesq=seq(0,21)^2)
results2.w2 <- data.frame(arm = "Polymyxin monotherapy",
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
  geom_line(aes(y=risk0, color="Ceftazidime/avibactam based")) +
  geom_line(aes(y=risk1, color="Polymyxin combination")) +
  geom_line(aes(y=risk2, color="Polymyxin monotherapy")) +
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

