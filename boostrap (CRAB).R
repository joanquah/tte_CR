library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggpubr)

# Bootstrap

################################################
### Obtain percentile-based bootstrapped 95% CIs for each quantity ###

# Resampling individuals
# Extract individual recordid, resampled into bootstrap samples
crab_ids <-data.frame(id = unique(crab.surv$recordid))

# Create a function to obtain risks, RD, and RR from each bootstrap sample
risk.boot <- function(data, indices) {
  # Select individuals into each boostrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the boostrapped sample
  d <- left_join(boot.ids, crab.surv, by = c("id" = "recordid"))
  # At this point, d is the person-time dataset limited to the selected individuals
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      arm + time + timesq,
                    weights = d$w2,
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  #Include all time points under each treatment level
  results0 <- data.frame(arm = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(arm = "Polymyxin combination",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results2 <- data.frame(arm = "Polymyxin Sulbactam",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results3 <- data.frame(arm = "Sulbactam based",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("arm", "time", "timesq")
  colnames(results1) <- c("arm", "time", "timesq")
  colnames(results2) <- c("arm", "time", "timesq")
  colnames(results3) <- c("arm", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  #Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  results2$p.event2 <- predict(pool.boot, results2, type="response")
  results3$p.event3 <- predict(pool.boot, results3, type="response")
  
  # Estimate survival probabilityes S(t) from hazads, h(t)
  # S (t) = cumulative product of (1 - h(t))
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  results2$surv2 <- cumprod(1-results2$p.event2)
  results3$surv3 <- cumprod(1-results3$p.event3)
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  results2$risk2 <- 1 - results2$surv2
  results3$risk3 <- 1 - results3$surv3
  
  # Merge data from four groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- merge(graph, results2, by=c("time", "timesq"))
  graph <- merge(graph, results3, by=c("time", "timesq"))
  return(c(graph$risk0[which(graph$time==20)], 
           graph$risk1[which(graph$time==20)], 
           graph$risk2[which(graph$time==20)],
           graph$risk3[which(graph$time==20)]
           ))
}

# Run 2 boostrap samples 
library(boot)
set.seed(123)
risk.results <- boot(data = crab_ids, statistic = risk.boot, R = 50)

# Print point estimates from the original data
head(risk.results$t0)


# 95% CI for risk in polmyxin monotherapy group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =1)

# 95% CI for risk in polmyxin combination group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =2)

# 95% CI for risk in polymyxin sulbactam group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =3)

# 95% CI for risk in sulbactam based group
boot.ci(risk.results,
        conf = 0.95,
        type = "perc",
        index =4)

# # 95% CI for risk difference
# boot.ci(risk.results,
#         conf = 0.95,
#         type = "perc",
#         index =3)
# 
# # 95% CI for risk ratio
# boot.ci(risk.results,
#         conf = 0.95,
#         type = "perc",
#         index =4)


### Create a parametric cumulative incidnce (risk) plot that includes 95% CIs

# Estimate CIs for plot using bootstrapping

### Polymyxin monotherapy group

# Create a function to obtain risk in polymyxin monotherapy group at each time t from each bootstrap sample
risk.boot.0 <- function(data, indices) {
  #Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, crab.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      arm + time + timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(arm = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(arm = "Polymyxin combination",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results2 <- data.frame(arm = "Polymyxin Sulbactam",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results3 <- data.frame(arm = "Sulbactam based",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("arm", "time", "timesq")
  colnames(results1) <- c("arm", "time", "timesq")
  colnames(results2) <- c("arm", "time", "timesq")
  colnames(results3) <- c("arm", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  results2$p.event2 <- predict(pool.boot, results2, type="response")
  results3$p.event3 <- predict(pool.boot, results3, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  results2$surv2 <- cumprod(1-results2$p.event2)
  results3$surv3 <- cumprod(1-results3$p.event3)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  results2$risk2 <- 1 - results2$surv2
  results3$risk3 <- 1 - results3$surv3
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- merge(graph, results2, by=c("time", "timesq"))
  graph <- merge(graph, results3, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk0)
}

# Run 2 boostrap samples
set.seed(123)
risk.results0 <- boot(data = crab_ids, statistic = risk.boot.0, R = 50)

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
  d <- left_join(boot.ids, crab.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      arm + time + timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(arm = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(arm = "Polymyxin combination",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results2 <- data.frame(arm = "Polymyxin Sulbactam",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results3 <- data.frame(arm = "Sulbactam based",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("arm", "time", "timesq")
  colnames(results1) <- c("arm", "time", "timesq")
  colnames(results2) <- c("arm", "time", "timesq")
  colnames(results3) <- c("arm", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  results2$p.event2 <- predict(pool.boot, results2, type="response")
  results3$p.event3 <- predict(pool.boot, results3, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  results2$surv2 <- cumprod(1-results2$p.event2)
  results3$surv3 <- cumprod(1-results3$p.event3)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  results2$risk2 <- 1 - results2$surv2
  results3$risk3 <- 1 - results3$surv3
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- merge(graph, results2, by=c("time", "timesq"))
  graph <- merge(graph, results3, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk1)
}

# Run 2 boostrap samples
set.seed(123)
risk.results1 <- boot(data = crab_ids, statistic = risk.boot.1, R = 50)

# Combine relevant boostrapped results into a dataframe
risk.boot.results.1 <- data.frame(cbind(risk1 = risk.results1$t0,
                                        t(risk.results1$t)))
# Format boostrapped results for plotting
risk.boot.graph.1 <- data.frame (cbind(time = seq(0,20),
                                       mean.1 = risk.boot.results.1$risk1),
                                 ll.1 = (apply((risk.boot.results.1)[,-1], 1, quantile, probs = 0.025)),
                                 ul.1 = (apply((risk.boot.results.1)[,-1], 1, quantile, probs = 0.975)))

### Polymyxin Sulbactam group
# Create a function to obtain risk in polymyxin sulbactam group at each time t from each bootstrap sample
risk.boot.2 <- function(data, indices) {
  #Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, crab.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      arm + time + timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(arm = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(arm = "Polymyxin combination",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results2 <- data.frame(arm = "Polymyxin Sulbactam",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results3 <- data.frame(arm = "Sulbactam based",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("arm", "time", "timesq")
  colnames(results1) <- c("arm", "time", "timesq")
  colnames(results2) <- c("arm", "time", "timesq")
  colnames(results3) <- c("arm", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  results2$p.event2 <- predict(pool.boot, results2, type="response")
  results3$p.event3 <- predict(pool.boot, results3, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  results2$surv2 <- cumprod(1-results2$p.event2)
  results3$surv3 <- cumprod(1-results3$p.event3)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  results2$risk2 <- 1 - results2$surv2
  results3$risk3 <- 1 - results3$surv3
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- merge(graph, results2, by=c("time", "timesq"))
  graph <- merge(graph, results3, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk2)
}

# Run 2 boostrap samples
set.seed(123)
risk.results2 <- boot(data = crab_ids, statistic = risk.boot.2, R = 50)
# Combine relevant boostrapped results into a dataframe
risk.boot.results.2 <- data.frame(cbind(risk2 = risk.results2$t0,
                                        t(risk.results2$t)))
# Format boostrapped results for plotting
risk.boot.graph.2 <- data.frame (cbind(time = seq(0,20),
                                       mean.2 = risk.boot.results.2$risk2),
                                 ll.2 = (apply((risk.boot.results.2)[,-1], 1, quantile, probs = 0.025)),
                                 ul.2 = (apply((risk.boot.results.2)[,-1], 1, quantile, probs = 0.975)))

### Sulbactam based group
# Create a function to obtain risk in sulbactam based group at each time t from each bootstrap sample
risk.boot.3 <- function(data, indices) {
  #Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <-data.frame(id = ids[indices])
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, crab.surv, by = c("id" = "recordid"))
  
  # Fit pooled logistic model to estimate discrete hazards
  pool.boot <- glm (event ~ 
                      arm + time + timesq, 
                    data = d, 
                    family = "binomial")
  
  # Create a dataset to store results
  # Include all time points under each treatment level
  results0 <- data.frame(arm = "Polymyxin monotherapy", 
                         time =seq(0,20), timesq=seq(0,20)^2)
  results1 <- data.frame(arm = "Polymyxin combination",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results2 <- data.frame(arm = "Polymyxin Sulbactam",
                         time =seq(0,20), timesq=seq(0,20)^2)
  results3 <- data.frame(arm = "Sulbactam based",
                         time =seq(0,20), timesq=seq(0,20)^2)
  
  # Set column names
  colnames(results0) <- c("arm", "time", "timesq")
  colnames(results1) <- c("arm", "time", "timesq")
  colnames(results2) <- c("arm", "time", "timesq")
  colnames(results3) <- c("arm", "time", "timesq")
  
  # Extract predicted values from pooled logistic regression model
  # Predicted values correspond to discrete-time hazards
  results0$p.event0  <- predict(pool.boot, results0, type="response")
  results1$p.event1 <- predict(pool.boot, results1, type="response")
  results2$p.event2 <- predict(pool.boot, results2, type="response")
  results3$p.event3 <- predict(pool.boot, results3, type="response")
  
  # Convert from discrete-time hazards to survival probabilities
  results0$surv0 <- cumprod(1-results0$p.event0)
  results1$surv1 <- cumprod(1-results1$p.event1)
  results2$surv2 <- cumprod(1-results2$p.event2)
  results3$surv3 <- cumprod(1-results3$p.event3)
  
  # Convert from survival probabilities to risks
  results0$risk0 <- 1 - results0$surv0
  results1$risk1 <- 1 - results1$surv1
  results2$risk2 <- 1 - results2$surv2
  results3$risk3 <- 1 - results3$surv3
  
  # Merge data from two groups and format
  graph <- merge(results0, results1, by=c("time", "timesq"))
  graph <- merge(graph, results2, by=c("time", "timesq"))
  graph <- merge(graph, results3, by=c("time", "timesq"))
  graph <- graph[order(graph$time),]
  return(graph$risk3)
}
  
# Run 2 boostrap samples
set.seed(123)
risk.results3 <- boot(data = crab_ids, statistic = risk.boot.3, R = 50)
# Combine relevant boostrapped results into a dataframe
risk.boot.results.3 <- data.frame(cbind(risk3 = risk.results3$t0,
                                        t(risk.results3$t)))
# Format boostrapped results for plotting
risk.boot.graph.3 <- data.frame (cbind(time = seq(0,20),
                                       mean.3 = risk.boot.results.3$risk3),
                                 ll.3 = (apply((risk.boot.results.3)[,-1], 1, quantile, probs = 0.025)),
                                 ul.3 = (apply((risk.boot.results.3)[,-1], 1, quantile, probs = 0.975)))

### Merge all risk estimates into one data frame

# Prepare data
risk.boot.graph.pred <- merge(risk.boot.graph.0, risk.boot.graph.1, by = "time")
risk.boot.graph.pred <- merge(risk.boot.graph.pred, risk.boot.graph.2, by = "time")
risk.boot.graph.pred <- merge(risk.boot.graph.pred, risk.boot.graph.3, by = "time")


# Edit data frame to reflect that risks are estimated at the END of each interval
risk.boot.graph.pred$time_0 <- risk.boot.graph.pred$time + 1
zero <-data.frame(cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
zero <-setNames(zero, names(risk.boot.graph.pred))
risk.boot.graph <- rbind(zero, risk.boot.graph.pred)

head(risk.boot.graph)

# Use color-blind palette (Color Universal Design)
line_colors <- c(
  "Polymyxin monotherapy" = "#E69F00",
  "Polymyxin combination" = "#56B4E9",
  "Polymyxin Sulbactam" = "#009E73",
  "Sulbactam based" = "#D55E00"
)

# STEP 1: Plot with legend for point estimates only
plot.plr.ci <- ggplot(risk.boot.graph, aes(x = time_0)) +
  geom_step(aes(y = mean.0, color = "Polymyxin monotherapy", linetype = "a"), size = 1.2, direction = "hv") +
  geom_step(aes(y = mean.1, color = "Polymyxin combination", linetype = "a"), size = 1.2, direction = "hv") +
  geom_step(aes(y = mean.2, color = "Polymyxin Sulbactam", linetype = "a"), size = 1.2, direction = "hv") +
  geom_step(aes(y = mean.3, color = "Sulbactam based", linetype = "a"), size = 1.2, direction = "hv") +
  
  geom_ribbon(aes(ymin = ll.0, ymax = ul.0, fill = "Polymyxin monotherapy"), alpha = 0.2, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ll.1, ymax = ul.1, fill = "Polymyxin combination"), alpha = 0.2, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ll.2, ymax = ul.2, fill = "Polymyxin Sulbactam"), alpha = 0.2, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ll.3, ymax = ul.3, fill = "Sulbactam based"), alpha = 0.2, show.legend = FALSE) +
  
  scale_color_manual(name = "Strata", values = line_colors) +
  scale_fill_manual(values = line_colors) +
  scale_linetype_manual(values = "solid", guide = "none") +
  
  xlab("Time (days)") +
  ylab("Cumulative incidence of mortality (%)") +
  scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 7)) +
  # scale_y_continuous(name = "Cumulative incidence of mortality", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  
  labs(title = "All-cause mortality") +
  
  theme_minimal(base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, face = "bold", family = "Times"),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 10, family = "Times"),
    axis.title = element_text(size = 10, family = "Times"),
    legend.title = element_text(size = 10, family = "Times"),
    legend.text = element_text(size = 10, family = "Times"),
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA))
  

ggsave("figure1_lancetstyle.tiff", plot = plot.plr.ci,
       width = 107, height = 90, units = "mm", dpi = 300, compression = "lzw")

library(cowplot)

# STEP 2: Prepare number-at-risk data for display at day 0, 7, 14, 21
# --- Risk/event table calculation ---
# Adjust time values
# crab.surv$time_updated <- crab.surv$time + 1
# 
# time_points <- c(0, 7, 14, 21)
# 
# risk_event_table <- lapply(time_points, function(t) {
#   crab.surv %>%
#     filter(time_updated <= t) %>%
#     group_by(recordid, arm) %>%
#     summarise(
#       has_event = any(event == 1),
#       first_event_time = ifelse(any(event == 1), min(time_updated[event == 1]), Inf),
#       .groups = "drop"
#     ) %>%
#     mutate(
#       at_risk = first_event_time > t,
#       event_by_t = ifelse(first_event_time <= t, 1, 0)
#     ) %>%
#     group_by(arm) %>%
#     summarise(
#       day = t,
#       at_risk = sum(at_risk),
#       cum_events = sum(event_by_t),
#       .groups = "drop"
#     )
# }) %>%
#   bind_rows()

# Step 1: Add time_updated
crab.surv <- crab.surv %>%
  mutate(time_updated = time + 1)

# Step 2: Calculate first event time per recordid
subject_event_info <- crab.surv %>%
  group_by(recordid, arm) %>%
  summarise(
    first_event_time = if (any(event == 1)) min(time_updated[event == 1]) else Inf,
    .groups = "drop"
  )

# Step 3: Compute risk and event counts at each timepoint
time_points <- c(0, 7, 14, 21)

risk_event_table <- lapply(time_points, function(t) {
  subject_event_info %>%
    mutate(
      at_risk = first_event_time > t,
      event_by_t = ifelse(first_event_time <= t, 1, 0)
    ) %>%
    group_by(arm) %>%
    summarise(
      day = t,
      at_risk = sum(at_risk),
      cum_events = sum(event_by_t),
      .groups = "drop"
    )
}) %>%
  bind_rows()

# Pivot table to long-wide format
table_plot_ready <- risk_event_table %>%
  pivot_wider(names_from = day, values_from = c(at_risk, cum_events)) %>%
  arrange(arm)

# Create separate rows for "Number at risk" and "Cumulative events"
risk_row <- table_plot_ready %>%
  select(arm, starts_with("at_risk")) %>%
  mutate(type = "Number at risk")

event_row <- table_plot_ready %>%
  select(arm, starts_with("cum_events")) %>%
  mutate(type = "Cumulative events")

# Combine into a table
final_table <- bind_rows(risk_row, event_row) %>%
  relocate(type, .before = arm)

# Clean column names
colnames(final_table) <- gsub(".*_(\\d+)", "\\1", colnames(final_table))
colnames(final_table)[1:2] <- c("", "Treatment")

# Convert to grob
table_grob <- tableGrob(
  final_table,
  rows = NULL,
  theme = ttheme_default(
    core = list(fg_params = list(fontfamily = "Times", fontsize = 10)),
    colhead = list(fg_params = list(fontfamily = "Times", fontface = "bold", fontsize = 10))
  )
)
plot.plr.ci <- plot.plr.ci +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha('white', 0.6), color = NA)
  )


final_plot <- ggarrange(
  plot.plr.ci,
  table_grob,
  ncol = 1,
  heights = c(3, 0.8)
)

print(final_plot)


ggsave("figure1_lancetstyle_withrisk.tiff", plot = final_plot,
       width = 107, height = 110, units = "mm", dpi = 300, compression = "lzw")
###########################
test <- risk.boot.graph %>%
  mutate(valid_CI = ll.1 <= mean.1 & mean.1 <= ul.1) %>%
  filter(!valid_CI)
############################
# STEP 2: Reformat table
# Create wide table for plotting
table_matrix <- matrix(
  nrow = 12, ncol = 5,
  dimnames = list(NULL, c("", "0", "7", "14", "21"))
)

# Row 1: label
table_matrix[1, 1] <- "Number at risk"

# Row 3–6: treatments (number at risk)
table_matrix[2, ] <- c("Polymyxin monotherapy", 
                       table_plot_ready %>% filter(arm == "Polymyxin monotherapy") %>% select(starts_with("at_risk")) %>% as.character())
table_matrix[3, ] <- c("Polymyxin combination", 
                       table_plot_ready %>% filter(arm == "Polymyxin combination") %>% select(starts_with("at_risk")) %>% as.character())
table_matrix[4, ] <- c("Polymyxin Sulbactam", 
                       table_plot_ready %>% filter(arm == "Polymyxin Sulbactam") %>% select(starts_with("at_risk")) %>% as.character())
table_matrix[5, ] <- c("Sulbactam based", 
                       table_plot_ready %>% filter(arm == "Sulbactam based") %>% select(starts_with("at_risk")) %>% as.character())

# Row 7: label
table_matrix[6, 1] <- "Cumulative number of events"

# Row 8–11: treatments (events)
table_matrix[7, ] <- c("Polymyxin monotherapy", 
                       table_plot_ready %>% filter(arm == "Polymyxin monotherapy") %>% select(starts_with("cum_events")) %>% as.character())
table_matrix[8, ] <- c("Polymyxin combination", 
                       table_plot_ready %>% filter(arm == "Polymyxin combination") %>% select(starts_with("cum_events")) %>% as.character())
table_matrix[9, ] <- c("Polymyxin Sulbactam", 
                        table_plot_ready %>% filter(arm == "Polymyxin Sulbactam") %>% select(starts_with("cum_events")) %>% as.character())
table_matrix[10, ] <- c("Sulbactam based", 
                        table_plot_ready %>% filter(arm == "Sulbactam based") %>% select(starts_with("cum_events")) %>% as.character())

# Convert to tableGrob
table_grob <- tableGrob(
  table_matrix,
  rows = NULL,
  theme = ttheme_default(
    core = list(fg_params = list(fontfamily = "Times", fontsize = 10)),
    colhead = list(fg_params = list(fontfamily = "Times", fontface = "bold", fontsize = 10))
  )
)

# Step 3: Final assembly 
spacer <- ggplot() + theme_void()  # creates an empty space

final_plot <- ggarrange(
  plot.plr.ci,
  spacer,
  table_grob,
  ncol = 1,
  heights = c(4, 0.3, 2)  # tweak the middle value to control spacing
)

print(final_plot)
