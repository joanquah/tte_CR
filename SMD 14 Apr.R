############ Standardized mean difference (CRAB) ###################
# Generate SMD before weighting
bal.tab(arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
          diabetes + malignancy + renal + liver +
          sofa_imp + infection_types + 
          los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
          delay + ent + pae,
        data = crab, 
        estimand = "ATE", 
        weights = NULL)  # No weights for before-weighting SMD

# Generate SMD after weighting (applying inverse probability weights)
bal.tab(arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
          diabetes + malignancy + renal + liver +
          sofa_imp + infection_types + 
          los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
          delay + ent + pae,
        data = crab,
        estimand = "ATE", 
        weights = crab$w2)  # Apply IPW weights

# Love plot for visualization
love.plot(bal.tab(arm ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                    diabetes + malignancy + renal + liver +
                    sofa_imp + infection_types + 
                    los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
                    delay + ent + pae,
                  data = crab, 
                  estimand = "ATE", 
                  weights = crab$w2),
          abs = TRUE, 
          threshold = 0.1, 
          title = "Standardized Mean Differences Before and After Weighting")


####################################################################################################
############ Standardized mean difference (CRE_CRPAE) ###################
# Generate SMD before weighting
bal.tab(arm3 ~ age_new + sex + country_income2 + comorbidities_Chalson + 
          diabetes + malignancy + renal + liver +
          sofa_imp + infection_types + 
          los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
          delay + aci,
        data = cre_crpae, 
        estimand = "ATE", 
        weights = NULL)  # No weights for before-weighting SMD

# Generate SMD after weighting (applying inverse probability weights)
bal.tab(arm3 ~ age_new + sex + country_income2 + comorbidities_Chalson + 
          diabetes + malignancy + renal + liver +
          sofa_imp + infection_types + 
          los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
          delay + aci,
        data = cre_crpae,
        estimand = "ATE", 
        weights = cre_crpae$w2)  # Apply IPW weights

# Love plot for visualization
love.plot(bal.tab(arm3 ~ age_new + sex + country_income2 + comorbidities_Chalson + 
                    diabetes + malignancy + renal + liver +
                    sofa_imp + infection_types + 
                    los_onset4 + hai_icu48days + hai_have_med_device___vent + monopoly +
                    delay + aci,
                  data = cre_crpae, 
                  estimand = "ATE", 
                  weights = cre_crpae$w2),
          abs = TRUE, 
          threshold = 0.1, 
          title = "Standardized Mean Differences Before and After Weighting")

