# Baseline table
# Stratify age 
CR_final$age_group <- cut(CR_final$age_new, breaks = c(17, 40, 60, 80, 100), labels = c("18-40", "41-60", "61-80", "81-100"))

# Baseline demographic table with p-value, for abx_outcome, intervention is Arm, outcome is mortality
table1 <- CR_final %>% 
  #filter(monopoly == "monomicrobial") %>%
  select(arm, 
         age_group,
         age_new, sex, country_income,  
         comorbidities_Chalson,
         malignancy,
         diabetes,
         liver,
         renal,
         cmb_comorbidities___onc,
         cmb_comorbidities___mst, 
         cmb_comorbidities___cpd,   
         cmb_comorbidities___cog,
         cmb_comorbidities___rheu,
         cmb_comorbidities___diab, 
         cmb_comorbidities___diad,  
         cmb_comorbidities___liv,
         cmb_comorbidities___renal,
         sofa_imp, qsofa,
         infection_types, monopoly, 
         delay,delay_group,
         icu_at_onset4, vent_at_onset4, los_onset4, iculos_onset4, mvdur_onset4, 
         hai_icu48days, hai_have_med_device___vent,
         crea_imp,
         aci, pae, ent, cr_aci, cr_pae, cr_ent)%>% 
  #sofa_score_sum_group,
  tbl_summary(
    by = arm,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})", # Show median (IQR) instead of mean (SD)
      # add mean and sd
      all_continuous() ~ "{mean} ({sd})", # Show mean (SD)
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany" # Show missing data only if present
  ) %>%
  add_overall() 

# %>%
#   add_p(test = list(
#     all_continuous() ~ "kruskal.test",  # Kruskal-Wallis test for continuous variables
#     all_categorical() ~ "chisq.test" # Chi-squared test for categorical variables
#   ))

print(table1)


library(gtsummary)
library(dplyr)

# Pre-process data
table1_data <- CR_final %>% 
  mutate(
    # Create new variables
    malignancy = ifelse(cmb_comorbidities___onc == 1 | cmb_comorbidities___mst == 1, 1, 0),
    diabetes = ifelse(cmb_comorbidities___diab == 1 | cmb_comorbidities___diad == 1, 1, 0),
    sofa_group = case_when(
      sofa_imp <= 2 ~ "0-2",
      sofa_imp <= 5 ~ "3-5",
      sofa_imp >= 6 ~ "6 and above"
    ),
    sex = factor(sex, levels = c("M", "F"), labels = c("Male", "Female"))
  ) %>%
  select(
    arm, 
    age_new, sex, country_income2,
    comorbidities_Chalson,
    malignancy, diabetes,
    cmb_comorbidities___cpd,
    cmb_comorbidities___cog,
    cmb_comorbidities___rheu,
    cmb_comorbidities___liv,
    cmb_comorbidities___renal,
    sofa_imp, sofa_group, qsofa,
    infection_types, monopoly,
    delay, delay_group,
    icu_at_onset4, vent_at_onset4,
    los_onset4, iculos_onset4, mvdur_onset4,
    hai_icu48days, hai_have_med_device___vent,
    crea_imp, aci, cr_aci,pae, cr_pae, ent, cr_ent
  )

# Create Lancet-style table
lancet_table <- table1_data %>%
  tbl_summary(
    by = arm,
    type = list(
      c(comorbidities_Chalson, sofa_imp, qsofa, delay, los_onset4, iculos_onset4, mvdur_onset4, crea_imp) ~ "continuous",
      c(malignancy, diabetes, cmb_comorbidities___cpd, cmb_comorbidities___cog, cmb_comorbidities___rheu, 
        cmb_comorbidities___liv, cmb_comorbidities___renal, delay_group, icu_at_onset4, vent_at_onset4,
        hai_have_med_device___vent,hai_icu48days,aci, pae, ent, cr_aci, cr_pae, cr_ent) ~ "dichotomous"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age_new ~ "Age, years",
      sex ~ "Sex",
      country_income2 ~ "Country income status",
      comorbidities_Chalson ~ "Charlson Comorbidity Index",
      malignancy ~ "Malignancy",
      diabetes ~ "Diabetes",
      cmb_comorbidities___cpd ~ "Chronic pulmonary disease",
      cmb_comorbidities___cog ~ "Congestive heart failure",
      cmb_comorbidities___rheu ~ "Connective tissue disease",
      cmb_comorbidities___liv ~ "Moderate to severe liver disease",
      cmb_comorbidities___renal ~ "Renal disease",
      sofa_imp ~ "SOFA score",
      sofa_group ~ "SOFA category",
      qsofa ~ "qSOFA score",
      infection_types ~ "Infection type",
      monopoly ~ "Number of pathogen isolated in index culture",
      delay ~ "Days to effective antibiotics",
      delay_group ~ "Delayed antibiotics (>1 day)",
      icu_at_onset4 ~ "ICU admission at index date",
      vent_at_onset4 ~ "Ventilatory support at index date",
      los_onset4 ~ "Hospital length of stay prior to infection, days",
      iculos_onset4 ~ "ICU length of stay prior to infection, days",
      mvdur_onset4 ~ "Mechanical ventilation duration prior to infection, days",
      hai_icu48days ~ "ICU admission >48h before infection",
      hai_have_med_device___vent ~ "Ventilatory support at infection diagnosis",
      crea_imp ~ "Serum creatinine at infection onset, mg/dL",
      aci ~ "A. calcoaceticus-baumannii complex",
      cr_aci ~ "Carbapenem-resistant A. calcoaceticus-baumannii",
      ent ~ "Enterobacterales",
      cr_ent ~ "Carbapenem-resistant Enterobacterales",
      pae ~ "P. aeruginosa",
      cr_pae ~ "Carbapenem-resistant P. aeruginosa"
    ),
    missing = "ifany",
    digits = list(all_continuous() ~ c(1, 1, 1))
  ) %>%
  modify_header(
    label ~ "**Characteristic**",
    all_stat_cols() ~ "**{level}**\n(N = {n})"
  ) %>%
  add_overall(
    last = TRUE,
    col_label = "**Overall**\n(N = {N})"
  ) %>%
  # add_p(
  #   test = list(
  #     all_continuous() ~ "kruskal.test",
  #     all_categorical() ~ "fisher.test"
  #   ),
  #   pvalue_fun = function(x) style_pvalue(x, digits = 2)
  # ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Treatment Group**") %>%
  bold_labels() %>%
  as_gt() %>%
  # Lancet-specific styling
  gt::tab_options(
    table.font.names = "Arial",
    table.font.size = "90%",
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    heading.border.bottom.style = "solid",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    row_group.border.bottom.style = "none"
  ) %>%
  gt::tab_style(
    style = cell_text(size = "small"),
    locations = cells_body()
  ) %>%
  gt::tab_footnote(
    footnote = "Data are median (IQR) for continuous variables and n (%) for categorical variables",
    locations = cells_column_labels(columns = 1)
  )

# Print table
lancet_table

###### Another format ###
# Pre-process data
table1_data <- CR_final %>% 
  mutate(
    # Create new variables
    malignancy = ifelse(cmb_comorbidities___onc == 1 | cmb_comorbidities___mst == 1, 1, 0),
    diabetes = ifelse(cmb_comorbidities___diab == 1 | cmb_comorbidities___diad == 1, 1, 0),
    sofa_group = case_when(
      sofa_imp <= 2 ~ "0-2",
      sofa_imp <= 5 ~ "3-5",
      sofa_imp >= 6 ~ "6 and above"
    ),
    sex = factor(sex, levels = c("M", "F"), labels = c("Male", "Female")),
    delay_cat = factor(delay, levels = 0:4)  # NEW: Factor for delay 0–4
  ) %>%
  select(
    arm, 
    age_new, sex, country, country_income, country_income2,
    comorbidities_Chalson,
    malignancy, diabetes,
    cmb_comorbidities___cpd,
    cmb_comorbidities___cog,
    cmb_comorbidities___rheu,
    cmb_comorbidities___liv,
    cmb_comorbidities___renal,
    sofa_imp, sofa_group, qsofa,
    infection_types, monopoly,
    delay_cat, delay_group,
    icu_at_onset4, vent_at_onset4,
    los_onset4, iculos_onset4, mvdur_onset4,
    hai_icu48days, hai_have_med_device___vent,
    crea_imp,
    aci, cr_aci, ent, cr_ent, pae, cr_pae  # reordered isolated pathogen vars
  )

# Create Lancet-style table
lancet_table <- table1_data %>%
  tbl_summary(
    by = arm,
    type = list(
      c(comorbidities_Chalson, sofa_imp, qsofa, los_onset4, iculos_onset4, mvdur_onset4, crea_imp) ~ "continuous",
      c(malignancy, diabetes, cmb_comorbidities___cpd, cmb_comorbidities___cog, cmb_comorbidities___rheu, 
        cmb_comorbidities___liv, cmb_comorbidities___renal, delay_group, delay_cat,
        icu_at_onset4, vent_at_onset4, hai_have_med_device___vent, hai_icu48days,
        aci, cr_aci, ent, cr_ent, pae, cr_pae) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age_new ~ "Age, years",
      sex ~ "Sex",
      country_income2 ~ "Country income status",
      comorbidities_Chalson ~ "Charlson Comorbidity Index",
      malignancy ~ "Malignancy",
      diabetes ~ "Diabetes",
      cmb_comorbidities___cpd ~ "Chronic pulmonary disease",
      cmb_comorbidities___cog ~ "Congestive heart failure",
      cmb_comorbidities___rheu ~ "Connective tissue disease",
      cmb_comorbidities___liv ~ "Moderate to severe liver disease",
      cmb_comorbidities___renal ~ "Renal disease",
      sofa_imp ~ "SOFA score",
      sofa_group ~ "SOFA category",
      qsofa ~ "qSOFA score",
      infection_types ~ "Infection type",
      monopoly ~ "Number of pathogens isolated in index culture",
      delay_cat ~ "Days to effective antibiotics (0–4)",
      delay_group ~ "Delayed antibiotics (>1 day)",
      icu_at_onset4 ~ "ICU admission at infection onset",
      vent_at_onset4 ~ "Ventilatory support at infection onset",
      los_onset4 ~ "Hospital length of stay prior to infection, days",
      iculos_onset4 ~ "ICU length of stay prior to infection, days",
      mvdur_onset4 ~ "Mechanical ventilation duration prior to infection, days",
      hai_icu48days ~ "ICU admission >48h before infection",
      hai_have_med_device___vent ~ "Ventilatory support at infection diagnosis",
      crea_imp ~ "Serum creatinine at infection onset, mg/dL",
      aci ~ "A. calcoaceticus-baumannii complex",
      cr_aci ~ "Carbapenem-resistant A. calcoaceticus-baumannii",
      ent ~ "Enterobacterales",
      cr_ent ~ "Carbapenem-resistant Enterobacterales",
      pae ~ "P. aeruginosa",
      cr_pae ~ "Carbapenem-resistant P. aeruginosa"
    ),
    missing = "ifany",
    digits = list(all_continuous() ~ c(1, 1, 1))
  ) %>%
  modify_header(
    label ~ "**Characteristic**",
    all_stat_cols() ~ "**{level}**\n(N = {n})"
  ) %>%
  add_overall(
    last = TRUE,
    col_label = "**Overall**\n(N = {N})"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Treatment Group**") %>%
  bold_labels() %>%
  as_gt() %>%
  gt::tab_options(
    table.font.names = "Arial",
    table.font.size = "90%",
    column_labels.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    heading.border.bottom.style = "solid",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    row_group.border.bottom.style = "none"
  ) %>%
  gt::tab_style(
    style = cell_text(size = "small"),
    locations = cells_body()
  ) %>%
  gt::tab_footnote(
    footnote = "Data are median (IQR) for continuous variables and n (%) for categorical variables",
    locations = cells_column_labels(columns = 1)
  )

# Print table
lancet_table


