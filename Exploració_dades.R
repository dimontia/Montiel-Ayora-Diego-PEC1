load("se_object.rda")
se
# Creem una taula per veure differents features per timepoint
table_timepoints <- table(rowData(se)$timepoint)
table_timepoints
# Creem una funció per trobar els numero de valors buits o infinits que hi ha en el dataset
check_SE_missing <- function(se) {
  cat("Assay matrix:\n")
  cat("  NA: ", sum(is.na(assay(se))), "\n")
  cat("  NaN:", sum(is.nan(assay(se))), "\n")
  cat("  Inf:", sum(is.infinite(assay(se))), "\n")
  cat("  % of missing or inifinite values:", (sum(is.na(assay(se)))+sum(is.nan(assay(se)))+sum(is.infinite(assay(se))))/(689*39)*100, "\n")
  cat("\nRowData NA count:", sum(is.na(as.data.frame(rowData(se)))), "\n")
  cat("ColData NA count:", sum(is.na(as.data.frame(colData(se)))), "\n")
}
check_SE_missing(se) # Cridem la funció
# Creem un dataset que inclogui el assay i a més le nom sense timepoint de cada feature
df_assay_feature <- as.data.frame(assay(se)) %>%
  tibble::rownames_to_column("RowID") %>%
  mutate(Feature = rowData(se)$feature)

# Creem un altre dataset on hem inclós un valor true o false en funció de si la row te un missing value o no
df_assay_feature_NA <- df_assay_feature %>%
  mutate(has_missing = apply(select(., -RowID, -Feature), 1, function(x) any(!is.finite(x))))

# Ara ajuntem tots els features iguals, independentment del timepoint, si en cap cas tenien un missing value, es mante el TRUE 
df_assay_feature_NA <- df_assay_feature_NA %>%
  group_by(Feature) %>%
  summarise(all_rows_good = all(!has_missing))

# Mostrem el numero de features no tenen cap missing value
sum(df_assay_feature_NA$all_rows_good)
# Creem df gran amb el id dels pacients, els diferents features i els seus valors i timepoints 
df_subject_feature_value_timepoint <- assay(se) %>%
  as.data.frame() %>%
  rownames_to_column("feature_id") %>%
  pivot_longer(-feature_id, names_to = "Subject", values_to = "Value") %>%
  left_join(
    rowData(se) %>%
      as.data.frame() %>%
      select(timepoint, feature) %>%
      tibble::rownames_to_column("feature_id"),
    by = "feature_id"
  )
# Fem el ANOVA per veure com influeixen els pacients sobre la variació en els metabòlits
ANOVA_subject <- df_subject_feature_value_timepoint %>%
  group_split(feature) %>%
  map_dfr(function(df) {
    feat <- unique(df$feature)
    if (length(unique(df$Subject)) > 1) {
      p <- tryCatch({
        summary(aov(Value ~ Subject, data = df))[[1]][["Pr(>F)"]][1]
      }, error = function(e) NA)
    } else {
      p <- NA
    }
    tibble(feature = feat, p_value = p)
  })
# Fem el ANOVA per veure com influeixe el temps sobre la variació en els metabòlits
ANOVA_timepoints <- df_subject_feature_value_timepoint %>%
  group_split(feature) %>%
  map_dfr(function(df) {
    feat <- unique(df$feature)
    if (length(unique(df$timepoint)) > 1) {
      p <- tryCatch({
        summary(aov(Value ~ timepoint, data = df))[[1]][["Pr(>F)"]][1]
      }, error = function(e) NA)
    } else {
      p <- NA
    }
    tibble(feature = feat, p_value = p)
  })
# Filtrem i mostrem el resultat del número de features que tenen variació significativa en funció dels pacients
significant_subject <- ANOVA_subject %>%
  filter(!is.na(p_value) & p_value < 0.05) %>%
  arrange(p_value)
cat("Number of features with significant variation amongst subjects, independent to time:", length(significant_subject$feature), "\n")
# Filtrem i mostrem el resultat del número de features que tenen variació significativa en funció dels timepoints
significant_timepoint <- ANOVA_timepoints %>%
  filter(!is.na(p_value) & p_value < 0.05) %>%
  arrange(p_value)
cat("Number of features with significant variation amongst timepoints, independent to the subjects:", length(significant_timepoint$feature), "\n")
time_only_features <- setdiff(significant_timepoint$feature, significant_subject$feature)
time_only_features
# ens quedem  només amb les features (sense timpeoint associat) del llistat que em mencionat
keep_rows <- rowData(se)$feature %in% time_only_features

# creem un df amb aquestes features i els seus valors associats
significant_features <- assay(se)[keep_rows, ]
significant_features_scaled <- t(scale(t(significant_features))) # escalem els valors per igualar diferents magnituds

# Afegim els timepoints de cada feature
df_significant_features_tiempoint <- as.data.frame(significant_features_scaled) %>%
  tibble::rownames_to_column("RowID") %>%
  pivot_longer(-RowID, names_to = "Sample", values_to = "Expression") %>%
  mutate(
    Timepoint = rowData(se)$timepoint[match(RowID, rownames(se))],
    Feature = rowData(se)$feature[match(RowID, rownames(se))]
  )  # Add proper feature names

# Centrem la variació al voltant del valor incial de T0 igualant a 0
df_significant_features_tiempoint_centered <- df_significant_features_tiempoint %>%
  group_by(Feature, Sample) %>%
  mutate(Expression = Expression - Expression[Timepoint == "T0"])

# Fem el plot
timepoint_plot <- ggplot(df_significant_features_tiempoint_centered, aes(x = Timepoint, y = Expression, group = Feature, color = Feature)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  labs(x = "Timepoint", y = "Mean Expression", title = "Timepoint-driven feature variation") +
  theme_minimal()
# 1. Seleccionem les features significatives que varien entre pacients però no pel temps
subject_only_features <- setdiff(significant_subject$feature, significant_timepoint$feature)

# 2. Ens quedem amb les files del SE que compleixin això
keep_rows_subject <- rowData(se)$feature %in% subject_only_features
subject_features <- assay(se)[keep_rows_subject, ]

# 3. Afegim la info de timepoint i feature
df_subject_features <- as.data.frame(subject_features) %>%
  tibble::rownames_to_column("RowID") %>%
  pivot_longer(-RowID, names_to = "Subject", values_to = "Expression") %>%
  mutate(
    Timepoint = rowData(se)$timepoint[match(RowID, rownames(se))],
    Feature = rowData(se)$feature[match(RowID, rownames(se))]
  )

# 4. Calculem la mitjana d’expressió per feature i pacient (ignorem el timepoint)
df_subject_features_avg <- df_subject_features %>%
  group_by(Feature, Subject) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Top 6 més variables entre pacients
top_features <- df_subject_features_avg %>%
  group_by(Feature) %>%
  summarise(sd = sd(Mean_Expression, na.rm = TRUE)) %>%
  slice_max(sd, n = 20) %>%
  pull(Feature)

df_top <- df_subject_features_avg %>% filter(Feature %in% top_features)

subject_plot <- ggplot(df_top, aes(x = Subject, y = Mean_Expression)) +
  geom_boxplot() +
  facet_wrap(~Feature, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))