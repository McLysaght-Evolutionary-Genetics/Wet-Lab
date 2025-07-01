# ──────────────────────────────────────────────────────────────────────────────
# 0) Install & load required packages
#    Checks for each package and installs missing ones, then loads them.
# ──────────────────────────────────────────────────────────────────────────────
packages <- c("lme4", "lmerTest", "emmeans", 
              "ggplot2", "dplyr", "tidyr", "tcltk")
# Identify packages not yet installed
to_install <- setdiff(packages, rownames(installed.packages()))
# Install any missing packages
if (length(to_install)) install.packages(to_install)
# Load all packages into the session
lapply(packages, library, character.only = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 1) Import & preprocess
#    Prompts user to select a CSV, reads it, and creates derived columns.
# ──────────────────────────────────────────────────────────────────────────────
file_path <- file.choose()  # Opens a dialog for file selection
# Validate that the chosen file is a CSV
if (tolower(tools::file_ext(file_path)) != "csv") 
  stop("The selected file is not a CSV file.")
data <- read.csv(file_path)  # Read data into a data.frame
head(data)  # Print the first few rows for inspection

# Add derived variables:
data <- data %>%
  mutate(
    Rep           = factor(rep),  # Convert replicate ID to factor
    Group         = factor(      # Create combined Strain×Condition factor
      paste(Strain, condition),
      levels = c(
        "N2 NoCS", "N2 10HCS", "N2 24HCS",
        "MCL2 NoCS", "MCL2 10HCS", "MCL2 24HCS"
      )
    ) %>% droplevels(),           # Drop unused factor levels
    combined_ie   = int.value + egg.value,        # Sum of int + egg fluorescence
    avg_egg_fluor = egg.value / number.of.eggs   # Per-egg average fluorescence
  )

# ──────────────────────────────────────────────────────────────────────────────
# 2) Per‐group t‐tests: worm.value vs combined_ie
#    Performs a paired t-test for each Group and assigns significance stars.
# ──────────────────────────────────────────────────────────────────────────────
# Function to map p-values to star labels
star_label <- function(p) {
  if      (is.na(p))      NA
  else if (p < 1e-4)      "****"
  else if (p < 1e-3)      "***"
  else if (p < 1e-2)      "**"
  else if (p < 0.05)      "*"
  else                    "ns"
}

# Compute p-values for each Group
ttest_df <- data %>%
  group_by(Group) %>%
  summarise(
    p.value = if (n() >= 2) t.test(worm.value, combined_ie)$p.value 
    else NA_real_,  # Only perform test if at least two observations
    .groups = "drop"
  ) %>%
  mutate(label = sapply(p.value, star_label))  # Add star labels

# ──────────────────────────────────────────────────────────────────────────────
# 3) Individual‐value LMM dot‐plots for worm/int/egg
#    For each metric, fit an LMM, extract emmeans and plot with replicate dots.
# ──────────────────────────────────────────────────────────────────────────────
# Define colors for each replicate
rep_colors  <- c("1"="red","2"="green","3"="blue",
                 "4"="orange","5"="yellow","6"="grey")
# Metrics to iterate over
value_types <- c("worm.value","int.value","egg.value")
# Pairs of groups for contrast filtering
comparisons <- list(
  c("N2 NoCS","N2 10HCS"), c("N2 NoCS","N2 24HCS"), c("N2 10HCS","N2 24HCS"),
  c("MCL2 NoCS","MCL2 10HCS"), c("MCL2 NoCS","MCL2 24HCS"), c("MCL2 10HCS","MCL2 24HCS"),
  c("N2 NoCS","MCL2 NoCS"),    c("N2 10HCS","MCL2 10HCS"),   c("N2 24HCS","MCL2 24HCS")
)

for (val_col in value_types) {
  message("Processing: ", val_col)
  
  # Subset and rename for modeling
  dat    <- data %>% select(Value = all_of(val_col), Group, Rep)
  # Fit the linear mixed model
  model  <- lmer(Value ~ Group + (1|Rep), data = dat)
  # Extract estimated marginal means per Group
  emm    <- emmeans(model, ~Group) %>% as.data.frame()
  
  # Compute all pairwise contrasts
  all_contrasts <- summary(pairs(emmeans(model, ~Group)))
  # Filter only the contrasts of interest
  contrast_df   <- subset(
    all_contrasts, 
    contrast %in% sapply(comparisons, paste, collapse = " - ")
  )
  # Add star labels for p-values
  contrast_df$stars <- sapply(contrast_df$p.value, star_label)
  cat(sprintf("\nPairwise comparisons for %s:\n", val_col))
  print(contrast_df)
  
  # Compute per-replicate means for overlay dots
  dot_df <- dat %>%
    group_by(Group, Rep) %>%
    summarise(DayMean = mean(Value), .groups = "drop")
  
  # Build the ggplot
  p <- ggplot() +
    geom_col(data = emm,
             aes(x = Group, y = emmean),
             # fill N2 groups white and MCL2 groups grey
             fill = ifelse(grepl("^MCL2", emm$Group), "grey80", "white"),
             color = "black", width = 0.7) +
    geom_errorbar(
      data = emm,
      aes(x = Group, ymin = emmean - SE, ymax = emmean + SE),
      width = 0.2
    ) +
    geom_point(
      data     = dot_df,
      aes(x = Group, y = DayMean, fill = Rep),
      position = position_jitter(width = 0.1),
      shape    = 21, size = 3.5, stroke = 0.2
    ) +
    scale_fill_manual(values = rep_colors) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste("LMM: Estimated Means for", val_col),
      x     = "Strain × Condition",
      y     = paste0(val_col, " (± SE)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid      = element_blank(),
      axis.line       = element_line(color = "black", size = 1.1),
      axis.ticks      = element_line(color = "black"),
      axis.text       = element_text(color = "black"),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  print(p)
}

# ──────────────────────────────────────────────────────────────────────────────
# 4) Combined bar: worm vs. stacked int+egg (with SEM & stars)
#    Summarizes per-replicate means, plots side-by-side worm and stacked int+egg
# ──────────────────────────────────────────────────────────────────────────────
rep_means_all <- data %>%
  group_by(Group, Rep) %>%
  summarise(
    worm_rep     = mean(worm.value, na.rm = TRUE),  # Mean worm fluorescence
    int_rep      = mean(int.value,  na.rm = TRUE),  # Mean int fluorescence
    egg_rep      = mean(egg.value,  na.rm = TRUE),  # Mean egg fluorescence
    combined_rep = int_rep + egg_rep,               # Combined int+egg
    .groups      = "drop"
  )

plot_summary <- rep_means_all %>%
  group_by(Group) %>%
  summarise(
    worm        = mean(worm_rep),
    int         = mean(int_rep),
    egg         = mean(egg_rep),
    combined    = mean(combined_rep),
    worm_se     = sd(worm_rep)     / sqrt(n()),
    combined_se = sd(combined_rep) / sqrt(n()),
    .groups     = "drop"
  ) %>%
  mutate(
    pos   = as.numeric(Group),
    y_pos = pmax(worm, combined) + 0.35 * max(pmax(worm, combined))
  ) %>%
  left_join(ttest_df, by = "Group")  # Attach t-test star labels

int_egg_df <- plot_summary %>%
  select(Group, pos, int, egg) %>%
  pivot_longer(c(int, egg), names_to = "Type", values_to = "Value")

type_colors <- c(int = "skyblue", egg = "tomato")

p2 <- ggplot() +
  geom_col(
    data = plot_summary,
    aes(x = pos - 0.15, y = worm),
    # fill N2 groups white and MCL2 groups grey
    fill = ifelse(grepl("^MCL2", emm$Group), "grey80", "white"),
    color = "black", width = 0.3) +
  geom_errorbar(
    data = plot_summary,
    aes(x = pos - 0.15, ymin = worm - worm_se, ymax = worm + worm_se),
    width = 0.1
  ) +
  geom_col(
    data     = int_egg_df,
    aes(x = pos + 0.15, y = Value, fill = Type),
    color    = "black",
    width    = 0.3,
    position = position_stack()
  ) +
  geom_errorbar(
    data = plot_summary,
    aes(x = pos + 0.15, ymin = combined - combined_se, ymax = combined + combined_se),
    width = 0.1
  ) +
  geom_text(
    data = plot_summary,
    aes(x = pos, y = y_pos, label = label),
    size = 5, vjust = 0
  ) +
  scale_x_continuous(
    breaks = plot_summary$pos,
    labels = plot_summary$Group,
    expand = c(0.01, 0)
  ) +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "Worm Value vs. Stacked Int and Egg Values",
    x     = "Strain × Condition",
    y     = "Mean Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black", size = 1.1),
    axis.ticks    = element_line(color = "black"),
    axis.text     = element_text(color = "black"),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_blank()
  )
print(p2)

# ──────────────────────────────────────────────────────────────────────────────
# 5) Third bar: avg fluorescence per egg (LMM)
#    Fit LMM on avg_egg_fluor and plot estimated means with replicate dots
# ──────────────────────────────────────────────────────────────────────────────
model_avg <- lmer(avg_egg_fluor ~ Group + (1|Rep), data = data)
emm_avg   <- emmeans(model_avg, ~Group) %>% as.data.frame()

dot_df_avg <- data %>%
  group_by(Group, Rep) %>%
  summarise(DayMean = mean(avg_egg_fluor), .groups = "drop")

p3 <- ggplot() +
  geom_col(
    data = emm_avg,
    aes(x = Group, y = emmean),
    # fill N2 groups white and MCL2 groups grey
    fill = ifelse(grepl("^MCL2", emm$Group), "grey80", "white"),
    color = "black", width = 0.7) +
  geom_errorbar(
    data = emm_avg,
    aes(x = Group, ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2
  ) +
  geom_point(
    data     = dot_df_avg,
    aes(x = Group, y = DayMean, fill = Rep),
    position = position_jitter(width = 0.1),
    shape    = 21, size = 3, stroke = 0.2
  ) +
  scale_fill_manual(values = rep_colors) +
  labs(
    title = "LMM: Avg Fluorescence per Egg",
    x     = "Strain × Condition",
    y     = "Avg Fluorescence per Egg (± SE)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black", size = 1.1),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
print(p3)

# ──────────────────────────────────────────────────────────────────────────────
# 6) Fourth bar: side‐by‐side egg.value vs avg fluorescence (dual‐axis, LMM)
#    Same treatment, new levels included via Group
# ──────────────────────────────────────────────────────────────────────────────
rep_means <- data %>%
  group_by(Group, Rep) %>%
  summarise(
    Egg    = mean(egg.value,     na.rm = TRUE),
    Fluor  = mean(avg_egg_fluor, na.rm = TRUE),
    .groups = "drop"
  )
model_e <- lmer(Egg   ~ Group + (1|Rep), data = rep_means)
model_f <- lmer(Fluor ~ Group + (1|Rep), data = rep_means)
emm_e <- emmeans(model_e, ~Group) %>% as.data.frame()
emm_f <- emmeans(model_f, ~Group) %>% as.data.frame()
scale_factor <- max(emm_e$emmean) / max(emm_f$emmean)

emm_e_df <- emm_e %>% transmute(Group, Metric = "Egg value", Mean = emmean, SE = SE)
emm_f_df <- emm_f %>% transmute(
  Group, Metric = "Avg fluorescence per egg",
  Mean = emmean * scale_factor, SE = SE * scale_factor
)
summary_long <- bind_rows(emm_e_df, emm_f_df) %>%
  mutate(Metric = factor(Metric, levels = c("Egg value","Avg fluorescence per egg")))
rep_long <- rep_means %>%
  pivot_longer(c(Egg, Fluor), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("Egg","Fluor"), labels = c("Egg value","Avg fluorescence per egg")),
    Value = if_else(Metric == "Egg value", Value, Value * scale_factor)
  )
bar_width <- 0.4
pd <- position_dodge(width = bar_width + 0.1)

# Prepare fill colors for Egg vs Avg-per-egg dual-axis
summary_long <- summary_long %>%
  mutate(
    fill_color = case_when(
      Metric == "Egg value" & grepl("^N2", Group)    ~ "white",
      Metric == "Egg value" & grepl("^MCL2", Group)  ~ "grey80",
      Metric == "Avg fluorescence per egg" & grepl("^N2", Group)    ~ "lightblue",
      Metric == "Avg fluorescence per egg" & grepl("^MCL2", Group)  ~ "deepskyblue"
    )
  )

# Plot dual-axis bar chart with replicate dots and custom fills
bar_width <- 0.4
pd        <- position_dodge(width = bar_width + 0.1)

p4 <- ggplot() +
  geom_col(
    data     = summary_long,
    aes(x = Group, y = Mean, fill = fill_color, group = Metric),
    color    = "black",
    width    = bar_width,
    position = pd
  ) +
  geom_errorbar(
    data     = summary_long,
    aes(x = Group, ymin = Mean - SE, ymax = Mean + SE, group = Metric),
    width    = 0.1,
    position = pd
  ) +
  geom_point(
    data     = rep_long,
    aes(x = Group, y = Value, color = factor(Rep), group = Metric),
    shape    = 19,
    size     = 2.5,
    position = pd
  ) +
  scale_fill_identity(guide = "none") +
  scale_color_manual(values = rep_colors) +
  scale_y_continuous(
    name     = "Egg value",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "Avg fluorescence per egg")
  ) +
  labs(
    title = "Egg Value vs. Average Fluorescence per Egg
(LMM‐estimated means)",
    x     = "Strain × Condition"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black", size = 1),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p4)



# ──────────────────────────────────────────────────────────────────────────────
# 7) LMM: Number of Eggs (Selective Contrasts Only)
#    LMM comparing timepoints within strains and matched conditions across strains
# ──────────────────────────────────────────────────────────────────────────────

# Fit LMM
model_eggs <- lmer(number.of.eggs ~ Group + (1|Rep), data = data)

# Get estimated marginal means
emm_eggs <- emmeans(model_eggs, ~ Group)

# Define desired comparisons
selected_contrasts <- list(
  c("N2 NoCS", "N2 10HCS"), c("N2 NoCS", "N2 24HCS"), c("N2 10HCS", "N2 24HCS"),
  c("MCL2 NoCS", "MCL2 10HCS"), c("MCL2 NoCS", "MCL2 24HCS"), c("MCL2 10HCS", "MCL2 24HCS"),
  c("N2 NoCS", "MCL2 NoCS"), c("N2 10HCS", "MCL2 10HCS"), c("N2 24HCS", "MCL2 24HCS")
)
contrast_labels <- sapply(selected_contrasts, \(x) paste(x, collapse = " - "))

# Get all pairwise contrasts
all_contrasts <- summary(pairs(emm_eggs))

# Filter to only selected contrasts
filtered_contrasts <- subset(all_contrasts, contrast %in% contrast_labels)

# Add stars
filtered_contrasts$stars <- sapply(filtered_contrasts$p.value, star_label)

# Show result
cat("\nSelected pairwise comparisons for number of eggs:\n")
print(filtered_contrasts)

# Plot estimated means and replicates
emm_df <- as.data.frame(emm_eggs)
dot_df <- data %>%
  group_by(Group, Rep) %>%
  summarise(DayMean = mean(number.of.eggs), .groups = "drop")

p_eggs <- ggplot() +
  geom_col(
    data = emm_df,
    aes(x = Group, y = emmean),
    fill = ifelse(grepl("^MCL2", emm_df$Group), "grey80", "white"),
    color = "black", width = 0.7
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = Group, ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2
  ) +
  geom_point(
    data = dot_df,
    aes(x = Group, y = DayMean, fill = Rep),
    position = position_jitter(width = 0.1),
    shape = 21, size = 3, stroke = 0.2
  ) +
  scale_fill_manual(values = rep_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "LMM: Number of Eggs per Worm (Selective Contrasts)",
    x = "Strain × Condition",
    y = "Number of Eggs (± SE)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1.1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_eggs)
