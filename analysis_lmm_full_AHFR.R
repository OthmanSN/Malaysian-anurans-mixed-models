# =============================================================================
# AHFR Full Dataset Analysis
# Linear Mixed Models + emmeans
# Kruskal–Wallis + Dunn (per species)
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(MASS)
  library(FSA)
})

# =============================================================================
# PATHS
# =============================================================================

DATA_DIR <- "data"
OUT_DIR  <- "outputs"

dir.create(OUT_DIR, showWarnings = FALSE)

INPUT_FILE <- file.path(DATA_DIR, "AHFR_GLMM_full.csv")

if (!file.exists(INPUT_FILE)) {
  stop("Place AHFR_GLMM_full.csv inside data/ folder.")
}

df <- read.csv(INPUT_FILE)

# =============================================================================
# PREPARE DATA
# =============================================================================

df <- df %>%
  mutate(
    Species = factor(Species),
    Season  = factor(Season),
    Zone    = factor(Zone,
                     levels = c("Controlled",
                                "Disturbed",
                                "Moderately disturbed",
                                "Protected"))
  )

# =============================================================================
# DERIVED VARIABLES
# =============================================================================

# TL size correction (standard residual method)
tl_lm <- lm(TL ~ SVL, data = df)
df$TL_resid <- residuals(tl_lm)

# Body condition index
df$BCI <- df$Weight / (df$SVL^3)

# Box-Cox transform BCI
min_bci <- min(df$BCI, na.rm = TRUE)
shift   <- ifelse(min_bci <= 0, abs(min_bci) + 1e-6, 0)

df$BCI_shift <- df$BCI + shift

bcx <- boxcox(BCI_shift ~ 1,
              data = df,
              lambda = seq(-2, 2, 0.1),
              plotit = FALSE)

lambda <- bcx$x[which.max(bcx$y)]

if (abs(lambda) < 1e-8) {
  df$BCI_trans <- log(df$BCI_shift)
} else {
  df$BCI_trans <- (df$BCI_shift^lambda - 1) / lambda
}

df$BCI_final <- log(df$BCI_trans - min(df$BCI_trans) + 1)

write.csv(df,
          file.path(OUT_DIR, "transformed_dataset.csv"),
          row.names = FALSE)

# =============================================================================
# FUNCTION 1: RUN LMM + EMMEANS + SAVE PLOTS
# =============================================================================

run_lmm <- function(response, label) {

  cat("\n==============================\n")
  cat("Running LMM for:", label, "\n")
  cat("==============================\n")

  formula <- as.formula(
    paste(response, "~ Zone * Season + (1 | Species)")
  )

  model <- lmer(formula, data = df)

  capture.output(summary(model),
    file = file.path(OUT_DIR, paste0(label, "_LMM_summary.txt"))
  )

  # emmeans
  emm <- emmeans(model, pairwise ~ Zone * Season, adjust = "tukey")

  write.csv(as.data.frame(emm$emmeans),
            file.path(OUT_DIR,
                      paste0(label, "_emmeans.csv")),
            row.names = FALSE)

  write.csv(as.data.frame(emm$contrasts),
            file.path(OUT_DIR,
                      paste0(label, "_pairwise.csv")),
            row.names = FALSE)

  # Boxplot (Zone)
  p1 <- ggplot(df,
               aes(x = Zone,
                   y = .data[[response]],
                   fill = Zone)) +
        geom_boxplot() +
        theme_minimal() +
        labs(title = paste(label, "Across Zones"),
             y = label)

  ggsave(file.path(OUT_DIR,
                   paste0(label, "_Zone.png")),
         p1)

  # Caterpillar
  emm_df <- as.data.frame(emm$emmeans)

  p2 <- ggplot(emm_df,
               aes(x = emmean,
                   y = interaction(Zone, Season),
                   xmin = lower.CL,
                   xmax = upper.CL,
                   color = Season)) +
        geom_point(size = 2.8) +
        geom_errorbarh(height = 0.3) +
        theme_minimal() +
        labs(title = paste("EMMeans:", label),
             x = "Estimated mean")

  ggsave(file.path(OUT_DIR,
                   paste0(label, "_EMMeans.png")),
         p2)

  return(model)
}

# =============================================================================
# FUNCTION 2: KRUSKAL + DUNN PER SPECIES
# =============================================================================

run_kruskal <- function(response, label) {

  results <- data.frame(
    Species = character(),
    KW_p    = numeric(),
    stringsAsFactors = FALSE
  )

  for (sp in levels(df$Species)) {

    sub <- df[df$Species == sp, ]

    if (length(unique(sub$Zone)) > 1) {

      kw <- kruskal.test(
        as.formula(paste(response, "~ Zone")),
        data = sub
      )

      results <- rbind(results,
                       data.frame(
                         Species = sp,
                         KW_p = kw$p.value
                       ))
    } else {

      results <- rbind(results,
                       data.frame(
                         Species = sp,
                         KW_p = NA
                       ))
    }
  }

  results$Significant <-
    ifelse(!is.na(results$KW_p) &
           results$KW_p < 0.05,
           "Yes", "No")

  write.csv(results,
            file.path(OUT_DIR,
                      paste0(label,
                             "_Kruskal_by_species.csv")),
            row.names = FALSE)

  # Plot KW results
  p <- ggplot(results,
              aes(x = Species,
                  y = KW_p,
                  color = Significant)) +
       geom_point(size = 3) +
       geom_hline(yintercept = 0.05,
                  linetype = "dashed") +
       theme_minimal() +
       theme(axis.text.x =
             element_text(angle = 45,
                          hjust = 1)) +
       labs(title = paste("Kruskal–Wallis:", label),
            y = "p-value")

  ggsave(file.path(OUT_DIR,
                   paste0(label,
                          "_Kruskal_plot.png")),
         p)

  # Dunn posthoc
  for (sp in levels(df$Species)) {

    sub <- df[df$Species == sp, ]

    if (length(unique(sub$Zone)) > 1) {

      dunn <- FSA::dunnTest(
        as.formula(paste(response, "~ Zone")),
        data = sub,
        method = "bonferroni"
      )$res

      write.csv(dunn,
        file.path(OUT_DIR,
          paste0(label, "_Dunn_", sp, ".csv")),
        row.names = FALSE)
    }
  }
}

# =============================================================================
# RUN ALL ANALYSES
# =============================================================================

model_BCI <- run_lmm("BCI_final", "BCI")
run_kruskal("BCI_final", "BCI")

model_SVL <- run_lmm("SVL", "SVL")
run_kruskal("SVL", "SVL")

model_TL  <- run_lmm("TL_resid", "TL_resid")
run_kruskal("TL_resid", "TL_resid")

cat("\nAll analyses completed successfully.\n")
