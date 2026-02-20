# ============================================================
# UPMBC mixed models (LMM) and non-parametric tests (Kruskal/Dunn)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
  library(FSA)      # dunnTest()
})

# -----------------------------
# 0) Paths (repo-friendly)
# -----------------------------
DATA_FILE <- file.path("data", "UPMBC_LMM.csv")   # put your csv in /data
OUT_DIR   <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load + clean
# -----------------------------
frog <- read.csv(DATA_FILE, header = TRUE, stringsAsFactors = FALSE)

# Required columns check (fails early with clear message)
required_cols <- c("TL", "SVL", "Forest_type", "Habitat_type", "Species")
missing_cols <- setdiff(required_cols, colnames(frog))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input CSV: ", paste(missing_cols, collapse = ", "))
}

# Clean numeric columns (handles accidental text)
frog <- frog %>%
  mutate(
    TL  = as.numeric(gsub("[^0-9.\\-]", "", TL)),
    SVL = as.numeric(gsub("[^0-9.\\-]", "", SVL)),
    Forest_type  = as.factor(Forest_type),
    Habitat_type = as.factor(Habitat_type),
    Species      = as.factor(Species)
  )

# Keep a clean subset for analyses needing TL + SVL
frog_tl <- frog %>% filter(!is.na(TL), !is.na(SVL))

cat("\nRows total:", nrow(frog), "\n")
cat("Rows with TL+SVL:", nrow(frog_tl), "\n")
cat("\nForest types (all data):\n"); print(table(frog$Forest_type, useNA = "ifany"))
cat("\nForest types (TL+SVL only):\n"); print(table(frog_tl$Forest_type, useNA = "ifany"))

# -----------------------------
# 2) Size correction: TL ~ SVL
#   -> TL_resid as size-corrected TL (BCI proxy in your workflow)
# -----------------------------
size_model <- lm(TL ~ SVL, data = frog_tl)
sink(file.path(OUT_DIR, "TL_SVL_regression_summary.txt"))
print(summary(size_model))
sink()

frog_tl <- frog_tl %>%
  mutate(
    TL_predicted = predict(size_model),
    TL_resid     = resid(size_model)   # size-corrected TL
  )

write.csv(frog_tl, file.path(OUT_DIR, "UPMBC_TL_corrected_dataset.csv"), row.names = FALSE)

# Diagnostic plot (optional, but useful)
p_tl_svl <- ggplot(frog_tl, aes(SVL, TL)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal(base_size = 13) +
  labs(title = "Raw TL vs SVL (UPMBC)", x = "SVL", y = "TL")
ggsave(file.path(OUT_DIR, "Fig_TL_vs_SVL.png"), p_tl_svl, width = 6.5, height = 5, dpi = 300)

# -----------------------------
# LMM + emmeans + save outputs
# -----------------------------
run_lmm <- function(df, formula, emm_specs = NULL, name = "model") {
  m <- lmer(formula, data = df)

  # Save summary + ANOVA
  sink(file.path(OUT_DIR, paste0(name, "_summary.txt")))
  cat("Formula: ", deparse(formula), "\n\n")
  print(summary(m))
  cat("\nANOVA:\n")
  print(anova(m))
  sink()

  # EMMeans (optional)
  if (!is.null(emm_specs)) {
    emm <- emmeans(m, emm_specs)
    emm_df <- as.data.frame(emm)

    write.csv(emm_df, file.path(OUT_DIR, paste0(name, "_emmeans.csv")), row.names = FALSE)

    sink(file.path(OUT_DIR, paste0(name, "_emmeans.txt")))
    print(emm)
    cat("\nPairwise:\n")
    print(pairs(emm))
    sink()

    return(list(model = m, emmeans = emm, emmeans_df = emm_df))
  }

  return(list(model = m))
}

# -----------------------------
# Kruskal–Wallis + Dunn (Bonferroni) + save the outputs
# -----------------------------
run_kw_dunn <- function(df, response, group, name = "kw") {
  f <- as.formula(paste(response, "~", group))
  kw <- kruskal.test(f, data = df)

  sink(file.path(OUT_DIR, paste0(name, "_kruskal.txt")))
  cat("Formula: ", deparse(f), "\n\n")
  print(kw)
  sink()

  # Dunn only if KW is significant AND >=2 groups
  if (!is.na(kw$p.value) && kw$p.value < 0.05) {
    dunn <- dunnTest(f, data = df, method = "bonferroni")
    dunn_df <- dunn$res %>%
      mutate(Significance = case_when(
        P.adj < 0.001 ~ "***",
        P.adj < 0.01  ~ "**",
        P.adj < 0.05  ~ "*",
        TRUE ~ ""
      )) %>% arrange(P.adj)

    write.csv(dunn_df, file.path(OUT_DIR, paste0(name, "_dunn.csv")), row.names = FALSE)
    return(list(kw = kw, dunn = dunn_df))
  } else {
    return(list(kw = kw, dunn = NULL))
  }
}

# ============================================================
# 3) LMM: SVL ~ Forest_type + Habitat_type + (1|Species)
# (uses all SVL rows; TL not required)
# ============================================================

frog_svl <- frog %>% filter(!is.na(SVL), !is.na(Forest_type), !is.na(Habitat_type), !is.na(Species))

res_svl <- run_lmm(
  df       = frog_svl,
  formula  = SVL ~ Forest_type + Habitat_type + (1 | Species),
  emm_specs = ~ Forest_type,
  name     = "SVL_LMM"
)

# Plot EMMeans (Forest_type) with CIs
p_svl_emm <- ggplot(res_svl$emmeans_df, aes(Forest_type, emmean)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  theme_minimal(base_size = 13) +
  labs(title = "Estimated marginal means of SVL by forest type (LMM)",
       x = "Forest type", y = "Estimated SVL (mm)") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(OUT_DIR, "Fig_SVL_emmeans_forest.png"), p_svl_emm, width = 7, height = 5, dpi = 300)

# Non-parametric (optional): SVL ~ Forest_type / Habitat_type / Species
run_kw_dunn(frog_svl, "SVL", "Forest_type",  name = "SVL_by_forest")
run_kw_dunn(frog_svl, "SVL", "Habitat_type", name = "SVL_by_habitat")
run_kw_dunn(frog_svl, "SVL", "Species",      name = "SVL_by_species")

# ============================================================
# 4) LMM: TL_resid ~ Forest_type + Habitat_type + (1|Species)
# (requires TL+SVL)
# ============================================================

frog_tl2 <- frog_tl %>% filter(!is.na(TL_resid), !is.na(Forest_type), !is.na(Habitat_type), !is.na(Species))

res_tlres <- run_lmm(
  df        = frog_tl2,
  formula   = TL_resid ~ Forest_type + Habitat_type + (1 | Species),
  emm_specs = ~ Habitat_type,
  name      = "TLresid_LMM"
)

# Plot TL_resid across habitats (community visual)
p_tlres_hab <- ggplot(frog_tl2, aes(Habitat_type, TL_resid)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 13) +
  labs(title = "Size-corrected TL (residuals) across microhabitats",
       x = "Habitat type", y = "TL residual (TL_resid)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "Fig_TLresid_by_habitat.png"), p_tlres_hab, width = 8, height = 5, dpi = 300)

# ============================================================
# 5) Optional: Planted forest only (BCI-type plots/tests)
# ============================================================

if ("Planted_forest" %in% levels(frog_tl2$Forest_type)) {
  planted <- frog_tl2 %>% filter(Forest_type == "Planted_forest") %>% droplevels()

  # Plot by species with habitat as fill
  p_bci <- ggplot(planted, aes(Species, TL_resid, fill = Habitat_type)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 1) +
    theme_minimal(base_size = 13) +
    labs(title = "Body condition proxy (TL_resid) across species (Planted forest)",
         x = "Species", y = "BCI proxy (TL residual)", fill = "Habitat") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "Fig_BCI_planted_species.png"), p_bci, width = 9, height = 5.5, dpi = 300)

  # Kruskal/Dunn in planted forest
  run_kw_dunn(planted, "TL_resid", "Habitat_type", name = "BCIplanted_by_habitat")
  run_kw_dunn(planted, "TL_resid", "Species",      name = "BCIplanted_by_species")

  write.csv(planted, file.path(OUT_DIR, "UPMBC_planted_BCI_dataset.csv"), row.names = FALSE)
}

cat("\nDONE. Outputs written to:", OUT_DIR, "\n")
