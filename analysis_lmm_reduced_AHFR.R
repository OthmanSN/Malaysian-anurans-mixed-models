# =============================================================================
# AHFR reduced dataset: Linear mixed models + nonparametric tests
# Outcomes: BoxCox-adjusted BCI, TL_corrected, SVL, BoxCox_Weight
# Random effect: Species
# Fixed effects: Zone * Season
#
# Author: Siti N. Othman
# Last updated: 2026-02-20
# =============================================================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(FSA)        # dunnTest
  library(dunn.test)  # dunn.test
})

# ---- User config ----
DATA_DIR   <- "data"     # put your CSVs in ./data
OUT_DIR    <- "outputs"  # results + plots
PLOT_DIR   <- file.path(OUT_DIR, "plots")
TABLE_DIR  <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

ZONE_LEVELS <- c("Controlled", "Disturbed", "Moderately disturbed", "Protected")

# ---- Helpers ----
assert_has_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

prep_ahfr <- function(df) {
  assert_has_cols(df, c("Species", "Zone", "Season"))
  df %>%
    mutate(
      Species = factor(Species),
      Zone    = factor(Zone, levels = ZONE_LEVELS),
      Season  = factor(Season)
    )
}

fit_lmm <- function(df, response) {
  fml <- as.formula(paste0(response, " ~ Zone * Season + (1 | Species)"))
  mod <- lmer(fml, data = df, REML = TRUE)
  list(
    model = mod,
    summary = summary(mod),
    singular = isSingular(mod)
  )
}

emm_zone_season <- function(model, adjust = "tukey") {
  emmeans(model, pairwise ~ Zone * Season, adjust = adjust)
}

save_emm_tables <- function(emm_obj, prefix) {
  emm_tbl  <- as.data.frame(emm_obj$emmeans)
  pair_tbl <- as.data.frame(emm_obj$contrasts)

  write.csv(emm_tbl,
            file.path(TABLE_DIR, paste0(prefix, "_emmeans.csv")),
            row.names = FALSE)
  write.csv(pair_tbl,
            file.path(TABLE_DIR, paste0(prefix, "_pairs.csv")),
            row.names = FALSE)

  invisible(list(emmeans = emm_tbl, pairs = pair_tbl))
}

plot_box_by_zone <- function(df, y, title, ylab, filename, fill_var = "Zone") {
  p <- ggplot(df, aes(x = Zone, y = .data[[y]], fill = .data[[fill_var]])) +
    geom_boxplot(outlier.colour = "red", linewidth = 0.6) +
    labs(title = title, x = "Zone", y = ylab) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(PLOT_DIR, filename), p, width = 7.5, height = 4.8, dpi = 300)
  p
}

plot_box_by_zone_season <- function(df, y, title, ylab, filename) {
  p <- ggplot(df, aes(x = Zone, y = .data[[y]], fill = Season)) +
    geom_boxplot(outlier.colour = "red", linewidth = 0.6) +
    labs(title = title, x = "Zone", y = ylab) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(PLOT_DIR, filename), p, width = 7.5, height = 4.8, dpi = 300)
  p
}

plot_caterpillar_emm <- function(emm_df, title, xlab, filename) {
  # expects columns: emmean, lower.CL, upper.CL, Zone, Season
  emm_df <- emm_df %>%
    mutate(Group = interaction(Zone, Season, sep = " · ", drop = TRUE))

  p <- ggplot(emm_df,
              aes(x = emmean, y = Group, xmin = lower.CL, xmax = upper.CL, color = Season)) +
    geom_point(size = 2.8) +
    geom_errorbarh(height = 0.25) +
    labs(title = title, x = xlab, y = "Zone · Season") +
    theme_minimal()

  ggsave(file.path(PLOT_DIR, filename), p, width = 8.0, height = 5.5, dpi = 300)
  p
}

kruskal_zone <- function(df, y) {
  kruskal.test(as.formula(paste0(y, " ~ Zone")), data = df)
}

kruskal_species <- function(df, y) {
  kruskal.test(as.formula(paste0(y, " ~ Species")), data = df)
}

kruskal_interaction <- function(df, y) {
  kruskal.test(as.formula(paste0(y, " ~ interaction(Species, Zone)")), data = df)
}

dunn_zone <- function(df, y, method = "bonferroni") {
  FSA::dunnTest(as.formula(paste0(y, " ~ Zone")), data = df, method = method)$res
}

dunn_by_species <- function(df, y, method = "bonferroni") {
  out <- lapply(split(df, df$Species), function(sub_df) {
    if (length(unique(sub_df$Zone)) < 2) return(NULL)
    FSA::dunnTest(as.formula(paste0(y, " ~ Zone")), data = sub_df, method = method)$res
  })
  out[!sapply(out, is.null)]
}

save_kruskal <- function(test_obj, filename) {
  df_out <- data.frame(
    statistic = unname(test_obj$statistic),
    df        = unname(test_obj$parameter),
    p_value   = unname(test_obj$p.value)
  )
  write.csv(df_out, file.path(TABLE_DIR, filename), row.names = FALSE)
  df_out
}

save_dunn <- function(dunn_df, filename) {
  write.csv(dunn_df, file.path(TABLE_DIR, filename), row.names = FALSE)
  dunn_df
}

# ---- Load data ----
df <- read.csv(file.path(DATA_DIR, "AHFR_GLMM_5SP.csv"))
df <- prep_ahfr(df)

# Optional: quick sanity check
message("Rows: ", nrow(df), " | Species: ", nlevels(df$Species))

# =============================================================================
# 1) BCI (already prepared in your CSV as log_adjusted_boxcox_BCI)
# =============================================================================
assert_has_cols(df, "log_adjusted_boxcox_BCI")

bci_fit <- fit_lmm(df, "log_adjusted_boxcox_BCI")
capture.output(bci_fit$summary,
               file = file.path(OUT_DIR, "model_BCI_summary.txt"))
writeLines(paste("Singular:", bci_fit$singular),
           con = file.path(OUT_DIR, "model_BCI_singularity.txt"))

bci_emm <- emm_zone_season(bci_fit$model, adjust = "tukey")
bci_tables <- save_emm_tables(bci_emm, prefix = "BCI")

plot_box_by_zone(df,
                 y = "log_adjusted_boxcox_BCI",
                 title = "BCI across zones",
                 ylab = "Box–Cox adjusted BCI (log scale)",
                 filename = "BCI_box_zone.png")

plot_box_by_zone_season(df,
                        y = "log_adjusted_boxcox_BCI",
                        title = "BCI across zones and seasons",
                        ylab = "Box–Cox adjusted BCI (log scale)",
                        filename = "BCI_box_zone_season.png")

plot_caterpillar_emm(bci_tables$emmeans,
                     title = "Estimated marginal means (BCI)",
                     xlab = "EMM (log-adjusted Box–Cox BCI)",
                     filename = "BCI_emm_caterpillar.png")

# Nonparametric checks
bci_kw_zone <- kruskal_zone(df, "log_adjusted_boxcox_BCI")
save_kruskal(bci_kw_zone, "BCI_kruskal_zone.csv")

bci_kw_species <- kruskal_species(df, "log_adjusted_boxcox_BCI")
save_kruskal(bci_kw_species, "BCI_kruskal_species.csv")

bci_kw_int <- kruskal_interaction(df, "log_adjusted_boxcox_BCI")
save_kruskal(bci_kw_int, "BCI_kruskal_species_zone_interaction.csv")

bci_dunn <- dunn_zone(df, "log_adjusted_boxcox_BCI")
save_dunn(bci_dunn, "BCI_dunn_zone_bonferroni.csv")

bci_dunn_sp <- dunn_by_species(df, "log_adjusted_boxcox_BCI")
# save each species to its own file
invisible(mapply(function(nm, dd) {
  write.csv(dd, file.path(TABLE_DIR, paste0("BCI_dunn_zone_", nm, ".csv")),
            row.names = FALSE)
}, nm = names(bci_dunn_sp), dd = bci_dunn_sp))

# =============================================================================
# 2) TL corrected (recommended: residual TL, NOT the contradictory lines)
# =============================================================================
assert_has_cols(df, c("TL", "SVL"))

# Residual-based "size-corrected TL" (standard approach)
tl_lm <- lm(TL ~ SVL, data = df)
df$TL_corrected <- residuals(tl_lm)  # <-- keep this, simple and correct

write.csv(df %>% select(Species, Zone, Season, TL, SVL, TL_corrected),
          file.path(TABLE_DIR, "TL_corrected_residuals.csv"),
          row.names = FALSE)

tl_fit <- fit_lmm(df, "TL_corrected")
capture.output(tl_fit$summary,
               file = file.path(OUT_DIR, "model_TLcorrected_summary.txt"))
writeLines(paste("Singular:", tl_fit$singular),
           con = file.path(OUT_DIR, "model_TLcorrected_singularity.txt"))

tl_emm <- emm_zone_season(tl_fit$model, adjust = "tukey")
tl_tables <- save_emm_tables(tl_emm, prefix = "TLcorrected")

plot_box_by_zone(df,
                 y = "TL_corrected",
                 title = "Size-corrected TL across zones",
                 ylab = "TL residuals (TL ~ SVL)",
                 filename = "TLcorrected_box_zone.png")

plot_box_by_zone_season(df,
                        y = "TL_corrected",
                        title = "Size-corrected TL across zones and seasons",
                        ylab = "TL residuals (TL ~ SVL)",
                        filename = "TLcorrected_box_zone_season.png")

plot_caterpillar_emm(tl_tables$emmeans,
                     title = "Estimated marginal means (TL corrected)",
                     xlab = "EMM (TL residuals)",
                     filename = "TLcorrected_emm_caterpillar.png")

# =============================================================================
# 3) SVL
# =============================================================================
svl_fit <- fit_lmm(df, "SVL")
capture.output(svl_fit$summary,
               file = file.path(OUT_DIR, "model_SVL_summary.txt"))
writeLines(paste("Singular:", svl_fit$singular),
           con = file.path(OUT_DIR, "model_SVL_singularity.txt"))

svl_emm <- emm_zone_season(svl_fit$model, adjust = "tukey")
svl_tables <- save_emm_tables(svl_emm, prefix = "SVL")

plot_box_by_zone(df,
                 y = "SVL",
                 title = "SVL across zones",
                 ylab = "SVL",
                 filename = "SVL_box_zone.png")

plot_box_by_zone_season(df,
                        y = "SVL",
                        title = "SVL across zones and seasons",
                        ylab = "SVL",
                        filename = "SVL_box_zone_season.png")

plot_caterpillar_emm(svl_tables$emmeans,
                     title = "Estimated marginal means (SVL)",
                     xlab = "EMM (SVL)",
                     filename = "SVL_emm_caterpillar.png")

# =============================================================================
# 4) Weight (BoxCox_Weight in separate file, as you had df2)
# =============================================================================
weight_path <- file.path(DATA_DIR, "updated_data_with_boxcox_weight3.csv")
if (file.exists(weight_path)) {
  df2 <- read.csv(weight_path) %>% prep_ahfr()
  assert_has_cols(df2, "BoxCox_Weight")

  # Kruskal + Dunn
  w_kw_zone <- kruskal_zone(df2, "BoxCox_Weight")
  save_kruskal(w_kw_zone, "Weight_kruskal_zone.csv")

  w_kw_species <- kruskal_species(df2, "BoxCox_Weight")
  save_kruskal(w_kw_species, "Weight_kruskal_species.csv")

  w_kw_int <- kruskal_interaction(df2, "BoxCox_Weight")
  save_kruskal(w_kw_int, "Weight_kruskal_species_zone_interaction.csv")

  w_dunn <- dunn_zone(df2, "BoxCox_Weight")
  save_dunn(w_dunn, "Weight_dunn_zone_bonferroni.csv")

  plot_box_by_zone(df2,
                   y = "BoxCox_Weight",
                   title = "Box–Cox weight across zones",
                   ylab = "Box–Cox weight",
                   filename = "Weight_box_zone.png")

  plot_box_by_zone_season(df2,
                          y = "BoxCox_Weight",
                          title = "Box–Cox weight across zones and seasons",
                          ylab = "Box–Cox weight",
                          filename = "Weight_box_zone_season.png")
} else {
  message("Weight file not found (skipping): ", weight_path)
}

message("Done. Outputs saved to: ", normalizePath(OUT_DIR))
