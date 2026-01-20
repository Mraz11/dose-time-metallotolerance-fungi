#!/usr/bin/env Rscript
# ============================================================
# Publication-ready figure pipeline (R) — GitHub-friendly
# Figures: 1–8 + tidy tables (CSV)
# Inputs : Excel sheet with columns
#          (fungi/heavy_metal/time/concentration/replication/diameter)
# Outputs: PDF + 600 dpi LZW-TIFF for ggplots; PNG/PDF for heatmap;
#          IC50 table + Fig 8 (grouped bars)
# ============================================================

# ----------------------------
# 0) Reproducibility & config
# ----------------------------
options(stringsAsFactors = FALSE)
set.seed(1)

# GitHub-friendly paths:
# - Put Excel in: data/Final_New_Combined.xlsx
# - Outputs to : outputs/figures
infile  <- file.path("data", "Final_New_Combined.xlsx")
out_dir <- file.path("outputs", "figures")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(file.exists(infile))

# ----------------------------
# 1) Packages
# ----------------------------
pkgs <- c(
  "readxl","dplyr","tidyr","stringr","janitor","ggplot2","forcats",
  "purrr","scales","pheatmap","drc","tibble",
  "emmeans","multcomp"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(scales)
  library(pheatmap)
  library(drc)
  library(tibble)
  library(emmeans)
  library(multcomp)
})

# ----------------------------
# 2) Helpers (style + save)
# ----------------------------
theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "right"
    )
}

save_pub <- function(p, name, w = 190, h = 140, units = "mm", dpi = 800) {
  # TIFF (LZW) + PDF
  ggsave(
    filename = file.path(out_dir, paste0(name, ".tiff")),
    plot = p, width = w, height = h, units = units,
    dpi = dpi, compression = "lzw"
  )
  ggsave(
    filename = file.path(out_dir, paste0(name, ".pdf")),
    plot = p, width = w, height = h, units = units
  )
}

# Tukey compact letters (stable: emmeans + multcomp::cld)
tukey_letters <- function(dat) {
  dat <- dat %>% filter(!is.na(diameter), !is.na(concentration_f))
  if (n_distinct(dat$concentration_f) < 2) {
    return(tibble(concentration_f = unique(dat$concentration_f), letters = "a"))
  }
  fit <- aov(diameter ~ concentration_f, data = dat)
  em  <- emmeans::emmeans(fit, ~ concentration_f)
  
  multcomp::cld(em, Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    transmute(
      concentration_f = concentration_f,
      letters = stringr::str_trim(.group)
    )
}

# ----------------------------
# 3) Read + clean data
# ----------------------------
df <- readxl::read_excel(infile, sheet = 1) %>%
  clean_names() %>%
  mutate(
    fungi       = na_if(str_squish(as.character(fungi)), ""),
    heavy_metal = na_if(str_squish(as.character(heavy_metal)), "")
  ) %>%
  tidyr::fill(fungi, heavy_metal, .direction = "down") %>%
  rename(
    fungus = fungi,
    metal  = heavy_metal,
    time   = time,
    concentration = concentration,
    replicate = replication,
    diameter  = diameter
  ) %>%
  mutate(
    fungus = str_squish(as.character(fungus)),
    metal  = str_squish(as.character(metal)),
    
    fungus = case_when(
      str_detect(str_to_lower(fungus), "aspergillus|a\\.?\\s*niger") ~ "Aspergillus niger",
      str_detect(str_to_lower(fungus), "rhizopus|r\\.?\\s*stolonifer") ~ "Rhizopus stolonifer",
      str_detect(str_to_lower(fungus), "trichoderma|t\\.?\\s*asperellum") ~ "Trichoderma asperellum",
      TRUE ~ fungus
    ),
    metal = case_when(
      str_detect(str_to_lower(metal), "ars|\\bas\\b") ~ "As",
      str_detect(str_to_lower(metal), "lead|\\bpb\\b") ~ "Pb",
      str_detect(str_to_lower(metal), "merc|\\bhg\\b") ~ "Hg",
      TRUE ~ metal
    ),
    
    time = as.numeric(time),
    concentration = as.numeric(concentration),
    diameter = as.numeric(diameter),
    replicate = as.character(replicate)
  ) %>%
  filter(
    !is.na(fungus), !is.na(metal), !is.na(time),
    !is.na(concentration), !is.na(diameter)
  ) %>%
  mutate(
    plate_id = paste(fungus, metal, concentration, replicate, sep = "_")
  )

conc_levels <- sort(unique(df$concentration))
message("Concentrations found (mg L-1): ", paste(conc_levels, collapse = ", "))

fungi_order <- c("Aspergillus niger","Rhizopus stolonifer","Trichoderma asperellum")
metal_order <- c("As","Pb","Hg")

df <- df %>%
  mutate(
    fungus = fct_relevel(factor(fungus), fungi_order),
    metal  = fct_relevel(factor(metal), metal_order),
    concentration_f = factor(concentration, levels = sort(unique(concentration)))
  )

maxdose <- max(conc_levels, na.rm = TRUE)

# ----------------------------
# 4) Summaries
# ----------------------------
summ <- df %>%
  group_by(fungus, metal, concentration, concentration_f, time) %>%
  summarise(
    mean_diam = mean(diameter, na.rm = TRUE),
    sd_diam   = sd(diameter, na.rm = TRUE),
    n = sum(!is.na(diameter)),
    .groups = "drop"
  )

# ---- Tolerance Index (TI) replicate-first (for SD/errorbars) ----
ctrl_means <- df %>%
  filter(concentration == 0) %>%
  group_by(fungus, metal, time) %>%
  summarise(ctrl_mean = mean(diameter, na.rm = TRUE), .groups = "drop")

ti_raw <- df %>%
  left_join(ctrl_means, by = c("fungus","metal","time")) %>%
  mutate(
    ti_val = ifelse(ctrl_mean > 0, diameter / ctrl_mean, NA_real_),
    inhibition = 100 * (1 - ti_val)
  )

ti <- ti_raw %>%
  group_by(fungus, metal, concentration, concentration_f, time) %>%
  summarise(
    ti    = mean(ti_val, na.rm = TRUE),
    sd_ti = sd(ti_val, na.rm = TRUE),
    n = sum(!is.na(ti_val)),
    .groups = "drop"
  )

# ---- Growth rate slope per plate (mm h-1) ----
gr <- df %>%
  group_by(plate_id, fungus, metal, concentration, concentration_f) %>%
  filter(sum(!is.na(diameter)) >= 2) %>%
  group_modify(~{
    m <- lm(diameter ~ time, data = .x)
    tibble(gr = unname(coef(m)[["time"]]))
  }) %>%
  ungroup()

gr_summ <- gr %>%
  group_by(fungus, metal, concentration, concentration_f) %>%
  summarise(
    mean_gr = mean(gr, na.rm = TRUE),
    sd_gr   = sd(gr, na.rm = TRUE),
    n = sum(!is.na(gr)),
    .groups = "drop"
  )

# ============================================================
# FIG 1: Time-course growth (mean ± SD)
# ============================================================
p1 <- summ %>%
  ggplot(aes(x = time, y = mean_diam, group = concentration_f, colour = concentration_f)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.7) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_diam - sd_diam), ymax = mean_diam + sd_diam),
    width = 1.5, linewidth = 0.35
  ) +
  facet_grid(fungus ~ metal) +
  scale_x_continuous(breaks = sort(unique(summ$time))) +
  labs(
    x = "Time (h)",
    y = "Colony diameter (mm)",
    colour = "Concentration (mg L\u207B\u00B9)",
    title = "Colony expansion under metal stress over time"
  ) +
  theme_pub()

save_pub(p1, "Fig1_timecourse_growth", w = 200, h = 145, dpi = 800)

# ============================================================
# FIG 2: Endpoint diameter at 48 h + Tukey letters
# ============================================================
endpoint_h <- 48

end_df <- df %>%
  filter(time == endpoint_h) %>%
  mutate(concentration_f = factor(concentration, levels = sort(unique(concentration))))

letters_df <- end_df %>%
  group_by(fungus, metal) %>%
  group_modify(~ tukey_letters(.x)) %>%
  ungroup()

end_summ <- end_df %>%
  group_by(fungus, metal, concentration_f) %>%
  summarise(
    mean_diam = mean(diameter, na.rm = TRUE),
    sd_diam   = sd(diameter, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(letters_df, by = c("fungus","metal","concentration_f")) %>%
  mutate(label_y = mean_diam + ifelse(is.na(sd_diam), 0, sd_diam) + 2)

p2 <- ggplot(end_summ, aes(x = concentration_f, y = mean_diam)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_diam - sd_diam), ymax = mean_diam + sd_diam),
    width = 0.15, linewidth = 0.35
  ) +
  geom_text(aes(y = label_y, label = letters), size = 3.0, vjust = 0) +
  facet_grid(fungus ~ metal) +
  labs(
    x = "Concentration (mg L\u207B\u00B9)",
    y = paste0("Colony diameter at ", endpoint_h, " h (mm)"),
    title = paste0("Endpoint growth at ", endpoint_h, " h (Tukey HSD within fungus × metal)")
  ) +
  theme_pub()

save_pub(p2, "Fig2_endpoint_48h_letters", w = 200, h = 145, dpi = 800)

# ============================================================
# FIG 3: Growth-rate slope (mm h-1)
# ============================================================
p3 <- gr_summ %>%
  ggplot(aes(x = concentration_f, y = mean_gr)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_gr - sd_gr), ymax = mean_gr + sd_gr),
    width = 0.15, linewidth = 0.35
  ) +
  facet_grid(fungus ~ metal) +
  labs(
    x = "Concentration (mg L\u207B\u00B9)",
    y = "Colony expansion rate (mm h\u207B\u00B9)",
    title = "Colony expansion rate estimated as slope of diameter vs time"
  ) +
  theme_pub()

save_pub(p3, "Fig3_growth_rate", w = 200, h = 145, dpi = 800)

# ============================================================
# FIG 4: Tolerance index (TI) at 48 h
# ============================================================
ti_h <- 48
ti_end <- ti %>% filter(time == ti_h)

p4 <- ti_end %>%
  ggplot(aes(x = concentration_f, y = ti)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pmax(0, ti - sd_ti), ymax = ti + sd_ti),
    width = 0.15, linewidth = 0.35
  ) +
  facet_grid(fungus ~ metal) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) +
  labs(
    x = "Concentration (mg L\u207B\u00B9)",
    y = paste0("Tolerance index at ", ti_h, " h"),
    title = paste0("Tolerance index (TI) relative to control at ", ti_h, " h")
  ) +
  theme_pub()

save_pub(p4, "Fig4_tolerance_index_48h", w = 200, h = 145, dpi = 800)

# ============================================================
# FIG 5: Heatmap of TI at 48 h (PNG + PDF)
# ============================================================
heat_dat <- ti_end %>%
  mutate(col = paste0(as.character(metal), "_", concentration)) %>%
  dplyr::select(fungus, col, ti) %>%
  group_by(fungus, col) %>%
  summarise(ti = mean(ti, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = col, values_from = ti)

heat_mat <- heat_dat %>%
  tibble::column_to_rownames("fungus") %>%
  as.matrix()

png(file.path(out_dir, "Fig5_heatmap_TI_48h.png"), width = 2200, height = 900, res = 300)
pheatmap(
  heat_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  fontsize = 11,
  main = "Heatmap of tolerance index (48 h)",
  angle_col = 45
)
dev.off()

pdf(file.path(out_dir, "Fig5_heatmap_TI_48h.pdf"), width = 12, height = 5)
pheatmap(
  heat_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  fontsize = 11,
  main = "Heatmap of tolerance index (48 h)",
  angle_col = 45
)
dev.off()

# ============================================================
# FIG 6: Mean tolerance (%) at 48 h
# ============================================================
mean_tol <- ti_end %>%
  filter(concentration > 0) %>%
  group_by(fungus, metal) %>%
  summarise(
    mean_ti = mean(ti, na.rm = TRUE),
    sd_ti   = sd(ti, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_tol_pct = 100 * mean_ti,
    sd_tol_pct   = 100 * sd_ti
  )

p6 <- mean_tol %>%
  ggplot(aes(x = metal, y = mean_tol_pct, fill = fungus)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75, colour = "black", linewidth = 0.25
  ) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_tol_pct - sd_tol_pct), ymax = mean_tol_pct + sd_tol_pct),
    position = position_dodge(width = 0.8),
    width = 0.15, linewidth = 0.35
  ) +
  labs(
    x = "Metal",
    y = "Mean tolerance (%) at 48 h",
    fill = "Fungus",
    title = "Mean tolerance percentage (average TI × 100) at 48 h"
  ) +
  theme_pub()

save_pub(p6, "Fig6_mean_tolerance_percent", w = 170, h = 120, dpi = 800)

# ============================================================
# FIG 7: Dose–response curves (LL.4) for TI at 48 h
# - Includes 0 using pseudo-log transform
# ============================================================
dose_means <- ti_end %>%
  transmute(fungus, metal, concentration, mean_ti = ti, sd_ti = sd_ti)

pos_conc <- sort(unique(dose_means$concentration[dose_means$concentration > 0]))
min_pos  <- if (length(pos_conc) > 0) min(pos_conc) else 1

safe_fit <- function(dat) {
  dat <- dat %>% filter(!is.na(mean_ti))
  if (n_distinct(dat$concentration) < 4) return(NULL)
  
  dat <- dat %>% mutate(dose = ifelse(concentration == 0, 1e-6, concentration))
  
  m1 <- try(
    drm(
      mean_ti ~ dose, data = dat,
      fct = LL.4(fixed = c(NA, 0, 1, NA), names = c("b","c","d","e")),
      control = drmc(errorm = FALSE, noMessage = TRUE, maxIt = 200)
    ),
    silent = TRUE
  )
  if (!inherits(m1, "try-error")) return(m1)
  
  m2 <- try(
    drm(
      mean_ti ~ dose, data = dat,
      fct = LL.4(names = c("b","c","d","e")),
      control = drmc(errorm = FALSE, noMessage = TRUE, maxIt = 200)
    ),
    silent = TRUE
  )
  if (inherits(m2, "try-error")) return(NULL)
  m2
}

fit_tbl <- dose_means %>%
  group_by(fungus, metal) %>%
  group_split() %>%
  map(~{
    m <- safe_fit(.x)
    if (is.null(m)) return(NULL)
    tibble(fungus = as.character(.x$fungus[1]), metal = as.character(.x$metal[1]), model = list(m))
  }) %>%
  compact() %>%
  bind_rows()

pred_df <- fit_tbl %>%
  mutate(pred = pmap(list(model, fungus, metal), function(model, fungus, metal) {
    xseq <- exp(seq(log(min_pos), log(maxdose), length.out = 200))
    yhat <- predict(model, newdata = data.frame(dose = xseq))
    yhat <- pmin(pmax(yhat, 0), 1.05)
    tibble(fungus = fungus, metal = metal, concentration = xseq, ti_pred = yhat)
  })) %>%
  dplyr::select(pred) %>%
  unnest(pred)

p7 <- ggplot() +
  geom_errorbar(
    data = dose_means,
    aes(
      x = concentration, y = mean_ti,
      ymin = pmax(0, mean_ti - sd_ti),
      ymax = pmin(1.05, mean_ti + sd_ti)
    ),
    width = 0, linewidth = 0.3
  ) +
  geom_point(data = dose_means, aes(x = concentration, y = mean_ti), size = 1.8) +
  geom_line(data = pred_df, aes(x = concentration, y = ti_pred), linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  facet_grid(fungus ~ metal) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 50),
    breaks = conc_levels,
    labels = conc_levels
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Concentration (mg L\u207B\u00B9)",
    y = "Tolerance index (48 h)",
    title = "Dose–response curves (LL.4) with tolerance index at 48 h"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_pub(p7, "Fig7_dose_response_TI48h_symlog", w = 200, h = 145, dpi = 800)

# ============================================================
# IC50 (TI50) at 48 h + FIG 8 (GROUPED like your example)
# ============================================================
df48 <- df %>% filter(time == 48)

ctrl48 <- df48 %>%
  filter(concentration == 0) %>%
  dplyr::select(fungus, metal, replicate, ctrl_diam = diameter)

rep_ti <- df48 %>%
  left_join(ctrl48, by = c("fungus","metal","replicate")) %>%
  mutate(TI = diameter / ctrl_diam) %>%
  filter(!is.na(TI), is.finite(TI)) %>%
  mutate(TI = pmin(pmax(TI, 0), 1.2)) # mild clamp for numerical stability

ic50_interp <- function(dat) {
  dd <- dat %>%
    group_by(concentration) %>%
    summarise(ti = median(TI, na.rm = TRUE), .groups = "drop") %>%
    arrange(concentration)
  
  if (nrow(dd) < 2) return(NA_real_)
  
  idx <- which(dd$concentration > 0 & dd$ti <= 0.5)
  if (length(idx) == 0) return(Inf)
  
  j <- idx[1]
  if (j == 1) return(dd$concentration[j])
  
  x1 <- dd$concentration[j - 1]; y1 <- dd$ti[j - 1]
  x2 <- dd$concentration[j];     y2 <- dd$ti[j]
  
  if (!is.finite(y2 - y1) || (y2 - y1) == 0) return(dd$concentration[j])
  x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
}

fit_ic50_once <- function(dat) {
  dat2 <- dat %>%
    mutate(dose = ifelse(concentration == 0, 1e-6, concentration)) %>%
    filter(is.finite(TI), is.finite(dose))
  
  if (n_distinct(dat2$dose) < 4 || n_distinct(dat2$TI) < 3) return(ic50_interp(dat2))
  
  m <- suppressWarnings(try(
    drm(
      TI ~ dose, data = dat2,
      fct = LL.4(fixed = c(NA, 0, 1, NA), names = c("b","c","d","e")),
      control = drmc(errorm = FALSE, noMessage = TRUE, maxIt = 200)
    ),
    silent = TRUE
  ))
  
  if (inherits(m, "try-error")) {
    m <- suppressWarnings(try(
      drm(
        TI ~ dose, data = dat2,
        fct = LL.4(names = c("b","c","d","e")),
        control = drmc(errorm = FALSE, noMessage = TRUE, maxIt = 200)
      ),
      silent = TRUE
    ))
  }
  
  if (inherits(m, "try-error")) return(ic50_interp(dat2))
  
  ed <- suppressWarnings(try(ED(m, 50, interval = "none", display = FALSE), silent = TRUE))
  if (inherits(ed, "try-error")) return(ic50_interp(dat2))
  
  val <- as.numeric(ed[1, 1])
  if (!is.finite(val) || val <= 0) return(ic50_interp(dat2))
  val
}

bootstrap_ic50 <- function(dat, B = 200, seed = 1) {
  set.seed(seed)
  vals <- numeric(B)
  
  for (i in seq_len(B)) {
    samp <- dat %>% slice_sample(prop = 1, replace = TRUE)
    vals[i] <- fit_ic50_once(samp)
  }
  
  n_inf <- sum(is.infinite(vals), na.rm = TRUE)
  clean_vals <- vals[!is.na(vals) & is.finite(vals)]
  
  if (length(clean_vals) < (B * 0.15) && n_inf > (B * 0.5)) {
    return(c(med = Inf, lo = Inf, hi = Inf))
  }
  if (length(clean_vals) < (B * 0.15)) {
    return(c(med = NA, lo = NA, hi = NA))
  }
  
  c(
    med = median(clean_vals),
    lo  = unname(quantile(clean_vals, 0.025, na.rm = TRUE)),
    hi  = unname(quantile(clean_vals, 0.975, na.rm = TRUE))
  )
}

ic50_rows <- rep_ti %>%
  group_by(fungus, metal) %>%
  group_modify(~{
    est  <- fit_ic50_once(.x)
    boot <- bootstrap_ic50(.x, B = 200, seed = 1)
    tibble(IC50 = est, CI_low = boot["lo"], CI_high = boot["hi"], n = nrow(.x))
  }) %>%
  ungroup() %>%
  mutate(
    censored = is.infinite(IC50) | IC50 > maxdose,
    IC50_mg_L = case_when(
      is.na(IC50) ~ NA_character_,
      censored ~ paste0(">", maxdose),
      TRUE ~ sprintf("%.1f", IC50)
    )
  )

write.csv(ic50_rows, file.path(out_dir, "IC50_48h_TI_bootstrap.csv"), row.names = FALSE)

# ----------------------------
# FIG 8: grouped bars + metal colours (match your example)
# ----------------------------
ic50_plot <- ic50_rows %>%
  mutate(
    IC50_num     = suppressWarnings(as.numeric(IC50)),
    CI_low_num   = suppressWarnings(as.numeric(CI_low)),
    CI_high_num  = suppressWarnings(as.numeric(CI_high)),
    IC50_plot    = ifelse(censored, maxdose, IC50_num),
    CI_low_plot  = ifelse(is.finite(CI_low_num),  pmin(CI_low_num,  maxdose), NA_real_),
    CI_high_plot = ifelse(is.finite(CI_high_num), pmin(CI_high_num, maxdose), NA_real_)
  ) %>%
  filter(!is.na(IC50_plot), is.finite(IC50_plot)) %>%
  mutate(
    fungus = factor(fungus, levels = fungi_order),
    metal  = factor(metal,  levels = metal_order)
  )

p8 <- ggplot(ic50_plot, aes(x = fungus, y = IC50_plot, fill = metal)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(ymin = CI_low_plot, ymax = CI_high_plot),
    position = position_dodge(width = 0.8),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_text(
    data = dplyr::filter(ic50_plot, censored),
    aes(label = paste0(">", maxdose), y = IC50_plot),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3
  ) +
  scale_fill_manual(values = c(As = "#1f77b4", Pb = "#ff7f0e", Hg = "#2ca02c")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Dose causing 50% growth reduction (IC50) across fungi and metals",
    x = NULL,
    y = expression(IC[50]~"(mg L"^{-1}*") at 48 h (TI = 0.5)"),
    fill = "Metal"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_pub(p8, "Fig8_IC50_TI05_summary", w = 200, h = 130, dpi = 800)

# ----------------------------
# 5) Export tidy datasets
# ----------------------------
write.csv(summ,    file.path(out_dir, "summary_growth_by_group.csv"), row.names = FALSE)
write.csv(ti,      file.path(out_dir, "tolerance_index_by_group.csv"), row.names = FALSE)
write.csv(gr_summ, file.path(out_dir, "growth_rate_summary.csv"), row.names = FALSE)

# Save session info for GitHub reproducibility
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("All figures + tables saved to: ", normalizePath(out_dir))
