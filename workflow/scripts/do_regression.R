# OLS vs PGLS (Brownian, Grafen) using reg_vars.csv + Kuderna tree

setwd("~/PycharmProjects/PrimateDFE/")
set.seed(42)

library(tidyverse)
library(ape)
library(nlme)
library(car)

theme_set(theme_minimal(base_size = 14))

# -----------------------
# Inputs
# -----------------------
reg_file  <- "scratch/reg_vars.discrete.full.csv"
tree_file <- "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"

# -----------------------
# Data
# -----------------------
df <- read_csv(reg_file, show_col_types = FALSE)

names(df)  # check the column name is actually "population"

df <- df %>%
  rename(species = .data$population) %>%
  mutate(
    species = ifelse(
      str_count(species, "_") >= 2,
      sub("(_[^_]+)$", "", species),
      species
    ), # remove subspecies postfix if present
    logNe = log10(Ne)
  ) %>%
  group_by(species) %>%
  summarise(
    Ne = mean(Ne, na.rm = TRUE), # average Ne across subspecies
    range_inf_m10 = mean(.data$`range_inf_-10`, na.rm = TRUE),
    logNe = mean(logNe, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------
# Phylogeny
# -----------------------
primphylo <- read.tree(tree_file)

# keep only overlapping taxa and align rows to tree tip order
common <- intersect(primphylo$tip.label, df$species)

# drop tips not in common and filter/arrange df to match tree tip order
primphylo <- drop.tip(primphylo, setdiff(primphylo$tip.label, common))
df <- df %>%
  filter(species %in% primphylo$tip.label) %>%
  mutate(species = factor(species, levels = primphylo$tip.label)) %>%
  arrange(species)

plot(primphylo, show.tip.label = FALSE)

# -----------------------
# Models
# -----------------------
form <- range_inf_m10 ~ logNe

# OLS
model_OLS <- lm(form, data = df)

# PGLS Brownian
model_PGLS_brownian <- gls(
  form, data = df,
  correlation = corBrownian(phy = primphylo, form = ~species),
  method = "ML"
)

# PGLS Grafen
model_PGLS_grafen <- gls(
  form, data = df,
  correlation = corGrafen(phy = primphylo, form = ~species, value = 0.8, fixed = FALSE),
  method = "ML"
)

# -----------------------
# Summary table (AIC, pred R^2, coefficients)
# -----------------------
predR2 <- function(y, fit) cor(y, fit, use = "complete.obs")^2

r2s <- c(
  predR2(df$range_inf_m10, fitted(model_OLS)),
  predR2(df$range_inf_m10, fitted(model_PGLS_brownian)),
  predR2(df$range_inf_m10, fitted(model_PGLS_grafen))
)

pvals <- rbind(
  summary(model_OLS)$coefficients[, "Pr(>|t|)"],
  summary(model_PGLS_brownian)$tTable[, "p-value"],
  summary(model_PGLS_grafen)$tTable[, "p-value"]
)

tbl_sum <- rbind(
  round(coef(model_OLS), 3),
  round(coef(model_PGLS_brownian), 3),
  round(coef(model_PGLS_grafen), 3)
)

models_spec <- c("OLS", "PGLS Brownian", "PGLS Grafen")
models_aic <- c(
  round(AIC(model_OLS), 2),
  round(AIC(model_PGLS_brownian), 2),
  round(AIC(model_PGLS_grafen), 2)
)

out_tbl <- cbind(
  `Regression Model` = models_spec,
  AIC = models_aic,
  `R^2` = round(r2s, 3),
  Intercept = tbl_sum[, 1],
  `p(Intercept)` = signif(pvals[, 1], 3),
  log10Ne = tbl_sum[, 2],
  `p(log10Ne)` = signif(pvals[, 2], 3)
)

write_csv(as.data.frame(out_tbl), "scratch/pgls_summary.csv")
