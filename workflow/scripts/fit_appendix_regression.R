# Fit OLS, PGLS Brownian, and PGLS Grafen regressions of the strongly
# deleterious DFE bin (S in (-inf, -10]) on log10(Ne) + generation_time,
# and emit a LaTeX table fragment matching the appendix format.
#
# Mirrors the regression chunks of
# https://github.com/tbata/PLS_DFE/blob/main/Report.Rmd
# but skips the visuals and writes a single .tex fragment per DFE model.

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

if (exists("snakemake")) {
  csv_file   <- snakemake@input[["csv"]]
  tree_file  <- snakemake@input[["tree"]]
  out_tex    <- snakemake@output[[1]]
  caption    <- snakemake@params[["caption"]]
  label      <- snakemake@params[["label"]]
} else {
  # local testing
  csv_file  <- "results/tables/regression/reg_vars.unfolded.8.gamma.full.noeps.csv"
  tree_file <- "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
  out_tex   <- "scratch/regression.gamma.full.tex"
  caption   <- "DFE model: gamma + exponential allowing for a proportion of beneficial mutations."
  label     <- "tab:regression_gamma_full"
}

# ---- load + filter + aggregate by species (mirror Report.Rmd) ----
raw <- read.csv(csv_file, check.names = FALSE)
raw <- raw[!is.na(raw$generation_time), ]
raw$species <- ifelse(
  sapply(gregexpr("_", raw$population), function(x) sum(x > 0)) >= 2,
  sub("(_[^_]+)$", "", raw$population),
  raw$population
)
raw$logNe <- log10(raw$Ne)

agg <- aggregate(
  cbind(Ne, logNe, generation_time,
        range_inf_m10 = raw$`range_inf_-10`,
        range_m1_0    = raw$`range_-1_0`) ~ species,
  data = raw,
  FUN  = function(x) mean(x, na.rm = TRUE)
)
# aggregate returns matrices for multi-column LHS — flatten
agg <- data.frame(species = agg$species,
                  Ne              = agg[, "Ne"],
                  logNe           = agg[, "logNe"],
                  generation_time = agg[, "generation_time"],
                  range_inf_m10   = agg[, "range_inf_m10"],
                  range_m1_0      = agg[, "range_m1_0"])

# ---- intersect with tree tips ----
phy    <- read.tree(tree_file)
common <- intersect(phy$tip.label, agg$species)
phy    <- drop.tip(phy, setdiff(phy$tip.label, common))
df     <- agg[agg$species %in% phy$tip.label, ]
df$species <- factor(df$species, levels = phy$tip.label)
df     <- df[order(df$species), ]

# ---- fit three models ----
form <- range_inf_m10 ~ logNe + generation_time

m_ols     <- lm(form, data = df)
m_brown   <- gls(form, data = df,
                 correlation = corBrownian(phy = phy, form = ~species),
                 method = "ML")
m_grafen  <- gls(form, data = df,
                 correlation = corGrafen(phy = phy, form = ~species,
                                         value = 0.8, fixed = FALSE),
                 method = "ML")

pred_r2 <- function(y, fit) cor(y, fit, use = "complete.obs")^2

r2s <- c(
  pred_r2(df$range_inf_m10, fitted(m_ols)),
  pred_r2(df$range_inf_m10, fitted(m_brown)),
  pred_r2(df$range_inf_m10, fitted(m_grafen))
)

aics <- c(AIC(m_ols), AIC(m_brown), AIC(m_grafen))

coefs <- rbind(
  coef(m_ols),
  coef(m_brown),
  coef(m_grafen)
)

pvals <- rbind(
  summary(m_ols)$coefficients[, "Pr(>|t|)"],
  summary(m_brown)$tTable[, "p-value"],
  summary(m_grafen)$tTable[, "p-value"]
)

# Grafen rho is stored in corStruct as a transformed param — extract & back-transform.
# corGrafen() stores log(rho) when fixed = FALSE, so exp() recovers rho.
grafen_rho <- as.numeric(coef(m_grafen$modelStruct$corStruct, unconstrained = FALSE))

# ---- format helpers ----
fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

fmt_p <- function(p) {
  # very small p in scientific LaTeX; otherwise plain
  if (is.na(p))           return("NA")
  if (p == 0)             return("$<10^{-300}$")
  if (p < 1e-3) {
    e <- floor(log10(p))
    m <- p / 10^e
    return(sprintf("$%.2f\\times10^{%d}$", m, e))
  }
  formatC(signif(p, 3), format = "g", digits = 3)
}

row_tex <- function(model_label, i) {
  paste(
    model_label,
    fmt_num(aics[i], 2),
    fmt_num(r2s[i],  3),
    fmt_num(coefs[i, 1], 3),
    fmt_p(pvals[i, 1]),
    fmt_num(coefs[i, 2], 3),
    fmt_p(pvals[i, 2]),
    fmt_num(coefs[i, 3], 3),
    fmt_p(pvals[i, 3]),
    sep = " & "
  )
}

rows <- c(
  paste(row_tex("OLS", 1), "\\\\"),
  paste(row_tex("PGLS Brownian", 2), "\\\\"),
  paste(row_tex(sprintf("PGLS Grafen ($\\rho=%.2f$)", grafen_rho), 3), "\\\\")
)

tex <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\begin{tabular}{lrrrrrrrr}",
  "\\toprule",
  "Model & AIC & $R^2$ & Intercept & p(Int) & $\\log_{10}N_e$ & p($\\log_{10}N_e$) & gen\\_time & p(gen\\_time) \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("\\caption{%s}", caption),
  sprintf("\\label{%s}", label),
  "\\end{table}",
  ""
)

writeLines(tex, out_tex)
cat(sprintf("Wrote %s (n_species = %d, Grafen rho = %.3f)\n",
            out_tex, nrow(df), grafen_rho))
