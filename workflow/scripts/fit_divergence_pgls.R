# Phylogenetically-controlled alpha-vs-Ne regression for the divergence
# experiment. The plain OLS in assemble_divergence_experiment.py treats the
# (up to) 38 populations as independent, but they collapse onto ~14 references
# and divergence is largely a per-reference quantity, so the OLS p-values are
# anticonservative. Here we aggregate to species and refit alpha ~ log10(Ne)
# under OLS, PGLS Brownian, and PGLS Grafen against the Kuderna tree, exactly
# as in fit_appendix_regression.R / do_regression.R, for both the
# polymorphism-only (alpha_off) and the +divergence (alpha_on) estimates.

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

if (exists("snakemake")) {
  table_file <- snakemake@input[["table"]]
  tree_file  <- snakemake@input[["tree"]]
  out_csv    <- snakemake@output[["csv"]]
  out_tex    <- snakemake@output[["tex"]]
  param      <- snakemake@wildcards[["parametrization"]]
} else {
  # local testing
  table_file <- "results/divergence/experiment/catarrhini/gamma/table.8.csv"
  tree_file  <- "resources/Kuderna/supplementary_files/science.abn7829_data_s4.nex.tree"
  out_csv    <- "scratch/divergence/pgls.gamma.8.csv"
  out_tex    <- "scratch/divergence/pgls.gamma.8.tex"
  param      <- "gamma"
}

# ---- load + aggregate population -> species (mirror fit_appendix_regression.R) ----
raw <- read.csv(table_file, check.names = FALSE)
raw$species <- ifelse(
  sapply(gregexpr("_", raw$population), function(x) sum(x > 0)) >= 2,
  sub("(_[^_]+)$", "", raw$population),
  raw$population
)
raw$logNe <- log10(raw$Ne)

agg <- aggregate(
  cbind(logNe, alpha_off, alpha_on) ~ species,
  data = raw,
  FUN  = function(x) mean(x, na.rm = TRUE),
  na.action = na.pass
)
agg <- data.frame(
  species   = agg$species,
  logNe     = agg[, "logNe"],
  alpha_off = agg[, "alpha_off"],
  alpha_on  = agg[, "alpha_on"]
)

# ---- intersect with tree tips ----
phy    <- read.tree(tree_file)
common <- intersect(phy$tip.label, agg$species)
phy    <- drop.tip(phy, setdiff(phy$tip.label, common))
df     <- agg[agg$species %in% phy$tip.label, ]
df$species <- factor(df$species, levels = phy$tip.label)
df     <- df[order(df$species), ]

cat(sprintf("Aggregated %d populations -> %d species (%d on tree)\n",
            nrow(raw), nrow(agg), nrow(df)))

# ---- fit OLS / PGLS Brownian / PGLS Grafen for one response ----
pred_r2 <- function(y, fit) cor(y, fit, use = "complete.obs")^2

fit_one <- function(response) {
  sub  <- df[is.finite(df[[response]]) & is.finite(df$logNe), ]
  form <- as.formula(sprintf("%s ~ logNe", response))
  tip  <- as.character(sub$species)
  phy_sub <- drop.tip(phy, setdiff(phy$tip.label, tip))

  m_ols <- lm(form, data = sub)
  m_brown <- tryCatch(
    gls(form, data = sub, method = "ML",
        correlation = corBrownian(phy = phy_sub, form = ~species)),
    error = function(e) NULL)
  m_grafen <- tryCatch(
    gls(form, data = sub, method = "ML",
        correlation = corGrafen(phy = phy_sub, form = ~species,
                                value = 0.8, fixed = FALSE)),
    error = function(e) NULL)

  grab <- function(model, label) {
    if (is.null(model)) {
      return(data.frame(response = response, model = label, n = nrow(sub),
                        AIC = NA, R2 = NA, slope = NA, p_slope = NA,
                        intercept = NA, p_intercept = NA, grafen_rho = NA))
    }
    if (inherits(model, "lm")) {
      co <- summary(model)$coefficients
      p  <- co[, "Pr(>|t|)"]
      rho <- NA
    } else {
      co <- summary(model)$tTable
      p  <- co[, "p-value"]
      rho <- if (label == "PGLS Grafen")
        as.numeric(coef(model$modelStruct$corStruct, unconstrained = FALSE))
      else NA
    }
    data.frame(
      response = response, model = label, n = nrow(sub),
      AIC = AIC(model), R2 = pred_r2(sub[[response]], fitted(model)),
      slope = coef(model)[["logNe"]], p_slope = p[["logNe"]],
      intercept = coef(model)[[1]], p_intercept = p[[1]],
      grafen_rho = rho
    )
  }

  rbind(grab(m_ols, "OLS"),
        grab(m_brown, "PGLS Brownian"),
        grab(m_grafen, "PGLS Grafen"))
}

res <- rbind(fit_one("alpha_off"), fit_one("alpha_on"))
res_csv <- res
res_csv[] <- lapply(res_csv, function(c) if (is.numeric(c)) round(c, 4) else c)
write.csv(res_csv, out_csv, row.names = FALSE)

# ---- LaTeX fragment (appendix style) ----
fmt_num <- function(x, d = 3) if (is.na(x)) "NA" else formatC(x, format = "f", digits = d)
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p <= 0) return("$<10^{-15}$")
  if (p < 1e-3) {
    e <- as.integer(floor(log10(p))); m <- p / 10^e
    return(sprintf("$%.2f\\times10^{%d}$", m, e))
  }
  formatC(signif(p, 3), format = "g", digits = 3)
}
mode_label <- function(r) if (r == "alpha_off") "Polymorphism only" else "With divergence"

rows <- c()
for (r in c("alpha_off", "alpha_on")) {
  block <- res[res$response == r, ]
  rows <- c(rows, sprintf("\\multicolumn{6}{l}{\\emph{%s}} \\\\", mode_label(r)))
  for (i in seq_len(nrow(block))) {
    rows <- c(rows, paste(
      block$model[i],
      fmt_num(block$AIC[i], 2),
      fmt_num(block$R2[i], 3),
      fmt_num(block$slope[i], 3),
      fmt_p(block$p_slope[i]),
      if (is.na(block$grafen_rho[i])) "--" else fmt_num(block$grafen_rho[i], 2),
      sep = " & "), "\\\\")
  }
  rows <- c(rows, "\\midrule")
}
rows <- head(rows, -1)  # drop trailing midrule

tex <- c(
  "\\begin{table}[H]", "\\centering",
  "\\begin{tabular}{lrrrrr}", "\\toprule",
  "Model & AIC & $R^2$ & $\\log_{10}N_e$ & p($\\log_{10}N_e$) & Grafen $\\rho$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
  sprintf("\\caption{Phylogenetically-controlled regression of $\\alpha$ on $\\log_{10}N_e$ (%s parametrization), aggregated to %d species. OLS treats species as independent; PGLS Brownian and PGLS Grafen account for shared ancestry via the Kuderna tree.}",
          param, nrow(df)),
  sprintf("\\label{tab:divergence_pgls_%s}", param),
  "\\end{table}", "")
writeLines(tex, out_tex)

cat(sprintf("Wrote %s and %s\n", out_csv, out_tex))
print(res)
