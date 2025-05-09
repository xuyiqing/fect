#############################################################
## 1) Load packages and data
#############################################################
library(dplyr)
library(fixest)
library(did)
library(PanelMatch)
library(DIDmultiplegtDYN)
library(ggplot2)
library(panelView)

# If your updated did_wrapper() and esplot() functions are saved in a file,
# uncomment and adjust the line below:
# source("did_wrapper_esplot.R")

data(fect)
df <- hh2019
head(df)

# Visualize the panel data
panelview(nat_rate_ord ~ indirect, data = df, index = c("bfs", "year"),
          xlab = "Year", ylab = "Unit", display.all = TRUE,
          gridOff = TRUE, by.timing = TRUE)

# Main variables
Y <- "nat_rate_ord"
D <- "indirect"
index <- c("bfs", "year")
#############################################################
## 2) TWFE example
#############################################################
res_twfe <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "twfe",
  se     = "default"
)
cat("\n>>> TWFE results:\n")
cat("ATT:", res_twfe$ATT, "SE:", res_twfe$ATT_se,
    "CI:", res_twfe$CI_lower, "to", res_twfe$CI_upper, "\n")

p_twfe <- esplot(data = res_twfe,
                 main = "TWFE event‐study", xlim = c(-12,10), show.count =TRUE)
print(p_twfe)

#############################################################
## 3) Stacked DID ("st") example
#############################################################
res_st <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "st",
  se     = "default"
)
cat("\n>>> Stacked DID results:\n")
cat("ATT:", res_st$ATT, "SE:", res_st$ATT_se, "\n")

p_st <- esplot(data = res_st,
               main = "Stacked DID event‐study", xlim = c(-12,10))
print(p_st)

#############################################################
## 4) Sun & Abraham's IW ("iw") example
#############################################################
res_iw <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "iw",
  se     = "default"
)
cat("\n>>> IW (Sun & Abraham) results:\n")
cat("ATT:", res_iw$ATT, "SE:", res_iw$ATT_se, "\n")

p_iw <- esplot(data = res_iw,
               main = "IW enumerated expansions", xlim = c(-12,10))
print(p_iw)

#############################################################
## 5) CSDID examples: "cs_never" and "cs_notyet"
#############################################################
### 5.1) cs_never
res_csnever <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "cs_never",
  se     = "default"
)
cat("\n>>> CSDID (never-treated) results:\n")
cat("ATT:", res_csnever$ATT, "SE:", res_csnever$ATT_se,
    "CI:", res_csnever$CI_lower, "to", res_csnever$CI_upper, "\n")

p_csnotyet <- esplot(data = res_csnever,
                     main = "CSDID not-yet-treated ES", xlim = c(-12,10))
print(p_csnotyet)

### 5.2) cs_notyet
res_csnotyet <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "cs_notyet",
  se     = "default"
)
cat("\n>>> CSDID (not-yet-treated) results:\n")
cat("ATT:", res_csnotyet$ATT, "SE:", res_csnotyet$ATT_se,
    "CI:", res_csnotyet$CI_lower, "to", res_csnotyet$CI_upper, "\n")

p_csnotyet <- esplot(data = res_csnotyet,
                     main = "CSDID not-yet-treated ES", xlim = c(-12,10))
print(p_csnotyet)

#############################################################
## 6) DIDmultiplegtDYN ("didm") example
#############################################################
res_didm <- did_wrapper(
  data         = df,
  Y            = Y,
  D            = D,
  index        = index,
  method       = "didm",
  didm.effects = 12,
  didm.placebo = 9,
  se           = "default"
)
cat("\n>>> DIDmultiplegtDYN results:\n")
cat("ATT:", res_didm$ATT, "SE:", res_didm$ATT_se,
    "CI:", res_didm$CI_lower, "to", res_didm$CI_upper, "\n")

p_didm <- esplot(data = res_didm,
                 main = "did_multiplegt_dyn ES", xlim = c(-12,9))
print(p_didm)

#############################################################
## 7) TWFE with cluster bootstrap example
#############################################################
res_twfe_boot <- did_wrapper(
  data   = df,
  Y      = Y,
  D      = D,
  index  = index,
  method = "twfe",
  se     = "boot",
  nboots = 50
)
cat("\n>>> TWFE with cluster bootstrap:\n")
cat("ATT:", res_twfe_boot$ATT, "\n")
cat("SE (boot):", res_twfe_boot$ATT_se, "\n")
cat("CI:", res_twfe_boot$CI_lower, "to", res_twfe_boot$CI_upper, "\n")



rmarkdown::render("/Users/rivka/Dropbox/fectbook/02-fect.Rmd", output_format = "html_document")
rmarkdown::render("/Users/rivka/Dropbox/fectbook/03-plots.Rmd", output_format = "html_document")
rmarkdown::render("/Users/rivka/Dropbox/fectbook/04-gsynth.Rmd", output_format = "html_document")
rmarkdown::render("/Users/rivka/Dropbox/fectbook/05-panel.Rmd", output_format = "html_document")
