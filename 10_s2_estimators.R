#### S^2 novērtējuma tests

# Packages ####

require(foreach)
require(ggplot2)
require(doMC)
require(data.table)
require(sampling)

require(xlsx)

# require(openxlsx)
# Sys.setenv(R_ZIPCMD = "zip")
# Sys.getenv("R_ZIPCMD")


# Multi core ####

registerDoMC(cores = 8)
options(stringsAsFactors = F)


# Reset ####
rm(list = ls())
gc()



# Folders ####

dir.root <- "/home/math/T/SMKD/MND/R/s2/Results"


# S^2 estimation
s2 <- function(y, w = rep(1, length(y))) {
  N <- sum(w)
  n <- length(y)
  s2 <- (sum(y ^ 2 * w) - sum(y * w) ^ 2 / N) / (N - 1)
  return(s2)
}

s2_alt <- function(y, w = rep(1, length(y))) {
  N <- sum(w)
  n <- length(y)
  s2 <- (N - 1) / N * n / (n - 1) *
    (sum(y ^ 2 * w) - sum(y * w) ^ 2 / N) / (N - 1)
  return(s2)
}

s2_lap_SRS <- function(y, w = rep(1, length(y))) {
  N <- sum(w)
  n <- length(y)

  Var_Y_vid <- (1 - n / N) / n * var(y)

  Y_vid_sq <- (sum(y * w) / N) ^ 2 - Var_Y_vid

  s2 <- (sum(y ^ 2 * w) - Y_vid_sq * N) / (N - 1)
  return(s2)
}

s2_lap_UPS <- function(y, w = rep(1, length(y)), pikl) {
  N <- sum(w)
  n <- length(y)

  Var_Y_vid <- varHT(y, pikl, method = 1) / N ^ 2

  Y_vid_sq <- (sum(y * w) / N) ^ 2 - Var_Y_vid

  s2 <- (sum(y ^ 2 * w) - Y_vid_sq * N) / (N - 1)
  return(s2)
}


# Data sim ####

N <- 1000

# set.seed(234)

pop <- data.table(y_norm = rnorm(N),
                  y_chisq = rchisq(N, 2),
                  y_binom = rbinom(N, 1, .4))

ynames <- grep("^y", names(pop), value = T)

pop[, x_chisq := y_chisq + rnorm(.N)][, x_chisq := x_chisq - min(x_chisq) + 1]

ggplot(pop, aes(x_chisq, y_chisq)) + geom_point() + theme_bw()

lapply(pop, hist)

S2 <- pop[, lapply(.SD, var), .SDcols = ynames]
S2



# Izlases apjoms
n <- 100

# Izasē iekļūšanas varbūtības
pop[, pik_SRS := n / N]
pop[, pik_UPS := inclusionprobabilities(x_chisq, n)]
pop

pikl <- UPmaxentropypi2(pop$pik_UPS)
dim(pikl)
pikl[1:5, 1:5]
all.equal(diag(pikl), pop$pik_UPS)


ggplot(pop, aes(pik_SRS, y_chisq)) + geom_point() + theme_bw()
ggplot(pop, aes(pik_UPS, y_chisq)) + geom_point() + theme_bw()



# Simulācija ####

K <- 80

res <- foreach(i = 1:K,
               .combine = function(x, y) rbindlist(list(x, y))) %dopar% {
  # SRS
  sampl_SRS <- sample(N, n)

  # UPS
  s <- pop[, UPmaxentropy(pik_UPS)] * (1:N)
  sampl_UPS <- s[s > 0]

  pikl_sam <- pikl[sampl_UPS, sampl_UPS]

  res_SRS_wgh <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2,
                                       w = 1 / pik_SRS)]
  res_SRS_wgh[, sample := "SRS"]
  res_SRS_wgh[, s2_est := "wgh"]

  res_SRS_alt <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2_alt,
                                       w = 1 / pik_SRS)]
  res_SRS_alt[, sample := "SRS"]
  res_SRS_alt[, s2_est := "alt"]

  res_SRS_lap <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2_lap_SRS,
                                       w = 1 / pik_SRS)]
  res_SRS_lap[, sample := "SRS"]
  res_SRS_lap[, s2_est := "lap"]

  res_UPS_wgh <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2,
                                       w = 1 / pik_UPS)]
  res_UPS_wgh[, sample := "UPS"]
  res_UPS_wgh[, s2_est := "wgh"]

  res_UPS_alt <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2_alt,
                                       w = 1 / pik_UPS)]
  res_UPS_alt[, sample := "UPS"]
  res_UPS_alt[, s2_est := "alt"]

  res_UPS_lap <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2_lap_UPS,
                                       w = 1 / pik_UPS, pikl = pikl_sam)]
  res_UPS_lap[, sample := "UPS"]
  res_UPS_lap[, s2_est := "lap"]

  res <- rbindlist(list(res_SRS_wgh, res_SRS_alt, res_SRS_lap,
                        res_UPS_wgh, res_UPS_alt, res_UPS_lap))
  res[, id := i]
  return(res)
}

res

res2 <- melt(res, measure.vars = ynames, variable.factor = F)
res2

sapply(res2, class)



# Īstās vērtības
S2

S2_tab <- melt(S2, measure.vars = ynames, value.name = "S2",
               variable.factor = F)
S2_tab



# Vidējais un sd
tab <- res2[, .(mean = mean(value), sd = sd(value)),
            keyby = list(variable, sample, s2_est)]
tab

tab <- merge(tab, S2_tab, by = "variable")
tab

setkey(tab, sample, variable, s2_est)
tab


# Novirze
tab[, bias := mean - S2]
tab


# MSE
tab[, MSE := bias ^ 2 + sd ^ 2]


# N un n
tab[, N := N]
tab[, n := n]
tab[, K := K]


# args(data.table:::print.data.table)
# args(print.data.frame)

tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 6) else x)]



# Plot ####

res2

gg <- function(sam, variab) {
  ggplot(res2[sample == sam & variable == variab]) +
    geom_density(aes(value, colour = s2_est, linetype = s2_est)) +
    geom_vline(aes(xintercept = S2),
               data = S2_tab[variable == variab], size = 1) +
    geom_vline(aes(xintercept = mean, colour = s2_est, linetype = s2_est),
               data = tab[sample == sam & variable == variab]) +
    ggtitle(paste(sam, variab)) +
    theme_bw()
}

tmp <- tab[, .N, keyby = list(sample, variable)]

grafiki <- mapply(gg, tmp$sample, tmp$variable, SIMPLIFY = F)

grafiki



# Save ####

fname <- gsub("-|:| ", "_", Sys.time())
fname

fname.tab <- file.path(dir.root, paste0("tabl_", fname, ".csv"))
fname.xls <- file.path(dir.root, paste0("tabl_", fname, ".xlsx"))
fname.ggp <- file.path(dir.root, paste0("plot_", fname, ".pdf"))

write.table(tab, file = fname.tab,
            quote = T, sep = ";", row.names = F, qmethod = "double")

# write.xlsx(tab, file = fname.xls, row.names = F)

names(tab)

wb <- createWorkbook()
sh <- createSheet(wb, sheetName = "tab")
csh <- CellStyle(wb) + Font(wb, isBold = TRUE) + Alignment(h = "ALIGN_CENTER")
addDataFrame(tab, sh, row.names = F, colnamesStyle = csh)
saveWorkbook(wb, file = fname.xls)


# wb <- loadWorkbook(fname.xls)
# sheets <- getSheets(wb)
# autoSizeColumn(sheets[[1]], colIndex = 1:ncol(tab))
# saveWorkbook(wb, "Final.xlsx")


cairo_pdf(fname.ggp, width = 16, height = 9, onefile = T)
grafiki
dev.off()
