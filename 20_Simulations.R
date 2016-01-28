#### S^2 novērtējuma tests

# Packages ####
require(foreach)
require(ggplot2)
require(doMC)
require(data.table)
require(sampling)
require(xlsx)


# Multi core ####
registerDoMC(cores = 4)

options(stringsAsFactors = F)


# Reset ####
rm(list = ls())
gc()


# Number of iterations
K <- 100e3


# Load functions
source("10_s2_estimators.R")


# Local folder for results ####
dir.resu <- "/home/math/T/SMKD/MND/R/s2/Results"



# Population data ####

N <- 10e3

set.seed(234)

pop <- data.table(rnd = runif(N),
                  y_norm = sort(round(rnorm(N) * 10)),
                  y_chisq = sort(round(rchisq(N, 2) * 10)),
                  y_binom = sort(rbinom(N, 1, .4)))

ynames <- grep("^y", names(pop), value = T)

pop[, x_chisq := y_chisq + rnorm(.N)]
pop[, x_chisq := round(x_chisq - min(x_chisq) + 1)]

pop

ggplot(pop, aes(x_chisq, y_chisq)) + geom_point() + theme_bw()
ggplot(pop, aes(x_chisq, y_norm)) + geom_point() + theme_bw()


# Order in random order
setorder(pop, rnd)
pop


# lapply(pop, hist)

S2 <- pop[, lapply(.SD, var), .SDcols = ynames]
S2

# S2 <- pop[, lapply(.SD, s2a), .SDcols = ynames]
# S2



# Izlases apjoms
n <- 10


# Izasē iekļūšanas varbūtības
# SRS
pop[, pik_SRS := n / N]

# UPS
pop[, pik_UPS := inclusionprobabilities(x_chisq, n)]
pop[, summary(pik_UPS)]


# ggplot(pop, aes(pik_SRS, y_chisq)) + geom_point() + theme_bw()
# ggplot(pop, aes(pik_UPS, y_chisq)) + geom_point() + theme_bw()
# ggplot(pop, aes(pik_UPS, x_chisq)) + geom_point() + theme_bw()



# Simulācija ####

t1 <- Sys.time()

res <- foreach(i = 1:K,
               .combine = function(x, y) rbindlist(list(x, y))) %dopar% {
  # SRS
  sampl_SRS <- sample(N, n)

  # UPS
  s <- pop[, UPrandomsystematic(pik_UPS)]
  sampl_UPS <- (1:N)[s == 1]

  res_SRS_a <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2a)]
  res_SRS_a[, sample := "SRS"]
  res_SRS_a[, s2_est := "s2a"]

  res_SRS_b <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2b,
                                     w = 1 / pik_SRS)]
  res_SRS_b[, sample := "SRS"]
  res_SRS_b[, s2_est := "s2b"]

  res_SRS_c <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2c,
                                     w = 1 / pik_SRS)]
  res_SRS_c[, sample := "SRS"]
  res_SRS_c[, s2_est := "s2c"]

  res_SRS_d <- pop[sampl_SRS, lapply(.SD[, ynames, with = F], s2d,
                                     w = 1 / pik_SRS)]
  res_SRS_d[, sample := "SRS"]
  res_SRS_d[, s2_est := "s2d"]

  res_UPS_a <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2a)]
  res_UPS_a[, sample := "UPS"]
  res_UPS_a[, s2_est := "s2a"]

  res_UPS_b <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2b,
                                       w = 1 / pik_UPS)]
  res_UPS_b[, sample := "UPS"]
  res_UPS_b[, s2_est := "s2b"]

  res_UPS_c <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2c,
                                       w = 1 / pik_UPS)]
  res_UPS_c[, sample := "UPS"]
  res_UPS_c[, s2_est := "s2c"]

  res_UPS_d <- pop[sampl_UPS, lapply(.SD[, ynames, with = F], s2d,
                                       w = 1 / pik_UPS)]
  res_UPS_d[, sample := "UPS"]
  res_UPS_d[, s2_est := "s2d"]

  res <- rbindlist(list(res_SRS_a, res_SRS_b, res_SRS_c, res_SRS_d,
                        res_UPS_a, res_UPS_b, res_UPS_c, res_UPS_d))
  res[, id := i]
  return(res)
}

t2 <- Sys.time()

run.time <- as.integer(difftime(t2, t1, units = "secs"))
run.time
run.time / K

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
tab[, mse := sqrt(bias ^ 2 + sd ^ 2)]

tab[, min := ifelse(mse == min(mse), "*", ""), by = .(variable, sample)]
tab


# N un n
tab[, N := N]
tab[, n := n]
tab[, K := K]
tab[, tt := run.time]
tab[, at := run.time / K]

tab


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
tmp

grafiki <- mapply(gg, tmp$sample, tmp$variable, SIMPLIFY = F)

grafiki



# Save ####

fname <- gsub("-|:| ", "_", Sys.time())
fname

fname.csv <- file.path(dir.resu, paste0("tabl_", fname, ".csv"))
fname.xls <- file.path(dir.resu, paste0("tabl_", fname, ".xlsx"))
fname.pdf <- file.path(dir.resu, paste0("plot_", fname, ".pdf"))


# CSV
write.table(tab, file = fname.csv,
            quote = T, sep = ";", row.names = F, qmethod = "double")


# XLSX
names(tab)

wb <- createWorkbook()
sh <- createSheet(wb, sheetName = "tab")
csh <- CellStyle(wb) + Font(wb, isBold = TRUE) + Alignment(h = "ALIGN_CENTER")
addDataFrame(tab, sh, row.names = F, colnamesStyle = csh)
saveWorkbook(wb, file = fname.xls)


# PDF
cairo_pdf(fname.pdf, width = 16, height = 9, onefile = T)
grafiki
dev.off()
