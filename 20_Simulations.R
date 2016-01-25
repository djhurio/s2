#### S^2 novērtējuma tests

# Packages ####

require(foreach)
require(ggplot2)
require(doMC)
require(data.table)
require(sampling)
# require(xlsx)


# Multi core ####

registerDoMC(cores = 1)
options(stringsAsFactors = F)


# Reset ####
rm(list = ls())
gc()


# Load functions
source("10_s2_estimators.R")


# Folders ####
dir.resu <- "/home/math/T/SMKD/MND/R/s2/Results"



# Data sim ####

N <- 1000

set.seed(234)

pop <- data.table(y_norm = rnorm(N),
                  y_chisq = rchisq(N, 2),
                  y_binom = rbinom(N, 1, .4))

ynames <- grep("^y", names(pop), value = T)

pop[, x_chisq := y_chisq + rnorm(.N)][, x_chisq := x_chisq - min(x_chisq) + 1]

ggplot(pop, aes(x_chisq, y_chisq)) + geom_point() + theme_bw()

lapply(pop, hist)

S2 <- pop[, lapply(.SD, var), .SDcols = ynames]
S2

S2 <- pop[, lapply(.SD, s2a), .SDcols = ynames]
S2



# Izlases apjoms
n <- 100

# Izasē iekļūšanas varbūtības
pop[, pik_SRS := n / N]
pop[, pik_UPS := inclusionprobabilities(x_chisq, n)]
pop

ggplot(pop, aes(pik_SRS, y_chisq)) + geom_point() + theme_bw()
ggplot(pop, aes(pik_UPS, y_chisq)) + geom_point() + theme_bw()



# Simulācija ####

K <- 10

res <- foreach(i = 1:K,
               .combine = function(x, y) rbindlist(list(x, y))) %dopar% {
  # SRS
  sampl_SRS <- sample(N, n)

  # UPS
  s <- pop[, UPmaxentropy(pik_UPS)] * (1:N)
  sampl_UPS <- s[s > 0]

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
tmp

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
