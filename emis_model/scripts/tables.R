# Tables for report, for most likely conditions

# Mixing
mxe <- data.frame(hours = 1:5)
x <- subset(emis, mtc == 0.01 & film == 0.0001 & dur == 18000 & treated == 1)
mxe$e.NH3.rel <- signif(approx(x$time.h, x$e.NH3.rel, xout = 1:5)$y, 2)
mxe$e.NH3 <- round(14.0067 * approx(x$time.h, x$e.NH3, xout = 1:5)$y, 0)
mxe$j.NH3 <- signif(diff(c(0, mxe$e.NH3))/1, 2)

write.table(mxe, '../output/mix_summ.csv', sep = '|')

# Mixing high
mxm <- data.frame(hours = 1:5)
x <- subset(emis, mtc == 0.01 & film == 0.00002 & dur == 18000 & treated == 1)
mxm$e.NH3.rel <- signif(approx(x$time.h, x$e.NH3.rel, xout = 1:5)$y, 2)
mxm$e.NH3 <- round(14.0067 * approx(x$time.h, x$e.NH3, xout = 1:5)$y, 0)
mxm$j.NH3 <- signif(diff(c(0, mxm$e.NH3))/1, 2)

write.table(mxm, '../output/mix_summ_high.csv', sep = '|')

# Long-term
lte <- data.frame(months = 1:5)
x <- subset(emis, mtc == 0.005 & film == 0.05 & dur == 12960000 & treated == 1)
lte$e.NH3.rel <- signif(approx(x$time.d, x$e.NH3.rel, xout = 1:5 * 30)$y, 2)
lte$e.NH3 <- round(14.0067 * approx(x$time.d, x$e.NH3, xout = 1:5 * 30)$y, 0)
lte$j.NH3 <- signif(diff(c(0, lte$e.NH3))/(30 * 24), 2)

write.table(lte, '../output/store_summ.csv', sep = '|')

# High long-term (1 cm film)
ltm <- data.frame(months = 1:5)
x <- subset(emis, mtc == 0.01 & film == 0.01 & dur == 12960000 & treated == 1)
ltm$e.NH3.rel <- signif(approx(x$time.d, x$e.NH3.rel, xout = 1:5 * 30)$y, 2)
ltm$e.NH3 <- round(14.0067 * approx(x$time.d, x$e.NH3, xout = 1:5 * 30)$y, 0)
ltm$j.NH3 <- signif(diff(c(0, ltm$e.NH3))/(30 * 24), 2)

write.table(ltm, '../output/store_high_summ.csv', sep = '|')
