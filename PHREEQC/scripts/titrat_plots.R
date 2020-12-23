# Plots titration simulations

##system('phreeqc titrat_sim1.phrq')

dat <- read.table('../output/titrat1.sel', header = TRUE)

# Ca in mmol/kgw
dat$lime <- 1000 * (dat$Ca + dat$Aragonite) 

# Without aragonite
d1 <- subset(dat, sim %in% c(1, 4))
# With aragonite
d2 <- subset(dat, sim %in% c(1, 5))

pdf('../../figs/titrat1.pdf', height = 4, width = 4)

  plot(pH ~ lime, data = d1, type = 'l', 
       xlab = expression(Ca(OH)[2]~('mmol'/'kg'[w])), ylab = 'Slurry pH')
  lines(pH ~ lime, data = d2, col = 'gray45')
  grid()

dev.off()
