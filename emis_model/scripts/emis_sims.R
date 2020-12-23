# NH3 emission from mink slurry
# S. Hafner

temp.c <- 10

# Use eqSpec to calculate equilibrium speciation
# Before Ca(OH)2 addition
equ <- eqSpec(tot = c(H. = 0.227, NH3 = 0.470, H2CO3 = 0.250, K. = 0, Na. = 0, Cl. = 0.227, HAc = 0), temp.c = temp.c, pH = 7.5, adjpH = 'HCl', of = 'all')

# After Ca(OH)2 addition
# But note that KOH is used for pH
eqt <- eqSpec(tot = c(H. = 0.227, NH3 = 0.470, H2CO3 = 0.250, K. = 0, Na. = 0, Cl. = 0.227, HAc = 0), temp.c = temp.c, pH = 12, adjpH = 'NaOH', of = 'all')

# Some info on arguments: tot = molalities (concentration in mol/kg water) of "components" (elements here, but not below, each represented by a "master species"), temp.c = temperature in C, of = output format ('all' is the most complete).

# Here we will extract total component concentrations for use in kinetic model (only difference from tot argument in call above is partitioning of H2CO3)
itot.u <- equ$totk
itot.t <- eqt$totk
itot.te <- eqt$tot
# itot is initial total concentration of components (including distribution of IC between H2CO3 and CO2, which are not necessarily in equilibrium in simulation below)

# Now for the emission model
# 1 m total thickness, but well-mixed below 1 mm
layers <- list(l20um = rep(5E-6, 4),
               l100um = c(rep(1E-5, 10)),
               #l1mm = c(rep(1E-5, 5), rep(1E-4, 5), 0.00045),
               l1cm = c(rep(1E-5, 5), rep(1E-4, 5), rep(9.45E-4, 10)),
               #l3cm = c(rep(1E-5, 5), rep(1E-4, 5), rep(9.45E-4, 10), rep(0.01, 2)),
               l5cm = c(rep(1E-5, 5), rep(1E-4, 5), rep(9.45E-4, 10), rep(0.01, 4)),
               l3m = c(rep(1E-5, 5), rep(1E-4, 5), rep(9.45E-4, 10), rep(0.01, 9), rep(0.1, 29)))

# Total thickness
thk <- 3

times <- list(t.mix = seq(0, 5 * 3600, length.out = 100),
              t.store = seq(0, 5 * 30 * 86400, length.out = 100))

mtc <- c(straw = 0.005, 
         typical = 0.01)

predsu <- list()
predst <- list()

emis <- data.frame()

# NTS: Could have run for only one set of times, as long as all required durations were present

cc <- 0
for (i in layers) {
  for (j in mtc) {
    for (k in times) {

    cc <- cc + 1
    cat(cc, '\n')

    lb <- if (sum(i) >= 3) 'nf' else 'mx'

    #k <- k[1:3]
    cat('Untreated simulation . . . .\n')
    pp <- kinEmisDiffMod(c.thk = i, 
                         h.m = j, p.CO2.a = 4E-6, times = k,
                         temp.c = temp.c, tot = itot.u,
                         lb = lb, thk = thk)

    predsu[[cc]] <- pp
    ee <- pp$emis
    ee$id <- paste(cc, 'untreated')
    ee$treated <- 0
    ee$film <- sum(i)
    ee$mtc <- j
    ee$dur = max(k)
    emis <- rbind(emis, ee)

    cat('Treated simulation . . . .\n')
    pp <- eqEmisDiffMod(c.thk = i, 
                        h.m = j, p.CO2.a = 4E-6, times = k,
                        temp.c = temp.c, tot = itot.te,
                        lb = lb, thk = thk)

    predst[[cc]] <- pp
    ee <- pp$emis
    ee$id <- paste(cc, 'treated')
    ee$treated <- 1
    ee$film <- sum(i)
    ee$mtc <- j
    ee$dur = max(k)
    emis <- rbind(emis, ee)

    }
  }
}

# Add/modify some variables
emis$j.NH3.gm2h <- 14 * 3600 * emis$j.NH3
emis$time.h <- emis$t / 3600
emis$time.d <- emis$t / 86400
emis$e.NH3.rel <- 100 * emis$e.NH3.rel

emis$condition <- ifelse(emis$treated == 1, 'Treated (pH 12)', 'Untreated (pH 7.5)')
emis$condition <- factor(emis$condition, levels = c('Untreated (pH 7.5)', 'Treated (pH 12)'))

emis$env <- factor(substr(emis$id, 1, 2))
