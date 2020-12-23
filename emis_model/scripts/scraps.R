

    head(ee)

    # Untreated
    predsu <- kinEmisDiffMod(c.thk = i, 
                            h.m = j, p.CO2.a = 4E-6, times = t.mix,
                            temp.c = temp.c, tot = itot.u,
                            lb = 'mx', thk = thk)



# Reference check with thinest film 100 um and mixing
pred1 <- kinEmisDiffMod(c.thk = layers[['l100um']], 
                        h.m = mtc[[1]], p.CO2.a = 4E-6, times = t.mix,
                        temp.c = temp.c, tot = itot.u,
                        lb = 'mx', thk = thk)

# Treated check with 100 um film and mixing
pred2 <- eqEmisDiffMod(c.thk = layers[['l100um']], 
                        h.m = mtc[[1]], p.CO2.a = 4E-6, times = t.mix,
                        temp.c = temp.c, tot = itot.te,
                        lb = 'mx', thk = thk)

# Treated check with full depth film and typical wind
pred3 <- eqEmisDiffMod(c.thk = layers[['l3m']], 
                        h.m = mtc[[1]], p.CO2.a = 4E-6, times = t.store,
                        temp.c = temp.c, tot = itot.te,
                        lb = 'nf', thk = thk)


tail(pred3$emis)
pred3$ph
x
(-j.CO2*c(1, rep(0, n.cells-1)))  
( diff(c(j.diff.out[, 'H2CO3'], j.diff.lb['H2CO3'])) )/(c.thk*1000)

traceback()

#pred2 <- eqEmisDiffMod(c.thk = layers1cm, 
#                        h.m = 1E-2, p.CO2.a = 4E-6, times = times,
#                        temp.c = 20, tot = itot2e,
#                        lb = 'mx', thk = thk)
#
#pred4 <- eqEmisDiffMod(c.thk = layers1mm, 
#                        h.m = 1E-2, p.CO2.a = 4E-6, times = times,
#                        temp.c = 20, tot = itot2e,
#                        lb = 'mx', thk = thk)
#

# Extract times
tt <- pred1$times
# position (center of cells) converted to mm
z <- pred1$pos*1000
# and pH
pH1 <- pred1$ph
pH2 <- pred2$ph
pH4 <- pred4$ph

# Plot pH
matplot(t(pH1), z, type = 'l', ylim = c(10, 0))

# "Surface" pH over time
plot(tt, pH1[, 1], type = 'o')
plot(tt, pH2, type = 'o')

# Get fluxes
emis1 <- pred1$emis
emis2 <- pred2$emis
emis4 <- pred4$emis

# Flux in g N/m2-h
emis1$j.NH3.gm2h <- 3600 * emis1$j.NH3
emis2$j.NH3.gm2h <- 3600 * emis2$j.NH3
emis4$j.NH3.gm2h <- 3600 * emis4$j.NH3

head(emis1)
emis1$time.h <- emis1$t / 3600
emis2$time.h <- emis2$t / 3600
emis4$time.h <- emis4$t / 3600

# Plot NH3 flux
pdf('../../figs/flux1.pdf', height = 4, width = 4)
  plot(j.NH3.gm2h ~ time.h, data = emis1, type = 'l',
       xlab = 'Time (h)', ylab = expression(NH[3]~'flux'~(g~m^'-2'~s^'-1')))
dev.off()

pdf('../../figs/flux2.pdf', height = 4, width = 4)
  plot(j.NH3.gm2h ~ time.h, data = emis2, type = 'l',
       xlab = 'Time (h)', ylab = expression(NH[3]~'flux'~(g~m^'-2'~s^'-1')))
dev.off()


plot(j.NH3.gm2h ~ time.h, data = emis2[-1, ])
# and cumulative emission
plot(e.NH3 ~ t, data = pred1$emis)






