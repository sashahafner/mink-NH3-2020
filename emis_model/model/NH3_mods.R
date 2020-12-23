# Equilibrium and kinetic (CO2) models for NH3 and CO2 emission from slurry and changes in speciation
# Author: S. Hafner
#
# History of revisions
# 2014 OCT 31  nh3_mods.R          Combined code from two separate files for equilibrium and kinetic models
#                                  See comments above each section for details on earlier revisions.
# 2014 NOV 02  nh3_mods.R          Added pH argument to eqSpec, if specified HAc or KOH will be added to match 
#                                  pH based on H balance. Changed while check for ionic strength. Added tot
#                                  check for kinEmisDiffMod.
# 2014 NOV 24  nh3_mods.R          Found problem with ul setting in kinSpec (just that it was too low, -6). Set to
#                                  -5 for now, and added error messages when charge balance is off, but this should
#                                  be improved.
# 2014 NOV 25  nh3_mods.R          Added error check on both optimize calls for when objective is too high.
#                                  Added acid argument to eqSpec so user can select which acid will be used to 
#                                  match entered pH. Choices are "HAc" of "HCl".
# 2014 DEC 01  nh3_mods.R          Made acid/base pH adjustment approach more flexible and replaced "acid" argument with 
#                                  adjpH. Can handle addition or subtraction of acid or base now
# 2015 JAN 03  nh3_mods.R          Changed calculation of relative emission, added imass.NH3 and imass.CO2 calculation
#                                  before emission calculations. Still is a problem for emission for lb = 'cc'.
# 2015 JAN 04  nh3_mods.R          Added 'NaOH' option to adjpH argument. Added ph.results argument to eqEmisDiffMod
#                                  for returning pH values
# 2015 JAN 05  nh3_mods.R          Added signif rounding to 5 digits for time in output matrix row names. Need to
#                                  revisit.
# 2015 MAR 07  nh3_mods.R          Note: made some changes to ch.cells, added Dn, plus other changes to a.globe
#                                  and tot.glob, and related, that decreased speed and were discarded. See
#                                  archived versions.
# 2015 MAR 07  nh3_mods.R          Added variable D with depth. Specify Dd for deep D and zD for depth of change.
#
# 2015 MAR 14  nh3_mods.R          Added mol = molalities of all species to output of kinEmisDiffMod().
#                                  Had to make m.glob for this.
# 2015 APR 20  nh3_mods.R          kinEmisDiffMod only, added h.m.constant logical argument. Set to TRUE to use 
#                                  single h.m value
# 2018 MAR 02  NH3_mods.R          File name change only
# 2019 APRIL 10  NH3_mods.R        Correct global assignment symbols, had added spaces between them.



# Equilibrium-only model
# History of revisions
# Date         File                Description
# 2010 OCT 08  SpecMod_v1.R        First version that really works. Cacluates speciation of NH3 and CO2.
# 2010 OCT 11  SpecMod_v2.R        Trying to incorporate temperature response of equilibrium constants, 
#                                  kinetics, and emission. Seems to work, but is slow.
# 2010 OCT 12  SpecMod_v3.R        Switching to H2CO3 as master C species.
# 2010 OCT 13  SpecMod_v4.R        Trying to speed things up. Kinetic component now works pretty well, but is only single layer.
# 2010 OCT 13  SpecMod_v5.R        Includes transport. Now uses extended Debye-Huckel equation for activity coefficients
# 2010 OCT 14  SpecMod_v6.R        Modifying so CO2 is returned by spec function. Everything seems to work.
# 2010 OCT 18  SpecMod_v7.R        Adding a function for running the complete transport model (was manual)
#                                  Corrected a 1000-fold error in derivative calculations
# 2010 OCT 19  SpecMod_v8.R        Trying to speed up execution
# 2010 OCT 20  SpecMod_v8.R        Corrected an error in dx
# 2010 OCT 27  SpecMod_v12.R       Between this and v8, tried several things for improving speed--all failed.
#                                  This version uses a slightly simplified version of spec.ph for opitimization
# 2010 OCT 28  SpecMod_v13.R       Uses global variables to track tot and a. Only calculates speciation when needed.
#                                  Slightly faster than v12, but results are sensitive to the change in tot that is 
#                                  considered insignificant. If it is too small, ode1D makes many more calls to der.calc, 
#                                  and results in slow runs. Saw this by writing out time at each der.calc call to a file.
# 2010 NOV 03  SpecMod_v14.R       Adding K+, Na+, Cl-, and acetic acid.
# 2010 NOV 04  SpecMod_v15.R       Trying to pull out non-volatile solutes from ODE solver to speed things up
# 2010 NOV 08-09SpecMod_v16.R      Tried to write my own solver that uses simple fixed relative change approach. 
#                                  As feared, it was very slow. Incorporated a new approach where CO2(aq) in the upper layer
#                                  is set to zero. Speeds things up slightly, and is optional.
# 2010 NOV 09  SpecMod_v17.R       Adding speciation-only option, with alternate speciation for simulations without kinetics
# 2010 NOV 09  SpecTransMod_v1.R   Separating from original model. This version does not include kinetics
#                                  Note that species order is the same as in original model. This makes for a slightly different
#                                  approach to the stoichiometry matrix.
# 2010 DEC 13  SpecTransMod_v1.R   Added line to deal with n.calls.spec when only spec is called up without transport
# 2010 DEC 24  SpecTransMod_v1.R   Added effect of ambient CO2
# 2010 DEC 24  SpecTransMod_v2.R   Made some changes to output. Adding option for constant concentration at lower surface.
#                                  Required changes in diffusion calculations--need to check!!!
# 2011 JAN 04  SpecTransMod_v3.R   Improved der.calc code a bit.
# 2011 APR 05  SpecTransMod_v4.R   Making h.m length of 2, for NH3 and CO2, and calculate h.m for CO2 from h.m for NH3.
# 2011 APR 07  SpecTransMod_v4.R   Added if statement so summ is in output only if relevant times are requested.
# 2011 APR 14  SpecTransMod_v5.R   Changing equations for equilibrium constants and kinetic constants to latest version, 
#                                  based on Plummer and Bunsenberg (1982).
# 2011 SEP 26  SpecTransMod_v7.R   Adding CO2.emis argument so CO2 emission can be shut off. Note that the changes in version 6 
#                                  were abandoned. Added pos to output (had been missing).
# 2011 SEP 27   SpecTransMod_v7.R  Replaced dielectric constant calculation and Debye-Huckel A and B with equations from PHREEQC and WATEQ.
# 2013 NOV 26   eqmod.R            Changed file name, replaced tabs with spaces
# 2014 MAR 21   eqmod.R            Changed names of functions, nested pHSpec and HBalErr in eqSpec, nested derivative functions inside model functions, renamed model functions
# 2014 MAY 30   eqmod.R            Added totk to output, with separate CO2 and other TIC species

# Required package
library(deSolve)

# Calculates speciation for a given set of total concentrations, with equilibrium control of CO2 hydration
eqSpec <- function(tot = tot, temp.c = 20, of = 'a', ll = -14, ul = 0, pH, adjpH = 'H2CO3') {

  # Check arguments
  if(length(tot) != 7 || any(!names(tot)%in%c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'))) stop('Check tot argument. Should be a numeric vector with elements ', "'H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'")

  # Define other functions
  pHSpec <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
  # Iteratively solve, updating ionic strength each time
    i <- sum(0.5*tot[tot>0]) # Very crude guess
    b <- 999
    k <- 10^l.k
    l.a <- 0*l.k # Just for length and names
    l.a['H.'] <- l.a.h
    a <- 10^l.a
    j <- 0  
    di <- 99
    while (di/i>1E-4){ #abs(log10(i/b))>log10(1.001)) {
      j <- j + 1
      b <- i
      l.g<- -a.dh*z^2*sqrt(i)/(1+b.dh*a.par[1, ]*sqrt(i)) + a.par[2, ]*i
      g <- 10^l.g
  
      l.a['K.']   <-log10(tot['K.']) + l.g['K.']
      l.a['Na.']  <-log10(tot['Na.']) + l.g['Na.']
      l.a['Cl.']  <-log10(tot['Cl.']) + l.g['Cl.']
      l.a['NH4.'] <-log10(tot['NH3']*k['NH4.']*g['NH3']*a['H.']*g['NH4.']/(g['NH4.'] + k['NH4.']*g['NH3']*a['H.']) )
      l.a['NH3']  <-l.a['NH4.'] - l.a['H.'] - l.k['NH4.']
      l.a['H2CO3'] <- log10(tot['H2CO3']*a['H.']*g['H2CO3']/(k['HCO3.']*g['H2CO3']/g['HCO3.'] + a['H.'] + k['CO3.2']*g['H2CO3']/(g['CO3.2']*a['H.']) + a['H.']*k['CO2']*g['H2CO3']/g['CO2'] ) )
      l.a['HCO3.'] <- l.k['HCO3.'] + l.a['H2CO3'] - l.a['H.']
      l.a['CO3.2'] <- l.k['CO3.2'] + l.a['H2CO3'] - 2*l.a['H.']
      l.a['CO2']  <-l.k['CO2'] + l.a['H2CO3']
      l.a['OH.']  <- 0 -l.a['H.'] + l.k['OH.']
      l.a['Ac.']  <- log10(k['Ac.']*g['HAc']*tot['HAc']*g['Ac.']/(a['H.']*g['Ac.'] + k['Ac.']*g['HAc']))
      l.a['HAc']  <- l.a['Ac.'] + l.a['H.'] - l.k['Ac.']
  
      l.m <- l.a - l.g
      m <- 10^l.m
      i <- sum(0.5*m*z^2)
      di <- abs(i - b)
    }
    a <- 10^l.a
    cb <- sum(z*m)
	totk <- c(tot[1:2], m[3:4], tot[4:7])
	totk[3] <- totk[3] + m['HCO3.'] + m['CO3.2']
    list(m = m, a = a, g = g, i = i, l.m = l.m, l.a = l.a, l.g = l.g, tot = tot, totk = totk, cb = cb, i.its = j)
  }

  # Calculates error in total H for a given pH guess
  HBalErr <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
    m <- pHSpec(l.a.h = l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)$m
    abs(tot['H.'] - sum(s[, 1]*m)) # All calculated totals given by t(s)%*%m. Alternative would be charge balance.
  }

  # Now code for speciation
  if(!exists('n.calls.spec')) n.calls.spec <<- 0
  n.calls.spec <<- n.calls.spec + 1
  # Temperature
  temp.k <- 273.15+temp.c

  # Henry's law constants
  kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
  kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
 
  # Equilibrium constants
  temp.k <- 273.15+temp.c
  l.k <- c('H.' = 0, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 2.778151, 'K.' = 0, 'Na.' = 0, 'Cl.' = 0, 'HAc' = 0, 
      'OH.'= -4.2195 -2915.16/temp.k, 
      'NH4.'= 0.0905 + 2729.31/temp.k, 
      'HCO3.'= -353.5305 -0.06092*temp.k + 21834.37/temp.k + 126.8339*log10(temp.k) -1684915/temp.k^2, 
      'CO3.2'= -461.4176 -0.093448*temp.k + 26986.16/temp.k + 165.7595*log10(temp.k) -2248629/temp.k^2, 
      'Ac.' = -4.8288 + 21.42/temp.k)
  # Matrix of stoichiometric coefficients. Not used for much, but is good for keeping track of reactions
  s <- matrix(c(1, 0,0, 0,0, 0,0, 0,-1, 1,-1, -2, -1,  0, 1,0, 0,0, 0,0, 0,0, 1,0, 0,0,  0, 0,1, 1,0, 0,0, 0,0, 0,1, 1,0, 
              0, 0,0, 0,1, 0,0, 0,0, 0,0, 0,0,      0, 0,0, 0,0, 1,0, 0,0, 0,0, 0,0,  0, 0,0, 0,0, 0,1, 0,0, 0,0, 0,0,  
	      0, 0,0, 0,0, 0,0, 1,0, 0,0, 0,1), 
            nrow = 13, dimnames = list(c("H.", "NH3", 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', "OH.", "NH4.", 'HCO3.', 'CO3.2', 'Ac.'), 
                                  c("H.", 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc'))
  )
  # Species charge
  z <- c('H.' = 1, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 0, 'K.' = 1, 'Na.' = 1, 'Cl.' = -1, 'HAc' = 0, 'OH.' = -1, 'NH4.' = 1, 'HCO3.' = -1, 'CO3.2' = -2, 'Ac.' = -1)
  # Parameters for calculation of activity coefficients, using extended Debye-Huckel equation (see PHREEQC manual)
  a.par <- matrix(c(9, 0,0, 0.1, 0,0.1, 0,0.1, 3,0, 4.25, 0,3, 0,0, 0.1, 3.5, 0,2.5, 0,5.4, 0,4.5, 0,4.5, 0), nrow = 2, dimnames = list(c('a', 'b'), c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))

  # Calculate a.dh and b.dh parameters for Debye-Huckel equation. Dielectric constant, density rearranged from PHREEQC code
  # Dielectric constant
  de <- 2727.586 + 0.6224107*temp.k - 1075.112*log10(temp.k) - 52000.87/temp.k
  # Water density
  c.d <- 647.26 - temp.k
  d.H2O <- (1 + .1342489*c.d^(1/3) - 3.946263E-3*c.d)/(3.1975 - 0.3151548*c.d^(1/3) - 1.203374E-3*c.d + 7.48908E-13*c.d^4)
  # a.dh and b.dh are A and B in Debye-Huckel equation. Equations are from Tuesdell & Jones (1974)
  a.dh <- 1.82483E6*d.H2O^0.5/(de*temp.k)^(3/2)
  b.dh <- 50.2916*d.H2O^0.5/(de*temp.k)^0.5

  # If pH not specified (assumed to be typical case) it is calculated, if specified, KOH and HAc is added to reach specified value, based on proton balance
  if(missing(pH)) {
    sol <- optimize(f = HBalErr, interval = c(ll, ul), a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z, tol = 1E-15)
    if(sol$objective>5E-7) {
      print(sol)
      stop('Around line 160, optimize didn\'t converge: complete results above, specified limits: ', ll, ' ', ul)
    }
    l.a.h <- sol$minimum
  } else { # pH is specified, acid or base adjusted to match it
    l.a.h<- -pH
    dhb <- 999
    hb <- tot['H.'] - sum(s[, 1]*pHSpec(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z)$m)
    while(dhb>1E-10) { # Must be solved iteratively because i will change
      if(adjpH == 'HAc') {
        tot['HAc'] <- tot['HAc'] - hb  
      } else if(adjpH == 'HCl') {
        tot['Cl.'] <- tot['Cl.'] - hb
        tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'H2CO3') {
        tot['H2CO3'] <- tot['H2CO3'] - hb
      } else if(adjpH == 'KOH') {
          tot['K.'] <- tot['K.'] + hb 
          tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'NaOH') {
          tot['Na.'] <- tot['Na.'] + hb 
          tot['H.'] <- tot['H.'] - hb
      } else if(adjpH == 'NH3') {
        tot['NH3'] <- tot['NH3'] + hb
       } else stop('adjpH argument must be "HAc", "H2CO3", "KOH", or "HCl" but is ', adjpH)
      hb2 <- tot['H.'] - sum(s[, 1]*pHSpec(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z)$m)
      dhb <- abs(hb2-hb)
      hb <- hb2
    }
  }

  out <- pHSpec(l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)
  if(abs(out$cb)>5E-8) warning('Charge balance off in eqSpec. Check results. cb =', out$cb)
  tot.g <- t(s)%*%out$m
  p.CO2 <- out$a['CO2']/kh.CO2
  p.NH3 <- out$a['NH3']/kh.NH3
  names(p.CO2) <- NULL
  if(of == 'a') return(out$a)
  if(of == 'm') return(out$m)
  if(of == 'k') return(l.k)
  if(of == 'all') return(c(out, p.CO2 = p.CO2, p.NH3 = p.NH3))
}


### To test model spec function
#tt <- c("H." = 0.57, 'NH3' = 0.57, 'H2CO3' = 0, 'K.' = 0, 'Na.' = 0, 'Cl.' = 0.57, 'HAc' = 0)
#out <- eqSpec(tot = tt, temp.c = 20, of = 'all', ll = -14, ul = 0)
#-log10(out$a)
#out$g
#tt <- c("H." = -0.101, 'NH3' = 0.1, 'H2CO3' = 0.1, 'K.' = 0.1, 'Na.' = 0.1, 'Cl.' = 0.1, 'HAc' = 0.1)
#out <- eqSpec(tot = tt, temp.c = 20, of = 'all', ll = -14, ul = 0)
#out$m
#out$l.a

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Interface function for emission without diffusive transport
# Emission for single solution--no diffusive transport
eqEmisMod <- function(
  tot,  # Vector of master species totals (mol/kgw), e.g., c("H." = 0, 'NH3' = 0.1, 'H2CO3' = 0.1, 'K.' = 0., 'Na.' = 0., 'Cl.' = 0., 'HAc' = 0.). Must be in this order.
  times,  # Vector of times, in seconds
  h.m,  # Mass transfer coefficient for NH3, m/s
  p.CO2.a,  # Ambient CO2 partial pressure (atm)
  temp.c,  # Temperature in °C
  thk,  # Layer thickness in m
  of = 'all'
  ) {

  # Define function for derivatives
  derCalc <- function(t, y,parms) {
  # Calculate rate constants
    R <- 0.000082057  # Universal gas constant (atm m3/mol-K)
    h.m <- parms$h.m
    temp.c <- parms$temp.c
    p.CO2.a <- parms$p.CO2.a
    thk <- parms$thk
    temp.k <- temp.c + 273.15
    k1 <- 10^(12.066 - 4018.0925/temp.k)
    kr1 <- 10^(14.844 - 4018.0925/temp.k)
    k2 <- 10^(14.354 - 3087.545/temp.k)
    kr2 <- 10^(14.881 - 5524.2258/temp.k)
  
    kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
    h.CO2<- R*temp.k*kh.CO2  # Henry's law constant aq:g (m3/kgw)
    kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
    h.NH3<- R*temp.k*kh.NH3
  
  # Calculate species concentrations
    y[y<0 & names(y) != 'H.'] <- 0
    #y[y<0] <- 0
    out <- eqSpec(tot = y, temp.c = temp.c, of = 'all')
    a <- out$a
    m <- out$m
    cb <- out$cb
  
  # Surface fluxes, mol/m2-s
    j.CO2<- h.m['CO2']*(a['CO2']/h.CO2 - p.CO2.a/(R*temp.k))
    j.NH3<- h.m['NH3']*(a['NH3']/h.NH3 - 0.0/(R*temp.k))
  
  # Derivatives
    dH2CO3.dt  <- -j.CO2/(thk*1000)
    dNH3.dt<- -j.NH3/(thk*1000)
  
    names(j.CO2) <- names(j.NH3) <- NULL
  
    list(tot = c(dH.dt = 0, dNH3.dt = dNH3.dt, dH2CO3.dt = dH2CO3.dt, dK.dt = 0, dNa.dt = 0, dCl.dt = 0, dHAc.dt = 0), act = a, mol = m, cb = cb, 
    pH = -log10(a['H.']), j.NH3 = j.NH3, j.CO2 = j.CO2)
  }

  # Model calculations
  h.m[2] <- h.m[1]*sqrt(17.0305/44.0095);names(h.m) <- c('NH3', 'CO2') # Calculate CO2 h.m based on that for NH3, based on Eq. 6 in Lee et al 2004
  out <- lsoda(y = tot, times = times, func = derCalc, parms = list(h.m = h.m, temp.c = temp.c, p.CO2.a = p.CO2.a, thk = thk))
  tot.out <- matrix(as.vector(out[, colnames(out) %in% c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')]), nrow = dim(out)[1], ncol = length(tot), dimnames = list(t = out[, 'time'], msp = c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')))
  act.out <- matrix(as.vector(out[, grep('act', colnames(out))]), nrow = dim(out)[1], ncol = 13, dimnames = list(t = out[, 'time'], sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  mol.out <- matrix(as.vector(out[, grep('mol', colnames(out))]), nrow = dim(out)[1], ncol = 13, dimnames = list(t = out[, 'time'], sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  cb.out <- matrix(as.vector(out[, grep('cb', colnames(out))]), nrow = dim(out)[1], ncol = 1, dimnames = list(t = out[, 'time'], cb = 'cb'))
  ph.out <- matrix(as.vector(out[, grep('pH', colnames(out))]), nrow = dim(out)[1], ncol = 1, dimnames = list(t = out[, 'time'], 'pH'))
  emis.out <- data.frame(t = out[, 'time'], j.NH3 = out[, 'j.NH3'], j.CO2 = out[, 'j.CO2'], e.NH3 = 1000*(tot['NH3'] - tot.out[, 'NH3'])*thk, 
            e.CO2 = 1000*(tot['H2CO3'] - tot.out[, 'H2CO3'])*thk)
  if(of == 'all') list(times = times, tot = tot.out, act = act.out, mol = mol.out, cb = cb.out, ph = ph.out, emis = emis.out, pars = c(h.m = h.m, p.CO2.a = p.CO2.a, temp.c = temp.c, thk = thk))
}

## To test kin.mod
#eqEmisMod(tot = c("H." = -0.101, 'NH3' = 0.1, 'H2CO3' = 0.1, 'K.' = 0.1, 'Na.' = 0.1, 'Cl.' = 0.1, 'HAc' = 0.1), times = c(0, 1,3600), h.m = 0.001, temp.c = 20, p.CO2.a = 400E-6, thk = 0.01)

# Interface function for complete model with transport
eqEmisDiffMod <- function(
  c.thk,          # Vector of cell thicknesses (m) or single value
  CO2.emis = TRUE,  # Include CO2 emission?
  D = 1E-9,         # Diffusion coefficient (m2/s)
  h.m,            # Mass transfer coefficient for NH3, m/s
  lb = 'nf',        # Lower boundary condition, no flow by default, set to 'cc' for constant concentration
  n.cells = length(c.thk),  # Total number of cells in simulation
  of = 'all',       # Output format
  p.CO2.a,        # Ambient CO2 partial pressure (atm)
  times,          # Vector of times (s)
  temp.c,         # Temperature (C)
  thk = 0,          # Total thickness (m), only used for lb = 'mx'
  tot,            # Vector of master species totals (mol/kgw), e.g., c("H." = 0, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0., 'Na.' = 0., 'Cl.' = 0., 'HAc' = 0.). Must be in this order.
  ph.results = FALSE # Set to true to calculate and return all pH values, otherwise just surface
  ) {

  if(!(lb %in% c('nf', 'cc', 'mx'))) stop('lb error')

  # Calculate derivatives
  
  # Calculate derivatives of emission + transport equations
  derCalc <- function(t, y,parms) {
    n.calls.der <<- n.calls.der + 1
    R <- 0.000082057  # Universal gas constant (atm m3/mol-K)
    c.thk <- parms$c.thk
    D <- parms$D
    dx <- parms$dx
    h.m <- parms$h.m
    lb <- parms$lb
    n.cells <- parms$n.cells
    p.CO2.a <- parms$p.CO2.a
    temp.c <- parms$temp.c
    temp.k <- temp.c + 273.15
    tot.lb <- parms$tot.lb
  #write.table(t, '14.dmp', row.names = F, col.names = F, append = T)
  # Put total concentrations in matrix for working with them
    y[y<0 & names(y) != 'H.'] <- 0 # CRUDE
    #y[y<0] <- 0 # CRUDE
    tot <- matrix(y, nrow = n.cells, dimnames = list(1:n.cells, c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')))
  
  # Calculate constants
    kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
    h.CO2<- R*temp.k*kh.CO2  # Henry's law constant aq:g (m3/kgw)
  
    kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
    h.NH3<- R*temp.k*kh.NH3
  
  # Diffusion/dispersion. Note single D for all solutes, and that totals diffuse, not species
    j.diff.out <- D*1000*rbind(c(0, 0,0, 0,0, 0,0), diff(tot))/dx # Diffusion out, toward cell 1 (g/m2-s). Note 1000 kgw/m3 conversion to make concentrations volumetric
    j.diff.lb <- as.vector(D*1000*diff(rbind(tot[n.cells, ],tot.lb))/(0.5*c.thk[n.cells])) # Diffusion into bottom cell. Used below only if lb = "c"
    names(j.diff.lb) <- c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')
  
  # Calculate species activities, only for upper-most layer
    ll<- -13  
    ul<- -4
    a <- eqSpec(tot = tot[1, ],temp.c = temp.c, of = 'a', ll = ll, ul = ul)
  
  # Derivatives
    # Surface fluxes, mol/m2-s
    j.CO2 <- h.m['CO2']*(a['CO2']/h.CO2 - p.CO2.a/(R*temp.k))
    j.NH3 <- h.m['NH3']*(a['NH3']/h.NH3 - 0.0/(R*temp.k))
  
    # Derivatives for concentration
    if(lb != 'cc') j.diff.lb <- 0*j.diff.lb
    dH2CO3.dt<- (-j.CO2*c(1, rep(0, n.cells-1)) + diff(c(j.diff.out[, 'H2CO3'], j.diff.lb['H2CO3'])) )/(c.thk*1000)
    dNH3.dt  <- (-j.NH3*c(1, rep(0, n.cells-1)) + diff(c(j.diff.out[, 'NH3'], j.diff.lb['NH3'])) )/(c.thk*1000)
  
    names(j.CO2) <- names(j.NH3) <- NULL
  
    list(tot = c(dH.dt = rep(0, n.cells), dNH3.dt, dH2CO3.dt, dK.dt = rep(0, n.cells), dNa.dt = rep(0, n.cells), dCl.dt = rep(0, n.cells), dHAc.dt = rep(0, n.cells)), act = a, pH = -log10(a['H.']), j.NH3 = j.NH3, j.CO2 = j.CO2, p.CO2.s = a['CO2']/h.CO2*R*temp.k, p.NH3.s = a['NH3']/h.NH3*R*temp.k)
  }

  if(lb == 'mx') {
    n.cells <- n.cells + 1 # For bottom, tracked cell
    c.thk <- c(c.thk, thk - sum(c.thk))
  }

  if (any(c.thk <= 0)) stop('c.thk error Auq1995')

  # Initial total masses, for calculating relative emission
  # 1000 goes from kg to m3
  imass.NH3 <- sum(1000*tot['NH3']*c.thk)
  imass.CO2 <- sum(1000*(tot['CO2'] + tot['H2CO3'])*c.thk) 
 
  time.start <- Sys.time()
  n.calls.der <<- 0
  n.calls.spec <<- 0

  tot.lb <- tot # Lower boundary constant concentrations
  if(length(c.thk) == 1) c.thk <- rep(c.thk, n.cells)
  if(length(tot) == 7) tot <- rep(tot, each = n.cells)
  pos <- cumsum(c.thk) - 0.5*c.thk  # position (m), based on cell centers
  dx <- c(c.thk[n.cells], diff(pos))
  h.m[2] <- ifelse(CO2.emis, h.m[1]*sqrt(17.0305/44.0095), 0);names(h.m) <- c('NH3', 'CO2') # Calculate CO2 h.m based on that for NH3, based on Eq. 6 in Lee et al 2004

  if(lb == 'mx') dx[n.cells] <- 0.5*c.thk[n.cells-1] # Eliminates resistance within bottom cell

  # Call up ODE solver. Note that only NH3 and IC solutes change over time
  out <- ode.1D(y = tot, times = times, func = derCalc, 
        parms = list(h.m = h.m, n.cells = n.cells, p.CO2.a = p.CO2.a, temp.c = temp.c, c.thk = c.thk, dx = dx, D = D, lb = lb, p.CO2.a = p.CO2.a, tot.lb = tot.lb), 
        nspec = 7, dimens = n.cells) #, method = 'lsodes') #, hini = 0.01)

  tot.out <- array(as.vector(out[, colnames(out) %in% c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')]), dim = c(dim(out)[1], n.cells, 7), dimnames = list(t = signif(times, 5), pos = pos, msp = c('H.', 'NH3', 'H2CO3', 'K.', 'Na.', 'Cl.', 'HAc')))
  act.out <- matrix(as.vector(out[, grep('act', colnames(out))]), nrow = dim(out)[1], ncol = 13, dimnames = list(t = signif(times, 5), sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  if(ph.results) {
    ph.out <- matrix(NA, nrow = length(times), ncol = length(pos), dimnames = list(time = signif(times, 5), pos = pos))
    tot.out[, ,-1][tot.out[, ,-1]<0] <- 0
    for(tt in 1:length(times)) {
      for(pp in 1:length(pos)) {
        ph.out[tt, pp]<- -as.numeric(eqSpec(tot = tot.out[tt, pp, ],temp.c = temp.c, of = 'all')$l.a['H.'])
      }
    }
  } else ph.out <- as.vector(out[, grep('pH', colnames(out))]) 
  names(ph.out) <- out[, 'time']
  emis.NH3<- sum(1000*tot['NH3']*c.thk) - apply(1000*tot.out[, ,'NH3'], 1,function(x) sum(x*c.thk)) 
  emis.CO2<- sum(1000*(tot['H2CO3'])*c.thk) - apply(1000*(tot.out[, ,'H2CO3']), 1,function(x) sum(x*c.thk)) 
  emis.NH3.rel<- emis.NH3/imass.NH3
  emis.CO2.rel<- emis.CO2/imass.CO2
  p.CO2.s.out <- out[, grep('p.CO2.s', colnames(out))]
  p.NH3.s.out <- out[, grep('p.NH3.s', colnames(out))]
  emis.out <- data.frame(t = signif(out[, 'time'], 5), j.NH3 = out[, 'j.NH3'], j.CO2 = out[, 'j.CO2'], e.NH3 = emis.NH3, e.CO2 = emis.CO2, e.NH3.rel = emis.NH3.rel, e.CO2.rel = emis.CO2.rel, p.CO2.s = p.CO2.s.out, p.NH3.s = p.NH3.s.out)
  pars <- c(D = D, h.m = h.m, kin.mult = NA, lb = lb, n.cells = n.cells, p.CO2.a = p.CO2.a, thk = sum(c.thk), temp.c = temp.c)
  if (sum(times %in% c('0', '3600', '86400')) == 3) {
    summ.0 <- emis.out[times == 0, -1];names(summ.0) <- paste(names(summ.0), '.0', sep = '')
    summ.1h <- emis.out[times == 3600, -1];names(summ.1h) <- paste(names(summ.1h), '.1h', sep = '')
    summ.1d <- emis.out[times == 86400, -1];names(summ.1d) <- paste(names(summ.1d), '.1d', sep = '')
    ph.summ <- ph.out[times%in%c('0', '3600', '86400')];names(ph.summ) <- c('ph.0', 'ph.1h', 'ph.1d')
  } else summ.0 <- summ.1h <- summ.1d <- ph.summ <- NA
  summ <- unlist(c(summ.0, summ.1h, summ.1d, ph.summ))

  exe.time <- Sys.time() - time.start
  if(of == 'all') list(times = times, pos = pos, tot = tot.out, act = act.out, ph = ph.out, emis = emis.out, exe.time = exe.time, 
  n.calls = c(der = n.calls.der, spec = n.calls.spec), pars = pars, summ = summ)
}

# Kinetic functions
# Model for NH3 and CO2 speciation and transport, including kinetically-limited reactions for H2CO3/CO2
# Author: Sasha Hafner
# History of revisions
# Date         File                   Description
# 2010 OCT 08  SpecMod_v1.R           First version that really works. Cacluates speciation of NH3 and CO2.
# 2010 OCT 11  SpecMod_v2.R           Trying to incorporate temperature response of equilibrium constants, 
#                                     kinetics, and emission. Seems to work, but is slow.
# 2010 OCT 12  SpecMod_v3.R           Switching to H2CO3 as master C species.
# 2010 OCT 13  SpecMod_v4.R           Trying to speed things up. Kinetic component now works pretty well, but is only single layer.
# 2010 OCT 13  SpecMod_v5.R           Includes transport. Now uses extended Debye-Huckel equation for activity coefficients
# 2010 OCT 14  SpecMod_v6.R           Modifying so CO2 is returned by spec function. Everything seems to work.
# 2010 OCT 18  SpecMod_v7.R           Adding a function for running the complete transport model (was manual)
#                                     Corrected a 1000-fold error in derivative calculations
# 2010 OCT 19  SpecMod_v8.R           Trying to speed up execution
# 2010 OCT 20  SpecMod_v8.R           Corrected an error in dx
# 2010 OCT 27  SpecMod_v12.R          Between this and v8, tried several things for improving speed--all failed.
#                                     This version uses a slightly simplified version of spec.ph for opitimization
# 2010 OCT 28  SpecMod_v13.R          Uses global variables to track tot and a. Only calculates speciation when needed.
#                                     Slightly faster than v12, but results are sensitive to the change in tot that is 
#                                     considered insignificant. If it is too small, ode1D makes many more calls to der.calc, 
#                                     and results in slow runs. Saw this by writing out time at each der.calc call to a file.
# 2010 NOV 03  SpecMod_v14.R          Adding K+, Na+, Cl-, and acetic acid.
# 2010 NOV 04  SpecMod_v15.R          Trying to pull out non-volatile solutes from ODE solver to speed things up
# 2010 NOV 08-09SpecMod_v16.R         Tried to write my own solver that uses simple fixed relative change approach. 
#                                     As feared, it was very slow. Incorporated a new approach where CO2(aq) in the upper layer
#                                     is set to zero. Speeds things up slightly, and is optional.
# 2010 NOV 10  SpecKinTransMod_v1.R   After version 16, I made an equilibrium-only version (SpecTransMod_v1.R).
# 2010 NOV 16  SpecKinTransMod_v1.R   Adding more detail to output.
# 2010 NOV 17  SpecKinTransMod_v2.R   Adding two-film-type option, with constant composition layer below bottom cell.
#                                     I changed the diffusion rate calculations slightly to do this.
# 2010 NOV 19  SpecKinTransMod_v2.R   Added multiplier for kinetic rates
# 2010 DEC 24  SpecKinTransMod_v2.R   Added effect of ambient CO2
# 2010 DEC 28  SpecKinTransMod_v2.R   Added line to deal with n.calls.spec when only speciation is called up without transport
# 2010 DEC 29  SpecKinTransMod_v3.R   Trying to improve code for output a bit, also adding some output for a concise summary
# 2010 JAN 04  SpecKinTransMod_v4.R   Improving the der.calc code a bit.
# 2010 JAN 20  SpecKinTransMod_v4.R   Corrected small error.
# 2011 JAN 21  SpecKinTransMod_v5.R   Added lower boundary option--essentially two-film with finite bulk layer, specified by lb = 'mx'.
# 2011 APR 05  SpecKinTransMod_v6.R   Making h.m length of 2, for NH3 and CO2, and calculate h.m for CO2 from h.m for NH3.
# 2011 APR 07  SpecKinTransMod_v7.R   Added kin.mult option to kin.mod.
# 2011 APR 08  SpecKinTransMod_v7.R   Added if statement so summ is in output only if relevant times are requested.
# 2011 APR 14  SpecKinTransMod_v8.R   Changing equations for equilibrium constants and kinetic constants to latest version, 
#                                     based on Plummer and Bunsenberg (1982).
# 2011 APR 14  SpecKinTransMod_v9.R   Same as 8 (made a change and changed back)
# 2011 SEP 26  SpecKinTransMod_v13.R  Added CO2.emis argument, which can be set to FALSE for turning off CO2 emission.
#                                     This change doesn't make much sense here, since an equilibrium model can be used if
#                                     there is no CO2 emission.
#                                     Note that changes made for simulating Ni and jar trials are skipped--see version
#                                     12 for these.
# 2011 SEP 27  SpecKinTransMod_v13.R  Replaced dielectric constant calculation and Debye-Huckel A and B with equations from PHREEQC and WATEQ.
# 2012 JAN 10  SpecKinTransMod_v14.R  Simplified the equation for calculating l.a['H2CO3'] to match the one in the equations document (just moved g['H2CO3']
#                                     around). Model simulations in paper were run using version 13, but the two should be identical.
# 2013 NOV 26  kinmod.R               Changed file name
# 2014 MAR 21  kinmod.R               Changing function names and nesting small functions within others
# 2014 OCT 31  kinmod.R               Removed degree sign from comments, was causing a warning in grepl() which was called when source() is called

# Required packages
library(deSolve)

# Calculates speciation for a given set of total concentrations, with kinetic control of CO2 hydration
kinSpec <- function(tot, temp.c = 20, of = 'a', ll = -14, ul = 0) {
  if(!exists('n.calls.spec')) n.calls.spec <<- 0
  n.calls.spec <<- n.calls.spec + 1

  # Define functions
  # Calculates speciation for a given pH, assuming kinetic control of CO2 hydration/dehydration
  pHSpec <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
  # Iteratively solve, updating ionic strength each time
    i <- sum(0.5*tot[tot>0]) # Very crude guess
    b <- 999
    k <- 10^l.k
    l.a <- 0*l.k # Just for length and names
    l.a['H.'] <- l.a.h
    a <- 10^l.a
    j <- 0  
    di <- 99
    while (di/i>1E-4){ #abs(log10(i/b))>log10(1.001)) {
      j <- j + 1
      b <- i
      l.g<- -a.dh*z^2*sqrt(i)/(1+b.dh*a.par[1, ]*sqrt(i)) + a.par[2, ]*i
      g <- 10^l.g
  
      l.a['CO2']  <-log10(tot['CO2']) + l.g['CO2']
      l.a['K.']  <-log10(tot['K.']) + l.g['K.']
      l.a['Na.']  <-log10(tot['Na.']) + l.g['Na.']
      l.a['Cl.']  <-log10(tot['Cl.']) + l.g['Cl.']
      l.a['NH4.'] <-log10(tot['NH3']*k['NH4.']*g['NH3']*a['H.']*g['NH4.']/(g['NH4.'] + k['NH4.']*g['NH3']*a['H.']) )
      l.a['NH3']  <-l.a['NH4.'] - l.a['H.'] - l.k['NH4.']
      l.a['H2CO3'] <- log10(tot['H2CO3']*a['H.']/(k['HCO3.']/g['HCO3.'] + a['H.']/g['H2CO3'] + k['CO3.2']/(g['CO3.2']*a['H.'])) )
      l.a['HCO3.'] <- l.k['HCO3.'] + l.a['H2CO3'] - l.a['H.']
      l.a['CO3.2'] <- l.k['CO3.2'] + l.a['H2CO3'] - 2*l.a['H.']
      l.a['OH.']  <- 0 -l.a['H.'] + l.k['OH.']
      l.a['Ac.']  <- log10(k['Ac.']*g['HAc']*tot['HAc']*g['Ac.']/(a['H.']*g['Ac.'] + k['Ac.']*g['HAc']))
      l.a['HAc']  <- l.a['Ac.'] + l.a['H.'] - l.k['Ac.']
  
      l.m <- l.a - l.g
      m <- 10^l.m
      i <- sum(0.5*m*z^2)
      di <- abs(i - b)
    }
    a <- 10^l.a
    cb <- sum(z*m)
    #if(abs(cb)>1E-10) stop('Charge balance error (', cb, '), maybe there is a problem in pHSpec with ul or ll. See lines 707 & 708, or add a browser() call on this line.')
    list(m = m, a = a, g = g, i = i, l.m = l.m, l.a = l.a, l.g = l.g, tot = tot, cb = cb, i.its = j)
  }

  # Calculates error in total H for a given pH guess
  HBalErr <- function(l.a.h, a.dh, b.dh, a.par, l.k, s,tot, z) {
    m <- pHSpec(l.a.h = l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)$m
    abs(tot['H.'] - sum(s[, 1]*m)) # All calculated totals given by t(s)%*%m. Alternative would be charge balance.
  }

  # Temperature
  temp.k <- 273.15+temp.c

  # Henry's law constants
  kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
  kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
 
  # Equilibrium constants
  l.k <- c('H.' = 0, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 0, 'K.' = 0, 'Na.' = 0, 'Cl.' = 0, 'HAc' = 0, 
      'OH.'= -4.2195 -2915.16/temp.k, 
      'NH4.'= 0.0905 + 2729.31/temp.k, 
      'HCO3.'= -353.5305 -0.06092*temp.k + 21834.37/temp.k + 126.8339*log10(temp.k) -1684915/temp.k^2, 
      'CO3.2'= -461.4176 -0.093448*temp.k + 26986.16/temp.k + 165.7595*log10(temp.k) -2248629/temp.k^2, 
      'Ac.' = -4.8288 + 21.42/temp.k)
  # Matrix of stoichiometric coefficients. Not used for much, but is good for keeping track of reactions
  s <- matrix(c(1, 0,0, 0,0, 0,0, 0,-1, 1,-1, -2, -1,  0, 1,0, 0,0, 0,0, 0,0, 1,0, 0,0,  0, 0,1, 0,0, 0,0, 0,0, 0,1, 1,0,  0, 0,0, 1,0, 0,0, 0,0, 0,0, 0,0, 
         0, 0,0, 0,1, 0,0, 0,0, 0,0, 0,0,     0, 0,0, 0,0, 1,0, 0,0, 0,0, 0,0,  0, 0,0, 0,0, 0,1, 0,0, 0,0, 0,0,  0, 0,0, 0,0, 0,0, 1,0, 0,0, 0,1), 
        nrow = 13, dimnames = list(c("H.", "NH3", 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', "OH.", "NH4.", 'HCO3.', 'CO3.2', 'Ac.'), 
                  c("H.", 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')))
  # Species charge
  z <- c('H.' = 1, 'NH3' = 0, 'H2CO3' = 0, 'CO2' = 0, 'K.' = 1, 'Na.' = 1, 'Cl.' = -1, 'HAc' = 0, 'OH.' = -1, 'NH4.' = 1, 'HCO3.' = -1, 'CO3.2' = -2, 'Ac.' = -1)
  # Parameters for calculation of activity coefficients, using extended Debye-Huckel equation (see PHREEQC manual)
  a.par <- matrix(c(9, 0,0, 0.1, 0,0.1, 0,0.1, 3,0, 4.25, 0,3, 0,0, 0.1, 3.5, 0,2.5, 0,5.4, 0,4.5, 0,4.5, 0), nrow = 2, dimnames = list(c('a', 'b'), c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))

  # Calculate a.dh and b.dh parameters for Debye-Huckel equation. Dielectric constant, density rearranged from PHREEQC code
  # Dielectric constant
  de <- 2727.586 + 0.6224107*temp.k - 1075.112*log10(temp.k) - 52000.87/temp.k
  # Water density
  c.d <- 647.26 - temp.k
  d.H2O <- (1 + .1342489*c.d^(1/3) - 3.946263E-3*c.d)/(3.1975 - 0.3151548*c.d^(1/3) - 1.203374E-3*c.d + 7.48908E-13*c.d^4)
  # a.dh and b.dh are A and B in Debye-Huckel equation. Equations are from Truesdell & Jones (1974)
  a.dh <- 1.82483E6*d.H2O^0.5/(de*temp.k)^(3/2)
  b.dh <- 50.2916*d.H2O^0.5/(de*temp.k)^0.5
  
  sol <- optimize(f = HBalErr, interval = c(ll, ul), a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z, tol = 1E-15)
  if(sol$objective>5E-7) {
    print(sol)
    stop('Around line 540, optimize didn\'t converge: complete results above, specified limits: ', ll, ' ', ul)
  }
  l.a.h <- sol$minimum
  out <- pHSpec(l.a.h, a.dh = a.dh, b.dh = b.dh, a.par = a.par, l.k = l.k, s = s, tot = tot, z = z)
  p.CO2 <- out$a['co2']/kh.CO2
  p.NH3 <- out$a['NH3']/kh.NH3
  if(of == 'a') return(out$a)
  if(of == 'm') return(out$m)
  if(of == 'k') return(l.k)
  if(of == 'all') return(c(out, p = c(p.CO2, p.NH3)))
}

## To test model spec function
#n.calls.spec <<- 0
#tt <- c("H." = -0.101, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0.1, 'Na.' = 0.1, 'Cl.' = 0.1, 'HAc' = 0.1)
#out <- kinSpec(tot = tt, temp.c = 25, of = 'all', ll = -14, ul = 0)
#out$m
#out$l.a
#tt <- c("H." = 0, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0, 'Na.' = 0, 'Cl.' = 0, 'HAc' = 0)
#system.time(for(i in 1:100) kinSpec(tot = tt, temp.c = 25, of = 'all'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Interface function for kinetics without transport
kinEmisMod <- function(
  tot,  # Vector of master species totals (mol/kgw), e.g., c("H." = 0, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0., 'Na.' = 0., 'Cl.' = 0., 'HAc' = 0.). Must be in this order.
  times,  # Vector of times, in seconds
  h.m,  # Mass transfer coefficient for NH3, m/s
  kin.mult = 1,        # Multiplier for kinetic rates
  p.CO2.a = 400E-6,  # Ambient CO2 partial pressure (atm)
  temp.c = 20,  # Temperature in C
  thk = 0.01,  # Layer thickness in m
  of = 'all'
  ) {

  # Check tot argument
  if(length(tot) != 8 || any(!names(tot)%in%c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc'))) warning('Is tot argument correct? Must be in the order', "'H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc'")

  # Calculate derivatives
  # Kinetics and emission for single solution--no diffusive transport
  derCalc <- function(t, y,parms) {
  # Calculate rate constants
    R <- 0.000082057  # Universal gas constant (atm m3/mol-K)
    h.m <- parms$h.m
    kin.mult <- parms$kin.mult
    temp.c <- parms$temp.c
    thk <- parms$thk
    p.CO2.a <- parms$p.CO2.a
    temp.k <- temp.c + 273.15
    k1 <-kin.mult*10^(12.0657 - 4018.09/temp.k)
    kr1 <- kin.mult*10^(14.8438 - 4018.09/temp.k)
    k2 <-kin.mult*10^(-337.2083 + -0.06092*temp.k + 19225.3/temp.k + 126.8339*log10(temp.k) + -1684915/temp.k^2)
    kr2 <- kin.mult*10^(14.8809 + -5524.23/temp.k)
  
    kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
    h.CO2<- R*temp.k*kh.CO2  # Henry's law constant aq:g (m3/kgw)
    kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
    h.NH3<- R*temp.k*kh.NH3
  
  # Calculate species concentrations
    y[y<0 & names(y) != 'H.'] <- 0
    #y[y<0] <- 0
    out <- kinSpec(tot = y, temp.c = temp.c, of = 'all')
    a <- out$a
    m <- out$m
    cb <- out$cb
    if(any(abs(cb)>5E-8)) stop('Charge balance error, maybe there is a problem in pHSpec with ul or ll.') 
  
  # Surface fluxes, mol/m2-s
    j.CO2<- h.m['CO2']*(a['CO2']/h.CO2 - p.CO2.a/(R*temp.k))
    j.NH3<- h.m['NH3']*(a['NH3']/h.NH3 - 0.0/(R*temp.k))
  
  # Derivatives
    dCO2.dt  <- -j.CO2/(thk*1000) - k1*a['CO2']*1.0 - k2*a['CO2']*a['OH.'] + kr1*a['H2CO3'] + kr2*a['HCO3.']
    dH2CO3.dt<-           k1*a['CO2']*1.0 + k2*a['CO2']*a['OH.'] - kr1*a['H2CO3'] - kr2*a['HCO3.']
    dNH3.dt<- -j.NH3/(thk*1000)
  #tt <- c("H." = -0.101, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0.1, 'Na.' = 0.1, 'Cl.' = 0.1, 'HAc' = 0.1)
    names(j.CO2) <- names(j.NH3) <- NULL
  
    list(tot = c(dH.dt = 0, dNH3.dt = dNH3.dt, dH2CO3.dt = dH2CO3.dt, dCO2.dt = dCO2.dt, dK.dt = 0, dNa.dt = 0, dCl.dt = 0, dHAc.dt = 0), 
      act = a, mol = m, cb = cb, pH = -log10(a['H.']), j.NH3 = j.NH3, j.CO2 = j.CO2, p.CO2.s = a['CO2']/h.CO2*R*temp.k)
  }
  
  time.start <- Sys.time()
  n.calls.der <<- 0
  n.calls.spec <<- 0

  h.m[2] <- h.m[1]*sqrt(17.0305/44.0095);names(h.m) <- c('NH3', 'CO2') # Calculate CO2 h.m based on that for NH3, based on Eq. 6 in Lee et al 2004
  out <- lsoda(y = tot, times = times, func = derCalc, parms = list(h.m = h.m, kin.mult = kin.mult, p.CO2.a = p.CO2.a, temp.c = temp.c, thk = thk))
  tot.out <- matrix(as.vector(out[, colnames(out) %in% c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')]), nrow = dim(out)[1], ncol = 8, dimnames = list(t = out[, 'time'], msp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')))
  act.out <- matrix(as.vector(out[, grep('act', colnames(out))]), nrow = dim(out)[1], ncol = 13, dimnames = list(t = out[, 'time'], sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  mol.out <- matrix(as.vector(out[, grep('mol', colnames(out))]), nrow = dim(out)[1], ncol = 13, dimnames = list(t = out[, 'time'], sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  cb.out <- matrix(as.vector(out[, grep('cb', colnames(out))]), nrow = dim(out)[1], ncol = 1, dimnames = list(t = out[, 'time'], cb = 'cb'))
  ph.out <- matrix(as.vector(out[, grep('pH', colnames(out))]), nrow = dim(out)[1], ncol = 1, dimnames = list(t = out[, 'time'], 'pH'))
  emis.out <- data.frame(t = out[, 'time'], j.NH3 = out[, 'j.NH3'], j.CO2 = out[, 'j.CO2'], e.NH3 = 1000*(tot['NH3'] - tot.out[, 'NH3'])*thk, 
            e.CO2 = 1000*(tot['H2CO3'] + tot['CO2'] - tot.out[, 'H2CO3'] - tot.out[, 'CO2'])*thk)
  emis.out$e.NH3.rel<- emis.out$e.NH3/(1000*tot['NH3']*thk)
  emis.out$e.CO2.rel<- emis.out$e.CO2/(1000*(tot['H2CO3']+tot['CO2'])*thk)
  if(of == 'all') list(times = times, tot = tot.out, act = act.out, mol = mol.out, cb = cb.out, ph = ph.out, emis = emis.out, pars = c(h.m = h.m, p.CO2.a = p.CO2.a, temp.c = temp.c, thk = thk))
}

### To test kin.mod
#kinEmisMod(tot = c("H." = -0.101, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0, 'K.' = 0.1, 'Na.' = 0.1, 'Cl.' = 0.1, 'HAc' = 0.1), times = c(0, 1,3600), h.m = 0.0, temp.c = 20, p.CO2.a = 400E-6, thk = 0.01)
#kinEmisMod(tot = c("H." = -0.1, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0., 'Na.' = 0., 'Cl.' = 0., 'HAc' = 0.1), times = c(0:60), h.m = 0.001, temp.c = 25, thk = 0.01)

# Interface function for complete model with transport
kinEmisDiffMod <- function(
  c.thk,          # Vector of cell thicknesses (m) or single value
  CO2.emis = TRUE,  # Include CO2 emission?
  D = 1E-9,         # Diffusion coefficient (m2/s)
  Dd = NULL,        # Deep diffusion coefficient (m2/s)
  zD = NULL,        # Location where diffusion coefficient switches
  e1 = 1,           # Set to 1 to explicitly model CO2 emission from the upper-most layer, 0 to base it on diffusion from below
  h.m,            # Mass transfer coefficient for NH3 (or both NH3 and CO2 if h.m.constant = TRUE), m/s
  h.m.constant = FALSE, # Set to TRUE to use single h.m value for CO2 and NH3
  kin.mult = 1,     # Multiplier for kinetic rates
  lb = 'nf',        # Lower boundary condition, no flow by default, set to 'cc' for constant concentration
  n.cells = length(c.thk),  # Total number of cells in simulation
  of = 'all',       # Output format
  p.CO2.a,        # Ambient CO2 partial pressure (atm)
  times,          # Vector of times (s)
  temp.c,         # Temperature (C)
  thk = 0,          # Total thickness (m), only used for lb = 'mx'
  tot,            # Vector of master species totals (mol/kgw), e.g., c("H." = 0, 'NH3' = 0.1, 'H2CO3' = 0.1, 'CO2' = 0.0, 'K.' = 0., 'Na.' = 0., 'Cl.' = 0., 'HAc' = 0.). Must be in this order.
  parallel = FALSE  # TRUE for parallel processing
  ) {

  if(!(lb %in% c('nf', 'cc', 'mx'))) stop('lb error')
  
  if(parallel) {
    library(foreach)
    library(parallel)
    library(doParallel)
    if(!exists('clstr')) {
      clstr <- makeCluster(detectCores()-1)
      registerDoParallel(clstr)
    }
  }

  # Define functions
  # Calculate derivatives of kinetic + transport equations
  derCalc <- function(t, y,parms) {
    n.calls.der <<- n.calls.der + 1
    R <- 0.000082057  # Universal gas constant (atm m3/mol-K)
    c.thk <- parms$c.thk
    D <- parms$D
    dx <- parms$dx
    e1 <- parms$e1
    h.m <- parms$h.m
    kin.mult <- parms$kin.mult
    n.cells <- parms$n.cells
    p.CO2.a <- parms$p.CO2.a 
    temp.c <- parms$temp.c
    lb <- parms$lb
    tot.lb <- parms$tot.lb
    temp.k <- temp.c + 273.15

    # Put total concentrations in matrix for working with them
    y[y<0] <- 0 # CRUDE
    tot.k <- matrix(y, nrow = n.cells, dimnames = list(1:n.cells, c('NH3', 'H2CO3', 'CO2')))
    tot <- cbind('H.' = tot.glob[, 'H.'], tot.k, tot.glob[, c('K.', 'Na.', 'Cl.', 'HAc')])
  
  # Calculate constants
    k1 <-kin.mult*10^(12.0657 - 4018.09/temp.k)
    kr1 <- kin.mult*10^(14.8438 - 4018.09/temp.k)
    k2 <-kin.mult*10^(-337.2083 + -0.06092*temp.k + 19225.3/temp.k + 126.8339*log10(temp.k) + -1684915/temp.k^2)
    kr2 <- kin.mult*10^(14.8809 + -5524.23/temp.k)
  
    kh.CO2<- 10^(108.38578 +0.01985076*temp.k - 6919.53/temp.k - 40.45154*log10(temp.k) + 669365/temp.k^2)
    h.CO2<- R*temp.k*kh.CO2  # Henry's law constant aq:g (m3/kgw)
  
    kh.NH3<- 10^(-3.51645 -0.0013637*temp.k + 1701.35/temp.k)
    h.NH3<- R*temp.k*kh.NH3
  
  # Diffusion/dispersion. Note single D for all solutes, and that totals diffuse, not species
    j.diff.out <- D*1000*rbind(c(0, 0,0), diff(tot.k))/dx # Diffusion out, toward cell 1 (mol/m2-s). Note 1000 kgw/m3 conversion to make concentrations volumetric
    j.diff.lb <- as.vector(D[1]*1000*diff(rbind(tot.k[n.cells, ],tot.lb))/(0.5*c.thk[n.cells])) # Diffusion into bottom cell. Used below only if lb = "c"
    names(j.diff.lb) <- c('NH3', 'H2CO3', 'CO2')
  
    # Find location of cells with changed tot
    d.tot <- abs(tot.k - tot.glob[, c('NH3', 'H2CO3', 'CO2')])
    ch.cells <- c(1:n.cells)[apply(d.tot, 1,max)>1E-15]
  
    # Update tot.glob, but only for those cells that changed (to prevent significant cumulative error)
    tot.glob[ch.cells, c('NH3', 'H2CO3', 'CO2')] <<- tot.k[ch.cells, ]
  
    # Set a to a.glob for all cells. Those with tot changes are updated below
    a <- a.glob
    m <- m.glob
  
    ll<- -13  
    ul<- -5.0
    
    n.calls.spec <<- n.calls.spec + length(ch.cells)
    if(parallel) {
      a[ch.cells, ] <- foreach(i = ch.cells, .combine = 'rbind', .export = ls(envir = globalenv())) %dopar% {
        tot.x <- as.numeric(tot[i, ])
        names(tot.x) <- c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')
        kinSpec(tot = tot.x, temp.c = temp.c, of = 'a', ll = ll, ul = ul)
      }
    } else {
      for (i in ch.cells) {
        tot.x <- as.numeric(tot[i, ])
        #ll<- log10(a[i, 'H.']) - 0.05
        #ul<- log10(a[i, 'H.']) + 0.05
        names(tot.x) <- c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')
        pred <- kinSpec(tot = tot.x, temp.c = temp.c, of = 'all', ll = ll, ul = ul)
        a[i, ] <- pred$a
        m[i, ] <- pred$m
      }
    }
  
    # Update a.glob
    a.glob[ch.cells, ] <<- a[ch.cells, ]
    m.glob[ch.cells, ] <<- m[ch.cells, ]
  
    # Derivatives
    # Surface fluxes, mol/m2-s
    if(e1 == 1) {
      j.CO2 <- h.m['CO2']*(a[1, 'CO2']/h.CO2 - p.CO2.a/(R*temp.k))
    }
    j.NH3 <- h.m['NH3']*(a[1, 'NH3']/h.NH3 - 0.0/(R*temp.k))
  
    # CO2 hydration, mol/kgw-s
    if(e1 == 0 & h.m['CO2']>0) a[1, 'CO2'] <- 0
    kin.CO2  <- k1*a[, 'CO2']*1.0 + k2*a[, 'CO2']*a[, 'OH.'] - kr1*a[, 'H2CO3'] - kr2*a[, 'HCO3.']
  
    # Derivatives for concentration
    if(lb != 'cc') j.diff.lb <- 0*j.diff.lb
    if(e1 == 1) {
      dCO2.dt <- (-j.CO2*c(1, rep(0, n.cells-1)) + diff(c(j.diff.out[, 'CO2'], j.diff.lb['CO2'])) )/(c.thk*1000) - kin.CO2
    }
    else if(e1 == 0) {
      dCO2.dt  <- diff(c(j.diff.out[, 'CO2'], j.diff.lb['CO2']))/(c.thk*1000) - kin.CO2  # Diffusion and kinetics only
      j.CO2 <- dCO2.dt[1]
      dCO2.dt[1] <- 0
    } else stop('e1 error', e1)
    dH2CO3.dt<- diff(c(j.diff.out[, 'H2CO3'], j.diff.lb['H2CO3']))/(c.thk*1000) + kin.CO2
    dNH3.dt  <- (-j.NH3*c(1, rep(0, n.cells-1)) + diff(c(j.diff.out[, 'NH3'], j.diff.lb['NH3'])) )/(c.thk*1000)
  
    names(j.CO2) <- names(j.NH3) <- names(kin.CO2) <- NULL
  
    list(tot = c(dNH3.dt, dH2CO3.dt, dCO2.dt), act = a, mol = m, pH = -log10(a[, 'H.']), j.NH3 = j.NH3, j.CO2 = j.CO2, kin.CO2 = kin.CO2, p.CO2.s = a[1, 'CO2']/h.CO2*R*temp.k, p.NH3.s = a[1, 'NH3']/h.NH3*R*temp.k) # ,
      #dCO2.dt = dCO2.dt, dH2CO3.dt = dH2CO3.dt, dNH3.dt = dNH3.dt) # Works but slow things down a bit
  }

  if(lb == 'mx') {
    n.cells <- n.cells + 1 # For bottom, tracked cell
    c.thk <- c(c.thk, thk - sum(c.thk))
  }

  # Diffusion coefficient, if it varies with depth
  if(!is.null(Dd) & !is.null(zD)) {
    Ds <- D
    D <- rep(D, n.cells)
    if(zD<sum(c.thk)) {
      di <- which(cumsum(c.thk)>=zD)[1]
      D[(di+1):n.cells] <- Dd
      ffs <- (sum(c.thk[1:di]) - zD)/c.thk[di]
      # If zD doesn't fall exactly on a break between cells, weight it
      D[di] <- 10^(ffs*log10(Dd) + (1 - ffs)*log10(Ds))
    }
  }

  # NTS: Why was this duplicated????
  #if(lb == 'mx') {
  #  n.cells <- n.cells + 1 # For bottom, tracked cell
  #  c.thk <- c(c.thk, thk - sum(c.thk))
  #}

  # Initial total masses, for calculating relative emission
  imass.NH3 <- sum(1000*tot['NH3']*c.thk)
  imass.CO2 <- sum(1000*(tot['CO2'] + tot['H2CO3'])*c.thk) 

  time.start <- Sys.time()
  n.calls.der <<- 0
  n.calls.spec <<- 0
  
  tot.lb <- tot[2:4] # Lower boundary constant concentrations, only for "kinetic" species NH3, H2CO3, CO2
  if(length(c.thk) == 1) c.thk <- rep(c.thk, n.cells)
  if(length(tot) == 8) tot <- rep(tot, each = n.cells)
  pos <- cumsum(c.thk) - 0.5*c.thk  # position (m), based on cell centers
  dx <- c(c.thk[n.cells], diff(pos))
  if(h.m.constant) {
    h.m[2] <- h.m[1] 
  } else {
    h.m[2] <- h.m[1]*sqrt(17.0305/44.0095)
  }
  if(!CO2.emis) h.m[2] <- 0
  names(h.m) <- c('NH3', 'CO2') # Calculate CO2 h.m based on that for NH3, based on Eq. 6 in Lee et al 2004

  if(lb == 'mx') dx[n.cells] <- 0.5*c.thk[n.cells-1] # Eliminates resistance within bottom cell

  # Global variables for tracking changes and avoiding speciation calculations if possible
  tot.glob <<- matrix(tot, nrow = n.cells, dimnames = list(1:n.cells, c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')))
  a.glob <<- m.glob <<- matrix(rep(999, 13*n.cells), nrow = n.cells)
  colnames(a.glob) <<- colnames(m.glob) <<- c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')

  # Initial speciation
  ll<- -14
  ul<- -0
  for (i in 1:n.cells) {
    tot.x <- as.numeric(tot.glob[i, ])
    names(tot.x) <- c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc')
    pred <- kinSpec(tot = tot.x, temp.c = temp.c, of = 'all', ll = ll, ul = ul)
    a.glob[i, ] <<- pred$a
    m.glob[i, ] <<- pred$m
  }

  y <- tot[names(tot) %in% c('NH3', 'H2CO3', 'CO2')]
  if(e1 == 0) y['CO2'][1] <- 0 # Set surface CO2 to zero

  # Call up ODE solver. Note that only NH3 and IC solutes change over time
  out <- ode.1D(y = y, times = times, func = derCalc, 
        parms = list(e1 = e1, h.m = h.m, n.cells = n.cells, temp.c = temp.c, c.thk = c.thk, dx = dx, D = D, kin.mult = kin.mult, lb = lb, p.CO2.a = p.CO2.a, tot.lb = tot.lb), 
        nspec = 3, dimens = n.cells) #, method = 'lsodes') #, hini = 0.01)

  tot.out <- array(as.vector(out[, colnames(out) %in% c('NH3', 'H2CO3', 'CO2')]), dim = c(dim(out)[1], n.cells, 3), dimnames = list(t = out[, 'time'], pos = pos, msp = c('NH3', 'H2CO3', 'CO2')))
  tot.f.out <- matrix(tot[names(tot) %in% c('H.', 'K.', 'Na.', 'Cl.', 'HAc')], nrow = n.cells, dimnames = list(pos = pos, msp = c('H.', 'K.', 'Na.', 'Cl.', 'HAc')))
  act.out <- array(as.vector(out[, grep('^act', colnames(out))]), dim = c(dim(out)[1], n.cells, 13), dimnames = list(t = out[, 'time'], pos = pos, sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  mol.out <- array(as.vector(out[, grep('^mol', colnames(out))]), dim = c(dim(out)[1], n.cells, 13), dimnames = list(t = out[, 'time'], pos = pos, sp = c('H.', 'NH3', 'H2CO3', 'CO2', 'K.', 'Na.', 'Cl.', 'HAc', 'OH.', 'NH4.', 'HCO3.', 'CO3.2', 'Ac.')))
  ph.out <- matrix(as.vector(out[, grep('pH', colnames(out))]), nrow = dim(out)[1], ncol = n.cells, dimnames = list(t = out[, 'time'], pos = pos))
  emis.NH3<- sum(1000*tot['NH3']*c.thk) - apply(1000*tot.out[, ,'NH3'], 1,function(x) sum(x*c.thk)) 
  emis.CO2<- sum(1000*(tot['CO2'] + tot['H2CO3'])*c.thk) - apply(1000*(tot.out[, ,'CO2'] + tot.out[, ,'H2CO3']), 1,function(x) sum(x*c.thk)) 
  emis.NH3.rel<- emis.NH3/imass.NH3
  emis.CO2.rel<- emis.CO2/imass.CO2
  p.CO2.s.out <- out[, grep('p.CO2.s', colnames(out))]
  p.NH3.s.out <- out[, grep('p.NH3.s', colnames(out))]
  emis.out <- data.frame(t = out[, 'time'], j.NH3 = out[, 'j.NH3'], j.CO2 = out[, 'j.CO2'], e.NH3 = emis.NH3, e.CO2 = emis.CO2, e.NH3.rel = emis.NH3.rel, e.CO2.rel = emis.CO2.rel, p.CO2.s = p.CO2.s.out, p.NH3.s = p.NH3.s.out)
  kin.CO2.out <- matrix(as.vector(out[, grep('kin.CO2', colnames(out))]), nrow = dim(out)[1], ncol = n.cells, dimnames = list(t = out[, 'time'], pos = pos))
  dCO2.out <- matrix(as.vector(out[, grep('dCO2.dt', colnames(out))]), nrow = dim(out)[1], ncol = n.cells, dimnames = list(t = out[, 'time'], pos = pos))
  dH2CO3.out <- matrix(as.vector(out[, grep('dH2CO3.dt', colnames(out))]), nrow = dim(out)[1], ncol = n.cells, dimnames = list(t = out[, 'time'], pos = pos))
  dNH3.out <- matrix(as.vector(out[, grep('dNH3.dt', colnames(out))]), nrow = dim(out)[1], ncol = n.cells, dimnames = list(t = out[, 'time'], pos = pos))
  pars <- c(D = D, h.m = h.m, kin.mult = kin.mult, lb = lb, n.cells = n.cells, p.CO2.a = p.CO2.a, thk = sum(c.thk), temp.c = temp.c)
  if (sum(times %in% c('0', '3600', '86400')) == 3) {
    summ.0 <- emis.out[times == 0, -1];names(summ.0) <- paste(names(summ.0), '.0', sep = '')
    summ.1h <- emis.out[times == 3600, -1];names(summ.1h) <- paste(names(summ.1h), '.1h', sep = '')
    summ.1d <- emis.out[times == 86400, -1];names(summ.1d) <- paste(names(summ.1d), '.1d', sep = '')
    ph.summ <- ph.out[times%in%c('0', '3600', '86400'), 1];names(ph.summ) <- c('ph.0', 'ph.1h', 'ph.1d')
  } else summ.0 <- summ.1h <- summ.1d <- ph.summ <- NA
  summ <- unlist(c(summ.0, summ.1h, summ.1d, ph.summ))

  exe.time <- Sys.time() - time.start
  if(of == 'all') list(times = times, pos = pos, tot.k = tot.out, tot.f = tot.f.out, act = act.out, mol = mol.out, ph = ph.out, kin.CO2 = kin.CO2.out, 
            emis = emis.out, exe.time = exe.time, n.calls = c(der = n.calls.der, spec = n.calls.spec), pars = pars, summ = summ)
            #dCO2 = dCO2.out, dH2CO3 = dH2CO3.out, dNH3 = dNH3.out, 
}
