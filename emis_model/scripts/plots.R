# Plots

# Names
emis$mtc.nm <- factor(emis$mtc, levels = c(0.01, 0.005), 
                      labels = c('Uncovered (0.01 m/s)', 'Straw covered (0.005 m/s)'))
emis$film.nm <- factor(emis$film, levels = c(2E-5, 1E-4, 0.01, 0.05, 3), 
                       labels = c('Very intense mixing\n20 um', 'Intense mixing\n100 um',
                                  'Moderate mixing\n1 cm', 'Limited/natural mixing\n5 cm', 'No mixing\n3 m'))

# Mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dd <- subset(emis, dur < 86400 & !film %in% c(3, 0.05))

# Flux
ggplot(dd, aes(time.h, j.NH3.gm2h, group = id, colour = film.nm)) +
  geom_line() +
  geom_point(data = subset(dd, t ==0), pch = 1) +
  geom_hline(yintercept = 0, colour = 'gray55') +
  #scale_y_continuous(trans = 'log10') +
  facet_grid(condition ~ mtc.nm, scale = 'free') +
  labs(x = 'Time (h)', y = expression(NH[3]~'flux'~(g/h-m^'2')), colour = '') +
  theme(legend.position = 'top')
ggsave('../../figs/flux_mixing.pdf', height = 6, width = 6.5, scale = 0.8)

# Cumulative emission
ggplot(dd, aes(time.h, e.NH3.rel, group = env, colour = film.nm)) +
  geom_line() +
  facet_grid(condition ~ mtc.nm, scale = 'free') +
  labs(x = 'Time (h)', y = expression('Cumulative'~NH[3]~'emission (% of initial TAN)'), colour = '') +
  theme(legend.position = 'top')
ggsave('../../figs/emis_mixing.pdf', height = 6, width = 6.5, scale = 0.8)

# Stored ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dd <- subset(emis, dur > 86400 & film %in% c(0.01, 0.05, 3))

ggplot(dd, aes(time.d, j.NH3.gm2h, group = id, colour = film.nm)) +
  geom_line() +
  geom_point(data = subset(dd, t ==0), pch = 1) +
  geom_hline(yintercept = 0, colour = 'gray55') +
  #scale_y_continuous(trans = 'log10') +
  #annotation_logticks(side = 'l') +  
  facet_grid(condition ~ mtc.nm, scale = 'free') +
  labs(x = 'Time (d)', y = expression(NH[3]~'flux'~(g/h-m^'2')), colour = '') +
  theme(legend.position = 'top')
ggsave('../../figs/flux_store.pdf', height = 6, width = 6.5, scale = 0.8)

ggplot(dd, aes(time.d, e.NH3.rel, group = env, colour = film.nm)) +
  geom_line() +
  facet_grid(condition ~ mtc.nm, scale = 'free') +
  labs(x = 'Time (d)', y = expression('Cumulative'~NH[3]~'emission (% of initial TAN)'), colour = '') +
  theme(legend.position = 'top')
ggsave('../../figs/emis_store.pdf', height = 6, width = 6.5, scale = 0.8)

table(emis$env, emis$mtc.nm : emis$film.nm)


x$time.d
