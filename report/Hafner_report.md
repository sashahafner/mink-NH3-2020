---
title: 'Prediction of ammonia emission from alkalinized mink slurry'
output: pdf_document
author: Sasha D. Hafner[^1]
date: December 22, 2020
---

[^1]:Hafner Consulting LLC, Reston, VA, USA *and* Aarhus University Department of Engineering, Aarhus, Denmark | `sasha@hafnerconsulting.com` 

# 1. Overview
This document addresses the problem of ammonia (NH~3~) volatilization (emission) from liquid mink manure (slurry) that has been treated with lime (Ca(OH)~2~) to inactivate SARS-CoV-2 particles. 
In this work, a reaction-transport model was applied to predict volatilization rates during mixing and storage without mixing.
Results provide estimates of upper and lower limits on likely emission from mink slurry.

# 2. Mink slurry composition
Slurry composition affects the lime (Ca(OH)~2~) dose required for pH elevation as well as NH~3~ emission rate and dynamics. 
The most important parameters are total ammoniacal N (TAN) concentration, which represents the pool of available NH~3~ but also controls pH; total inorganic carbon (TIC) concentration, which also contributes buffering; and pH.
In this work untreated slurry is assumed to have a TAN concentration of 6.6 g/kg (per kg total mass).
Very few measurements of TIC are available.
Mink slurry from a single source used in two Danish field trials at AU Foulum in 2015 had an average concentration of 240 mmol/kg, and this value was used here. 
Other solutes that contribute to buffering are present in mink slurry, including organic acids (e.g., acetic acid) and orthophosphate (PO~4~$^{-3}$ and related species), but they are less important.

\newpage

:Mink slurry composition assumed in this work.

Analyte      | Value  | Unit   
:----:       | :----: | :----: 
TAN          | 470    | mmol/kg
TIC          | 240    | mmol/kg
pH           | 7.5    | -      
pH (treated) | 12     | -      

# 3. Titration
Assuming the composition of mink slurry matches the values given in Table 1, charge can be balanced with addition of 230 mmol/kg chloride Cl$^-$ (effectively as hydrochloric acid, HCl).
Titration of this solution with Ca(OH)~2~ was simulated using the PHREEQC reaction-transport model from USGS, with the Amm.dat database (Parkhurst, 2020).
Results show that about 360 mmol/kg Ca(OH)~2~ (almost 30 g/kg, or 30 kg/t) is required to reach pH 12 (Fig. 1).
Both calcite and aragonite (both forms of CaCO~3~) are over-saturated immediately even with small additions of Ca(OH)~2~, but the rate at which precipitates form probably limits their importance in the short term.
Dolomite (CaMg(CO~3~)~2~) may also form, but will be limited by what are probably low Mg concentrations.
Regardless, the effect on pH change is probably minimal, as in Fig. 1.
The presense of other solutes ignored here will influence buffering, and probably tend to increase the required dose of lime.

![Simulation of mink slurry titration (Table 1) with calcium hydroxide. Black line excludes solid phase formation, gray line includes aragonite precipitation, assuming equilibrium.](../figs/titrat1.pdf)

\newpage

# 4. Volatilization predictions
## 4.1. Methods
The model described in Hafner et al. (2012) was used to predict NH~3~ emission. 
This model is a relatively complete (i.e., few simplifications, relatively accurate description of the physical/chemical system) numerical process-based model, including pH-dependent equilibrium, kinetically-limited inorganic carbon reactions, and diffusive transport within the slurry.
Emission from the exposed upper surface is predicted using a mass transfer coefficient approach.
Slurry may be assumed to be well-mixed below a liquid film, similar to the two-film model.
However, the numerical model is dynamic, so does not require any assumptions about steady state, and can include non-equilibrium effects, along with pH gradients and resulting nonlinear chemical species gradients within the liquid film.

For untreated slurry, the model predicts an elevation of surface pH due to CO~2~ loss, and in some cases a *reduction* in emission with mixing due to more transport of HCO~3~$^-$ to the surface, where it will moderate the pH gradient.
For treated slurry with a pH of 12, CO~2~ emission is nonexistent--actually, slight *absorption* of CO~2~ from the atmosphere will occur.
Therefore, a somewhat simpler equilibrium version can be used, where all inorganic carbon species are in equilibrium.
Loss of NH~3~ from the surface will result in a decline in surface pH, reducing emission.

Simulations of emission during both mixing and storage were carried out.
Untreated slurry was included as a reference.
Both treated and untreated slurry had the composition given in Table 1.

Four liquid conditions, with a liquid film in three (thickness varying from 20 $\mu$m to 5 cm) and one with no mixing (effectively a 3 m film) were considered (Table 2).
Film thickness (more accurately, liquid-side mass transfer resistance, since the use of films is a simplification) is difficult to accurately determine for slurry.
The lowest of the values used (20 $\mu$m) is close to conditions for the surface of the sea (Tsunogai and Tanka, 1980), or during mechanical surface aeration (calculated based on results in Roberts and Daendliker, 1983).
This value is unlikely to be exceeded or even met in mechanically mixed slurry.
At the other extreme, completely still slurry is effectively impossible, so it provides a convenient upper limit on liquid-side resistance.
Even with no mechanical mixing, slurry addition, or gas bubble formation, wind and temperature change at the surface will cause some mixing.

The air-side mass transfer coefficient estimated from the correlation given by US EPA (1994, Table 5-1) is 0.01 m/s with a wind speed of 3-4 m/s (close to mean wind speed for Denmark) and a 10-20 m tank diameter.
Assuming the resistance associated with a straw cover is approximately 100 s/m (Olesen and Sommer, 1993, Table 4), this value would be halved by adding a layer of straw (0.005 m/s).
These two values were paired with liquid-side values for film thickness in six combinations of parameter values (Tables 2 and 3), which were applied to both untreated (pH 7.5) and treated slurry (pH 12) (Table 1).

: Liquid-side film thickness used for emission prediction (all in m). For the still condition, there is no well-mixed layer (all solute transport within the slurry is by diffusion).

Condition          | Value 
:----:             | :----
Very intense mixing| 0.00002
Intense mixing     | 0.0001
Moderate mixing    | 0.01
Limited mixing     | 0.05
No mixing          | 3


: Air-side mass transfer coefficient values used for emission predictions (all in m/s). Values are based on wind speeds between 3 and 4 m/s.

Air-side condition     | Value 
:----:                 | :----
Straw cover            | 0.006 
Uncovered              | 0.01  

: Other variable or parameter values fixed for all simulations.

Variable              | Value  | Unit 
:----:                | :----: | :----:
Temperature           | 10     | &deg;C
Storage depth         | 3      | m
CO~2~ partial pressure| 0.0004 | atm

## 4.2. Results
Results are presented for treated slurry but also untreated slurry as a reference condition.

During mixing, the predicted NH~3~ flux from treated slurry tends to quickly decline as TAN is depleted from near the surface (Fig. 2). 
This change becomes less important as film thickness decreases.
Emission rate depends strongly on the film thickness within the slurry, and for thinner films (less liquid-side resistance), gas-side resistance has a moderate effect (Figs. 2 and 3).
Even with low liquid-side resistance (20 $\mu$m film), predicted total emission over 5 hours does not exceed 2.5% of initial TAN.
A more likely value is between 1 and 2%.
Shorter mixing periods will, of course, reduce emission (Fig. 3, Table 5).

Untreated slurry shows the highest predicted emission rates for the thickest film (1 cm). 
This surprising effect is due to the elevation of surface pH caused by volatilization of CO~2~.
With thinner films, transport of TIC from below reduces the pH change.
This response has been shown in laboratory measurements (Hafner et al., 2012, Figure S-2 in supplementary material) and is explained in detail in Hafner et al. (2012).

![Predicted NH~3~ flux during slurry mixing for both untreated (pH 7.5, top) and treated (pH 12, bottom) mink slurry over 5 hours.](../figs/flux_mixing.pdf)

![Predicted NH~3~ cumulative emission as a percentage of inital TAN mass during slurry mixing for both untreated (pH 7.5, top) and treated (pH 12, bottom) mink slurry over 5 hours.](../figs/emis_mixing.pdf)

Predicted NH~3~ flux from untreated slurry ranges from below 0.1 to above 1 g/h-m$^2$ (Fig. 4). 
Measured values from storage tanks for cattle and pig slurry from a recent literature review are similar to the lowest predicted values; measured values ranged from around 0.01 to almost 1 g/h-m$^2$, with median values of 0.05 and 0.09 g/h-m$^2$ (Sommer et al., in progress).
Kupper et al. (2020) reviewed NH~3~ emission rates from cattle and pig slurry in storage tanks and lagoons.
Median values were between 0.06 and 0.1 g/h-m$^2$ (Kupper et al., 2020, Table 7).
Considering the higher TAN concentration in mink slurry compared to cattle and pig, and the relatively low temperature used for these mink slurry predictions, measurements suggest that predictions based on limited mixing (5 cm film) may be the most accurate.

During long-term storage of treated slurry (5 months), predicted emission rate rapidly decreases as TAN is depleted from near the surface (Figs. 4 & 5).
With no mixing, predicted emission is only 4% of initial TAN over 5 months (Fig. 5).
However, this value probably represents an untainable lower limit, since mixing cannot be completely eliminated (See Section 4.1).
Without significant mixing (film thickness 1 cm or larger), predicted emission remains below 35% of initial TAN mass over 5 months, or 22% over 3 months (Fig. 5, Table 6).
Emission could be much higher--even approaching 100% of initial TAN--with intense mixing, highlighting the importance of minimizing slurry mixing.
Assuming a film thickness of 5 cm (limited/natural mixing), NH~3~ loss would be about 9% over 5 months, or under 6% over 3 months of storage (Fig. 5, Table 6). 


![Predicted NH~3~ flux during long-term slurry storage for both untreated (pH 7.5, top) and treated (pH 12, bottom) mink slurry over 5 months.](../figs/flux_store.pdf)

![Predicted NH~3~ cumulative emission as a percentage of initial TAN mass during long-term slurry storage for both untreated (pH 7.5, top) and treated (pH 12, bottom) mink slurry over 5 months.](../figs/emis_store.pdf)

Despite uncertainty in the degree of mixing within slurry storage tanks, the predictions presented here provide limits on the likely loss of NH~3~ during mixing and long-term storage of slurry.
Cumulative NH~3~ loss and average flux for what are thought to be the most likely scenarios are summarized below for mixing (Table 5) and storage (Table 6).
In addition to effects of mixing, the depth of slurry within the tank is important.
In all these predictions, depth was fixed at 3 m. 
Storage at shallower depths will increase emission, by increasing the area available for emission (per kg TAN) and reducing mass transfer resistance.
For very shallow depths, relative losses could be very high even for completely still slurry.

: Predicted cumulative emission and average flux by hour for alkalinized mink slurry (pH 12) undergoing mixing. For the predictions, the air-side mass transfer coefficient was taken as 0.01 m/s (uncovered slurry), and liquid-side film thickness 100 $\mu$m (intense mixing) (Table 2). Average flux is presented for each 1 hour interval (0-1 h, etc.).

Period (hour)   | Cum. emission (% init. TAN) | Ave. flux (g/h-m$^2$)
:----:       | :----:                 | :----: 
1            |0.36                    |71
2            |0.71                    |70
3            |1.1                     |70
4            |1.4                     |70 
5            |1.8                     |69

: Predicted cumulative emission and average flux by month for stored alkalinized mink slurry (pH 12). For the predictions, the air-side mass transfer coefficient was taken as 0.005 m/s (slurry covered with straw), and liquid-side film thickness 5 cm (limited/natural mixing) (Table 2). Average flux is presented for each 1 month interval (0-1 month, etc.).

Period (month)   | Cum. emission (% init. TAN) | Ave. flux (g/h-m$^2$)
:----:       | :----:                 | :----: 
1            |2.2                     |0.62
2            |3.9                     |0.46
3            |5.6                     |0.45
4            |7.2                     |0.44 
5            |8.8                     |0.44

\newpage 

# Conclusions
Predictions from a reaction-transport model show the following for mink slurry with pH adjusted to 12:

1. Ammonia volatilization losses during 5 hours of mixing is unlikely to exceed 2.5% of initial total ammoniacal nitrogen, and will probably be between 1 and 2%.
2. Ammonia volatilization losses during 3 months of storage is unlikley to exceed 22% of initial total ammoniacal nitrogen, and will probably be closer to 6%.

There are some uncertainties associated with these predictions.
In particular, mixing conditions and storage depth (assumed to be 3 m) are important.
If shallower storage is anticipated, these estimates should be revised.

# References
Hafner, S.D., Montes, F., Alan Rotz, C., 2012. The role of carbon dioxide in emission of ammonia from manure. Atmospheric Environment 66, 63–71. https://doi.org/10.1016/j.atmosenv.2012.01.026

Kupper, T., Häni, C., Neftel, A., Kincaid, C., Bühler, M., Amon, B., VanderZaag, A., 2020. Ammonia and greenhouse gas emissions from slurry storage - A review. Agriculture, Ecosystems & Environment 300, 106963. https://doi.org/10.1016/j.agee.2020.106963

Parkhurst, D., 2020. PHREEQC Version 3. https://www.usgs.gov/software/phreeqc-version-3. 

Roberts, P.V., Daendliker, P.G., 1983. Mass transfer of volatile organic contaminants from aqueous solution to the atmosphere during surface aeration. Environmental Science & Technology 17, 484–489. https://doi.org/10.1021/es00114a009

Tsunogai, S., Tanaka, N., 1980. Flux of oxygen across the air-sea interface as determined by the analysis of dissolved components in sea water. Geochemical Journal 14, 227–234. https://doi.org/10.2343/geochemj.14.227

US EPA, 1994. Air Emission Models for Waste and Wastewater. EPA-453/R-94-080A. US Environmental Protection Agency, Office of Air Quality Planning and Standards, Research Triangle Park, NC.
