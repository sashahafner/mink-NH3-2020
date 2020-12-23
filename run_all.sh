# Shell script to run complete mink slurry analysis
# S. Hafner
# 17 December 2020

# Titration
cd PHREEQC/scripts/
phreeqc titrat_sim1.phrq
cd ../..

# Volatilization/emission simulations
cd emis_model/scripts
R CMD BATCH main.R
cd ../..

# Create report (figures updated by two calls above)
cd report
pandoc Hafner_report.md -o Hafner_report.pdf
cd ..

# See report.pdf for results


