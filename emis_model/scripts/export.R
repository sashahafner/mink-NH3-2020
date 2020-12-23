# Export results

# Emission curves
write.csv(emis, '../output/emission.csv', row.names = FALSE)

# All objects
save.image('../R_images/workspace.RData')

