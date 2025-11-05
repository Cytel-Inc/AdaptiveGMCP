# ThisPkgDepends.R
library(pkgdepends)

deps <- pkgdepends::new_pkg_deps("local::.") # Use local for packages under dev

# Gather and interpret metadata about all packages mentioned (directly or indirectly)
deps$resolve()

# Optional step to inspect the resolved metadata
# deps$get_resolution()

# Use that metadata to find a consistent set of package versions that satisfy
# all dependency constraints simultaneously.
deps$solve()

sol <- deps$get_solution() # Full solution, contains lots of columns and data

# Focusing on relevant columns
sol_min <- sol$data[, c("package", "version", "license")]

str(sol_min) # Solution table

# Visualizing dependencies as a pretty ASCII tree of direct + transitive deps
deps$draw()

# Producing a Software Bill of Materials (SBOM)
# Keep the most useful columns; keep dep_types/direct so you can filter later
bom <- sol$data[, c("package", "version", "sources", "ref", "license", "dep_types", "direct")]

library(dplyr)

bom <- bom %>% mutate(across(where(is.list), ~ vapply(.x,
                function(x) if(length(x) == 0) "" else paste(x, collapse = ";"),
                character(1))))

write.csv(bom, "bom_pkgdepends.csv", row.names = FALSE)
