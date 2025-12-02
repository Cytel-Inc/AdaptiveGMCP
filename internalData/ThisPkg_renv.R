# ThisPkg_renv.R
# renv workflow for this package

library(renv)

# # Get current list of dependencies
# renv::dependencies()
#
# # Creating renv folder and lockfile
# renv::init(settings = list(use.cache = TRUE))  # accept prompts
#
# # Snapshot the current working state
# # renv::snapshot() # Using implicit snapshot type to capture any package
#                    # actually used in the codebase
#
# # Since DESCRIPTION file is now curated, we are switching to explicit snapshot
# # type
# renv::settings$snapshot.type("explicit")
# renv::snapshot()
#
# # To test if renv is properly set up
# renv::restore()
# devtools::check()

### FUTURE renv WORKFLOW ###################################
# 05 Nov 2025: Now that renv has been set up and DESCRIPTION file has been
# properly modified, in future do the following:

## INSTALLING A NEW PACKAGE >>>>>>>>>>>>>>>>
renv::install("<pkg name>")
renv::snapshot()





