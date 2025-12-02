# ThisPkg_renv.R
# renv workflow for this package

# install.packages("renv")
library(renv)

### ONE TIME STEPS FOR THE FIRST DEVELOPER ###################
### One of the package developers needs to create the renv folder and lockfile
### (renv.lock) using the following steps:
  # Get current list of dependencies
  renv::dependencies()

  # Creating renv folder and lockfile
  renv::init(settings = list(use.cache = TRUE))

  # Snapshot the current working state
  renv::snapshot() # Using implicit snapshot type to capture any package
                   # actually used in the codebase

  # To test if renv is properly set up
  renv::restore()
  devtools::check()

  # # Since DESCRIPTION file is now curated, we are switching to explicit snapshot
  # # type
  # renv::settings$snapshot.type("explicit")
  # renv::snapshot()

### ONE TIME WORKFLOW FOR OTHER DEVELOPERS ##################
### Once one of the developers has performed the above steps and published the
### package changes, other developers need to do the following after getting the
### changes:
  renv::restore()
  # This will read the lockfile and install the same versions of each package
  # into the local environment of the developer
  # The developer can then optionally run the following:
  devtools::check()

### FUTURE renv WORKFLOW FOR ALL DEVELOPERS ###################
### Whenever any developers takes an action that may impact the package
### dependencies, they need to do the following:

  ## CHECKING IF ANYTHING HAS CHANGED W.R.T LOCKFILE
  renv::status()
  # This shows packages that differ from the lockfile. Resolve any differences
  # that you may see.

  ## INSTALLING A NEW PACKAGE >>>>>>>>>>>>>>>>
  # If you want to install a new package:
  renv::install("<pkg name>") # Use this instead of install.packages()
  renv::snapshot()

  ## PINNING A SPECIFIC VERSION OF A DEPENDENT PACKAGE
  # If your package depends on a specific version of a package, you can pin it
  # in the lockfile:
  renv::record("dplyr@1.1.4")
  renv::snapshot()

  ## UPATING A PACKAGE TO USE ITS LATEST VERSION
  renv::update("dplyr")
  renv::snapshot()



##########################################################
# PRODUCTING packrat.lock FILE FOR BLACK DUCK ############
  # As of Nov 2025, Black Duck detector scan works only if the R package
  # contains a packrat.lock file. Black Duck currently does not work with
  # the lockfile produced by renv.

  # save as dev/make_packrat_lock_from_renv.R and run from project root
  json <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
  pkgs <- json$Packages

  # helper to map renv "Source" to Packrat "Source"
  map_source <- function(rec) {
    src <- rec$Source %||% "Repository"
    if (src %in% c("Repository","CRAN")) "CRAN" else if (src == "Bioconductor") "Bioconductor" else if (src == "GitHub") "GitHub" else src
  }

  # header (Packrat lockfile metadata)
  lines <- c(
    "PackratFormat: 1.4",
    paste0("RVersion: ", getRversion()),
    # Repos line is optional but nice to include if in lock
    if (!is.null(json$R$Repositories)) {
      repos <- vapply(json$R$Repositories, function(r) paste0(r$Name, "=", r$URL), character(1))
      paste0("Repos: ", paste(repos, collapse = ", "))
    } else NULL,
    "" # blank line before packages
  )

  # package entries
  for (name in sort(names(pkgs))) {
    rec <- pkgs[[name]]
    src <- map_source(rec)
    lines <- c(lines,
               paste0("Package: ", name),
               paste0("Source: ", src),
               paste0("Version: ", rec$Version),
               # include GitHub metadata if present
               if (!is.null(rec$RemoteType)) paste0("RemoteType: ", rec$RemoteType) else NULL,
               if (!is.null(rec$RemoteHost)) paste0("RemoteHost: ", rec$RemoteHost) else NULL,
               if (!is.null(rec$RemoteUsername)) paste0("RemoteUsername: ", rec$RemoteUsername) else NULL,
               if (!is.null(rec$RemoteRepo)) paste0("RemoteRepo: ", rec$RemoteRepo) else NULL,
               if (!is.null(rec$RemoteRef)) paste0("RemoteRef: ", rec$RemoteRef) else NULL,
               if (!is.null(rec$RemoteSha)) paste0("RemoteSha: ", rec$RemoteSha) else NULL,
               "" # blank line between package stanzas
    )
  }

  writeLines(lines, "packrat.lock")
  cat("Wrote packrat.lock with ", length(pkgs), " packages\n")
##########################################################
