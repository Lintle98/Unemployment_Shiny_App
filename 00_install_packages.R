# 00_install_packages.R - run ONCE (click "Source")
# Installs any packages the dashboard needs that are not installed yet.
pkgs <- c("shiny", "bslib", "sf", "dplyr", "spdep", "spatialreg",
          "leaflet", "plotly", "DT", "rsconnect")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing) else message("All packages already installed.")
