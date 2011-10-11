'.onAttach' <- function(lib, pkg="gof")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n\n")
  }
