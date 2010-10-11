'.onLoad' <- function(lib, pkg="gof")
  {
    desc <- packageDescription(pkg)
    cat("\nLoading '", desc$Package, "' package...\n", sep="")
    cat("Version    : ", desc$Version, "\n\n", sep="")
  }
