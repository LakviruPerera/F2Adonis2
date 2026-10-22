## The functions below are internal helpers of the 'vegan' package that
## F2_adonis2() needs in order to reproduce the behavior of vegan::adonis2()
## (formula parsing, distance adjustment, permutation setup). They are not
## exported by vegan, so they are retrieved at load time with
## utils::getFromNamespace(). The package requires vegan (>= 2.7-3) and will
## be updated if these functions change.


# 1. Initialize hidden placeholders at the package level
addCailliez      <- NULL
addLingoes       <- NULL
getPermuteMatrix <- NULL
initDBRDA        <- NULL
ordiParseFormula <- NULL
ordiTerminfo     <- NULL

# 2. This runs automatically when the package loads
.onLoad <- function(libname, pkgname) {
  addCailliez      <<- utils::getFromNamespace("addCailliez", "vegan")
  addLingoes       <<- utils::getFromNamespace("addLingoes", "vegan")
  getPermuteMatrix <<- utils::getFromNamespace("getPermuteMatrix", "vegan")
  initDBRDA        <<- utils::getFromNamespace("initDBRDA", "vegan")
  ordiParseFormula <<- utils::getFromNamespace("ordiParseFormula", "vegan")
  ordiTerminfo     <<- utils::getFromNamespace("ordiTerminfo", "vegan")
}
