# 1. Initialize hidden placeholders at the package level
addCailliez      <- NULL
addLingoes       <- NULL
getPermuteMatrix <- NULL
initDBRDA        <- NULL
ordiParseFormula <- NULL
ordiTerminfo     <- NULL

# 2. This runs automatically the microsecond your package loads
.onLoad <- function(libname, pkgname) {
  addCailliez      <<- utils::getFromNamespace("addCailliez", "vegan")
  addLingoes       <<- utils::getFromNamespace("addLingoes", "vegan")
  getPermuteMatrix <<- utils::getFromNamespace("getPermuteMatrix", "vegan")
  initDBRDA        <<- utils::getFromNamespace("initDBRDA", "vegan")
  ordiParseFormula <<- utils::getFromNamespace("ordiParseFormula", "vegan")
  ordiTerminfo     <<- utils::getFromNamespace("ordiTerminfo", "vegan")
}
