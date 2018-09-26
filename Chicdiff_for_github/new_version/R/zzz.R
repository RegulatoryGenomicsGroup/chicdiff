
.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("")
  packageStartupMessage("Welcome to Chicdiff - version ", utils::packageDescription("Chicdiff", fields="Version"))
  packageStartupMessage('If you are new to Chicdiff, please consider reading the vignette through the command: vignette("Chicdiff").')
  
}