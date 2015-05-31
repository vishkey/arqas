#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
#--------------------------------------------------------------------
#   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
    pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "arqas"),
                      fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      paste("\n Loading package arqas:", pkg.info["Title"], "\n"),
      paste(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n", sep=""),
      " Copyright Borja Varela 2015.\n")
}
