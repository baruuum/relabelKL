.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
      paste0("relabelKL version ", "0.0.0.9000")
  )
}

.onUnload <- function (libpath) {

  library.dynam.unload("relabelKL", libpath)

}
