.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
      paste0("relabelKL version ", "0.0.1")
  )
}

.onUnload <- function (libpath) {

  library.dynam.unload("relabelKL", libpath)

}
