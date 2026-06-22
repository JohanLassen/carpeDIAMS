# Package-internal cache for Python module proxies, populated by .onLoad.
# Access from other package files as .pkg_py$apricot / .pkg_py$np.
.pkg_py <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    .pkg_py$apricot <- reticulate::import("apricot", delay_load = TRUE)
    .pkg_py$np      <- reticulate::import("numpy",   delay_load = TRUE)
}