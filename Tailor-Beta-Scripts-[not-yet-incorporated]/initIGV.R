#' Initiate the Java IGV library
#'
#' @return The pathname to the igvtools.jar or igv.jar file used.
#'
#' @aliases findIGV findIGVTools
#' @export
#' @importFrom utils file_test
initIGV <- function() {
  pathnames <- c(findIGVTools(), findIGV())
  pathname <- pathnames[!is.na(pathnames)][1]
  if (!file_test("-f", pathname)) {
    stop("Failed to located igvtools.jar or igv.jar. See help('initIGV').")
  }
  rJava::.jaddClassPath(pathname)

  invisible(pathname)
}

#' @export
#' @importFrom utils file_test
findIGVTools <- function() {
  ## Locate igvtools.jar
  paths <- c(Sys.getenv("IGVTOOLS_HOME"),
             dirname(Sys.which(c("igvtools.jar", "igvtools", "igvtools.bat"))))
  paths <- unique(paths[file_test("-d", paths)])
  pathnames <- file.path(paths, "igvtools.jar")
  pathnames <- pathnames[file_test("-f", pathnames)]
  if (length(pathnames) > 0) return(pathnames[1])

  NA_character_
}

#' @export
#' @importFrom utils file_test
findIGV <- function() {
  paths <- c(Sys.getenv("IGV_HOME"),
             dirname(Sys.which(c("igv.jar", "igv.sh", "igv.bat", "igv.command"))))
  paths <- unique(paths[file_test("-d", paths)])
  pathnames <- file.path(paths, "igv.jar")
  pathnames <- pathnames[file_test("-f", pathnames)]
  if (length(pathnames) > 0) return(pathnames[1])

  NA_character_
}
