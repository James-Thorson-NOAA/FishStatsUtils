#' Determine latest version of VAST
#'
#' @param version The default is \code{NULL}, which will cause the function
#' to look for the latest version of the \code{.cpp} file to compile.
#' If a version is supplied, then \code{R} will look for that version name
#' within the appropriate folder the disk.
#'
#' @return A full file path as a character value supplying the location of
#' the latest, or supplied, version of the VAST code to compile.
#'
#' @author Kelli Faye Johnson
#'
#' @export
get_latest_version <- function(version = NULL, package = "VAST") {

  # Determine location of files on machine
  thedir <- system.file("executables", package = package)

  # Determine list of available files
  if (!is.null(version)) {
    thefile <- dir(thedir, pattern = version, ignore.case = TRUE, full.names = FALSE)
    if( length(thefile)==0 ){
      stop("The file ", version, " was not found in the dir, ", thedir, ".")
    }
  }else{
    thefile <- dir(thedir, full.names = FALSE, pattern = "\\.cpp")
    if( length(thefile)==0 ){
      stop("cpp files were not found in the dir, ", thedir, ".")
    }
  }

  # Determine which is latest version
  for(i in 1:length(thefile)){
    if(i==1) semantic_version = convert_version_name(thefile[1])
    if(i>=2) semantic_version = c( semantic_version, convert_version_name(thefile[i]) )
  }
  thefile = thefile[which(semantic_version==max(semantic_version))]

  # Remove .cpp from end and return
  thefile <- strsplit( thefile, ".", fixed=TRUE )[[1]][1]
  return(thefile)
}
