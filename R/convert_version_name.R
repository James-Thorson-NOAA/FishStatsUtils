#' Convert string to semantic version number
#'
#' @param char
#'
#' @return A semantic version numvber
#'
#'
#' @export
convert_version_name <- function( char, split="_" ) {
  convert = function(part) paste0(na.omit(as.numeric(strsplit(part,"")[[1]])),collapse="")

  Split = strsplit( char, split=split )[[1]]
  Codes = suppressWarnings(sapply( Split, FUN=convert))
  Codes = paste0( na.omit(ifelse(Codes=="",NA,Codes)), collapse="." )
  semantic_version = numeric_version(Codes)

  return(semantic_version)
}
