
#### Auxiliary functions ####

#### String joining for patterns ####
##
#' #### Function for string combinations
#'
#' Combine patterns for search pattern
#' Or and AND for search pattern combinations
#'
#' @param patterns a vector of strings to combine into search pattern
#' @param logic used in combination
#' @param custom used in combination
#'
#' @examples
#' patterns = paste0("x",seq(1:4))
#' combine_patterns(patterns, "OR" )
#' combine_patterns(patterns, "AND" )
#' combine_patterns(patterns, custom = "_")
#'
#' @export
combine_patterns <- function(patterns = "", logic = c("OR","AND"), custom = NULL, exacts = FALSE) {

  # choices
  logic = match.arg(logic, c("OR","AND"))
  if (!is.null(custom)) {
    return(paste(patterns, collapse=custom))
  }
  logic <- match.arg(logic)
  if (logic == "OR" & !exacts) {
    return(paste(patterns, collapse="|"))
  }
  if (logic == "AND" & !exacts) {
    return(paste(patterns, collapse="&"))
  }

  if (logic == "OR" & exacts) {
    patterns <- unlist(strsplit(patterns, "\\|") )
    patterns <- gsub('[\\^\\$]', '', patterns)
    patterns <- paste0("^",patterns,"$")
    patterns <- paste(patterns, collapse="|")
    return(patterns)
  }
  if (logic == "AND" & exacts) {
    print("Not yet implemented!")
  }

  print("No logic or custom found.")
}
