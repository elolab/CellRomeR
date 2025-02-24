
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

#'
#'
#'
excld.pattern.process <- function(excld.pattern, StdP = NULL) {
  
  # Parse exclusion pattern with  "standard_varnames" option
  if (is.null(StdP)) {
    StdP = "^LABEL|^label|ID$|ID[1-9]$|id$|INDEX|^TIME|QUALITY|LOCATION|START|STOP|GAP$|^NUMBER|DURATION|^POSITION|^FRAME|VISIBILITY|LINK_COST|^EDGE_TIME$"
  }
  
  if (is.null(excld.pattern)) {
    excld.pattern = "NO_EXCLUSION"
  } else if (any(excld.pattern %in% c("standard_varnames","standard_variables","standard_columns") )) {
    excld.pattern <- excld.pattern[-grep("standard_", excld.pattern)]
    excld.pattern <- combine_patterns(c(excld.pattern,StdP), logic = "OR")
  } else if (!is.null(excld.pattern)) {
    excld.pattern = combine_patterns(excld.pattern)
  } else excld.pattern = "NO_EXCLUSION"
  return(excld.pattern)
}

