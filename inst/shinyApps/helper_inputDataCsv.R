# This function converts an R object to a string so it can be sent to a csv.
# Currently, conversion of only list, dataframe and matrix is supported.

convertToString <- function(x) {
  if (is.list(x)) {
    elements <- paste0(names(x), " = ", sapply(x, deparse))
    result <- paste("list(", paste(elements, collapse = ", "), ")")
  } else if (is.data.frame(x)) {
    columns <- paste0(names(x), " = ", sapply(x, function(col) {
      if (is.character(col)) {
        paste0("c(", paste(shQuote(col, type = "cmd"), collapse = ", "), ")")
      } else {
        paste0("c(", paste(col, collapse = ", "), ")")
      }
    }))
    result <- paste("data.frame(", paste(columns, collapse = ", "), ")")
  } else if (is.matrix(x)) {
    values <- paste(c(x), collapse = ", ")
    result <- paste("matrix(c(", values, "), nrow = ", nrow(x),
                    ", ncol = ", ncol(x), ", byrow = TRUE)")
  } else {
    result <- deparse(x)
  }
  return(result)
}


ifElse <- function(test, yes, no) {
  if (test) {
    yes
  } else {
    no
  }
}
