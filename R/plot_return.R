print.fect_plot_return <- function(x, ...) {
  # When return.data=TRUE we return a list; printing it should behave like a plot:
  # draw the figure, but don't dump $data to console.
  if (!is.null(x$p)) {
    suppressWarnings(print(x$p))
  }
  invisible(x)
}

