## Internal modern-theme helpers for plot.fect() and esplot().
## Not exported. Loaded by esplot() / plot.R when theme.bw = TRUE.

# Color constants for the modern recipe
.MODERN_EST_COLOR    <- "grey20"   # connected line / point default
.MODERN_REF_COLOR    <- "grey75"   # y=0 hline / treatment-onset vline
.MODERN_PLACEBO_PT   <- "#E07A2B"  # orange — placebo highlight point
.MODERN_PLACEBO_FILL <- "#F5C8A8"  # peach — placebo window background
.MODERN_CARRYOVER_PT   <- "#2B6CB0" # blue — carryover highlight point
.MODERN_CARRYOVER_FILL <- "#C7DEF0" # light-blue — carryover window background

.modern_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = base_size, face = "plain",
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size - 2, color = "grey40",
                                                 margin = ggplot2::margin(b = 8)),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(8, 12, 8, 8),
      ## Compact legend so it doesn't dominate the panel in subfigure use:
      ## tight top margin against the panel, small text, small keys.
      legend.position    = "bottom",
      legend.margin      = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.text        = ggplot2::element_text(size = base_size - 3),
      legend.title       = ggplot2::element_text(size = base_size - 3),
      legend.key.size    = grid::unit(10, "pt"),
      legend.spacing.x   = grid::unit(4, "pt"),
      legend.spacing.y   = grid::unit(0, "pt")
    )
}

## Lighten a color toward white. Used to auto-derive rect-fill from a point color
## (e.g., orange #E07A2B -> peach #F5C8A8).
.lighten_color <- function(col, amount = 0.75) {
  if (is.null(col) || length(col) == 0) return(NA_character_)
  pal <- grDevices::colorRampPalette(c(col, "white"))(100)
  pal[max(1L, min(100L, round(amount * 100)))]
}
