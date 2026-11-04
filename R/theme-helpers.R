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

## Optional manuscript recipe. All sizes are points; cex.* multipliers are
## resolved by plot.fect before this theme is applied. This is visual only.
.publication_theme <- function(label_size = 18, axis_size = 12, legend_size = 11) {
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "grey45", linewidth = 0.35),
    axis.title = ggplot2::element_text(size = label_size, colour = "black"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
    axis.text = ggplot2::element_text(size = axis_size, colour = "black"),
    axis.text.x = ggplot2::element_text(size = axis_size),
    axis.text.y = ggplot2::element_text(size = axis_size),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = legend_size),
    legend.key.height = grid::unit(9, "pt"),
    legend.key.width = grid::unit(16, "pt"),
    legend.margin = ggplot2::margin(0, 0, 2, 0),
    legend.box.spacing = grid::unit(2, "pt"),
    plot.margin = ggplot2::margin(14, 18, 14, 18)
  )
}

## Lighten a color toward white. Used to auto-derive rect-fill from a point color
## (e.g., orange #E07A2B -> peach #F5C8A8).
.lighten_color <- function(col, amount = 0.75) {
  if (is.null(col) || length(col) == 0) return(NA_character_)
  pal <- grDevices::colorRampPalette(c(col, "white"))(100)
  pal[max(1L, min(100L, round(amount * 100)))]
}
