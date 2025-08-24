## new plot
# x: a fect object
# type of the plot; axes limits; axes labels;
# main: whether to show the title;
# id: plot a part of units

plot.fect <- function(
    x,
    type = NULL, # gap, equiv, status, exit, factors, loadings, calendar, counterfactual, heterogeneous
    restrict = "rm",
    loo = FALSE,
    highlight = NULL, ## for carryover test and placebo test
    plot.ci = NULL, ## "0.9", "0.95", "none"
    show.points = TRUE,
    show.group = NULL,
    bound = NULL, # "none", "min", "equiv", "both"
    show.count = TRUE,
    proportion = 0.3, # control the xlim
    pre.periods = NULL, # for testing
    f.threshold = NULL, # equiv f
    tost.threshold = NULL, # pre-trend placebo carryover
    effect.bound.ratio = FALSE,
    stats = NULL, ## "none", "F.p", "F.equiv.p", "placebo.p", "carryover.p", "equiv.p"
    stats.labs = NULL,
    raw = "none", ## "none", "band", "all"
    vis = NULL,
    main = NULL,
    xlim = NULL,
    ylim = NULL,
    xlab = NULL,
    ylab = NULL,
    xangle = NULL,
    yangle = NULL,
    xbreaks = NULL,
    ybreaks = NULL,
    xticklabels = NULL,
    yticklabels = NULL,
    gridOff = NULL,
    legendOff = FALSE,
    legend.pos = NULL,
    legend.nrow = NULL,
    legend.labs = NULL,
    stats.pos = NULL,
    show.stats = TRUE,
    theme.bw = TRUE,
    nfactors = NULL,
    include.FE = TRUE,
    id = NULL,
    cex.main = NULL,
    cex.main.sub = NULL,
    cex.axis = NULL,
    cex.lab = NULL,
    cex.legend = NULL,
    cex.text = NULL,
    axis.adjust = FALSE,
    axis.lab = "both",
    axis.lab.gap = c(0, 0),
    shade.post = FALSE,
    start0 = FALSE,
    return.test = FALSE,
    balance = NULL,
    weight = NULL,
    lcolor = NULL,
    lwidth = NULL,
    ltype = NULL,
    line.color = NULL,
    line.width = NULL,
    count = NULL,
    preset = NULL,
    connected = NULL,
    ci.outline = FALSE,
    color = NULL,
    est.lwidth = NULL,
    est.pointsize = NULL,
    count.color = NULL,
    count.alpha = NULL,
    count.outline.color = NULL,
    placebo.color = NULL,
    carryover.color = NULL,
    carryover.rm.color = NULL,
    sens.original.color = NULL,
    sens.colors = NULL,
    counterfactual.color = NULL,
    counterfactual.raw.controls.color = NULL,
    counterfactual.raw.treated.color = NULL,
    counterfactual.linetype = NULL,
    box.control = NULL,
    box.treat = NULL,
    calendar.color = NULL,
    calendar.lcolor = NULL,
    calendar.cicolor = NULL,
    heterogeneous.color = NULL,
    heterogeneous.lcolor = NULL,
    heterogeneous.cicolor = NULL,
    equiv.color = NULL,
    status.treat.color = NULL,
    status.control.color = NULL,
    status.missing.color = NULL,
    status.removed.color = NULL,
    status.placebo.color = NULL,
    status.carryover.color = NULL,
    status.carryover.rm.color = NULL,
    status.balanced.post.color = NULL,
    status.balanced.pre.color = NULL,
    status.background.color = NULL,
    covariate = NULL,
    covariate.labels = NULL,
    ...) {

  if (!missing(vis)) {
    warning("'vis' is deprecated and will be removed in future versions.", call. = FALSE)
  }
  if (!missing(line.width)) {
    warning("'line.width' is deprecated. For gap/equiv/exit plots, use 'est.lwidth' for estimate lines/points. For factor/counterfactual plots, hline/vline widths are now typically controlled by 'lwidth'.", call. = FALSE)
  }
  if (!missing(line.color)) {
    warning("'line.color' is deprecated. For gap/equiv/exit plots, use 'color' for estimate lines/points. For factor/counterfactual plots, hline/vline colors are now typically controlled by 'lcolor'.", call. = FALSE)
  }
  if (!missing(count)) {
    warning("'count' is deprecated. Use 'show.count'.", call. = FALSE)
    if (is.logical(count) && missing(show.count)) show.count <- count
  }
  if (is.null(preset)) {
    if (is.null(connected)) connected <- FALSE
    if (is.null(ltype)) ltype <- c("solid", "solid")
    if (is.null(gridOff)) gridOff <- FALSE
    if (is.null(color)) color <- "black"
    if (is.null(count.color)) count.color <- "grey70"
    if (is.null(count.alpha)) count.alpha <- 0.4
    if (is.null(count.outline.color)) count.outline.color <- "grey69"
    if (is.null(placebo.color)) placebo.color <- "blue"
    if (is.null(carryover.color)) carryover.color <- "red"
    if (is.null(carryover.rm.color)) carryover.rm.color <- "blue"
    if (is.null(sens.original.color)) sens.original.color <- "darkblue"
    if (is.null(sens.colors)) sens.colors <- c("#218C23", "#FF34B4", "#FF521B", "#2B59C3")

    if (is.null(counterfactual.color)) counterfactual.color <- "steelblue"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#4682B420"
    if (is.null(counterfactual.raw.treated.color)) counterfactual.raw.treated.color <- "#77777750"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "longdash"

    if (is.null(box.control)) box.control <- "skyblue"
    if (is.null(box.treat)) box.treat <- "pink"

    if (is.null(calendar.color)) calendar.color <- "#2C7FB8"
    if (is.null(calendar.cicolor)) calendar.cicolor <- "skyblue"
    if (is.null(calendar.lcolor)) calendar.lcolor <- "red"

    if (is.null(heterogeneous.color)) heterogeneous.color <- "#2C7FB8"
    if (is.null(heterogeneous.cicolor)) heterogeneous.cicolor <- "skyblue"
    if (is.null(heterogeneous.lcolor)) heterogeneous.lcolor <- "red"

    if (is.null(equiv.color)) equiv.color <- "red"

    if (is.null(status.treat.color)) status.treat.color <- "#06266F"
    if (is.null(status.control.color)) status.control.color <- "#B0C4DE"
    if (is.null(status.missing.color)) status.missing.color <- "#FFFFFF"
    if (is.null(status.removed.color)) status.removed.color <- "#A9A9A9"
    if (is.null(status.placebo.color)) status.placebo.color <- "#66C2A5"
    if (is.null(status.carryover.color)) status.carryover.color <- "#E78AC3"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#ffc425"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#00852B"
    if (is.null(status.balanced.pre.color)) status.balanced.pre.color <- "#A5CA18"
    if (is.null(status.background.color)) status.background.color <- "gray90"
  } else if (preset == "vibrant") {
    if (is.null(connected)) connected <- TRUE
    if (is.null(color)) color <- "#054A91"
    if (is.null(count.color)) count.color <- "#E6AF2E"
    if (is.null(count.alpha)) count.alpha <- 1
    if (is.null(count.outline.color)) count.outline.color <- "#E6AF2E"
    if (is.null(lwidth)) lwidth <- c(0.8, 1)
    if (is.null(ltype)) ltype <- c("solid", "dashed")
    if (is.null(placebo.color)) placebo.color <- "#386641"
    if (is.null(carryover.color)) carryover.color <- "#A40E4C"
    if (is.null(carryover.rm.color)) carryover.rm.color <- "#FF5400"
    if (is.null(sens.original.color)) sens.original.color <- "#054A91"
    if (is.null(sens.colors)) sens.colors <- c("#A40E4C", "#FF5400", "#E6AF2E", "#386641", "#ACC3DA")
    if (is.null(counterfactual.color)) counterfactual.color <- "#777777"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#D5E1ED"
    if (is.null(counterfactual.raw.treated.color)) counterfactual.raw.treated.color <- "#77777750"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "dashed"
    if (is.null(box.control)) box.control <- "#ACC3DA"
    if (is.null(box.treat)) box.treat <- "#E1AFC3"
    if (is.null(calendar.color)) calendar.color <- "#ACC3DA"
    if (is.null(calendar.lcolor)) calendar.lcolor <- "#054A91"
    if (is.null(equiv.color)) equiv.color <- "#A40E4C"
    if (is.null(status.treat.color)) status.treat.color <- "#054A91"
    if (is.null(status.control.color)) status.control.color <- "#ACC3DA"
    if (is.null(status.missing.color)) status.missing.color <- "#dddddd"
    if (is.null(status.removed.color)) status.removed.color <- "#D7E8E0"
    if (is.null(status.placebo.color)) status.placebo.color <- "#386641"
    if (is.null(status.carryover.color)) status.carryover.color <- "#A40E4C"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#FF5400"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#E6AF2E"
    if (is.null(status.balanced.pre.color)) status.balanced.pre.color <- "#777777"
    if (is.null(status.background.color)) status.background.color <- "#FFFFFF"
  } else if (preset %in% c("grayscale", "greyscale")) {
    if (is.null(connected)) connected <- FALSE
    if (is.null(color)) color <- "black"
    if (is.null(count.color)) count.color <- "gray80"
    if (is.null(count.alpha)) count.alpha <- 0.5
    if (is.null(count.outline.color)) count.outline.color <- "black"
    if (is.null(lwidth)) lwidth <- c(1, 1)
    if (is.null(ltype)) ltype <- c("solid", "dashed")
    if (is.null(placebo.color)) placebo.color <- "gray40"
    if (is.null(carryover.color)) carryover.color <- "gray70"
    if (is.null(carryover.rm.color)) carryover.rm.color <- "gray40"
    if (is.null(sens.original.color)) sens.original.color <- "gray20"
    if (is.null(sens.colors)) sens.colors <- c("gray80", "gray50", "black")
    if (is.null(counterfactual.color)) counterfactual.color <- "#777777"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#eeeeee"
    if (is.null(counterfactual.raw.treated.color)) counterfactual.raw.treated.color <- "#666666"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "dashed"
    if (is.null(box.control)) box.control <- "gray90"
    if (is.null(box.treat)) box.treat <- "gray50"
    if (is.null(calendar.color)) calendar.color <- "gray80"
    if (is.null(calendar.lcolor)) calendar.lcolor <- "black"
    if (is.null(equiv.color)) equiv.color <- "gray80"
    if (is.null(status.treat.color)) status.treat.color <- "black"
    if (is.null(status.control.color)) status.control.color <- "gray90"
    if (is.null(status.missing.color)) status.missing.color <- "white"
    if (is.null(status.removed.color)) status.removed.color <- "white"
    if (is.null(status.placebo.color)) status.placebo.color <- "gray50"
    if (is.null(status.carryover.color)) status.carryover.color <- "gray50"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#ffffff"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#aaaaaa"
    if (is.null(status.balanced.pre.color)) status.balanced.pre.color <- "#eeeeee"
    if (is.null(status.background.color)) status.background.color <- "#FFFFFF"
  }
  if (is.null(est.lwidth) || is.null(est.pointsize)) {
    default_est_lwidth <- .8
    default_est_pointsize <- 3

    if (!connected) {
      default_est_lwidth <- 0.6
      default_est_pointsize <- 2
    } else {
      if (show.points) {
        default_est_lwidth <- 0.7
        default_est_pointsize <- 1.2
      } else {
        default_est_lwidth <- 1.2
        default_est_pointsize <- 3
      }
    }
    if (is.null(est.lwidth)) est.lwidth <- default_est_lwidth
    if (is.null(est.pointsize)) est.pointsize <- default_est_pointsize
  }


  # come needed variables
  equiv.p <- NULL
  ATT <- ATT2 <- ATT3 <- ATT4 <- NULL
  CI.lower3 <- CI.lower4 <- CI.upper3 <- CI.upper4 <- NULL
  p <- NULL
  outcome <- NULL ## global variable
  labels1 <- labels2 <- labels3 <- NULL
  ATT.OFF <- ATT.ON <- CI.lower <- CI.upper <- NULL
  xmax <- xmin <- ymax <- ymin <- NULL

  # check the class
  if (class(x)[1] != "fect") {
    stop("Not a \"fect\" object.")
  }
  loo.test.out <- test.out <- x$test.out

  if (is.null(balance)) {
    if (!is.null(x$balance.period)) {
      balance <- TRUE
    } else {
      balance <- FALSE
    }
  }

  if (is.null(weight)) {
    if (!is.null(x$W)) {
      weight <- TRUE
    } else {
      weight <- FALSE
    }
  }

  if (is.null(lwidth) == TRUE) {
    lwidth <- 2
    if (theme.bw == TRUE) {
      lwidth <- 1.5
    }
  }
  if (length(as.vector(lwidth)) == 1) {
    lwidth <- rep(lwidth, 2)
  } else if (length(as.vector(lwidth)) != 2) {
    stop("\"lwidth\" must be a numeric vector of length 1 or 2.")
  }



  # check if the input has the loo results
  pequiv <- !is.null(x$pre.est.att) ## if have leave one out pre-treatment results
  if (loo == TRUE && pequiv == FALSE) {
    stop("No leave one out results for pre-treatment periods.")
  } else if (loo == TRUE && pequiv == TRUE) {
    loo <- 1
  } else {
    loo <- 0
  }
  ## y=0 line type
  if (is.null(lcolor) == TRUE) {
    lcolor <- "white"
    if (theme.bw == TRUE) {
      lcolor <- "#AAAAAA70"
    }
  }
  if (length(as.vector(lcolor)) == 1) {
    lcolor <- rep(lcolor, 2)
  } else if (length(as.vector(lcolor)) != 2) {
    stop("\"lcolor\" must be a numeric vector of length 1 or 2.")
  }


  placeboTest <- x$placeboTest
  placebo.period <- x$placebo.period
  carryoverTest <- x$carryoverTest

  carryover.period <- x$carryover.period
  binary <- x$binary
  Ftest <- !is.null(x$pre.test)

  if (!is.null(show.group)) {
    if (length(show.group) > 1) {
      stop("\"show.group\" should contain only one group.\n")
    }
    all.group.name <- names(x$g.level)
    if (!show.group %in% all.group.name) {
      message("The specified group does not exist or its treatment effects cannot be estimated.\n")
      return(0)
    }
  }
  if (!(restrict %in% c("rm", "sm"))) {
    stop("\"restrict\" option misspecified. Must be either \"rm\" or \"sm\".")
  }

  # check the key option type
  if (!is.null(type)) {
    if (type == "ct") {
      type <- "counterfactual"
    }
    if (type == "es") {
      type <- "gap"
    }
    if (type == "hte") {
      type <- "heterogeneous"
    }
    if (!type %in% c("status", "gap", "equiv", "exit", "factors", "loadings", "calendar", "box", "counterfactual", "sens", "sens_es", "cumul", "heterogeneous")) {
      stop("\"type\" option misspecified. Must be one of the following:\"status\",\"gap\",\"equiv\",\"exit\",\"calendar\",\"box\",\"counterfactual\",\"equiv\",\"sens\",\"sens_es\",\"cumul\",\"heterogeneous\".")
    }
    if (type == "exit" && is.null(x$att.off)) {
      stop("No exiting treatment effect to be plotted.")
    }
  } else { # default setting
    if (loo == 1) {
      type <- "equiv"
    } else if (placeboTest == 1) {
      type <- "gap"
    } else if (carryoverTest == 1) {
      type <- "exit"
    } else {
      type <- "gap"
    }
  }
  if (!is.null(x$effect.est.att)) {
    type <- "cumul"
  }

  type.old <- type

  if (is.null(gridOff)) {
    if (type == "status") {
      gridOff <- FALSE
    } else {
      gridOff <- TRUE
    }
  }
  provided_xlim <- xlim
  # add spacing
  if (!is.null(xlim) && !(type %in% c("gap", "equiv", "exit", "sens", "sens_es"))) {
    if (length(xlim) == 2) {
      xlim <- c(xlim[1] - 0.2, xlim[2] + 0.2)
    }
  }


  # factors, loadings, fe
  if (type %in% c("loadings", "factors")) {
    if (type == "loadings") {
      if (!x$method %in% c("gsynth", "ife", "fe")) {
        stop("Can't Visualize the Loadings.\n")
      }
      if (x$r.cv == 0) {
        stop("No factors are included in the model.\n")
      } else {
        ## number of loadings to be plotted
        if (is.null(nfactors) == TRUE) {
          nfactors <- min(x$r.cv, 4)
        } else if (nfactors > x$r.cv) {
          message("Too many factors specified. ")
          nfactors <- min(x$r.cv, 4)
        }
        if (nfactors == 1) {
          message("Loadings for the first factor are shown...\n")
        } else if (nfactors < x$r.cv) {
          message(paste("Loadings for the first", nfactors, "factors are shown...\n"))
        }

        ## title
        if (is.null(main) == TRUE) {
          main <- "Factor Loadings"
        } else if (main == "") {
          main <- NULL
        }

        ## prepare data
        L.hat <- rbind(x$lambda.tr, x$lambda.co)
        Lname <- Llabel <- c()
        r <- x$r.cv
        for (i in 1:r) {
          Lname <- c(Lname, paste("L", i, sep = ""))
          Llabel <- c(Llabel, paste("Factor", i))
        }
        colnames(L.hat) <- Lname
        rownames(L.hat) <- c()

        if (x$force %in% c(1, 3) && include.FE == TRUE) {
          L.hat <- cbind(c(x$alpha.tr, x$alpha.co), L.hat)
          colnames(L.hat) <- c(paste("Factor", 0), Lname)
          rownames(L.hat) <- c()
          nfactors <- nfactors + 1
          Llabel <- c("FE", Llabel)
        }

        data <- cbind.data.frame(L.hat,
          "id" = c(x$tr, x$co),
          "group" = as.factor(c(
            rep("Treated", x$Ntr),
            rep("Control", x$Nco)
          ))
        )

        if (nfactors == 1) {
          p <- ggplot(data, aes(x = .data$group, y = .data$L1, fill = .data$group)) +
            geom_boxplot(alpha = 0.7) +
            coord_flip() +
            guides(fill = FALSE) +
            xlab("") +
            ylab("Factor Loading")
        } else {
          if (x$Ntr >= 5) {
            my_dens <- function(data, mapping, ...) {
              ggplot(data = data, mapping = mapping) +
                geom_density(..., alpha = 0.7, color = NA)
            }
            p <- GGally::ggpairs(data,
              mapping = aes(color = .data$group, fill = .data$group),
              columns = 1:nfactors,
              columnLabels = Llabel[1:nfactors],
              diag = list(continuous = my_dens),
              title = main
            ) +
              theme(plot.title = element_text(hjust = 0.5))
          } else if (x$Ntr > 1) {
            my_dens <- function(data, mapping, ...) {
              ggplot(data = data, mapping = mapping) +
                geom_density(..., fill = "gray", alpha = 0.7, color = "gray50")
            }
            p <- GGally::ggpairs(data,
              mapping = aes(color = .data$group),
              columns = 1:nfactors,
              columnLabels = Llabel[1:nfactors],
              diag = list(continuous = my_dens),
              title = main
            )
          } else {
            my_dens <- function(data, mapping, ...) {
              ggplot(data = data, mapping = mapping) +
                geom_density(..., fill = "gray", alpha = 0.7, color = "gray50")
            }
            p <- GGally::ggpairs(data,
              mapping = aes(color = .data$group),
              columns = 1:nfactors, upper = "blank",
              columnLabels = Llabel[1:nfactors],
              diag = list(continuous = my_dens),
              title = main
            )
          }
        }
        # suppressWarnings(print(p))
        return(p)
      }
    }
    if (type == "factors") {
      if (is.null(legend.pos) == TRUE) {
        legend.pos <- "bottom"
      }

      if (axis.adjust == TRUE) {
        angle <- 45
        x.v <- 1
        x.h <- 1
      } else {
        angle <- 0
        x.v <- 0
        if (type == "missing") {
          x.h <- 0.5
        } else {
          x.h <- 0
        }
      }
      r <- x$r.cv
      if (!x$method %in% c("gsynth", "ife")) {
        stop("Can't Visualize the Loadings.\n")
      }
      if (x$r.cv == 0) {
        stop("No factors are included in the model.\n")
      }

      time <- x$rawtime
      if (!is.numeric(time[1])) {
        time <- 1:x$T
      }
      if (length(provided_xlim) != 0) {
        show <- which(time >= provided_xlim[1] & time <= provided_xlim[2])
        proportion <- 0
      } else {
        show <- 1:length(time)
      }
      nT <- length(show)
      time.label <- x$rawtime[show]
      F.hat <- x$factor

      if (x$force %in% c(2, 3) & include.FE == TRUE) {
        F.hat <- cbind(x$xi, F.hat)
        r <- r + 1
      }

      if (x$r.cv == 0) {
        message("No factors included in the model.\n")
      } else {
        ## axes labels
        if (is.null(xlab) == TRUE) {
          xlab <- x$index[2]
        } else if (xlab == "") {
          xlab <- NULL
        }
        if (is.null(ylab) == TRUE) {
          ylab <- "Estimate"
        } else if (ylab == "") {
          ylab <- NULL
        }
        ## title
        if (is.null(main) == TRUE) {
          main <- "Latent Factors"
        } else if (main == "") {
          main <- NULL
        }
        ## prepare data
        L.co <- x$lambda.co

        if (x$force %in% c(2:3) & include.FE == TRUE) {
          L.co <- cbind(rep(1, dim(L.co)[1]), L.co)
          r.use <- c(0:(r - 1))
        } else {
          r.use <- c(1:r)
        }
        norm <- sqrt(diag(t(L.co) %*% L.co) / (x$N - x$Ntr))
        data <- cbind.data.frame(
          "time" = rep(time[show], r),
          "factor" = c(F.hat[show, ]) * rep(norm, each = nT),
          "group" = as.factor(c(rep(r.use, each = nT)))
        )
        ## theme
        p <- ggplot(data)
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        if (gridOff == TRUE) {
          p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
        p <- p + xlab(xlab) + ylab(ylab) + ggtitle(main) +
          geom_hline(yintercept = 0, colour = lcolor[1], size = lwidth[1], linetype = ltype[1]) +
          theme(
            legend.position = legend.pos,
            axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
            plot.title = element_text(
              size = 20,
              hjust = 0.5,
              face = "bold",
              margin = margin(10, 0, 10, 0)
            )
          )
        ## main plot
        p <- p + geom_line(aes(time, factor,
          colour = .data$group,
          group = .data$group
        ), size = 1.2)


        brew.colors <- c("black", "steelblue", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
        set.colors <- brew.colors[1:r]
        p <- p + scale_colour_manual(values = set.colors)

        ## legend
        p <- p + guides(colour = guide_legend(title = "Factor(s)", ncol = 4))

        if (!is.numeric(time.label)) {
          p <- p +
            scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
        }

        ## ylim
        if (is.null(ylim) == FALSE) {
          p <- p + coord_cartesian(ylim = ylim)
        }
        # suppressWarnings(print(p))
        return(p)
      }
    }
  }


  if (!is.null(show.group)) {
    if (is.null(x$group.output[[show.group]]$att.on) & type != "status") {
      stop(paste0(show.group, " doesn't contain treated units. Can't plot.\n"))
    }
  }

  if (type == "status" && !is.null(show.group)) {
    if (!is.null(id)) {
      stop("\"show.group\" can't be used with \"id\" in the status plot.\n")
    }
  }

  if (type == "gap" | type == "equiv") {
    carryoverTest <- FALSE
  }

  if (type == "exit") {
    placeboTest <- FALSE
  }

  if (!is.null(plot.ci)) {
    if (!plot.ci %in% c("0.9", "0.95", "none")) {
      stop("\"plot.ci\" must be one of \"0.95\", \"0.9\" or \"none\".")
    }
    if (plot.ci %in% c("0.90", "0.95") && is.null(x$est.att)) {
      stop("No uncertainty estimates.")
    }
    if (plot.ci == "0.90" && type %in% c("gap", "exit")) {
      warning("90% CI in gap plots or exiting treatment plots.\n")
    }
    if (plot.ci == "0.95" && type == "equiv") {
      warning("95% CI in equivalence test plots.\n")
    }
  } else { # default settings for plot.ci
    if (is.null(x$est.att)) {
      plot.ci <- "none"
    } else if (type == "equiv") {
      plot.ci <- "0.9"
    } else { # gap plot or exiting plot
      plot.ci <- "0.95"
    }
  }

  if (plot.ci == "0.95") {
    plot.ci <- "95"
  }
  if (plot.ci == "0.9") {
    plot.ci <- "90"
  }

  if (type == "equiv" && plot.ci == "none") {
    stop("No uncertainty estimates. Can't perform equivalence tests.\n")
  }

  if ((placeboTest | carryoverTest) && plot.ci == "none" && type != "status") {
    stop("No uncertainty estimates. Can't perform placebo test or carryover test.\n")
  }

  # if (!connected) {
  #     if (is.null(plot.ci)) {
  #         plot.ci.point <- "both"
  #     } else {
  #         plot.ci.point <- plot.ci
  #     }
  #     plot.ci <- "both"
  # }

  if (is.null(highlight)) {
    if (placeboTest | carryoverTest) {
      highlight <- TRUE
    } else {
      highlight <- FALSE
    }
  }

  if (is.null(bound) == FALSE) {
    if (!bound %in% c("none", "min", "equiv", "both")) {
      stop("\"bound\" option misspecified.")
    }
  } else { # default settings for bound
    if (type == "equiv") {
      bound <- "both"
    } else {
      bound <- "none"
    }
  }

  if (is.null(xlim) == FALSE) {
    if (is.numeric(xlim) == FALSE) {
      stop("Some element in \"xlim\" is not numeric.")
    } else {
      if (length(xlim) != 2) {
        stop("xlim must be of length 2.")
      }
    }
  }
  if (is.null(ylim) == FALSE) {
    if (is.numeric(ylim) == FALSE) {
      stop("Some element in \"ylim\" is not numeric.")
    } else {
      if (length(ylim) != 2) {
        stop("ylim must be of length 2.")
      }
    }
  }

  if (is.null(xlab) == FALSE) {
    if (is.character(xlab) == FALSE) {
      stop("\"xlab\" is not a string.")
    } else {
      xlab <- xlab[1]
    }
  }
  if (is.null(ylab) == FALSE) {
    if (is.character(ylab) == FALSE) {
      stop("\"ylab\" is not a string.")
    } else {
      ylab <- ylab[1]
    }
  }

  if (!is.null(stats)) {
    if (!placeboTest && !carryoverTest) {
      for (i in 1:length(stats)) {
        if (!stats[i] %in% c("none", "F.p", "F.equiv.p", "F.stat", "equiv.p")) {
          stop("Choose \"stats\" from c(\"none\", \"F.stat\", \"F.p\", \"F.equiv.p\", \"equiv.p\").")
        }
      }
    } else if (placeboTest) {
      for (i in 1:length(stats)) {
        if (!stats[i] %in% c("none", "placebo.p", "equiv.p")) {
          stop("Choose \"stats\" from c(\"none\", \"placebo.p\", \"equiv.p\").")
        }
      }
    } else if (carryoverTest) { # carry over test
      for (i in 1:length(stats)) {
        if (!stats[i] %in% c("none", "carryover.p", "equiv.p")) {
          stop("Choose \"stats\" from c(\"none\", \"carryover.p\", \"equiv.p\").")
        }
      }
    }
    if ("none" %in% stats) {
      stats <- "none"
    }
  } else { # default settings for the option stats
    if (type == "gap") {
      if (placeboTest == TRUE) {
        stats <- c("placebo.p", "equiv.p")
      } else {
        stats <- c("none")
      }
    } else if (type == "equiv") {
      stats <- c("F.p", "equiv.p")
      if (placeboTest == TRUE) {
        stats <- c("placebo.p", "equiv.p")
      }
    } else if (type == "exit") {
      if (carryoverTest == TRUE) {
        stats <- c("carryover.p", "equiv.p")
      } else {
        stats <- "none"
      }
    } else {
      stats <- "none"
    }
  }

  if (type == "calendar") {
    stats <- "none"
  }

  # names for all statistics
  if (!("none" %in% stats)) {
    if (is.null(stats.labs) == FALSE) {
      if (length(stats.labs) != length(stats)) {
        stop("\"stats.lab\" should have the same length as \"stats\".")
      }
    } else {
      stats.labs <- rep(NA, length(stats))
      for (i in 1:length(stats)) {
        if (stats[i] == "F.p") {
          stats.labs[i] <- "F test p-value"
        }
        if (stats[i] == "F.equiv.p") {
          stats.labs[i] <- "F equivalence test p-value"
        }
        if (stats[i] == "F.stat") {
          stats.labs[i] <- "F statistics"
        }
        if (stats[i] == "placebo.p") {
          stats.labs[i] <- "Placebo test p-value"
        }
        if (stats[i] == "carryover.p") {
          stats.labs[i] <- "Carryover effect test p-value"
        }
        if (stats[i] == "equiv.p") {
          if (placeboTest) {
            stats.labs[i] <- "Placebo equivalence test p-value"
          } else if (carryoverTest) {
            stats.labs[i] <- "Carryover effect equivalence test p-value"
          } else {
            stats.labs[i] <- "Equivalence test p-value"
          }
        }
      }
    }
  }

  # titles; xlim and ylim
  ytitle <- NULL
  bound.old <- bound
  if (type == "gap") { # classic plots or placebo tests
    maintext <- "Estimated ATT"
    ytitle <- paste("Effect on", x$Y)
    if (placeboTest == 1) {
      maintext <- "Placebo Test"
    }
  } else if (type == "equiv") { # equiv plot is a gap plot (with some options changes)
    if (is.null(x$est.att)) {
      stop("No uncertainty estimates.\n")
    }
    # classic equivalence test:
    if (length(xlim) != 0 && xlim[2] > 0) {
      xlim[2] <- 0.2
    }
    if (loo == 0) {
      maintext <- "Equivalence Test"
      # ytitle <- paste("Residual Average of",x$Y)
      ytitle <- paste("Effect on", x$Y)
    } else {
      # loo equivalence test
      maintext <- "Leave-one-out Equivalence Test"
      ytitle <- paste("Effect on", x$Y)
    }
  } else if (type == "exit") {
    maintext <- "Estimated ATT"
    ytitle <- paste("Effect on", x$Y)
    if (carryoverTest == 1) {
      maintext <- "Carryover Effects"
    }
  } else if (type == "calendar") {
    maintext <- "CATT by Calendar Time"
    ytitle <- paste("Effect on", x$Y)
  } else if (type == "heterogeneous") {
    maintext <- paste("CATT by", covariate)
    ytitle <- paste("Effect on", x$Y)
  } else if (type == "box") {
    maintext <- "Individual Treatment Effects"
    ytitle <- paste("Effect on", x$Y)
  }

  if (is.logical(legendOff) == FALSE & is.numeric(legendOff) == FALSE) {
    stop("\"legendOff\" is not a logical flag.")
  }
  if (is.logical(gridOff) == FALSE & is.numeric(gridOff) == FALSE) {
    stop("\"gridOff\" is not a logical flag.")
  }

  if (is.null(main) == FALSE) {
    if (is.character(main) == FALSE) {
      stop("\"main\" is not a string.")
    } else {
      main <- main[1]
    }
  }
  if (axis.adjust == TRUE) {
    angle <- 45
    x.v <- 1
    x.h <- 1
  } else {
    angle <- 0
    x.v <- 0
    if (type == "status" || type == "heterogeneous" || type == "calendar") {
      x.h <- 0.5
    } else {
      x.h <- 0
    }
  }

  #### font size
  ## title
  if (is.null(cex.main) == FALSE) {
    if (is.numeric(cex.main) == FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
    cex.main <- 16 * cex.main
  } else {
    cex.main <- 16
  }
  ## subtitle
  if (is.null(cex.main.sub) == FALSE) {
    if (is.numeric(cex.main.sub) == FALSE) {
      stop("\"cex.main.sub\" is not numeric.")
    }
    cex.main.sub <- 16 * cex.main.sub
  } else {
    cex.main.sub <- 16
  }
  ## axis label
  if (is.null(cex.lab) == FALSE) {
    if (is.numeric(cex.lab) == FALSE) {
      stop("\"cex.lab\" is not numeric.")
    }
    cex.lab <- 15 * cex.lab
  } else {
    cex.lab <- 15
  }
  ## axis number
  if (is.null(cex.axis) == FALSE) {
    if (is.numeric(cex.axis) == FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }
    cex.axis <- 15 * cex.axis
  } else {
    cex.axis <- 15
  }
  ## legend
  if (is.null(cex.legend) == FALSE) {
    if (is.numeric(cex.legend) == FALSE) {
      stop("\"cex.legend\" is not numeric.")
    }
    cex.legend <- 15 * cex.legend
  } else {
    cex.legend <- 15
  }
  ## text
  if (is.null(cex.text) == FALSE) {
    if (is.numeric(cex.text) == FALSE) {
      stop("\"cex.text\" is not numeric.")
    }
    cex.text <- 5 * cex.text
  } else {
    cex.text <- 5
  }

  ## text label position
  if (!is.null(stats.pos)) {
    if (length(stats.pos) != 2) {
      stop(" \"stats.pos\" must be of length 2. ")
    }
  }

  # key function
  # generate a data frame contains the results for the plots
  # colnames of est.att: ATT S.E. CI.lower CI.upper p.value count CI.lower.90 CI.upper.90
  est.bound <- est.att <- NULL
  if (!is.null(show.group)) {
    target.group <- x$est.group.output[[show.group]]
    info.group <- x$group.output[[show.group]]

    x$att <- info.group$att.on
    x$time <- info.group$time.on
    x$count <- info.group$count.on
    x$time.off <- info.group$time.off
    x$count.off <- info.group$count.off

    x$est.att <- target.group$att.on
    x$att.bound <- target.group$att.on.bound
    x$att.boot <- target.group$att.on.boot

    x$est.att.off <- target.group$att.off
    x$att.off.bound <- target.group$att.off.bound

    x$est.placebo <- target.group$att.placebo
    x$est.carryover <- target.group$att.carryover

    if (loo == 1) {
      loo.group <- x$pre.est.group.output[[show.group]]
      x$pre.est.att <- loo.group$pre.est.att
      x$pre.att.bound <- loo.group$pre.att.bound
      x$pre.att.boot <- loo.group$pre.att.boot
    }

    if (type == "status") {
      NN <- dim(x$obs.missing)[2]
      TT <- dim(x$obs.missing)[1]
      T.name <- rownames(x$obs.missing)
      N.name <- colnames(x$obs.missing)
      use.obs.missing <- x$obs.missing[which(x$G == x$g.level[show.group])]
      use.id <- colnames(x$obs.missing)
      use.index <- apply(x$G, 2, mean)
      use.id <- use.id[which(use.index == x$g.level[show.group])]
      use.obs.missing <- matrix(use.obs.missing, nrow = TT)
      rownames(use.obs.missing) <- T.name
      colnames(use.obs.missing) <- use.id
      x$obs.missing <- use.obs.missing
    }
  }

  use.balance <- FALSE
  if (!is.null(x$balance.period) & balance == TRUE) {
    x$att <- x$balance.att
    x$time <- x$balance.time
    x$count <- x$balance.count
    x$att.avg <- x$balance.att.avg
    x$est.att <- x$est.balance.att
    x$att.bound <- x$balance.att.bound
    x$att.boot <- x$balance.att.boot
    x$est.placebo <- x$est.balance.placebo
    x$obs.missing <- x$obs.missing.balance
    use.balance <- TRUE
  }

  use.weight <- FALSE
  if (!is.null(x$W) & weight == TRUE) {
    x$att <- x$att.on.W
    x$time <- x$time.on.W
    x$count <- x$count.on.W
    x$att.avg <- x$att.avg.W
    x$est.att <- x$est.att.W
    x$att.bound <- x$att.W.bound
    x$att.boot <- x$att.W.boot
    x$est.placebo <- x$est.placebo.W

    x$att.off <- x$att.off.W
    x$time.off <- x$time.off.W
    x$count.off <- x$count.off.W
    x$att.off.bound <- x$att.off.W.bound
    x$est.att.off <- x$est.att.off.W
    x$est.carryover <- x$est.carryover.W

    use.weight <- TRUE
  }


  if (!is.null(x$est.att)) { # have uncertainty estimation
    est.att <- x$est.att
    est.bound <- x$att.bound

    colnames(est.bound) <- c("CI.lower.90", "CI.upper.90") ## 90% ci
    est.att <- cbind(est.att, est.bound)
    pre.est.att <- pre.att.bound <- NULL
    if (loo == 1) { ## replace pre-treatment period with loo results
      pre.est.att <- x$pre.est.att
      pre.att.bound <- x$pre.att.bound
      colnames(pre.att.bound) <- c("CI.lower.90", "CI.upper.90") ## 90% ci
      pre.est.att <- cbind(pre.est.att, pre.att.bound)
      t0 <- t1 <- NULL
      t.s <- t.e <- NULL
      t0 <- rownames(pre.est.att)
      t1 <- rownames(est.att)
      t.s <- which(t1 == t0[1])
      t.e <- which(t1 == t0[length(t0)])
      est.att[t.s:t.e, ] <- pre.est.att
    }
  }

  # est.att.off for the exit plot
  est.bound.off <- est.att.off <- NULL
  if (!is.null(x$est.att.off)) {
    est.att.off <- x$est.att.off
    est.bound.off <- x$att.off.bound
    colnames(est.bound.off) <- c("CI.lower.90", "CI.upper.90") ## 90% ci
    est.att.off <- cbind(est.att.off, est.bound.off)
  }

  ## -------------------------------##
  ## Plotting
  ## -------------------------------##
  align_time_series <- function(Y, T0_counts, D_full, tr_idx_logical_full) {
    # Y: T × N matrix; T0_counts: length-N vector
    T_rows <- nrow(Y)
    N_tr <- ncol(Y)
    # subset out just the treated columns of D_full
    D_tr <- D_full[, tr_idx_logical_full, drop = FALSE]
    if (ncol(D_tr) != N_tr) {
      stop("Number of columns in Y and in D_full[,tr_idx] must match.")
    }

    # compute timeline
    T0_valid <- T0_counts[!is.na(T0_counts) & is.finite(T0_counts)]
    if (length(T0_valid) == 0) {
      return(list(
        timeline = integer(0),
        Y_aug    = matrix(NA_real_, 0, N_tr)
      ))
    }
    min_rel <- 1 - max(T0_valid)
    max_rel <- T_rows - min(T0_valid)
    timeline <- if (min_rel > max_rel) integer(0) else seq(min_rel, max_rel)

    # build aligned matrix
    Y_aug <- matrix(NA_real_, nrow = length(timeline), ncol = N_tr)
    if (length(timeline) > 0) rownames(Y_aug) <- as.character(timeline)

    for (i in seq_len(N_tr)) {
      t0 <- T0_counts[i]
      if (!is.finite(t0)) next
      for (t_abs in seq_len(T_rows)) {
        rel_time <- t_abs - t0
        row_idx <- match(as.character(rel_time), rownames(Y_aug))
        if (is.na(row_idx)) next
        include <- (rel_time <= 0) ||
          (t_abs <= nrow(D_tr) && D_tr[t_abs, i] == 1 && !is.na(D_tr[t_abs, i]))
        if (include) Y_aug[row_idx, i] <- Y[t_abs, i]
      }
    }

    list(timeline = timeline, Y_aug = Y_aug)
  }

  ct.adjust <- function(Y.tr, Y.ct, T0_counts, D_full, tr_idx_logical_full) {
    # align both
    res.tr <- align_time_series(Y.tr, T0_counts, D_full, tr_idx_logical_full)
    res.ct <- align_time_series(Y.ct, T0_counts, D_full, tr_idx_logical_full)

    # mask treated where control missing
    res.tr$Y_aug[is.na(res.ct$Y_aug)] <- NA_real_

    # compute row‐means
    if (length(res.tr$timeline) > 0) {
      Y.tr.bar <- apply(res.tr$Y_aug, 1, mean, na.rm = TRUE)
      Y.ct.bar <- apply(res.ct$Y_aug, 1, mean, na.rm = TRUE)
      Y.tr.bar[is.nan(Y.tr.bar)] <- NA_real_
      Y.ct.bar[is.nan(Y.ct.bar)] <- NA_real_
      Yb <- cbind(Y.tr.bar, Y.ct.bar)
      colnames(Yb) <- c("Y.tr.bar", "Y.ct.bar")
    } else {
      Yb <- matrix(NA_real_, 0, 2,
        dimnames = list(NULL, c("Y.tr.bar", "Y.ct.bar"))
      )
    }

    list(
      timeline = res.tr$timeline,
      Y.tr.aug = res.tr$Y_aug,
      Y.ct.aug = res.ct$Y_aug,
      Yb       = Yb
    )
  }

  if (type %in% c("counterfactual", "ct")) {
    # --- Basic Setup ---
    scaleFUN <- function(x) sprintf("%.f", x)

    if (!raw %in% c("none", "band", "all")) {
      cat("\"raw\" option misspecifed. Reset to \"none\".\n")
      raw <- "none"
    }
    if (axis.adjust == TRUE) {
      angle <- 45
      x.v <- 1
      x.h <- 1
    } else {
      angle <- 0
      x.v <- 0
      x.h <- 0.5
    }

    plot_xlim <- xlim
    if (!is.null(xlim)) {
      if (is.numeric(xlim) && length(xlim) == 2) {
        proportion <- 0 # When xlim is given, don't filter by proportion
      } else {
        warning("Invalid xlim provided. It must be a numeric vector of length 2. Ignoring xlim.")
      }
    }


    # --- Data Extraction and Cleaning ---
    subset <- 1:x$N
    subset.tr <- x$tr
    if (!is.null(id)) {
      subset.tr <- which(x$id %in% as.character(id))
      if (!all(subset.tr %in% x$tr) || length(subset.tr) < 1) {
        stop("One or more of the units given in \"id\" is not treated.\n")
      }
      x$Y.avg <- NULL
      subset <- sort(c(subset.tr, x$co))
    }
    tr <- x$tr
    co <- x$co
    I_orig <- x$I
    II_orig <- x$II
    Y_orig <- x$Y.dat
    Y.ct_orig <- x$Y.ct
    D_orig <- x$D.dat
    rawid_orig <- x$id
    time_orig_vec <- x$rawtime
    TT_dim <- dim(Y_orig)[1]
    if (!is.numeric(time_orig_vec[1])) {
      time_orig_vec <- 1:TT_dim
    }
    I_mat <- as.matrix(I_orig)[, subset]
    colnames(I_mat) <- rawid_orig[subset]
    rownames(I_mat) <- time_orig_vec
    II_mat <- as.matrix(II_orig)[, subset]
    colnames(II_mat) <- rawid_orig[subset]
    rownames(II_mat) <- time_orig_vec
    Y_mat <- as.matrix(Y_orig)[, subset]
    colnames(Y_mat) <- rawid_orig[subset]
    rownames(Y_mat) <- time_orig_vec
    Y.ct_mat <- as.matrix(Y.ct_orig)[, subset]
    colnames(Y.ct_mat) <- rawid_orig[subset]
    rownames(Y.ct_mat) <- time_orig_vec
    D_mat <- as.matrix(D_orig)[, subset]
    colnames(D_mat) <- rawid_orig[subset]
    rownames(D_mat) <- time_orig_vec
    if (!is.null(id)) {
      tr_idx_logical <- (1:ncol(Y_mat)) %in% which(x$id[subset] %in% as.character(id))
    } else {
      tr_idx_logical <- (1:ncol(Y_mat)) %in% tr
    }
    I.tr <- I_mat[, tr_idx_logical, drop = FALSE]
    II.tr <- II_mat[, tr_idx_logical, drop = FALSE]
    D.tr <- D_mat[, tr_idx_logical, drop = FALSE]
    Y.tr <- Y_mat[, tr_idx_logical, drop = FALSE]
    Y.ct <- Y.ct_mat[, tr_idx_logical, drop = FALSE]
    co_idx_logical <- (1:ncol(Y_mat)) %in% co
    Y.co <- Y_mat[, co_idx_logical, drop = FALSE]

    if (!0 %in% I.tr) {
      pre <- as.matrix(D.tr == 0 & II.tr == 1)
    } else {
      pre <- as.matrix(D.tr == 0 & I.tr == 1 & II.tr == 1)
    }
    T0 <- apply(x$D.dat[, subset.tr, drop = FALSE], 2, function(col) {
      first_one <- which(col == 1)[1]
      if (is.na(first_one)) length(col) else first_one - 1
    })
    id.tr_names <- colnames(Y.tr)
    id.co_names <- colnames(Y.co)
    sameT0 <- length(unique(T0)) == 1
    is_case1_scenario <- length(id.tr_names) == 1 || sameT0 == TRUE

    if (is_case1_scenario) {
      num_treated_units_for_plot <- ncol(D.tr)
      num_time_periods_for_plot <- nrow(D.tr)
      if (num_treated_units_for_plot > 0 && num_time_periods_for_plot > 0) {
        for (j_unit_idx in 1:num_treated_units_for_plot) {
          first_treated_period_indices_for_unit <- which(D.tr[, j_unit_idx] == 1)
          if (length(first_treated_period_indices_for_unit) > 0) {
            first_ever_treated_period_row_idx <- min(first_treated_period_indices_for_unit)
            for (t_row_idx in first_ever_treated_period_row_idx:num_time_periods_for_plot) {
              if (is.na(D.tr[t_row_idx, j_unit_idx]) || D.tr[t_row_idx, j_unit_idx] == 0) {
                Y.tr[t_row_idx, j_unit_idx] <- NA
                Y.ct[t_row_idx, j_unit_idx] <- NA
              }
            }
          }
        }
      }
    }
    Y.ct[is.na(Y.tr)] <- NA_real_
    Yb <- cbind(apply(Y.tr, 1, mean, na.rm = TRUE), apply(Y.ct, 1, mean, na.rm = TRUE))
    colnames(Yb) <- c("Tr_Avg", "Ct_Avg")

    if (is.null(shade.post)) {
      shade.post <- TRUE
    } else if (!is.logical(shade.post)) {
      stop("Option \"shade.post\" must be logical (TRUE/FALSE).")
    }
    if (legendOff == TRUE) {
      legend.pos <- "none"
    } else {
      if (is.null(legend.pos)) {
        legend.pos <- "bottom"
      }
    }
    y_data_for_range_calc <- c()
    Y.tr[is.na(Y.ct)] <- NA_real_
    # Case 1: Single Treated Unit or All Treated Units have Same T0 (Absolute Time)
    if (is_case1_scenario) {
      plot_time_abs <- time_orig_vec
      time_bf_abs_val <- NA
      time_step_abs <- if (length(plot_time_abs) > 1) min(diff(sort(unique(plot_time_abs))), na.rm = TRUE) else 1
      if (!is.finite(time_step_abs) || time_step_abs <= 0) time_step_abs <- 1
      if (sameT0) {
        unique_t0_val_count <- unique(T0)[1]
        if (!is.na(unique_t0_val_count) && unique_t0_val_count >= 0 && (unique_t0_val_count + 1) <= length(plot_time_abs)) {
          time_bf_abs_val <- plot_time_abs[unique_t0_val_count + 1]
        } else if (!is.na(unique_t0_val_count) && unique_t0_val_count == length(plot_time_abs)) {
          time_bf_abs_val <- plot_time_abs[length(plot_time_abs)] + time_step_abs
        }
      }
      vline_pos_abs <- if (!is.na(time_bf_abs_val)) time_bf_abs_val - (time_step_abs / 2) else NA

      # Determine which points to show based on provided_xlim (un-padded)
      show_abs <- 1:length(plot_time_abs)
      if (!is.null(provided_xlim)) {
        show_abs_check <- which(plot_time_abs >= provided_xlim[1] & plot_time_abs <= provided_xlim[2])
        if (length(show_abs_check) == 0) {
          warning("No data points in xlim range.")
        } else {
          show_abs <- show_abs_check
        }
      }

      counts_for_filtering_abs <- apply(Y.tr, 1, function(row) sum(!is.na(row)))
      max_count_for_filtering_abs <- max(counts_for_filtering_abs, na.rm = TRUE)
      if (is.finite(max_count_for_filtering_abs) && max_count_for_filtering_abs > 0) {
        indices_from_proportion_abs <- which(counts_for_filtering_abs >= proportion * max_count_for_filtering_abs)
        show_abs <- intersect(show_abs, indices_from_proportion_abs)
        if (length(show_abs) == 0) {
          stop("No data points remain after applying the 'proportion' filter. Try a smaller value for 'proportion'.")
        }
      }

      if (length(show_abs) == 0) {
        stop("Cannot determine time range for absolute time plot.")
      }

      nT_abs <- length(show_abs)
      time_label_abs <- plot_time_abs[show_abs]
      T.b_abs <- integer(0)
      if (nT_abs > 0) {
        if (is.numeric(time_label_abs) && length(time_label_abs) > 1) {
          T.b_breaks_abs <- pretty(time_label_abs)
          T.b_breaks_abs <- T.b_breaks_abs[T.b_breaks_abs >= min(time_label_abs) & T.b_breaks_abs <= max(time_label_abs)]
          if (length(T.b_breaks_abs) > 0) {
            T.b_abs <- sapply(T.b_breaks_abs, function(br) which.min(abs(time_label_abs - br)))
            T.b_abs <- unique(T.b_abs)
          }
        }
        if (length(T.b_abs) == 0) {
          max_labels <- 10
          step <- max(1, floor(length(time_label_abs) / max_labels))
          T.b_abs <- seq(1, length(time_label_abs), by = step)
        }
      }
      xlab_final <- if (is.null(xlab)) x$index[2] else if (xlab == "") NULL else xlab
      ylab_final <- if (is.null(ylab)) x$Yname else if (ylab == "") NULL else ylab
      plot_single_unit_flag <- FALSE
      unit_to_plot_name <- NULL
      plot_single_unit_flag <- length(id.tr_names) == 1
      unit_to_plot_name <- id.tr_names[1]
      if (plot_single_unit_flag) {
        maintext <- paste("Treated and Estimated Counterfactual (", unit_to_plot_name, ")", sep = "")
        unit_col_idx_in_Y.tr <- which(colnames(Y.tr) == unit_to_plot_name)
        if (length(unit_col_idx_in_Y.tr) != 1) stop(paste("Could not find unique column for unit", unit_to_plot_name, "in Y.tr."))
        tr.info_unit <- Y.tr[show_abs, unit_col_idx_in_Y.tr]
        ct.info_unit <- Y.ct[show_abs, unit_col_idx_in_Y.tr] # Subset by show_abs here
        y_data_for_range_calc <- c(tr.info_unit, ct.info_unit) # Already subsetted
        if (raw == "none") {
          if (x$vartype == "parametric" & !is.null(id) & !is.null(x$eff.boot)) {
            if (plot.ci == "95") {
              para.ci <- basic_ci_alpha(
                rowMeans(x$eff.boot[, which(tr %in% subset.tr), ]),
                x$eff.boot[, which(tr %in% subset.tr), ],
                alpha = 0.05
              )
              x$Y.avg <- data.frame(
                period   = x$rawtime,
                lower.tr = Yb[, "Tr_Avg"],
                upper.tr = Yb[, "Tr_Avg"],
                lower.ct = Yb[, "Ct_Avg"] + para.ci[, "upper"],
                upper.ct = Yb[, "Ct_Avg"] + para.ci[, "lower"]
              )
            } else if (plot.ci == "90") {
              para.ci <- basic_ci_alpha(
                rowMeans(x$eff.boot[, which(tr %in% subset.tr), ]),
                x$eff.boot[, which(tr %in% subset.tr), ],
                alpha = 0.1
              )
              x$Y.avg <- data.frame(
                period = x$rawtime,
                lower90.tr = Yb[, "Tr_Avg"],
                upper90.tr = Yb[, "Tr_Avg"],
                lower90.ct = Yb[, "Ct_Avg"] + para.ci[, "upper"],
                upper90.ct = Yb[, "Ct_Avg"] + para.ci[, "lower"]
              )
            } else {
              warning("Invalid plot.ci provided.")
            }
          } else if (x$vartype == "parametric" & raw == "none" & is.null(id)) {
            if (plot.ci == "95") {
              x$Y.avg <- data.frame(
                period   = x$rawtime,
                lower.tr = Yb[, "Tr_Avg"],
                upper.tr = Yb[, "Tr_Avg"],
                lower.ct = Yb[, "Ct_Avg"] - (x$est.att[, "CI.upper"] - x$est.att[, "ATT"]),
                upper.ct = Yb[, "Ct_Avg"] - (x$est.att[, "CI.lower"] - x$est.att[, "ATT"])
              )
            } else if (plot.ci == "90") {
              x$Y.avg <- data.frame(
                period = x$rawtime,
                lower90.tr = Yb[, "Tr_Avg"],
                upper90.tr = Yb[, "Tr_Avg"],
                lower90.ct = Yb[, "Ct_Avg"] - (x$est.att90[, "CI.upper"] - x$est.att90[, "ATT"]),
                upper90.ct = Yb[, "Ct_Avg"] - (x$est.att90[, "CI.lower"] - x$est.att90[, "ATT"])
              )
            } else {
              warning("Invalid plot.ci provided.")
            }
          }
          data_plot_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(tr.info_unit, ct.info_unit),
            type = factor(
              rep(c("tr", "ct"), each = nT_abs),
              levels = c("tr", "ct")
            )
          )
          p <- ggplot()
          if (plot.ci %in% c("90", "95") && !is.null(x$Y.avg)) {
            tr_lo_col <- if (plot.ci == "95") "lower.tr" else "lower90.tr"
            tr_hi_col <- if (plot.ci == "95") "upper.tr" else "upper90.tr"
            cf_lo_col <- if (plot.ci == "95") "lower.ct" else "lower90.ct"
            cf_hi_col <- if (plot.ci == "95") "upper.ct" else "upper90.ct"
            ci_data <- x$Y.avg[x$Y.avg$period %in% plot_time_abs[show_abs], ]

            if (nrow(ci_data) > 0) {
              p <- p +
                geom_ribbon(
                  data = ci_data,
                  aes(x = period, ymin = .data[[cf_lo_col]], ymax = .data[[cf_hi_col]], fill = "ct"),
                  alpha = 0.2, inherit.aes = FALSE
                )
            }
          }
          p <- p +
            geom_line(
              data = data_plot_abs,
              aes(
                x = time, y = outcome,
                colour = type,
                linetype = type,
                linewidth = type
              )
            )
          set.limits <- c("tr", "ct")
          set.labels <- c(
            paste0("Treated (", unit_to_plot_name, ")"),
            paste0("Est. Y(0) (", unit_to_plot_name, ")")
          )
          set.colors <- c(color, counterfactual.color)
          set.linetypes <- c("solid", counterfactual.linetype)
          set.linewidth <- c(est.lwidth, est.lwidth)
          set.fill <- c(color, counterfactual.color) # for the ribbon legend
        } else if (raw == "band") {
          Y.co.quantiles <- t(apply(Y.co[show_abs, , drop = FALSE], 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE))
          main_lines_data_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(tr.info_unit, ct.info_unit), # Use already subsetted data
            type = factor(c(rep("tr", nT_abs), rep("ct", nT_abs)), levels = c("tr", "ct", "co.band"))
          )
          data.band_abs <- data.frame(time = plot_time_abs[show_abs], co5 = Y.co.quantiles[, 1], co95 = Y.co.quantiles[, 2], type = "co.band") # Add type for legend
          y_data_for_range_calc <- c(y_data_for_range_calc, data.band_abs$co5, data.band_abs$co95)
          p <- ggplot() +
            geom_ribbon(data = data.band_abs, aes(x = time, ymin = .data$co5, ymax = .data$co95, fill = type), alpha = 0.15, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
            geom_line(data = main_lines_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          set.limits <- c("tr", "ct", "co.band")
          set.labels <- c(paste0("Treated (", unit_to_plot_name, ")"), paste0("Est. Y(0) (", unit_to_plot_name, ")"), "Controls (5-95% Quantiles)")
          set.colors <- c(color, counterfactual.color, NA)
          set.linetypes <- c("solid", counterfactual.linetype, "blank")
          set.linewidth <- c(est.lwidth, est.lwidth, 0)
          set.fill <- c(NA, NA, counterfactual.color)
        } else if (raw == "all") {
          # Data for Raw Control Lines
          Y.co.subset <- Y.co[show_abs, , drop = FALSE]
          raw_co_data_abs <- NULL
          if (ncol(Y.co.subset) > 0 && nrow(Y.co.subset) > 0) {
            melt_temp_co <- reshape2::melt(Y.co.subset, varnames = c("time_idx_abs", "id_co_idx_abs"), value.name = "outcome")
            if (nrow(melt_temp_co) > 0) {
              raw_co_data_abs <- data.frame(time = as.numeric(melt_temp_co$time_idx_abs), id = as.character(melt_temp_co$id_co_idx_abs), outcome = melt_temp_co$outcome, type = "raw.co")
            }
          }
          # Data for Main Treated Line
          main_tr_line_data_abs <- data.frame(time = plot_time_abs[show_abs], outcome = tr.info_unit, type = "tr", id = "_MAIN_TREATED_LINE_")
          # Data for Main Counterfactual Line
          main_ct_line_data_abs <- data.frame(time = plot_time_abs[show_abs], outcome = ct.info_unit, type = "ct", id = "_MAIN_COUNTERFACTUAL_LINE_")
          y_data_for_range_calc <- c(y_data_for_range_calc, if (!is.null(raw_co_data_abs)) raw_co_data_abs$outcome else NULL)

          p <- ggplot()
          if (!is.null(raw_co_data_abs)) {
            p <- p + geom_line(data = raw_co_data_abs, aes(x = time, y = outcome, group = id, colour = type, linetype = type, linewidth = type))
          }
          p <- p + geom_line(data = main_tr_line_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          p <- p + geom_line(data = main_ct_line_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          set.limits <- c("tr", "ct", "raw.co")
          set.labels <- c(paste0("Treated (", unit_to_plot_name, ")"), paste0("Est. Y(0) (", unit_to_plot_name, ")"), "Controls")
          set.colors <- c(color, counterfactual.color, counterfactual.raw.controls.color)
          set.linetypes <- c("solid", counterfactual.linetype, "solid")
          lw <- c(est.lwidth, est.lwidth / 2)
          set.linewidth <- c(lw[1], lw[1], lw[2])
        }
      } else { # Multiple treated units, all with the same T0
        maintext <- "Treated vs Estimated Counterfactuals"
        Yb_show_abs <- Yb[show_abs, , drop = FALSE]
        y_data_for_range_calc <- c(Yb_show_abs[, 1], Yb_show_abs[, 2])
        if (raw == "none") {
          if (x$vartype == "parametric" & !is.null(id) & !is.null(x$eff.boot)) {
            subset.eff.boot <- sapply(seq_len(dim(x$eff.boot)[3]), function(j) {
              rowMeans(x$eff.boot[, which(tr %in% subset.tr), j], na.rm = TRUE)
            })
            if (plot.ci == "95") {
              para.ci <- basic_ci_alpha(
                rowMeans(subset.eff.boot),
                subset.eff.boot,
                alpha = 0.05
              )
              x$Y.avg <- data.frame(
                period   = x$rawtime,
                lower.tr = Yb[, "Tr_Avg"],
                upper.tr = Yb[, "Tr_Avg"],
                lower.ct = Yb[, "Ct_Avg"] + para.ci[, "upper"],
                upper.ct = Yb[, "Ct_Avg"] + para.ci[, "lower"]
              )
            } else if (plot.ci == "90") {
              para.ci <- basic_ci_alpha(
                rowMeans(subset.eff.boot),
                subset.eff.boot,
                alpha = 0.1
              )
              x$Y.avg <- data.frame(
                period = x$rawtime,
                lower90.tr = Yb[, "Tr_Avg"],
                upper90.tr = Yb[, "Tr_Avg"],
                lower90.ct = Yb[, "Ct_Avg"] + para.ci[, "upper"],
                upper90.ct = Yb[, "Ct_Avg"] + para.ci[, "lower"]
              )
            } else {
              warning("Invalid plot.ci provided.")
            }
          } else if (x$vartype == "parametric" & raw == "none" & is.null(id)) {
            if (plot.ci == "95") {
              x$Y.avg <- data.frame(
                period   = x$rawtime,
                lower.tr = Yb[, "Tr_Avg"],
                upper.tr = Yb[, "Tr_Avg"],
                lower.ct = Yb[, "Ct_Avg"] - (x$est.att[, "CI.upper"] - x$est.att[, "ATT"]),
                upper.ct = Yb[, "Ct_Avg"] - (x$est.att[, "CI.lower"] - x$est.att[, "ATT"])
              )
            } else if (plot.ci == "90") {
              x$Y.avg <- data.frame(
                period = x$rawtime,
                lower90.tr = Yb[, "Tr_Avg"],
                upper90.tr = Yb[, "Tr_Avg"],
                lower90.ct = Yb[, "Ct_Avg"] - (x$est.att90[, "CI.upper"] - x$est.att90[, "ATT"]),
                upper90.ct = Yb[, "Ct_Avg"] - (x$est.att90[, "CI.lower"] - x$est.att90[, "ATT"])
              )
            } else {
              warning("Invalid plot.ci provided.")
            }
          }
          data_plot_abs <- data.frame(time = rep(plot_time_abs[show_abs], 2), outcome = c(Yb_show_abs[, 1], Yb_show_abs[, 2]), type = factor(c(rep("tr", nT_abs), rep("co", nT_abs)), levels = c("tr", "co")))
          p <- ggplot() # Start with empty ggplot for CI
          if (plot.ci %in% c("90", "95") && !is.null(x$Y.avg)) {
            tr_lo_col <- if (plot.ci == "95") "lower.tr" else "lower90.tr"
            tr_hi_col <- if (plot.ci == "95") "upper.tr" else "upper90.tr"
            cf_lo_col <- if (plot.ci == "95") "lower.ct" else "lower90.ct"
            cf_hi_col <- if (plot.ci == "95") "upper.ct" else "upper90.ct"
            required_cols_ci <- c("period", tr_lo_col, tr_hi_col, cf_lo_col, cf_hi_col)
            if (all(required_cols_ci %in% names(x$Y.avg))) {
              ci_data_abs <- x$Y.avg
              ci_data_filtered_abs <- NULL
              # Filter CI data based on the points being shown
              ci_data_filtered_abs <- ci_data_abs[show_abs, ]
              if (!is.null(ci_data_filtered_abs) && nrow(ci_data_filtered_abs) > 0) {
                p <- p + geom_ribbon(data = ci_data_filtered_abs, aes(x = plot_time_abs[show_abs], ymin = .data[[tr_lo_col]], ymax = .data[[tr_hi_col]], fill = "tr"), alpha = 0.2, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE) +
                  geom_ribbon(data = ci_data_filtered_abs, aes(x = plot_time_abs[show_abs], ymin = .data[[cf_lo_col]], ymax = .data[[cf_hi_col]], fill = "co"), alpha = 0.2, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE)
                y_data_for_range_calc <- c(y_data_for_range_calc, ci_data_filtered_abs[[tr_lo_col]], ci_data_filtered_abs[[tr_hi_col]], ci_data_filtered_abs[[cf_lo_col]], ci_data_filtered_abs[[cf_hi_col]])
              } else {
                warning("No CI data for treated/counterfactual averages.")
              }
            } else {
              warning("CI columns not found for treated/counterfactual averages.")
            }
          }
          p <- p + geom_line(data = data_plot_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          set.limits <- c("tr", "co")
          set.labels <- c("Treated Average", "Estimated Y(0) Average")
          set.colors <- c(color, counterfactual.color)
          set.linetypes <- c("solid", counterfactual.linetype)
          set.linewidth <- c(est.lwidth, est.lwidth)
          set.fill <- c(color, counterfactual.color) # fill for CI ribbons
        } else if (raw == "band") {
          Y.tr.quantiles <- t(apply(Y.tr[show_abs, , drop = FALSE], 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE))
          Y.co.quantiles <- t(apply(Y.co[show_abs, , drop = FALSE], 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE))
          main_lines_data_abs <- data.frame(time = rep(plot_time_abs[show_abs], 2), outcome = c(Yb_show_abs[, 1], Yb_show_abs[, 2]), type = factor(c(rep("tr", nT_abs), rep("co", nT_abs)), levels = c("tr", "co", "tr_band", "co_band")))
          data.band_tr_abs <- data.frame(time = plot_time_abs[show_abs], tr5 = Y.tr.quantiles[, 1], tr95 = Y.tr.quantiles[, 2], type = "tr_band")
          data.band_co_abs <- data.frame(time = plot_time_abs[show_abs], co5 = Y.co.quantiles[, 1], co95 = Y.co.quantiles[, 2], type = "co_band")
          y_data_for_range_calc <- c(y_data_for_range_calc, data.band_tr_abs$tr5, data.band_tr_abs$tr95, data.band_co_abs$co5, data.band_co_abs$co95)
          p <- ggplot() +
            geom_ribbon(data = data.band_tr_abs, aes(x = time, ymin = .data$tr5, ymax = .data$tr95, fill = type), alpha = 0.15, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
            geom_ribbon(data = data.band_co_abs, aes(x = time, ymin = .data$co5, ymax = .data$co95, fill = type), alpha = 0.15, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
            geom_line(data = main_lines_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          set.limits <- c("tr", "co", "tr_band", "co_band")
          set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated (5-95% Quantiles)", "Controls (5-95% Quantiles)")
          set.colors <- c(color, counterfactual.color, NA, NA)
          set.linetypes <- c("solid", counterfactual.linetype, "blank", "blank")
          set.linewidth <- c(est.lwidth, est.lwidth, 0, 0)
          set.fill <- c(NA, NA, color, counterfactual.color)
        } else if (raw == "all") {
          Y.tr.subset <- Y.tr[show_abs, , drop = FALSE]
          Y.co.subset <- Y.co[show_abs, , drop = FALSE]
          raw_tr_data_abs <- NULL
          if (ncol(Y.tr.subset) > 0 && nrow(Y.tr.subset) > 0) {
            melt_temp_tr <- reshape2::melt(Y.tr.subset, varnames = c("t", "id"), value.name = "o")
            if (nrow(melt_temp_tr) > 0) raw_tr_data_abs <- data.frame(time = as.numeric(melt_temp_tr$t), id = as.character(melt_temp_tr$id), outcome = melt_temp_tr$o, type = "raw.tr")
          }
          raw_co_data_abs <- NULL
          if (ncol(Y.co.subset) > 0 && nrow(Y.co.subset) > 0) {
            melt_temp_co <- reshape2::melt(Y.co.subset, varnames = c("t", "id"), value.name = "o")
            if (nrow(melt_temp_co) > 0) raw_co_data_abs <- data.frame(time = as.numeric(melt_temp_co$t), id = as.character(melt_temp_co$id), outcome = melt_temp_co$o, type = "raw.co")
          }
          avg_tr_data_abs <- data.frame(time = plot_time_abs[show_abs], outcome = Yb_show_abs[, 1], type = "tr", id = "tr_avg_line")
          avg_co_data_abs <- data.frame(time = plot_time_abs[show_abs], outcome = Yb_show_abs[, 2], type = "co", id = "co_avg_line")
          y_data_for_range_calc <- c(y_data_for_range_calc, if (!is.null(raw_tr_data_abs)) raw_tr_data_abs$outcome else NULL, if (!is.null(raw_co_data_abs)) raw_co_data_abs$outcome else NULL)
          p <- ggplot()
          if (!is.null(raw_tr_data_abs)) {
            p <- p + geom_line(data = raw_tr_data_abs, aes(x = time, y = outcome, group = id, colour = type, linetype = type, linewidth = type))
          }
          if (!is.null(raw_co_data_abs)) {
            p <- p + geom_line(data = raw_co_data_abs, aes(x = time, y = outcome, group = id, colour = type, linetype = type, linewidth = type))
          }
          p <- p + geom_line(data = avg_tr_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          p <- p + geom_line(data = avg_co_data_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
          set.limits <- c("tr", "co", "raw.tr", "raw.co")
          set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated Raw Data", "Controls Raw Data")
          set.colors <- c(color, counterfactual.color, counterfactual.raw.treated.color, counterfactual.raw.controls.color)
          set.linetypes <- c("solid", counterfactual.linetype, "solid", "solid")
          lw <- c(est.lwidth, est.lwidth / 2)
          set.linewidth <- c(lw[1], lw[1], lw[2], lw[2])
        }
      }
      p <- p + xlab(xlab_final) + ylab(ylab_final) + theme(legend.position = legend.pos, axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v), axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab), plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0)))
      if (theme.bw == TRUE) {
        p <- p + theme_bw()
      }
      p <- p + theme(
        legend.position = legend.pos,
        axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
        axis.text = element_text(size = cex.axis), # This sets axis text generally
        axis.title = element_text(size = cex.lab),
        plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
      )
      if (!is.na(vline_pos_abs)) {
        p <- p + geom_vline(xintercept = vline_pos_abs, colour = lcolor[2], linewidth = lwidth[2], linetype = ltype[2])
        if (shade.post == TRUE) {
          p <- p + annotate("rect", xmin = vline_pos_abs, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = .3)
        }
      }
      p <- p + scale_colour_manual(name = NULL, limits = set.limits, values = set.colors, labels = set.labels, na.value = NA, drop = FALSE) +
        scale_linetype_manual(name = NULL, limits = set.limits, values = set.linetypes, labels = set.labels, na.value = "blank", drop = FALSE) +
        scale_linewidth_manual(name = NULL, limits = set.limits, values = set.linewidth, labels = set.labels, na.value = 0, drop = FALSE)
      guide_obj_abs <- guide_legend(title = NULL, ncol = if (length(set.limits) > 2) ceiling(length(set.limits) / 2) else length(set.limits))
      if (exists("set.fill") && !is.null(set.fill) && any(!is.na(set.fill))) { # Check if set.fill has actual values
        p <- p + scale_fill_manual(name = NULL, limits = set.limits, values = set.fill, labels = set.labels, na.value = NA, drop = FALSE)
        p <- p + guides(colour = guide_obj_abs, linetype = guide_obj_abs, linewidth = guide_obj_abs, fill = guide_obj_abs)
      } else {
        p <- p + guides(colour = guide_obj_abs, linetype = guide_obj_abs, linewidth = guide_obj_abs)
      }
      if (length(T.b_abs) > 0) {
        if (!is.numeric(time_label_abs)) {
          p <- p + scale_x_continuous(breaks = plot_time_abs[show_abs][T.b_abs], labels = time_label_abs[T.b_abs])
        } else {
          p <- p + scale_x_continuous(labels = scaleFUN, breaks = plot_time_abs[show_abs][T.b_abs])
        }
      } else {
        p <- p + scale_x_continuous()
      }
      final_main_text <- if (is.null(main)) maintext else if (main == "") NULL else main
      if (!is.null(final_main_text)) p <- p + ggtitle(final_main_text)
      if (show.count == TRUE) {
        counts_values_abs <- NULL
        current_times_abs <- plot_time_abs[show_abs]
        if (plot_single_unit_flag) {
          unit_col_idx_for_count <- which(colnames(Y.tr) == unit_to_plot_name)
          counts_values_abs <- ifelse(!is.na(Y.tr[show_abs, unit_col_idx_for_count, drop = FALSE]), 1, 0)
        } else {
          counts_values_abs <- apply(Y.tr[show_abs, , drop = FALSE], 1, function(row) sum(!is.na(row)))
        }
        counts_for_plot_df_abs <- data.frame(time = current_times_abs, count = counts_values_abs)
        counts_for_plot_df_abs <- counts_for_plot_df_abs[counts_for_plot_df_abs$count > 0 & !is.na(counts_for_plot_df_abs$count), ]
        if (nrow(counts_for_plot_df_abs) > 0) {
          max_count_val_abs <- max(counts_for_plot_df_abs$count, na.rm = TRUE)
          if (max_count_val_abs > 0 && !is.na(max_count_val_abs)) {
            current_plot_yrange <- NULL
            if (!is.null(ylim)) {
              current_plot_yrange <- ylim
            } else {
              if (length(y_data_for_range_calc) == 0 || all(is.na(y_data_for_range_calc))) {
                gb <- ggplot_build(p)
                current_plot_yrange <- gb$layout$panel_scales_y[[1]]$range$range
              } else {
                current_plot_yrange <- range(y_data_for_range_calc, na.rm = TRUE)
              }
            }
            count_bar_space_prop <- 0.20
            count_bar_space_height_abs <- (current_plot_yrange[2] - current_plot_yrange[1]) * count_bar_space_prop
            actual_rect_length_abs <- count_bar_space_height_abs * 0.8
            rect_min_val_abs <- if (!is.null(ylim)) ylim[1] else current_plot_yrange[1] - count_bar_space_height_abs
            counts_for_plot_df_abs$ymin <- rect_min_val_abs
            counts_for_plot_df_abs$ymax <- rect_min_val_abs + (counts_for_plot_df_abs$count / max_count_val_abs) * actual_rect_length_abs
            time_step_for_bars_abs <- if (length(unique(counts_for_plot_df_abs$time)) > 1) min(diff(sort(unique(counts_for_plot_df_abs$time))), na.rm = TRUE) else 1
            if (!is.finite(time_step_for_bars_abs) || time_step_for_bars_abs <= 0) time_step_for_bars_abs <- 1
            bar_width_half_abs <- time_step_for_bars_abs * 0.20
            counts_for_plot_df_abs$xmin <- counts_for_plot_df_abs$time - bar_width_half_abs
            counts_for_plot_df_abs$xmax <- counts_for_plot_df_abs$time + bar_width_half_abs
            max_count_time_pos_abs <- counts_for_plot_df_abs$time[which.max(counts_for_plot_df_abs$count)[1]]
            text_y_pos_abs <- rect_min_val_abs + actual_rect_length_abs + (count_bar_space_height_abs - actual_rect_length_abs) * 0.5
            p <- p + geom_rect(data = counts_for_plot_df_abs, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = count.color, color = count.outline.color, inherit.aes = FALSE, alpha = count.alpha, linewidth = 0.2) + annotate("text", x = max_count_time_pos_abs, y = text_y_pos_abs, label = max_count_val_abs, size = cex.text * 0.8, hjust = 0.5, vjust = 0.5)
          }
        }
      }
      coord_args_abs <- list(clip = "on", expand = TRUE)
      if (!is.null(ylim)) {
        coord_args_abs$ylim <- ylim
      } else if (show.count == TRUE && exists("rect_min_val_abs") && exists("counts_for_plot_df_abs") && nrow(counts_for_plot_df_abs) > 0 && exists("current_plot_yrange") && !is.null(current_plot_yrange)) {
        final_yrange_min <- min(current_plot_yrange[1], rect_min_val_abs, na.rm = TRUE)
        final_yrange_max <- current_plot_yrange[2]
        coord_args_abs$ylim <- c(final_yrange_min, final_yrange_max)
      }
      # Use the padded plot_xlim for the final axis range
      if (!is.null(plot_xlim)) {
        coord_args_abs$xlim <- plot_xlim
      }
      p <- p + do.call(coord_cartesian, coord_args_abs)
      p <- p + theme(legend.position = legend.pos)
    } else { # Case 2: Staggered Adoption (Different T0) -> Relative Time Plot
      xx <- ct.adjust(Y.tr, Y.ct, T0, D_mat, tr_idx_logical)
      event_time_full_series <- xx$timeline
      if (length(event_time_full_series) == 0 || all(is.na(event_time_full_series))) {
        stop("Relative timeline calculation resulted in zero length or all NA.")
      }

      maintext <- "Treated vs Estimated Counterfactuals"
      xlab_final <- if (is.null(xlab)) "Time Since the Treatment's Onset" else if (xlab == "") NULL else xlab
      ylab_final <- if (is.null(ylab)) x$Yname else if (ylab == "") NULL else ylab

      # Determine indices of event_time_full_series to show, based on provided_xlim (un-padded)
      indices_to_show <- 1:length(event_time_full_series)
      if (!is.null(provided_xlim)) {
        xlim_to_check <- provided_xlim
        if (exists("start0", inherits = FALSE) && isTRUE(start0)) {
          xlim_to_check <- provided_xlim + 1
        }
        show_check <- which(event_time_full_series >= xlim_to_check[1] & event_time_full_series <= xlim_to_check[2])
        if (length(show_check) == 0) {
          warning("User xlim for relative time contains no points.")
        } else {
          indices_to_show <- show_check
        }
      }

      counts_for_filtering_rel <- rowSums(!is.na(xx$Y.tr.aug), na.rm = TRUE)
      max_count_for_filtering_rel <- max(counts_for_filtering_rel, na.rm = TRUE)
      if (is.finite(max_count_for_filtering_rel) && max_count_for_filtering_rel > 0) {
        indices_from_proportion_rel <- which(counts_for_filtering_rel >= proportion * max_count_for_filtering_rel)
        indices_to_show <- intersect(indices_to_show, indices_from_proportion_rel)
        if (length(indices_to_show) == 0) {
          stop("No data points remain after applying the 'proportion' filter. Try a smaller value for 'proportion'.")
        }
      }

      event_times_for_data_subset <- event_time_full_series[indices_to_show]
      Yb_data_subset <- xx$Yb[indices_to_show, , drop = FALSE]
      Ytr_aug_data_subset <- xx$Y.tr.aug[indices_to_show, , drop = FALSE]

      time_for_plot_axis <- event_times_for_data_subset
      vline_pos_for_plot_axis <- 0.5
      apply_start0_shift <- exists("start0", inherits = TRUE) && isTRUE(start0)

      if (apply_start0_shift) {
        time_for_plot_axis <- time_for_plot_axis - 1
        vline_pos_for_plot_axis <- -0.5
      }

      nT_on_axis <- length(time_for_plot_axis)
      T_b_axis_ticks_indices <- integer(0)
      if (nT_on_axis > 0) {
        if (is.numeric(time_for_plot_axis) && length(time_for_plot_axis) > 1) {
          pretty_breaks_axis <- pretty(time_for_plot_axis)
          pretty_breaks_axis <- pretty_breaks_axis[pretty_breaks_axis >= min(time_for_plot_axis, na.rm = TRUE) & pretty_breaks_axis <= max(time_for_plot_axis, na.rm = TRUE)]
          if (length(pretty_breaks_axis) > 0) {
            T_b_axis_ticks_indices <- sapply(pretty_breaks_axis, function(br) which.min(abs(time_for_plot_axis - br)))
            T_b_axis_ticks_indices <- unique(T_b_axis_ticks_indices)
          }
        }
        if (length(T_b_axis_ticks_indices) == 0) {
          max_labels <- 10
          step <- max(1, floor(length(time_for_plot_axis) / max_labels))
          T_b_axis_ticks_indices <- seq(1, length(time_for_plot_axis), by = step)
        }
      }

      y_data_for_range_calc <- c(Yb_data_subset[, 1], Yb_data_subset[, 2])

      if (raw == "none") {
        if (x$vartype == "parametric" & !is.null(id) & !is.null(x$eff.boot)) {
          subset.eff.boot <- sapply(seq_len(dim(x$eff.boot)[3]), function(j) {
            rowMeans(align_time_series(
              x$eff.boot[, which(tr %in% subset.tr), j],
              T0,
              D_mat, tr_idx_logical
            )$Y_aug, na.rm = TRUE)
          })
          if (plot.ci == "95") {
            para.ci <- basic_ci_alpha(
              rowMeans(subset.eff.boot),
              subset.eff.boot,
              alpha = 0.05
            )
            x$Y.avg <- data.frame(
              period   = xx$timeline,
              lower.tr = xx$Yb[, "Y.tr.bar"],
              upper.tr = xx$Yb[, "Y.tr.bar"],
              lower.ct = xx$Yb[, "Y.ct.bar"] + para.ci[, "upper"],
              upper.ct = xx$Yb[, "Y.ct.bar"] + para.ci[, "lower"]
            )
          } else if (plot.ci == "90") {
            para.ci <- basic_ci_alpha(
              rowMeans(subset.eff.boot),
              subset.eff.boot,
              alpha = 0.1
            )
            x$Y.avg <- data.frame(
              period = xx$timeline,
              lower90.tr = xx$Yb[, "Y.tr.bar"],
              upper90.tr = xx$Yb[, "Y.tr.bar"],
              lower90.ct = xx$Yb[, "Y.ct.bar"] + para.ci[, "upper"],
              upper90.ct = xx$Yb[, "Y.ct.bar"] + para.ci[, "lower"]
            )
          } else {
            warning("Invalid plot.ci provided.")
          }
        } else if (x$vartype == "parametric" & raw == "none" & is.null(id)) {
          if (plot.ci == "95") {
            x$Y.avg <- data.frame(
              period   = xx$timeline,
              lower.tr = xx$Yb[, "Y.tr.bar"],
              upper.tr = xx$Yb[, "Y.tr.bar"],
              lower.ct = xx$Yb[, "Y.ct.bar"] - (x$est.att[, "CI.upper"] - x$est.att[, "ATT"]),
              upper.ct = xx$Yb[, "Y.ct.bar"] - (x$est.att[, "CI.lower"] - x$est.att[, "ATT"])
            )
          } else if (plot.ci == "90") {
            x$Y.avg <- data.frame(
              period = xx$timeline,
              lower90.tr = xx$Yb[, "Y.tr.bar"],
              upper90.tr = xx$Yb[, "Y.tr.bar"],
              lower90.ct = xx$Yb[, "Y.ct.bar"] - (x$est.att90[, "CI.upper"] - x$est.att90[, "ATT"]),
              upper90.ct = xx$Yb[, "Y.ct.bar"] - (x$est.att90[, "CI.lower"] - x$est.att90[, "ATT"])
            )
          } else {
            warning("Invalid plot.ci provided.")
          }
        }
        data_plot_main <- data.frame(
          time = rep(time_for_plot_axis, 2),
          outcome = c(Yb_data_subset[, 1], Yb_data_subset[, 2]),
          type = factor(c(rep("tr", nT_on_axis), rep("co", nT_on_axis)), levels = c("tr", "co"))
        )
        p <- ggplot()
        if (plot.ci %in% c("90", "95") && !is.null(x$Y.avg)) {
          tr_lo_col <- if (plot.ci == "95") "lower.tr" else "lower90.tr"
          tr_hi_col <- if (plot.ci == "95") "upper.tr" else "upper90.tr"
          cf_lo_col <- if (plot.ci == "95") "lower.ct" else "lower90.ct"
          cf_hi_col <- if (plot.ci == "95") "upper.ct" else "upper90.ct"
          required_cols_ci <- c("period", tr_lo_col, tr_hi_col, cf_lo_col, cf_hi_col)

          if (all(required_cols_ci %in% names(x$Y.avg))) {
            ci_data_for_plot <- x$Y.avg
            ci_data_for_plot$time_axis_period <- ci_data_for_plot$period
            if (apply_start0_shift) {
              ci_data_for_plot$time_axis_period <- ci_data_for_plot$time_axis_period - 1
            }
            ci_data_filtered <- ci_data_for_plot[ci_data_for_plot$period %in% event_times_for_data_subset, ]

            if (!is.null(ci_data_filtered) && nrow(ci_data_filtered) > 0) {
              p <- p + geom_ribbon(data = ci_data_filtered, aes(x = .data$time_axis_period, ymin = .data[[tr_lo_col]], ymax = .data[[tr_hi_col]], fill = "tr"), alpha = 0.2, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE) +
                geom_ribbon(data = ci_data_filtered, aes(x = .data$time_axis_period, ymin = .data[[cf_lo_col]], ymax = .data[[cf_hi_col]], fill = "co"), alpha = 0.2, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE)
              y_data_for_range_calc <- c(y_data_for_range_calc, ci_data_filtered[[tr_lo_col]], ci_data_filtered[[tr_hi_col]], ci_data_filtered[[cf_lo_col]], ci_data_filtered[[cf_hi_col]])
            } else {
              warning("No CI data for treated/counterfactual averages.")
            }
          } else {
            warning("CI columns not found for treated/counterfactual averages.")
          }
        }
        p <- p + geom_line(data = data_plot_main, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
        set.limits <- c("tr", "co")
        set.labels <- c("Treated Average", "Estimated Y(0) Average")
        set.colors <- c(color, counterfactual.color)
        set.linetypes <- c("solid", counterfactual.linetype)
        set.linewidth <- c(est.lwidth, est.lwidth)
        set.fill <- c(color, counterfactual.color)
      } else if (raw == "band") {
        Ytr_aug_quantiles <- t(apply(Ytr_aug_data_subset, 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE))
        main_lines_data_plot <- data.frame(
          time = rep(time_for_plot_axis, 2),
          outcome = c(Yb_data_subset[, 1], Yb_data_subset[, 2]),
          type = factor(c(rep("tr", nT_on_axis), rep("co", nT_on_axis)), levels = c("tr", "co", "tr_band"))
        )
        data_band_plot <- data.frame(
          time = time_for_plot_axis,
          tr5 = Ytr_aug_quantiles[, 1], tr95 = Ytr_aug_quantiles[, 2], type = "tr_band"
        )
        y_data_for_range_calc <- c(y_data_for_range_calc, data_band_plot$tr5, data_band_plot$tr95)
        p <- ggplot() +
          geom_ribbon(data = data_band_plot, aes(x = time, ymin = .data$tr5, ymax = .data$tr95, fill = type), alpha = 0.15, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
          geom_line(data = main_lines_data_plot, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
        set.limits <- c("tr", "co", "tr_band")
        set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated (5-95% Quantiles)")
        set.colors <- c(color, counterfactual.color, NA)
        set.linetypes <- c("solid", counterfactual.linetype, "blank")
        set.linewidth <- c(est.lwidth, est.lwidth, 0)
        set.fill <- c(NA, NA, color)
      } else if (raw == "all") {
        raw_tr_data_plot <- NULL
        if (ncol(Ytr_aug_data_subset) > 0 && nrow(Ytr_aug_data_subset) > 0) {
          melt_temp_tr <- reshape2::melt(Ytr_aug_data_subset, varnames = c("event_t_char", "id"), value.name = "o")
          if (nrow(melt_temp_tr) > 0) {
            plot_time_for_melted_data <- as.numeric(as.character(melt_temp_tr$event_t_char))
            if (apply_start0_shift) {
              plot_time_for_melted_data <- plot_time_for_melted_data - 1
            }
            raw_tr_data_plot <- data.frame(
              time = plot_time_for_melted_data,
              id = as.character(melt_temp_tr$id),
              outcome = melt_temp_tr$o,
              type = "raw.tr"
            )
          }
        }
        avg_tr_data_plot <- data.frame(time = time_for_plot_axis, outcome = Yb_data_subset[, 1], type = "tr", id = "tr_avg_line")
        avg_co_data_plot <- data.frame(time = time_for_plot_axis, outcome = Yb_data_subset[, 2], type = "co", id = "co_avg_line")
        y_data_for_range_calc <- c(y_data_for_range_calc, if (!is.null(raw_tr_data_plot)) raw_tr_data_plot$outcome else NULL)
        p <- ggplot()
        if (!is.null(raw_tr_data_plot)) {
          p <- p + geom_line(data = raw_tr_data_plot, aes(x = time, y = outcome, group = id, colour = type, linetype = type, linewidth = type))
        }
        p <- p + geom_line(data = avg_tr_data_plot, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
        p <- p + geom_line(data = avg_co_data_plot, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))
        set.limits <- c("tr", "co", "raw.tr")
        set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated Raw Data")
        set.colors <- c(color, counterfactual.color, counterfactual.raw.treated.color)
        set.linetypes <- c("solid", counterfactual.linetype, "solid")
        lw <- c(est.lwidth, est.lwidth / 2)
        set.linewidth <- c(lw[1], lw[1], lw[2])
      }

      p <- p + xlab(xlab_final) + ylab(ylab_final) +
        geom_vline(xintercept = vline_pos_for_plot_axis, colour = lcolor[2], linewidth = lwidth[2], linetype = ltype[2]) +
        theme(
          legend.position = legend.pos,
          axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
          axis.text = element_text(size = cex.axis),
          axis.title = element_text(size = cex.lab),
          plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
        )
      if (theme.bw == TRUE) {
        p <- p + theme_bw()
      }
      p <- p + theme(
        legend.position = legend.pos,
        axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
        axis.text = element_text(size = cex.axis),
        axis.title = element_text(size = cex.lab),
        plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
      )
      if (shade.post == TRUE) {
        p <- p + annotate("rect", xmin = vline_pos_for_plot_axis, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = .3)
      }
      p <- p + scale_colour_manual(name = NULL, limits = set.limits, values = set.colors, labels = set.labels, na.value = NA, drop = FALSE) +
        scale_linetype_manual(name = NULL, limits = set.limits, values = set.linetypes, labels = set.labels, na.value = "blank", drop = FALSE) +
        scale_linewidth_manual(name = NULL, limits = set.limits, values = set.linewidth, labels = set.labels, na.value = 0, drop = FALSE)

      guide_obj_staggered <- guide_legend(title = NULL, ncol = if (length(set.limits) > 2) ceiling(length(set.limits) / 2) else length(set.limits))
      if (exists("set.fill") && !is.null(set.fill) && any(!is.na(set.fill))) {
        p <- p + scale_fill_manual(name = NULL, limits = set.limits, values = set.fill, labels = set.labels, na.value = NA, drop = FALSE)
        p <- p + guides(colour = guide_obj_staggered, linetype = guide_obj_staggered, linewidth = guide_obj_staggered, fill = guide_obj_staggered)
      } else {
        p <- p + guides(colour = guide_obj_staggered, linetype = guide_obj_staggered, linewidth = guide_obj_staggered)
      }

      if (length(T_b_axis_ticks_indices) > 0 && nT_on_axis > 0) {
        p <- p + scale_x_continuous(labels = scaleFUN, breaks = time_for_plot_axis[T_b_axis_ticks_indices])
      } else {
        p <- p + scale_x_continuous(labels = scaleFUN)
      }

      final_main_text <- if (is.null(main)) maintext else if (main == "") NULL else main
      if (!is.null(final_main_text)) p <- p + ggtitle(final_main_text)

      if (show.count == TRUE) {
        time_for_count_bars <- time_for_plot_axis
        counts_values <- NULL
        if (nrow(Ytr_aug_data_subset) > 0 && ncol(Ytr_aug_data_subset) > 0) {
          counts_values <- rowSums(!is.na(Ytr_aug_data_subset), na.rm = TRUE)
        } else {
          counts_values <- rep(0, length(time_for_count_bars))
        }

        counts_for_plot_df <- data.frame(time = time_for_count_bars, count = counts_values)
        counts_for_plot_df <- counts_for_plot_df[counts_for_plot_df$count > 0 & !is.na(counts_for_plot_df$count), ]

        if (nrow(counts_for_plot_df) > 0) {
          max_count_val <- max(counts_for_plot_df$count, na.rm = TRUE)
          if (max_count_val > 0 && !is.na(max_count_val)) {
            current_plot_yrange <- NULL
            if (!is.null(ylim)) {
              current_plot_yrange <- ylim
            } else {
              if (length(y_data_for_range_calc) == 0 || all(is.na(y_data_for_range_calc))) {
                gb <- ggplot_build(p)
                current_plot_yrange <- gb$layout$panel_scales_y[[1]]$range$range
              } else {
                current_plot_yrange <- range(y_data_for_range_calc, na.rm = TRUE)
              }
            }
            count_bar_space_prop <- 0.20
            count_bar_space_height <- (current_plot_yrange[2] - current_plot_yrange[1]) * count_bar_space_prop
            actual_rect_length <- count_bar_space_height * 0.8
            rect_min_val <- if (!is.null(ylim)) ylim[1] else current_plot_yrange[1] - count_bar_space_height

            counts_for_plot_df$ymin <- rect_min_val
            counts_for_plot_df$ymax <- rect_min_val + (counts_for_plot_df$count / max_count_val) * actual_rect_length

            time_step_for_bars <- if (length(unique(counts_for_plot_df$time)) > 1) min(diff(sort(unique(counts_for_plot_df$time))), na.rm = TRUE) else 1
            if (!is.finite(time_step_for_bars) || time_step_for_bars <= 0) time_step_for_bars <- 1
            bar_width_half <- time_step_for_bars * 0.20
            counts_for_plot_df$xmin <- counts_for_plot_df$time - bar_width_half
            counts_for_plot_df$xmax <- counts_for_plot_df$time + bar_width_half

            max_count_time_pos <- counts_for_plot_df$time[which.max(counts_for_plot_df$count)[1]]
            text_y_pos <- rect_min_val + actual_rect_length + (count_bar_space_height - actual_rect_length) * 0.5

            p <- p + geom_rect(data = counts_for_plot_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = count.color, inherit.aes = FALSE, alpha = count.alpha, color = count.outline.color, linewidth = 0.2) +
              annotate("text", x = max_count_time_pos, y = text_y_pos, label = max_count_val, size = cex.text * 0.8, hjust = 0.5, vjust = 0.5)
          }
        }
      }

      coord_args_staggered <- list(clip = "on", expand = TRUE)
      if (!is.null(ylim)) {
        coord_args_staggered$ylim <- ylim
      } else if (show.count == TRUE && exists("rect_min_val") && exists("counts_for_plot_df") && nrow(counts_for_plot_df) > 0 && exists("current_plot_yrange") && !is.null(current_plot_yrange)) {
        final_yrange_min <- min(current_plot_yrange[1], rect_min_val, na.rm = TRUE)
        final_yrange_max <- current_plot_yrange[2]
        coord_args_staggered$ylim <- c(final_yrange_min, final_yrange_max)
      }
      # Use the padded plot_xlim for the final axis range
      if (!is.null(plot_xlim)) {
        coord_args_staggered$xlim <- plot_xlim
      }
      p <- p + do.call(coord_cartesian, coord_args_staggered)
      p <- p + theme(legend.position = legend.pos)
    }

    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    if (exists("set.fill")) rm(set.fill)
    return(p)
  }

  show.T0 <- which(x$time == 0)
  if (type == "exit") {
    show.T0 <- which(x$time.off == 0)
  }
  Y <- x$Y.dat
  D <- x$D.dat
  I <- x$I.dat
  Yname <- x$Y
  index <- x$index
  unit.type <- x$unit.type
  obs.missing <- x$obs.missing
  tname <- time <- x$rawtime
  TT <- dim(Y)[1]
  N <- dim(Y)[2]

  if (type == "status") {
    if (is.null(id) == TRUE) {
      id <- colnames(obs.missing)
    }
    m.l <- length(id)
    for (i in 1:m.l) {
      if (!id[i] %in% colnames(obs.missing)) {
        stop("Some specified units are not in the data.")
      }
    }
  } else { ## raw plot
    if (is.null(id) == TRUE) {
      id <- x$id
    }
    m.l <- length(id)
    id.pos <- rep(NA, m.l)
    for (i in 1:m.l) {
      if (!id[i] %in% x$id) {
        stop("Some specified units are not in the data.")
      } else {
        id.pos[i] <- which(x$id == id[i])
      }
    }
    Y <- as.matrix(Y[, id.pos])
    I <- as.matrix(I[, id.pos])
    D <- as.matrix(D[, id.pos])
    unit.type <- unit.type[id.pos]
  }

  ## type of plots
  if (type == "gap" | type == "equiv") {
    if (loo == TRUE & x$loo == TRUE) {
      time.loo <- as.numeric(rownames(x$pre.est.att))
      time.post <- as.numeric(rownames(x$est.att))
      time.post <- time.post[which(time.post > 0)]
      time <- sort(c(time.loo, time.post))
      count.num <- c(x$pre.est.att[, "count.on"], x$est.att[which(as.numeric(rownames(x$est.att)) > 0), "count"])
      best.pos <- 1
      max.count <- max(count.num)
    } else {
      time <- x$time
      count.num <- x$count
      best.pos <- 1
      max.count <- max(count.num)
    }
  } else if (type == "exit") {
    time <- x$time.off
    count.num <- x$count.off
    best.pos <- 0
    max.count <- max(count.num)
  } else {
    if (!is.numeric(time[1])) {
      time <- 1:TT
    }
  }

  ## periods to show
  time.end <- length(time)
  show.time <- 1:time.end
  if (length(xlim) != 0) {
    show.time <- which(time >= xlim[1] & time <= xlim[2])
    proportion <- 0
  }
  if (type %in% c("gap", "equiv", "exit")) {
    if (is.null(proportion) == TRUE) {
      show.c <- 1:time.end
    } else {
      show.c <- which(count.num >= max.count * proportion)
    }
    # which periods to be shown
    show <- intersect(show.c, show.time)

    # maximum number of cases to be shown
    max.count <- max(count.num[show])

    # where on x-axis to show the number
    max.count.pos <- time[intersect(show, which(count.num == max.count))]

    if (length(max.count.pos) > 1) {
      if (best.pos %in% max.count.pos) {
        max.count.pos <- best.pos
      } else if ((1 - best.pos) %in% max.count.pos) {
        max.count.pos <- 1 - best.pos
      } else {
        max.count.pos <- max.count.pos[1]
      }
    }
  } else {
    show <- show.time
  }

  if (length(show) < 2 & type %in% c("gap", "equiv", "exit")) {
    stop("Cannot plot.\n")
  }

  nT <- length(show)
  time.label <- tname[show]
  T.b <- 1:length(show)

  ## legend on/off
  if (legendOff == TRUE) {
    legend.pos <- "none"
  } else {
    if (is.null(legend.pos)) {
      legend.pos <- "bottom"
    }
  }

  ## adjust x-axis
  scaleFUN <- function(x) sprintf("%.f", x) ## integer value at x axis

  # get equivalence p values for two-one-sided-t tests
  tost <- function(coef, se, range) {
    z <- coef / se
    p1 <- 1 - pnorm((-range[1] + coef) / se) # left bound
    p2 <- 1 - pnorm((range[2] - coef) / se) # right bound
    tost.p <- max(p1, p2)
    return(tost.p)
  }


  if (type %in% c("gap", "equiv", "exit")) {
    if (type == "exit") {
      switch.on <- FALSE
    } else {
      switch.on <- TRUE
    }

    ## axes labels
    if (is.null(xlab) == TRUE) {
      if (switch.on == TRUE) {
        xlab <- paste("Time Since the Treatment's Onset")
      } else {
        xlab <- paste("Time Relative to Exiting the Treatment")
      }
    } else if (xlab == "") {
      xlab <- NULL
    }

    if (is.null(ylab) == TRUE) {
      ylab <- ytitle
    } else if (ylab == "") {
      ylab <- NULL
    }




    ## equivalence range

    if (is.null(f.threshold) == TRUE) {
      f.threshold <- x$test.out$f.threshold
      change.f.threshold <- 0
    } else {
      change.f.threshold <- 1
    }

    if (is.null(tost.threshold) == TRUE) {
      if (placeboTest == 1 | x$carryoverTest == 1) {
        tost.threshold <- x$tost.threshold
      } else {
        tost.threshold <- x$test.out$tost.threshold
      }
      change.tost.threshold <- 0
    } else {
      change.tost.threshold <- 1
    }

    if (proportion == x$proportion) {
      change.proportion <- 0
    } else {
      change.proportion <- 1
    }

    if (is.null(pre.periods)) {
      pre.periods <- x$pre.periods
    } else {
      max.count.test <- max(x$count)
      max.pre.periods <- x$time[which(x$count >= max.count.test * proportion & x$time <= 0)]
      pre.periods <- intersect(pre.periods[1]:pre.periods[length(pre.periods)], max.pre.periods)
    }

    if (length(pre.periods) != length(x$pre.periods)) {
      change.pre.periods <- 1
    } else if (all(pre.periods == x$pre.periods)) {
      change.pre.periods <- 0
    } else {
      change.pre.periods <- 1
    }

    ## bound
    data2 <- NULL
    if (bound != "none") {
      if (is.null(x$est.att)) {
        message("No uncertainty estimates.\n")
        bound <- "none"
      }
    }
    if (bound != "none" || "equiv.p" %in% stats) {
      time0 <- NULL
      if (switch.on == TRUE) {
        if (sum(time[show] <= 0) == 0) {
          message("No pretreatment periods are to be plotted.\n")
          time0 <- 1:length(time[show])
        } else {
          time0 <- which(time[show] <= 0)
        }
        att.sub <- as.matrix(est.att[show, c("CI.lower.90", "CI.upper.90")])
        minBound <- max(abs(att.sub[time0, c("CI.lower.90", "CI.upper.90")]), na.rm = TRUE)
      } else {
        if (sum(time[show] > 0) == 0) {
          message("No non-treatment periods are to be plotted.\n")
          time0 <- 1:length(time[show])
        } else {
          time0 <- which(time[show] >= 1)
        }
        att.sub <- as.matrix(x$att.off.bound[show, ])
        minBound <- max(abs(att.sub[time0, c("CI.lower", "CI.upper")]), na.rm = TRUE)
      }
      # Min range
      minbound <- c(-minBound, minBound)
      equiv.range <- c(-1 * tost.threshold, tost.threshold)

      # bound period
      bound.time <- time[show]
      if (switch.on == TRUE) {
        bound.time <- bound.time[which(bound.time <= 0)]
      } else {
        bound.time <- bound.time[which(bound.time >= 1)]
      }

      ## add legend for 95\% CI
      set.limits <- "ci"
      if (is.null(legend.labs) == TRUE) {
        if (plot.ci == "90") {
          set.labels <- "Residual Average (w/ 90% CI)"
        } else if (plot.ci == "95") {
          set.labels <- "ATT (w/ 95% CI)"
        } else {
          set.labels <- "ATT"
        }
      } else {
        set.labels <- legend.labs
      }
      set.colors <- "#000000FF"
      set.linetypes <- "solid"
      set.size <- 1

      ## create a dataset for bound

      if (bound.old == "equiv") {
        use2 <- c(rep(equiv.range, each = length(bound.time)))
        data2 <- cbind.data.frame(bound = use2)
        names(data2) <- "bound"
        data2$time <- rep(bound.time, 2)
        data2$type <- rep(c("equiv"), 2 * length(bound.time))
        data2$id <- rep(1:2, each = length(bound.time))

        set.limits <- c(set.limits, "equiv")
        if (is.null(legend.labs) == TRUE) {
          set.labels <- c(set.labels, "equiv. Bound")
        } else {
          set.labels <- legend.labs
        }
        set.colors <- c(set.colors, equiv.color)
        set.linetypes <- c(set.linetypes, "dashed")
        set.size <- c(set.size, 0.7)
      } else if (bound.old == "min") {
        data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time))))
        names(data2) <- "bound"
        data2$time <- rep(bound.time, 2)
        data2$type <- rep(c("min"), 2 * length(bound.time))
        data2$id <- rep(1:2, each = length(bound.time))
        set.limits <- c(set.limits, "min")
        if (is.null(legend.labs) == TRUE) {
          set.labels <- c(set.labels, "Min. Range")
        } else {
          set.labels <- legend.labs
        }
        set.colors <- c(set.colors, "gray50")
        set.linetypes <- c(set.linetypes, "dashed")
        set.size <- c(set.size, 0.7)
      } else if (bound.old == "both") {
        data2 <- cbind.data.frame(c(rep(minbound, each = length(bound.time)), rep(equiv.range, each = length(bound.time))))
        names(data2) <- "bound"
        data2$time <- rep(bound.time, 4)
        data2$type <- rep(c("min", "equiv"), each = 2 * length(bound.time))
        data2$id <- rep(1:4, each = length(bound.time))
        set.limits <- c(set.limits, "min", "equiv")
        if (is.null(legend.labs) == TRUE) {
          set.labels <- c(set.labels, "Min. Range", "Equiv. Range")
        } else {
          set.labels <- legend.labs
        }
        set.colors <- c(set.colors, "gray50", equiv.color)
        set.linetypes <- c(set.linetypes, "dashed", "dashed")
        set.size <- c(set.size, 0.7, 0.7)
      }
    }

    CI <- NULL
    if (switch.on == TRUE) {
      if (is.null(x$est.att) == TRUE) {
        CI <- FALSE
      } else {
        CI <- TRUE
      }
    } else if (switch.on == FALSE) {
      if (is.null(x$est.att.off) == TRUE) {
        CI <- FALSE
      } else {
        CI <- TRUE
      }
    }

    if (plot.ci == "none") {
      CI <- FALSE
    }
    ## data frame for main estimates
    if (switch.on == TRUE) {
      ## switch-on effect
      if (CI == FALSE) {
        data <- cbind.data.frame(time, ATT = x$att, count = count.num)[show, ]
      } else {
        tb <- est.att
        data <- cbind.data.frame(time, tb)[show, ]
        colnames(data)[2] <- "ATT" # rename 2nd column to "ATT"

        if (plot.ci %in% c("90", "95")) {
          # Optionally set the names of the CI columns:
          ci.name <- if (plot.ci == "95") {
            c("CI.lower", "CI.upper")
          } else {
            c("CI.lower.90", "CI.upper.90")
          }
        }
      }
    } else {
      ## exit treatment plot
      if (CI == FALSE) {
        data <- cbind.data.frame(time, ATT = x$att.off, count = count.num)[show, ]
      } else {
        tb <- est.att.off
        data <- cbind.data.frame(time, tb)[show, ]

        # second column is the estimated ATT, rename it
        colnames(data)[2] <- "ATT"

        # rename 7th column to "count"
        colnames(data)[7] <- "count"

        if (plot.ci %in% c("90", "95")) {
          # Optionally set the names of the CI columns:
          ci.name <- if (plot.ci == "95") {
            c("CI.lower", "CI.upper")
          } else {
            c("CI.lower.90", "CI.upper.90")
          }
        }
      }
    }


    # height of the histogram
    if (CI == FALSE) {
      message("Uncertainty estimates not available.\n")
      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(data[, "ATT"], na.rm = TRUE) - min(data[, "ATT"], na.rm = TRUE)) / 2
        rect.min <- min(data[, "ATT"], na.rm = TRUE) - rect.length
      }
    } else {
      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(data[, "CI.upper"], na.rm = TRUE) - min(data[, "CI.lower"], na.rm = TRUE)) / 2
        rect.min <- min(data[, "CI.lower"], na.rm = TRUE) - rect.length
      }
    }

    ## plotting
    ## line

    if (start0 == TRUE) {
      data$time <- data$time - 1
      data2$time <- data2$time - 1
      if (!is.null(placebo.period)) {
        placebo.period <- placebo.period - 1
      }
      if (!is.null(carryover.period)) {
        carryover.period <- carryover.period - 1
      }
    }

    ## theme
    if (theme.bw == TRUE) {
      p <- p + theme_bw()
    }


    ## add ATT point estimates
    classic <- 0
    if (highlight == FALSE) {
      classic <- 1
    }
    if (placeboTest == TRUE && length(placebo.period) == 1 && plot.ci %in% c("90", "95")) {
      classic <- 1
    }
    if (carryoverTest == TRUE && length(carryover.period) == 1 && plot.ci %in% c("90", "95")) {
      classic <- 1
    }
    if (carryoverTest == FALSE && placeboTest == FALSE) {
      classic <- 1
    }


    # Initialize containers for stats
    stats_values <- numeric(0) # Will store numeric statistic values
    stats_labels <- character(0) # Will store the labels for the statistics
    p.label.lines <- character(0) # To reconstruct p.label string if needed elsewhere

    # --- Cache for diagtest results ---
    # These will store the output of diagtest to avoid redundant calls.
    # We need separate caches because diagtest might be called with f.threshold or tost.threshold,
    # and also in a LOO (Leave-One-Out) context.

    # For standard (non-LOO like) calls
    cached_test_out_f <- NULL # Results when f.threshold is primary
    cached_test_out_tost <- NULL # Results when tost.threshold is primary

    # For LOO calls
    cached_loo_test_out_f <- NULL # LOO results when f.threshold is primary
    cached_loo_test_out_tost <- NULL # LOO results when tost.threshold is primary

    # Determine the effective LOO setting for the current context
    # This flag indicates if the overall scenario is LOO for 'equiv' or 'gap' types.
    is_current_scenario_loo <- (x$loo == TRUE && loo == TRUE && type %in% c("equiv", "gap"))

    # Loop through each requested statistic
    for (i in 1:length(stats)) {
      stat_name_from_user <- stats[i]
      stat_label_text_from_user <- stats.labs[i] # Assuming stats.labs corresponds by index
      value_to_store <- NULL

      # --- Path 1: Non-LOO calculations or types other than 'equiv'/'gap' LOO ---
      if (!is_current_scenario_loo) {
        original_x_loo_state <- x$loo # Preserve original x$loo
        x$loo <- FALSE # Ensure diagtest runs in non-LOO mode

        # Category 1.1: F-statistics (F.p, F.stat, F.equiv.p)
        if (stat_name_from_user %in% c("F.p", "F.stat")) {
          # Recalculate if basic parameters changed or if not cached
          if (change.proportion || change.pre.periods || !is.null(show.group) || use.balance || is.null(cached_test_out_f)) {
            cached_test_out_f <- diagtest(x,
              proportion = proportion,
              pre.periods = pre.periods,
              f.threshold = f.threshold
            )
          }
          if (stat_name_from_user == "F.p") value_to_store <- cached_test_out_f$f.p
          if (stat_name_from_user == "F.stat") value_to_store <- cached_test_out_f$f.stat
        } else if (stat_name_from_user == "F.equiv.p") {
          # Recalculate also if f.threshold changed
          if (change.f.threshold || change.proportion || change.pre.periods || !is.null(show.group) || use.balance || is.null(cached_test_out_f)) {
            cached_test_out_f <- diagtest(x,
              proportion = proportion,
              pre.periods = pre.periods,
              f.threshold = f.threshold
            )
          }
          value_to_store <- cached_test_out_f$f.equiv.p
        }

        # Category 1.2: Statistics related to TOST, Placebo, Carryover (equiv.p, placebo.p, carryover.p)
        # These typically depend on tost.threshold
        else if (stat_name_from_user == "equiv.p") {
          if (type %in% c("equiv", "gap") && placeboTest == 0) { # Standard TOST equivalence
            if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || use.balance || is.null(cached_test_out_tost)) {
              cached_test_out_tost <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            }
            value_to_store <- cached_test_out_tost$tost.equiv.p
          } else if ((type %in% c("equiv", "gap") || type == "gap") && placeboTest == 1) { # Placebo equivalence p-value
            # Note: original logic for placebo.equiv.p did not include 'use.balance' in its specific condition block
            if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_test_out_tost)) {
              cached_test_out_tost <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            }
            value_to_store <- cached_test_out_tost$placebo.equiv.p
          } else if (type == "exit" && carryoverTest == TRUE) { # Carryover equivalence p-value
            if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_test_out_tost)) {
              cached_test_out_tost <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            }
            value_to_store <- cached_test_out_tost$carryover.equiv.p
          }
        } else if (stat_name_from_user == "placebo.p" && type == "gap" && placeboTest == TRUE) { # Standard placebo p-value
          if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_test_out_tost)) {
            cached_test_out_tost <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
          }
          value_to_store <- cached_test_out_tost$placebo.p
        } else if (stat_name_from_user == "carryover.p" && type == "exit" && carryoverTest == TRUE) { # Standard carryover p-value
          if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_test_out_tost)) {
            cached_test_out_tost <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
          }
          value_to_store <- cached_test_out_tost$carryover.p
        }
        x$loo <- original_x_loo_state # Restore x$loo state
      }

      # --- Path 2: LOO calculations (only for type 'equiv' or 'gap') ---
      else { # This means is_current_scenario_loo is TRUE
        original_x_loo_state <- x$loo # Preserve original x$loo
        x$loo <- TRUE # Ensure diagtest runs in LOO mode

        # Category 2.1: LOO F-statistics
        if (stat_name_from_user %in% c("F.p", "F.stat")) {
          # Recalculate for LOO if basic LOO parameters changed or not cached
          # Note: Original LOO F.p/F.stat fell back to non-LOO x$test.out if no specific changes.
          # This version consistently uses a LOO-mode diagtest call.
          if (change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_loo_test_out_f)) {
            cached_loo_test_out_f <- diagtest(x,
              proportion = proportion,
              pre.periods = pre.periods,
              f.threshold = f.threshold
            )
          }
          # Assuming diagtest output fields (f.p, f.stat) are the LOO versions when x$loo = TRUE
          if (stat_name_from_user == "F.p") value_to_store <- cached_loo_test_out_f$f.p
          if (stat_name_from_user == "F.stat") value_to_store <- cached_loo_test_out_f$f.stat
        } else if (stat_name_from_user == "F.equiv.p") { # LOO F-equivalence p-value
          if (change.f.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_loo_test_out_f)) {
            cached_loo_test_out_f <- diagtest(x,
              proportion = proportion,
              pre.periods = pre.periods,
              f.threshold = f.threshold
            )
          }
          # Assuming diagtest output field f.equiv.p is the LOO version when x$loo = TRUE
          value_to_store <- cached_loo_test_out_f$f.equiv.p
        }

        # Category 2.2: LOO TOST equivalence (equiv.p)
        else if (stat_name_from_user == "equiv.p") {
          if (change.tost.threshold || change.proportion || change.pre.periods || !is.null(show.group) || is.null(cached_loo_test_out_tost)) {
            cached_loo_test_out_tost <- diagtest(x,
              proportion = proportion,
              pre.periods = pre.periods,
              tost.threshold = tost.threshold
            )
          }
          # Assuming diagtest output field tost.equiv.p is the LOO version when x$loo = TRUE
          value_to_store <- cached_loo_test_out_tost$tost.equiv.p
        }
        x$loo <- original_x_loo_state # Restore x$loo state
      }

      # If a value was determined for the current statistic, store it
      if (!is.null(value_to_store)) {
        stats_labels <- c(stats_labels, stat_label_text_from_user)
        stats_values <- c(stats_values, as.numeric(value_to_store)) # Ensure it's stored as numeric
        p.label.lines <- c(p.label.lines, paste0(stat_label_text_from_user, ": ", sprintf("%.3f", value_to_store)))
      }
    } # End of loop over stats

    # Reconstruct the p.label string if it's used elsewhere (e.g., for a single annotation object)
    p.label <- if (length(p.label.lines) > 0) paste(p.label.lines, collapse = "\n") else NULL

    if (type == "equiv" && is.null(ylim)) {
      ylim <- c(-max(abs(data2$bound)) * 1.4, max(abs(data2$bound)) * 1.4)
    }
    ## point estimates
    if (classic == 1) {
      #
      # --- REGULAR EVENT-STUDY PLOT (no highlights) ---
      #

      # 1) Build a data frame that esplot() expects
      data_es <- data.frame(
        Period = data$time,
        ATT = data$ATT,
        `S.E.` = if (CI) data$S.E. else NA, # if you have standard errors
        CI.lower = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else {
          NA
        },
        CI.upper = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else {
          NA
        },
        count = if (!is.null(data$count)) data$count else NA
      )
      # 2) Call esplot()
      p <- esplot(
        data_es,
        Period = "Period",
        Estimate = "ATT",
        SE = "S.E.",
        CI.lower = "CI.lower",
        CI.upper = "CI.upper",
        Count = "count",
        show.count = show.count,
        show.points = show.points,
        ci.outline = ci.outline,
        connected = connected,
        color = color,
        count.color = count.color,
        count.alpha = count.alpha,
        count.outline.color = count.outline.color,
        xlab = xlab,
        ylab = ylab,
        main = main,
        xlim = xlim,
        ylim = ylim,
        gridOff = gridOff,
        start0 = start0,
        cex.main = cex.main / 16,
        cex.axis = cex.axis / 15,
        cex.lab = cex.lab / 15,
        cex.text = cex.text / 5,
        proportion = proportion,
        est.lwidth = est.lwidth,
        est.pointsize = est.pointsize,
        fill.gap = FALSE,
        lcolor = lcolor,
        lwidth = lwidth,
        ltype = ltype,
        axis.adjust = axis.adjust,
        stats = if (exists("stats_values")) stats_values else NULL,
        stats.labs = if (exists("stats_labels")) stats_labels else NULL,
        stats.pos = stats.pos,
        theme.bw = theme.bw,
        only.pre = type == "equiv"
      )
    } else if (classic == 0 && switch.on == TRUE) {
      #
      # --- PLACEBO TEST PLOT ---
      #

      data_es <- data.frame(
        Period = data$time,
        ATT = data$ATT,
        `S.E.` = if (CI) data$S.E. else NA,
        CI.lower = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else {
          NA
        },
        CI.upper = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else {
          NA
        },
        count = if (!is.null(data$count)) data$count else NA
      )

      # highlight the placebo interval [placebo.period[1], placebo.period[2]]
      placebo_seq <- seq(placebo.period[1], placebo.period[2])
      n_placebo <- length(placebo_seq)
      p <- esplot(
        data_es,
        Period = "Period",
        Estimate = "ATT",
        SE = "S.E.",
        CI.lower = "CI.lower",
        CI.upper = "CI.upper",
        Count = "count",
        show.count = show.count,
        connected = connected,
        color = color,
        count.color = count.color,
        count.alpha = count.alpha,
        count.outline.color = count.outline.color,
        show.points = show.points,
        ci.outline = ci.outline,
        highlight.periods = placebo_seq,
        highlight.colors = rep(placebo.color, n_placebo),
        xlab = xlab,
        ylab = ylab,
        main = main,
        xlim = xlim,
        ylim = ylim,
        gridOff = gridOff,
        start0 = start0,
        proportion = proportion,
        fill.gap = FALSE,
        lcolor = lcolor,
        lwidth = lwidth,
        ltype = ltype,
        cex.main = cex.main / 16,
        cex.axis = cex.axis / 15,
        cex.lab = cex.lab / 15,
        cex.text = cex.text / 5,
        axis.adjust = axis.adjust,
        est.lwidth = est.lwidth,
        est.pointsize = est.pointsize,
        stats = if (exists("stats_values")) stats_values else NULL,
        stats.labs = if (exists("stats_labels")) stats_labels else NULL,
        theme.bw = theme.bw,
        stats.pos = stats.pos
      )
    } else if (classic == 0 && switch.on == FALSE) {
      #
      # --- CARRYOVER TEST OR EXITING TREATMENT ---
      #
      if (is.null(x$est.carry.att)) {
        placebo_seq <- c()
        n_placebo <- 0
        shift <- 0
      } else {
        placebo_seq <- seq(carryover.period[1], carryover.period[1] - 1 + dim(x$est.carry.att)[1])
        n_placebo <- length(placebo_seq)
        shift <- dim(x$est.carry.att)[1]
      }

      data_es <- data.frame(
        Period = data$time + shift,
        ATT = data$ATT,
        `S.E.` = if (CI) data$S.E. else NA,
        CI.lower = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else {
          NA
        },
        CI.upper = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else {
          NA
        },
        count = if (!is.null(data$count)) data$count else NA
      )


      carry_seq <- seq(carryover.period[1] + shift, carryover.period[2] + shift)
      n_carry <- length(carry_seq)

      p <- esplot(
        data_es,
        Period = "Period",
        Estimate = "ATT",
        SE = "S.E.",
        CI.lower = "CI.lower",
        CI.upper = "CI.upper",
        Count = "count",
        show.count = show.count,
        show.points = show.points,
        ci.outline = ci.outline,
        connected = connected,
        color = color,
        count.color = count.color,
        count.alpha = count.alpha,
        count.outline.color = count.outline.color,
        highlight.periods = c(placebo_seq, carry_seq),
        highlight.colors = c(rep(carryover.rm.color, n_placebo), rep(carryover.color, n_carry)),
        xlab = xlab,
        ylab = ylab,
        main = main,
        xlim = xlim,
        ylim = ylim,
        gridOff = gridOff,
        fill.gap = FALSE,
        lcolor = lcolor,
        lwidth = lwidth,
        ltype = ltype,
        cex.main = cex.main / 16,
        cex.axis = cex.axis / 15,
        cex.lab = cex.lab / 15,
        cex.text = cex.text / 5,
        axis.adjust = axis.adjust,
        start0 = start0,
        proportion = proportion,
        ## newly added options:
        est.lwidth = est.lwidth,
        est.pointsize = est.pointsize,
        stats = if (exists("stats_values")) stats_values else NULL,
        stats.labs = if (exists("stats_labels")) stats_labels else NULL,
        theme.bw = theme.bw,
        stats.pos = stats.pos
      )
    }

    # plot bound
    if (bound.old != "none") { ## with bounds
      p <- p + geom_line(data = data2, aes(time, bound, colour = type, linetype = type, size = type, group = id))
      ## legends for bounds
      if (is.null(legend.nrow) == TRUE) {
        legend.nrow <- ifelse(length(set.limits) <= 3, 1, 2)
      }
      p <- p + scale_colour_manual(limits = set.limits, labels = set.labels, values = set.colors) +
        scale_size_manual(limits = set.limits, labels = set.labels, values = set.size) +
        scale_linetype_manual(limits = set.limits, labels = set.labels, values = set.linetypes) +
        guides(
          linetype = guide_legend(title = NULL, nrow = legend.nrow), colour = guide_legend(title = NULL, nrow = legend.nrow),
          size = guide_legend(title = NULL, nrow = legend.nrow)
        )

      p <- p + theme(legend.position = "bottom")
    }
  }

  if (type == "calendar") {
    CI <- NULL
    if (is.null(x$est.eff.calendar) == TRUE) {
      CI <- FALSE
    } else {
      CI <- TRUE
    }
    if (!is.null(provided_xlim)) {
      x$est.eff.calendar <- x$est.eff.calendar[which(rownames(x$est.eff.calendar) >= min(provided_xlim) & rownames(x$est.eff.calendar) <= max(provided_xlim)), ]
      x$est.eff.calendar.fit <- x$est.eff.calendar.fit[which(rownames(x$est.eff.calendar.fit) >= min(provided_xlim) & rownames(x$est.eff.calendar.fit) <= max(provided_xlim)), ]
    }
    if (plot.ci == "none") {
      CI <- FALSE
    }
    ## axes labels
    if (is.null(xlab) == TRUE) {
      xlab <- "Calendar Time"
    } else if (xlab == "") {
      xlab <- NULL
    }

    if (is.null(ylab) == TRUE) {
      ylab <- ytitle
    } else if (ylab == "") {
      ylab <- NULL
    }


    if (CI == FALSE) {
      message("Uncertainty estimates not available.\n")
      data.1 <- x$eff.calendar
      data.2 <- x$eff.calendar.fit
      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(c(data.1, data.2), na.rm = TRUE) - min(c(data.1, data.2), na.rm = TRUE)) / 2
        rect.min <- min(c(data.1, data.2), na.rm = TRUE) - rect.length
      }
      d1 <- data.1 <- as.matrix(x$eff.calendar[which(!is.na(x$eff.calendar[, 1])), ])
      d2 <- data.2 <- as.matrix(x$eff.calendar.fit[which(!is.na(x$eff.calendar.fit[, 1])), ])
      if (dim(d1)[2] == 1) {
        d1 <- data.1 <- t(d1)
        rownames(d1) <- rownames(data.1) <- rownames(x$eff.calendar)[which(!is.na(x$est.eff.calendar[, 1]))]
      }
      if (dim(d2)[2] == 1) {
        d2 <- data.2 <- t(d2)
        rownames(d2) <- rownames(data.2) <- rownames(x$eff.calendar.fit)[which(!is.na(x$est.eff.calendar.fit[, 1]))]
      }
    } else {
      if (is.null(x$est.eff.calendar)) {
        stop("Uncertainty estimates not available.\n")
      }
      d1 <- data.1 <- as.matrix(x$est.eff.calendar[which(!is.na(x$est.eff.calendar[, 1])), ])
      d2 <- data.2 <- as.matrix(x$est.eff.calendar.fit[which(!is.na(x$est.eff.calendar.fit[, 1])), ])
      if (dim(d1)[2] == 1) {
        d1 <- data.1 <- t(d1)
        rownames(d1) <- rownames(data.1) <- rownames(x$est.eff.calendar)[which(!is.na(x$est.eff.calendar[, 1]))]
      }
      if (dim(d2)[2] == 1) {
        d2 <- data.2 <- t(d2)
        rownames(d2) <- rownames(data.2) <- rownames(x$est.eff.calendar.fit)[which(!is.na(x$est.eff.calendar.fit[, 1]))]
      }

      if (length(ylim) != 0) {
        rect.length <- (ylim[2] - ylim[1]) / 5
        rect.min <- ylim[1]
      } else {
        rect.length <- (max(c(data.1[, 4], data.2[, 4]), na.rm = TRUE) - min(c(data.1[, 3], data.2[, 3]), na.rm = TRUE)) / 2
        rect.min <- min(c(data.1[, 3], data.2[, 3]), na.rm = TRUE) - rect.length
      }
    }
    p <- ggplot()
    ## xlab and ylab
    p <- p + xlab(xlab) + ylab(ylab)

    ## theme
    if (theme.bw == TRUE) {
      p <- p + theme_bw()
    }

    ## grid
    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    # horizontal 0 line
    p <- p + geom_hline(yintercept = 0, colour = lcolor[1], size = lwidth[1], linetype = ltype[1])

    TTT <- as.numeric(rownames(data.1))
    TTT.2 <- as.numeric(rownames(data.2))

    if (CI == FALSE) {
      p <- p + geom_hline(yintercept = x$att.avg, color = calendar.lcolor, size = 0.8, linetype = "dashed")
      p <- p + geom_line(aes(x = TTT.2, y = d2[, 1]), color = calendar.color, size = 1.1)
      p <- p + geom_point(aes(x = TTT, y = d1[, 1]), color = "gray50", fill = "gray50", alpha = 1, size = 1.2)
    } else {
      p <- p + geom_ribbon(aes(x = TTT.2, ymin = d2[, 3], ymax = d2[, 4]), color = calendar.cicolor, fill = calendar.cicolor, alpha = 0.5, size = 0)
      p <- p + geom_hline(yintercept = x$att.avg, color = calendar.lcolor, size = 0.8, linetype = "dashed")
      p <- p + geom_line(aes(x = TTT.2, y = d2[, 1]), color = calendar.color, size = 1.1)
      p <- p + geom_pointrange(aes(x = TTT, y = d1[, 1], ymin = d1[, 3], ymax = d1[, 4]), color = "gray50", fill = "gray50", alpha = 1, size = 0.6)
    }

    if (show.count == TRUE & !(type == "gap" | type == "equiv")) {
      T.start <- c()
      T.end <- c()
      ymin <- c()
      ymax <- c()
      T.gap <- (max(TTT) - min(TTT)) / length(TTT)
      for (i in c(1:dim(d1)[1])) {
        T.start <- c(T.start, TTT[i] - 0.25 * T.gap)
        T.end <- c(T.end, TTT[i] + 0.25 * T.gap)
        ymin <- c(ymin, rect.min)
        ymax <- c(ymax, rect.min + rect.length * d1[i, "count"] / max(d1[, "count"]))
      }
      data.toplot <- cbind.data.frame(
        xmin = T.start,
        xmax = T.end,
        ymin = ymin,
        ymax = ymax
      )
      max.count.pos <- mean(TTT[which.max(d1[, "count"])])
      p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = data.toplot, fill = count.color, size = 0.3, alpha = count.alpha, color = count.outline.color, linewidth = 0.2)
      p <- p + annotate("text",
        x = max.count.pos - 0.02 * T.gap,
        y = max(data.toplot$ymax) + 0.2 * rect.length,
        label = max(x$N.calendar), size = cex.text * 0.8, hjust = 0.5
      )
    }

    ## title
    if (is.null(main) == TRUE) {
      p <- p + ggtitle(maintext)
    } else if (main != "") {
      p <- p + ggtitle(main)
    }

    ## ylim
    if (is.null(ylim) == FALSE) {
      p <- p + coord_cartesian(ylim = ylim)
    }

    if (length(TTT) <= 10) {
      p <- p + scale_x_continuous(breaks = TTT)
    } else {
      p <- p + scale_x_continuous(labels = scaleFUN)
    }

    # ## xlim
    # if(is.null(xlim)){
    #     if(is.na(d1[1,1])){
    #         ## drop all periods before first non-missing
    #         for(j in c(2:dim(d1)[1])){
    #             if(!is.na(d1[j,1])){
    #                 xlim <- c(TTT[j],max(TTT))
    #                 break
    #             }
    #         }
    #     }
    # }
    ## xlim
    if (is.null(xlim) == FALSE) {
      p <- p + coord_cartesian(xlim = xlim)
    }


    p <- p + theme(
      legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
      legend.position = legend.pos,
      legend.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_text(size = cex.lab),
      axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text = element_text(color = "black", size = cex.axis),
      axis.text.x = element_text(size = cex.axis, angle = angle, hjust = x.h, vjust = x.v),
      axis.text.y = element_text(size = cex.axis),
      plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
    )
  }

  if (type == "heterogeneous") {
    Xs <- x[names(x) == "X"][[2]]

    if ((is.null(covariate) == TRUE) && ((dim(Xs) > 2) && (dim(Xs)[3] > 1))) {
      stop("Please provide a covariate to plot heterogeneous effects.\n")
    }

    if ((is.null(covariate) == FALSE) && (!covariate %in% x$X)) {
      stop("Please provide a valid covariate to plot heterogeneous effects.\n")
    }

    if (is.null(xlab) == TRUE) {
      xlab <- covariate
    } else if (xlab == "") {
      xlab <- NULL
    }

    if (is.null(ylab) == TRUE) {
      ylab <- ytitle
    } else if (ylab == "") {
      ylab <- NULL
    }

    D.missing <- x$D.dat
    D.missing[which(D == 0)] <- NA
    D.missing.vec <- as.vector(D.missing)

    eff.vec <- as.vector(x$eff)
    X.vec <- as.vector(x[names(x) == "X"][2]$X[, , which(x$X == covariate)])

    eff.vec <- eff.vec[which(!is.na(D.missing.vec))]
    X.vec <- X.vec[which(!is.na(D.missing.vec))]

    j <- order(X.vec)
    eff.vec <- eff.vec[j]
    X.vec <- X.vec[j]
    
    if (length(unique(X.vec)) <= 4) {
      p <- ggplot()
      ## xlab and ylab
      p <- p + xlab(xlab) + ylab(ylab)

      ## theme
      if (theme.bw == TRUE) {
        p <- p + theme_bw()
      }

      ## grid
      if (gridOff == TRUE) {
        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }

      ## title
      if (is.null(main) == TRUE) {
        p <- p + ggtitle(maintext)
      } else if (main != "") {
        p <- p + ggtitle(main)
      }

      ## Build discrete stats at evenly spaced factor positions
      x_levels_sorted <- sort(unique(X.vec))
      x_factor <- factor(X.vec, levels = x_levels_sorted, ordered = TRUE)
      # stats per level
      level_means <- tapply(eff.vec, x_factor, function(v) mean(v, na.rm = TRUE))
      level_ns    <- tapply(eff.vec, x_factor, function(v) sum(!is.na(v)))
      level_sds   <- tapply(eff.vec, x_factor, function(v) sd(v, na.rm = TRUE))
      level_se    <- level_sds / sqrt(pmax(level_ns, 1))
      level_qt    <- stats::qt(0.975, pmax(level_ns - 1, 1))
      level_ci    <- level_qt * level_se
      level_lower <- as.numeric(level_means) - level_ci
      level_upper <- as.numeric(level_means) + level_ci
      data_disc <- cbind.data.frame(
        x = factor(names(level_means), levels = names(level_means), ordered = TRUE),
        mean = as.numeric(level_means),
        lower = level_lower,
        upper = level_upper,
        count = as.numeric(level_ns)
      )

      ## core geoms (even spacing because x is factor)
      p <- p + geom_hline(yintercept = 0, colour = lcolor[1], size = lwidth[1], linetype = ltype[1])
      p <- p + geom_hline(yintercept = x$att.avg, color = heterogeneous.lcolor, size = 0.8, linetype = "dashed")
      # nicer CI + mean: thick error bars + solid point
      p <- p + geom_linerange(
        aes(x = x, ymin = .data$lower, ymax = .data$upper),
        data = data_disc, color = heterogeneous.color, linewidth = 0.9
      )
      p <- p + geom_point(
        aes(x = x, y = mean), data = data_disc,
        shape = 21, fill = "white", color = heterogeneous.color, stroke = 1, size = 2.6
      )

      ## bottom rectangles for counts under each discrete level (stick to very bottom)
      if (show.count == TRUE && any(!is.na(data_disc$count))) {
        # obtain current plot y-range from built plot if ylim not set
        current_plot_yrange <- NULL
        if (!is.null(ylim)) {
          current_plot_yrange <- ylim
        } else {
          gb <- ggplot_build(p)
          current_plot_yrange <- gb$layout$panel_scales_y[[1]]$range$range
        }
        count_bar_space_prop <- 0.20
        count_bar_space_height <- (current_plot_yrange[2] - current_plot_yrange[1]) * count_bar_space_prop
        actual_rect_length <- count_bar_space_height * 0.8
        rect_min_val <- if (!is.null(ylim)) ylim[1] else current_plot_yrange[1] - count_bar_space_height

        # Even spacing: rectangles centered on factor positions with fixed half-width
        x_idx <- as.numeric(data_disc$x)
        bar_half <- 0.1
        counts <- data_disc$count
        ymax_scaled <- rect_min_val + actual_rect_length * counts / max(counts, na.rm = TRUE)
        data_counts <- cbind.data.frame(
          xmin = x_idx - bar_half,
          xmax = x_idx + bar_half,
          ymin = rep(rect_min_val, length(counts)),
          ymax = ymax_scaled,
          xcenter = x_idx,
          count = counts
        )
        p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          data = data_counts, fill = count.color, color = count.outline.color, alpha = count.alpha, linewidth = 0.2
        )
        p <- p + geom_text(aes(x = .data$xcenter, y = .data$ymax + 0.12 * count_bar_space_height, label = .data$count),
                           data = data_counts, size = cex.text * 0.85, hjust = 0.5, vjust = 0.5, color = "#444444")
        if (!is.null(covariate.labels)) {
          p <- p + scale_x_discrete(labels = covariate.labels)
        } else {
          p <- p + scale_x_discrete(labels = as.character(x_levels_sorted))
        }
      }

      if (is.null(ylim) == FALSE) {
        p <- p + coord_cartesian(ylim = ylim)
      }

      if (is.null(xlim) == FALSE) {
        p <- p + coord_cartesian(xlim = xlim)
      }

    } else {
      plx <- predict(loess(eff.vec ~ X.vec), se = T)
      se <- stats::qt(0.975, plx$df) * plx$se
      y_hat <- plx$fit
      y_hat_lower <- y_hat - se
      y_hat_upper <- y_hat + se

      p <- ggplot()
      ## xlab and ylab
      p <- p + xlab(xlab) + ylab(ylab)

      ## theme
      if (theme.bw == TRUE) {
        p <- p + theme_bw()
      }

      ## grid
      if (gridOff == TRUE) {
        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }

      p <- p + geom_hline(yintercept = 0, colour = lcolor[1], size = lwidth[1], linetype = ltype[1])
      p <- p + geom_ribbon(aes(x = X.vec, ymin = y_hat_lower, ymax = y_hat_upper), color = heterogeneous.cicolor, fill = heterogeneous.cicolor, alpha = 0.5, size = 0)
      p <- p + geom_hline(yintercept = x$att.avg, color = heterogeneous.lcolor, size = 0.8, linetype = "dashed")
      p <- p + geom_line(aes(x = X.vec, y = y_hat), color = heterogeneous.color, size = 1.1)

      ## title
      if (is.null(main) == TRUE) {
        p <- p + ggtitle(maintext)
      } else if (main != "") {
        p <- p + ggtitle(main)
      }

      ## ylim
      if (is.null(ylim) == FALSE) {
        p <- p + coord_cartesian(ylim = ylim)
      }

      if (is.null(xlim) == FALSE) {
        p <- p + coord_cartesian(xlim = xlim)
      }

      if (show.count == TRUE) {
        current_plot_yrange <- NULL
        if (!is.null(ylim)) {
          current_plot_yrange <- ylim
        } else {
          gb <- ggplot_build(p)
          current_plot_yrange <- gb$layout$panel_scales_y[[1]]$range$range
        }
        count_bar_space_prop <- 0.20
        count_bar_space_height <- (current_plot_yrange[2] - current_plot_yrange[1]) * count_bar_space_prop
        actual_rect_length <- count_bar_space_height * 0.8
        rect_min_val <- if (!is.null(ylim)) ylim[1] else current_plot_yrange[1] - count_bar_space_height

        X.vec.for.hist <- X.vec
        if (!is.null(xlim)) {
          X.vec.for.hist <- X.vec.for.hist[which(X.vec.for.hist >= min(xlim) & X.vec.for.hist <= max(xlim))]
        }
        if (length(na.omit(X.vec.for.hist)) > 0) {
          breaks <- pretty(range(X.vec.for.hist, na.rm = TRUE), n = 50)
          h <- hist(X.vec.for.hist, breaks = breaks, plot = FALSE)
          bin_xmin <- h$breaks[-length(h$breaks)]
          bin_xmax <- h$breaks[-1]
          counts <- h$counts
          ymax_scaled <- rect_min_val + actual_rect_length * counts / max(counts, na.rm = TRUE)
          data.toplot <- cbind.data.frame(
            xmin = bin_xmin,
            xmax = bin_xmax,
            ymin = rep(rect_min_val, length(counts)),
            ymax = ymax_scaled
          )
          max_idx <- which.max(counts)
          max_count_pos <- (bin_xmin[max_idx] + bin_xmax[max_idx]) / 2
          p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            data = data.toplot, fill = count.color, size = 0.3, alpha = count.alpha, color = count.outline.color, linewidth = 0.2
          )
          p <- p + annotate("text",
            x = max_count_pos,
            y = max(data.toplot$ymax) + 0.12 * count_bar_space_height,
            label = max(counts, na.rm = TRUE), size = cex.text * 0.8, hjust = 0.5
          )
          if (is.null(ylim)) {
            final_yrange_min <- min(current_plot_yrange[1], rect_min_val)
            final_yrange_max <- current_plot_yrange[2]
            p <- p + coord_cartesian(ylim = c(final_yrange_min, final_yrange_max))
          }
        }
      }
    }
 
    p <- p + theme(
      legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
      legend.position = legend.pos,
      legend.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_text(size = cex.lab),
      axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text = element_text(color = "black", size = cex.axis),
      axis.text.x = element_text(size = cex.axis, angle = angle, hjust = x.h, vjust = x.v),
      axis.text.y = element_text(size = cex.axis),
      plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
    )
  }

  if (type == "box") {
    if (is.null(xlab) == TRUE) {
      xlab <- index[2]
    } else if (xlab == "") {
      xlab <- NULL
    }
    if (is.null(ylab) == TRUE) {
      ylab <- "Estimated Treatment Effects"
    } else if (ylab == "") {
      ylab <- NULL
    }
    if (is.null(main) == TRUE) {
      main <- maintext
    } else if (main == "") {
      main <- NULL
    }

    p <- ggplot()
    ## xlab and ylab
    p <- p + xlab(xlab) + ylab(ylab)

    ## theme
    if (theme.bw == TRUE) {
      p <- p + theme_bw()
    }

    ## title
    if (!is.null(main)) {
      p <- p + ggtitle(main)
    }

    ## grid
    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    # horizontal 0 line
    p <- p + geom_hline(yintercept = 0, colour = lcolor[1], size = lwidth[1], linetype = ltype[1])

    complete.index.eff <- which(!is.na(x$eff))
    complete.index.time <- which(!is.na(x$T.on))
    complete.index <- intersect(complete.index.eff, complete.index.time)
    eff.use <- x$eff[complete.index]
    time.use <- x$T.on[complete.index]
    id.mat <- rep(colnames(x$eff), each = dim(x$eff)[1])
    id.use <- id.mat[complete.index]
    data.toplot <- cbind.data.frame(time = time.use, id = id.use, eff = eff.use)
    data.count <- cbind.data.frame(time = x$time, count = x$count)

    if (start0 == TRUE) {
      data.toplot$time <- data.toplot$time - 1
      data.count$time <- data.count$time - 1
    }

    if (!is.null(xlim)) {
      data.count <- data.count[which(data.count[, "time"] >= min(xlim) & data.count[, "time"] <= max(xlim)), ]
      data.toplot <- data.toplot[which(data.toplot[, "time"] >= min(xlim) & data.toplot[, "time"] <= max(xlim)), ]
    }


    data.use <- merge(data.toplot, data.count, by = "time")
    # print(data.use)

    if (length(ylim) != 0) {
      rect.length <- (ylim[2] - ylim[1]) / 5
      rect.min <- ylim[1]
    } else {
      rect.length <- (max(data.use[, "eff"], na.rm = TRUE) - min(data.use[, "eff"], na.rm = TRUE)) / 3
      rect.min <- min(data.use[, "eff"], na.rm = TRUE) - rect.length
    }

    if (start0 == FALSE) {
      data.pre.1 <- data.use[which(data.use$time <= 0 & data.use$count >= 10), ]
      data.pre.2 <- data.use[which(data.use$time <= 0 & data.use$count < 10), ]
      data.post.1 <- data.use[which(data.use$time > 0 & data.use$count >= 10), ]
      data.post.2 <- data.use[which(data.use$time > 0 & data.use$count < 10), ]
    } else {
      data.pre.1 <- data.use[which(data.use$time < 0 & data.use$count >= 10), ]
      data.pre.2 <- data.use[which(data.use$time < 0 & data.use$count < 10), ]
      data.post.1 <- data.use[which(data.use$time >= 0 & data.use$count >= 10), ]
      data.post.2 <- data.use[which(data.use$time >= 0 & data.use$count < 10), ]
    }


    levels <- as.factor(as.character(data.count[, 1]))
    data.pre.1$time <- as.character(data.pre.1$time)
    data.pre.2$time <- as.character(data.pre.2$time)
    data.post.1$time <- as.character(data.post.1$time)
    data.post.2$time <- as.character(data.post.2$time)
    data.pre.1$time <- factor(data.pre.1$time, levels = levels)
    data.pre.2$time <- factor(data.pre.2$time, levels = levels)
    data.post.1$time <- factor(data.post.1$time, levels = levels)
    data.post.2$time <- factor(data.post.2$time, levels = levels)

    p <- p + geom_boxplot(aes(x = .data$time, y = .data$eff),
      position = "dodge", alpha = 0.5,
      data = data.pre.1, fill = box.control,
      outlier.fill = box.control, outlier.size = 1.25,
      outlier.color = box.control,
    )
    p <- p + geom_boxplot(aes(x = .data$time, y = .data$eff),
      position = "dodge", alpha = 0.5,
      data = data.post.1, fill = box.treat, outlier.fill = box.treat,
      outlier.size = 1.25, outlier.color = box.treat,
    )

    p <- p + geom_point(aes(x = .data$time, y = .data$eff),
      data = data.post.2,
      color = box.treat, size = 1.25, alpha = 0.8
    )
    p <- p + geom_point(aes(x = .data$time, y = .data$eff),
      data = data.pre.2,
      color = box.control, size = 1.25, alpha = 0.8
    )

    if (show.count == TRUE & !(type == "gap" | type == "equiv")) {
      T.start <- c()
      T.end <- c()
      ymin <- c()
      ymax <- c()
      T.gap <- 1
      for (i in c(1:dim(data.count)[1])) {
        T.start <- c(T.start, data.count[i, 1] - 0.25 * T.gap - min(data.count[, 1]) + 1)
        T.end <- c(T.end, data.count[i, 1] + 0.25 * T.gap - min(data.count[, 1]) + 1)
        ymin <- c(ymin, rect.min)
        ymax <- c(ymax, rect.min + rect.length * data.count[i, 2] / max(data.count[, 2]))
      }
      data.toplot <- cbind.data.frame(
        xmin = T.start,
        xmax = T.end,
        ymin = ymin,
        ymax = ymax
      )
      max.count.pos <- data.count[which.max(data.count[, 2]), 1][1] - min(data.count[, 1]) + 1
      p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = data.toplot, fill = count.color, size = 0.3, alpha = count.alpha, color = count.outline.color, linewidth = 0.2)
      p <- p + annotate("text",
        x = max.count.pos - 0.02 * T.gap,
        y = max(data.toplot$ymax) + 0.1 * rect.length,
        label = max(data.count[, 2]), size = cex.text * 0.7, hjust = 0.5
      )
    }

    ## ylim
    if (is.null(ylim) == FALSE) {
      p <- p + coord_cartesian(ylim = ylim)
    }

    p <- p + theme(
      legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
      legend.position = legend.pos,
      legend.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_text(size = cex.lab),
      axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text = element_text(color = "black", size = cex.axis),
      axis.text.x = element_text(size = cex.axis, angle = angle, hjust = x.h, vjust = x.v),
      axis.text.y = element_text(size = cex.axis),
      plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
    )

    if (is.null(xticklabels) == FALSE) {
      xticklabels.all <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
      labels <- c()
      for (xticklabel.all in xticklabels.all) {
        if (xticklabel.all %in% xticklabels) {
          labels <- c(labels, xticklabel.all)
        } else {
          labels <- c(labels, "")
        }
      }
      p <- p + scale_x_discrete(limits = levels, breaks = xticklabels.all, labels = labels)
    } else {
      p <- p + scale_x_discrete(limits = levels)
    }
  }

  if (type == "status") {
    if (is.null(xlab) == TRUE) {
      xlab <- index[2]
    } else if (xlab == "") {
      xlab <- NULL
    }
    if (is.null(ylab) == TRUE) {
      ylab <- index[1]
    } else if (ylab == "") {
      ylab <- NULL
    }
    if (is.null(main) == TRUE) {
      main <- "Treatment Status"
    } else if (main == "") {
      main <- NULL
    }

    m <- obs.missing
    m <- as.matrix(m[show, which(colnames(m) %in% id)])

    all <- unique(c(m))
    col <- col2 <- breaks <- label <- NULL

    if (1 %in% all) {
      col <- c(col, status.treat.color)
      col2 <- c(col2, "1" = NA)
      breaks <- c(breaks, 1)
      label <- c(label, "Under Treatment")
    }
    if (2 %in% all) {
      col <- c(col, status.control.color)
      col2 <- c(col2, "2" = NA)
      breaks <- c(breaks, 2)
      label <- c(label, "Under Control")
    }
    if (3 %in% all) {
      col <- c(col, status.missing.color)
      col2 <- c(col2, "3" = NA)
      breaks <- c(breaks, 3)
      label <- c(label, "Missing")
    }
    if (4 %in% all) {
      col <- c(col, status.removed.color)
      col2 <- c(col2, "4" = NA)
      breaks <- c(breaks, 4)
      label <- c(label, "Removed")
    }
    if (5 %in% all) {
      col2 <- c(col2, "5" = NA)
      breaks <- c(breaks, 5)
      if (placeboTest == TRUE) {
        col <- c(col, status.placebo.color)
        label <- c(label, "Placebo Tests")
      } else if (carryoverTest == TRUE) {
        col <- c(col, status.carryover.color)
        label <- c(label, "Carryover Tests")
      }
    }
    if (6 %in% all) {
      col <- c(col, status.carryover.rm.color)
      col2 <- c(col2, "6" = NA)
      breaks <- c(breaks, 6)
      label <- c(label, "Carryover Removed")
    }
    if (7 %in% all) {
      col <- c(col, status.balanced.post.color)
      col2 <- c(col2, "7" = NA)
      breaks <- c(breaks, 7)
      label <- c(label, "Balanced Sample: Post")
    }
    if (8 %in% all) {
      col <- c(col, status.balanced.pre.color)
      col2 <- c(col2, "8" = NA)
      breaks <- c(breaks, 8)
      label <- c(label, "Balanced Sample: Pre")
    }

    TT <- dim(m)[1]
    N <- dim(m)[2]
    units <- rep(rev(1:N), each = TT)
    period <- rep(1:TT, N)
    res <- c(m)
    data <- cbind.data.frame(units = units, period = period, res = res)
    data[, "res"] <- as.factor(data[, "res"])

    ## check if N >= 200
    if (dim(m)[2] >= 200) {
      if (axis.lab == "both") {
        axis.lab <- "time"
      } else if (axis.lab == "unit") {
        axis.lab <- "off"
      }
    }

    ## labels
    N.b <- 1:N
    if (axis.lab == "both") {
      if (length(axis.lab.gap) == 2) {
        x.gap <- axis.lab.gap[1]
        y.gap <- axis.lab.gap[2]
      } else {
        x.gap <- y.gap <- axis.lab.gap[1]
      }
    } else {
      x.gap <- y.gap <- axis.lab.gap[1]
    }

    if (x.gap != 0) {
      T.b <- seq(from = 1, to = length(show), by = (x.gap + 1))
    }
    if (y.gap != 0) {
      N.b <- seq(from = N, to = 1, by = -(y.gap + 1))
    }
    id <- rev(id)

    p <- ggplot(data, aes(
      x = period, y = units,
      fill = res
    ), position = "identity")

    if (gridOff == FALSE) {
      p <- p + geom_tile(colour = status.background.color, linewidth = 0.05, stat = "identity")
    } else {
      p <- p + geom_tile(stat = "identity")
    }

    p <- p +
      labs(
        x = xlab, y = ylab,
        title = main
      ) +
      theme_bw() +
      scale_fill_manual(NA, breaks = breaks, values = col, labels = label)

    # if(4%in%all) {
    #    p <- p + geom_point(aes(colour=res),size=0.5)
    #    p <- p + scale_color_manual(NA, breaks=breaks,
    #                                values=col2, labels=label)
    # }

    p <- p +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = status.background.color, linewidth = 0.5, linetype = "solid"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = cex.lab),
        axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.text = element_text(color = "black", size = cex.axis),
        axis.text.x = element_text(size = cex.axis, angle = angle, hjust = x.h, vjust = x.v),
        axis.text.y = element_text(size = cex.axis),
        plot.background = element_rect(fill = status.background.color),
        legend.background = element_rect(fill = status.background.color),
        legend.position = legend.pos,
        legend.margin = margin(c(0, 5, 5, 0)),
        legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
        legend.title = element_blank(),
        plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(8, 0, 8, 0))
      )

    if (is.null(xticklabels) == FALSE) {
      xticks <- c()
      for (xticklabel in xticklabels) {
        xticks <- c(xticks, which(time.label == xticklabel))
      }
    } else {
      xticks <- T.b
    }

    if (is.null(yticklabels) == FALSE) {
      yticks <- c()
      for (yticklabel in yticklabels) {
        yticks <- c(yticks, which(id == yticklabel))
      }
    } else {
      yticks <- N.b
    }

    if (axis.lab == "both") {
      p <- p + scale_x_continuous(expand = c(0, 0), breaks = xticks, labels = time.label[xticks]) +
        scale_y_continuous(expand = c(0, 0), breaks = yticks, labels = id[yticks])
    } else if (axis.lab == "unit") {
      p <- p + scale_x_continuous(expand = c(0, 0), breaks = xticks, labels = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = yticks, labels = id[yticks])
    } else if (axis.lab == "time") {
      p <- p + scale_x_continuous(expand = c(0, 0), breaks = xticks, labels = time.label[xticks]) +
        scale_y_continuous(expand = c(0, 0), breaks = yticks, labels = NULL)
    } else if (axis.lab == "off") {
      p <- p + scale_x_continuous(expand = c(0, 0), breaks = 1:length(show), labels = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = 1:N, labels = NULL)
    }

    if (length(all) >= 3) {
      p <- p + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
    }
  } else if (type == "sens") {
    # Determine plot-specific variables based on the 'restrict' option
    sens_data_obj <- NULL
    x_var_name <- NULL
    default_xlab_val <- NULL

    if (restrict == "rm") {
      if (is.null(x$sensitivity.rm) || is.null(x$sensitivity.rm$results) || is.null(x$sensitivity.rm$original)) {
        stop("No sensitivity results found in x$sensitivity.rm. Required for restrict = 'rm'.")
      }
      sens_data_obj <- x$sensitivity.rm
      x_var_name <- "Mbar"
      default_xlab_val <- expression(bar(M))
    } else if (restrict == "sm") {
      if (is.null(x$sensitivity.smooth) || is.null(x$sensitivity.smooth$results) || is.null(x$sensitivity.smooth$original)) {
        stop("No sensitivity results found in x$sensitivity.smooth. Required for restrict = 'sm'.")
      }
      sens_data_obj <- x$sensitivity.smooth
      x_var_name <- "M"
      default_xlab_val <- "M"
    }

    # Prepare data for plotting
    data_original <- sens_data_obj$original
    data_results <- sens_data_obj$results

    # Assign a default x-axis position for the original estimate to plot it left of zero
    if (!(x_var_name %in% colnames(data_original))) {
      data_original[[x_var_name]] <- -0.05
    }

    # Validate required columns
    if (!all(c("lb", "ub") %in% colnames(data_original))) {
      stop("Original sensitivity data requires 'lb' and 'ub' columns.")
    }
    if (!all(c(x_var_name, "lb", "ub") %in% colnames(data_results))) {
      stop(paste("Sensitivity results data requires '", x_var_name, "', 'lb', and 'ub' columns."))
    }

    # Assign color groups and combine into a single data frame
    data_original$color_group <- "Original"
    data_results$color_group <- "Robust Confidence Set"
    data_combined <- rbind(data_original, data_results)

    # --- Handle Plot Labels and Title ---
    # Use user-provided 'main' or set a default title
    maintext <- if (is.null(main)) "Smoothness Restriction Sensitivity Analysis" else main

    # Use user-provided 'xlab' or set a default based on 'restrict'
    xlab_final <- if (is.null(xlab)) default_xlab_val else xlab

    # Use user-provided 'ylab' or set a default
    ylab_final <- if (is.null(ylab)) "Treatment Effect" else ylab

    # --- Create the Plot ---
    p <- ggplot(data_combined, aes_string(x = x_var_name))

    if (theme.bw == TRUE) {
      p <- p + theme_bw()
    }
    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    p <- p + geom_hline(yintercept = 0, color = lcolor[1], size = lwidth[1], linetype = ltype[1]) +
      geom_errorbar(
        aes(ymin = .data$lb, ymax = .data$ub, color = .data$color_group),
        width = 0.02, # Width of error bar caps
        linewidth = 1
      ) +
      scale_color_manual(values = c("Original" = sens.original.color, "Robust Confidence Set" = sens.colors[1])) +
      labs(
        x = xlab_final,
        y = ylab_final
      )

    # Add title only if it's not an empty string
    if (!is.null(maintext) && maintext != "") {
      p <- p + ggtitle(maintext)
    }

    # --- Apply Consistent Theming ---
    p <- p + theme(
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(size = cex.legend),
      axis.title = element_text(size = cex.lab),
      axis.text = element_text(size = cex.axis),
      plot.title = element_text(
        size = cex.main,
        hjust = 0.5,
        face = "bold",
      )
    )
  } else if (type == "sens_es") {
    if (restrict == "rm") {
      if (is.null(x$sensitivity.rm) || is.null(x$sensitivity.rm$periods)) {
        stop("No period-by-period Smoothness data found in x$sensitivity.rm$periods.")
      }

      dte_output <- x$sensitivity.rm$periods
      required_cols <- c("lb", "ub", "postPeriod", "Mbar")
      if (!all(required_cols %in% names(dte_output))) {
        stop("sensitivity.rm$periods is missing one of: 'lb','ub','postPeriod','Mbar'.")
      }

      # Reorder the data so that higher values of Mbar are plotted first.
      dte_output <- dte_output[order(dte_output$Mbar, decreasing = TRUE), ]

      # Build base ES plot
      fect.output.p <- as.data.frame(x$est.att)
      fect.output.p$Time <- as.numeric(rownames(fect.output.p))
      if (is.null(ylim)) {
        ylim <- c(
          min(c(fect.output.p$CI.lower, dte_output$lb)) * 1.3,
          max(c(fect.output.p$CI.upper, dte_output$ub)) * 1.3
        )
      }
      p <- esplot(
        fect.output.p,
        Period = "Time",
        Estimate = "ATT",
        SE = "S.E.",
        CI.lower = "CI.lower",
        CI.upper = "CI.upper",
        Count = "count",
        show.count = show.count,
        proportion = proportion,
        show.points = show.points,
        ci.outline = ci.outline,
        connected = connected,
        color = color,
        count.color = count.color,
        count.alpha = count.alpha,
        count.outline.color = count.outline.color,
        highlight.periods = x$placebo.period[1]:x$placebo.period[2],
        highlight.colors = rep(placebo.color, x$placebo.period[2] - x$placebo.period[1] + 1),
        xlab = xlab,
        ylab = ylab,
        main = main,
        fill.gap = FALSE,
        lcolor = lcolor,
        lwidth = lwidth,
        ltype = ltype,
        cex.main = cex.main / 16,
        cex.axis = cex.axis / 15,
        cex.lab = cex.lab / 15,
        cex.text = cex.text / 5,
        axis.adjust = axis.adjust,
        xlim = xlim,
        ylim = ylim,
        gridOff = gridOff,
        start0 = start0,
        est.lwidth = est.lwidth,
        theme.bw = theme.bw,
        est.pointsize = est.pointsize
      )

      mbar_levels <- sort(unique(dte_output$Mbar))
      n_colors <- length(mbar_levels)
      final_palette <- sens.colors[1:min(n_colors, length(sens.colors))]
      p <- p +
        geom_linerange(
          aes(x = .data$postPeriod + 0.2, ymin = .data$lb, ymax = .data$ub, color = factor(.data$Mbar, levels = mbar_levels)),
          data = dte_output,
          linewidth = 1, inherit.aes = FALSE
        ) +
        scale_color_manual(name = "Mbar", values = setNames(final_palette, mbar_levels)) +
        guides(color = guide_legend(title = expression(bar(M)))) +
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c("left", "top")
        )
    } else if (restrict == "sm") {
      if (is.null(x$sensitivity.smooth) || is.null(x$sensitivity.smooth$periods)) {
        stop("No period-by-period Smoothness data found in x$sensitivity.smooth$periods.")
      }

      dte_output <- x$sensitivity.smooth$periods
      required_cols <- c("lb", "ub", "postPeriod", "M")
      if (!all(required_cols %in% names(dte_output))) {
        stop("sensitivity.smooth$periods is missing one of: 'lb','ub','postPeriod','M'.")
      }

      # Reorder the data so that higher values of M are plotted first.
      dte_output <- dte_output[order(dte_output$M, decreasing = TRUE), ]

      # Build base ES plot
      fect.output.p <- as.data.frame(x$est.att)
      fect.output.p$Time <- as.numeric(rownames(fect.output.p))
      if (is.null(ylim)) {
        ylim <- c(
          min(c(fect.output.p$CI.lower, dte_output$lb)) * 1.3,
          max(c(fect.output.p$CI.upper, dte_output$ub)) * 1.3
        )
      }
      p <- esplot(
        fect.output.p,
        Period = "Time",
        Estimate = "ATT",
        SE = "S.E.",
        CI.lower = "CI.lower",
        CI.upper = "CI.upper",
        Count = "count",
        show.count = show.count,
        proportion = proportion,
        show.points = show.points,
        ci.outline = ci.outline,
        connected = connected,
        color = color,
        count.color = count.color,
        count.alpha = count.alpha,
        count.outline.color = count.outline.color,
        highlight.periods = x$placebo.period[1]:x$placebo.period[2],
        highlight.colors = rep(placebo.color, x$placebo.period[2] - x$placebo.period[1] + 1),
        xlab = xlab,
        ylab = ylab,
        main = main,
        xlim = xlim,
        ylim = ylim,
        gridOff = gridOff,
        fill.gap = FALSE,
        lcolor = lcolor,
        lwidth = lwidth,
        ltype = ltype,
        cex.main = cex.main / 16,
        cex.axis = cex.axis / 15,
        cex.lab = cex.lab / 15,
        cex.text = cex.text / 5,
        axis.adjust = axis.adjust,
        start0 = start0,
        theme.bw = theme.bw,
        est.lwidth = est.lwidth,
        est.pointsize = est.pointsize
      )

      m_levels <- sort(unique(dte_output$M))
      n_colors <- length(m_levels)
      final_palette <- sens.colors[1:min(n_colors, length(sens.colors))]
      p <- p +
        geom_linerange(
          aes(x = .data$postPeriod + 0.2, ymin = .data$lb, ymax = .data$ub, color = factor(.data$M, levels = m_levels)),
          data = dte_output,
          linewidth = 1, inherit.aes = FALSE
        ) +
        scale_color_manual(name = "M", values = setNames(final_palette, m_levels)) +
        guides(color = guide_legend(title = "M")) +
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c("left", "top")
        )
    }
  } else if (type == "cumul") {
    if (is.null(x$effect.est.att)) {
      stop("No period-by-period cumulative ATT data found in x$est.eff.")
    }
    if (is.null(main)) {
      main <- "Estimated Cumulative Treatment Effects"
    }
    p <- esplot(
      x$effect.est.att,
      Estimate = "ATT",
      SE = "S.E.",
      CI.lower = "CI.lower",
      CI.upper = "CI.upper",
      main = main,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      fill.gap = FALSE,
      lcolor = lcolor,
      lwidth = lwidth,
      ltype = ltype,
      cex.main = cex.main / 16,
      cex.axis = cex.axis / 15,
      cex.lab = cex.lab / 15,
      cex.text = cex.text / 5,
      show.points = show.points,
      ci.outline = ci.outline,
      axis.adjust = axis.adjust,
      Count = "count",
      show.count = FALSE,
      proportion = proportion,
      est.pointsize = est.pointsize,
      est.lwidth = est.lwidth,
      color = color,
      count.color = count.color,
      count.alpha = count.alpha,
      count.outline.color = count.outline.color,
      theme.bw = theme.bw,
      connected = connected,
      only.post = TRUE,
    )
  }

  if (!is.null(xbreaks)) {
    p <- p + scale_x_continuous(breaks = xbreaks)
  }

  if (!is.null(ybreaks)) {
    p <- p + scale_y_continuous(breaks = ybreaks)
  }

  if (!is.null(xangle)) {
    p <- p + theme(axis.text.x = element_text(angle = xangle))
  }

  if (!is.null(yangle)) {
    p <- p + theme(axis.text.y = element_text(angle = yangle))
  }

  if (return.test == TRUE) {
    return(list(p = p, test.out = test.out))
  } else {
    return(p)
  }

  # suppressWarnings(print(p))
  if (return.test == TRUE) {
    return(list(p = p, test.out = test.out))
  } else {
    return(p)
  }
}
