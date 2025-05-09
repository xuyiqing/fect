## new plot
# x: a fect object
# type of the plot; axes limits; axes labels;
# main: whether to show the title;
# id: plot a part of units

plot.fect <- function(
    x,
    type = NULL, # gap, equiv, status, exit, factors, loadings, calendar, counterfactual
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
    gridOff = TRUE,
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
    line.color = NULL,
    line.width = c(1.2, 0.5),
    lcolor = NULL,
    lwidth = NULL,
    theme = NULL,
    connected = NULL,
    ci.outline = FALSE,
    color =  NULL,
    est.linewidth = 1.2,
    est.pointsize = 3,
    count.color =  NULL,
    placebo.color =  NULL,
    carryover.color =  NULL,
    carryover.rm.color = NULL,
    sens.original.color = NULL,
    sens.colors =  NULL,
    counterfactual.color =  NULL,
    counterfactual.raw.controls.color =  NULL,
    counterfactual.raw.treated.color =  NULL,
    counterfactual.linetype =  NULL,
    box.control =  NULL,
    box.treat =  NULL,
    calendar.color =  NULL,
    calendar.line.color =  NULL,
    equiv.color =  NULL,
    status.treat.color =  NULL,
    status.control.color = NULL,
    status.missing.color =  NULL,
    status.removed.color =  NULL,
    status.placebo.color = NULL,
    status.carryover.color =  NULL,
    status.carryover.rm.color =  NULL,
    status.balanced.post.color = NULL,
    status.balanced.pre.color = NULL,
    status.background.color = NULL,
    ...) {
  group <- ATT5 <- ATT6 <- CI.lower.90 <- CI.lower6 <- CI.upper.90 <- CI.upper6 <- L1 <- eff <- NULL

  if (!missing(vis)) {
    warning("'vis' is deprecated and will be removed in future versions.", call. = FALSE)
  }
  if (is.null(theme)) {
    if (is.null(connected))              connected              <- FALSE
    if (is.null(color))                  color                  <- "black"
    if (is.null(count.color))            count.color            <- "gray80"
    if (is.null(placebo.color))          placebo.color          <- "blue"
    if (is.null(carryover.color))        carryover.color        <- "red"
    if (is.null(carryover.rm.color))     carryover.rm.color     <- "blue"
    if (is.null(sens.original.color))    sens.original.color    <- "darkblue"
    if (is.null(sens.colors))            sens.colors            <- c("#218C23","#FF34B4","#FF521B","#2B59C3")
    if (is.null(counterfactual.color))   counterfactual.color   <- "blue"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#D5E1ED"
    if (is.null(counterfactual.raw.treated.color))  counterfactual.raw.treated.color  <- "#77777750"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "dashed"
    if (is.null(box.control))            box.control            <- "lightblue"
    if (is.null(box.treat))              box.treat              <- "lightpink"
    if (is.null(calendar.color))         calendar.color         <- "lightblue"
    if (is.null(calendar.line.color))    calendar.line.color    <- "red"
    if (is.null(equiv.color))            equiv.color            <- "red"
    if (is.null(status.treat.color))     status.treat.color     <- "#06266F"
    if (is.null(status.control.color))   status.control.color   <- "#B0C4DE"
    if (is.null(status.missing.color))   status.missing.color   <- "#eeeeee"
    if (is.null(status.removed.color))   status.removed.color   <- "#A9A9A9"
    if (is.null(status.placebo.color))   status.placebo.color   <- "#66C2A5"
    if (is.null(status.carryover.color)) status.carryover.color <- "#E78AC3"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#ffc425"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#00852B"
    if (is.null(status.balanced.pre.color))  status.balanced.pre.color  <- "#A5CA18"
    if (is.null(status.background.color)) status.background.color <- "#FFFFFF"
  }
  else if (theme == "vibrant") {
    if (is.null(connected))              connected              <- TRUE
    if (is.null(color))                  color                  <- "#054A91"
    if (is.null(count.color))            count.color            <- "#E6AF2E"
    if (is.null(placebo.color))          placebo.color          <- "#386641"
    if (is.null(carryover.color))        carryover.color        <- "#A40E4C"
    if (is.null(carryover.rm.color))     carryover.rm.color     <- "#FF5400"
    if (is.null(sens.original.color))    sens.original.color    <- "#054A91"
    if (is.null(sens.colors))            sens.colors            <- c("#A40E4C", "#FF5400", "#E6AF2E", "#386641", "#ACC3DA")
    if (is.null(counterfactual.color))       counterfactual.color       <- "#777777"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#D5E1ED"
    if (is.null(counterfactual.raw.treated.color))  counterfactual.raw.treated.color  <- "#77777750"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "dashed"
    if (is.null(box.control))            box.control            <- "#ACC3DA"
    if (is.null(box.treat))              box.treat              <- "#E1AFC3"
    if (is.null(calendar.color))         calendar.color         <- "#ACC3DA"
    if (is.null(calendar.line.color))    calendar.line.color    <- "#054A91"
    if (is.null(equiv.color))            equiv.color            <- "#A40E4C"
    if (is.null(status.treat.color))        status.treat.color        <- "#054A91"
    if (is.null(status.control.color))      status.control.color      <- "#ACC3DA"
    if (is.null(status.missing.color))      status.missing.color      <- "#dddddd"
    if (is.null(status.removed.color))      status.removed.color      <- "#D7E8E0"
    if (is.null(status.placebo.color))      status.placebo.color      <- "#386641"
    if (is.null(status.carryover.color))    status.carryover.color    <- "#A40E4C"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#FF5400"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#E6AF2E"
    if (is.null(status.balanced.pre.color))  status.balanced.pre.color  <- "#777777"
    if (is.null(status.background.color))    status.background.color    <- "#FFFFFF"

    }
  else if (theme %in% c("grayscale","greyscale")) {
    if (is.null(connected))              connected              <- FALSE
    if (is.null(color))                  color                  <- "black"
    if (is.null(count.color))            count.color            <- "gray80"
    if (is.null(placebo.color))          placebo.color          <- "gray40"
    if (is.null(carryover.color))        carryover.color        <- "gray70"
    if (is.null(carryover.rm.color))     carryover.rm.color     <- "gray40"
    if (is.null(sens.original.color))    sens.original.color    <- "gray20"
    if (is.null(sens.colors))            sens.colors            <- c("gray80", "gray50", "black")
    if (is.null(counterfactual.color))       counterfactual.color       <- "#777777"
    if (is.null(counterfactual.raw.controls.color)) counterfactual.raw.controls.color <- "#eeeeee"
    if (is.null(counterfactual.raw.treated.color))  counterfactual.raw.treated.color  <- "#666666"
    if (is.null(counterfactual.linetype)) counterfactual.linetype <- "dashed"
    if (is.null(box.control))            box.control            <- "gray90"
    if (is.null(box.treat))              box.treat              <- "gray50"
    if (is.null(calendar.color))         calendar.color         <- "gray80"
    if (is.null(calendar.line.color))    calendar.line.color    <- "black"
    if (is.null(equiv.color))            equiv.color            <- "gray80"
    if (is.null(status.treat.color))        status.treat.color        <- "black"
    if (is.null(status.control.color))      status.control.color      <- "gray90"
    if (is.null(status.missing.color))      status.missing.color      <- "white"
    if (is.null(status.removed.color))      status.removed.color      <- "white"
    if (is.null(status.placebo.color))      status.placebo.color      <- "gray50"
    if (is.null(status.carryover.color))    status.carryover.color    <- "gray50"
    if (is.null(status.carryover.rm.color)) status.carryover.rm.color <- "#ffffff"
    if (is.null(status.balanced.post.color)) status.balanced.post.color <- "#aaaaaa"
    if (is.null(status.balanced.pre.color))  status.balanced.pre.color  <- "#eeeeee"
    if (is.null(status.background.color))    status.background.color    <- "#FFFFFF"

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
    lwidth <- 1.5
    if (theme.bw == TRUE) {
      lwidth <- 1
    }
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
  if (!(restrict %in% c("rm","sm"))) {
    stop("\"restrict\" option misspecified. Must be either \"rm\" or \"sm\".")
  }

  # check the key option type
  if (!is.null(type)) {
    if (type == "ct") {
      type <- "counterfactual"
    }
    if (!type %in% c("status", "gap", "equiv", "exit", "factors", "loadings", "calendar", "box", "counterfactual", "sens","sens_es", "cumul")) {
      stop("\"type\" option misspecified. Must be one of the following:\"status\",\"gap\",\"equiv\",\"exit\",\"calendar\",\"box\",\"counterfactual\",\"equiv\",\"sens\",\"sens_es\",\"cumul\".")
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

  if (type == "es") {
    type = "gap"
  }
  type.old <- type

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

        if (x$force %in% c(1, 3) & include.FE == TRUE) {
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
          p <- ggplot(data, aes(x = group, y = L1, fill = group)) +
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
                                 mapping = aes(color = group, fill = group),
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
                                 mapping = aes(color = group),
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
                                 mapping = aes(color = group),
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

      if (is.null(line.color) == TRUE) {
        if (theme.bw == FALSE) {
          line.color <- "#AAAAAA70"
        } else {
          line.color <- "white"
        }
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
      if (length(xlim) != 0) {
        show <- which(time >= xlim[1] & time <= xlim[2])
        proportion = 0
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
        p <- p + xlab(xlab) + ylab(ylab) + ggtitle(main) +
          geom_hline(yintercept = 0, colour = line.color, size = 2) +
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
                               colour = group,
                               group = group
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
      xlim[2] <- 0
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
    maintext <- "ATT by Calendar Time"
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
    if (type == "status") {
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
  # Main plotting logic for counterfactual type
  if (type %in% c("counterfactual", "ct")) {

    # --- Basic Setup ---
    scaleFUN <- function(x) sprintf("%.f", x) # integer value at x axis

    # Validate 'raw' option
    if (!raw %in% c("none", "band", "all")) {
      cat("\"raw\" option misspecifed. Reset to \"none\".\n")
      raw <- "none"
    }
    # Validate 'id' option
    if (!is.null(id) && length(id) > 1) {
      stop("More than 1 element in \"id\".")
    }

    # Axis text angle adjustment
    if (axis.adjust == TRUE) {
      angle <- 45; x.v <- 1; x.h <- 1
    } else {
      angle <- 0; x.v <- 0; x.h <- 0.5 # Default h just for text
    }

    # --- Data Extraction and Cleaning ---
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

    I_mat <- as.matrix(I_orig); colnames(I_mat) <- rawid_orig; rownames(I_mat) <- time_orig_vec
    II_mat <- as.matrix(II_orig); colnames(II_mat) <- rawid_orig; rownames(II_mat) <- time_orig_vec
    Y_mat <- as.matrix(Y_orig); colnames(Y_mat) <- rawid_orig; rownames(Y_mat) <- time_orig_vec
    Y.ct_mat <- as.matrix(Y.ct_orig); colnames(Y.ct_mat) <- rawid_orig; rownames(Y.ct_mat) <- time_orig_vec
    D_mat <- as.matrix(D_orig); colnames(D_mat) <- rawid_orig; rownames(D_mat) <- time_orig_vec

    tr_idx_logical <- (1:ncol(Y_mat)) %in% tr

    I.tr <- I_mat[, tr_idx_logical, drop = FALSE]
    II.tr <- II_mat[, tr_idx_logical, drop = FALSE]
    D.tr <- D_mat[, tr_idx_logical, drop = FALSE] # This D.tr is for the selected treated units
    Y.tr <- Y_mat[, tr_idx_logical, drop = FALSE] # This Y.tr is for the selected treated units
    Y.ct <- Y.ct_mat[, tr_idx_logical, drop = FALSE] # This Y.ct is for the selected treated units

    co_idx_logical <- (1:ncol(Y_mat)) %in% co
    Y.co <- Y_mat[, co_idx_logical, drop = FALSE]

    # T0 calculation must use the original D.tr state.
    if (!0 %in% I.tr) {
      pre <- as.matrix(D.tr == 0 & II.tr == 1)
    } else {
      pre <- as.matrix(D.tr == 0 & I.tr == 1 & II.tr == 1)
    }
    T0 <- apply(pre, 2, sum, na.rm = TRUE) # Count of pre-treatment periods for each treated unit

    Ntr <- x$Ntr
    Nco <- x$Nco
    id.tr_names <- colnames(Y.tr) # Names of treated units being plotted (could be empty if tr is empty)
    id.co_names <- colnames(Y.co)

    sameT0 <- length(unique(T0)) == 1 # Critical for Case 1 vs Case 2 decision

    # --- MODIFICATION START: Filter Y.tr and Y.ct for nonstaggered case (Case 1) ---
    # This modification applies if the plot will be in absolute time (Case 1).
    # It ensures that for any treated unit, after its first treatment period,
    # only periods where it is actively treated (D.tr=1) contribute to Y.tr and Y.ct.
    # Pre-treatment periods are always included.
    # The condition here matches the entry condition for Case 1 plotting logic below.
    is_case1_scenario <- ( !is.null(id) && length(id) == 1 ) || length(id.tr_names) == 1 || sameT0 == TRUE

    if (is_case1_scenario) {
      # We are modifying Y.tr and Y.ct in place.
      # These matrices are already subsetted by tr_idx_logical.
      num_treated_units_for_plot <- ncol(D.tr)
      num_time_periods_for_plot <- nrow(D.tr)

      if (num_treated_units_for_plot > 0 && num_time_periods_for_plot > 0) {
        for (j_unit_idx in 1:num_treated_units_for_plot) {
          # For each treated unit that will be part of the plot:
          # Find its first treatment period based on D.tr for this unit.
          first_treated_period_indices_for_unit <- which(D.tr[, j_unit_idx] == 1)

          if (length(first_treated_period_indices_for_unit) > 0) {
            first_ever_treated_period_row_idx <- min(first_treated_period_indices_for_unit)

            # Iterate over time periods from the first treatment onwards for this unit
            for (t_row_idx in first_ever_treated_period_row_idx:num_time_periods_for_plot) {
              # If, in this post-first-treatment period, the unit is untreated (D.tr=0 or NA)
              if (is.na(D.tr[t_row_idx, j_unit_idx]) || D.tr[t_row_idx, j_unit_idx] == 0) {
                # Set its outcome and counterfactual outcome to NA for this period
                Y.tr[t_row_idx, j_unit_idx] <- NA
                Y.ct[t_row_idx, j_unit_idx] <- NA
              }
            }
          }
          # If a unit in D.tr (i.e., a "treated" unit selected by tr_idx_logical) was never actually treated,
          # all its periods are considered pre-treatment relative to this logic, so no NA-ing here.
        }
      }
    }
    # --- MODIFICATION END ---

    # Yb (averages) must be calculated *after* Y.tr and Y.ct might have been modified.
    Yb <- cbind(apply(Y.tr, 1, mean, na.rm = TRUE), apply(Y.ct, 1, mean, na.rm = TRUE))
    colnames(Yb) <- c("Tr_Avg", "Ct_Avg")


    if (is.null(line.color)) { line.color <- if (theme.bw) "#AAAAAA70" else "white" }
    if (is.null(shade.post)) { shade.post <- TRUE }
    else if (!is.logical(shade.post)) { stop("Option \"shade.post\" must be logical (TRUE/FALSE).") }

    legend.pos <- if (legendOff) "none" else "bottom"
    user_xlim <- xlim
    validated_user_xlim <- NULL
    if (length(user_xlim) != 0) {
      if (!is.numeric(user_xlim) || length(user_xlim) != 2) {
        warning("Invalid xlim provided. It should be a numeric vector of length 2. Ignoring xlim.")
      } else {
        validated_user_xlim <- sort(user_xlim)
      }
    }

    y_data_for_range_calc <- c()

    # ================================================================
    # Case 1: Single Treated Unit or All Treated Units have Same T0 (Absolute Time)
    # ================================================================
    if (is_case1_scenario) { # Use the already determined flag

      plot_time_abs <- time_orig_vec
      time_bf_abs_val <- NA
      time_step_abs <- if(length(plot_time_abs) > 1) min(diff(sort(unique(plot_time_abs))), na.rm=TRUE) else 1
      if(!is.finite(time_step_abs) || time_step_abs <= 0) time_step_abs <- 1


      if (!is.null(id) && length(id) == 1) {
        if (!id[1] %in% id.tr_names) stop(paste(id[1], "not in the set of plotted treated units or not a valid treated unit name."))
        t0_unit_col_idx <- which(id.tr_names == id[1])
        if(length(t0_unit_col_idx) == 1 && t0_unit_col_idx <= length(T0)) {
          time_bf_period_count <- T0[t0_unit_col_idx]
          if (!is.na(time_bf_period_count) && time_bf_period_count >= 0 && (time_bf_period_count + 1) <= length(plot_time_abs)) {
            time_bf_abs_val <- plot_time_abs[time_bf_period_count + 1]
          } else if (!is.na(time_bf_period_count) && time_bf_period_count == length(plot_time_abs)) {
            time_bf_abs_val <- plot_time_abs[length(plot_time_abs)] + time_step_abs
          }
        }
      } else if (sameT0) { # This covers (is.null(id) && sameT0) including length(id.tr_names)==1
        unique_t0_val_count <- unique(T0)[1] # T0 could be empty if id.tr_names is empty, leading to issues.
        # However, if id.tr_names is empty, sameT0 is likely FALSE or T0 is length 0.
        # If T0 is length 0, unique(T0)[1] is NA.
        if(!is.na(unique_t0_val_count) && unique_t0_val_count >= 0 && (unique_t0_val_count + 1) <= length(plot_time_abs)) {
          time_bf_abs_val <- plot_time_abs[unique_t0_val_count + 1]
        } else if (!is.na(unique_t0_val_count) && unique_t0_val_count == length(plot_time_abs)) {
          time_bf_abs_val <- plot_time_abs[length(plot_time_abs)] + time_step_abs
        }
      }

      vline_pos_abs <- if (!is.na(time_bf_abs_val)) time_bf_abs_val - (time_step_abs / 2) else NA

      plot_xlim_abs <- NULL
      if (!is.null(validated_user_xlim)) {
        plot_xlim_abs <- validated_user_xlim
      } else {
        if (length(plot_time_abs) > 0) {
          plot_xlim_abs <- range(plot_time_abs, na.rm = TRUE)
          if(any(!is.finite(plot_xlim_abs))) plot_xlim_abs <- NULL
        }
      }

      show_abs <- 1:length(plot_time_abs)
      if (!is.null(plot_xlim_abs)) {
        show_abs_check <- which(plot_time_abs >= plot_xlim_abs[1] & plot_time_abs <= plot_xlim_abs[2])
        if (length(show_abs_check) == 0) {
          warning(paste0("No data points fall within the specified xlim range for absolute time. Plotting all available data."))
        } else {
          show_abs <- show_abs_check
        }
      } else if (length(plot_time_abs) > 0) {
        plot_xlim_abs <- range(plot_time_abs[show_abs], na.rm = TRUE)
        if(any(!is.finite(plot_xlim_abs))) plot_xlim_abs <- NULL
      } else {
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

      if (!is.null(id) && length(id) == 1) { # id[1] is already validated to be in id.tr_names
        plot_single_unit_flag <- TRUE
        unit_to_plot_name <- id[1]
      } else if (is.null(id) && length(id.tr_names) == 1) {
        plot_single_unit_flag <- TRUE
        unit_to_plot_name <- id.tr_names[1]
      }

      if (plot_single_unit_flag) {
        maintext <- paste("Treated and Counterfactual (", unit_to_plot_name, ")", sep = "")
        unit_col_idx_in_Y.tr <- which(colnames(Y.tr) == unit_to_plot_name) # Y.tr is potentially modified
        # This if condition should ideally not be needed if id[1] %in% id.tr_names is robustly checked
        if(length(unit_col_idx_in_Y.tr) != 1) stop(paste("Could not find unique column for unit", unit_to_plot_name, "in Y.tr (this should not happen)."))

        tr.info_unit <- Y.tr[, unit_col_idx_in_Y.tr] # Uses (potentially) modified Y.tr
        ct.info_unit <- Y.ct[, unit_col_idx_in_Y.tr] # Uses (potentially) modified Y.ct
        y_data_for_range_calc <- c(tr.info_unit[show_abs], ct.info_unit[show_abs])

        if (raw == "none") {
          data_plot_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(tr.info_unit[show_abs], ct.info_unit[show_abs]), # Will use NAs from modification
            type = factor(c(rep("tr", nT_abs), rep("ct", nT_abs)), levels = c("tr", "ct"))
          )
          p <- ggplot(data_plot_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type)) +
            geom_line()
          set.limits <- c("tr", "ct")
          set.labels <- c("Treated", "Estimated Y(0)")
          set.colors <- c(color, counterfactual.color)
          set.linetypes <- c("solid", counterfactual.linetype)
          set.linewidth <- c(line.width[1], line.width[1])

        } else if (raw == "band") {
          Y.co.quantiles <- t(apply(Y.co, 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE)) # Controls unchanged
          data_plot_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(tr.info_unit[show_abs], ct.info_unit[show_abs]), # Uses modified data
            type = factor(c(rep("tr", nT_abs), rep("ct", nT_abs)), levels = c("tr", "ct"))
          )
          data.band_abs <- data.frame(time = plot_time_abs[show_abs], co5 = Y.co.quantiles[show_abs, 1], co95 = Y.co.quantiles[show_abs, 2])
          y_data_for_range_calc <- c(y_data_for_range_calc, data.band_abs$co5, data.band_abs$co95)

          p <- ggplot() +
            geom_ribbon(data = data.band_abs, aes(x = time, ymin = co5, ymax = co95, fill = "co.band"), alpha = 0.15, color = if (ci.outline) adjustcolor(counterfactual.raw.controls.color, offset = c(0.3, 0.3, 0.3, 0)) else NA, show.legend = TRUE) +
            geom_line(data = data_plot_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))

          set.limits <- c("tr", "ct", "co.band")
          set.labels <- c("Treated", "Estimated Y(0)", "Controls (5-95% Quantiles)")
          set.colors <- c(color, counterfactual.color, NA)
          set.linetypes <- c("solid", counterfactual.linetype, "blank")
          set.linewidth <- c(line.width[1], line.width[1], 0)
          set.fill <- c(NA, NA, counterfactual.raw.controls.color)

        } else if (raw == "all") {
          Y.co.subset <- Y.co[show_abs, , drop = FALSE] # Controls unchanged
          raw_co_data_abs <- NULL
          if (ncol(Y.co.subset) > 0) {
            melt_temp_co <- reshape2::melt(Y.co.subset, varnames = c("time_idx_abs", "id_co_idx_abs"), value.name = "outcome")
            if(nrow(melt_temp_co) > 0) {
              raw_co_data_abs <- data.frame(
                time = plot_time_abs[show_abs][melt_temp_co$time_idx_abs],
                id = id.co_names[melt_temp_co$id_co_idx_abs],
                outcome = melt_temp_co$outcome,
                type = "raw.co"
              )
            }
          }

          main_lines_data_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(tr.info_unit[show_abs], ct.info_unit[show_abs]), # Uses modified data
            type = factor(c(rep("tr", nT_abs), rep("ct", nT_abs)), levels=c("tr", "ct", "raw.co")),
            id = c(rep(paste0("tr_", unit_to_plot_name), nT_abs), rep(paste0("ct_", unit_to_plot_name), nT_abs))
          )

          plot_data_list_abs <- list(main_lines_data_abs[, c("time", "outcome", "type", "id")])
          if (!is.null(raw_co_data_abs)) plot_data_list_abs[[2]] <- raw_co_data_abs[, c("time", "outcome", "type", "id")]
          plot_data_combined_abs <- do.call(rbind, plot_data_list_abs)
          y_data_for_range_calc <- plot_data_combined_abs$outcome

          p <- ggplot(plot_data_combined_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type, group = id)) +
            geom_line()

          set.limits <- c("tr", "ct", "raw.co")
          set.labels <- c(paste0("Treated (",unit_to_plot_name,")"), paste0("Est. Y(0) (",unit_to_plot_name,")"), "Controls")
          set.colors <- c(color, counterfactual.color, counterfactual.raw.controls.color)
          set.linetypes <- c("solid", counterfactual.linetype, "solid")
          lw <- if(length(line.width) >= 2) line.width else rep(line.width[1], 2)
          set.linewidth <- c(lw[1], lw[1], lw[2])
        }
      } else { # Multiple treated units, all with the same T0 (plot_single_unit_flag is FALSE)
        maintext <- "Treated and Counterfactual Averages"
        Yb_show_abs <- Yb[show_abs, , drop = FALSE] # Yb is already based on modified Y.tr, Y.ct
        y_data_for_range_calc <- c(Yb_show_abs[,1], Yb_show_abs[,2])

        if (raw == "none") {
          data_plot_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(Yb_show_abs[, 1], Yb_show_abs[, 2]), # Uses modified averages
            type = factor(c(rep("tr", nT_abs), rep("co", nT_abs)), levels = c("tr", "co"))
          )
          p <- ggplot(data_plot_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))

          if (plot.ci %in% c("90", "95") && !is.null(x$Y.avg)) {
            # CI data from x$Y.avg is based on original model estimates, not this plot-specific filtering.
            # This might lead to CIs appearing around lines that have NAs due to filtering. This is a known behavior.
            tr_lo_col <- if (plot.ci == "95") "lower.tr" else "lower90.tr"
            tr_hi_col <- if (plot.ci == "95") "upper.tr" else "upper90.tr"
            cf_lo_col <- if (plot.ci == "95") "lower.ct" else "lower90.ct"
            cf_hi_col <- if (plot.ci == "95") "upper.ct" else "upper90.ct"
            required_cols_ci <- c("period", tr_lo_col, tr_hi_col, cf_lo_col, cf_hi_col)

            if (all(required_cols_ci %in% names(x$Y.avg))) {
              ci_data_abs <- x$Y.avg
              ci_data_filtered_abs <- NULL
              if (!is.null(plot_xlim_abs) && is.numeric(ci_data_abs$period)) {
                ci_data_filtered_abs <- ci_data_abs[ci_data_abs$period >= plot_xlim_abs[1] & ci_data_abs$period <= plot_xlim_abs[2], ]
              } else {
                ci_data_filtered_abs <- ci_data_abs[ci_data_abs$period %in% plot_time_abs[show_abs],]
              }

              if (!is.null(ci_data_filtered_abs) && nrow(ci_data_filtered_abs) > 0) {
                p <- p +
                  geom_ribbon(data = ci_data_filtered_abs, aes(x = period, ymin = .data[[tr_lo_col]], ymax = .data[[tr_hi_col]], fill = "tr"), alpha = 0.2, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE) +
                  geom_ribbon(data = ci_data_filtered_abs, aes(x = period, ymin = .data[[cf_lo_col]], ymax = .data[[cf_hi_col]], fill = "co"), alpha = 0.2, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE)
                y_data_for_range_calc <- c(y_data_for_range_calc, ci_data_filtered_abs[[tr_lo_col]], ci_data_filtered_abs[[tr_hi_col]], ci_data_filtered_abs[[cf_lo_col]], ci_data_filtered_abs[[cf_hi_col]])
              } else { warning("No CI data points fall within the specified xlim range for absolute time CI.") }
            } else { warning("CI columns not found in x$Y.avg for absolute time CI. Skipping CI plotting.") }
          }
          p <- p + geom_line()

          set.limits <- c("tr", "co")
          set.labels <- c("Treated Average", "Estimated Y(0) Average")
          set.colors <- c(color, counterfactual.color)
          set.linetypes <- c("solid", counterfactual.linetype)
          set.linewidth <- c(line.width[1], line.width[1])
          set.fill <- c(color, counterfactual.color)

        } else if (raw == "band") {
          # Y.tr.quantiles must use the (potentially) modified Y.tr
          Y.tr.quantiles <- t(apply(Y.tr, 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE)) # Y.tr is modified
          Y.co.quantiles <- t(apply(Y.co, 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE)) # Controls unchanged
          quantiles_all_abs <- cbind(Y.tr.quantiles, Y.co.quantiles)
          data.band_abs <- data.frame(time = plot_time_abs[show_abs], tr5 = quantiles_all_abs[show_abs, 1], tr95 = quantiles_all_abs[show_abs, 2], co5 = quantiles_all_abs[show_abs, 3], co95 = quantiles_all_abs[show_abs, 4])
          data_plot_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(Yb_show_abs[, 1], Yb_show_abs[, 2]), # Uses modified averages
            type = factor(c(rep("tr", nT_abs), rep("co", nT_abs)), levels = c("tr", "co"))
          )
          y_data_for_range_calc <- c(y_data_for_range_calc, data.band_abs$tr5, data.band_abs$tr95, data.band_abs$co5, data.band_abs$co95)

          p <- ggplot() +
            geom_ribbon(data = data.band_abs, aes(x = time, ymin = tr5, ymax = tr95, fill = "tr_band"), alpha = 0.15, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
            geom_ribbon(data = data.band_abs, aes(x = time, ymin = co5, ymax = co95, fill = "co_band"), alpha = 0.15, color = if (ci.outline) adjustcolor(counterfactual.raw.controls.color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
            geom_line(data = data_plot_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))

          set.limits <- c("tr", "co", "tr_band", "co_band")
          set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated (5-95% Quantiles)", "Controls (5-95% Quantiles)")
          set.colors <- c(color, counterfactual.color, NA, NA)
          set.linetypes <- c("solid", counterfactual.linetype, "blank", "blank")
          set.linewidth <- c(line.width[1], line.width[1], 0, 0)
          set.fill <- c(NA, NA, color, counterfactual.raw.controls.color)

        } else if (raw == "all") {
          # Y.tr.subset must use the (potentially) modified Y.tr
          Y.tr.subset <- Y.tr[show_abs, , drop = FALSE] # Y.tr is modified
          Y.co.subset <- Y.co[show_abs, , drop = FALSE] # Controls unchanged
          plot_data_list_abs <- list()

          avg_lines_data_abs <- data.frame(
            time = rep(plot_time_abs[show_abs], 2),
            outcome = c(Yb_show_abs[, 1], Yb_show_abs[, 2]), # Uses modified averages
            type = factor(c(rep("tr", nT_abs), rep("co", nT_abs)), levels=c("tr", "co", "raw.tr", "raw.co")),
            id = c(rep("tr_avg_line", nT_abs), rep("co_avg_line", nT_abs))
          )
          plot_data_list_abs[[1]] <- avg_lines_data_abs[, c("time", "outcome", "type", "id")]

          raw_tr_data_abs <- NULL
          if (ncol(Y.tr.subset) > 0) { # Y.tr.subset is from modified Y.tr
            melt_temp_tr <- reshape2::melt(Y.tr.subset, varnames = c("time_idx_abs", "id_tr_idx_abs"), value.name = "outcome")
            if(nrow(melt_temp_tr) > 0) {
              raw_tr_data_abs <- data.frame(
                time = plot_time_abs[show_abs][melt_temp_tr$time_idx_abs],
                id = id.tr_names[melt_temp_tr$id_tr_idx_abs], # id.tr_names are colnames of original Y.tr
                outcome = melt_temp_tr$outcome,
                type = "raw.tr"
              )
              plot_data_list_abs[[length(plot_data_list_abs) + 1]] <- raw_tr_data_abs[, c("time", "outcome", "type", "id")]
            }
          }
          raw_co_data_abs <- NULL
          if (ncol(Y.co.subset) > 0) {
            melt_temp_co <- reshape2::melt(Y.co.subset, varnames = c("time_idx_abs", "id_co_idx_abs"), value.name = "outcome")
            if(nrow(melt_temp_co) > 0) {
              raw_co_data_abs <- data.frame(
                time = plot_time_abs[show_abs][melt_temp_co$time_idx_abs],
                id = id.co_names[melt_temp_co$id_co_idx_abs],
                outcome = melt_temp_co$outcome,
                type = "raw.co"
              )
              plot_data_list_abs[[length(plot_data_list_abs) + 1]] <- raw_co_data_abs[, c("time", "outcome", "type", "id")]
            }
          }
          plot_data_combined_abs <- do.call(rbind, plot_data_list_abs)
          y_data_for_range_calc <- plot_data_combined_abs$outcome

p <- ggplot(plot_data_combined_abs, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type, group = id)) +
            geom_line()

          set.limits <- c("tr", "co", "raw.tr", "raw.co")
          set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated Raw Data", "Controls Raw Data")
          set.colors <- c(color, counterfactual.color, counterfactual.raw.treated.color, counterfactual.raw.controls.color)
          set.linetypes <- c("solid", counterfactual.linetype, "solid", "solid")
          lw <- if(length(line.width) >= 2) line.width else rep(line.width[1], 2)
          set.linewidth <- c(lw[1], lw[1], lw[2], lw[2])
        }
      }

      p <- p + xlab(xlab_final) + ylab(ylab_final) +
        theme(
          legend.position = legend.pos,
          axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
          axis.text = element_text(size = cex.axis),
          axis.title = element_text(size = cex.lab),
          plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
        )
      if (theme.bw == TRUE) { p <- p + theme_bw() }

      if (!is.na(vline_pos_abs)) {
        p <- p + geom_vline(xintercept = vline_pos_abs, colour = line.color, linewidth = lwidth, linetype = "dashed")
        if (shade.post == TRUE) {
          p <- p + annotate("rect", xmin = vline_pos_abs, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = .3)
        }
      }

      p <- p + scale_colour_manual(name = NULL, limits = set.limits, values = set.colors, labels = set.labels, na.value = NA) +
        scale_linetype_manual(name = NULL, limits = set.limits, values = set.linetypes, labels = set.labels, na.value = "blank") +
        scale_linewidth_manual(name = NULL, limits = set.limits, values = set.linewidth, labels = set.labels, na.value = 0)

      guide_obj_abs <- guide_legend(title = NULL, ncol = if(length(set.limits) > 2) ceiling(length(set.limits)/2) else length(set.limits))

      if (exists("set.fill") && !is.null(set.fill)) {
        p <- p + scale_fill_manual(name = NULL, limits = set.limits, values = set.fill, labels = set.labels, na.value = NA)
        p <- p + guides(colour = guide_obj_abs, linetype = guide_obj_abs, size = guide_obj_abs, fill = guide_obj_abs)
      } else {
        p <- p + guides(colour = guide_obj_abs, linetype = guide_obj_abs, size = guide_obj_abs)
      }

      if (length(T.b_abs) > 0) {
        if (!is.numeric(time_label_abs)) {
          p <- p + scale_x_continuous(expand = c(0.01, 0.01), breaks = plot_time_abs[show_abs][T.b_abs], labels = time_label_abs[T.b_abs])
        } else {
          p <- p + scale_x_continuous(expand = c(0.01, 0.01), labels = scaleFUN, breaks = plot_time_abs[show_abs][T.b_abs])
        }
      } else {
        p <- p + scale_x_continuous(expand = c(0.01, 0.01))
      }

      final_main_text <- if (is.null(main)) maintext else if (main == "") NULL else main
      if (!is.null(final_main_text)) p <- p + ggtitle(final_main_text)

      if (show.count == TRUE) {
        counts_values_abs <- NULL
        current_times_abs <- plot_time_abs[show_abs]

        if (plot_single_unit_flag) {
          # Count based on non-NA values in the (potentially) modified Y.tr for the specific unit
          unit_col_idx_for_count <- which(colnames(Y.tr) == unit_to_plot_name) # Y.tr is (potentially) modified
          counts_values_abs <- ifelse(!is.na(Y.tr[show_abs, unit_col_idx_for_count, drop = FALSE]), 1, 0)
        } else {
          # Count based on non-NA values in the (potentially) modified Y.tr for all units in this group
          counts_values_abs <- apply(Y.tr[show_abs, , drop=FALSE], 1, function(row) sum(!is.na(row))) # Y.tr is (potentially) modified
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

            rect_min_val_abs <- NULL
            if (!is.null(ylim)) {
              rect_min_val_abs <- ylim[1]
            } else {
              rect_min_val_abs <- current_plot_yrange[1] - count_bar_space_height_abs
            }

            counts_for_plot_df_abs$ymin <- rect_min_val_abs
            counts_for_plot_df_abs$ymax <- rect_min_val_abs + (counts_for_plot_df_abs$count / max_count_val_abs) * actual_rect_length_abs

            time_step_for_bars_abs <- if(length(unique(counts_for_plot_df_abs$time)) > 1) min(diff(sort(unique(counts_for_plot_df_abs$time))),na.rm=TRUE) else 1
            if(!is.finite(time_step_for_bars_abs) || time_step_for_bars_abs <=0) time_step_for_bars_abs <- 1
            bar_width_half_abs <- time_step_for_bars_abs * 0.20

            counts_for_plot_df_abs$xmin <- counts_for_plot_df_abs$time - bar_width_half_abs
            counts_for_plot_df_abs$xmax <- counts_for_plot_df_abs$time + bar_width_half_abs

            max_count_time_pos_abs <- counts_for_plot_df_abs$time[which.max(counts_for_plot_df_abs$count)[1]]
            text_y_pos_abs <- rect_min_val_abs + actual_rect_length_abs + (count_bar_space_height_abs - actual_rect_length_abs) * 0.5

            p <- p + geom_rect(data = counts_for_plot_df_abs,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                               fill = count.color, inherit.aes = FALSE)
            p <- p + annotate("text", x = max_count_time_pos_abs, y = text_y_pos_abs,
                              label = max_count_val_abs, size = cex.text * 0.8, hjust = 0.5, vjust = 0.5)
          }
        }
      }

      coord_args_abs <- list(clip = "on")
      if (!is.null(ylim)) { coord_args_abs$ylim <- ylim }
      else if (show.count == TRUE && exists("rect_min_val_abs") && exists("counts_for_plot_df_abs") && nrow(counts_for_plot_df_abs) > 0 && !is.null(current_plot_yrange)) { # Added !is.null(current_plot_yrange)
        final_yrange_min <- min(current_plot_yrange[1], rect_min_val_abs, na.rm=TRUE)
        final_yrange_max <- current_plot_yrange[2]
        coord_args_abs$ylim <- c(final_yrange_min, final_yrange_max)
      }
      if (!is.null(plot_xlim_abs)) { coord_args_abs$xlim <- plot_xlim_abs }
      if (length(coord_args_abs) > 1) {
        p <- p + do.call(coord_cartesian, coord_args_abs)
      }
      p <- p + theme(legend.position = legend.pos)

      # ================================================================
      # Case 2: Staggered Adoption (Different T0) -> Relative Time Plot
      # ================================================================
    } else { # This is the staggered case (relative time plot, is_case1_scenario is FALSE)
      # Y.tr and Y.ct passed to ct.adjust are original, as modification block was skipped.
      maintext <- "Treated and Counterfactual Averages (Relative Time)"
      xlab_final <- if (is.null(xlab)) "Time relative to Treatment" else if (xlab == "") NULL else xlab
      ylab_final <- if (is.null(ylab)) x$Yname else if (ylab == "") NULL else ylab

      ct.adjust <- function(Y.tr, Y.ct, T0_counts, D_full, tr_idx_logical_full) {
        if (!is.matrix(Y.tr) || !is.matrix(Y.ct) || !is.vector(T0_counts)) stop("Y.tr, Y.ct must be matrices, T0_counts must be a vector.")
        if (!all(dim(Y.tr) == dim(Y.ct))) stop("Dimensions of Y.tr and Y.ct must match.")

        T_rows <- nrow(Y.tr)
        N_tr <- ncol(Y.tr)

        if (length(T0_counts) != N_tr) stop("Length of T0_counts must match N_tr.")
        if (T_rows == 0 || N_tr == 0) {
          warning("Input matrices for ct.adjust have 0 rows or 0 columns.")
          return(list(timeline = integer(0), Y.tr.aug = matrix(NA, 0, N_tr), Y.ct.aug = matrix(NA, 0, N_tr), Yb = matrix(NA, 0, 2, dimnames=list(NULL, c("Y.tr.bar", "Y.ct.bar")))))
        }

        D_tr_local <- D_full[, tr_idx_logical_full, drop = FALSE]
        if (ncol(D_tr_local) != N_tr) stop("Dimension mismatch for D_tr_local in ct.adjust.")


        T0_counts_valid <- T0_counts[!is.na(T0_counts) & is.finite(T0_counts)]
        if (length(T0_counts_valid) == 0) {
          warning("No valid T0_counts values found.")
          return(list(timeline = integer(0), Y.tr.aug = matrix(NA,0,N_tr), Y.ct.aug = matrix(NA,0,N_tr), Yb = matrix(NA,0,2)))
        }

        min_rel_time <- 1 - max(T0_counts_valid)
        max_rel_time <- T_rows - min(T0_counts_valid)
        if (min_rel_time > max_rel_time) {
          timeline <- integer(0)
        } else {
          timeline <- min_rel_time:max_rel_time
        }


        Y.tr.aug <- matrix(NA_real_, nrow = length(timeline), ncol = N_tr)
        Y.ct.aug <- matrix(NA_real_, nrow = length(timeline), ncol = N_tr)
        if(length(timeline)>0) rownames(Y.tr.aug) <- rownames(Y.ct.aug) <- timeline

        for (i in seq_len(N_tr)) {
          unit_t0_count_val <- T0_counts[i]
          if (is.na(unit_t0_count_val) || !is.finite(unit_t0_count_val)) next

          for (t_abs_idx in seq_len(T_rows)) {
            rel_time <- t_abs_idx - unit_t0_count_val
            row_idx_aug <- match(rel_time, timeline)
            if (is.na(row_idx_aug)) next

            include_value <- FALSE
            if (rel_time <= 0) {
              include_value <- TRUE
            } else {
              if (t_abs_idx <= nrow(D_tr_local) && i <= ncol(D_tr_local)) {
                if (!is.na(D_tr_local[t_abs_idx, i]) && D_tr_local[t_abs_idx, i] == 1) {
                  include_value <- TRUE
                }
              }
            }

            if (include_value) {
              if (t_abs_idx <= nrow(Y.tr) && t_abs_idx <= nrow(Y.ct)) {
                Y.tr.aug[row_idx_aug, i] <- Y.tr[t_abs_idx, i]
                Y.ct.aug[row_idx_aug, i] <- Y.ct[t_abs_idx, i]
              }
            }
          }
        }

        Y.tr.bar <- if(length(timeline)>0) apply(Y.tr.aug, 1, mean, na.rm = TRUE) else numeric(0)
        Y.ct.bar <- if(length(timeline)>0) apply(Y.ct.aug, 1, mean, na.rm = TRUE) else numeric(0)
        if(length(Y.tr.bar)>0) Y.tr.bar[is.nan(Y.tr.bar)] <- NA_real_
        if(length(Y.ct.bar)>0) Y.ct.bar[is.nan(Y.ct.bar)] <- NA_real_

        Yb_rel <- cbind(Y.tr.bar, Y.ct.bar)
        if(length(timeline)>0 && ncol(Yb_rel)>0) colnames(Yb_rel) <- c("Y.tr.bar", "Y.ct.bar")

        return(list(timeline = timeline, Y.tr.aug = Y.tr.aug, Y.ct.aug = Y.ct.aug, Yb = Yb_rel))
      }
      # --- End ct.adjust definition ---

      xx <- ct.adjust(Y.tr, Y.ct, T0, D_mat, tr_idx_logical)  # D_mat is the full D matrix, tr_idx_logical for full D

      plot_time_rel <- xx$timeline
      if (length(plot_time_rel) == 0 || all(is.na(plot_time_rel))) {
        stop("Relative timeline calculation resulted in zero length or all NA. Cannot plot.")
      }

      plot_xlim_rel <- NULL
      if (!is.null(validated_user_xlim)) {
        plot_xlim_rel <- validated_user_xlim
        show_check_rel <- which(plot_time_rel >= plot_xlim_rel[1] & plot_time_rel <= plot_xlim_rel[2])
        if (length(show_check_rel) == 0) {
          warning(paste0("User-provided xlim for relative time contains no points. Plotting full relative time range."))
          plot_xlim_rel <- range(plot_time_rel, na.rm = TRUE)
          if(any(!is.finite(plot_xlim_rel))) plot_xlim_rel <- NULL
        }
      } else {
        plot_xlim_rel <- range(plot_time_rel, na.rm = TRUE)
        if(any(!is.finite(plot_xlim_rel))) plot_xlim_rel <- NULL
      }

      show_rel <- 1:length(plot_time_rel)
      if (!is.null(plot_xlim_rel)) {
        show_rel_check <- which(plot_time_rel >= plot_xlim_rel[1] & plot_time_rel <= plot_xlim_rel[2])
        if (length(show_rel_check) == 0) {
          warning("No relative time data points fall within the calculated xlim range. Plotting all available relative time data.")
        } else {
          show_rel <- show_rel_check
        }
      } else if (length(plot_time_rel) > 0) {
        plot_xlim_rel <- range(plot_time_rel[show_rel], na.rm = TRUE)
        if(any(!is.finite(plot_xlim_rel))) plot_xlim_rel <- NULL
      } else {
        stop("Cannot determine relative time range for plot.")
      }

      nT_rel <- length(show_rel)
      time_label_rel <- plot_time_rel[show_rel]
      time_bf_rel_vline <- 0.5

      T.b_rel <- integer(0)
      if (nT_rel > 0) {
        if (is.numeric(time_label_rel) && length(time_label_rel) > 1) {
          T.b_breaks_rel <- pretty(time_label_rel)
          T.b_breaks_rel <- T.b_breaks_rel[T.b_breaks_rel >= min(time_label_rel) & T.b_breaks_rel <= max(time_label_rel)]
          if (length(T.b_breaks_rel) > 0) {
            T.b_rel <- sapply(T.b_breaks_rel, function(br) which.min(abs(time_label_rel - br)))
            T.b_rel <- unique(T.b_rel)
          }
        }
        if (length(T.b_rel) == 0) {
          max_labels <- 10
          step <- max(1, floor(length(time_label_rel) / max_labels))
          T.b_rel <- seq(1, length(time_label_rel), by = step)
        }
      }

      Yb_show_rel <- xx$Yb[show_rel, , drop = FALSE]
      Ytr_aug_show_rel <- xx$Y.tr.aug[show_rel, , drop = FALSE]
      y_data_for_range_calc <- c(Yb_show_rel[,1], Yb_show_rel[,2])

      if (raw == "none") {
        data_plot_rel <- data.frame(
          time = rep(plot_time_rel[show_rel], 2),
          outcome = c(Yb_show_rel[, 1], Yb_show_rel[, 2]),
          type = factor(c(rep("tr", nT_rel), rep("co", nT_rel)), levels = c("tr", "co"))
        )
        p <- ggplot(data_plot_rel, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))

        if (plot.ci %in% c("90", "95") && !is.null(x$Y.avg)) {
          tr_lo_col <- if (plot.ci == "95") "lower.tr" else "lower90.tr"
          tr_hi_col <- if (plot.ci == "95") "upper.tr" else "upper90.tr"
          cf_lo_col <- if (plot.ci == "95") "lower.ct" else "lower90.ct"
          cf_hi_col <- if (plot.ci == "95") "upper.ct" else "upper90.ct"
          required_cols_ci <- c("period", tr_lo_col, tr_hi_col, cf_lo_col, cf_hi_col)

          if (all(required_cols_ci %in% names(x$Y.avg))) {
            ci_data_rel <- x$Y.avg
            ci_data_filtered_rel <- NULL
            if (!is.null(plot_xlim_rel) && is.numeric(ci_data_rel$period)) {
              ci_data_filtered_rel <- ci_data_rel[ci_data_rel$period >= plot_xlim_rel[1] & ci_data_rel$period <= plot_xlim_rel[2], ]
            } else {
              ci_data_filtered_rel <- ci_data_rel[ci_data_rel$period %in% plot_time_rel[show_rel],]
            }

            if (!is.null(ci_data_filtered_rel) && nrow(ci_data_filtered_rel) > 0) {
              p <- p +
                geom_ribbon(data = ci_data_filtered_rel, aes(x = period, ymin = .data[[tr_lo_col]], ymax = .data[[tr_hi_col]], fill = "tr"), alpha = 0.2, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE) +
                geom_ribbon(data = ci_data_filtered_rel, aes(x = period, ymin = .data[[cf_lo_col]], ymax = .data[[cf_hi_col]], fill = "co"), alpha = 0.2, color = if (ci.outline) adjustcolor(counterfactual.color, offset = c(0.3, 0.3, 0.3, 0)) else NA, inherit.aes = FALSE)
              y_data_for_range_calc <- c(y_data_for_range_calc, ci_data_filtered_rel[[tr_lo_col]], ci_data_filtered_rel[[tr_hi_col]], ci_data_filtered_rel[[cf_lo_col]], ci_data_filtered_rel[[cf_hi_col]])
            } else { warning("No CI data points fall within the specified relative xlim range for CI.") }
          } else { warning("Relative time CI columns not found in x$Y.avg. Skipping CI plotting.") }
        }
        p <- p + geom_line()

        set.limits <- c("tr", "co")
        set.labels <- c("Treated Average", "Estimated Y(0) Average")
        set.colors <- c(color, counterfactual.color)
        set.linetypes <- c("solid", counterfactual.linetype)
        set.linewidth <- c(line.width[1], line.width[1])
        set.fill <- c(color, counterfactual.color)

      } else if (raw == "band") {
        Y.tr.aug.quantiles_rel <- t(apply(Ytr_aug_show_rel, 1, quantile, prob = c(0.05, 0.95), na.rm = TRUE))
        data_plot_rel <- data.frame(
          time = rep(plot_time_rel[show_rel], 2),
          outcome = c(Yb_show_rel[, 1], Yb_show_rel[, 2]),
          type = factor(c(rep("tr", nT_rel), rep("co", nT_rel)), levels = c("tr", "co"))
        )
        data.band_rel <- data.frame(time = plot_time_rel[show_rel], tr5 = Y.tr.aug.quantiles_rel[, 1], tr95 = Y.tr.aug.quantiles_rel[, 2])
        y_data_for_range_calc <- c(y_data_for_range_calc, data.band_rel$tr5, data.band_rel$tr95)


        p <- ggplot() +
          geom_ribbon(data = data.band_rel, aes(x = time, ymin = tr5, ymax = tr95, fill = "tr_band"), alpha = 0.15, color = if (ci.outline) adjustcolor(color, offset = c(0.3, 0.3, 0.3, 0)) else NA) +
          geom_line(data = data_plot_rel, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type))

        set.limits <- c("tr", "co", "tr_band")
        set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated (5-95% Quantiles)")
        set.colors <- c(color, counterfactual.color, NA)
        set.linetypes <- c("solid", counterfactual.linetype, "blank")
        set.linewidth <- c(line.width[1], line.width[1], 0)
        set.fill <- c(NA, NA, color)

      } else if (raw == "all") {
        plot_data_list_rel <- list()
        avg_lines_data_rel <- data.frame(
          time = rep(plot_time_rel[show_rel], 2),
          outcome = c(Yb_show_rel[, 1], Yb_show_rel[, 2]),
          type = factor(c(rep("tr", nT_rel), rep("co", nT_rel)), levels=c("tr", "co", "raw.tr")),
          id = c(rep("tr_avg_line", nT_rel), rep("co_avg_line", nT_rel))
        )
        plot_data_list_rel[[1]] <- avg_lines_data_rel[, c("time", "outcome", "type", "id")]

        if (ncol(Ytr_aug_show_rel) > 0 && nrow(Ytr_aug_show_rel) > 0) {
          melt_temp_tr_rel <- reshape2::melt(Ytr_aug_show_rel, varnames = c("time_idx_rel", "id_tr_idx_rel"), value.name = "outcome")
          if(nrow(melt_temp_tr_rel) > 0) {
            raw_tr_data_rel <- data.frame(
              time = plot_time_rel[show_rel][melt_temp_tr_rel$time_idx_rel],
              id = colnames(Y.tr)[melt_temp_tr_rel$id_tr_idx_rel], # Y.tr colnames are original treated unit names
              outcome = melt_temp_tr_rel$outcome,
              type = "raw.tr"
            )
            plot_data_list_rel[[length(plot_data_list_rel) + 1]] <- raw_tr_data_rel[, c("time", "outcome", "type", "id")]
          }
        }
        plot_data_combined_rel <- do.call(rbind, plot_data_list_rel)
        y_data_for_range_calc <- plot_data_combined_rel$outcome


        p <- ggplot(plot_data_combined_rel, aes(x = time, y = outcome, colour = type, linetype = type, linewidth = type, group = id)) +
          geom_line()

        set.limits <- c("tr", "co", "raw.tr")
        set.labels <- c("Treated Average", "Estimated Y(0) Average", "Treated Raw Data")
        set.colors <- c(color, counterfactual.color, counterfactual.raw.treated.color)
        set.linetypes <- c("solid", counterfactual.linetype, "solid")
        lw <- if(length(line.width) >= 2) line.width else rep(line.width[1], 2)
        set.linewidth <- c(lw[1], lw[1], lw[2])
      }


      p <- p + xlab(xlab_final) + ylab(ylab_final) +
        geom_vline(xintercept = time_bf_rel_vline, colour = line.color, linewidth = lwidth, linetype = "dashed") +
        theme(
          legend.position = legend.pos,
          axis.text.x = element_text(angle = angle, hjust = x.h, vjust = x.v),
          axis.text = element_text(size = cex.axis),
          axis.title = element_text(size = cex.lab),
          plot.title = element_text(size = cex.main, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0))
        )
      if (theme.bw == TRUE) { p <- p + theme_bw() }
      if (shade.post == TRUE) {
        p <- p + annotate("rect", xmin = time_bf_rel_vline, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = .3)
      }

      p <- p + scale_colour_manual(name = NULL, limits = set.limits, values = set.colors, labels = set.labels, na.value = NA) +
        scale_linetype_manual(name = NULL, limits = set.limits, values = set.linetypes, labels = set.labels, na.value = "blank") +
        scale_linewidth_manual(name = NULL, limits = set.limits, values = set.linewidth, labels = set.labels, na.value = 0)

      guide_obj_rel <- guide_legend(title = NULL, ncol = if(length(set.limits) > 2) ceiling(length(set.limits)/2) else length(set.limits))

      if (exists("set.fill") && !is.null(set.fill)) {
        p <- p + scale_fill_manual(name = NULL, limits = set.limits, values = set.fill, labels = set.labels, na.value = NA)
        p <- p + guides(colour = guide_obj_rel, linetype = guide_obj_rel, size = guide_obj_rel, fill = guide_obj_rel)
      } else {
        p <- p + guides(colour = guide_obj_rel, linetype = guide_obj_rel, size = guide_obj_rel)
      }

      if (length(T.b_rel) > 0) {
        p <- p + scale_x_continuous(expand = c(0.01, 0.01), labels = scaleFUN, breaks = plot_time_rel[show_rel][T.b_rel])
      } else {
        p <- p + scale_x_continuous(expand = c(0.01, 0.01), labels = scaleFUN)
      }

      final_main_text <- if (is.null(main)) maintext else if (main == "") NULL else main
      if (!is.null(final_main_text)) p <- p + ggtitle(final_main_text)

      if (show.count == TRUE) {
        current_times_rel <- plot_time_rel[show_rel]
        counts_values_rel <- NULL
        if (nrow(xx$Y.tr.aug) > 0 && ncol(xx$Y.tr.aug) > 0) {
          counts_values_rel <- rowSums(!is.na(xx$Y.tr.aug[show_rel, , drop = FALSE]), na.rm = TRUE)
        } else {
          counts_values_rel <- rep(0, length(current_times_rel))
        }

        counts_for_plot_df_rel <- data.frame(time = current_times_rel, count = counts_values_rel)
        counts_for_plot_df_rel <- counts_for_plot_df_rel[counts_for_plot_df_rel$count > 0 & !is.na(counts_for_plot_df_rel$count), ]

        if (nrow(counts_for_plot_df_rel) > 0) {
          max_count_val_rel <- max(counts_for_plot_df_rel$count, na.rm = TRUE)
          if (max_count_val_rel > 0 && !is.na(max_count_val_rel)) {
            current_plot_yrange_rel <- NULL
            if (!is.null(ylim)) {
              current_plot_yrange_rel <- ylim
            } else {
              if (length(y_data_for_range_calc) == 0 || all(is.na(y_data_for_range_calc))) {
                gb_rel <- ggplot_build(p)
                current_plot_yrange_rel <- gb_rel$layout$panel_scales_y[[1]]$range$range
              } else {
                current_plot_yrange_rel <- range(y_data_for_range_calc, na.rm = TRUE)
              }
            }

            count_bar_space_prop <- 0.20
            count_bar_space_height_rel <- (current_plot_yrange_rel[2] - current_plot_yrange_rel[1]) * count_bar_space_prop
            actual_rect_length_rel <- count_bar_space_height_rel * 0.8

            rect_min_val_rel <- NULL
            if (!is.null(ylim)) {
              rect_min_val_rel <- ylim[1]
            } else {
              rect_min_val_rel <- current_plot_yrange_rel[1] - count_bar_space_height_rel
            }

            counts_for_plot_df_rel$ymin <- rect_min_val_rel
            counts_for_plot_df_rel$ymax <- rect_min_val_rel + (counts_for_plot_df_rel$count / max_count_val_rel) * actual_rect_length_rel

            time_step_for_bars_rel <- if(length(unique(counts_for_plot_df_rel$time)) > 1) min(diff(sort(unique(counts_for_plot_df_rel$time))), na.rm=TRUE) else 1
            if(!is.finite(time_step_for_bars_rel) || time_step_for_bars_rel <=0) time_step_for_bars_rel <- 1
            bar_width_half_rel <- time_step_for_bars_rel * 0.20

            counts_for_plot_df_rel$xmin <- counts_for_plot_df_rel$time - bar_width_half_rel
            counts_for_plot_df_rel$xmax <- counts_for_plot_df_rel$time + bar_width_half_rel

            max_count_time_pos_rel <- counts_for_plot_df_rel$time[which.max(counts_for_plot_df_rel$count)[1]]
            text_y_pos_rel <- rect_min_val_rel + actual_rect_length_rel + (count_bar_space_height_rel - actual_rect_length_rel) * 0.5

            p <- p + geom_rect(data = counts_for_plot_df_rel,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                               fill = count.color, inherit.aes = FALSE)
            p <- p + annotate("text", x = max_count_time_pos_rel, y = text_y_pos_rel,
                              label = max_count_val_rel, size = cex.text * 0.8, hjust = 0.5, vjust = 0.5)
          }
        }
      }

      coord_args_rel <- list(clip = "on")
      if (!is.null(ylim)) { coord_args_rel$ylim <- ylim }
      else if (show.count == TRUE && exists("rect_min_val_rel") && exists("counts_for_plot_df_rel") && nrow(counts_for_plot_df_rel) > 0 && !is.null(current_plot_yrange_rel)) { # Added !is.null
        final_yrange_min_rel <- min(current_plot_yrange_rel[1], rect_min_val_rel, na.rm=TRUE)
        final_yrange_max_rel <- current_plot_yrange_rel[2]
        coord_args_rel$ylim <- c(final_yrange_min_rel, final_yrange_max_rel)
      }
      if (!is.null(plot_xlim_rel)) { coord_args_rel$xlim <- plot_xlim_rel }
      if (length(coord_args_rel) > 1) {
        p <- p + do.call(coord_cartesian, coord_args_rel)
      }
      p <- p + theme(legend.position = legend.pos)

    }

    if (gridOff == TRUE) {
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    if (exists("set.fill")) rm(set.fill) # clean up set.fill if it was created

    return(p)

  } # End of main if statement for type == "counterfactual"

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
    if(loo==TRUE & x$loo==TRUE){
      time.loo <- as.numeric(rownames(x$pre.est.att))
      time.post <- as.numeric(rownames(x$est.att))
      time.post <- time.post[which(time.post>0)]
      time <- sort(c(time.loo,time.post))
      count.num <- c(x$pre.est.att[,'count.on'],x$est.att[which(as.numeric(rownames(x$est.att))>0),'count'])
      best.pos <- 1
      max.count <- max(count.num)
    }else{
      time <- x$time
      count.num <- x$count
      best.pos <- 1
      max.count <- max(count.num)
    }
  }
  else if (type == "exit") {
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
        xlab <- paste("Time Since the Treatment’s Onset")
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

    ## y=0 line type
    if (is.null(lcolor) == TRUE) {
      lcolor <- "white"
      if (theme.bw == TRUE) {
        lcolor <- "#aaaaaa"
      }
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
        tb   <- est.att
        data <- cbind.data.frame(time, tb)[show, ]
        colnames(data)[2] <- "ATT"    # rename 2nd column to "ATT"

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
        tb   <- est.att.off
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



    ## print stats
    p.label <- NULL

    if (x$loo == TRUE && loo == TRUE) {
      # recalculate p value and f value
      loo.equiv <- 1
    } else {
      loo.equiv <- 0
    }

    if (type %in% c("equiv", "gap") && loo.equiv == 0) {
      for (i in 1:length(stats)) {
        if ("F.p" %in% stats[i]) {
          if (change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
            x$loo <- FALSE
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 f.threshold = f.threshold
            )
            f.p <- test.out$f.p
          } else {
            f.p <- x$test.out$f.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("F.stat" %in% stats[i]) {
          if (change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
            x$loo <- FALSE
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 f.threshold = f.threshold
            )
            f.stat <- test.out$f.stat
          } else {
            f.stat <- x$test.out$f.stat
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.stat))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("F.equiv.p" %in% stats[i]) {
          # calculate new p value (ziyi: re-add this)
          if (change.f.threshold | change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
            x$loo <- FALSE
            # some problems here, should change to change.f.threshold; change.proportion; change.pre.periods
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 f.threshold = f.threshold
            )
            f.equiv.p <- test.out$f.equiv.p
          } else {
            f.equiv.p <- x$test.out$f.equiv.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("equiv.p" %in% stats[i] && placeboTest == 0) {
          # calculate new p value (ziyi: re-add this)
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group) | use.balance) {
            x$loo <- FALSE
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 tost.threshold = tost.threshold
            )
            tost.equiv.p <- test.out$tost.equiv.p
          } else {
            tost.equiv.p <- x$test.out$tost.equiv.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", tost.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("placebo.p" %in% stats[i]) {
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 tost.threshold = tost.threshold
            )
            placebo.p <- test.out$placebo.p
          } else {
            placebo.p <- x$test.out$placebo.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", placebo.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("equiv.p" %in% stats[i] && placeboTest == 1) {
          p.label1 <- NULL
          # calculate new p value (ziyi: re-add this)
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            placebo.equiv.p <- test.out$placebo.equiv.p
          } else {
            placebo.equiv.p <- x$test.out$placebo.equiv.p
          }
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", placebo.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
      }
    } else if (type %in% c("equiv", "gap") && loo.equiv == 1) { # loo
      for (i in 1:length(stats)) {
        if ("F.p" %in% stats[i]) {
          if (change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 f.threshold = f.threshold
            )
            f.p <- test.out$f.p
          } else {
            f.p <- x$test.out$f.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("F.stat" %in% stats[i]) {
          if (change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 f.threshold = f.threshold
            )
            f.stat <- test.out$f.stat
          } else {
            f.stat <- x$test.out$f.stat
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.stat))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("F.equiv.p" %in% stats[i]) {
          # calculate new p value (ziyi: re-add this)
          if (change.f.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            loo.test.out <- diagtest(x,
                                     proportion = proportion,
                                     pre.periods = pre.periods,
                                     f.threshold = f.threshold
            )
            f.equiv.p <- loo.test.out$f.equiv.p
          } else {
            f.equiv.p <- x$loo.test.out$f.equiv.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", f.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("equiv.p" %in% stats[i]) {
          # calculate new p value (ziyi: re-add this)
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            loo.test.out <- diagtest(x,
                                     proportion = proportion,
                                     pre.periods = pre.periods,
                                     tost.threshold = tost.threshold
            )
            tost.equiv.p <- loo.test.out$tost.equiv.p
          } else {
            tost.equiv.p <- x$loo.test.out$tost.equiv.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", tost.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
      }
    } else if (type == "gap" && placeboTest == TRUE) {
      ## stats
      for (i in 1:length(stats)) {
        if ("placebo.p" %in% stats[i]) {
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 tost.threshold = tost.threshold
            )
            placebo.p <- test.out$placebo.p
          } else {
            placebo.p <- x$test.out$placebo.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", placebo.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("equiv.p" %in% stats[i]) {
          p.label1 <- NULL
          # calculate new p value (ziyi: re-add this)
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            placebo.equiv.p <- test.out$placebo.equiv.p
          } else {
            placebo.equiv.p <- x$test.out$placebo.equiv.p
          }
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", placebo.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
      }
    } else if (type == "exit" && carryoverTest == TRUE) {
      ## stats
      for (i in 1:length(stats)) {
        if ("carryover.p" %in% stats[i]) {
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x,
                                 proportion = proportion,
                                 pre.periods = pre.periods,
                                 tost.threshold = tost.threshold
            )
            carryover.p <- test.out$carryover.p
          } else {
            carryover.p <- x$test.out$carryover.p
          }
          p.label1 <- NULL
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", carryover.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
        if ("equiv.p" %in% stats[i]) {
          p.label1 <- NULL
          # calculate new p value (ziyi: re-add this)
          if (change.tost.threshold | change.proportion | change.pre.periods | !is.null(show.group)) {
            test.out <- diagtest(x, proportion = proportion, pre.periods = pre.periods, tost.threshold = tost.threshold)
            carryover.equiv.p <- test.out$carryover.equiv.p
          } else {
            carryover.equiv.p <- x$test.out$carryover.equiv.p
          }
          p.label1 <- paste0(stats.labs[i], ": ", sprintf("%.3f", carryover.equiv.p))
          p.label <- paste0(p.label, p.label1, "\n")
        }
      }
    }


    if (("none" %in% stats == FALSE) && (show.stats == TRUE)) {
      if (is.null(stats.pos)) {
        if (switch.on == TRUE) {
          stats.pos[1] <- min(data[, "time"], na.rm = 1)
        } else {
          stats.pos[1] <- min(data[, "time"], na.rm = 1)
        }
        ci.top <- max(data[, "CI.upper"], na.rm = 1)
        stats.pos[2] <- ifelse(is.null(ylim), ci.top, ylim[2])
      }

      # Instead of annotating directly, prepare stats to pass to esplot
      # Format the stats values for display
      stats_values <- c()
      # Keep the labels to pass to esplot
      stats_labels <- c()

      # Extract the statistics that were set up earlier
      if (!is.null(p.label)) {
        # Split the existing p.label by newlines to get individual stats
        stats_lines <- strsplit(p.label, "\n")[[1]]
        stats_lines <- stats_lines[stats_lines != ""] # Remove empty lines

        for (line in stats_lines) {
          # Extract the value part from each line (after the colon)
          parts <- strsplit(line, ": ")[[1]]
          if (length(parts) == 2) {
            stats_labels <- c(stats_labels, parts[1])
            stats_values <- c(stats_values, parts[2])
          }
        }
      }
    }
    ## point estimates
    if (classic == 1) {

      #
      # --- REGULAR EVENT-STUDY PLOT (no highlights) ---
      #

      # 1) Build a data frame that esplot() expects
      data_es <- data.frame(
        Period     = data$time,
        ATT        = data$ATT,
        `S.E.`     = if (CI) data$S.E. else NA, # if you have standard errors
        CI.lower   = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else NA,
        CI.upper   = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else NA,
        count      = if (!is.null(data$count)) data$count else NA
      )
      # 2) Call esplot()
      p <- esplot(
        data_es,
        Period       = "Period",
        Estimate     = "ATT",
        SE           = "S.E.",
        CI.lower     = "CI.lower",
        CI.upper     = "CI.upper",
        Count        = "count",
        show.count   = show.count,      # user’s choice
        show.points = show.points,
        ci.outline = ci.outline,
        connected    = connected,  # line+CI or point-range
        color        = color,
        count.color  = count.color,
        xlab         = xlab,
        ylab         = ylab,
        main         = main,
        xlim         = xlim,
        ylim         = ylim,
        gridOff      = gridOff,
        start0       = start0,
        proportion = proportion,
        est.linewidth = est.linewidth,
        est.pointsize  = est.pointsize,
        stats        = if(exists("stats_values")) as.numeric(stats_values) else NULL,
        stats.labs   = if(exists("stats_labels")) stats_labels else NULL,
        stats.pos    = if(show.stats == TRUE && exists("stats_pos")) stats.pos else NULL,
        theme.bw = theme.bw,
        only.pre = type == "equiv")

    } else if (classic == 0 && switch.on == TRUE) {

      #
      # --- PLACEBO TEST PLOT ---
      #

      data_es <- data.frame(
        Period     = data$time,
        ATT        = data$ATT,
        `S.E.`     = if (CI) data$S.E. else NA,
        CI.lower   = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else NA,
        CI.upper   = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else NA,
        count      = if (!is.null(data$count)) data$count else NA
      )

      # highlight the placebo interval [placebo.period[1], placebo.period[2]]
      placebo_seq <- seq(placebo.period[1], placebo.period[2])
      n_placebo   <- length(placebo_seq)

      p <- esplot(
        data_es,
        Period       = "Period",
        Estimate     = "ATT",
        SE           = "S.E.",
        CI.lower     = "CI.lower",
        CI.upper     = "CI.upper",
        Count        = "count",
        show.count   = show.count,
        connected    = connected,
        color        = color,
        count.color  = count.color,
        show.points = show.points,
        ci.outline = ci.outline,
        highlight.periods = placebo_seq,
        highlight.colors  = rep(placebo.color, n_placebo),
        xlab         = xlab,
        ylab         = ylab,
        main         = main,
        xlim         = xlim,
        ylim         = ylim,
        gridOff      = gridOff,
        start0       = start0,
        proportion = proportion,
        est.linewidth = est.linewidth,
        est.pointsize  = est.pointsize,
        stats        = if(exists("stats_values")) as.numeric(stats_values) else NULL,
        stats.labs   = if(exists("stats_labels")) stats_labels else NULL,
        theme.bw = theme.bw,
        stats.pos    = if(show.stats == TRUE && exists("stats_pos")) stats.pos else NULL)

    } else if (classic == 0 && switch.on == FALSE) {
      #
      # --- CARRYOVER TEST OR EXITING TREATMENT ---
      #
      if (is.null(x$est.carry.att)) {
        placebo_seq <- c()
        n_placebo   <- 0
        shift  = 0
      } else{
        placebo_seq <- seq(carryover.period[1], carryover.period[1]-1+dim(x$est.carry.att)[1])
        n_placebo   <- length(placebo_seq)
        shift = dim(x$est.carry.att)[1]
      }

      data_es <- data.frame(
        Period     = data$time + shift,
        ATT        = data$ATT,
        `S.E.`     = if (CI) data$S.E. else NA,
        CI.lower   = if (CI) {
          if (plot.ci == "95") data$CI.lower else data$CI.lower.90
        } else NA,
        CI.upper   = if (CI) {
          if (plot.ci == "95") data$CI.upper else data$CI.upper.90
        } else NA,
        count      = if (!is.null(data$count)) data$count else NA
      )


      carry_seq <- seq(carryover.period[1] + shift, carryover.period[2] + shift)
      n_carry   <- length(carry_seq)

      p <- esplot(
        data_es,
        Period       = "Period",
        Estimate     = "ATT",
        SE           = "S.E.",
        CI.lower     = "CI.lower",
        CI.upper     = "CI.upper",
        Count        = "count",
        show.count   = show.count,
        show.points = show.points,
        ci.outline = ci.outline,
        connected    = connected,
        color        = color,
        count.color  = count.color,
        highlight.periods = c(placebo_seq,carry_seq),
        highlight.colors  = c(rep(carryover.rm.color, n_placebo),rep(carryover.color, n_carry)),
        xlab         = xlab,
        ylab         = ylab,
        main         = main,
        xlim         = xlim,
        ylim         = ylim,
        gridOff      = gridOff,
        start0       = start0,
        proportion = proportion,
        ## newly added options:
        est.linewidth = est.linewidth,
        est.pointsize  = est.pointsize,
        stats        = if(exists("stats_values")) as.numeric(stats_values) else NULL,
        stats.labs   = if(exists("stats_labels")) stats_labels else NULL,
        theme.bw = theme.bw,
        stats.pos    = if(show.stats == TRUE && exists("stats_pos")) stats.pos else NULL)
    }

    # plot bound
    if (bound.old != "none")  { ## with bounds
      p <- p + geom_line(data= data2,aes(time, bound, colour = type, linetype = type, size = type, group = id))
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

      if (effect.bound.ratio == TRUE) {
        if (is.null(stats.pos)) {
          stats.pos[1] <- min(data[, "time"], na.rm = 1)
          stats.pos[2] <- ifelse(is.null(ylim), max(data[, "CI.upper"], na.rm = 1), ylim[1])
        }
        p.label <- paste("ATT / Min. Range = ", sprintf("%.3f", x$att.avg / minBound), sep = "")
        p <- p + annotate("text",
                          x = stats.pos[1], y = stats.pos[2],
                          label = p.label, size = cex.text, hjust = 0
        )
      }
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

    ## y=0 line type
    if (is.null(lcolor) == TRUE) {
      lcolor <- "white"
      if (theme.bw == TRUE) {
        lcolor <- "#aaaaaa"
      }
    }

    if (is.null(lwidth) == TRUE) {
      lwidth <- 1.5
      if (theme.bw == TRUE) {
        lwidth <- 1
      }
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
    p <- p + geom_hline(yintercept = 0, colour = lcolor, size = lwidth)

    TTT <- as.numeric(rownames(data.1))
    TTT.2 <- as.numeric(rownames(data.2))

    if (CI == FALSE) {
      p <- p + geom_hline(yintercept = x$att.avg, color = calendar.line.color, size = 0.8, linetype = "dashed")
      p <- p + geom_point(aes(x = TTT, y = d1[, 1]), color = "gray50", fill = "gray50", alpha = 1, size = 1.2)
      p <- p + geom_line(aes(x = TTT.2, y = d2[, 1]), color = calendar.color, size = 1.1)
    } else {
      p <- p + geom_line(aes(x = TTT.2, y = d2[, 1]), color = calendar.color, size = 1.1)
      p <- p + geom_ribbon(aes(x = TTT.2, ymin = d2[, 3], ymax = d2[, 4]), color = calendar.color, fill = calendar.color, alpha = 0.5, size = 0)
      p <- p + geom_hline(yintercept = x$att.avg, color = calendar.line.color, size = 0.8, linetype = "dashed")
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
      p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = data.toplot, fill = count.color, size = 0.3)
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

    ## y=0 line type
    if (is.null(lcolor) == TRUE) {
      lcolor <- "white"
      if (theme.bw == TRUE) {
        lcolor <- "#aaaaaa"
      }
    }

    if (is.null(lwidth) == TRUE) {
      lwidth <- 1.5
      if (theme.bw == TRUE) {
        lwidth <- 1
      }
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
    p <- p + geom_hline(yintercept = 0, colour = lcolor, size = lwidth)

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

    p <- p + geom_boxplot(aes(x = time, y = eff),
                          position = "dodge", alpha = 0.5,
                          data = data.pre.1, fill =box.control,
                          outlier.fill = box.control, outlier.size = 1.25,
                          outlier.color = box.control,
    )
    p <- p + geom_boxplot(aes(x = time, y = eff),
                          position = "dodge", alpha = 0.5,
                          data = data.post.1, fill = box.treat, outlier.fill =box.treat,
                          outlier.size = 1.25, outlier.color = box.treat,
    )

    p <- p + geom_point(aes(x = time, y = eff),
                        data = data.post.2,
                        color = box.treat, size = 1.25, alpha = 0.8
    )
    p <- p + geom_point(aes(x = time, y = eff),
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
      p <- p + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = data.toplot, fill = count.color, size = 0.3)
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
      p <- p + geom_tile(colour = status.background.color, size = 0.05, stat = "identity")
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
        legend.background = element_rect(fill =status.background.color),
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
  }
  else if (type == "sens") {

    if (restrict == "rm") {
      # Check for the existence of sensitivity results and original data
      if (is.null(x$sensitivity.rm) ||
          is.null(x$sensitivity.rm$results) ||
          is.null(x$sensitivity.rm$original)) {
        stop("No sensitivity results found in x$sensitivity.rm$results or x$sensitivity.rm$original.")
      }

      # Extract the two sets of data
      data_original <- x$sensitivity.rm$original
      data_results  <- x$sensitivity.rm$results

      # For the original data, if the 'Mbar' column is missing, assign a default value to plot it left of zero.
      if (!("Mbar" %in% colnames(data_original))) {
        data_original$Mbar <- -0.05  # Adjust this value if needed.
      }

      # Ensure that the necessary columns are present in each dataset.
      # For original, we require 'lb' and 'ub'; for results, we require 'Mbar', 'lb', and 'ub'.
      if (!all(c("lb", "ub") %in% colnames(data_original))) {
        stop("Original sensitivity results require 'lb' and 'ub'. Check your data.")
      }
      if (!all(c("Mbar", "lb", "ub") %in% colnames(data_results))) {
        stop("Results sensitivity require 'Mbar', 'lb', and 'ub'. Check your data.")
      }

      # Assign color groups: original data will be blue, and results will be red.
      data_original$color_group <- "Original"
      data_results$color_group  <- "Robust Confidence Set"

      # Combine the data into one data frame
      data <- rbind(data_original, data_results)

      p <- ggplot(data, aes(x = Mbar)) +
        geom_hline(yintercept = 0, color = "slategray") +
        geom_errorbar(
          aes(ymin = lb, ymax = ub, color = color_group),
          width = 0.02,  # Width of error bar caps
          linewidth = 1
        ) +
        scale_color_manual(values = c("Original" = sens.original.color, "Robust Confidence Set" = sens.colors[1])) +
        labs(
          x = "M",
          y = "Treatment Effect",
          title = "Smoothness Restriction Sensitivity Analysis"
        ) +
        theme_bw() +
        theme(
          legend.title = element_blank(),
          legend.position = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
        )
    }


  else if (restrict == "sm") {
    # Check for the existence of sensitivity results and original data
    if (is.null(x$sensitivity.smooth) ||
        is.null(x$sensitivity.smooth$results) ||
        is.null(x$sensitivity.smooth$original)) {
      stop("No sensitivity results found in x$sensitivity.smooth$results or x$sensitivity.smooth$original.")
    }

    # Extract the two sets of data
    data_original <- x$sensitivity.smooth$original
    data_results  <- x$sensitivity.smooth$results

    # For the original data, if the 'M' column is missing, assign a default value to plot it left of zero.
    if (!("M" %in% colnames(data_original))) {
      data_original$M <- -0.05  # Adjust this value if needed.
    }

    # Ensure that the necessary columns are present in each dataset.
    # For original, we require 'lb' and 'ub'; for results, we require 'M', 'lb', and 'ub'.
    if (!all(c("lb", "ub") %in% colnames(data_original))) {
      stop("Original sensitivity results require 'lb' and 'ub'. Check your data.")
    }
    if (!all(c("M", "lb", "ub") %in% colnames(data_results))) {
      stop("Results sensitivity require 'M', 'lb', and 'ub'. Check your data.")
    }

    # Assign color groups: original data will be blue, and results will be red.
    data_original$color_group <- "Original"
    data_results$color_group  <- "Robust Confidence Set"

    # Combine the data into one data frame
    data <- rbind(data_original, data_results)

    p <- ggplot(data, aes(x = M)) +
      geom_hline(yintercept = 0, color = "slategray") +
      geom_errorbar(
        aes(ymin = lb, ymax = ub, color = color_group),
        width = 0.02,  # Width of error bar caps
        linewidth = 1
      ) +
      scale_color_manual(values = c("Original" = sens.original.color, "Robust Confidence Set" = sens.colors[1])) +
      labs(
        x = "M",
        y = "Treatment Effect",
        title = "Smoothness Restriction Sensitivity Analysis"
      ) +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
      )
  }
  }

  else if (type == "sens_es") {
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
        ylim <- c(min(c(fect.output.p$CI.lower, dte_output$lb)) * 1.3,
                  max(c(fect.output.p$CI.upper, dte_output$ub)) * 1.3)
      }
      p <- esplot(
        data_es,
        Period       = "Time",
        Estimate     = "ATT",
        SE           = "S.E.",
        CI.lower     = "CI.lower",
        CI.upper     = "CI.upper",
        Count        = "count",
        show.count   = show.count,
        proportion = proportion,
        show.points = show.points,
        ci.outline = ci.outline,
        connected    = connected,
        color        = color,
        count.color  = count.color,
        highlight.periods = x$placebo.period[1]:x$placebo.period[2],
        highlight.colors = rep(placebo.color,x$placebo.period[2]-x$placebo.period[1]+1),
        xlab         = xlab,
        ylab         = ylab,
        main         = main,
        xlim         = xlim,
        ylim         = ylim,
        gridOff      = gridOff,
        start0       = start0,
        est.linewidth = est.linewidth,
        theme.bw = theme.bw,
        est.pointsize  = est.pointsize)

      mbar_levels <- sort(unique(dte_output$Mbar))
      n_colors <- length(mbar_levels)
      final_palette <- sens.colors[1:min(n_colors, length(sens.colors))]
      p <- p +
        geom_linerange(
          aes(x = postPeriod + 0.2, ymin = lb, ymax = ub, color = factor(Mbar, levels = mbar_levels)),
          data = dte_output,
          linewidth = 1
        ) +
        scale_color_manual(name = "Mbar", values = setNames(final_palette, mbar_levels)) +
        guides(color = guide_legend(title = "Mbar")) +
        theme(
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
        )
    }


    else if (restrict == "sm") {

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
        ylim <- c(min(c(fect.output.p$CI.lower, dte_output$lb)) * 1.3,
                  max(c(fect.output.p$CI.upper, dte_output$ub)) * 1.3)
      }
      p <- esplot(
        data_es,
        Period       = "Time",
        Estimate     = "ATT",
        SE           = "S.E.",
        CI.lower     = "CI.lower",
        CI.upper     = "CI.upper",
        Count        = "count",
        show.count   = show.count,
        proportion = proportion,
        show.points = show.points,
        ci.outline = ci.outline,
        connected    = connected,
        color        = color,
        count.color  = count.color,
        highlight.periods = x$placebo.period[1]:x$placebo.period[2],
        highlight.colors = rep(placebo.color,x$placebo.period[2]-x$placebo.period[1]+1),
        xlab         = xlab,
        ylab         = ylab,
        main         = main,
        xlim         = xlim,
        ylim         = ylim,
        gridOff      = gridOff,
        start0       = start0,
        theme.bw = theme.bw,
        est.linewidth = est.linewidth,
        est.pointsize  = est.pointsize)

      m_levels <- sort(unique(dte_output$M))
      n_colors <- length(m_levels)
      final_palette <- sens.colors[1:min(n_colors, length(sens.colors))]
      p <- p +
        geom_linerange(
          aes(x = postPeriod + 0.2, ymin = lb, ymax = ub, color = factor(M, levels = m_levels)),
          data = dte_output,
          linewidth = 1
        ) +
        scale_color_manual(name = "M", values = setNames(final_palette, m_levels)) +
        guides(color = guide_legend(title = "M")) +
        theme(
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
        )

    }
  }    else if (type == "cumul") {
    if (is.null(x$est.eff)) {
      stop("No cumulative ATT data found in x$est.eff.")
    }
    p <- esplot(
      x$est.eff,
      Period    = "Time",
      Estimate  = "ATT",
      SE        = "S.E.",
      CI.lower  = "CI.lower",
      CI.upper  = "CI.upper",
      main      = main,
      xlab      = xlab,
      ylab      = ylab,
      xlim      = xlim,
      ylim      = ylim,
      Count     = "count",
      show.count = show.count,
      proportion = proportion,
      color = color,
      count.color = count.color,
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
