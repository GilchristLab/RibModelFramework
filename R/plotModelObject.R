#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#' @param samples The number of samples in the trace
#' @param mixture The mixture for which to graph values.
#' @param simulated A boolean value that determines whether to use the simulated genome.
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'  
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.

plot.Rcpp_ROCModel <- function(x, genome = NULL, samples = 100, mixture = 1,
                               simulated = FALSE, layout = "original",
                               show.gene.hist = FALSE, show.date = TRUE,
                               color.codon.groups = FALSE, aa.include = NULL,
                               panels = NULL, subtitle = NULL,
                               options = NULL, ...)
{
  model <- x
  opar <- par(no.readonly = T)

  # Apply options list (overrides individual parameters when provided)
  if(!is.null(options)) {
    if(!is.null(options$layout))             layout             <- options$layout
    if(!is.null(options$show.gene.hist))     show.gene.hist     <- options$show.gene.hist
    if(!is.null(options$show.date))          show.date          <- options$show.date
    if(!is.null(options$color.codon.groups)) color.codon.groups <- options$color.codon.groups
    if(!is.null(options$aa.include))         aa.include         <- options$aa.include
    if(!is.null(options$panels))             panels             <- options$panels
    if(!is.null(options$subtitle))           subtitle           <- options$subtitle
  }

  input_list <- as.list(list(...))
  # Backward compat: legacy.layout=FALSE -> layout="compact-v1"
  if("legacy.layout" %in% names(input_list)) {
    warning("'legacy.layout' is deprecated; use layout = \"compact-v1\" instead")
    if(isFALSE(input_list$legacy.layout)) layout <- "compact-v1"
    input_list$legacy.layout <- NULL
  }
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }

  num.genes <- length(genome)
  parameter <- model$getParameter()
  mixtureAssignment <- unlist(lapply(1:num.genes, function(geneIndex){
      parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
      parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)}))
  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  names.aa <- aminoAcids()

  if(layout %in% c("compact", "compact-v1")) {
    # -- compact layout: interleaved AA + marginal columns, n-per-bin strip --
    # Pre-compute bins (shared across all AA panels).
    quantiles  <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE)
    n.bins     <- length(quantiles)
    xlimit.global <- range(expressionValues, na.rm = TRUE)
    # Bin counts computed from the first AA panel result; initialise here.
    bin.counts <- NULL
    bin.mids   <- NULL

    lay <- .buildModelLayout(20L)
    par(oma = c(0, 2, 0, 4), mgp = c(3, 0.56, 0))
    layout(lay$mat, lay$widths, lay$heights, respect = FALSE)

    # title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.92, main, font = 2)
    if(show.date) text(0.5, 0.72, date(), cex = 0.5)

    total.aa.bin.counts <- NULL
    phi.breaks.global   <- NULL
    total.aa.all        <- 0L
    n.aa.drawn <- 0L
    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa, model.type = "ROC")
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability,
                                precomputed.quantiles = quantiles,
                                show.gene.hist = show.gene.hist)
      xlimit <- result$xlimit
      if(is.null(bin.counts)) {
        bin.counts <- result$bin.counts
        bin.mids   <- result$bin.mids
      }
      if(!is.null(result$aa.bin.totals)) {
        if(is.null(total.aa.bin.counts)) {
          total.aa.bin.counts <- result$aa.bin.totals
          phi.breaks.global   <- result$phi.breaks
        } else {
          total.aa.bin.counts <- total.aa.bin.counts + result$aa.bin.totals
        }
        total.aa.all <- total.aa.all + result$total.aa.count
      }
      box()
      aa.x   <- mean(xlimit)
      cex.aa <- 1.2
      usr    <- par("usr")
      pin    <- par("pin")
      tick.y <- 0.02 * min(pin) / pin[2L] * (usr[4L] - usr[3L])
      aa.y   <- usr[4L] - 4 * tick.y - strheight(aa, cex = cex.aa, font = 2) / 2
      rect(aa.x - strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y - strheight(aa, cex = cex.aa, font = 2) * 0.65,
           aa.x + strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y + strheight(aa, cex = cex.aa, font = 2) * 0.65,
           col = adjustcolor("white", alpha.f = 0.75), border = NA)
      text(aa.x, aa.y, aa, cex = cex.aa, font = 2)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1, at = seq(0, 1, by = 0.2))
      if(aa %in% c("V", "Y", "Z"))             axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1, at = seq(0, 1, by = 0.2))
      if(show.gene.hist) {
        axis(1, tck = 0.02, labels = FALSE, pos = 0, lwd = 0, lwd.ticks = 0.5)
        axis(2, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2))
        axis(3, tck = 0.02, labels = FALSE)
        axis(4, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2))
      } else {
        axis(1, tck = 0.02, labels = FALSE)
        axis(2, tck = 0.02, labels = FALSE)
        axis(3, tck = 0.02, labels = FALSE)
        axis(4, tck = 0.02, labels = FALSE)
      }
      axis(2, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2))
      axis(4, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2))
      n.aa.drawn <- n.aa.drawn + 1L
    }

    # phi histogram (20th grid slot): ylim matches AA panels; sub-zero grey strip for alignment
    h.phi <- hist(expressionValues,
                  breaks = seq(xlimit.global[1L], xlimit.global[2L], length.out = 21L),
                  plot = FALSE)
    max.phi        <- max(h.phi$counts)
    hist.space.val <- if(show.gene.hist) 0.15 else 0.05
    ylim.top.val   <- if(show.gene.hist) 1.02 else 1.05
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL,
         xlim = xlimit.global, ylim = c(-(hist.space.val + 0.02), ylim.top.val),
         axes = FALSE, xlab = "", ylab = "")
    if(show.gene.hist) {
      phi.usr <- par("usr")
      rect(phi.usr[1L], phi.usr[3L], phi.usr[2L], 0, col = "grey80", border = NA)
      if(!is.null(total.aa.bin.counts) && max(total.aa.bin.counts) > 0) {
        n.b.ta   <- length(phi.breaks.global)
        scale.ta <- (-0.02 - phi.usr[3L]) / max(total.aa.bin.counts)
        rect(phi.breaks.global[-n.b.ta], phi.usr[3L],
             phi.breaks.global[-1L],     phi.usr[3L] + total.aa.bin.counts * scale.ta,
             col = "white", border = "black", lwd = 0.6)
        text(phi.usr[1L] + 0.02 * (phi.usr[2L] - phi.usr[1L]), -0.01,
             paste0("Count: ", format(total.aa.all, big.mark = ",")),
             adj = c(0, 1), cex = 0.45, col = "black")
      }
      segments(phi.usr[1L], 0, phi.usr[2L], 0, lwd = 0.25, col = "grey50")
    }
    n.b <- length(h.phi$breaks)
    rect(h.phi$breaks[-n.b], 0,
         h.phi$breaks[-1L],  h.phi$counts / max.phi,
         col = "grey80", border = "grey40", lwd = 0.5)
    at.phi <- seq(0, 1, by = 0.2)
    axis(4, las = 1, at = at.phi, labels = round(at.phi * max.phi), cex.axis = 0.8)
    par(xpd = NA)
    gc.usr <- par("usr")
    one.line.x <- par("csi") / par("pin")[1L] * (gc.usr[2L] - gc.usr[1L])
    text(gc.usr[2L] + 2.5 * one.line.x, mean(gc.usr[3:4]),
         "Gene Count", srt = 270, cex = 0.65, adj = 0.5, font = 2)
    par(xpd = FALSE)
    axis(1, cex.axis = 0.8)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    box()

    # x-label (two lines: formula + axis name)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.55, "Gene Expression", cex = 1.4, font = 2)
    text(0.5, 0.28, expression(bold(log[10]~"(Protein Synthesis Rate"~phi~")")), cex = 0.9)

    # y-label (bold, larger, shifted left)
    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.2, 0.5, "Proportion", srt = 90, cex = 1.4, font = 2)

  } else if(!is.null(panels) ||
            layout %in% c("split-ag", "split-ct", "split-ry", "split-6codon")) {
    # -- split wobble layout (preset name, or user-supplied panels=) --
    .runSplitLayout(layout, main, show.date, show.gene.hist, color.codon.groups,
                    aa.include,
                    parameter, model, genome, expressionValues, samples, mixture,
                    prob.fn = function(aa)
                        calculateProbabilityVector(parameter, model, expressionValues,
                            mixture, samples, aa, model.type = "ROC"),
                    panels = panels, subtitle = subtitle)

  } else {
    # -- original layout (default) --
    mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                  nrow = 7, ncol = 4, byrow = TRUE)
    mat <- cbind(rep(23, 7), mat, rep(24, 7))
    nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.75, main, font = 2)
    if(show.date) text(0.5, 0.3, date(), cex = 0.5)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa, model.type = "ROC")
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability)
      xlimit <- result$xlimit
      box()
      aa.x   <- mean(xlimit)
      cex.aa <- 1.2
      usr    <- par("usr")
      pin    <- par("pin")
      tick.y <- 0.02 * min(pin) / pin[2L] * (usr[4L] - usr[3L])
      aa.y   <- usr[4L] - 4 * tick.y - strheight(aa, cex = cex.aa, font = 2) / 2
      rect(aa.x - strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y - strheight(aa, cex = cex.aa, font = 2) * 0.65,
           aa.x + strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y + strheight(aa, cex = cex.aa, font = 2) * 0.65,
           col = adjustcolor("white", alpha.f = 0.75), border = NA)
      text(aa.x, aa.y, aa, cex = cex.aa, font = 2)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("V", "Y", "Z"))             axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1); axis(4, las = 1)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.5, "Propotion", srt = 90)
  }

  par(opar)
}

#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#'
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#'
#' @param samples The number of samples in the trace
#'
#' @param mixture The mixture for which to graph values.
#'
#' @param simulated A boolean value that determines whether to use the simulated genome.
#'
#' @param codon.window A boolean value that determines the codon window to use for calculating codon frequencies. If NULL (the default), use complete sequences.
#'
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.
#'
plot.Rcpp_FONSEModel <- function(x, genome, samples = 100, mixture = 1,
                               simulated = FALSE, codon.window = NULL,
                               layout = "original", show.gene.hist = FALSE,
                               show.date = TRUE, color.codon.groups = FALSE,
                               aa.include = NULL, panels = NULL, subtitle = NULL,
                               options = NULL, ...)
{
  model <- x
  opar <- par(no.readonly = T)

  # Apply options list (overrides individual parameters when provided)
  if(!is.null(options)) {
    if(!is.null(options$layout))             layout             <- options$layout
    if(!is.null(options$show.gene.hist))     show.gene.hist     <- options$show.gene.hist
    if(!is.null(options$show.date))          show.date          <- options$show.date
    if(!is.null(options$codon.window))       codon.window       <- options$codon.window
    if(!is.null(options$color.codon.groups)) color.codon.groups <- options$color.codon.groups
    if(!is.null(options$aa.include))         aa.include         <- options$aa.include
    if(!is.null(options$panels))             panels             <- options$panels
    if(!is.null(options$subtitle))           subtitle           <- options$subtitle
  }

  input_list <- as.list(list(...))
  # Backward compat: legacy.layout=FALSE -> layout="compact-v1"
  if("legacy.layout" %in% names(input_list)) {
    warning("'legacy.layout' is deprecated; use layout = \"compact-v1\" instead")
    if(isFALSE(input_list$legacy.layout)) layout <- "compact-v1"
    input_list$legacy.layout <- NULL
  }
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }

  num.genes <- length(genome)
  parameter <- model$getParameter()
  mixtureAssignment <- unlist(lapply(1:num.genes, function(geneIndex){
      parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
      parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)}))
  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  genes <- genome$getGenes(simulated)
  genome$clear()

  if (is.null(codon.window)) {
    codon.window <- seq(1, 100000)
  } else if (length(codon.window) == 2) {
    codon.window <- seq(codon.window[1], codon.window[2])
  }
  for (i in seq_along(genes)) {
    dna   <- genes[[i]]$seq
    start <- seq(1, nchar(dna), 3)
    stop  <- pmin(start + 2, nchar(dna))
    cods  <- substring(dna, start, stop)
    cods  <- cods[codon.window]
    cods  <- cods[!is.na(cods)]
    genes[[i]]$seq <- paste(cods, collapse = '')
    genome$addGene(genes[[i]], simulated)
  }
  names.aa <- aminoAcids()

  if(layout %in% c("compact", "compact-v1")) {
    # -- compact layout: interleaved AA + marginal columns, n-per-bin strip --
    quantiles     <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE)
    xlimit.global <- range(expressionValues, na.rm = TRUE)
    bin.counts    <- NULL
    bin.mids      <- NULL

    lay <- .buildModelLayout(20L)
    par(oma = c(0, 2, 0, 4), mgp = c(3, 0.56, 0))
    layout(lay$mat, lay$widths, lay$heights, respect = FALSE)

    # title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.92, main, font = 2)
    if(show.date) text(0.5, 0.72, date(), cex = 0.5)

    total.aa.bin.counts <- NULL
    phi.breaks.global   <- NULL
    total.aa.all        <- 0L
    n.aa.drawn <- 0L
    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa,
          model.type = "FONSE", codon.window = codon.window)
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability,
                                precomputed.quantiles = quantiles,
                                show.gene.hist = show.gene.hist)
      xlimit <- result$xlimit
      if(is.null(bin.counts)) {
        bin.counts <- result$bin.counts
        bin.mids   <- result$bin.mids
      }
      if(!is.null(result$aa.bin.totals)) {
        if(is.null(total.aa.bin.counts)) {
          total.aa.bin.counts <- result$aa.bin.totals
          phi.breaks.global   <- result$phi.breaks
        } else {
          total.aa.bin.counts <- total.aa.bin.counts + result$aa.bin.totals
        }
        total.aa.all <- total.aa.all + result$total.aa.count
      }
      box()
      aa.x   <- mean(xlimit)
      cex.aa <- 1.2
      usr    <- par("usr")
      pin    <- par("pin")
      tick.y <- 0.02 * min(pin) / pin[2L] * (usr[4L] - usr[3L])
      aa.y   <- usr[4L] - 4 * tick.y - strheight(aa, cex = cex.aa, font = 2) / 2
      rect(aa.x - strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y - strheight(aa, cex = cex.aa, font = 2) * 0.65,
           aa.x + strwidth(aa,  cex = cex.aa, font = 2) * 0.65,
           aa.y + strheight(aa, cex = cex.aa, font = 2) * 0.65,
           col = adjustcolor("white", alpha.f = 0.75), border = NA)
      text(aa.x, aa.y, aa, cex = cex.aa, font = 2)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1, at = seq(0, 1, by = 0.2))
      if(aa %in% c("V", "Y", "Z"))             axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1, at = seq(0, 1, by = 0.2))
      if(show.gene.hist) {
        axis(1, tck = 0.02, labels = FALSE, pos = 0, lwd = 0, lwd.ticks = 0.5)
        axis(2, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2))
        axis(3, tck = 0.02, labels = FALSE)
        axis(4, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2))
      } else {
        axis(1, tck = 0.02, labels = FALSE)
        axis(2, tck = 0.02, labels = FALSE)
        axis(3, tck = 0.02, labels = FALSE)
        axis(4, tck = 0.02, labels = FALSE)
      }
      axis(2, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2))
      axis(4, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2))
      n.aa.drawn <- n.aa.drawn + 1L
    }

    # phi histogram (20th grid slot): ylim matches AA panels; sub-zero grey strip for alignment
    h.phi <- hist(expressionValues,
                  breaks = seq(xlimit.global[1L], xlimit.global[2L], length.out = 21L),
                  plot = FALSE)
    max.phi        <- max(h.phi$counts)
    hist.space.val <- if(show.gene.hist) 0.15 else 0.05
    ylim.top.val   <- if(show.gene.hist) 1.02 else 1.05
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL,
         xlim = xlimit.global, ylim = c(-(hist.space.val + 0.02), ylim.top.val),
         axes = FALSE, xlab = "", ylab = "")
    if(show.gene.hist) {
      phi.usr <- par("usr")
      rect(phi.usr[1L], phi.usr[3L], phi.usr[2L], 0, col = "grey80", border = NA)
      if(!is.null(total.aa.bin.counts) && max(total.aa.bin.counts) > 0) {
        n.b.ta   <- length(phi.breaks.global)
        scale.ta <- (-0.02 - phi.usr[3L]) / max(total.aa.bin.counts)
        rect(phi.breaks.global[-n.b.ta], phi.usr[3L],
             phi.breaks.global[-1L],     phi.usr[3L] + total.aa.bin.counts * scale.ta,
             col = "white", border = "black", lwd = 0.6)
        text(phi.usr[1L] + 0.02 * (phi.usr[2L] - phi.usr[1L]), -0.01,
             paste0("Count: ", format(total.aa.all, big.mark = ",")),
             adj = c(0, 1), cex = 0.45, col = "black")
      }
      segments(phi.usr[1L], 0, phi.usr[2L], 0, lwd = 0.25, col = "grey50")
    }
    n.b <- length(h.phi$breaks)
    rect(h.phi$breaks[-n.b], 0,
         h.phi$breaks[-1L],  h.phi$counts / max.phi,
         col = "grey80", border = "grey40", lwd = 0.5)
    at.phi <- seq(0, 1, by = 0.2)
    axis(4, las = 1, at = at.phi, labels = round(at.phi * max.phi), cex.axis = 0.8)
    par(xpd = NA)
    gc.usr <- par("usr")
    one.line.x <- par("csi") / par("pin")[1L] * (gc.usr[2L] - gc.usr[1L])
    text(gc.usr[2L] + 2.5 * one.line.x, mean(gc.usr[3:4]),
         "Gene Count", srt = 270, cex = 0.65, adj = 0.5, font = 2)
    par(xpd = FALSE)
    axis(1, cex.axis = 0.8)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    box()

    # x-label (two lines: formula + axis name)
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.55, "Gene Expression", cex = 1.4, font = 2)
    text(0.5, 0.28, expression(bold(log[10]~"(Protein Synthesis Rate"~phi~")")), cex = 0.9)

    # y-label (bold, larger, shifted left)
    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.2, 0.5, "Proportion", srt = 90, cex = 1.4, font = 2)

  } else if(!is.null(panels) ||
            layout %in% c("split-ag", "split-ct", "split-ry", "split-6codon")) {
    # -- split wobble layout (preset name, or user-supplied panels=) --
    .runSplitLayout(layout, main, show.date, show.gene.hist, color.codon.groups,
                    aa.include,
                    parameter, model, genome, expressionValues, samples, mixture,
                    prob.fn = function(aa)
                        calculateProbabilityVector(parameter, model, expressionValues,
                            mixture, samples, aa,
                            model.type = "FONSE", codon.window = codon.window),
                    panels = panels, subtitle = subtitle)

  } else {
    # -- original layout (default) --
    mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                  nrow = 7, ncol = 4, byrow = TRUE)
    mat <- cbind(rep(23, 7), mat, rep(24, 7))
    nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.75, main, font = 2)
    if(show.date) text(0.5, 0.3, date(), cex = 0.5)

    for(aa in names.aa) {
      if(aa == "M" || aa == "W" || aa == "X") next
      codon.probability <- calculateProbabilityVector(
          parameter, model, expressionValues, mixture, samples, aa,
          model.type = "FONSE", codon.window = codon.window)
      result <- plotSinglePanel(parameter, model, genome, expressionValues,
                                samples, mixture, aa, codon.probability)
      xlimit <- result$xlimit
      box()
      main.aa <- aa #TODO map to three letter code
      text(mean(xlimit), 1, main.aa, cex = 1.5)
      if(aa %in% c("A", "F", "K", "Q", "V")) axis(2, las = 1)
      if(aa %in% c("V", "Y", "Z"))             axis(1)
      if(aa %in% c("A", "C", "D", "E"))       axis(3)
      if(aa %in% c("E", "I", "P", "T"))       axis(4, las = 1)
      axis(1, tck = 0.02, labels = FALSE)
      axis(2, tck = 0.02, labels = FALSE)
      axis(3, tck = 0.02, labels = FALSE)
      axis(4, tck = 0.02, labels = FALSE)
    }

    hist.values <- hist(expressionValues, plot = FALSE, nclass = 30)
    plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
    axis(1); axis(4, las = 1)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.5, "Propotion", srt = 90)
  }

  par(opar)
}

calculateProbabilityVector <- function(parameter,model,expressionValues,mixture,samples,aa,model.type="ROC",codon.window = c(1,300))
{
  codons <- AAToCodon(aa, T)
  
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T,log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T,log_scale = F)
  }
  
  # calculate codon probabilities with respect to phi
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  if (model.type == "ROC")
  {
    codonProbability <- lapply(10^phis,  
                               function(phi){
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi)
                               })
    codonProbability <- do.call("rbind",codonProbability)
  } else if (model.type == "FONSE")
  {
    codonProbability <- vector(mode = "list", length = length(codon.window))
    for (i in 1:length(codon.window))
    {
      codonProbability.tmp <- lapply(10^phis,
                               function(phi)
                                 {
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi,codon.window[i])
                               })
      codonProbability[[i]] <- do.call("rbind",codonProbability.tmp)
    }
    codonProbability <- Reduce("+",codonProbability)/length(codon.window)
  }
  return(codonProbability)
}



## Current plot layout version. Bump when a layout changes incompatibly.
.PLOT_LAYOUT_VERSION <- 1L

#' Query the current plot layout version
#' @return Integer version number.
#' @export
plotLayoutVersion <- function() .PLOT_LAYOUT_VERSION

## Internal: common plot options shared by ROC and FONSE plot methods.
.plotModelOptions.common <- function(
    layout             = "original",
    show.gene.hist     = FALSE,
    show.date          = TRUE,
    color.codon.groups = FALSE,
    aa.include         = NULL,
    panels             = NULL,
    subtitle           = NULL
) {
    list(layout = layout, show.gene.hist = show.gene.hist, show.date = show.date,
         color.codon.groups = color.codon.groups, aa.include = aa.include,
         panels = panels, subtitle = subtitle)
}

#' Plot options for ROC model
#'
#' Returns a named list of plotting options for \code{plot.Rcpp_ROCModel}.
#' Pass the result to the \code{options} argument of \code{plot()}.
#'
#' @param layout Character. Layout to use. One of:
#'   \code{"original"} (default, classic layout),
#'   \code{"compact"} / \code{"compact-v1"} (compact 4x5 grid),
#'   \code{"split-ag"} (A/G wobble-pair panels),
#'   \code{"split-ct"} (T/C wobble-pair panels),
#'   \code{"split-ry"} (R vs Y aggregate comparison panels).
#' @param show.gene.hist Logical. When \code{TRUE} and layout supports it,
#'   draws per-AA gene-count histogram strips below each panel.
#' @param show.date Logical. Whether to include a timestamp in the title panel.
#' @return A named list of plot options to pass to \code{options} in \code{plot()}.
#' @examples
#' \dontrun{
#'   opts <- plotROCOptions(layout = "compact-v1", show.gene.hist = TRUE)
#'   plot(model, genome, samples = 500, main = "My fit", options = opts)
#' }
#' @export
plotROCOptions <- function(...) {
    .plotModelOptions.common(...)
}

#' Plot options for FONSE model
#'
#' Returns a named list of plotting options for \code{plot.Rcpp_FONSEModel}.
#' Pass the result to the \code{options} argument of \code{plot()}.
#'
#' @inheritParams plotROCOptions
#' @param codon.window Integer vector of codon positions to include.
#'   \code{NULL} (default) uses all positions.
#' @return A named list of plot options to pass to \code{options} in \code{plot()}.
#' @examples
#' \dontrun{
#'   opts <- plotFONSEOptions(layout = "compact-v1", codon.window = c(1, 50))
#'   plot(model, genome, samples = 500, options = opts)
#' }
#' @export
plotFONSEOptions <- function(..., codon.window = NULL) {
    opts <- .plotModelOptions.common(...)
    c(opts, list(codon.window = codon.window))
}

## Internal: build layout matrix for compact codon-freq vs phi plot.
## n.slots: number of data panel slots (20 for ROC). 19 slots show AA codon-freq
## panels; the 20th slot (bottom-right of the grid) shows the phi histogram.
## Skipped AAs (M/W/X) are not plotted, so the phi histogram naturally lands in
## slot 20 after 19 AA plot.new() calls -- no blank fill needed.
## Panel numbering (in drawing order):
##   1          = title
##   2 .. n+1   = n data slots (19 AA panels + phi histogram as slot n)
##   n+2        = x-label
##   n+3        = y-label (left column)
# n.rows: optional minimum row count; pads with 0-cells so panel height stays
# consistent across layouts with different panel counts.  0-cells are skipped
# automatically by R's layout() device cursor -- no explicit fill calls needed.
.buildModelLayout <- function(n.slots = 20L, n.rows = NULL) {
    # Natural rows = rows actually filled by panels.  Total rows may be floored
    # higher (n.rows) so panel HEIGHT stays consistent across layouts with
    # different panel counts.  The extra (pad) rows are placed at the BOTTOM,
    # below the x-label row, so the x-label sits directly under the last used
    # row and the y-label spans only the used rows.
    n.rows.natural <- ceiling(n.slots / 4L)
    n.rows.total   <- max(n.rows.natural,
                          if(!is.null(n.rows)) as.integer(n.rows) else 0L)
    n.rows.pad     <- n.rows.total - n.rows.natural
    base  <- 1L
    x.lbl <- base + n.slots + 1L
    y.lbl <- base + n.slots + 2L

    # 4-column data block over the used (natural) rows only.
    data.mat <- matrix(0L, nrow = n.rows.natural, ncol = 4L)
    pnum <- 2L
    for(r in seq_len(n.rows.natural)) {
        for(cp in seq_len(4L)) {
            if(pnum <= base + n.slots) {
                data.mat[r, cp] <- pnum
                pnum <- pnum + 1L
            }
        }
    }

    # Assemble: title, used data rows (y-label spans these in col 1), x-label
    # row immediately below the last used row, then blank padding rows at the
    # bottom.  Col 1 of the title and x-label rows is left blank (0) so the
    # y-label is centred on the used data rows only.
    title.row <- c(0L,     rep(1L,    4L))
    data.rows <- cbind(y.lbl, data.mat)
    xlbl.row  <- c(0L,     rep(x.lbl, 4L))
    mat       <- rbind(title.row, data.rows, xlbl.row)
    heights   <- c(3, rep(8, n.rows.natural), 4)
    if(n.rows.pad > 0L) {
        for(i in seq_len(n.rows.pad)) mat <- rbind(mat, rep(0L, 5L))
        heights <- c(heights, rep(8, n.rows.pad))
    }

    list(
        mat            = mat,
        widths         = c(3, 8, 8, 8, 8),
        heights        = heights,
        x.lbl          = x.lbl,
        y.lbl          = y.lbl,
        n.rows.data    = n.rows.total,
        n.rows.natural = n.rows.natural
    )
}

# NOT EXPOSED
# precomputed.quantiles: if provided (compact layout), use instead of recomputing.
# show.gene.hist: if TRUE, draw a phi-binned gene-count histogram at the bottom
#   of the panel (bars scaled to ~20% of panel height, same x-axis).
# Returns a list(xlimit, codonCounts, codons, bin.mids, bin.counts) so the caller
# can draw the n-per-bin strip without re-reading the genome.
plotSinglePanel <- function(parameter, model, genome, expressionValues, samples,
                            mixture, aa, codon.probability,
                            precomputed.quantiles = NULL,
                            show.gene.hist = FALSE)
{
  codons <- AAToCodon(aa, T)
  #
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T, log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T, log_scale = F)
  }
  #
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  codonProbability <- codon.probability
  #get codon counts
  codons <- AAToCodon(aa, F)
  codonCounts <- vector("list", length(codons))
  for(i in 1:length(codons))
  {
    codonCounts[[i]] <- genome$getCodonCountsPerGene(codons[i])
  }
  codonCounts <- do.call("cbind", codonCounts)
  aa.count.per.gene <- rowSums(codonCounts)  # raw AA occurrences per gene (before normalizing)
  # codon proportions
  codonCounts <- codonCounts / rowSums(codonCounts)
  codonCounts[is.nan(codonCounts)] <- NA # necessary if AA does not appear in gene

  # make empty plot
  xlimit    <- range(expressionValues, na.rm = T)
  hist.space <- if(show.gene.hist) 0.15 else 0.05
  ylim.top   <- if(show.gene.hist) 1.02 else 1.05
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim=xlimit, ylim=c(-(hist.space + 0.02), ylim.top),
       xlab = "", ylab="", axes = FALSE)
  # bin expression values of genes
  quantiles <- if(!is.null(precomputed.quantiles)) precomputed.quantiles else
                   quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)
  n.bins    <- length(quantiles)
  bin.counts <- integer(n.bins)
  bin.mids   <- numeric(n.bins)

  # AA-count histogram: equal-width phi bins; drawn first so colored data appears on top
  aa.bin.totals.ret <- NULL
  phi.breaks.ret    <- NULL
  total.aa.count    <- 0L
  if(show.gene.hist) {
    # par("usr") gives actual plot extents after axis expansion -- use for all geometry
    usr <- par("usr")
    rect(usr[1L], usr[3L], usr[2L], 0, col = "grey80", border = NA)
    aa.present  <- !apply(is.na(codonCounts), 1L, all)
    phi.aa      <- expressionValues[aa.present]
    counts.aa   <- aa.count.per.gene[aa.present]
    if(length(phi.aa) > 1L) {
      n.hist.bins   <- 20L
      breaks        <- seq(xlimit[1L], xlimit[2L], length.out = n.hist.bins + 1L)
      bin.idx       <- findInterval(phi.aa, breaks, rightmost.closed = TRUE)
      bin.idx       <- pmax(1L, pmin(bin.idx, n.hist.bins))
      aa.bin.totals <- vapply(seq_len(n.hist.bins), function(b) {
                         sum(counts.aa[bin.idx == b])
                       }, numeric(1L))
      aa.bin.totals.ret <- aa.bin.totals
      phi.breaks.ret    <- breaks
      total.aa.count    <- as.integer(sum(counts.aa))
      bar.max <- max(aa.bin.totals)
      if(bar.max > 0L) {
        scale <- (-0.02 - usr[3L]) / bar.max  # tallest bar stops 0.02 below separator
        rect(breaks[-length(breaks)], usr[3L],
             breaks[-1L],             usr[3L] + aa.bin.totals * scale,
             col = "white", border = "black", lwd = 0.6)
      }
      # total observed AA count label: left of grey strip
      total.aa <- sum(counts.aa)
      text(usr[1L] + 0.02 * (usr[2L] - usr[1L]), -0.01,
           paste0("Count: ", format(total.aa, big.mark = ",")),
           adj = c(0, 1), cex = 0.45, col = "black")
    }
    # separator line at y=0: lighter than panel frame so strip reads as part of panel above
    segments(usr[1L], 0, usr[2L], 0, lwd = 0.25, col = "grey50")
  }

  # reference line at uniform codon usage (1/n synonymous codons)
  abline(h = 1 / ncol(codonCounts), col = "grey70", lty = 2, lwd = 0.8)

  for(i in 1:n.bins)
  {
    if(i == 1){
      tmp.id <- expressionValues < quantiles[i]
    }else if(i == n.bins){
      tmp.id <- expressionValues > quantiles[i]
    }else{
      tmp.id <- expressionValues > quantiles[i] & expressionValues < quantiles[i + 1]
    }
    bin.counts[i] <- sum(tmp.id, na.rm = TRUE)
    bin.mids[i]   <- median(expressionValues[tmp.id], na.rm = TRUE)

    # plot quantiles
    means <- colMeans(codonCounts[tmp.id,], na.rm = T)
    std <- apply(codonCounts[tmp.id,], 2, sd, na.rm = T)
    for(k in 1:length(codons))
    {
      points(bin.mids[i], means[k],
             col=.codonColors[[ codons[k] ]] , pch=19, cex = 0.5)
      lines(rep(bin.mids[i], 2), pmax(0, pmin(1, c(means[k]-std[k], means[k]+std[k]))),
            col=.codonColors[[ codons[k] ]], lwd=0.8)
    }
  }

  # draw model fit
  for(i in 1:length(codons))
  {
    lines(phis, codonProbability[, i], col=.codonColors[[ codons[i] ]])
  }
  colors <- unlist(.codonColors[codons])

  # add indicator to optimal codon
  optim.codon.index <- which(min(c(selection, 0)) == c(selection, 0))
  codons[optim.codon.index] <- paste0(codons[optim.codon.index], "*")
  legend("topleft", legend = codons, col=colors, lty=1, cex=0.75,
         bty = "o", bg = adjustcolor("white", alpha.f = 0.7), box.lwd = 0.5)

  invisible(list(xlimit = xlimit, codonCounts = codonCounts,
                 codons = codons, bin.mids = bin.mids, bin.counts = bin.counts,
                 aa.bin.totals = aa.bin.totals.ret, phi.breaks = phi.breaks.ret,
                 total.aa.count = total.aa.count))
}

## ============================================================
## Wobble-split layout helpers
## ============================================================

# Classify codons by 3rd-position (wobble) nucleotide.
# Returns: r.idx = A/G-ending indices, y.idx = C/T-ending indices.
# Build a codon-family label without a type prefix, e.g. c("CTA","CTC","CTG","CTT") -> "CTN",
# c("TTA","TTG") -> "TTA/G".
.blockLabel <- function(codon.names) {
    if(length(codon.names) == 0L) return("--")
    prefix  <- substr(codon.names[1L], 1L, 2L)
    wobbles <- sort(substr(codon.names, 3L, 3L))
    if(identical(wobbles, c("A","C","G","T"))) paste0(prefix, "N")
    else paste0(prefix, paste(wobbles, collapse = "/"))
}

# Build a legend label from actual codon names, e.g. c("GCA","GCG") -> "R: GCA/G".
# Codons sharing a 2-nt prefix are collapsed to "PFX<w1>/<w2>/...";
# multi-prefix groups (e.g. ry-all6) join the per-prefix tokens with "+".
.codonGroupLabel <- function(type.prefix, codon.names) {
    if(length(codon.names) == 0L) return(paste0(type.prefix, ": --"))
    if(length(codon.names) == 1L) return(paste0(type.prefix, ": ", codon.names))
    pfx <- substr(codon.names, 1L, 2L)
    make.token <- function(cods) {
        p <- substr(cods[1L], 1L, 2L)
        w <- substr(cods, 3L, 3L)
        paste0(p, paste(w, collapse = "/"))
    }
    token <- if(length(unique(pfx)) == 1L) {
        make.token(codon.names)
    } else {
        paste(vapply(split(codon.names, pfx), make.token, character(1L)),
              collapse = "+")
    }
    paste0(type.prefix, ": ", token)
}

# Build a position-pattern label from a set of codons.
# Each position: fixed nt shown as-is; variable shown as "(A/G)" etc. (sorted).
# E.g. c("CTA","CTG","TTA","TTG") -> "(C/T)T(A/G)"
.codonPattern <- function(codon.names) {
    if(length(codon.names) == 0L) return("--")
    m <- do.call(rbind, strsplit(codon.names, ""))
    paste(apply(m, 2L, function(nts) {
        u <- sort(unique(nts))
        if(length(u) == 1L) u else paste0("(", paste(u, collapse = "/"), ")")
    }), collapse = "")
}

.classifyWobble <- function(codons) {
    wobble <- substr(codons, 3L, 3L)
    list(r.idx  = which(wobble %in% c("A", "G")),
         y.idx  = which(wobble %in% c("C", "T")),
         wobble = wobble)
}

# For 6-codon AAs (L=Leu, R=Arg), return 4-codon and 2-codon block indices.
.get6CodonBlocks <- function(aa, codons) {
    prefix4 <- switch(aa, L = "CT", R = "CG", NULL)
    if(is.null(prefix4)) return(NULL)
    prefix2 <- switch(aa, L = "TT", R = "AG")
    list(block4.idx = which(startsWith(codons, prefix4)),
         block2.idx = which(startsWith(codons, prefix2)))
}

# Return ordered list of panel descriptors for a split layout page.
# Each descriptor: list(aa, block, type, label).
#   aa:    single-letter AA code
#   block: "4", "2", "all", or NULL
#   type:  "ag" | "tc" | "ry" | "4v2" | "ry-all6"
#   label: display label placed inside the panel
# ---------------------------------------------------------------------------
# Grouping-driven panel construction
#
# A split plot is described as a list of GROUPS, each = a grouping description
# applied to a list of amino acids (aaGroup(grouping, aa)).  The grouping names
# are user-facing aliases for the internal panel `type`/`block` fields consumed
# by plotWobbleSplitPanel; .expandGrouping() is the single place that encodes
# how each grouping expands one AA into its panel(s).  Pre-designed layouts
# (.SPLIT.PRESETS) and a user-supplied `panels=` argument both flow through the
# same builder (.buildSplitPanels) -- no per-layout switch.
# ---------------------------------------------------------------------------

# Recognised grouping descriptions (user-facing names).
.SPLIT.GROUPINGS <- c("wobble-purine", "wobble-pyrimidine", "ry-aggregate",
                      "block-4v2", "r-block", "6codon-detail")

# Codon-count class of an AA: "2","3","4","6" (19-AA convention; Ser split S/Z).
.aaCodonGroup <- function(aa) {
    if(aa %in% c("C","D","E","F","H","K","N","Q","Y","Z")) "2"
    else if(aa == "I")                                     "3"
    else if(aa %in% c("A","G","P","S","T","V"))            "4"
    else                                                    "6"
}

# One panel descriptor (fields consumed downstream by plotWobbleSplitPanel).
.panelDescriptor <- function(aa, block = NULL, type, label = aa, sublabel = NULL) {
    list(aa = aa, block = block, type = type, label = label,
         sublabel = sublabel, codon.group = .aaCodonGroup(aa))
}

# Expand one (grouping, aa) into its panel descriptor(s).
.expandGrouping <- function(grouping, aa) {
    cc <- .aaCodonGroup(aa)
    P  <- .panelDescriptor
    switch(grouping,
        "wobble-purine" =
            if(cc == "6")
                list(P(aa, "4", "ag", sublabel = "4 Codon Block"),
                     P(aa, "2", "ag", sublabel = "2 Codon Block"))
            else list(P(aa, type = "ag")),
        "wobble-pyrimidine" =
            if(cc == "6")
                list(P(aa, "4", "tc", sublabel = "4 Codon Block"))
            else list(P(aa, type = "tc")),
        "ry-aggregate" =
            if(cc == "6")
                list(P(aa, "4", "ry"), P(aa, "all", "ry-all6"))
            else list(P(aa, type = "ry")),
        "block-4v2" = list(P(aa, "all", "4v2")),
        "r-block"   = list(P(aa, NULL,  "r-block")),
        "6codon-detail" = list(
            P(aa, "all", "4v2",     sublabel = "Block Proportion"),
            P(aa, NULL,  "r-block", sublabel = "R-Wobble Block"),
            P(aa, "4",   "tc",      sublabel = "Y-Wobble")),
        stop("unknown grouping: ", grouping))
}

#' Describe a group of amino-acid panels for a split codon plot
#'
#' Building block for the \code{panels} argument of \code{plotROCOptions}: a
#' grouping description applied to a list of amino acids.  The same constructor
#' defines the built-in layouts and any user-defined plot.
#'
#' @param grouping Character. One of \code{"wobble-purine"},
#'   \code{"wobble-pyrimidine"}, \code{"ry-aggregate"}, \code{"block-4v2"},
#'   \code{"r-block"}, \code{"6codon-detail"}.
#' @param aa Character vector of amino-acid letters (19-AA convention; Ser
#'   split into \code{S}=TCN and \code{Z}=AGY).
#' @return A group spec consumed by the split-layout builder.
#' @examples
#' \dontrun{
#'   plotROCOptions(panels = list(
#'       aaGroup("wobble-purine", c("E","K","Q")),
#'       aaGroup("block-4v2",     c("L","R"))))
#' }
#' @export
aaGroup <- function(grouping, aa) {
    grouping <- match.arg(grouping, .SPLIT.GROUPINGS)
    list(grouping = grouping, aa = aa)
}

# Flatten a list of aaGroup() specs into a flat list of panel descriptors.
.buildSplitPanels <- function(group.list) {
    do.call(c, lapply(group.list, function(g)
        do.call(c, lapply(g$aa, function(a) .expandGrouping(g$grouping, a)))))
}

# Pre-designed layouts expressed as group lists (the former hard-coded switch).
.SPLIT.PRESETS <- list(
    # 2-codon (pure R) -> 4-codon -> 6-codon (each 4-block + 2-block).
    "split-ag" = list(
        aaGroup("wobble-purine", c("E","K","Q")),
        aaGroup("wobble-purine", c("A","G","P","S","T","V")),
        aaGroup("wobble-purine", c("L","R"))),
    # 2-codon (pure Y) -> 3-codon (Ile) -> 4-codon -> 6-codon (4-block only).
    "split-ct" = list(
        aaGroup("wobble-pyrimidine", c("C","D","N","H","F","Y","Z")),
        aaGroup("wobble-pyrimidine", "I"),
        aaGroup("wobble-pyrimidine", c("A","G","P","S","T","V")),
        aaGroup("wobble-pyrimidine", c("L","R"))),
    # 3-codon (Ile) -> 4-codon -> 6-codon (4-block + all-6 R/Y).
    "split-ry" = list(
        aaGroup("ry-aggregate", "I"),
        aaGroup("ry-aggregate", c("A","G","P","S","T","V")),
        aaGroup("ry-aggregate", c("L","R"))),
    # 6-codon AA (Leu, Arg) block-structure page: 4v2 / r-block / Y-wobble per AA.
    "split-6codon" = list(
        aaGroup("6codon-detail", c("L","R")))
)

# Resolve a preset layout name to its flat panel-descriptor list.
.getSplitPanels <- function(layout.type) {
    preset <- .SPLIT.PRESETS[[layout.type]]
    if(is.null(preset)) return(list())
    .buildSplitPanels(preset)
}

# Pad each codon-count group to a row boundary (multiple of n.cols) by
# appending NULL entries.  NULL slots are drawn as empty panels (plot.new()).
# Groups are ordered 2->3->4->6 so they appear on consecutive row bands.
.padGroupsToRows <- function(panels, n.cols = 4L) {
    if(length(panels) == 0L) return(list())
    grp.order <- c("2", "3", "4", "6")
    present   <- sapply(panels, `[[`, "codon.group")
    groups    <- grp.order[grp.order %in% present]
    result    <- list()
    for(g in groups) {
        grp.panels <- Filter(function(pd) pd$codon.group == g, panels)
        n   <- length(grp.panels)
        pad <- (n.cols - n %% n.cols) %% n.cols
        result <- c(result, grp.panels, vector("list", pad))
    }
    result
}

# Draw a single panel for split-ag / split-ct / split-ry layouts.
# panel.desc: one list from .getSplitPanels().
# codon.probability: full probability matrix from calculateProbabilityVector().
plotWobbleSplitPanel <- function(parameter, model, genome, expressionValues,
                                  samples, mixture, panel.desc, codon.probability,
                                  precomputed.quantiles = NULL,
                                  show.gene.hist = FALSE)
{
    aa          <- panel.desc$aa
    codon.block <- panel.desc$block
    split.type  <- panel.desc$type

    codons.par <- AAToCodon(aa, TRUE)
    codons     <- AAToCodon(aa, FALSE)

    expression.range <- range(expressionValues)
    phis <- seq(from = expression.range[1L], to = expression.range[2L], by = 0.01)

    # Codon counts per gene
    codonCounts.raw <- do.call("cbind", lapply(codons, function(c)
        genome$getCodonCountsPerGene(c)))
    aa.count.per.gene <- rowSums(codonCounts.raw)
    codonCounts.prop  <- codonCounts.raw / aa.count.per.gene
    codonCounts.prop[is.nan(codonCounts.prop)] <- NA

    wob    <- .classifyWobble(codons)
    blocks <- .get6CodonBlocks(aa, codons)

    # ---- Determine display mode ----
    if(split.type %in% c("ag", "tc")) {
        # Individual codons filtered by wobble class
        base.idx <- if(split.type == "ag")
            (if(aa == "I") wob$y.idx else wob$r.idx)
        else wob$y.idx
        if(!is.null(codon.block) && codon.block != "all" && !is.null(blocks)) {
            bi       <- if(codon.block == "4") blocks$block4.idx else blocks$block2.idx
            base.idx <- intersect(base.idx, bi)
        }
        show.codons    <- codons[base.idx]
        show.probs.raw <- codon.probability[, base.idx, drop = FALSE]

        # For >2-codon AAs, rescale to pair-relative (0-1 within just these codons).
        # pair.overall holds the pair's combined share of all AA codons (for background).
        should.rescale <- length(base.idx) < length(codons)
        if(should.rescale) {
            raw.pair     <- codonCounts.raw[, base.idx, drop = FALSE]
            pair.total   <- rowSums(raw.pair)
            pair.total[pair.total == 0L] <- NA_real_
            show.counts  <- raw.pair / pair.total
            show.counts[is.nan(show.counts)] <- NA_real_
            prob.sum     <- pmax(rowSums(show.probs.raw), 1e-300)
            show.probs   <- show.probs.raw / prob.sum
            # Histogram counts: just the codons in this pair
            hist.counts.pair <- pair.total
            hist.counts.pair[is.na(hist.counts.pair)] <- 0L
        } else {
            show.counts      <- codonCounts.prop[, base.idx, drop = FALSE]
            show.probs       <- show.probs.raw
            hist.counts.pair <- aa.count.per.gene
        }
        show.colors <- unlist(.codonColors[show.codons])

        # Optimal codon (mark with *)
        sel.all <- vapply(seq_along(codons.par), function(i)
            parameter$getCodonSpecificPosteriorMean(
                mixture, samples, codons.par[i], 1L, TRUE, log_scale = FALSE),
            numeric(1L))
        sel.sub  <- sel.all[base.idx]
        opt.i    <- which(min(c(sel.sub, 0)) == c(sel.sub, 0))
        show.lbl <- show.codons
        if(length(opt.i) && opt.i[1L] <= length(show.lbl))
            show.lbl[opt.i[1L]] <- paste0(show.lbl[opt.i[1L]], "*")

    } else {
        should.rescale <- FALSE
        hist.counts.pair <- aa.count.per.gene
        # Aggregate groups
        if(split.type == "ry") {
            if(aa == "I") {
                # Ile: ATA (A-ending, purine) vs ATC+ATT (C/T-ending, pyrimidine)
                g1.idx <- wob$r.idx; g2.idx <- wob$y.idx
            } else {
                if(!is.null(codon.block) && codon.block == "4" && !is.null(blocks)) {
                    g1.idx <- intersect(wob$r.idx, blocks$block4.idx)
                    g2.idx <- intersect(wob$y.idx, blocks$block4.idx)
                } else {
                    g1.idx <- wob$r.idx; g2.idx <- wob$y.idx
                }
            }
            g1.lbl <- .codonPattern(codons[g1.idx])
            g2.lbl <- .codonPattern(codons[g2.idx])
        } else if(split.type == "4v2") {
            g1.idx <- blocks$block4.idx; g2.idx <- blocks$block2.idx
            g1.lbl <- .blockLabel(codons[blocks$block4.idx])
            g2.lbl <- .blockLabel(codons[blocks$block2.idx])
        } else if(split.type == "r-block") {
            # R-wobble codons split by block: 4-block R vs 2-block R.
            g1.idx <- intersect(wob$r.idx, blocks$block4.idx)
            g2.idx <- intersect(wob$r.idx, blocks$block2.idx)
            g1.lbl <- .codonPattern(codons[g1.idx])
            g2.lbl <- .codonPattern(codons[g2.idx])
            # Histogram counts: only R-wobble codons (not full AA).
            hist.counts.pair <- rowSums(codonCounts.raw[, c(g1.idx, g2.idx), drop = FALSE])
        } else {  # ry-all6
            g1.idx <- wob$r.idx; g2.idx <- wob$y.idx
            g1.lbl <- .codonPattern(codons[g1.idx])
            g2.lbl <- .codonPattern(codons[g2.idx])
        }
        # 4v2 and r-block (block comparisons) share vermillion/green.
        # ry / ry-all6 keep orange/blue for purine/pyrimidine.
        if(split.type %in% c("4v2", "r-block")) {
            g1.col <- "#D55E00"   # vermillion (4-codon block)
            g2.col <- "#009E73"   # bluish-green (2-codon block)
        } else {
            g1.col <- "#E69F00"   # orange (purine / R wobble)
            g2.col <- "#56B4E9"   # blue   (pyrimidine / Y wobble)
        }
        # Denominator = sum over only the shown codons (g1 + g2).
        # For ry block-4, this is the CG*/CT* 4-codon total, not all-AA total.
        grp.denom <- rowSums(codonCounts.raw[, c(g1.idx, g2.idx), drop = FALSE])
        grp.denom[grp.denom == 0L] <- NA_real_
        g1.prop  <- rowSums(codonCounts.raw[, g1.idx, drop = FALSE]) / grp.denom
        g2.prop  <- rowSums(codonCounts.raw[, g2.idx, drop = FALSE]) / grp.denom
        g1.prop[is.nan(g1.prop)] <- NA
        g2.prop[is.nan(g2.prop)] <- NA
        grp.model.tot <- rowSums(codon.probability[, c(g1.idx, g2.idx), drop = FALSE])
        grp.model.tot[grp.model.tot == 0] <- NA_real_
        g1.model <- rowSums(codon.probability[, g1.idx, drop = FALSE]) / grp.model.tot
        g2.model <- rowSums(codon.probability[, g2.idx, drop = FALSE]) / grp.model.tot
    }

    # ---- Set up plot ----
    xlimit     <- range(expressionValues, na.rm = TRUE)
    hist.space <- if(show.gene.hist) 0.15 else 0.05
    ylim.top   <- if(show.gene.hist) 1.02 else 1.05
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = xlimit, ylim = c(-(hist.space + 0.02), ylim.top),
         xlab = "", ylab = "", axes = FALSE)

    # ---- Histogram strip ----
    aa.bin.totals.ret <- NULL
    phi.breaks.ret    <- NULL
    total.aa.count    <- 0L
    if(show.gene.hist) {
        usr        <- par("usr")
        aa.present <- hist.counts.pair > 0
        phi.aa     <- expressionValues[aa.present]
        counts.aa  <- hist.counts.pair[aa.present]
        rect(usr[1L], usr[3L], usr[2L], 0, col = "grey80", border = NA)
        if(length(phi.aa) > 1L) {
            n.hb   <- 20L
            breaks <- seq(xlimit[1L], xlimit[2L], length.out = n.hb + 1L)
            bidx   <- pmax(1L, pmin(findInterval(phi.aa, breaks,
                                                  rightmost.closed = TRUE), n.hb))
            abt    <- vapply(seq_len(n.hb), function(b) sum(counts.aa[bidx == b]),
                             numeric(1L))
            aa.bin.totals.ret <- abt
            phi.breaks.ret    <- breaks
            total.aa.count    <- as.integer(sum(counts.aa))
            bmax <- max(abt)
            if(bmax > 0L) {
                sc <- (-0.02 - usr[3L]) / bmax
                rect(breaks[-length(breaks)], usr[3L],
                     breaks[-1L], usr[3L] + abt * sc,
                     col = "white", border = "black", lwd = 0.6)
            }
            text(usr[1L] + 0.02 * (usr[2L] - usr[1L]), -0.01,
                 paste0("Count: ", format(total.aa.count, big.mark = ",")),
                 adj = c(0, 1), cex = 0.45, col = "black")
        }
        segments(usr[1L], 0, usr[2L], 0, lwd = 0.25, col = "grey50")
    }

    # ---- Reference line (ag/tc only) ----
    quantiles <- if(!is.null(precomputed.quantiles)) precomputed.quantiles else
                 quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = TRUE)
    n.bins     <- length(quantiles)
    if(split.type %in% c("ag", "tc")) {
        # Rescaled panels: reference at 0.5 (equal usage within pair).
        # 2-codon AAs (no rescaling): reference at 1/n_synonymous.
        ref.h <- if(should.rescale) 0.5 else 1 / length(codons)
        abline(h = ref.h, col = "grey70", lty = 2, lwd = 0.8)
    } else {
        # ry / 4v2 / ry-all6: reference lines at null expectation (uniform
        # per-codon usage).  n1/(n1+n2) for g1, n2/(n1+n2) for g2.
        n1 <- length(g1.idx); n2 <- length(g2.idx)
        abline(h = n1 / (n1 + n2), col = adjustcolor(g1.col, 0.5), lty = 2, lwd = 0.8)
        abline(h = n2 / (n1 + n2), col = adjustcolor(g2.col, 0.5), lty = 2, lwd = 0.8)
    }

    # ---- Bin data and draw points + error bars ----
    bin.mids   <- numeric(n.bins)
    bin.counts <- integer(n.bins)

    for(i in seq_len(n.bins)) {
        if(i == 1L)          tmp.id <- expressionValues < quantiles[i]
        else if(i == n.bins) tmp.id <- expressionValues > quantiles[i]
        else                 tmp.id <- expressionValues > quantiles[i] &
                                       expressionValues < quantiles[i + 1L]
        bin.counts[i] <- sum(tmp.id, na.rm = TRUE)
        bin.mids[i]   <- median(expressionValues[tmp.id], na.rm = TRUE)

        if(split.type %in% c("ag", "tc")) {
            means <- colMeans(show.counts[tmp.id, , drop = FALSE], na.rm = TRUE)
            stds  <- apply(show.counts[tmp.id, , drop = FALSE], 2, sd, na.rm = TRUE)
            for(k in seq_along(show.codons)) {
                points(bin.mids[i], means[k], col = show.colors[k], pch = 19, cex = 0.5)
                lines(rep(bin.mids[i], 2),
                      pmax(0, pmin(1, c(means[k] - stds[k], means[k] + stds[k]))),
                      col = show.colors[k], lwd = 0.8)
            }
        } else {
            for(grp in 1:2) {
                gp   <- if(grp == 1L) g1.prop  else g2.prop
                gcol <- if(grp == 1L) g1.col   else g2.col
                gm   <- mean(gp[tmp.id], na.rm = TRUE)
                gsd  <- sd(gp[tmp.id], na.rm = TRUE)
                if(!is.na(gm)) {
                    points(bin.mids[i], gm, col = gcol, pch = 19, cex = 0.5)
                    lines(rep(bin.mids[i], 2),
                          pmax(0, pmin(1, c(gm - gsd, gm + gsd))),
                          col = gcol, lwd = 0.8)
                }
            }
        }
    }

    # ---- Model fit lines ----
    if(split.type %in% c("ag", "tc")) {
        for(k in seq_along(show.codons))
            lines(phis, show.probs[, k], col = show.colors[k])
    } else {
        lines(phis, g1.model, col = g1.col, lwd = 1.2)
        lines(phis, g2.model, col = g2.col, lwd = 1.2)
    }

    # ---- Legend ----
    if(split.type %in% c("ag", "tc")) {
        legend("topleft", legend = show.lbl, col = show.colors, lty = 1, cex = 0.75,
               bty = "o", bg = adjustcolor("white", alpha.f = 0.7), box.lwd = 0.5)
    } else {
        legend("topleft", legend = c(g1.lbl, g2.lbl),
               col = c(g1.col, g2.col), lty = 1, cex = 0.75,
               bty = "o", bg = adjustcolor("white", alpha.f = 0.7), box.lwd = 0.5)
    }

    invisible(list(xlimit = xlimit, bin.mids = bin.mids, bin.counts = bin.counts,
                   aa.bin.totals = aa.bin.totals.ret, phi.breaks = phi.breaks.ret,
                   total.aa.count = total.aa.count))
}

# Internal: shared split layout loop body for ROC and FONSE.
# prob.fn: function(aa) -> codon.probability matrix
.runSplitLayout <- function(layout, main, show.date, show.gene.hist, color.codon.groups,
                             aa.include,
                             parameter, model, genome, expressionValues,
                             samples, mixture, prob.fn,
                             panels = NULL, subtitle = NULL)
{
    # Panels come from a user-supplied list of aaGroup() specs, else the named
    # preset for `layout`.  Both flow through the same builder.
    panels <- if(!is.null(panels)) .buildSplitPanels(panels) else .getSplitPanels(layout)

    # Filter and reorder by aa.include (order of aa.include = panel order).
    if(!is.null(aa.include)) {
        panels <- do.call(c, lapply(aa.include, function(aa)
            Filter(function(pd) pd$aa == aa, panels)))
    }

    # Pad each codon-count group to a row boundary so groups align to rows.
    # NULL entries in slots are drawn as empty panels (plot.new()).
    slots         <- .padGroupsToRows(panels, n.cols = 4L)
    n.slots       <- length(slots)
    last.real     <- max(c(0L, which(!vapply(slots, is.null, logical(1L)))))
    # Floor at 5 rows to match compact-v1 proportions; grows naturally for
    # layouts with enough groups to require more (e.g. split-ct = 6 rows).
    lay           <- .buildModelLayout(as.integer(n.slots), n.rows = 5L)
    n.rows.data   <- lay$n.rows.data
    last.data.row <- if(last.real > 0L) ((last.real - 1L) %/% 4L) + 1L else n.rows.data
    # Page subtitle: explicit user value wins, else the preset's default.
    subtitle <- if(!is.null(subtitle)) subtitle else switch(layout,
        "split-ag"      = "Purine (A/G) Wobble Codons",
        "split-ct"      = "Pyrimidine (C/T) Wobble Codons",
        "split-ry"      = "Purine vs Pyrimidine (R/Y) Wobble Codons",
        "split-6codon"  = "6-Codon AA Block Structure (Leu, Arg)",
        NULL)

    # Color map: codon-count group -> border/axis color.
    grp.colors <- c("2" = "steelblue",
                    "3" = "darkorange",
                    "4" = "forestgreen",
                    "6" = "purple4")

    # Right oma must hold, for the full col-4 rows: tick numbers + the offset
    # gap + the rotated group label, without clipping at the device edge.
    par(oma = c(0, 2, 0, 7), mgp = c(3, 0.56, 0))
    layout(lay$mat, lay$widths, lay$heights, respect = FALSE)

    # title
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    # Vertical gaps between title / subtitle / timestamp:
    # title->subtitle 0.21, subtitle->date 0.13.
    text(0.5, 0.92, main, font = 2)
    if(!is.null(subtitle)) text(0.5, 0.71, subtitle, cex = 0.75, font = 3)
    if(show.date) text(0.5, 0.58, date(), cex = 0.5)

    quantiles     <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05),
                              na.rm = TRUE)
    xlimit.global <- range(expressionValues, na.rm = TRUE)
    total.aa.bin.counts <- NULL
    phi.breaks.global   <- NULL
    total.aa.all        <- 0L

    grp.labels <- c("2" = "2 Codon AAs",
                    "3" = "3 Codon AAs",
                    "4" = "4 Codon AAs",
                    "6" = "6 Codon AAs")

    # Group membership per slot (NA for blank slots).
    slot.grps <- vapply(slots, function(s)
        if(is.null(s)) NA_character_ else s$codon.group, character(1L))

    # Index of the last real panel in each data row: label is drawn from there,
    # so par("fig")[2] gives the right edge of whichever column it falls in.
    n.total.rows <- ceiling(length(slots) / 4L)
    last.real.per.row <- vapply(seq_len(n.total.rows), function(r) {
        idxs <- ((r - 1L) * 4L + 1L) : min(r * 4L, length(slots))
        real <- idxs[!is.na(slot.grps[idxs])]
        if(length(real) > 0L) real[length(real)] else NA_integer_
    }, integer(1L))
    last.real.per.row <- last.real.per.row[!is.na(last.real.per.row)]

    # Labels are collected here and drawn from the y-label panel at the end,
    # where the coordinate context is unambiguous.
    row.label.list <- list()

    for(j in seq_along(slots)) {
        pd      <- slots[[j]]
        col.pos <- ((j - 1L) %% 4L) + 1L
        row.pos <- ((j - 1L) %/% 4L) + 1L

        if(is.null(pd)) {
            par(mar = c(0, 0, 0, 0))
            plot.new()
            next
        }

        # Per-panel color: group color when color.codon.groups, else black.
        grp.col <- if(color.codon.groups) grp.colors[[pd$codon.group]] else "black"

        codon.probability <- prob.fn(pd$aa)
        result <- plotWobbleSplitPanel(parameter, model, genome, expressionValues,
                                       samples, mixture, pd, codon.probability,
                                       precomputed.quantiles = quantiles,
                                       show.gene.hist = show.gene.hist)
        if(!is.null(result$aa.bin.totals)) {
            if(is.null(total.aa.bin.counts)) {
                total.aa.bin.counts <- result$aa.bin.totals
                phi.breaks.global   <- result$phi.breaks
            } else {
                total.aa.bin.counts <- total.aa.bin.counts + result$aa.bin.totals
            }
            total.aa.all <- total.aa.all + result$total.aa.count
        }

        # Box and axis color reflects codon-count group.
        box(col = grp.col)

        # Panel label (+ optional sublabel)
        cex.lbl <- 1.0; cex.sub <- 0.6
        usr <- par("usr"); pin <- par("pin")
        tick.y  <- 0.02 * min(pin) / pin[2L] * (usr[4L] - usr[3L])
        lbl.h   <- strheight(pd$label, cex = cex.lbl, font = 2)
        sub.h   <- if(!is.null(pd$sublabel)) strheight(pd$sublabel, cex = cex.sub) else 0
        gap     <- if(!is.null(pd$sublabel)) lbl.h * 0.15 else 0
        total.h <- lbl.h + gap + sub.h
        lbl.y   <- usr[4L] - 4 * tick.y - total.h / 2 + (total.h - lbl.h) / 2
        sub.y   <- lbl.y - lbl.h / 2 - gap - sub.h / 2
        lbl.x   <- mean(xlimit.global)
        pad.x   <- max(strwidth(pd$label, cex = cex.lbl, font = 2),
                       if(!is.null(pd$sublabel)) strwidth(pd$sublabel, cex = cex.sub) else 0) * 0.6
        pad.y   <- total.h * 0.35
        rect(lbl.x - pad.x, lbl.y - lbl.h / 2 - gap - sub.h - pad.y,
             lbl.x + pad.x, lbl.y + lbl.h / 2 + pad.y,
             col = adjustcolor("white", alpha.f = 0.75), border = NA)
        text(lbl.x, lbl.y, pd$label, cex = cex.lbl, font = 2)
        if(!is.null(pd$sublabel))
            text(lbl.x, sub.y, pd$sublabel, cex = cex.sub, font = 1)

        # Labeled axes: left col, right col, top row, bottom row.
        if(col.pos == 1L)
            axis(2, las = 1, at = seq(0, 1, by = 0.2),
                 col = grp.col, col.ticks = grp.col, col.axis = grp.col)
        # Labeled right axis on the LAST real panel of EVERY row (not just the
        # full col-4 rows).  Partial rows always have blank slots to the right,
        # so the tick numbers have room.  Giving every row's last panel the same
        # numbered axis means the group label sits a constant distance past
        # identical tick numbers regardless of how many columns the row has --
        # this is what makes the label placement look consistent.
        if(j %in% last.real.per.row)
            axis(4, las = 1, at = seq(0, 1, by = 0.2),
                 col = grp.col, col.ticks = grp.col, col.axis = grp.col)
        if(row.pos == 1L)
            axis(3, col = grp.col, col.ticks = grp.col, col.axis = grp.col)
        if(row.pos == last.data.row)
            axis(1, col = grp.col, col.ticks = grp.col, col.axis = grp.col)

        # Inner ticks (color coded, no labels).
        if(show.gene.hist) {
            axis(1, tck = 0.02, labels = FALSE, pos = 0, lwd = 0, lwd.ticks = 0.5,
                 col.ticks = grp.col)
            axis(2, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2),
                 col = grp.col, col.ticks = grp.col)
            axis(3, tck = 0.02, labels = FALSE,
                 col = grp.col, col.ticks = grp.col)
            axis(4, tck = 0.02, labels = FALSE, at = seq(0, 1, by = 0.2),
                 col = grp.col, col.ticks = grp.col)
        } else {
            axis(1, tck = 0.02, labels = FALSE, col = grp.col, col.ticks = grp.col)
            axis(2, tck = 0.02, labels = FALSE, col = grp.col, col.ticks = grp.col)
            axis(3, tck = 0.02, labels = FALSE, col = grp.col, col.ticks = grp.col)
            axis(4, tck = 0.02, labels = FALSE, col = grp.col, col.ticks = grp.col)
        }
        axis(2, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2),
             col = grp.col, col.ticks = grp.col)
        axis(4, tck = 0.01, labels = FALSE, at = seq(0.1, 0.9, by = 0.2),
             col = grp.col, col.ticks = grp.col)

        # Collect row label info; drawn from the y-label panel after the loop.
        # BOTH coordinates come straight from layout geometry, never from
        # per-panel par("fig")/par("usr") (which vary panel-to-panel and made
        # placement inconsistent).
        #   x.ndc: right edge of data column `col.pos` + a fixed offset that
        #          clears the (now always present) right-axis tick numbers.
        #          Column weights are y-label=3 then 8 per data col, total 35.
        #   y.ndc: vertical centre of data row `row.pos`.  Layout heights are
        #          title=3, each data row=8, x-label=4 => total 7+8*n.rows.data.
        #          oma top/bottom are 0 lines, so no outer-margin correction.
        if(j %in% last.real.per.row) {
            omi   <- par("omi"); din <- par("din")
            l.ndc <- omi[2L] / din[1L]
            r.ndc <- omi[4L] / din[1L]
            col.right.ndc <- l.ndc + (3L + 8L * col.pos) / 35L *
                                     (1.0 - l.ndc - r.ndc)
            offset.ndc <- 4.0 * par("csi") / din[1L]
            y.ndc <- (8L + 8L * (n.rows.data - row.pos)) /
                     (7L + 8L * n.rows.data)
            row.label.list[[length(row.label.list) + 1L]] <- list(
                x.ndc = col.right.ndc + offset.ndc,
                y.ndc = y.ndc,
                grp   = pd$codon.group,
                col   = grp.col
            )
        }
    }

    # x-label
    par(mar = c(0, 0, 0, 0), xpd = FALSE)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.5, 0.55, "Gene Expression", cex = 1.4, font = 2)
    text(0.5, 0.28,
         expression(bold(log[10]~"(Protein Synthesis Rate"~phi~")")), cex = 0.9)

    # y-label (drawn last; also used as drawing surface for row group labels)
    par(mar = c(0, 0, 0, 0), xpd = NA)
    plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "")
    text(0.2, 0.5, "Proportion", srt = 90, cex = 1.4, font = 2)

    # Draw group labels collected during the panel loop.  Converting stored NDC
    # coords to this panel's user coords keeps placement device-independent.
    for(info in row.label.list) {
        x.usr <- grconvertX(info$x.ndc, "ndc", "user")
        y.usr <- grconvertY(info$y.ndc, "ndc", "user")
        text(x.usr, y.usr, grp.labels[[info$grp]],
             srt = 270, font = 2L, col = info$col, cex = 1.0,
             adj = c(0.5, 0.5), xpd = NA)
    }
}

