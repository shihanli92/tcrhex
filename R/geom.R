#' Hexagon polygons with an optional texture
#'
#' A [ggplot2::geom_polygon()] that also understands a `pattern` aesthetic:
#' `"none"`, `"stripe"`, `"crosshatch"` or `"circle"` (dots; `"dots"` is an
#' alias), drawn as near-white lines / dots clipped to each (convex) polygon,
#' and drawn in the legend keys as well. No dependency on ggpattern; its
#' `scale_pattern_manual()` can be used with it too.
#'
#' @inheritParams ggplot2::geom_polygon
#' @param pattern_colour,pattern_linewidth,pattern_spacing Texture line
#'   colour, width (pt) and spacing (data units for the hexagons).
#' @return A ggplot2 layer.
#' @export
geom_hextile <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", ..., pattern_colour = "#fcfcfb",
                         pattern_linewidth = 1.1, pattern_spacing = 0.62,
                         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomHexTile,
                 position = position, show.legend = show.legend,
                 inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm, pattern_colour = pattern_colour,
                               pattern_linewidth = pattern_linewidth,
                               pattern_spacing = pattern_spacing, ...))
}

#' @rdname geom_hextile
#' @format NULL
#' @usage NULL
#' @export
GeomHexTile <- ggplot2::ggproto("GeomHexTile", ggplot2::GeomPolygon,
  default_aes = ggplot2::aes(colour = "#f3f2ee", fill = "#2a78d6", linewidth = 0.28,
                             linetype = 1, alpha = NA, subgroup = NULL,
                             pattern = "none"),
  extra_params = c("na.rm", "pattern_colour", "pattern_linewidth", "pattern_spacing"),

  draw_panel = function(self, data, panel_params, coord, pattern_colour = "#fcfcfb",
                        pattern_linewidth = 1.1, pattern_spacing = 0.62,
                        lineend = "butt", linejoin = "round", linemitre = 10) {
    poly <- ggplot2::GeomPolygon$draw_panel(data, panel_params, coord, lineend = lineend,
                                            linejoin = linejoin, linemitre = linemitre)
    pat <- if (is.null(data$pattern)) rep("none", nrow(data)) else as.character(data$pattern)
    data$pattern <- .norm_pattern(pat)
    pd <- data[data$pattern != "none", , drop = FALSE]
    if (!nrow(pd)) return(poly)
    segs <- list()
    dots <- list()
    for (g in unique(pd$group)) {
      d <- pd[pd$group == g, ]
      pt <- d$pattern[1]
      s <- .pattern_geom(d$x, d$y, pt, pattern_spacing)
      if (!is.null(s$segments)) segs[[length(segs) + 1]] <- s$segments
      if (!is.null(s$dots)) dots[[length(dots) + 1]] <- s$dots
    }
    gp <- grid::gpar(col = pattern_colour, lwd = pattern_linewidth / 0.75, lineend = "butt")
    kids <- list(poly)
    if (length(segs)) {
      s <- do.call(rbind, segs)
      a <- coord$transform(data.frame(x = s$x0, y = s$y0), panel_params)
      b <- coord$transform(data.frame(x = s$x1, y = s$y1), panel_params)
      kids[[length(kids) + 1]] <- grid::segmentsGrob(a$x, a$y, b$x, b$y,
                                                     default.units = "native", gp = gp)
    }
    if (length(dots)) {
      dd <- do.call(rbind, dots)
      a <- coord$transform(data.frame(x = dd$x, y = dd$y), panel_params)
      kids[[length(kids) + 1]] <- grid::pointsGrob(a$x, a$y, pch = 16,
        size = grid::unit(pattern_linewidth * 1.6, "pt"),
        default.units = "native", gp = grid::gpar(col = pattern_colour))
    }
    grid::gTree(children = do.call(grid::gList, kids))
  },

  draw_key = function(data, params, size) {
    r <- 0.5
    ang <- pi / 2 - (0:5) * pi / 3
    hx <- 0.5 + r * cos(ang) * 0.95
    hy <- 0.5 + r * sin(ang)
    fill <- if (is.na(data$alpha %||% NA)) data$fill else grDevices::adjustcolor(data$fill, data$alpha)
    poly <- grid::polygonGrob(hx, hy, gp = grid::gpar(fill = fill, col = data$colour %||% NA,
                                                      lwd = (data$linewidth %||% 0.28) * ggplot2::.pt))
    pt <- .norm_pattern(as.character(data$pattern %||% "none"))
    if (pt == "none") return(poly)
    s <- .pattern_geom(hx, hy, pt, 0.3)
    pc <- params$pattern_colour %||% "#fcfcfb"
    lw <- (params$pattern_linewidth %||% 1.1) * 0.6
    kids <- list(poly)
    if (!is.null(s$segments)) kids[[2]] <- grid::segmentsGrob(s$segments$x0, s$segments$y0,
      s$segments$x1, s$segments$y1, gp = grid::gpar(col = pc, lwd = lw / 0.75))
    if (!is.null(s$dots)) kids[[2]] <- grid::pointsGrob(s$dots$x, s$dots$y, pch = 16,
      size = grid::unit(lw * 1.6, "pt"), gp = grid::gpar(col = pc))
    grid::gTree(children = do.call(grid::gList, kids))
  }
)

.norm_pattern <- function(p) {
  p <- tolower(p)
  p[is.na(p)] <- "none"
  p[p %in% c("dots", "dot", "points")] <- "circle"
  p[p %in% c("stripes", "lines", "hatch")] <- "stripe"
  p[!p %in% c("none", "stripe", "crosshatch", "circle")] <- "none"
  p
}

# Texture of one convex polygon (vertices x, y): stripes at 45 degrees,
# crosshatch at +-45, or a dot lattice; clipped to the polygon shrunk by
# `inset` (fraction of the spacing).
.pattern_geom <- function(x, y, pattern, spacing, inset = 0.45) {
  cx <- mean(x)
  cy <- mean(y)
  n <- length(x)
  # half-planes of the polygon: for edge i -> i+1, outward normal
  x2 <- c(x[-1], x[1])
  y2 <- c(y[-1], y[1])
  ex <- x2 - x
  ey <- y2 - y
  len <- sqrt(ex^2 + ey^2)
  nx <- ey / len
  ny <- -ex / len
  # orient outward
  sgn <- sign(nx * (x - cx) + ny * (y - cy))
  nx <- nx * sgn
  ny <- ny * sgn
  off <- nx * (x - cx) + ny * (y - cy) - inset * spacing
  rmax <- max(sqrt((x - cx)^2 + (y - cy)^2))
  clip_line <- function(px, py, ux, uy) {
    lo <- -Inf
    hi <- Inf
    for (i in seq_len(n)) {
      a <- nx[i] * ux + ny[i] * uy
      b <- off[i] - (nx[i] * (px - cx) + ny[i] * (py - cy))
      if (abs(a) < 1e-12) {
        if (b < 0) return(NULL)
      } else if (a > 0) hi <- min(hi, b / a) else lo <- max(lo, b / a)
    }
    if (hi <= lo) return(NULL)
    c(px + lo * ux, py + lo * uy, px + hi * ux, py + hi * uy)
  }
  stripes <- function(theta) {
    ux <- cos(theta)
    uy <- sin(theta)
    vx <- -uy
    vy <- ux
    ts <- seq(-rmax, rmax, by = spacing)
    ts <- ts - mean(ts)
    out <- lapply(ts, function(t) clip_line(cx + t * vx, cy + t * vy, ux, uy))
    out <- do.call(rbind, out[!vapply(out, is.null, TRUE)])
    if (is.null(out)) return(NULL)
    data.frame(x0 = out[, 1], y0 = out[, 2], x1 = out[, 3], y1 = out[, 4])
  }
  if (pattern == "stripe") return(list(segments = stripes(pi / 4)))
  if (pattern == "crosshatch") return(list(segments = rbind(stripes(pi / 4), stripes(-pi / 4))))
  if (pattern == "circle") {
    s <- spacing * 1.05
    g <- expand.grid(i = -6:6, j = -6:6)
    gx <- cx + s * (g$i + g$j / 2)
    gy <- cy + s * g$j * sqrt(3) / 2
    inside <- rep(TRUE, length(gx))
    for (i in seq_len(n)) inside <- inside & (nx[i] * (gx - cx) + ny[i] * (gy - cy) <= off[i])
    return(list(dots = data.frame(x = gx[inside], y = gy[inside])))
  }
  list()
}
