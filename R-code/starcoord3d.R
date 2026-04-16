## starcoords3D functions
## anchor points coordinates
##
## function to create the anchor points, following Zhu et al (JCGS, 2022)
##

anchors_sphere <- function(p) {
  phi <- (sqrt(5) + 1)/2
  anchor_coord <- list()
  anchor_coord[[4]] <- matrix(c(1, 1, 1, 1, -1, -1, -1, 1, -1, -1, -1, 1), byrow = T, 
                              ncol = 3)/sqrt(3)
  anchor_coord[[6]] <- matrix(c(1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, 
                                -1), byrow = T, ncol = 3)
  anchor_coord[[8]] <- matrix(c(1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 
                                1, -1, 1, -1, -1, -1, -1, -1), byrow = T, ncol = 3)/sqrt(3)
  anchor_coord[[12]] <- matrix(c(0, 1, phi, 0, 1, -phi, 0, -1, phi, 0, -1, -phi, 1, 
                                 phi, 0, -1, phi, 0, 1, -phi, 0, -1, -phi, 0, phi, 0, 1, -phi, 0, 1, phi, 0, -1, 
                                 -phi, 0, -1), byrow = T, ncol = 3)/sqrt(1 + phi^2)
  anchor_coord[[20]] <- matrix(c(1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, 
                                 -1, 1, -1, 1, -1, -1, -1, -1, -1, 0, 1/phi, phi, 0, 1/phi, -phi, 0, -1/phi, phi, 
                                 0, -1/phi, -phi, 1/phi, phi, 0, -1/phi, phi, 0, 1/phi, -phi, 0, -1/phi, -phi, 
                                 0, phi, 0, 1/phi, -phi, 0, 1/phi, phi, 0, -1/phi, -phi, 0, -1/phi), byrow = T, 
                               ncol = 3)/sqrt(3)
  
  Fibonacci_set <- function(p) {
    anchor <- matrix(nrow = p, ncol = 3)
    for (i in 1:p) {
      anchor[i, 3] <- (2 * i - 1)/p - 1
      anchor[i, 1] <- sqrt(1 - anchor[i, 3]^2) * cos(2 * pi * i * phi^{
        -1
      })
      anchor[i, 2] <- sqrt(1 - anchor[i, 3]^2) * sin(2 * pi * i * phi^{
        -1
      })
    }
    return(anchor)
  }
  
  if (p %in% c(4, 6, 8, 12, 20)) {
    anchors <- anchor_coord[[p]]
  } else {
    anchors <- Fibonacci_set(p)
  }
  return(anchors)
}


simple.starcoords3D <- function(data, class) {
  x <- as.matrix(data)
  x.range <- apply(x, 2, range)
  z <- t((t(x) - x.range[1,]) / (x.range[2,] - x.range[1,]))
  d <- dim(z)[2] # Dimensions
  prj <- anchors_sphere(d)
  z %*% prj
}


plot.starcoords3D <- function(data, class, pradius = 0.02, lwd = 1, 
                              colors = RColorBrewer::brewer.pal(name = "Dark2", n = 8),
                              pch = letters[1:26],
                              #pch = c(9818:9823, 9830, 9824, 9827, 9829),
                              cex = 1, axes.cex = 1, utf.chars = FALSE,...) {
  x <- data
  x.range <- apply(x, 2, range)
  z <- t((t(x) - x.range[1,]) / (x.range[2,] - x.range[1,]))
  d <- dim(z)[2] # Dimensions
  
  prj <- anchors_sphere(d)
  axes.prj <- diag(d) %*% prj
  
  data_trans <- z %*% prj
  rgl::open3d()
  apply(axes.prj, 1, function(v) (rgl::rgl.lines(rbind(rep(0, 3), v), # "grey40"
                                                 col = "black", lwd = lwd, ...)))
  
  ## ---- ADD AXIS LABELS ----
  var.names <- colnames(data)
  if (is.null(var.names))
    var.names <- paste0("V", 1:d)
  
  for (i in 1:d) {
    rgl::text3d(axes.prj[i, ], texts = var.names[i], 
                col = "black", cex = axes.cex, ...)
  }
  
  class.labels <- unlist(lapply(class[,drop=T], as.numeric))
  if (!utf.chars) {
      rgl::pch3d(x = data_trans, y = NULL, z = NULL, radius = pradius, 
                 color = colors[class.labels], 
                 pch = pch[class.labels], 
                 cex = cex,...)
  }
  else {
      for (i in 1:max(class.labels))
          rgl::text3d(x = data_trans[class.labels == i, ],
                      col=colors[i], text = intToUtf8(pch[i]),
                      pos = 0, cex = cex, usePlotmath = TRUE,...)
  }
  apply(axes.prj, 1, function(v) (rgl::rgl.lines(rbind(rep(0, 3), v), 
                                                 col = "gray40", lwd = lwd, ...))) ## axes redrawn, perhaps not needed
}


simple.starcoords2D <- function(data) {
  x <- as.matrix(data)
  x.range <- apply(x, 2, range)
  z <- t((t(x) - x.range[1,]) / (x.range[2,] - x.range[1,]))
  d <- dim(z)[2] # Dimensions
  prj <- t(sapply((1:d)/d, function(i) c(cos(2*pi*i + pi/4), sin(2*pi*i + pi/4))))
  z %*% prj
}

plot.starcoords2D <- function(data, class, colors = RColorBrewer::brewer.pal(name = "Dark2", n = 8), 
                              cex = 1, pch = -1*c(9818:9823, 9830, 9824, 9827, 9829),
                              utf.chars = FALSE,...) {
  star <- simple.starcoords2D(data)
  d <- dim(data)[2] # Dimensions
  prj <- t(sapply((1:d)/d, function(i) c(cos(2*pi*i + pi/4), sin(2*pi*i + pi/4))))
  plot(rbind(t(apply(star, 2, range)), apply(prj*1.25, 2, range)), type="n", bty="n", xaxt="n", yaxt="n", main="", xlab="", ylab="")
  tmp <- apply(prj, 1, function(v) lines(rbind(c(0,0), v)))
  text(prj * 1.1, labels=colnames(data), cex=0.8, col="Gray")
  if (!utf.chars) {
    points(star, pch=letters[as.numeric(class)], col= colors[unclass(class)], cex = cex)
  }else {
    points(star, pch=pch[as.numeric(class)], col= colors[unclass(class)], cex = cex)
  }
}

#plot.starcoords3D(data = iris[,-5], class = iris[,5])
