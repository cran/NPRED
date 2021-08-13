## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>",
  
  out.width='100%'
)

## ----setup--------------------------------------------------------------------
library(NPRED)

op <- par()
require(zoo)
require(ggplot2)

## ----dat, fig.cap='Example of AR models', fig.height=5, fig.width=9-----------
set.seed(2020)
# AR1 model from paper with 9 dummy variables
data.ar1 <- data.gen.ar1(500)


# AR4 model from paper with total 9 dimensions
data.ar4 <- data.gen.ar4(500)


# AR9 model from paper with total 9 dimensions
data.ar9 <- data.gen.ar9(500)

plot.zoo(cbind(data.ar1$x, data.ar4$x, data.ar9$x),
  xlab = NA, main = "Example of AR models",
  ylab = c("AR1", "AR4", "AR9"))

## ----pic----------------------------------------------------------------------
# R version of pic
pic.calc(data.ar9$x, data.ar9$dp)


# fortran version of pic
calc.PIC(data.ar9$x, data.ar9$dp)


## ----fig, fig.cap='Example of PIC selection implemented in NPRED', fig.height=5, fig.width=9----
# pic <- stepwise.PIC(data.ar1$x, data.ar1$dp)
# pic <- stepwise.PIC(data.ar4$x, data.ar4$dp)
pic <- stepwise.PIC(data.ar9$x, data.ar9$dp)
pic

## ----pw-----------------------------------------------------------------------
# R version of pw
pw.calc(data.ar9$x, data.ar9$dp, pic$cpy, pic$cpyPIC)

# fortran version of pw
calc.PW(data.ar9$x, data.ar9$dp, pic$cpy, pic$cpyPIC)


## ----knn, fig.cap='Example of KNN implemented in NPRED', fig.height=5, fig.width=9----
data("data3")
x <- ts(data3[, 1]) # response
z <- ts(data3[, -1]) # possible predictors
zout <- ts(data.gen.ar1(500, ndim = 15)$dp) # new input

xhat1 <- xhat2 <- x
# xhat1 <- NPRED::knn(x,z,zout,k=5,reg=T,extrap=F)
# xhat2 <- NPRED::knn(x,z,zout,k=5,reg=T,extrap=T)

for (i in 1:500) {
  xhat1[i] <- NPRED::knn(x[-i], z[-i, ], z[i, ], extrap = F)
  xhat2[i] <- NPRED::knn(x[-i], z[-i, ], z[i, ], extrap = T)
}

if (TRUE) {
  par(mfrow = c(1, 2), pty = c("s"))

  ts.plot(x, xhat1, xhat2, col = c("black","red","blue"), ylim = c(-10, 10), lwd = c(1, 1, 1))
  legend("topleft",
    bty = "n", lwd = 3, cex = 1, lty = 1,
    # inset = c(-0.5, 0),
    legend = c("OBS", "Pred", "Pred(extrap=T)"),
    x.intersp = 0, xjust = 0, yjust = 0, text.width = c(0, 50, 50), horiz = T,
    col = c("black","red","blue")
  )

  plot(xhat1, xhat2, xlim = c(-10, 10), ylim = c(-10, 10))
  abline(coef = c(0, 1), lwd = 1, col = 2)
}

## ----weights, fig.cap='Illustration of the usage of partial weights', fig.height=6, fig.width=9, out.width='100%'----
par(op)
sample <- 500
k <- 0
u <- runif(sample, 0, 5 * pi)
z <- sin(u) + rnorm(sample, sd = 0.2)

u1 <- cbind(u, runif(sample, 0, 5 * pi), runif(sample, 0, 5 * pi), runif(sample, 0, 5 * pi))
# zhat1 <- knnregl1cv(x=z, z=u1, k=k)
zhat1 <- sapply(1:sample, function(i) knn(x = z[-i], z = u1[-i, ], zout = u1[i, ], k = k))

sel <- stepwise.PIC(x = z, py = u1)
sel
# zhat2 <- knnregl1cv(x=z, z=u1[,sel$cpy], k=k)
zhat2 <- sapply(1:sample, function(i) knn(x = z[-i], z = u1[-i, sel$cpy], zout = u1[i, sel$cpy], k = k))

if (TRUE) {
  plot(u, z, pch = 16)
  lines(sort(u), zhat1[order(u)], col = "green")
  lines(sort(u), zhat2[order(u)], col = "red")
  abline(a = 0, b = 0)
}
