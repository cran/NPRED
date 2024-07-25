## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, collapse = TRUE, comment = "#>", warning = FALSE,
  
  out.width='100%',
  fig.align = "center", fig.pos = "h!"
)

## ----setup--------------------------------------------------------------------
library(NPRED)
require(zoo)
require(ggplot2)

## ----dat, fig.cap='Example of AR models', fig.height=5, fig.width=9-----------
# Other statistical models for generating synthetic data see Package - synthesis
require(synthesis)

set.seed(2020) # set seed for reproducible results
# AR1 model from paper with 9 dummy variables
data.ar1 <- synthesis::data.gen.ar1(500)

# AR4 model from paper with total 9 dimensions
data.ar4 <- synthesis::data.gen.ar4(500)

# AR9 model from paper with total 9 dimensions
data.ar9 <- synthesis::data.gen.ar9(500)

# plot.zoo(cbind(data.ar1$x, data.ar4$x, data.ar9$x),
#   xlab = NA, main = "Example of AR models",
#   ylab = c("AR1", "AR4", "AR9"))

## ----pic----------------------------------------------------------------------
# calculate PIC
pic.calc(data.ar1$x, data.ar1$dp)

pic.calc(data.ar4$x, data.ar4$dp)

pic.calc(data.ar9$x, data.ar9$dp)


## ----fig, fig.cap='Example of PIC selection implemented in NPRED', fig.height=5, fig.width=9----
pic1 <- stepwise.PIC(data.ar1$x, data.ar1$dp)

pic4 <- stepwise.PIC(data.ar4$x, data.ar4$dp)

pic9 <- stepwise.PIC(data.ar9$x, data.ar9$dp)


## ----pw-----------------------------------------------------------------------
# calculate partial weights
pw.calc(data.ar1$x, data.ar1$dp, pic1$cpy, pic1$cpyPIC)

pw.calc(data.ar4$x, data.ar4$dp, pic4$cpy, pic4$cpyPIC)

pw.calc(data.ar9$x, data.ar9$dp, pic9$cpy, pic9$cpyPIC)


## ----knn, fig.cap='Example of KNN implemented in NPRED', fig.height=5,fig.width=9----
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
  par(mfrow = c(1, 1),
      mar=c(3,2,1,1), mgp=c(1.5,0.5,0),
      pty = c("m"))

  ts.plot(x, xhat1, xhat2, col = c("black","red","blue"), ylim = c(-10, 10), lwd = c(1, 1, 1))
  legend("topleft",
    bty = "n", lwd = 3, cex = 1, lty = 1,
    # inset = c(-0.5, 0),
    legend = c("OBS", "Pred", "Pred(extrap=T)"),
    x.intersp = 0, xjust = 0, yjust = 0, text.width = c(0, 50, 50), horiz = T,
    col = c("black","red","blue"))

  par(mfrow = c(1, 1), pty = c("s"))
  plot(xhat1, xhat2, xlim = c(-10, 10), ylim = c(-10, 10))
  abline(coef = c(0, 1), lwd = 1, col = 2)

}

## ----weights, fig.cap='Illustration of the usage of partial weights', fig.height=6, fig.width=9----
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
  par(mfrow = c(1, 1), pty = c("m"))
  
  plot(u, z, pch = 16)
  lines(sort(u), zhat1[order(u)], col = "green")
  lines(sort(u), zhat2[order(u)], col = "red")
  abline(a = 0, b = 0)
}

## ----exp, warning=FALSE, fig.cap='Application to predicting drought anomalies', fig.height=9, fig.width=7----
# An demo example used in Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020)
# Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2021). A wavelet-based tool to modulate variance in predictors: An application to predicting drought anomalies. Environmental Modelling and Software, 135, 104907.

require(WASP)
#load response and predictor variables
data("SPI.12", package="WASP"); 
data("data.CI", package="WASP")
data("Ind_AWAP.2.5", package="WASP")

#study grids and period
Grid <- sample(Ind_AWAP.2.5,1)
Grid = 149 #Grid used in the SI
SPI.12.ts <- window(SPI.12, start=c(1910,1),end=c(2009,12))
data.CI.ts <- window(data.CI, start=c(1910,1),end=c(2009,12))
#partition into two folds
folds <- cut(seq(1,nrow(SPI.12.ts)),breaks=2,labels=FALSE)
sub.cali <- which(folds==1, arr.ind=TRUE); sub.vali <- which(folds==2, arr.ind=TRUE)
#-------------------------------------------------------------------------------------
###calibration and selection
data <- list(x=SPI.12.ts[sub.cali,Grid],dp=data.CI.ts[sub.cali,])

#variance transformation - calibration
modwt <- modwt.vt(data, wf="haar", J=8, boundary="periodic")

#stepwise PIC selection
sel <- NPRED::stepwise.PIC(modwt$x, modwt$dp.n)
#-------------------------------------------------------------------------------------
###validation and prediction
data.val <- list(x=SPI.12.ts[sub.vali,Grid],dp=data.CI.ts[sub.vali,])

#variance transformation - validation
modwt.val <- modwt.vt.val(data.val, J=8, modwt)

#knn prediction
cpy <- sel$cpy; wt <- sel$wt
x=data$x; z=modwt$dp.n[,cpy]; zout=modwt.val$dp.n[,cpy]
mod <- knn(x, z, zout, k=5, pw=wt, extrap=T)

#-------------------------------------------------------------------------------------
###plot
start.cal <- c(1910,1); start.val <- c(1960,1)
ndim = ncol(data.CI.ts); CI.names = colnames(data.CI.ts)
par(mfcol=c(ndim+1,2),mar=c(2,4,2,2),mgp=c(1.5,0.5,0),
    bg = "white",pty="m",cex.lab=1.5,ps=8)
#----------------------------------------------
#plot before and after vt - calibration
if(TRUE){

  x <- ts(modwt$x, start=start.cal, frequency = 12)
  dp <- ts(modwt$dp, start=start.cal, frequency = 12)
  dp.n <- ts(modwt$dp.n, start=start.cal, frequency = 12)

  ts.plot(x,xlab=NA, main=paste0("Sampled Grid: ", Grid), ylab=paste0("SPI",12), col=c("black"),lwd=c(1))
  for(nc in 1:ndim)
    ts.plot(dp[,nc],dp.n[,nc],xlab=NA,ylab=paste0(CI.names[nc]),
            col=c("red","blue"),lwd=c(1,1),lty=c(1,2))

}
#----------------------------------------------
#plot before and after vt - validation
if(TRUE){

  x <- ts(modwt.val$x, start=start.val, frequency = 12)
  dp <- ts(modwt.val$dp, start=start.val, frequency = 12)
  dp.n <- ts(modwt.val$dp.n, start=start.val, frequency = 12)

  ts.plot(x, xlab=NA, main=paste0("Sampled Grid: ",Grid), ylab=paste0("SPI",12), col=c("black"),lwd=c(1))
  for(nc in 1:ndim)
    ts.plot(dp[,nc],dp.n[,nc],xlab=NA,ylab=paste0(CI.names[nc]),
            col=c("red","blue"),lwd=c(1,1),lty=c(1,2))

}


