% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.PICPW.R
\name{stepwise.PIC}
\alias{stepwise.PIC}
\title{Calculate stepwise PIC}
\usage{
stepwise.PIC(x, py, nvarmax = 100, alpha = 0.1)
}
\arguments{
\item{x}{A vector of response.}

\item{py}{A matrix containing possible predictors of x.}

\item{nvarmax}{The maximum number of variables to be selected.}

\item{alpha}{The significance level used to judge whether the sample estimate in Equation \deqn{\hat{PIC} = sqrt(1-exp(-2\hat{PI})} is significant or not. A default alpha value is 0.1.}
}
\value{
A list of 2 elements: the column numbers of the meaningful predictors (cpy), and partial informational correlation (cpyPIC).
}
\description{
Calculate stepwise PIC
}
\examples{
\donttest{
data(data1) # AR9 model   x(i)=0.3*x(i-1)-0.6*x(i-4)-0.5*x(i-9)+eps
x <- data1[, 1] # response
py <- data1[, -1] # possible predictors
stepwise.PIC(x, py)

data(data2) # AR4 model:  x(i)=0.6*x(i-1)-0.4*x(i-4)+eps
x <- data2[, 1] # response
py <- data2[, -1] # possible predictors
stepwise.PIC(x, py)

data(data3) # AR1 model  x(i)=0.9*x(i-1)+0.866*eps
x <- data3[, 1] # response
py <- data3[, -1] # possible predictors
stepwise.PIC(x, py)
}
}
\references{
Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
}
