\name{PH2EWMA}
\alias{PH2EWMA}
\title{Build Phase II EWMA control chart with a corrected charting constant}
\usage{
PH2EWMA(
  X2,
  X1,
  cc = NULL,
  ARL0 = 370,
  interval = c(1, 4),
  lambda = 0.1,
  EPC.p0 = 0.05,
  EPC.epstilda = 0,
  cc.option = c('CUC', 'EPC'),
  t.interval = c(20, 50),
  state = 'steady',
  tol = 1e-2,
  ubCons.option = TRUE,
  plot.option = TRUE)
}
\arguments{
  \item{X2}{ Phase II sample }
  \item{X1}{ Phase I sample }
  \item{cc}{ given charting constant.  If cc = NULL, it will be calculated using the CUC and EPC methods. }
  \item{ARL0}{ the nominal value of average run length }
  \item{interval}{ the searching range of charting constants.  The default setting is between 1 and 4. }
  \item{lambda}{ smoothing constant }
  \item{EPC.p0}{ the nominal value of probability for the EPC method.  It is not recommended using values greater than 0.5.  The default setting is 0.05 }
  \item{EPC.epstilda}{ the nominal value of tolerance factor for the EPC method.  The default setting is 0.}
  \item{cc.option}{ the charting constant option.  In the default setting, the charting constant is calculated using the EPC method}
  \item{t.interval}{ searching range of number of states in Markov chain.  The default setting is between 20 and 50. }
  \item{state}{ steady state or zero state.  The default setting is "steady".  Note that using zero state reduces efficiency significantly }
  \item{tol}{ tolerance level for the expectation of the conditional average run length.  The default setting is 1e-2. }
  \item{ubCons.option}{ the unbiasing constant option.  For the default setting, the unbiasing constant is invovlved and c4 is applied. }
  \item{plot.option}{ the plot option.  The default setting is to graph the Phase II control chart. }
}
\description{
Build a Phase II Shewhart EWMA control chart.
}
\examples{
%\dontrun{
# Build a Phase II basic Shewhart control chart

data(Ph1data)
data(Ph2data)

X1 <- as.matrix(Ph1data[, 2:4]) ^ (1/3)
X2 <- as.matrix(Ph2data[, 2:4]) ^ (1/3)
X2[which(is.na(X2))] <- mean(X1)

PH2EWMA(
  X2 = X2,
  X1 = X1,
  cc = NULL,
  ARL0 = 370,
  interval = c(1, 4),
  lambda = 0.1,
  EPC.p0 = 0.05,
  EPC.epstilda = 0,
  cc.option = c('CUC', 'EPC'),
  t.interval = c(20, 50),
  state = 'steady',
  tol = 1e-2,
  ubCons.option = TRUE,
  plot.option = TRUE) 
%}
}
