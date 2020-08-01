\name{getCC}
\alias{getCC}
\title{Get Phase II EWMA corrected charting constant using the CUC or EPC method}
\usage{
getCC(
    m,
    nu = m - 1,
	n = 5,
    ARL0 = 370,
    interval = c(1, 4),
	lambda = 0.1,
    EPC.p0 = 0.05,
    EPC.epstilda = 0,
    cc.option = c('EPC'),
	ubCons = 1, 
	t.interval = c(20, 50),
	state = 'steady',
	tol = 1e-2
)
}
\arguments{
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{n}{ sample size }
  \item{ARL0}{ the nominal value of average run length }
  \item{interval}{ the searching range of charting constants.  The default setting is between 1 and 4. }
  \item{lambda}{ smoothing constant }
  \item{EPC.p0}{ the nominal value of probability for the EPC method.  It is not recommended using values greater than 0.5.  The default setting is 0.05 }
  \item{EPC.epstilda}{ the nominal value of tolerance factor for the EPC method.  The default setting is 0.}
  \item{cc.option}{ the charting constant option.  In the default setting, the charting constant is calculated using the EPC method}
  \item{ubCons}{ the unbiasing constant.  The default setting is 1. }
  \item{t.interval}{ searching range of number of states in Markov chain.  The default setting is between 20 and 50. }
  \item{state}{ steady state or zero state.  The default setting is "steady".  Note that using zero state reduces efficiency significantly }
  \item{tol}{ tolerance level for the expectation of the conditional average run length.  The default setting is 1e-2. }
}
\description{
Obtain a Phase II EWMA corrected charting constant using the CUC or EPC method
}
\examples{
%\dontrun{

ub <- c4.f(49)

# get the charting constant using the CUC method
getCC(m = 50, nu = 49, cc.option = 'CUC', ubCons = ub)

# get the charting constant using the EPC method
getCC(m = 50, nu = 49, cc.option = 'EPC', ubCons = ub)

%}
}
