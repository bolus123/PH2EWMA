\name{getCC.EPC}
\alias{getCC.EPC}
\title{Get Phase II EWMA corrected charting constant using the EPC method}
\usage{
getCC.EPC(
  p0,
  interval = c(1, 4),
  ARL0,
  epstilda,
  lambda,
  m,
  nu = m - 1,
  n = 5,
  ubc = 1,
  t.interval = c(20, 50),
  state = "steady",
  tol = 1e-2
)
}
\arguments{

\item{p0}{nominal value of probability.  It is not recommended using values greater than 0.5}

\item{interval}{searching range of charting constants for the CUC method.  The default setting is between 1 and 4.}

\item{ARL0}{nominal value of average run length}

\item{epstilda}{nominal value of tolerance factor}

\item{lambda}{smoothing constant}

\item{m}{number of subgroups}

\item{nu}{degrees of freedom}

\item{n}{sample size}

\item{ubc}{unbiasing constant.  The default setting is 1.}

\item{t.interval}{searching range of number of states in Markov chain.  The default setting is between 20 and 50.}

\item{state}{steady state or zero state.  The default setting is "steady".  Note that using zero state reduces efficiency significantly}

\item{tol}{tolerance level for the expectation of the conditional average run length.  The default setting is 1e-2.}
}
\value{
The sum of \code{x} and \code{y}
}
\description{
Obtain a Phase II EWMA corrected charting constant using the EPC method
}
\examples{
getCC.EPC(p0 = 0.05, interval = c(1, 5), ARL0 = 370, epstilda = 0, lambda = 0.1, m = 50)
}
