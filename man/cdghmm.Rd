\name{cdghmm}
\alias{cdghmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hidden Markov Models for Multivariate Panel Data}
\description{
Estimates hidden Markov models from the CDGHMM family under various missingness schemes.
}
\usage{
cdghmm(x,m,id,mu=NULL,sigma=NULL,gamma=NULL,delta=NULL,alpha=NULL,beta=NULL,
                        maxiter=10000,tol=1e-6,type="s",covtype="VVA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
Data frame or matrix to perform variable selection on
}
\item{m}{
Number to indicate the number of states to fit.
}
\item{id}{
A vector of indicators to indicate observational unit.
}
\item{mu}{
An \code{m}-length list of matrices to be used as an initial estimate for \code{mu}. If no initial \code{mu} is provided the algorithm will initialize via k-means.
}
  \item{sigma}{
An \code{m}-length list of matrices to be used as an initial estimate for sigma. If no initial sigma is provided the algorithm will initialize.
}
\item{gamma}{
A \code{m}x\code{m} matrix to be used as an initial estimate for \code{gamma}. If no initial gamma is provided the algorithm will initialize.
}
\item{delta}{
A vector to be used as an initial estimate for \code{delta}. If no initial \code{delta} is provided the algorithm will initialize.
}
\item{alpha}{
A \code{m}x\code{p}x\code{t} array to be used as an initial estimate for \code{alpha}. If no initial \code{alpha} is provided the algorithm will initialize.
}
\item{beta}{
A \code{m}x\code{p}x\code{t} array to be used as an initial estimate for \code{beta}. If no initial \code{beta} is provided the algorithm will initialize.
}
\item{maxiter}{
A number to indicate the maximum number of iterations allowed, default is \code{10000}.
}
\item{tol}{
A number to indicate the tolerance value, default is \code{1e-6}.
}
\item{type}{
A character to indicate which type of missingness mechanism to use. The allowed values are:
\code{"mar"} (missing at random), \code{"s"} (3.2.1 in cited paper), \code{"sv"} (3.2.2), \code{"st"} (3.2.3 model without beta), \code{"svt"} (3.2.4 model without beta), \code{"st2"} (3.2.3 model with beta), \code{"svt2"} (3.2.4 model with beta). The default is \code{type="s"}.
}
\item{covtype}{
A string to indicate which covariance estimate to use. The allowed values are:
\code{EEA}, \code{VVA}, \code{VEA}, \code{EVA}, \code{VVI}, \code{VEI}, \code{EVI}, \code{EEI}. The default is \code{covtype="VVA"}.
}
}
\value{
\item{mu}{
Estimated mean matrices.
}
\item{sigma}{
Estimated covariance matrices.
}
\item{gamma}{
Estimated gamma matrix.
}
\item{delta}{
Estimated delta vector.
}
\item{alpha}{
Estimated alpha missingness parameters.
}
\item{beta}{
Estimated beta missingness parameters.
}
\item{llk}{
Estimated log-likelihood.
}
\item{AIC}{
The value of the Akaike information criterion.
}
\item{BIC}{
The value of the Bayes information criterion.
}
\item{ICL}{
The value of the integrated completed likelihood.
}
\item{Avg_Silouette}{
The value of the average silhouette score.
}
\item{probs}{
A matrix whose entries correspond to the probability of belonging to a state.
}
\item{states}{
Estimated states via map(probs).
}
\item{mod}{
The CDGHMM family member fit.
}
}
\references{
See \code{citation("CDGHMM")}.
}
\author{
Mackenzie R. Neal, Alexa A. Sochaniwsky, Paul D. McNicholas
}
\examples{
\donttest{
data("simulated_data")
id=simulated_data$V5
x <- simulated_data[,1:4]
EEI_mod=cdghmm(x,2,id=id,covtype="EEI",tol=1e-4)
table(simulated_data$V7,EEI_mod$states)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
