\name{spsl}
\alias{spsl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  Function to exploit Independent SPSL model to estimate the dynamic gene coexpression among all the gene pairs systematically.
Fitting SPSL model
}
\description{
spsl is used to estimate the dynamic gene coexpression among all gene pairs with independent spike-and-slab prior in a Bayesian variable selection framework.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spsl(x, z, iterations, burnIn, chains, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The matrix of gene expression data, where each row represents the genes of a patient and each column represents a gene.
%%     ~~Describe \code{x} here~~
}
\item{z}{
A vector of independent variable, where each element represents the cancer recurrence-free survival probability of a patient.
}
\item{iterations}{Number of iterations in MCMC}
\item{burnIn}{Length of burn-in period}
\item{chains}{Number of the MCMC chains}
\item{startvalue}{Initial values of MCMC. The order of the array of initial value is (group1, group2, ...). In each group, the values are in the order of tau0, tau1, phi and w. Since transformations are implemented in the sampling process, the initial values for all parameters can be set in the range from -Inf to Inf.}
\item{t0}{History tracing size of MCMC samples used for calculating the covariance matrix.}
\item{v0}{Shrinkage parameter used in the mixture of inverse gamma distribution. It controls the shrinkage effect in spike-and-slab prior.}
\item{a1}{Shape parameter used in the inverse gamma distribution}
\item{a2}{Scale parameter used in the inverse gamma distribution}
\item{e}{Step size used in the initial stage where independent normal proposal distributions are applied.}
\item{Sd}{Step size used in the adaptive stage to control the acceptance rate.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list object. Within each chain, the parameters are collected in list with labels of "tau0", "tau1", "phi", "w", respectively.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wenda Zhang
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# set tau0 and tau1 to generate data
tau0<-c(0,0,0,0,0,0,0,0,0,0)
tau1<-c(0,0,0,0,0,0,0.1,0.2,0.3,0.5)
data <- get_data(n=100, nb=5, tau0=tau0, tau1=tau1)
x <- data[[1]]
z <- data[[2]]
# Implement SPSL model
result <- spsl(x, z, iterations=5000, burnIn=1000, chains=1)
summary(result[[1]]$tau1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
