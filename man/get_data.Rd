\name{get_data}
\alias{get_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Get simulation data of dynamic gene coexpressions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get_data is used to generate the simulation data for dynamic gene coexpressions. It generates a matrix of gene expression data x and a vector of uniform independent varialbe. n is the number of rows of x as well as the observations in z. nb is the number of columns of x.
}
\usage{
get_data(n, nb, tau0, tau1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{The number of rows in x as well as the number of observations in z}
\item{nb}{The number of columns in x}
\item{tau0}{The fisher's Z-transformed correlation coefficient between every two genes when z is 0}
\item{tau1}{The association between z and fisher's Z-transformed correlation coefficient between every two genes}

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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wenda Zhang
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
