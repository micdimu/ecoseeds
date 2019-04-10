\name{Germination.Analysis}
\alias{Germination.Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Germination.Analysis
%%  ~~function to do ... ~~
}
\description{The function analyses germination tests automatically
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Germination.Analysis(germ_dyn, Nv.seed=NULL, cv=1.5, n.seed=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{germ_dyn}{germination dynamic matrix, row=tests, col=days, frist column=type of test }
  \item{n.seed}{number of seeds for each petri}
  \item{Nv.seed}{Number of not empty seeds, must be a vector or the column number where this values are located}
  \item{cv}{Weibull curve shape parameter}

%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{test_type}{type of test}
  \item{Dyn_perc_germ_SP}{Germination percentage for each day - each petri}
  \item{M_SP}{Final germination percentage each petri}
  \item{Z_SP}{Germination delay each petri}
  \item{T50_SP}{T50 each petri}
  \item{Dyn_perc_germ}{Germination percentage for each day - each test}
  \item{M}{Final germination percentage each test}
  \item{M_SD}{Standard deviation of final germination percentage each test}
  \item{Z}{Germination delay each test}
  \item{Z_SD}{Standard deviation of Germination delay for each test}
  \item{T50}{T50 each test}
  \item{T50_SD}{Standard deviation of T50 for each test}
  \item{k}{Relative germination rate}
  \item{y}{germination percentage from weibull function for each days)}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Michele Di Musciano (michele.dimusciano@graduate.univaq.it)
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
## load file and run the function

file_germ<-read.csv("your_germination_file.csv", sep=";", header=T)

Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=10, cv=1.5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
