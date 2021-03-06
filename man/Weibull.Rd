\name{Weibull}
\alias{Weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weibull}

\description{The function performs a Weibull dynamic for each test. The figure is automatically saved in 16:9 at 300 dpi.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Weibull(Germination.Analysis.output, Test.int)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Germination.Analysis.output}{output from Germination.Analysis function}
  %\item{colour}{can be "yes" (coloured by test.type) or "no" (BAW output) or a vector that specify the groups of test (length of vector must be equal to n° of petri)}
  \item{Test.int}{character vector where are indicated the types of test that would be compared}
%%     ~~Describe \code{x} here~~

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Weibull}{germination percentage from weibull function for each days}
  \item{Test }{Tests compared}
  \item{plot_Weibull}{ggplot output}

%% ...
}
\references{

Weibull, W. (1951). A statistical distribution function of wide applicability. Journal of applied mechanics, 18(3), 293-297.

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

Germination.Analysis.output<-Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=20, cv=1.5)

Weibull(Germination.Analysis.output)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line\name{Weibull}
