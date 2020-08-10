\name{Germination.niche}
\alias{Germination.niche}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Germination.niche}

\description{The function performs linear model with quadratic term to identify the value of the ideal environmental variable for germination (germination niche) from a series of tests. Indeed perform ANOVA and Tukey test. The results are plotted using boxplots, different letters indicate significant differences among tests. The figure is automatically saved in 16:9 at 300 dpi.}

\usage{Germination.niche(germ_dyn, Nv.seed=NULL, n.seed=20,Env_Var=null)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{germ_dyn}{germination dynamic matrix, row=tests, col=days, frist column=type of test }
  \item{n.seed}{number of seeds for each petri}
  \item{Nv.seed}{Number of not empty seeds, must be a vector or the column number where these values are     located}
  \item{Env.var}{Environmental variable, the column number where these values are located}

%%     ~~Describe \code{x} here~~
}

\details{
  \item{M_env}{Final germination percentage each petri, with associated environmental varable}
  \item{Z_env}{Germination delay each petri, with associated environmental varable}
  \item{T50_env}{T50 each petri, with associated environmental varable}
  \item{Plot_M_env}{Smoothed fit for final germination percentage, with boxplot for each test}
  \item{Plot_Z_env}{Smoothed fit for germination delay, with boxplot for each test}
  \item{Plot_T50_env}{Smoothed fit for germination speed T50, with boxplot for each test}
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
file_germ<-read.csv("your_germination_file.csv", sep=";", header=T)

Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=10, Env.var=3)

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
