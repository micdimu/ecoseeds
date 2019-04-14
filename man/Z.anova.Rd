\name{Z.anova}
\alias{Z.anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Z.anova}

\description{The function performs ANOVA and Tukey test. The results are plotted using boxplots, different letters indicate significant differences among tests. The figure is automatically saved in 16:9 at 300 dpi.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Z.anova(Germination.Analysis.output, colour = "yes", Test.int)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Germination.Analysis.output}{output from Germination.Analysis function}
  \item{colour}{can be "yes" (coloured by test.type) or "no" (BAW output) or a vector that specify the groups of testS (length of vector must be equal to n° of petri)}
  \item{Test.int}{character vector where are indicated the types of testS that would be compared}
%%     ~~Describe \code{x} here~~

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Anova }{anova results, summary(aov)}
  \item{Tukey }{Post-hoc Tukey test output}
  \item{Test }{TestS compared}
  \item{Z }{Z for each test analysed}
  \item{Boxplot_Z }{ggplot output}

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

Germination.Analysis.output<-Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=20, cv=1.5)

Z.anova(Germination.Analysis.output, colour = "yes")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line\name{Z.anova}
