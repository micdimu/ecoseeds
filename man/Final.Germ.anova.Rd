\name{Final.Germ.anova}
\alias{Final.Germ.anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Final.Germ.anova}

\description{The function performs anova and Tukey test than plot results using barplots, different letters indicate significant differences among test. The figure is automatically saved in 16:9 at 300 dpi.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Final.Germ.anova(Germination.Analysis.output, colour = "yes", Test.int)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Germination.Analysis.output}{}
  \item{colour}{can be "yes" (coloured by test.type) or "no" (BAW output) or a vector that specify the groups of test (length of vector must be equal to n° of petri)}
  \item{Test.int}{character vector where are indicated the types of test that would be compared}
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
  \item{Test }{Test compared}
  \item{Fin_germ }{Final germination percentage for each test analysed}
  \item{Boxplot_Fin.Germ }{ggplot output}

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

Final.Germ.anova(Germination.Analysis.output, colour = "yes")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
