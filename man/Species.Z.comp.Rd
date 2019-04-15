\name{Species.Z.comp}
\alias{Species.Z.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Species.Z.comp
%%  ~~function to do ... ~~
}
\description{Compare germination delay percentage for two species. The results are plotted using boxplots, different letters indicate significant differences among testS and species. The figure is automatically saved in 16:9 at 300 dpi.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Species.Z.comp(Germ.Analysis.exp_sp1, Germ.Analysis.exp_sp2, sp_name=NULL , colour="yes", Test.int=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Germination.Analysis.output_sp1}{output from Germination.Analysis function for species 1}
  \item{Germination.Analysis.output_sp2}{output from Germination.Analysis function for species 2}
  \item{colour}{can be "yes" (coloured by test.type) or "no" (B/W output) or a vector that specify the groups of tests (length of vector must be equal to n° of petri)}
  \item{Test.int}{character vector where are indicated the types of tests that would be compared}
%%     ~~Describe \code{x} here~~

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Anova }{anova results, summary(aov)}
  \item{Tukey }{Post-hoc Tukey test output}
  \item{Test }{Test compared}
  \item{Germination_delay}{Germination delay percentage for each test analysed}
  \item{Boxplot_Z }{ggplot output}
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
file_germ_sp1<-read.csv("your_germination_file_sp1.csv", sep=";", header=T)
file_germ_sp2<-read.csv("your_germination_file_sp2.csv", sep=";", header=T)

Germination.Analysis.output_sp1<-Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=20, cv=1.5)
Germination.Analysis.output_sp2<-Germination.Analysis(file_germ, Nv.seed = NULL, n.seed=20, cv=1.5)

Species.Z.comp(Germination.Analysis.output_sp1, Germination.Analysis.output_sp2, colour = "yes")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
