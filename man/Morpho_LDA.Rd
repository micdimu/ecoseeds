\name{Morpho_LDA}
\alias{Morpho_LDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Morpho_LDA
%%  ~~function to do ... ~~
}
\description{The function perform Linear Disciminant Analysis with morphometric data to discriminate taxa (at least 3 groups). The results of LDA are the permormance of discrimination with test-training data, and plotted both with 2 and 3 dimension. The 2 dimension figure is automatically saved in 16:9 at 300 dpi.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Morpho_LDA(Morpho_df, N.iter=100, split_size=0.80)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Morpho_df}{Morphometric matrix, row=replicates, col=morphometric variable, frist column=names of taxa (e.g. Genus, species, population...)}
  \item{N.iter}{Number of iterations to assest the LDA performance}
  \item{split_size}{Number from 0 to 1 to split training test data}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{morpho.lda}{output from LDA function}
  \item{accuracy_train}{accurancy discrimination with train data}
  \item{accuracy_test}{accurancy discrimination with test data}
  \item{confusion_train}{confusion matrix with train data}
  \item{confusion_test}{confusion matrix with test data}
  \item{plot_lda}{ggplot output}
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
Morpho_data<-read.csv("your_morpho_file.csv", sep=";", header=T)

Morpho_LDA_output<-Morpho_LDA(Morpho_data, N.iter=100, split_size=0.80)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
