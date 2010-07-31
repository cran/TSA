\name{fitted.Arimax}
\alias{fitted.Arimax}
\title{Fitted values of an arimax model.}
\description{
Computes the fitted values of an arimax model.}
\usage{
\method{fitted}{Arimax}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model from the arimax function.}
  \item{...}{other arguments; not used here but kept to be consistent with the 
generic method}
}
\value{
fitted values}
\author{Kung-Sik Chan}
\seealso{ \code{\link{arimax}}}
\examples{
data(airmiles)
air.m1=arimax(log(airmiles),order=c(0,1,1),seasonal=list(order=c(0,1,1),
period=12),xtransf=data.frame(I911=1*(seq(airmiles)==69),
I911=1*(seq(airmiles)==69)),
transfer=list(c(0,0),c(1,0)),xreg=data.frame(Dec96=1*(seq(airmiles)==12),
Jan97=1*(seq(airmiles)==13),Dec02=1*(seq(airmiles)==84)),method='ML')
plot(log(airmiles),ylab="log(airmiles)")
points(fitted(air.m1))
}
\keyword{methods}
