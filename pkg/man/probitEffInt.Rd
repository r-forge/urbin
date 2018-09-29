\name{probitEffInt}
\alias{probitEffInt}

\title{Probit Models: Effects of Changes between Intervals}

\usage{
probitEffInt( allCoef, allXVal, xPos, refBound, intBound, 
  allCoefVcov = NULL, xMeanSd = NULL )
}

\description{
  This function calculates 
  the effects of continuous variables in a probit model,
  when these variables change between discrete intervals.
  The calculated effect indicates
  by how many percentage points the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  changes from the the 'reference' interval 
  to a selected interval of interest.
  This function also applies the Delta-method
  to calculate the approximate standard error 
  of the calculated effect.
}

\arguments{
  \item{allCoef}{a vector with all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector with the values of all explanatory variables
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xPos}{a posivite integer 
    indicating the position of the explanatory variable of interest
    or a vector of two positive integers 
    indicating the positions of the linear term and the quadratic term
    of the explanatory variable of interest
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefSE}.}
  \item{refBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the 'reference' interval, respectively.}
  \item{intBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the interval of interest, respectively.}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xMeanSd}{a vector with two elements: 
    the mean value and the standard deviation of the explanatory variable 
    of interest.
    If argument \code{allCoefVcov} is a vector of standard errors
    and the explanatory variables include a linear and a quadratic term
    of the variable of interest,
    the information provided in argument \code{xMeanSd} 
    is used to approximate the covariance 
    between the coefficient of the linear term 
    and the coefficient of the quadratic term.}
}

\value{
  A numeric vector with two elements: 
  the semi-elasticity and its approximate standard error.
}

\examples{
# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

# estimate probit model with linear and quadratic age terms
estProbitQuad <- glm( lfp ~ kids + age + I(age^2) + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitQuad )

# mean values of the 'other' explanatory variables
xMeanQuad <- c( 1, mean( Mroz87$kids ), NA, NA, mean( Mroz87$educ ) )

# effect of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors of coefficients)
probitEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), sqrt( diag( vcov( estProbitQuad ) ) ) )
  
# effect of age changing from the 30-40 interval to the 50-60 interval
# (with standard errors of coefficients as well as 
# mean and standard deviation of age) 
probitEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
}

\keyword{models}
\keyword{regression}
