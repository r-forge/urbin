\name{uProbitEla}
\alias{uProbitEla}

\title{Semi-Elasticities from Probit Models}

\usage{
uProbitEla( allCoef, allXVal, 
  allCoefSE = rep( NA, length( allCoef )), xPos )
}

\description{
  Calculate Semi-Elasticities Based on Results from Probit Models
}

\arguments{
  \item{allCoef}{a vector with all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector with the values of all explanatory variables
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{allCoefSE}{an optional vector 
    with the values of the standard errors of all estimated coefficients
    (including intercept; the order of these standard errors must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xPos}{a posivite integer 
    indicating the position of the explanatory variable of interest
    or a vector of two positive integers 
    indicating the positions of the linear term and the quadratic term
    of the explanatory variable of interest
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefSE}.}
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

# mean values of the explanatory variables
xMeanQuad <- c( 1, mean( Mroz87$kids ), mean( Mroz87$age ), 
  mean( Mroz87$age )^2, mean( Mroz87$educ ) )

# semi-elasticity of age (including standard error)
uProbitEla( coef( estProbitQuad ), xMeanQuad, 
  sqrt( diag( vcov( estProbitQuad ) ) ), c( 3, 4 ) )
}

\keyword{models}
\keyword{regression}
