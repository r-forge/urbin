\name{probitEla}
\alias{lpmEla}
\alias{probitEla}
\alias{logitEla}

\title{Semi-Elasticities of Linear and Quadratic Explanatory Variables}

\usage{
lpmEla( xCoef, xVal, xCoefVcov = NULL, xMeanSd = NULL )
probitEla( allCoef, allXVal, xPos, allCoefVcov = NULL,
  seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL )
logitEla( allCoef, allXVal, xPos, allCoefVcov = NULL, 
  allCoefBra = NULL, allXValBra = NULL, yCat = NULL, yCatBra = NULL, 
  lambda = NULL, method =  "binary",
  seSimplify = !is.matrix( allCoefVcov ) )
}

\description{
  These functions calculate 
  the semi-elasticity of a continuous explanatory variable
  that is used as a linear explanatory variable
  or as a linear and quadratic explanatory variable
  in a linear probability model, a probit model, or a logit model.
  The calculated semi-elasticity indicates
  by how many percentage points the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  increases by one percent.
  This function also applies the Delta-method
  to calculate the approximate standard error 
  of the calculated semi-elasticity.
}

\arguments{
  \item{xCoef}{a numeric scalar that indicates 
    the estimated coefficient of the explanatory variable of interest
    OR a numeric vector with two elements that indicates
    the estimated coefficients of the linear term and the quadratic term
    of the explanatory variable of interest.}
  \item{xVal}{a numeric scalar that indicates 
    the value of the explanatory variable of interest,
    at which the semi-elasticities should be calculated.}
  \item{xCoefVcov}{an optional argument
    that can be a numeric scalar equal to the standard error 
    of the coefficient of the explanatory variable of interest 
    (if there is only a linear term of the explanatory variable of interest)
    OR a numeric vector with two elements equal to the standard errors
    of the coefficients of the linear and quadratic terms
    of the explanatory variable of interest
    OR a numeric 2x2 matrix with the elements of the variance covariance
    matrix of the estimated coefficients
    that correspond to the linear and quadratic terms
    of the explanatory variable of interest.}
  \item{allCoef}{a vector with all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector with the values of the explanatory variables,
    at which the semi-elasticity should be calculated
    (including intercept; the order of its elements must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xPos}{a posivite integer 
    indicating the position of the explanatory variable of interest
    or a vector of two positive integers 
    indicating the positions of the linear term and the quadratic term
    of the explanatory variable of interest
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefVcov}.}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{seSimplify}{logical value indicating whether the standard errors
    of the semi-elasticities should be calculated by a simplified equation
    (i.e., assuming that \eqn{\phi( \beta' x )} does not depend on \eqn{\beta}),
    which may be more appropriate 
    if the off-diagonal elements of the variance-covariance matrix are unknown.}
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
  \item{allCoefBra,allXValBra,yCat,yCatBra,lambda,method}{
    will be documented later}
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

# semi-elasticity of age (full covariance matrix of coefficients)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  vcov( estProbitQuad ) )

# semi-elasticity of age (only standard errors of coefficients)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ), seSimplify = FALSE )

# semi-elasticity of age (only standard errors of coefficients, 
# simplified equation for calculating standard error of the semi-elasticity)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ) )

# semi-elasticity of age (standard errors of coefficients, 
# mean and standard deviation of age provided)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ), seSimplify = FALSE,
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

# semi-elasticity of age (standard errors of coefficients,
# mean and standard deviation of age provided, 
# simplified equation for calculating standard error of the semi-elasticity)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
}

\keyword{models}
\keyword{regression}
