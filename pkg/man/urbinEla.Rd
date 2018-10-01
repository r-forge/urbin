\name{urbinEla}
\alias{urbinEla}

\title{Semi-Elasticities of Linear and Quadratic Explanatory Variables}

\description{
  This function calculates 
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
  This function applies the Delta-method
  to calculate the approximate standard error 
  of the calculated semi-elasticity.
}

\usage{
urbinEla( allCoef, allXVal, xPos, model, allCoefVcov = NULL, 
  allCoefBra = NULL, allXValBra = NULL, yCat = NULL, yCatBra = NULL, 
  lambda = NULL,
  seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL )
}

\arguments{
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
  \item{model}{a character string that indicates the type of model 
    that was estimated to obtain the parameter estimates:
    \code{"lpm"} = binary probit model (see \bold{Details} section below), 
    \code{"probit"} = binary probit model, 
    \code{"logit"} = binary logit model, 
    \code{"MNL"} = multinomial logit model,
    \code{"CondL"} = conditional logit model, or
    \code{"NestedL"} = nested logit model.}
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
  \item{allCoefBra,allXValBra,yCat,yCatBra,lambda}{
    will be documented later}
}

\details{
If \code{urbin} calculates the semi-elasticity 
based on a linear probability model (\code{model = "lpm"}),
argument \code{allCoef} can only include the coefficient(s)
of the explanatory variable of interest.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly,
whereas argument \code{allXVal} can be a scalar 
with the value of the variable of interest
even if the model includes both a linear and a quadratic term.
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
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit",
  vcov( estProbitQuad ) )

# semi-elasticity of age (only standard errors of coefficients)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ), seSimplify = FALSE )

# semi-elasticity of age (only standard errors of coefficients, 
# simplified equation for calculating standard error of the semi-elasticity)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ) )

# semi-elasticity of age (standard errors of coefficients, 
# mean and standard deviation of age provided)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ), seSimplify = FALSE,
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

# semi-elasticity of age (standard errors of coefficients,
# mean and standard deviation of age provided, 
# simplified equation for calculating standard error of the semi-elasticity)
urbinEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
}

\keyword{models}
\keyword{regression}
