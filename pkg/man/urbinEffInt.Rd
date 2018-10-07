\name{urbinEffInt}
\alias{urbinEffInt}

\title{Effects of Changes between Intervals}

\description{
  This function calculates 
  the effect of a continuous variable
  that is used as a linear explanatory variable
  or as a linear and quadratic explanatory variable
  in a linear probability model, a probit model, or a logit model
  if this variable changes between discrete intervals.
  The calculated effect indicates
  by how much the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  changes from the 'reference' interval 
  to a selected interval of interest
  (this effect multiplied by 100
  indicates the increase in percentage points).
  This function applies the Delta-method
  to calculate the approximate standard error 
  of the calculated effect.
}

\usage{
urbinEffInt( allCoef, allXVal = NULL, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL, yCat = NULL,
  allCoefBra = NULL, allXValBra = NULL, yCatBra = NULL, lambda = NULL )
}

\arguments{
  \item{allCoef}{a vector of all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector of the values of the explanatory variables,
    at which the effect should be calculated
    (including intercept; the order of its elements must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).
    In case of linear probability models,
    this argument is ignored.
    In case of all other types of models, 
    the element(s) that correspond to the variable of interest
    is (are) ignored and can be set to \code{NA}.}
  \item{xPos}{a posivite integer 
    indicating the position of the explanatory variable of interest
    or a vector of two positive integers 
    indicating the positions of the linear term and the quadratic term
    of the explanatory variable of interest
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefVcov}.}
  \item{refBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the 'reference' interval, respectively.}
  \item{intBound}{a numeric vector of two elements
    that indicate the lower boundary and the upper boundary
    of the interval of interest, respectively.}
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
  \item{yCat,allCoefBra,allXValBra,yCatBra,lambda}{not yet documented.}
}

\details{
If \code{urbinEffInt()} calculates the effect of a variable of interest 
based on a linear probability model (\code{model = "lpm"}),
argument \code{allCoef} can only include the coefficient(s)
of the explanatory variable of interest.
In this case, argument \code{allCoefVcov}
must be adjusted accordingly.
}

\value{
  A numeric vector with two elements: 
  the effect and its approximate standard error.
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
urbinEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit",
  sqrt( diag( vcov( estProbitQuad ) ) ) )
  
# effect of age changing from the 30-40 interval to the 50-60 interval
# (with standard errors of coefficients as well as 
# mean and standard deviation of age) 
urbinEffInt( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
}

\keyword{models}
\keyword{regression}
