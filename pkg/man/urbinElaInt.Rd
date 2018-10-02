\name{urbinElaInt}
\alias{urbinElaInt}

\title{Semi-Elasticities of Interval-Coded Explanatory Variables}

\description{
  This function calculates 
  the semi-elasticity of an interval-coded variable
  that is used as an explanatory variable
  in a linear probability model or a probit model.
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
urbinElaInt( allCoef, allXVal, xPos, xBound, model, allCoefVcov = NULL )
}

\arguments{
  \item{allCoef}{a vector of all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector of the values of all explanatory variables
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}),
    whereas the values that correspond to the intervals 
    of the explanatory variable of interested must indicate 
    the shares of these intervals in the sample of interest.}
  \item{xPos}{a vector of positive integers
    with length equal to the number of intervals.
    Each element of this vector must refer to one interval,
    whereas these intervals must be in ascending order.
    Each element of this vector indicates 
    the position of the coefficient of the respective dummy variable 
    (i.e., the dummy variable that indicates 
    whether the values of the explanatory variable of interest lie
    in the corresponding interval or not)
    in arguments \code{allCoef} 
    and the position of the value (=share) of this dummy variable 
    in argument \code{allXVal}.
    The position of the reference interval
    (i.e., the interval without a corresponding dummy variable
    as explanatory variable)
    must be indicated by a zero.}
  \item{xBound}{a numeric vector
    indicating the boundaries of the intervals.}
  \item{model}{a character string that indicates the type of model 
    that was estimated to obtain the parameter estimates:
    \code{"lpm"} = binary probit model (see \bold{Details} section below), 
    \code{"probit"} = binary probit model.}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
}

\details{
If \code{urbinElaInt()} calculates the semi-elasticity 
based on a linear probability model (\code{model = "lpm"}),
argument \code{allCoef} can only include the coefficient(s)
of the explanatory variable of interest.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly.
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

# create dummy variables for age intervals
Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60

# probit estimation
estProbitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitInt )

# mean values of the explanatory variables
xMeanInt <- c( 1, colMeans( Mroz87[ , 
  c( "kids", "age30.37", "age38.44", "age53.60", "educ" ) ] ) )

# semi-elasticity of age (full covariance matrix of coefficients)
urbinElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  model = "probit", vcov( estProbitInt ) )

# semi-elasticity of age (only standard errors of coefficients)
urbinElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  model = "probit", sqrt( diag( vcov( estProbitInt ) ) ) )
}

\keyword{models}
\keyword{regression}
