\name{probitElaInt}
\alias{probitElaInt}

\title{Semi-Elasticities of Interval-Coded Variables from Probit Models}

\usage{
probitElaInt( allCoef, allXVal, xPos, xBound, allCoefVcov = NULL )
}

\description{
  This function calculates 
  the semi-elasticity of an interval-coded explanatory variable
  in a probit model.
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
  \item{allCoef}{a vector with all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector with the values of all explanatory variables
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{xPos}{a vector of posivite integers
    with length equal to the number of intervals.
    Each element of this vector refers to one interval,
    whereas these intervals need to be in ascending order.
    Each element indicates 
    the position of the coefficient of the respective dummy variable 
    (i.e., the dummy variable that indicates 
    whether the values of the explanatory variable of interest lie
    in the corresponding interval or not)
    in arguments \code{allCoef} 
    and the position of the mean value (=share) of this dummy variable 
    in argument \code{allXVal}.
    The position of the reference interval
    (i.e., the interval without a corresponding dummy variable
    as explanatory variable)
    must be indicated by a zero.}
  \item{xBound}{a numeric vector
    indicating the boundaries of the intervals.}
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

# estimation
estProbitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitInt )

# mean values of the explanatory variables
xMeanInt <- c( 1, colMeans( Mroz87[ , 
  c( "kids", "age30.37", "age38.44", "age53.60", "educ" ) ] ) )

# semi-elasticity of age (full covariance matrix)
probitElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  vcov( estProbitInt ) )

# semi-elasticity of age with standard errors (only standard errors)
probitElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  sqrt( diag( vcov( estProbitInt ) ) ) )
}

\keyword{models}
\keyword{regression}
