\name{probitEffGroup}
\alias{probitEffGroup}

\title{Probit Models: Grouping and Re-basing Categorical Variables}

\usage{
probitEffGroup( allCoef, allXVal, xPos, xGroups, allCoefVcov = NULL )
}

\description{
  This function calculates 
  the effect of a categorical variable in a probit model,
  particularly if one wants to change the reference category and/or
  wants to group some of the categories together to a new category.
  The calculated effect indicates
  by how many percentage points the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  changes from the 'reference' category 
  to a selected category of interest.
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
  \item{xPos}{a vector of posivite integers 
    indicating the positions of the coefficients of the categorical variable
    in argument \code{allCoef}
    and, equally, the positions of the shares of the categories 
    of the categorical variable in argument \code{allXVal}.}
  \item{xGroups}{a vector consisting of \eqn{-1}s, \eqn{0}s, and \eqn{1}s, 
    where a \eqn{-1} indicates 
    that the category belongs to the new reference category, 
    a \eqn{1} indicates 
    that the category belongs to the new category of interest,
    and a zero indicates 
    that the category belongs to neither.}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
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

# mean values of the explanatory variables (including the reference category) 
xMeanGroup <- c( 1, colMeans( Mroz87[ , 
  c( "kids", "age30.37", "age38.44", "age45.52", "age53.60", "educ" ) ] ) )

# coefficients including the coefficient of the reference category
allCoefGroup <- c( coef( estProbitInt )[1:4], 0, coef( estProbitInt )[5:6] )

# variance covariance matrix of the coefficients 
# including the coefficient of the reference category
allCoefVcovGroup <- vcov( estProbitInt )
allCoefVcovGroup <- rbind( 
  allCoefVcovGroup[ 1:4, ], rep( 0, 6 ), allCoefVcovGroup[ 5:6, ] )
allCoefVcovGroup <- cbind( 
  allCoefVcovGroup[ , 1:4 ], rep( 0, 7 ), allCoefVcovGroup[ , 5:6 ] )

# effects of age changing from 30-44 years to 53-60 years
# (with full covariance matrix)
probitEffGroup( allCoefGroup, xMeanGroup, c( 3:6 ), c( -1, -1, 0, 1 ),
  allCoefVcovGroup )

# effects of age changing from 30-44 years to 53-60 years
# (with standard errors only)
probitEffGroup( allCoefGroup, xMeanGroup, c( 3:6 ), c( -1, -1, 0, 1 ),
  sqrt( diag( allCoefVcovGroup ) ) )
}

\keyword{models}
\keyword{regression}
