\name{urbinEffCat}
\alias{urbinEffCat}

\title{Grouping and Re-basing Categorical Variables}

\description{
  This function calculates 
  the effect of a categorical variable
  that is used as an explanatory variable
  in a linear probability model or a probit model,
  particularly if one wants to change the reference category and/or
  wants to group some of the categories together to a new category.
  The calculated effect indicates
  by how much the probability
  that the dependent variable has a value of one
  increases
  if the explanatory variable of interest
  changes from the 'reference' category 
  to a selected category of interest
  (this effect multiplied by 100
  indicates the increase in percentage points).
  This function applies the Delta-method
  to calculate the approximate standard error 
  of the calculated effect.
}

\usage{
urbinEffCat( allCoef, allXVal, xPos, xGroups, model, allCoefVcov = NULL )
}

\arguments{
  \item{allCoef}{a vector of all estimated coefficients
    (including intercept).}
  \item{allXVal}{a vector of the values of the explanatory variables,
    at which the effect should be calculated
    (including intercept; the order of its elements must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}),
    whereas the values that correspond to the categories 
    of the explanatory variable of interest must indicate 
    the shares of these categories in the sample of interest.}
  \item{xPos}{a vector of posivite integers 
    indicating the positions of the coefficients of the categorical variable
    in argument \code{allCoef}
    and, equally, the positions of the shares of the non-reference categories 
    in argument \code{allXVal};
    this vector must have one element less than the number of categories,
    because argument \code{allCoef} does not include a coefficient 
    of the reference category.}
  \item{xGroups}{a vector consisting of \eqn{-1}s, \eqn{0}s, and \eqn{1}s, 
    where a \eqn{-1} indicates 
    that the category should belong to the new reference category, 
    a \eqn{1} indicates 
    that the category should belong to the new category of interest,
    and a zero indicates 
    that the category belongs to neither;
    the length of this vector must be equal to the number of categories;
    the last element of this vector corresponds to the reference category,
    while all other elements correspond to the categories 
    as indicated by argument \code{xPos}.}
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
If \code{urbinEffCat()} calculates the effect of a variable of interest 
based on a linear probability model (\code{model = "lpm"}),
argument \code{allCoef} can only include the coefficient(s)
of the explanatory variable of interest.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
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

# effect of age changing from 38-52 years (2nd category + reference category) 
# to 53-60 years (3rd category) (with full covariance matrix)
urbinEffCat( coef( estProbitInt ), xMeanInt, c( 3:5 ), c( 0, -1, 1, -1 ),
  model = "probit", vcov( estProbitInt ) )

# effect of age changing from 38-52 years (2nd category + reference category) 
# to 53-60 years (3rd category) (with standard errors only)
urbinEffCat( coef( estProbitInt ), xMeanInt, c( 3:5 ), c( 0, -1, 1, -1 ),
  model = "probit", sqrt( diag( vcov( estProbitInt ) ) ) )
}

\keyword{models}
\keyword{regression}
