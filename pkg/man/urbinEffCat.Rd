\name{urbinEffCat}
\alias{urbinEffCat}

\title{Grouping and Re-basing Categorical Variables}

\description{
  This function calculates 
  the effect of a categorical variable
  that is used as an explanatory variable
  in a linear probability model, a binary probit model,
  an ordered probit model, or a binary logit model,
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
urbinEffCat( allCoef, allXVal, xPos, xGroups, model, 
  allCoefVcov = NULL, iPos = ifelse( all( xPos != 1 ), 1, 0 ), 
  yCat = NULL )
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
    \code{"probit"} = binary probit model,
    \code{"oprobit"} = ordered probit model (see \bold{Details} section below), 
    \code{"logit"} = binary logit model, 
    \code{"MNL"} = multinomial logit model, or
    \code{"CondL"} = conditional logit model.}
  \item{allCoefVcov}{an optional argument
    that can be the variance-covariance matrix of all estimated coefficients
    OR the vector of the standard errors of all estimated coefficients
    (including intercept; the order of these values must be the same
    as the order of the corresponding coefficients in argument 
    \code{allCoef}).}
  \item{iPos}{a posivite integer 
    indicating the position of the intercept
    in arguments \code{allCoef}, \code{allXVal}, 
    and eventually \code{allCoefVcov}.
    A value of zero indicates that the model does not have an intercept.}
  \item{yCat}{not documented yet.}
}

\details{
Linear Probability Model:
If a user wants to calculate the effect of a variable of interest
based on a linear probability model 
using \code{urbinEffCat( ..., model = "lpm")},
he/she can include only the coefficient(s)
of the explanatory variable of interest in argument \code{allCoef}
and omit all other coefficients.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly.
This simplified use of \code{urbinEffCat()}
is only unsuitable for linear probability models
and unsuitable for all other types of models.

Ordered Probit Model:
The dependent variable has \eqn{P>2} distinct and strictly ordered categories.
A user needs to choose a number \eqn{p \in \{1,\ldots, P-1\}}
that separates the \eqn{P} ordered categories 
into two mutually exclusive combined categories,
where the lower \eqn{p} categories indicate a binary variable with value zero
and the upper \eqn{P-p} categories indicate a binary variable with value one.
An ordered probit model does not have an intercept
but it has \eqn{P-1} thresholds between the \eqn{P} ordered categories.
If a user wants to calculate the semi-elasticity 
based on an ordered probit model 
using \code{urbinEla( ..., model = "oprobit")},
argument \code{allCoef} must include 
both the coefficients and the thresholds.
Consequently, argument \code{allCoefVcov} 
must also include the variences and covariances
or the standard arrors of both the coefficients and the thresholds.
The vector specified by argument \code{allXVal} 
must include elements that correspond to the thresholds,
where the element that corresponds to the \eqn{p}'s threshold
must be equal to \eqn{-1}, 
while the elements that correspond to the other thresholds
must be equal to \eqn{0}.
Argument \code{iPos} must specify the position of the \eqn{p}'s threshold
in arguments \code{allCoef}, \code{allXVal}, 
and eventually \code{allCoefVcov}.
While the \eqn{p}'s threshold must be specified in argument \code{allCoef}, 
all other thresholds can be omitted.
In this case, arguments \code{allXVal} and \code{allCoefVcov}
must be adjusted accordingly.
}

\value{
\code{urbinElaInt()} returns a list of class \code{"urbin"} 
that contains the following components:
\item{effect}{the calculated effect.}
\item{stdEr}{the approximate standard error of the calculated effect.}
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
