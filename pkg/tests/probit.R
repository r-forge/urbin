library( "urbin" )
library( "maxLik" )
library( "mfx" )

# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

### linear in age
estProbitLin <- glm( lfp ~ kids + age + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitLin )
# mean values of the explanatory variables
xMeanLin <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# semi-elasticity of age without standard errors
probitEla( coef( estProbitLin ), xMeanLin, xPos = 3 )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estProbitLin, 
  newdata = as.data.frame( t( xMeanLin * c( 1, 1, 1.005, 1 ) ) ), 
  type = "response" ) -
    predict( estProbitLin, 
      newdata = as.data.frame( t( xMeanLin * c( 1, 1, 0.995, 1 ) ) ), 
      type = "response" ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitElaDeriv( coef( estProbitLin ), xMeanLin, 3 )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( probitEla, t0 = coef( estProbitLin ), 
  allXVal = xMeanLin, xPos = 3 )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitElaDeriv( coef( estProbitLin ), xMeanLin, 3,
  simplified = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
probitEla( coef( estProbitLin ), xMeanLin, 3, vcov( estProbitLin ) )
# semi-elasticity of age with standard errors (only standard errors)
probitEla( coef( estProbitLin ), xMeanLin, 3,
  sqrt( diag( vcov( estProbitLin ) ) ), seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
probitEla( coef( estProbitLin ), xMeanLin, 3, 
  sqrt( diag( vcov( estProbitLin ) ) ) )
# semi-elasticity of age based on partial derivative calculated by the mfx package
estProbitLinMfx <- probitmfx( lfp ~ kids + age + educ, data = Mroz87 )
estProbitLinMfx$mfxest[ "age", 1:2 ] * xMeanLin[ "age" ]

### quadratic in age
estProbitQuad <- glm( lfp ~ kids + age + I(age^2) + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitQuad )
# mean values of the explanatory variables
xMeanQuad <- c( xMeanLin[ 1:3], xMeanLin[3]^2, xMeanLin[4] )
# semi-elasticity of age without standard errors
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ) )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estProbitQuad, 
  newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 1.005, 1.005^2, 1 ) ) ), 
  type = "response" ) -
    predict( estProbitQuad, 
      newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 0.995, 0.995^2, 1 ) ) ), 
      type = "response" ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitElaDeriv( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ) )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( probitEla, t0 = coef( estProbitQuad ), 
  allXVal = xMeanQuad, xPos = c( 3, 4 ) )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitElaDeriv( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ),
  simplified = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), vcov( estProbitQuad ) )
# semi-elasticity of age with standard errors (only standard errors)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ), seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ) )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on the original data
se <- sqrt( diag( vcov( estProbitQuad ) ) )
X <- cbind( Mroz87$age, Mroz87$age^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), vcovApp )
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), vcovApp,
  seSimplify = TRUE )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on simulated data
se <- sqrt( diag( vcov( estProbitQuad ) ) )
set.seed( 123 )
x <- rnorm( 1000, xMeanQuad[ "age" ], sd( Mroz87$age ) )
X <- cbind( x, x^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), vcovApp )
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), vcovApp,
  seSimplify = TRUE )
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ), 
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
probitEla( coef( estProbitQuad ), xMeanQuad, c( 3, 4 ), 
  sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
# semi-elasticity of age based on partial derivatives calculated by the mfx package
# (differs from the above, because mean(age)^2 is not the same as mean(age^2))
estProbitQuadMfx <- probitmfx( lfp ~ kids + age + I(age^2) + educ, data = Mroz87 )
estProbitQuadMfx$mfxest[ "age", 1:2 ] * xMeanQuad[ "age" ] +
  2 * estProbitQuadMfx$mfxest[ "I(age^2)", 1:2 ] * xMeanQuad[ "age" ]^2

### age is interval-coded (age is in the range 30-60)
# create dummy variables for age intervals
Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60
all.equal( 
  Mroz87$age30.37 + Mroz87$age38.44 + Mroz87$age45.52 + Mroz87$age53.60,
  rep( 1, nrow( Mroz87 ) ) )
# estimation
estProbitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitInt )
# mean values of the explanatory variables
xMeanInt <- c( xMeanLin[1:2], mean( Mroz87$age30.37 ), 
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMeanLin[4] )
# semi-elasticity of age without standard errors
probitElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitElaIntDeriv( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ) )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( probitElaInt, t0 = coef( estProbitInt ), allXVal = xMeanInt, 
  xPos = c( 3, 4, 0, 5 ), xBound = c( 30, 37.5, 44.5, 52.5, 60 ) )
# semi-elasticity of age with standard errors (full covariance matrix)
probitElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  vcov( estProbitInt ) )
# semi-elasticity of age with standard errors (only standard errors)
probitElaInt( coef( estProbitInt ), xMeanInt, 
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), 
  sqrt( diag( vcov( estProbitInt ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear explanatory variable 
# mean values of the 'other' explanatory variables
xMeanLinInt <- c( xMeanLin[ 1:2 ], NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
probitEffInt( coef( estProbitLin ), xMeanLinInt, 3,
  c( 30, 40 ), c( 50, 60 ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitEffIntDeriv( coef( estProbitLin ), xMeanLinInt, 3, 
  c( 30, 40 ), c( 50, 60 ) )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( probitEffInt, t0 = coef( estProbitLin ), 
  allXVal = xMeanLinInt, xPos = 3, 
  refBound = c( 30, 40 ), intBound = c( 50, 60 ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
probitEffInt( coef( estProbitLin ), xMeanLinInt, 3,
  c( 30, 40 ), c( 50, 60 ), vcov( estProbitLin ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
probitEffInt( coef( estProbitLin ), xMeanLinInt, 3,
  c( 30, 40 ), c( 50, 60 ), sqrt( diag( vcov( estProbitLin ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear and quadratic explanatory variable 
# mean values of the 'other' explanatory variables
xMeanQuadInt <- c( xMeanLin[ 1:2 ], NA, NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
probitEffInt( coef( estProbitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::probitEffIntDeriv( coef( estProbitQuad ), xMeanQuadInt, c( 3, 4 ), 
  c( 30, 40 ), c( 50, 60 ) )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( probitEffInt, t0 = coef( estProbitQuad ), 
  allXVal = xMeanQuadInt, xPos = c( 3, 4 ), 
  refBound = c( 30, 40 ), intBound = c( 50, 60 ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
probitEffInt( coef( estProbitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), vcov( estProbitQuad ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
probitEffInt( coef( estProbitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), sqrt( diag( vcov( estProbitQuad ) ) ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (standard errors + mean value and standard deviation of age)
probitEffInt( coef( estProbitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), sqrt( diag( vcov( estProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
