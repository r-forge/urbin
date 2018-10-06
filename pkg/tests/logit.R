library( "urbin" )
library( "maxLik" )
library( "mfx" )

# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

### linear in age
estLogitLin <- glm( lfp ~ kids + age + educ, 
  family = binomial(link = "logit"), 
  data = Mroz87 )
summary( estLogitLin )
# mean values of the explanatory variables
xMeanLin <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# semi-elasticity of age without standard errors
urbinEla( coef( estLogitLin ), xMeanLin, xPos = 3, model = "logit" )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estLogitLin, 
  newdata = as.data.frame( t( xMeanLin * c( 1, 1, 1.005, 1 ) ) ), 
  type = "response" ) -
    predict( estLogitLin, 
      newdata = as.data.frame( t( xMeanLin * c( 1, 1, 0.995, 1 ) ) ), 
      type = "response" ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinAttr <- xMeanLin
attr( xMeanLinAttr, "derivOnly" ) <- 1 
urbinEla( coef( estLogitLin ), xMeanLinAttr, 3, seSimplify = FALSE, 
  model = "logit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( estLogitLin ), 
  allXVal = xMeanLin, xPos = 3, model = "logit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coef( estLogitLin ), xMeanLinAttr, 3, model = "logit", 
  seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coef( estLogitLin ), xMeanLin, 3, model = "logit", 
  vcov( estLogitLin ) )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( estLogitLin ), xMeanLin, 3, model = "logit",
  sqrt( diag( vcov( estLogitLin ) ) ), seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coef( estLogitLin ), xMeanLin, 3, model = "logit", 
  sqrt( diag( vcov( estLogitLin ) ) ) )
# semi-elasticity of age based on partial derivative calculated by the mfx package
estLogitLinMfx <- logitmfx( lfp ~ kids + age + educ, data = Mroz87 )
estLogitLinMfx$mfxest[ "age", 1:2 ] * xMeanLin[ "age" ]

### quadratic in age
estLogitQuad <- glm( lfp ~ kids + age + I(age^2) + educ, 
  family = binomial(link = "logit"), 
  data = Mroz87 )
summary( estLogitQuad )
# mean values of the explanatory variables
xMeanQuad <- c( xMeanLin[ 1:3 ], xMeanLin[3]^2, xMeanLin[4] )
# semi-elasticity of age without standard errors
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit" )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estLogitQuad, 
  newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 1.005, 1.005^2, 1 ) ) ), 
  type = "response" ) -
    predict( estLogitQuad, 
      newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 0.995, 0.995^2, 1 ) ) ), 
      type = "response" ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanQuadAttr <- xMeanQuad
attr( xMeanQuadAttr, "derivOnly" ) <- 1 
urbinEla( coef( estLogitQuad ), xMeanQuadAttr, c( 3, 4 ), model = "logit",
  seSimplify = FALSE )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( estLogitQuad ), 
  allXVal = xMeanQuad, xPos = c( 3, 4 ), model = "logit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coef( estLogitQuad ), xMeanQuadAttr, c( 3, 4 ), model = "logit",
  seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  vcov( estLogitQuad ) )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  sqrt( diag( vcov( estLogitQuad ) ) ), seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit",
  sqrt( diag( vcov( estLogitQuad ) ) ) )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on the original data
se <- sqrt( diag( vcov( estLogitQuad ) ) )
X <- cbind( Mroz87$age, Mroz87$age^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  vcovApp )
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  vcovApp, seSimplify = TRUE )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on simulated data
se <- sqrt( diag( vcov( estLogitQuad ) ) )
set.seed( 123 )
x <- rnorm( 1000, xMeanQuad[ "age" ], sd( Mroz87$age ) )
X <- cbind( x, x^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  vcovApp )
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit", 
  vcovApp, seSimplify = TRUE )
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit",
  sqrt( diag( vcov( estLogitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
urbinEla( coef( estLogitQuad ), xMeanQuad, c( 3, 4 ), model = "logit",
  sqrt( diag( vcov( estLogitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
# semi-elasticity of age based on partial derivatives calculated by the mfx package
# (differs from the above, because mean(age)^2 is not the same as mean(age^2))
estLogitQuadMfx <- logitmfx( lfp ~ kids + age + I(age^2) + educ, data = Mroz87 )
estLogitQuadMfx$mfxest[ "age", 1:2 ] * xMeanQuad[ "age" ] +
  2 * estLogitQuadMfx$mfxest[ "I(age^2)", 1:2 ] * xMeanQuad[ "age" ]^2

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
estLogitInt <- glm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  family = binomial(link = "logit"), 
  data = Mroz87 )
summary( estLogitInt )
# mean values of the explanatory variables
xMeanInt <- c( xMeanLin[1:2], mean( Mroz87$age30.37 ), 
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMeanLin[4] )
# semi-elasticity of age without standard errors
urbinElaInt( coef( estLogitInt ), xMeanInt,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit" )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanIntAttr <- xMeanInt
attr( xMeanIntAttr, "derivOnly" ) <- 1 
urbinElaInt( coef( estLogitInt ), xMeanIntAttr,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinElaInt, t0 = coef( estLogitInt ), allXVal = xMeanInt,
  xPos = c( 3, 4, 0, 5 ), xBound = c( 30, 37.5, 44.5, 52.5, 60 ),
  model = "logit" )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinElaInt( coef( estLogitInt ), xMeanInt,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit",
  allCoefVcov = vcov( estLogitInt ) )
# semi-elasticity of age with standard errors (only standard errors)
urbinElaInt( coef( estLogitInt ), xMeanInt,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "logit",
  allCoefVcov = sqrt( diag( vcov( estLogitInt ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear explanatory variable 
# mean values of the 'other' explanatory variables
xMeanLinInt <- c( xMeanLin[ 1:2 ], NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( estLogitLin ), allXVal = xMeanLinInt, xPos = 3,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "logit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
predict( estLogitLin, 
  newdata = as.data.frame( t( replace( xMeanLin, 3, 55 ) ) ), 
  type = "response" ) -
  predict( estLogitLin, 
    newdata = as.data.frame( t( replace( xMeanLin, 3, 35 ) ) ), 
    type = "response" )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinIntAttr <- xMeanLinInt
attr( xMeanLinIntAttr, "derivOnly" ) <- 1 
urbinEffInt( coef( estLogitLin ), xMeanLinIntAttr, 3,
  c( 30, 40 ), c( 50, 60 ), model = "logit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( estLogitLin ),
  allXVal = xMeanLinInt, xPos = 3,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "logit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( estLogitLin ), xMeanLinInt, 3,
  c( 30, 40 ), c( 50, 60 ), model = "logit", allCoefVcov = vcov( estLogitLin ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( estLogitLin ), allXVal = xMeanLinInt, xPos = 3,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "logit", 
  allCoefVcov = sqrt( diag( vcov( estLogitLin ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear and quadratic explanatory variable 
# mean values of the 'other' explanatory variables
xMeanQuadInt <- c( xMeanLin[ 1:2 ], NA, NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( estLogitQuad ), allXVal = xMeanQuadInt, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "logit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
predict( estLogitQuad, 
  newdata = as.data.frame( t( replace( xMeanQuad, 3:4, c( 55, 55^2 ) ) ) ), 
  type = "response" ) -
  predict( estLogitQuad, 
    newdata = as.data.frame( t( replace( xMeanQuad, 3:4, c( 35, 35^2 ) ) ) ), 
    type = "response" )
# partial derivatives of the effect wrt the coefficients
xMeanQuadIntAttr <- xMeanQuadInt
attr( xMeanQuadIntAttr, "derivOnly" ) <- 1 
urbinEffInt( coef( estLogitQuad ), xMeanQuadIntAttr, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "logit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( estLogitQuad ),
  allXVal = xMeanQuadInt, xPos = c( 3, 4 ),
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "logit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( estLogitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "logit", 
  allCoefVcov = vcov( estLogitQuad ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( estLogitQuad ), allXVal = xMeanQuadInt, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "logit", allCoefVcov = sqrt( diag( vcov( estLogitQuad ) ) ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (standard errors + mean value and standard deviation of age)
urbinEffInt( coef( estLogitQuad ), xMeanQuadInt, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "logit", 
  allCoefVcov = sqrt( diag( vcov( estLogitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

### grouping and re-basing categorical variables
### effects of age changing from the 30-44 category to the 53-60 category
# without standard errors
coefLogitInt <- coef( estLogitInt )
coefLogitInt <- c( coefLogitInt[1:5], 0, coefLogitInt[6] )
xMeanIntShares <- c( xMeanInt[1:5], 1-sum( xMeanInt[3:5] ), xMeanInt[6] )
urbin:::logitEffCat( coefLogitInt, xMeanIntShares, c( 3:6 ), c( -1, -1, 1, 0 ) )
# partial derivatives of the effect wrt the coefficients
# urbin:::logitEffCat( coef( estLogitInt ), xMeanIntAttr, c( 3:5 ), c( -1, -1, 1, 0 ) )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbin:::logitEffCat, t0 = coefLogitInt,
  allXVal = xMeanIntShares, xPos = c( 3:6 ), Group = c( -1, -1, 1, 0 ) )
# with full covariance matrix
vcovLogitInt <- vcov( estLogitInt )
vcovLogitInt <- rbind( vcovLogitInt[ 1:5, ], 0, vcovLogitInt[ 6, ] )
vcovLogitInt <- cbind( vcovLogitInt[ , 1:5 ], 0, vcovLogitInt[ , 6 ] )
# urbin:::logitEffCat( coef( estLogitInt ), xMeanInt, c( 3:5 ), c( -1, -1, 1, 0 ),
#   vcov( estLogitInt ) )
# with standard errors only
urbin:::logitEffCat( coefLogitInt, xMeanIntShares, c( 3:6 ), c( -1, -1, 1, 0 ),
  allCoefSE = sqrt( diag( vcovLogitInt ) ) )
### effects of age changing from the 53-60 category to the 38-52 category
# without standard errors
urbin:::logitEffCat( coefLogitInt, xMeanIntShares, c( 3:6 ), c( 0, 1, -1, 1 ) )
# partial derivatives of the effect wrt the coefficients
# urbin:::logitEffCat( coef( estLogitInt ), xMeanIntAttr, c( 3:5 ), c( 0, 1, -1, 1 ) )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbin:::logitEffCat, t0 = coefLogitInt,
  allXVal = xMeanIntShares, xPos = c( 3:6 ), Group = c( 0, 1, -1, 1 ) )
# with full covariance matrix
# urbin:::logitEffCat( coef( estLogitInt ), xMeanInt, c( 3:5 ), c( 0, 1, -1, 1 ),
#   vcov( estLogitInt ) )
# with standard errors only
urbin:::logitEffCat( coefLogitInt, xMeanIntShares, c( 3:6 ), c( 0, 1, -1, 1 ),
  allCoefSE = sqrt( diag( vcovLogitInt ) ) )

