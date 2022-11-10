########## R function: OmegaOSull ##########

# For creation of O'Sullivan-type Omega penalty matrices.

# Last changed: 20 MAR 2020

OmegaOSull <- function(a, b, intKnots)
{
   allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
   K <- length(intKnots) ; L <- 3*(K+8)
   xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
              rep(allKnots,each=3)[-c(1,2,L)])/2
   wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
   Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                     outer.ok=TRUE)$design  
   Omega     <- crossprod(Bdd*wts,Bdd)
   
   return(Omega)
}

########## End of OmegaOSull ##########

