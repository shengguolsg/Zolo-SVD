  
    DOUBLE PRECISION FUNCTION hfunc( x )
!
    implicit none
!
    DOUBLE PRECISION :: x
!
!   Purpose
!  ============
!  This function is used to compute the parameter ak in QDWH algorithm.
!
!  =====================
!  Written by Shengguo Li, NUDT, China.
!  =============================
!
   DOUBLE PRECISION, PARAMETER :: FOUR=4.0D+0, EIGHT=8.0D+0, TWO=2.0D+0, &
                    ONE = 1.0D+0, THREE=3.0D+0, ZERO=0.0D+0

   DOUBLE PRECISION :: tmp, ww, sum1, x2, eps, tol
!
   INTRINSIC   :: SQRT
   DOUBLE PRECISION, EXTERNAL :: DLAMCH
!
   
    
    EPS  = DLAMCH('Precision');
    TOL = 30*EPS
    x2 = x*x
    if ( abs(x-one) .le. Tol ) then
       ww = ZERO
    else
        tmp = four*(ONE - x2) / (x2*x2)  
        ww  = EXP( LOG( tmp )/ THREE )  
    endif

    sum1 = SQRT( ONE+ ww ) 

    hfunc = SQRT( EIGHT-FOUR*ww + EIGHT*(TWO-x2)/(x2*sum1) )  
    hfunc = sum1+ hfunc / TWO
    

    END FUNCTION hfunc
