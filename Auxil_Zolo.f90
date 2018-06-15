!*******************************************************************************
!> author: Shengguo Li
!>         College of Computer Science, National University of Defense Technology,
!>         Changsha 410073, China
!  license: BSD
!  date: 08/31/2016
!
!  Compute the complete elliptic integral K and the function values of the Jacobi elliptic
!  functions sn and cn. The standard approach is based on the Arithmetic Geometric Mean (AGM)
!  method. 
!
!  The following two routines are based Matlab function 'ellipke' (for the computation of complete
!  elliptic integrals) and 'ellipj' (for the computation of function values of the Jacobi elliptic
!  function). These routines have been modified by Yuji Nakatsukasa and Roland W. Freund. 
 
   module Auxil_Zolo

     implicit none

!     use iso_fortran_env, only: error_unit, wp => real64

     public  :: mellipke, mellipj, choosem, compute_coeff

   contains
!*******************************************************************************

!*******************************************************************************
!>
! MELLIPKE computes complete elliptic integral, and this version is modified from 
! Matlab's built-in code ELLIPKE. Rewritten from Yuji's matlab codes for Zolotarev
! functions. 
!
! *** This routine is only designed for ZOLO-PD, and it may not suitable for general purpose ***
!
     
   SUBROUTINE mellipke( alpha, tol, k, e )
!
!
     double precision, intent(in)  :: alpha, tol
     double precision, intent(out) :: k, e
!
     double precision, parameter :: ONE=1.0D+0, ZERO=0.0D+0, TWO=2.0D+0, pi=3.141592653589793D+0
!
     double precision :: m, m1, mm, a0, b0, s0, a1, b1, w1, c1
     integer  :: i1
!     
     m  = sin(alpha) * sin(alpha)
     m1 = cos(alpha) * cos(alpha)

     a0 = ONE
     b0 = cos(alpha)
     s0 = m
     i1 = 0
     mm = ONE

     DO WHILE (mm > tol)
        a1 = (a0+b0)/TWO
        b1 = sqrt( a0*b0 )
        c1 = (a0 -b0 )/TWO
        i1 = i1 + 1
        w1 = (2**i1)*(c1**2)
        mm = w1  ! mm = max( w1(:) )
        s0 = s0 + w1
        a0 = a1
        b0 = b1
     END DO

     k = pi / (TWO*a1)
     e = k*( ONE-s0/TWO )
     
   END SUBROUTINE MELLIPKE
   
   SUBROUTINE MELLIPJ( u,alpha,tol,sn,cn,dn )
!
!  This subroutine returns the values of the Jacobi elliptic functions     
!  sn, cn and dn, evaluated for corresponding elements of argument 'U' and 
!  parameter M. U and M must be arrays of the same size or either can be scalar.
!  As currently implemented, M is limited to 0 <= M <= 1. 
! 
     
     double precision, intent(in)  :: u, alpha, tol
     double precision, intent(out) :: sn, cn, dn
!
     integer,  parameter  ::  mmax=1, mchunk=1000
     double precision, parameter  ::  ONE = 1.0D+0, ZERO=0.0D+0, HALF=0.5D+0
!     
     double precision   ::  m, m1, nin, i, in
!
     double precision   ::  a(mchunk),b(mchunk),c(mchunk), phin(mchunk)
     
     m  = sin(alpha) * sin(alpha)
     m1 = cos(alpha) * cos(alpha)
     cn = u
     sn = cn
     dn = sn

     c(1) = sin(alpha)
     b(1) = cos(alpha)
     a(1) = ONE
     i = 1
     nin = 0
     DO WHILE ( (abs(c(i))>tol ) .AND. (i<mchunk) )
        i = i + 1
        a(i) = ( a(i-1)+b(i-1) )* HALF
        b(i) = sqrt( a(i-1)*b(i-1) )
        c(i) = ( a(i-1)-b(i-1) )* HALF

        IF( ( abs(c(i)).le.tol) .AND. (abs(c(i-1))> tol) ) THEN
           nin = i-1
        END IF
     END DO !(WHILE)

     phin(1) = u
     phin(i) = (2**nin)*a(i)*u

     DO WHILE ( i>1 )
        i = i -1
        phin(i) = phin(i+1)
        IF( nin.GE. i ) THEN
           phin(i) = HALF*( asin( c(i+1) *sin( phin(i+1) )/a(i+1) ) +phin(i+1) )
        END IF
     END DO !(WHILE2)

     sn = sin( phin(1) )
     cn = cos( phin(1) )
     dn = sqrt( ONE-m*sn*sn )
     
   END SUBROUTINE MELLIPJ

   FUNCTION choosem( con )  result(m)
  
     double precision  :: con
     integer   :: m

      if ( con<1.001 ) then
         m = 2
      elseif (con<=1.01) then
         m = 3
      elseif (con<=1.1) then
         m = 4
      elseif (con<=1.2) then
         m = 5 
      elseif (con<=1.5 ) then
         m = 6
      elseif (con<=2) then
         m = 8     ! one-step convergence till here
      elseif (con<6.5) then
         m = 2
      elseif (con<180) then
         m = 3
      elseif (con<1.5*1.0D+4) then
         m = 4
      elseif (con<2*1.0D+6) then
         m = 5
      elseif (con<1*1.0D+9) then
         m = 6
      elseif (con< 3*1.D+12) then
         m = 7 
      else
         m = 8
      end if
   
    END FUNCTION choosem

    SUBROUTINE Compute_Coeff( CON, Tol, rk, coeff)
!
      integer, intent(in)  :: rk
      double precision, intent(in) :: CON, Tol
      double precision, intent(inout) :: Coeff(:)
      !
      double precision, parameter :: ONE = 1.0D+0
      !
      integer           :: i, rk2, j
      double precision  :: kappa, alpha, tmp1, KK, sn, cn, tn
!
      kappa = ONE/ con
      alpha = ACOS( kappa )
      CALL mellipke( alpha,tol,KK,tmp1 )
      rk2 = 2*rk
!$OMP PARALLEL PRIVATE( I,tmp1,sn,cn,tn )
!$OMP DO SCHEDULE(dynamic)
      DO i = 1, 2*rk
         tmp1 = dble(i) * KK / (rk2+ONE)
         call mellipj( tmp1, alpha, tol, sn, cn, tn)
         coeff(i) = sn**2 / cn**2 
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      ! parameter aj
      DO I = 1, rk
         coeff(rk2+I) = ONE
         DO J = 1, rk
            coeff(rk2+I) = coeff(rk2+I) * (coeff(2*I-1)-coeff(2*J) )
            IF ( I.NE.J ) THEN
               coeff(rk2+I) = coeff(rk2+I) / (coeff(2*I-1)-coeff(2*J-1) )
            END IF
         END DO
      END DO

    END SUBROUTINE Compute_Coeff
   
   end module Auxil_Zolo
