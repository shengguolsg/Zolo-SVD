
   SUBROUTINE PDNORMEST( M,N,A,IA,JA,DESCA,Tol,Cond,CNT )
!
    IMPLICIT NONE
!
    INTEGER, INTENT(IN) :: M,N,IA,JA
    DOUBLE PRECISION, INTENT(IN) :: Tol
    INTEGER, INTENT(INOUT) :: CNT
    DOUBLE PRECISION, INTENT(INOUT) :: Cond
!  
    INTEGER             :: DESCA(*)
    DOUBLE PRECISION    :: A(*)
!
! ==========
!           
!   This routine simulates MATLAB routine normest on the distributed
!   paralel platform. The input matrix is A, and we estimate the 2-norm
!   of matrix A.  
! 
!   NORMEST Estimates the matrix 2-norm.
!   NORMEST(S) is an estimate of the 2-norm of the matrix S.
!   NORMEST(S,tol) uses relative error tol instead of 1.e-6.
!   [nrm,cnt] = NORMEST(..) also gives the number of iterations used.
!
!   This function is intended primarily for sparse matrices,
!   although it works correctly and may be useful for large, full
!   matrices as well.  Use NORMEST when your problem is large
!   enough that NORM takes too long to compute and an approximate
!   norm is acceptable.

! ==================
! Written by Shengguo Li on 2018-05-23, Changsha
!
! ===============================================

    ! Local Variables 
    INTEGER  :: IT, MaxIter, Info; 
    DOUBLE PRECISION :: E0, Alpha, Beta, NormX, NormSx;

    ! Local Workspace
    DOUBLE PRECISION, ALLOCATABLE :: x(:), Ax(:)

    ! x = sum(abs(A),1)';
    
    //e0 = pdlange_ ( "1", &M, &N, A, &i1, &i1, descA, W); 
    double *w  = (double *)malloc(1*sizeof(double)) ;
    
    *e = pdlange_ ( "f", &N, &i1, W, &i1, &i1, descW, w);
    //pdnrm2_( &N, e, W, &i1 , &i1, descW , &i1);
    
    if (*e == 0){ return;}
    //x = x/e;
    alpha = 1.0;
    pdlascl_( "G", e, &alpha, &N, &i1, W, &i1, &i1, descW, &info);

    e0 = 0;
    while ( (cnt < maxiter) &&
           (fabs((*e) - e0) > (tol * (*e))) )
   {
        e0 = *e; alpha = 1.0; beta = 0.0;
        pdgemv_ ("N", &M, &N, &alpha, A, &i1, &i1, descA, W, &i1, &i1, descW, &i1, &beta, Sx, &i1, &i1, descSx, &i1);
        normSx = pdlange_ ( "f", &N, &i1, Sx, &i1, &i1, descSx, w);

        //if nnz(Sx) == 0
        //    Sx = rand(size(Sx),class(Sx));
        //end

        pdgemv_ ("N", &M, &N, &alpha, A, &i1, &i1, descA, Sx, &i1, &i1, descSx, &i1, &beta, W, &i1, &i1, descW, &i1);
        normx = pdlange_ ( "f", &N, &i1, W, &i1, &i1, descW, w);
   
        *e = normx/normSx;
        pdlascl_( "G", &normx, &alpha, &N, &i1, W, &i1, &i1, descW, &info);
        cnt = cnt+1;
        if ( (cnt >= maxiter) &&
            (fabs((*e) - e0) > (tol * (*e))) ) {
            fprintf(stderr, "normest: didn't converge\n");
        }
    }
    return;
}
