!
  SUBROUTINE ZOLOCHOL( SYM,N,NBLK,FONE,A,Z,CONTXT,WORK,COEFF,Rk,IP,INFO  )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   CHARACTER:: SYM
   INTEGER  :: N, NBLK, CONTXT, INFO, Rk, IP
   REAL*8   :: FONE
!   
   REAL*8   :: WORK(*), A(*), Z(*), COEFF(*)
!
! ============
!   
!  This subroutine computes the Cholesky factorization of matrix X1*X1**T + c(2j-1)I, where 
!  matrix X1 is computed from the first iteration. At the beginning of this routine, we
!  need to compute X1 first.    
!   
! .. Arguments ..
! ================

! N      INTEGER  The global column dimension of matrix A
!
! nblk   INTEGER  Blocking factor in block cyclic distribution  
!   
! FONE   DOUBLE PRECISION Scalar 
!        It is used to scale matrix (A+Z) to get X1 = (A+Z)/FONE. 
!
! A      DOUBLE PRECISION Array. On entry, it stores the initial matrix X0.
!        On exit, it contains the first iteration matrix X1. 
!
! Z      DOUBLE PRECISION Array, stores the summation of matrices gamma*Q_j1*Q_j2**T. 
!        On exit, it contains the result X2.
!
! WORK  DOUBLE PRECISION Array, used as workspace for Cholesky factorization
!
! CONTXT INTEGER The BLACS context.    
!
! COEFF  DOUBLE PRECISION Array
!        Coeff(1:2*rk) stores parameters cj, coeff(2*rk+1:3*rk) stores the parameters aj. 
! ========================
! Written by Shengguo Li, 09/14/2016, at Beijing
! NUDT, Changsha 410073, China
!
! 2018-01-03
! Fix one bug for nonsymmetric matrices
! 
! ==================================================
!
   INTEGER  :: np_rows, np_cols, myprow, mypcol, na, na_rows, na_cols, &
        I,J,IA,JA,ID,IID,JJD,IDROW,IDCOL,IPQ, myid, nprocs
   DOUBLE PRECISION  :: akk
!   
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0,MONE=-1.0D+0
!    
   INTEGER  :: DESCA(9)
!
!  .. Local Arrays ..
   DOUBLE PRECISION, ALLOCATABLE :: Tmat(:,:) 
!
!  .. External Functions ..
   INTEGER, EXTERNAL :: INDXL2G, NUMROC
   LOGICAL, EXTERNAL :: LSAME
!
   CALL BLACS_PINFO( myid,nprocs )
   CALL BLACS_Gridinfo( CONTXT,np_rows,np_cols,myprow,mypcol )

   na = N
   IA = 1
   JA = 1
   na_rows = NUMROC( na, nblk, myprow, 0, np_rows )
   na_cols = NUMROC( na, nblk, mypcol, 0, np_cols )

   ! Allocate temp workspace for solving linear equations
   ALLOCATE( Tmat(na_rows,na_cols) )

   ! Set up a scalapack descriptor
   CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )

   ! Compute X1 from A, Z and FONE
   ! CALL PDGEADD( 'N',N,N,ONE,Z,1,1,DESCA,ONE,A,IA,JA,DESCA )  ! A = A + Z
!   J = na_rows*na_cols
!   A(1:J)=Z(1:J)+A(1:J)   
!   CALL DLASCL( 'G',IA,JA,FONE,ONE,na_rows,na_cols,A,na_rows,info ) ! A = A / FONE

   ! Compute Zj = X1* X1**T + cj*I
   CALL DLACPY( 'All',na_rows,na_cols,A,na_rows,Z,na_rows )
   CALL PDGEMM( 'T','N',N,N,N,ONE,A,IA,JA,DESCA,Z,IA,JA,DESCA,&
        ZERO,WORK,IA,JA,DESCA )
   DO ID=1, N, 1
      CALL INFOG2L( IA-1+ID,JA-1+ID,DESCA,NP_ROWS,NP_COLS,MYPROW,MYPCOL,&
           IID,JJD,IDROW,IDCOL )
      IF( MYPROW.EQ.IDROW .AND. MYPCOL.EQ.IDCOL ) THEN
         IPQ = IID + ( JJD-1 )*na_rows 
         WORK( IPQ ) = Coeff(2*IP-1) + WORK( IPQ )
      END IF
   END DO

   IF( LSAME(SYM,'S') ) THEN
       CALL PDPOSV( 'U',N,N,WORK,IA,JA,DESCA,Z,IA,IA,DESCA,INFO )
   ELSE
       ! Compute Xk*Inv(Zk)
       CALL PDTRAN( N,N,ONE,Z,1,1,DESCA,ZERO,Tmat,1,1,DESCA )  ! Tmat = Xk^T
       CALL PDPOSV( 'U',N,N,WORK,IA,JA,DESCA,Tmat,IA,JA,DESCA,INFO )
    
       CALL PDTRAN( N,N,ONE,Tmat,1,1,DESCA,ZERO,Z,1,1,DESCA )  ! 
   
   ENDIF

   ! Scale the matrix Z by aj, Z = aj*Z
   akk = MONE*coeff(2*rk+IP)
   !CALL PDLASCL( 'General',ONE,akk,N,N,Z,IA,JA,DESCA,info )  ! Z = aj * Z
   CALL DLASCL( 'G',na_rows,na_cols,ONE,akk,na_rows,na_cols,Z,na_rows,info  )

   DEALLOCATE( Tmat )
   
 END SUBROUTINE ZOLOCHOL
