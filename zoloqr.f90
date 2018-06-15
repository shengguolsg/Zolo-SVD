!
  SUBROUTINE ZOLOQR( M,N,NBLK,A,Z,CONTXT,WORKQ,TAU,COEFF,Rk,IP,INFO,IT  )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   INTEGER  :: M, N, NBLK, CONTXT, INFO, Rk, IP,IT
!   
   REAL*8   :: WORKQ(*), TAU(*), A(*), Z(*), COEFF(*)
!
!  This subroutine computes the QR factorization of matrix [A; SIGMA*I], where A is a sparse
!  matrix, defined by (INDX, JNDX, RVAL). Right now this routine only works for an square matrix
!  A, M == N. 
!   
! .. Arguments ..
! ================
! M      INTERGER The global row dimension of matrix A
!
! N      INTEGER  The global column dimension of matrix A
!
! nblk   INTEGER  Blocking factor in block cyclic distribution  
!   
! LDA    INTEGER  The local leading dimension of matrix A
!
! LDQ    INTEGER  The local leading dimension of matrix [A; I]
!
! Z      DOUBLE PRECISION Array, workspace
!        Used as a copy of matrix A. On exit, it contains the result Q1*Q2**T. 
!
! WORKQ  DOUBLE PRECISION Array, used to stored the distributed matrix 
!        [A; I], and it is large enough to store two matrix A
!
! TAU    DOUBLE PRECISION Array with dimension (N) 
!        It is used as a workspace.
! 
! CONTXT INTEGER The BLACS context.    
!
! COEFF  DOUBLE PRECISION Array
!        Coeff(1:2*rk) stores parameters cj, coeff(2*rk+1:3*rk) stores the parameters aj. 
!
! ========================
!
! Written by Shengguo Li, 09/13/2016, at Beijing
! NUDT, Changsha 410073, China
! ==================================================
!
   INTEGER  :: np_rows, np_cols, my_prow, my_pcol, na, na_rows, na_cols, &
        I,J,IROW,IID,JJD, IDROW,IDCOL,JCOL,nq,IPQ,IA,JA,LWORK2, &
        nq_rows, nii, njj, ID, myid, nprocs
   DOUBLE PRECISION :: akk
!
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0, MONE=-1.0D+0
!    
   INTEGER  :: DESCA(9), DESCQ(9)
   DOUBLE PRECISION, ALLOCATABLE :: WORK2(:)
!
!  .. External Functions ..
   INTEGER, EXTERNAL :: INDXL2G, NUMROC
!  ..
!  .. External Routines ..
   EXTERNAL :: INFOG2L
!   
   CALL BLACS_PINFO( myid,nprocs )
   CALL BLACS_Gridinfo( CONTXT,np_rows,np_cols,my_prow,my_pcol )

   na = N
   nq = 2*N
   IA = 1
   JA = 1
   na_rows = NUMROC( na, nblk, my_prow, 0, np_rows )
   na_cols = NUMROC( na, nblk, my_pcol, 0, np_cols )
   nq_rows = NUMROC( nq, nblk, my_prow, 0, np_rows )

   ! Set up a scalapack descriptor
   CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )
   CALL DESCINIT( DESCQ, nq, na, nblk, nblk, 0, 0, CONTXT, nq_rows, info )

   ! Print out the scaled matrix A
   !IF( IT==0 .and. ip==2 ) CALL PDLAPRNT( na,na,A,IA,JA,DESCA,1,1,'A3',6,WORKQ )

   ! Form matrix [Xk; sck*I], prepare for QR
   CALL DLASET( 'A',nq_rows,na_cols,zero,zero,WORKQ,nq_rows )
   CALL PDLACPY( 'Full',na,na,A,IA,JA,DESCA,WORKQ,1,1,DESCQ )

   ! Append the identity matrix
!$OMP PARALLEL PRIVATE( ID,IID,JJD,IDROW,IDCOL,IPQ )
!$OMP DO SCHEDULE(dynamic)
   DO ID = 1, na, 1
      CALL INFOG2L( ID+na,ID,DESCQ,np_rows,np_cols,my_prow,my_pcol,IID,JJD,IDROW,IDCOL )
      IF( my_prow.EQ.IDROW .AND. my_pcol.EQ.IDCOL ) THEN
         IPQ = IID + ( JJD-1 )*nq_rows
         WORKQ( IPQ ) = ONE*SQRT( coeff(2*IP-1) )
      END IF
   END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   LWORK2 = -1 
   CALL MPDGEQRF( NQ,na,WORKQ,1,1,DESCQ,TAU,Z,LWORK2,INFO )
   LWORK2 = int(Z(1) ) + 1

   ALLOCATE( WORK2( LWORK2)  )

   CALL MPDGEQRF( NQ,na,WORKQ,1,1,DESCQ,TAU,WORK2,LWORK2,INFO )
   CALL MPDORGQR( NQ,na,na,WORKQ,1,1,DESCQ,TAU,WORK2,LWORK2,INFO )

   DEALLOCATE( WORK2 )
   !
   ! Compute Tmat = Q1*Q2**T, stored in WORK2.
   akk = MONE*coeff(2*rk+IP) / sqrt( coeff(2*IP-1) )
   CALL PDGEMM( 'N','T',na,na,na,akk,WORKQ,1,1,DESCQ,WORKQ,na+1,1,DESCQ, &
        ZERO,Z,1,1,DESCA  )
   !

 END SUBROUTINE ZOLOQR
