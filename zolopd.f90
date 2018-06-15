!
      SUBROUTINE ZOLOPD( M,N,SIGMA,A,IA,JA,DESCA,Q,IQ,JQ,DESCQ,WORK,LWORK,DESCW, &
          TAU,WORK2,LWORK2,INFO,POLAR )
!
        IMPLICIT NONE
        include 'mpif.h'
!
!  -- FASTPACK routine (version 0.1) --
!     National University of Defense Technology
!     Sep. 1, 2016
!
!     .. Scalar Arguments ..
      CHARACTER          POLAR  
      INTEGER            M, N, IA, JA, IQ, JQ, LWORK, LWORK2, INFO
      DOUBLE PRECISION   SIGMA
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCQ( * ), DESCW( * )
      DOUBLE PRECISION   A( * ), Q( * ), WORK( * ), TAU( * ), WORK2( * )
!     ..
!
!  Purpose
!  =======
!
!  PQDWHFAC computes a polar factorization of a real distributed M-by-N
!  matrix sub( A )-SIGMA = A(IA:IA+M-1,JA:JA+N-1) -SIGMA*i = Q*H. QDWH denotes QR-based
!  dynamically weighted Halley iteration. This algorithm is based on the
!  Matlab routine written by Yuji. 
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!
!  Arguments
!  =========
!
!  M       (global input) INTEGER
!          The number of rows to be operated on, i.e. the number of rows
!          of the distributed submatrix sub( A ). M >= 0. It is assumed that
!          M >= N. 
!
!  N       (global input) INTEGER
!          The number of columns to be operated on, i.e. the number of
!          columns of the distributed submatrix sub( A ). N >= 0.
!
!  SIGMA   (global input) DOUBLE PRECISION
!          The shift, an constant. 
!
!  A       (local input) DOUBLE PRECISION pointer into the
!          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
!          On entry, the local pieces of the M-by-N distributed matrix
!          sub( A ) which is to be factored.  On exit, the elements are
!          unchanged. 
!
!  IA      (global input) INTEGER
!          The row index in the global array A indicating the first
!          row of sub( A ).
!
!  JA      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( A ).
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!
!  Q       (local output) DOUBLE PRECISION pointer into the
!          local memory to an array of dimension (LLD_Q, LOCc(JQ+N-1)).
!          Q is originally copied from A, and it iteratively computes
!          the polar factor, the current one, and A would be the previous one.  
!
!  IQ      (global input) INTEGER
!          The row index in the global array Q indicating the first
!          row of sub( Q ). It is useless right now. 
!
!  JQ      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( Q ). It is useless right now. 
!
!  WORK    (local workspace/local output) DOUBLE PRECISION array,
!                                                     dimension (LWORK)
!          On exit, WORK(1) returns the minimal and optimal LWORK.
!
!  DESCW   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix W. DESCW is for 
!          QR factorization of an 2n x n matrix.
!
!  LWORK   (global or local input) INTEGER
!          The dimension of the array WORK. It is due to PDGECON (for an N-by-N matrix)
!          The least size is computed by calling PDGECON with LWORK=-1. 
!           
!          If LWORK = -1, then LWORK is global input and a workspace query is assumed;
!          the routine only calculate the minimum and optimal size of workspace. 
!          LWORK is return in WORK(1). 
!
!  TAU     (local output)  DOUBLE PRECISION array, dimension
!          LOCc(JA+MIN(M,N)-1). This array contains the scalar factors
!          TAU of the elementary reflectors. TAU is used for PDEGEQR and
!          PDORGQR. 
!      
!  WORK2    (local workspace/local output) DOUBLE PRECISION array,
!                                                     dimension (LWORK2)
!          On exit, WORK(2) returns the minimal and optimal LWORK2.
!      
!          If POLAR = 'True', WORK2 stores the computed polar factor H. Otherwise,
!          It is used as a workspace. 
!
!  LWORK2   (global or local input) INTEGER
!          The dimension of the array WORK2. It is due to PDGEQRF for an
!          (M+N)-by-N matrix (plus the workspace for TAU). It is computed by PDGEQRF. 
!           
!          If LWORK2 = -1, then LWORK is global input and a workspace query is assumed;
!          the routine only calculate the minimum and optimal size of workspace. 
!          LWORK2 is returned in WORK(2). 
!
!  INFO    (global output) INTEGER
!          = 0:  successful exit
!          < 0:  If the i-th argument is an array and the j-entry had
!                an illegal value, then INFO = -(i*100+j), if the i-th
!                argument is a scalar and had an illegal value, then
!                INFO = -i.
!
!  POLAR   (global input) CHARACTER 
!          = 'Not'  Do not compute the polar factor H
!          = 'True' Compute the polar factor H and store it in WORK2. 
!
!  Further Details
!  ===============
!  2016.07.29
!  
!  1) This routine only works when M == N.
!  2) We do not check whether A is symmetric or not. We assume it is.
!  3) We do not allow row or column pivoting. 
!
!  This routine is written by Shengguo Li.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                         LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                           CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                           RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION, PARAMETER :: THREE = 3.0D+0, ZERO=0.0D+0, ONE=1.0D+0, &
           TWO = 2.0D+0, TEN = 1.0D+1, FOUR=4.0D+0, NEONE=-1.0D+0, HALF=0.5D+0
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, ICTXT, J, JN, K,ID,IPQ,LDQ,LWIN,LWIN2, JJD,IDCOL,&
                         MYCOL, MYROW, NPCOL, narows, nacols, IID,IDROW,&
                         NPROW, NQ, IT, NP, LIWORK, ierr, nqrows, nqcols, &
                         nii, njj
      DOUBLE PRECISION   TOL1, TOL2, EPS, ALPHA, BETA, t0,t1,ak,bk, &
                         ck,sck,lowk,akk,bkk,err,err2,lowk2
      ! TOL2 is used to check for stopping
      ! ALPHA is an upper bound of sigma_max
      ! BETA is an lower bound of sigma_min
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION     :: rtmp(4)
      INTEGER, ALLOCATABLE :: IPIV(:)
!     ..       
!     External Subroutines ..
      EXTERNAL     PDLACPY, PDGETRF, PDLAPRNT
!     ..
!     External Functions ..
      INTEGER, EXTERNAL :: NUMROC, INDXL2G
      LOGICAL, EXTERNAL :: LSAME
      DOUBLE PRECISION, EXTERNAL :: DLAMCH, PDLANGE, hfunc
!     ..
!     .. Intrinsic ..
      INTRINSIC   SQRT
!     .. Executable Statements ..
!
!     Get grid parameters       
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

!     Test the input parameters
!
      INFO = 0
      IF( M .LT. N ) THEN
         INFO = -1
      END IF
!      
      LQUERY = ( LWORK.EQ.-1 .OR. LWORK2.EQ.-1 )
      IF( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PQDWHFAC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
        RETURN      
!      
      NQ = 2*N
      narows = NUMROC( M,  DESCA( MB_ ), MYROW, 0, NPROW )
      nacols = NUMROC( N,  DESCA( NB_ ), MYCOL, 0, NPCOL )
      nqrows = NUMROC( NQ, DESCQ( MB_ ), MYROW, 0, NPROW )
      nqcols = NUMROC( N,  DESCQ( NB_ ), MYCOL, 0, NPCOL )
!
!    Calculate the least size of WORK
      NP = narows + DESCA( MB_ )
      LIWORK = NP
      ALLOCATE( IPIV(NP), STAT=ierr )
      IF( ierr .ne. 0 ) THEN
         write(*,*) 'Allocate IPIV failed, Processor', MYROW, MYCOL
         return
      ENDIF

      LWIN = -1
      rtmp(:) = ZERO
      CALL PDGECON('1',N,Q,IQ,JQ,DESCQ,ONE,BETA,rtmp,LWIN,IPIV,LIWORK,INFO)
      LWIN = MAX( int(rtmp(1))+1, nqrows*nqcols)
      LWIN = MAX( LWIN, 2*narows*nacols )

      LWIN2 = -1
      call PDGEQRF( NQ, N, WORK, 1, 1, descw, Q, rtmp, LWIN2, INFO )
      IF( info .ne. 0 ) THEN
         PRINT *, 'pdgeqrf returns error during 1st call', info
      END IF
      LWIN2=MAX( int(rtmp(1))+1, narows*nacols )
      WORK(1) = LWIN
      WORK(2) = LWIN2
      
!     Quick return the size of WORK       
      IF( LQUERY ) THEN
         RETURN
      END IF
      IF( LWORK.LT.LWIN ) THEN
         WRITE(*,*) 'WORK is too small, at least', LWIN
         INFO = -12
         RETURN
      ELSEIF( LWORK2.LT.LWIN2 ) THEN
         WRITE(*,*) 'WORK2 is too small, at least', LWIN2
         INFO = -16
         RETURN
      END IF

      EPS  = DLAMCH('Precision'); 
      TOL1 = TEN*EPS/TWO
      TOL2 = EXP( LOG(TOL1)/THREE ) 
!
!     Estimate the lower and upper bounds
      ! Copy A to Q
      CALL PDLACPY('FULL',M,N,A,IA,JA,DESCA,Q,IQ,JQ,DESCQ )
      IF( SIGMA.NE.ZERO ) THEN
!         WRITE(*,*) 'shift is performing', myrow, mycol
         DO I = 1, narows
            DO J = 1, nacols
               nii = INDXL2G( I, DESCA( MB_ ), MYROW, 0, NPROW  )
               njj = INDXL2G( J, DESCA( NB_ ), MYCOL, 0, NPCOL  )
               IF ( nii == njj ) THEN
                  IPQ = I +(J-1)*narows
                  Q(IPQ) = Q(IPQ) - SIGMA
               END IF
            END DO
         END DO
      END IF

      ALPHA = PDLANGE('Fro',M,N,Q,IA,JA,DESCQ,WORK )
      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
          WRITE(*,*) 'ALPHA = ', ALPHA
       END IF
!     Scale the original matrix to form X0
      CALL PDLASCL( 'General',ALPHA,ONE,N,N,Q,IQ,JQ,DESCQ,info )
!
      ! WORK2 is used as a workspace, and it also stores Q of the previous step in
      ! the iterative loop. 
      t0 = MPI_Wtime()
      CALL PDLACPY('FULL',M,N,Q,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ )
      CALL PDGETRF( N,N,WORK2,IQ,JQ,DESCQ,IPIV,INFO )
      t1 = MPI_Wtime()
!     Return if INFO is non-zero.
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
!
      t0 = MPI_Wtime()
      CALL PDGECON( '1', N, WORK2, IQ, JQ, DESCQ, ONE, BETA, WORK, LWORK, &
           IPIV, LIWORK, INFO )
      BETA = BETA / SQRT( DBLE(N) )
      t1 = MPI_Wtime()
      DEALLOCATE(IPIV)
!!$      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!!$          WRITE(*,*) 'PDGECON finishes', BETA, TOL1, TOL2
!!$       END IF

!      CALL PDLAPRNT( N,N,A,IA,JA,DESCA,1,1,'A0',6, WORK )
      IT = 0
      lowk = BETA     ! a lower bound of sigma_min
      err  = ONE
      err2 = ABS( ONE-lowk )
      LDQ = DESCW( LLD_ )
!      LDQ = nqrows

      DO WHILE ( IT.EQ.0 .OR. err.GT.TOL2 .OR. err2.GT.TOL1 )
!
         ak = hfunc(lowk) 
         bk = (ak-ONE)**2 / FOUR
         ck = ak+bk-ONE
         sck = SQRT( ck )
!        Update the lower bound 
         lowk2 = lowk*lowk
         lowk = lowk * (ak+bk*lowk2) / (ONE + ck*lowk2)
         
!!$      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!!$         WRITE(*,*) 'IT =', IT, 'ck=',ck, 'ak', ak, 'bk', bk, ldq
!!$         WRITE(*,*) 'DESCA=', DESCA(1:9)
!!$         WRITE(*,*) 'DESCW=', DESCW(1:9)
!!$      END IF
      
      IF( ck > 100 ) THEN  ! use QR
!        Form matrix [Xk; 1/sck * I], prepare for QR
         CALL PDLASET( 'Full',NQ,N,ZERO,ZERO,WORK,1,1,DESCW )  ! set WORK to ZERO
         CALL PDLACPY( 'Full',N,N,Q,IQ,JQ,DESCQ,WORK,1,1,DESCW )

         ! Append the identity matrix
         DO ID=1, N, 1
            CALL INFOG2L( IQ-1+ID+N,IQ-1+ID,DESCW,NPROW,NPCOL,MYROW,MYCOL,&
                 IID,JJD,IDROW,IDCOL )
            IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
               IPQ = IID + ( JJD-1 )*LDQ 
               WORK( IPQ ) = ONE/sck
            END IF
         END DO
!         IF( IT==0 ) CALL PDLAPRNT( NQ,N,WORK,IA,JA,DESCW,1,1,'Q0',6, WORK2 )

       CALL PDGEQRF( NQ,N,WORK,IA,JA,DESCW,TAU,WORK2, LWORK2, INFO ) 
       CALL PDORGQR( NQ,N,N,WORK,IA,JA,DESCW,TAU,WORK2,LWORK2,INFO )
!       IF(IT ==0) CALL PDLAPRNT( NQ,N,WORK,IA,JA,DESCW,1,1,'Q1',6, WORK2 )

!       CALL PDLACPY( 'Full',N,N,WORK,IA,JA,DESCW,WORK2,1,1,DESCA )
!       CALL PDGEMM( 'N','T',N,N,N,ONE,WORK2,1,1,DESCA,WORK,N+1,1,DESCW, &
       CALL PDGEMM( 'N','T',N,N,N,ONE,WORK,1,1,DESCW,WORK,N+1,1,DESCW, &
	              ZERO,WORK2,IQ,JQ,DESCQ  )
!       IF(IT==0 ) CALL PDLAPRNT( N,N,Q,IA,JA,DESCQ,1,1,'Qm',6, WORK )
        
!       Update Xk  
        akk = ( ak-bk/ck )/ sck
        bkk = bk / ck
        CALL PDGEADD( 'N',N,N,bkk,Q,IQ,JQ,DESCA,akk,WORK2,IQ,JQ,DESCQ )
!        IF(IT==0) CALL PDLAPRNT( N,N,Q,IA,JA,DESCQ,1,1,'X1',6, WORK )

      ELSE  ! use Cholesky factorization when X is well conditioned 
!
!       Compute matrix Zk        
        CALL PDLACPY( 'Full',N,N,Q,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ )
        CALL PDGEMM( 'T','N',N,N,N,ck,Q,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ,&
             ZERO,WORK,IQ,JQ,DESCQ)
      	DO ID=1, N, 1
          CALL INFOG2L( IA-1+ID,JA-1+ID,DESCA,NPROW,NPCOL,MYROW,MYCOL,&
               IID,JJD,IDROW,IDCOL )
          IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
             IPQ = IID + ( JJD-1 )*narows 
             WORK( IPQ ) = ONE + WORK( IPQ )
          END IF
      	END DO 
       
        CALL PDPOSV( 'U',N,N,WORK,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ,INFO)    

!       Update Xk, stored in Q
        akk = ak - bk/ck
        bkk = bk / ck
        CALL PDGEADD( 'N',N,N,bkk,Q,IQ,JQ,DESCQ,akk,WORK2,IQ,JQ,DESCQ )

      END IF  ! (ck > 100)

!     Compute the error and store it in the matrix A
      CALL PDGEADD( 'N',N,N,NEONE,WORK2,IQ,JQ,DESCQ,ONE,Q,IA,JA,DESCA  )
      err = PDLANGE('Fro',N,N,Q,IA,JA,DESCA,WORK )

      IT = IT + 1
      err2 = ABS( ONE -lowk)
!!$      IF(MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
!!$         WRITE(*,*) 'err', err, 'err2', err2
!!$      END IF

      CALL PDLACPY( 'Full',N,N,WORK2,IQ,JQ,DESCQ,Q,IA,JA,DESCA )
 
     END DO ! WHILE

     ! Compute the polar factor H and store it in WORK2
     IF( LSAME(POLAR,'T') ) THEN
        ! H = Q'*A; H = (H'+H)/2
        CALL PDGEMM( 'T','N', N,N,N,ONE,Q,IQ,JQ,DESCQ,A,IA,JA,DESCA,&
             ZERO,WORK,IQ,JQ,DESCQ)
        CALL PDLACPY( 'Full',N,N,WORK,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ )
        CALL PDGEADD( 'T',N,N,HALF,WORK,IQ,JQ,DESCQ,HALF,WORK2,IQ,JQ,DESCQ )
     ENDIF
     
      RETURN 
 
    END SUBROUTINE ZOLOPD
