!
   SUBROUTINE PQDWHSUBIT( N,A,IA,JA,DESCA,NORMA,Q,IQ,JQ,DESCQ,WORK,& 
     LWORK,TAU,WORK2,LWORK2,Rk,INFO )
!
        IMPLICIT NONE
        include 'mpif.h'
!
!  -- FASTPACK routine (version 0.1) --
!     National University of Defense Technology
!     July 28, 2016
!
!     .. Scalar Arguments ..
      INTEGER            N, IA, JA, IQ, JQ, LWORK, LWORK2, INFO, Rk
      DOUBLE PRECISION   NORMA
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCQ( * )
      DOUBLE PRECISION   A( * ), Q( * ), WORK( * ), TAU( * ), WORK2( * )
!     ..
!
!  Purpose
!  =======
!
!  PQDWHSUBIT uses subspace iteration to compute an invariant subspace of
!  matrix C = (1/2)*(Q+I) which is an N-by-N matrix. Q is assumed to be
!  an orthogonal polar factor of a matrix A which is a square matrix.
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
!  N       (global input) INTEGER
!          The number of columns to be operated on, i.e. the number of
!          columns of the distributed submatrix sub( A ). N >= 0.
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
!  NORMA   (global input) DOUBLE PRECISION
!          The norm of matrix, usually Frobenius norm.
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
!          On exit, the first Rk columns of WORK(1) contains V1; the first Rk rows and
!          Rk columns of WORK(IW) contains A1 = V1**T*A*V1, and the last (N-Rk) rows and
!          (N-Rk) columns of WORK(IW) contain A2=V2**T*A*V2. 
!
!  LWORK   (global or local input) INTEGER
!          The dimension of the array WORK. It is due to PDGEQRF (for an N-by-N matrix),
!          and it should be large enough to contain two A's. The least size is returned from
!          pqdwhfacs.f90. 
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
!  Further Details
!  ===============
!  2016.08.15
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
      INTEGER            I, ICTXT, J, JN, K,ID,IPQ,LWIN,LWIN2, JJD,IDCOL,&
                         MYCOL, MYROW, NPCOL, narows, nacols, IID,IDROW, &
                         NPROW,NP,LIWORK,ierr,nii,njj,pp,myid,MB,NB,nxcols, &
                         EM,EN, it,kk,IW
      DOUBLE PRECISION   TOL1, EPS, t0, t1, err, err2, ALPHA
!     ..
!     .. Local Arrays ..
!     ..       
!     External Subroutines ..
      EXTERNAL     PDLACPY, PDLAPRNT, PDMATGEN2
!     ..
!     External Functions ..
      INTEGER, EXTERNAL :: NUMROC, INDXL2G
      LOGICAL, EXTERNAL :: LSAME
      DOUBLE PRECISION, EXTERNAL :: DLAMCH, PDLANGE
!     ..
!     .. Intrinsic ..
      INTRINSIC   SQRT, MIN
!     .. Executable Statements ..
!
!     Get grid parameters       
      ICTXT = DESCA( CTXT_ )
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

!     Test the input parameters
!
      INFO = 0
!      
      IF( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICTXT, 'PQDWHSUBIT', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
      IF( N.EQ.0 ) &
        RETURN      
!      
      narows = NUMROC( N,  MB, MYROW, 0, NPROW )
      nacols = NUMROC( N,  NB, MYCOL, 0, NPCOL )
      EPS  = DLAMCH('Precision')
      TOL1 = TEN*EPS/TWO
      pp = 3              ! the oversampling parameter
!
!     Construct the matrix C = (Q+I)/2, and it is stored in C.
      DO I = 1, narows
         DO J = 1, nacols
            nii = INDXL2G( I, MB, MYROW, 0, NPROW  )
            njj = INDXL2G( J, NB, MYCOL, 0, NPCOL  )
            IF ( nii == njj ) THEN
               IPQ = I +(J-1)*narows
               ! Q(IPQ) = Q(IPQ) + ONE
               Q(IPQ) = Q(IPQ) - ONE
            END IF
         END DO
      END DO
      CALL PDLASCL( 'General',TWO,ONE,N,N,Q,IQ,JQ,DESCQ,info )
!     Estimate the rank of C 
      ALPHA = PDLANGE( 'Fro',N,N,Q,IQ,JQ,DESCQ,WORK )
      ALPHA = ALPHA*ALPHA
      Rk = NINT( ALPHA )
      kk = MIN( Rk  + pp, N )
      nxcols = NUMROC( kk,  NB, MYCOL, 0, NPCOL )
      IT = 1

!      WRITE(*,*) 'alpha in subit', alpha, 'Rk=', Rk, myrow, mycol, nprow, npcol
!      
!     Generate an intial random matrix X \in R^{N\times kk}, and store it 
!     in WORK, which is assumed to be an N-by-N matrix. 
!
! *** Here the value myid may make trouble when using multiple communicators ***
! *** For the problem of Zolo-SVD, the computation of myid is more difficult.***
!      myid = MYROW*NPROW+MYCOL+1
!       CALL PDMATGEN2( ICTXT,'R','NoDiag',N,kk,MB,NB,WORK,narows,0,0,myid,0,narows,&
!           0,nxcols,MYROW,MYCOL,NPROW,NPCOL )
!
     call random_number( work( 1:narows*nacols ))
!         
!     Use the matrix C to construct an intial matrix X        
!     CALL PDLACPY( 'Full',N,kk,Q,1,1,DESCQ,WORK,1,1,DESCA )
!
!     Compute the QR factorization of X, and generate V=[V1 V2]
!
 100  CALL PDGEQRF( N,kk,WORK,IA,JA,DESCA,TAU,WORK2,LWORK2,INFO ) 
      CALL PDORGQR( N,kk,kk,WORK,IA,JA,DESCA,TAU,WORK2,LWORK2,INFO )
      ! X = C*X
      CALL PDGEMM( 'N','N',N,kk,N,ONE,Q,IQ,JQ,DESCQ,WORK,1,1,DESCA,ZERO,&
           WORK2,1,1,DESCA )
      CALL PDGEQRF( N,kk,WORK2,1,1,DESCA,TAU,WORK,LWORK,INFO )
      CALL PDORGQR( N,N,kk,WORK2,IA,JA,DESCA,TAU,WORK,LWORK,INFO )
!
!     Check the error
      IW = narows*nacols+1
      CALL PDGEMM( 'N','N',N,N,N,ONE,A,IA,JA,DESCA,WORK2,1,1,DESCA,ZERO,&
           WORK(IW),1,1,DESCA )  ! E=A*V
      CALL PDGEMM( 'T','N',N,N,N,ONE,WORK2,1,1,DESCA,WORK(IW),1,1,DESCA, &
           ZERO,WORK,1,1,DESCA ) ! E=V**T *A*V.
      ! It is assumed that WORK is big enough to store two matrix A.
      ! WORK stores A1 and A2, and WORK2 stores V1 and V2. 
      EM = Rk
      EN = N-EM
      err = PDLANGE( 'Fro',EM,EN,WORK,1,EM+1,DESCA,WORK(IW) )
      err = err / NORMA
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         WRITE(*,*) 'Iteration ', IT, err, ALPHA, Rk, nprow, npcol
      END IF

      ! This code implements at most twice subspace. If it doesn't converge, it
      ! returns error.
      IF( err .LE. FOUR*TOL1 ) THEN
         ! WRITE(*,*) 'Iteration ', IT
         CALL PDLACPY( 'Full',N,N,WORK2,1,1,DESCA,Q,IQ,JQ,DESCQ )
         CALL PDLACPY( 'Full',N,N,WORK,1,1,DESCA,A,IA,JA,DESCA )
         RETURN
      ELSEIF( IT .EQ. 2 ) THEN
         IF(MYROW.EQ.0 .AND. MYCOL.EQ. 0 ) THEN
            WRITE(*,*) 'The backward error of subspace iteration is', err
         END IF
         CALL PDLACPY( 'Full',N,N,WORK2,1,1,DESCA,Q,IQ,JQ,DESCQ )
         CALL PDLACPY( 'Full',N,N,WORK,1,1,DESCA,A,IA,JA,DESCA )
      ELSE
         ! Copy V1 to WORK2
         CALL PDLACPY( 'Full',N,Rk,WORK2,1,1,DESCA,WORK(IW),1,1,DESCA ) 
         ! X = C*V1
         CALL PDGEMM( 'N','N',N,Rk,N,ONE,Q,1,1,DESCQ,WORK(IW),1,1,DESCA,&
              ZERO,WORK,1,1,DESCA )
         kk = Rk
         IT = IT + 1 
         GO TO 100
      END IF
      
      RETURN 
 
    END SUBROUTINE PQDWHSUBIT
