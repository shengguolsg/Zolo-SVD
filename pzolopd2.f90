!
    SUBROUTINE PZOLOPD2( SYM,RK,M,N,ALPHA,BETA,As,IA,JA,DESCA,Q,IQ,JQ,DESCQ,INFO,POLAR )
!
      use auxil_zolo
      use mod_check_routines
!
      IMPLICIT NONE
      include 'mpif.h' 
!
!  -- FASTPACK routine (version 0.1) --
!  National University of Defense Technology
!  Feb 16, 2017
!
! .. Scalar Arguments ..
    CHARACTER         SYM,POLAR
    INTEGER           RK, M, N, IA, JA, IQ, JQ, INFO
    DOUBLE PRECISION  ALPHA, BETA
! ..
! .. Array Arguments ..
    INTEGER           DESCA(*), DESCQ(*)
    DOUBLE PRECISION  As(*), Q(*)
!
! ..
! PURPOSE
! ========
!
! This routine computes the poloar decomposition of matrix A-SIGMA*I by using Zolo-PD. 
! It uses the distributed parallel version and uses rk communicators. On entry,
! Matrix A is distributed over all processes. On exit, matrix Q stores the
! computed polar factor over all processes. 
! 
! This routine is rewritten from pzolopd1.f90, and it is aimed for
! optimizing it. 
! 
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
!  RK      (global input) INTEGER 
!          The number of subgroups to communicators. 
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
!  ALPHA   (global input) DOUBLE PRECISION
!          The upper bound of the singular values of As.
!
!  BETA    (global input) DOUBLE PRECISION
!          The lower bound of the singular values of As.
!
!  As      (local input) DOUBLE PRECISION pointer into the
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
! Further Details
! ===============
! 2017.02.19
!
! 1) This routine only works for M == N. 
!
!
! This routine is written by Shengguo Li
!
! 2) 2017-12-27
!    
! The difference between this routine and pzolopd.f90 is that
! this routine use alpha and beta as input parameters, while they are computed in
! pzolopd.f90. 
!
! =============================================

   integer :: nblk, na
   ! 
   ! .. Parameters ..
   INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                      LLD_, MB_, M_, NB_, N_, RSRC_
   PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                        CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                        RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
   real*8, parameter :: ZERO = 0.d0, ONE = 1.d0, TWO=2.0D+0, TEN=1.0D+01
   !-------------------------------------------------------------------------------
   !  Local Variables
   INTEGER :: IT 
   INTEGER :: NP_ROWS, NP_COLS, NPROCS, sep_nprows, sep_npcols, sep_nprocs, myid, &
              MY_PROW, MY_PCOL, mpierr, na_rows, na_cols, & 
              na_rows1, na_cols1, nq_rows1, i, j, LWORK, LWORK1, liwork, IP, &
              nq, npr_dummy, npc_dummy
   INTEGER :: tp_nprows, tp_npcols, tp_myprow, tp_mypcol
   REAL*8  :: err, errmax,ttt0, ttt1, ttt2, CON, TOL, TOL1, TOL2, eps, FCON, FONE, RRk, &
              ERR2, lowk, scale
   !
   INTEGER :: ALL_CONTXT, TOP_CONTXT
   INTEGER :: CONTXT(8), myprow(8), mypcol(8), Sep_DESC(9)
   !
   ! There are rk+2 BLACS_CONTXT or communicators. At the top level, all the processes are
   ! organized into two grids, one is DESCA(CTXT_) and another one is for 
   ! for implementing the communications among rk sub-communicators. The other rk 
   ! sub-communicators are used for parallel computations. 
   ! 
   ! Vectors for computing zolo-function
   REAL*8  :: coeff(24), MaxCoeff
   ! coeff(1:2*rk) stores parameters cj, coeff(2*rk+1:3*rk) stores aj. 
   ! 
   REAL*8, ALLOCATABLE  :: A(:), Z(:), TA(:), WORK(:), WORK1(:), TAU(:)
   INTEGER, ALLOCATABLE :: IPIV(:), UMAP(:,:,:)
   !
   ! Z is used as workspace, usually a temple copy of A;
   ! TA is the matrix A distributed among fewer processes; 
   ! work1 and TAU are used as workspaces. 
   ! umap is used for splitting the communicators. 
   !
   DOUBLE PRECISION     :: rtmp(4)
   ! ..
   ! .. External Functions ..
   INTEGER, EXTERNAL :: NUMROC
   DOUBLE PRECISION, EXTERNAL :: DLAMCH, PDLANGE

   integer :: STATUS
!-------------------------------------------------------------------------------
!
!  .. Executable Statements ..
!
!  Get grid parameters
   ALL_CONTXT = DESCA( CTXT_ )
   NBLK  = DESCA( MB_ ) 
   CALL BLACS_GRIDINFO( ALL_CONTXT, NP_ROWS, NP_COLS, MY_PROW, MY_PCOL )
   CALL BLACS_PINFO(MYID, NPROCS)

   na     = N
   STATUS = 0
   !SIGMA  = ZERO   ! the shift
!
   ! Determine the necessary size of the distributed matrices,
   na_rows = NUMROC( M, NBLK, MY_PROW, 0, NP_ROWS )
   na_cols = NUMROC( N, NBLK, MY_PCOL, 0, NP_COLS )

   !-------------------------------------------------------------------------------
   ! Allocate workspace 
   ALLOCATE( TA(na_rows*na_cols),IPIV(N) )

   ! **********************************************************
   !       Initialize Matrix TA and estimate alpha and beta    *
   ! **********************************************************
   
   ! Use all the processes to estimate the lower and upper bounds of singular values
   CALL PDLACPY( 'FULL',M,N,As,1,1,DESCA,Q,1,1,DESCQ )

   !CON = ALPHA / BETA
   !CON = ONE / BETA
   !scale = ALPHA

   scale = ALPHA*BETA
   CON = ONE / BETA

   if( myid==0 ) write(*,*) 'scale=', scale

   CALL PDLASCL( 'General',scale,ONE,M,N,Q,1,1,DESCQ,info )

!   CALL PDLAPRNT( na,na,As,1,1,DESCA,1,1,'A1',6, Q )

   ! **********************************************************
   !       Initialize another TOP level grid process TOP2     *
   ! **********************************************************
   sep_nprocs = NPROCS / rk   

   CALL BLACS_GET( -1,0,TOP_CONTXT )
   tp_nprows = rk
   tp_npcols = sep_nprocs
   CALL BLACS_GRIDINIT( TOP_CONTXT,'R',tp_nprows,tp_npcols )
   CALL BLACS_GRIDINFO( TOP_CONTXT,tp_nprows,tp_npcols,tp_myprow,tp_mypcol )

! **************************************************************
!          Split the processes into rk groups                  * 
!          r=2 or 3. Choose nprocs carefully.                  *
! **************************************************************
   DO sep_npcols = NINT(SQRT(REAL(sep_nprocs))),2,-1
      IF( mod(sep_nprocs,sep_npcols)==0 ) EXIT
   END DO
   sep_nprows = sep_nprocs/sep_npcols
   ALLOCATE( UMAP(rk,sep_nprows,sep_npcols ) )

   !IF( myid == 0 ) WRITE(*,*) 'The subgrid is sepnrow', sep_nprows, sep_npcols

   ! We need (rk+1) grids. top_contxt is a (rk)-by-(sep_nprocs)
   DO IP = 1, rk
      DO I = 1, sep_nprows
         DO J = 1, sep_npcols
            UMAP(IP,I,J) = (IP-1)*sep_nprocs+(I-1)*sep_npcols+(J-1) 
         END DO
      END DO
      CALL BLACS_GET( ALL_CONTXT,0,CONTXT(IP) )
      CALL BLACS_GRIDMAP( CONTXT(IP),UMAP(IP,:,:),sep_nprows,sep_nprows,sep_npcols )
   END DO
   DEALLOCATE( UMAP )
  
   myprow(1:rk) = -1
   mypcol(1:rk) = -1
   DO IP = 1, rk
      CALL BLACS_GRIDINFO( CONTXT(IP),npr_dummy,npc_dummy,myprow(IP),mypcol(IP) )
   END DO

   !Allocate the workspace for each sub-communicator
   nq = 2*na
   DO IP = 1,rk 
      IF( myprow(IP).GT.-1 .AND. mypcol(IP).GT.-1 ) THEN
         na_rows1 = NUMROC( M, nblk, myprow(IP), 0, sep_nprows )
         na_cols1 = NUMROC( N, nblk, mypcol(IP), 0, sep_npcols )
         nq_rows1 = NUMROC( nq,nblk, myprow(IP), 0, sep_nprows )
         LWORK1 = MAX( nq_rows1*na_cols1, 2*na_rows1*na_cols1 )+ 8*na
      END IF
   END DO

   ! Some process may only use partial of the workspace 
   ALLOCATE( A(na_rows1*na_cols1),Z(na_rows1*na_cols1),WORK1(LWORK1),TAU(na) )

! **************************************************************
!       Each process read the data and construct matrix A      * 
!                                                              *
! **************************************************************

   CALL DLASET( 'A',na_rows1,na_cols1,ZERO,ZERO,A,na_rows1 )
   IP = 1
   IF( myprow(IP).GT.-1 .AND. mypcol(IP).GT.-1 ) THEN
!      print *, 'myid===', myid, na_rows1
      CALL DESCINIT( Sep_DESC,na,na,nblk,nblk,0,0,CONTXT(IP),na_rows1,info )
   ELSE
      Sep_DESC(CTXT_) = -1
   END IF

   !  Redistribute the data to each process in each sub-communicator
   ttt1= MPI_Wtime()
   CALL PDGEMR2D( na,na,Q,1,1,DESCQ,A,1,1,Sep_DESC,ALL_CONTXT )
   ttt2 = MPI_Wtime()
   IF(myid==0) write(*,*) 'Data redistribution1 costs', ttt2-ttt1

   DO IP = 0, tp_npcols-1
      IF( tp_mypcol == IP ) THEN
         CALL DGSUM2D( TOP_CONTXT,'COLUMN',' ',na_rows1,na_cols1,A,na_rows1,-1,IP )
      END IF
   END DO


! **************************************************************
!                      Main iterative steps                    * 
!                                                              *
! **************************************************************
   EPS  = DLAMCH('Precision'); 
   Tol  = EPS
   TOL1 = TEN*EPS/TWO
   RRk  = dble(rk)*2+1 
   TOL2 = EXP( LOG(TOL1)/RRk ) 
   
   IT = 0
   ERR = ONE
   !lowk = BETA/ALPHA           ! A lower bound of sigma_min 

   lowk = BETA                  ! A lower bound of sigma_min 
   ERR2 = ABS( ONE-lowk) 

   ttt0 = MPI_Wtime() 
   DO WHILE ( (ERR.GT.TOL2 .AND. IT.LT.12) .OR. CON.EQ.ONE  )
   !DO WHILE ( IT.LT.3  )

     ! We use explicit restart strategy, and check whether it can make
     ! LF10000 matrix converge.  
     IF( CON .EQ. ONE ) THEN
       ! Recompute the upper and lower bound, and the condition number
       IF( myid ==0 ) write(*,*) 'Use formXtwo'
       CON = ONE/BETA
       lowk = BETA
       ERR2 = ABS( ONE-lowk )
     END IF 

     IF( myid==0 ) write(*,*) 'Iteration No.', IT, 'TOL2=', TOL2, ERR, &
          'lowk=',lowk, 'ERR2=', ERR2

     !Compute the coefficient cj and aj
     CALL Compute_Coeff( CON, Tol, rk, Coeff )
     MaxCoeff = MAXVAL( coeff(1:2*rk) )

     IF( myid==0 ) write(*,*) 'Iteration No.', IT, CON, 'Coeff', MaxCoeff, &
            'ERR=',ERR

     IF ( (IT.EQ.0 .OR. CON.GT.1.0D+03) .AND. MaxCoeff .GT. 3.0D+02 ) THEN

        ! Each group compute the QR factorization independently. 
        ! The matrix X0 = A/SCALE is stored in A, SCALE=ALPHA*BETA.  
        ttt1 = MPI_Wtime()
        DO IP = 1, rk
           IF( myprow(IP).GT.-1 .AND. mypcol(IP).GT.-1 ) THEN
              CALL ZOLOQR( na,na,nblk,A,Z,CONTXT(IP),WORK1,TAU, &
                   coeff,rk,IP,INFO,IT )
           END IF
        END DO ! (IP for QR)
        ttt2 = MPI_Wtime()
        IF( myid==0 ) write(*,*) 'QR factorization finishes.', ttt2-ttt1
!
!       Use BLACS routine DGSUM2D to compute the sum of WORK2, the processes in the same
!       column of Top_contxt. All the processes will obtain a copy of the data. 
!
!       Compute sum_{j=1}^rk gamma*Q_j1*Q_j2**T. 
        ttt1 = MPI_Wtime() 
        DO IP = 0, tp_npcols-1
           IF( tp_mypcol == IP ) THEN
              !write(*,*) 'myid=', myid, 'my_prow', my_prow, my_pcol, 'IP', IP
              CALL DGSUM2D( TOP_CONTXT,'COLUMN',' ',na_rows1,na_cols1,Z,na_rows1,-1,IP )
           END IF
        END DO
        ttt2 = MPI_Wtime() 
        IF( myid==0 ) write(*,*) 'The first combination finishes.',ttt2-ttt1,CON
  
     ELSE

!       Each group compute the Cholesky factorization independently. 
        ttt1 = MPI_Wtime() 
        DO IP = 1, rk
           IF( myprow(IP).GT.-1 .AND. mypcol(IP).GT.-1 ) THEN
              CALL ZOLOCHOL( SYM,na,nblk,FONE,A,Z,CONTXT(IP),WORK1,coeff,rk,IP,INFO )
           END IF
        END DO
        ttt2 = MPI_Wtime() 
        IF( myid==0 ) write(*,*) 'The Cholesky factorization finishes.', ttt2-ttt1

        ! Compute the summation of iternal matrices
        ttt1 = MPI_Wtime() 
        DO IP = 0, tp_npcols-1
           IF( tp_mypcol == IP ) THEN
              CALL DGSUM2D( TOP_CONTXT,'COLUMN',' ',na_rows1,na_cols1,Z,na_rows1,-1,IP )
           END IF
        END DO
        ttt2 = MPI_Wtime() 

        IF( myid==0 ) write(*,*) 'The second combination finishes.',ttt2-ttt1

     END IF ! (QR)

     !Compute FONE
     FONE = ONE
     DO I = 1, rk
        FONE = FONE * ( ONE+Coeff(2*I) )/( ONE+Coeff(2*I-1) )
     END DO

     !---------------------------------------------------
     ! Compute CON for next iteration
     !---------------------------------------------------
     ! Update CON and recompute cj and aj
     FCON = CON
     DO I = 1, rk
        FCON = FCON * ( CON*CON+Coeff(2*I) ) / ( CON*CON+Coeff(2*I-1) )
     END DO
     CON = MAX( FCON/FONE, ONE )
    
     ! Form X2 and check the distance between X1 and X2 
     ! After FormX2, X2 is stored in both A and Z
     ttt1 = MPI_Wtime() 
     DO IP = 1, rk
        IF( myprow(IP).GT.-1 .AND. mypcol(IP).GT.-1 ) THEN
           IF( CON .NE. ONE ) THEN
              CALL FormX2( na,nblk,FONE,ERR,A,Z,CONTXT(IP),WORK1,INFO,IT,IP )
           ELSE
              IF( myid ==0 ) write(*,*) 'Use formXtwo'
              CALL FormXTwo( na,nblk,FONE,ERR,A,Z,CONTXT(IP),WORK1,INFO,IT,IP,&
                   TOL2,ALPHA,BETA )
           END IF
        END IF
     END DO ! (IP for FormX2)
     ttt2 = MPI_Wtime() 
     IF(myid ==0 ) write(*,*) 'FormX2 costs', ttt2-ttt1, ERR, ttt2-ttt0
  
     IT = IT + 1
     
  END DO ! (WHILE)
!
!   After calling FormX2, A stores the polar factor. We only need to redistribute it
!   to the whole processes, and the polar factor would be stored in TA. 
    IP = 1
    ttt1 = MPI_Wtime()
    CALL PDGEMR2D( na,na,A,1,1,Sep_DESC,Q,1,1,DESCA,ALL_CONTXT )
    ttt2 = MPI_Wtime()
    IF( myid == 0 ) write(*,*) 'Data redistribution2 costs', ttt2-ttt1
 
!    CALL PDLAPRNT( na,na,TA,1,1,DESCA,1,1,'A1',6, WORK1 )
  
!  Let one group check the correctness of computed PD. 
!   CALL ZoloCHK( na,nblk,TA,Z,ALL_CONTXT,MPI_COMM_WORLD,INFO,As,WORK1,LWORK1,TAU )
!   ttt1 = MPI_Wtime()
   IF( myid==0 ) write(*,*) 'zolopd time costs', ttt1-ttt0
   
    ! Compute the polar factor H and store it in WORK2
!     IF( LSAME(POLAR,'T') ) THEN
!       !H = Q'*A; H = (H'+H)/2
!       CALL PDGEMM( 'T','N', N,N,N,ONE,Q,IQ,JQ,DESCQ,As,IA,JA,DESCA,&
!            ZERO,WORK,IQ,JQ,DESCQ )
!       CALL PDLACPY( 'Full',N,N,WORK,IQ,JQ,DESCQ,WORK2,IQ,JQ,DESCQ )
!       CALL PDGEADD( 'T',N,N,HALF,WORK,IQ,JQ,DESCQ,HALF,WORK2,IQ,JQ,DESCQ )
!     ENDIF
   
   ! Check the orthogonality of the polar factor
    call CheckOrth( na,na,myid,my_prow,my_pcol,np_rows,np_cols, &
                       DESCA,Q )

   deallocate(TA)
   deallocate(A)
   deallocate(Z)
   deallocate(work1)
   deallocate(Tau)

!   
 END SUBROUTINE PZOLOPD2

!-------------------------------------------------------------------------------
