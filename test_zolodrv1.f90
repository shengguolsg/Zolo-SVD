!
     program test_zolodrv1
  
  ! 
  !-------------------------------------------------------------------------------
  ! This routine tests the polar factorization algorithm based on Zolotarev's function.
  ! It uses rk times more compute resource, and rk equals to 2 or 3 in this routine. 
  ! The number of iterations may be larger than 2, different from the approach in
  ! Nakatsukasa's SIAM Review paper.      
  !
  ! This routine uses PDGEMR2D to redistribute the matrix, and therefore it works for
  ! dense matrix. This would be our final version of the code. 
  ! 
  ! The difference with other test_zolosvd's is that this one works for dense matrix.
  ! 
  !
  !-------------------------------------------------------------------------------
     use mod_prepare_matrix
!       
     implicit none
     include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:     System size
   ! sigma:  the eigenvalues below alpha are computed
   ! nblk:   Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   character :: SYM
   integer :: nblk, na, rk, ncols
   real*8  :: sigma, alpha, beta, temp_beta
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
   INTEGER :: np_rows, np_cols, nprocs, sep_nprows, sep_npcols, sep_nprocs, myid, &
              my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols, mpierr, na_rows, na_cols, & 
              na_rows1, na_cols1, nq_rows1, i, j, lwork, lwork1, liwork, IP, &
              nq, info, nbw, NNZ, iunit
   integer :: STATUS
   REAL*8  :: ttt0, ttt1,ttt2
   !
   INTEGER :: ALL_CONTXT
   INTEGER :: ALL_DESC(9)
   real*8  :: rtmp(4)
   !
   CHARACTER(128) :: filename, arg
   !
   ! There are rk+1 BLACS_CONTXT or communicators. At the top level, all the processes are
   ! organized into two grids, one for estimating the bounds of singular values, and another
   ! for implementing the communications among rk sub-communicators. The other rk 
   ! sub-communicators are used for parallel computations. 
   ! 
   REAL*8, ALLOCATABLE  :: Q(:), As(:), TA(:), WORK(:), TAU(:), WORK1(:)
   INTEGER, ALLOCATABLE :: IPIV(:)
   !
   ! As stores the original matrix among all the processes;
   ! Q is used as workspace, usually a temple copy of A;
   ! ..
   ! .. External Functions ..
   INTEGER, EXTERNAL :: NUMROC
   DOUBLE PRECISION, EXTERNAL :: DLAMCH, PDLANGE
   !-------------------------------------------------------------------------------

!   i = 1
!   call get_command_argument(i, arg)
!   filename = arg
 
!   iunit = 13
!   open( unit=iunit, file=filename )
!   call mmread_nnz( iunit, na, ncols, nnz)
!   close( iunit )

!  random
   na = 100
!   i = 1
!   call get_command_argument(i, na) 

   nblk  = 2 
   rk    = 1 
   !SYM = 'Gen'
   sigma = ZERO   ! the shift
   !-------------------------------------------------------------------------------
   ! MPI Initialization
   call mpi_init( mpierr )
   call mpi_comm_rank( mpi_comm_world,myid,mpierr )
   call mpi_comm_size( mpi_comm_world,nprocs,mpierr )

   STATUS = 0

   ! Initialize the BLACS: refer to user's guide or lbll manual
   CALL BLACS_PINFO( myid,nprocs )

   ! ************************************************
   !       Initialize the TOP process grid          * 
   !          Used for estimate bounds              *
   ! ************************************************
   
   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns. 
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs 
   ! and decrement it until a divisor is found. 
   do np_cols = NINT(SQRT(REAL(nprocs))),2,-1
      if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   np_rows = nprocs/np_cols

   if(myid==0) then
      print *
      print '(a)','Test routine for the Zolo-Pd algorithm'
      print *
      print '(3(a,i0))','Matrix size=', na, 'mpi_comm_world=', mpi_comm_world
      print '(3(a,i0))','Number of processor rows=',np_rows,',cols=',np_cols,',total=',nprocs
      write(*,*) 'sigma=', sigma, 'Rk=', rk
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context for all the processes. TOP1
   ! ALL_CONTXT = MPI_COMM_WORLD
   CALL BLACS_GET( -1,0,ALL_CONTXT )
   CALL BLACS_Gridinit( ALL_CONTXT, 'R', np_rows, np_cols )
   call BLACS_Gridinfo( ALL_CONTXT, np_rows, np_cols, my_prow, my_pcol )

   ! Determine the necessary size of the distributed matrices,
   ! We use the Scalapack tools routine NUMROC for that.
   na_rows = NUMROC( na, nblk, my_prow, 0, np_rows )
   na_cols = NUMROC( na, nblk, my_pcol, 0, np_cols )

   ! Set up a scalapack descriptor
   CALL DESCINIT( ALL_DESC, na, na, nblk, nblk, 0, 0, ALL_CONTXT, na_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem
   ! The symmetric matrix in MM form is only stroed its upper or lower triangular part.
   LWORK = na_rows*na_cols
   ALLOCATE( As(na_rows*na_cols),Q(na_rows*na_cols), TA(na_rows*na_cols),WORK(na_rows*na_cols),&
        TAU(na), IPIV(na) )

   ! **********************************************************
   !                    Initialize Matrix A                   *
   ! **********************************************************
   ! For random matrix.
   call prepare_real_rand( 'N',na,myid,my_prow,my_pcol,np_rows,np_cols,ALL_DESC,As,Q )
   
   ! For sparse matrices
!!   CALL PDLASET( 'A',na,na,ZERO,ZERO,As,1,1,ALL_DESC )
!!   CALL PDLASET( 'A',na,na,ZERO,ZERO,Q,1,1,ALL_DESC )
   !CALL prepare_real_symm_MM( filename,na,nnz,myid,my_prow,my_pcol, &
   !         np_rows,np_cols,all_desc,As,Q,sigma )

!!   CALL prepare_real_Gem( SYM,filename,na,nnz,myid,my_prow,my_pcol, &
!!`            np_rows,np_cols,all_desc,As,Q,sigma )

   IF( myid == 0) write(*,*) "The matrix has been initialized"
   !write(*,*) "The matrix has been initialized", myid

   ! Use all the processes to estimate the lower and upper bounds of singular values
   CALL PDLACPY( 'FULL',na,na,As,1,1,ALL_DESC,Q,1,1,ALL_DESC )

   ! Allocate the workspace for PDGECON
   LWORK1 = -1
   LIWORK = na 
   rtmp(:) = ZERO
   CALL PDGECON( '1',na,Q,1,1,ALL_DESC,ONE,BETA,rtmp,LWORK1,IPIV,LIWORK,INFO )
   LWORK1 = INT( rtmp(1) )+1
   ALLOCATE( WORK1(LWORK1) )

   ttt0 = MPI_Wtime() 
   ! The Frobenius norm of matrix A-sigma*I 
   ALPHA = PDLANGE( 'Fro',na,na,Q,1,1,ALL_DESC,WORK1 )
   IF( myid.eq.0 ) WRITE(*,*) 'Alpha =', alpha

   CALL PDLASCL( 'General',ALPHA,ONE,na,na,Q,1,1,ALL_DESC,info )
   !TEMP_BETA = PDLANGE('1',na,na,Q,1,1,ALL_DESC,WORK )  

   CALL PDGETRF( Na,Na,Q,1,1,ALL_DESC,IPIV,INFO )
   CALL PDGECON( '1',na,Q,1,1,ALL_DESC,ONE,BETA,WORK1,LWORK1, &
           IPIV,LIWORK,INFO )
   !BETA = TEMP_BETA*BETA/ SQRT( DBLE(na) ) 
   BETA = BETA/ SQRT( DBLE(na) ) 

   DEALLOCATE( IPIV,WORK1 )
   !Beta = Beta*alpha

   !IF( myid.eq.0 ) WRITE(*,*) 'beta =', beta, 'temp_beta=', temp_beta
   IF( myid.eq.0 ) WRITE(*,*) 'beta =', beta

   ttt0 = MPI_Wtime()
   call pzolopd1( SYM,rk,na,na,alpha,beta,As,1,1,ALL_DESC,Q,1,1,ALL_DESC,INFO,'T' )
   !call pzolopd2( SYM,rk,na,na,alpha,beta,As,1,1,ALL_DESC,Q,1,1,ALL_DESC,INFO,'T' )
  
!  Let one group check the correctness of computed PD. 
   CALL ZoloCHK1( na,nblk,Q,TA,ALL_CONTXT,MPI_COMM_WORLD,INFO,As,WORK,LWORK,TAU )
   ttt1 = MPI_Wtime()
   IF( myid==0 ) write(*,*) 'Total time costs', ttt1-ttt0
   
   deallocate(Q)
   deallocate(As)
   deallocate(TA)
   deallocate(TAU)
   deallocate(WORK)

   call blacs_gridexit(all_contxt)
   call mpi_finalize(mpierr)
   call EXIT(STATUS)
!   
 end program test_zolodrv1

!-------------------------------------------------------------------------------
