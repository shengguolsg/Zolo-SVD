!
   program pqdwh_caller2

!-------------------------------------------------------------------------------
! This routine calls pqdwhfac.f90 to compute a polar factorization of an symmetric
! n x n real matrix A. It use the QDWH algorithm. 
!  
!-------------------------------------------------------------------------------
   use mod_prepare_matrix

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:  the size of matrix A
   ! nq:  the size of matrix Q, nq = 2*na  
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer :: nblk
   integer na, nq
!     .. Parameters ..
   INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                      LLD_, MB_, M_, NB_, N_, RSRC_
   PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                        CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                        RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols, nq_rows, nq_cols, nbw, nnz, &
           IA, JA

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, desca(9), descw(9),info, nprow, npcol, &
           lwork,iptau,ipwork2, lwork2, lworkt

   integer, external :: numroc
   character(len=128)   :: filename

   real*8  err, errmax, ttt0, ttt1, ttt2, sigma, xc
   real*8  rrtmp(4)
   real*8  rtmp(4)
   real*8, allocatable :: a(:,:), q(:,:), as(:,:), work(:)

   real*8, parameter :: ZERO = 0.d0, ONE = 1.d0, NEONE=-1.0D+0

   integer :: iseed(4096)     ! Random seed, size should be sufficient for every generator

   DOUBLE PRECISION, EXTERNAL :: PDLANGE

   integer :: STATUS
   !-------------------------------------------------------------------------------
   !  Parse command line argumnents, if given

   IA = 1
   JA = 1
   nblk = 64 
  
!  toep
!   na = 16
!   nnz = 31

!  linverse
!   na = 11999
!   nnz = 53988

!  fv1
!   na = 9604
!   nnz = 47434

!  nemeth03
!   na    = 9506 
!   nnz = 202157

!   bcsstk18
!   na  = 11948
!   nnz = 80519 

   nq  = 2*na
   sigma = zero
   !-------------------------------------------------------------------------------
   !  MPI Initialization
   call mpi_init(mpierr)

   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

   STATUS = 0

!  Initialize the BLACS: refer to user's guide or lbll manual
   call BLACS_PINFO(myid,nprocs)

   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs
   ! and decrement it until a divisor is found.

   do np_cols = NINT(SQRT(REAL(nprocs))),2,-1
      if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_rows = nprocs/np_cols

   if(myid==0) then
      print *
      print '(a)','Test polar factorization - DOUBLE PRECISION version'
      print *
      print '(2(a,i0))','Matrix size=',na, 'Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !

!  Initialize a single BLACS context

!   my_blacs_ctxt = mpi_comm_world
   call BLACS_GET(-1,0,my_blacs_ctxt)
   call BLACS_Gridinit( my_blacs_ctxt, 'R', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   if(myid == 0 ) then
      print '(2(a,i3))', 'np_rows=', np_rows, 'nprow=', nprow
      print '(2(a,i3))', 'np_cols=', np_cols, 'nprow=', npcol
   end if

   !print '(3(a,i3))', 'my_prow=', my_prow, 'mypcol=', my_pcol, 'myid=', myid

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
   nq_rows = numroc(nq, nblk, my_prow, 0, np_rows)
   nq_cols = numroc(na, nblk, my_pcol, 0, np_cols)

!   write(*,*) 'narows', na_rows, na_cols, my_prow, my_pcol, myid
!   write(*,*) 'nqrows', nq_rows, nq_cols, my_prow, my_pcol, myid

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( desca, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )
   call descinit( descw, nq, na, nblk, nblk, 0, 0, my_blacs_ctxt, nq_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem

   allocate( a (na_rows,na_cols) )
   allocate( q (na_rows,na_cols) )
   allocate( as(na_rows,na_cols) )

   ! For getting a hermitian test matrix A we get a random matrix Z
   ! and calculate A = Z + Z**H

   ! We want different random numbers on every process
   ! (otherways A might get rank deficient):

   !iseed(:) = myid
   !call RANDOM_SEED(put=iseed)

   !call RANDOM_NUMBER(q)

   !a(:,:) = q(:,:)
   !call pdtran(na, na, ONE, q, 1, 1, desca, ONE, a, 1, 1, desca) ! A = Z + Z**T

   !filename = 'linverse/linverse.mtx'
   !filename = 'bcsstk06.mtx'
   !filename = 'bcsstk18/bcsstk18.mtx'
   !filename = 'nemeth03/nemeth03.mtx'
   filename = 'fv1/fv1.mtx'
   !filename = 'toep.mtx'
   CALL PDLASET( 'A',na,na,ZERO,ZERO,A,1,1,desca )
   CALL PDLASET( 'A',na,na,ZERO,ZERO,q,1,1,desca )
   CALL prepare_real_symm_MM( filename,na,nnz,myid,my_prow,my_pcol, &
        np_rows,np_cols,desca,A,Q,sigma )

   ! Save original matrix A for later accuracy checks
   as = a

   !-------------------------------------------------------------------------------
   ! set up workspace and use pqdwh to compute the polar factor of A 
   ! The problem left is to estimate the size of WORK
   LWORK = -1
   CALL pqdwhfacs( na,na,zero,a,1,1,desca,q,1,1,desca,rtmp,lwork,descw,rrtmp,rtmp,lwork2,info,'N' )
   LWORK = int(rtmp(1)) + 1 +8*na 
   LWORK2 = int(rtmp(2))+1
   LWORKT= LWORK + LWORK2 + na

!   WRITE(*,*) 'The total workspace', lworkt
!   LWORK = nq_rows*nq_cols+4*na + na_rows*na_cols

!   write(*,*) 'nqrows', nq_rows, nq_cols, lwork, my_prow, my_pcol
!   write(*,*) 'narows', na_rows, na_cols, my_prow, my_pcol

   allocate( work(lworkt) )
   IPTAU = lwork +1 
   IPWORK2 = na + IPTAU

   !CALL PDLAPRNT( na,na,A,1,1,DESCA,1,1,'A0',6, WORK )

   ttt0 = MPI_Wtime()
   call pqdwhfacs(na,na,ZERO,a,1,1,desca,q,1,1,desca,WORK,lwork,descw,WORK(IPTAU),WORK(IPWORK2),&
        lwork2,info,'T' )
   ttt1 = MPI_Wtime()
   if(myid == 0) print *,'QDWH polar decomposition costs  :',ttt1-ttt0, 'info', info

   !CALL PDLAPRNT( na,na,Q,1,1,DESCA,1,1,'Q0',6,WORK )

   ! If the last parameter of pqdwhfacs equals to 'T', WORK(IPWORK2) stores the polar factor. 
   ! Compute the eigendecomposition of polar factor, i.e., the SVD.

   ! Call ELPA routines to compute the eigendecomposition of H
   CALL elpa_pdsyevd2( 'A',na,na,WORK(IPWORK2),na_rows,na_cols,WORK(IPTAU),Work,nblk,DESCA,na_rows, &
                      na_cols,my_prow,my_pcol,myid,mpi_comm_world,info,5 )
   CALL PDGEMM( 'N','N',na,na,na,ONE,Q,IA,JA,DESCA,WORK,IA,JA,DESCA, &
                 ZERO,WORK(IPWORK2),IA,JA,DESCA )
   ttt2 = MPI_Wtime()
   ! WORK(IPWORK2) stores the right singular vectors
   ! WORK stores the left singular vectors
   ! WORK(IPTAU) stores the singular values.

   if(myid == 0 ) print *, 'Total SVD costs :', ttt2-ttt0
  
   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)

   ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)
   ! tmp1 =  A * v
   call pdgemm('N','N',na,na,na,ONE,a,1,1,desca, &
                work(ipwork2),1,1,desca,ZERO,q,1,1,desca )

   ! tmp2 = Zi*EVi
   do i=1,na
      xc = work(iptau+i-1)
      call pdscal(na,xc,work,1,i,desca,1)
   enddo

   !  tmp1 = A*Zi - Zi*EVi
   CALL PDGEADD( 'N',na,na,NEONE,WORK,IA,JA,DESCA,ONE,Q,IA,JA,DESCA )

   ! Get maximum norm of columns of tmp1
   errmax = 0
   do i=1,na
      xc = 0
      call pddot(na,xc,Q,1,i,desca,1,Q,1,i,desca,1)
      errmax = max(errmax, sqrt(real(xc,8)))
   enddo

   ! Get maximum error norm over all processors
   err = errmax
   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
   if(myid==0) print *
   if(myid==0) print *,'Error Residual     :',errmax
   
   deallocate(a)
   deallocate(q)
   deallocate(work)
   deallocate(as)

   call blacs_gridexit(my_blacs_ctxt)
!   call mpi_finalize(mpierr)
   call EXIT(STATUS)

end

!-------------------------------------------------------------------------------
