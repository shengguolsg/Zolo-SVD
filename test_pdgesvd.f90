!
   program test_pdgesvd2

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - DOUBLE PRECISION version
!-------------------------------------------------------------------------------

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:   System size
   ! nev:  Number of eigenvectors to be calculated
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer :: nblk
   integer na, nev

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols, iam, lwork

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol,itmp(4)
   integer  nnz

   integer, external :: numroc

   real*8 err, errmax, ttt0, ttt1
   real*8  rrtmp(4)
   real*8, allocatable :: ev(:), vt(:,:)
   real*8 :: xc,rtmp(4)
   real*8, allocatable :: a(:,:), u(:,:), as(:,:), work(:)

   real*8, parameter :: ZERO = 0.d0, ONE = 1.d0
   

   integer :: iseed(4096)     ! Random seed, size should be sufficient for every generator

   integer :: STATUS
   logical :: write_to_file
   !-------------------------------------------------------------------------------
   !  Parse command line argumnents, if given

   write_to_file = .false.

   nblk = 64 
   na  = 19998 
   nev = na 
   !-------------------------------------------------------------------------------
   !  MPI Initialization
!   call mpi_init(mpierr)

!   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
!   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

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
      print '(a)','Standard eigenvalue problem - DOUBLE PRECISION version'
      print *
      print '(3(a,i0))','Matrix size=',na,', Number of eigenvectors=',nev,', Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !
   ! The BLACS context is only necessary for using Scalapack.
   !
   ! For ELPA, the MPI communicators along rows/cols are sufficient,
   ! and the grid setup may be done in an arbitrary way as long as it is
   ! consistent (i.e. 0<=my_prow<np_rows, 0<=my_pcol<np_cols and every
   ! process has a unique (my_prow,my_pcol) pair).

!  Initialize a single BLACS context

!   my_blacs_ctxt = mpi_comm_world
   call BLACS_GET(-1,0,my_blacs_ctxt)
   call BLACS_Gridinit( my_blacs_ctxt, 'R', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem

   allocate( a (na_rows,na_cols) )
   allocate( u (na_rows,na_cols) )
   allocate( as(na_rows,na_cols) )

   allocate(ev(na))

   ! For getting a hermitian test matrix A we get a random matrix Z
   ! and calculate A = Z + Z**H

   ! We want different random numbers on every process
   ! (otherways A might get rank deficient):

   iseed(:) = myid
   call RANDOM_SEED(put=iseed)

   allocate(vt(na_rows,na_cols))
   call RANDOM_NUMBER(vt)
   u(:,:) = vt(:,:)
   call RANDOM_NUMBER(vt)
   u(:,:) = u(:,:) + 1.d0*vt(:,:)

   a(:,:) = u(:,:)
   call pdtran(na, na, ONE, u, 1, 1, sc_desc, ONE, a, 1, 1, sc_desc) ! A = Z + Z**T

   ! Save original matrix A for later accuracy checks
   as = a

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors by calling pzheevd

   lwork=-1
   rtmp(:) =ZERO
   rrtmp(:)=0.0d0

   call pdgesvd('V','V',na,na,a,1,1,sc_desc,ev,u,1,1,sc_desc,vt,1,1,sc_desc,rtmp,lwork,&
        info)

   if( info .ne. 0 ) then
      print *, 'pzheevd returns error during 1st call', info
   end if
   lwork=2*int(rtmp(1))+1

   if(myid ==0 ) print*, 'lwork, ',lwork

   allocate( work(lwork) )

   ttt0 = MPI_Wtime()
   call pdgesvd('V','V',na,na,a,1,1,sc_desc,ev,u,1,1,sc_desc,vt,1,1,sc_desc,WORK,lwork, &
        info)
   ttt1 = MPI_Wtime()
   nnz = 0
   nnz = count(ev > ZERO)
   if(myid == 0) print *,'Total time (sum above)  :',ttt1-ttt0, 'nnz',nnz 
   

   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)

   ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)

   ! tmp1 =  A * v
!   call pdgemm('N','T',na,na,na,ONE,as,1,1,sc_desc, &
!                vt,1,1,sc_desc,ZERO,a,1,1,sc_desc)

   ! tmp2 = Zi*EVi
!   as(:,:) = u(:,:)
!   do i=1,nev
!      xc = ev(i)
!      call pdscal(na,xc,as,1,i,sc_desc,1)
!   enddo

   !  tmp1 = A*Zi - Zi*EVi
!   a(:,:) =  a(:,:) - as(:,:)

   ! Get maximum norm of columns of tmp1
!   errmax = 0
!   do i=1,nev
!      xc = 0
!      call pddot(na,xc,a,1,i,sc_desc,1,a,1,i,sc_desc,1)
!      errmax = max(errmax, sqrt(real(xc,8)))
!   enddo

   ! Get maximum error norm over all processors
!   err = errmax
!   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
!   if(myid==0) print *
!   if(myid==0) print *,'Error Residual     :',errmax
!

   deallocate(a)
   deallocate(as)
   deallocate(u)
   deallocate(vt)
   deallocate(ev)
   deallocate(work)

   call blacs_gridexit(my_blacs_ctxt)
!   call mpi_finalize(mpierr)
   call EXIT(STATUS)
end

!-------------------------------------------------------------------------------
