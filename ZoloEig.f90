!
  SUBROUTINE ZoloEig( N,NBLK,Q,As,CONTXT,MPI_newcomm,WORK,LWORK,WORK2,LWORK2,TAU,INFO )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   INTEGER  :: N, NBLK, CONTXT, INFO, LWORK, LWORK2, MPI_newcomm
!   
   REAL*8   :: Q(*), WORK2(*), As(*), WORK(*), TAU(*)
!
! ============
!   
!  This subroutine computes the partial eigendecomposition of matrix As, which
!  is the original matrix. Matrix Q contains the polar decomposition of 
!  matrix As -sigma*I. It should further call 
!    a) pqdwhsubit.f90 to compute the invariant subspace V1; 
!    b) Compute the compressed matrix AA1 = V1^T * As * V1;
!    c) Compute the eigendecomposition of AA1 = VV1^T * Lambda * VV1;
!    d) Compute the eigenvectors of As, V1*VV1;
!   
! ..Arguments..
! =============
!
! N      INTEGER  The global column dimension of matrix A
!
! NBLK   INTEGER  Blocking factor in block cyclic distribution  
!   
! Q      DOUBLE PRECISION Array. On entry, it stores the polar orthogonal factor.
!
! As     DOUBLE PRECISION Array. On entry, it stores the original matrix.
!
! CONTXT INTEGER The BLACS context.    
!
! WORK   DOUBLE PRECISION Array, workspace
!        It should be large enough to contain at least TWO matrix A. 
!        On exit, it contains the computed matrix AA1. 
!
! LWORK  INTEGER The size of WORK 
!
! WORK2  DOUBLE PRECISION Array, WORKSPACE
!        It is used for PDGEQRF, and we let it be as large as matrix A.    
!
! ========================
! Written by Shengguo Li, 01/20/2017, at Changsha
! NUDT, Changsha 410073, China
! ==================================================
!
   INTEGER  :: np_rows, np_cols, myprow, mypcol, na, na_rows, na_cols, &
        I, J, IA, JA, myid, nprocs, LIWORK, Rk, LWORK1
   DOUBLE PRECISION :: err,normA, tt0, tt1
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0, HALF=0.5D+0, MONE=-1.0D+0
!    
   INTEGER :: DESCA(9)
   REAL*8  :: rrtmp(4)
   REAL*8  :: rtmp(4)
   INTEGER :: itmp(4)
   integer,allocatable :: iwork(:)   
!  ..   
!  ..External Functions..
   INTEGER, EXTERNAL :: INDXL2G, NUMROC
   DOUBLE PRECISION, EXTERNAL :: PDLANGE
!  ..   
!  ..External Subroutine..
   EXTERNAL :: PDLAPRNT
!
   CALL BLACS_PINFO( myid,nprocs )
   CALL BLACS_Gridinfo( CONTXT,np_rows,np_cols,myprow,mypcol )

   na = N
   IA = 1
   JA = 1
   na_rows = NUMROC( na, nblk, myprow, 0, np_rows )
   na_cols = NUMROC( na, nblk, mypcol, 0, np_cols )

   ! Set up a scalapack descriptor
   CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )

   NormA = PDLANGE('Fro',na,na,As,1,1,DESCA,WORK )

   CALL PQDWHSUBIT( Na,As,1,1,DESCA,NormA,Q,1,1,DESCA,WORK,LWORK,TAU,WORK2, &
        LWORK2,Rk,info )

   LIWORK = -1
   LWORK1 = -1
   rtmp(:) =ZERO
   rrtmp(:)=0.0d0
   itmp(:) =0

   CALL PDSYEVD( 'V','U',Na,As,1,1,DESCA,TAU,Q,1,1,DESCA,rtmp,LWORK1, &
        itmp,LIWORK,info )
   LIWORK=itmp(1)+1

!   if(myid ==0 ) print*, 'lwork, liwork ', lwork,liwork, int(rtmp)
   if(myid == 0 ) print*, 'Rk=',Rk 

   ALLOCATE( iwork(liwork) )

   CALL PDSYEVD( 'V','U',Rk,As,1,1,DESCA,TAU,WORK2,1,1,DESCA,WORK,LWORK, &
        IWORK,LIWORK,info )
   
   call pdgemm( 'N','N',na,Rk,Rk,ONE,Q,1,1,DESCA,WORK2,1,1,DESCA,ZERO,&
        WORK,1,1,DESCA )

 !  If(myid ==0 ) print *, 'eigens', TAU(1:Rk)
   
   !CALL PDLAPRNT( Na,Rk,WORK,1,1,DESCA,1,1,'QQ',6,Q )

   deallocate(iwork)   
   
 END SUBROUTINE ZoloEig
