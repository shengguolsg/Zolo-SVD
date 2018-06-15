!
  SUBROUTINE ZoloCHK( N,NBLK,A,Z,CONTXT,MPI_newcomm,INFO,As,WORK,LWORK,TAU )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   INTEGER  :: N, NBLK, CONTXT, INFO, LWORK, MPI_newcomm
   REAL*8   :: FONE
!   
   REAL*8   :: A(*), Z(*), As(*), WORK(*), TAU(*)
!
! ============
!   
!  This subroutine checks the correctness of Zolo-PD. 
!   
! ..Arguments..
! ================
!
! N      INTEGER  The global column dimension of matrix A
!
! nblk   INTEGER  Blocking factor in block cyclic distribution  
!   
! FONE   DOUBLE PRECISION Scalar 
!        It is used to scale matrix (A+Z) to get X1 = (A+Z)/FONE. 
!
! A      DOUBLE PRECISION Array. On entry, it stores the initial matrix X0.
!        On exit, it contains the polar orthogonal factor. 
!
! Z      DOUBLE PRECISION Array. On entry, it stores the summation of matrices gamma*Q_j1*Q_j2**T. 
!        On exit, it contains the result X2, the polar orthogonal factor. 
!
! CONTXT INTEGER The BLACS context.    
!
! COEFF  DOUBLE PRECISION Array
!        Coeff(1:2*rk) stores parameters cj, coeff(2*rk+1:3*rk) stores the parameters aj. 
!
! As     DOUBLE PRECISION Array, Optional
!        It stores the original matrix, used for accuracy test. 
!
! WORK  DOUBLE PRECISION Array, Optional
!       It is used as workspace for accuracy test. 
!
! ========================
! Written by Shengguo Li, 12/22/2016, at Changsha
! NUDT, Changsha 410073, China
! ==================================================
!
   INTEGER  :: np_rows, np_cols, myprow, mypcol, na, na_rows, na_cols, &
        I, J, IA, JA, myid, nprocs, &
        LWORK2, LIWORK
   DOUBLE PRECISION :: err,normA, tt0, tt1
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0, HALF=0.5D+0, MONE=-1.0D+0
!    
   INTEGER  :: DESCA(9)
!  ..   
!  .. External Functions ..
   INTEGER, EXTERNAL :: INDXL2G, NUMROC
   DOUBLE PRECISION, EXTERNAL :: PDLANGE
!  ..   
!  .. External Subroutine ..
   EXTERNAL :: PDLAPRNT
!
   CALL BLACS_PINFO( myid,nprocs )
   CALL BLACS_Gridinfo( CONTXT,np_rows,np_cols,myprow,mypcol )

   na = N
   IA = 1
   JA = 1
   na_rows = NUMROC( na, nblk, myprow, 0, np_rows )
   na_cols = NUMROC( na, nblk, mypcol, 0, np_cols )

!   write(*,*) 'myid=', myid, myprow, mypcol, np_rows, np_cols, IP

   ! Set up a scalapack descriptor
   CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )

   ! The following codes are testing the accuracy. A stores Q and Z will store H.
   ! H = Q'*A=(A'*As); H=(H'+H)/2
   CALL PDGEMM( 'T','N', N,N,N,ONE,A,IA,JA,DESCA,As,IA,JA,DESCA,&
        ZERO,Z,IA,JA,DESCA )
   CALL PDLACPY( 'Full',N,N,Z,IA,JA,DESCA,WORK,IA,JA,DESCA )
   CALL PDGEADD( 'T',N,N,HALF,WORK,IA,JA,DESCA,HALF,Z,IA,JA,DESCA )

   ! Check the accuracy of computed Polar Decomposition
   ! TA = A*Z; Err = As-TA
   CALL PDGEMM( 'N','N', N,N,N,ONE,A,IA,JA,DESCA,Z,IA,JA,DESCA,&
        ZERO,WORK,IA,JA,DESCA )
   CALL PDGEADD( 'N',N,N,MONE,As,IA,JA,DESCA,ONE,WORK,IA,JA,DESCA )
   err = PDLANGE( 'Fro',N,N,WORK,1,1,DESCA,Z )
   if(myid ==0 ) write(*,*) 'err= ', err

   if(myid ==0 ) write(*,*) 'Zolo-PD is finished'

   ! Call ELPA routines to compute the eigendecomposition of H
   CALL elpa_pdsyevd2( 'A',N,N,Z,na_rows,na_cols,TAU,As,nblk,DESCA,np_rows,&
                      np_cols,myprow,mypcol,myid,mpi_newcomm,info,5 ) ! 12 
!   CALL elpa_pdsyevd( 'A',N,N,Z,na_rows,na_cols,TAU,As,nblk,DESCA,np_rows, &
!                      np_cols,myprow,mypcol,myid,mpi_newcomm,info )  
   CALL PDGEMM( 'N','N',N,N,N,ONE,A,IA,JA,DESCA,Z,IA,JA,DESCA, & 
                 ZERO,WORK,IA,JA,DESCA ) 
   ! WORK stores the final left singular vecotrs; The right ones are stored in
   ! As. 
   
 END SUBROUTINE ZoloCHK
