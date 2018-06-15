!
  SUBROUTINE FormX2( N,NBLK,FONE,ERR,A,Z,CONTXT,WORK,INFO,IT,IP )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   INTEGER  :: N, NBLK,CONTXT, INFO, IT, IP
   REAL*8   :: FONE, ERR
!   
   REAL*8   :: A(*), Z(*),WORK(*)
!
! ============
!   
!  This subroutine forms X2 for next iteration and X1 is stored in A.
!  X2 = A + Z : = Z, Z would first store X2 and then copy back to A.  
!  Finally, A stores X2, and Z would be destroyed. 
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
! ========================
! Written by Shengguo Li, 12/21/2016, at Changsha
! NUDT, Changsha 410073, China
! ==================================================
!
   INTEGER  :: np_rows, np_cols, myprow, mypcol, na, na_rows, na_cols, &
        I, J, IA, JA, myid, nprocs
   DOUBLE PRECISION :: tt0, tt1
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0,MONE=-1.0D+00
!    
   INTEGER  :: DESCA(9)
!  ..   
!  .. External Functions ..
   INTEGER, EXTERNAL :: NUMROC
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

   ! Compute X1 from A, Z and FONE
   CALL PDGEADD( 'N',N,N,ONE,A,1,1,DESCA,ONE,Z,IA,JA,DESCA )  ! Z = A + Z
   CALL PDLASCL( 'General',FONE,ONE,N,N,Z,IA,JA,DESCA,info )  ! X2 = X2 / FONE
!    J = na_rows*na_cols
!    Z(1:J)=Z(1:J)+A(1:J)   
!    CALL DLASCL( 'G',IA,JA,FONE,ONE,na_rows,na_cols,Z,na_rows,info ) ! X2 = X2 / FONE

   ! Print out the matrix X2
!   IF(IT==0 .AND. IP==1) CALL PDLAPRNT(na,na,Z,IA,JA,DESCA,1,1,'X1',6,WORK ) 

    ! Compute X2-X1 and store it in A 
!    A(1:J) = Z(1:J)-A(1:J)
    CALL PDGEADD( 'N',N,N,MONE,Z,1,1,DESCA,ONE,A,IA,JA,DESCA )  ! A = A - Z
    ERR = PDLANGE( 'Fro',na,na,A,IA,JA,DESCA,WORK )

    ! Copy Z back to A, A = X2
    !A(1:J) = Z(1:J) 
    CALL PDLACPY( 'Full',N,N,Z,IA,JA,DESCA,A,IA,JA,DESCA )
   
 END SUBROUTINE FormX2
