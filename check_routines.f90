!
module mod_check_routines

  interface check_routines
    module procedure CheckResidual_SVD
    module procedure CheckOrth
	module procedure CheckOrth1
  end interface

  contains

    SUBROUTINE CheckResidual_SVD( na,ma,myid,myrow,mycol,nprow,npcol, &
                       DESCA,As,U,VT,Tau )
!
      implicit none
!
      INTEGER,INTENT(IN)          :: myid,na,ma,myrow,mycol,nprow,npcol
      INTEGER,INTENT(IN)          :: DESCA(:)
      DOUBLE PRECISION,INTENT(IN) :: As(*),U(*),VT(*),Tau(*)
!
! ===============
!     This routine tests the residual of the computed SVD of As = U*Tau*VT. 
!     Afterward, As is unchanged. 
!     
! ========
!     As :: The original matrix with dimensions na-by-ma
!           Right now, we assume na == ma. 
!
!     U  :: The left singular vectors of computed
!
!     VT :: The right singular vectors of computed
!
!     Tau:: The computed singular values
! 
! ===================
!
! Written by Shengguo Li, on Jan. 18th, 2018. 
!
!====================================
!
      INTEGER ::   i,j,NB,na_rows,na_cols,N,IA,JA
      DOUBLE PRECISION  :: err,normA,xc
!
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+0,ONE=1.0D+0,&
             MONE=-1.0D+0
!
      DOUBLE PRECISION,ALLOCATABLE :: WORK1(:),WORK2(:) 
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      DOUBLE PRECISION, EXTERNAL :: PDLANGE
!     ..
!
      N  = na
      IA = 1
      JA = 1
      NB = DESCA( 6 )
      na_rows = numroc( na, NB, myrow, 0, nprow )
      na_cols = numroc( na, NB, mycol, 0, npcol )
   
      ALLOCATE( WORK1(na_rows*na_cols),WORK2(na_rows*na_cols) )
      ! WORK1 and WORK2 are used as workspace. 
      !
      !-------------------------------------------------------------------------------
      ! Check the correctness of computed SVD (using plain scalapack routines)
      ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)

      ! WORK1 =  A * V
      CALL PDGEMM('N','N',N,N,N,ONE,As,1,1,DESCA,VT,1,1,DESCA,ZERO,WORK1,1,1,DESCA )

      ! WORK2 = Ui*EVi
      CALL PDLACPY( 'Full',N,N,U,IA,JA,DESCA,WORK2,IA,JA,DESCA )
      DO i=1, N
         xc = TAU(i)
         call pdscal( N,xc,WORK2,1,i,DESCA,1 )
      ENDDO

      ! WORK2 = A*V - U*EVi
      CALL PDGEADD( 'N',N,N,MONE,WORK1,IA,JA,DESCA,ONE,WORK2,IA,JA,DESCA )

      ! Compute the Frobenius norm of the error matrix
      normA = maxval( TAU(1:na) )
      err = PDLANGE( 'Fro',N,N,WORK2,1,1,DESCA,WORK1 )
      if( myid==0 ) print *, 'NormA is', normA
      if( myid==0 ) print *,'Fro Error Residual     :',err, err/(normA)

      DEALLOCATE( WORK1,WORK2 )

    END SUBROUTINE CheckResidual_SVD


    SUBROUTINE CheckOrth( na,ma,myid,myrow,mycol,nprow,npcol, &
                       DESCA,A )
!
      implicit none
!
      INTEGER,INTENT(IN)          :: myid,na,ma,myrow,mycol,nprow,npcol
      INTEGER,INTENT(IN)          :: DESCA(*)
      DOUBLE PRECISION,INTENT(IN) :: A(*)
!
! ===============
!     This routine tests the orthogonality of A, and afterward, As is unchanged. 
!     
! ========
!     A :: The original matrix with dimensions na-by-ma
!          Right now we assume that na == ma. 
!
! ===================
!
! Written by Shengguo Li, on Jan. 18th, 2018. 
!
!====================================
!
      INTEGER ::   i,j,NB,na_rows,na_cols,N,IA,JA
      DOUBLE PRECISION :: err 
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+0,ONE=1.0D+0,&
             MONE=-1.0D+0
!
      DOUBLE PRECISION,ALLOCATABLE :: WORK1(:),WORK2(:) 
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      DOUBLE PRECISION, EXTERNAL :: PDLANGE
!     ..
!
      N  = na
      IA = 1
      JA = 1
      NB = DESCA( 6 )
      na_rows = NUMROC( na, NB, myrow, 0, nprow )
      na_cols = NUMROC( na, NB, mycol, 0, npcol )
   
      ALLOCATE( WORK1(na_rows*na_cols),WORK2(na_rows*na_cols) )
      ! WORK1 and WORK2 are used as workspace. 
      !
      !--------------------------------------------------------------
      ! Check the orthogonality of computed singular vectors

      ! tmp1 = U**T * U
      CALL PDGEMM('T','N',na,na,na,ONE,A,1,1,DESCA, &
               A,1,1,DESCA,ZERO,WORK1,1,1,DESCA)
      ! Initialize Z to unit matrix
      CALL PDLASET('A',na,na,ZERO,ONE,WORK2,1,1,DESCA )

      ! Q = U**T * U - Unit Matrix, Q = Q -Z
      CALL PDGEADD( 'N',N,N,MONE,WORK1,IA,JA,DESCA,ONE,WORK2,IA,JA,DESCA )

      ! The Frobenius norm error
      err = PDLANGE( 'Fro',N,N,WORK2,1,1,DESCA,WORK1 )
      if(myid==0) print *,'Fro Error Orthogonality of input matrix:',err, &
                  err/(na)

      DEALLOCATE( WORK1,WORK2 )

   END SUBROUTINE CheckOrth
   
   SUBROUTINE CheckOrth1( na,ma,nblk,A,CONTXT )
!
      implicit none
!
      INTEGER  :: na, ma, nblk
      INTEGER  :: CONTXT
      DOUBLE PRECISION,INTENT(IN) :: A(*)
!
! ===============
!     This routine tests the orthogonality of A, and afterward, As is unchanged. 
!     
! ========
!     A :: The original matrix with dimensions na-by-ma
!          Right now we assume that na == ma. 
!
! ===================
!
! Written by Shengguo Li, on Jan. 18th, 2018. 
!
!====================================
!
      INTEGER :: myid,myrow,mycol,nprow,npcol,nprocs
      INTEGER ::   i,j,NB,na_rows,na_cols,N,IA,JA,info
      DOUBLE PRECISION :: err 
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+0,ONE=1.0D+0,&
             MONE=-1.0D+0
!
      INTEGER  :: DESCA(9)
      DOUBLE PRECISION,ALLOCATABLE :: WORK1(:),WORK2(:) 
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      DOUBLE PRECISION, EXTERNAL :: PDLANGE
!     ..
!
      CALL BLACS_PINFO( myid,nprocs )
      CALL BLACS_Gridinfo( CONTXT,nprow,npcol,myrow,mycol )
	  
	  N  = na
      NB = nblk
      IA = 1
      JA = 1
      na_rows = NUMROC( na, NB, myrow, 0, nprow )
      na_cols = NUMROC( na, NB, mycol, 0, npcol )
	  
	  ! Set up a scalapack descriptor
	  CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )
	  
   
      ALLOCATE( WORK1(na_rows*na_cols),WORK2(na_rows*na_cols) )
      ! WORK1 and WORK2 are used as workspace. 
      !
      !--------------------------------------------------------------
      ! Check the orthogonality of computed singular vectors

      ! tmp1 = U**T * U
      CALL PDGEMM('T','N',na,na,na,ONE,A,1,1,DESCA, &
               A,1,1,DESCA,ZERO,WORK1,1,1,DESCA)
      ! Initialize Z to unit matrix
      CALL PDLASET('A',na,na,ZERO,ONE,WORK2,1,1,DESCA )

      ! Q = U**T * U - Unit Matrix, Q = Q -Z
      CALL PDGEADD( 'N',N,N,MONE,WORK1,IA,JA,DESCA,ONE,WORK2,IA,JA,DESCA )

      ! The Frobenius norm error
      err = PDLANGE( 'Fro',N,N,WORK2,1,1,DESCA,WORK1 )
      if(myid==0) print *,'Fro Error Orthogonality of input matrix:',err, &
                  err/(na)

      DEALLOCATE( WORK1,WORK2 )

   END SUBROUTINE CheckOrth1

end module
