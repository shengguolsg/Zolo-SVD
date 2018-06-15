!
  SUBROUTINE ZoloCHK1( N,NBLK,Q,Z,CONTXT,MPI_newcomm,INFO,As,WORK,LWORK,TAU )
!  
   IMPLICIT NONE
   include 'mpif.h'
!
   INTEGER  :: N, NBLK, CONTXT, INFO, LWORK, MPI_newcomm
   REAL*8   :: FONE
!   
   REAL*8   :: Q(*), Z(*), As(*), WORK(*), TAU(*)
!
! ============
!   
!  This subroutine checks the correctness of the eigendecomposition computed by ELPA. 
!
! ..Arguments..
! ==============
!
! N      INTEGER  The global column dimension of matrix A
!
! nblk   INTEGER  Blocking factor in block cyclic distribution  
!   
! FONE   DOUBLE PRECISION Scalar 
!        It is used to scale matrix (A+Z) to get X1 = (A+Z)/FONE. 
!
! Q      DOUBLE PRECISION Array. On entry, it stores the polar matrix Q.
!        On exit, it contains the left singular matrix. 
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
! WORK   DOUBLE PRECISION Array, Optional
!        It is used as workspace for accuracy test. 
!
! ========================
! Written by Shengguo Li, 01/19/2017, at Changsha
! NUDT, Changsha 410073, China
! ==================================================
!
   INTEGER  :: np_rows, np_cols, myprow, mypcol, na, na_rows, na_cols, &
        I, J, IA, JA, myid, nprocs, mpierr, LWORK2, LIWORK
   DOUBLE PRECISION :: err,normA, tt0, tt1, errmax, xc 
   DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0, HALF=0.5D+0, MONE=-1.0D+0
!    
   INTEGER  :: DESCA(9)
   DOUBLE PRECISION, ALLOCATABLE :: WORK2(:)
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

   Allocate( WORK2(na_rows*na_cols) ) 

   ! Set up a scalapack descriptor
   CALL DESCINIT( DESCA, na, na, nblk, nblk, 0, 0, CONTXT, na_rows, info )

   ! The following codes are testing the accuracy. Q stores Q and Z will store H.
   ! H = Q'*A=(Q'*As); H=(H'+H)/2
   tt0 = MPI_Wtime() 
   CALL PDGEMM( 'T','N', N,N,N,ONE,Q,IA,JA,DESCA,As,IA,JA,DESCA,&
        ZERO,Z,IA,JA,DESCA )
   CALL PDLACPY( 'Full',N,N,Z,IA,JA,DESCA,WORK,IA,JA,DESCA )
   CALL PDGEADD( 'T',N,N,HALF,WORK,IA,JA,DESCA,HALF,Z,IA,JA,DESCA )

   ! Call ELPA routines to compute the eigendecomposition of H
   CALL elpa_pdsyevd2( 'A',N,N,Z,na_rows,na_cols,TAU,Work2,nblk,DESCA,np_rows,&
                          np_cols,myprow,mypcol,myid,mpi_newcomm,info,12 ) !5,12
!   CALL elpa_pdsyevd('A',N,N,Z,na_rows,na_cols,TAU,Work2,nblk,DESCA,np_rows, &
!                      np_cols,myprow,mypcol,myid,mpi_newcomm,info )  
   CALL PDGEMM( 'N','N',N,N,N,ONE,Q,IA,JA,DESCA,WORK2,IA,JA,DESCA, & 
                 ZERO,WORK,IA,JA,DESCA ) 
   ! WORK stores the final left singular vecotrs; 
   ! The right ones are stored in Work2. 
   ! The singular values are stored in TAU. 
   tt1 = MPI_Wtime() 
   if(myid==0) print *,'Compute SVD costs    :', tt1 - tt0

   !-------------------------------------------------------------------------------
   ! Check the correctness of computed SVD (using plain scalapack routines)
   ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)

   ! tmp1 =  A * V
   call pdgemm('N','N',N,N,N,ONE,As,1,1,DESCA,WORK2,1,1,DESCA,ZERO,Z,1,1,DESCA )

   ! tmp2 = Zi*EVi
   CALL PDLACPY( 'Full',N,N,WORK,IA,JA,DESCA,As,IA,JA,DESCA )
   do i=1, N
      xc = TAU(i)
      call pdscal( N,xc,As,1,i,DESCA,1 )
   enddo

   !  tmp1 = A*Zi - Zi*EVi
   CALL PDGEADD( 'N',N,N,MONE,As,IA,JA,DESCA,ONE,Z,IA,JA,DESCA )

   ! Get maximum norm of columns of tmp1
!   errmax = 0
!   do i=1,N
!      xc = 0
!      call pddot( N,xc,Z,1,i,DESCA,1,Z,1,i,DESCA,1 )
!      errmax = max( errmax, sqrt(real(xc,8)) )
!   enddo

!   ! Get maximum error norm over all processors
!   normA = maxval(TAU(1:na)) 
!   err = errmax
!   CALL mpi_allreduce( err,errmax,1,MPI_REAL8,MPI_MAX,MPI_NEWCOMM,mpierr )
!   if(myid==0) print *
!   if(myid==0) print *,'Error Residual     :',errmax, errmax/(na*normA)
   
   ! Compute the Frobenius norm of the error matrix
   normA = maxval(TAU(1:na)) 
   err = PDLANGE( 'Fro',N,N,Z,1,1,DESCA,As )
   if(myid==0) print *, 'NormA is', normA
   if(myid==0) print *,'Fro Error Residual     :',err, err/(normA)
   
   !--------------------------------------------------------------
   ! Check the orthogonality of computed singular vectors

   ! tmp1 = U**T * U
   call pdgemm('T','N',na,na,na,ONE,WORK,1,1,DESCA, &
               WORK,1,1,DESCA,ZERO,Q,1,1,DESCA)
   ! Initialize Z to unit matrix
   call pdlaset('A',na,na,ZERO,ONE,Z,1,1,DESCA )

   ! Q = U**T * U - Unit Matrix, Q = Q -Z
   CALL PDGEADD( 'N',N,N,MONE,Z,IA,JA,DESCA,ONE,Q,IA,JA,DESCA )

   ! Get maximum error (max abs value in tmp1)
!   err = maxval(abs(Q(1:na_rows*na_cols)))
!   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_NEWCOMM,mpierr)
!   if(myid==0) print *,'Error Orthogonality of left SVecs:',errmax, &
!               errmax/(na), normA

   ! The Frobenius norm error  
   err = PDLANGE( 'Fro',N,N,Q,1,1,DESCA,Work )
   if(myid==0) print *,'Fro Error Orthogonality of left SVecs:',err, &
                  err/(na), normA

   ! Right singular Svecs  tmp1 = V**T * V
   call pdgemm('C','N',na,na,na,ONE,WORK2,1,1,DESCA, &
               WORK2,1,1,DESCA,ZERO,Q,1,1,DESCA)
   ! Initialize Z to unit matrix
   call pdlaset('A',na,na,ZERO,ONE,Z,1,1,DESCA )

   ! Q = U**T * U - Unit Matrix, Q = Q -Z
   CALL PDGEADD( 'N',N,N,MONE,Z,IA,JA,DESCA,ONE,Q,IA,JA,DESCA )

!   ! Get maximum error (max abs value in tmp1)
!   err = maxval(abs(Q(1:na_rows*na_cols)))
!   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_NEWCOMM,mpierr)
!   if(myid==0) print *,'Error Orthogonality of Right SVecs:',errmax, &
!               errmax/(na), normA

   ! The Frobenius norm error  
   err = PDLANGE( 'Fro',N,N,Q,1,1,DESCA,Work )
   if(myid==0) print *,'Fro Error Orthogonality of Right SVecs:',err, &
                  err/(na), normA

   DEALLOCATE( WORK2 )

 END SUBROUTINE ZoloCHK1
