!
!
! Written by Shengguo Li, National University of Defense Technology, Changsha 410073, China
! Some routines are modified from or obtained from ELPA. 
! If you have any questions, please tell me: nudtlsg@gmail.com. 
!
! ==================================================
! 
module mod_prepare_matrix

  interface prepare_matrix
    module procedure prepare_band_real
    module procedure prepare_band_complex
    module procedure prepare_real_symm_MM
    module procedure prepare_readMat
    module procedure prepare_real_GEM
    module procedure prepare_real_symm
    module procedure prepare_real_rand
    module procedure mmread
    module procedure mminfo
    module procedure mmwrite
    module procedure lowerc
    module procedure getwd
    module procedure countwd
  end interface

  contains

    subroutine prepare_band_real( na,bw,myid,myrow,mycol,nprow,npcol, &
                       sc_desc,iseed,a,z,as )
!
      implicit none

      integer, intent(in)       :: myid, na, sc_desc(:),myrow,mycol,nprow, &
                                   npcol, bw
      integer, intent(inout)    :: iseed(:)
      real*8, intent(inout)     :: z(:,:), a(:,:), as(:,:)
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
      integer i, j, nii, njj, NB, na_rows, na_cols
      double precision, parameter :: ZERO = 0.0D+0
!
      ! for getting a hermitian test matrix A we get a random matrix Z
      ! and calculate A = Z + Z**H

      ! we want different random numbers on every process
      ! (otherwise A might get rank deficient):

      NB = sc_desc( 6 )
      iseed(:) = myid
      call RANDOM_SEED(put=iseed)
      call RANDOM_NUMBER(Z)

      a(:,:) = z(:,:)

      if (myid == 0) then
        print '(a)','| Random matrix block has been set up. (only processor 0 confirms this step)'
        write(*,*) "NB=", NB, bw, na, nprow, npcol
      endif

      call pdtran(na, na, 1.d0, z, 1, 1, sc_desc, 1.d0, a, 1, 1, sc_desc) ! A = Z + Z**T

      if (myid == 0) then
        print '(a)','| Random matrix block has been symmetrized'
      endif

       na_rows = numroc( na, NB, myrow, 0, nprow )
       na_cols = numroc( na, NB, mycol, 0, npcol )

!      if(myid .lt. nprow*npcol) then
        DO I = 1, na_rows
           DO J = 1, na_cols
     	      nii = indxl2g( i, nb, myrow, 0, nprow  )
	      njj = indxl2g( j, nb, mycol, 0, npcol  )
	      if ( abs(nii-njj) .gt. bw ) then
	         A(i,j) = ZERO
	      end if
           END DO
        END DO
 !     end if

      ! save original matrix A for later accuracy checks
      as = a

    end subroutine prepare_band_real

    subroutine prepare_band_complex(na,bw,myid,myrow,mycol,nprow,npcol, &
    	       sc_desc, iseed, xr, a, z, as)

      implicit none

      integer, intent(in)       :: myid, na, sc_desc(:),myrow,mycol,nprow, &
      	       			   npcol, bw
      integer, intent(inout)    :: iseed(:)
      real*8, intent(inout)     :: xr(:,:)
      complex*16, intent(inout) :: z(:,:), a(:,:), as(:,:)
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
      complex*16, parameter     :: CZERO = (0.d0, 0.d0), CONE = (1.d0, 0.d0)
      integer i, j, nii, njj, NB, na_rows, na_cols

      ! for getting a hermitian test matrix A we get a random matrix Z
      ! and calculate A = Z + Z**H

      ! we want different random numbers on every process
      ! (otherwise A might get rank deficient):

      NB = sc_desc( 6 )
      iseed(:) = myid
      call RANDOM_SEED(put=iseed)
      call RANDOM_NUMBER(xr)
      z(:,:) = xr(:,:)
      call RANDOM_NUMBER(xr)
      z(:,:) = z(:,:) + (0.d0,1.d0)*xr(:,:)

      a(:,:) = z(:,:)

      if (myid == 0) then
        print '(a)','| Random matrix block has been set up. (only processor 0 confirms this step)'
      endif

      call pztranc(na, na, CONE, z, 1, 1, sc_desc, CONE, a, 1, 1, sc_desc) ! A = Z + Z**H

      if (myid == 0) then
        print '(a)','| Random matrix block has been symmetrized'
      endif

       na_rows = numroc( na, NB, myrow, 0, nprow )
       na_cols = numroc( na, NB, mycol, 0, npcol )

!      if(myid .lt. nprow*npcol) then
        DO I = 1, na_rows
           DO J = 1, na_cols
     	      nii = indxl2g( i, nb, myrow, 0, nprow  )
	      njj = indxl2g( j, nb, mycol, 0, npcol  )
	      if ( abs(nii-njj) .gt. bw ) then
	         A(i,j) = CZERO
	      end if
           END DO
        END DO
 !     end if

      ! save original matrix A for later accuracy checks
      as = a

    end subroutine prepare_band_complex

    subroutine prepare_real_GEM( SYM,filename,na,nnz,myid,myrow,mycol,nprow,npcol, &
                       sc_desc,a,z,sigma )
!
      implicit none

      CHARACTER                 :: SYM
      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(inout)    :: nnz
      character(len=128)        :: filename
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: sigma
      real*8, intent(inout)     :: z(*), a(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================

      integer nntmp, nnzmax, bw, mtype
      character rep*10
      character field*7
      character symm*19
      parameter (nntmp=3 )
      integer ival(nntmp)
      complex cval(nntmp)

      integer, allocatable :: indx(:), jndx(:), tindx(:)
      real*8, allocatable  :: rval(:)
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D, INFOG2L, PDLACPY, PDTRAN
      
      integer  i, j, nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      LOGICAL, EXTERNAL :: LSAME
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
      mtype = -1
      nnzmax = nnz
      sigma = zero

      ALLOCATE( indx(nnz),jndx(nnz),rval(nnz),tindx(nnz) )
!
      if(myid ==0 ) then
         print *,'Reading header and data...'
         open( unit=iunit, file=filename )
         call mmread(iunit,rep,field,symm,nrows,ncols,nnz,nnzmax, &
              indx,jndx,ival,rval,cval )
         tindx = indx(1:nnz)-jndx(1:nnz)
         bw = MAXVAL( ABS(tindx(1:nnz) ) ) +1
         print *,'  Matrix is type: ',rep,' ',field,' ',symm
         print *,'  Matrix size: ',nrows,' by ',ncols,' with ', &
              nnz,' nonzeros', 'bandwidth', bw, nnzmax
         close(iunit)

         IF( LSAME(symm,'symmetric') )  THEN
            mtype = 0
            write(*,*) 'symm in prepare_matrix is ', symm
         ELSE
           IF( LSAME(symm,'general') ) THEN
              write(*,*) 'symm in prepare_matrix is ', symm
              mtype = 1
           ELSE
              write(*,*) 'symm in prepare_matrix is not defined'
              mtype = 2
           END IF
         END IF

!#if defined(OUTPUT_Mat)
!        Write(*,*) "The values of matrix A" 
!        DO I =1, NNZ, 1
!          Write(*,*) rval(I)  
!        END DO 
!#endif

         ! broadcast indx, jndx and rval
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, indx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, jndx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, rval, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, nnz, 1 )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, mtype, 1 )
      else
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,indx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,jndx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,rval,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,nnz,1,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,mtype,1,IQROW,IQCOL )
      endif

      narows = numroc( na, NB, myrow, 0, nprow )
      nacols = numroc( na, NB, mycol, 0, npcol )

      ! Use indx, jndx and rval to inital matrix A
      DO I = 1, NNZ, 1
         irow = indx( I )
         jcol = jndx( I )
         CALL INFOG2L( irow,jcol,sc_desc,nprow,npcol,myrow,mycol,IID,JJD,IDROW,IDCOL )
          IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
             IPQ = IID + ( JJD-1 )*narows
             A( IPQ ) = rval( I )
          END IF
!         if(irow .eq. jcol) sigma = sigma+rval(I)
      END DO

      IF( mtype .eq. 0 )  THEN

         if (myid == 0) then
           print '(a)','| Symmetric sparse matrix starts. '
         endif
   
         SYM = 'S'
         call pdlacpy( 'Full',na,na,A,1,1,sc_desc,Z,1,1,sc_desc )
         call pdtran( na,na,ONE,Z,1,1,sc_desc,ONE,A,1,1,sc_desc ) ! A = Z+Z**T

         DO I = 1, narows
            DO J = 1, nacols
               nii = INDXL2G( i, nb, myrow, 0, nprow  )
               njj = INDXL2G( j, nb, mycol, 0, npcol  )
               if ( nii == njj ) then
                  IPQ = I +(J-1)*narows
                  A(IPQ) = A(IPQ) / TWO - SIGMA
               end if
            END DO
         END DO

         if (myid == 0) then
           print '(a)','| Symmetric sparse matrix is initialized. '
         endif

      ELSE
        IF( mtype .eq. 1 ) THEN
          SYM = 'Gen'
          if (myid == 0) then
            print '(a)','| General sparse matrix is initialized. '
          endif
        ELSE
          SYM = 'Other'
          if (myid == 0) then
            print '(a)','| The matrix type is not symm either general.'
          endif
        END IF 

      END IF

     deallocate( indx,jndx,tindx,rval ) 

   end subroutine prepare_real_GEM

    subroutine prepare_real_symm(filename,na,nnz,indx,jndx,rval,myid,myrow,mycol,nprow,npcol, &
                       sc_desc,a,z,sigma )
!
      implicit none

      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(inout)    :: nnz,indx(*),jndx(*)
      character(len=128)        :: filename
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: sigma
      real*8, intent(inout)     :: z(*), a(*), rval(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================

      integer nntmp, nnzmax, bw
      character rep*10
      character field*7
      character symm*19
      parameter (nntmp=3 )
      integer ival(nntmp)
      complex cval(nntmp)

      integer, allocatable :: tindx(:)
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D, INFOG2L, PDLACPY, PDTRAN
      
      integer  i, j, nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
      nnzmax = nnz
      sigma = zero

      ALLOCATE( tindx(nnz) )
!
      if(myid ==0 ) then
         print *,'Reading header and data...'
         open( unit=iunit, file=filename )
         call mmread(iunit,rep,field,symm,nrows,ncols,nnz,nnzmax, &
              indx,jndx,ival,rval,cval )
         tindx = indx(1:nnz)-jndx(1:nnz)
         bw = MAXVAL( ABS(tindx(1:nnz) ) ) +1
         print *,'  Matrix is type: ',rep,' ',field,' ',symm
         print *,'  Matrix size: ',nrows,' by ',ncols,' with ', &
              nnz,' nonzeros', 'bandwidth', bw
         !close(iunit)

         ! broadcast indx, jndx and rval
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, indx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, jndx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, rval, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, nnz, 1 )
      else
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,indx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,jndx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,rval,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,nnz,1,IQROW,IQCOL )
      endif

      narows = numroc( na, NB, myrow, 0, nprow )
      nacols = numroc( na, NB, mycol, 0, npcol )

      if(myid ==0 ) then
        write(*,*) 'Broadcast finishes.'
      endif

      ! Use indx, jndx and rval to inital matrix A
      DO I = 1, NNZ, 1
         irow = indx( I )
         jcol = jndx( I )
         CALL INFOG2L( irow,jcol,sc_desc,nprow,npcol,myrow,mycol,IID,JJD,IDROW,IDCOL )
          IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
             IPQ = IID + ( JJD-1 )*narows
             A( IPQ ) = rval( I )
          END IF
!         if(irow .eq. jcol) sigma = sigma+rval(I)
      END DO
!      sigma = sigma / (two*na) 

      if(myid ==0 ) then
        write(*,*) 'Start to copy'
      endif
      call pdlacpy( 'Full',na,na,A,1,1,sc_desc,Z,1,1,sc_desc )
      call pdtran(na, na, ONE, z, 1, 1, sc_desc, ONE, a, 1, 1, sc_desc) ! A = Z+Z**T

      DO I = 1, narows
         DO J = 1, nacols
            nii = INDXL2G( i, nb, myrow, 0, nprow  )
            njj = INDXL2G( j, nb, mycol, 0, npcol  )
            if ( nii == njj ) then
               IPQ = I +(J-1)*narows
               A(IPQ) = A(IPQ) / TWO - SIGMA
            end if
         END DO
      END DO

      if (myid == 0) then
        print '(a)','| MM sparse matrix is initialized. '
     endif

     deallocate( tindx ) 

   end subroutine prepare_real_symm

    subroutine prepare_real_rand( SYM,na,myid,myrow,mycol,nprow,npcol,sc_desc,a,z )
!
      implicit none

      character                 :: SYM
      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: z(*), a(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!      
      integer :: iseed(4096)     ! Random seed, size should be sufficient for every generator      
!
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D, INFOG2L, PDLACPY, PDTRAN
      LOGICAL, EXTERNAL :: LSAME
      
      integer  j, nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
!
      iseed(:) = myid
      call RANDOM_SEED(put=iseed)
      
      narows = numroc( na, NB, myrow, 0, nprow )
      nacols = numroc( na, NB, mycol, 0, npcol )

      call RANDOM_NUMBER( A( 1:narows*nacols ) )

      call pdlacpy( 'Full',na,na,A,1,1,sc_desc,Z,1,1,sc_desc )

      IF ( LSAME(SYM,'S') ) THEN
         call pdtran(na, na, ONE, z, 1, 1, sc_desc, ONE, a, 1, 1, sc_desc) ! A = Z+Z**T
      ENDIF

      if (myid == 0) then
        print '(a)','| Randomized matrix is initialized. '
     endif

   end subroutine prepare_real_rand

   
    subroutine prepare_real_nonsym(filename,na,nnz,indx,jndx,rval,myid,myrow,mycol,nprow,npcol, &
                       sc_desc,a,z,sigma )
!
      implicit none

      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(inout)    :: nnz,indx(*),jndx(*)
      character(len=128)        :: filename
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: sigma
      real*8, intent(inout)     :: z(*), a(*), rval(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================

      integer nntmp, nnzmax, bw
      character rep*10
      character field*7
      character symm*19
      parameter (nntmp=3 )
      integer ival(nntmp)
      complex cval(nntmp)

      integer, allocatable :: tindx(:)
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D, INFOG2L, PDLACPY, PDTRAN
      
      integer  i, nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
      nnzmax = nnz
      sigma = zero

      ALLOCATE( tindx(nnz) )
!
      if(myid ==0 ) then
         print *,'Reading header and data...'
         open( unit=iunit, file=filename )
         call mmread(iunit,rep,field,symm,nrows,ncols,nnz,nnzmax, &
              indx,jndx,ival,rval,cval )
         tindx = indx(1:nnz)-jndx(1:nnz)
         bw = MAXVAL( ABS(tindx(1:nnz) ) ) +1
         print *,'  Matrix is type: ',rep,' ',field,' ',symm
         print *,'  Matrix size: ',nrows,' by ',ncols,' with ', &
              nnz,' nonzeros', 'bandwidth', bw
         close(iunit)

         ! broadcast indx, jndx and rval
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, indx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, jndx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, rval, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, nnz, 1 )
      else
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,indx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,jndx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,rval,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,nnz,1,IQROW,IQCOL )
      endif

      narows = numroc( na, NB, myrow, 0, nprow )
      nacols = numroc( na, NB, mycol, 0, npcol )

      ! Use indx, jndx and rval to inital matrix A
      DO I = 1, NNZ, 1
         irow = indx( I )
         jcol = jndx( I )
         CALL INFOG2L( irow,jcol,sc_desc,nprow,npcol,myrow,mycol,IID,JJD,IDROW,IDCOL )
          IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
             IPQ = IID + ( JJD-1 )*narows
             A( IPQ ) = rval( I )
          END IF
!         if(irow .eq. jcol) sigma = sigma+rval(I)
      END DO
!      sigma = sigma / (two*na) 

!      DO I = 1, narows
!         DO J = 1, nacols
!            nii = INDXL2G( i, nb, myrow, 0, nprow  )
!            njj = INDXL2G( j, nb, mycol, 0, npcol  )
!            if ( nii == njj ) then
!               IPQ = I +(J-1)*narows
!               A(IPQ) = A(IPQ) - SIGMA
!            end if
!         END DO
!      END DO

      if (myid == 0) then
        print '(a)','| MM sparse matrix is initialized. '
     endif

     deallocate( tindx ) 

   end subroutine prepare_real_nonsym

    subroutine prepare_real_symm_MM(filename,na,nnz,myid,myrow,mycol,nprow,npcol, &
                       sc_desc,a,z,sigma )
!
      implicit none

      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(inout)    :: nnz
      character(len=128)        :: filename
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: sigma
      real*8, intent(inout)     :: z(*), a(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================

      integer nntmp, nnzmax, bw
      character rep*10
      character field*7
      character symm*19
      parameter (nntmp=3 )
      integer ival(nntmp)
      complex cval(nntmp)

      integer, allocatable :: indx(:), jndx(:), tindx(:)
      real*8, allocatable  :: rval(:)
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!
!     ..
!     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D, INFOG2L, PDLACPY, PDTRAN
      
      integer  i, j, nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
      nnzmax = nnz
      sigma = zero

      ALLOCATE( indx(nnz), jndx(nnz), rval(nnz),tindx(nnz) )
!
      if(myid ==0 ) then
         print *,'Reading header and data...'
         open( unit=iunit, file=filename )
         call mmread(iunit,rep,field,symm,nrows,ncols,nnz,nnzmax, &
              indx,jndx,ival,rval,cval )
         tindx = indx(1:nnz)-jndx(1:nnz)
         bw = MAXVAL( ABS(tindx(1:nnz) ) ) +1
         print *,'  Matrix is type: ',rep,' ',field,' ',symm
         print *,'  Matrix size: ',nrows,' by ',ncols,' with ', &
              nnz,' nonzeros', 'bandwidth', bw
         close(iunit)

         ! broadcast indx, jndx and rval
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, indx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, jndx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, rval, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, nnz, 1 )
      else
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,indx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,jndx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,rval,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,nnz,1,IQROW,IQCOL )
      endif

      narows = numroc( na, NB, myrow, 0, nprow )
      nacols = numroc( na, NB, mycol, 0, npcol )

      ! Use indx, jndx and rval to inital matrix A
      DO I = 1, NNZ, 1
         irow = indx( I )
         jcol = jndx( I )
         CALL INFOG2L( irow,jcol,sc_desc,nprow,npcol,myrow,mycol,IID,JJD,IDROW,IDCOL )
          IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
             IPQ = IID + ( JJD-1 )*narows
             A( IPQ ) = rval( I )
          END IF
!         if(irow .eq. jcol) sigma = sigma+rval(I)
      END DO
!      sigma = sigma / (two*na) 

      call pdlacpy( 'Full',na,na,A,1,1,sc_desc,Z,1,1,sc_desc )
      call pdtran(na, na, ONE, z, 1, 1, sc_desc, ONE, a, 1, 1, sc_desc) ! A = Z+Z**T

      DO I = 1, narows
         DO J = 1, nacols
            nii = INDXL2G( i, nb, myrow, 0, nprow  )
            njj = INDXL2G( j, nb, mycol, 0, npcol  )
            if ( nii == njj ) then
               IPQ = I +(J-1)*narows
               A(IPQ) = A(IPQ) / TWO - SIGMA
            end if
         END DO
      END DO

      if (myid == 0) then
        print '(a)','| MM sparse matrix is initialized. '
     endif

     deallocate( indx,jndx,tindx,rval ) 

   end subroutine prepare_real_symm_MM

    subroutine prepare_readMat(filename,na,nnz,indx,jndx,rval,myid,myrow,mycol,nprow,npcol, &
                       sc_desc )
!
      implicit none

      integer, intent(in)       :: myid,myrow,mycol,nprow,npcol,na
      integer, intent(inout)    :: nnz,indx(*),jndx(*)
      character(len=128)        :: filename
      integer, intent(in)       :: sc_desc(*)
      real*8, intent(inout)     :: rval(*)
!
!     This file reads a sparse matrix in MM format to A, which is assumed to be real
!     and symmetric. Z is used as workspace to obtain an symmetric matrix. Both lower
!     upper parts of a symmetrix matrix are stored in ELPA, different from SCALAPACK or
!     LAPACK routines. A is calculated A = 0.5*A + 0.5*Z**T, and A=Z. 
!
! =================
!  Written by Shengguo Li, 08/11/2016 
!
! =============================================

      integer nntmp, nnzmax, bw
      character rep*10
      character field*7
      character symm*19
      parameter (nntmp=3 )
      integer ival(nntmp)
      complex cval(nntmp)
!      
      integer  :: nrows,ncols,iunit,maxbw,irow,jcol,IID,JJD,IDROW,IDCOL, &
           narows, nacols, IPQ
      double precision, parameter :: ONE = 1.0D+0, TWO=2.0D+0
      integer, parameter :: CTXT_ = 2
!
!     ..
!     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
!     ..
!     .. External Subroutines ..      
      EXTERNAL           DGEBS2D, DGEBR2D
      
      integer  nii, njj, NB, IQROW, IQCOL
      double precision, parameter :: ZERO = 0.0D+0
!
!      
      INTRINSIC  :: MAXVAL, ABS, CEILING
      
      NB = sc_desc( 6 )
      iunit = 8
      IQROW = 0
      IQCOL = 0
      nnzmax = nnz

      if(myid ==0 ) then
         print *,'Reading header and data...'
         open( unit=iunit, file=filename )
         call mmread(iunit,rep,field,symm,nrows,ncols,nnz,nnzmax, &
              indx,jndx,ival,rval,cval )
         print *,'  Matrix is type: ',rep,' ',field,' ',symm
         print *,'  Matrix size: ',nrows,' by ',ncols,' with ', &
              nnz,' nonzeros'
         close(iunit)

         ! broadcast indx, jndx and rval
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, indx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, jndx, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax, 1, rval, nnzmax )
         CALL DGEBS2D( sc_desc( CTXT_ ), 'A', ' ', 1, 1, nnz, 1 )
      else
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,indx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,jndx,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', nnzmax,1,rval,nnzmax,IQROW,IQCOL )
         CALL DGEBR2D( sc_desc( CTXT_ ), 'A', ' ', 1,1,nnz,1,IQROW,IQCOL )
      endif

      if (myid == 0) then
        print '(a)','| MM sparse matrix is initialized. '
     endif

   end subroutine prepare_readMat

!  The following subroutine reads the number of rows, columns, and nnzeros of a MTX matrix. 

      subroutine mmread_nnz(iunit,rows,cols,nnz )
!        
        implicit none
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will read data from a matrix market formatted file.
! It is derivated from mmread, and only reads the number of rows, columns, and nnz.  
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   iunit    integer     in   Unit identifier for the file
!                             containing the data to be read.
!                             Must be open prior to call.
!                             Will be rewound on return.
!
!   rep     character*10 out  Matrix Market 'representation'
!                             indicator. On return:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                                elemental    (to be added)
!
!   field   character*7  out  Matrix Market 'field'. On return:
!
!                                real
!                                complex
!                                integer
!                                pattern
!
!   symm    character*19 out  Matrix Market 'field'. On return:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     out  Number of rows in matrix.
!
!   cols     integer     out  Number of columns in matrix.
!
!   nnz      integer     out  Number of nonzero entries required to
!                             store matrix.
!
!   nnzmax   integer     in   Maximum dimension of data arrays.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer i, rows, cols, nnz, nnzreq, nnzmax, iunit
      integer count,next
      character mmhead*15
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
!
! Read header line and check validity:
!
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
!
! Convert type code to lower case for easier comparisons:
!
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
         stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' ) &
         go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' .and. &
          field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and. &
          symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric') &
         go to 9000
!
! Read through comment lines, ignoring content:
!
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
!
! Just read a non-comment.
!   Now, back up a line, and read for first int, and back up
!   again. This will set pointer to just before apparent size
!   info line.
!   Before continuing with free form input, count the number of
!   words on the size info line to ensure there is the right amount
!   of info (2 words for array matrices, 3 for coordinate matrices).
!
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
!
!   Correct number of words are present, now back up and read them.
!
      backspace (iunit)
!
      if ( rep .eq. 'coordinate' ) then
!
! Read matrix in sparse coordinate format
!
        read (iunit,fmt=*) rows,cols,nnz
        return
!
      elseif ( rep .eq. 'array' ) then
!
! Read matrix in dense column-oriented array format
!
        read (iunit,fmt=*) rows,cols
!
! Check to ensure adequate storage is available
!
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnzreq = (rows*cols - rows)/2 + rows
          nnz = nnzreq
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnzreq = (rows*cols - rows)/2
          nnz = nnzreq
        else
          nnzreq = rows*cols
          nnz = nnzreq
        endif
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
!
! Various error conditions:
!
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data lines found.'
      stop
 3000 print *,'Size info inconsistant with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistant with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 4000 print *,'Premature end-of-file.'
      print *,'Check that the data file contains ',nnz, &
              ' lines of  i,j,[val] data.'
      print *,'(it appears there are only ',i,' such lines.)'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    end subroutine mmread_nnz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmread_nnz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!  The following files are from Matrix-Market website, which can read and write MM format
!  sparse matrices.    

      subroutine mmread(iunit,rep,field,symm,rows,cols,nnz,nnzmax, &
           indx,jndx,ival,rval,cval)
!        
        implicit none
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will read data from a matrix market formatted file.
! The data may be either sparse coordinate format, or dense array format.
!
! The unit iunit must be open, and the file will be rewound on return.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
!             Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
! 15-Oct-08  fixed illegal attempt of mimicking "do while" construct
!            by redifing limits inside loop. (lines 443-450)
!            (Thanks to Geraldo Veiga for his comments.)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   iunit    integer     in   Unit identifier for the file
!                             containing the data to be read.
!                             Must be open prior to call.
!                             Will be rewound on return.
!
!   rep     character*10 out  Matrix Market 'representation'
!                             indicator. On return:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                                elemental    (to be added)
!
!   field   character*7  out  Matrix Market 'field'. On return:
!
!                                real
!                                complex
!                                integer
!                                pattern
!
!   symm    character*19 out  Matrix Market 'field'. On return:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     out  Number of rows in matrix.
!
!   cols     integer     out  Number of columns in matrix.
!
!   nnz      integer     out  Number of nonzero entries required to
!                             store matrix.
!
!   nnzmax   integer     in   Maximum dimension of data arrays.
!
!   indx     integer(nnz)out  Row indices for coordinate format.
!                             Undefined for array format.
!
!   jndx     integer(nnz)out  Column indices for coordinate format.
!                             Undefined for array format.
!
!   ival     integer(nnz) out Integer data (if applicable, see 'field')
!
!   rval     double(nnz) out  Real data (if applicable, see 'field')
!
!   cval     complex(nnz)out  Complex data (if applicable, see 'field')
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      double precision rpart,ipart
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, nnzmax, iunit
      integer count,next
      character mmhead*15
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
!
! Read header line and check validity:
!
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
!
! Convert type code to lower case for easier comparisons:
!
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
         stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' ) &
         go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' .and. &
          field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and. &
          symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric') &
         go to 9000
!
! Read through comment lines, ignoring content:
!
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
!
! Just read a non-comment.
!   Now, back up a line, and read for first int, and back up
!   again. This will set pointer to just before apparent size
!   info line.
!   Before continuing with free form input, count the number of
!   words on the size info line to ensure there is the right amount
!   of info (2 words for array matrices, 3 for coordinate matrices).
!
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
!
!   Correct number of words are present, now back up and read them.
!
      backspace (iunit)
!
      if ( rep .eq. 'coordinate' ) then
!
! Read matrix in sparse coordinate format
!
        read (iunit,fmt=*) rows,cols,nnz
!
! Check to ensure adequate storage is available
!
        if ( nnz .gt. nnzmax ) then
          print *,'insufficent array lengths for matrix of ',nnz, &
                  ' nonzeros.'
          print *,'resize nnzmax to at least ',nnz,'. (currently ', &
                  nnzmax,')'
          stop
        endif
!
! Read data according to data type (real,integer,complex, or pattern)
!
        if ( field .eq. 'integer' ) then
          do 30 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),ival(i)
 30       continue
        elseif ( field .eq. 'real' ) then
          do 35 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rval(i)
 35       continue
        elseif ( field .eq. 'complex' ) then
          do 40 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i),rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 40       continue
        elseif ( field .eq. 'pattern' ) then
          do 50 i=1,nnz
            read (iunit,fmt=*,end=4000) indx(i),jndx(i)
            rval(i) = 1.0E0
 50       continue
        else
           print *,'''',field,''' data type not recognized.'
           stop
        endif
        rewind(iunit)
        return
!
      elseif ( rep .eq. 'array' ) then
!
! Read matrix in dense column-oriented array format
!
        read (iunit,fmt=*) rows,cols
!
! Check to ensure adequate storage is available
!
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnzreq = (rows*cols - rows)/2 + rows
          nnz = nnzreq
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnzreq = (rows*cols - rows)/2
          nnz = nnzreq
        else
          nnzreq = rows*cols
          nnz = nnzreq
        endif
        if ( nnzreq .gt. nnzmax ) then
          print *,'insufficent array length for ',rows, ' by ', &
                   cols,' dense ',symm,' matrix.'
          print *,'resize nnzmax to at least ',nnzreq,'. (currently ', &
                   nnzmax,')'
          stop
        endif
!
! Read data according to data type (real,integer,complex, or pattern)
!
        if ( field .eq. 'integer' ) then
          do 60 i=1,nnzreq
            read (iunit,fmt=*,end=4000) ival(i)
 60      continue
        elseif ( field .eq. 'real' ) then
          do 65 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rval(i)
 65      continue
        elseif ( field .eq. 'complex' ) then
          do 70 i=1,nnzreq
            read (iunit,fmt=*,end=4000) rpart,ipart
            cval(i) = cmplx(rpart,ipart)
 70      continue
        else
           print *,'''pattern'' data not consistant with type ''array'''
           stop
        endif
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
!
! Various error conditions:
!
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data lines found.'
      stop
 3000 print *,'Size info inconsistant with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistant with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 4000 print *,'Premature end-of-file.'
      print *,'Check that the data file contains ',nnz, &
              ' lines of  i,j,[val] data.'
      print *,'(it appears there are only ',i,' such lines.)'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    end subroutine mmread
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmread
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine mminfo(iunit,rep,field,symm,rows,cols,nnz)
!
      implicit none
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will read header information from a Matrix Market
! formatted file.
!
! The unit iunit must be open, and the file will be rewound on return.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
!             Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   iunit  integer     in   Unit identifier for the open file
!                             containing the data to be read.
!
!   rep     character*10 out  Matrix Market 'representation'
!                             indicator. On return:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                                elemental    (to be added)
!
!   field   character*7  out  Matrix Market 'field'. On return:
!
!                                real
!                                complex
!                                integer
!                                pattern
!
!   symm    character*19 out  Matrix Market 'field'. On return:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     out  Number of rows in matrix.
!
!   cols     integer     out  Number of columns in matrix.
!
!   nnz      integer     out  Number of nonzero entries required to store
!                             the matrix.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer rows, cols, nnz, iunit
      integer count,next
      character mmhead*14
      character mmtype*6
      character rep*10
      character field*7
      character symm*19
      character tmp1*1024
      character tmp2*2
!
! Read header line and check validity:
!
      read (iunit,end=1000,fmt=5) tmp1
 5    format(1024A)
!
! Parse words from header line:
!
      call getwd(mmhead,tmp1,1024,1,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(mmtype,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(rep,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(field,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      call getwd(symm,tmp1,1024,next,next,count)
      if ( count .eq. 0 ) go to 5000
      if ( mmhead .ne. '%%MatrixMarket' ) go to 5000
!
! Convert type code to upper case for easier comparisons:
!
      call lowerc(mmtype,1,6)
      if ( mmtype .ne. 'matrix' ) then
         print *,'Invalid matrix type: ',mmtype
         print *,'This reader only understands type ''matrix''.'
        stop
      else
         call lowerc(rep,1,10)
         call lowerc(field,1,7)
         call lowerc(symm,1,19)
      endif
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' ) &
         go to 6000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' .and. &
          field .ne. 'pattern') go to 7000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' ) go to 8000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and. &
          symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric') &
         go to 9000
!
! Read through comment lines, ignoring content:
!
      read (iunit,end=2000,fmt=200) tmp2
 200  format(1a)
 10   continue
        if ( tmp2(1:1) .ne. '%' ) then
           go to 20
        endif
        read (iunit,end=2000,fmt=200) tmp2
        go to 10
 20   continue
!
! Just read a non-comment.
!   Now, back up a line, and read for first int, and back up
!   again. This will set pointer to just before apparent size
!   info line.
!   Before continuing with free form input, count the number of
!   words on the size info line to ensure there is the right amount
!   of info (2 words for array matrices, 3 for coordinate matrices).
!
      backspace (iunit)
      read (iunit,end=1000,fmt=5) tmp1
      call countwd(tmp1,1024,1,count)
      if ( rep .eq. 'array' .and. count .ne. 2 ) go to 3000
      if ( rep .eq. 'coordinate' .and. count .ne. 3 ) go to 3500
!
!   Correct number of words are present, now back up and read them.
!
      backspace (iunit)
!
      if ( rep .eq. 'coordinate' ) then
!
! Read matrix in sparse coordinate format
!
        read (iunit,fmt=*) rows,cols,nnz
!
! Rewind before returning
!
        rewind(iunit)
        return
!
      elseif ( rep .eq. 'array' ) then
!
! Read matrix in dense column-oriented array format
!
        read (iunit,fmt=*) rows,cols
        if ( symm .eq. 'symmetric' .or. symm .eq. 'hermitian' ) then
          nnz = (rows*cols - rows)/2 + rows
        elseif ( symm .eq. 'skew-symmetric' ) then
          nnz = (rows*cols - rows)/2
        else
          nnz = rows*cols
        endif
!
! Rewind before returning
!
        rewind(iunit)
        return
      else
        print *,'''',rep,''' representation not recognized.'
        print *, 'Recognized representations:'
        print *, '   array'
        print *, '   coordinate'
        stop
      endif
!
! Various error conditions:
!
 1000 print *,'Premature end-of-file.'
      print *,'No lines found.'
      stop
 2000 print *,'Premature end-of-file.'
      print *,'No data found.'
      stop
 3000 print *,'Size info inconsistant with representation.'
      print *,'Array matrices need exactly 2 size descriptors.'
      print *, count,' were found.'
      stop
 3500 print *,'Size info inconsistant with representation.'
      print *,'Coordinate matrices need exactly 3 size descriptors.'
      print *, count,' were found.'
      stop
 5000 print *,'Invalid matrix header: ',tmp1
      print *,'Correct header format:'
      print *,'%%MatrixMarket type representation field symmetry'
      print *
      print *,'Check specification and try again.'
      stop
 6000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 7000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 8000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 9000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine mminfo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mminfo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mmwrite(ounit,rep,field,symm,rows,cols,nnz, &
           indx,jndx,ival,rval,cval)
!
        implicit none
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This routine will write data to a matrix market formatted file.
! The data may be either sparse coordinate format, or dense array format.
!
! The unit ounit must be open.
!
! 20-Sept-96  Karin A. Remington, NIST ACMD (karin@cam.nist.gov)
! 18-Oct-96   Change in routine name to match C and Matlab routines.
! 30-Oct-96   Bug fixes in mmio.f:
!                  -looping for comment lines
!                  -fixed non-ansi zero stringlength
!                  -incorrect size calculation for skew-symmetric arrays
!             Other changes in mmio.f:
!                  -added integer value parameter to calling sequences
!                  -enforced proper count in size info line
!                  -added routine to count words in string (countwd)
!            (Thanks to G.P.Leendetse and H.Oudshoom for their review
!             of the initial version and suggested fixes.)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Arguments:
!
!   name     type      in/out description
!   ---------------------------------------------------------------
!
!   ounit  integer     in   Unit identifier for the file
!                             to which the data will be written.
!                             Must be open prior to call.
!
!   rep     character*   in   Matrix Market 'representation'
!                             indicator. Valid inputs:
!
!                                coordinate   (for sparse data)
!                                array        (for dense data)
!                               *elemental*    (to be added)
!
!   field   character*   in   Matrix Market 'field'. Valid inputs:
!
!                                real
!                                complex
!                                integer
!                                pattern (not valid for dense arrays)
!
!   symm    character*   in   Matrix Market 'field'. Valid inputs:
!
!                                symmetric
!                                hermitian
!                                skew-symmetric
!                                general
!
!   rows     integer     in   Number of rows in matrix.
!
!   cols     integer     in   Number of columns in matrix.
!
!   nnz      integer     in   Number of nonzero entries in matrix.
!                             (rows*cols for array matrices)
!
!   indx     integer(nnz)in   Row indices for coordinate format.
!                             Undefined for array format.
!
!   jndx     integer(nnz)in   Column indices for coordinate format.
!                             Undefined for array format.
!
!   ival     integer(nnz) in  Integer data (if applicable, see 'field')
!
!   rval     double(nnz) in   Real data (if applicable, see 'field')
!
!   cval     complex(nnz)in   Complex data (if applicable, see 'field')
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Declarations:
!
      integer ival(*)
      double precision rval(*)
      complex cval(*)
      integer indx(*)
      integer jndx(*)
      integer i, rows, cols, nnz, nnzreq, ounit
      character*(*)rep,field,symm
!
! Test input qualifiers:
!
      if (rep .ne. 'coordinate' .and. rep .ne. 'array' ) &
         go to 1000
      if (rep .eq. 'coordinate' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' .and. &
          field .ne. 'pattern') go to 2000
      if (rep .eq. 'array' .and. field .ne. 'integer' .and. &
          field .ne. 'real' .and. field .ne. 'complex' ) go to 3000
      if (symm .ne. 'general' .and. symm .ne. 'symmetric' .and. &
          symm .ne. 'hermitian' .and. symm .ne. 'skew-symmetric') &
         go to 4000
!
! Write header line:
!
      write(unit=ounit,fmt=5)rep,' ',field,' ',symm
 5    format('%%MatrixMarket matrix ',11A,1A,8A,1A,20A)
!
! Write size information:
!
      if ( rep .eq. 'coordinate' ) then
         nnzreq=nnz
         write(unit=ounit,fmt=*) rows,cols,nnz
         if ( field .eq. 'integer' ) then
            do 10 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),ival(i)
 10         continue
         elseif ( field .eq. 'real' ) then
            do 20 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i),rval(i)
 20         continue
         elseif ( field .eq. 'complex' ) then
            do 30 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i), &
                                      real(cval(i)),aimag(cval(i))
 30         continue
         else
!        field .eq. 'pattern'
            do 40 i=1,nnzreq
               write(unit=ounit,fmt=*)indx(i),jndx(i)
 40         continue
         endif
      else
!        rep .eq. 'array'
         if ( symm .eq. 'general' ) then
           nnzreq = rows*cols
         elseif ( symm .eq. 'symmetric' .or. &
                  symm .eq. 'hermitian' ) then
           nnzreq = (rows*cols - rows)/2 + rows
         else
!        symm .eq. 'skew-symmetric'
           nnzreq = (rows*cols - rows)/2
         endif
         write(unit=ounit,fmt=*)rows,cols
         if ( field .eq. 'integer' ) then
            do 50 i=1,nnzreq
               write(unit=ounit,fmt=*)ival(i)
 50         continue
         elseif ( field .eq. 'real' ) then
            do 60 i=1,nnzreq
               write(unit=ounit,fmt=*)rval(i)
 60         continue
         else
!        field .eq. 'complex'
            do 70 i=1,nnzreq
               write(unit=ounit,fmt=*)real(cval(i)),aimag(cval(i))
 70         continue
         endif
      endif
      return
!
! Various errors
!
 1000 print *,'''',rep,''' representation not recognized.'
      print *, 'Recognized representations:'
      print *, '   array'
      print *, '   coordinate'
      stop
 2000 print *,'''',field,''' field is not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      print *, '   pattern'
      stop
 3000 print *,'''',field,''' arrays are not recognized.'
      print *, 'Recognized fields:'
      print *, '   real'
      print *, '   complex'
      print *, '   integer'
      stop
 4000 print *,'''',symm,''' symmetry is not recognized.'
      print *, 'Recognized symmetries:'
      print *, '   general'
      print *, '   symmetric'
      print *, '   hermitian'
      print *, '   skew-symmetric'
      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine mmwrite
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! End of subroutine mmwrite
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lowerc(string,pos,len)
!
        implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Convert uppercase letters to lowercase letters in string with
! starting postion pos and length len.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer pos, len, i, k
      character*(*) string

      character*26 lcase, ucase
      save lcase,ucase
      data lcase/'abcdefghijklmnopqrstuvwxyz'/
      data ucase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do 10 i=pos,len
        k = index(ucase,string(i:i))
        if (k.ne.0) string(i:i) = lcase(k:k)
 10   continue
      return
      end subroutine lowerc

      subroutine getwd(word,string,slen,start,next,wlen)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Getwd extracts the first  word from string starting
!     at position start.  On return, next is the position
!     of the blank which terminates the word in string.
!     If the found word is longer than the allocated space
!     for the word in the calling program, the word will be
!     truncated to fit.
!     Count is set to the length of the word found.
!
! 30-Oct-96   Bug fix: fixed non-ansi zero stringlength
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer i, slen, start, next, begin, space, wlen
      character*(*) word
      character*(*) string

      begin = start
      do 5 i=start,slen
         space = index(string(i:slen),' ')
         if ( space .gt. 1) then
            next = i+space-1
            go to 100
         endif
         begin=begin+1
 5    continue
 100  continue
      wlen=next-begin
      if ( wlen .le. 0 ) then
        wlen = 0
        word = ' '
        return
      endif
      word=string(begin:begin+wlen)
      return
    end subroutine getwd

    subroutine countwd(string,slen,start,count)
!
      implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Countwd counts the number of words in string starting
!     at position start.  On return, count is the number of words.
! 30-Oct-96   Routine added
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*(*) string
      integer slen, start, next, wordlength, count
      character tmp2*2

      count = 0
      next = 1
 10   call getwd(tmp2,string,1024,next,next,wordlength)
      if ( wordlength .gt. 0 ) then
         count = count + 1
         go to 10
      endif
      return
    end subroutine countwd
!!
!! The six routines above are for reading MM format sparse matrices. 

end module
