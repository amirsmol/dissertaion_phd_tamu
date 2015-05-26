module linearsolvers
!use mkl95_lapack

   implicit none
      logical,parameter::updated_lagrangian=.false.

      integer:: i_calibration
      integer::gauss_point_number
      integer::ngauss
      integer::curved_node

      integer,parameter :: ndf=4 !number of degrees of freedom
      integer,parameter :: dimen=3 !dimention of problem
      integer,parameter::iwp=selected_real_kind(15)

      integer,parameter::out=2 !output file
      integer,parameter::csv=11 !comma seperated file
      integer,parameter::msh=3 ! mesh file

      integer,parameter::vtu=13 !output file
      integer,parameter::OUTPUT_UNIT=7

   ! the default value for the smallest pivot that will be accepted
   ! using the linearsolvers subroutines.  pivots smaller than this
   ! threshold will cause premature termination of the linear equation
   ! solver and return false as the return value of the function.

   real(iwp), parameter :: default_smallest_pivot = 1.0e-6
   integer::iter


contains

function absh(x)
real(iwp)::x,absh
absh=x*tanh(10*x)
end function absh


function steph(x)
real(iwp)::x,steph
steph=(tanh(10*x)+1)/2
end function steph

function macaulayh(x)
real(iwp)::x,macaulayh
macaulayh=(absh(x)+x)/2
end function macaulayh

function sech(x)
real(iwp)::x,sech
sech=2.0d0*exp(x)/(1.0d0+exp(2.0d0*x))
end function sech


function macaulay_brackets(x)
real(iwp)::macaulay_brackets
real(iwp)::x
macaulay_brackets=x
if (x.lt.0) macaulay_brackets=0.0d0
endfunction macaulay_brackets

function signum(x)
real(iwp)::signum
real(iwp)::x
signum=1.0
if (x.lt.0) signum=-1.0d0
if (x.eq.0) signum= 0.0d0
end function signum


function signumh(x)
real(iwp)::signumh
real(iwp)::x
signumh=tanh(100*x)
end function signumh

function sign_function(x)
real(iwp)::sign_function
real(iwp)::x
sign_function=1.0d0
if (x.lt.0) sign_function=0.0d0
endfunction sign_function

function sawtooth(time)
implicit none
real(iwp)::time,sawtooth
sawtooth=4*(abs(time+0.25-floor(time+0.75))-0.25)
end function sawtooth

!     ================================================================
!                 norm of a vector
!     ================================================================
function norm_vect(a)
real(iwp)::a(:),norm_vect
norm_vect=dsqrt(dot_product(a,a))
end function norm_vect

subroutine lapack_gesv(a,b)
implicit none
! ============================================================================
! name        : lapack solver compact
! author      : amir
! version     :
! copyright   : your copyright notice
! description : this solves a * x=b and put the solution
! into the b vector
! ============================================================================


real(iwp)::a(:,:),b(:)
integer::n, nrhs, lda, ldb, info
integer,allocatable::ipiv(:)

n=size(b)
nrhs=1
lda=n
ldb=n;
allocate(ipiv(n))


call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

!call gesv(a,b)


end subroutine lapack_gesv


function lapack_matmul(a,b)
      implicit none
      integer, parameter :: nmax=1000

! ============================================================================
! name        : lapack solver compact
! author      : amir
! version     :
! copyright   : your copyright notice
! description : this solves a * x=b and put the solution
! into the b vector
! ============================================================================
!  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

! NAME
!      DGEMM - perform one of the matrix-matrix operations   C :=
!      alpha*op( A )*op( B ) + beta*C,

real(iwp)::a(:,:),b(:,:)
real(iwp),allocatable::lapack_matmul(:,:),c(:,:)
real(iwp),parameter::alpha=1.0d0,beta=1.0d0
integer::m, n, k
integer::lda, ldb, ldc

!! lapack_matmul_{m*n} = a_{m*k} * b_{k*n}


m=size(a,dim=1)
k=size(a,dim=2)
n=size(b,dim=2)




allocate(lapack_matmul(m,n),c(m,n))

lda=size(a,dim=1)
ldb=size(b,dim=1)
ldc=size(c,dim=1)

lapack_matmul=0.0d0
c=0.0d0
! call dgemm('n'   ,'n'   ,m,n,k,myone  ,Z,n  ,Z,n  ,myone ,X,n)

call dgemm ( 'n', 'n', m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

lapack_matmul=c

end function lapack_matmul


subroutine parafem_cg_solver(glk,glr)
implicit none
integer, parameter :: nmax=1000

! ============================================================================
! name        : conjugate gradient solver taken from parafem
! author      : amir
! version     :
! copyright   :
! description : this solves a * x=b and put the solution
! into the b vector
! ============================================================================
real(iwp)::glk(:,:),glr(:)
!!    ============================pcg solver variables
integer::cg_iters,cg_limit,neq
real(iwp)::alpha,beta,cg_tol  , zero=0.0_iwp ! ,one=1.0_iwp
real(iwp)::up
logical::cg_converged ! ,solid=.true.
integer::i
real(iwp),allocatable::diag_precon(:),x(:),xnew(:)
real(iwp),allocatable::p(:),d(:),u(:),loads(:)

neq=size(glr)

 allocate(p(neq),loads(neq),x(neq),xnew(neq),u(neq),            &
          diag_precon(neq),d(neq))


diag_precon=0.0D0

do i=1,neq; diag_precon(i)=1/glk(i,i);enddo
 loads=glr
 d=diag_precon*loads;
 p=d;
 x=zero;
 cg_iters=0

cg_limit=200
pcg: DO; cg_iters=cg_iters+1;
u=zero
u=u+MATMUL(glk,p)
up=DOT_PRODUCT(loads,d)
alpha=up/DOT_PRODUCT(p,u);
xnew=x+p*alpha;
loads=loads-u*alpha
d=diag_precon*loads;
beta=DOT_PRODUCT(loads,d)/up;
p=d+p*beta
if(norm_vect(xnew-x).le.norm_vect(x)*cg_tol) cg_converged=.true.
IF(cg_converged.OR.cg_iters==cg_limit)EXIT
END DO pcg
glr=xnew


end subroutine parafem_cg_solver



    !     ================================================================
    !                          show matrix subroutine
    !     ================================================================
    subroutine show_matrix(a, namea)
        implicit none
        !     ================================================================
        !                           the variables
        !     ================================================================
        real(iwp), intent(in)         :: a(:, :) ! this is the matrix
!        integer         :: i1, i2 !these are dimentions
        character (len=*), intent(in) :: namea ! the name of the matrix
        integer            :: i ,idim, jdim,output_unit
        integer,save::counter
!        character (len=2)  :: jdim_string
!        character (len=16) :: fmtstring
        !     ================================================================
        !                           the printers
        !     ================================================================
        counter=counter+1
        output_unit=10
        open (output_unit,file='fem_matrices_output_file.txt')
        !      write(output_unit,*)'counter=',counter
        write(output_unit,*)namea
        idim = size(a, 1)
        jdim = size(a, 2)
        !            "(f10.6,a2)"

!        write(fmtstring,fmt="(a1,I3,a8)")"(",jdim,"(f10.6))"
        do i=1,idim
            !      write(output_unit,fmt=(fmtstring))a(i,:)
        write(output_unit,*)a(i,:)
        enddo
        write(output_unit,*)
    end subroutine show_matrix


    !     ================================================================
    !                          show matrix subroutine
    !     ================================================================
    subroutine show_matrix_int(a, namea)
        implicit none
        !     ================================================================
        !                           the variables
        !     ================================================================
        integer, intent(in)         :: a(:, :) ! this is the matrix
        character (len=*), intent(in) :: namea ! the name of the matrix
        integer            :: i, idim,jdim,output_unit
        integer,save::counter


        !     ================================================================
        !                           the printers
        !     ================================================================
        counter=counter+1
        output_unit=10
        open (output_unit,file='fem_matrices_output_file.txt')
        write(output_unit,*)namea
        idim = size(a, 1)
        jdim = size(a, 2)
        do i=1,idim
            write(output_unit,*)a(i,:)
        enddo
        write(output_unit,*)
    end subroutine show_matrix_int
    !     ================================================================
    !                          show vector subroutine
    !     ================================================================
    subroutine show_vector(a, namea)
        implicit none
        !     ================================================================
        !                           the variables
        !     ================================================================
        real(iwp), intent(in)         :: a(:) ! this is thevector
        character (len=*), intent(in) :: namea ! the name of the matrix
        integer            :: output_unit
        integer::idim
        integer,save::counter

        !     ================================================================
        !                           the printers
        !     ================================================================
        counter=counter+1
        output_unit=10
        idim = size(a, 1)
        open (output_unit,file='fem_matrices_output_file.txt')
        !      write(output_unit,*)'counter=',counter
        write(output_unit,*)namea

        do idim=1,size(a, 1)
            write(output_unit,*)a(idim);enddo
            write(output_unit,*)
        end subroutine show_vector



subroutine Gaussian_Elimination_Solver(a,y)
    integer :: k,i
    real(iwp)::a(:,:),y(:)
    integer::n
    n=size(a,dim=1)

    do k = 1, n-1
        a(k+1: n, k) = a(k+1: n, k) / a(k, k)

        a(k+1: n, k+1: n) = a(k+1: n, k+1: n) - &
                matmul(a(k+1: n, k: k), a(k: k, k+1: n))
    end do


    ! L x = f  =>  x = L \ f
    do i = 1, n
        y(i) = y(i) - dot_product(a(i, 1: i-1), y(1: i-1))
    end do

    ! U y = x  =>  y = U \ x
    do i = n, 1, -1
        y(i) = y(i) - dot_product(a(i, i+1: n), y(i+1: n))
        y(i) = y(i) / a(i, i)
    end do

end subroutine Gaussian_Elimination_Solver


!     ===============================================================
!                       TRACE OF A 3X3 TENSOR
!     ===============================================================

      FUNCTION TRACE(A)
      real(iwp)::A(:,:),TRACE
      INTEGER::I
      TRACE=0.0D0;
      DO I=1,SIZE(A,DIM=1);TRACE=TRACE+A(I,I);ENDDO
      END FUNCTION TRACE

      FUNCTION MATMUL_EPZ(EPZ,ELECT)
      real(iwp)::ELECT(3),EPZ(3,3,3),MATMUL_EPZ(3,3)
      INTEGER::I,J,K
      MATMUL_EPZ=0.0D0;
      DO I=1,3;
      DO J=1,3;
      DO K=1,3;
      MATMUL_EPZ(I,J)=MATMUL_EPZ(I,J)+EPZ(K,I,J)*ELECT(K)
      ENDDO;ENDDO;ENDDO
      END FUNCTION MATMUL_EPZ

      FUNCTION MATMUL_D_EPZ(EPZ,STRAN)
      real(iwp)::MATMUL_D_EPZ(3),EPZ(3,3,3),STRAN(3,3)
      INTEGER::I,J,K
      MATMUL_D_EPZ=0.0D0;
      DO I=1,3;
      DO J=1,3;
      DO K=1,3;
      MATMUL_D_EPZ(K)=MATMUL_D_EPZ(K)+EPZ(K,I,J)*STRAN(I,J)
      ENDDO;ENDDO;ENDDO
      END FUNCTION MATMUL_D_EPZ

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV,DET, OK_FLAG)

      IMPLICIT NONE

      real(iwp), DIMENSION(3,3), INTENT(IN)  :: A
      real(iwp), DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      real(iwp), PARAMETER :: EPS = 1.0D-20
      real(iwp) :: DET
      real(iwp), DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)


      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET
      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV

!     ================================================================
!                 NORM OF A VECTOR
!     ================================================================
      FUNCTION NORM_MATRIX(A)
      INTEGER::I,J
      real(iwp)::A(:,:),NORM_MATRIX


      NORM_MATRIX=0.0D0
      DO I=1,SIZE(A,DIM=1);DO J=1,SIZE(A,DIM=2);
      NORM_MATRIX=NORM_MATRIX+A(I,J)*A(I,J)
      ENDDO;ENDDO


      END FUNCTION NORM_MATRIX

!     ================================================================
!                 MACAULAY BRACKETS
!     ================================================================
      FUNCTION MACAULAY(A)
      real(iwp)::MACAULAY
      real(iwp)::A
      MACAULAY=0
      IF (A>=0)THEN
       MACAULAY=1
      ENDIF
      END FUNCTION MACAULAY

!     ================================================================
!                 MACAULAY BRACKETS
!     ================================================================
      FUNCTION MACAULAY_s(A)
      real(iwp)::MACAULAY_s
      real(iwp)::A
      MACAULAY_s=0.5*(abs(a)+a)

      END FUNCTION MACAULAY_s
!     ================================================================
!                 MACAULAY BRACKETS
!     ================================================================
      FUNCTION MACAULAY_H(A)
      real(iwp)::MACAULAY_H
      real(iwp)::A
      MACAULAY_h=0.5*(tanh(a)+1)

      END FUNCTION MACAULAY_H

!     ================================================================
!                 UNIT VECTOR IN THE DIRECTION OF A VECTOR
!     ================================================================
      FUNCTION UNIT_VECT(A)
      real(iwp)::A(:)
      real(iwp),ALLOCATABLE :: UNIT_VECT(:)

      ALLOCATE(UNIT_VECT(SIZE(A)))
      UNIT_VECT=0.0D0

      UNIT_VECT=A/NORM_VECT(A)

      IF(NORM_VECT(A).EQ.0)UNIT_VECT=0.0D0


      END FUNCTION UNIT_VECT


!
!      subroutine elapsedtime (END_STAMP_V,START_STAMP_V)
!      integer ( kind = 4 )START_STAMP_V(8),END_STAMP_V(8),values(8),SECONDS,MINUTES,HOURS,MILISECONDS,START_STAMP,END_STAMP,ELAPSED_TIME
!
!      MILISECONDS=mod(END_STAMP_V(8)-START_STAMP_V(8),100)
!      SECONDS=mod(END_STAMP_V(7)-START_STAMP_V(7),60)
!      MINUTES=mod(END_STAMP_V(6)-START_STAMP_V(6),60)
!      HOURS=END_STAMP_V(5)-START_STAMP_V(5)
!      write ( *,'(1X,a,2x,i2,a1,i4.4,a1,i4.4,a1,i4.4)' )'Elapsed Time=', HOURS,':',MINUTES,':',SECONDS,'.',ABS(MILISECONDS)
!      end subroutine elapsedtime


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end SUBROUTINE

      SUBROUTINE M22INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR


      DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +A(2,2)
      COFACTOR(1,2) = -A(2,1)
      COFACTOR(2,1) = -A(1,2)
      COFACTOR(2,2) = +A(1,1)

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M22INV

function rotate_deg_rz(thetadeg)
    real(iwp)::rotate_deg_rz(3,3),theta,thetadeg,pi
    pi=datan(1.0d0)*4

    theta=(thetadeg/180.0)*pi


    rotate_deg_rz=0.0d0;

    rotate_deg_rz(1,1)= cos(theta)
    rotate_deg_rz(2,2)= cos(theta)
    rotate_deg_rz(1,2)=-sin(theta)
    rotate_deg_rz(2,1)= sin(theta)
    rotate_deg_rz(3,3)= 1.0d0



end function rotate_deg_rz

!     ===============================================================
!   function for cross product normal
!     ===============================================================

function cross3n(a,b) result(axbn)

    implicit none
    integer,parameter :: wp=selected_real_kind(15,307) !double precision
    real(wp)::euc_norm
    real(wp),dimension(3) :: axb,axbn
    real(wp),dimension(3),intent(in) :: a
    real(wp),dimension(3),intent(in) :: b

    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

    euc_norm=sqrt(dot_product(axb,axb))

    if(euc_norm.eq.0)then
        axbn=0;
    else
        axbn=axb/euc_norm;
    endif

end function cross3n

!====================================================================
subroutine generate_trans_tensor(two_point,qtran)
!====================================================================
!      code related to the phd thesis of amir sohrabi mollayousef
!      programmer(s):amir sohrabi
!      objective: to find the rotation tensor for a given point with respect
!      to the origin of the problem. this transformation matrix is to be used
!      in non linear electromechancial analysis of telescopic actoator.
!      qtran transforms a vector e1 to the element local frame
!=====================================================================
      implicit none
!=====================================================================
      real(iwp),intent(in) :: two_point(:,:) ! full local material matrix local
      real(iwp),intent(out) :: qtran(:,:) ! full global material matrix local
!     =================================================================
      real(iwp):: lent_prog(dimen)!transformation_matrix
      real(iwp):: lent_ele !lenght of element and its projection on x1, x2
      real(iwp) :: aiden(ndf,ndf) ! the identity matrix
      real(iwp) :: truss_local_system(dimen,dimen) ! the identity matrix
      real(iwp) :: refrence_system(dimen,dimen) ! the identity matrix
      real(iwp) :: element_direction_normal(dimen) ! the identity matrix

      integer::i,j
!     ================================================================
!     the transformation matrix
!     ================================================================
      aiden = 0.0; do i = 1,ndf; aiden(i,i) = 1.0d0;enddo
      refrence_system=aiden(1:dimen,1:dimen);

      lent_prog(:)= two_point(2,:)-two_point(1,:)
      lent_ele=dsqrt(dot_product(lent_prog,lent_prog))

      lent_ele=norm_vect(lent_prog)

      element_direction_normal=lent_prog/lent_ele


      truss_local_system=0.0d0
      truss_local_system(1,:)=lent_prog/lent_ele
      truss_local_system(3,:)= cross3n(truss_local_system(1,:),refrence_system(2,:))

            if(norm_vect(cross3n(truss_local_system(1,:),refrence_system(2,:))).lt.1e-3)then
            truss_local_system(3,:)=refrence_system(3,:)
            endif

      truss_local_system(2,:)= cross3n(truss_local_system(3,:),truss_local_system(1,:))
      do 10 i=1,dimen
      do 10 j=1,dimen
10    qtran(i,j)=dot_product(refrence_system(i,:),truss_local_system(j,:))

      end subroutine generate_trans_tensor

    end module linearsolvers
