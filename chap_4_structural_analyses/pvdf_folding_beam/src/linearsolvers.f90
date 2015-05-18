module linearsolvers
!use mkl95_lapack

   implicit none
      logical,parameter::updated_lagrangian=.false.


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
! copyright   : your copyright notice
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



    SUBROUTINE CONJGRAD_SOLVER(A,B,X)
    IMPLICIT NONE
    real(iwp),DIMENSION(:),INTENT(IN)::B
    real(iwp),DIMENSION(:)::X
    real(iwp),DIMENSION(:,:),INTENT(IN)::A
    real(iwp),DIMENSION(:),ALLOCATABLE :: P,R,AP
    real(iwp)::RSOLD,RSNEW,ALPHA
    INTEGER::DIMEN
! ========================CONVERGENCE VARIABLES
    INTEGER::CG_ITERS,CG_LIMIT
    LOGICAL::CG_CONVERGED
    real(iwp)::CG_TOL
! ================================================================

    DIMEN=SIZE(B)
    CG_TOL=1E-5
    CG_LIMIT=500
    CG_ITERS=0

    ALLOCATE(R(DIMEN),P(DIMEN),AP(DIMEN))
    R=B-MATMUL(A,X)
    P=R
    RSOLD=DOT_PRODUCT(R,R)

 PCG: DO; CG_ITERS=CG_ITERS+1;
    AP=MATMUL(A,P)
    ALPHA=RSOLD/DOT_PRODUCT(P,AP)
    X=X+ALPHA*P
    R=R-ALPHA*AP
    RSNEW=DOT_PRODUCT(R,R)
    CG_CONVERGED=(SQRT(RSNEW).LT.CG_TOL)
    IF(CG_CONVERGED.OR.CG_ITERS==CG_LIMIT)EXIT
       P=R+(RSNEW/RSOLD)*P;
       RSOLD=RSNEW;
     END DO PCG
     WRITE(1,*) CG_ITERS
    end subroutine CONJGRAD_SOLVER



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

   ! Output A in Matlab format, using name in the Matlab assignment statement.
   subroutine printMatrix( A, name )
      implicit none
      real(iwp), dimension( :, : ) :: A   ! Assume the shape of A.
      character name  ! Name for use in assignment, ie, name = ......

      integer n, m, i, j


      OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      n = size( A, 1 )
      m = size( A, 2 )

      write(OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no `;'.
      do i = 1, n-1

         ! Output current row.
         do j = 1, m-1
            write( OUTPUT_UNIT, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
         end do

         ! Output last element in row and end current row.
         write( OUTPUT_UNIT, fmt="(f10.6,a1)" ) A( i, m ), ';'

      end do

      ! Output the last row.
      do j = 1, m-1
         write( OUTPUT_UNIT, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
      end do

      ! Output last element in row and end.
      write( OUTPUT_UNIT, fmt="(f10.6,a1)" ) A( i, m ), ']'

   end subroutine printMatrix


   subroutine printMatrix_int( A, name )
      implicit none
      integer, dimension( :, : ) :: A   ! Assume the shape of A.
      character name  ! Name for use in assignment, ie, name = ......

      integer n, m, i, j

      n = size( A, 1 )
      m = size( A, 2 )
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no `;'.
      do i = 1, n-1

         ! Output current row.
         do j = 1, m-1
            write( OUTPUT_UNIT, fmt="(I3,a2)", advance = "no" ) A( i, j ), ', '
         end do

         ! Output last element in row and end current row.
         write( OUTPUT_UNIT, fmt="(I3,a2)" ) A( i, m ), ';'

      end do

      ! Output the last row.
      do j = 1, m-1
         write( OUTPUT_UNIT, fmt="(I3,a2)", advance = "no" ) A( i, j ), ', '
      end do

      ! Output last element in row and end.
      write( 7, fmt="(I3,a2)" ) A( i, m ), ']'

   end subroutine printMatrix_int


   ! Output b in Matlab format, using name in the Matlab assignment statement.
   subroutine printVector_int( b, name )
      implicit none
      integer, dimension( : ) :: b   ! Assume the shape of b.
      character name   ! Name for use in assignment, ie, name = ......

      integer n, i
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      n = size( b )

      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      do i = 1, n-1
         write( OUTPUT_UNIT, fmt = "(i3,a2)", advance = "no" ) b( i ), ', '
      end do

      write( OUTPUT_UNIT, fmt = "(i3,a2)" ) b( n ), ']'''

   end subroutine printVector_int

     subroutine printVector( b, name )
      implicit none
      real(iwp), dimension( : ) :: b   ! Assume the shape of b.
      character name   ! Name for use in assignment, ie, name = ......

      integer n, i

      n = size( b )
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')

      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      do i = 1, n-1
         write( OUTPUT_UNIT, fmt = "(E14.5,a2)", advance = "no" ) b( i ), ', '
      end do

      write( OUTPUT_UNIT, fmt = "(E14.5,a2)" ) b( n ), ']'''

   end subroutine printVector



      SUBROUTINE SLVUNSYM(A,NRMAX,NCMAX,N,ITERM)
!C     _________________________________________________________________
!C
!C       Solver for BANDED UNSYMMETRIC system of algebraic equations
!C     ______________________________________________________________
!C
      IMPLICIT real(iwp) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER::NRMAX,NCMAX,N,ITERM
      real(iwp)::A(NRMAX,NCMAX)
      CERO=1.0D-15
      PARE=CERO**2
      NBND=2*ITERM
      NBM=NBND-1
!C
!C     Begin elimination of the lower left
!C
      DO 80 I=1,N
      IF (DABS(A(I,ITERM)).LT.CERO) GO TO 10
      GO TO 20
   10 IF (DABS(A(I,ITERM)).LT.PARE) GO TO 110
   20 JLAST=MIN0(I+ITERM-1,N)
      L=ITERM+1
      DO 40 J=I,JLAST
      L=L-1
      IF (DABS(A(J,L)).LT.PARE) GO TO 40
      B=A(J,L)
      DO 30 K=L,NBND
   30 A(J,K)=A(J,K)/B
      IF (I.EQ.N) GO TO 90
   40 CONTINUE
      L=0
      JFIRST=I+1
      IF (JLAST.LE.I) GO TO 80
      DO 70 J=JFIRST,JLAST
      L=L+1
      IF (DABS(A(J,ITERM-L)).LT.PARE) GO TO 70
      DO 50 K=ITERM,NBM
   50 A(J,K-L)=A(J-L,K)-A(J,K-L)
      A(J,NBND)=A(J-L,NBND)-A(J,NBND)
      IF (I.GE.N-ITERM+1) GO TO 70
      DO 60 K=1,L
   60 A(J,NBND-K)=-A(J,NBND-K)
   70 CONTINUE
   80 CONTINUE
   90 L=ITERM-1
      DO 100 I=2,N
      DO 100 J=1,L
      IF (N+1-I+J.GT.N) GO TO 100
      A(N+1-I,NBND)=A(N+1-I,NBND)-A(N+1-I+J,NBND)*A(N+1-I,ITERM+J)
  100 CONTINUE
      RETURN
  110 WRITE (*,140) I,A(I,ITERM)
  140 FORMAT (/,2X,'Computation stopped in SLVUNSYM because zero appears on the main diagonal *** Eqn no. and value:',I5,E12.4)
      end subroutine SLVUNSYM

! -------------------------- MODULE cg.f90 ----------------------------

!************************************************************************
!*                                                                      *
!* Conjugate Gradient Method (CG Method)                                *
!* -------------------------------------                                *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Compiler:             Turbo C 2.0                                    *
!* Computer:             IBM PS/2 70 with 80387                         *
!* Sources:              [BUNS85], [SCHW], [MAES84]                     *
!* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
!* Date:                 7.31.1992                                      *
!*                                                                      *
!*             F90 version by J-P Moreau (without dynamic allocations). *
!*                                (www.jpmoreau.fr)                     *
!************************************************************************
Subroutine cg_method (    &     ! Conjugate Gradient Method
                       n, &     ! Size of the linear system
                       a, &     ! System matrix
                       y, &     ! right hand side
                       x, &     ! solution vector
                       fehler & ! error code
                     )
! original name: cg_verfahren()
integer::n
real(iwp), parameter :: ZERO=0.d0, MACH_EPS=2.d-16
real(iwp)  a(0:n-1,0:n-1),x(0:n-1),y(0:n-1)
integer fehler
!************************************************************************
!* CG solves the linear system                                          *
!*                         A * X = Y                                    *
!* for a symmetric, positive definite matrix A via the conjugate        *
!* gradient method.                                                     *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* n  Size of the linear system                                         *
!* a  [0..n-1,0..n-1] system matrix A. Only the upper triangle of A is  *
!*    used.                                                             *
!* y  [0..n-1] vector of the right hand side                            *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* x  [0..n-1] vector giving the solution                               *
!*                                                                      *
!* Return value:                                                        *
!* =============                                                        *
!* = 0: all is ok                                                       *
!* = 1: n < 2 or other disallowed input parameters                      *
!* = 2: memory exceeded                                                 *
!*                                                                      *
!************************************************************************
  real(iwp) d(0:n-1), &   ! (0..n-1) auxiliary vectors d and g
         g(0:n-1), &
         AmalD(0:n-1)  ! (0..n-1) auxiliary vector A * d
  real(iwp) alpha,    &   ! coefficient
         beta,     &   ! coefficient
         dividend, &   ! numerator and denominator of a fraction
         divisor,  &   ! respectively, used to compute alpha, beta
         hilf,     &   ! auxiliary variables
         hilf2,    &
         abstand,  &   ! distance of two successive approximations
                       ! for the solution vector x (taken in the
                       ! euclidean norm)
         xnorm         ! euklidean norm of x
  integer k, i, j      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
    return
  end if

  !------------------------------------------------------------------
  ! start with x at the origin
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    x(i) = ZERO
  end do

  !------------------------------------------------------------------
  ! initialize  d and g :
  ! d = -g = -(a*x - y) = y (since x = 0)
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    hilf = y(i)
    d(i) = hilf
    g(i) = -hilf
  end do


  !------------------------------------------------------------------
  ! perform at most n steps of the CG Method
  !------------------------------------------------------------------
  do k = n, 0, -1

    !----------------------------------------------------------------
    ! compute new alpha:
    ! alpha = -(d(transp) * g) / (d(transp) * (a * d))
    !----------------------------------------------------------------

    dividend = ZERO
    divisor  = ZERO

    do i = n - 1, 0, -1
      dividend = dividend + d(i) * g(i)
      hilf = ZERO
      do j = 0, i-1
        hilf = hilf + a(j,i) * d(j)
      end do
      do j = i, n-1
        hilf = hilf + a(i,j) * d(j)
      end do
      AmalD(i) = hilf
      divisor = divisor + d(i) * hilf
    end do

    if (divisor.eq.ZERO) then
      fehler=0
      return
    end if

    alpha = -dividend / divisor

    !----------------------------------------------------------------
    ! compute the norm of x und  alpha * d  and find a new x:
    ! x  =  x + alpha * d, then check whether x is close enough,
    ! in order to stop the process before n complete steps
    !----------------------------------------------------------------
    xnorm   = ZERO
    abstand = ZERO

    do i = n - 1, 0, -1
      hilf =  x(i)
      xnorm   = xnorm + hilf*hilf
      hilf2   =  alpha * d(i)
      abstand = abstand + hilf2*hilf2
      x(i)    =  hilf + hilf2
    end do

    if (abstand < MACH_EPS * xnorm) then
      fehler=0
      return
    end if


    !----------------------------------------------------------------
    ! compute new g:   g  =  g + alpha * (a * d)
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      g(i) = g(i) + alpha * AmalD(i)
    end do

    !----------------------------------------------------------------
    ! compute new beta :
    ! beta = (g(transp) * (a * d)) / (d(transp) * (a * d))
    !----------------------------------------------------------------
    dividend = ZERO
    do i = n - 1, 0, -1
      dividend = dividend + g(i) * AmalD(i)
    end do

    beta = dividend / divisor

    !----------------------------------------------------------------
    ! compute new d :   d  =  - g + beta * d
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      d(i) = -g(i) + beta * d(i)
    end do

  end do  !k loop

  fehler=0
  return
end subroutine cg_method



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



      SUBROUTINE TENSOR_TRANSFORM_RANK_1(TRANSFORMATION,TENSOR_LOCAL_COORD,TENSOR_GLOBAL_COORD)
!-------------------------------------------------------------------
!This will transform the tensor from local into the global coordinate
!it used the transformation tensor to present the tensor that is
!represented in local coordinate in the local coordinate
!-------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------
      real(iwp),INTENT(IN) :: TRANSFORMATION(3,3) ! TRANSFORMATION TENSOR
      real(iwp),INTENT(IN) :: TENSOR_LOCAL_COORD(3) ! LOCAL COORDINATE TENSOR
      real(iwp),INTENT(OUT) :: TENSOR_GLOBAL_COORD(3) ! GLOBAL COORDINATE TENSOR
      INTEGER :: I,J ! The counter integers
! ----------------------------------------------------------
      TENSOR_GLOBAL_COORD=0.0D0 !INITIALIZING TENSOR
! ----------------------------------------------------------
!   IMPLEMENTING THE SYMMETRY
    DO  I=1,3
    DO  J=1,3
    TENSOR_GLOBAL_COORD(I)=TENSOR_GLOBAL_COORD(I)+TRANSFORMATION(I,J)*TENSOR_LOCAL_COORD(J)
    ENDDO;ENDDO
      end subroutine TENSOR_TRANSFORM_RANK_1
!===================================================================
!
      SUBROUTINE TENSOR_TRANSFORM_RANK_2(TRANSFORMATION,TENSOR_LOCAL_COORD,TENSOR_GLOBAL_COORD)
!-------------------------------------------------------------------
!This will transform the tensor from local into the global coordinate
!it used the transformation tensor to present the tensor that is
!represented in local coordinate in the local coordinate
!-------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------
      real(iwp),INTENT(IN) :: TRANSFORMATION(3,3) ! TRANSFORMATION TENSOR
      real(iwp),INTENT(IN) :: TENSOR_LOCAL_COORD(3,3) ! LOCAL COORDINATE TENSOR
      real(iwp),INTENT(OUT) :: TENSOR_GLOBAL_COORD(3,3) ! GLOBAL COORDINATE TENSOR
      INTEGER :: I,J,K,L! The counter integers
! ----------------------------------------------------------
      TENSOR_GLOBAL_COORD=0.0D0 !INITIALIZING TENSOR
! ----------------------------------------------------------
    DO  I=1,3
    DO  J=1,3
    DO  K=1,3
    DO  L=1,3
    TENSOR_GLOBAL_COORD(I,J)=TENSOR_GLOBAL_COORD(I,J)+TRANSFORMATION(I,K)*TRANSFORMATION(J,L)*TENSOR_LOCAL_COORD(K,L)
    ENDDO;ENDDO;ENDDO;ENDDO

      end subroutine TENSOR_TRANSFORM_RANK_2

!
!!===================================================================
!      SUBROUTINE ROTATION(XLOCAL,QTRANSFORM)
!!====================================================================
!      ! THE PhD THESIS OF AMIR SOHRABI MOLLAYOUSEF
!      ! Programmer(s):Amir Sohrabi
!      ! Objective: to find the rotation tensor for a given point with respect
!      ! to the origin of the problem. This transformation matrix is to be used
!      ! in non linear electromechancial analysis of telescopic actoator
!!=====================================================================
!      IMPLICIT NONE
!!=====================================================================
!      real(iwp) :: XLOCAL(3) ! THIS IS THE COORDINATE OF THE POINT OF INTEREST
!      real(iwp) :: QTRANSFORM(3,3) ! THIS IS THE TRANSFORMTION MATRIX
!      real(iwp) :: RADIOUS,EUCLI_NORM ! THIS IS THE TRANSFORMTION MATRIX
!      INTEGER::DIMENSIONV,INPUT_UNIT
!!=====================================================================
!      INPUT_UNIT=2
!!      OPEN(INPUT_UNIT,file='INPUT_UNIT.TXT')
!!      READ(INPUT_UNIT,*)XLOCAL(:)
!      QTRANSFORM=0.0D0
!      DIMENSIONV=3
!      RADIOUS=SQRT(DOT_PRODUCT(XLOCAL,XLOCAL))
!
!      IF (RADIOUS>0) THEN
!      QTRANSFORM(1,1)=XLOCAL(1)/RADIOUS
!      QTRANSFORM(1,2)=XLOCAL(2)/RADIOUS
!      QTRANSFORM(2,1)=-XLOCAL(2)/RADIOUS
!      QTRANSFORM(2,2)=XLOCAL(1)/RADIOUS
!      QTRANSFORM(3,3)=1
!      ENDIF
!
!
!!     ================================================================
!!                           FORMATS
!!     ================================================================
!  100 FORMAT(3E14.5)
!      end subroutine ROTATION

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
    end module linearsolvers
