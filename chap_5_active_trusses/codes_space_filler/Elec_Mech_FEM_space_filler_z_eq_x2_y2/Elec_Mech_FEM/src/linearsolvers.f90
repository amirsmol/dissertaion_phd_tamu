 !!  \brief     Electro Mechanical Finite Element Model
 !!  \author    Amir Sohrabi Mollayousef
 !!  \version   0.1a
 !!  \date      2011 - 2014
 !!  \copyright     This program is free software: you can redistribute it and/or modify
 !!    it under the terms of the GNU General Public License as published by
 !!    the Free Software Foundation, either version 3 of the License, or
 !!    (at your option) any later version.
 !!
 !!    This program is distributed in the hope that it will be useful,
 !!    but WITHOUT ANY WARRANTY; without even the implied warranty of
 !!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !!    GNU General Public License for more details.
 !!
 !!    You should have received a copy of the GNU General Public License
 !!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

module linearsolvers
use fem_functions_and_parameters
!use mkl95_lapack

contains


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



    end module linearsolvers
