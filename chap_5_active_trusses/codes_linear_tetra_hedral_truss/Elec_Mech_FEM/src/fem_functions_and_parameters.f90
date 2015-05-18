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

module fem_functions_and_parameters
    implicit none
      logical,parameter::updated_lagrangian=.false.
      integer:: i_calibration
      integer::gauss_point_number
      integer::ngauss
      integer::curved_node

      integer,parameter :: ndf=4 !! number of degrees of freedom
      integer,parameter :: dimen=3 !! dimention of problem
      integer,parameter::iwp=selected_real_kind(15)

      integer,parameter::out=2 !! output file
      integer,parameter::csv=11 !! comma seperated file
      integer,parameter::msh=3 !!  mesh file

      integer,parameter::vtu=13 !! output file
      integer,parameter::OUTPUT_UNIT=7
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

!     ===============================================================
!                       TRACE OF A 3X3 TENSOR
!     ===============================================================

      FUNCTION TRACE(A)
      real(iwp)::A(:,:),TRACE
      INTEGER::I
      TRACE=0.0D0;
      DO I=1,SIZE(A,DIM=1);TRACE=TRACE+A(I,I);ENDDO
      END FUNCTION TRACE


!! This is taken from http://www.netlib.org/napack/sort.f
!!      ________________________________________________________
!!     |                                                        |
!!     |            sort an array in increasing order           |
!!     |                                                        |
!!     |    input:                                              |
!!     |                                                        |
!!     |         x     --array of numbers                       |
!!     |                                                        |
!!     |         y     --working array (length  at least n)     |
!!     |                                                        |
!!     |         n     --number of array elements to sort       |
!!     |                                                        |
!!     |    output:                                             |
!!     |                                                        |
!!     |         x     --sorted array                           |
!!     |________________________________________________________|
!!
      subroutine sort(x)
      integer:: x(:),s,t
      integer,allocatable::y(:)
      integer i,j,k,l,m,n
      n=size(x,dim=1)
      y=x
      i = 1
10    k = i
20    j = i
      i = i + 1
      if ( j .eq. n ) goto 30
      if ( x(i) .ge. x(j) ) goto 20
      y(k) = i
      goto 10
30    if ( k .eq. 1 ) return
      y(k) = n + 1
40    m = 1
      l = 1
50    i = l
      if ( i .gt. n ) goto 120
      s = x(i)
      j = y(i)
      k = j
      if ( j .gt. n ) goto 100
      t = x(j)
      l = y(j)
      x(i) = l
60    if ( s .gt. t ) goto 70
      y(m) = s
      m = m + 1
      i = i + 1
      if ( i .eq. k ) goto 80
      s = x(i)
      goto 60
70    y(m)= t
      m = m + 1
      j = j + 1
      if ( j .eq. l ) goto 110
      t = x(j)
      goto 60
80    y(m) = t
      k = m + l - j
      i = j - m
90    m = m + 1
      if ( m .eq. k ) goto 50
      y(m) = x(m+i)
      goto 90
100   x(i) = j
      l = j
110   y(m) = s
      k = m + k - i
      i = i - m
      goto 90
120   i = 1
130   k = i
      j = x(i)
140   x(i) = y(i)
      i = i + 1
      if ( i .lt. j ) goto 140
      y(k) = i
      if ( i .le. n ) goto 130
      if ( k .eq. 1 ) return
      goto 40
      end subroutine sort


!! This is taken from http://www.netlib.org/napack/sort.f
!!      ________________________________________________________
!!     |                                                        |
!!     |            sort an array in increasing order           |
!!     |                                                        |
!!     |    input:                                              |
!!     |                                                        |
!!     |         x     --array of numbers                       |
!!     |                                                        |
!!     |         y     --working array (length  at least n)     |
!!     |                                                        |
!!     |         n     --number of array elements to sort       |
!!     |                                                        |
!!     |    output:                                             |
!!     |                                                        |
!!     |         x     --sorted array                           |
!!     |________________________________________________________|
!!
      subroutine sort_real(x)
      real(iwp):: x(:),s,t
      real(iwp),allocatable::y(:)
      integer i,j,k,l,m,n
      n=size(x,dim=1)
      y=x
      i = 1
10    k = i
20    j = i
      i = i + 1
      if ( j .eq. n ) goto 30
      if ( x(i) .ge. x(j) ) goto 20
      y(k) = i
      goto 10
30    if ( k .eq. 1 ) return
      y(k) = n + 1
40    m = 1
      l = 1
50    i = l
      if ( i .gt. n ) goto 120
      s = x(i)
      j = y(i)
      k = j
      if ( j .gt. n ) goto 100
      t = x(j)
      l = y(j)
      x(i) = l
60    if ( s .gt. t ) goto 70
      y(m) = s
      m = m + 1
      i = i + 1
      if ( i .eq. k ) goto 80
      s = x(i)
      goto 60
70    y(m)= t
      m = m + 1
      j = j + 1
      if ( j .eq. l ) goto 110
      t = x(j)
      goto 60
80    y(m) = t
      k = m + l - j
      i = j - m
90    m = m + 1
      if ( m .eq. k ) goto 50
      y(m) = x(m+i)
      goto 90
100   x(i) = j
      l = j
110   y(m) = s
      k = m + k - i
      i = i - m
      goto 90
120   i = 1
130   k = i
      j = x(i)
140   x(i) = y(i)
      i = i + 1
      if ( i .lt. j ) goto 140
      y(k) = i
      if ( i .le. n ) goto 130
      if ( k .eq. 1 ) return
      goto 40
      end subroutine sort_real


function identity_matrix(n_size)
real(iwp),allocatable::identity_matrix(:,:)
integer::n_size,i_counter

allocate(identity_matrix(n_size,n_size))
identity_matrix=0.0d0

do i_counter=1,n_size
    identity_matrix(i_counter,i_counter)=1.0d0
enddo

end function identity_matrix


function dyad_of_two_vector(x,y)
real(iwp),allocatable::dyad_of_two_vector(:,:)
real(iwp)::x(:),y(:)
integer::i_size,j_size,i_counter,j_counter

i_size=size(x,dim=1)
j_size=size(y,dim=1)

allocate(dyad_of_two_vector(i_size,j_size))

do i_counter=1,i_size
    do j_counter=1,i_size
    dyad_of_two_vector(i_counter,j_counter)=x(i_counter)*y(j_counter)
    enddo
enddo

end function dyad_of_two_vector







!!>This function finds the derivative of derivative of norm of a vector with respect to the vector

function derivative_of_direction(x)
real(iwp),allocatable::derivative_of_direction(:,:)
real(iwp)::x(:)
integer::i_size,i_counter,j_counter

i_size=size(x,dim=1)
allocate(derivative_of_direction(i_size,i_size))
derivative_of_direction=0.0d0

do i_counter=1,i_size
    do j_counter=1,i_size
    derivative_of_direction(i_counter,j_counter)=-x(i_counter)*x(j_counter)
    enddo
enddo



derivative_of_direction=derivative_of_direction/(norm_vect(x)**3.0d0)




    do j_counter=1,i_size
    derivative_of_direction(j_counter,j_counter)=derivative_of_direction(j_counter,j_counter)+1.0d0/norm_vect(x)
    enddo


if(norm_vect(x).eq.0)derivative_of_direction=0.0d0

end function derivative_of_direction


! ================================================================
!   polarization swithing function
!this will find the polarization with the given state of material
! ================================================================
subroutine direction_polarization(pr_vec,a_direc)
real(iwp)::a_direc(dimen),normed,pr_vec(dimen)
a_direc=0.0d0
normed=norm_vect(pr_vec)
if (normed.gt.0.0d0)then
a_direc=pr_vec/normed
endif



end subroutine direction_polarization


! ================================================================
!   polarization swithing function
!this will find the polarization with the given state of material
! ================================================================
function direction_of_vector(a)
real(iwp)::a(:)
real(iwp),allocatable::direction_of_vector(:)
real(iwp)::normed
allocate( direction_of_vector( size(a) ) );


direction_of_vector=0.0d0
normed=norm_vect(a)

if (normed.gt.0.0d0)then
direction_of_vector=a/normed
endif

end function direction_of_vector

!     ================================================================
!                 norm of a vector
!     ================================================================
      function norm_matrix(a)
      integer::i,j
      real(iwp)::a(:,:),norm_matrix


      norm_matrix=0.0d0
      do i=1,size(a,dim=1);do j=1,size(a,dim=2);
      norm_matrix=norm_matrix+a(i,j)*a(i,j)
      enddo;enddo


      end function norm_matrix

!     ================================================================
!                 macaulay brackets
!     ================================================================
      function macaulay(a)
      real(iwp)::macaulay
      real(iwp)::a
      macaulay=0
      if (a>=0)then
       macaulay=1
      endif
      end function macaulay

!     ================================================================
!                 macaulay brackets
!     ================================================================
      function macaulay_s(a)
      real(iwp)::macaulay_s
      real(iwp)::a
      macaulay_s=0.5*(abs(a)+a)

      end function macaulay_s
!     ================================================================
!                 macaulay brackets
!     ================================================================
      function macaulay_h(a)
      real(iwp)::macaulay_h
      real(iwp)::a
      macaulay_h=0.5*(tanh(a)+1)

      end function macaulay_h

!     ================================================================
!                 unit vector in the direction of a vector
!     ================================================================
      function unit_vect(a)
      real(iwp)::a(:)
      real(iwp),allocatable :: unit_vect(:)

      allocate(unit_vect(size(a)))
      unit_vect=0.0d0

      unit_vect=a/norm_vect(a)

      if(norm_vect(a).eq.0)unit_vect=0.0d0


      end function unit_vect

end module fem_functions_and_parameters
