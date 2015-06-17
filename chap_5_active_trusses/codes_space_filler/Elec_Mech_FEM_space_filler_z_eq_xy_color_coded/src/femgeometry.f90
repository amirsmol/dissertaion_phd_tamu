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
module fem_geometry
 use fem_functions_and_parameters
    implicit none
    real(iwp),allocatable::each_truss_strain(:)
    integer::npe !<number of nodes in each element
    integer::nem !<number of elements in the domain
    integer::nnm !<number of nodes in the domain
    integer::iel   !< element type tag that can be linear=1 and quadratic=2
    integer::ipdf  !< bumber of gauss pints in each direction
    integer::neq !< number of equations to be solved
    integer::nn !<number of degrees of freedom in an element
    integer::nhbw !<half band width
    integer::nw !<full band width
    integer::nbw !<full band width
    real(iwp),parameter  :: pi = 3.14159265 !<the pie number
    !  ================================================================
    ! geormetry variables
    !  ================================================================
    real(iwp),allocatable  ::coords(:,:) !<the nodes geometry array
    integer,allocatable ::nod(:,:) !<the connectivity matrix
    !  ================================================================
    ! boundary condition variables
    !  ================================================================
    integer,allocatable::ispv(:,:)!<the primary boundary condition array
    real(iwp),allocatable::vspv(:),vspvt(:)!<the primary boundary condition array
    integer,allocatable::issv(:,:)!<the secondry boundary condition
    real(iwp),allocatable::vssv(:),vssvt(:)!<the secondry boundary conditi

    integer,allocatable::bnd_no_pr_vec(:) !<the secondry boundary conditif
    real(iwp),allocatable::bnd_va_pr_vec(:) !<the secondry boundary conditi

    integer,allocatable::bnd_no_se_vec(:) !<the secondry boundary conditi
    real(iwp),allocatable::bnd_va_se_vec(:) !<the secondry boundary conditi

    integer::nspv,np,nb  !<number of specified primary boundary condition
    integer::nssv  !<number of specified secondry boundary condition

    integer,allocatable::afc_boundary_x_periodic(:) !<this holds the nodes on the periodic boundary conditions in x direction
    integer,allocatable::afc_boundary_y_periodic(:) !<this holds the nodes on the periodic boundary conditions in y direction
    integer,allocatable::afc_boundary_z_periodic(:) !<this holds the nodes on the periodic boundary conditions in z direction
    integer,allocatable::afc_negative_electr_z0_periodic(:) !<this holds the nodes on the periodic boundary conditions in y direction
    integer,allocatable::afc_negative_electr_zL_periodic(:) !<this holds the nodes on the periodic boundary conditions in z direction

    logical,allocatable::is_constrained(:) !<this holds the nodes on the periodic boundary conditions in z direction
    real(iwp),allocatable::constraint(:,:)
    logical::periodic_on


contains

subroutine truss_tetrahedral_space_filler(neldirectional,length,num_tetra_units)
        implicit none
        !     ================================================================
        !                           the variables
        !     ================================================================
        integer,parameter:: one=1,two=2,three=3 !number of points in the element
        real(iwp),parameter::  pi = 3.14159265 !the pie number
        !     ================================================================
        real(iwp)  ::length(dimen)
        integer::neldirectional(dimen)


        !     ============================geometrical variables
        real(iwp)::side_lenght
        real(iwp),dimension(3)::v1,v2,v3,vhnormal,planecenterpont, nextpoint

        real(iwp),dimension(3)::basecoodinate_1,basecoodinate_2 ,basecoodinate_3
        integer::base_node_num(4)
        integer,dimension(3)::base_elem,side_elem
        integer,dimension(4,3)::surface_nodes
        integer,dimension(6,2)::vertex_nodes
        integer::num_tetra_units,index_tetra_units

        !     =========================first and last node
        real(iwp),dimension(3)::firspoint,lastpoint
        real(iwp)::lenght_beam
        real(iwp),allocatable::coords_refed(:,:),total_coords(:,:),qtran(:,:)
        integer::ix,iy,iz
        integer::nx,ny,nz
        integer::i,j,k

        integer::ine
        integer::i_firstnode,i_secondnode

  integer :: numrows
  integer :: numcols

  integer, dimension(:,:), allocatable :: arraya
  logical, dimension(:) ,  allocatable ::  mask
  integer, dimension(:,:), allocatable :: arrayb
!  integer :: ix
  integer, dimension(:), allocatable :: index_vector


        !     ================================================================
        common/geometry_beam/firspoint,lastpoint
        !     ================================================================
       nem=3
       nnm=3
       npe=2



nx=7
ny=5
nz=3

nnm=nx*ny*nz

allocate(coords(nnm,dimen))

        coords=0.0d0
        do iz=0,nz-1
        do iy=0,ny-1
        do ix=0,nx-1
        coords(1+ix+iy*(nx)+iz*(ny*nx),:)=[length(1)*ix,length(1)*2.0d0*sqrt(3.0d0/4.0d0)*iy,length(1)*iz]
        end do
        end do
        end do


!!>fomring the total coordinate array

allocate( total_coords(nnm*2,dimen) )

total_coords(1:nnm,:)=coords(:,:)

FORALL(i=nnm+1:2*nnm) total_coords(i,:)=coords(i-nnm,:)+[length(1)*0.5d0,length(1)*sqrt(3.0d0/4.0d0),0.0d0]

deallocate(coords)
coords=total_coords
call show_matrix(total_coords,"total_coords")

nnm=2*nnm

nem=0
do i_firstnode=1,nnm
do i_secondnode=1,nnm
if ( (sqrt(norm_vect( coords(i_firstnode,:) -coords(i_secondnode,:) ))-length(1) ).le.length(1)*default_smallest_pivot) then
nem=nem+1
endif
end do
end do

write(*,*)nem

!nem=nem/2
allocate(nod(nem,npe))


ine=0
do i_firstnode=1,nnm
    do i_secondnode=1,nnm
        if ( (sqrt(norm_vect( coords(i_firstnode,:) -coords(i_secondnode,:) ))-length(1) ).le.length(1)*default_smallest_pivot) then
            ine=ine+1
            nod(ine,:)=[i_firstnode,i_secondnode]
        endif
    end do
end do


do ine=1,nem
    call sort(nod(ine,:))
end do


arraya=nod

numcols=size(arraya,dim=1)
numrows=size(arraya,dim=2)

ALLOCATE(mask(numcols))
mask = .TRUE.


arrayb=nod
arrayb=0
arrayb(1,:)=nod(1,:)

!mask=( arrayb(:,1)== nod(:,1) )


!write(*,*)
  ! make an index vector
!  allocate(index_vector, ountmask))
!
!  ! now copy the unique elements of a into b
!  allocate(arrayb, source=arraya(:,index_vector))

arrayb(1,:) = nod(1,:)
k = 1
  outer: do i=2,nem
     do j=1,k
        if ( arrayb(j,1) == nod(i,1).and.arrayb(j,2) == nod(i,2)) then
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     arrayb(k,:) = nod(i,:)
  end do outer

deallocate(nod)
allocate(nod(k,npe))
nod=arrayb(1:k,:)

end subroutine truss_tetrahedral_space_filler




subroutine truss_3d_space_filler_connectivity(length,neldirectional)
implicit none
!     ================================================================
!                           the variables
!     ================================================================
        integer,parameter:: one=1,two=2,three=3 !number of points in the element
        real(iwp),parameter::  pi = 3.14159265 !the pie number
!     ================================================================
real(iwp)::length(dimen)
real(iwp)::deltal(dimen)
real(iwp)::corner_point(dimen),length_1(dimen)
integer::nnddirectional(dimen),neldirectional(dimen)

integer::i,j,k
integer::inode
integer::nblock
integer::ielem

!integer::nblock,nnplane
!integer::inode,ix,iy,iz,i,j,k,id,in,iel,idn,ibl,nhbw,nw,n
!integer::nel_teuss,nel_teuss_tr,nel_teuss_pl,nroddx ! total number of truss element
!     ================================================================
       nnddirectional=neldirectional+1
       deltal=length/neldirectional
       npe=2
       nnm=1
       nblock=1
       do 7 i=1,dimen
       nblock=nblock*neldirectional(i)
7      nnm=nnm*nnddirectional(i)
!
       nem=0
       do 8 i=1,dimen
8      nem=nem+nnm*neldirectional(i)/nnddirectional(i)
       nem=nem+nblock*4

!       write(*,*)nem,npe,nnm,dimen

      allocate(nod(nem,npe),coords(nnm,dimen))

      coords=0.0d0
      nod=1
      select case(dimen)
      case(two)
!     ================================================================
      do 9 i=0,neldirectional(1)
      do 9 j=0,neldirectional(2)
      inode=(j)*nnddirectional(1)+i+1
      coords(inode,1)=i*deltal(1)
9     coords(inode,2)=j*deltal(2)

!     ============================x1 direction nod
      ielem=0
      do 12 j=0,nnddirectional(2)-1
      do 12 i=0,neldirectional(1)-1
      ielem=ielem+1
      nod(ielem,1)=1+i+j*nnddirectional(1)
12    nod(ielem,2)=2+i+j*nnddirectional(1)
!     ============================x2 direction nod
      do 13 j=0,nnddirectional(1)-1
      do 13 i=0,neldirectional(2)-1
      ielem=ielem+1
      nod(ielem,1)=1+i*nnddirectional(1)+j
13    nod(ielem,2)=1+(i+1)*nnddirectional(1)+j


!     ============================diagonal direction nod
      do 18 j=0,neldirectional(2)-1
      do 18 i=0,neldirectional(1)-1
      ielem=ielem+1
      nod(ielem,1)=1+i+j*nnddirectional(1)
18    nod(ielem,2)=2+i+nnddirectional(1)+j*nnddirectional(1)

      case(three)
!     ================================================================
      do 10 i=0,neldirectional(1)
      do 10 j=0,neldirectional(2)
      do 10 k=0,neldirectional(3)
      inode=k*nnddirectional(1)*nnddirectional(2)+j*nnddirectional(1)+i+1
      coords(inode,1)=i*deltal(1)
      coords(inode,2)=j*deltal(2)
10    coords(inode,3)=k*deltal(3)

!     ============================x1 direction nod
      ielem=0
      do 15 k=0,neldirectional(3)
      do 15 j=0,nnddirectional(2)-1
      do 15 i=0,neldirectional(1)-1
      ielem=ielem+1
      nod(ielem,1)=1+i+j*nnddirectional(1) +k*nnddirectional(1)*nnddirectional(2)
15    nod(ielem,2)=2+i+j*nnddirectional(1)+k*nnddirectional(1)*nnddirectional(2)
!     ============================x2 direction nod
      do 16 k=0,neldirectional(3)
      do 16 j=0,nnddirectional(1)-1
      do 16 i=0,neldirectional(2)-1
      ielem=ielem+1
      nod(ielem,1)=1+i*nnddirectional(1)+j+k*nnddirectional(1)*nnddirectional(2)
16    nod(ielem,2)=1+(i+1)*nnddirectional(1)+j+k*nnddirectional(1)*nnddirectional(2)
!     ============================x3 direction nod
      do 17 j=0,nnddirectional(2)-1
      do 17 i=0,nnddirectional(1)-1
      do 17 k=0,neldirectional(3)-1
      ielem=ielem+1
      nod(ielem,1)=1+k*nnddirectional(1)*nnddirectional(2)+i+j*nnddirectional(1)
17    nod(ielem,2)=1+(k+1)*nnddirectional(1)*nnddirectional(2)+i+j*nnddirectional(1)

!     ============================diagonal direction nod
      do 19 k=0,neldirectional(3)-1
      do 19 j=0,neldirectional(2)-1
      do 19 i=0,neldirectional(1)-1
      ielem=ielem+1
      nod(ielem,1)=1+i+j*nnddirectional(1)+k*nnddirectional(1)*nnddirectional(2)
19    nod(ielem,2)=2+i+nnddirectional(1)*nnddirectional(2) +nnddirectional(1)+&
       j*nnddirectional(1)&
      +k*nnddirectional(1)*nnddirectional(2)

      do 21 k=0,neldirectional(3)-1
      do 21 j=0,neldirectional(2)-1
      do 21 i=0,neldirectional(1)-1
      ielem=ielem+1
      nod(ielem,1)=2+i+j*nnddirectional(1) +k*nnddirectional(1)*nnddirectional(2)
21    nod(ielem,2)=1+i+nnddirectional(1)*nnddirectional(2) +nnddirectional(1)+j*nnddirectional(1)&
      +k*nnddirectional(1)*nnddirectional(2)

!
      do 22 k=0,neldirectional(3)-1
      do 22 j=0,neldirectional(2)-1
      do 22 i=0,neldirectional(1)-1
      ielem=ielem+1

      nod(ielem,1)=1+i+(j+1)*nnddirectional(1) +k*nnddirectional(1)*nnddirectional(2)
22    nod(ielem,2)=2+i+(j)*nnddirectional(1) +(k+1)*nnddirectional(1)*nnddirectional(2)

      do 23 k=0,neldirectional(3)-1
      do 23 j=0,neldirectional(2)-1
      do 23 i=0,neldirectional(1)-1
      ielem=ielem+1

      nod(ielem,1)=2+i+(j+1)*nnddirectional(1) +k*nnddirectional(1)*nnddirectional(2)
23    nod(ielem,2)=1+i+(j)*nnddirectional(1) +(k+1)*nnddirectional(1)*nnddirectional(2)
      end select

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length_1(i)=maxval(coords(:,i))-minval(coords(:,i)); enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length_1(i);  enddo

end subroutine truss_3d_space_filler_connectivity



!< This subroutine define and producd the boundary conditions for space filler tsuss !>
!! @param nspv number of specified primary bounday conditions
!! @param nsvv number of specified secondry bounday conditions
!! @param ispv(i boundady,2) this array contains the information for the boundary conditions
!!        i boundary the number of boundary conditions and ispv(i boundady,1) is the node number
!!        i boundary the number of boundary conditions and ispv(i boundady,2) is the degree of freedom number
!! @todo apply the periodic boundary conditions
subroutine truss_3d_space_filler_boundary_z_eq_xy(length)
real(iwp)::length(dimen)
integer::corner_node
integer::inode,ibond

corner_node=0
!!< the space filler truss is going to be bounded
!! if x_1=0,x_2=0,x_3=-L_3/2 then u_1=0, u_2=0, u_3=0, this will restrain it from movment in x_1, x_2, x_3
!! if x_1=0,x_2=0,x_3=+L_3/2 then u_1=0, u_2=0,  this will restrain it from rotation around x_1, x_2 axis
!! if x_1=0,x_2=L_2/2,x_3=+L_3/2 then u_1=0, this will restrain it from rotation around x_3 axis


nspv=6

! write(*,*) nspv
!
!write(*,*)'nspv'
!write(*,*)nspv
!
!write(*,*)'ispv'
!write(*,*)ispv

!!< Manually forming the bounrary condition index arrays
!!  this will rotate trhough all the nodes and if the nodes
!!  matches the bounday condition
!!  it will add it to the index arrays and assigns values to the
!!  boundary condition

nssv=1

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))

vspv=0.0d0
vssv=0.0d0

ibond=0
do inode=1,nnm

!! pinning the x_1=-l/2 direction

if(coords(inode,1)==0 .and. coords(inode,2)==0 .and.coords(inode,3)==-0.5*length(3) )then
 ibond=ibond+1;

 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;

ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;

 ibond=ibond+1;

 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;

endif


if(coords(inode,1)==0 .and. coords(inode,2)==0 .and.coords(inode,3)==0.5*length(3) )then

 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;

 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;

endif


if(coords(inode,1)==0.and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
corner_node=inode
endif

enddo

call show_matrix_int((ispv),'ispv')

!! // issv(1,:)=[node number , degree of freedom number];

issv(1,:)=[corner_node,3];
vssv(1)=0.0e0

vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb

    enddo



end subroutine truss_3d_space_filler_boundary_z_eq_xy


subroutine truss_3d_space_filler_boundary(length)
real(iwp)::length(dimen)
integer::corner_node
integer::inode,ibond

corner_node=0
!!< the space filler truss is going to be bounded
!! if the x_1=-l/2 then u_1=0 this will restrain it from movment in x_1 direction and rotation in x_3 and x_2
!! if the x_1=-l/2 and x_3=-l/2 then u_3=0 this will restrain it from movment in x_2 direction and also rotation in x_1 axis
!! if the x_3=-l/2  then u_2=0 this will restrain it from movment in x_3
!! if the fourth degrees of freedom for all truss

nspv=count(coords(:,1)==-0.5*length(1))+ &
     count(coords(:,2)==-0.5*length(2)) + &
     count((coords(:,1)==-0.5*length(1)).and.(coords(:,3)==-0.5*length(3))) + &
     nnm

! write(*,*) nspv
!
!write(*,*)'nspv'
!write(*,*)nspv
!
!write(*,*)'ispv'
!write(*,*)ispv

!!< Manually forming the bounrary condition index arrays
!!  this will rotate trhough all the nodes and if the nodes
!!  matches the bounday condition
!!  it will add it to the index arrays and assigns values to the
!!  boundary condition

nssv=1

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))

vspv=0.0d0
vssv=0.0d0

ibond=0
do inode=1,nnm

!! pinning the x_1=-l/2 direction

if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif


!! pinning the x_2=-l/2

if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
endif


!! pinning the intersection of x_1=-l/2 and x_3=-l/2 planes in x_3 direction

if( coords(inode,1)==-0.5*length(1) .and. coords(inode,3)==-0.5*length(3)  )then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;
endif
!
!! putting the electric potential for all the nodes to be zero

 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0;


if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
corner_node=inode
endif

enddo




!write(*,*)ispv
!! // issv(1,:)=[node number , degree of freedom number];

issv(1,:)=[corner_node,3];
vssv(1)=0.0e0

vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb

    enddo

end subroutine truss_3d_space_filler_boundary


!    subroutine make_element_direction_normal(two_point)
!        implicit none
!        real(iwp)::two_point(:,:),lent_ele
!        real(iwp)::lent_prog(dimen)
!
!        lent_prog(:)= two_point(2,:)-two_point(1,:)
!        lent_ele=dsqrt(dot_product(lent_prog,lent_prog))
!
!        lent_ele=norm_vect(lent_prog)
!        element_direction_normal=lent_prog/lent_ele
!    end subroutine make_element_direction_normal


subroutine truss_actua_3d_connectivity(length,num_tetra_units)
        implicit none
        !     ================================================================
        !                           the variables
        !     ================================================================
        integer,parameter:: one=1,two=2,three=3 !number of points in the element
        real(iwp),parameter::  pi = 3.14159265 !the pie number
        !     ================================================================
        real(iwp)  ::length(dimen)


        !     ============================geometrical variables
        real(iwp)::side_lenght
        real(iwp),dimension(3)::v1,v2,v3,vhnormal,planecenterpont, nextpoint

        real(iwp),dimension(3)::basecoodinate_1,basecoodinate_2 ,basecoodinate_3
        integer::base_node_num(4)
        integer,dimension(3)::base_elem,side_elem
        integer,dimension(4,3)::surface_nodes
        integer,dimension(6,2)::vertex_nodes
        integer::num_tetra_units,index_tetra_units

        !     =========================first and last node
        real(iwp),dimension(3)::firspoint,lastpoint
        real(iwp)::lenght_beam
        real(iwp),allocatable::coords_refed(:,:),total_coords(:,:),qtran(:,:)
        integer::i
        !     ================================================================
        common/geometry_beam/firspoint,lastpoint
        !     ================================================================
       nem=num_tetra_units*3+3
       nnm=num_tetra_units+3
       npe=2


        allocate(nod(nem,npe),coords(nnm,dimen))
        allocate(coords_refed(nnm,dimen),total_coords(2,dimen), qtran(dimen,dimen))
!
        coords=0.0d0
        nod=0

        coords(1,:)=[length(1),0.0d0,0.0d0]
        coords(2,:)=matmul(rotate_deg_rz(120.0d0),coords(1,:))
        coords(3,:)=matmul(rotate_deg_rz(120.0d0),coords(2,:))

        base_node_num(1)=1
        base_node_num(2)=2
        base_node_num(3)=3
        base_node_num(4)=base_node_num(3)+1

        base_elem=[1,2,3]
        side_elem=[1,2,3]+3
!!
        surface_nodes(1,:)=base_node_num(1:3)
        surface_nodes(2,:)=base_node_num(2:4)
        surface_nodes(3,:)=[base_node_num(1),base_node_num(3:4)]

        vertex_nodes(1,:)=base_node_num([1,2])
        vertex_nodes(2,:)=base_node_num([2,3])
        vertex_nodes(3,:)=base_node_num([3,1]) ! [base_node_num(3),base_node_num(1)]
        vertex_nodes(4,:)=base_node_num([1,4]) ![base_node_num(1),base_node_num(4)]
        vertex_nodes(5,:)=base_node_num([2,4]) ! [base_node_num(2),base_node_num(4)]
        vertex_nodes(6,:)=base_node_num([3,4]) !base_node_num(3:4)
!!
!!
        nod(1:3,:)= vertex_nodes(1:3,:)

        do index_tetra_units=1,num_tetra_units

            !     ==========================connectivity
            base_node_num(4)= index_tetra_units+3
            base_node_num(3)=base_node_num(4)-1
            base_node_num(2)=base_node_num(4)-2
            base_node_num(1)=base_node_num(4)-3


            surface_nodes(1,:)=base_node_num(1:3)
            surface_nodes(2,:)=base_node_num(2:4)
            surface_nodes(3,:)=[base_node_num(1),base_node_num(3:4)]


            vertex_nodes(1,[1,2])=base_node_num([1,2])
            vertex_nodes(2,[1,2])=base_node_num([2,3]) ! base_node_num(2:3)
            vertex_nodes(3,[1,2])=base_node_num([3,1]) ! [base_node_num(3),base_node_num(1)]
            vertex_nodes(4,[1,2])=base_node_num([1,4]) ! [base_node_num(1),base_node_num(4)]
            vertex_nodes(5,[1,2])=base_node_num([2,4]) ! [base_node_num(2),base_node_num(4)];
            vertex_nodes(6,[1,2])=base_node_num([3,4]) ! [base_node_num(3),base_node_num(4)];

             nod(index_tetra_units*3+[1,2,3],:)=vertex_nodes([4,5,6],:)
!
            basecoodinate_1=coords(base_node_num(1),:);
            basecoodinate_2=coords(base_node_num(2),:);
            basecoodinate_3=coords(base_node_num(3),:);


            v1= basecoodinate_2-basecoodinate_1
            v2= basecoodinate_3-basecoodinate_2
            v3= basecoodinate_1-basecoodinate_3

            side_lenght=sqrt(dot_product(v1,v1))
            vhnormal=cross3n(v2,v3)

            planecenterpont=       &
                (basecoodinate_1+basecoodinate_2+basecoodinate_3)/3.0
            nextpoint= planecenterpont+vhnormal*side_lenght*sqrt(2.0/3.0)

            coords(base_node_num(1),:)=basecoodinate_1;
            coords(base_node_num(2),:)=basecoodinate_2;
            coords(base_node_num(3),:)=basecoodinate_3;
            coords(base_node_num(4),:)=nextpoint;


        enddo

        firspoint=sum(coords(1:4,:), dim = 1)
        firspoint=firspoint/4.0


        lastpoint=sum(coords(base_node_num(1:4),:), dim = 1)
        lastpoint=lastpoint/4.0

        total_coords(1,:)=firspoint
        total_coords(2,:)=lastpoint

        lenght_beam=sqrt(dot_product(firspoint-lastpoint,firspoint-lastpoint))
        call generate_trans_tensor(total_coords,qtran)
        !
        do i=1,nnm
            coords_refed(i,:)= matmul(coords(i,:),qtran)
        enddo


        coords= coords_refed/lenght_beam
!        height_of_beam=maxval(coords(:,2))-minval(coords(:,2))

!       write(2,*)height_of_beam,lenght_beam,maxval(coords(:,1))

    end subroutine truss_actua_3d_connectivity

!     ================================================================
!                           the para view writer subroutine
!     ================================================================
  subroutine truss_paraview_3d_vtu_xml_writer(glu)
      implicit none
!     ================================================================
!                           the variables
!     ================================================================
      real(iwp),intent(in)::glu(:)
      real(iwp),allocatable::coords_deformed(:,:),coordst(:,:)
!     real(iwp)::ele_pol_vector(:,:)  !ele_pol_vector(nem,dimen)
      integer::ipdf
      integer::per_element_gauss

!     ================================================================
!                      internal variables
!     ================================================================
      integer,save::counter
!     ===============================================================
      integer :: i,cell_types
      character(len=5)::x1
      character(len=20) :: fmt,filename ! format descriptor
      integer::inode,inem
      real(iwp)::scale
      real(iwp) :: each_truss_voltage(nem)
!     ===============================================================
      fmt = '(i4.4)' ! an integer of width 5 with zeros at the left
!     ===============================================================
      allocate(coords_deformed(nnm,dimen),coordst(nnm,dimen))


      if(npe.eq.20)then;cell_types=25;ipdf=2;endif
      if(npe.eq.8)then;cell_types= 12;ipdf=1;endif
      if(npe.eq.2)then;cell_types= 3 ;ipdf=1;endif

      per_element_gauss =ipdf**dimen
      each_truss_voltage=each_truss_strain/0.340E-9;

!      write(*,*)'maximum voltage=',MAXVAL(each_truss_voltage),MAXLOC(each_truss_voltage)
!      write(*,*)'minumum voltage=',minval(each_truss_voltage),minloc(each_truss_voltage)
!
!
!      write(*,*)'maximum strain=',MAXVAL(each_truss_strain),MAXLOC(each_truss_strain)
!      write(*,*)'minumum strain=',minval(each_truss_strain),minloc(each_truss_strain)
!      write(*,*)'maximum deflection=',MAXVAL(abs(glu)),MAXLOC(glu)


      counter=counter+1
      write(x1,fmt)counter
      coordst=0.0d0
      coords_deformed=0.0d0
      scale=1.0

      coordst=coords

      do i=1,nnm
      coords_deformed(i,:)=coords(i,:)+glu((i-1)*ndf+[1,2,3]);
      enddo !i=1,nnm

!      coords_deformed=coords;

!      coords_deformed

!      forall(i=1:nnm) coords_deformed([1,2,3],i)=coordst([1,2,3],i)+glu((i-1)*ndf+[1,2,3])

      filename='vtu/outpar'//trim(x1)//'.vtu'
      write(*,*)filename
      !write(*,*)nnm
      open (vtu,file=filename)
      write(vtu,*)'<VTKFile  type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
      write(vtu,*)'<UnstructuredGrid>'
      write(vtu,*)'<Piece NumberOfPoints="',nnm,'" NumberOfCells="',nem,'">'

    WRITE(vtu,*)('<CellData Scalars="strain">')

        WRITE(vtu,*)'<DataArray type="Float32" Name="strain" NumberOfComponents="1" format="ascii">'
        WRITE(vtu,105)each_truss_strain;
        WRITE(vtu,*)'</DataArray>'
!
      WRITE(vtu,*)'<DataArray type="Float32" Name="Electric Field (M V/m)" NumberOfComponents="1" format="ascii">'
      WRITE(vtu,105)(each_truss_voltage/1e9);
      WRITE(vtu,*)'</DataArray>'
!
    WRITE(vtu,*)('</CellData>')


    write(vtu,*)'<PointData Scalars="elec_pot" Vectors="displacement">'
    write(vtu,*)' <DataArray type="Float32" Name="elec_pot" NumberOfComponents="1" Format="ascii">'
    write(vtu,*)(glu((i-1)*ndf+4),i=1,nnm);
    write(vtu,*)' </DataArray >'
    write(vtu,*)' <DataArray  type="Float32" Name="displacement" NumberOfComponents="3" Format="ascii">'
    write(vtu,fmt='(a3,3(e16.5))')'  ',(glu((i-1)*ndf+1:(i-1)*ndf+3),i=1,nnm);
    write(vtu,*)' </DataArray>'
    write(vtu,*)'</PointData>'
!
!
      write(vtu,*)'<Points>'
      write(vtu,*)' <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
      write(vtu,fmt='(a3,3(e16.5))')'  ',(coords_deformed(inode,:),inode=1,nnm)
      write(vtu,*)' </DataArray>'
      write(vtu,*)'</Points>'

      write(vtu,*)'<Cells>'
      write(vtu,*)' <DataArray type="Int32" Name="connectivity" Format="ascii">'
      write(vtu,fmt='(2(i5))')(nod(inem,:)-1,inem=1,nem)
      write(vtu,*)' </DataArray>'

      write(vtu,*)'<DataArray type="Int32" Name="offsets" Format="ascii">'
      write(vtu,*)(inem*npe,inem=1,nem)
      write(vtu,*)'</DataArray>'

      write(vtu,*)'<DataArray type="Int32" Name="types" Format="ascii">'
      write(vtu,*)(cell_types,inem=1,nem)

      write(vtu,*)'</DataArray>'
      write(vtu,*)'</Cells>'
      write(vtu,*)'</Piece>'
      write(vtu,*)'</UnstructuredGrid>'
      write(vtu,*)'</VTKFile>'



      close (vtu)

105   format(3(e16.5))

      end subroutine truss_paraview_3d_vtu_xml_writer




subroutine truss_tetra_hedral_bounday()

nspv=6; !number of nodes, number of elements
nssv=0

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
ispv(:,1)=[1,2,3,1,1,2]
ispv(1:3,2)=3
ispv(4,2)=1
ispv(5,2)=2
ispv(6,2)=1




vspv=0.0d0
vssv=0.0d0
vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb
    enddo

!write(*,*)bnd_va_pr_vec
!write(*,*)vssv

end subroutine truss_tetra_hedral_bounday


subroutine linear_truss_bending_boundary()
integer::i

nspv=nnm+3; !number of nodes, number of elements
nssv=1

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))

ispv(1:nnm,1)=(/(i, i=1, nnm)/)
ispv(1:nnm,2)=2

ispv(nnm+1,:)=[1,1]
ispv(nnm+2,:)=[1,3]
ispv(nnm+3,:)=[nnm,3]



issv(1,:)=[1,1]

vspv=0.0d0
vssv=0.0d0

vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb
    enddo

!write(*,*)bnd_va_pr_vec
!write(*,*)vssv

end subroutine linear_truss_bending_boundary



subroutine linear_truss_bending_boundary_2()
integer::i

nspv=10; !number of nodes, number of elements
nssv=1

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))


ispv(1,:)=[1,1]
ispv(2,:)=[1,2]
ispv(3,:)=[1,3]
ispv(4,:)=[2,1]
ispv(5,:)=[2,2]
ispv(6,:)=[2,3]


ispv(7,:)=[nnm,2]
ispv(8,:)=[nnm,3]

ispv(9,:)=[nnm-1,2]
ispv(10,:)=[nnm-1,3]





issv(1,:)=[1,1]

vspv=0.0d0
vssv=0.0d0

vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb
    enddo

!write(*,*)bnd_va_pr_vec
!write(*,*)vssv

end subroutine linear_truss_bending_boundary_2

subroutine truss_tetra_hedral_displacement_control()
integer::i

nspv=nnm+2; !number of nodes, number of elements
nssv=2

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))

ispv(1:nnm,1)=(/(i, i=1, nnm)/)
ispv(1:nnm,2)=3

ispv(nnm+1,:)=1

ispv(nnm+2,:)=[1,2]


vspv=0.0d0
vssv=0.0d0

issv(1,:)=[nnm,2]
issv(2,:)=[nnm,1]



vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb
    enddo

!write(*,*)bnd_va_pr_vec
!write(*,*)vssv

end subroutine truss_tetra_hedral_displacement_control














subroutine truss_single_element_boundary()

nspv=5+2; !number of nodes, number of elements
nssv=1

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))

ispv(1,:)=[1,1]
ispv(2,:)=[1,2]
ispv(3,:)=[1,3]
ispv(4,:)=[1,4]

ispv(5,:)=[2,2]
ispv(6,:)=[2,3]
ispv(7,:)=[2,4]


issv(1,:)=[2,1];
vssv(1)=1.0

vspv=0.0d0

vspvt=vspv
vssvt=vssv

    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

    do np=1,nspv
        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
        bnd_va_pr_vec(np)=vspv(np)
        bnd_no_pr_vec(np)=nb
    enddo

    do np=1,nssv
        nb=(issv(np,1)-1)*ndf+issv(np,2)
        bnd_va_se_vec(np)=vssv(np)
        bnd_no_se_vec(np)=nb

    enddo

!write(*,*) bnd_va_pr_vec
!write(*,*) bnd_no_pr_vec

end subroutine truss_single_element_boundary








subroutine form_constraint_bimorph()

real(iwp),allocatable::constraint_x(:,:),constraint_y(:,:),constraint_z(:,:) !<the constraint matrix
real(iwp),allocatable::full_constraint(:,:)
logical,allocatable::fixed_boundary_mask(:)
integer::n_constraint !<number of constraintd
integer::i_const,j_const,const_num,n_fixed_boundary
integer::i_boundary



n_fixed_boundary=count(vspv(:)==0)


fixed_boundary_mask=vspv(:)==0

!
!write(*,*)'n_fixed_boundary',n_fixed_boundary
!write(1,*)'vspv',vspv
!write(1,*)'fixed_boundary_mask',fixed_boundary_mask

n_constraint=n_fixed_boundary

write(*,*)'n_constraint=',n_constraint
!! we will add three sets of periodic boundary condistions
!! so we have to add one to the dimension of constraint for each set
allocate(constraint(neq,neq-n_constraint),is_constrained(neq))
allocate(full_constraint(neq,neq))


full_constraint=0
do j_const=1,neq; full_constraint(j_const,j_const)=1; enddo


is_constrained=.false.

!write(*,*)afc_boundary_x_periodic



do i_boundary=1,size(vspv(:))
    if (vspv(i_boundary)==0) then
    i_const=(ispv(i_boundary,1)-1)*ndf+ispv(i_boundary,2)
    full_constraint(i_const,:)=0;
    endif
enddo

constraint=0
i_const=0

do j_const=1,neq
if (   sum(  full_constraint(:,j_const) ) .gt. 0 ) then
 i_const=i_const+1
 constraint(:,i_const)=full_constraint(:,j_const);
endif
enddo

end subroutine form_constraint_bimorph




!  ================================================================
!!>  This subroutine will apply the  boundary conditions for bimorph beam
!!   The bimorph beam is clamped at x=0
!!   The middle section is under 1 volt potential
!!   The upped and lower surfaces of the beam are under zero potential
!  ================================================================
subroutine bimorph_beam_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode
! integer::curved_node
integer::np,nb

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction


!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! body of the subroutine
!  ================================================================
coords=coords*1.0;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

!! At x_1=0 we will have u_1=0,
!! At x_2=0 we will have u_2=0;
!! At x_3=0 and x_1= we will have u_3=0

write(*,*)'length=',length
nspv=count( coords(:,1) ==- 0.5*length(1) )+ &
     nnm+ &
     count( (coords(:,3) ==- 0.5*length(3) ).and.( coords(:,1)==-0.5*length(1) )  )+ &
     count( coords(:,3) == -0.5*length(3) )+ &
     count( coords(:,3) == +0.5*length(3) )+ &
     count( coords(:,3) ==  0.0 )

nssv=1
!
write(*,*)'nspv',nspv
allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0

ix_periodic=0
iy_periodic=0
iz_periodic=0

ibond=0
do inode=1,nnm

if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif
!
!if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
! endif
!
if( (coords(inode,3) ==-0.5*length(3) ).and.( coords(inode,1)==-0.5*length(1) )  )then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;
endif
!
if(coords(inode,3)==0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0;
endif

if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0;
endif

if(coords(inode,3)==-0.0d0*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=1.0;
endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.0*length(2).and.coords(inode,3)==0.0*length(3))then
curved_node=inode
endif
enddo

! ispv(nspv,1)=curved_node;
! ispv(nspv,2)=4;
! vspv(nspv)=1.0;


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo
write(1,*)'vspv',vspv
end subroutine bimorph_beam_boundary




! ================================================================
      subroutine beam_mesh_file_reader()
      implicit none
! ================================================================
!                          input variables
! ================================================================

! ================================================================
!                          output variables
! ================================================================
!      integer,allocatable::nod(:,:)  !nod(nem,npe) !nod(mxelm,mxnod)
!      real*8,allocatable::coords(:,:) !coords(nnm,dimen) the nodes geometry array
! ================================================================
!                          trivial variables
! ================================================================
      integer::n_node,n_elem,n_npe,node_num,n_dimen
      integer::nem,nnm
      integer::nele,i
      real(iwp)::corner_point(dimen),length_1(dimen)
! ================================================================
iel=2;npe=20

     open (msh,file='meshes/3d_bimorh_quadratic_mesh.msh')
     read(msh,*)nnm,nem; !number of nodes, number of elements


     allocate(nod(nem,npe),coords(nnm,dimen))



      nod=0
      coords=0.0d0
      write(*,*)'nem=',nem
      write(*,*)'nnm=',nnm

      do n_node=1,nnm
        read(msh,*)node_num,(coords(n_node,n_dimen),n_dimen=1,dimen);
      enddo;

!     reading the connectivity array from the mesh file
      read(msh,*)
      do n_elem=1,nem
        read(msh,*)nele,(nod(n_elem,n_npe),n_npe=1,npe);
      enddo


coords=0.0010d0*coords

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length_1(i)=maxval(coords(:,i))-minval(coords(:,i)); enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5d0*length_1(i);  enddo

write(*,*)'length_1',length_1


      end subroutine beam_mesh_file_reader

      subroutine quadratic_c3d20_3d_mesh_geometry(neldirectional,length)
!     ================================================================
!                          input variables
!     ================================================================
      real*8::length(:)
!     ===================================================================
!                       output variables
!     ===================================================================

      integer::neldirectional(:)
!     ===================================================================
!                       trivial variables
!     ===================================================================
      integer::inode
      real*8,allocatable::deltal(:)
      integer,allocatable::nnddirectional(:)
!     ===================================================================

!     ===================================================================
!                       trivial variables
!     ===================================================================
      integer::ix,iy,iz,i,j,k
      integer::node_1,i_ele,nblock
!     ===================================================================
       iel=2
       npe=20

      allocate(deltal(dimen),nnddirectional(dimen))

       nnddirectional=neldirectional+1
       nnm=1
       nblock=product(neldirectional)
       nnm=product(nnddirectional)
       nem=nblock

       nem=0
       do  i=1,dimen
       nem=nem+nnm*neldirectional(i)/nnddirectional(i)
       enddo
       nnm=nem+nnm
       nem=nblock


       nnddirectional=neldirectional+1
       deltal=length/neldirectional
!c
!c        cubical elements
!c
      allocate(nod(nem,npe),coords(nnm,dimen))

       coords=0.0d0
       inode=0
 !     ===================================nodes in the corner
          do  k = 1, nnddirectional(3)
          do  j = 1, nnddirectional(2)
          do  i = 1, nnddirectional(1)
         inode=inode+1
         coords(inode,:) =   [deltal(1)*(i-1),deltal(2)*(j-1),deltal(3)*(k-1)]
          enddo
          enddo
          enddo

 !     ===================================nodes on x direction bars
          do  k = 1, nnddirectional(3)
          do  j = 1, nnddirectional(2)
          do  i = 1, neldirectional(1)
         inode=inode+1
         coords(inode,:) =   [deltal(1)*(i-1)+deltal(1)*0.5,deltal(2)*(j-1),deltal(3)*(k-1)]

          enddo
          enddo
          enddo

 !     ===================================nodes on bars on y direction bars
          do  k = 1, nnddirectional(3)
          do  i = 1, nnddirectional(1)
          do  j = 1, neldirectional(2)
         inode=inode+1
         coords(inode,:) =   [deltal(1)*(i-1),deltal(2)*(j-1)+deltal(2)*0.5,deltal(3)*(k-1)]

          enddo
          enddo
          enddo

 !     ===================================nodes on on bars on z direction bars
          do  k = 1, neldirectional(3)
          do  j = 1, nnddirectional(2)
          do  i = 1, nnddirectional(1)
         inode=inode+1
         coords(inode,:) =    [deltal(1)*(i-1),deltal(2)*(j-1),deltal(3)*(k-1)+deltal(3)*0.5]

          enddo
          enddo
          enddo
!     ===================================the connectivity matrix
!     ===================================corner nodes
          nod(1,1)=1
          nod(1,2)= nod(1,1)+1
          nod(1,3)= nod(1,2)+nnddirectional(1)
          nod(1,4)= nod(1,3)-1
          nod(1,5:8)=nod(1,1:4)+nnddirectional(1)*nnddirectional(2)


       node_1=0
       nod=0
      do ix=0,neldirectional(1) -1
      do iy=0,neldirectional(2) -1
      do iz=0,neldirectional(3) -1
      i_ele=1+ix+iy*neldirectional(1)+ iz*neldirectional(1)*neldirectional(2)

      nod(i_ele,1)=1 +ix+iy*nnddirectional(1)+           iz*nnddirectional(1)*nnddirectional(2)
      nod(i_ele,2)=nod(i_ele,1)+1
      nod(i_ele,3)=nod(i_ele,2)+nnddirectional(1)
      nod(i_ele,4)=nod(i_ele,3)-1
      nod(i_ele,5:8)=nod(i_ele,1:4)+nnddirectional(1)*nnddirectional(2)
 !     ===================================nodes on x and y direction bars
      nod(i_ele,9)=product(nnddirectional)+ix+iy*neldirectional(1) +iz*nnddirectional(2)*neldirectional(1)+1

      nod(i_ele,11)=nod(i_ele,9)+neldirectional(1)



      nod(i_ele,12)=product(nnddirectional)+ nnddirectional(1)*neldirectional(2)*nnddirectional(3)+ &
      iy+ix*neldirectional(2)+iz*neldirectional(2)*nnddirectional(3)  +1

       nod(i_ele,10)=nod(i_ele,12)+neldirectional(2)

 !     ===================================nodes on x and y direction bars  ob top shelf
       nod(i_ele,13)=nod(i_ele,9)+ nnddirectional(2)*neldirectional(1)

       nod(i_ele,16)=nod(i_ele,12)+ neldirectional(1)*nnddirectional(2)

       nod(i_ele,15)=nod(i_ele,13)+neldirectional(1)
       nod(i_ele,14)=nod(i_ele,16)+neldirectional(2)


!     !number of bars in xy plane

 !     ===================================nodes on z and bars
      nod(i_ele,17)=product(nnddirectional)+                          &
         neldirectional(1) * nnddirectional(2) *nnddirectional(3)    & !bars in x direction
      +  nnddirectional(1) * neldirectional(2) *nnddirectional(3)    & !bars in y direction
      +  ix                                                          &
      +  iy*nnddirectional(1)                                       &
      +  iz*nnddirectional(2)*nnddirectional(1)+1

      nod(i_ele,18)= nod(i_ele,17)+1

      nod(i_ele,19)=nod(i_ele,18)+nnddirectional(1)
      nod(i_ele,20)=nod(i_ele,17)+nnddirectional(1)

      enddo
      enddo
      enddo




!
      end subroutine quadratic_c3d20_3d_mesh_geometry
    !  ================================================================
    !   boundary conditions
    !  ================================================================
    subroutine linder_miehe_boundary(length)
        !  ================================================================
        ! input geormetry variables
        !  ================================================================
        implicit none
        !  ================================================================
        ! input geormetry variables
        !  ================================================================
        real*8:: length(:)
        integer::ibond,inode
        ! integer::curved_node
        integer::np,nb
        !  ================================================================
        real*8::corner_point(dimen)
        integer::i ! ,j,k
        !  ================================================================
        ! body of the subroutine
        !  ================================================================
        coords=coords*1.0;

        do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
            do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
                do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
                enddo
                do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

                    nspv=count(coords(:,1)==-0.5*length(1))+ &
                        count(coords(:,2)==-0.5*length(2))+ &
                        2*count(coords(:,3)==-0.5*length(3))+ &
                        count(coords(:,3)== 0.5*length(3))

                    nssv=1
                    !
                    !write(*,*)nspv
                    allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
                    allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
                    vssv=0.0d0
                    issv(1,1)=1
                    issv(1,2)=1
                    !
                    !!manually forming the bounrary condition
                    !
                    vspv=0.0d0

                    ibond=0
                    do inode=1,nnm

                        if(coords(inode,1)==-0.5*length(1))then
                            ibond=ibond+1;
                            ispv(ibond,1)=inode;
                            ispv(ibond,2)=1;
                            vspv(ibond)=0.0d0;
                        endif
                        !
                        if(coords(inode,2)==-0.5*length(2))then
                            ibond=ibond+1;
                            ispv(ibond,1)=inode;
                            ispv(ibond,2)=2;
                            vspv(ibond)=0.0d0;
                        endif
                        !
                        if(coords(inode,3)==-0.5*length(3))then
                            ibond=ibond+1;
                            ispv(ibond,1)=inode;
                            ispv(ibond,2)=3;
                            vspv(ibond)=0.0d0;

                            ibond=ibond+1;
                            ispv(ibond,1)=inode;
                            ispv(ibond,2)=4;
                            vspv(ibond)=-1.0d0;

                        endif
                        !
                        if(coords(inode,3)==0.5*length(3))then
                            ibond=ibond+1;
                            ispv(ibond,1)=inode;
                            ispv(ibond,2)=4;
                            vspv(ibond)=+1.0;
                        endif

                        if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
                            curved_node=inode
                        endif
                    enddo

                    allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
                    allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

                    do np=1,nspv
                        nb=(ispv(np,1)-1)*ndf+ispv(np,2)
                        bnd_va_pr_vec(np)=vspv(np)
                        bnd_no_pr_vec(np)=nb
                    enddo

                    do np=1,nssv
                        nb=(issv(np,1)-1)*ndf+issv(np,2)
                        bnd_va_se_vec(np)=vssv(np)
                        bnd_no_se_vec(np)=nb
                    enddo

end subroutine linder_miehe_boundary

subroutine linear_c3d8_3d_fem_geometry(neldirectional,length)
    ! ================================================================
    !the variables
    ! ================================================================
    ! ================================================================
    real(iwp),intent(in)::length(:)
    integer,intent(in)::neldirectional(:)


    integer,allocatable::nnddirectional(:)
    real(iwp),allocatable::deltal(:)

    integer::nblock,i,j,k,inode
    integer::ix,iy,iz,i_ele
    ! ================================================================
    npe=8

    allocate(deltal(dimen),nnddirectional(dimen))
    nnddirectional=neldirectional+1
    nnm=1

    nblock=1
    iel=1
    do i=1,dimen
        nblock=nblock*neldirectional(i)
        nnm=nnm*nnddirectional(i); enddo
        !
        nem=nblock

        allocate(nod(nem,npe),coords(nnm,dimen))
        nod=0
        coords=0


        nnddirectional=neldirectional+1
        deltal=length/neldirectional
        !c
        !c  cubical elements
        !c
        coords=0.0d0
        inode=0

        do  k = 1, nnddirectional(3)
            do  j = 1, nnddirectional(2)
                do  i = 1, nnddirectional(1)

                    inode=inode+1
                    coords(inode,:) =[deltal(1)*(i-1),deltal(2)*(j-1),deltal(3)*(k-1)]

                enddo
            enddo
        enddo

        do ix=0,neldirectional(1) -1
            do iy=0,neldirectional(2) -1
                do iz=0,neldirectional(3) -1
                    i_ele=1+ix+iy*neldirectional(1)+iz*neldirectional(1)*neldirectional(2)

                    nod(i_ele,1)=1 +ix+iy*nnddirectional(1)+iz*nnddirectional(1)*nnddirectional(2)
                    nod(i_ele,2)=nod(i_ele,1)+1
                    nod(i_ele,3)=nod(i_ele,2)+nnddirectional(1)
                    nod(i_ele,4)=nod(i_ele,3)-1
                    nod(i_ele,5:8)=nod(i_ele,1:4)+nnddirectional(1)*nnddirectional(2)

                enddo
            enddo
        enddo

    end subroutine linear_c3d8_3d_fem_geometry


    !  ================================================================
    !   the para view writer subroutine
    !  ================================================================
    subroutine paraview_3d_vtu_xml_writer(glu)
        implicit none
        !  ================================================================
        !   the variables
        !  ================================================================
        real(iwp),intent(in)::glu(:)
        real(iwp),allocatable::coords_deformed(:,:),coordst(:,:)

        !  ================================================================
        ! internal variables
        !  ================================================================
        integer,save::counter
        !  ===============================================================
        integer :: i,cell_types
        character(len=5)::x1
        character(len=40) :: fmt,filename ! format descriptor
        integer::inode,inem
        real(iwp)::scale
        !  ===============================================================
        fmt = '(i4.4)' ! an integer of width 5 with zeros at the left
        !  ===============================================================




        allocate(coords_deformed(nnm,dimen),coordst(nnm,dimen))


        if(npe.eq.20)cell_types= 25
        if(npe.eq.8)cell_types= 12


        counter=counter+1
        write(x1,fmt)counter
        coordst=0.0d0
        coords_deformed=0.0d0
        scale=1.0

        coordst=coords
        coords_deformed=coords
        !do j=1,dimen
        !forall(i=1:nnm) coords_deformed(i,j)=coordst(i,j)+glu((i-1)*ndf+j)
        !enddo


        filename='vtu/geom_test'//trim(x1)//'.vtu'
        write(*,*)filename
        open (vtu,file=filename)

        write(vtu,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
        write(vtu,*)'<UnstructuredGrid>'
        write(vtu,*)'<Piece NumberOfPoints="',nnm,'" NumberOfCells="',nem,'">'


        write(vtu,*)('<PointData Scalars="Elec_Pot" Vectors="displacement">')
        write(vtu,*)('<DataArray type="Float32" Name="Elec_Pot" NumberOfComponents="1" format="ascii">')

        write(vtu,*)(glu((i-1)*ndf+4),i=1,nnm);

        write(vtu,*)'</DataArray>'
        write(vtu,*)'<DataArray type="Float32" Name="displacement" NumberOfComponents="3" format="ascii">'
        write(vtu,105)(glu((i-1)*ndf+1:(i-1)*ndf+3),i=1,nnm);
        write(vtu,*)'</DataArray>'
        write(vtu,*)'</PointData>'


        write(vtu,*)'<Points>'
        write(vtu,*)'<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
        write(vtu,105)(coords_deformed(inode,:),inode=1,nnm)
        write(vtu,*)'</DataArray>'
        write(vtu,*)'</Points>'

        write(vtu,*)'<Cells>'
        write(vtu,*)'<DataArray type="Int32" Name="connectivity" Format="ascii">'
        write(vtu,fmt='(8(I5))')(NOD(INEM,:)-1,INEM=1,NEM)
        write(vtu,*)'</DataArray>'

        write(vtu,*)'<DataArray type="Int32" Name="offsets" Format="ascii">'
        write(vtu,*)(INEM*NPE,INEM=1,NEM)
        write(vtu,*)'</DataArray>'

        write(vtu,*)'<DataArray type="Int32" Name="types" Format="ascii">'
        write(vtu,*)(CELL_TYPES,INEM=1,NEM)

        write(vtu,*)'</DataArray>'
        write(vtu,*)'</Cells>'
        write(vtu,*)'</Piece>'
        write(vtu,*)'</UnstructuredGrid>'
        write(vtu,*)'</VTKFile>'

        close (vtu)

        !100   format(1x,a120)
105     format(3(e16.5))

    end subroutine paraview_3d_vtu_xml_writer


!  ================================================================
!   the para view writer subroutine
!  ================================================================
subroutine paraview_3d_vtu_xml_writer_vector(glu,elements_electric_field,elements_electric_polar)
implicit none
!  ================================================================
!   the variables
!  ================================================================
real(iwp),intent(in)::glu(:)
real(iwp),intent(in)::elements_electric_field(:,:),elements_electric_polar(:,:)
real(iwp),allocatable::coords_deformed(:,:),coordst(:,:)

!  ================================================================
! internal variables
!  ================================================================
integer,save::counter
!  ===============================================================
integer :: i,cell_types
character(len=5)::x1
character(len=52) :: fmt,filename ! format descriptor
integer::inode,inem
real(iwp)::scale

!  ===============================================================
fmt = '(i4.4)' ! an integer of width 5 with zeros at the left
!  ===============================================================

allocate(coords_deformed(nnm,dimen),coordst(nnm,dimen))

if(npe.eq.20)cell_types= 25
if(npe.eq.8)cell_types= 12

counter=counter+1
write(x1,fmt)counter
coordst=0.0d0
coords_deformed=0.0d0
scale=1.0

coordst=coords
coords_deformed=coords

!do j=1,dimen
!forall(i=1:nnm) coords_deformed(i,j)=coordst(i,j)+glu((i-1)*ndf+j)
!enddo


filename='./vtu/paraview_output_'//trim(x1)//'.vtu'
write(*,*)filename
open (vtu,file=filename)

   write(vtu,*)'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
   write(vtu,*)'<UnstructuredGrid>'
   write(vtu,*)'<Piece NumberOfPoints="',nnm,'" NumberOfCells="',nem,'">'
!  ===============================================================
 write(vtu,*)('<CellData Vectors="Electric_Field">')

 write(vtu,100)'<DataArray type="Float32" Name="Electric_Field" NumberOfComponents="3" Format="ascii">'
 do i=1,nem; write(vtu,105)elements_electric_field(i,:);enddo
 write(vtu,*)'</DataArray>'
! write(vtu,*)('</CellData>')
!!  =====================
! write(vtu,*)('<CellData Vectors="Electric_Polarization">')
 write(vtu,100)'<DataArray type="Float32" Name="Electric_Polarization" NumberOfComponents="3" Format="ascii">'
 do i=1,nem; write(vtu,105)elements_electric_polar(i,:);enddo
 write(vtu,*)'</DataArray>'
 write(vtu,*)('</CellData>')
!  ===============================================================


 write(vtu,*)('<PointData Scalars="Elec_Pot" Vectors="displacement">')
 write(vtu,100)('<DataArray type="Float32" Name="Elec_Pot" NumberOfComponents="1" Format="ascii">')

 write(vtu,*)(glu((i-1)*ndf+4),i=1,nnm);

 write(vtu,*)'</DataArray>'
 write(vtu,100)'<DataArray type="Float32" Name="displacement" NumberOfComponents="3" Format="ascii">'
 write(vtu,105)(glu((i-1)*ndf+1:(i-1)*ndf+3),i=1,nnm);
 write(vtu,*)'</DataArray>'
 write(vtu,*)'</PointData>'


   write(vtu,*)'<Points>'
   write(vtu,*)'<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
   write(vtu,105)(coords_deformed(inode,:),inode=1,nnm)
   write(vtu,*)'</DataArray>'
   write(vtu,*)'</Points>'

   write(vtu,*)'<Cells>'
   write(vtu,*)'<DataArray type="Int32" Name="connectivity" Format="ascii">'
   write(vtu,fmt='(8(I5))')(NOD(INEM,:)-1,INEM=1,NEM)
   write(vtu,*)'</DataArray>'

   write(vtu,*)'<DataArray type="Int32" Name="offsets" Format="ascii">'
   write(vtu,*)(INEM*NPE,INEM=1,NEM)
   write(vtu,*)'</DataArray>'

   write(vtu,*)'<DataArray type="Int32" Name="types" Format="ascii">'
   write(vtu,*)(CELL_TYPES,INEM=1,NEM)

   write(vtu,*)'</DataArray>'
   write(vtu,*)'</Cells>'
   write(vtu,*)'</Piece>'
   write(vtu,*)'</UnstructuredGrid>'
   write(vtu,*)'</VTKFile>'

close (vtu)

100   format(1x,a120)
105   format(3(e16.5))

end subroutine paraview_3d_vtu_xml_writer_vector


!  ================================================================
!   boundary conditions
!  ================================================================
subroutine crawley_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode
! integer::curved_node
integer::np,nb

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction


!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! body of the subroutine
!  ================================================================
coords=coords*1.0;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

nspv=count(coords(:,1)==-0.5*length(1))+ &
  count(coords(:,2)==-0.5*length(2))+ &
   2*count(coords(:,3)==-0.5*length(3))+ &
  count(coords(:,3)== 0.5*length(3))

ix_periodic=count(coords(:,1)==0.5*length(1))
iy_periodic=count(coords(:,2)==0.5*length(2))
iz_periodic=count(coords(:,3)==0.5*length(3))



allocate(afc_boundary_x_periodic(ix_periodic) , &
   afc_boundary_y_periodic(iy_periodic) , &
   afc_boundary_z_periodic(iz_periodic) )
!
nssv=1
!
!write(*,*)nspv
allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0

ix_periodic=0
iy_periodic=0
iz_periodic=0

ibond=0
do inode=1,nnm

if(coords(inode,1)==0.5*length(1))then
ix_periodic=ix_periodic+1
afc_boundary_x_periodic(ix_periodic)=inode
endif

if(coords(inode,2)==0.5*length(2))then
iy_periodic=iy_periodic+1
afc_boundary_y_periodic(iy_periodic)=inode
endif

if(coords(inode,3)==0.5*length(3))then
iz_periodic=iz_periodic+1
afc_boundary_z_periodic(iz_periodic)=inode
endif




if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif
!
!
!
if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;

 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0d0;
! write(1,*)inode

endif
!
if(coords(inode,3)==0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=1.0;
endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif
enddo

! ispv(nspv,1)=curved_node;
! ispv(nspv,2)=4;
! vspv(nspv)=1.0;


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo
! write(1,*)bnd_va_pr_vec
! write(1,*)bnd_no_pr_vec

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo





end subroutine crawley_boundary


!! This subroutine define and producd constraint matrix for periodic boundary conditions
!! @param constraint
!! >The essential boundary conditions are applied through constrainted matrix
!! First the number of fixed boundary conditions are added to constrained matrix
!! Number of essential boundary conditions are substracted from number of column
!! Then rows related to the zero or fixed boundary conditions are put to zero

subroutine afc_form_constraint()

real(iwp),allocatable::constraint_x(:,:),constraint_y(:,:),constraint_z(:,:) !<the constraint matrix
real(iwp),allocatable::full_constraint(:,:)
logical,allocatable::fixed_boundary_mask(:)
integer::n_constraint !<number of constraintd
integer::i_const,j_const,const_num,n_fixed_boundary
integer::i_boundary



n_fixed_boundary=count(vspv(:)==0)
fixed_boundary_mask=vspv(:)==0
write(*,*)'n_fixed_boundary',n_fixed_boundary
write(1,*)'vspv',vspv
write(1,*)'fixed_boundary_mask',fixed_boundary_mask

n_constraint=size(afc_boundary_x_periodic,dim=1)+  &
             size(afc_boundary_y_periodic,dim=1)+  &
             size(afc_boundary_z_periodic,dim=1) + &
             n_fixed_boundary

write(*,*)'n_constraint=',n_constraint
!! we will add three sets of periodic boundary condistions
!! so we have to add one to the dimension of constraint for each set
allocate(constraint(neq,neq-n_constraint+3),is_constrained(neq))
allocate(full_constraint(neq,neq))


full_constraint=0
do j_const=1,neq; full_constraint(j_const,j_const)=1; enddo


!constraint_x=constraint;
!constraint_y=constraint;
!constraint_z=constraint;

constraint=0
is_constrained=.false.

!write(*,*)afc_boundary_x_periodic

do j_const=1,size(afc_boundary_x_periodic,dim=1)

    const_num=(afc_boundary_x_periodic(j_const)-1)*ndf+1
    constraint(const_num,(afc_boundary_x_periodic(1)-1)*ndf+1)=1

    full_constraint(const_num,:)=0;
    full_constraint(const_num,(afc_boundary_x_periodic(1)-1)*ndf+1)=1.0d0;

    is_constrained(const_num)=.true.
enddo;




do j_const=1,size(afc_boundary_y_periodic,dim=1)

    const_num=(afc_boundary_y_periodic(j_const)-1)*ndf+2
    full_constraint(const_num,:)=0;
    full_constraint(const_num,(afc_boundary_y_periodic(1)-1)*ndf+2)=1.0d0;
    is_constrained(const_num)=.true.
enddo;



do j_const=1,size(afc_boundary_z_periodic,dim=1)
    const_num=(afc_boundary_z_periodic(j_const)-1)*ndf+3
    full_constraint(const_num,:)=0;
    full_constraint(const_num,(afc_boundary_z_periodic(1)-1)*ndf+3)=1;

    is_constrained(const_num)=.true.
enddo;



do i_boundary=1,size(vspv(:))
    if (vspv(i_boundary)==0) then
    i_const=(ispv(i_boundary,1)-1)*ndf+ispv(i_boundary,2)
    full_constraint(i_const,:)=0;
    endif
enddo

constraint=0
i_const=0

do j_const=1,neq
if (   sum(  full_constraint(:,j_const) ) .gt. 0 ) then
 i_const=i_const+1
 constraint(:,i_const)=full_constraint(:,j_const);
endif
enddo

end subroutine afc_form_constraint


!  ================================================================
!   boundary conditions
!  ================================================================
subroutine wedge_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode
! integer::curved_node
integer::np,nb

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction


!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! body of the subroutine
!  ================================================================
coords=coords*1.0;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

nspv=count(coords(:,1)==-0.5*length(1))+ &
  count(coords(:,2)==-0.5*length(2))+ &
  count(coords(:,3)==-0.5*length(3))+ &
  count(coords(:,1)== 0.5*length(1))+ &
  count(coords(:,2)== 0.5*length(2))- &
 2* count(   coords(:,1)== 0.5*length(1)  .and.  coords(:,2)== 0.5*length(2))

write(*,*)nspv
!ix_periodic=count(coords(:,1)==0.5*length(1))
!iy_periodic=count(coords(:,2)==0.5*length(2))
!iz_periodic=count(coords(:,3)==0.5*length(3))



!allocate(afc_boundary_x_periodic(ix_periodic) , &
!   afc_boundary_y_periodic(iy_periodic) , &
!   afc_boundary_z_periodic(iz_periodic) )
!
nssv=1
!
!write(*,*)nspv
allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0

ix_periodic=0
iy_periodic=0
iz_periodic=0

ibond=0
do inode=1,nnm

if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif
!
!
!
if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;

endif
!
if(coords(inode,1)==0.5*length(1))then
if(.not.coords(inode,2)==0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=1.0;
endif
endif

if(coords(inode,2)==0.5*length(2))then
if(.not.coords(inode,1)==0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=-1.0;
endif
endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif
enddo

allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))



do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo

write(*,*)bnd_no_pr_vec

end subroutine wedge_boundary



!< This subroutine define and producd the boundary conditions for afc unit cell !>
!! This subroutine define and producd the boundary conditions for afc unit cell
!! @param lentgh of the domain in each direction
!! @param curved_node End point of the domain. This is the node that plots will be bases on
!! @todo apply the periodic boundary conditions
subroutine afc_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode,jbond
! integer::curved_node
integer::np,nb
!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! afc specific boundary
!  ================================================================
integer,allocatable::afc_positive_electr(:)
integer,allocatable::afc_negative_electr(:)

integer::number_afc_positive_electr
integer::number_afc_negative_electr

integer::number_afc_positive_z0_electr
integer::number_afc_negative_zL_electr

integer,allocatable::afc_negative_electr_z0(:)
integer,allocatable::afc_positive_electr_zL(:)

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction

integer::phi_0_periodic !<number of periodic boundary condition nodes on electrode z0 direction
integer::phi_L_periodic !<number of periodic boundary condition nodes on electrode zL  direction

!  ================================================================
! body of the subroutine
!  ================================================================
!coords=coords*1.0e6;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo


do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

!  ================================================================
! afc specific boundary
!  ================================================================
number_afc_positive_electr=32
number_afc_negative_electr=32


allocate(afc_positive_electr(number_afc_positive_electr),afc_negative_electr(number_afc_negative_electr))

afc_positive_electr=[   1,   2,   3,   4,  11,  12,  17,  18,  29,  30,  31,  32,  33,  34,  35,  36, &
                        37,  38,  39,  40,  71,  72,  73,  74,  75,  76, 103, 104, 105, 106, 107, 108]

afc_negative_electr=[   5,   6,   7,   8,   9,  10,  19,  20,  41,  42,  43,  44,  45,  46,  47,  48 , &
                        49,  50,  51,  52,  65,  66,  67,  68,  69,  70, 109, 110, 111, 112, 113, 114]



number_afc_positive_z0_electr=16
number_afc_negative_zL_electr=16

allocate(afc_negative_electr_z0(number_afc_negative_zL_electr),afc_positive_electr_zL(number_afc_negative_zL_electr))

afc_negative_electr_z0=[ &
   9,  10,  19,  20,  65,  66,  67,  68,  69,  70, 109, 110, 111, 112, 113, 114]
afc_positive_electr_zL=[ &
  11,  12,  17,  18,  71,  72,  73,  74,  75,  76, 103, 104, 105, 106, 107, 108]


!write(*,*)length

!nspv=count(coords(:,1)==-0.5*length(1))+ &
!     count(coords(:,2)==-0.5*length(2))+ &
!     count(coords(:,3)==-0.5*length(3))+ &
!     number_afc_positive_electr+number_afc_negative_electr

nspv=count(coords(:,1)==-0.5*length(1))+ &
     count(coords(:,2)==-0.5*length(2))+ &
     count(coords(:,3)==-0.5*length(3))+ &
     number_afc_positive_z0_electr+number_afc_negative_zL_electr



ix_periodic=count(coords(:,1)==0.5*length(1))
iy_periodic=count(coords(:,2)==0.5*length(2))
iz_periodic=count(coords(:,3)==0.5*length(3))

phi_0_periodic=number_afc_positive_z0_electr
phi_L_periodic=number_afc_negative_zL_electr



allocate(afc_boundary_x_periodic(ix_periodic) , &
         afc_boundary_y_periodic(iy_periodic) , &
         afc_boundary_z_periodic(iz_periodic) , &
         afc_negative_electr_z0_periodic(phi_0_periodic) , &
         afc_negative_electr_zL_periodic(phi_L_periodic) )

nssv=1
!
!write(*,*)nspv
!write(*,*)nssv

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0
ibond=0


ix_periodic=0
iy_periodic=0
iz_periodic=0

do inode=1,nnm

if(coords(inode,1)==0.5*length(1))then
ix_periodic=ix_periodic+1
afc_boundary_x_periodic(ix_periodic)=inode
endif

if(coords(inode,2)==0.5*length(2))then
iy_periodic=iy_periodic+1
afc_boundary_y_periodic(iy_periodic)=inode
endif

if(coords(inode,3)==0.5*length(3))then
iz_periodic=iz_periodic+1
afc_boundary_z_periodic(iz_periodic)=inode
endif



if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif


if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;
endif


if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif

enddo

!  ================================================================
! afc specific boundary
!  ================================================================
do jbond=1,number_afc_positive_z0_electr
 ibond=ibond+1;
 ispv(ibond,1)=afc_negative_electr_z0(jbond);
 ispv(ibond,2)=4;
 vspv(ibond)=1.0;
enddo

do jbond=1,number_afc_negative_zL_electr
 ibond=ibond+1;
 ispv(ibond,1)=afc_positive_electr_zL(jbond);
 ispv(ibond,2)=4;
 vspv(ibond)=-1.0;
enddo


!write(*,*)ix_periodic,iy_periodic,iz_periodic
!write(*,*)afc_positive_electr


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
    nb=(ispv(np,1)-1)*ndf+ispv(np,2)
    bnd_va_pr_vec(np)=vspv(np)
    bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo

end subroutine afc_boundary



!  ================================================================
!   boundary conditions
!  ================================================================
subroutine pull_z_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode
! integer::curved_node
integer::np,nb

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction


!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! body of the subroutine
!  ================================================================
coords=coords*1.0;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

nspv=count(coords(:,1)==-0.5*length(1))+ &
  count(coords(:,2)==-0.5*length(2))+ &
   2*count(coords(:,3)==-0.5*length(3))+ &
  count(coords(:,3)== 0.5*length(3))

ix_periodic=count(coords(:,1)==0.5*length(1))
iy_periodic=count(coords(:,2)==0.5*length(2))
iz_periodic=count(coords(:,3)==0.5*length(3))



allocate(afc_boundary_x_periodic(ix_periodic) , &
   afc_boundary_y_periodic(iy_periodic) , &
   afc_boundary_z_periodic(iz_periodic) )
!
nssv=1
!
!write(*,*)nspv
allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0

ix_periodic=0
iy_periodic=0
iz_periodic=0

ibond=0
do inode=1,nnm

if(coords(inode,1)==0.5*length(1))then
ix_periodic=ix_periodic+1
afc_boundary_x_periodic(ix_periodic)=inode
endif

if(coords(inode,2)==0.5*length(2))then
iy_periodic=iy_periodic+1
afc_boundary_y_periodic(iy_periodic)=inode
endif

if(coords(inode,3)==0.5*length(3))then
iz_periodic=iz_periodic+1
afc_boundary_z_periodic(iz_periodic)=inode
endif




if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif
!
!
!
if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;

 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0d0;

endif
!
if(coords(inode,3)==0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=1.0;
endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif
enddo

! ispv(nspv,1)=curved_node;
! ispv(nspv,2)=4;
! vspv(nspv)=1.0;


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo

end subroutine pull_z_boundary






! ================================================================
!! This subroutine reads mesh file containing geometry of unit cell.
!! @param constraint

subroutine afc_mesh_file_reader()

!the variables
! ================================================================
! ================================================================


real(iwp)::corner_point(dimen),length_1(dimen)

integer::i,j,inode
integer::ix,i_ele
 ! ================================================================
 npe=8
 iel=1


open (msh,file='meshes/3d_afc_mesh.msh')
read(msh,*)nnm,nem; !number of nodes, number of elements

allocate(nod(nem,npe),coords(nnm,dimen))
nod=0
coords=0

do inode=1,nnm
read(msh,*)ix,(coords(inode,i),i=1,dimen);enddo;


!  reading the connectivity array from the mesh file
 read(msh,*)
 do i=1,nem
   read(msh,*)i_ele,(nod(i,j),j=1,npe);enddo
close(msh)


do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length_1(i)=maxval(coords(:,i))-minval(coords(:,i)); enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length_1(i);  enddo


end subroutine afc_mesh_file_reader


!  ================================================================
!   boundary conditions
!  ================================================================
subroutine dealii_test_boundary(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode
! integer::curved_node
integer::np,nb

integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction


!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! body of the subroutine
!  ================================================================
coords=coords*1.0;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo
do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

!nspv=count(coords(:,1)==-0.5*length(1))+ &
!  count(coords(:,2)==-0.5*length(2))+ &
!   2*count(coords(:,3)==-0.5*length(3))+ &
!  nnm

nspv=count(coords(:,1)==-0.5*length(1))+ &
  count(coords(:,1)== 0.5*length(2))+ &
  count(coords(:,2)==-0.5*length(2))+ &
  count(coords(:,3)==-0.5*length(3))+ &
  nnm
!  count(coords(:,3)== 0.5*length(3))

!! This will put all the potential in all the domain to be zero.

ix_periodic=count(coords(:,1)==0.5*length(1))
iy_periodic=count(coords(:,2)==0.5*length(2))
iz_periodic=count(coords(:,3)==0.5*length(3))



allocate(afc_boundary_x_periodic(ix_periodic) , &
   afc_boundary_y_periodic(iy_periodic) , &
   afc_boundary_z_periodic(iz_periodic) )
!
nssv=1
!
!write(*,*)nspv
allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0

ix_periodic=0
iy_periodic=0
iz_periodic=0

ibond=0
do inode=1,nnm

if(coords(inode,1)==0.5*length(1))then
ix_periodic=ix_periodic+1
afc_boundary_x_periodic(ix_periodic)=inode
endif

if(coords(inode,2)==0.5*length(2))then
iy_periodic=iy_periodic+1
afc_boundary_y_periodic(iy_periodic)=inode
endif

if(coords(inode,3)==0.5*length(3))then
iz_periodic=iz_periodic+1
afc_boundary_z_periodic(iz_periodic)=inode
endif




if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif
!
!
!
if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;
endif

if(coords(inode,1)==0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=1.0d0;
endif

!
!if(coords(inode,3)==0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=4;
 vspv(ibond)=0.0;
!endif

if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif
enddo

! ispv(nspv,1)=curved_node;
! ispv(nspv,2)=4;
! vspv(nspv)=1.0;


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo

end subroutine dealii_test_boundary



!< This subroutine define and producd the boundary conditions for afc unit cell !>
!! This subroutine define and producd the boundary conditions for afc unit cell
!! @param lentgh of the domain in each direction
!! @param curved_node End point of the domain. This is the node that plots will be bases on
!! @todo apply the periodic boundary conditions
subroutine afc_boundary_full_electrode(length)
!  ================================================================
! input geormetry variables
!  ================================================================
implicit none
!  ================================================================
! input geormetry variables
!  ================================================================
real*8:: length(:)
integer::ibond,inode,jbond
! integer::curved_node
integer::np,nb
!  ================================================================
real*8::corner_point(dimen)
integer::i ! ,j,k
!  ================================================================
! afc specific boundary
!  ================================================================
integer,allocatable::afc_positive_electr(:)
integer,allocatable::afc_negative_electr(:)
integer::number_afc_positive_electr
integer::number_afc_negative_electr



integer::ix_periodic !<number of periodic boundary condition nodes on x direction
integer::iy_periodic !<number of periodic boundary condition nodes on y direction
integer::iz_periodic !<number of periodic boundary condition nodes on z direction

!  ================================================================
! body of the subroutine
!  ================================================================
!coords=coords*1.0e6;

do i=1,dimen; corner_point(i)=minval(coords(:,i));enddo
do i=1,dimen; coords(:,i)=coords(:,i)-corner_point(i);enddo
do i=1,dimen; length(i)=maxval(coords(:,i))-minval(coords(:,i));
enddo


do i=1,dimen; coords(:,i)= coords(:,i)-0.5*length(i);  enddo

!  ================================================================
! afc specific boundary
!  ================================================================
number_afc_positive_electr=32
number_afc_negative_electr=32


allocate(afc_positive_electr(number_afc_positive_electr),afc_negative_electr(number_afc_negative_electr))

afc_positive_electr=[   1,   2,   3,   4,  11,  12,  17,  18,  29,  30,  31,  32,  33,  34,  35,  36, &
                        37,  38,  39,  40,  71,  72,  73,  74,  75,  76, 103, 104, 105, 106, 107, 108]

afc_negative_electr=[   5,   6,   7,   8,   9,  10,  19,  20,  41,  42,  43,  44,  45,  46,  47,  48 , &
                        49,  50,  51,  52,  65,  66,  67,  68,  69,  70, 109, 110, 111, 112, 113, 114]

!write(*,*)length

nspv=count(coords(:,1)==-0.5*length(1))+ &
     count(coords(:,2)==-0.5*length(2))+ &
   1*count(coords(:,3)==-0.5*length(3))+ &
   0*count(coords(:,3)== 0.5*length(3))+ &
     number_afc_positive_electr+number_afc_negative_electr



ix_periodic=count(coords(:,1)==0.5*length(1))
iy_periodic=count(coords(:,2)==0.5*length(2))
iz_periodic=count(coords(:,3)==0.5*length(3))



allocate(afc_boundary_x_periodic(ix_periodic) , &
         afc_boundary_y_periodic(iy_periodic) , &
         afc_boundary_z_periodic(iz_periodic) )

nssv=1
!
!write(*,*)nspv
!write(*,*)nssv

allocate(ispv(nspv,2),vspv(nspv),vspvt(nspv))
allocate(issv(nssv,2),vssv(nssv),vssvt(nssv))
vssv=0.0d0
issv(1,1)=1
issv(1,2)=1
!
!!manually forming the bounrary condition
!
vspv=0.0d0
ibond=0


ix_periodic=0
iy_periodic=0
iz_periodic=0

do inode=1,nnm

if(coords(inode,1)==0.5*length(1))then
ix_periodic=ix_periodic+1
afc_boundary_x_periodic(ix_periodic)=inode
endif

if(coords(inode,2)==0.5*length(2))then
iy_periodic=iy_periodic+1
afc_boundary_y_periodic(iy_periodic)=inode
endif

if(coords(inode,3)==0.5*length(3))then
iz_periodic=iz_periodic+1
afc_boundary_z_periodic(iz_periodic)=inode
endif



if(coords(inode,1)==-0.5*length(1))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=1;
 vspv(ibond)=0.0d0;
 endif


if(coords(inode,2)==-0.5*length(2))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=2;
 vspv(ibond)=0.0d0;
 endif
!
!
if(coords(inode,3)==-0.5*length(3))then
 ibond=ibond+1;
 ispv(ibond,1)=inode;
 ispv(ibond,2)=3;
 vspv(ibond)=0.0d0;
endif


if(coords(inode,1)==0.5*length(1).and.coords(inode,2)==0.5*length(2).and.coords(inode,3)==0.5*length(3))then
curved_node=inode
endif

enddo

!  ================================================================
! afc specific boundary
!  ================================================================
do jbond=1,number_afc_positive_electr
 ibond=ibond+1;
 ispv(ibond,1)=afc_positive_electr(jbond);
 ispv(ibond,2)=4;
 vspv(ibond)=-1.0;
enddo

do jbond=1,number_afc_negative_electr
 ibond=ibond+1;
 ispv(ibond,1)=afc_negative_electr(jbond);
 ispv(ibond,2)=4;
 vspv(ibond)=1.0;
enddo


!write(*,*)ix_periodic,iy_periodic,iz_periodic
!write(*,*)afc_positive_electr


allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

do np=1,nspv
nb=(ispv(np,1)-1)*ndf+ispv(np,2)
bnd_va_pr_vec(np)=vspv(np)
bnd_no_pr_vec(np)=nb
enddo

do np=1,nssv
  nb=(issv(np,1)-1)*ndf+issv(np,2)
  bnd_va_se_vec(np)=vssv(np)
  bnd_no_se_vec(np)=nb
enddo

end subroutine afc_boundary_full_electrode

end module fem_geometry
