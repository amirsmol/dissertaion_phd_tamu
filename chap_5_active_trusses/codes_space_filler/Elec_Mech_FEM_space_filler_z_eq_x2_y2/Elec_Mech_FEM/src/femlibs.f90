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

module fem_libs

use fem_geometry
use material_behavior
use linearsolvers

implicit none

real(iwp),allocatable:: elcrds(:,:) !the geometries of nodes of the element

contains
!     ================================================================
!     forming the global martixes and vetors
!     ================================================================
subroutine glbmatrcs(glu,glk,glq,glp,glb)
implicit none
!     ================================================================
!  input variables
!     ================================================================
real(iwp),intent(in)::glu(:),glp(:),glb(:) ! the dof vector, current time
!     ================================================================
!  output variables
!     ================================================================
real(iwp),intent(out)::glq(:) ! residual vector
real(iwp),intent(out)::glk(:,:) ! the global stiffness matrix
!     ================================================================
!  trivial variables
!     ================================================================
real(iwp),allocatable:: elq(:),elk(:,:) !the force and coeffi
real(iwp),allocatable:: elu(:),elp(:),elb(:) ! degrees of freedom vector in current time
integer::elneq
!     ================================================================
integer::inpe,iglobal,idof,jdof
!     =========================element global node and dof map
integer,allocatable::el_nod_vec(:) !el_nod_vec(npe)
integer,allocatable::el_dof_vec(:)  !el_dof_vec(npe*ndf)
!     ================================================================

elneq=npe*ndf;
allocate(elk(elneq,elneq),elq(elneq),elu(elneq),elp(elneq),elb(elneq))
!     ================================================================
elu=0.0d0;elq=0.0d0;elk=0.0d0;

gauss_point_number=0

!>     starting the element (i,:)

allocate(el_nod_vec(npe),el_dof_vec(elneq),elcrds(npe,dimen))

do noelem=1,nem
el_nod_vec=nod(noelem,:)
elcrds=coords(el_nod_vec,:)

do inpe=1,npe
iglobal=nod(noelem,inpe)
idof=inpe*ndf
do jdof=0,ndf-1
el_dof_vec(idof-ndf+1+jdof)=(iglobal*ndf-ndf+1+jdof)

enddo !j=1,ndf
enddo !inpe=1,npe

elu=glu(el_dof_vec)
elp=glp(el_dof_vec)
elb=glb(el_dof_vec)

!     ===============================================================
!     forming the coefficitne matrix and source vector for each element
!     ===============================================================
    if(npe.eq.20.or.npe.eq.8) call elmatrcs_3d(elu,elp,elb,elk,elq)
    if(npe.eq.2)              call truss_elmatrcs_3d(elu,elp,elb,elk,elq)
!     ===============================================================
!     assembling the global source vector and coefficient matrix
!     ===============================================================

glq(el_dof_vec)=glq(el_dof_vec)+elq
glk(el_dof_vec,el_dof_vec)=glk(el_dof_vec,el_dof_vec)+elk


enddo ! noelem=1,nem

deallocate(el_nod_vec,el_dof_vec,elcrds)

end subroutine glbmatrcs


!     ===============================================================
!                       residula vector subroutine
!     ===============================================================
      subroutine truss_residual3d(inode,gdsf,ur,der,sigma,res_vect)
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      integer,intent(in):: inode !node number
      real(iwp),intent(in)::gdsf(dimen,npe) !  the shape functions
      real(iwp),intent(in),dimension(ndf,dimen):: ur
      real(iwp),intent(in),dimension(dimen)::der ! electric displacement tensor
      real(iwp),intent(in),dimension(dimen,dimen)::sigma !stress tensor
!     ================================================================
!                          output variables
!     ================================================================
      real(iwp),intent(out)::res_vect(ndf) !the coefficient matri
!     ================================================================
      integer::i,k,j,m
      integer::del(dimen,dimen),eye(ndf,ndf)
!     ================================================================
!     deformmation and strain tensors
!     ================================================================

!     ===========================shape function tensors
      real(iwp)::bepsiloni(dimen,dimen,ndf),beleci(dimen,ndf)
!     ================================================================
      del = 0; do i = 1,dimen; del(i,i) = 1;enddo
      eye = 0; do i = 1,ndf; eye(i,i) = 1;enddo
!!     ================================================================
!!            c     shape function tensor
!!     ================================================================
      beleci=0.0d0;

      bepsiloni=0.0d0;


      beleci(:,dimen+1) = -gdsf(:, inode)

      do 10 k=1,ndf
      do 10 i=1,dimen
      do 10 j=1,dimen

      bepsiloni(i,j,k)= 0.5*gdsf(j,inode)*eye(k,i)+0.5*gdsf(i,inode)*eye(k,j)

      do 10 m=1,dimen

      bepsiloni(i,j,k)=bepsiloni(i,j,k)                      &
                   + 0.5*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
                   + 0.5*gdsf(j,inode)*eye(k,m)*ur(m,i)

 10   continue

!write(1,*)'gdsf=',gdsf

!     ================================================================
      res_vect=0.0d0


      do 40 k=1,ndf
      do 40 i=1,dimen
      res_vect(k)=res_vect(k)-der(i)*beleci(i,k)

      do 40 j=1,dimen

40    res_vect(k)=res_vect(k)+sigma(i,j)*bepsiloni(i,j,k)

      end subroutine truss_residual3d
!     ===============================================================
!                       the element coefficient matrix
!     ===============================================================
subroutine truss_elmatrcs_3d(elu,elp,elb,elk,elq) !(elcrds,elu,elk,elq)
!     ================================================================
!     element  calculations based on  linear and quadratic rectangular
!     relements with isoparametric  formulation.
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real(iwp),intent(in):: elu(:),elp(:),elb(:) !dof vectors
!      real(iwp),intent(in):: elcrds(npe,dimen)
!     ================================================================
!                          otput variables
!     ================================================================
      real(iwp),intent(out):: elq(:) !element source vector
      real(iwp),intent(out):: elk(:,:) !the element coefficient matrix
!     ================================================================
      integer::ni,ii,jj,li,lj,idof
      integer::inode,jnode
!     ================================================================
      real(iwp):: elu_tense(npe,ndf) !dof current time
      real(iwp):: elp_tense(npe,ndf) !dof current time
      real(iwp):: elb_tense(npe,ndf) !dof current time
      real(iwp)::coord(dimen) !the global coordinate of each node
      real(iwp)::xi(dimen)! transformed coordinates
      real(iwp)::sf(npe),gdsf(dimen,npe)
      real(iwp)::ur(ndf,dimen),up(ndf,dimen),ub(ndf,dimen) !value of the function and its derivative
      real(iwp)::k_coef(ndf,ndf) ! the auxilary matrix for cuefficitne matr
      real(iwp)::res_vect(ndf)
      real(iwp)::cnst,det
      real(iwp)::der(dimen) !electric displacement
      real(iwp)::sigma(dimen,dimen) !stress
!     ================================================================
      real(iwp):: gauspt(5,5),gauswt(5,5) ! gauss points and gauss weights
!     =========================element center point and transformation variables
      real(iwp):: el_center_coord(dimen)
      real(iwp):: qtran(dimen,dimen)
!     ================================================================
      real(iwp),dimension(dimen,dimen) ::sigma_shape
      real(iwp),dimension(dimen,dimen) ::R_rotation_tensor


 !     ================================================================
       data gauspt/5*0.0d0, -0.57735027d0, 0.57735027d0, 3*0.0d0,            &
       -0.77459667d0, 0.0d0, 0.77459667d0, 2*0.0d0, -0.86113631d0,           &
       -0.33998104d0, 0.33998104d0, 0.86113631d0, 0.0d0, -0.90617984d0,      &
       -0.53846931d0,0.0d0,0.53846931d0,0.90617984d0/
!
      data gauswt/2.0d0, 4*0.0d0, 2*1.0d0, 3*0.0d0, 0.55555555d0,        &
        0.88888888d0, 0.55555555d0, 2*0.0d0, 0.34785485d0,               &
      2*0.65214515d0, 0.34785485d0, 0.0d0, 0.23692688d0,                 &
        0.47862867d0, 0.56888888d0, 0.47862867d0, 0.23692688d0/
!     ===============================================================
      elq=0.0d0
      elk=0.0d0

!     ================================================================
    do ni=1,npe
        elu_tense(ni,:)=elu(1+(ni-1)*ndf:ni*ndf)
        elp_tense(ni,:)=elp(1+(ni-1)*ndf:ni*ndf)
        elb_tense(ni,:)=elb(1+(ni-1)*ndf:ni*ndf)
    enddo !i=1,npe
!     ================================================================
!     finding stress due to shape change
!     ================================================================
          element_direction_normal=unit_vect(elcrds(2,:)-elcrds(1,:))
!     ===============================================================
          do 200 ni = 1,ipdf
          xi(1) = gauspt(ni,ipdf)
          gauss_point_number=gauss_point_number+1

          call truss_3d_shape_3d_linea(elcrds,xi,det,sf,gdsf)
          cnst = det*gauswt(ni,ipdf)
          coord=0
!     ===============================================================
!                       coordinate in refrence and element frame
!     ===============================================================
          coord=matmul(sf,elcrds)
!     ===============================================================
!                       post process
!     ===============================================================
          ur       =transpose(       matmul( gdsf,elu_tense     ) )
          up       =transpose(       matmul( gdsf,elp_tense     ) )
          ub       =transpose(       matmul( gdsf,elb_tense     ) )
          call truss_material_properties()
          call truss_shape_change_stress_z_eq_xy(coord,sigma_shape)
          ! call shape_change_stress_beam_folding(coord,sigma_shape)
          call truss_stress_elect_displacement(ur,up,ub,der,sigma)

     ! do ii=1,dimen;do jj=1,dimen;
     !    each_truss_strain(noelem)=each_truss_strain(noelem)+ &
     !    sigma_shape(ii,jj)*element_direction_normal(ii)*element_direction_normal(jj)
     ! enddo;enddo


!    call show_matrix(sigma_shape,'sigma_shape')
der=0.0d0;
!     ===============================================================
!                       filling the elemt matrixes
!     ===============================================================

do inode=1,npe
    ii = ndf*(inode-1)+1
    li=(ii+ndf-1)
    call truss_residual3d(inode,gdsf,ur,der,sigma-sigma_shape,res_vect)

    elq(ii:li)=elq(ii:li)+res_vect*cnst
    do  jnode=1,npe
  call  truss_k_gen3d(jnode,inode,gdsf,ur,k_coef)
  jj = ndf*(jnode-1)+1
  lj=(jj+ndf-1)
    elk(ii:li,jj:lj)=elk(ii:li,jj:lj)+k_coef*cnst
    enddo;
enddo



200 continue

!call show_vector(elq,'elq')
!call show_matrix(elk,'elk')

end subroutine truss_elmatrcs_3d

!     ===============================================================
!                       THE SHAPE FUNCTION SUBROUTINE
!     ===============================================================
subroutine truss_3d_shape_3d_linea(elcrds,xi,det,sf,gdsf_3d)
      real(iwp)::elcrds(:,:)
      real(iwp)::det ! the jacobian determinant
      real(iwp),intent(out) ::sf(:),gdsf_3d(:,:) ! shape function sf(npe),gdsf(dimen,npe)

      real(iwp)::dsf(dimen,npe) !(dimen,npe) ! the differentiatian of shape function{


      real(iwp)::xnode(2) ! coordinate of nodes in refrence
      real(iwp)::xi(dimen) !locad coordinate in element
      integer::i
      real(iwp)::delXi_delzeta(dimen) ! local coordinate derivitives
      real(iwp)::delzeta_delXi(dimen) ! local coordinate derivitives
!     ================================================================
      data xnode/ -1,1 /!master element coordinates
!     ================================================================

      do 10 i=1,npe
      sf(i)=0.5d0*( 1+xi(1)*xnode(i) )
      dsf(1,i)=0.5d0*xnode(i)
10    continue

!     ===================================================================
!         compute the jacobian matrix [gj] and its inverse [gjinv]
!     ===================================================================
      det=norm_vect(elcrds(1,:)-elcrds(2,:))/abs(xnode(2)-xnode(1))

      delXi_delzeta=0.0d0

      do i=1,npe
      delXi_delzeta(:)=delXi_delzeta(:)+dsf(1,i)*elcrds(i,:)
      enddo

      delzeta_delXi=1/delXi_delzeta

      do i=1,dimen
      if(   (  abs(elcrds(1,i)-elcrds(2,i) ) .lt. default_smallest_pivot ) ) then
      delzeta_delXi(i)=0.0d0;
      endif
      enddo

      do i=1,npe
      gdsf_3d(:,i)=dsf(1,i)*delzeta_delXi(:)
      enddo

end subroutine truss_3d_shape_3d_linea
!     ===============================================================
!                       coefficient matrix subroutine
!     ===============================================================
      subroutine  truss_k_gen3d(inode,jnode,gdsf,ur,k_coef)
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      integer,intent(in)::inode,jnode !node number iterators
      real(iwp),intent(in)::gdsf(:,:) !  the shape functions
      real(iwp),intent(in)::ur(:,:)! displacement
!     ================================================================
!                          output variables
!     ================================================================
      real(iwp),intent(out)::k_coef(:,:) !the coefficient matrix
!     ================================================================
      integer::i,j,k,l,m,n
      integer,save::counter
      integer::del(dimen,dimen),eye(ndf,ndf)
!     ===========================shape function tensors
      real(iwp)::bepsiloni(dimen,dimen,ndf),beleci(dimen,ndf)
      real(iwp)::bepsilonj(dimen,dimen,ndf),belecj(dimen,ndf)
!!     ================================================================
      counter=counter+1
      del = 0; do i = 1,dimen; del(i,i) = 1;enddo
      eye = 0; do i = 1,ndf;   eye(i,i) = 1;enddo
!!     ================================================================
!!            c     material properties
!!     ================================================================
      beleci=0.0d0;
      belecj=0.0d0;
      bepsiloni=0.0d0;
      bepsilonj=0.0d0;

      beleci(:,dimen+1) = -gdsf(:, inode)
      belecj(:,dimen+1) = -gdsf(:, jnode)

      do 10 k=1,ndf
      do 10 i=1,dimen
      do 10 j=1,dimen

      bepsilonj(i,j,k)= 0.5*gdsf(j,jnode)*eye(k,i)+0.5*gdsf(i,jnode)*eye(k,j)

      bepsiloni(i,j,k)= 0.5*gdsf(j,inode)*eye(k,i)+0.5*gdsf(i,inode)*eye(k,j)

      do 10 m=1,dimen

      bepsilonj(i,j,k)=bepsilonj(i,j,k)                      &
                   + 0.5*gdsf(i,jnode)*eye(k,m)*ur(m,j)   &
                   + 0.5*gdsf(j,jnode)*eye(k,m)*ur(m,i)

      bepsiloni(i,j,k)=bepsiloni(i,j,k)                      &
                   + 0.5*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
                   + 0.5*gdsf(j,inode)*eye(k,m)*ur(m,i)

 10   continue

      call truss_material_properties()

!     =====================================time dependent materials properties

ctens=(k00+0.5d0*k01)*ctens
partial_sigma_to_partial_elec_t=(k_elec_00+0.5d0*k_elec_01)*partial_sigma_to_partial_elec_t
epz=(k_elec_00+0.5d0*k_elec_01)*epz

if(time(2).lt.time(1))then;
ctens=(k00+k01)*ctens;
partial_sigma_to_partial_elec_t=(k_elec_00+k_elec_01)*partial_sigma_to_partial_elec_t
epz=(k_elec_00+k_elec_01)*epz
endif


!!     ===========================componenst of tangent matrixe
      k_coef=0.0d0;
!!     ===========================the tangent matrix
      do 40 k=1,ndf;
      do 40 l=1,ndf;

      do 40 i=1,dimen;
      do 40 j=1,dimen;

      k_coef(k,l)=k_coef(k,l)-beleci(i,l)*ktense(i,j)*belecj(j,k)

      do 40 m=1,dimen;

      k_coef(k,l)=k_coef(k,l)-epz(m,i,j)*beleci(m,l)*bepsilonj(i,j,k)
      k_coef(k,l)=k_coef(k,l)-epz(m,i,j)*belecj(m,k)*bepsiloni(i,j,l)

      do 40 n=1,dimen;
40    k_coef(k,l)=k_coef(k,l)+ctens(i,j,m,n)*bepsilonj(i,j,k)*bepsiloni(m,n,l)

      end subroutine truss_k_gen3d



!     ===============================================================
!     the element coefficient matrix
!     ===============================================================
subroutine elmatrcs_3d(elu,elp,elb,elk,elq)
!     ================================================================
!     element  calculations based on  linear and quadratic rectangular
!     relements with isoparametric  formulation.
!     ================================================================
implicit none
!     ================================================================
!  input variables
!     ================================================================
real(iwp),intent(in):: elu(:),elp(:),elb(:) !dof vector
!     ================================================================
!  otput variables
!     ================================================================
real(iwp),intent(out):: elq(:) !element source vector
real(iwp),intent(out):: elk(:,:) !the element coefficient matrix
!     ================================================================
integer::ni,nj,nk,i,ii,jj,li,lj
integer::inode,jnode
!     ================================================================
real(iwp):: elu_tense(npe,ndf) !dof current time
real(iwp):: elp_tense(npe,ndf) !dof current time
real(iwp):: elb_tense(npe,ndf) !dof current time
real(iwp)::coord(dimen) !the global coordinate of each node
real(iwp)::xi(dimen)! transformed coordinates
real(iwp)::sf(npe),gdsf(dimen,npe)
real(iwp)::ur(ndf,dimen),up(ndf,dimen),ub(ndf,dimen) !value of the function and its derivative
real(iwp)::k_coef(ndf,ndf) ! the auxilary matrix for cuefficitne matr
real(iwp)::res_vect(ndf)
real(iwp)::cnst,det
real(iwp)::der(dimen) !electric displacement
real(iwp)::sigma(dimen,dimen) !stress

real(iwp)::elec_element(dimen) !electric field in element
real(iwp)::pola_element(dimen) !polarization field in element
!     ================================================================
real(iwp):: gauspt(5,5),gauswt(5,5) ! gauss points and gauss weights
!     ================================================================
 data gauspt/5*0.0d0, -0.57735027d0, 0.57735027d0, 3*0.0d0,&
 -0.77459667d0, 0.0d0, 0.77459667d0, 2*0.0d0, -0.86113631d0,     &
 -0.33998104d0, 0.33998104d0, 0.86113631d0, 0.0d0, -0.90617984d0,&
 -0.53846931d0,0.0d0,0.53846931d0,0.90617984d0/
!
data gauswt/2.0d0, 4*0.0d0, 2*1.0d0, 3*0.0d0, 0.55555555d0,  &
  0.88888888d0, 0.55555555d0, 2*0.0d0, 0.34785485d0,   &
2*0.65214515d0, 0.34785485d0, 0.0d0, 0.23692688d0,     &
  0.47862867d0, 0.56888888d0, 0.47862867d0, 0.23692688d0/
!     ===============================================================
elq=0.0d0
elk=0.0d0
!     ===============================================================
    do i=1,npe
    elu_tense(i,:)=elu(1+(i-1)*ndf:i*ndf)
    elp_tense(i,:)=elp(1+(i-1)*ndf:i*ndf)
    elb_tense(i,:)=elb(1+(i-1)*ndf:i*ndf)
    enddo !i=1,npe




     elec_element=0.0d0  !electric field in the element
     pola_element=0.0d0

    do ni = 1,ipdf
    do nj = 1,ipdf
    do nk = 1,ipdf

    xi(1) = gauspt(ni,ipdf)
    xi(2) = gauspt(nj,ipdf)
    xi(3) = gauspt(nk,ipdf)

    gauss_point_number=gauss_point_number+1

    if(npe.eq.20) call shap_3d_quad(elcrds,xi,det,sf,gdsf)
    if(npe.eq.8)  call shapefuncs_3d(xi,det,sf,gdsf)


    cnst = det*gauswt(ni,ipdf)
     coord=0
    do i=1,npe
    coord(:)=coord(:)+elcrds(i,:)*sf(i)
    enddo !i=1,npe
!     ===============================================================
!     post process
!     ===============================================================
    ur=transpose(matmul(gdsf,elu_tense))
    up=transpose(matmul(gdsf,elp_tense))
    ub=transpose(matmul(gdsf,elb_tense))




call stress_elect_displacement(ur,up,ub,der,sigma)
!     ===============================================================
!     filling the elemt matrixes
!     ===============================================================


do inode=1,npe
    ii = ndf*(inode-1)+1
    li=(ii+ndf-1)
    call residual(inode,gdsf,ur,der,sigma,res_vect)
    elq(ii:li)=elq(ii:li)+res_vect*cnst
    do  jnode=1,npe
  call  k_gen(inode,jnode,gdsf,ur,k_coef)
  jj = ndf*(jnode-1)+1
  lj=(jj+ndf-1)
  elk(ii:li,jj:lj)=elk(ii:li,jj:lj)+k_coef*cnst
    enddo;
enddo

!if(noelem.eq.5)then
!  write(10,*)'npe=',npe
!  write(10,*)'cnst=',cnst
!  call show_matrix(gdsf,'gdsf')
! call show_matrix(k_coef,'k_coef')
! call show_vector(res_vect,'res_vect')
! endif

elec_element=elec_element+curn_electric_field(gauss_point_number,:)
pola_element=pola_element+curn_polarization_function(gauss_point_number,:)

    enddo !ni = 1,ipdf
    enddo !nj = 1,ipdf
    enddo !nk = 1,ipdf

!write(1,*)elec_element

elements_electric_field(noelem,:)=elec_element
elements_electric_polar(noelem,:)=pola_element

!if(noelem.eq.5)then
! call show_matrix(elk,'elk')
! call show_vector(elq,'elq')
! endif

end subroutine elmatrcs_3d

!     ===============================================================
!                       the shape function subroutine
!     ===============================================================
      subroutine shap_3d_quad(elcrds,xi_l,det,sf,gdsf)
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real*8,intent(in):: elcrds(:,:) !nodes geometry array of
      real*8,intent(in):: xi_l(:) !gauss point integration point
      !     ================================================================
      !                          output variables
      !     ================================================================
      real(iwp)::det ! the jacobian determinant
      real(iwp),intent(out) ::sf(:),gdsf(:,:) ! shape function sf(npe),gdsf(dimen,npe)
      real(iwp)::dsf(dimen,npe) ! the differentiatian of shape function
      real(iwp)::gj(dimen,dimen)! global jacobian matrix
      real(iwp)::gjinv(dimen,dimen) ! inverse of global jacobian matrix
      real*8,allocatable::sf_master(:),gdsf_master(:,:) ! shape function gdsf(dimen,npe)
!     ================================================================
!                          trivial variables
!     ================================================================
      real*8,allocatable:: elcrds_master(:,:) !nodes geometry array of master element
      real*8::xnode(20,3) ! coordinate of nodes in refrence
      integer,allocatable::ipvt(:)
      integer::i,ni,np(20) !the counter integers
      logical::ok_flag
      real*8::xi,zeta,eta
      real*8::xi0,zeta0,eta0
      real*8::xi1,zeta1,eta1
      real*8::xi2,zeta2,eta2
      real*8::xp,yp,zp
      real*8::fncc,a,b,c
!     ===================================================================
!                              data
!     ===================================================================
    data xnode/ &
     -1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,-1,-1, 1,1,-1, 0,1,0,-1, &
     -1,-1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,-1,-1,1, 1,-1,0,1, 0, &
     -1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0,0, 0, 1,1,1, 1/

    data np/ &
     1 , 2, 3, 4, 5, 6, 7, 8, &
     9 ,11,17,19, &
     10,12,18,20, &
     13,14,15,16/
!     ================================================================




      fncc(a,b,c) = (a*b)*c
      allocate(ipvt(dimen),elcrds_master(npe,dimen), &
      sf_master(npe),gdsf_master(dimen,npe) )
      elcrds_master=elcrds
      elcrds_master(17:20,:)=elcrds(13:16,:)
      elcrds_master(13:16,:)=elcrds(17:20,:)

      xi=xi_l(1)
      eta=xi_l(2)
      zeta=xi_l(3)

      do 20 i = 1, npe


      ni = np(i)
    xp = xnode(ni,1)
    yp = xnode(ni,2)
    zp = xnode(ni,3)
!c
    xi0 = 1.0+xi*xp
    eta0 = 1.0+eta*yp
    zeta0 = 1.0+zeta*zp
!c
    xi1 = 1.0-xi*xi
    eta1 = 1.0-eta*eta
    zeta1 = 1.0-zeta*zeta
!c
    xi2 = xp*xi
    eta2 = yp*eta
    zeta2 = zp*zeta
!c
    if(i.le.8)then
    sf(ni) = 0.125*fncc(xi0,eta0,zeta0)*(xi2+eta2+zeta2-2.0)
!c
    dsf(1,ni) = 0.125*fncc(xp,eta0,zeta0)*(2*xi2+eta2+zeta2-1.0)
    dsf(2,ni) = 0.125*fncc(xi0,yp,zeta0)*(xi2+2*eta2+zeta2-1.0)
    dsf(3,ni) = 0.125*fncc(xi0,eta0,zp)*(xi2+eta2+2*zeta2-1.0)
!c   xi=0 nodes, 9,11,17,19
    else
    if(i.le.12)then
    sf(ni) = 0.25*fncc(xi1,eta0,zeta0)
!c
    dsf(1,ni) = -0.5*fncc(xi,eta0,zeta0)
    dsf(2,ni) = 0.25*fncc(xi1,yp,zeta0)
    dsf(3,ni) = 0.25*fncc(xi1,eta0,zp)
!c   eta=0 nodes, 10,12,18,20
    else
    if(i.le.16)then
    sf(ni) = 0.25*fncc(xi0,eta1,zeta0)
!c
    dsf(1,ni) = 0.25*fncc(xp,eta1,zeta0)
    dsf(2,ni) = -0.5*fncc(xi0,eta,zeta0)
    dsf(3,ni) = 0.25*fncc(xi0,eta1,zp)
!c   zeta=0 nodes, 13,14,15,16
    else
    sf(ni) = 0.25*fncc(xi0,eta0,zeta1)
!c
    dsf(1,ni) = 0.25*fncc(xp,eta0,zeta1)
    dsf(2,ni) = 0.25*fncc(xi0,yp,zeta1)
    dsf(3,ni) = -0.5*fncc(xi0,eta0,zeta)
!c
      endif
      endif
      endif
   20 continue




!     ===================================================================
!         compute the jacobian matrix [gj] and its inverse [gjinv]
!     ===================================================================
      gj=matmul(dsf,elcrds_master)
      gjinv=gj

      call m33inv (gj, gjinv,det, ok_flag)
!     ===================================================================
!                        finding the global derivitive
!     ===================================================================
      gdsf=matmul(gjinv,dsf)
      gdsf_master=gdsf
      sf_master=sf

      gdsf(:,13:16)=gdsf_master(:,17:20)
      gdsf(:,17:20)=gdsf_master(:,13:16)

      sf(13:16)=sf_master(17:20)
      sf(17:20)=sf_master(13:16)

!  920 format (5x,i3,3e12.4,2x,2i5,2x,4i5,2x,2i5)


 end subroutine shap_3d_quad

!     ===============================================================
!     THE SHAPE FUNCTION SUBROUTINE
!     ===============================================================
subroutine shapefuncs_3d(xi,det,sf,gdsf)
real(iwp)::det ! the jacobian determinant
real(iwp),intent(out) ::sf(:),gdsf(:,:) ! shape function sf(npe),gdsf(dimen,npe)
real(iwp)::dsf(dimen,npe) ! the differentiatian of shape function
real(iwp)::gj(dimen,dimen)! global jacobian matrix
real(iwp)::gjinv(dimen,dimen) ! inverse of global jacobian matrix
real(iwp)::xnode(20,3) ! coordinate of nodes in refrence
real(iwp)::xp(dimen),xi(dimen),xi0(dimen) !locad coordinate in element
logical::ok_flag
integer::i
!     ================================================================
data xnode/     &
     -1,1,1,-1,-1,1,1,-1,0,1,0,-1,-1,1,1,-1,0,1,0,-1,   &
     -1,-1,1,1,-1,-1,1,1,-1,0,1,0,-1,-1,1,1,-1,0,1,0,   &
     -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0,1,1,1,1/
     !     ================================================================
     do i=1,npe
         xp(:)=xnode(i,:)
         xi0(:)=1.0+xi(:)*xp(:)
         sf(i) = 0.125*xi0(1)*xi0(2)*xi0(3)
         dsf(1,i)= 0.125*xp(1)*xi0(2)*xi0(3)
         dsf(2,i)= 0.125*xi0(1)*xp(2)*xi0(3)
         dsf(3,i)= 0.125*xi0(1)*xi0(2)*xp(3)
     enddo ! i=1,npe

!     ===================================================================
!   compute the jacobian matrix [gj] and its inverse [gjinv]
!     ===================================================================
gj=matmul(dsf,elcrds)
gjinv=gj

call m33inv (gj, gjinv,det, ok_flag)

if(.not.ok_flag)then
    write(*,*)'determinant of jacobian is too small'
endif

gdsf=matmul(gjinv,dsf)

end subroutine shapefuncs_3d
!     ===============================================================
!     residula vector subroutine
!     ===============================================================
subroutine residual(inode,gdsf,ur,der,sigma,res_vect)
!     ================================================================
implicit none
!     ================================================================
!  input variables
!     ================================================================
integer,intent(in):: inode !node number
real(iwp),intent(in)::gdsf(dimen,npe) !  the shape functions
real(iwp),intent(in),dimension(ndf,dimen):: ur
real(iwp),intent(in),dimension(dimen)::der ! electric displacement tensor
real(iwp),intent(in),dimension(dimen,dimen)::sigma !stress tensor
!     ================================================================
!  output variables
!     ================================================================
real(iwp),intent(out)::res_vect(ndf) !the coefficient matri
!     ================================================================
integer::i,j,k,m
integer::del(dimen,dimen),eye(ndf,ndf)
!     ================================================================
!     deformmation and strain tensors
!     ================================================================
!real(iwp),dimension(dimen,dimen)::udsp_t,udsp_p,udsp_b
!     ===========================shape function tensors
real(iwp)::bepsilon(dimen,dimen,ndf),belec(dimen,ndf)
!     ================================================================
del = 0; do i = 1,dimen; del(i,i) = 1;enddo
eye = 0; do i = 1,ndf; eye(i,i) = 1;enddo
!!     ================================================================
!!c     shape function tensor
!!     ================================================================
belec=0.0d0;
bepsilon=0.0d0;
!

belec(:,dimen+1) = -gdsf(:,inode)
do 10 k=1,ndf
do 10 i=1,dimen
do 10 j=1,dimen
bepsilon(i,j,k)= 0.5d0*gdsf(j,inode)*eye(k,i)  &
   + 0.5d0*gdsf(i,inode)*eye(k,j)
do 10 m=1,dimen
bepsilon(i,j,k)=bepsilon(i,j,k)     &
 + 0.5d0*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
 + 0.5d0*gdsf(j,inode)*eye(k,m)*ur(m,i)
 10   continue

!     ================================================================
res_vect=0.0d0

do 40 k=1,ndf
do 40 i=1,dimen
      res_vect(k)=res_vect(k)-der(i)*belec(i,k)
do 40 j=1,dimen
40    res_vect(k)=res_vect(k)+sigma(i,j)*bepsilon(i,j,k)

end subroutine residual

!     ===============================================================
!     coefficient matrix subroutine
!     ===============================================================
subroutine  k_gen(jnode,inode,gdsf,ur,k_coef)
!     ================================================================
implicit none
!     ================================================================
!  input variables
!     ================================================================
integer,intent(in)::inode,jnode !node number iterators
real(iwp),intent(in)::gdsf(:,:) !  the shape functions
real(iwp),intent(in)::ur(:,:)! displacement
!     ================================================================
!  output variables
!     ================================================================
real(iwp),intent(out)::k_coef(:,:) !the coefficient matrix
!     ================================================================
integer::i,j,k,l,m,n
integer,save::counter
integer::del(dimen,dimen),eye(ndf,ndf)
!     ===========================shape function tensors
real(iwp)::bepsiloni(dimen,dimen,ndf),beleci(dimen,ndf)
real(iwp)::bepsilonj(dimen,dimen,ndf),belecj(dimen,ndf)
!!     ================================================================
counter=counter+1
del = 0; do i = 1,dimen; del(i,i) = 1;enddo
eye = 0; do i = 1,ndf;   eye(i,i) = 1;enddo
!!     ================================================================
!!c     material properties
!!     ================================================================
beleci=0.0d0;
belecj=0.0d0;
bepsiloni=0.0d0;
bepsilonj=0.0d0;

beleci(:,dimen+1) = -gdsf(:, inode)
belecj(:,dimen+1) = -gdsf(:, jnode)

do k=1,ndf
do i=1,dimen
do j=1,dimen

bepsilonj(i,j,k)= 0.5*gdsf(j,jnode)*eye(k,i)+0.5*gdsf(i,jnode)*eye(k,j)
bepsiloni(i,j,k)= 0.5*gdsf(j,inode)*eye(k,i)+0.5*gdsf(i,inode)*eye(k,j)

do m=1,dimen

bepsilonj(i,j,k)=bepsilonj(i,j,k) &
 + 0.5*gdsf(i,jnode)*eye(k,m)*ur(m,j)   &
 + 0.5*gdsf(j,jnode)*eye(k,m)*ur(m,i)
bepsiloni(i,j,k)=bepsiloni(i,j,k) &
 + 0.5*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
 + 0.5*gdsf(j,inode)*eye(k,m)*ur(m,i)

enddo ! m=1,dimen
enddo ! 10 k=1,ndf
enddo ! 10 i=1,dimen
enddo ! 10 j=1,dimen

call material_properties()

!     ================================================================

ctens=(k00+0.5d0*k01)*ctens
partial_sigma_to_partial_elec_t=(k_elec_00+0.5d0*k_elec_01)*partial_sigma_to_partial_elec_t
epz=(k_elec_00+0.5d0*k_elec_01)*epz

if(time(2).lt.time(1))then;
ctens=(k00+k01)*ctens;
partial_sigma_to_partial_elec_t=(k_elec_00+k_elec_01)*partial_sigma_to_partial_elec_t
epz=(k_elec_00+k_elec_01)*epz
endif
!     ===========================componenst of tangent matrixe
k_coef=0.0d0;
!     ===========================the tangent matrix
do 40 k=1,ndf;
do 40 l=1,ndf;

do 40 i=1,dimen;
do 40 j=1,dimen;

k_coef(k,l)=k_coef(k,l)-beleci(i,l)*ktense(i,j)*belecj(j,k)

do 40 m=1,dimen;

k_coef(k,l)=k_coef(k,l)-partial_sigma_to_partial_elec_t(m,i,j)*beleci(m,l)*bepsilonj(i,j,k)
k_coef(k,l)=k_coef(k,l)-partial_sigma_to_partial_elec_t(m,i,j)*belecj(m,k)*bepsiloni(i,j,l)

do 40 n=1,dimen;
40    k_coef(k,l)=k_coef(k,l)+ctens(i,j,m,n)*bepsilonj(i,j,k)*bepsiloni(m,n,l)

end subroutine k_gen

!     ===============================================================
!                  the boundary condition subroutine
!     ===============================================================
     subroutine symmetric_primary_bounday(glk,glr)
     implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real(iwp),intent(inout) :: glk(:,:),glr(:)! global coefficient
      integer::i !integer counters
      real(iwp),allocatable::gls(:)
!     ================================================================
!     ================================================================
!                          primary variables
!     for u(i)=alpha, it puts k(i,j)=1 and f(i)=alpha
!     ================================================================
      allocate(gls(size(glr)));
        gls=0.0d0
        gls(bnd_no_pr_vec)=vspv
        glr=glr-matmul(glk,gls)

        glk(bnd_no_pr_vec,:)=0.0d0;
        glk(:,bnd_no_pr_vec)=0.0d0;

        do i=1,nspv;
        glk(bnd_no_pr_vec(i),bnd_no_pr_vec(i))=1.0d0; enddo
        glr(bnd_no_pr_vec)=vspv

    end subroutine symmetric_primary_bounday

subroutine result_printer(iter,glu)
real(iwp),INTENT(IN)::glu(:)

integer::iter;
!     ================================================================
!  trivial variables
!     ================================================================
integer::i,j,pdf,k;
! integer::curved_node;
integer,parameter::gnuplot=5
integer,save::counter
character(len=5)::x1
character(len=30) :: fmt,filename ! format descriptor
fmt = '(i4.4)' ! an integer of width 5 with zeros at the left

! write(*,*)dimen
!write(*,*)maxval(glu)
!write(*,*)maxloc(glu)

curved_node=78
!curved_node=nnm
counter=counter+1


write(x1,fmt)i_calibration

if (time(2).eq.0)then

open (out,file='out_data_outpu.txt')

filename='./gnuplot/gnuplot'//trim(x1)//'.gnu'
open (gnuplot,file=filename)

filename='./excel/results_'//trim(x1)//'.csv'
open (csv,file=filename)



endif

if(counter.eq.1)then
!write(csv,672)
endif

write(out,910);write(out,*)'time',time
write(out,910);write(out,910);write(out,670);write(out,910)
do i=1,nnm
pdf=(i-1)*ndf
! write(out,950)i,(coords(i,j),j=1,dimen),(glu(pdf+k),k=1,ndf)
enddo

!write(*,*)'curved_node',curved_node
pdf=(curved_node-1)*ndf

write(csv,951)iter,time(2),loadfactor,glu(pdf+ndf-1) ! ,curn_polarization_function(3,:)
write(gnuplot,951)iter,time(2),loadfactor,glu(pdf+ndf-1) 

write(out,*)
!     ===============================================================
!                                 formats
!     ===============================================================
910   format (2x,80('_'),/)
670   format(5x,'node     x-coord. u-solut.')
950   format(5x,i5,10e14.5)
951   format(5x,i5,' , ',10(e14.5,' , '))
!672   format(5x,'itr, time ,u1-solut,  u2-solut,loadfactor')
672   format(5x,'itr, time , potential ,strain(%)')

end subroutine result_printer

end module fem_libs


