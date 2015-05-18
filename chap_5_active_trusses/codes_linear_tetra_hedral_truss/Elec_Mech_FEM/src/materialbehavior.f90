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

module material_behavior
use linearsolvers
! ================================================================
real(iwp)::k00
real(iwp)::k01
real(iwp)::lambda_01

real(iwp)::k_elec_00
real(iwp)::k_elec_01
real(iwp)::lambda_elec_01

! =============================pzt
real(iwp),parameter::k00_pzt=1.00d0
real(iwp),parameter::k01_pzt=0.0d0
real(iwp),parameter::lambda_01_pzt=1.0

real(iwp),parameter::k_elec_00_pzt=1.0d0
real(iwp),parameter::k_elec_01_pzt=-0.0d0
real(iwp),parameter::lambda_elec_01_pzt=1.5

! =============================epoxy
real(iwp),parameter::k00_epx=1.0d0
real(iwp),parameter::k01_epx=0.0d0
real(iwp),parameter::lambda_01_epx=0.8

real(iwp),parameter::k_elec_00_epx=1.0d0
real(iwp),parameter::k_elec_01_epx=0.0d0
real(iwp),parameter::lambda_elec_01_epx=1.0
! ================================================================
real(iwp),parameter::ec=1.2e6;
real(iwp),parameter::pr_sat=0.25;
real(iwp),parameter::m=4.0;

real(iwp),parameter::eta=0.1e0;
real(iwp),parameter::h0=1.0e0! 100.0;


real(iwp),parameter::eps_s_r=1e-3;
real(iwp),parameter::poarization_speed=5 !1/s
! ================================================================ pzt


real(iwp),parameter::e33 = -14.52*1.3d0;
real(iwp),parameter::e15 = -7.56*1.3d0;
real(iwp),parameter::e31 = 6.15*1.3d0;

real(iwp),parameter::eps_11 = 4.0e-9;
real(iwp),parameter::eps_33 = 2.0e-9;

real(iwp),parameter::ey_pzt =60.6e9 !Pa
real(iwp),parameter::nu_pzt = 0.3;

real(iwp),parameter::ey_epoxy =  1.5e9 ! pa
real(iwp),parameter::nu_epoxy =  0.35;
real(iwp),parameter::eps_epoxy = 8.854187e-9; !  with eps_r=1000 F/ micro m

! ================================================================ aluminum
real(iwp),parameter::ey_aluminum =69.0e9 !N/micro m !69.0e9  N/m!
real(iwp),parameter::nu_aluminum = 0.33;

real(iwp),parameter::k_00_aluminum=1.0d0
real(iwp),parameter::k_01_aluminum=0.0d0
real(iwp),parameter::lambda_01_aluminum=4.0
! ================================================================ aluminum

logical::is_polarized
real(iwp)::a
real(iwp)::direc_a(dimen)
real(iwp) :: element_direction_normal(dimen)


real(iwp)::lambda
real(iwp)::mu
real(iwp)::alpha1
real(iwp)::alpha2
real(iwp)::alpha3
real(iwp)::beta1
real(iwp)::beta2
real(iwp)::beta3
real(iwp)::gamma1
real(iwp)::gamma2

real(iwp)::partial_sigma_to_partial_elec_t(dimen,dimen,dimen)
real(iwp)::partial_sigma_to_partial_elec_p(dimen,dimen,dimen)
real(iwp)::b_tilt(dimen,dimen,dimen,dimen)

! =========================element center point
real(iwp) :: el_center_coord(dimen)
real(iwp) :: qtransform(dimen,dimen) ! this is the transformtion matrix
! =======================mechanical properties variables
real(iwp)::ctens(dimen,dimen,dimen,dimen)
! ===========================materials tensors
real(iwp),dimension(dimen,dimen,dimen):: epz
real(iwp),dimension(dimen,dimen) ::ktense
! ================================================================
!  time variables
! ================================================================
real(iwp)::time(2),dtime
integer::ntime,itime  ! time
! ================================================================
integer::noelem
! ================================================================
!  history variables
! ================================================================
real(iwp),allocatable::remanent_polarization_vector(:,:)

real(iwp),allocatable::hist_polarization_function(:,:)
real(iwp),allocatable::curn_polarization_function(:,:)

real(iwp),allocatable::mechanical_hist_curentt(:,:,:,:,:)
real(iwp),allocatable::mechanical_hist_previus(:,:,:,:,:)
real(iwp),allocatable::mechanical_hist_beforep(:,:,:,:,:)


real(iwp),allocatable::sigma_el_hist_curentt(:,:,:,:)
real(iwp),allocatable::sigma_el_hist_previus(:,:,:,:)
real(iwp),allocatable::sigma_el_hist_beforep(:,:,:,:)

real(iwp),allocatable::displ_el_hist_curentt(:,:,:,:)
real(iwp),allocatable::displ_el_hist_previus(:,:,:,:)
real(iwp),allocatable::displ_el_hist_beforep(:,:,:,:)

real(iwp),allocatable::mechanical_strain_curentt(:,:,:)
real(iwp),allocatable::mechanical_strain_previus(:,:,:)
real(iwp),allocatable::mechanical_strain_beforep(:,:,:)

real(iwp),allocatable::hist_electric_field(:,:)
real(iwp),allocatable::curn_electric_field(:,:)

real(iwp),allocatable::hist_electric_displ(:,:)
real(iwp),allocatable::curn_electric_displ(:,:)

real(iwp)::total_displacement(dimen)
real(iwp)::total_polarization(dimen)
real(iwp)::pr(dimen)


! the output parameters
real(iwp),allocatable::elements_electric_field(:,:) !(nem,dimen)
real(iwp),allocatable::elements_electric_polar(:,:) !(nem,dimen)


real(iwp),allocatable::hist_back_electric_field(:,:)
real(iwp),allocatable::curn_back_electric_field(:,:)

real(iwp)::electric_field(dimen),electric_field_p(dimen),electric_field_b(dimen);

contains


subroutine stress_elect_displacement(ur,up,ub,der,sigma)
implicit none
! ================================================================
!  input variables
! ================================================================
real(iwp),intent(in) :: ur(:,:),up(:,:),ub(:,:) ! values of field function n point
! ================================================================
! output variables
! ================================================================
integer :: i , j, k , l !integer counters
real(iwp)  :: der(:),sigma(:,:)
integer::eye(dimen,dimen)
! ================================================================
real(iwp)::strain(dimen,dimen),strain_p(dimen,dimen),strain_b(dimen,dimen);
real(iwp)::d_strain(dimen,dimen),d_strain_p(dimen,dimen)
real(iwp)::d_electric_field(3),d_electric_field_p(3)
real(iwp)::strain_r(dimen,dimen),strain_elastic(dimen,dimen);

! ================================================================
eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
! ================================================================
strain=0.0d0;electric_field=0.0d0;
strain_p=0.0d0;electric_field_p=0.0d0;
strain_b=0.0d0;electric_field_b=0.0d0;

strain(1:dimen,1:dimen)=0.5d0*(   ur(1:dimen,:)+transpose(ur(1:dimen,:))+  &
matmul(  transpose(   ur(1:dimen,:)  ),ur(1:dimen,:)));

strain_p(1:dimen,1:dimen)=0.5d0*(   up(1:dimen,:)+transpose(up(1:dimen,:))+  &
matmul(  transpose(   up(1:dimen,:)  ),up(1:dimen,:)));

strain_b(1:dimen,1:dimen)=0.5d0*(   ub(1:dimen,:)+transpose(ub(1:dimen,:))+  &
matmul(  transpose(   ub(1:dimen,:)  ),ub(1:dimen,:)));

electric_field(1:dimen)=-ur(dimen+1,:)
electric_field_p(1:dimen)=-up(dimen+1,:)
electric_field_b(1:dimen)=-ub(dimen+1,:)

d_electric_field=electric_field-electric_field_p
d_electric_field_p=electric_field_p-electric_field_b
! ================================================================
curn_electric_field(gauss_point_number,:)=electric_field
hist_electric_field(gauss_point_number,:)=electric_field_p
! ================================================================
!curn_polarization_function=0.0d0
call material_properties_afc()


! ===================================remanent strain
strain_r=0.0d0


do i=1,dimen
do j=1,dimen
strain_r(i,j)=3.0d0*0.5d0*eps_s_r*(direc_a(i)*direc_a(j)-eye(i,j)/3.0d0)*norm_vect(pr)/pr_sat
end do;
end do;

!write(1,*)strain_r
strain_r=0.0d0


strain_elastic=strain-strain_r
! ==============================================increment in strain
mechanical_strain_curentt(gauss_point_number,:,:)=strain_elastic
d_strain=mechanical_strain_curentt(gauss_point_number,:,:)-mechanical_strain_previus(gauss_point_number,:,:)
d_strain_p=mechanical_strain_previus(gauss_point_number,:,:)-mechanical_strain_beforep(gauss_point_number,:,:)

!!> derivative of stress with respect to electric field.
! ==============================the history variable
do i=1,dimen
do j=1,dimen
do k=1,dimen

sigma_el_hist_curentt(gauss_point_number,i,j,k)= &
sigma_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime)+ &
(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime)*d_electric_field_p(k)*partial_sigma_to_partial_elec_p(k,i,j) &
                          +d_electric_field(k)  *partial_sigma_to_partial_elec_t(k,i,j) &
)

displ_el_hist_curentt(gauss_point_number,i,j,k)= &
displ_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime)+ &
(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime)*d_strain_p(i,j)*partial_sigma_to_partial_elec_p(k,i,j)+ &
                           d_strain(i,j)  *partial_sigma_to_partial_elec_t(k,i,j) &
)


do l=1,dimen
mechanical_hist_curentt(gauss_point_number,i,j,k,l)= &
mechanical_hist_previus(gauss_point_number,i,j,k,l)*exp(-lambda_01*dtime)+ &
ctens(i,j,k,l)*(k01*0.5d0)* &
(  &
exp(-lambda_01*dtime)*d_strain_p(k,l)+d_strain(k,l) &
)
end do;
end do;
end do;
end do;



sigma =  0.0d0 ;
der   =  0.0d0 ;

do i=1,dimen
do j=1,dimen
der(i)=der(i)+ktense(i,j)*electric_field(j)
do k=1,dimen
   der(k)  =  der(k)    +  k_elec_00*partial_sigma_to_partial_elec_t(k,i,j)*( strain_elastic(i,j) ) &
                        +  displ_el_hist_curentt(gauss_point_number,i,j,k)

sigma(i,j) =  sigma(i,j)  -  k_elec_00*partial_sigma_to_partial_elec_t(k,i,j)*  electric_field(k) &
                          -  sigma_el_hist_curentt(gauss_point_number,i,j,k)
do l=1,dimen
sigma(i,j)=sigma(i,j)+k00*ctens(i,j,k,l)*(  strain_elastic(k,l) )  &
         + mechanical_hist_curentt(gauss_point_number,i,j,k,l)

!!>> Nonlinear parts of electric stress and electric displacement field


end do;
end do;
end do;
end do;

end subroutine stress_elect_displacement


! ================================================================
!   polarization swithing function
!this will find the polarization with the given state of material
! ================================================================
subroutine remanent_polarization()
implicit none
real(iwp)::d_electric_field(3)
real(iwp)::electric_field(3)
real(iwp)::electric_drive(dimen)
real(iwp)::electric_field_t(3)
real(iwp)::electric_field_p(3)
real(iwp)::back_electric_field_t(dimen)
real(iwp)::back_electric_field_p(dimen)
real(iwp)::normalized_polarization(dimen)

!>newton raphson parameters
real(iwp)::ktan_pol(6,6),del_u_n_pol(6)
real(iwp)::r_n_pol(6),u_n_pol(6),force_n_pol(6)
real(iwp)::u_n_r_pol(6),u_n_0_pol(6)
real(iwp)::ktan_inv_pol(6,6)
!!<Iteration parameters
integer::polarization_iteration_number
integer,parameter::max_iteration=50
integer::required_number_of_iteration
real(iwp),parameter::tolerance=1.0e-2
logical::converged
real(iwp):: polarization_error

!>polarization parameters
real(iwp):: electric_polarization_t(dimen)
real(iwp):: electric_polarization_p(dimen)


electric_field=curn_electric_field(gauss_point_number,:)
electric_field_p=hist_electric_field(gauss_point_number,:)


electric_field_t=electric_field
electric_field_p=electric_field-d_electric_field

electric_polarization_p=hist_polarization_function(gauss_point_number,:)
electric_polarization_t=curn_polarization_function(gauss_point_number,:)

back_electric_field_p=hist_back_electric_field(gauss_point_number,:)
back_electric_field_t=curn_back_electric_field(gauss_point_number,:)

u_n_pol=[electric_polarization_p/pr_sat,back_electric_field_p]


do polarization_iteration_number=1,max_iteration

back_electric_field_t=u_n_pol(4:6)
normalized_polarization=u_n_pol(1:3)
electric_drive=electric_field_t/ec-back_electric_field_t

r_n_pol(1:3)=normalized_polarization- electric_polarization_p/pr_sat &
-( dtime/eta)*( ( macaulay_brackets ( norm_vect(electric_drive ) -1 )  ) ** m ) * &
                                      unit_vect( electric_drive )
r_n_pol(4:6)=back_electric_field_t &
- h0*atanh( norm_vect( normalized_polarization ) )*unit_vect( normalized_polarization )


polarization_error=norm_vect(r_n_pol)

ktan_pol=identity_matrix(6)
ktan_pol(1:3,4:6)=  &
(m*dtime/eta)*( ( macaulay_brackets ( norm_vect(electric_drive ) -1 ) ) ** (m-1) ) &
*dyad_of_two_vector( unit_vect( electric_drive ), &
                     unit_vect( electric_drive ) ) &
+(dtime/eta)*( ( macaulay_brackets ( norm_vect(electric_drive ) -1 )  ) ** m )  &
* derivative_of_direction(electric_drive )

ktan_pol(4:6,1:3)=( h0 /(norm_vect( normalized_polarization )**2.0-1 ) ) &
*dyad_of_two_vector( unit_vect( normalized_polarization ), &
                     unit_vect( normalized_polarization ) ) &
+h0*atanh( norm_vect( normalized_polarization )) &
*derivative_of_direction(normalized_polarization)


call Gaussian_Elimination_Solver(ktan_pol,r_n_pol);del_u_n_pol=r_n_pol
u_n_pol=u_n_pol-del_u_n_pol


write(23,951)polarization_iteration_number,time(2),electric_field_t
write(22,951)polarization_iteration_number,time(2),u_n_pol,polarization_error

         if (polarization_error.le.tolerance*ec)then
             converged=.true.
             exit
         endif
enddo ! iteration_number=1,max_iteration



if(.not.converged)then
write(*,*)'********************************************'
write(*,*)'iteration did not converge at time',time(2)
stop
endif


curn_polarization_function(gauss_point_number,:)=u_n_pol(1:3)*pr_sat
curn_back_electric_field(gauss_point_number,:)=u_n_pol(4:6)

951   format(5x,i5,' , ',10 (e14.5,' , ') )
end subroutine remanent_polarization



! ===============================================================
!   forming history variables
! ===============================================================
subroutine  clear_history()


elements_electric_field=0.0d0
elements_electric_polar=0.0d0


mechanical_hist_curentt=0.0d0
mechanical_hist_previus=0.0d0
mechanical_hist_beforep=0.0d0


         displ_el_hist_curentt=0.0d0
         displ_el_hist_previus=0.0d0
         displ_el_hist_beforep=0.0d0

         sigma_el_hist_curentt=0.0d0
         sigma_el_hist_previus=0.0d0
         sigma_el_hist_beforep=0.0d0

mechanical_strain_curentt=0.0d0;
mechanical_strain_previus=0.0d0;
mechanical_strain_beforep=0.0d0;


remanent_polarization_vector=0.0d0

hist_polarization_function=0.0d0
curn_polarization_function=0.0d0

hist_electric_field=0.0d0
curn_electric_field=0.0d0

hist_electric_displ=0.0d0
curn_electric_displ=0.0d0

end subroutine  clear_history

! ===============================================================
!   forming history variables
! ===============================================================
subroutine  form_history(ngauss,nem)
integer::ngauss,nem
! ===================allocate history variables
allocate( remanent_polarization_vector(ngauss,dimen) ,&
hist_polarization_function(ngauss,dimen), &
curn_polarization_function(ngauss,dimen), &
hist_electric_field (ngauss,dimen),&
curn_electric_field (ngauss,dimen),&
hist_electric_displ (ngauss,dimen),&
curn_electric_displ (ngauss,dimen))

allocate(mechanical_hist_curentt(ngauss,dimen,dimen,dimen,dimen), &
 mechanical_hist_previus(ngauss,dimen,dimen,dimen,dimen), &
 mechanical_hist_beforep(ngauss,dimen,dimen,dimen,dimen)  )


allocate(sigma_el_hist_curentt(ngauss,dimen,dimen,dimen) ,&
         sigma_el_hist_previus(ngauss,dimen,dimen,dimen) ,&
         sigma_el_hist_beforep(ngauss,dimen,dimen,dimen) )

allocate(displ_el_hist_curentt(ngauss,dimen,dimen,dimen) ,&
         displ_el_hist_previus(ngauss,dimen,dimen,dimen) ,&
         displ_el_hist_beforep(ngauss,dimen,dimen,dimen) )

allocate(  &
mechanical_strain_curentt(ngauss,dimen,dimen) ,&
mechanical_strain_previus(ngauss,dimen,dimen) ,&
mechanical_strain_beforep(ngauss,dimen,dimen)  &
)

allocate(elements_electric_field(nem,dimen) ,& !(nem,dimen)
         elements_electric_polar(nem,dimen) ) !(nem,dimen))

elements_electric_field=0.0d0
elements_electric_polar=0.0d0


mechanical_hist_curentt=0.0d0
mechanical_hist_previus=0.0d0
mechanical_hist_beforep=0.0d0


         displ_el_hist_curentt=0.0d0
         displ_el_hist_previus=0.0d0
         displ_el_hist_beforep=0.0d0

         sigma_el_hist_curentt=0.0d0
         sigma_el_hist_previus=0.0d0
         sigma_el_hist_beforep=0.0d0

mechanical_strain_curentt=0.0d0;
mechanical_strain_previus=0.0d0;
mechanical_strain_beforep=0.0d0;


remanent_polarization_vector=0.0d0

hist_polarization_function=0.0d0
curn_polarization_function=0.0d0

hist_electric_field=0.0d0
curn_electric_field=0.0d0

hist_electric_displ=0.0d0
curn_electric_displ=0.0d0
is_polarized=.false.

allocate(hist_back_electric_field(ngauss,dimen) , &
         curn_back_electric_field(ngauss,dimen)   )

hist_back_electric_field=0.0d0
curn_back_electric_field=0.0d0


end subroutine  form_history

! ===============================================================
!   forming history variables
! ===============================================================
subroutine  update_history()


hist_polarization_function=curn_polarization_function
hist_electric_field =curn_electric_field
hist_electric_displ =curn_electric_displ


mechanical_hist_beforep=mechanical_hist_previus
mechanical_hist_previus=mechanical_hist_curentt

sigma_el_hist_beforep=sigma_el_hist_previus
sigma_el_hist_previus=sigma_el_hist_curentt

displ_el_hist_beforep=displ_el_hist_previus
displ_el_hist_previus=displ_el_hist_curentt

mechanical_strain_beforep=mechanical_strain_previus
mechanical_strain_previus=mechanical_strain_curentt

hist_back_electric_field=curn_back_electric_field

end subroutine  update_history

! ===============================================================
!   Material properties
! ===============================================================
subroutine material_properties()
implicit none
! ================================================================
!  material variables
! ================================================================

integer::eye(dimen,dimen)
integer::i,j,k,l ! ,m,n
real(iwp)::ey,nu

! ================================================================

eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
!

direc_a=unit_vect(pr)

!3D formulation
ctens=0;
epz=0.0d0;
b_tilt=0.0d0;
ktense=0.0d0

beta1   =-e31;
beta2   =-e33+2.0d0*e15+e31;
beta3   =-2.0d0*e15;
gamma1  =-eps_11/2.0d0;
gamma2  =(eps_11-eps_33)/2.0d0;


k00=k00_pzt
k01=k01_pzt
lambda_01=lambda_01_pzt

k_elec_00=k_elec_00_pzt
k_elec_01=k_elec_01_pzt
lambda_elec_01=lambda_elec_01_pzt

ey=ey_pzt;
nu=nu_pzt;

mu=ey/(1+nu)/2.0
lambda=ey*nu/(1+nu)/(1-2.0*nu)

do i = 1,dimen;do j = 1,dimen;do k = 1,dimen;do l = 1,dimen;
ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+ &
   mu*( eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)  )
enddo;enddo;enddo;enddo;

do i = 1,dimen; do k = 1,dimen; do l = 1,dimen;
epz(i,k,l)=( beta1*direc_a(i)*eye(k,l) &
            + beta2*direc_a(i)*direc_a(k)*direc_a(l) &
            + beta3*0.5d0*( eye(i,l)*direc_a(k) + eye(i,k)*direc_a(l) )  ) *norm_vect(pr)/pr_sat
enddo;enddo;enddo;

ktense(1,1)=eps_11
ktense(2,2)=eps_11 ;
ktense(3,3)=eps_33 ;


b_tilt=0.0d0

partial_sigma_to_partial_elec_t=epz
partial_sigma_to_partial_elec_p=epz


end subroutine material_properties

!     ================================================================
!     finding stress due to shape change
!     ================================================================
subroutine truss_stress_elect_displacement(ur,der,sigma)
!     ================================================================
!                          input variables
!     ================================================================
real(iwp),intent(in) :: ur(:,:)! global coefficient
!     ================================================================
!                         output variables
!     ================================================================
integer :: i,j,k,l !integer counters
real(iwp)  :: der(:),sigma(:,:)
!     ================================================================
real(iwp)::strain(3,3);
real(iwp)::electric_feild(3);
!     ================================================================
strain=0.0d0;electric_feild=0.0d0;
strain(1:dimen,1:dimen)=0.5d0*(   ur(1:dimen,:)+transpose(ur(1:dimen,:))+      &
                              matmul(  transpose(   ur(1:dimen,:)  ),ur(1:dimen,:))    );

electric_feild(1:dimen)=-ur(dimen+1,:)

!write(1,*)'strain=',strain(1,1)
!write(1,*)'electric_feild=',electric_feild(1)

sigma =  0.0d0 ;
der   =  0.0d0 ;

call truss_material_properties()

do i=1,dimen; do j=1,dimen;
der(i)=der(i)+ktense(i,j)*electric_feild(j)
do k=1,dimen;
sigma(i,j)=sigma(i,j)-epz(i,j,k)*electric_feild(k)
der(k)=der(k)+epz(i,j,k)*strain(i,j)
do l=1,dimen
sigma(i,j)=sigma(i,j)+ctens(i,j,k,l)*strain(k,l)
enddo;enddo;enddo;enddo;

der=0.0d0

end subroutine truss_stress_elect_displacement


      subroutine truss_shape_change_stress_bending(el_center_coord,ur,sigma_shape)
      implicit none
!     ================================================================
!     input variables
!     ================================================================
      real(iwp), intent(in) :: el_center_coord(:),ur(:,:)
!     ================================================================
!     output variables
!     ================================================================
      real(iwp),intent(out)::sigma_shape(dimen,dimen)
!     ================================================================
!     trivial variables
!     =========================element center point
      real(iwp)::x(dimen)
!     =========================the stress field in the point
      real(iwp)::lagrange_strain_global(dimen,dimen),e(dimen,dimen)
      integer::i,j,k,l
!     =========================geometrical variables
      real(iwp)::lenght_beam
      real(iwp)::pi
      real(iwp)::radius_of_curvature,curvature_in_time
      real(iwp)::x1,x2,x3,r,kappa,epsilon
real(iwp)::strain(3,3);
real(iwp)::electric_feild(3);
real(iwp)::volumetric_stress
!     ================================================================
strain=0.0d0;electric_feild=0.0d0;
strain(1:dimen,1:dimen)=0.5d0*(   ur(1:dimen,:)+transpose(ur(1:dimen,:))+      &
                              matmul(  transpose(   ur(1:dimen,:)  ),ur(1:dimen,:)) );

electric_feild(1:dimen)=-ur(dimen+1,:)
volumetric_stress=strain(1,1)+strain(2,2)+strain(3,3)

!call show_matrix(strain,'strain')
!write(1,*)volumetric_stress
!     ================================================================
      pi=datan(1.0d0)*4.0
!     =========================lagrange green strain tensor
      x=el_center_coord
      lenght_beam=1.0d0


      radius_of_curvature=lenght_beam /pi/2.0d0

      x1=x(1)
      x2=x(2)
      x3=x(3)
      R=radius_of_curvature
      e=0.0d0
      curvature_in_time=time(2)/r
      kappa=curvature_in_time
      epsilon=time(2)*1.0d0

e=0.0d0
e(1,1)=x3*time(2)
e=e-identity_matrix(dimen)*time(2)/12.0
!lagrange_strain_global=e


call truss_material_properties()

!write(10,*)'time',time
!call show_matrix(e,'e')

sigma_shape=0.0d0
do i=1,dimen;
do j=1,dimen;
do k=1,dimen;
do l=1,dimen;
sigma_shape(i,j)=sigma_shape(i,j)+ctens(i,j,k,l)*e(k,l)
enddo;enddo;enddo;enddo;

!
!volumetric_stress=sigma_shape(1,1)+sigma_shape(2,2)+sigma_shape(3,3)
!!call show_matrix(sigma_shape,'sigma_shape_before')
!
!!sigma_shape(3,:)=0.0d0
!!sigma_shape(:,3)=0.0d0
!
!sigma_shape=sigma_shape-identity_matrix(dimen)*volumetric_stress/3.0

!sigma_shape=sigma_shape-identity_matrix(dimen)*volumetric_stress*0.0001
!call show_matrix(sigma_shape,'sigma_shape')

end subroutine truss_shape_change_stress_bending


!     ================================================================
!     finding stress due to shape change
!     ================================================================
      subroutine shape_change_stress_beam_folding(el_center_coord,sigma_shape)
      implicit none
!     ================================================================
!     input variables
!     ================================================================
      real(iwp), intent(in) :: el_center_coord(:)
!     ================================================================
!     output variables
!     ================================================================
      real(iwp),intent(out)::sigma_shape(dimen,dimen)
!     ================================================================
!     trivial variables
!     =========================element center point
      real(iwp)::x(dimen)
!     =========================the stress field in the point
      real(iwp)::lagrange_strain_global(dimen,dimen),e(dimen,dimen)
      integer::i,j,k,l
!     =========================geometrical variables
      real(iwp)::lenght_beam
      real(iwp)::pi
      real(iwp)::radius_of_curvature,curvature_in_time
      real(iwp)::x1,x2,x3,r,kappa,epsilon
!     ================================================================
      pi=datan(1.0d0)*4.0
!     =========================lagrange green strain tensor
      x=el_center_coord
      lenght_beam=1.0d0


      radius_of_curvature=lenght_beam /pi/2.0d0

      x1=x(1)
      x2=x(2)
      x3=x(3)
      r=radius_of_curvature

       e=0.0d0
!      e(1,1) = cos(X1 / R) ** 2 * X2 / R
!      e(1,2) = -cos(X1 / R) * X2 / R * sin(X1 / R)
!      e(2,1) = -cos(X1 / R) * X2 / R * sin(X1 / R)
!      e(2,2) = sin(X1 / R) ** 2 * X2 / R

!      e(2,2)=0.3
!      e(2,2)=0.3

!       e(3,2)=R
!       e(2,3)=-e(3,2)

      curvature_in_time=time(2)/r
      kappa=curvature_in_time

!e(1,1)=-x3*kappa+0.5*x3*x3*kappa*kappa

!if(abs(x3**2).lt.height_of_beam*0.01)then
!e(1,1)=0
!endif


!e(1,1)=kappa*0.5
!e(1,1)=-x3*kappa+0.5*x3*x3*kappa*kappa
!e(3,3)=epsilon

!e(1,1) =epsilon+0.5d0*epsilon**2.0d0-x3*kappa-2*x3*kappa*epsilon-x3*kappa*epsilon**2.0d0+ &
!        0.5d0*x3**2.0d0*kappa**2.0d0+x3**2.0d0*kappa**2.0d0*epsilon+ &
!        0.5d0*x3**2.0d0*kappa**2.0d0*epsilon**2.0d0
!e(2,2)=-e(1,1)
!e(3,3)=-e(1,1)




e(1,1)=x3*kappa
!
e(2,2)=-x3*kappa/2
e(3,3)=-x3*kappa/2

e=e+identity_matrix(dimen)*time(2)*0.5
!e(2,2)=kappa
!e(3,3)=kappa

!each_truss_strain(noelem)=e(1,1)

!e=e+identity_matrix(dimen)*kappa

      lagrange_strain_global=e



!      lagrange_strain_global=time(2)* &
!      matmul( matmul( R_rotation_tensor ,e ), transpose(R_rotation_tensor))
!!     =========================test for extension

! call show_matrix(lagrange_strain_global,'lagrange_strain_global')

sigma_shape=0.0d0
do i=1,dimen;
do j=1,dimen;
do k=1,dimen;
do l=1,dimen;
sigma_shape(i,j)=sigma_shape(i,j)+ctens(i,j,k,l)*lagrange_strain_global(k,l)

enddo;enddo;enddo;enddo;

!each_truss_strain(noelem)=sigma_shape(1,1)



end subroutine shape_change_stress_beam_folding


      subroutine truss_shape_change_stress_z_eq_xy(el_center_coord,sigma_shape)
      implicit none
!     ================================================================
!     input variables
!     ================================================================
      real(iwp), intent(in) :: el_center_coord(:)
!     ================================================================
!     output variables
!     ================================================================
      real(iwp),intent(out)::sigma_shape(dimen,dimen)
!     ================================================================
!     trivial variables
!     =========================element center point
      real(iwp)::x(dimen)
!     =========================the stress field in the point
      real(iwp)::lagrange_strain_global(dimen,dimen),e(dimen,dimen)
      integer::i,j,k,l
!     =========================geometrical variables
      real(iwp)::lenght_beam
      real(iwp)::pi
      real(iwp)::radius_of_curvature,curvature_in_time
      real(iwp)::x1,x2,x3,r,kappa,epsilon
real(iwp)::strain(3,3);
real(iwp)::electric_feild(3);
real(iwp)::volumetric_stress
!     ================================================================

!     ================================================================
      pi=datan(1.0d0)*4.0
!     =========================lagrange green strain tensor
      x=el_center_coord

      e=0.0d0

      e(1,1) = x(2) ** 2 / 0.2d1
      e(1,2) = x(1) * x(2) / 0.2d1
      e(1,3) = x(2) / 0.2d1
      e(2,1) = x(1) * x(2) / 0.2d1
      e(2,2) = x(1) ** 2 / 0.2d1
      e(2,3) = x(1) / 0.2d1
      e(3,1) = x(2) / 0.2d1
      e(3,2) = x(1) / 0.2d1


e=e*time(2)

!e=e+identity_matrix(dimen)*time(2)*0.001

call truss_material_properties()

sigma_shape=0.0d0
do i=1,dimen;
do j=1,dimen;
do k=1,dimen;
do l=1,dimen;
sigma_shape(i,j)=sigma_shape(i,j)+ctens(i,j,k,l)*e(k,l)
enddo;enddo;enddo;enddo;

!if (abs(x(3)).le.0.001)sigma_shape=0.0d0

end subroutine truss_shape_change_stress_z_eq_xy


      subroutine truss_material_properties()
      implicit none
!     ================================================================
!                      material variables
!     ================================================================

!     ================================================================
      integer::eye(dimen,dimen)
      integer::i,j,k,l

!     ================================================================


eye = 0; do i = 1,dimen; eye(i,i) = 1;enddo
mu=ey_epoxy/(1+nu_epoxy)/2.0
lambda=ey_epoxy*nu_epoxy/(1+nu_epoxy)/(1-2.0*nu_epoxy)
ktense=0.0d0  ;ktense(1,1)=1;

      ctens=0.0d0
      do i=1,dimen;
          do j=1,dimen;
              do k=1,dimen;
                  do l=1,dimen;
                      ctens(i,j,k,l)  =ctens(i,j,k,l)+          &
                          element_direction_normal(i)*          &
                          element_direction_normal(j)*          &
                          element_direction_normal(k)*          &
                          element_direction_normal(l)* ey_epoxy
                  enddo;enddo;enddo;enddo;


k00=1.0
k01=0.0
lambda_01=1.0

k_elec_00=1.0
k_elec_01=0.0
lambda_elec_01=1.0


epz=0.0d0

      end subroutine truss_material_properties

! ===============================================================
!   Material properties
! ===============================================================
subroutine material_properties_afc()
implicit none
! ================================================================
!  material variables
! ================================================================

integer::eye(dimen,dimen)
integer::i,j,k,l ! ,m,n
real(iwp)::ey,nu

! ================================================================

eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
!
!3D formulation
ctens=0;
epz=0.0d0;
b_tilt=0.0d0;
ktense=0.0d0

! ================================================================
k00=k00_epx
k01=k01_epx
lambda_01=lambda_01_epx

k_elec_00=k_elec_00_epx
k_elec_01=k_elec_01_epx
lambda_elec_01=lambda_elec_01_epx
!
! ================================================================
!   E poxy
! ================================================================
mu=ey_epoxy/(1+nu_epoxy)/2.0
lambda=ey_epoxy*nu_epoxy/(1+nu_epoxy)/(1-2.0*nu_epoxy)
do i = 1,dimen;do j = 1,dimen;do k = 1,dimen;do l = 1,dimen;
ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+ &
   mu*( eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)  )
enddo;enddo;enddo;enddo;
epz=0.0d0;
b_tilt=0.0d0;

ktense=eps_epoxy*eye
partial_sigma_to_partial_elec_t=0
partial_sigma_to_partial_elec_p=0

! ================================================================
!   pzt
! ================================================================
if((noelem.ge.127).and.(noelem.le.336))then !it is pzt if 336>noelem>127

call remanent_polarization()
pr= curn_polarization_function(gauss_point_number,:)
direc_a=unit_vect(pr)

beta1   =-e31;
beta2   =-e33+2.0d0*e15+e31;
beta3   =-2.0d0*e15;
gamma1  =-eps_11/2.0d0;
gamma2  =(eps_11-eps_33)/2.0d0;


k00=k00_pzt
k01=k01_pzt
lambda_01=lambda_01_pzt

k_elec_00=k_elec_00_pzt
k_elec_01=k_elec_01_pzt
lambda_elec_01=lambda_elec_01_pzt

ey=ey_pzt;
nu=nu_pzt;

mu=ey/(1+nu)/2.0
lambda=ey*nu/(1+nu)/(1-2.0*nu)

do i = 1,dimen;do j = 1,dimen;do k = 1,dimen;do l = 1,dimen;
ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+ &
   mu*( eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)  )
enddo;enddo;enddo;enddo;

do i = 1,dimen; do k = 1,dimen; do l = 1,dimen;
epz(i,k,l)=( beta1*direc_a(i)*eye(k,l) &
            + beta2*direc_a(i)*direc_a(k)*direc_a(l) &
            + beta3*0.5d0*( eye(i,l)*direc_a(k) + eye(i,k)*direc_a(l) )  ) *norm_vect(pr)/pr_sat
enddo;enddo;enddo;

ktense(1,1)=eps_11
ktense(2,2)=eps_11 ;
ktense(3,3)=eps_33 ;


b_tilt=0.0d0
b_tilt(3,3,3,3)=  1.5e-5 ! * 2.0;

!b_tilt(3,3,2,2)=  1.8e-5  !*2.0;
!b_tilt(3,3,1,1)=  1.8e-5  !*2.0;


partial_sigma_to_partial_elec_t=0.0
partial_sigma_to_partial_elec_p=0.0

partial_sigma_to_partial_elec_t=epz
partial_sigma_to_partial_elec_p=epz



partial_sigma_to_partial_elec_t(3,3,3)=partial_sigma_to_partial_elec_t(3,3,3) &
        +b_tilt(3,3,3,3)*( abs(  electric_field(3) ) )

partial_sigma_to_partial_elec_p(3,3,3)=partial_sigma_to_partial_elec_t(3,3,3) &
        +b_tilt(3,3,3,3)*( abs(  electric_field_p(3) ) )
endif !if(noelem.ge.127.and.le.336)then !it is pzt if 336>noelem>127


! ================================================================
!   aluminum
! ================================================================
if((noelem.ge.113).and.(noelem.le.126))then !it is aluminum if 126>noelem>113

k00=k_00_aluminum
k01=k_01_aluminum
lambda_01=lambda_01_aluminum

ey=ey_aluminum
nu=nu_aluminum

mu=ey/(1+nu)/2.0
lambda=ey*nu/(1+nu)/(1-2.0*nu)

do i = 1,dimen;do j = 1,dimen;do k = 1,dimen;do l = 1,dimen;
ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+ &
   mu*( eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)  )
enddo;enddo;enddo;enddo;

epz=0.0d0;
b_tilt=0.0d0;
ktense=eye


partial_sigma_to_partial_elec_t=0
partial_sigma_to_partial_elec_p=0

endif ! (noelem.ge.113).and.(noelem.le.126))then !it is aluminum if 126>noelem>113

end subroutine material_properties_afc




end module material_behavior
