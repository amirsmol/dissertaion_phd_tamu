module material_behavior
use linearsolvers


!> The time scaling variables
real(iwp),parameter::gamma_e=1.75e-6;
real(iwp)::a_electric_shift 

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
real(iwp),parameter::k_elec_01_pzt=-0.3d0
real(iwp),parameter::lambda_elec_01_pzt=1.5

! =============================epoxy
real(iwp),parameter::k00_epx=1.0d0
real(iwp),parameter::k01_epx=0.6d0
real(iwp),parameter::lambda_01_epx=0.8

real(iwp),parameter::k_elec_00_epx=1.0d0
real(iwp),parameter::k_elec_01_epx=0.0d0
real(iwp),parameter::lambda_elec_01_epx=1.0
! ================================================================
real(iwp),parameter::ec=1;
real(iwp),parameter::pr_sat=26e-2;
real(iwp),parameter::eps_s_r=1e-3;
real(iwp),parameter::c=1.0;
real(iwp),parameter::eta=10e-2;
real(iwp),parameter::m=2;
! ================================================================ pzt


real(iwp),parameter::e33 = -14.52*1.3d0;
real(iwp),parameter::e15 = -7.56*1.3d0;
real(iwp),parameter::e31 = 6.15*1.3d0;

!real(iwp),parameter::e33 = -20.92*1.3d0;
!real(iwp),parameter::e31 =  6.15*1.3d0;
!real(iwp),parameter::e15 = -45.74*1.3d0;

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
real(iwp)::time(2),dtime,dtime_scaled
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



real(iwp)::electric_field(dimen),electric_field_p(dimen),electric_field_b(dimen);

contains

!< This subroutine defines the time shift and reduced time factor!>
!! @param electric_field this is the electric field vector
!! @param a_electric_shift this is the electric shift vector
!! @todo noting for now
subroutine get_time_scale_factor()
        implicit none
        real(iwp)::electric_norm


        electric_norm=norm_vect(electric_field)
        ! a_electric_shift=exp(- gamma_e*(electric_norm-electric_field_o) /electric_field_o)
        a_electric_shift=exp(- gamma_e * electric_norm)        

end subroutine get_time_scale_factor


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
! ================================================================
curn_polarization_function=0.0d0
!call remanent_polarization(d_electric_field,electric_field)
!pr= curn_polarization_function(gauss_point_number,:)
pr=0.0d0
pr(3)=pr_sat
call direction_polarization(pr,direc_a)
call material_properties()


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
sigma_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime_scaled)+ &
(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime_scaled)*d_electric_field_p(k)*partial_sigma_to_partial_elec_p(k,i,j) &
                          +d_electric_field(k)  *partial_sigma_to_partial_elec_t(k,i,j) &
)

displ_el_hist_curentt(gauss_point_number,i,j,k)= &
displ_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime_scaled)+ &
(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime_scaled)*d_strain_p(i,j)*partial_sigma_to_partial_elec_p(k,i,j)+ &
                           d_strain(i,j)  *partial_sigma_to_partial_elec_t(k,i,j) &
)


do l=1,dimen
mechanical_hist_curentt(gauss_point_number,i,j,k,l)= &
mechanical_hist_previus(gauss_point_number,i,j,k,l)*exp(-lambda_01*dtime_scaled)+ &
ctens(i,j,k,l)*(k01*0.5d0)* &
(  &
exp(-lambda_01*dtime_scaled)*d_strain_p(k,l)+d_strain(k,l) &
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
subroutine remanent_polarization(d_electric_field,electric_field)
implicit none
real(iwp)::d_electric_field(3)
real(iwp)::electric_field(3)
real(iwp)::delta_time

real(iwp)::remanent_polarization_vector_1(dimen)
real(iwp)::electric_drive(dimen)
real(iwp)::elect_force

delta_time=time(1)

remanent_polarization_vector_1=0.0d0

electric_drive=electric_field-3*ec*direction_of_vector(d_electric_field)
elect_force=norm_vect(electric_field)-3.0d0*ec

remanent_polarization_vector_1=tanh( elect_force ) * &
direction_of_vector(electric_field) * pr_sat


if(norm_vect(electric_field).gt.8*ec)then
is_polarized=.true.
endif

if(.not.is_polarized)then
remanent_polarization_vector_1=0.5d0*pr_sat*(tanh(elect_force)+1)*direction_of_vector(electric_field)
endif



if(is_polarized)then
if(norm_vect(electric_field).lt.8*ec)then

if( dot_product(electric_field,d_electric_field).lt.0 )then
remanent_polarization_vector_1=hist_polarization_function(gauss_point_number,:)
endif

endif
endif



!if(norm_vect(electric_field).lt.ec)then
!remanent_polarization_vector_1=hist_polarization_function(gauss_point_number,:)
!endif



curn_polarization_function(gauss_point_number,:)=remanent_polarization_vector_1
curn_electric_field(gauss_point_number,:)=electric_field

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

end subroutine  update_history

! ===============================================================
!   Material properties
! ===============================================================
subroutine material_properties()
implicit none
! ================================================================
!  material variables
! ================================================================

! ================================================================
integer::eye(dimen,dimen)
integer::i,j,k,l ! ,m,n
real(iwp)::ey,nu
! ================================================================
! ================================================================
! ================================================================

eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
!

direc_a=direction_of_vector(pr)

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
dtime_scaled=dtime
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

dtime_scaled=dtime
call get_time_scale_factor()
dtime_scaled=dtime/a_electric_shift

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
! b_tilt(3,3,3,3)=  1.5e-5 ! * 2.0;

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
dtime_scaled=dtime

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

end subroutine material_properties

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


      subroutine truss_shape_change_stress(el_center_coord,ur,sigma_shape)
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
      epsilon=time(2)*1.0d0

e=0.0d0
e(1,1)=x3*time(2)
!e(2,2)=kappa
!e(3,3)=kappa

lagrange_strain_global=e

sigma_shape=0.0d0
do i=1,dimen;
do j=1,dimen;
do k=1,dimen;
do l=1,dimen;
sigma_shape(i,j)=sigma_shape(i,j)+ctens(i,j,k,l)*lagrange_strain_global(k,l)

enddo;enddo;enddo;enddo;

!each_truss_strain(noelem)=sigma_shape(1,1)

end subroutine truss_shape_change_stress


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
                          element_direction_normal(l)* 1.0d0
                  enddo;enddo;enddo;enddo;

epz=0.0d0

      end subroutine truss_material_properties




end module material_behavior
