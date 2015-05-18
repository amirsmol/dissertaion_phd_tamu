module material_behavior
use linearsolvers
implicit none
! ================================================================
real(iwp),parameter::k00=1.0d0
real(iwp),parameter::k01=0.0d0
real(iwp),parameter::lambda_01=0.0

real(iwp),parameter::k_elec_00=1.0d0
real(iwp),parameter::k_elec_01=0.0d0
real(iwp),parameter::lambda_elec_01=0.0d0
! ================================================================
real(iwp),parameter::ec=1.0;
real(iwp),parameter::pr_sat=0.26;
real(iwp),parameter::eps_s_r=1.0e-10;
real(iwp),parameter::c=1.0;
real(iwp),parameter::eta=10e-8;
real(iwp),parameter::m=2.0;
real(iwp),parameter::poarization_speed=10.0
real(iwp),parameter::poarization_tolerance=1e-3
! ================================================================
real(iwp),parameter::c11 = 13.9e10;
real(iwp),parameter::c12 = 7.78e10;
real(iwp),parameter::c13 = 7.43e10;
real(iwp),parameter::c33 = 11.3e10;
real(iwp),parameter::c44 = 2.56e10;

real(iwp),parameter::e31 = -6.98 ! e-6;
real(iwp),parameter::e33 = 13.84 ! 2e-6;
real(iwp),parameter::e15 = 13.44 ! e-6;

real(iwp),parameter::eps_11 = 6.0e-9 ! ;
real(iwp),parameter::eps_33 = 5.47e-9 ;

logical::is_polarized

real(iwp)::a
real(iwp)::direc_a(dimen)


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

! =========================element center point
real(iwp) :: el_center_coord(dimen)
real(iwp) :: qtransform(dimen,dimen) ! this is the transformtion matrix
! =======================mechanical properties variables
real(iwp)::ctens(3,3,3,3)
! ===========================materials tensors
real(iwp),dimension(3,3,3):: epz
real(iwp),dimension(3,3)::ktense
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


!! vecotr is polarized is defined for all the gauss points
!! it is .true. if the material is polarized in the specific point
!! and it is false otherwise.
!! this is used to predict polarization switching behaviour of
!! material

logical,allocatable::vector_is_polarized(:)


real(iwp)::total_displacement(dimen)
real(iwp)::total_polarization(dimen)
real(iwp)::pr(dimen),prn(dimen)


! the output parameters
real(iwp),allocatable::elements_electric_field(:,:) !(nem,dimen)
real(iwp),allocatable::elements_electric_polar(:,:) !(nem,dimen)

! the invariants of electric field and strain



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
real(iwp)::d_electric_field(dimen),d_electric_field_p(dimen)
real(iwp)::strain_r(dimen,dimen),strain_elastic(dimen,dimen);
real(iwp)::electric_field(dimen),electric_field_p(dimen),electric_field_b(dimen);
! ================================================================
eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
! ================================================================

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
is_polarized=vector_is_polarized(gauss_point_number)
! ================================================================

call remanent_polarization(electric_field,d_electric_field,strain)

pr= curn_polarization_function(gauss_point_number,:)
prn=hist_polarization_function(gauss_point_number,:)
vector_is_polarized(gauss_point_number)=is_polarized

!pr=0.0d0
!pr(3)=pr_sat
!pr=prn

pr=prn

call direction_polarization(pr,direc_a)
call material_properties()
! ===================================remanent strain
strain_r=0.0d0

!pr= curn_polarization_function(gauss_point_number,:)
do i=1,dimen
do j=1,dimen
strain_r(i,j)=3.0d0*0.5d0*eps_s_r*(direc_a(i)*direc_a(j)-eye(i,j)/3.0d0)*norm_vect(pr)/pr_sat
end do;
end do;

!write(1,*)strain_r
!strain_r=0.0d0

strain_elastic=strain-strain_r
! ==============================================increment in strain
mechanical_strain_curentt(gauss_point_number,:,:)=strain_elastic
d_strain=mechanical_strain_curentt(gauss_point_number,:,:)-mechanical_strain_previus(gauss_point_number,:,:)
d_strain_p=mechanical_strain_previus(gauss_point_number,:,:)-mechanical_strain_beforep(gauss_point_number,:,:)
! ==============================the history variable

do i=1,dimen
do j=1,dimen
do k=1,dimen

sigma_el_hist_curentt(gauss_point_number,i,j,k)= &
sigma_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime)+ &
epz(k,i,j)*(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime)*d_electric_field_p(k)+d_electric_field(k) &
)


displ_el_hist_curentt(gauss_point_number,i,j,k)= &
displ_el_hist_previus(gauss_point_number,i,j,k)*exp(-lambda_elec_01*dtime)+ &
epz(k,i,j)*(k_elec_01*0.5d0)* &
(  &
exp(-lambda_elec_01*dtime)*d_strain_p(i,j)+d_strain(i,j) &
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


!> Here we difined a subroutine that gets the electric field, strain and the direction of polarization
! and outputs the energy invariants ii_inv
!

!call invariants_electric_entalpy(electric_field,direc_a,strain_elastic,pr)

sigma =  0.0d0 ;
der   =  0.0d0 ;

do i=1,dimen
do j=1,dimen
der(i)=der(i)+ktense(i,j)*electric_field(j)
do k=1,dimen
   der(k)  =  der(k)    +  k_elec_00*epz(k,i,j)*( strain_elastic(i,j) ) &
                        +  displ_el_hist_curentt(gauss_point_number,i,j,k)

sigma(i,j) =  sigma(i,j)  -  k_elec_00*epz(k,i,j)*  electric_field(k) &
                          -  sigma_el_hist_curentt(gauss_point_number,i,j,k)
do l=1,dimen
sigma(i,j)=sigma(i,j)+k00*ctens(i,j,k,l)*(  strain_elastic(k,l) )  &
         + mechanical_hist_curentt(gauss_point_number,i,j,k,l)
end do;
end do;
end do;
end do;


end subroutine stress_elect_displacement


!> Here we difined a subroutine that gets the electric field, strain and the direction of polarization
!< This subroutine define and producd the boundary conditions for afc unit cell !>
!! This subroutine define and producd the boundary conditions for afc unit cell
!! @param electric_field the electric field
!! @param strain_elastic the elastic strain
!! @param direc_a direction of polarization

!! @todo write the code




! ================================================================
!   polarization swithing function
!this will find the polarization with the given state of material
! ================================================================
subroutine remanent_polarization(electric_field,d_electric_field,strain)
implicit none
real(iwp)::d_electric_field(dimen)
real(iwp)::electric_field(3)
real(iwp)::delta_time
real(iwp)::strain(dimen,dimen)

real(iwp)::remanent_polarization_vector(dimen)
real(iwp)::remhisto_polarization_vector(dimen)
real(iwp)::electric_drive(dimen)
real(iwp)::elect_force
real(iwp)::er,e
real(iwp)::ps
real(iwp)::pr
real(iwp)::pn(dimen)
real(iwp)::del_p_del_er(dimen,dimen)
real(iwp)::del_er_del_e(dimen,dimen)
real(iwp)::del_e_del_t(dimen,dimen)

ps=pr_sat
!real(iwp)::elect_force

delta_time=time(1)
remanent_polarization_vector=0.0d0
remhisto_polarization_vector=hist_polarization_function(gauss_point_number,:)

electric_drive=electric_field-3*ec*direction_of_vector(d_electric_field)
elect_force=norm_vect(electric_field)-3.0d0*ec

remanent_polarization_vector=tanh( elect_force ) * &
direction_of_vector(electric_field) * pr_sat


!if(norm_vect(electric_field).gt.8*ec)then
!is_polarized=.true.
!endif

if( ( 1-norm_vect(remhisto_polarization_vector)/(pr_sat) ).lt. poarization_tolerance )then
    is_polarized=.true.
endif


if(.not.is_polarized)then
remanent_polarization_vector=0.5d0*pr_sat*(tanh(elect_force)+1)*direction_of_vector(electric_field)
endif



if(is_polarized)then
if(norm_vect(electric_field).lt.8*ec)then

    if( dot_product(electric_field,d_electric_field).lt.0 )then
        remanent_polarization_vector=hist_polarization_function(gauss_point_number,:)
    endif

endif
endif



curn_polarization_function(gauss_point_number,:)=remanent_polarization_vector
curn_electric_field(gauss_point_number,:)=electric_field

end subroutine remanent_polarization


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

allocate( vector_is_polarized (ngauss) )

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

vector_is_polarized=.false.



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
!real(iwp)::ey,nu

! ================================================================
! ================================================================
! ================================================================
lambda  =c12;
mu      =(c11-c12)/2.0d0;
alpha1  =2*c44+c12-c11;
alpha2  =(c11+c33)/2.0d0-2.0d0*c44-c
alpha3  =c13-c12;
beta1   =-e31;
beta2   =-e33+2.0d0*e15+e31;
beta3   =-2.0d0*e15;
gamma1  =-eps_11/2.0d0;
gamma2  =(eps_11-eps_33)/2.0d0;

eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo
!

direc_a=direction_of_vector(pr)

!3D formulation
ctens=0;
epz=0.0d0;
ktense=0.0d0


do i = 1,dimen;do j = 1,dimen;do k = 1,dimen;do l = 1,dimen;
ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+ &
   mu*( eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)  )
enddo;enddo;enddo;enddo;

do i = 1,dimen; do k = 1,dimen; do l = 1,dimen;
epz(i,k,l)=(- beta1*direc_a(i)*eye(k,l) &
            - beta2*direc_a(i)*direc_a(k)*direc_a(l) &
            - beta3*0.5d0*( eye(i,l)*direc_a(k) + eye(i,k)*direc_a(l) )  ) *norm_vect(pr)/pr_sat
enddo;enddo;enddo;

do i = 1,dimen;do j = 1,dimen;
ktense(i,j)=-(2*gamma2*direc_a(i)*direc_a(j)+gamma1*eye(i,j))
enddo;enddo;

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




end module material_behavior
