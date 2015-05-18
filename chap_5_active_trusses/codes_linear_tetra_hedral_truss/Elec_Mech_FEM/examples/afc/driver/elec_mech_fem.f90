!!> ============================================================================
!! Name        : Elec_Mech_FEM.f90
!! Author      : Amir Sohrabi
!! Version     : V0.00
!! Copyright   : Electro Mechanical Finite Element Model Copyright (C) 2015  Amir Sohrabi <amirsmol@gmail.com>
!!               This code is distributed under the GNU LGPL license.
!! Description :
!!> ============================================================================


!>  \mainpage Finite Element Model for Active Fiber Composite
 !!
 !! \section intro_sec Introduction
 !!
 !! This is the introduction.
 !!
 !! \section install_sec Installation
 !!
 !! \subsection step1 Step 1: Opening the box
 !!
 !!
 !!
!>
 !!  \brief     Electro Mechanical Finite Element Model
 !!  \details   This is the driver file
 !!  \author    Amir Sohrabi Mollayousef
 !!  \version   0.1a
 !!  \date      2011 - 2014
 !!  \pre       Pay attention to the mesh file.
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

program fem_piezeo_3d_visco_polarization_switching
use linearsolvers
use fem_geometry
use material_behavior
use fem_libs
implicit none
!     ================================================================
!                      solution variables
!     ================================================================
      real(iwp),allocatable::glk(:,:)
      real(iwp),allocatable::glu(:) !the global solution
      real(iwp),allocatable::glp(:) !the previus degrees of freedom vector
      real(iwp),allocatable::glb(:) !the before previus degrees of freedom
!     ================================================================
!                           the variables
!     ================================================================
      real(iwp),allocatable::glk_constrain(:,:) !< the global constrainted constant matrix
      real(iwp),allocatable::glr_constrain(:) !< the global constraint solution vector
!      real(iwp),allocatable::aux_glk_constrain(:,:) !< the global constraint solution vector
      real(iwp),allocatable::aux_glk_constrain_t(:,:) !< the global constraint solution vector
      real(iwp),allocatable::transpose_constrain(:,:) !< the global constraint solution vector
!      real(iwp),allocatable::constraint_real(:,:)
!     ================================================================
!                      solver variables
!     ================================================================
      real(iwp),allocatable::glq(:) !internal force vector
      real(iwp),allocatable::glt(:) !external fource vector
      real(iwp),allocatable::glr(:) !total residual vector
      real(iwp),allocatable::glu_ref(:) !displacement in refrence coordinate
!     ============================iteration variables
      real(iwp)::  error ! ,eps!the error tolerance and error

!      integer:: itmax !maxumim number of iteration
      logical:: converged

      integer::iteration_number
      integer::time_step_number

      integer::max_iteration
      integer::max_time_numb

      real(iwp)::tolerance
      real(iwp)::normalforce
!     ================================================================
!                      time variables
!     ================================================================
      real(iwp)::loadfactor    ! time
      real(iwp)::freq
      integer :: clck_counts_beg, clck_counts_end, clck_rate
      real(iwp) :: beg_cpu_time, end_cpu_time
      real(iwp):: timer_begin, timer_end
!     ================================================================
!      integer::igauss
      integer::i,j
      real(iwp),parameter::myone=1.0d0
!     =======================trivial meshing arrays the meshing seed parameters
      integer :: neldirectional(dimen)
      real(iwp) :: length(dimen)
            real(iwp):: load_factors(2)
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
!       length=[1.4000001800000002E-004,1.7499999400000000E-004,7.4999988999999980E-004]
!       length=1
!       neldirectional=[1,1,1]
!       call linear_c3d8_3d_fem_geometry(neldirectional,length)
      call afc_mesh_file_reader()
!     ===============================================================
!     define the solution parameters
!     ===============================================================
      nnm=size(coords,dim=1)
      nem=size(nod,dim=1)
      neq  = nnm*ndf;     nn=npe*ndf;
      ipdf = iel+1;    ngauss=(ipdf**dimen)*nem
      write(*,*)'Number to equations to solve=',neq
      write(*,*)'Number to nodes=',nnm
!     ===============================================================
!     reading the boundary conditions
!     ===============================================================
      call form_history(ngauss,nem)
      allocate(glk(neq,neq),glu(neq),glq(neq),glt(neq),glr(neq),glp(neq),glb(neq),glu_ref(neq))
      glu=0.0d0;glt=0.0d0;glr=0.0d0;glp=0.0d0;glb=0.0d0;glu_ref=0.0d0
!!     ===============================================================
      call afc_boundary_full_electrode(length)
!      call crawley_boundary(length)

      write(*,*)'length',length

      call afc_form_constraint()
      allocate( glk_constrain( size(constraint,dim=2), size(constraint,dim=2) ) )
      allocate( glr_constrain( size(constraint,dim=2) ) )
!
      allocate(transpose_constrain(   size(constraint,dim=2), size(constraint,dim=1)  )      )
      allocate(aux_glk_constrain_t(   size(transpose_constrain,dim=1), size(glk,dim=2)  )      )
!     write(*,*)'check======================='
      transpose_constrain=transpose(constraint)
!     ===============================================================
!                        time increment starts here
!     ===============================================================
      vspvt=vspv
      vssvt=vssv
!>
!!    reading iteration and convergence critria's
!<
      tolerance=1.0e-3
      max_iteration=30;
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      dtime=0.01;      freq=1.0d0;
      max_time_numb= int(2.0e0/dtime)
!     ===============================================================
    write(*,*)'number of time steps',max_time_numb
    load_factors(1)=83.0
    load_factors(2)=375.0
     do i_calibration=1,2
     do time_step_number=0, max_time_numb
     call cpu_time(timer_begin)


         time(1)=dtime;time(2)=time_step_number*dtime
                  loadfactor=sin(2*3.14515*freq*time(2))*load_factors(i_calibration)
!        loadfactor=-sin(2*3.14515*freq*time(2))*375
!        loadfactor=375 !175.0d0
!        loadfactor=time(2)
         vspv=loadfactor*vspvt
         glu(bnd_no_pr_vec)=0.0d0;
!        glu=0.0d0;
!!     ===============================================================
!!                 nonlinear solution iteration starts here
!!     ===============================================================
     normalforce=norm_vect(vspv)+norm_vect(vssv)
     converged=.false.
     do iteration_number=0,max_iteration
!!     ===============================================================
!!                        forming the global matrices
!!     ===============================================================
      glk=0.0d0;glq=0.0d0;!
      call glbmatrcs(glu,glk,glq,glp,glb)
!!     ===============================================================
!!                             newton raphson sprocedure
!!     ===============================================================
      glt(bnd_no_se_vec) = vssv ;
      glr=glt-glq ;
!!     ===============================================================
!!                        solving the linear system of equations
!!     ===============================================================
      call symmetric_primary_bounday(glk,glr)
      error=norm_vect(glr)
      vspv=0.0d0
!     ===============================================================
!                        solving constraint
!     ===============================================================
!     call cpu_time(timer_begin)

      aux_glk_constrain_t=lapack_matmul(transpose_constrain,glk)
      glk_constrain=lapack_matmul (aux_glk_constrain_t, constraint)
      glr_constrain=matmul(transpose_constrain,glr)
!     call cpu_time(timer_end); write(*,*)'time to apply constrainst',timer_end- timer_begin
!     call cpu_time(timer_begin)
!     write(*,*)'size of the coefficient matrix=',size(glk_constrain,dim=1)

     call lapack_gesv(glk_constrain,glr_constrain)

!     call cpu_time(timer_end); write(*,*)'time to solve linear system',timer_end- timer_begin
     glr=matmul (constraint,glr_constrain)
!!     ===============================================================
!!                        updating the solution
!!     ===============================================================
      glu=glu+glr

      write(out,*)'iteration_number=', iteration_number, 'time=',time,'error=',error,'normalforce=',normalforce
      write(*,*)'iteration_number=', iteration_number, 'time=',time,'error=',error,'normalforce=',normalforce
      write(*,*)'loadfactor amplitude=',load_factors(i_calibration)
      write(*,*)

      if (error.le.tolerance*(normalforce))then;
          converged=.true.;exit;
      endif

     enddo ! time_step_number=0,max_time_numb
!     ======================updatng coordinates  updated lagrangian
    if(.not.converged)then
     write(*,*)"iteration did not converge"
     stop
    endif

if(updated_lagrangian)then
glu_ref=glu+glu_ref
do j=1,dimen;forall(i=1:nnm) coords(i,j)=coords(i,j)+glu((i-1)*ndf+j);enddo
endif

    call update_history()
    glb=glp;glp=glu

    call result_printer(iter,glu,loadfactor)
!    call paraview_3d_vtu_xml_writer(glu)
!    call paraview_3d_vtu_xml_writer_vector(glu,elements_electric_field,elements_electric_polar)
!    write(*,*)'max_time_iteration',max_time_numb
    enddo !itime=0,ntime

    call clear_history()

    glb=0.0d0
    glp=0.0d0
    glp=0.0d0
    glu=0.0d0

    end do ! i_calibration=1,2


      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()
      close (csv)

end program fem_piezeo_3d_visco_polarization_switching
