!!> ============================================================================
!! Name        : Elec_Mech_FEM.f90
!! Author      : Amir Sohrabi
!! Version     : V0.00
!! Copyright   : Electro Mechanical Finite Element Model Copyright (C) 2015  Amir Sohrabi <amirsmol@gmail.com>
!!               This code is distributed under the GNU LGPL license.
!! Description :
!!> ============================================================================


!>  \mainpage tetra_hedral_space_filler_truss
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

program tetra_hedral_linear_truss
use fem_functions_and_parameters
use fem_geometry
use material_behavior
use fem_libs

implicit none
!!     ================================================================
!!                      solution variables
!!     ================================================================
      real(iwp),allocatable::glk(:,:) !!> the global stiffness matrix solution
      real(iwp),allocatable::glu(:) !!>   the global solution vector
      real(iwp),allocatable::glp(:) !!>   the previus degrees of freedom vector
      real(iwp),allocatable::glb(:) !!>  the before previus degrees of freedom
!     ================================================================
!                      solver variables
!     ================================================================
      real(iwp),allocatable::glq(:) !!>internal force vector
      real(iwp),allocatable::glt(:) !!>external fource vector
      real(iwp),allocatable::glr(:) !!>total residual vector
!     ================================================================
!                      time variables
!     ================================================================
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
      integer::num_tetra_units
!     ============================iteration variables
      real(iwp)::  error ! ,eps!the error tolerance and error
      logical:: converged

      integer::iteration_number
      integer::time_step_number

      integer::max_iteration
      integer::max_time_numb

      real(iwp)::tolerance
      real(iwp)::normalforce
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
       length=1
       num_tetra_units=150

      call truss_actua_3d_connectivity(length,num_tetra_units)
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
      allocate(glk(neq,neq),glu(neq),glq(neq),glt(neq),glr(neq),glp(neq),glb(neq))
      glu=0.0d0;glt=0.0d0;glr=0.0d0;glp=0.0d0;glb=0.0d0;
!!     ===============================================================

!      call linear_truss_bending_boundary()
      call linear_truss_bending_boundary_2()
!     ===============================================================
!                        time increment starts here
!     ===============================================================
      vspvt=vspv
      vssvt=vssv
!>
!!    reading iteration and convergence critria's
!<
      tolerance=1.0e-4
      max_iteration=50;
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      dtime=0.1;      ! freq=1.0d0;
      max_time_numb= int(20.0e0/dtime)
!     ===============================================================
do time_step_number=0, max_time_numb
     call cpu_time(timer_begin)
         time(1)=dtime;time(2)=time_step_number*dtime

         vspv=time(2)*vspvt
         glu(bnd_no_pr_vec)=0.0d0;

!         glu=0.0d0;
!!     ===============================================================
!!                 nonlinear solution iteration starts here
!!     ===============================================================
     normalforce=norm_vect(vspv)+norm_vect(vssv)+dtime

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
write(*,*)'check'
write(*,*)'check',bnd_no_se_vec
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
      call lapack_gesv(glk,glr)
!!     ===============================================================
!!                        updating the solution
!!     ===============================================================
      glu=glu+glr

      write(out,*)'iteration_number=', iteration_number, 'time=',time,'error=',error,'normalforce=',normalforce
      write(*,*)'iteration_number=', iteration_number, 'time=',time,'error=',error,'normalforce=',normalforce

      if (error.le.tolerance*(normalforce))then;
          converged=.true.;exit;
      endif

     enddo ! time_step_number=0,max_time_numb
!     ======================updatng coordinates  updated lagrangian
    if(.not.converged)then
     write(*,*)"iteration did not converge"
     stop
    endif


    glb=glp;glp=glu
   call truss_paraview_3d_vtu_xml_writer(glu)
    write(*,*)'max_time_iteration',max_time_numb
    enddo !itime=0,ntime

      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()
      close (csv)

end program tetra_hedral_linear_truss




