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

program tetra_hedral_space_filler_truss
use fem_functions_and_parameters
use fem_geometry
implicit none
!!     ================================================================
!!                      solution variables
!!     ================================================================
      real(iwp),allocatable::glu(:) !the global solution
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
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
       length=1
       num_tetra_units=10
       neldirectional=0
       neldirectional(1)=10
       neldirectional(2)=2
       neldirectional(3)=1

      call truss_tetrahedral_space_filler(neldirectional,length,num_tetra_units)

call show_matrix(coords,"coords")
call show_matrix_int(nod,"nod")

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
      allocate(glu(neq))
      glu=0.0d0;
!!     ===============================================================
      call truss_paraview_3d_vtu_xml_writer(glu)

      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()
      close (csv)

end program tetra_hedral_space_filler_truss




