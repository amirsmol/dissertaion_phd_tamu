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
 !!  \brief     Finite Element Model for Active Fiber Composite
 !!  \details   This class is used to demonstrate a number of section commands.
 !!  \author    Amir Sohrabi Mollayousef
 !!  \version   0.1a
 !!  \date      2011 - 2014
 !!  \pre       Pay attention to the mesh file.
 !!  \bug       Not all memory is freed when deleting an object of this class.
 !!  \warning   Improper use can crash your application
 !!  \copyright This is free software; you can use it, redistribute
 !! it, and/or modify it under the terms of the GNU Lesser General

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
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
      length=1.0d0;neldirectional=[2,2,2]
      call linear_c3d8_3d_fem_geometry(neldirectional,length)
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
      call crawley_boundary(length)

      write(*,*)'length',length
!     ===============================================================
!                        time increment starts here
!     ===============================================================
      vspvt=vspv
      vssvt=vssv
!>
!!    reading iteration and convergence critria's
!<
      tolerance=1.0e-3
      max_iteration=20;
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      dtime=0.01;      freq=1.0d0;
      max_time_numb= int(5.0e0/dtime)
!     ===============================================================
!    write(*,*)'number of time steps',max_time_numb

     do time_step_number=0,max_time_numb
     call cpu_time(timer_begin)
         time(1)=dtime;time(2)=time_step_number*dtime
         loadfactor=-sin(2*3.14515*freq*time(2))*5.0*ec
         vspv=loadfactor*vspvt
         glu(bnd_no_pr_vec)=0.0*vspv;
         glu=0.0d0
!!     ===============================================================
!!                 nonlinear solution iteration starts here
!!     ===============================================================
     normalforce=norm_vect(vspv)+norm_vect(vssv)
     converged=.false.
     do iteration_number=0,max_iteration
     write(22,*)'global iteration_number=',iteration_number
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


    write(24,951)time_step_number,time(2),loadfactor,curn_polarization_function(3,:)

    call update_history()
    glb=glp;glp=glu

    call result_printer(iter,glu,loadfactor)
    call paraview_3d_vtu_xml_writer(glu)
    enddo !itime=0,ntime


      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()
      close (csv)

951   format(5x,i5,' , ',10 (e14.5,' , ') )

end program fem_piezeo_3d_visco_polarization_switching
