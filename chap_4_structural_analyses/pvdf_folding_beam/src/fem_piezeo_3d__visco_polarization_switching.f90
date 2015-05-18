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
!!>  Reading the mesh for the solution
!     ===============================================================
!      neldirectional=4
!      length=1.0
      call beam_mesh_file_reader()
!      call linear_c3d8_3d_fem_geometry(neldirectional,length)
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
      call beam_bimorph_boundary(length)

!      call crawley_boundary(length)
!      call wedge_boundary(length)
!      write(*,*)'length',length

!      call afc_form_constraint()
      call beam_bimorph_form_constraint()

      allocate( glk_constrain( size(constraint,dim=2), size(constraint,dim=2) ) )
      allocate( glr_constrain( size(constraint,dim=2) ) )
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
      max_iteration=50;
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      max_time_numb=5
      dtime=455.0e3/max_time_numb;  freq=1.0d0;
!      max_time_numb= int(500.0e3/dtime)
!     ===============================================================
    write(*,*)'number of time steps',max_time_numb
     do time_step_number=0, max_time_numb
     call cpu_time(timer_begin)
         time(1)=dtime;time(2)=time_step_number*dtime
!         loadfactor=-sin(2*3.14515*freq*time(2))*10*ec
         loadfactor=time(2) !175.0d0
!         loadfactor=time(2)
         vspv=loadfactor*vspvt
         glu(bnd_no_pr_vec)=0.0d0;
         glu=0.0d0;
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
     call cpu_time(timer_begin)
      aux_glk_constrain_t=lapack_matmul(transpose_constrain,glk)
      glk_constrain=lapack_matmul (aux_glk_constrain_t, constraint)
      glr_constrain=matmul(transpose_constrain,glr)

!     call cpu_time(timer_end); write(*,*)'time to apply constrainst',timer_end- timer_begin
!     call cpu_time(timer_begin)
!     write(*,*)'size of the coefficient matrix=',size(glk_constrain,dim=1)

     call lapack_gesv(glk_constrain,glr_constrain)

!    call cpu_time(timer_end); write(*,*)'time to solve linear system',timer_end- timer_begin
     glr=matmul (constraint,glr_constrain)
!!     ===============================================================
!!                        solving non constraint
!!     ===============================================================
!    call gaussian_elimination_solver(glk,glr)
!    call parafem_cg_solver(glk,glr)
!     call cpu_time(timer_end); write(*,*)'time to prepare the system',timer_end- timer_begin
!     call cpu_time(timer_begin)
!     call lapack_gesv(glk,glr)
!     call cpu_time(timer_end); write(*,*)'time to solve the system',timer_end- timer_begin
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

if(updated_lagrangian)then
glu_ref=glu+glu_ref
do j=1,dimen;forall(i=1:nnm) coords(i,j)=coords(i,j)+glu((i-1)*ndf+j);enddo
endif

    call update_history()
    glb=glp;glp=glu

    call result_printer(iter,glu,loadfactor)
!    call paraview_3d_vtu_xml_writer(glu)
    call paraview_3d_vtu_xml_writer_vector(glu,elements_electric_field,elements_electric_polar)
!    write(*,*)'max_time_iteration',max_time_numb
    enddo !itime=0,ntime


      call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()
      close (csv)

end program fem_piezeo_3d_visco_polarization_switching
