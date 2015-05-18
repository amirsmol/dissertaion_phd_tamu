module fem_geometry
 use linearsolvers
    implicit none
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

     open (msh,file='src/3d_bimorh_quadratic_mesh.msh')
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
        character(len=20) :: fmt,filename ! format descriptor
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


        filename='geom_test'//trim(x1)//'.vtu'
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
character(len=20) :: fmt,filename ! format descriptor
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


filename='outpar'//trim(x1)//'.vtu'
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


open (msh,file='src/3d_afc_mesh.msh')
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
