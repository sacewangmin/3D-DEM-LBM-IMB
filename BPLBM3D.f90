!*****************************************************************************************************************************************
!
!  PROGRAM: BPLBM3D_NEW
!
!  DEVELOPED BY:  Min Wang
!               Zienkiewicz Centre for Computational Engineering
!	              College of Engineering
!	              Swansea University, UK
!               June 2015   
!  
!  PURPOSE:  This program is aimed at modelling 1) the mechanical behaviour of geomaterials (continuous, discontinuous and granular media),
!            2) fluid flow using Lattice Boltzmann Equation and 3) the fluid-solid interactions.
!            It includes DEM, BPM, LBM and Their coupling schemes.
!
!
!  Modified by: 
!            Anders Bahrami            
!            1) Velocity and Periodic Boundary conditions, May 2026
!
!            Ryan Nguyen
!            1) MRT D3Q15 model, Jun 2026
!            2) Improved paraview output, Jun 2026
!            3) Body force subroutine added, Jun 2026
! 
!            Min Wang
!            1) modified the B to be nonlinear function, Jun 2026
!            2) modified body force subroutine for LBM, July 2026
!
!*****************************************************************************************************************************************
!
!****************************************************Instruction for DELBM****************************************************************
!
!Program BPLBM3D
!  Subroutine read_parameters
!  initialisation and allocate space
!----------------------------------------------
!  Subroutine boundary_boundingbox
!  Subroutine contact_search_first
!  Subroutine install_bond
!----------------------------------------------
!  Subroutine boundary_profile
!  Subroutine init_density
!  Subroutine generate_stationary_particles
!  Subroutine generate_moving_particles
!  Subroutine generate_fluid_boundary_nodes
!-----------------MAIN LOOP--------------------
!  Subroutine contact_main
!  Subroutine update_moving_particle_position
!  Subroutine update_moving_particle_nodes
!  Subroutine update_fluid_boundary_nodes
!  Subroutine relaxation
!  Subroutine bounceback
!  Subroutine propagate
!  Subroutine boundary_treatment_regular
!----------------------------------------------
!  Subroutine output_result  
!End program
!
!????incompleted functions:
! 1--friction for DEM
! 2--IBM
!*****************************************************************************************************************************************
!    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Definition of variables and parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------system data---------------------------------------------------------------------
    Module system
    implicit none
!.....number of vertexes, number of edges/walls/surfaces
      integer nv,ne
!.....number of stationary and moving edges
      integer ne0,ne1
!.....ith timestep and maximum number of timestep, screen output frequency, file output frequency
      integer istep,msteps, screen_steps, outsteps
!.....flag for the first loop
      logical first,output
!.....flag for lbm,lubrication force calculation , uniaxial test
      logical LBM,lubrication,uniaxial
!.....flags for Immersed Boundary method and Immersed Moving Boundary
      logical IBM,IMB
      logical fix_packing      
!===== for coarse-grid LBM-DEM coupled by empirical equations
      logical empirical
!=====
!.....flags for periodic boundary (only used if bctype=3)
      logical periodic_flag
!.....flag for MRT
      logical use_MRT
!.....Relaxation values for MRT
      real(8) S1, S2, S4, S14
!.....simulation time    time:real time used, factor:time factor
      real(8) start_time,end_time,time,factor
!.....boundary vertexes coordinate
	  real(8),allocatable::vertex(:,:)
!.....maximum coordinate used for LBM domain
      real(8) xmax,ymax,zmax
!.....i,j,k are used for loop;   is:start point;  ie:end point ip:ith particle
	  integer i,j,k,is,ie,ip,ix,iy,iz
!.....parameter used in pointers err!=0:out of memory;  flag for pointer code check
      integer err,pointer_check
!     implementation mthods for periodic
      integer implement
!.....used for dealing with moving wall nodes  (integeration)
      real(8),allocatable::vertex_mw(:,:)
!.....used for dealing with stationary nodes under conditions of moving wall-fluid interaction 
      real(8),allocatable::vertex_sw(:,:)
      real(8), parameter::PI=3.14159265385 
!
!....................................................Definition of type wall property.....................................................
      type wall_prop
!.......wall stiffness,contact force
        real(8) kn,FC(3)
!.......1:start point;   2:end point;  
	    integer point(2)
!.......flag   0:stationary   1:moving;   rate-dependent flag
        integer flag,rate
!.......the velocity of each wall
!.......coor(:):1 x of start point 2 y of start point  3 z of start point  4 x of end point 5 y of end point 6 z of end point
        real(8) vel(3),coor(6)
!.......number of each moving wall nodes, number of moving wall_fluid nodes corresponding to each wall
        integer nmwn,nmwfn
!.......wfluid:flag for interaction between specific moving wall and fluid boundary nodes    
        integer wfluid
!????????????????????????????????????????????????????Maybe useful in the future
!.......list of each moving wall nodes coordinate, list of moving wall_fluid nodes corresponding to each wall coordinate
        integer, allocatable :: lmwnx(:),lmwfnx(:),lmwny(:),lmwfny(:)              
      end type      
!
!.................................................Definition of type pointers.............................................................
!.....this one is for wall and particle list
      type datalink
        integer num
        type(datalink),pointer::next
      end type
!.....this one is for wall nodes,particle nodes list
      type datalink2
        integer coord(3)
        type(datalink2),pointer::next
      end type
!
      type(wall_prop),allocatable::wall(:)
!.....p1,p2 are temporary pointers,  head_le0,head_le1 are heads for stationary wall and moving wall lists
      type(datalink),pointer::p1,p2,head_le0,head_le1
!.....head_lps,head_lpm are heads for stationary and moving particles lists
      type(datalink),pointer::head_lps,head_lpm
!.....coordinate list of stationary particle nodes,moving particle interior nodes,moving boundary nodes and fluid boundary nodes
!.....fixed wall boundary nodes,non-fixed wall boundary nodes, p3 is temporary pointer
      type(datalink2),pointer::head_lns,p3
!
    End module system
!
!----------------------------------------------------------solid data---------------------------------------------------------------------
    Module solid
    implicit none
!.....number of particles, number of stationary and moving particles
      integer np,ns,nm      
!.....number of maximum contact, contact number for each particle ,  it used for loop structure    
      integer mcnts,ntar,it
!.....solid densities
      real(8) ds
!.....gravational acceleration,total maximum displacement of particles in DEM contact,reduced_gravity
      real(8) gacce,sum_disp,gacce1
!.....damp_option
      integer damp_option
!.....damp_coef and timestep for DEM
      real(8) damp_coef,dtime
!.....number of subcycles for dem 
      integer nsub,isub
!.....number of particle stationary,interior and boundary nodes
      integer MS,MI,MB
!.....particle IP list of stationary and moving particle nodes)
!      integer,allocatable::lps(:),lpm(:)
!.....number of walls
      integer nwall
!.....Boundary type for DEM, 0:stationary wall; 1:periodic 2:special wall 3:used for tri-axial moving mall 4:others
      integer boundary_dem
!.....Particle contact list :List of total number of contacted bodies,  List of pointer to contacted bodies
	  integer,allocatable::laccoc(:),laccop(:)
!.....laccot:List of targetor numbers 
	  integer,allocatable::laccot(:)
!.....iccor: coordination of bounding boxes  xs: X-axis start;  xe: X end ys ye zs ze    
	  integer,allocatable::icoor(:,:)   
!.....buff for particle contact, minimum radius, maximum particle radius, maximum velocity,
!.....distance between two particles centres, gap between two particles, critical gap for bond installation
      real(8)  buff,minrad,maxrad,vmax,distance,gap,b_install
!.....bonded particle contact list :List of total number of contacted bodies for each particle,  
!.....List of pointer to contacted bodies for each particle,  List of IP of targetor particles
	  integer,allocatable::b_laccoc(:),b_laccot(:,:)
!.....initial gap between particles used for sampling        
      real(8),allocatable::gap_ini(:,:)
!.....output parameters for shperes
 !     real(8),allocatable::sphere_mesh(:,:)
 !     integer,allocatable::sphere_elem(:,:)
!.....b_exist: check if the current contact belongs to bond contact;   b_ip: pointer for current bond particle
      integer b_exist,b_ip
!.....True:use bond model;  False: without bond model   
      logical bond,sam_prep2
!.....computational parameters for bond model:stiffness, normal force, critical normal force,critical shear force for bond
      real(8)  cb_kn,cb_nf,cb_cnf,b_ftmax
!.....physical parameters for bond model:stiffness, normal force, critical normal force, critical shear force         
      real(8)  phy_kn,phy_nf,phy_cnf,phy_ftmax
!.....dimension of b_laccot
      integer b_ntar
!.....Flag for friction force calculation    
      integer fric_option
!.......normal stiffness,
      real(8) kn
!.....shear stiffness,friction coefficient,bond shear stiffness,friction coefficient for bond
      real(8) kt,fric_coef,b_kt,b_fric_coef
!.....tangential force, tangential forces for bond
      real(8) ft(0:3),b_ft(0:3)
!.....output parameters for spheres
      integer nodes, nelem
!.....previous tengential force for particles with and without bond; pre_ft1:periodic boundary particles
      real(8),allocatable::pre_ft(:,:,:),b_pre_ft(:,:,:),pre_ft1(:,:,:)
!
!...........................Definition of boundary points(for IBM).......................
      type bound_pts
!.......coordinate of boundary points,hydrodynamic forces of fluid nodes at interface
        real(8) coorp(3),hf(3)
!.......angle coordinate
        real(8) angle(3)
!.......(reference or desired) velocity of boundary nodes
        real(8) vel_des(3)
      end type
!
!..............................................Definition of type particle property.......................................................
      type par_prop
!.......coordinate, radius and mass
        real(8) coor(3),radius,mass
!.......flag for the existence of particle,  flag:particle type(0:stationary,  1:moving)
        integer active,flag
!.......particle velocities, hydrodynamic forces and particle-particle contact forces
        real(8) U(6),F(6),FC(6)
!.....................................................................................
!.......angular velocity, moment of inertia, angle of rotation
        real(8) av(3),moi,aor(3)
!.....................................................................................
!.......nb:number of boundary nodes for each particle   ni: number of interior nodes for each particle
!.......nf:number of fluid boundary nodes for each particle
        integer nb,ni,nf
!.......pb: pointer pointing to the boundary node list, 
!.......pn: pointer pointing to the interior node list  pf(:):pointer to the moving particle    
        integer pb,pn,pf 
      end type
!---------------------------------
!.....for periodic boundary output
      integer num_vir
      type part_virtual
!.......coordinate, active flag to use
        real(8) coor(3)
        integer active
      end type
!
      type(bound_pts),allocatable::bpoint(:,:)
      type(par_prop),allocatable::particle(:)
      type(part_virtual),allocatable::part_vir(:)
!
    End module solid
!
!----------------------------------------------------------fluid data---------------------------------------------------------------------
    Module fluid
    implicit none
!.....grid size in x- and y-dimension, z-dimention
      integer nx,ny,nz
!.....fluid densities
      real(8) d0
!.....relaxation parameter & fluid viscosity  tao2=tao*tao   rtao=1.d0/tao
      real(8) tao,visco,tao2,rtao
!.....type of boundary conditions 1=periodic; 2=specified pressure(density) at inlet and outlet
!.............................................3=specified inlet velocities
!.....bc_mode:type of periodic boundary   1,horizontal   2,vertical
      integer bctype,bc_mode
!.....number of fluid boundary edges & edge number
      integer neb
      integer,allocatable::eb(:)
!.....LBM boundary conditions like density or velocity. 
	  real(8),allocatable::bc(:)
!.....coordinate list of fixed wall boundary nodes,non-fixed wall boundary nodes, stationary particle nodes,
!.....moving particle interior nodes,moving boundary nodes and fluid boundary nodes 
	  integer,allocatable::lnw0(:,:),lnw1(:,:),lns(:,:),&
                             lni(:,:),lnb(:,:),lnf(:,:)
!.....wall nodes, pointers to each wall
      integer,allocatable::nw0(:),nw1(:),pw0(:),pw1(:)
!.....number of staionary and moving wall nodes 
      integer now0,now1                             
!.....average velocity, computed by subroutine 'average_velocity',hydraulic radius ratio
      real(8)  vel,hyrad_ratio
!.....xb(1:):bottom boundary;  xb(2:):top boundary;  yb(1:):left boundary;  yb(2:):right boundary, zb(1):front zb(2):back
	    integer,allocatable::xb(:,:),yb(:,:),zb(:,:)
!.....dx: lattice spacing;  dt: time step for LBM;  umax: maximum velocity  cc:lattice speed
      real(8)  dx,dt,umax,cc
!.....fixed edge(wall) list; non-fixed edge(wall) list
!	  integer,allocatable::le0(:),le1(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....fb:logical symbol for fluid boundary  1:true;   mturb:turbulence flag,body force density
!.....fs_body density for IBM method like fi_body
      integer fb,mturb,BODYF
      real(8) fi_body(0:14),w(0:14),fs_body(0:14)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....lmv:list of moving boundary
!	  integer,allocatable::lmv(:)
!.....temporary fluid density distribution function 
      real(8),allocatable::temp(:,:,:,:)
!
!....................................................Definition of lattice node property..................................................
      type node_prop
!.......obstacle type    0:fluid nodes, 1:wall nodes or stationary nodes, 2:moving solid boundary nodes,
!........................3:interior particle nodes, 4:fluid boundary nodes
        integer obst
!.......relazation parameter:tao for turbulence model
        real(8) taostar
!===== for coarsE-grid LBM-DEM coupled by empirical equations
        real(8) velocity(3)
!=====
!.......fluid density distribution function(fdd)
        real(8) fdd(0:14)  
!.......fluid density distribution function(fdd)        
        real(8) body_fs(3)                 
      end type
!===== for coarsE-grid LBM-DEM coupled by empirical equations
      type cell_prop
      ! mean velocity
        real(8) mean_vel(3)
      ! particle volume faction  
        real(8) concentration
        integer num_p
      end type
      type(cell_prop),allocatable::cell(:,:,:)
!=====
      type(node_prop),allocatable::node(:,:,:)


!
    End module fluid
!
!----------------------------------------------------------fluid-solid interaction--------------------------------------------------------
    Module fs_inter
    implicit none
!.....nos: number of stationary nodes,  noi: number of interior nodes
!.....nob: number of boundary nodes,  nof: number of fluid boundary nodes
      integer nos,noi,nob,nof
!.....nnb,nni: previous number of solid noundary nodes and previous number of interior noundary nodes; 
!.....nnf:previous number of fluid noundary nodes
      integer nnb,nni,nnf
!.....number of boundary points of each paricle for IBM
      integer nbp
!.....imwp:logical symbol for interaction between moving wall and particles
	  integer imwp
!.....imwf:logical symbol for interaction between moving wall and fluid nodes
      integer imwf
!.....list of moving wall_fluid nodes corresponding to each wall coordinate
      integer,allocatable::lmwfn1(:,:),lmwfn2(:,:),lmwfn3(:,:),lmwfn4(:,:)
    End module fs_inter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*****************************************************************Main program************************************************************
    Program BPLBM3D_NEW
    use system
    use solid
    use fluid
    use fs_inter
    use MRT, only : initialize_MRT_lattice, initialize_MRT_matrices, MRT_collision
    implicit none
!.....function cputime daclaration and deactive particles
      real(8),external::cputime      
      integer,external::inactive_particle
      character(20) date,start_date,end_date,time1
      character(4) start_h,end_h,start_m,end_m,start_s,end_s
      real(8) gatherStart, gatherEnd

      real(8) ux2,uy2,uz2,den2
      integer ix2, iy2, iz2

      integer :: c0,c1,c2,c3,c4

      !MRT varibles
      integer :: cx_mrt(0:14)
      integer :: cy_mrt(0:14)
      integer :: cz_mrt(0:14)

!.....Read parameter file
      call read_parameters
!.....initial CPU time 
      call date_and_time(start_date,time1)
      start_h=time1(1:2)
      start_m=time1(3:4)
      start_s=time1(5:6)
!.....if necessary, convert some variables from physical to lattice unit.
      if(LBM)then
        call boundary_max_coor
        nx=xmax/dx
	      ny=ymax/dx
        nz=zmax/dx
        xmax=nx
        ymax=ny
        zmax=nz
!        write(*,*) nx,ny,vertex(1,1),vertex(2,1)
!.......convert bond parameters
        if(bond)then
          cb_kn=phy_kn/cc/cc
          cb_cnf=phy_cnf/cc/cc/dx
          b_ftmax=phy_ftmax/cc/cc/dx
          b_kt=b_kt/cc/cc
        end if
!
!.......Scaled gravity or unit convertion
      if(BODYF.eq.1)then
	      gacce=gacce*dt/cc
      else
        gacce=0.d0
      endif
        gacce1=(1.0d0-d0/ds)*9.81d0*dt/cc
!        gacce1=0.d0
!.......physical time
        time=dt
!.......subcycles of particle contact loop
        if(nm.gt.0)then
	        nsub=1.0/dtime+0.999999
!...........reset particle timestep and penalty
	        dtime=1.0d0/nsub
	        kn=kn/cc/cc
          kt=kt/cc/cc
	      else
	        nsub=1
	        dtime=1.0d0
        end if
!=========== coarse-grid        
        ! if(empirical)then
        !   nsub=100
        !   dtime=1.0d0/nsub
        !   write(*,*) nsub,dtime
        ! endif
!.......Startup information message
	    write(*,1000) nx,ny,nz,dt,cc,nsub        
      else
        if(bond)then
          cb_kn=phy_kn
          cb_cnf=phy_cnf
          b_ftmax=phy_ftmax
        end if
!        gacce=9.81d0
!        gacce=0.0d0
        time=dtime
      end if     
!
!----------------------------------------------------------------------
!.....Initialization and allocate space for allocatable array
      buff=minrad*0.1d0
      vmax=0.0d0
      umax=0.0d0
      sum_disp=0.0d0
!.....initializtion simulation time and first time step
!      start_time=cputime()
!
	    first=.true.
      mcnts=np*20
!.....Particle active
	  do i=1,np
	    particle(i)%active=1
      particle(i)%FC=0.0d0
	  enddo
!  contact force on walls    
    do i=1,ne
      wall(i)%FC=0.0d0
    end do       
!
!.....Allocate spaces
    if(LBM)then
!.......weighting coefficient
      w(0)=2.0d0/9.0d0
      do i=1,6
        w(i)=1.0d0/9.0d0
      enddo
      
      do i=7,14
        w(i)=1.0d0/72.0d0
      enddo

!----------------------------------------------------------------------
!.....Initialize MRT lattice once, before the main loop
!----------------------------------------------------------------------
      if(use_MRT)then
        ! MRT needs information from d3q15
        cx_mrt = (/ &
             0, &
             1, -1,  0,  0,  0,  0, &
             1, -1,  1, -1,  1, -1,  1, -1 /)

        cy_mrt = (/ &
             0, &
             0,  0,  1, -1,  0,  0, &
             1, -1,  1, -1, -1,  1, -1,  1 /)

        cz_mrt = (/ &
             0, &
             0,  0,  0,  0,  1, -1, &
             1, -1, -1,  1,  1, -1, -1,  1 /)

        call initialize_MRT_lattice(cx_mrt, cy_mrt, cz_mrt, w )
             
        call initialize_MRT_matrices(tao, S1, S2, S4, S14, d0)

        ! write(*,*) "MRT check: tau = ", tao
        ! write(*,*) "MRT check: S1  = ", S1
        ! write(*,*) "MRT check: stress S = ", 1.0d0/tao
        ! write(*,*) "MRT check: S = ", S

      endif


!.......fluid boundary for moving particles.  check it later????
	    fb=1   
      allocate(vertex_mw(3,nv),vertex_sw(3,nv))
      allocate(node(0:nx,0:ny,0:nz),lnw0(3,3*(nx*ny+ny*nz+nx*nz)),&
      lnw1(3,3*(nx*ny+ny*nz+nx*nz)),temp(0:14,-1:nx+1,-1:ny+1,-1:nz+1))
      if(ne0.ge.1) allocate(nw0(ne0),pw0(ne0))
      if(ne1.ge.1) allocate(nw1(ne1),pw1(ne1))
      if(IMB) allocate(lni(3,MI),lnb(3,MB),lnf(3,3*MB))
!....................................................
        if(IBM)then
          if(nm.gt.0) allocate(bpoint(nbp,np))
        end if    
! coarse-grid LBM-DEM        
        if(empirical)then
          allocate(cell(1:nx,1:ny,1:nz))
          do i=1,nz
            do j=1,ny
              do k=1,nx
              cell(k,j,i)%num_p=0
              cell(k,j,i)%mean_vel=0.0d0
              enddo
            enddo
          enddo
        endif     
!....................................................
      end if
!    only consider ghost or virtual particles for moving particles
      if(boundary_dem==1) allocate(part_vir(nm))
      allocate(laccoc(np),laccop(np),icoor(6,np),xb(3,0:nx),yb(3,0:ny),zb(3,0:nz))
      allocate(pre_ft(0:3,np,np),pre_ft1(0:3,np,np),laccot(mcnts))
!.....pre_ft=0
      do i=1,np
        do j=1,np
          do k=0,3
            pre_ft(k,j,i)=0.0d0
          enddo
        enddo
      enddo 
!.....pre_ft1=0
      do i=1,np
        do j=1,np
          do k=0,3
            pre_ft1(k,j,i)=0.0d0
          enddo
        enddo
      enddo
!  PBC for moving particles
      if(boundary_dem==1)then
      do i=1,nm
        part_vir(i)%active=0
      end do
      end if
!       
      do i=1,np
        icoor(:,i)=0
      end do
!     
      if (bond) then
!.....Sample preparation method by gap=gap-gap_initial
        allocate(b_pre_ft(0:3,np,np),b_laccoc(np),b_laccot(10,np))    
        do i=1,np
          do j=1,np
            do k=0,3
              b_pre_ft(k,j,i)=0.0d0
            enddo
          enddo
        enddo   
        if (sam_prep2) then
          allocate(gap_ini(np,np))
          gap_ini(:,:)=0.0d0
        end if
        b_laccot=0
      endif
!----------------------------------------------------------------------
!
!-----------------------------DEM part---------------------------------
!.....Boundary bounding box
!      call boundary_boundingbox
!.....used for check subroutine boundary_boundingbox
!      do j=1,ne0+ne1
!        write(*,*) icoor(1,np+j),icoor(2,np+j),icoor(3,np+j),icoor(4,np+j)
!      end do
!
	  if(nm.gt.0)then        
!.......Initialise moving particle velocities & forces
        p1=>head_lpm
        do while(.true.)
          ip= p1%num
!          write(*,*) ip
          do j=1,6
	        particle(ip)%U(j)=0.0d0
	        particle(ip)%F(j)=0.0d0
	      enddo
! deal with rotation in the future ??
          particle(ip)%av=0.0d0
          particle(ip)%aor=0.0d0
!.........calculate moment of inertia for particles
          particle(ip)%moi=0.4*particle(ip)%mass*(particle(ip)%radius**2)
!                 
          if(.not. associated(p1%next)) exit
          p1=>p1%next
        end do
!
!.......initial contact detection for bond installation
	    if(bond)then
!          call contact_search_first        
!.......Sort out Bond list
!          call intall_bond
        end if
      end if
!------------------------------------------------------------------------
!
!--------------------------LBM part--------------------------------------
      if(LBM)then
!.....Open output files
!	      open(11,file='BPLBM_model.plt')  
	      open(12,file='model_info.dat')
	      open(13,file='veloc.dat')
        open(14,file='aveloc.dat')
        open(20,file='particle_force.txt')
        open(21,file='particle_posi.txt')
!        open(34,file='spheres.plt')
        open(unit=98, file='paraview/Fluids.pvd', status='replace')
        open(unit=99, file='paraview/particles.pvd', status='replace')
        !open(unit=199, file='paraview/particles_.csv', status='replace')

        write(98,'(A)') '<?xml version="1.0"?>'
        write(98,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
        write(98,'(A)') '  <Collection>'

        write(99,'(A)') '<?xml version="1.0"?>'
        write(99,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
        write(99,'(A)') '  <Collection>'


!.......Generate boundary_profile.  deal with stationary walls ???
!        call boundary_profile
!.......Initialise dentisies
        call init_density
! 
!.....Generate wall boundary
      call generate_wall_boundary
!.......Generate stationary particles
        if(ns.gt.0)call generate_stationary_particles
        if(nm.gt.0)then
!.................................................................................
          if(IMB)then
!...........Generate solid nodes within moving particles & fluid nodes outside moving particles
            call generate_moving_particles
!...........Generate fluid nodes outside moving particles
!	          call generate_fluid_boundary_nodes
!.........generate boundary points of each particles for Immersed Boundary Method
          else if(IBM)then
            call generate_boundary_nodes
          end if 
!.................................................................................
        end if
!.......For turbulence model
        if(mturb.eq.1)then
	      tao2=tao*tao
          do iz=0,nz
	          do iy=0,ny
              do ix=0,nx
    	          node(ix,iy,iz)%taostar=tao
              enddo
            end do
          end do
        else
	      rtao=1.d0/tao
        endif
      else
        ! open(11,file='BPM_model.plt') 
        open(unit=99, file='paraview/particles.pvd', status='replace')
      end if
!	read sphere-output-template
!	  open(19, file='fe_sphere_block.dat')
!	  read(19,*)nodes,nelem
!	  allocate(sphere_mesh(3,nodes))
!	  allocate(sphere_elem(4,nelem))
!	  call sphere_output_temp
!	  nodes=nodes-1
!	  close(19)
!--------------------------------------------------------------------------
!
!--------------------------Main Loop---------------------------------------
      loopmain: do istep=1, msteps
!.......Particle interactions
        if(nm.gt.0.and.istep.ge.1)then   
!          write(*,*) "Loop", istep     
	        if(LBM)then
!...........................DEM part.......................................
!.........Subcycling
            loopinner: do isub=1,nsub
!.............Check if all the moving particles are deactivated
!	          i=inactive_particle()
              if(inactive_particle().eq.nm)then
	            write(*,*) 'ALL the moving particles have been deactivated!'
	            exit loopmain
	            endif
!.............Perform contact detection and compute contact forces
	            if((.not.fix_packing) .and. nm.ge.1)call contact_main
!.............Update moving particle velocity and position
	            call update_moving_particle_position(nm,dtime,gacce1,vmax,isub)
!.............update boundary bounding box
!             call boundary_boundingbox(ne0,le0,le1,ne1,edge,vertex,np,dx,icoor)
!              if(empirical .and. istep==1) exit
	          enddo loopinner
!..........................................................................
            if(IMB)then
!.............Update moving particle boundary and interior node list
!              write(*,*) "IMB", istep 
              call update_moving_particle_nodes
!.............Update moving particle fluid boundary nodes list
!!!	          call update_fluid_boundary_nodes
!              write(*,*) "nodes", istep 
            end if
!...........IBM. need to modify in the future ???
            if(IBM) call update_boundary_nodes
!===== coarse-grid LBM-DEM
            if(empirical) call update_cell_particle_num
  
          else
!----------------
            if(nm.ge.1)call contact_main
!...........Update moving particle velocity and position
	        call update_moving_particle_position(nm,dtime,gacce,vmax,isub)
!...........update boundary bounding box
!            call update_wall_position    
!----------------   
          end if
        endif
!
        if(LBM)then
!.........Introduce periodic boundary condition
          if(bctype.eq.1) call redistribute(bc(1)) 
!          write(*,*) "boundary", istep 
!.........Relaxation
!.........initialise body forces in IBM
          if(IBM)then
            do i=0,nz
            do j=0,ny
              do k=0,nx
                node(k,j,i)%body_fs=0.0d0
              enddo
            enddo
            end do
          end if
!===== coarse-grid LBM-DEM: to get hydro force on particles and body force at nodes
          if(empirical) call empirical_force
!=====       
! mturb=1 will be revised in the future
!          if(mturb.eq.1)then
! in the future
!        	 call relaxation_turbulence
!	        else
            call relaxation
!	        endif
!write(*,*) "relaxation", istep 
!.........Bounce back from obstacles:for fixed wall and stationary particles
          call bounceback
!
!.........Propagation
          call propagate
!
!.........Boundary treatment: apply different boundary conditions Periodic Boundary and Pressure Boundary. 
!         Other boundaries like Velocity Boundary etc. are not included and should be added
	        if(bctype.gt.0) then
!            if(ne==6)then
              call boundary_treatment_regular
!            else
!              call boundary_treatment
!            end if
          end if
!.......The integral fluid density is checked regularly.

          if(mod(istep,screen_steps).eq.0) call check_density(nx,ny,nz,istep,time)
        end if
                   
!.......Output results
!           write(*,1001) istep, particle(1)%F(1:3)!, particle(1)%coor(1:3)   
        if(LBM)then
          if(istep.eq.1)then
           call write_results1
           call write_spheres
           write(20,1006) istep, particle(1)%coor(1:3),particle(1)%U(1:3),particle(1)%F(1:3)
          endif

          if(mod(istep,outsteps).eq.0)then
           call write_results1
           call write_spheres
           write(20,1006) istep, particle(1)%coor(1:3),particle(1)%U(1:3),particle(1)%F(1:3)
           !.....calculate FLUID velocity at fluid nodes
           call write_velocity(nx,ny,nz)
          endif

          write(21,1001) istep, particle(1)%coor(1:3)  
!          write(20,1001) istep, particle(1)%U(1:3)
        else
          if(istep.eq.1) call write_results0
          if(istep.ge.0.and.(mod(istep,outsteps).eq.0)) call write_results0
        end if
!
!
        if(LBM)then  
          time=time+dt
        else
          time=time+dtime  
          if(mod(istep,screen_steps).eq.0) write(*,1002) istep,time     
        end if
!...................................................................
!.......initailise angle of rotation for next time step
!        write(*,*)particle(10)%aor
        do i=1,np
          particle(i)%aor=0.0d0
        end do
!        write(*,*) "Check end!", istep
!...................................................................
      end do loopmain
!
!.....simulation time
!      end_time=cputime()
      call date_and_time(end_date,time1)
      write(12,*) '----------------Simulation time:----------------'
      write(12,"('***  Start Time:',A10,'  ',A2,':',A2,':',A2,'  ***')") start_date,start_h,start_m,start_s
      write(12,"('***  End   Time:',A10,'  ',A2,':',A2,':',A2,'  ***')") end_date,time1(1:2),time1(3:4),time1(5:6)
      write(12,*) '---------------Outputting result----------------'
!.....generate sample
      output=.false.
      if(output)then
        open (15,file='particleprop.dat')
        do i=1,np
          write(15,*) particle(i)%flag,particle(i)%coor(1),particle(i)%coor(2),particle(i)%coor(3),particle(i)%radius
        end do
        close(15)
      end if
!.......................................................
      if(LBM)then
!.....calculate reynold number
        if(np.gt.0)then
          call comp_rey(visco,minrad,umax,cc,dx,dt)
        else
          call comp_rey(visco,ny*1.0d0,umax,cc,dx,dt)
        end if
      end if
!      write(*,*) "LBM ouput", istep
!.......................................................
! for closing paraview output
      write(98,'(A)') '  </Collection>'
      write(98,'(A)') '</VTKFile>'
      close(98)
      write(99,'(A)') '  </Collection>'
      write(99,'(A)') '</VTKFile>'
      close(99)

!      close(11)
      if(LBM)then
	    close(12)
	    close(13)
        close(14)
        close(20)
        close(21)
!        close(34)
      end if
!
!-----------------------------------------------------------------------------------------------
!.....Deallocate spaces 
     deallocate(particle)
     deallocate(laccoc,laccop,icoor,xb,yb,pre_ft,pre_ft1,laccot)
     if(boundary_dem==1) deallocate(part_vir)

     if(LBM)then
       deallocate(vertex_mw,vertex_sw,node,lnw0,lnw1,temp)
       if(ne0.ge.1) deallocate(nw0,pw0)
      if(ne1.ge.1) deallocate(nw1,pw1)
       if(bctype.eq.1) deallocate(bc)
       if(bctype.eq.2) deallocate(eb)
       if(IMB) deallocate(lni,lnb,lnf)
     end if
     if(bond)then
       deallocate(b_pre_ft,b_laccoc,b_laccot)
       if(sam_prep2) deallocate(gap_ini)
     end if
     if(ne0.gt.0) deallocate(head_le0)
     if(ne1.gt.0) deallocate(head_le1)
     if(ns.gt.0) deallocate(head_lps)
     if(nm.gt.0) deallocate(head_lpm)
     if(nos.gt.0) deallocate(head_lns)
!------------------------------------------------------------------------------------------------
!
!.....simulation end
      write(6,*) '********************    end     ********************'
!
    1001    format(i6,1x,g12.4,1x,g12.4,1x,g12.4)
    1006    format(i6,9(1x,g12.4))
    1002    format('** Step =',i7, ' Time =',g11.4, ' **')
    1000    format('****************************************************'/,&
                    '***          Lattice size  nx = ',i8, /&
                    '***                        ny = ',i8, /&
                    '***                        nz = ',i8, /&
                    '***              TimeStep  dt = ',g12.4, /&
                    '***          Lattice Speed  C = ',g12.4, /&
                    '***          No. SubCycles ns = ',i8, /&
                    '****************************************************'//) 
!
      stop
    End program BPLBM3D_NEW
!*****************************************************************************************************************************************
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine read_parameters
    use system
    use solid
    use fluid
    use fs_inter
    implicit none
!
!      real(8), parameter::PI=3.14159265385 
!.....origin of coordinates
      real(8) x0,y0,z0
!.....initialisation
!      MS=0
	    MI=0
	    MB=0
      ns=0
	    nm=0
      ne0=0
  	  ne1=0
	    maxrad=0.0d0
	    minrad=1000.d0
!
!.....open system parameter file
      open(10,file='parameter3D.dat')
!.....read DEM flag,logical symbol for turbulence model
      read(10,*) LBM,mturb,BODYF,hyrad_ratio
!.....read IBM flag,IMB flag,number of boundary points of each particle
      read(10,*) IBM,IMB,nbp,empirical
      read(10,*) sam_prep2,lubrication,uniaxial
!.....read if MRT (else is SRT/BGK), S1, 2, 4, 14 are relaxation value
      read(10,*) use_MRT
      if(use_MRT) read(10,*) S1, S2, S4, S14
!.....open particle parameter file      
      open(15,file='ballproperty.dat') 
!
!.....boundary information
!
      read(10,*) xmax,ymax,zmax
!      read(10,*) nv,ne
! current box for cuiboid only ???
      nv=8
      ne=6
!.....allocate memory for vertex(:,:)
      allocate(vertex(3,nv),wall(ne))
      do i=1,ne
        wall(i)%vel=0.0d0
      end do
!.....read vertex data and convert the first vertex to origin point
!	  do i=1,nv
!	    read(10,*) vertex(1,i),vertex(2,i)
!vertex 1
vertex(1,1)=0.0d0
vertex(2,1)=0.0d0
vertex(3,1)=0.0d0
!vetex 2
vertex(1,2)=xmax
vertex(2,2)=0.0d0
vertex(3,2)=0.0d0
!vetex 3
vertex(1,3)=xmax
vertex(2,3)=0.0d0
vertex(3,3)=zmax
!vetex 4
vertex(1,4)=0.0d0
vertex(2,4)=0.0d0
vertex(3,4)=zmax
!vertex 5
vertex(1,5)=0.0d0
vertex(2,5)=ymax
vertex(3,5)=0.0d0
!vertex 6
vertex(1,6)=xmax
vertex(2,6)=ymax
vertex(3,6)=0.0d0
!vertex 7
vertex(1,7)=xmax
vertex(2,7)=ymax
vertex(3,7)=zmax
!vertex 8
vertex(1,8)=0.0d0
vertex(2,8)=ymax
vertex(3,8)=zmax
!	    if(i.eq.1)then
!	      x0=vertex(1,1)
!	      y0=vertex(2,1)
!	    endif
!	    vertex(1,i)=vertex(1,i)-x0
!	    vertex(2,i)=vertex(2,i)-y0
!	  enddo
!.....read start, end points of walls and wall type
!	  do i=1,ne
!	    read(10,*) wall(i)%point(1),wall(i)%point(2),wall(i)%flag, wall(i)%rate
!        wall(i)%coor(1:2)=vertex(:,wall(i)%point(1))
!        wall(i)%coor(3:4)=vertex(:,wall(i)%point(2))
!        write(*,*) wall(i)%coor
!        if(wall(i)%rate.eq.1)then
!          read(10,*) wall(i)%vel(1),wall(i)%vel(2)
!          write(*,*) wall(i)%vel(1),wall(i)%vel(2)
!        end if
!        if(wall(i)%flag.eq.0)then
!	      ne0=ne0+1
!........................................................................
!          if(ne0.eq.1)then
!            allocate(head_le0)
!            head_le0%num=i
!            nullify(head_le0%next)
!            p1=>head_le0
!          else
!            allocate(p1%next,stat=err)
!            if(err/=0) then
!              write(*,*) "Out of memory!"
!              stop
!            end if
!            p1=>p1%next
!            p1%num=i
!            nullify(p1%next)
!          end if         
!........................................................................
!	      le0(ne0)=i
!	    else if(wall(i)%flag.eq.1)then
!	      ne1=ne1+1
!........................................................................
!          if(ne1.eq.1)then
!            allocate(head_le1)
!            head_le1%num=i
!            nullify(head_le1%next)
!            p2=>head_le1
 !         else
!            allocate(p2%next,stat=err)
 !           if(err/=0) then
!              write(*,*) "Out of memory!"
 !             stop
!            end if
!            p2=>p2%next
!            p2%num=i
!            nullify(p2%next)
!          end if
          
!.......................................................................
!	      le1(ne1)=i
!	    endif
!	  enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....ouput array to check pointer code
!      pointer_check=0
!      if(pointer_check)then
!        p1=>head_le0
!        do while(.true.)
!          write(*,*) p1%num
!          if(.not. associated(p1%next)) exit
!          p1=>p1%next
!        end do
!
!        p1=>head_le1
!        do while(.true.)
!          write(*,*) p1%num
!          if(.not. associated(p1%next)) exit
!          p1=>p1%next
!        end do
!      end if      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!.......lattice spacing
	    read(10,*) dx
!.......initial value for fluid density, and solid density
        read(10,*) d0,ds
!.......relaxation parameter & fluid voscosity
        read(10,*) tao,visco
!.......number of particles, and the flag to make all moving particles fixed enabling hydrodynamic force calculation using IMB
        read(10,*) np,fix_packing
!.......allocate memory for particle array
        allocate(particle(np))
!.......particle nromal penalty, damping ratio and timestep factor
        read(10,*) kn,kt,damp_option,damp_coef,factor
!.......friction option, friction coefficient
        read(10,*) fric_option,fric_coef
!
        read(10,*) imwp
        read(10,*) imwf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!        do i=1,ne
!          read(10,*) wall(i)%wfluid,wall(i)%vel(1),wall(i)%vel(2)
!        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!.......parameters for bond model
        read(10,*) bond
        read(10,*) phy_kn, phy_cnf
        read(10,*) b_kt, phy_ftmax, b_fric_coef
!
!        if(LBM)then
!          do i=1,ne
!          do j=1,4
!            wall(i)%coor(j)=wall(i)%coor(j)/dx
!          end do
!          end do
!        end if
!.......read particle parameter
        if(np.gt.0)then
          do 20 i=1,np
!...........particle info: flag (0=stationary, 1=move), center coordinates and radius
            read(15,*) particle(i)%flag,particle(i)%coor(1),particle(i)%coor(2),particle(i)%coor(3),particle(i)%radius
!...........Convert to the lattice space
	        if(LBM)then
            particle(i)%coor(1)=(particle(i)%coor(1))/dx
	          particle(i)%coor(2)=(particle(i)%coor(2))/dx
            particle(i)%coor(3)=(particle(i)%coor(3))/dx
	          particle(i)%radius=particle(i)%radius/dx          
	        end if
!...........find maximum and minimum radius
            if(maxrad.lt.particle(i)%radius) maxrad=particle(i)%radius
            if(particle(i)%radius.lt.minrad) minrad=particle(i)%radius
!
!...........get stationary and moving partcle number and list
	        if(particle(i)%flag.eq.0)then
!...............Estimate maximum no of nodes inside partilces
!              if(LBM)then
!	            MS=MS+4*(particle(i)%radius+1)*(particle(i)%radius+1)
!              end if
		      ns=ns+1
              if(ns.eq.1)then
                allocate(head_lps)
                head_lps%num=i
                nullify(head_lps%next)
                p1=>head_lps
              else
                allocate(p1%next,stat=err)
                if(err/=0) then
                  write(*,*) "Out of memory!"
                  stop
                end if
                p1=>p1%next
                p1%num=i
                nullify(p1%next)
              end if             
!          
!	          lps(ns)=i
	        else if(particle(i)%flag.eq.1)then
              if(LBM.and.IMB)then
!...............Estimate maximum no of nodes inside moving partilces
                MI=MI+8*(particle(i)%radius+1)**3              
!...............Estimate maximum no of boundary nodes of moving partilces
  	            MB=MB+24*(particle(i)%radius+1)**2
              end if 
		      nm=nm+1
              if(nm.eq.1)then
                allocate(head_lpm)
                head_lpm%num=i
                nullify(head_lpm%next)
                p2=>head_lpm
              else
                allocate(p2%next,stat=err)
                if(err/=0) then
                  write(*,*) "Out of memory!"
                  stop
                end if
                p2=>p2%next
                p2%num=i
                nullify(p2%next)
              end if
!         
!              lpm(nm)=i
              if(LBM)then
	            particle(i)%mass=4.0d0/3.0d0*ds*PI*(particle(i)%radius**3)
              else
                particle(i)%mass=4.0d0/3.0d0*ds*PI*(particle(i)%radius**3)
              end if
	        endif
     20   continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.........ouput array to check pointer code
          pointer_check=0
          if(pointer_check.eq.1)then
            p1=>head_lps
            do while(.true.)
              write(*,*) p1%num
              if(.not. associated(p1%next)) exit
              p1=>p1%next
            end do
!
            p1=>head_lpm
            do while(.true.)
              write(*,*) p1%num
              if(.not. associated(p1%next)) exit
              p1=>p1%next
            end do
          end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
!.........Estimate critical time step for particle contact  should we convert kn before timestep calculation???
	      dtime=factor*2.0*minrad*sqrt(PI*ds*minrad/3.0/kn)
!          write(*,*) dtime
	    endif
!
!.......number of iterations, screen output frequency, file output frequency
        read(10,*) msteps,screen_steps,outsteps
!.......read DEM boundary
        read(10,*) boundary_dem
!
!..........................................................................................................
!.......If consider fluid read the following data and convert them
        if(LBM)then
!.........type of boundary condition
          read(10,*) bctype
!.........type 0: all wall fixed boundary condition
          if(bctype.eq.0)then
!.........type 1: periodic boundary condition
          else if(bctype.eq.1)then
!.........acceleration
            allocate(bc(1))
!            read(10,*) bc(1),bc_mode
              read(10,*) periodic_flag,implement
              read(10,*) gacce
!.........type 2: specified pressure/density
          else if(bctype.eq.2)then
!...........number of pressure edges 
            read(10,*)neb
!...........allocate eb(:) memory
            allocate(eb(neb),bc(neb))
	        do i=1,neb
!.............edge number and density
              read(10,*) eb(i),bc(i)
	        enddo
!.........type 3: specified inlet velocities
          else if(bctype.eq.3) then
            allocate(bc(2))
            read(10,*) periodic_flag,implement
            read(10,*) bc(1),bc(2)
	      endif          
!
!.........LB time step based on relaxation time and fluid vicosity
	      dt=(tao-0.5d0)*dx*dx/visco/3.d0
!           write(*,*) dt
!.........Lattice speed
          cc=dx/dt

!       unit convertion
!
        if(bctype.eq.3) then
            bc(1)=bc(1)/cc
            bc(2)=bc(2)/cc
	      endif   
!.........Estimate critical time step for particle contact
	      dtime=factor*4.0*sqrt(PI*ds*(minrad*dx)**3.0/kn/3.0)/dt
!           write(*,*) dtime,factor,PI,ds,minrad*dx,kn,dt,tao,visco
!           stop
        end if
!..........................................................................................................
!    
        close(10)
        close(15) 
!                            
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine boundary_max_coor
    use system
    use fluid
    implicit none
    real(8) xs,ys,zs
!
      xmax=0.d0
	    ymax=0.d0
      zmax=0.0d0
      do 10 i=1,nv
        xs=vertex(1,i)
        ys=vertex(2,i)
        zs=vertex(3,i)
	      vertex(1,i)=vertex(1,i)/dx
	      vertex(2,i)=vertex(2,i)/dx
        vertex(3,i)=vertex(3,i)/dx

	    if(xs.gt.xmax)xmax=xs
	    if(ys.gt.ymax)ymax=ys
      if(zs.gt.zmax)zmax=zs
   10 continue
      return
	End subroutine
!


!============================================================================
!=======for coarse-grid LBM-DEM
  Subroutine update_cell_particle_num
    use system
    use solid
    use fluid
    implicit none
    !
    do k=1,nz
      do j=1,ny
        do i=1,nx
           cell(i,j,k)%num_p=0
        enddo
      enddo
    enddo
    do ip=1,np
       ix=floor(particle(ip)%coor(1))+1
       iy=floor(particle(ip)%coor(2))+1
       iz=floor(particle(ip)%coor(3))+1
       cell(ix,iy,iz)%num_p=cell(ix,iy,iz)%num_p+1
       ! comment later
       write(*,*) istep,ix,iy,iz,cell(ix,iy,iz)%num_p
    end do


    return
	End subroutine
!

  Subroutine empirical_force
    use system
    use solid
    use fluid
    implicit none
    real(8) ux,uy,uz,den,relative_vel(3)
    real(8) f1, f2, f3, Rem,concentration
    real(8) B1,B2, PFP
    !
    do ip=1,np
       ix=floor(particle(ip)%coor(1))
       iy=floor(particle(ip)%coor(2))
       iz=floor(particle(ip)%coor(3))
       ! get nodal fluid velocity
       cell(ix+1,iy+1,iz+1)%mean_vel(:)=0.0d0
       do k=iz,iz+1
         do j=iy,iy+1
           do i=ix,ix+1
              call nodal_velocity(i,j,k,nx,ny,nz,ux,uy,uz,den)
              ! node(i,j,k)%velocity(1)=ux
              ! node(i,j,k)%velocity(2)=uy
              ! node(i,j,k)%velocity(3)=uz
              cell(ix+1,iy+1,iz+1)%mean_vel(1)=cell(ix+1,iy+1,iz+1)%mean_vel(1)+ux
              cell(ix+1,iy+1,iz+1)%mean_vel(2)=cell(ix+1,iy+1,iz+1)%mean_vel(2)+uy
              cell(ix+1,iy+1,iz+1)%mean_vel(3)=cell(ix+1,iy+1,iz+1)%mean_vel(3)+uz
           enddo
         enddo
       enddo
!  !      real mean velocity
!        cell(ix+1,iy+1,iz+1)%mean_vel(1)=cell(ix+1,iy+1,iz+1)%mean_vel(1)/8.0
!        cell(ix+1,iy+1,iz+1)%mean_vel(2)=cell(ix+1,iy+1,iz+1)%mean_vel(2)/8.0
!        cell(ix+1,iy+1,iz+1)%mean_vel(3)=cell(ix+1,iy+1,iz+1)%mean_vel(3)/8.0

       !particle volume fraction
       cell(ix+1,iy+1,iz+1)%concentration=cell(ix+1,iy+1,iz+1)%num_p*&
                 4.0d0/3.0d0*PI*(particle(ip)%radius**3)
      concentration=cell(ix+1,iy+1,iz+1)%concentration
!       concentration=0.001
     if(concentration>=0.999d0)then
       write(*,*) "Simulation crushed, because the fluid grid size is too small!!! "
       stop
     endif
!     comment later
!      write(*,*) "Concentration for particle ID: ", concentration, ip,4.0d0/3.0d0*PI*(particle(ip)%radius**3)

!calculate relative velocity
       do k=1,3
        !      real mean velocity
         cell(ix+1,iy+1,iz+1)%mean_vel(k)=cell(ix+1,iy+1,iz+1)%mean_vel(k)/8.0
         relative_vel(k)=particle(ip)%U(k)-cell(ix+1,iy+1,iz+1)%mean_vel(k)
!         relative_vel(k)=0.18428
        !    get hydroforce and Mean velocity Re number
        !------- Re should be calculated by velocity magnitude, for single particle sedimentation the current way is ok
         Rem=abs(relative_vel(k)*2.0*particle(ip)%radius/(visco/cc/dx) * (1-concentration))
         if(k==2) write(*,*) "relative vel and Rem: ", relative_vel(k), Rem, concentration
        ! Tenneti model 2011
        !Rem=3.564592895
         f1 = (1+0.15*Rem**0.687)/ (1-concentration)**3.0
!         write(*,*) "num ",(1+0.15*Rem**0.687)
         f2 = 5.81*concentration/(1-concentration)**3.0 + (0.48 * concentration**(1.0/3.0))/(1-concentration)**4.0
         f3 = concentration**3.0 * Rem * (0.95 + 0.61*concentration**3.0/(1-concentration)**2.0)
!         write(*,*) "f1,f2,f3: ", f1,f2,f3
!         stop
          !particle(ip)%F(k) = -(f1+f2+f3) * 6.0*PI * (particle(ip)%radius*dx) * (visco*d0) * (relative_vel(k)*cc)/cc/cc/dx/dx
         !particle(ip)%F(k) = -(f1+f2+f3) * (1-concentration) * 6.0*PI * (particle(ip)%radius*dx) * (visco*d0) * (relative_vel(k)*cc)/cc/cc/dx/dx
         particle(ip)%F(k) = -(f1+f2+f3) * (1-concentration) * 6.0*PI * (particle(ip)%radius) * (visco*d0) /cc/dx * (relative_vel(k))
! calculate PFP stress  ! only for single sedimentation inflow direction
!          if (k==2) then
! !          Rem=Rem
!            B1=0.95*(Rem/(1-concentration))**0.02-0.01*log(concentration)-1.0915
!            B2=0.4046*(Rem/(1-concentration))**(-0.3)-0.0412
!            PFP=d0*concentration*(B1*relative_vel(k)**2.0+B2*relative_vel(k)**2.0)
!            if (Rem<10) then
!               particle(ip)%F(k)=particle(ip)%F(k)+PFP*particle(ip)%mass/dx
!            else
!               particle(ip)%F(k)=particle(ip)%F(k)+PFP*particle(ip)%mass
!            endif
!          endif
       enddo
        write(*,*) "hydro forces total and PFP: ",particle(ip)%F(1:3),PFP*particle(ip)%mass/ds
!         stop

! distribute to fluid grid nodes

    end do
!minwang

    return
  End subroutine
!==================================================================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine boundary_boundingbox
    use system
    use solid
    use fluid
	implicit none
      integer::MAGNIFY=10000
      real(8) xs,xe,ys,ye,tmp,buff1
      
	  buff1=dx
!.....Loop over fixed wall boundaries
	    if(ne0>0)then
        p1=>head_le0
        do while(.true.)
          j= p1%num
          is=wall(j)%point(1)
          ie=wall(j)%point(2)
!.........xs,xe are the start and end x coordinate; ys,ye are for y coordinate
          xs=vertex(1,is)
	      ys=vertex(2,is)
	      xe=vertex(1,ie)
	      ye=vertex(2,ie)
          if(xs.gt.xe)then
	        tmp=xs
	        xs=xe
	        xe=tmp
	      endif
	      if(ys.gt.ye)then
	        tmp=ys
	        ys=ye
	        ye=tmp
	      endif
          icoor(1,np+j)=(xs-buff1)*MAGNIFY
	      icoor(2,np+j)=(ys-buff1)*MAGNIFY
	      icoor(3,np+j)=(xe+buff1)*MAGNIFY
	      icoor(4,np+j)=(ye+buff1)*MAGNIFY
!
          if(.not. associated(p1%next)) exit
          p1=>p1%next       
        end do
        end if
!        
!.....  Loop over moving wall boundaries
	    if(ne1>0)then
        p1=>head_le1
        do while(.true.)
          j= p1%num
          is=wall(j)%point(1)
          ie=wall(j)%point(2)
!.........xs,xe are the start and end x coordinate; ys,ye are for y coordinate
          xs=vertex(1,is)
	      ys=vertex(2,is)
	      xe=vertex(1,ie)
	      ye=vertex(2,ie)
          if(xs.gt.xe)then
	        tmp=xs
	        xs=xe
	        xe=tmp
	      endif
	      if(ys.gt.ye)then
	        tmp=ys
	        ys=ye
	        ye=tmp
	      endif
          icoor(1,np+j)=(xs-buff1)*MAGNIFY
	      icoor(2,np+j)=(ys-buff1)*MAGNIFY
	      icoor(3,np+j)=(xe+buff1)*MAGNIFY
	      icoor(4,np+j)=(ye+buff1)*MAGNIFY
!
          if(.not. associated(p1%next)) exit
          p1=>p1%next       
        end do
        end if
!
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine intall_bond
    use system
    use solid
    implicit none          
!.....loop each particle
      loop1: do i=1,np
        if(particle(i)%active.eq.0) cycle loop1
!.......initialise target number for each bonded particle
        b_ntar=0
!.......target number in contact
        ntar=laccoc(i)
        if(ntar.eq.0) cycle loop1
        ip=laccop(i)-1
!        write(*,*) ip
!.......loop the target particles
        loop2: do it=1,ntar
          j=laccot(ip+it)
          if(j.ge.i .or. particle(j)%active.eq.0) cycle loop2
!.........critical bond gap for bond installation 
          b_install=0.0001d0*(particle(i)%radius+particle(j)%radius)
!.........distance tween particles
          distance=sqrt((particle(i)%coor(1)-particle(j)%coor(1))*(particle(i)%coor(1)-particle(j)%coor(1))+&
          (particle(i)%coor(2)-particle(j)%coor(2))*(particle(i)%coor(2)-particle(j)%coor(2)))
!.........gap between particles
          gap=distance-(particle(i)%radius+particle(j)%radius)
          if(gap<=b_install)then
            b_ntar=b_ntar+1
            b_laccot(b_ntar,i)=j
            if(sam_prep2)then
            gap_ini(j,i)=gap
            end if
!            write(*,*) i,j,gap_ini(j,i)
          end if  
        end do loop2
        b_laccoc(i)=b_ntar
      end do loop1                          
!...........Check bond list
!            open(101,file='Bond_list.txt') 
!            do i=2,np
!            write(101,*)  b_laccoc(i)
!            write(101,*)  i,b_laccot(1:b_laccoc(i),i),  gap_ini(1,i)
!            end do         
!            close(101)  

      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***************************************************************
	subroutine generate_wall_boundary
!***************************************************************
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
	  use system
    use solid
    use fluid
    use fs_inter
    implicit none
!
!.....Initialisation
      do iz=0,nz
        do iy=0, ny
          do ix=0, nx
          node(ix,iy,iz)%obst=0
          end do
        end do
      enddo

! body force driven periodic BCs
    if (bctype.eq.1 .and. periodic_flag) then
          print *, "body foce with periodic BC activated!"
          return
    end if

    if (bctype.eq.3) then

        if (periodic_flag) then
          print *, "Periodic activated!"
          return
        end if
        print *, "NOT Periodic!"

        do iy=0, ny
        do ix=0, nx
          node(ix,iy,nz)%obst=1
          now0=now0+1
          lnw0(1,now0)=ix
          lnw0(2,now0)=iy
          lnw0(3,now0)=nz
        enddo
        enddo
  !.....Set bottom wall boundary
        do iy=0, ny
        do ix=0, nx
          node(ix,iy,0)%obst=1
          now0=now0+1
          lnw0(1,now0)=ix
          lnw0(2,now0)=iy
          lnw0(3,now0)=0
        enddo
        enddo
  !.....Set front wall boundary
        do iz=1, nz-1
        do ix=0, nx
          node(ix,0,iz)%obst=1
          now0=now0+1
          lnw0(1,now0)=ix
          lnw0(2,now0)=0
          lnw0(3,now0)=iz
        enddo
        enddo
  !.....Set back wall boundary
        do iz=1, nz-1
        do ix=0, nx
          node(ix,ny,iz)%obst=1
          now0=now0+1
          lnw0(1,now0)=ix
          lnw0(2,now0)=ny
          lnw0(3,now0)=iz
        enddo
        enddo

    else
      now0=0
!.....Set top wall boundary
      do iy=0, ny
	    do ix=0, nx
   	    node(ix,iy,nz)%obst=1
        now0=now0+1
	      lnw0(1,now0)=ix
	      lnw0(2,now0)=iy
   	    lnw0(3,now0)=nz
      enddo
      enddo
!.....Set bottom wall boundary
      do iy=0, ny
	    do ix=0, nx
   	    node(ix,iy,0)%obst=1
	      now0=now0+1
	      lnw0(1,now0)=ix
	      lnw0(2,now0)=iy
   	    lnw0(3,now0)=0
      enddo
      enddo
!.....Set front wall boundary
      do iz=1, nz-1
	    do ix=0, nx
   	    node(ix,0,iz)%obst=1
	      now0=now0+1
	      lnw0(1,now0)=ix
	      lnw0(2,now0)=0
   	    lnw0(3,now0)=iz
      enddo
      enddo
!.....Set back wall boundary
      do iz=1, nz-1
	    do ix=0, nx
   	    node(ix,ny,iz)%obst=1
	      now0=now0+1
	      lnw0(1,now0)=ix
	      lnw0(2,now0)=ny
        lnw0(3,now0)=iz
      enddo
      enddo
!.....Set left wall boundary
      do iz=1, nz-1
	    do iy=1, ny-1
   	    node(0,iy,iz)%obst=1
	      now0=now0+1
	      lnw0(1,now0)=0
	      lnw0(2,now0)=iy
        lnw0(3,now0)=iz
      enddo
      enddo
!.....Set right wall boundary
      do iz=1, nz-1
	    do iy=1, ny-1
   	    node(nx,iy,iz)%obst=1
	      now0=now0+1
	      lnw0(1,now0)=nx
	      lnw0(2,now0)=iy
        lnw0(3,now0)=iz
      enddo
      enddo  
    endif
!      open(10001,file='wall_nodes_list.txt')      
!      do i=1,now0
!        write(10001,*) i, lnw0(1:3,i)
!       enddo 
!       close(10001)
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine boundary_profile
	use system
    use solid
    use fluid
    use fs_inter
    implicit none
!.....nn0,nn1 are total nodes before current wall and they are used to calculate nodes of current wall
	  integer nn0,nn1
!	  common /wall/nw0(100),nw1(100),pw0(100),pw1(100),lnf(200)
!
!.....Initialisation
      do iz=0,nz
        do iy=0, ny
          do ix=0, nx
          node(ix,iy,iz)%obst=0
          end do
        end do
      enddo
!
	  do i=0,nx
!.......top wall
	    xb(1,i)=ny
!.......bottom wall
	    xb(2,i)=0
	  enddo
!
	  do i=0,ny
!.......right wall
	    yb(1,i)=nx
!.......left wall
	    yb(2,i)=0
	  enddo
!
	  do i=0,nz
!.......front wall
	    zb(1,i)=nx
!.......back wall
	    zb(2,i)=0
	  enddo
!.....Loop over fixed boundary edges/surfaces/walls
    now0=0
	  nn0=0
    k=0
!   

!  in the future deal with walls ????
	  if(ne0>1)then
      p1=>head_le0
      do while(.true.)
        i=p1%num       
        is=wall(i)%point(1)
	      ie=wall(i)%point(2)
        vertex_sw(1,is)=vertex(1,is)
        vertex_sw(2,is)=vertex(2,is)
        vertex_sw(1,ie)=vertex(1,ie)
        vertex_sw(2,ie)=vertex(2,ie)
!	    call edge_nodes(1,nx,ny,is,ie,vertex,xb,yb,now0,lnw0)
!.......number of boundary nodes for current edge
	      k=k+1
        nw0(k)=now0-nn0
!.......pointer pointing to the boundary node list
        pw0(k)=nn0+1
	      nn0=now0
!
        if(.not. associated(p1%next)) exit
        p1=>p1%next       
      end do	
      end if    
!
!write(*,*) now0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        do k=1,ne0
!          write(*,*) nw0(k)
!        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....Loop over non-fixed boundary edges
      now1=0
	  nn1=0
      k=0
!
      if(ne1>0)then
      p1=>head_le1
      do while(.true.)
        i= p1%num
	    is=wall(i)%point(1)
	    ie=wall(i)%point(2)
!	    call edge_nodes(0,nx,ny,is,ie,vertex,xb,yb,now1,lnw1)
!.......number of boundary nodes for current edge
	    k=k+1
        nw1(k)=now1-nn1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        wall(i)%nmwn=nw1(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.......pointer pointing to the boundary node list
	    pw1(k)=nn1+1
	    nn1=now1
!
        if(.not. associated(p1%next)) exit
        p1=>p1%next       
      end do	
      end if
!  
!write(*,*) now1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        do k=1,ne1
!          write(*,*) nw1(k)
!        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!       
      return
      end	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine init_density
    use fluid
    use system, only: use_MRT
    implicit none
!.....Weighting factors (depending on lattice geometry)
      real(8) dd,w2,w1,w0
      integer iy,ix,iz,i
!------------------------------------------------
!.....dd is used for reducing computing time 
      dd=d0
      if(bctype.eq.2)then
	    if(neb.eq.2)then
	      dd=(bc(1)+bc(2))/2.d0
	    else if(neb.eq.3)then
	      dd=(bc(1)+bc(2)+bc(3))/3.d0
	    endif
	  endif
!------------------------------------------------
! Initilize the weights for MRT
    if(use_MRT)then
      ! MRT D3Q15 weights are consistent with (from MRT module):
      ! e_eq       = -rho
      ! epsilon_eq = -rho
      w0 = (2.0d0/15.0d0)  * dd
      w1 = (2.0d0/15.0d0)  * dd
      w2 = (1.0d0/120.0d0) * dd
    else  
  !      w2=dd/36.d0
      w2=dd/72.d0
      w1=8.0d0*w2
      w0=2.d0*w1
    end if
!
!.....Loop over rows
    do iz=-1,nz+1
      do iy=-1,ny+1
!.......Get start and end nodes of the row
!        ixs=yb(1,iy)
!	       ixe=yb(2,iy)
!.......Loop over nodes of current row
        do ix=-1,nx+1
          temp(0,ix,iy,iz)=w0
          do i=1,6
            temp(i,ix,iy,iz)=w1
	        enddo
          do i=7,14
            temp(i,ix,iy,iz)=w2
	        enddo
        enddo
      enddo
    enddo

!   write(*,*) temp(:,10,10)
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine generate_stationary_particles
    use system
    use solid
    use fluid
    use fs_inter
    implicit none  
!.....local variables
      integer  xl,xr,yl,yr,zl,zr,deltax,deltay,deltaz
      real(8) radius2
!
      nos=0
	    p1=>head_lps
      do while(.true.)
        ip= p1%num
	      radius2=particle(ip)%radius*particle(ip)%radius
	      xl=particle(ip)%coor(1)-particle(ip)%radius
	      xr=particle(ip)%coor(1)+particle(ip)%radius
	      yl=particle(ip)%coor(2)-particle(ip)%radius
	      yr=particle(ip)%coor(2)+particle(ip)%radius
        zl=particle(ip)%coor(3)-particle(ip)%radius
	      zr=particle(ip)%coor(3)+particle(ip)%radius
        do iz=zl,zr
          do iy=yl,yr
          do ix=xl,xr
	          deltax=ix-particle(ip)%coor(1)
	          deltay=iy-particle(ip)%coor(2)
            deltaz=iz-particle(ip)%coor(3)
	if((deltax*deltax+deltay*deltay+deltaz*deltaz).le.radius2)then
	          node(ix,iy,iz)%obst=1
	          nos=nos+1
!.............form the list using pointer..............	          
            if(nos.eq.1)then
                allocate(head_lns)
                head_lns%coord(1)=ix
                head_lns%coord(2)=iy
                head_lns%coord(3)=iz
                nullify(head_lns%next)
                p3=>head_lns
            else
                allocate(p3%next,stat=err)
                if(err/=0) then
                  write(*,*) "Out of memory!"
                  stop
                end if
                p3=>p3%next
                p3%coord(1)=ix
                p3%coord(2)=iy
                p3%coord(3)=iz
                nullify(p3%next)
            end if                       
!.......................................................              
	endif
          enddo
          enddo
        end do
!    
        ip= p1%num
        if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.........ouput array to check pointer code
          pointer_check=0
          if(pointer_check.eq.1)then
            p3=>head_lns
            do while(.true.)
              write(*,*) p3%coord(1),p3%coord(2),p3%coord(3)
              if(.not. associated(p3%next)) exit
              p3=>p3%next
            end do
!
          end if      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
    End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine generate_moving_particles
    use system
    use solid
    use fluid
    use fs_inter
    implicit none       
!      integer nnb,nni
!
      noi=0
	    nni=0
	    nob=0
	    nnb=0
      nof=0
      nnf=0
!
      p1=>head_lpm
      do while(.true.)
        ip=p1%num
        if(particle(ip)%active==0)goto 100
        call nodes_of_moving_particles(ip)
!.......number of boundary nodes for current particle
	    particle(ip)%nb=nob-nnb
!        write(*,*) ip,particle(ip)%nb
!.......pointer pointing to the boundary node list
	    particle(ip)%pb=nnb+1
!.......number of interior nodes for current particle
	    particle(ip)%ni=noi-nni
!.......pointer pointing to the interior node list
	    particle(ip)%pn=nni+1
!        write(*,*) ip,particle(ip)%ni,particle(ip)%pn
	    nnb=nob
	    nni=noi
!
100        if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do

!      write(*,*) 	    nni, nob, nnb, nof, nnf
! stop
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine nodes_of_moving_particles(ip)
    use system,only:istep
    use solid
    use fluid,only:nx,ny,nz,lni,lnb,lnf,node
    use fs_inter,only:noi,nob,nof
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)    
!      
!      real(8) xc,yc,rad,rad2
!.....local variables
    integer lb(3,50000),l0(3,14),xl,xr,yf,yb,zu,zd,x,y,z
    integer  ip,direction,current,inside_particle,NNN
    integer ixx,iyy,izz
!
!.....Particle paramters
	  xc=particle(ip)%coor(1)
	  yc=particle(ip)%coor(2)
    zc=particle(ip)%coor(3)
	  rad=particle(ip)%radius
	  rad2=rad*rad
!.....min and max of the outer cubiod of the sphere
	  ixmin=xc-rad
	  dx=xc+rad
	  ixmax=dx
	  if(ixmax*1.d0.lt.dx)ixmax=ixmax+1

	  iymin=yc-rad
	  dy=yc+rad
	  iymax=dy
	  if(iymax*1.d0.lt.dy)iymax=iymax+1

	  izmin=zc-rad
    dz=zc+rad
	  izmax=dz
	  if(izmax*1.d0.lt.dz)izmax=izmax+1

    if(boundary_dem==0)then
!.....modification
      if(ixmin<0) ixmin=0
!.....modification
      if(ixmax>nx) ixmax=nx
!.....modification
      if(iymin<0) iymin=0
!.....modification
      if(iymax>ny) iymax=ny
!.....modification
      if(izmin<0) izmin=0
!.....modification
      if(izmax>nz) izmax=nz
!    elseif(boundary_dem==1)then 
    endif

!.....min and max of the inner cuboid of the sphere
    a=0.577d0*rad
	  ixmi=xc-a
	  iymi=yc-a
	  izmi=zc-a
	  ixma=xc+a
	  iyma=yc+a
	  izma=zc+a
    if(boundary_dem==0)then
      if(ixmi<0) ixmi=0
      if(iymi<0) iymi=0
      if(izmi<0) izmi=0
      if(ixma>nx) ixma=nx
      if(iyma>ny) iyma=ny
      if(izma>nz) izma=nz
    endif
!.....Find the interior nodes
    do iz=izmi+1, izma-1
	    do iy=iymi+1, iyma-1
        do ix=ixmi+1, ixma-1
          !record the interior nodes for virtual particle? 
          ixx=ix
          iyy=iy
          izz=iz
          if(boundary_dem==1)then
            if(ixx>nx) ixx=ix-nx-1
            if(ixx<0) ixx=ix+nx+1
            if(iyy>ny) iyy=iy-ny-1
            if(iyy<0) iyy=iy+ny+1
            if(izz>nz) izz=iz-nz-1
            if(izz<0) izz=iz+nz+1
          endif
          ! mark pure solid nodes including those belonging to virtual particles
          node(ixx,iyy,izz)%obst=3
	        noi=noi+1
	        lni(1,noi)=ixx
	        lni(2,noi)=iyy
   	      lni(3,noi)=izz

        enddo
      enddo
    enddo
!
!.....Rough list of boundary nodes
  nb=0
	do iz=izmin, izmax
	  do iy=iymin, iymax
	    do ix=ixmin, ixmax
        if (boundary_dem==1)then
          ! check node array for PBC but not link them to virtual particle now
          ixx=ix
          iyy=iy
          izz=iz
          if(ix>nx) ixx=ix-nx-1
          if(ix<0) ixx=ix+nx+1
          if(iy>ny) iyy=iy-ny-1
          if(iy<0) iyy=iy+ny+1
          if(iz>nz) izz=iz-nz-1
          if(iz<0) izz=iz+nz+1
          if(node(ixx,iyy,izz)%obst.eq.3) cycle
        else
	        if(node(ix,iy,iz)%obst.eq.3) cycle
        endif
        !
	      if(inside_particle(ix,iy,iz,xc,yc,zc,rad2)==1)then
          if (boundary_dem==1)then
            node(ixx,iyy,izz)%obst=2
          else
	          node(ix,iy,iz)%obst=2
          endif
	        nb=nb+1
	        lb(1,nb)=ix
	        lb(2,nb)=iy
   	      lb(3,nb)=iz
	       endif
      enddo
    enddo
  end do
!      open(10002,file='Particle_nodes_list.txt')      
!      do i=1,nb
!        write(10002,*) i, lb(1:3,i)
!      enddo 
!      close(10002)
! if(istep==2)      stop

   NNN=0
!.....Detail check
  do i=1,nb
	  ix=lb(1,i)
	  iy=lb(2,i)
	  iz=lb(3,i)

    ixx=ix
    iyy=iy
    izz=iz

    if (boundary_dem==1)then
          ! check node array for PBC
          if(ix>nx) ixx=ix-nx-1
          if(ix<0) ixx=ix+nx+1
          if(iy>ny) iyy=iy-ny-1
          if(iy<0) iyy=iy+ny+1
          if(iz>nz) izz=iz-nz-1
          if(iz<0) izz=iz+nz+1
      if(node(ixx,iyy,izz)%obst.eq.4.or.node(ixx,iyy,izz)%obst.eq.3)goto 30
    else
	    if(node(ix,iy,iz)%obst.eq.4.or.node(ix,iy,iz)%obst.eq.3)goto 30
    endif

!	  if(node(ix,iy,iz)%obst.eq.4.or.node(ix,iy,iz)%obst.eq.3)goto 30
!.......neighbour nodes check
      xr=ix+1
      xl=ix-1
      yb=iy+1
      yf=iy-1
	    zu=iz+1
      zd=iz-1
	    n0=0
	    n4=0
      call check_neighbour_nodes(xr,iy,iz,n0,n4,l0,NNN)
      call check_neighbour_nodes(xl,iy,iz,n0,n4,l0,NNN)
      call check_neighbour_nodes(ix,yf,iz,n0,n4,l0,NNN)
      call check_neighbour_nodes(ix,yb,iz,n0,n4,l0,NNN)
      call check_neighbour_nodes(ix,iy,zu,n0,n4,l0,NNN)
      call check_neighbour_nodes(ix,iy,zd,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xr,yb,zu,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xl,yf,zd,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xr,yb,zd,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xl,yf,zu,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xr,yf,zu,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xl,yb,zd,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xr,yf,zd,n0,n4,l0,NNN)
	    call check_neighbour_nodes(xl,yb,zu,n0,n4,l0,NNN)

! finalize solid boundary nodes using lnb array considering PBC/virtual particles
	    if(n0.gt.0)then
      	nob=nob+1
        if(boundary_dem==1)then
          node(ixx,iyy,izz)%obst=2
          lnb(1,nob)=ixx
	        lnb(2,nob)=iyy
	        lnb(3,nob)=izz
        else
	        node(ix,iy,iz)%obst=2
          lnb(1,nob)=ix
	        lnb(2,nob)=iy
	        lnb(3,nob)=iz
        endif
!find fluid boundary nodes
	      do j=1,n0
	        x=l0(1,j)
	        y=l0(2,j)
	        z=l0(3,j)
          if (boundary_dem==1)then
          ! check node belonging to virtual particles 
            if(x>nx) x=x-nx-1
            if(x<0)  x=x+nx+1
            if(y>ny) y=y-ny-1
            if(y<0)  y=y+ny+1
            if(z>nz) z=z-nz-1
            if(z<0)  z=z+nz+1
          end if
	          node(x,y,z)%obst=4
	          nof=nof+1
	          lnf(1,nof)=x
	          lnf(2,nof)=y
	          lnf(3,nof)=z
!          endif
	      enddo
	      goto 30
	    endif

	    if(n4.gt.0)then

	      node(ixx,iyy,izz)%obst=2
	      nob=nob+1
	      lnb(1,nob)=ixx
	      lnb(2,nob)=iyy
	      lnb(3,nob)=izz
  
        goto 30
	    endif
	    node(ixx,iyy,izz)%obst=3
	    noi=noi+1
	    lni(1,noi)=ixx
	    lni(2,noi)=iyy
	    lni(3,noi)=izz
  30 continue 
    enddo
!    write(*,*) nof,nob,noi
!  stop
     
    return
  End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!**********************************************************************
  Subroutine check_neighbour_nodes(ix,iy,iz,n0,n4,l0,NNN)
!**********************************************************************
    use fluid,only:nx,ny,nz,node
    use solid,only:boundary_dem
    implicit none     
    integer  l0(3,14),n0,n4,NNN,ix,iy,iz
    integer ixx,iyy,izz
!     use temporal variable for node obst check
        ixx=ix
        iyy=iy
        izz=iz
      if(boundary_dem==1)then
        if(ix>nx) ixx=ix-nx-1
        if(ix<0) ixx=ix+nx+1
        if(iy>ny) iyy=iy-ny-1
        if(iy<0) iyy=iy+ny+1
        if(iz>nz) izz=iz-nz-1
        if(iz<0) izz=iz+nz+1
      else 
        if(ix<0.or.ix>nx .or. iy<0.or.iy>ny .or. iz<0.or.iz>nz) goto 30
      endif
! check the current neighbour node is fluid (boundary) node or not
      NNN=NNN+1
	    if(node(ixx,iyy,izz)%obst.eq.0)then
	      n0=n0+1
        ! for PBC only record the nodes outside domain, we will correct it later for virtual particles
	      l0(1,n0)=ix
	      l0(2,n0)=iy
	      l0(3,n0)=iz
!          write(*,*) ix,iy,iz,node(ix,iy,iz)%obst
	    else if(node(ixx,iyy,izz)%obst.eq.4)then
	      n4=n4+1
!          write(*,*) ix,iy,iz,node(ix,iy,iz)%obst
	    endif
  30  return

	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine node_incr(direction,i,j)
    implicit none
	  integer direction,i,j
!
	  select case (direction)
!       west
	    case (1,5) 
	      i=-1
	      j=0
!       north
	    case (2,6)
	    i=0
	    j=1
!       east
	    case (3)
	    i=1
	    j=0
!       south
	    case (4)
	      i=0
	      j=-1
	  end select
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Integer function inside_particle(ix,iy,iz,xc,yc,zc,rad2)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
	  inside_particle=1
	  dx=xc-ix
	  dy=yc-iy
    dz=zc-iz
	  if((dx*dx+dy*dy+dz*dz).gt.rad2) inside_particle=0
	  return
	End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine contact_main
    use system
    use solid
    use fs_inter,only:imwp
    implicit none
!.....used for reducing computing time
      integer::MAGNIFY=10000

!.....Initialise contact force
      do i=1,np
	    do j=1,6
	      particle(i)%FC(j)=0.d0
        enddo
	  enddo
!
      sum_disp=sum_disp+vmax*dtime
!      write(*,*) istep
!
!      if(istep.eq.1.or.sum_disp.gt.buff)then
        sum_disp=0.d0
!.......Particle bounding box (with buffer) coordinates
        do i=1,np	  
	      do j=1,3
	        icoor(j,i)  =(particle(i)%coor(j)-particle(i)%radius-buff)*MAGNIFY
	        icoor(j+3,i)=(particle(i)%coor(j)+particle(i)%radius+buff)*MAGNIFY
	      enddo
	    enddo
!.......Perform global search with NBS algorithm
	    call nbsw3d(icoor,laccoc,laccop,laccot,mcnts,3,6,np) 
!	  endif
!
!.....Compute contact force
      call contact_force
!
      if(bond) call bond_force
!            write(*,*) "i LOVE mIN"
!.....Boundary treatment
      if(boundary_dem==0)then
        call wall_particle 
      elseif(boundary_dem==1)then
        call periodic_boundary
      elseif(boundary_dem==2)then
        call wall_particle_special
      elseif(boundary_dem==3)then
!       call triaxial_boundary
      else

      end if   
!
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine contact_force
    use system
    use solid
    use fluid,only:dx
    use fs_inter,only:imwp
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  real(8) dn(3)
      integer kkk
! ....temporal array used for tangential force calculation
      real(8),allocatable::save_ft(:,:,:)
!      
      allocate(save_ft(0:3,np,np))
      kkk=0
!.....save_ft=0
      do i=1,np
        do j=1,np
          do k=0,3
            save_ft(k,j,i)=0
          enddo
        enddo
      enddo 
!     Assume a small thichness of a fluid compressive layer
      if(lubrication)then
	    ad=0.1d0*dx
      else
        ad=0.0d0
      end if
!
!.....Particle-Particle contact
      loop1: do i=1,np
	    if(particle(i)%active.eq.0) cycle loop1
	    ntar=laccoc(i)
        
	    if(ntar<=0) cycle loop1
	    x1=particle(i)%coor(1)
	    y1=particle(i)%coor(2)
        z1=particle(i)%coor(3)
	    r1=particle(i)%radius
	    ip=laccop(i)-1
!
	    loop2: do it=1,ntar
	      j=laccot(ip+it)
	      if(j.ge.i.or.particle(j)%active.eq.0) cycle loop2
!---------Check if j is in bond list-----------
          if(bond)then
!            nnn=b_laccoc(i)
            b_exist=0
            loop3: do ik=1,b_laccoc(i)
              kkk=b_laccot(ik,i)
              if(j-kkk==0)then
                b_exist=1
                exit loop3
              end if
            end do loop3
!...
            if(b_exist==1) cycle loop2
          end if
!----------------------------------------------

	      x2=particle(j)%coor(1)
	      y2=particle(j)%coor(2)
          z2=particle(j)%coor(3)
	      r2=particle(j)%radius
	      deltax=x1-x2
	      deltay=y1-y2
          deltaz=z1-z2
	      rs=r1+r2
!.........normal gap
	      gap=deltax*deltax+deltay*deltay+deltaz*deltaz
!. if no contact initionalise previous tangential force==0
	      if(gap.ge.(rs+ad)*(rs+ad)) cycle loop2
	      dnorm=sqrt(gap)
	      gap=dnorm-rs-ad
!.........contact direction
	      dn(1)=deltax/dnorm
	      dn(2)=deltay/dnorm
          dn(3)=deltaz/dnorm
!.........normal contact force
          fn=kn*gap
!.........viscous damping
          fnd=fn
          if(damp_option.eq.1)then
	        veln1=particle(i)%U(1)*dn(1)+particle(i)%U(2)*dn(2)+particle(i)%U(3)*dn(3)
	        veln2=particle(j)%U(1)*dn(1)+particle(j)%U(2)*dn(2)+particle(j)%U(3)*dn(3)
	        velnr=veln1-veln2
	        if(abs(velnr).gt.0.d0)then
		      fmin=sqrt(particle(i)%mass*particle(j)%mass*kn/(particle(i)%mass+particle(j)%mass))
		      fnd=fn+2.d0*damp_coef*velnr*fmin
	        endif
	      endif	
!.........assemble the normal component of contact force
          do k=1,3
	        comp=fnd*dn(k)
	        particle(i)%FC(k)=particle(i)%FC(k)-comp
	        particle(j)%FC(k)=particle(j)%FC(k)+comp
	      enddo
!--------------------------------------Tangential Force Calculation---------------------------------------
        if(fric_option==1) then
				
!.........Extract tangential force at previous step
          do k=0,3
            ft(k)=pre_ft(k,j,i)
          enddo
!            
!.........Compute tangential force
		  call tangential_force_calculation(i,j,damp_coef,damp_option,dn,r1,r2,fn,ft,dtime,kt,fric_coef) 
!    
!........ Save current tangential force for next timestep
           do k=0,3
             save_ft(k,j,i)=ft(k)
           end do                     
!...				
        end if
!-----------------------------------------------------------------------------------------------------------
        end do loop2
      end do loop1
!
!..............................................................................
!......... Store current tangential force as previous for next calculation
        do i=1,np
          do j=1,np
            do k=0,3
              pre_ft(k,j,i)=save_ft(k,j,i)
            enddo
          enddo
        enddo        
        deallocate(save_ft)	
!..............................................................................         			

      return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine tangential_force_calculation(i,j,damping_coef,damp_option,dn,rc,rt,fn,ft,dtime,penalty,fric_coef)                                                            
    use solid,only:particle
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer i,j,k,damp_option
    real(8) dn(3),rc,rt,dtime
    real(8) penalty,fric_coef,damping_coef,mass1,mass2,vc(6),&
             vt(6),vcc(3),vct(3),vr(3),vrn,vrt(3),vrrt,&
             ftd(0:3),dt(3),orn,gapt,ftn,ftt,fto(0:3),ft1,ft2,dd,sita,&
             ft11,ft22,ftmax,fmin,fm(3),ft(0:3)
!
!.....Extract variables
      mass1=particle(i)%mass
      mass2=particle(j)%mass
!  in the future deal with the abgular velocity here???????      
      do k=1,6
        vc(k)=particle(i)%U(k)
        vt(k)=particle(i)%U(k)
      end do 
!      vc(4:6)=particle(i)%av
!      vt(4:6)=particle(j)%av    
!.....Components of velocities at contact point for contactor and target
      vcc(1)=vc(1)+(vc(5)*dn(3)-vc(6)*dn(2))*rc
      vcc(2)=vc(2)+(vc(6)*dn(1)-vc(4)*dn(3))*rc
      vcc(3)=vc(2)+(vc(3)*dn(2)-vc(5)*dn(1))*rc
      vct(1)=vt(1)-(vt(5)*dn(3)-vt(6)*dn(2))*rt
      vct(2)=vt(2)-(vt(6)*dn(1)-vt(4)*dn(3))*rt
      vct(3)=vt(2)-(vt(3)*dn(1)-vt(4)*dn(1))*rt
!.....Componenets of relative velocities
      do k=1,3
        vr(k)=vct(k)-vcc(k)
      end do
!.... Projection of relative velocity on the normal direction
      vrn=0.00
      do k=1,3
        vrn=vrn+vr(k)*dn(k)
      end do
!.... Relative tangential velocities
      vrrt=0.00
      do k=1,3
        vrt(k)=vr(k)-vrn*dn(k)
        vrrt=vrrt+vrt(k)*vrt(k)
      end do
        vrrt=sqrt(vrrt)
!......Get previous tangential force  0:total vector    1,3: components 
      do k=0,3
        ftd(k)=ft(k)
      end do
!...         
      if(ABS(vrrt).gt. 0.00) then
!.....Direction cosines of relative tangential velocity
	    do k=1,3
          dt(k)=vrt(k)/vrrt
        end do
!...... Sliding distance
        gapt=vrrt*dtime
          
!.......Co-rotate ftold to current tangential plane
!..... 
!.......Projection of old tangential force on normal direction
        ftn=ft(1)*dn(1)+ft(2)*dn(2)+ft(3)*dn(3)
!.......projection of old tangential force on tangential plane
        ftt=ft(0)*ft(0)-ftn*ftn
!          write(*,*)  ftt        
!..
          if(ftt <=0)then 
            ftt=0.0d0
            ft11=0.0d0
            ft22=0.0d0
          else
            ftt=sqrt(ftt)
!..... Percentage of ftt to ftold
            dd=ft(0)/ftt
!..... Components of ftt vector (enlagrged by dd to make |ftt|=|ftold|)
			  do k=1,3
                fto(k)=dd*(ft(k)-ftn*dn(k))
              end do
              
!..... Decompose ftt in two orthogonal directions (nt and nt')
              ft1=fto(1)*dt(1)+fto(2)*dt(2)+fto(3)*dt(3)
              ft2=ft(0)*ft(0)-ft1*ft1
!..
              if(ft2 > 0)then
                ft2=sqrt(ft2)
              else 
                ft2=0.0d0
              end if
!..
!..... Relative rotational angle of previous and current tangential planes			   
              orn=0.0d0
              do k=1,3
              orn=orn+(vt(k+3)-vc(k+3))*dn(k) 
              end do        
!..                                                                       
			  sita=orn*dtime
!...
!...... Rotate ftt to current local coordinates (normal and tangential   sita~~sine(sita))
              ft11=ft1+ft2*sita
              ft22=-ft1*sita+ft2
          end if
!                      write(*,*)ft11,ft22
!........ Trial tangential fore	  	    
          ft(0)=ft11+penalty*gapt    
          ft(0)=ft(0)*ft(0)+ft22*ft22
          if(ft(0)> 0.0) ft(0)=sqrt(ft(0))
    
!....... Check limiting friction force
		  ftmax=ABS(fric_coef*fn)
		  if(ft(0)>ftmax)then
            ft(0)=ftmax
		  else if(ft(0) < -ftmax)then
            ft(0)=-ftmax 
          end if
           
!...... Components of tangential force
		  do k=1,3
            ft(k)=ft(0)*dt(k)
          enddo
!...
!..... Damping tangetial component of contact force
		  do k=0,3
            ftd(k)=ft(k)
          enddo
!....
		  if(damping_option==1) then 		
			 fmin=sqrt(mass1*mass2*penalty/(mass1+mass2))
			 ftd(0)=ft(0)- 2.0d0*damping_coef*vrrt*fmin
			 if(ftd(0) > ftmax)then
               ftd(0)=ftmax
			 else if(ftd(0) < -ftmax)then
               ftd(0)=-ftmax
             end if
		   end if
!...
		   do k=1,3
             ftd(k)=ftd(0)*dt(k)
           enddo
!.
!..... Assemble tangential component at the contact point to calculate moment 
		   fm(1)=ftd(3)*dn(2)-ftd(2)*dn(3)
           fm(2)=ftd(1)*dn(3)-ftd(3)*dn(1)
           fm(3)=ftd(2)*dn(1)-ftd(1)*dn(2)
!           write(*,*) fm,ftd(1)
    
!..... Set contact force on contactor
		   do k=1,3
           particle(i)%FC(k)=particle(i)%FC(k)-ftd(k)
           particle(i)%FC(k+3)=particle(i)%FC(k+3)-rc*fm(k)
           end do
!.... Set contact force on target		   
		   do k=1,3
           particle(j)%FC(k)=particle(j)%FC(k)+ftd(k)
           particle(j)%FC(k+3)=particle(j)%FC(k+3)-rc*fm(k)
           end do
        end if
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine bond_force
    use system
    use solid
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  real(8) dn(2)
! ....temporal array used for tangential force calculation
      real(8),allocatable::b_save_ft(:,:,:)
!      
      allocate(b_save_ft(0:2,np,np))
!.....save_ft=0
      do i=1,np
        do j=1,np
          do k=0,2
            b_save_ft(k,j,i)=0
          enddo
        enddo
      enddo 
!    
!     calculate the critical gap (negative overlap) corresponding to critical cohesive forces
	  b_gap=cb_cnf/cb_kn
!      write(*,*) b_gap
!      ad=0.125
!.....Particle-Particle contact using bond model
      do 100 i=2,np
	    if(particle(i)%active.eq.0)goto 100
	    ntar=b_laccoc(i)
	    if(ntar.eq.0)goto 100
	    x1=particle(i)%coor(1)
	    y1=particle(i)%coor(2)
	    r1=particle(i)%radius
!
	    do 50 it=1,ntar
	      j=b_laccot(it,i)
	      if(j.eq.0 )goto 50
          if( j.ge.i .or. particle(j)%active.eq.0)goto 50
!
!
	      x2=particle(j)%coor(1)
	      y2=particle(j)%coor(2)
	      r2=particle(j)%radius
	      deltax=x1-x2
	      deltay=y1-y2
	      rs=r1+r2
!.........normal gap
	      gap=deltax*deltax+deltay*deltay
!-------------------------Damage criterion-------------------------------
	      if(gap-(rs+b_gap+gap_ini(j,i))*(rs+b_gap+gap_ini(j,i))>0.0d0)then
              b_laccot(it,i)=0
              goto 50
	      end if
!------------------------------------------------------------------------
          dnorm=sqrt(gap)	      
          gap=dnorm-rs
!.........contact direction
	      dn(1)=deltax/dnorm
	      dn(2)=deltay/dnorm
!.........normal cohesion force
          if (sam_prep2) then
            fn=cb_kn*(gap-gap_ini(j,i))
!            write(*,*) fn
          else
            fn=cb_kn*gap
          end if
!.........viscous damping
          fnd=fn
!          
          if(damp_option.eq.1)then
	        veln1=particle(i)%U(1)*dn(1)+particle(i)%U(2)*dn(2)
	        veln2=particle(j)%U(1)*dn(1)+particle(j)%U(2)*dn(2)
	        velnr=veln1-veln2;
	        if(abs(velnr).gt.0.d0)then
		      fmin=sqrt(particle(i)%mass*particle(j)%mass*cb_kn/(particle(i)%mass+particle(j)%mass))
		      fnd=fn+2.d0*damp_coef*velnr*fmin
	        endif
	      endif	
!.........assemble the normal component of cohesion force
          do k=1,2
	        comp=fnd*dn(k)
!            write(*,*) comp
	        particle(i)%FC(k)=particle(i)%FC(k)-comp
	        particle(j)%FC(k)=particle(j)%FC(k)+comp
	      enddo
!          write(*,*) i,j
!
!----------------------Tangential Force Calculation--------------------------------
        if(fric_option==1) then
				
!......... Extract tangential force at previous step
            do k=0,2
              b_ft(k)=b_pre_ft(k,j,i)
            enddo
!            write(*,*) ft(0)
!......... Compute tangential force..................
			call tangential_force_calculation_bond(i,j,damp_coef,damp_option,dn,r1,r2,fn,b_ft,dtime,b_kt,b_fric_coef,it)
!   write(*,*) b_kt,b_fric_coef,b_ft

!....... 
      
!........ Save current tangential force for next timestep
           do k=0,2
             b_save_ft(k,j,i)=b_ft(k)
           end do                     
!...				
        end if
!------------------------------------------------------------------------------------
!
   50   continue
  100 continue  
!------------------------------------------------------------------------------------
!......... Store current tangential force as previous for next calculation
        do i=1,np
          do j=1,np
            do k=0,2
              b_pre_ft(k,j,i)=b_save_ft(k,j,i)
            enddo
          enddo
        enddo        
        deallocate(b_save_ft)	
!------------------------------------------------------------------------------------     
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine tangential_force_calculation_bond(i,j,damping_coef,damp_option,dn,rc,rt,fn,ft,dtime,penalty,fric_coef,it)   
    use solid,only:particle,b_laccot
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer i,j,k,damp_option,it
      real(8) dn(2),rc,rt,dtime
      real(8) penalty,fric_coef,damping_coef,mass1,mass2,vc(3),&
               vt(3),vcc(2),vct(2),vr(2),vrn,vrt(2),vrrt,&
               ftd(0:2),dt(2),orn,gapt,ftn,ftt,fto(0:2),ft1,ft2,dd,sita,&
               ft11,ft22,ftmax,fmin,fm,ft(0:2)
!........... Extract variables
        mass1=particle(i)%mass
        mass2=particle(j)%mass
!        write(*,*)  fric_coef,penalty,dtime
        do k=1,3
          vc(k)=particle(i)%U(k)
          vt(k)=particle(j)%U(k)
        end do 
!        vc(3)=particle(i)%av
!        vt(3)=particle(j)%av    
!..... Components of velocities at contact point for contactor and target
        vcc(1)=vc(1)-vc(3)*dn(2)*rc
        vcc(2)=vc(2)+vc(3)*dn(1)*rc
        vct(1)=vt(1)+vt(3)*dn(2)*rt
        vct(2)=vt(2)-vt(3)*dn(1)*rt
      
!..... Componenets of relative velocities
        do k=1,2
          vr(k)=vct(k)-vcc(k)
        end do
!.... Projection of relative velocity on the normal direction
        vrn=0.00
        do k=1,2
          vrn=vrn+vr(k)*dn(k)
        end do
!.... Relative tangential velocities
        vrrt=0.00
        do k=1,2
          vrt(k)=vr(k)-vrn*dn(k)
		  vrrt=vrrt+vrt(k)*vrt(k)
        end do
        vrrt=sqrt(vrrt)
!......Get previous tangential force  0:total vector    1,2: components 
        do k=0,2
          ftd(k)=ft(k)
        end do
!...         
        if(ABS(vrrt).gt. 0.00) then
!......Direction cosines of relative tangential velocity
		  do k=1,2
            dt(k)=vrt(k)/vrrt
          end do
!..... Sliding distance
          gapt=vrrt*dtime
          
!..... Co-rotate ftold to current tangential plane
!..... 
!..... Projection of old tangential force on normal direction
          ftn=ft(1)*dn(1)+ft(2)*dn(2)
!..... projection of old tangential force on tangential plane
          ftt=ft(0)*ft(0)-ftn*ftn
!          write(*,*)  ftt        
!..
            if(ftt <=0)then 
              ftt=0.0d0
              ft11=0.0d0
              ft22=0.0d0
            else
              ftt=sqrt(ftt)
!..... Percentage of ftt to ftold
              dd=ft(0)/ftt
!..... Components of ftt vector (enlagrged by dd to make |ftt|=|ftold|)
			  do k=1,2
                fto(k)=dd*(ft(k)-ftn*dn(k))
              end do
              
!..... Decompose ftt in two orthogonal directions (nt and nt')
              ft1=fto(1)*dt(1)+fto(2)*dt(2)
              ft2=ft(0)*ft(0)-ft1*ft1
!..
              if(ft2 > 0)then
                ft2=sqrt(ft2)
              else 
                ft2=0.0d0
              end if
!..
!..... Relative rotational angle of previous and current tangential planes			   
              orn=vt(3)-vc(3)            
!..                                                                       
			  sita=orn*dtime
!...
!...... Rotate ftt to current local coordinates (normal and tangential   sita~~sine(sita))
              ft11=ft1+ft2*sita
              ft22=-ft1*sita+ft2
            end if
!                      write(*,*)ft11,ft22
!........ Trial tangential fore
	  	    
          ft(0)=ft11+penalty*gapt    
          ft(0)=ft(0)*ft(0)+ft22*ft22
          if(ft(0)> 0) ft(0)=sqrt(ft(0))
    
!....... Check limiting friction force
		  ftmax=ABS(fric_coef*fn)
		  if(ft(0)>ftmax)then
            ft(0)=ftmax
		  else if(ft(0) < -ftmax)then
            ft(0)=-ftmax 
          end if
!
!---------------------------Damage criterion-------------------------------
!	      if(ABS(ft(0))>=ftmax)then
!              b_laccot(it,i)=0
!              goto 50
!	      end if
!--------------------------------------------------------------------------
!       write(*,*) ft(0),b_ftmax
!           
!...... Components of tangential force
		  do k=1,2
            ft(k)=ft(0)*dt(k)
          enddo
!...
!.
!..... Assemble tangential component at the contact point to calculate moment 
		   fm=ft(1)*dn(2)-ft(2)*dn(1)
    
!..... Set contact force on contactor
		   particle(i)%FC(1)=particle(i)%FC(1)-ft(1)
           particle(i)%FC(2)=particle(i)%FC(2)-ft(2)
           particle(i)%FC(3)=particle(i)%FC(3)-rc*fm

!.... Set contact force on target		   
           particle(j)%FC(1)=particle(j)%FC(1)+ft(1)
           particle(j)%FC(2)=particle(j)%FC(2)+ft(2)
           particle(j)%FC(3)=particle(j)%FC(3)-rt*fm 
        end if
 50     continue
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine wall_particle
    use system
    use solid
    use fluid,only:dx
    use fs_inter,only:imwp
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!.....wall - partcile contact
!
!write(*,*) "Wall subroutine is entered=="
!.....Fixed boundary edge - partcile contact
      do i=1,nm
        xc=particle(i)%coor(1)
	      yc=particle(i)%coor(2)
        zc=particle(i)%coor(3)
	      rad=particle(i)%radius
!.......YZ walls
        if(xc-rad<0.d0 .or. xc+rad>xmax)then
          if(xc-rad<0.d0)then
!.........calculate overlap
            gap=rad-xc
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(1)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn-2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(1)=particle(i)%FC(1)+fnd
          else
!.........calculate overlap
            gap=rad+xc-xmax
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(1)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn+2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(1)=particle(i)%FC(1)-fnd
          end if
        end if
!.......XZ walls
        if(yc-rad<0.d0 .or. yc+rad>ymax)then
!        write(*,*) "Calculate force=="
          if(yc-rad<0.d0)then
!.........calculate overlap
            gap=rad-yc
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(2)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn-2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(2)=particle(i)%FC(2)+fnd
!===            write(*,*) "Calculate force==",particle(i)%FC(2),kn,fn,2.d0*damp_coef*velnr*fmin
          else
!.........calculate overlap
            gap=rad+yc-ymax
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(2)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn+2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(2)=particle(i)%FC(2)-fnd
          end if
        end if 
!.......XY walls
        if(zc-rad<0.d0 .or. zc+rad>zmax)then
          if(zc-rad<0.d0)then
!.........calculate overlap
            gap=rad-zc
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(3)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn-2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(3)=particle(i)%FC(3)+fnd
          else
!.........calculate overlap
            gap=rad+zc-zmax
!.........normal contact force
            fn=kn*gap
!.........viscous damping
            fnd=fn
            if(damp_option.eq.1)then
	          velnr=particle(i)%U(3)
	          if(abs(velnr).gt.0.d0)then
		        fmin=1.0
		        fnd=fn+2.d0*damp_coef*velnr*fmin
	          endif
	        endif	
!.........assemble the normal component of contact force
            particle(i)%FC(3)=particle(i)%FC(3)-fnd
          end if
        end if        

      end do
!

      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine periodic_boundary
    use system
    use solid
    use fluid,only:dx
    use fs_inter,only:imwp
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      real(8) dn(3)
! ....temporal array used for tangential force calculation
      real(8),allocatable::save_ft(:,:,:)
!      
	  mmm=0
      if(mmm.eq.1)then
      allocate(save_ft(0:2,np,np))
!.....save_ft=0
      do i=1,np
        do j=1,np
          do k=0,2
            save_ft(k,j,i)=0
          enddo
        enddo
      enddo 
      end if
!.....wall - partcile contact
!
!.....Fixed boundary edge - partcile contact
      if(ne0>1)then
      p1=>head_le0
!
      do while(.true.)
        i=p1%num
!      do 200 ii=np+1,np+ne0
	    ntar=laccoc(np+i)
	    if(ntar.eq.0)goto 200
!	    i=p1%num
	    is=wall(i)%point(1)
	    ie=wall(i)%point(2)
	    x1=vertex(1,is)
	    y1=vertex(2,is)
	    x2=vertex(1,ie)
	    y2=vertex(2,ie)
!.......line parameters
        ac=y2-y1
        bc=x1-x2 
	    dlen=sqrt(ac*ac+bc*bc)
	    ac=ac/dlen
	    bc=bc/dlen
        cc=-(ac*x1+bc*y1)
!.......direction cosines
        cosa=ac
        sina=bc
!.......Loop over targets
	    ip=laccop(i+np)-1
	    do 150 it=1,ntar
	      j=laccot(ip+it)
	      if(j.gt.np)goto 150
	      xc=particle(j)%coor(1)
	      yc=particle(j)%coor(2)
	      rad=particle(j)%radius
!.........check if the centre of the particle projected outside of the segment
	      dist=-bc*(xc-x1)+ac*(yc-y1)
	      if(dist.gt.(dlen+rad).or.dist.lt.-rad)goto 150
!.........square distance between disk centre and segment line
	      if(dist.gt.dlen)then
	        if((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2).gt.rad*rad)goto 150
	      else if(dist.lt.0.0)then 
	        if((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1).gt.rad*rad)goto 150
	      endif
!.........distance from the centre to the segment
          dist=-(ac*xc+bc*yc+cc)
          if(dist.le.0.d0.or.dist.ge.rad)goto 150
!.........contact gap (penetration)
          gapnm=dist-rad
	      if(gapnm.ge.0.d0)goto 150
!.........normal contact force
          fn=-kn*gapnm
!.........assemble the normal component of contact force
          particle(j)%FC(1)=particle(j)%FC(1)-fn*cosa
	      particle(j)%FC(2)=particle(j)%FC(2)-fn*sina
      150   continue

      200 continue
!
        if(.not. associated(p1%next)) exit
        p1=>p1%next  
!    
        end do	
      end if 
!
      ad=0.0d0
!.....periodic boundary for right

!.....search particles exceeding left domain
      loop1:do i=1,np
        if(particle(i)%coor(1).lt.(2*maxrad))then
          x1=particle(i)%coor(1)+vertex(1,2)
          y1=particle(i)%coor(2)
          r1=particle(i)%radius
!.....search particles on the right
          loop2:do j=1,np
            if(particle(j)%coor(1).gt.(vertex(1,2)-2*maxrad))then
              	      x2=particle(j)%coor(1)
	                  y2=particle(j)%coor(2)
	                  r2=particle(j)%radius
	                  deltax=x1-x2
	                  deltay=y1-y2
	                  rs=r1+r2
!.........normal gap
	                  gap=deltax*deltax+deltay*deltay
!. if no contact initionalise previous tangential force==0
	                  if(gap.ge.(rs)*(rs)) cycle loop2
	                  dnorm=sqrt(gap)
	                  gap=dnorm-rs
!.........contact direction
	                  dn(1)=deltax/dnorm
	                  dn(2)=deltay/dnorm
!.........normal contact force
                      fn=kn*gap
!.........viscous damping
                      fnd=fn
                      if(damp_option.eq.1)then
	                    veln1=particle(i)%U(1)*dn(1)+particle(i)%U(2)*dn(2)
	                    veln2=particle(j)%U(1)*dn(1)+particle(j)%U(2)*dn(2)
	                    velnr=veln1-veln2
	                    if(abs(velnr).gt.0.d0)then
		                  fmin=sqrt(particle(i)%mass*particle(j)%mass*kn/(particle(i)%mass+particle(j)%mass))
		                  fnd=fn+2.d0*damp_coef*velnr*fmin
	                    endif
	                  endif	
!.........assemble the normal component of contact force
                      do k=1,2
	                    comp=fnd*dn(k)
	                    particle(i)%FC(k)=particle(i)%FC(k)-comp
	                    particle(j)%FC(k)=particle(j)%FC(k)+comp
	                  enddo
!--------------------------------------Tangential Force Calculation---------------------------------------
!
                      if(mmm.eq.1)then
                      if(fric_option==1) then
				
!.........Extract tangential force at previous step
                        do k=0,2
                          if(i>j)then
                          ft(k)=pre_ft1(k,j,i)
                          else
                          ft(k)=pre_ft1(k,i,j)
                          end if
                        enddo
!            
!.........Compute tangential force
		                call tangential_force_calculation(i,j,damp_coef,damp_option,dn,r1,r2,fn,ft,dtime,kt,fric_coef) 
!    
!........ Save current tangential force for next timestep
                        do k=0,2
                          if(i>j)then
                          save_ft(k,j,i)=ft(k)
                          else
                          save_ft(k,i,j)=ft(k)
                          end if
                        end do                     
!...				
                      end if
                      end if 
            end if
          end do loop2
        end if
      end do loop1
!-----------------------------------------------------------------------------------------------------------
!..............................................................................
!......... Store current tangential force as previous for next calculation
        if(mmm.eq.1)then
        do i=1,np
          do j=1,np
            do k=0,2
              pre_ft1(k,j,i)=save_ft(k,j,i)
            enddo
          enddo
        enddo        
        deallocate(save_ft)
        end if	
!..............................................................................
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine wall_particle_special
    use system
    use solid
    use fluid,only:dx
    use fs_inter,only:imwp
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!.....wall - partcile contact
!
!.....Fixed boundary edge - partcile contact
      if(ne0>1)then
      p1=>head_le0
!
      do while(.true.)
        i=p1%num
!      do 200 ii=np+1,np+ne0
	    ntar=laccoc(np+i)
	    if(ntar.eq.0)goto 200
!	    i=p1%num
	    is=wall(i)%point(1)
	    ie=wall(i)%point(2)
	    x1=vertex(1,is)
	    y1=vertex(2,is)
	    x2=vertex(1,ie)
	    y2=vertex(2,ie)
!.......line parameters
        ac=y2-y1
        bc=x1-x2 
	    dlen=sqrt(ac*ac+bc*bc)
	    ac=ac/dlen
	    bc=bc/dlen
        cc=-(ac*x1+bc*y1)
!.......direction cosines
        cosa=ac
        sina=bc
!.......Loop over targets
	    ip=laccop(i+np)-1
	    do 150 it=1,ntar
	      j=laccot(ip+it)
	      if(j.gt.np)goto 150
	      xc=particle(j)%coor(1)
	      yc=particle(j)%coor(2)
	      rad=particle(j)%radius
!.........check if the centre of the particle projected outside of the segment
	      dist=-bc*(xc-x1)+ac*(yc-y1)
	      if(dist.gt.(dlen+rad).or.dist.lt.-rad)goto 150
!.........square distance between disk centre and segment line
	      if(dist.gt.dlen)then
	        if((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2).gt.rad*rad)goto 150
	      else if(dist.lt.0.0)then 
	        if((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1).gt.rad*rad)goto 150
	      endif
!.........distance from the centre to the segment
          dist=-(ac*xc+bc*yc+cc)
          if(dist.le.0.d0.or.dist.ge.rad)goto 150
!.........contact gap (penetration)
          gapnm=dist-rad
	      if(gapnm.ge.0.d0)goto 150
!.........normal contact force
          fn=-kn*gapnm
!.........assemble the normal component of contact force
          particle(j)%FC(1)=particle(j)%FC(1)-fn*cosa
	      particle(j)%FC(2)=particle(j)%FC(2)-fn*sina
      150   continue

      200 continue
!
        if(.not. associated(p1%next)) exit
        p1=>p1%next  
!    
      end do	
    end if 
!
!.....left wall
    do i=1,np
      dist1=particle(i)%coor(1)
      if(dist1<particle(i)%radius)then
        fn=kn*(dist1-particle(i)%radius)	        
	    particle(i)%FC(1)=particle(i)%FC(1)-fn
      end if
    end do

      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer function inactive_particle()
	use system
    use solid,only:particle,boundary_dem
    use fluid,only:ny,nx
    implicit none
!.....number of inactive particle
      integer nia
!.....Loop over moving particles
      nia=0
!
      if(boundary_dem.ne.1)then
      p1=>head_lpm
      do while(.true.)
        ip=p1%num
	    if(particle(ip)%coor(2).gt.(ny) .or. particle(ip)%coor(1).gt. nx)then
	      particle(ip)%active=0
	      nia=nia+1
	    endif


!        write(*,*)  ip,particle(ip)%nf, particle(ip)%pf
        if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do 
      end if
!
      inactive_particle=nia
	  return
	End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine update_moving_particle_position(nm,dt,ga,vmax,isub)
    use system,only:p1,head_lpm,vertex,xmax,ymax,zmax,fix_packing,LBM
    use solid,only:particle,boundary_dem,part_vir,num_vir
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    integer isub
	  real(8) G(3),lx,ly,lz
!
	  G(1)=0.d0
!    
	  G(3)=0.d0
	  G(2)=-ga
	  vmax=0.0d0
!...initialization for virtual particles with periodic boundary 
    if(boundary_dem==1)then
      num_vir=0
      do i=1,nm
        part_vir(i)%active=0
      end do
    end if
!
!.....Loop over moving particles
      p1=>head_lpm
      do while(.true.)
        ip=p1%num      
!
	    if(particle(ip)%active.eq.0) goto 100
!      write(*,*) "particle force when updating motion: ",particle(ip)%F(2),particle(ip)%FC(1:3)
!.......Update particle velocities
        do j=1,3
          particle(ip)%U(j)=particle(ip)%U(j)+(G(j)+(particle(ip)%FC(j)+particle(ip)%F(j))/particle(ip)%mass)*dt
	        if(particle(ip)%U(j).gt.vmax) vmax=particle(ip)%U(j)
          if(fix_packing) particle(ip)%U(j)=0.0d0
        end do
!        write(*,*) "Body force is: ", particle(ip)%mass*G(2)
!.....................................................................................................
!.......update angle velocity
!        particle(ip)%av=particle(ip)%av+((particle(ip)%FC(3)+particle(ip)%F(3))/particle(ip)%moi)*dt
!.......update angle of rotation
!        particle(ip)%aor=particle(ip)%aor+particle(ip)%av*dt
!        write(*,*) ip,particle(ip)%av,particle(ip)%aor
!         write(*,*) particle(ip)%FC(3)
!.....................................................................................................
!.......Update particle position
          particle(ip)%coor(1)=particle(ip)%coor(1)+particle(ip)%U(1)*dt
          particle(ip)%coor(2)=particle(ip)%coor(2)+particle(ip)%U(2)*dt
          particle(ip)%coor(3)=particle(ip)%coor(3)+particle(ip)%U(3)*dt

! periodic BCs: 
          if(boundary_dem==1)then
          ! update physical position
          lx=xmax
          ly=ymax
          lz=zmax
          if(LBM)then
          !  the lx in LBM depends on the treatment of PBC used
            lx=nx+1
            ly=ny+1
            lz=nz+1
          endif
            ! X direction
            if(particle(ip)%coor(1)>xmax) particle(ip)%coor(1)=particle(ip)%coor(1)-lx
            if(particle(ip)%coor(1)<0)  particle(ip)%coor(1)=particle(ip)%coor(1)+lx
            ! Y direction
            if(particle(ip)%coor(2)>ymax) particle(ip)%coor(2)=particle(ip)%coor(2)-ly
            if(particle(ip)%coor(2)<0)  particle(ip)%coor(2)=particle(ip)%coor(2)+ly
            ! Z direction
            if(particle(ip)%coor(3)>zmax) particle(ip)%coor(3)=particle(ip)%coor(3)-lz
            if(particle(ip)%coor(3)<0)  particle(ip)%coor(3)=particle(ip)%coor(3)+lz

          ! generate virtual particles
            if((particle(ip)%coor(1)+particle(ip)%radius)>xmax)then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(1)=particle(ip)%coor(1)-lx
            else if( (particle(ip)%coor(1)-particle(ip)%radius)<0.0d0 )then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(1)=particle(ip)%coor(1)+lx
            endif
!
            if((particle(ip)%coor(2)+particle(ip)%radius)>ymax)then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(2)=particle(ip)%coor(2)-ly
            else if( (particle(ip)%coor(2)-particle(ip)%radius)<0.0d0 )then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(2)=particle(ip)%coor(2)+ly
            endif
!
            if((particle(ip)%coor(3)+particle(ip)%radius)>zmax)then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(3)=particle(ip)%coor(3)-lz
            else if( (particle(ip)%coor(3)-particle(ip)%radius)<0.0d0 )then
              if(part_vir(ip)%active==0) call generate_virtual_particles(ip)
              part_vir(ip)%coor(3)=particle(ip)%coor(3)+lz
            endif

          end if          
!        enddo   
 100       if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do 
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!====================
              Subroutine  generate_virtual_particles(ip)
                use solid,only:particle,boundary_dem,part_vir,num_vir
                implicit none
                integer ip,i
                part_vir(ip)%active=1
                num_vir=num_vir+1
                do i=1,3
                  part_vir(ip)%coor(i)=particle(ip)%coor(i)
                enddo
                return
              end Subroutine
!====================
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine update_wall_position
    use system
    use solid
    implicit none
!      
      do i=1,ne
      if(wall(i)%flag.eq.1)then
        if(wall(i)%rate.eq.1)then
          wall(i)%coor(1)=wall(i)%coor(1)+wall(i)%vel(1)*dtime
          wall(i)%coor(3)=wall(i)%coor(1)
          wall(i)%coor(2)=wall(i)%coor(2)+wall(i)%vel(2)*dtime
          wall(i)%coor(4)=wall(i)%coor(2)
        else
          
        end if

      end if
      end do
!      write(*,*) wall(3)%vel(2),wall(3)%coor(1:2)
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine  write_results0
    use system
    use solid
    use fluid
    use fs_inter
    implicit none
      integer no_output,n0
      real(8) x1,y1,z1,u2,rad
	  save    no_output
!

!.....Write header for postprocessing with TECPLOT software
  if(first)then
        no_output=1
        write(11,*) 'TITLE = DE3D' 
        write(11,*) 'VARIABLES = X, Y, Z, VV, Diameter' 
!        write(11,*) 'ZONE T=',istep
	      write(11,*) 'ZONE T="Step:',istep,'", I=',np, ', J=1, K=1'
!       write(11,*) 'I=',np,'J=1, K=1'
        write(11,*) "ZONETYPE=ORDERED, DATAPACKING=POINT"
	else
!	  write(11,*) 'ZONE T=',istep
	      write(11,*) 'ZONE T="Step:',istep,'", I=',np, ', J=1, K=1'
!        write(11,*) 'I=',np,'J=1, K=1'
        write(11,*) "ZONETYPE=ORDERED, DATAPACKING=POINT"
	endif
!
	do i=1,np
		  u2=sqrt(particle(i)%U(1)**2+particle(i)%U(2)**2+particle(i)%U(3)**2)
          write(11,110) particle(i)%coor(1:3),u2,2.0*particle(i)%radius
	enddo
!

	first=.false.
!
  110 format(8(1x,g10.4)) 
  120 format(4(1x,i4)) 
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine  write_spheres
  use system
  use solid
  use fluid
  use fs_inter
  implicit none
  integer no_output,n0
  real(8) x1,y1,z1,u2,rad
  save no_output
  !
  character(len=100) :: output_name,output_name2
  character(len=20) :: step_string
  character(len=100) :: output_path,output_path2
  ! for paraview files
  write(step_string,'(I8.8)') istep
  output_name = 'particles_' // trim(step_string) // '.vtp'
  output_name2 = 'particles_' // trim(step_string) // '.csv'
  output_path = 'paraview/' // trim(output_name)
  output_path2 = 'paraview/' // trim(output_name2)
  ! output each istep/file
  write(99,'(A,I0,A,A,A)') ' <DataSet timestep="', istep, &
  '" group="" part="0" file="', trim(output_name), '"/>'

    open(unit=999,file=output_path,status='replace')
    !output csv file
!    open(unit=9999,file=output_path2,status='replace')
!    write(9999,'(A)') ' ID, X, Y, Z, Diameter, VX, VY, VZ, FX, FY, FZ'


  ! ---- XML VTK PolyData file for particles ----
  write(999,'(A)') '<?xml version="1.0"?>'
  write(999,'(A)') '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
  write(999,'(A)') ' <PolyData>'
  write(999,'(A,I0,A,I0,A)') ' <Piece NumberOfPoints="', np, '" NumberOfVerts="', np, &
  '" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'

  ! ---- Point data: radius, velocity, hydrodynamic force ----
  write(999,'(A)') ' <PointData Scalars="radius" Vectors="velocity">'

  ! ---- Scalar: radius ----
  write(999,'(A)') ' <DataArray type="Float64" Name="radius" NumberOfComponents="1" format="ascii">'
  do i=1,np
    write(999,'(ES24.16)') particle(i)%radius
!    write(9999,*) i, ",",particle(i)%coor(1),",",particle(i)%coor(2),",",particle(i)%coor(3),",",2.0*particle(i)%radius,",",particle(i)%U(1),",",particle(i)%U(2), &
!    ",",particle(i)%U(3),",",particle(i)%F(1),",",particle(i)%F(2),",",particle(i)%F(3)
  end do
  write(999,'(A)') ' </DataArray>'

  ! ---- Vector: Velocity ----
  write(999,'(A)') ' <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">'
    do i=1,np
        write(999,'(ES24.16)') particle(i)%U(1:3)
!        write(9999,*) i, ",",particle(i)%coor(1),",",particle(i)%coor(2),",",particle(i)%coor(3),",",2.0*particle(i)%radius,",",particle(i)%U(1),",",particle(i)%U(2), &
!        ",",particle(i)%U(3),",",particle(i)%F(1),",",particle(i)%F(2),",",particle(i)%F(3)
    end do
  write(999,'(A)') ' </DataArray>'

  ! ---- Vector: hydrodynamic force ----
  write(999,'(A)') ' <DataArray type="Float64" Name="hydrodynamic_force" NumberOfComponents="3" format="ascii">'
    do i=1,np
    write(999,'(3ES24.16)') particle(i)%F(1:3)
    end do
  write(999,'(A)') ' </DataArray>'

  write(999,'(A)') ' </PointData>'

  ! ---- Points: particle centers ----
  write(999,'(A)') ' <Points>'
  write(999,'(A)') ' <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i=1,np
    write(999,'(3ES24.16)') particle(i)%coor(1:3)
    end do
  write(999,'(A)') ' </DataArray>'
  write(999,'(A)') ' </Points>'

  ! ---- Vertices: one vertex cell per particle ----
  write(999,'(A)') ' <Verts>'

  write(999,'(A)') ' <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i=0,np-1
    write(999,'(I0)') i
    end do
  write(999,'(A)') ' </DataArray>'

  write(999,'(A)') ' <DataArray type="Int32" Name="offsets" format="ascii">'
    do i=1,np
    write(999,'(I0)') i
    end do
  write(999,'(A)') ' </DataArray>'

  write(999,'(A)') ' </Verts>'

  write(999,'(A)') ' </Piece>'
  write(999,'(A)') ' </PolyData>'
  write(999,'(A)') '</VTKFile>'

    close(999)
!    close(9999)
  return
End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine  write_results1
    use system
    use solid
    use fluid
    use fs_inter
    implicit none
      integer no_output,ixs,ixe,iob,npoints
      real(8) x1,y1,z1,ux,uy,uz,press,c_squ,den,velocity
    character(len=100) :: output_name
    character(len=20)  :: step_string
    character(len=100) :: output_path
	  save    no_output
!
    npoints = (nx+1) * (ny+1) * (nz+1)

!.....square fluid's speed of sound
    c_squ=3.d0/8.d0

! for paraview files
    write(step_string,'(I8.8)') istep
    output_name = 'Fluids_' // trim(step_string) // '.vts'
    output_path = 'paraview/' // trim(output_name)

! output each istep/file
    write(98,'(A,I0,A,A,A)') '    <DataSet timestep="', istep, &
         '" group="" part="0" file="', trim(output_name), '"/>'

    open(unit=999,file=output_path,status='replace')

! ---- XML VTK StructuredGrid file for fluid ----
    write(999,'(A)') '<?xml version="1.0"?>'
    write(999,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(999,'(A,I0,A,I0,A,I0,A)') '  <StructuredGrid WholeExtent="0 ', nx, ' 0 ', ny, ' 0 ', nz, '">'
    write(999,'(A,I0,A,I0,A,I0,A)') '    <Piece Extent="0 ', nx, ' 0 ', ny, ' 0 ', nz, '">'

! -------------------------
! Write nodal velocity data
! -------------------------
    write(999,'(A)') '      <PointData Vectors="velocity">'
    write(999,'(A)') '        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">'

    do k = 0, nz
        do j = 0, ny
            do i = 0, nx
                call nodal_velocity(i,j,k,nx,ny,nz,ux,uy,uz,den)

                if(BODYF.eq.1)then
                  uy=uy-0.5d0*gacce
                end if

                velocity=sqrt(ux*ux+uy*uy+uz*uz)
                if(umax.lt.velocity) umax=velocity
                write(999,'(3ES24.16)') ux, uy, uz
            end do
        end do
    end do
    write(999,'(A)') '        </DataArray>'

!   ! ---- Scalar: node flag ----
    write(999,'(A)') '        <DataArray type="Int8" Name="obstacle" NumberOfComponents="1" format="ascii">'

    do k = 0, nz
        do j = 0, ny
            do i = 0, nx
                write(999,'(I8)') node(i,j,k)%obst
            end do
        end do
    end do
    write(999,'(A)') '        </DataArray>'
!    
!   ! ---- Scalar: density ----
    write(999,'(A)') '        <DataArray type="Float64" Name="density" NumberOfComponents="1" format="ascii">'

    do k = 0, nz
        do j = 0, ny
            do i = 0, nx
              call nodal_velocity(i,j,k,nx,ny,nz,ux,uy,uz,den)
              write(999,'(ES24.16)') den
            end do
        end do
    end do
    write(999,'(A)') '        </DataArray>'

    write(999,'(A)') '      </PointData>'

    write(999,'(A)') '      <CellData>'
    write(999,'(A)') '      </CellData>'

! ---- Write coordinates ----
    write(999,'(A)') '      <Points>'
    write(999,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'

    do k = 0, nz
        do j = 0, ny
            do i = 0, nx
                write(999,'(3ES24.16)') dble(i), dble(j), dble(k)
            end do
        end do
    end do

    write(999,'(A)') '        </DataArray>'
    write(999,'(A)') '      </Points>'

    write(999,'(A)') '    </Piece>'
    write(999,'(A)') '  </StructuredGrid>'
    write(999,'(A)') '</VTKFile>'

    close(999)

    return

! !

! !.....Write header for postprocessing with TECPLOT software
! !.....uncomment following line, if this header should be printed
!       if(istep.eq.1)then
! 	  no_output=1
!         write(11,*) 'TITLE = LB3D' 
!         write(11,*) 'VARIABLES = X, Y, Z, VX, VY, VZ, VV,  PRESS' 
! 	  write(11,*)'ZONE T="Step:',istep,'",I=',nx,', J=',ny, ', &
!                                            K=',nz, ', F=POINT'
! 	else
! 	  write(11,*)'ZONE T="Step:',istep,'",I=',nx,', J=',ny, ',&
!                  K=',nz, ', F=POINT', ' VARSHARELIST=([1,2,3]=1)'
! 	endif
! ! Loop over all nodes
!       do iz=0,nz
! 	     do iy=0,ny
!           do ix =0,nx
!           if(node(ix,iy,iz)%obst.ne.0) then 
! !
! !...........obstacle indicator
!             iob=1
! !
! !...........velocity components = 0
!             ux=0.d0
!             uy=0.d0
! 	          uz=0.d0
! !...........pressure = average pressure
!             press=d0*c_squ
!           else
!             call nodal_velocity(ix,iy,iz,nx,ny,nz,ux,uy,uz,den)
! !...........pressure
!             press=den*c_squ
! ! 
!             iob=0
!           end if
!           velocity=sqrt(ux*ux+uy*uy+uz*uz)
!           if(umax.lt.velocity) umax=velocity
! !
! !      press=taostar(ix,iy)-tao
! !.......write results to file
!           if(first)then
! 	          write(11,100) ix, iy, iz, ux, uy, uz, velocity, press
! 	        else
!             write(11,110) ux, uy, uz, velocity,press
! 	        endif
!         enddo
!       enddo
!     end do
! !
! !.....wall boundary
! !	if(first)then
! !	    write(11,120) vertex(1,1),vertex(2,1),1,ne+1
! !	  do i=1,ne
! !          write(11,130)(vertex(1,i)-vertex(1,1)),&
! !     &                 (vertex(2,i)-vertex(2,1))
! !	  enddo
! !	  write(11,130) vertex(1,1),vertex(2,1)
! !	endif
! !
! !.....write stationary particles
! 	  if(first.and.ns.gt.0)then
! 	    do i=1,ns
!           write(11,140) particle(i)%coor(1),particle(i)%coor(2),particle(i)%coor(3)
! 	      write(11,*) particle(i)%radius
! 	    enddo
! 	  endif	
! !.....write moving particles
! !
! 	  first=.false.
! 	  no_output=no_output+1
!   100 format(1x,i5,i5,i5,1x,g10.4,1x,g10.4,1x,g10.4,1x,g10.4,1x,g10.4) 
!   110 format(            1x,g10.4,1x,g10.4,1x,g10.4,1x,g10.4,1x,g10.4) 
!   140 format(            1x,g10.4,1x,g10.4,1x,g10.4) 
      ! return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine update_moving_particle_nodes
    use system
    use solid
    use fluid
    use fs_inter
    implicit none       
!      integer nnb,nni
!      write(*,*) noi, nob, nof
!.....Release old nodes in obst
      do i=1,noi
	    ix=lni(1,i)
	    iy=lni(2,i)
      iz=lni(3,i)
	    node(ix,iy,iz)%obst=0
      end do
!      write(*,*) "noi"
      do i=1,nob
	    ix=lnb(1,i)
	    iy=lnb(2,i)
      iz=lnb(3,i)
	    node(ix,iy,iz)%obst=0
      end do
!            write(*,*) "nob",nof,3*MB
      do i=1,nof
	      ix=lnf(1,i)
	      iy=lnf(2,i)
        iz=lnf(3,i)
	      node(ix,iy,iz)%obst=0
      end do
!      write(*,*) "nof"
!.....Get new nodes
      noi=0
	    nni=0
	    nob=0
	    nnb=0
	    nof=0
	    nnf=0
!
!      write(*,*) ip, nob
      p1=>head_lpm
      do while(.true.)
        ip=p1%num
        if(particle(ip)%active==0)goto 100
        call nodes_of_moving_particles(ip)
!        write(*,*) ip, nob
!.......number of boundary nodes for current particle
	    particle(ip)%nb=nob-nnb
!        write(*,*) ip,particle(ip)%nb
!.......pointer pointing to the boundary node list
	    particle(ip)%pb=nnb+1
!.......number of interior nodes for current particle
	    particle(ip)%ni=noi-nni
!.......pointer pointing to the interior node list
	    particle(ip)%pn=nni+1
!.......number of fluid boundary nodes for current particle
	    particle(ip)%nf=nof-nnf
!.......pointer pointing to the interior node list
	    particle(ip)%pf=nnf+1
!        write(*,*) ip,particle(ip)%ni,particle(ip)%pn
	    nnb=nob
	    nni=noi
      nnf=nof
    !  open(10002,file='Particle_nodes_list.txt')      
   !   do i=1,nob
   !     write(10002,*) i, lnb(1:3,i)
  !    enddo 
  !    close(10002)
 !if(istep==2)      stop
!
100        if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do
!  write(*,*) "MOVING--nnb,nni,nnf: ",nnb,nni,nnf
    return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine redistribute(accel)
    use fluid
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!.....Weighting factors
      w2=d0*accel/72.d0
	    w1=8.d0*w2
!
      if(bc_mode.eq.1)then
! check this setting later?? 2 or 1?
        do iz=0,nz
        do iy=0, ny
!
!.......Accelerate flow only on non-occupied nodes
        if (node(0,iy,iz)%obst.eq.0 .and.&
!.........check to avoid negative densities
          temp(2,1,iy,iz)-w1 .gt.0..and. &
          temp(8,1,iy,iz)-w2 .gt.0..and. &
          temp(10,1,iy,iz)-w2 .gt.0..and. &
          temp(12,1,iy,iz)-w2 .gt.0.)then 
          temp(1,1,iy,iz)=temp(1,1,iy,iz)+w1
          temp(2,1,iy,iz)=temp(2,1,iy,iz)-w1
          temp(7,1,iy,iz)=temp(7,1,iy,iz)+w2
          temp(8,1,iy,iz)=temp(8,1,iy,iz)-w2
          temp(9,1,iy,iz)=temp(9,1,iy,iz)+w2
          temp(10,1,iy,iz)=temp(10,1,iy,iz)-w2
          temp(11,1,iy,iz)=temp(11,1,iy,iz)+w2
          temp(12,1,iy,iz)=temp(12,1,iy,iz)-w2
          temp(13,1,iy,iz)=temp(13,1,iy,iz)+w2
          temp(14,1,iy,iz)=temp(14,1,iy,iz)-w2
        end if
        enddo
        end do
      else
! in Y direction need to revise in the future
!
      end if
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine nodal_velocity(ix,iy,iz,nx,ny,nz,ux,uy,uz,den)        
    use fluid,only:temp
    implicit none
!      
      integer i,ix,iy,iz,nx,ny,nz
      real(8) den,ux,uy,uz
!
!.....Integrate local density
      den=0.d0 
      do i=0,14 
         den=den+temp(i,ix,iy,iz)
	  enddo
!
!.....x-, and y- and z- velocity components
      ux=(temp(1,ix,iy,iz)+temp(7,ix,iy,iz)+temp(9,ix,iy,iz) &
       +temp(11,ix,iy,iz)+temp(13,ix,iy,iz) &
       -(temp(2,ix,iy,iz)+temp(8,ix,iy,iz)+temp(10,ix,iy,iz) &
       +temp(12,ix,iy,iz)+temp(14,ix,iy,iz)))/den 
	    uy=(temp(3,ix,iy,iz)+temp(7,ix,iy,iz)+temp(9,ix,iy,iz) &
       +temp(12,ix,iy,iz)+temp(14,ix,iy,iz) &
       -(temp(4,ix,iy,iz)+temp(8,ix,iy,iz)+temp(10,ix,iy,iz) &
       +temp(11,ix,iy,iz)+temp(13,ix,iy,iz)))/den
      uz=(temp(5,ix,iy,iz)+temp(7,ix,iy,iz)+temp(10,ix,iy,iz) &
       +temp(11,ix,iy,iz)+temp(14,ix,iy,iz) &
       -(temp(6,ix,iy,iz)+temp(8,ix,iy,iz)+temp(9,ix,iy,iz) &
       +temp(12,ix,iy,iz)+temp(13,ix,iy,iz)))/den
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine equilibrium_function(ux,uy,uz,den,neq)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!.....neq=feq:equilibrium distribution functions
      real(8) neq(0:14)
!
!.....square velocity
      u_squ=ux*ux+uy*uy+uz*uz
!.....equilibrium densities
      w2=den/72.d0
	    w1=8.d0*w2
	    w0=2.d0*w1
	    tmp=1.d0-(3.0d0/2.0d0)*u_squ
!
!.....zero velocity density
      neq(0)=w0*tmp
!
!.....axis speeds (factor: w1)
      neq(1)=w1*(tmp+(3.0d0+9.0d0/2.0d0*ux)*ux)
      neq(2)=neq(1)-6.0d0*w1*ux
      neq(3)=w1*(tmp+(3.0d0+9.0d0/2.0d0*uy)*uy)
      neq(4)=neq(3)-6.0d0*w1*uy
      neq(5)=w1*(tmp+(3.0d0+9.0d0/2.0d0*uz)*uz)
      neq(6)=neq(5)-6.0d0*w1*uz
!
!.....diagonal speeds (factor: w2)
!.....ei*u=(eix*ux)+(eiy*uy)
      uxyz=ux+uy+uz
      neq(7)=w2*(tmp+(3.0d0+9.0d0/2.0d0*uxyz)*uxyz)
      neq(8)=neq(7)-6.0d0*w2*uxyz
      uxyz=ux+uy-uz
      neq(9)=w2*(tmp+(3.0d0+9.0d0/2.0d0*uxyz)*uxyz)
      neq(10)=neq(9)-6.0d0*w2*uxyz
      uxyz=ux-uy+uz
      neq(11)=w2*(tmp+(3.0d0+9.0d0/2.0d0*uxyz)*uxyz)
      neq(12)=neq(11)-6.0d0*w2*uxyz
      uxyz=ux-uy-uz
      neq(13)=w2*(tmp+(3.0d0+9.0d0/2.0d0*uxyz)*uxyz)
      neq(14)=neq(13)-6.0d0*w2*uxyz
!
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine nodal_volume_covered_by_particles(ix,iy,iz,ipar,vol)
    use fluid,only:nx
    use solid,only:boundary_dem
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
	  real(8)   ipar(4),rand
    integer seed

!
    xc=ipar(1)
	  yc=ipar(2)
	  zc=ipar(3)
	  rad=ipar(4)
	  rad2=rad*rad
	  xmin=ix-0.5
	  ymin=iy-0.5
	  zmin=iz-0.5
	  N=100
    M=0
  seed=1234567
  call RANDOM_SEED(seed)
	do i=1,N
    call random_number(rand)
	  x=xmin+rand
!        write(*,*) rand,x
    call random_number(rand)
	  y=ymin+rand
!       write(*,*) rand,y
    call random_number(rand)
	  z=zmin+rand
!        write(*,*) rand,z
!        stop
	  dx=x-xc
	  dy=y-yc
	  dz=z-zc  
	  if(dx*dx+dy*dy+dz*dz.le.rad2) M=M+1
	enddo
	vol=M*1.d0/N
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine relaxation
    use system
    use solid
    use fluid
    use fs_inter 
    use MRT, only : MRT_collision, MRT_collision_body_forcing
    implicit none
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!    
!.....local variables: feq:equilibrium distribution function;  peq:equilibrium distribution function for particles 
!.....ux:velocity of current node in x coordinate, uy:velocity of current node in y coordinate
      real(8) feq(0:14),peq(0:14),ux,uy,uz,den
      real(8)::ex(14)=(/1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1/)
      real(8)::ey(14)=(/0,0,1,-1,0,0,1,-1,1,-1,-1,1,-1,1/)
      real(8)::ez(14)=(/0,0,0,0,1,-1,1,-1,-1,1,1,-1,-1,1/)
      real(8) Fx, Fy, Fz
      integer ixs,ixe,nn,ipr
!
! write(*,*) "relaxation is working at the beginning"
!.....Initialise force vector
      if(IMB .OR. IBM)then
      if(nm>0)then
      p1=>head_lpm
      do while(.true.)
        ip=p1%num
	    do j=1,6
          particle(ip)%F(j)=0.0d0
	    enddo
!
        if(.not. associated(p1%next)) exit
        p1=>p1%next
      end do
      endif
      end if
!
!.....Loop over fixed wall boundary nodes
      do i=1,now0
!      write(*,*) i,now0
	    ix=lnw0(1,i)
	    iy=lnw0(2,i)
      iz=lnw0(3,i)
	      do j=0,14
            node(ix,iy,iz)%fdd(j)=temp(j,ix,iy,iz)
	      enddo
      end do
!      write(*,*) "now=",now0
!      stop
!.....Loop over stationary particle nodes
      if (ns.gt.0)then
        p3=>head_lns
        do while(.true.)
          ix=p3%coord(1)
          iy=p3%coord(2)
          iz=p3%coord(3)   
          do j=0,14
             node(ix,iy,iz)%fdd(j)=temp(j,ix,iy,iz)
	      enddo
          if(.not. associated(p3%next)) exit
          p3=>p3%next
        end do
      end if
!.....Loop over free nodes
      do iz=0,nz
      do iy=0,ny
!	    ixs=yb(1,iy)
!	    ixe=yb(2,iy)
	      do ix=0,nx
	      if(node(ix,iy,iz)%obst.eq.0)then
	        call nodal_velocity(ix,iy,iz,nx,ny,nz,ux,uy,uz,den)

          if(use_MRT)then

            if(BODYF.eq.1)then
              !stop "ERROR: BODYF=1 is not implemented for MRT Relaxation yet"
              ! body force density and gravity in y-dir
              Fx = 0.0d0
              Fy = -den*gacce
              Fz = 0.0d0

              ! Note that guo forcing is already within this subroutine
              call MRT_collision_body_forcing(temp(:,ix,iy,iz), &
                                              node(ix,iy,iz)%fdd, &
                                              ux,uy,uz,den, &
                                              Fx,Fy,Fz)

            else
              ! Normal collision
              call MRT_collision(temp(:, ix, iy, iz), node(ix, iy, iz)%fdd, ux, uy, uz, den)
              
            end if

          else

	          call equilibrium_function(ux,uy,uz,den,feq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(BODYF.eq.1)then
              uy=uy-0.5*gacce
            end if

	        call equilibrium_function(ux,uy,uz,den,feq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	        if(BODYF.eq.1) call body_force_density(ix,iy,iz,ux,uy,uz,den)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) rtao,feq
!stop
	          do j=0,14
              if(BODYF.eq.1)then
                node(ix,iy,iz)%fdd(j)=temp(j,ix,iy,iz)+rtao*(feq(j)-temp(j,ix,iy,iz))+fi_body(j)
              else
                node(ix,iy,iz)%fdd(j)=temp(j,ix,iy,iz)+rtao*(feq(j)-temp(j,ix,iy,iz))
              end if
	          enddo
          end if

	      endif
        enddo

      enddo
      end do
!

!.....Loop over moving particles
    if(nm>0)then
      if(IMB)then
        p1=>head_lpm
        loop001:do while(.true.)
          ip=p1%num
!
	      if(particle(ip)%active.eq.0) goto 100
!.........Loop over particle boundary nodes
          nn=particle(ip)%nb
	      ipr=particle(ip)%pb-1
	      do j=1,nn
	        ix=lnb(1,ipr+j)
	        iy=lnb(2,ipr+j)
          iz=lnb(3,ipr+j)
	        call moving_particle_relaxation1(0,ip,ix,iy,iz,nx,ny,nz,rtao,istep,ex,ey,ez)
	      enddo
!    if(istep==2)     write(*,*) "nn", nn,particle(ip)%F(1:3)

!.........Loop over particle fluid boundary nodes
          nn=particle(ip)%nf
	      ipr=particle(ip)%pf-1
	      do j=1,nn
	        ix=lnf(1,ipr+j)
	        iy=lnf(2,ipr+j)
          iz=lnf(3,ipr+j)
	        call moving_particle_relaxation1(0,ip,ix,iy,iz,nx,ny,nz,rtao,istep,ex,ey,ez)
	      enddo
!    if(istep==2)     write(*,*) "nn", nn,  particle(ip)%F(1:3)
!.........Loop over particle interior nodes
          nn=particle(ip)%ni
	      ipr=particle(ip)%pn-1
	      do j=1,nn
	        ix=lni(1,ipr+j)
	        iy=lni(2,ipr+j)
          iz=lni(3,ipr+j)
	        call moving_particle_relaxation1(1,ip,ix,iy,iz,nx,ny,nz,rtao,istep,ex,ey,ez)
	      enddo
 !               write(*,*) particle(ip)%F(1:3)
! if(istep==2)                stop
!
   100    if(.not. associated(p1%next)) exit loop001
          p1=>p1%next
        end do loop001
      else if(IBM)then



      end if
    end if
!
    return
  End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine moving_particle_relaxation1(mode,ip,ix,iy,iz,nx,ny,nz,rtao,istep,ex,ey,ez)
    use solid,only:particle,gacce,part_vir
    use fluid,only:node,temp,BODYF,fi_body,tao
    use system, only: use_MRT
    use MRT, only : MRT_collision, MRT_IMB_collision, MRT_collision_body_forcing, MRT_IMB_collision_body_forcing
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real(8) ex(14),ey(14),ez(14),u1,u2,u3
!.....local variables
      integer  opposite_direction
      real(8)  feq(0:14),peq(0:14),par(4),vol,rad2
      ! real(8)  f_mrt_post(0:14), delta_f
      real(8)  os_mrt(0:14)
      real(8) Fx, Fy, Fz
!
!.....Compute nodal density and velocities
      call nodal_velocity(ix,iy,iz,nx,ny,nz,ux,uy,uz,den)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            write(*,*) uy
      ! Get body force density
      if(BODYF.eq.1)then

        Fx=0.0d0
        Fy=-den*gacce
        Fz=0.0d0

      else

        Fx=0.0d0
        Fy=0.0d0
        Fz=0.0d0

      end if


      if(BODYF.eq.1 .and. .not.use_MRT)then
        uy=uy-0.5*gacce
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....Compute equilibrium function
      if(.not. use_MRT)then
        call equilibrium_function(ux,uy,uz,den,feq)
      end if
!      write(*,*) den,feq
!      stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....FOR MRT
      ! Note that currently feq and peq are from BGK equilibrium
      ! Perhaps a MRT version of equilibrium in vel space is needed to be fully consistnent for IMB
      ! This is done now


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  if(BODYF.eq.1 .and. .not.use_MRT) call body_force_density(ix,iy,iz,ux,uy,uz,den)
!      write(*,*) BODYF,uy,fi_body(7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....Compute nodal area covered by particle
!.....for interior node
      if(mode.eq.1)then
	    B=1.d0
	  elseif(mode.eq.2)then
      B=0.0d0
    else
        par(1:3)=particle(ip)%coor(1:3)
        par(4)=particle(ip)%radius
        rad2=(par(4)+2)*(par(4)+2)
        ! check the nodes for virtual particles in PBC
        if(inside_particle(ix,iy,iz,par(1),par(2),par(3),rad2)==0)then
          par(1:3)=part_vir(ip)%coor(1:3)
        endif
	    call nodal_volume_covered_by_particles(ix,iy,iz,par,vol)
! B is a function of solid volume fraction
!	     B=vol
      B=vol*(tao-0.5)/((1.0-vol)+(tao-0.5))
! used the constant for debugging
!      B=0.2
!      write(*,*) "volume is ", B
!      stop
	  endif
!
	  if(B.gt.0.0d0)then
!.......Particle equilibrium function
        u1=particle(ip)%U(1)
        u2=particle(ip)%U(2)
        u3=particle(ip)%U(3)
!        write(*,*) u1,u2,u3,den,peq
!        stop
! 
        fx=0.d0
        fy=0.d0
        fz=0.0d0
        fm=0.d0
        rx=ix-par(1)
        ry=iy-par(2)
        rz=iz-par(3)
          
          if(use_MRT)then
              ! Fully MRT-space IMB:
              ! m_eq_f = m_eq(rho, fluid velocity)
              ! m_eq_p = m_eq(rho, particle velocity)
              ! m_new  = m_post + B * (m_eq_p - m_eq_f)

                  if(BODYF.eq.1)then
                    !.....MRT-IMB with body forcing
                    call MRT_IMB_collision_body_forcing(              &
                        temp(:,ix,iy,iz),                             &
                        node(ix,iy,iz)%fdd,                           &
                        os_mrt,                                       &
                        ux,uy,uz,den,                                 &
                        u1,u2,u3,B,                                   &
                        Fx,Fy,Fz)

                  
                  else
                    !.....MRT-IMB without body forcing
                    call MRT_IMB_collision(                           &
                        temp(:,ix,iy,iz),                             &
                        node(ix,iy,iz)%fdd,                           &
                        os_mrt,                                       &
                        ux,uy,uz,den,                                 &
                        u1,u2,u3,B)

                  end if

              do i = 1, 14
                fx = fx + os_mrt(i) * ex(i)
                fy = fy + os_mrt(i) * ey(i)
                fz = fz + os_mrt(i) * ez(i)
              end do
          else
            ! BGK IMB
            call equilibrium_function(u1,u2,u3,den,peq)
            do i=0,14

              os=peq(i)-temp(i,ix,iy,iz)+(1.d0-rtao)*(temp(i,ix,iy,iz)-feq(i))
               if(BODYF.eq.1)then
                  node(ix,iy,iz)%fdd(i)=temp(i,ix,iy,iz)+rtao*(1.d0-B)*(feq(i)-temp(i,ix,iy,iz))+B*os+(1.0-B)*fi_body(i)
                else
                  node(ix,iy,iz)%fdd(i)=temp(i,ix,iy,iz)+rtao*(1.d0-B)*(feq(i)-temp(i,ix,iy,iz))+B*os
               end if
    !.........Compute forces on particle
              if(i.gt.0)then
              fx=fx+os*ex(i)
              fy=fy+os*ey(i)
              fz=fz+os*ez(i)
    ! in the future do this torque	        fm=fm+fy*rx-fx*ry
              endif
            enddo
      end if
	      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	    if(mode.eq.1)then
      particle(ip)%F(1)=particle(ip)%F(1)-B*fx
	    particle(ip)%F(2)=particle(ip)%F(2)-B*fy
	    particle(ip)%F(3)=particle(ip)%F(3)-B*fz
!      write(*,*) fx,fy,fz
!      stop
      ! in the future consider torque
!        else
!        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     zero volume
	  else
        if ( use_MRT ) then
          if(BODYF.eq.1)then

            call MRT_collision_body_forcing(              &
                 temp(:,ix,iy,iz),                        &
                 node(ix,iy,iz)%fdd,                      &
                 ux,uy,uz,den,                            &
                 Fx,Fy,Fz)

          else

            call MRT_collision(                           &
                 temp(:,ix,iy,iz),                        &
                 node(ix,iy,iz)%fdd,                      &
                 ux,uy,uz,den)

          end if
          ! node(ix,iy,iz)%fdd(i) = f_mrt_post(i)
        else
          do i=0,14
            if(BODYF.eq.1)then
              node(ix,iy,iz)%fdd(i)=temp(i,ix,iy,iz)+rtao*(feq(i)-temp(i,ix,iy,iz))+fi_body(i)
            else
                node(ix,iy,iz)%fdd(i)=temp(i,ix,iy,iz)+rtao*(feq(i)-temp(i,ix,iy,iz))
            end if
          enddo
        end if
	    endif
	  return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine bounceback
    use system
    use solid
    use fluid
    use fs_inter
    implicit none
!
     real(8) w2,w1,w0
    
!.....Weighting factors
      w2=1.d0/72.d0
	    w1=8.d0*w2
!.....Loop over fixed wall boundary nodes
      do i=1, now0
        ix=lnw0(1,i)
	      iy=lnw0(2,i)
        iz=lnw0(3,i)
        call bounceback_stationary_particles(ix,iy,iz,nx,ny,nz)
      end do
!
!.....Loop over stationary particle nodes
      if(ns .gt. 0)then
        p3=>head_lns
        do while(.true.)
          ix=p3%coord(1)
          iy=p3%coord(2)    
          iz=p3%coord(3)
          call bounceback_stationary_particles(ix,iy,iz,nx,ny,nz)
 !
          if(.not. associated(p3%next)) exit
          p3=>p3%next
        end do
      end if
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine bounceback_stationary_particles(ix,iy,iz,nx,ny,nz)
    use fluid,only:node,temp
    implicit none
      integer ix,iy,iz,nx,ny,nz
!
!.....Rotate all ensities and write back to node
      node(ix,iy,iz)%fdd(1) = temp(2,ix,iy,iz)
	    node(ix,iy,iz)%fdd(2) = temp(1,ix,iy,iz)
	    node(ix,iy,iz)%fdd(3) = temp(4,ix,iy,iz)
	    node(ix,iy,iz)%fdd(4) = temp(3,ix,iy,iz)
	    node(ix,iy,iz)%fdd(5) = temp(6,ix,iy,iz)
	    node(ix,iy,iz)%fdd(6) = temp(5,ix,iy,iz)
	    node(ix,iy,iz)%fdd(7) = temp(8,ix,iy,iz)
	    node(ix,iy,iz)%fdd(8) = temp(7,ix,iy,iz)
	    node(ix,iy,iz)%fdd(9) = temp(10,ix,iy,iz)
	    node(ix,iy,iz)%fdd(10) = temp(9,ix,iy,iz)
	    node(ix,iy,iz)%fdd(11) = temp(12,ix,iy,iz)
	    node(ix,iy,iz)%fdd(12) = temp(11,ix,iy,iz)
	    node(ix,iy,iz)%fdd(13) = temp(14,ix,iy,iz)
	    node(ix,iy,iz)%fdd(14) = temp(13,ix,iy,iz)
!
!
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine propagate
!========
    use system,only:istep,periodic_flag,implement
!========
    use fluid,only:node,temp,nx,ny,nz
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!.....local variables
      integer xr,xl,yb,yf,zu,zd
!
!.....Loop over all nodes
      do iz = 0, nz
        do iy = 0, ny
          do ix = 0, nx
!
!...........neighbour nodes 
            xr = ix + 1
            xl = ix - 1
            yb = iy + 1
            yf = iy - 1
	          zu = iz + 1
            zd = iz - 1
!
!=======
          if(periodic_flag .and. implement==2)then
            if (xr>nx) xr=0
            if (xl<0)  xl=nx
            if (yb>ny) yb=0
            if (yf<0)  yf=ny
            if (zu>nz) zu=0
            if (zd<0)  zd=nz
!            write(*,*) "PBC 2 is used"
          endif
!========
!...........zero
            temp(0,ix,iy,iz) = node(ix,iy,iz)%fdd(0)
!
!...........axis
            temp(1,xr,iy,iz) = node(ix,iy,iz)%fdd(1)
	          temp(2,xl,iy,iz) = node(ix,iy,iz)%fdd(2)
	          temp(3,ix,yb,iz) = node(ix,iy,iz)%fdd(3)
	          temp(4,ix,yf,iz) = node(ix,iy,iz)%fdd(4)
	          temp(5,ix,iy,zu) = node(ix,iy,iz)%fdd(5)
	          temp(6,ix,iy,zd) = node(ix,iy,iz)%fdd(6)
!
!...........diagonal
            temp(7,xr,yb,zu) = node(ix,iy,iz)%fdd(7)
	          temp(8,xl,yf,zd) = node(ix,iy,iz)%fdd(8)
	          temp(9,xr,yb,zd) = node(ix,iy,iz)%fdd(9)
	          temp(10,xl,yf,zu)= node(ix,iy,iz)%fdd(10)
	          temp(11,xr,yf,zu)= node(ix,iy,iz)%fdd(11)
	          temp(12,xl,yb,zd)= node(ix,iy,iz)%fdd(12)
	          temp(13,xr,yf,zd)= node(ix,iy,iz)%fdd(13)
	          temp(14,xl,yb,zu)= node(ix,iy,iz)%fdd(14)
          enddo
        enddo
      end do
!
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine boundary_treatment_regular
    use system,only:wall,vertex,periodic_flag,istep,implement
    use fluid,only:bctype,bc,neb,eb,nx,ny,nz,temp,bc_mode,node
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!.....Periodic boundaries at inlet and outlet
!       if(bctype.eq.1)then
!         if(bc_mode.eq.1)then
! 	        do iz=0,nz
! 	        do iy=0,ny
!             temp(1,0,iy,iz)=temp(1,nx+1,iy,iz)
! 	          temp(7,0,iy,iz)=temp(7,nx+1,iy,iz)
! 	          temp(9,0,iy,iz)=temp(9,nx+1,iy,iz)
! 	          temp(11,0,iy,iz)=temp(11,nx+1,iy,iz)
! 	          temp(13,0,iy,iz)=temp(13,nx+1,iy,iz)
! 	          temp(2,nx,iy,iz)=temp(2,-1,iy,iz)
! 	          temp(8,nx,iy,iz)=temp(8,-1,iy,iz)
! 	          temp(10,nx,iy,iz)=temp(10,-1,iy,iz)
! 	          temp(12,nx,iy,iz)=temp(12,-1,iy,iz)
! 	          temp(14,nx,iy,iz)=temp(14,-1,iy,iz)
!           enddo
!           enddo
!         else
! !          do ix=1,nx-1
! !            temp(4, ix,ny )=temp(4,ix,-1)
! !            temp(7, ix,ny )=temp(7,ix,-1)
! !            temp(8, ix,ny )=temp(8,ix,-1)
! !            temp(2,ix,0 )=temp(2,ix,ny+1)
! !            temp(5,ix,0 )=temp(5,ix,ny+1)
! !            temp(6,ix,0 )=temp(6,ix,ny+1) 
! !          end do
!         end if

      ! Periodic boundaries in y and z directions
      if(bctype.eq.3)then

        if(periodic_flag .and. implement==1)then
	        do ix=0,nx
	        do iz=0,nz
            ! front
            temp(3,ix,0,iz)=temp(3,ix,ny+1,iz)
	          temp(7,ix,0,iz)=temp(7,ix,ny+1,iz)
	          temp(9,ix,0,iz)=temp(9,ix,ny+1,iz)
	          temp(12,ix,0,iz)=temp(12,ix,ny+1,iz)
	          temp(14,ix,0,iz)=temp(14,ix,ny+1,iz)

            ! back
	          temp(4,ix,ny,iz)=temp(4,ix,-1,iz)
	          temp(8,ix,ny,iz)=temp(8,ix,-1,iz)
	          temp(10,ix,ny,iz)=temp(10,ix,-1,iz)
	          temp(11,ix,ny,iz)=temp(11,ix,-1,iz)
	          temp(13,ix,ny,iz)=temp(13,ix,-1,iz)
          enddo
          enddo
          do ix=0,nx
          do iy=0,ny
            ! bottom
            temp(5,ix,iy,0)=temp(5,ix,iy,nz+1)
	          temp(7,ix,iy,0)=temp(7,ix,iy,nz+1)
	          temp(10,ix,iy,0)=temp(10,ix,iy,nz+1)
	          temp(11,ix,iy,0)=temp(11,ix,iy,nz+1)
	          temp(14,ix,iy,0)=temp(14,ix,iy,nz+1)

            ! top
	          temp(6,ix,iy,nz)=temp(6,ix,iy,-1)
	          temp(8,ix,iy,nz)=temp(8,ix,iy,-1)
	          temp(9,ix,iy,nz)=temp(9,ix,iy,-1)
	          temp(12,ix,iy,nz)=temp(12,ix,iy,-1)
	          temp(13,ix,iy,nz)=temp(13,ix,iy,-1)
          end do
          end do
! ========          
!            call nodal_velocity(0,0,0,nx,ny,nz,ux,uy,uz,den)
!            write(*,*) "BC1 before edge: ",istep, ux,uy,uz,den
!            write(*,*) temp(:,0,0,0)
!========    
          do ix=0,nx
            ! bottom front edge
            temp(7,ix,0,0)=temp(7,ix,ny+1,nz+1)
            temp(14,ix,0,0)=temp(14,ix,ny+1,nz+1)

            ! temp(3,ix,ny,0)=temp(3,ix,-1,0)
            ! temp(9,ix,ny,0)=temp(9,ix,-1,0)
            ! temp(12,ix,ny,0)=temp(12,ix,-1,0)

            ! temp(5,ix,ny,0)=temp(5,ix,ny,nz+1)
            ! temp(10,ix,ny,0)=temp(10,ix,ny,nz+1)
            ! temp(11,ix,ny,0)=temp(11,ix,ny,nz+1) ! from top
            

            ! top front edge
            temp(9,ix,0,nz)=temp(9,ix,ny+1,-1)
            temp(12,ix,0,nz)=temp(12,ix,ny+1,-1)

            ! temp(3,ix,ny,nz)=temp(3,ix,-1,nz)
            ! temp(7,ix,ny,nz)=temp(7,ix,-1,nz) ! back
            ! temp(14,ix,ny,nz)=temp(14,ix,-1,nz)

            ! temp(6,ix,ny,nz)=temp(6,ix,ny,-1)
            ! temp(8,ix,ny,nz)=temp(8,ix,ny,-1)
            ! temp(13,ix,ny,nz)=temp(13,ix,ny,-1) ! bottom
            

            ! bottom back edge
            temp(10,ix,ny,0)=temp(10,ix,-1,nz+1)
            temp(11,ix,ny,0)=temp(11,ix,-1,nz+1)

            ! temp(4,ix,0,0)=temp(4,ix,ny+1,0)
            ! temp(8,ix,0,0)=temp(8,ix,ny+1,0)
            ! temp(13,ix,0,0)=temp(13,ix,ny+1,0) ! front 

            ! temp(5,ix,0,0)=temp(5,ix,0,nz+1)
            ! temp(7,ix,0,0)=temp(7,ix,0,nz+1)
            ! temp(14,ix,0,0)=temp(14,ix,0,nz+1) ! top


            ! top back edge
            temp(8,ix,ny,nz)=temp(8,ix,-1,-1)
            temp(13,ix,ny,nz)=temp(13,ix,-1,-1)

            ! temp(4,ix,0,nz)=temp(4,ix,ny+1,nz)
            ! temp(10,ix,0,nz)=temp(10,ix,ny+1,nz) ! front
            ! temp(11,ix,0,nz)=temp(11,ix,ny+1,nz)

            ! temp(6,ix,0,nz)=temp(6,ix,0,-1)
            ! temp(9,ix,0,nz)=temp(9,ix,0,-1) ! bottom
            ! temp(12,ix,0,nz)=temp(12,ix,0,-1)
          enddo
!=========
            ! call nodal_velocity(0,0,0,nx,ny,nz,ux,uy,uz,den)
            ! write(*,*) "BC1 after edge: ",istep, ux,uy,uz,den
            ! write(*,*) temp(:,0,0,0)
!=========
        end if

        ! print *, istep, temp(0:14,nx/2,ny/2,0)

        ux_in=bc(1)
        ux_out=bc(2)
        ! print *, ux_in, ux_out
!.......inlet
        do iz=-1,nz+1
          do iy=-1,ny+1
        din2=-(temp(0,0,iy,iz)+temp(3,0,iy,iz)+temp(4,0,iy,iz) &
             +temp(5,0,iy,iz)+temp(6,0,iy,iz)+2.d0*(temp(2,0,iy,iz) &
             +temp(8,0,iy,iz)+temp(10,0,iy,iz)+temp(12,0,iy,iz) &
             +temp(14,0,iy,iz)))/(ux_in-1.d0)
        tmp=1.d0/12.d0*din2*ux_in
	      temp(1,0,iy,iz)=temp(2,0,iy,iz)+8.d0*tmp
	      tmp1=temp(3,0,iy,iz)-temp(4,0,iy,iz)
	      tmp2=temp(5,0,iy,iz)-temp(6,0,iy,iz)
        ! if (iy == 0 .and. iz == 0) print *, "compare to temp13", 0, iy, iz, temp(14,0,iy,iz), tmp, -0.25d0*(-tmp1-tmp2)
	      temp(7,0,iy,iz)=temp(8,0,iy,iz)+tmp-0.25d0*(tmp1+tmp2)
	      temp(9,0,iy,iz)=temp(10,0,iy,iz)+tmp-0.25d0*(tmp1-tmp2)
	      temp(11,0,iy,iz)=temp(12,0,iy,iz)+tmp-0.25d0*(-tmp1+tmp2)
	      temp(13,0,iy,iz)=temp(14,0,iy,iz)+tmp-0.25d0*(-tmp1-tmp2)
          enddo
        enddo
!.......outlet
        do iz=-1,nz+1
          do iy=-1,ny+1
        dout2=(temp(0,nx,iy,iz)+temp(3,nx,iy,iz)+temp(4,nx,iy,iz) &
             +temp(5,nx,iy,iz)+temp(6,nx,iy,iz)+2.d0*(temp(1,nx,iy,iz) &
             +temp(7,nx,iy,iz)+temp(9,nx,iy,iz)+temp(11,nx,iy,iz) &
             +temp(13,nx,iy,iz)))/(ux_out+1.d0)
        tmp=1.d0/12.d0*dout2*ux_out
	      temp(2,nx,iy,iz)=temp(1,nx,iy,iz)-8.d0*tmp
	      tmp1=temp(3,nx,iy,iz)-temp(4,nx,iy,iz)
	      tmp2=temp(5,nx,iy,iz)-temp(6,nx,iy,iz)
	      temp(8,nx,iy,iz)=temp(7,nx,iy,iz)-tmp-0.25d0*(-tmp1-tmp2)
	      temp(10,nx,iy,iz)=temp(9,nx,iy,iz)-tmp-0.25d0*(-tmp1+tmp2)
	      temp(12,nx,iy,iz)=temp(11,nx,iy,iz)-tmp-0.25d0*(tmp1-tmp2)
	      temp(14,nx,iy,iz)=temp(13,nx,iy,iz)-tmp-0.25d0*(tmp1+tmp2)

        if (istep==1 .and. periodic_flag)then
           

        endif
                  ! if(istep==1 .and. periodic_flag)then
          !   ux=bc(1)
          !   uy=0.0d0
          !   uz=0.0d0
          ! endif
 !         write(*,*) istep, ix,iy,iz, ux,uy,uz
          enddo
        end do
!========
            ! call nodal_velocity(0,0,0,nx,ny,nz,ux,uy,uz,den)
            ! write(*,*) "Ver == after velocity treatment: ",istep, ux,uy,uz,den
            ! write(*,*) temp(:,0,0,0)
            ! call nodal_velocity(0,0,1,nx,ny,nz,ux,uy,uz,den)
            ! write(*,*) "Not == after velocity treatment: ",istep, ux,uy,uz,den
            ! write(*,*) temp(:,0,0,1)
            ! write(*,*) "--------------------"
!========
      end if ! bctype.eq.3
!
!.....Specified densities
!       if(bctype.eq.2)then
! 	      din=bc(1)
! 	      dout=bc(2)
! !.......inlet
!         do iz=0,nz
!           do iy=0,ny
! 	      ux=1.d0-(temp(0,0,iy,iz)+temp(3,0,iy,iz)+temp(4,0,iy,iz) &
!              +temp(5,0,iy,iz)+temp(6,0,iy,iz)+2.d0*(temp(2,0,iy,iz) &
!              +temp(8,0,iy,iz)+temp(10,0,iy,iz)+temp(12,0,iy,iz) &
!              +temp(14,0,iy,iz)))/din
!         tmp=1.d0/12.d0*din*ux
! 	      temp(1,0,iy,iz)=temp(2,0,iy,iz)+8.d0*tmp
! 	      tmp1=temp(3,0,iy,iz)-temp(4,0,iy,iz)
! 	      tmp2=temp(5,0,iy,iz)-temp(6,0,iy,iz)
! 	      temp(7,0,iy,iz)=temp(8,0,iy,iz)+tmp-0.25d0*(tmp1+tmp2)
! 	      temp(9,0,iy,iz)=temp(10,0,iy,iz)+tmp-0.25d0*(tmp1-tmp2)
! 	      temp(11,0,iy,iz)=temp(12,0,iy,iz)+tmp-0.25d0*(-tmp1+tmp2)
! 	      temp(13,0,iy,iz)=temp(14,0,iy,iz)+tmp-0.25d0*(-tmp1-tmp2)
!           enddo
!         enddo
! !.......outlet
!         do iz=0,nz
!           do iy=0,ny
! 	      ux=-1.d0+(temp(0,nx,iy,iz)+temp(3,nx,iy,iz)+temp(4,nx,iy,iz) &
!              +temp(5,nx,iy,iz)+temp(6,nx,iy,iz)+2.d0*(temp(1,nx,iy,iz) &
!              +temp(7,nx,iy,iz)+temp(9,nx,iy,iz)+temp(11,nx,iy,iz) &
!              +temp(13,nx,iy,iz)))/dout
!             tmp=1.d0/12.d0*dout*ux
! 	      temp(2,nx,iy,iz)=temp(1,nx,iy,iz)-8.d0*tmp
! 	      tmp1=temp(3,nx,iy,iz)-temp(4,nx,iy,iz)
! 	      tmp2=temp(5,nx,iy,iz)-temp(6,nx,iy,iz)
! 	      temp(8,nx,iy,iz)=temp(7,nx,iy,iz)-tmp-0.25d0*(-tmp1-tmp2)
! 	      temp(10,nx,iy,iz)=temp(9,nx,iy,iz)-tmp-0.25d0*(-tmp1+tmp2)
! 	      temp(12,nx,iy,iz)=temp(11,nx,iy,iz)-tmp-0.25d0*(tmp1-tmp2)
! 	      temp(14,nx,iy,iz)=temp(13,nx,iy,iz)-tmp-0.25d0*(tmp1+tmp2)
!           enddo
!         end do
! !.......Special treatment for corner nodes
! !       left bottom inlet segment
!         do iy=0,ny
!           temp(1,0,iy,0)=temp(2,0,iy,0)
! 	    temp(5,0,iy,0)=temp(6,0,iy,0)
! 	    tmp1=temp(3,0,iy,0)-temp(4,0,iy,0)
! 	    temp(7,0,iy,0)=temp(8,0,iy,0)-0.5d0*tmp1
! 	    temp(11,0,iy,0)=temp(12,0,iy,0)+0.5d0*tmp1
! 	    temp(9,0,iy,0)=0.25d0*(din-temp(0,0,iy,0)-temp(1,0,iy,0) &
!                        -temp(2,0,iy,0)-temp(3,0,iy,0)-temp(4,0,iy,0) &
!                        -temp(5,0,iy,0)-temp(6,0,iy,0)-temp(7,0,iy,0) &
!                        -temp(8,0,iy,0)-temp(11,0,iy,0)-temp(12,0,iy,0))
! 	    temp(10,0,iy,0)=temp(9,0,iy,0)
! 	    temp(13,0,iy,0)=temp(9,0,iy,0)
! 	    temp(14,0,iy,0)=temp(9,0,iy,0)
!          end do
! !       left top inlet segment
!         do iy=0,ny
!           temp(1,0,iy,nz)=temp(2,0,iy,nz)
! 	    temp(6,0,iy,nz)=temp(5,0,iy,nz)
! 	    tmp1=temp(3,0,iy,nz)-temp(4,0,iy,nz)
! 	    temp(9,0,iy,nz)=temp(10,0,iy,nz)-0.5d0*tmp1
! 	    temp(13,0,iy,nz)=temp(14,0,iy,nz)+0.5d0*tmp1
! 	    temp(7,0,iy,nz)=0.25d0*(din-temp(0,0,iy,nz)-temp(1,0,iy,nz) &
!                    -temp(2,0,iy,nz)-temp(3,0,iy,nz)-temp(4,0,iy,nz) &
!                    -temp(5,0,iy,nz)-temp(6,0,iy,nz)-temp(9,0,iy,nz) &
!                    -temp(10,0,iy,nz)-temp(13,0,iy,nz)-temp(14,0,iy,nz))
! 	    temp(8,0,iy,nz)=temp(7,0,iy,nz)
! 	    temp(11,0,iy,nz)=temp(7,0,iy,nz)
! 	    temp(12,0,iy,nz)=temp(7,0,iy,nz)
!         enddo
! !       right bottom outlet segment
!         do iy=0,ny
!           temp(2,nx,iy,0)=temp(1,nx,iy,0)
! 	    temp(5,nx,iy,0)=temp(6,nx,iy,0)
! 	    tmp1=temp(3,nx,iy,0)-temp(4,nx,iy,0)
! 	    temp(10,nx,iy,0)=temp(9,nx,iy,0)+0.5d0*tmp1
! 	    temp(14,nx,iy,0)=temp(13,nx,iy,0)-0.5d0*tmp1
! 	    temp(7,nx,iy,0)=0.25d0*(dout-temp(0,nx,iy,0)-temp(1,nx,iy,0) &
!                    -temp(2,nx,iy,0)-temp(3,nx,iy,0)-temp(4,nx,iy,0) &
!                    -temp(5,nx,iy,0)-temp(6,nx,iy,0)-temp(9,nx,iy,0) &
!                    -temp(10,nx,iy,0)-temp(13,nx,iy,0)-temp(14,nx,iy,0))
! 	    temp(8,nx,iy,0)=temp(7,nx,iy,0)
! 	    temp(11,nx,iy,0)=temp(7,nx,iy,0)
! 	    temp(12,nx,iy,0)=temp(7,nx,iy,0)
!         enddo
! !       right top outlet segment
!         do iy=0,ny
!           temp(2,nx,iy,nz)=temp(1,nx,iy,nz)
! 	    temp(6,nx,iy,nz)=temp(5,nx,iy,nz)
! 	    tmp1=temp(3,nx,iy,nz)-temp(4,nx,iy,nz)
! 	    temp(8,nx,iy,nz)=temp(7,nx,iy,nz)+0.5d0*tmp1
! 	    temp(12,nx,iy,nz)=temp(11,nx,iy,nz)-0.5d0*tmp1
! 	   temp(9,nx,iy,nz)=0.25d0*(dout-temp(0,nx,iy,nz)-temp(1,nx,iy,nz) &
!                 -temp(2,nx,iy,nz)-temp(3,nx,iy,nz)-temp(4,nx,iy,nz) &
!                 -temp(5,nx,iy,nz)-temp(6,nx,iy,nz)-temp(7,nx,iy,nz) &
!                 -temp(8,nx,iy,nz)-temp(11,nx,iy,nz)-temp(12,nx,iy,nz))
! 	    temp(10,nx,iy,nz)=temp(9,nx,iy,nz)
! 	    temp(13,nx,iy,nz)=temp(9,nx,iy,nz)
! 	    temp(14,nx,iy,nz)=temp(9,nx,iy,nz)
!         enddo
! 	    endif
       return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine check_density(nx,ny,nz,istep,time)  
    use fluid,only:node,temp,yb
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      sum=0.d0
!.....Loop over rows
      do iz=0,nz
       do iy=0, ny       
        do ix=0, nx
!
          if(node(ix,iy,iz)%obst.eq.0.or.node(ix,iy,iz)%obst.eq.4)then
            do i=0, 14
              sum=sum+temp(i,ix,iy,iz)
            end do

              if(sum.lt.0.d0)then
		      write(*,*)'negative sum',istep,ix,iy,iz
		      stop
	          endif
          endif
        end do
       enddo
      end do
!      cpu_time=cputime()
!      write(*,*) cpu_time
!
      write(6,100) istep,time,sum
  100 format('** Step =',i7, '  Time =',g11.4, '  Density =',g12.5, ' **')
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Real(8) function cputime()
!      real trray(1)
!      real elatim

!      elatim=etime(trray)
!      cputime=trray(1)
!      return
!    End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine write_velocity(nx,ny,nz)
	use fluid,only:node,temp,yb,BODYF
    use solid,only:gacce
    use system,only:istep
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
    ix=nx/2!test
	  iy=ny/2
	  iz=nz/2
      call nodal_velocity(ix,iy,iz,nx,ny,nz,ux,uy,uz,den)
      write(13,100) istep,ix,iy,iz,ux,uy,uz
!          if(umax.lt.uvelocity) umax=uvelocity
!   	  write(13,100) ix,iy,iz,ux,uy,uz
  100     format(1x,i6,1x,i6,1x,i6,1x,i6,3(1x,g10.4))
      return
	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine comp_rey(visc,rad,umax,cc,dx,dt)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      umax=umax*cc
	  rey=2.d0*umax*rad*dx/visc
      write (6,*) '*** Calculations finished, results:'
      write (6,*) '***'
      write (6,*) '*** viscosity = ', visc
      write (6,*) '*** Max velocity = ', umax
      write (6,*) '*** Lattice speed = ', cc
      write (6,*) '*** Reynolds number = ', rey
      write (6,*) '***'
      write (12,*) '*** Reynolds number = ', rey
      write (12,*) '*** Max velocity = ', umax
      write (12,*) '*** Lattice speed = ', cc
      write (12,*) '*** Time step = ', dt
      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine generate_boundary_nodes
    use system,only:PI,p1,head_lpm
    use fs_inter,only:nbp
    use solid,only:bpoint,particle
    implicit none
      integer ip,j
      real(8) dangle

      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine update_boundary_nodes
    use system,only:PI,p1,head_lpm
    use fs_inter,only:nbp
    use solid,only:bpoint,particle
    implicit none
      integer ip,j

      return
    End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(8) function dirac(r)
    use system,only:PI
    implicit none
      real(8) r
      if(abs(r).gt.2.0d0)then
        dirac=0.0d0
      else
        dirac=0.25*(1.0+cos(PI*abs(r)/2.0))
      end if
      return
    End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*******************************************
      integer function opposite_direction(i)
!*******************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!	            
	select case (i)
        case (0)
	    j=0
	  case (1)
	    j=2
	  case (2)
	    j=1
	  case (3)
	    j=4
	  case (4)
	    j=3
	  case (5)
	    j=6
	  case (6)
	    j=5
	  case (7)
	    j=8
	  case (8)
	    j=7
	  case (9)
	    j=10
	  case (10)
	    j=9
	  case (11)
	    j=12
	  case (12)
	    j=11
	  case (13)
	    j=14
	  case (14)
	    j=13
	end select
	opposite_direction=j
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Subroutine sphere_output_temp
!    use solid
!    implicit none
 !     integer i,j
!	  real*8 sphere_mesh(3,nodes)
!	  integer sphere_elem(4,nelem)
!	  read(19,*)(sphere_mesh(1,i),i=1,nodes)
!	  read(19,*)(sphere_mesh(2,i),i=1,nodes)
!	  read(19,*)(sphere_mesh(3,i),i=1,nodes)
!	  do i=1,nelem
!	   read(19,*)(sphere_elem(j,i),j=1,4)
!	  enddo
!	return
!	End subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine body_force_density(ix,iy,iz,ux,uy,uz,den)
    use fluid, only : rtao,fi_body,w
    use solid, only : gacce
    implicit none

    integer ix,iy,iz,i
    real(8) den,ux,uy,uz
    real(8) eiu(0:14),ei_u_f(0:14),eif(0:14)

!.....calculate Ei dot U
    eiu(0)  = 0.0d0

!.....axis directions
    eiu(1)  =  ux
    eiu(2)  = -ux
    eiu(3)  =  uy
    eiu(4)  = -uy
    eiu(5)  =  uz
    eiu(6)  = -uz

!.....diagonal directions
    eiu(7)  =  ux + uy + uz
    eiu(8)  = -ux - uy - uz
    eiu(9)  =  ux + uy - uz
    eiu(10) = -ux - uy + uz
    eiu(11) =  ux - uy + uz
    eiu(12) = -ux + uy - uz
    eiu(13) =  ux - uy - uz
    eiu(14) = -ux + uy + uz

!.....calculate (Ei-U) dot F
!.....F = (0,-den*gacce,0)

!.....directions with Ey = 0
    ei_u_f(0) = den*uy*gacce
    ei_u_f(1) = den*uy*gacce
    ei_u_f(2) = den*uy*gacce
    ei_u_f(5) = den*uy*gacce
    ei_u_f(6) = den*uy*gacce

!.....directions with Ey = +1
    ei_u_f(3)  = (1.0d0-uy)*(-den*gacce)
    ei_u_f(7)  = (1.0d0-uy)*(-den*gacce)
    ei_u_f(9)  = (1.0d0-uy)*(-den*gacce)
    ei_u_f(12) = (1.0d0-uy)*(-den*gacce)
    ei_u_f(14) = (1.0d0-uy)*(-den*gacce)

!.....directions with Ey = -1
    ei_u_f(4)  = (1.0d0+uy)*(den*gacce)
    ei_u_f(8)  = (1.0d0+uy)*(den*gacce)
    ei_u_f(10) = (1.0d0+uy)*(den*gacce)
    ei_u_f(11) = (1.0d0+uy)*(den*gacce)
    ei_u_f(13) = (1.0d0+uy)*(den*gacce)

!.....calculate Ei dot F
    eif(0) = 0.0d0

!.....directions with Ey = 0
    eif(1) = 0.0d0
    eif(2) = 0.0d0
    eif(5) = 0.0d0
    eif(6) = 0.0d0

!.....directions with Ey = +1
    eif(3)  = -den*gacce
    eif(7)  = -den*gacce
    eif(9)  = -den*gacce
    eif(12) = -den*gacce
    eif(14) = -den*gacce

!.....directions with Ey = -1
    eif(4)  = den*gacce
    eif(8)  = den*gacce
    eif(10) = den*gacce
    eif(11) = den*gacce
    eif(13) = den*gacce

!.....calculate Guo body-force contribution
!.....For D3Q15: cs^2 = 1/3, so 1/cs^2 = 3 and 1/cs^4 = 9
    do i=0,14
        fi_body(i) = w(i)*(1.0d0-0.5d0*rtao) * &
                     (3.0d0*ei_u_f(i) + 9.0d0*eiu(i)*eif(i))
    end do

    return
End subroutine body_force_density
