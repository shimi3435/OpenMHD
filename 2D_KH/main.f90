program main
!-----------------------------------------------------------------------
!     OpenMHD  Riemann solver (serial version)
!-----------------------------------------------------------------------
!     2010/09/27  S. Zenitani  K-H instability
!-----------------------------------------------------------------------
  use run_config_mod, only: sim_config, load_run_config, set_output_dir, ensure_output_directory
  implicit none
  include 'param.h'
  integer :: ix, jx
  integer :: loop_max
  real(8) :: tend
  real(8) :: dtout
  real(8) :: cfl
! Slope limiter, numerical flux, and time integration settings
  integer :: lm_type
  integer :: flux_type
  integer :: time_type
!-----------------------------------------------------------------------
! See also model.f90
!-----------------------------------------------------------------------
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
  character(len=256) :: params_arg, out_arg
!-----------------------------------------------------------------------
  real(8), allocatable :: x(:), y(:)
  real(8) :: dx
  real(8), allocatable :: U(:,:,:)
  real(8), allocatable :: U1(:,:,:)
  real(8), allocatable :: V(:,:,:)
  real(8), allocatable :: VL(:,:,:), VR(:,:,:)
  real(8), allocatable :: F(:,:,:), G(:,:,:)
!-----------------------------------------------------------------------

  call get_command_argument(1, params_arg)
  if (len_trim(params_arg) > 0) then
     call load_run_config(trim(params_arg))
  else
     call load_run_config()
  endif
  call get_command_argument(2, out_arg)
  if (len_trim(out_arg) > 0) call set_output_dir(trim(out_arg))
  call ensure_output_directory()

  ix = sim_config%nx + 2
  jx = sim_config%ny + 2
  loop_max = sim_config%loop_max
  tend     = sim_config%tend
  dtout    = sim_config%dtout
  cfl      = sim_config%cfl
  lm_type  = sim_config%lm_type
  flux_type = sim_config%flux_type
  time_type = sim_config%time_type

  allocate(x(ix), y(jx))
  allocate(U(ix,jx,var1))
  allocate(U1(ix,jx,var1))
  allocate(V(ix,jx,var2))
  allocate(VL(ix,jx,var1), VR(ix,jx,var1))
  allocate(F(ix,jx,var1), G(ix,jx,var1))

  t    =  0.d0
  dt   =  0.d0
  call model(U,V,x,y,dx,ix,jx)
  call bc_for_U(U,ix,jx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
  t_output = -dt/3.d0
  n_output =  0

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  write(6,*) '[Params]'
  write(6,998) dt, dtout, ix, jx
  write(6,999) lm_type, flux_type, time_type
998 format (' dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:', i5, i5 )
999 format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
  write(6,*) '== start =='

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     write(6,*) ' t = ', t
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------
!    [ output ]
     if ( t >= t_output ) then
        write(6,*) 'data output   t = ', t
        write(filename,'(A,"/field-",i5.5,".dat")') trim(sim_config%output_dir), n_output
        call fileio_output(filename,ix,jx,t,x,y,U,V)
        n_output = n_output + 1
        t_output = t_output + dtout
     endif
!    [ end? ]
     if ( t >= tend )  exit
     if ( n_loop >= loop_max ) then
        write(6,*) 'max loop'
        exit
     endif
!   -----------------
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
!    GLM solver for the first half timestep
!    This should be done after set_dt()
!     write(6,*) 'U --> SS'
     call glm_ss2(U,ch,dt,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,1,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,1,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(U(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    fix VL/VR for periodic bc (F)
     call bc_for_F(VL,VR,ix,jx)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,1,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,2,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,2,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,2,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,2,lm_type)
     call limiter(U(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,2,lm_type)
     call limiter(U(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,2,lm_type)
     call limiter(U(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,2,lm_type)
     call limiter(U(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,2,lm_type)
     call limiter(U(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,2,lm_type)
!    fix VL/VR for wall bc (G)
     call bc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,2,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U* = U + (dt/dx) (F-F)'
        call rk_tvd21(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U*(n+1/2) = U + (0.5 dt/dx) (F-F)'
        call rk_std21(U,U1,F,G,dt,dx,ix,jx)
     endif
!    boundary condition
     call bc_for_U(U1,ix,jx)
!     write(6,*) 'U* --> V'
     call u2v(U1,V,ix,jx)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (F)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,1,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,1,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,1,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,1,lm_type)
     call limiter(U1(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,1,lm_type)
     call limiter(U1(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,1,lm_type)
     call limiter(U1(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,1,lm_type)
     call limiter(U1(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,1,lm_type)
     call limiter(U1(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,1,lm_type)
!    fix VL/VR for periodic bc (F)
     call bc_for_F(VL,VR,ix,jx)
!    Numerical flux in the X direction (F)
!     write(6,*) 'VL, VR --> F'
     call flux_solver(F,VL,VR,ch,ix,jx,1,flux_type)

!    Slope limiters on primitive variables
!     write(6,*) 'V --> VL, VR (G)'
     call limiter(V(:,:,vx),VL(:,:,vx),VR(:,:,vx),ix,jx,2,lm_type)
     call limiter(V(:,:,vy),VL(:,:,vy),VR(:,:,vy),ix,jx,2,lm_type)
     call limiter(V(:,:,vz),VL(:,:,vz),VR(:,:,vz),ix,jx,2,lm_type)
     call limiter(V(:,:,pr),VL(:,:,pr),VR(:,:,pr),ix,jx,2,lm_type)
     call limiter(U1(:,:,ro),VL(:,:,ro),VR(:,:,ro),ix,jx,2,lm_type)
     call limiter(U1(:,:,bx),VL(:,:,bx),VR(:,:,bx),ix,jx,2,lm_type)
     call limiter(U1(:,:,by),VL(:,:,by),VR(:,:,by),ix,jx,2,lm_type)
     call limiter(U1(:,:,bz),VL(:,:,bz),VR(:,:,bz),ix,jx,2,lm_type)
     call limiter(U1(:,:,ps),VL(:,:,ps),VR(:,:,ps),ix,jx,2,lm_type)
!    fix VL/VR for wall bc (G)
     call bc_for_G(VL,VR,ix,jx)
!    Numerical flux in the Y direction (G)
!     write(6,*) 'VL, VR --> G'
     call flux_solver(G,VL,VR,ch,ix,jx,2,flux_type)

     if( time_type == 0 ) then
!       write(6,*) 'U_new = 0.5( U_old + U* + F dt )'
        call rk_tvd22(U,U1,F,G,dt,dx,ix,jx)
     elseif( time_type == 1 ) then
!       write(6,*) 'U_new = U + (dt/dx) (F-F) (n+1/2)'
        call rk_std22(U,F,G,dt,dx,ix,jx)
     endif
!    boundary condition
     call bc_for_U(U,ix,jx)

!    GLM solver for the second half timestep
     call glm_ss2(U,ch,dt,ix,jx)

     t=t+dt

  enddo
!-----------------------------------------------------------------------

  write(6,*) '== end =='
end program main
!-----------------------------------------------------------------------
