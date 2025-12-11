program main
!-----------------------------------------------------------------------
!     OpenMHD  Riemann solver (parallel version)
!-----------------------------------------------------------------------
!     2010/09/27  S. Zenitani  K-H instability (MPI version)
!     2015/04/05  S. Zenitani  MPI-IO
!     2018/05/02  S. Zenitani  parallel module (MPI-3 version)
!-----------------------------------------------------------------------
  use run_config_mod, only: sim_config, load_run_config, broadcast_run_config, &
       set_output_dir, ensure_output_directory
  use mpi
  use parallel
  implicit none
  include 'param.h'
  integer :: ix, jx
  integer :: mpi_nums(2)
  logical :: bc_periodicity(2)
  integer :: loop_max
  real(8) :: tend
  real(8) :: dtout
  real(8) :: cfl
! If non-zero, restart from a previous file. If negative, find a restart file backword in time.
  integer :: n_start = 0
! Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  integer :: lm_type
! Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  integer :: flux_type
! Time marching  (0: TVD RK2, 1: RK2)
  integer :: time_type
! File I/O  (0: Standard, 1: MPI-IO)
  integer :: io_type
!-----------------------------------------------------------------------
! See also modelp.f90
!-----------------------------------------------------------------------
  integer :: k
  integer :: n_loop,n_output
  real(8) :: t, dt, t_output
  real(8) :: ch
  character*256 :: filename
  integer :: merr, myrank, mreq(2)   ! for MPI
  integer :: world_rank
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
  call get_command_argument(2, out_arg)

  call mpi_init(merr)
  call mpi_comm_rank(MPI_COMM_WORLD, world_rank, merr)
  if (world_rank == 0) then
     if (len_trim(params_arg) > 0) then
        call load_run_config(trim(params_arg))
     else
        call load_run_config()
     endif
     if (len_trim(out_arg) > 0) call set_output_dir(trim(out_arg))
     call ensure_output_directory()
  endif
  call broadcast_run_config(MPI_COMM_WORLD, 0)
  call mpi_barrier(MPI_COMM_WORLD, merr)

  ix = sim_config%nx + 2
  jx = sim_config%ny + 2
  mpi_nums = sim_config%mpi_nums
  bc_periodicity = sim_config%bc_periodicity
  loop_max = sim_config%loop_max
  tend     = sim_config%tend
  dtout    = sim_config%dtout
  cfl      = sim_config%cfl
  lm_type  = sim_config%lm_type
  flux_type = sim_config%flux_type
  time_type = sim_config%time_type
  io_type  = sim_config%io_type
  n_start  = sim_config%n_start

  allocate(x(ix), y(jx))
  allocate(U(ix,jx,var1))
  allocate(U1(ix,jx,var1))
  allocate(V(ix,jx,var2))
  allocate(VL(ix,jx,var1), VR(ix,jx,var1))
  allocate(F(ix,jx,var1), G(ix,jx,var1))

!-----------------------------------------------------------------------
! for MPI
  call parallel_init(mpi_nums,bc_periodicity,ix,jx)
  myrank = ranks%myrank
!-----------------------------------------------------------------------

  t    =  0.d0
  dt   =  0.d0
  call modelp(U,V,x,y,dx,ix,jx)
! boundary conditions
  call parallel_exchange(U,ix,jx,1)
  call parallel_exchange(U,ix,jx,2)
  call mpibc_for_U(U,ix,jx)
  call set_dt(U,V,ch,dt,dx,cfl,ix,jx)

  call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
  call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
  call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)
  call mpi_barrier(cart2d%comm,merr)

  if ( dt > dtout ) then
     write(6,*) 'error: ', dt, '>', dtout
     stop
  endif
  if( myrank == 0 ) then
     write(6,*) '[Params]'
     write(6,997) version, ranks%size, cart2d%sizes
     write(6,998) dt, dtout, cart2d%sizes(1)*(ix-2)+2, ix, cart2d%sizes(2)*(jx-2)+2, jx
     write(6,999) lm_type, flux_type, time_type
997  format ('Code version: ', i8, '  MPI node # : ', i5,' (',i4,' x ',i4,')' )
998  format (' dt:', 1p, e10.3, ' dtout:', 1p, e10.3, ' grids:',i6,' (',i5,') x ',i6,' (',i5,') ')
999  format (' limiter: ', i1, '  flux: ', i1, '  time-marching: ', i1 )
     write(6,*) '== start =='
  endif

  ! If n_start is negative, look for a latest restart file.
  if ( n_start < 0 ) then
     if ( myrank == 0 ) then
        do k = floor(tend/dtout),0,-1
           n_start = k
           if ( io_type == 0 ) then
              write(filename,'(A,"/field-rank",i5.5,"-",i5.5,".dat")') &
                   trim(sim_config%output_dir), myrank, n_start
           else
              write(filename,'(A,"/field-",i5.5,".dat")') trim(sim_config%output_dir), n_start
           endif
           open(15,file=filename,form='unformatted',access='stream',status='old',err=100)
           close(15)
           exit
100        continue
        enddo
     endif
     call mpi_bcast(n_start,1,mpi_integer,0,cart2d%comm,merr)
  endif
  call mpi_barrier(cart2d%comm,merr)

  ! If n_start is non-zero, restart from a previous file.
  if ( n_start == 0 ) then
     n_output = 0
     t_output = -dt/3.d0
  else
           if ( io_type == 0 ) then
              write(6,*) 'reading data ...   rank = ', myrank
              write(filename,'(A,"/field-rank",i5.5,"-",i5.5,".dat")') &
                   trim(sim_config%output_dir), myrank, n_start
              call fileio_input(filename,ix,jx,t,x,y,U)
           else
              if( myrank == 0 )  write(6,*) 'reading data ...'
              write(filename,'(A,"/field-",i5.5,".dat")') trim(sim_config%output_dir), n_start
              call mpiio_input(filename,ix,jx,t,x,y,U)
           endif
     n_output = n_start + 1
     t_output = t + dtout
  endif

!-----------------------------------------------------------------------
  do n_loop=1,loop_max

     if( myrank == 0 ) then
        write(6,*) ' t = ', t
     endif
!    Recovering primitive variables
!     write(6,*) 'U --> V'
     call u2v(U,V,ix,jx)
!   -----------------
!    [ output ]
     if ( t >= t_output ) then
        if (( n_loop > 1 ).or.( n_start == 0 )) then
           if ( io_type == 0 ) then
              write(6,*) 'writing data ...   t = ', t, ' rank = ', myrank
              write(filename,'(A,"/field-rank",i5.5,"-",i5.5,".dat")') &
                   trim(sim_config%output_dir), myrank, n_output
              call fileio_output(filename,ix,jx,t,x,y,U,V)
           else
              if( myrank == 0 )  write(6,*) 'writing data ...   t = ', t
              write(filename,'(A,"/field-",i5.5,".dat")') trim(sim_config%output_dir), n_output
              call mpiio_output(filename,ix,jx,t,x,y,U,V)
           endif
        endif
        n_output = n_output + 1
        t_output = t_output + dtout
        call mpi_barrier(cart2d%comm,merr)
     endif
!    [ end? ]
     if ( t >= tend )  exit
     if ( n_loop >= loop_max ) then
        if( myrank == 0 )  write(6,*) 'max loop'
        exit
     endif
!   -----------------
!    CFL condition
     call set_dt(U,V,ch,dt,dx,cfl,ix,jx)
     call mpi_iallreduce(mpi_in_place,ch,1,mpi_real8,mpi_max,cart2d%comm,mreq(1),merr)
     call mpi_iallreduce(mpi_in_place,dt,1,mpi_real8,mpi_min,cart2d%comm,mreq(2),merr)
     call mpi_waitall(2,mreq,mpi_statuses_ignore,merr)

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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (F)'
     call parallel_exchange2(VL,VR,ix,jx,1)
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (G)'
     call parallel_exchange2(VL,VR,ix,jx,2)
     call mpibc_for_G(VL,VR,ix,jx)
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
!    boundary conditions
     call parallel_exchange(U1,ix,jx,1)
     call parallel_exchange(U1,ix,jx,2)
     call mpibc_for_U(U1,ix,jx)
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (F)'
     call parallel_exchange2(VL,VR,ix,jx,1)
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
!    Interpolated primitive variables at MPI boundaries
!     write(6,*) 'adjusting VL/VR (G)'
     call parallel_exchange2(VL,VR,ix,jx,2)
     call mpibc_for_G(VL,VR,ix,jx)
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
!    boundary conditions
     call parallel_exchange(U,ix,jx,1)
     call parallel_exchange(U,ix,jx,2)
     call mpibc_for_U(U,ix,jx)

!    GLM solver for the second half timestep
     call glm_ss2(U,ch,dt,ix,jx)

     t=t+dt

  enddo
!-----------------------------------------------------------------------

  call parallel_finalize()
  if( myrank == 0 )  write(6,*) '== end =='

end program main
!-----------------------------------------------------------------------
