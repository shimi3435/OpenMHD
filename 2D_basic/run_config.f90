module run_config_mod
  implicit none

  type :: run_config_t
     integer :: nx = 200                 ! interior cells in X
     integer :: ny = 200                 ! interior cells in Y
     integer :: loop_max = 30000
     real(8) :: tend  = 4.0d0
     real(8) :: dtout = 0.1d0
     real(8) :: cfl   = 0.4d0
     integer :: lm_type   = 1
     integer :: flux_type = 3
     integer :: time_type = 0
     real(8) :: domain_x(2)  = (/0.d0, 6.2831853071795864769d0/)
     real(8) :: domain_y_min = 0.d0
     real(8) :: density0      = 2.7777777777777777d0
     real(8) :: pressure0     = 1.6666666666666667d0
     real(8) :: vx_amplitude  = -1.d0
     real(8) :: vx_wavenumber = 1.d0
     real(8) :: vy_amplitude  = 1.d0
     real(8) :: vy_wavenumber = 1.d0
     real(8) :: bx_amplitude  = -1.d0
     real(8) :: bx_wavenumber = 1.d0
     real(8) :: by_amplitude  = 1.d0
     real(8) :: by_wavenumber = 2.d0
     real(8) :: bz0 = 0.d0
     real(8) :: ps0 = 0.d0
     integer :: mpi_nums(2)       = (/2, 2/)
     logical :: bc_periodicity(2) = (/.true., .true./)
     integer :: io_type = 1
     integer :: n_start = 0
     character(len=256) :: param_filename = 'params.nml'
     character(len=256) :: output_dir     = 'data'
  end type run_config_t

  type(run_config_t) :: sim_config

contains

  subroutine load_run_config(filename)
    character(len=*), intent(in), optional :: filename
    character(len=256) :: fname
    logical :: has_file
    integer :: ios
    integer :: nx, ny, loop_max, lm_type, flux_type, time_type
    integer :: mpi_nums(2), io_type, n_start
    real(8) :: tend, dtout, cfl, domain_x(2), domain_y_min
    real(8) :: density0, pressure0
    real(8) :: vx_amplitude, vx_wavenumber
    real(8) :: vy_amplitude, vy_wavenumber
    real(8) :: bx_amplitude, bx_wavenumber
    real(8) :: by_amplitude, by_wavenumber
    real(8) :: bz0, ps0
    logical :: bc_periodicity(2)
    character(len=256) :: output_dir

    namelist /run_params/ nx, ny, loop_max, tend, dtout, cfl, lm_type, flux_type, time_type, &
&       domain_x, domain_y_min, density0, pressure0, vx_amplitude, vx_wavenumber,            &
&       vy_amplitude, vy_wavenumber, bx_amplitude, bx_wavenumber, by_amplitude,              &
&       by_wavenumber, bz0, ps0, mpi_nums, bc_periodicity, io_type, n_start, output_dir

    fname = sim_config%param_filename
    if (present(filename)) then
       fname = filename
    endif

    nx             = sim_config%nx
    ny             = sim_config%ny
    loop_max       = sim_config%loop_max
    tend           = sim_config%tend
    dtout          = sim_config%dtout
    cfl            = sim_config%cfl
    lm_type        = sim_config%lm_type
    flux_type      = sim_config%flux_type
    time_type      = sim_config%time_type
    domain_x       = sim_config%domain_x
    domain_y_min   = sim_config%domain_y_min
    density0       = sim_config%density0
    pressure0      = sim_config%pressure0
    vx_amplitude   = sim_config%vx_amplitude
    vx_wavenumber  = sim_config%vx_wavenumber
    vy_amplitude   = sim_config%vy_amplitude
    vy_wavenumber  = sim_config%vy_wavenumber
    bx_amplitude   = sim_config%bx_amplitude
    bx_wavenumber  = sim_config%bx_wavenumber
    by_amplitude   = sim_config%by_amplitude
    by_wavenumber  = sim_config%by_wavenumber
    bz0            = sim_config%bz0
    ps0            = sim_config%ps0
    mpi_nums       = sim_config%mpi_nums
    bc_periodicity = sim_config%bc_periodicity
    io_type        = sim_config%io_type
    n_start        = sim_config%n_start
    output_dir     = sim_config%output_dir

    inquire(file=trim(fname), exist=has_file)
    if (.not. has_file) then
       write(6,*) 'run_config: "', trim(fname), '" not found. Using defaults.'
       sim_config%param_filename = trim(fname)
       return
    endif

    open(unit=99,file=trim(fname),status='old',action='read',iostat=ios)
    if (ios /= 0) then
       write(6,*) 'run_config: unable to open "', trim(fname), '". Using defaults.'
       sim_config%param_filename = trim(fname)
       return
    endif

    read(99, nml=run_params, iostat=ios)
    close(99)
    if (ios /= 0) then
        write(6,*) 'run_config: error reading "', trim(fname), '". Check syntax.'
        sim_config%param_filename = trim(fname)
        return
    endif

    sim_config%nx             = max(4, nx)
    sim_config%ny             = max(4, ny)
    sim_config%loop_max       = max(1, loop_max)
    sim_config%tend           = max(0.d0, tend)
    sim_config%dtout          = max(1.d-12, dtout)
    sim_config%cfl            = cfl
    sim_config%lm_type        = lm_type
    sim_config%flux_type      = flux_type
    sim_config%time_type      = time_type
    sim_config%domain_x       = domain_x
    sim_config%domain_y_min   = domain_y_min
    sim_config%density0       = density0
    sim_config%pressure0      = pressure0
    sim_config%vx_amplitude   = vx_amplitude
    sim_config%vx_wavenumber  = vx_wavenumber
    sim_config%vy_amplitude   = vy_amplitude
    sim_config%vy_wavenumber  = vy_wavenumber
    sim_config%bx_amplitude   = bx_amplitude
    sim_config%bx_wavenumber  = bx_wavenumber
    sim_config%by_amplitude   = by_amplitude
    sim_config%by_wavenumber  = by_wavenumber
    sim_config%bz0            = bz0
    sim_config%ps0            = ps0
    sim_config%mpi_nums       = mpi_nums
    sim_config%bc_periodicity = bc_periodicity
    sim_config%io_type        = io_type
    sim_config%n_start        = n_start
    sim_config%param_filename = trim(fname)
    if (len_trim(output_dir) > 0) then
       sim_config%output_dir = adjustl(output_dir)
    endif

    write(6,*) 'run_config: loaded parameters from "', trim(fname), '"'
  end subroutine load_run_config

  subroutine set_output_dir(dir)
    character(len=*), intent(in) :: dir
    if (len_trim(dir) > 0) then
      sim_config%output_dir = adjustl(dir)
    endif
  end subroutine set_output_dir

  subroutine ensure_output_directory()
    character(len=256) :: dir
    logical :: exists
    integer :: istat
    character(len=512) :: cmd

    dir = adjustl(sim_config%output_dir)
    if (len_trim(dir) == 0) dir = 'data'
    inquire(file=trim(dir), exist=exists)
    if (.not. exists) then
       cmd = 'mkdir -p ' // trim(dir)
       call execute_command_line(trim(cmd), exitstat=istat)
       if (istat /= 0) then
          write(6,*) 'run_config: warning - failed to create directory ', trim(dir)
       endif
    endif
  end subroutine ensure_output_directory

  subroutine broadcast_run_config(comm, root)
    use mpi
    integer, intent(in) :: comm, root
    integer :: ierr
    integer :: int_vals(10)
    real(8) :: real_vals(18)
    integer :: logic_vals(2)

    int_vals = (/ sim_config%nx, sim_config%ny, sim_config%loop_max, sim_config%lm_type, &
                  sim_config%flux_type, sim_config%time_type, sim_config%mpi_nums(1),     &
                  sim_config%mpi_nums(2), sim_config%io_type, sim_config%n_start /)
    real_vals = (/ sim_config%tend, sim_config%dtout, sim_config%cfl, sim_config%domain_x(1), &
                   sim_config%domain_x(2), sim_config%domain_y_min, sim_config%density0,      &
                   sim_config%pressure0, sim_config%vx_amplitude, sim_config%vx_wavenumber,   &
                   sim_config%vy_amplitude, sim_config%vy_wavenumber, sim_config%bx_amplitude,&
                   sim_config%bx_wavenumber, sim_config%by_amplitude, sim_config%by_wavenumber,&
                   sim_config%bz0, sim_config%ps0 /)
    logic_vals = 0
    if (sim_config%bc_periodicity(1)) logic_vals(1) = 1
    if (sim_config%bc_periodicity(2)) logic_vals(2) = 1

    call mpi_bcast(int_vals, size(int_vals), mpi_integer, root, comm, ierr)
    call mpi_bcast(real_vals, size(real_vals), mpi_real8, root, comm, ierr)
    call mpi_bcast(logic_vals, size(logic_vals), mpi_integer, root, comm, ierr)
    call mpi_bcast(sim_config%output_dir, len(sim_config%output_dir), mpi_character, root, comm, ierr)

    sim_config%nx        = int_vals(1)
    sim_config%ny        = int_vals(2)
    sim_config%loop_max  = int_vals(3)
    sim_config%lm_type   = int_vals(4)
    sim_config%flux_type = int_vals(5)
    sim_config%time_type = int_vals(6)
    sim_config%mpi_nums(1) = int_vals(7)
    sim_config%mpi_nums(2) = int_vals(8)
    sim_config%io_type     = int_vals(9)
    sim_config%n_start     = int_vals(10)

    sim_config%tend            = real_vals(1)
    sim_config%dtout           = real_vals(2)
    sim_config%cfl             = real_vals(3)
    sim_config%domain_x(1)     = real_vals(4)
    sim_config%domain_x(2)     = real_vals(5)
    sim_config%domain_y_min    = real_vals(6)
    sim_config%density0        = real_vals(7)
    sim_config%pressure0       = real_vals(8)
    sim_config%vx_amplitude    = real_vals(9)
    sim_config%vx_wavenumber   = real_vals(10)
    sim_config%vy_amplitude    = real_vals(11)
    sim_config%vy_wavenumber   = real_vals(12)
    sim_config%bx_amplitude    = real_vals(13)
    sim_config%bx_wavenumber   = real_vals(14)
    sim_config%by_amplitude    = real_vals(15)
    sim_config%by_wavenumber   = real_vals(16)
    sim_config%bz0             = real_vals(17)
    sim_config%ps0             = real_vals(18)

    sim_config%bc_periodicity(1) = (logic_vals(1) == 1)
    sim_config%bc_periodicity(2) = (logic_vals(2) == 1)
  end subroutine broadcast_run_config

end module run_config_mod
