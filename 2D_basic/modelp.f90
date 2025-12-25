subroutine modelp(U,V,x,y,dx,ix,jx)
  use parallel
  use run_config_mod, only: sim_config
  implicit none
  include 'param.h'
  real(8), intent(out) :: U(ix,jx,var1)
  real(8), intent(out) :: V(ix,jx,var2)
  real(8), intent(out) :: x(ix), y(jx), dx
  integer, intent(in)  :: ix, jx
! ---------------------------------------------------
  real(8) :: domain_x(2)
  real(8) :: domain_y_min
  real(8) :: density0, pressure0
  real(8) :: vx_amp, vx_k
  real(8) :: vy_amp, vy_k
  real(8) :: bx_amp, bx_k
  real(8) :: by_amp, by_k
  real(8) :: bz0, ps0
  integer :: i, j
  integer :: iix, jjx
  real(8) :: B2, v2, f1
  real(8) :: tmpx(cart2d%sizes(1)*(ix-2) + 2)
  real(8) :: tmpy(cart2d%sizes(2)*(jx-2) + 2)
! ---------------------------------------------------

  domain_x     = sim_config%domain_x
  domain_y_min = sim_config%domain_y_min
  density0     = sim_config%density0
  pressure0    = sim_config%pressure0
  vx_amp       = sim_config%vx_amplitude
  vx_k         = sim_config%vx_wavenumber
  vy_amp       = sim_config%vy_amplitude
  vy_k         = sim_config%vy_wavenumber
  bx_amp       = sim_config%bx_amplitude
  bx_k         = sim_config%bx_wavenumber
  by_amp       = sim_config%by_amplitude
  by_k         = sim_config%by_wavenumber
  bz0          = sim_config%bz0
  ps0          = sim_config%ps0

! grid in the X direction
  iix = cart2d%sizes(1)*(ix-2) + 2
  dx = ( domain_x(2) - domain_x(1) ) / dble( iix-2 )
  tmpx(1)   = domain_x(1) - dx/2
! tmpx(iix) = domain_x(2) + dx/2
  do i=2,iix
     tmpx(i) = tmpx(i-1) + dx
  enddo
  do i=1,ix
     x(i) = tmpx(cart2d%coords(1)*(ix-2) + i)
  enddo

! grid in the Y direction
  jjx = cart2d%sizes(2)*(jx-2) + 2
  tmpy(1) = domain_y_min - dx/2
  do j=2,jjx
     tmpy(j) = tmpy(j-1) + dx
  enddo
  do j=1,jx
     y(j) = tmpy(cart2d%coords(2)*(jx-2) + j)
  enddo
! ---------------------------------------------------

  do j=1,jx
  do i=1,ix

     U(i,j,ro) = density0

     V(i,j,pr) = pressure0
     V(i,j,vx) = vx_amp * sin(vx_k * y(j))
     V(i,j,vy) = vy_amp * sin(vy_k * x(i))
     V(i,j,vz) = 0.d0

     U(i,j,bx) = bx_amp * sin(bx_k * y(j))
     U(i,j,by) = by_amp * sin(by_k * x(i))
     U(i,j,bz) = bz0
     U(i,j,ps) = ps0

     f1 = 1.d0 / ( gamma - 1 )
     v2 = dot_product( V(i,j,vx:vz), V(i,j,vx:vz) )
     B2 = dot_product( U(i,j,bx:bz), U(i,j,bx:bz) )
     U(i,j,mx:mz) = U(i,j,ro) * V(i,j,vx:vz)
     U(i,j,en)    = 0.5d0 * ( U(i,j,ro)*v2 + B2 ) + f1*V(i,j,pr)

  enddo
  enddo

  return
end subroutine modelp
