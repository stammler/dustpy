subroutine cs_adiabatic(gamma, mu, T, cs, Nr)
   ! Subroutine calculates the adiabatic sound speed.
   !
   ! Parameters
   ! ----------
   ! gamma(Nr) : heat capacity ratios
   ! mu(Nr) : mean molecular masses of gas
   ! T(Nr) : Gas temperature
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! cs(Nr) : adiabatic sound speed

   use constants, only: k_B

   implicit none

   double precision, intent(in)  :: gamma(Nr)
   double precision, intent(in)  :: mu(Nr)
   double precision, intent(in)  :: T(Nr)
   double precision, intent(out) :: cs(Nr)
   integer,          intent(in)  :: Nr

   cs(:) = SQRT(gamma(:) * k_B * T(:) / mu(:))

end subroutine cs_adiabatic


subroutine enforce_floor(Sigma, SigmaFloor, R, Nr)
   ! Subroutine enforces the floor value on the surface density.
   !
   ! Parameters
   ! ----------
   ! Sigma(Nr) : Gas surface density
   ! SigmaFloor(Nr) : Gas surface density floor value
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! R(Nr) : Gas surface density with enforced floor value

   implicit none

   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(in)  :: SigmaFloor(Nr)
   double precision, intent(out) :: R(Nr)
   integer, intent(in) :: Nr

   integer :: ir

   do ir=1, Nr
      R(ir) = MAX(SigmaFloor(ir), Sigma(ir))
   end do

end subroutine enforce_floor


subroutine eta_midplane(Hp, P, r, ri, eta, Nr)
   ! Subroutine calculates the eta pressure gradient parameter in the midplane.
   ! Pressure is linearly interpolated on radial grid cell interfaces.
   !
   ! Parameters
   ! ----------
   ! Hp(Nr) : Pressure scale height
   ! P(Nr) : Pressure
   ! r(Nr) : Radial grid cell centers
   ! ri(Nr+1) : Radial grid cell interfaces
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! eta(Nr) : eta pressure gradient parameter

   implicit none

   double precision, intent(in)  :: Hp(Nr)
   double precision, intent(in)  :: P(Nr)
   double precision, intent(in)  :: r(Nr)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(out) :: eta(Nr)
   integer,          intent(in)  :: Nr

   double precision :: grad(Nr)
   double precision :: Pi(Nr+1)
   integer :: ir

   call interp1d(ri, r, P, Pi, Nr)

   do ir=1, Nr
      grad(ir) = (Pi(ir+1) - Pi(ir)) / (ri(ir+1)-ri(ir))
   end do

   eta(:)  = -0.5d0 * Hp(:)**2 / (r(:) * P(:)) * grad(:)

end subroutine eta_midplane


subroutine fi(Sigma, v, r, ri, F_i, Nr)
   ! Subroutine calculates the mass fluxes through the grid cell interfaces.
   ! Velocity v is linearly interpolated on grid cell interfaces with
   ! vi(1) = vi(2) and vi(Nr+1) = vi(Nr).
   !
   ! Parameters
   ! ----------
   ! Sigma(Nr) : Surface density
   ! v(Nr) : Radial velocity at grid cell centers
   ! r(Nr) : Radial grid cell centers
   ! ri(Nr+1) : Radial grid cell interfaces
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! F_i(Nr+1) : Flux through grid cell interfaces.

   implicit none

   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(in)  :: v(Nr)
   double precision, intent(in)  :: r(Nr)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(out) :: F_i(Nr+1)
   integer,          intent(in)  :: Nr

   double precision :: vi(Nr+1)
   double precision :: vim(Nr+1)
   double precision :: vip(Nr+1)
   integer :: ir

   call interp1d(ri, r, v, vi, Nr)
   vi(1)    = v(1)
   vi(Nr+1) = v(Nr)

   ! Positive/Negative velocity contributions
   vip(:)    = max(0.d0, vi(:))
   vim(:)    = min(vi(:), 0.d0)
   ! No inflow through inner boundary
   vip(1)    = 0.d0
   ! No inflow through outer boundary
   vim(Nr+1) = 0.d0

   do ir=2, Nr
      F_i(ir) = (Sigma(ir-1)*vip(ir) + Sigma(ir)*vim(ir))
   end do
   F_i(1)    = Sigma(1)  * vim(1)
   F_i(Nr+1) = Sigma(Nr) * vip(Nr+1)

end subroutine fi


subroutine hp(cs, gam, OmegaK, h_p, Nr)
   ! Subroutine calculates the gas scale heights.
   !
   ! Parameters
   ! ----------
   ! cs(Nr) : sound speed
   ! gam(Nr) : heat capacity ratios
   ! OmegaK(Nr) : Keplerian frequency
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! h_p(Nr) : gas scale heights

   implicit none

   double precision, intent(in)  :: cs(Nr)
   double precision, intent(in)  :: gam(Nr)
   double precision, intent(in)  :: OmegaK(Nr)
   double precision, intent(out) :: h_p(Nr)
   integer,          intent(in)  :: Nr

   h_p(:) = cs(:) / ( SQRT(gam(:)) * OmegaK(:) )

end subroutine hp


subroutine implicit_boundaries(dt, Fi, ri, Sigma, SigmaOld, ret, Nr)
   ! Subroutine calculates the implicit source terms and fluxes at the boundaries
   ! after an integration step.
   !
   ! Parameters
   ! ----------
   ! dt(Nr) : Time step that has been taken in the previous step
   ! Fi(Nr+1) : Fluxes at grid cell interfaces
   ! ri(Nr+1) : location of grid cell interfaces
   ! Sigma(Nr) : Gas surface density
   ! SigmaOld(Nr+1) : Old gas surface density from previous time step
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! ret(4) : Array with source terms and fluxes at the grid boundaries

   implicit none

   double precision, intent(in)  :: dt
   double precision, intent(in)  :: Fi(Nr+1)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(in)  :: SigmaOld(Nr)
   double precision, intent(out) :: ret(4)
   integer,          intent(in)  :: Nr

   ret(1) = ( Sigma(1)  - SigmaOld(1)  ) / (dt+1.d-100)
   ret(2) = ( Sigma(Nr) - SigmaOld(Nr) ) / (dt+1.d-100)

   ret(3) = ( 0.5d0*ret(1)*(ri(2)**2-ri(1)**2) + ri(2)*Fi(2)                          ) / ri(1)
   ret(4) = ( Fi(Nr)*ri(Nr)                    - 0.5d0*ret(2)*(ri(Nr+1)**2-ri(Nr)**2) ) / ri(Nr+1)

end subroutine implicit_boundaries


subroutine jac_abc(area, nu, r, ri, v, A, B, C, Nr)
   ! Subroutine calculates the diagonals of the gas Jacobian for advection.
   !
   ! Parameters
   ! ----------
   ! area(Nr) : radial grid annulus area
   ! nu(Nr) : gas viscosity
   ! r(Nr) : Radial grid cell centers
   ! ri(Nr+1) : Radial grid cell interfaces
   ! v(Nr) : Radial velocity (only for backreaction purposes)
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! A(Nr) : sub-diagonal, A(1) not used
   ! B(Nr) : diagonal
   ! C(Nr) : super-diagoanl, C(Nr) not used

   use constants, only: twopi

   implicit none

   double precision, intent(in)  :: area(Nr)
   double precision, intent(in)  :: nu(Nr)
   double precision, intent(in)  :: r(Nr)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(in)  :: v(Nr)
   double precision, intent(out) :: A(Nr)
   double precision, intent(out) :: B(Nr)
   double precision, intent(out) :: C(Nr)
   integer,          intent(in)  :: Nr

   integer :: ir
   double precision :: Di(Nr+1)
   double precision :: g(Nr)
   double precision :: vi(Nr+1)
   double precision :: vim(Nr+1)
   double precision :: vip(Nr+1)
   double precision :: Vinv(Nr)
   double precision :: w(Nr)

   ! Velocity at grid cell interfaces
   call interp1d(ri, r, v, vi, Nr)
   vim(:)  = min(vi(:), 0.d0)
   vip(:)  = max(0.d0, vi(:))

   ! Helper quantities
   g(:) = nu(:) / sqrt(r(:))
   Di(:) = 3.d0 * sqrt(ri(:))

   ! Initialization
   A(:)    = 0.d0
   B(:)    = 0.d0
   C(:)    = 0.d0

   ! Grid cell volumes and distances
   do ir=1, Nr
      Vinv(ir) = twopi / area(ir)
      w(ir) = r(ir+1) - r(ir)
   end do

   do ir=2, Nr-1

      A(ir) = A(ir) + vip(ir)   * r(ir-1)
      B(ir) = B(ir) - vip(ir+1) * r(ir)   + vim(ir)   * r(ir)
      C(ir) = C(ir)                       - vim(ir+1) * r(ir+1)

      A(ir) = A(ir) + Di(ir)   * g(ir-1) / w(ir-1) * r(ir-1)
      B(ir) = B(ir) - Di(ir)   * g(ir)   / w(ir-1) * r(ir)
      B(ir) = B(ir) - Di(ir+1) * g(ir)   / w(ir)   * r(ir)
      C(ir) = C(ir) + Di(ir+1) * g(ir+1) / w(ir)   * r(ir+1)

   end do

   ! Normailization
   A(:) = A(:) * Vinv(:)
   B(:) = B(:) * Vinv(:)
   C(:) = C(:) * Vinv(:)

end subroutine jac_abc


subroutine mfp_midplane(n, mfp, Nr)
   ! Subroutine calculates the mean free path in the midplane.
   !
   ! Parameters
   ! ----------
   ! n(Nr) : Gas number density
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! mfp(Nr) : mean free paths

   use constants, only: sigma_H2

   implicit none

   double precision, intent(in)  :: n(Nr)
   double precision, intent(out) :: mfp(Nr)
   integer,          intent(in)  :: Nr

   mfp(:) = 1.d0 / ( SQRT(2.d0) * sigma_H2 * n(:) )

end subroutine mfp_midplane


subroutine modified_jacobian(dt, dat, ind, indptr, A, Nd, Nc)
   ! Subroutine modifies the gas Jacobian and calculates
   ! id - dt*J
   ! for a sparse matrix in csc format.
   ! The diagonals have to be completely in the sparse matrix.
   ! Otherwise, the algorithm will fail.
   !
   ! Parameters
   ! ----------
   ! dt : Time step to be taken
   ! dat(Nd) : Data array in csc format
   ! ind(Nd) : Row indices of data values
   ! indptr(Nc+1) : Index pointers in csc format
   ! Nd : Number of elements in dat
   ! Nc : Number of columns in Jacobian
   !
   ! Returns
   ! -------
   ! A(Nd) : New values of data array

   implicit none

   double precision, intent(in)  :: dt
   double precision, intent(in)  :: dat(Nd)
   integer,          intent(in)  :: ind(Nd)
   integer,          intent(in)  :: indptr(Nc+1)
   double precision, intent(out) :: A(Nd)
   integer,          intent(in)  :: Nd
   integer,          intent(in)  :: Nc

   integer :: ic
   integer :: j

   A(:) = -dt*dat(:)

   do ic=1, Nc
      do j=indptr(ic)+1, indptr(ic+1)
         if(ind(j)+1.eq.ic) then
            A(j) = A(j) + 1.d0
         end if
      end do
   end do

end subroutine modified_jacobian


subroutine modified_rhs(dt, rhs, S, r, Nr)
   ! Subroutine modifies the right-hand side of the matrix equation
   ! by adding the source terms.
   !
   ! Parameters
   ! ----------
   ! dt : Time step to be taken
   ! rhs(Nr) : right-hand side
   ! S(Nr) : Gas source terms
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! r(Nr) : New values of right-hand side

   implicit none

   double precision, intent(in)  :: dt
   double precision, intent(in)  :: rhs(Nr)
   double precision, intent(in)  :: S(Nr)
   double precision, intent(out) :: r(Nr)
   integer,          intent(in)  :: Nr

   integer :: ir

   r(:) = 0.d0
   do ir=2, Nr-1
      r(ir) = rhs(ir) + dt*S(ir)
   end do

end subroutine modified_rhs


subroutine n_midplane(mu, rho, n, Nr)
   ! Subroutine calculates the gas number density in the midplane
   !
   ! Parameters
   ! ----------
   ! mu(Nr) : Mean molecular weight of the gas
   ! rho(Nr) : Midplane gas mass density
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! n(Nr) : Midplane gas number density

   implicit none

   double precision, intent(in)  :: mu(Nr)
   double precision, intent(in)  :: rho(Nr)
   double precision, intent(out) :: n(Nr)
   integer,          intent(in)  :: Nr

   n(:) = rho(:) / mu(:)

end subroutine n_midplane


subroutine p_midplane(cs, gam, rho, P, Nr)
   ! Subroutine calculates the midplane gas pressure.
   !
   ! Parameters
   ! ----------
   ! cs(Nr) : Sound speed
   ! gam(Nr) : Heat capacity ratios
   ! rho(Nr) : Midplane gas mass density
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! P(Nr) : Midplane gas pressure

   implicit none

   double precision, intent(in)  :: cs(Nr)
   double precision, intent(in)  :: gam(Nr)
   double precision, intent(in)  :: rho(Nr)
   double precision, intent(out) :: P(Nr)
   integer,          intent(in)  :: Nr

   P(:) = rho(:) * cs(:)**2 / gam(:)

end subroutine p_midplane


subroutine rho_midplane(Hp, Sigma, rho, Nr)
   ! Subroutine calculates the midplane gas mass density.
   !
   ! Parameters
   ! ----------
   ! Hp(Nr) : Gas scale heights
   ! Sigma(Nr) : Gas surface density
   ! Nr : Number of radial gridd cells
   !
   ! Returns
   ! rho(Nr) : Midplane gas mass density

   use constants, only: pi

   implicit none

   double precision, intent(in)  :: Hp(Nr)
   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(out) :: rho(Nr)
   integer,          intent(in)  :: Nr

   rho(:) = Sigma(:) / (SQRT(2.d0*pi) * Hp(:))

end subroutine rho_midplane


subroutine s_hyd(Fi, ri, Shyd, Nr)
   ! Subroutine calculates the hydrodynamic sources from the interface fluxes.
   !
   ! Parameters
   ! ---------
   ! Fi(Nr+1) : Mass fluxes through grid interfaces
   ! ri(Nr+1) : Grid itnerfaces
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! Shyd(Nr) : Hydrodynamic source terms

   implicit none

   double precision, intent(in)  :: Fi(Nr+1)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(out) :: Shyd(Nr)
   integer,          intent(in)  :: Nr

   integer :: ir

   do ir=1, Nr
      Shyd(ir) = 2.d0 * (Fi(ir)*ri(ir) - Fi(ir+1)*ri(ir+1)) / (ri(ir+1)**2 - ri(ir)**2)
   end do

end subroutine s_hyd


subroutine s_tot(s_ext, s_hyd, s, Nr)
   ! Subrountine calculates the total gas source terms.
   !
   ! Parameters
   ! ----------
   ! s_ext(Nr) : External source terms
   ! s_hyd(Nr) : Hydrodynamic source terms
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! s(Nr) : Total gas source terms

   implicit none

   double precision, intent(in)  :: s_ext(Nr)
   double precision, intent(in)  :: s_hyd(Nr)
   double precision, intent(out) :: s(Nr)
   integer,          intent(in)  :: Nr

   integer :: ir

   s(:) = s_hyd(:)
   do ir=2, Nr-1
      s(ir) = s(ir) + s_ext(ir)
   end do


end subroutine s_tot


subroutine t_pass(L, r, T, Nr)
   ! Function calculates the gas temperature in a passively irradiated disk.
   !
   ! Parameters
   ! ----------
   ! L : Stellar luminosity
   ! r(Nr) : Radial grid cell centers
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! T(Nr) : Gas temperature

   use constants, only: pi, sigma_sb

   implicit none

   double precision, intent(in)  :: L
   double precision, intent(in)  :: r(Nr)
   double precision, intent(out) :: T(Nr)
   integer,          intent(in)  :: Nr

   ! The prefactor is
   ! 0.5*phi/4. = 0.5*0.05/4. = 0.00625
   T(:) = ( 6.25d-3 * L / (pi * r(:)**2 * sigma_sb) )**2.5d-1

end subroutine t_pass


subroutine timestep(S, Sigma, SigmaFloor, dt, Nr)
   ! Subroutine calculates the maximum possible time step for the gas.
   !
   ! Parameters
   ! ----------
   ! S(Nr) : Gas source terms
   ! Sigma(Nr) : Gas surface density
   ! SigmaFloor(Nr) : Gas surface density floor value
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! dt : Maximum possible time step for the gas

   implicit none

   double precision, intent(in)  :: S(Nr)
   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(in)  :: SigmaFloor(Nr)
   double precision, intent(out) :: dt
   integer,          intent(in)  :: Nr

   integer :: ir

   dt = 1.d100

   do ir=2, Nr-1

      if( S(ir).LT.0.d0 .and. Sigma(ir).GT.SigmaFloor(ir)) then
         dt = MIN(dt, -Sigma(ir)/S(ir))
      end if

   end do

end subroutine timestep


subroutine v_rad(A, B, eta, OmegaK, r, vv, v, Nr)
   ! Function calculates the radial gas velocity.
   !
   ! Parameters
   ! ----------
   ! A(Nr) : A backreaction coefficient
   ! B(Nr) : B backreaction coefficient
   ! eta(Nr) : eta pressure gradient parameter
   ! OmegaK(Nr) : Keplerian frequency
   ! r(Nr) : Radial grid cell centers
   ! vv(Nr) : Viscous gas velocity
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! v(Nr) : Radial gas velocity

   implicit none

   double precision, intent(in)  :: A(Nr)
   double precision, intent(in)  :: B(Nr)
   double precision, intent(in)  :: eta(Nr)
   double precision, intent(in)  :: OmegaK(Nr)
   double precision, intent(in)  :: r(Nr)
   double precision, intent(in)  :: vv(Nr)
   double precision, intent(out) :: v(Nr)
   integer,          intent(in)  :: Nr

   double precision :: vb(Nr)

   vb(:) = 2.d0 * eta(:) * r(:) * OmegaK(:)
   v(:) = A(:)*vv(:) + B(:)*vb(:)

end subroutine v_rad


subroutine v_visc(Sigma, nu, r, ri, vvisc, Nr)
   ! Function calculates the radial viscous gas velocity.
   ! Mass flux is linearly interpolated on grid cell interfaces.
   !
   ! Parameters
   ! ----------
   ! Sigma(Nr) : Surface density
   ! nu(Nr) : Kinematic viscosity
   ! r(Nr) : Radial grid cell centers
   ! ri(Nr+1) : Radial grid cell interfaces
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! vvisc(Nr) : Radial viscous gas velocity

   implicit none

   double precision, intent(in)  :: Sigma(Nr)
   double precision, intent(in)  :: nu(Nr)
   double precision, intent(in)  :: r(Nr)
   double precision, intent(in)  :: ri(Nr+1)
   double precision, intent(out) :: vvisc(Nr)
   integer,          intent(in)  :: Nr

   double precision :: arg(Nr)
   double precision :: argi(Nr+1)
   double precision :: grad(Nr)

   integer :: ir

   arg(:) = Sigma(:) * nu(:) * sqrt(r(:))
   call interp1d(ri, r, arg, argi, Nr)

   do ir=1, Nr
      grad(ir) = (argi(ir+1) - argi(ir)) / (ri(ir+1) - ri(ir))
   end do

   vvisc(:)  = -3.d0 * grad(:) / (Sigma(:) * sqrt(r(:)))
   vvisc(1)  = (vvisc(3) - vvisc(2)) * (r(1) - r(2)) / (r(3) - r(2)) + vvisc(2)
   vvisc(Nr) = (vvisc(Nr-1) - vvisc(Nr-2)) * (r(Nr) - r(Nr-2)) / (r(Nr-1) - r(Nr-2)) + vvisc(Nr-2)

end subroutine v_visc


subroutine viscosity(alpha, cs, Hp, nu, Nr)
   ! Subroutine calculates the kinematic viscosity.
   !
   ! Parameters
   ! ----------
   ! alpha(Nr) : Alpha viscosity parameter
   ! cs(Nr) : sound speed
   ! Hp(Nr) : Gas scale height
   ! Nr : Number of radial grid cells
   !
   ! Returns
   ! -------
   ! nu : Kinematic viscosity

   implicit none

   double precision, intent(in)  :: alpha(Nr)
   double precision, intent(in)  :: cs(Nr)
   double precision, intent(in)  :: Hp(Nr)
   double precision, intent(out) :: nu(Nr)
   integer,          intent(in)  :: Nr

   nu(:) = alpha(:) * cs(:) * Hp(:)

end subroutine viscosity
