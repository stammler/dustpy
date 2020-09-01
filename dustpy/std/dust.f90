subroutine a(m, rho, sizes, Nr, Nm)
  ! Subroutine calculates the particle sizes.
  !
  ! Parameters
  ! ----------
  ! m(Nm) : Mass grid
  ! rho(Nr, Nm) : Particle bulk densities
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! sizes(Nr, Nm) : Particle sizes

  use constants, only: pi

  implicit None

  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: rho(Nr, Nm)
  double precision, intent(out) :: sizes(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i
  double precision :: onethird = 1.d0/3.d0

  do i=1, Nm
    sizes(:, i) = ( 3.d0 * m(i) / (4.d0 * pi * rho(:, i)) )**onethird
  end do

end subroutine a


subroutine d(v2, OmegaK, St, Diff, Nr, Nm)
  ! Subroutine calculates the dust diffusivity.
  !
  ! Parameters
  ! ----------
  ! v2(Nr) : turbulent gas RMS velocity
  ! OmegaK(Nr) : Keplerian frequency
  ! St(Nr, Nm) : Stokes number
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! Diff(Nr, Nm) : Dust diffusivity

  implicit none

  double precision, intent(in)  :: v2(Nr)
  double precision, intent(in)  :: OmegaK(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(out) :: Diff(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i

  do i=1, Nm
    Diff(:, i) = v2(:) / OmegaK(:) / (1.d0 + St(:, i)**2)
  end do

end subroutine d


subroutine fi_adv(Sigma, v, r, ri, Fi, Nr, Nm)
  ! Function calculates the advective mass fluxes through the grid cell interfaces.
  ! Velocity v is linearly interpolated on grid cell interfaces with
  ! vi(1, :) = vi(2, :) and vi(Nr+1, :) = vi(Nr, :).
  !
  ! Parameters
  ! ----------
  ! Sigma(Nr, Nm) : Surface density
  ! v(Nr, Nm) : Radial velocity at grid cell centers
  ! r(Nr) : Radial grid cell centers
  ! ri(Nr+1) : Radial grid cell interfaces
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! Fi(Nr+1, Nm) : Flux through grid cell interfaces.

  use constants, only: pi

  implicit none

  double precision, intent(in)  :: Sigma(Nr, Nm)
  double precision, intent(in)  :: v(Nr, Nm)
  double precision, intent(in)  :: r(Nr)
  double precision, intent(in)  :: ri(Nr+1)
  double precision, intent(out) :: Fi(Nr+1, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: vi(Nr+1)
  integer :: i
  integer :: ir

  do i=1, Nm
    call interp1d(ri, r, v(:, i), vi, Nr)
    do ir=2, Nr
      Fi(ir, i) = 2.d0*pi*ri(ir) * (Sigma(ir-1, i)*max(0.0, vi(ir)) + Sigma(ir, i)*min(vi(ir), 0.d0))
    end do
    Fi(1, i) = 2.d0*pi*ri(1) * Sigma(1, i)*min(vi(2), 0.d0)
    Fi(Nr+1, i) = 2.d0*pi*ri(Nr+1) * Sigma(Nr, i)*max(0.0, vi(Nr))
  end do

end subroutine fi_adv


subroutine fi_diff(D, eps, Hd, Hg, OmegaK, rhod, rhog, r, ri, SigmaD, SigmaG, St, v2, Fi, Nr, Nm)
  ! Subroutine calculates the diffusive mass flux through the grid cell interfaces according
  ! Klahr et al. (in prep.). The fluxes at the boundaries are set to zero.
  !
  ! Parameters
  ! ----------
  ! D(Nr, Nm) : Dust diffusivity
  ! eps(Nr) : Vertically integrated dust-to-gas ratio
  ! Hd(Nr, Nm) : Dust scale heights
  ! Hg(Nr) : Gas pressure scale height
  ! OmegaK(Nr) : Keplerian frequency
  ! rhod(Nr, Nm) : Dust midplane mass density
  ! rhog(Nr) : Gas midplane mas density
  ! r(Nr) : Radial grid cell centers
  ! ri(Nr+1) : Radial grid cell interfaces
  ! SigmaD(Nr, Nm) : Dust surface density
  ! SigmaG(Nr) : Gas surface density
  ! St(Nr, Nm) : Stokes number
  ! v2(Nr) : Gas turbulent RMS velocity
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! Fi(Nr+1, Nm) : Diffusive mass flux

  implicit none

  double precision, intent(in)  :: D(Nr, Nm)
  double precision, intent(in)  :: eps(Nr)
  double precision, intent(in)  :: Hd(Nr, Nm)
  double precision, intent(in)  :: Hg(Nr)
  double precision, intent(in)  :: OmegaK(Nr)
  double precision, intent(in)  :: rhod(Nr, Nm)
  double precision, intent(in)  :: rhog(Nr)
  double precision, intent(in)  :: r(Nr)
  double precision, intent(in)  :: ri(Nr+1)
  double precision, intent(in)  :: SigmaD(Nr, Nm)
  double precision, intent(in)  :: SigmaG(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(in)  :: v2(Nr)
  double precision, intent(out) :: Fi(Nr+1, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: A(Nr)
  double precision :: gradA(Nr+1)
  double precision :: B(Nr)
  double precision :: gradB(Nr+1)
  double precision :: C(Nr)
  double precision :: gradC(Nr+1)

  double precision :: chi(Nr)
  double precision :: epsp1(Nr)
  double precision :: f(Nr)
  double precision :: Hdi(Nr+1)
  double precision :: K(Nr)
  double precision :: Ki(Nr+1)
  double precision :: tau(Nr)
  double precision :: taui(Nr+1)

  double precision :: dum(Nr+1)
  integer :: i
  integer :: ir

  epsp1(:) = 1.d0 + eps(:)

  do ir=2, Nr+1
    dum(ir) = (r(ir) - r(ir-1)) * ri(ir)
  end do

  do i=1, Nm
    chi(:) = Hd(:, i) / Hg(:)
    call f_st(St(:, i), f, Nr)
    tau(:) = St(:, i) / OmegaK(:)
    call interp1d(ri, r, tau, taui, Nr)

    K(:) = D(:, i) * SigmaG(:) * chi(:) / epsp1
    call interp1d(ri, r, K, Ki, Nr)

    call interp1d(ri, r, Hd(:, i), Hdi, Nr)

    A(:) = (1.d0 - chi(:)**2) / Hd(:, i)**2
    B(:) = rhod(:, i) / rhog(:)
    C(:) = f(:) * v2(:) * Sigmad(:, i) / epsp1

    do ir=2, Nr
      gradA(ir) = (r(ir)*A(ir) - r(ir-1)*A(ir-1)) / dum(ir)
      gradB(ir) = (r(ir)*B(ir) - r(ir-1)*B(ir-1)) / dum(ir)
      gradC(ir) = (r(ir)*C(ir) - r(ir-1)*C(ir-1)) / dum(ir)
    end do

    Fi(:, i) = -Ki(:)*(0.5d0*Hdi**2*gradA(:) + gradB(:)) - taui(:)*gradC(:)
  end do

  Fi(1, :) = 0.d0
  Fi(Nr+1, :) = 0.d0


end subroutine fi_diff


subroutine f_st(St, f, Nr)
  ! Helper function that distinguishes two diffusion regimes by the particles
  ! Stokes numbers.
  !
  ! Parameters
  ! ----------
  ! St(Nr) : Stokes number
  ! Nr : Number of radial grid cells
  !
  ! Returns
  ! -------
  ! f(Nr) : Modification factor

  implicit none

  double precision, intent(in)  :: St(Nr)
  double precision, intent(out) :: f(Nr)
  integer :: Nr

  f(:) = 2.d0 * St(:) / (1.d0 + St(:)**3)

end subroutine f_st


subroutine h_dubrulle1995(Hp, St, delta, h, Nr, Nm)
  ! Subroutine calculates the particle scale height according Dubrulle et al. (1995).
  !
  ! Parameters
  ! ----------
  ! Hp(Nr) : Gas pressure scale heights
  ! St(Nr, Nm) : Stokes number
  ! delta(Nr) : vertical mixing parameter
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! h(Nr, Nm) : Dust scale heights

  implicit none

  double precision, intent(in)  :: Hp(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(in)  :: delta(Nr)
  double precision, intent(out) :: h(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i

  do i=1, Nm
    h(:, i) = min(Hp(:) / sqrt(1.d0 + St(:, i)/delta(:)), Hp(:))
  end do

end subroutine h_dubrulle1995


subroutine s_hyd(Fi, ri, Shyd, Nr, Nm)
  ! Subroutine calculates the hydrodynamic sources from the interface fluxes.
  !
  ! Parameters
  ! ---------
  ! Fi(Nr+1, Nm) : Mass fluxes through grid interfaces
  ! ri(Nr+1) : Grid itnerfaces
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! Shyd(Nr, Nm) : Hydrodynamic source terms

  implicit none

  double precision, intent(in)  :: Fi(Nr+1, Nm)
  double precision, intent(in)  :: ri(Nr+1)
  double precision, intent(out) :: Shyd(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i
  integer :: ir

  do i=1, Nm
    do ir=1, Nr
        Shyd(ir, i) = 2.d0 * (Fi(ir, i) - Fi(ir+1, i)) / (ri(ir+1)**2 - ri(ir)**2)
    end do
  end do

end subroutine s_hyd


subroutine st_epstein_stokes1(a, mfp, rho, Sigma, St, Nr, Nm)
  ! Subroutine calculates the Stokes number using the Epstein and the
  ! Stokes I drag regimes.
  !
  ! Parameters
  ! ----------
  ! a(Nr, Nm) : Particle sizes
  ! mfp(Nr) : Mean free path of gas
  ! rho(Nr, Nm) : Particle bulk density
  ! Sigma(Nr) : Gas surface density
  ! Nr : Number of radial grid points
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! St(Nr, Nm) : Stokes number

  use constants, only: pi

  implicit none

  double precision, intent(in)  :: a(Nr, Nm)
  double precision, intent(in)  :: mfp(Nr)
  double precision, intent(in)  :: rho(Nr, Nm)
  double precision, intent(in)  :: Sigma(Nr)
  double precision, intent(out) :: St(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i
  integer :: ir
  double precision :: twoninth = 2.d0/9.d0

  do ir=1, Nr
    do i=1, Nm
      if (a(ir, i) .LT. 2.25d0 * mfp(ir)) then
        St(ir, i) = 0.5d0 * pi * a(ir, i) * rho(ir, i) / Sigma(ir)
      else
        St(ir, i) = twoninth * pi * a(ir, i)**2 * rho(ir, i) / (mfp(ir) * Sigma(ir))
      end if
    end do
  end do

end subroutine st_epstein_stokes1


subroutine vrad(St, vdm, vrg, vr, Nr, Nm)
  ! Subroutine calculates the radial dust velocity with consists of a drift component
  ! and a gas drag component.
  !
  ! Parameters
  ! ----------
  ! St(Nr, Nm) : Stokes number
  ! vdm(Nr) : Maximum drift velocity
  ! vrg(Nr) : Radial gas velocity
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vr(Nr, Nm) : Radial dust velocities

  implicit none

  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(in)  :: vdm(Nr)
  double precision, intent(in)  :: vrg(Nr)
  double precision, intent(out) :: vr(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: i

  do i=1, Nm
    vr(:, i) = (vrg(:) + 2.d0 * vdm(:) * St(:, i)) / (St(:, i)**2 + 1.d0)
  end do

end subroutine vrad


subroutine vrel_azimuthal_drift(vdriftmax, St, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to azimuthal drift.
  !
  ! Parameters
  ! ----------
  ! vdriftmax(Nr) : Maximum drift velocity
  ! St(Nr, Nm) : Stokes number
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vrel(Nr, Nm) : Relative velocities

  implicit none

  double precision, intent(in)  :: vdriftmax(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(out) :: vrel(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: dum
  double precision :: St2p1(Nr, Nm)
  integer          :: ir, i, j

  St2p1(:, :) = St(:, :)**2 + 1.d0

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        dum = vdriftmax(ir) * abs(St2p1(ir, i) - St2p1(ir, j)) / (St2p1(ir, i) + St2p1(ir, j))
        vrel(ir, j, i) = dum
        vrel(ir, i, j) = dum
      end do
    end do
  end do

end subroutine vrel_azimuthal_drift


subroutine vrel_brownian_motion(cs, m, T, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to Brownian motion.
  ! Its maximum value is the sound speed.
  !
  ! Parameters
  ! ----------
  ! cs(Nr) : Sound speed
  ! m(Nm) : Mass grid
  ! T(Nr) : Temperature
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vrel(Nr, Nm) : Relative velocities

  use constants, only: k_B, pi

  implicit none

  double precision, intent(in)  :: cs(Nr)
  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: T(Nr)
  double precision, intent(out) :: vrel(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: ir, i, j
  double precision :: fac1
  double precision :: fac2(Nr)

  double precision :: dum

  fac1 = 8.d0 * k_B / pi

  do ir=1, Nr
    fac2(ir) = fac1 * T(ir)
  end do

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        dum = min(sqrt(fac2(ir)*(m(j)+m(i))/(m(j)*m(i))), cs(ir))
        vrel(ir, j, i) = dum
        vrel(ir, i, j) = dum
      end do
    end do
  end do

end subroutine vrel_brownian_motion


subroutine vrel_cuzzi_ormel_2007(alpha, cs, mump, OmegaK, SigmaGas, &
  & St, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to turbulent motion
  ! accourding the prescription of Cuzzi & Ormel (2007).
  !
  ! Parameters
  ! ----------
  ! alpha(Nr) : Turbulent alpha parameters
  ! cs(Nr) : Sound speed
  ! mump(Nr) : Mean molecular weight of the gas
  ! OmegaK(Nr) : Keplerian frequency
  ! SigmaGas(Nr) : Gas surface density
  ! St(Nr, Nm) : Stokes number
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vrel(Nr, Nm) : Relative velocities

  use constants, only: sigma_H2

  implicit none

  double precision, intent(in)  :: alpha(Nr)
  double precision, intent(in)  :: cs(Nr)
  double precision, intent(in)  :: mump(Nr)
  double precision, intent(in)  :: OmegaK(Nr)
  double precision, intent(in)  :: SigmaGas(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(out) :: vRel(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: eps(Nr, Nm, Nm)
  double precision :: OmKinv(Nr)
  double precision :: Re
  double precision :: ReInvSqrt(Nr)
  double precision :: StL(Nr, Nm, Nm)
  double precision :: StS(Nr, Nm, Nm)
  double precision :: tauL(Nr, Nm, Nm)
  double precision :: tauS(Nr, Nm, Nm)
  double precision :: ts(Nr)
  double precision :: vg2(Nr)
  double precision :: vn
  double precision :: vs(Nr)

  double precision :: c0, c1, c2, c3, ya, yap1inv
  double precision :: h1(Nr, Nm, Nm)
  double precision :: h2(Nr, Nm, Nm)
  double precision :: ys(Nr, Nm, Nm)

  integer :: ir, i, j
  double precision :: dum

  c0      =  1.6015125d0
  c1      = -0.63119577d0
  c2      =  0.32938936d0
  c3      = -0.29847604d0
  ya      =  1.6d0
  yap1inv =  1.d0/ (1.d0 + ya)

  do ir=1, Nr
    OmKinv(ir)    = 1.d0 / OmegaK(ir)
    Re            = 0.5d0 * alpha(ir) * SigmaGas(ir) * sigma_H2 / mump(ir)
    ReInvSqrt(ir) = sqrt(1.d0 / Re)
    vn            = sqrt(alpha(ir)) * cs(ir)
    vs(ir)        = Re**(-0.25) * vn
    ts(ir)        = OmKinv(ir) * ReInvSqrt(ir)
    vg2(ir)       = 1.5d0 * vn**2
  end do

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        StL(ir, j, i) = max(St(ir, j), St(ir, i))
        StS(ir, j, i) = min(St(ir, j), St(ir, i))
        eps(ir, j, i) = StS(ir, j, i) / StL(ir, j, i)

        tauL(ir, j, i) = StL(ir, j, i) * OmKinv(ir)
        tauS(ir, j, i) = StS(ir, j, i) * OmKinv(ir)

        ys(ir, j, i) = c0 + c1*StL(ir, j, i) + c2*StL(ir, j, i)**2 &
          & + c3*StL(ir, j, i)**3

        h1(ir, j, i) = (StL(ir, j, i) - StS(ir, j, i))                 &
          & / (StL(ir, j, i) + StS(ir, j, i))                          &
          & * (StL(ir, j, i) * yap1inv                                 &
          & - StS(ir, j, i)**2 / (StS(ir, j, i) + ya * StL(ir, j, i)))
        h2(ir, j, i) = 2.d0 * (ya * StL(ir, j, i) - ReInvSqrt(ir))     &
          & + StL(ir, j, i) * yap1inv                                     &
          & - StL(ir, j, i)**2 / (StL(ir, j, i) + ReInvSqrt(ir))       &
          & + StS(ir, j, i)**2 / (ya * StL(ir, j, i) + StS(ir, j, i))  &
          & - StS(ir, j, i)**2 / (StS(ir, j, i) + ReInvSqrt(ir))
      end do
    end do
  end do

  do i=1, Nm
    do j=1, i

      where(tauL(:, j, i) .LT. 0.2d0 * ts(:))

        vRel(:, j, i) = 1.5d0 * (vs(:) / ts(:) &
          & * (tauL(:, j, i) - tauS(:, j, i)))**2

      elsewhere(tauL(:, j, i)*ya .LT. ts(:))

        vRel(:, j, i) = vg2(:) * (StL(:, j, i) - StS(:, j, i)) &
          & / (StL(:, j, i) + StS(:, j, i)) * (StL(:, j, i)**2 &
          & / (StL(:, j, i) + ReInvSqrt(:))                    &
          & - StS(:, j, i)**2 / (StS(:, j, i) + ReInvSqrt(:)))

      elsewhere(tauL(:, j, i) .LT. 5.d0 * ts(:))

        vRel(:, j, i) = vg2(:) * (h1(:, j, i) + h2(:, j, i))

      elsewhere(tauL(:, j, i) .LT. 0.2d0 * OmKinv(:))

        vRel(:, j, i) = vg2(:) * StL(:, j, i)                             &
          & * (2.d0*ya - 1.d0 - eps(:, j, i) + 2.d0/(1.d0 + eps(:, j, i)) &
          & * (yap1inv + eps(:, j, i)**3/(ya+eps(:, j, i))))

      elsewhere(tauL(:, j, i) .LT. OmKinv(:))

        vRel(:, j, i) =  vg2(:) * StL(:, j, i)                  &
          & * (2.d0*ys(:, j, i) - 1.d0 - eps(:, j, i)           &
          & + 2.d0/(1.d0+eps(:, j, i))*(1.d0/(1.d0+ys(:, j, i)) &
          & + eps(:, j, i)**3/(ys(:, j, i)+eps(:, j, i))))

      elsewhere(tauL(:, j, i) .GE. OmKinv(:))

        vRel(:, j, i) = vg2(:) &
          & * (2.d0 + StL(:, j, i) + StS(:, j, i)) &
          & / (1.d0 + StL(:, j, i) + StS(:, j, i) + StL(:, j, i)*StS(:, j, i))

      end where

    end do
  end do

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        dum = sqrt(vRel(ir, j, i))
        vRel(ir, j, i) = dum
        vRel(ir, i, j) = dum
      end do
    end do
  end do

end subroutine vrel_cuzzi_ormel_2007


subroutine vrel_radial_drift(vrad, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to radial drift.
  !
  ! Parameters
  ! ----------
  ! vrad(Nr, Nm) : Radial velocities
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vrel(Nr, Nm) : Relative velocities

  implicit none

  double precision, intent(in)  :: vrad(Nr, Nm)
  double precision, intent(out) :: vrel(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm
  
  integer :: ir, i, j
  double precision :: dum

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        dum = abs(vrad(ir, j) - vrad(ir, j))
        vrel(ir, i, j) = abs(vrad(ir, j) - vrad(ir, j))
        vrel(ir, j, i) = abs(vrad(ir, j) - vrad(ir, j))
      end do
    end do
  end do

end subroutine vrel_radial_drift


subroutine vrel_vertical_settling(h, OmK, St, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to vertical settling.
  !
  ! Parameters
  ! ----------
  ! h(Nr, Nm) : Particle scale heights
  ! Omk(Nr) : Keplerian frequency
  ! St(Nr, Nm) : Stokes number
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! vrel(Nr, Nm) : Relative velocities

  implicit none

  double precision, intent(in)  :: h(Nr, Nm)
  double precision, intent(in)  :: OmK(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(out) :: vrel(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: ir, i, j

  double precision :: dum
  double precision :: OmMinSt(Nr, Nm)

  do i=1, Nm
    do ir=1, Nr
      OmMinSt(ir, i) = OmK(ir) * h(ir, i) * min(St(ir, i), 0.5d0)
    end do
  end do

  do i=1, Nm
    do j=1, i
      do ir=1, Nr
        dum = abs(OmMinSt(ir, j) - OmMinSt(ir, i))
        vrel(ir, j, i) = dum
        vrel(ir, i, j) = dum
      end do
    end do
  end do

end subroutine vrel_vertical_settling