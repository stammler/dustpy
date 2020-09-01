subroutine eta_midplane(Hp, P, r, ri, eta, Nr)
    ! Function calculates the eta pressure gradient parameter in the midplane.
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

    eta(:) = -0.5d0 * Hp(:) / (r(:) * P(:)) * grad(:)

end subroutine eta_midplane


subroutine fi(Sigma, v, r, ri, F_i, Nr)
    ! Function calculates the mass fluxes through the grid cell interfaces.
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

    use constants, only: pi

    implicit none

    double precision, intent(in)  :: Sigma(Nr)
    double precision, intent(in)  :: v(Nr)
    double precision, intent(in)  :: r(Nr)
    double precision, intent(in)  :: ri(Nr+1)
    double precision, intent(out) :: F_i(Nr+1)
    integer,          intent(in)  :: Nr

    double precision :: vi(Nr+1)
    integer :: ir

    call interp1d(ri, r, v, vi, Nr)

    do ir=2, Nr
        F_i(ir) = 2.d0*pi*ri(ir) * (Sigma(ir-1)*max(0.0, vi(ir)) + Sigma(ir)*min(vi(ir), 0.d0))
    end do
    F_i(1) = 2.d0*pi*ri(1)*Sigma(1)*min(vi(2), 0.d0)
    F_i(Nr+1) = 2.d0*pi*ri(Nr+1)*Sigma(Nr)*max(0.d0, vi(Nr))

end subroutine fi

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
        Shyd(ir) = 2.d0 * (Fi(ir) - Fi(ir+1)) / (ri(ir+1)**2 - ri(ir)**2)
    end do

end subroutine s_hyd

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

    arg = Sigma * nu * sqrt(r)
    call interp1d(ri, r, arg, argi, Nr)
    
    do ir=1, Nr
        grad(ir) = (argi(ir+1) - argi(ir)) / (ri(ir+1) - ri(ir))
    end do

    vvisc(:) = -3 * grad(:) / (Sigma(:) * sqrt(r(:)))

end subroutine v_visc