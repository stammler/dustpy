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

    eta(:)  = -0.5d0 * Hp(:)**2 / (r(:) * P(:)) * grad(:)

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