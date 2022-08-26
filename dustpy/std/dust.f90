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


subroutine check_mass_conservation_erosion(AFrag, epsFrag, klf, krm, m, phiFrag, Nm, errmax, imax, jmax)
  ! Function calculates the maximum relative mass error in erosive collisions.
  !
  ! Parameters
  ! ----------
  ! AFrag(Nm, Nm) : Normalization factor of fragment distribution
  ! epsFrag(Nm, Nm) : Fraction of remnant mass that get distributed into the
  !                   lower mass bin
  ! klf(Nm, Nm) : Index of largest fragment
  ! krm(Nm, Nm) : Smaller index of remnant mass
  ! m(Nm) : Mass grid
  ! phiFrag(Nm, Nm) : The distribution of fragments
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! errmax : maximum relative mass error
  ! imax, jmax : indices of maximum mass error

  implicit none

  double precision, intent(in)  :: AFrag(Nm, Nm)
  double precision, intent(in)  :: epsFrag(Nm, Nm)
  integer,          intent(in)  :: klf(Nm, Nm)
  integer,          intent(in)  :: krm(Nm, Nm)
  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: phiFrag(Nm, Nm)
  integer,          intent(in)  :: Nm
  double precision, intent(out) :: errmax
  integer,          intent(out) :: imax
  integer,          intent(out) :: jmax

  integer :: i
  integer :: j
  integer :: k
  double precision :: merr
  double precision :: mfrag
  double precision :: mrm
  double precision :: mtot

  ! Initialization
  errmax = 0.d0

  do i=1, Nm
    do j=1, i

      if(krm(j, i)+1 .EQ. 0) cycle

      mtot = m(i) + m(j)

      mfrag = 0.d0
      do k=1, Nm
        mfrag = mfrag + AFrag(j, i) * phiFrag(klf(j, i)+1, k)
      end do
      mrm = epsFrag(j, i)*m(krm(j, i)+1) + (1.d0-epsFrag(j, i))*m(krm(j, i)+2)
      
      merr = ABS(1.d0 - mfrag/mtot - mrm/mtot)
      if(merr .GT. errmax) then
        errmax = merr
        imax = i - 1  ! Conversion to Python indexing
        jmax = j - 1
      end if

    end do
  end do

end subroutine check_mass_conservation_erosion


subroutine check_mass_conservation_full_fragmentation(AFrag, klf, m, phiFrag, Nm, errmax, imax, jmax)
  ! Function calculates the maximum relative mass error in full fragmenting collisions.
  !
  ! Parameters
  ! ----------
  ! AFrag(Nm, Nm) : Normalization factor of fragment distribution
  ! klf(Nm, Nm) : Index of largest fragment
  ! phiFrag(Nm, Nm) : The distribution of fragments
  ! m(Nm) : Mass grid
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! errmax : maximum relative mass error
  ! imax, jmax : indices of maximum mass error

  implicit none

  double precision, intent(in)  :: AFrag(Nm, Nm)
  integer,          intent(in)  :: klf(Nm, Nm)
  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: phiFrag(Nm, Nm)
  integer,          intent(in)  :: Nm
  double precision, intent(out) :: errmax
  integer,          intent(out) :: imax
  integer,          intent(out) :: jmax

  integer :: i
  integer :: j
  integer :: k
  double precision :: merr
  double precision :: mfrag
  double precision :: mtot

  ! Initialization
  errmax = 0.d0

  do i=1, Nm
    do j=1, i

      if(klf(j, i)+1 .NE. i) cycle

      mfrag = 0.d0
      mtot = m(i) + m(j)

      do k=1, Nm
        mfrag = mfrag + AFrag(j, i) * phiFrag(i, k)
      end do

      merr = ABS(1.d0 - mfrag/mtot)
      if(merr .GT. errmax) then
        errmax = merr
        imax = i - 1  ! Conversion to Python indexing
        jmax = j - 1
      end if

    end do
  end do

end subroutine check_mass_conservation_full_fragmentation


subroutine check_mass_conservation_sticking(cstick, cstick_ind, m, Nm, errmax, imax, jmax)
  ! Function calculates the maximum relative mass error in sticking collisions.
  !
  ! Parameters
  ! ----------
  ! cstick(4, Nm, Nm) : Non-zero elements of sticking matrix
  ! cstick_ind(4, Nm, Nm) : location of non-zero elements in third dimension
  !                         in Python indexing starting from 0
  ! m(Nm) : Mass grid
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! errmax : maximum relative mass error
  ! imax, jmax : indices of maximum mass error

  implicit none

  double precision, intent(in)  :: cstick(4, Nm, Nm)
  integer,          intent(in)  :: cstick_ind(4, Nm, Nm)
  double precision, intent(in)  :: m(Nm)
  integer,          intent(in)  :: Nm
  double precision, intent(out) :: errmax
  integer,          intent(out) :: imax
  integer,          intent(out) :: jmax

  integer :: i
  integer :: j
  integer :: k
  integer :: l
  double precision :: merr
  double precision :: msum
  double precision :: mtot

  ! Initializing
  errmax = 0.d0

  do i=1, Nm
    do j=1, i

      mtot = m(i) + m(j)
      msum = 0.d0

      do l=1, 4
        k = cstick_ind(l, j, i) + 1 ! Conversion to Fortran indexing
        if(k .EQ. 0) cycle
        msum = msum + cstick(l, j, i)*m(k)
      end do

      merr = ABS(msum/mtot)
      if(merr .GT. errmax) then
        errmax = merr
        imax = i - 1  ! Conversion to Python indexing
        jmax = j - 1
      end if

    end do
  end do

end subroutine check_mass_conservation_sticking


subroutine coagulation_parameters(cratRatio, fExcav, fragSlope, m, cstick, cstick_ind, AFrag, epsFrag, klf, krm, phiFrag, Nr, Nm)
  ! Subroutine calculates the coagulation parameters needed to calculate the
  ! coagulation sources. The sticking matrix is calculated with the method
  ! described in appendices A.1. and A.2. of Brauer et al. (2008).
  ! The fragmentation follows the prescription of Rafikov et al. (2020) with
  ! the difference that it is strictly mass conserving.
  !
  ! Parameters
  ! ----------
  ! cratRatio : Mass ratio below which particles fully fragment
  ! fExcav : Excavated erosive mass in units of smaller particle
  ! fragSlope : Power of fragment distribution
  ! m(Nm) : Mass grid
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! cstick(4, Nm, Nm) : Non-zero elements of sticking matrix
  ! cstick_ind(4, Nm, Nm) : location of non-zero elements in third dimension
  !                         in Python indexing starting from 0
  ! AFrag(Nm, Nm) : Normalization factor of fragment distribution
  ! epsFrag(Nm, Nm) : Fraction of remnant mass that get distributed into the
  !                   lower mass bin
  ! klf(Nm, Nm) : Index of largest fragment
  ! krm(Nm, Nm) : Smaller index of remnant mass
  ! phiFrag(Nm, Nm) : The distribution of fragments

  implicit none

  double precision, intent(in)  :: cratRatio
  double precision, intent(in)  :: fExcav
  double precision, intent(in)  :: fragSlope
  double precision, intent(in)  :: m(Nm)
  double precision, intent(out) :: AFrag(Nm, Nm)
  double precision, intent(out) :: epsFrag(Nm, Nm)
  integer,          intent(out) :: klf(Nm, Nm)
  integer,          intent(out) :: krm(Nm, Nm)
  double precision, intent(out) :: cstick(4, Nm, Nm)
  integer,          intent(out) :: cstick_ind(4, Nm, Nm)
  double precision, intent(out) :: phiFrag(Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: a
  double precision :: cpod(Nm, Nm, Nm)
  double precision :: cpodmod(Nm, Nm, Nm)
  double precision :: dum(Nm, Nm, Nm)
  double precision :: D(Nm, Nm)
  double precision :: E(Nm, Nm)
  double precision :: eps
  double precision :: kdelta
  double precision :: mi_m_mim1(Nm)
  double precision :: mip1_m_mi(Nm)
  double precision :: mip1_m_mim1(Nm)
  double precision :: mrm
  double precision :: mtot(Nm, Nm)
  double precision :: theta
  integer :: ce
  integer :: i
  integer :: j
  integer :: jmax
  integer :: jmin
  integer :: k
  integer :: inc
  integer :: p

  ! Initialization
  AFrag(:, :) = 0.d0
  cpod(:, :, :) = 0.d0
  cpodmod(:, :, :) = 0.d0
  cstick(:, :, :) = 0.d0
  cstick_ind(:, :, :) = -1
  D(:, :) = -1.d0
  epsFrag(:, :) = 0.d0
  klf(:, :) = 0
  krm(:, :) = 0
  mtot(:, :) = 0.d0
  phiFrag(:, :) = 0.d0

  ! Grid constants
  ! These grid constants are only valid for a regular logarithmic grid.
  ! For other grids this coagulation method will produce wrong results.

  ! a is defined as
  ! m(i+1) = m(i) * 10**a
  a = log10(m(1)/m(Nm) ) / (1.d0 - Nm)
  
  ! ce from Brauer et al. (2008) equation (A.6):
  ! m(k-1) + m(i) < m(k) for any i with i ≤ k - ce
  ce = aint( -1.d0/a * log10(1.d0 - 10.d0**(-a) ) ) + 1

  ! Since the grid is regular logarithmic particles with masses m(i) and m(j)
  ! fully fragment if j >= i-p.
  p = floor(log10(cratRatio)/a)

  ! COAGULATION

  ! Helper arrays to avoid out-of-bounds errors
  ! mi_m_mim1(i) = m(i) - m(i-1)
  mi_m_mim1(:)   = m(:) * ( 1.d0    - 10.d0**(-a))
  ! mip1_m_mi(i) = m(i+1) - m(i)
  mip1_m_mi(:)   = m(:) * (10.d0**a -  1.d0      )
  ! mip1_m_mim1(i) = m(i+1) - m(i-1)
  mip1_m_mim1(:) = m(:) * (10.d0**a - 10.d0**(-a))

  do i=1, Nm
    do j=1, Nm

      ! Total mass of collision partners
      mtot(j, i) = m(j) + m(i)

      ! Podolak coefficients according Brauer et al. (2008), A.1.
      ! Meaning of cpod(k, j, i):
      ! if mmass (i) collides with mass m(j), the fraction of eps(j, i)
      ! will be filled into mass m(k) and the fraction of 1 - eps(j, i)
      ! into mass m(k+1), since the total mass of both collision partners
      ! will lie between the two masses m(k) and m(k+1).
      k = minloc(m, 1, mtot(j, i) < m) - 1
      if(k .gt. 0) then
        eps             = ( m(k+1) - mtot(j, i) ) / mip1_m_mi(k)
        cpod(k, j, i)   = eps
        cpod(k+1, j, i) = 1.d0 - eps
      end if

      ! Modified Podolak coefficients according Brauer et al. (2008), A.2.
      ! If the two colliding masses differ by more than fifteen orders of 
      ! magnitude double precision is not precise enough, i.e.,
      ! m(i) + m(j) = m(i) from a numerical point of view. This violates
      ! mass conservation. To prevent this, the following procedure is applied.
      ! Please read Brauer et al. (2008), appendix A.2. for details.

      ! D matrix
      if(j .le. i+1-ce) then
        D(j, i) = - m(j) /  mip1_m_mi(i)
      end if

      ! E matrix
      if(j .le. i-ce) then
        E(j, i) = m(j) / mi_m_mim1(i)
      else
        E(j, i) = ( 1.d0 - ( m(j) - mi_m_mim1(i) ) / mip1_m_mi(i) )  *  theta( mip1_m_mim1(i) - m(j)  )
      end if
    end do
  end do

  ! Building modified Podolak coefficients
  ! The outer sum does not need to go all the way up to Nm,
  ! because sticking collisions involving m(Nm) will be outside the mass grid.
  ! The inner sum will only go to a maximum integer number jmax, such that
  ! the resulting mass is within the mass grid.
  do i=1, Nm-1
    jmax = MIN(Nm-1, INT(LOG10(10.d0**(a*Nm)-10.d0**(a*i))/a))
    do j=1, jmax
      ! Only if resulting mass inside of grid
      do k=1, Nm
        cpodmod(k, j, i) = 0.5d0 * kdelta(j, i) * cpod(k, j, i) &
                           & + cpod(k, j, i) * theta( (k-j)*1.d0 - 1.5d0 ) * theta( (j-i)*1.d0 - 0.5d0 ) &
                           & + kdelta(j, k) * d(i, j) &
                           & + kdelta(j, k-1) * e(i, k) * theta( (k-i)*1.d0 - 1.5d0 )
      end do
    end do
  end do

  ! cpodmod so far is not symmetric:
  ! cpodmod(k, j, i) + cpodmod(k, i, j) is the fraction that gets added to
  ! mass m(k), if masses m(i) and m(j) collide.
  ! Here we make it symmetric and the discard one half.
  ! So cpodmod(k, j, i) contains the full value if j <= i and is zero,
  ! if j > i.
  dum(:, :, :) = cpodmod(:, :, :)
  do i=1, Nm
    do j=1, i
      cpodmod(:, i, j) = 0.d0
      cpodmod(:, j, i) = dum(:, j, i) + dum(:, i, j)
    end do
  end do

  ! copodmod(:, j, i) has at most four non-zero elements for any combination of (j, i).
  ! We only store the non-zero elemtens.
  do i=1, Nm
    do j=1, i
      inc = 1
      do k=1, Nm
        if (cpodmod(k, j, i) .NE. 0.d0) then
          cstick_ind(inc, j, i) = k - 1 ! Minus one to convert to Python indexing
          cstick(inc, j, i) = cpodmod(k, j, i)
          inc = inc + 1
        end if
      end do
    end do
  end do

  ! FRAGMENTATION

  ! The fragment distribution
  do i=1, Nm
    phiFrag(i, 1:i) = m(1:i)**(2.d0 + fragSlope)
    phiFrag(i, 1:i) = phiFrag(i, 1:i) / SUM( phiFrag(i, 1:i) )
  end do

  do i=1, Nm
    ! Erosion
    do j=1, i-p-1
      ! FRAGMENT DISTRTIBUTION
      ! The largest fragment has the mass of the smaller
      ! collision partner

      ! Mass bin of largest fragment
      klf(j, i) = j

      ! Remnant particle mass is mass of larger collision partner minues
      ! excavated mass
      mrm = m(i) - fExcav*m(j)

      ! Normalization factor of fragment distribution (klf+1 for Fortran indexing)
      AFrag(j, i) = (1.d0+fExcav)*m(j)
      !             |________________|
      !                      |
      !                      Mass of fragments

      ! REMNANT PARTICLE
      ! The remnant particle mass is distributed between to adjacent
      ! mass bins

      ! Lower mass bin in which the remnant mass is distributed
      krm(j, i) = minloc(m, 1, mrm < m) - 1
      ! Fraction of remnant particle that gets distributed in the lower
      ! mass bin krm. The fraction that gets distributed into the larger
      ! mass bin krm+1 i given by 1-epsFrag
      epsFrag(j, i) = ( m( krm(j, i) + 1 ) - mrm ) / mip1_m_mi( krm(j, i) )
      ! ATTENTION: If the mass difference is large enough the remnant
      ! mass can be equal to the largest particle. In that case we deal
      ! with it separately.
      if(mrm .eq. m(i)) then
        krm(j, i) = i - 1
        epsFrag(j, i) = fExcav*m(j) / mi_m_mim1(i)
      end if

    end do

    ! Full fragmentation
    jmin = MAX(1, i-p)
    do j=jmin, i

      ! The largest fragment has the mass of the larger collison partner
      klf(j, i) = i

      ! Normalization factor of fragment distribution
      AFrag(j, i) = (m(i)+m(j))

    end do
  end do

  ! Converting to Python indexing
  klf(:, :) = klf(:, :) - 1
  krm(:, :) = krm(:, :) - 1
  

end subroutine coagulation_parameters


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
    Diff(:, i) = v2(:) / ( OmegaK(:) * (1.d0 + St(:, i)**2) )
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
      Fi(ir, i) = Sigma(ir-1, i)*max(0.0, vi(ir)) + Sigma(ir, i)*min(vi(ir), 0.d0)
    end do
    Fi(1, i) = Sigma(1, i)*min(vi(2), 0.d0)
    Fi(Nr+1, i) = Sigma(Nr, i)*max(0.0, vi(Nr))
  end do

end subroutine fi_adv


subroutine fi_diff(D, SigmaD, SigmaG, St, u, r, ri, Fi, Nr, Nm)
  ! Subroutine calculates the diffusive dust fluxes at the grid cell interfaces.
  ! The flux at the boundaries is assumed to be constant.
  !
  ! Parameters
  ! ----------
  ! D(Nr, Nm) : Dust diffusivity
  ! SigmaD(Nr, Nm) : Dust surface densities
  ! SigmaG(Nr) : Gas surface density
  ! St(Nr, Nm) : Stokes number
  ! u(Nr) : Gas turbulent RMS velocity
  ! r(Nr) : Radial grid cell centers
  ! ri(Nr+1) : Radial grid cell interfaces
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! Fi(Nr+1, Nm) : Diffusive fluxes at grid cell interfaces

  implicit none

  double precision, intent(in)  :: D(Nr, Nm)
  double precision, intent(in)  :: SigmaD(Nr, Nm)
  double precision, intent(in)  :: SigmaG(Nr)
  double precision, intent(in)  :: St(Nr, Nm)
  double precision, intent(in)  :: u(Nr)
  double precision, intent(in)  :: r(Nr)
  double precision, intent(in)  :: ri(Nr+1)
  double precision, intent(out) :: Fi(Nr+1, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: Di(Nr+1, Nm)
  double precision :: eps(Nr, Nm)
  double precision :: gradepsi(Nr+1, Nm)
  double precision :: lambda
  double precision :: P
  double precision :: SigDi(Nr+1, Nm)
  double precision :: SigGi(Nr+1)
  double precision :: Sti(Nr+1, Nm)
  double precision :: ui(Nr+1)
  double precision :: w
  integer :: ir
  integer :: i

  Fi(:, :) = 0.d0

  call interp1d(ri, r, SigmaG, SigGi, Nr)
  call interp1d(ri, r, u, ui, Nr)

  do i=1, Nm
    call interp1d(ri(:), r(:), D(:, i), Di(:, i), Nr)
    call interp1d(ri(:), r(:), SigmaD(:, i), SigDi(:, i), Nr)
    eps(:, i) = SigmaD(:, i) / SigmaG(:)
    call interp1d(ri(:), r(:), St(:, i), Sti(:, i), Nr)
  end do

  do ir=2, Nr
    gradepsi(ir, :) = ( eps(ir, :) - eps(ir-1, :) ) / ( r(ir) - r(ir-1) )
  end do

  do i=1, Nm
    do ir=2, Nr
      w = ui(ir) * SigDi(ir, i) / (1.d0 + Sti(ir, i)**2)
      Fi(ir, i) = -Di(ir, i) * SigGi(ir) * gradepsi(ir, i)
      P = abs( Fi(ir, i) / w )
      lambda = ( 1.d0 + P ) / ( 1.d0 + P + P**2 )
      if(lambda .GT. HUGE(lambda)) then
        Fi(ir, i) = w
      else
        Fi(ir, i) = lambda * Fi(ir, i)
      end if
    end do
  end do

  Fi(   1, :) = Fi( 2, :)
  Fi(Nr+1, :) = Fi(Nr, :)

end subroutine fi_diff


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


subroutine jacobian_coagulation_generator(A, cStick, eps, iLF, iRM, iStick, m, phi, Rf, Rs, Sigma, SigmaFloor, &
  & dat, row, col, Nr, Nm)
  ! Subroutine calculates the coagulation Jacobian at every radial grid cell except for the boundaries.
  !
  ! Parameters
  ! ----------
  ! A(Nm, Nm) : Normalization factor of fragment distribution
  ! cStick(4, Nm, Nm) : Non-zero elements of sticking matrix
  ! eps(Nm, Nm) : Remnant mass distribution
  ! iLF(Nm, Nm) : Index of largest fragment
  ! iRM(Nm, Nm) : Smaller index of remnant mass
  ! iStick(4, Nm, Nm) : Non-zero indices of sticking matrix
  ! m(Nm) : mass grid
  ! phi(Nm, Nm) : Fragment distribution
  ! Rf(Nr, Nm, Nm) : Fragmentation collision rates
  ! Rs(Nr, Nm, Nm) : Sticking collision rates
  ! Sigma(Nr, Nm) : Dust surface density
  ! SigmaFloor(Nr, Nm) : Floor value of surface density
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! jac(Nr, Nm, Nm) : Coagulation Jacobian
  
  implicit none

  double precision, intent(in)  :: A(Nm, Nm)
  double precision, intent(in)  :: cStick(4, Nm, Nm)
  double precision, intent(in)  :: eps(Nm, Nm)
  integer,          intent(in)  :: iLF(Nm, Nm)
  integer,          intent(in)  :: iRM(Nm, Nm)
  integer,          intent(in)  :: iStick(4, Nm, Nm)
  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: phi(Nm, Nm)
  double precision, intent(in)  :: Rf(Nr, Nm, Nm)
  double precision, intent(in)  :: Rs(Nr, Nm, Nm)
  double precision, intent(in)  :: Sigma(Nr, Nm)
  double precision, intent(in)  :: SigmaFloor(Nr, Nm)
  double precision, intent(out) :: dat((Nr-2)*Nm*Nm)
  integer,          intent(out) :: row((Nr-2)*Nm*Nm)
  integer,          intent(out) :: col((Nr-2)*Nm*Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: agrid
  double precision :: jac(Nr, Nm, Nm)
  double precision :: N(Nr, Nm)
  double precision :: D(Nm, Nm)
  double precision :: ratef
  double precision :: rates

  integer :: imax
  integer :: ir
  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: p
  integer :: q
  integer :: s
  integer :: start

  ! Determination of the mass bin range that contributes to
  ! full fragmentation
  p = Nm - minloc( m(:), 1, iRM(:, Nm) .EQ. -1 )

  ! Grid constant
  agrid = LOG10( m(1)/m(Nm) ) / (1.d0 - Nm)
  ! Number of non-zero off-diagonal
  q = INT( CEILING( LOG10(2.d0)/agrid ) )

  ! Initialization
  jac(:, :, :) = 0.d0
  dat(:) = 0.d0
  row(:) = 0
  col(:) = 0

  ! Conversion to number density
  do ir=1, Nr
    N(ir, :) = Sigma(ir, :) / m(:)
  end do

  ! Normalization factor for sticking
  do i=1, Nm
    do j=1, Nm
      D(j, i) = m(j)/m(i)
    end do
  end do

  do ir=2, Nr-1

    ! Largest mass bin that is above floor value
    imax = min( maxloc(m(:), 1, Sigma(ir, :) .GT. SigmaFloor(ir, :)), Nm )

    do i=1, imax
      do j=1, i

        ! STICKING

        ! Only if particles actually stick
        if(Rs(ir, j, i) .GT. 0.d0) then

          rates = N(ir, j) * Rs(ir, j, i)

          do l=1, 4
            k = iStick(l, j, i) + 1
            if(k .EQ. 0) cycle
            jac(ir, k, i) = jac(ir, k, i) + D(k, i) * cStick(l, j, i) * rates
          end do

        end if

        ! FRAGMENTATION

        ! Only if particles actually fragment
        if(Rf(ir, j ,i) .GT. 0.d0) then

          ! Fragmentation rate factor
          ratef = N(ir, j) * Rf(ir, j, i)

          ! Fragments distribution
          jac(ir, :, i) = jac(ir, :, i) + A(j, i) * phi(ilf(j, i)+1, :) / m(i) * ratef

          ! Negative terms and remnant masses
          ! Erosion
          if(j .LE. i-p-1) then

            k = irm(j, i) + 1
            ! Remnant mass falls between i-1 and i
            ! This distinction is beneficial for mass conservation
            if(k .EQ. i-1) then
              jac(ir, i-1, i) = jac(ir, i-1, i) + D(i-1, i) * eps(j, i)          * ratef
              jac(ir, i,   i) = jac(ir, i,   i) -             eps(j, i)          * ratef
              jac(ir, j,   i) = jac(ir, j,   i) - D(j,   i)                      * ratef
            else
              jac(ir, k,   i) = jac(ir, k,   i) + D(k,   i) * eps(j, i)          * ratef
              jac(ir, k+1, i) = jac(ir, k+1, i) + D(k+1, i) * (1.d0 - eps(j, i)) * ratef
              jac(ir, i,   i) = jac(ir, i,   i) -                                  ratef
              jac(ir, j,   i) = jac(ir, j,   i) - D(j,   i)                      * ratef
            end if

          ! Full fragmentation
          else
            jac(ir, i,   i) = jac(ir, i,   i) -             ratef
            jac(ir, j,   i) = jac(ir, j,   i) - D(j,   i) * ratef
          end if

        end if

      end do
    end do
  end do

  ! Filling the data array
  k = 1
  do ir=2, Nr-1
    start = (ir-1)*Nm - 1
    do i=1, Nm
      s = MAX(i-q, 1)
      do j=s, Nm
        dat(k) = jac(ir, i, j)
        row(k) = start + i
        col(k) = start + j
        k = k + 1
      end do
    end do
  end do

end subroutine jacobian_coagulation_generator


subroutine jacobian_hydrodynamic_generator(area, D, r, ri, SigmaGas, v, A, B, C, Nr, Nm)
  ! Subroutine calculates the diagonals of the hydrodynamic dust Jacobian.
  !
  ! Parameters
  ! ----------
  ! area(Nr) : Radial grid annulus area
  ! D(Nr, Nm) : Dust diffusivity
  ! r(Nr) : Radial grid cell centers
  ! ri(Nr+1) : Radial grid cell interfaces
  ! SigmaGas(Nr) : Gas surface density
  ! v(Nr) : Radial velocity
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! A(Nr) : sub-diagonal, A(1) not used
  ! B(Nr) : diagonal
  ! C(Nr) : super-diagoanl, C(Nr) not used

  use constants, only: twopi

  implicit none

  double precision, intent(in)  :: area(Nr)
  double precision, intent(in)  :: D(Nr, Nm)
  double precision, intent(in)  :: r(Nr)
  double precision, intent(in)  :: ri(Nr+1)
  double precision, intent(in)  :: v(Nr, Nm)
  double precision, intent(in)  :: SigmaGas(Nr)
  double precision, intent(out) :: A(Nr, Nm)
  double precision, intent(out) :: B(Nr, Nm)
  double precision, intent(out) :: C(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer :: ir
  integer :: i
  double precision :: Di(Nr+1)
  double precision :: Fa(Nr+1)
  double precision :: h(Nr)
  double precision :: hi(Nr+1)
  double precision :: vi(Nr+1)
  double precision :: vim(Nr+1)
  double precision :: vip(Nr+1)
  double precision :: Vinv(Nr)
  double precision :: w(Nr)

  ! Helper quantity
  h(:) = SigmaGas(:) * r(:)
  call interp1d(ri, r, h, hi, Nr)

  ! Initialization
  A(:, :)    = 0.d0
  B(:, :)    = 0.d0
  C(:, :)    = 0.d0

  ! Grid cell volumes and distances
  do ir=1, Nr
      Vinv(ir) = twopi / area(ir)
      w(ir) = r(ir+1) - r(ir)
  end do

  do i=1, Nm

    ! Interface velocities
    call interp1d(ri, r, v(:, i), vi, Nr)
    vim(:)  = min(vi(:), 0.d0)
    vip(:)  = max(0.d0, vi(:))
    ! Interface diffusivity
    call interp1d(ri, r, D(:, i), Di, Nr)

    do ir=2, Nr-1

      ! Advection terms
      A(ir, i) = A(ir, i) + vip(ir)   * r(ir-1)
      B(ir, i) = B(ir, i) - vip(ir+1) * r(ir)   + vim(ir)   * r(ir)
      C(ir, i) = C(ir, i)                       - vim(ir+1) * r(ir+1)

      ! Diffusion terms
      A(ir, i) = A(ir, i) + Di(ir)   * hi(ir)   / ( w(ir-1) * h(ir-1) ) * r(ir-1)
      B(ir, i) = B(ir, i) - Di(ir)   * hi(ir)   / ( w(ir-1) * h(ir)   ) * r(ir)
      B(ir, i) = B(ir, i) - Di(ir+1) * hi(ir+1) / ( w(ir)   * h(ir)   ) * r(ir)
      C(ir, i) = C(ir, i) + Di(ir+1) * hi(ir+1) / ( w(ir)   * h(ir+1) ) * r(ir+1)
        
    end do

    ! Dividing by grid cell volume
    A(:, i) = A(:, i) * Vinv(:)
    B(:, i) = B(:, i) * Vinv(:)
    C(:, i) = C(:, i) * Vinv(:)

  end do

end subroutine jacobian_hydrodynamic_generator


double precision function kdelta(i, j)
  ! Function returns the Kronecker delta.
  !
  ! Parameters
  ! ----------
  ! i : integer
  ! j : integer
  !
  ! Returns
  ! -------
  ! kdelta : float
  !   1. if i==j
  !   0. else

  if(i == j) then
    kdelta = 1.d0
  else
    kdelta = 0.d0
  end if

end function kdelta


subroutine kernel(a, H, Sigma, SigmaFloor, vrel, K, Nr, Nm)
  ! Subroutine calculates the vertically integrated collision kernel.
  ! Has to be multiplied with fragmentation/sticking probabilities
  ! for fragmentation/sticking rates.
  !
  ! Parameters
  ! ----------
  ! a(Nr, Nm) : Particle sizes
  ! H(Nr, Nm) : Particle scale heights
  ! Sigma(Nr, Nm) : Surface density
  ! SigmaFloor(Nr, Nm) : Floor value of surface density
  ! vfrag(Nr) : Fragmentation velocities
  ! vrel(Nr, Nm, Nm) : Relative collision velocities
  ! Nr : Number of radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! K(Nr, Nm, Nm) : Collision kernel

  use constants, only: pi

  implicit none

  double precision, intent(in)  :: a(Nr, Nm)
  double precision, intent(in)  :: H(Nr, Nm)
  double precision, intent(in)  :: Sigma(Nr, Nm)
  double precision, intent(in)  :: SigmaFloor(Nr, Nm)
  double precision, intent(in)  :: vrel(Nr, Nm, Nm)
  double precision, intent(out) :: K(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  integer          :: ir
  integer          :: i
  integer          :: j
  double precision :: kdelta

  ! Initialization
  K(:, :, :) = 0.d0

  do ir=2, Nr-1
    do i=1, Nm
      if(Sigma(ir, i) .LT. SigmaFloor(ir, i)) cycle
      do j=1, i
        if(Sigma(ir, j) .LT. SigmaFloor(ir, j)) cycle
        K(ir, j, i) = (1.d0 - 0.5d0*kdelta(j, i)) * pi * (a(ir, j) + a(ir, i))**2 * vrel(ir, j, i) &
          & / sqrt( 2.d0 * pi * ( H(ir, j)**2 + H(ir, i)**2 ) )
      end do
    end do
  end do

end subroutine kernel


subroutine pfrag(vrel, vfrag, pf, Nr, Nm)
  ! Subroutine calculates the fragmentation probability.
  ! It is assuming a Maxwell-Boltzmann velocity distribution.
  ! 
  ! Parameters
  ! ----------
  ! vrel(Nr, Nm, Nm) : Relative velocity
  ! vfrag(Nr) : Fragmentation velocity
  ! Nr : Number or radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! pf(Nr, Nm, Nm) : Fragmentation probability in [0, 1]
  !
  ! Notes
  ! -----
  ! The sticking probability is ps = 1 - pf

  implicit none

  double precision, intent(in)  :: vrel(Nr, Nm, Nm)
  double precision, intent(in)  :: vfrag(Nr)
  double precision, intent(out) :: pf(Nr, Nm, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm
  
  double precision :: dum
  integer :: ir
  integer :: i
  integer :: j

  do i=1, Nm
    do j=1, i
      do ir=2, Nr-1
        dum = (vfrag(ir)/vrel(ir, j, i))**2
        pf(ir, j, i) = (1.5d0*dum + 1.d0) * exp(-1.5d0*dum)
        pf(ir, i, j) = pf(ir, j, i)
      end do
    end do
  end do

end subroutine pfrag


subroutine s_coag(cstick, cstick_ind, A, eps, klf, krm, phi, Kf, Ks, m, Sigma, SigmaFloor, S, Nr, Nm)
  ! Subroutine calculates the coagulation source terms.
  ! 
  ! Parameters
  ! ----------
  ! cstick(4, Nm, Nm) : Non-zero elements of sticking matrix
  ! cstick_ind(4, Nm, Nm) : Indices of non-zero elements of sticking matrix
  ! A(Nm, Nm) : Normalization factor of fragment distribution
  ! eps(Nm, Nm) : Distribution of remnant mass
  ! klf(Nm, Nm) : Index of largest fragment
  ! krm(Nm, Nm) : Smaller index of remnant mass
  ! phi(Nm, Nm) : Fragment distribution
  ! Kf(Nr, Nm, Nm) : Fragmentation kernel
  ! Ks(Nr, Nm, Nm) : Sticking kernel
  ! m(Nm) : Mass grid
  ! Sigma(Nr, Nm) : Dust surface density
  ! SigmaFloor(Nr, Nm) Floor value of surface density
  ! Nr : Number or radial grid cells
  ! Nm : Number of mass bins
  !
  ! Returns
  ! -------
  ! S(Nr, Nm) : Coagulation source terms of dust surface density

  implicit none

  double precision, intent(in)  :: cstick(4, Nm, Nm)
  integer,          intent(in)  :: cstick_ind(4, Nm, Nm)
  double precision, intent(in)  :: A(Nm, Nm)
  double precision, intent(in)  :: eps(Nm, Nm)
  integer,          intent(in)  :: klf(Nm, Nm)
  integer,          intent(in)  :: krm(Nm, Nm)
  double precision, intent(in)  :: phi(Nm, Nm)
  double precision, intent(in)  :: Kf(Nr, Nm, Nm)
  double precision, intent(in)  :: Ks(Nr, Nm, Nm)
  double precision, intent(in)  :: m(Nm)
  double precision, intent(in)  :: Sigma(Nr, Nm)
  double precision, intent(in)  :: SigmaFloor(Nr, Nm)
  double precision, intent(out) :: S(Nr, Nm)
  integer,          intent(in)  :: Nr
  integer,          intent(in)  :: Nm

  double precision :: As(Nm)
  double precision :: n(Nm)
  double precision :: Rf(Nm, Nm)
  double precision :: Rs
  integer :: ir
  integer :: i
  integer :: imax
  integer :: j
  integer :: jmin
  integer :: k
  integer :: nz
  integer :: p

  ! Determination of the mass bin range that contributes to
  ! full fragmentation
  p = Nm - minloc( m(:), 1, krm(:, Nm) .EQ. -1 )

  ! Initialization
  As(:) = 0.d0
  S(:, :) = 0.d0

  do ir=2, Nr-1
    ! Conversion to number density
    n(:) = Sigma(ir, :) / m(:)

    ! Largest active bin
    imax = min( maxloc(m(:), 1, Sigma(ir, :) .gt. SigmaFloor(ir, :)), Nm )

    ! Coagulation
    do i=1, imax
      do j=1, i
        Rs = Ks(ir, j, i) * n(j) * n(i)
        do nz=1, 4
          k = cstick_ind(nz, j, i) + 1
          if(k .eq. 0) cycle
          S(ir, k) = S(ir, k) + cstick(nz, j, i) * Rs
        end do
      end do
    end do

    ! FRAGMENTATION

    ! Resetting
    As(:) = 0.d0

    ! Adding the collision rates to the fragments distribution
    do i=1, imax
      do j=1, i
        Rf(j, i) = Kf(ir, j, i) * Sigma(ir, j)/m(j) * Sigma(ir, i)/m(i)
        k = klf(j, i) + 1
        As(k) = As(k) + A(j, i)*Rf(j, i)
      end do
    end do

    ! Adding fragment distribution
    do i=1, nm
      do j=i, nm
        S(ir, i) = S(ir, i) + As(j)*phi(j, i)/m(i)
      end do
    end do

    ! Negative terms and erosive remannts
    do i=1, imax
      ! Erosion
      do j=1, i-p-1
        k = krm(j, i) + 1
        ! It's better for mass conservation to distinguish both cases.
        if(k .eq. i-1) then
          S(ir, k)   = S(ir, k)   + eps(j, i) * Rf(j, i)
          S(ir, k+1) = S(ir, k+1) - eps(j, i) * Rf(j, i) ! <- Both terms in one == better mass conservation
          S(ir, j)   = S(ir, j)   - Rf(j, i)
        else
          S(ir, k)   = S(ir, k)   + eps(j, i)          * Rf(j, i)
          S(ir, k+1) = S(ir, k+1) + (1.d0 - eps(j, i)) * Rf(j, i)
          S(ir, i)   = S(ir, i)   - Rf(j, i)
          S(ir, j)   = S(ir, j)   - Rf(j, i)
        end if
      end do
      ! Full fragmentation (only negative terms)
      jmin = MAX(1, i-p)
      do j=jmin, i
        S(ir, i) = S(ir, i) - Rf(j, i)
        S(ir, j) = S(ir, j) - Rf(j, i)
      end do
    end do

    S(ir, :) = S(ir, :) * m(:)

  end do

end subroutine s_coag


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
        Shyd(ir, i) = 2.d0 * (Fi(ir, i)*ri(ir) - Fi(ir+1, i)*ri(ir+1)) / (ri(ir+1)**2 - ri(ir)**2)
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


doubleprecision function theta(x)
  ! Heaviside step function.
  !
  ! Parameters
  ! ----------
  ! x : float
  !
  ! Returns
  ! -------
  ! theta : float
  !   0. if x < 0
  !   1. if x ≥ 0

  implicit none

  double precision :: x

  theta = 0.d0
  if(x .GE. 0.d0) then
    theta = 1.d0
  end if

end function theta


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
        dum = abs(vdriftmax(ir) * (St2p1(ir, i) - St2p1(ir, j)) / (St2p1(ir, i) * St2p1(ir, j)))
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


subroutine vrel_ormel_cuzzi_2007(alpha, cs, mump, OmegaK, SigmaGas, &
  & St, vrel, Nr, Nm)
  ! Subroutine calculates the relative particle velocities due to turbulent motion
  ! accourding the prescription of Ormel & Cuzzi (2007).
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

end subroutine vrel_ormel_cuzzi_2007


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
        dum = abs(vrad(ir, j) - vrad(ir, i))
        vrel(ir, i, j) = dum
        vrel(ir, j, i) = dum
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