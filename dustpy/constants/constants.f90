module constants

    implicit none

    double precision :: au = 1.4959787070000d13
    double precision :: G = 6.674299999999999d-08
    double precision :: k_B = 1.380649d-16
    double precision :: m_p = 1.67262192369d-24
    double precision :: M_sun = 1.988409870698051d33
    double precision :: pi = ACOS(-1.d0)
    double precision :: R_sun = 6.9570000000d10
    double precision :: sigma_H2 = 2.d-15
    double precision :: sigma_sb = 5.6703744191844314d-05
    double precision :: year = 3.1557600d7

    ! This is only used by the fragmentation algorithm
    ! and should not be changed by the user. It is not
    ! automatically passed to the Python module.
    ! Particles fully fragment if j >= i-p.
    integer :: p = 0 

end module constants