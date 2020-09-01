subroutine interp1d(xi, x, y, yi, N)
    ! Function linearly interpolates field from grid cell centers
    ! onto grid cell interfaces. First and last value are
    ! linearly extrapolated.
    !
    ! Parameters
    ! ----------
    ! xi(N+1) : Grid cell interfaces
    ! x(N) : Grid cell centers
    ! y(N) : Value at grid cell centers
    ! N : Number of grid cells
    !
    ! Returns
    ! -------
    ! yi : Value at grid cell interfaces

    implicit none

    double precision, intent(in)  :: x(N)
    double precision, intent(in)  :: xi(N+1)
    double precision, intent(in)  :: y(N)
    double precision, intent(out) :: yi(N+1)
    integer,          intent(in)  :: N

    double precision :: m(N-1)
    integer :: i

    do i=1, N-1
        m(i) = (y(i+1) - y(i)) / (x(i+1) - x(i))
    end do

    do i=2, N
        yi(i) = y(i-1) + m(i-1)*(xi(i)-x(i-1))
    end do

    yi(1) = y(1) + m(1)*(xi(1)-x(1))
    yi(N+1) = y(N) + m(N-1)*(xi(N+1)-x(N))

end subroutine interp1d