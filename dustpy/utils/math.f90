module math
   ! Module containing helper functions needed by other subroutines.

   implicit none
   
   contains
   
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
         !   1. if x >= 0
         implicit none
         double precision :: x
         theta = 0.d0
         if(x .GE. 0.d0) then
            theta = 1.d0
         end if
      end function theta
      
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
         implicit none
         integer :: i
         integer :: j
         kdelta = 0.d0
         if(i == j) then
           kdelta = 1.d0
         end if
      end function kdelta

end module math