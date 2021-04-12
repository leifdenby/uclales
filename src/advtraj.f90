!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
! Doxygen:
!> Lagrangian Particle using advected fields
!! Tracks massless trajectories using passive tracers on an Eulerian LES grid
!>
!! \author Leif Denby


module modadvtraj

  !--------------------------------------------------------------------------
  ! module modadvtraj
  !--------------------------------------------------------------------------
  use defs, only          : long, pi
  use mcrp, only          : lpartdrop                   !< Switch for rain drop like advtraj
  use grid, only          : a_advtrc_xr, a_advtrc_xi, a_advtrc_yr, a_advtrc_yi, &
                            a_advtrc_z, ladvtrc
  use grid, only          : nxp, nyp, nzp
  implicit none
  public :: reinit_advtraj

contains
  !
  !--------------------------------------------------------------------------
  ! Subroutine advtraj
  !> Set up the "position scalars" for the advected scalars trajectory method
  !> This should be every time these 3D fields have been saved to file, i.e.
  !> after the analysis file has been written
  !--------------------------------------------------------------------------
  !
  subroutine reinit_advtraj
    integer :: i, j, k
    if (ladvtrc) then
      do j=1,nyp
        do i=1,nxp
          do k=1,nzp
            a_advtrc_xr(k,i,j) = cos(2.0 * pi * i / nxp)
            a_advtrc_xi(k,i,j) = sin(2.0 * pi * i / nxp)
            a_advtrc_yr(k,i,j) = cos(2.0 * pi * j / nyp)
            a_advtrc_yi(k,i,j) = sin(2.0 * pi * j / nyp)
            a_advtrc_z(k,i,j) = k
          end do
        end do
      end do
    endif
  end subroutine

end module modadvtraj
