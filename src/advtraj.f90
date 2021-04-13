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
  use grid, only          : a_advtrc_xr, a_advtrc_xi, a_advtrc_yr, a_advtrc_yi, &
                            a_advtrc_z, ladvtrc
  use grid, only          : nxp, nyp, nzp, xt, yt, zt, deltax, deltay
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
    use mpi_interface, only : myid, appl_abort, nxg, nyg, nxprocs, nyprocs
    ! nxg, nyg are the actual number if grid points in the interior of the domain
    ! and so they should be used for setting up the position scalars

    integer :: i, j, k, i_g, j_g

    if (ladvtrc) then
      do j=1,nyp
        do i=1,nxp
          do k=1,nzp
            ! have to calculate the global grid position, each local
            ! decomposition has an extra of two ghost cells either side
            ! these indecies will be negative and larger than the number of
            ! cells in the x/y-direction, but this is ok since the sin/cos
            ! functions are periodic and will wrap
            i_g = xt(i) / deltax
            j_g = yt(j) / deltay
            a_advtrc_xr(k,i,j) = cos(2.0 * pi * i_g / nxg)
            a_advtrc_xi(k,i,j) = sin(2.0 * pi * i_g / nxg)
            a_advtrc_yr(k,i,j) = cos(2.0 * pi * j_g / nyg)
            a_advtrc_yi(k,i,j) = sin(2.0 * pi * j_g / nyg)
            a_advtrc_z(k,i,j) = k
          end do
        end do
      end do
    endif
  end subroutine

end module modadvtraj
