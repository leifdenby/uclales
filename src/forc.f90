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
!
module forc

  use defs, only      : cp
  use radiation, only : d4stream
  !irina
  use rad_gcss, only : gcss_rad
  implicit none

   !character (len=5), parameter :: case_name = 'xxxx'
   !irina
   !initializes the large scale forcing
   integer, parameter    :: nls = 500
   real, dimension(nls)  :: t_ls=0.
   real, dimension(nls)  :: div_ls=0.
   real, dimension(nls)  :: sst_ls=0.
   real, dimension(nls)  :: ugeo_ls=0.
   real, dimension(nls)  :: vgeo_ls=0.

    
contains
  !
  ! -------------------------------------------------------------------
  ! subroutine forcings:  calls the appropriate large-scale forcings
  !irina
  subroutine forcings(time_in, cntlat, sst,div, case_name)

!irina
    use grid, only: nxp, nyp, nzp, zm, zt, dzt, dzm, dn0, iradtyp, liquid  &
         , a_rflx, a_sflx, albedo, a_tt, a_tp, a_rt, a_rp, a_pexnr, a_scr1 &
         , vapor, a_rpp, CCN, pi0, pi1, level, a_ut, a_up, a_vt, a_vp,a_theta,&
          a_lflxu, a_lflxd, a_sflxu, a_sflxd

    use mpi_interface, only : myid, appl_abort

!irina
    real, optional, intent (in) :: time_in, cntlat, sst,div

    character (len=5), intent (in) :: case_name
!irina
    real :: xref1, xref2
    integer :: i, j, k, kp1


!irina
    select case(iradtyp)
    case (1)
        call case_forcing(nzp,nxp,nyp,case_name,zt,dzt,dzm,a_tp,a_rp,a_tt,a_rt)
    case (2)
       select case(level)
       case(1) 
          call smoke_rad(nzp, nxp, nyp, dn0, a_rflx, zm, dzt,a_tt,a_rp)
       case(2)
          call gcss_rad(nzp, nxp, nyp, cntlat, time_in, case_name, div, sst, liquid, dn0,   &
              a_rflx, a_sflx, zt, zm, dzt, a_tt, a_tp, a_rt, a_rp)
       case(3)
          call gcss_rad(nzp, nxp, nyp, cntlat, time_in, case_name, div, sst, liquid, dn0,   &
              a_rflx, a_sflx, zt, zm, dzt, a_tt, a_tp, a_rt, a_rp)
       end select
    case (3)
       call bellon(nzp, nxp, nyp, a_rflx, a_sflx, zt, dzt, dzm, a_tt, a_tp&
            ,a_rt, a_rp, a_ut, a_up, a_vt, a_vp)
    case (4)
       if (present(time_in) .and. present(cntlat) .and. present(sst)) then
          !irina 
          !a_scr1 = a_theta/a_pexnr 
          if (level == 3) then
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,   &
                  dn0, pi0, pi1, dzt, a_pexnr, a_theta, vapor, liquid, a_tt,&
                  a_rflx, a_sflx, a_lflxu, a_lflxd,a_sflxu,a_sflxd, albedo, rr=a_rpp)
            !old      
            ! call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,   &
            !      dn0, pi0, pi1, dzt, a_pexnr, a_scr1, vapor, liquid, a_tt,&
            !      a_rflx, a_sflx, albedo, rr=a_rpp)
          else
             xref1 = 0.
             xref2 = 0.
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,    &
                  dn0, pi0, pi1, dzt, a_pexnr, a_theta, vapor, liquid, a_tt, &
                  a_rflx, a_sflx,a_lflxu, a_lflxd,a_sflxu,a_sflxd,albedo)
             !call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,    &
             !     dn0, pi0, pi1, dzt, a_pexnr, a_scr1, vapor, liquid, a_tt, &
             !     a_rflx, a_sflx, albedo)
             xref1 = xref1 + a_sflx(nzp,3,3)/albedo(3,3)
             xref2 = xref2 + a_sflx(nzp,3,3)
             albedo(3,3) = xref2/xref1
          end if
!irina 
          !
          ! subsidence
          !
       do j=3,nyp-2
          do i=3,nxp-2
           if (div /= 0.) then
             do k=2,nzp-2
                kp1 = k+1
                a_tt(k,i,j) = a_tt(k,i,j) + &
                        div*zt(k)*(a_tp(kp1,i,j)-a_tp(k,i,j))*dzt(k)
                a_rt(k,i,j)=a_rt(k,i,j) + &
                        div*zt(k)*(a_rp(kp1,i,j)-a_rp(k,i,j))*dzt(k)
             end do
           end if
          enddo
       enddo
       else
          if (myid == 0) print *, '  ABORTING: inproper call to radiation'
          call appl_abort(0)
       end if

    end select 

  end subroutine forcings
  !
  ! -------------------------------------------------------------------
  ! subroutine smoke_rad:  call simple radiative parameterization for 
  ! the smoke cloud
  !
  subroutine smoke_rad(n1,n2,n3,dn0,flx,zm,dzt,tt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: zm(n1),dzt(n1),dn0(n1),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3)
    real, intent (out)  :: flx(n1,n2,n3)
    real, parameter     :: xka= 50.0, fr0=60.0

    integer :: i,j,k, km1, ki
    real    :: smoke(n2,n3)

    smoke=0.
    do j=3,n3-2
       do i=3,n2-2
          ki = n1
          do k=1,n1
             km1=max(1,k-1)
             smoke(i,j)=smoke(i,j)+max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
          enddo

          do k=2,n1
             km1=max(2,k-1)
             smoke(i,j)=smoke(i,j)-max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=fr0*exp(-1.*xka*smoke(i,j))
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzt(k)/(dn0(k)*cp)
          enddo
       enddo
    enddo

  end subroutine smoke_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine case_forcing: adjusts tendencies according to a specified
  ! large scale forcing.  Normally case (run) specific.
  !
  subroutine case_forcing(n1,n2,n3,case_name,zt,dzt,dzm,tl,rt,tt,rtt)

    use mpi_interface, only : pecount, double_scalar_par_sum,myid, appl_abort
    use stat, only : get_zi

    integer, intent (in):: n1,n2, n3
    real, dimension (n1), intent (in)          :: zt, dzt, dzm
    real, dimension (n1,n2,n3), intent (in)    :: tl, rt
    real, dimension (n1,n2,n3), intent (inout) :: tt, rtt
    character (len=5), intent (in) :: case_name

    integer :: i,j,k,kp1
    real, dimension (n1) :: sf
    real, parameter :: zmx_sub = 2260. ! originally 2260.

    real (kind=8) :: zig, zil
    real          :: zibar

    select case (trim(case_name))
    case('rico')
       
       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < zmx_sub) then
             sf(k) =  -0.005*zt(k)/zmx_sub
          else
             sf(k) =  -0.005 
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! subsidence
                ! 
                kp1 = k+1
                tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                !
                ! temperature advection and radiative cooling
                !
                tt(k,i,j) = tt(k,i,j)  - 2.5/86400.
                !
                ! moisture advection
                !
                if (zt(k) <= 2980.) then
                   rtt(k,i,j) = rtt(k,i,j)  - (1. -  1.3456*zt(k)/2980.)/8.64e7
                else
                   rtt(k,i,j) = rtt(k,i,j)  + .3456/8.64e7
                end if
             enddo
          enddo
       enddo

    case ('bomex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < 1500.) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065  + 0.0065*(zt(k)-1500.)/600.)
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < 1500.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5
                else if (zt(k) < 2000.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-1500.)*1.e-3)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < 300.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8
                elseif (zt(k) < 500.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8*(1.- (zt(k)-300.)/200.)
                end if
             enddo
          enddo
       enddo
    case ('atex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       zil = get_zi (n1, n2, n3, 2, rt, dzm, zt, 6.5e-3)
       call double_scalar_par_sum(zil,zig)
       zibar = real(zig/pecount)

       do k=2,n1-2
          if (zt(k) < zibar) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065*(1 - (zt(k)-zibar)/300.))
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < zibar) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1. + (1.- zt(k)/zibar)/2.)
                else if (zt(k) < zibar+300.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-zibar)/300.)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < zibar) rtt(k,i,j) = rtt(k,i,j)  - 1.5e-8
             enddo
          enddo
       enddo
    case default
       if (myid == 0) print *, '  ABORTING: inproper call to radiation'
       call appl_abort(0)
    end select

  end subroutine case_forcing
  !
  ! -------------------------------------------------------------------
  ! subroutine bellon_rad:  call simple radiative parameterization
  !
  subroutine bellon(n1,n2,n3,flx,sflx,zt,dzt,dzm,tt,tl,rtt,rt, ut,u,vt,v)

    integer, intent (in) :: n1,n2, n3

    real, dimension (n1), intent (in)            :: zt, dzt, dzm
    real, dimension (n1, n2, n3), intent (inout) :: tt, tl, rtt, rt, ut,u,vt,v
    real,  dimension (n1, n2, n3), intent (out)  :: flx, sflx
    real, parameter      :: w0= 7.5e-3, H=1000., Qrate = 2.5/86400.

    integer :: i,j,k,kp1
    real    :: grad,wk

    do j=3,n3-2
       do i=3,n2-2
          !
          ! subsidence
          !
          flx(1,i,j)  = 0.
          sflx(1,i,j) = 0.
          do k=2,n1-2
             kp1 = k+1
             wk = w0*(1.-exp(-zt(k)/H))
             grad = Qrate/wk
             flx(k,i,j)  = wk*((tl(kp1,i,j)-tl(k,i,j))*dzt(k)-grad)
             sflx(k,i,j) = wk*((rt(kp1,i,j)-rt(k,i,j))*dzt(k)-grad)
             tt(k,i,j) = tt(k,i,j) + flx(k,i,j)
             rtt(k,i,j)=rtt(k,i,j) + &
                  wk*(rt(kp1,i,j)-rt(k,i,j))*dzt(k)
             ut(k,i,j) =  ut(k,i,j) + &
                  wk*(u(kp1,i,j)-u(k,i,j))*dzm(k)
             vt(k,i,j) =  vt(k,i,j) + &
                  wk*(v(kp1,i,j)-v(k,i,j))*dzm(k)
          end do
          flx(n1,  i,j)  = 0.
          flx(n1-1,i,j)  = 0.
          sflx(n1,  i,j) = 0.
          sflx(n1-1,i,j) = 0.
       enddo
    enddo

  end subroutine bellon
  

end module forc