module advance
  use library
  use transform
  use fluxes
  use source
  use initialisation
  use boundary
  use recon
  use local
  use flow_operations
  use communications
  use implicit_time
  use implicit_fluxes
  use declaration
  use io
  use moodr
  implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!subroutine called to calculate cfl number depending on the scheme!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!employed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_cfl(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, l, kmaxe, j, ingtmax, ingtmin, whgu, whgl, srf
    real::suvi, suv3, maxu, minu
    real::ccfl, veln, agrt
    real, dimension(1:nof_variables)::leftv, rightv
    real, dimension(1:nof_variables)::srf_speed
    real::mp_pinfl, gammal
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    kmaxe = xmpielrank(n)
    ccfl = (cfl/3.0d0)
    dt = tolbig
    if (itestcase .lt. 3) then
      do i = 1, kmaxe
        veln = max(abs(lamx), abs(lamy), abs(lamz))
        if (dg .eq. 1) then
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1)))
        else
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln))))
        end if
      end do
    end if

    if (itestcase .eq. 3) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        if (multispecies .eq. 1) then
          agrt = sqrt((leftv(5) + mp_pinfl)*gammal/leftv(1))
        else
          agrt = sqrt(leftv(5)*gamma/leftv(1))
        end if
        if (rframe .eq. 0) then
          veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
        end if
        if (srfg .eq. 1) then
          pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
          poy(1:3) = srf_velocity
          srf_speed = zero
          srf_speed(2:4) = vect_function(pox, poy)
          veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
        end if
        if (mrf .eq. 1) then
          srf = ilocal_recon3(i)%mrf
          if (ilocal_recon3(i)%mrf .eq. 0) then
            veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
          else
            pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
            pox = pox - ilocal_recon3(i)%mrf_origin
            poy(1:3) = ilocal_recon3(i)%mrf_velocity
            srf_speed = zero
            srf_speed(2:4) = vect_function(pox, poy)
            veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
          end if
        end if
        if (dg .eq. 1) then
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1)))
        else
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln))))
        end if
      end do
    end if
    if (itestcase .eq. 4) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        rightv(1:nof_variables) = leftv(1:nof_variables)
        call sutherland(n, leftv, rightv, viscl, laml)
        agrt = sqrt(leftv(5)*gamma/leftv(1))
        if (rframe .eq. 0) then
          veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
        end if
        if (srfg .eq. 1) then
          pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
          poy(1:3) = srf_velocity
          srf_speed = zero
          srf_speed(2:4) = vect_function(pox, poy)
          veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
        end if
        if (mrf .eq. 1) then
          srf = ilocal_recon3(i)%mrf
          if (ilocal_recon3(i)%mrf .eq. 0) then
            veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
          else
            pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
            pox(1:3) = pox(1:3) - ilocal_recon3(i)%mrf_origin(1:3)
            poy(1:3) = ilocal_recon3(i)%mrf_velocity(1:3)
            srf_speed = zero
            srf_speed(2:4) = vect_function(pox, poy)
            veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
          end if
        end if
        if (turbulence .eq. 1) then
          if (turbulencemodel .eq. 1) then
            turbmv(1) = u_ct(i)%val(1, 1); turbmv(2) = u_ct(i)%val(1, 1); 
            eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
            call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            laml(1) = laml(1) + laml(3)
            viscl(1) = viscl(1) + viscl(3)
          end if
        end if
        if (dg .eq. 1) then
          dt=min(dt,(ccfl/(2*iorder+1))*(ielem(n,i)%minedge/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem(n,i)%minedge)))))
        else
          ielem(n, i)%viscx = viscl(1)/visc
          dt = min(dt, ccfl*(1.0d0/((abs(veln)/((ielem(n, i)%minedge))) + (0.5d0*(laml(1) + viscl(1))/((ielem(n, i)%minedge))**2))))
!         dt=min(dt,(ccfl)*(ielem(n,i)%minedge/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*(1.0d0/ielem(n,i)%minedge)))))
        end if
      end do
    end if

    return

  end subroutine calculate_cfl

  subroutine calculate_cfll(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, l, kmaxe, j, ingtmax, ingtmin, whgu, whgl, srf
    real::suvi, suv3, maxu, minu
    real::ccfl, veln, agrt
    real, dimension(1:nof_variables)::leftv, rightv
    real, dimension(1:nof_variables)::srf_speed
    real::mp_pinfl, gammal
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    kmaxe = xmpielrank(n)
    ccfl = (cfl/3.0d0)
    if (itestcase .lt. 3) then
      do i = 1, kmaxe
        veln = max(abs(lamx), abs(lamy), abs(lamz))
        ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))
      end do
    end if

    if (itestcase .eq. 3) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        if (multispecies .eq. 1) then
          agrt = sqrt((leftv(5) + mp_pinfl)*gammal/leftv(1))
        else
          agrt = sqrt(leftv(5)*gamma/leftv(1))
        end if
        if (rframe .eq. 0) then
          veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
        end if
        if (srfg .eq. 1) then
          pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
          poy(1:3) = srf_velocity
          srf_speed = zero
          srf_speed(2:4) = vect_function(pox, poy)
          veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
        end if
        if (mrf .eq. 1) then
          srf = ilocal_recon3(i)%mrf
          if (ilocal_recon3(i)%mrf .eq. 1) then
            veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
          else
            pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
            pox = pox - ilocal_recon3(i)%mrf_origin
            poy(1:3) = ilocal_recon3(i)%mrf_velocity
            srf_speed = zero
            srf_speed(2:4) = vect_function(pox, poy)
            veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
          end if
        end if
        if (dg .eq. 1) then
          ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1))
        else
          ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))
        end if
      end do
    end if

    if (itestcase .eq. 4) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        rightv(1:nof_variables) = leftv(1:nof_variables)
        call sutherland(n, leftv, rightv, viscl, laml)
        agrt = sqrt(leftv(5)*gamma/leftv(1))
        if (srfg .eq. 0 .and. mrf .eq. 0) then
          veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
        end if
        if (srfg .eq. 1) then
          pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
          poy(1:3) = srf_velocity
          srf_speed = zero
          srf_speed(2:4) = vect_function(pox, poy)
          veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
        end if
        if (mrf .eq. 1) then
          srf = ilocal_recon3(i)%mrf
          if (ilocal_recon3(i)%mrf .eq. 0) then
            veln = max(abs(leftv(2)), abs(leftv(3)), abs(leftv(4))) + agrt
          else
            pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
            pox(1:3) = pox(1:3) - ilocal_recon3(i)%mrf_origin(1:3)
            poy(1:3) = ilocal_recon3(i)%mrf_velocity
            srf_speed = zero
            srf_speed(2:4) = vect_function(pox, poy)
            veln = max(abs(leftv(2) - srf_speed(2)), abs(leftv(3) - srf_speed(3)), abs(leftv(4) - srf_speed(4))) + agrt
          end if
        end if
        if (turbulence .eq. 1) then
          if (turbulencemodel .eq. 1) then
            turbmv(1) = u_ct(i)%val(1, 1); turbmv(2) = u_ct(i)%val(1, 1); 
            eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
            call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            laml(1) = laml(1) + laml(3)
            viscl(1) = viscl(1) + viscl(3)
          end if
        end if
        if (dg .eq. 1) then
            ielem(n,i)%dtl=(ccfl/(2*iorder+1))*(ielem(n,i)%minedge/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem(n,i)%minedge))))
        else
      ielem(n, i)%dtl = ccfl*(1.0d0/((abs(veln)/((ielem(n, i)%minedge))) + (0.5d0*(laml(1) + viscl(1))/((ielem(n, i)%minedge))**2)))
        end if
      end do
    end if
    return
  end subroutine calculate_cfll

  subroutine calculate_cfl2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, l, kmaxe, j, ingtmax, ingtmin, whgu, whgl
    real::suvi, suv3, maxu, minu, sum_dt1, sum_dt2
    real::ccfl, veln, agrt, lamxl, lamyl
    real, dimension(1:nof_variables)::leftv, rightv
    real, dimension(1:nof_variables)::srf_speed
    real::mp_pinfl, gammal
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    kmaxe = xmpielrank(n)
    ccfl = (cfl/2.0d0)
    dt = tolbig
    if (itestcase .lt. 3) then
      do i = 1, kmaxe
        if (initcond .eq. 3) then
          lamxl = -ielem(n, i)%yyc + 0.5d0
          lamyl = ielem(n, i)%xxc - 0.5
        else
          lamxl = lamx
          lamyl = lamy
        end if
        veln = max(abs(lamxl), abs(lamyl))
        if (dg .eq. 1) then
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1)))
        else
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln))))
        end if
      end do
    end if
    if (itestcase .eq. 3) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        if (multispecies .eq. 1) then
          agrt = sqrt((leftv(4) + mp_pinfl)*gammal/leftv(1))
        else
          agrt = sqrt(leftv(4)*gamma/leftv(1))
        end if
        veln = max(abs(leftv(2)), abs(leftv(3))) + agrt
        if (dg .eq. 1) then
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1)))
        else
          dt = min(dt, ccfl*((ielem(n, i)%minedge)/(abs(veln))))
        end if
      end do
    end if
    if (itestcase .eq. 4) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        rightv(1:nof_variables) = leftv(1:nof_variables)
        call sutherland2d(n, leftv, rightv, viscl, laml)
        agrt = sqrt(leftv(4)*gamma/leftv(1))
        veln = max(abs(leftv(2)), abs(leftv(3))) + agrt
        if (turbulence .eq. 1) then
        if (turbulencemodel .eq. 1) then
          turbmv(1) = u_ct(i)%val(1, 1); turbmv(2) = u_ct(i)%val(1, 1); 
          eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
          call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
          laml(1) = laml(1) + laml(3)
          viscl(1) = viscl(1) + viscl(3)
        end if
        end if
        if (dg .eq. 1) then
        dt=min(dt,(ccfl/(2*iorder+1))*(ielem(n,i)%minedge/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem(n,i)%minedge)))))
        else
          dt = min(dt, ccfl*(1.0d0/((abs(veln)/((ielem(n, i)%minedge))) + (0.5d0*(laml(1) + viscl(1))/((ielem(n, i)%minedge))**2))))
        end if
      end do
    end if
    return

  end subroutine calculate_cfl2d
  subroutine calculate_cfll2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, l, kmaxe, j, ingtmax, ingtmin, whgu, whgl
    real::suvi, suv3, maxu, minu
    real::ccfl, veln, agrt
    real, dimension(1:nof_variables)::leftv, rightv
    real, dimension(1:nof_variables)::srf_speed
    real::mp_pinfl, gammal, mp_pinfr, gammar
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    kmaxe = xmpielrank(n)
    ccfl = (cfl/2.0d0)
    if (itestcase .lt. 3) then
      do i = 1, kmaxe
        veln = max(abs(lamx), abs(lamy))
        ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))
      end do
    end if
    if (itestcase .eq. 3) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        agrt = sqrt(leftv(4)*gamma/leftv(1))
        veln = max(abs(leftv(2)), abs(leftv(3))) + agrt
        if (dg .eq. 1) then
          ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))*(1.0d0/(2*iorder + 1))
        else
          ielem(n, i)%dtl = ccfl*((ielem(n, i)%minedge)/(abs(veln)))
        end if
      end do
    end if
    if (itestcase .eq. 4) then
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        rightv(1:nof_variables) = leftv(1:nof_variables)
        call sutherland2d(n, leftv, rightv, viscl, laml)
        agrt = sqrt(leftv(4)*gamma/leftv(1))
        veln = max(abs(leftv(2)), abs(leftv(3))) + agrt
        if (turbulence .eq. 1) then
          if (turbulencemodel .eq. 1) then
            turbmv(1) = u_ct(i)%val(1, 1); turbmv(2) = u_ct(i)%val(1, 1); 
            eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
            call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            laml(1) = laml(1) + laml(3)
            viscl(1) = viscl(1) + viscl(3)
          end if
        end if
        if (dg .eq. 1) then
          ielem(n,i)%dtl=(ccfl/(2*iorder+1))*(ielem(n,i)%minedge/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem(n,i)%minedge))))
        else
          ielem(n, i)%dtl = ccfl*(1.0d0/((abs(veln)/((ielem(n, i)%minedge))) + (0.5d0*(laml(1) + viscl(1))/((ielem(n, i)%minedge))**2)))
        end if
      end do
    end if
    return

  end subroutine calculate_cfll2d

! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!subroutine called to advance solution by one time step size!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine runge_kutta3_mood(n)
    implicit none
    integer, intent(in)::n
    integer::iavr, nvar, i, kmaxe, inds
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    if (mood .eq. 1) then
      inds = 4
    else
      inds = 1
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      if (mood .eq. 1) then
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      end if
      u_c(i)%val(inds, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (mood .eq. 1) then
      call mood_operator_2(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
          u_c(i)%val(inds, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
 u_c(i)%val(inds, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(to4*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo4*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((oo4))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    if (mood .eq. 1) then
      call mood_operator_2(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
 u_c(i)%val(inds, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
          u_c(i)%val(1, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
        call vortexcalc(n)
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (mood .eq. 1) then
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      end if
   u_c(i)%val(inds, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
    end do
    if (mood .eq. 1) then
      call mood_operator_2(n)
      !$omp do
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
   u_c(i)%val(inds, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(1, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(1,1:turbulenceequations+passivescalar)=((oo3)*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+((to3)*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(((to3))*&
                                                           ((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta3_mood
  subroutine runge_kutta3(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        u_c(i)%valdg(2, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(2,1:nof_variables,:) - dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        u_c(i)%valdg(3, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:)=to4*u_c(i)%valdg(2,1:nof_variables,:) + oo4*u_c(i)%valdg(3,1:nof_variables,:) - oo4*dt* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    u_c(i)%val(1, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(to4*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo4*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((oo4))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      if (dg .eq. 1) then
!         u_c(i)%valdg(1,1:nof_variables,:)=oo3*u_c(i)%valdg(2,1:nof_variables,:) + to3*u_c(i)%valdg(1,1:nof_variables,:) - to3*dt*transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
!         call dgemm('n','n',num_dg_dofs,nof_variables,num_dg_dofs,alpha,m_1(i)%val(1:num_dg_dofs,1:num_dg_dofs),num_dg_dofs,&
!           rhs(i)%valdg(1:num_dg_dofs,1:nof_variables),&
!         num_dg_dofs,beta,rhs(i)%sol_mm_dg,num_dg_dofs)
          rhs(i)%sol_mm_dg(1:num_dg_dofs,1:nof_variables)=matmul(m_1(i)%val(1:num_dg_dofs,1:num_dg_dofs),rhs(i)%valdg(1:num_dg_dofs,1:nof_variables))
          u_c(i)%valdg(1,1:nof_variables,:)=oo3*u_c(i)%valdg(2,1:nof_variables,:) + to3*u_c(i)%valdg(1,1:nof_variables,:) - to3*dt * transpose(rhs(i)%sol_mm_dg)
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(1, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_ct(i)%val(1,1:turbulenceequations+passivescalar)=((oo3)*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+((to3)*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(((to3))*&
                                                           ((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta3
  subroutine runge_kutta1(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, kx
    real::avrgs, oovolume, to4, oo4, to3, oo3
    real::t1, t2, t3
    kmaxe = xmpielrank(n)
    call call_flux_subroutines_3d
    do i = 1, kmaxe
    if (dg .eq. 1) then
    if ((u_c(i)%valdg(1, 1, 1) .ne. u_c(i)%valdg(1, 1, 1))) then
      if (n .eq. 0) print *, 'stopping because nans1'
      stop ! stop if nans
    end if
    end if
    end do
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(1,1:nof_variables,:) - dt* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))!*oovolume
      else
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(1,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta1
  subroutine runge_kutta2(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(1,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(1,1:nof_variables)=(oo2*u_c(i)%val(2,1:nof_variables))+(oo2*u_c(i)%val(1,1:nof_variables))-(dt*oo2*(rhs(i)%val(1:nof_variables)*oovolume))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(oo2*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo2*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(dt*oo2*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta2
  subroutine runge_kutta5(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(2, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(2,1:nof_variables,:) - ielem(n,i)%dtl * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) - (ielem(n, i)%dtl*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(1,1:turbulenceequations+passivescalar)-(ielem(n,i)%dtl*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:)=(oo2*u_c(i)%valdg(2,1:nof_variables,:)) +(oo2*u_c(i)%valdg(1,1:nof_variables,:))- (ielem(n,i)%dtl *oo2* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables))))
      else
  u_c(i)%val(1,1:nof_variables)=(oo2*u_c(i)%val(2,1:nof_variables))+(oo2*u_c(i)%val(1,1:nof_variables))-(ielem(n,i)%dtl*oo2*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(oo2*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo2*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(ielem(n,i)%dtl*oo2*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta5
  subroutine runge_kutta5_2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%valdg(2, :, :) = u_c(i)%valdg(1, :, :)
        u_c(i)%valdg(1, :, :) = u_c(i)%valdg(2, :, :) - (ielem(n, i)%dtl*transpose(matmul(m_1(i)%val(:, :), rhs(i)%valdg(:, :))))!*oovolume
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (ielem(n, i)%dtl*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(1,1:turbulenceequations+passivescalar)-(ielem(n,i)%dtl*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%valdg(1,:,:)=(oo2*u_c(i)%valdg(2,:,:))+(oo2*u_c(i)%valdg(1,:,:))-(ielem(n,i)%dtl* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,:))))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_c(i)%val(1,1:nof_variables)=(oo2*u_c(i)%val(2,1:nof_variables))+(oo2*u_c(i)%val(1,1:nof_variables))-(ielem(n,i)%dtl*oo2*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(oo2*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo2*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(ielem(n,i)%dtl*oo2*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta5_2d
  subroutine runge_kutta2_2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    call call_flux_subroutines_2d
    do i = 1, kmaxe
    if (dg .eq. 1) then
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%valdg(2, :, :) = u_c(i)%valdg(1, :, :)
      u_c(i)%valdg(1, :, :) = u_c(i)%valdg(2, :, :) - (dt*transpose(matmul(m_1(i)%val(:, :), rhs(i)%valdg(:, :))))!*oovolume
    else
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
    end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
    if (dg .eq. 1) then
      oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%valdg(1,:,:)=oo2*u_c(i)%valdg(2,:,:) + oo2*u_c(i)%valdg(1,:,:) - (oo2*dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,:))))!*oovolume
    else
      oovolume = 1.0d0/ielem(n, i)%totvolume
  u_c(i)%val(1,1:nof_variables)=(oo2*u_c(i)%val(2,1:nof_variables))+(oo2*u_c(i)%val(1,1:nof_variables))-(dt*oo2*(rhs(i)%val(1:nof_variables)*oovolume))
    end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(u_ct(i)%val(2,1:turbulenceequations+passivescalar))-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta2_2d
  subroutine sol_integ_dg(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real, dimension(1:nof_variables)::solution_integ2
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      call solution_integ(i, solution_integ2)
      u_c(i)%val(1, 1:nof_variables) = solution_integ2(1:nof_variables)
    end do

  end subroutine sol_integ_dg
  subroutine sol_integ_dgx(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, iconsidered
    real, dimension(1:nof_variables)::solution_integ2
    real, dimension(1:nof_variables)::solution_integ_weak
    real, dimension(1:nof_variables)::solution_integ_strong
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call solution_integ(i, solution_integ2)
      u_c(i)%val(1, :) = solution_integ2(1:nof_variables)
      if (filtering .eq. 1) then
        call solution_integ_s(i, solution_integ_strong)
        call solution_integ_w(i, solution_integ_weak)
        u_cs(i)%val(1, :) = solution_integ_strong(1:nof_variables)
        u_cw(i)%val(1, :) = solution_integ_weak(1:nof_variables)
      end if
    end do

  end subroutine sol_integ_dgx
  subroutine sol_integ_dg_init(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real, dimension(1:nof_variables)::solution_integ2
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      call solution_integ(i, solution_integ2)
      u_c(i)%val(1, :) = solution_integ2(1:nof_variables)
      u_e(i)%val(1, :) = u_c(i)%val(1, :)
    end do
  end subroutine sol_integ_dg_init
  subroutine runge_kutta1_2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, kx
    real::avrgs, oovolume, to4, oo4, to3, oo3
    real, dimension(nof_variables)::tempsol
    kmaxe = xmpielrank(n)
    call call_flux_subroutines_2d
    do i = 1, kmaxe
    if (dg .eq. 1) then
    if ((u_c(i)%valdg(1, 1, 1) .ne. u_c(i)%valdg(1, 1, 1))) then
      if (n .eq. 0) print *, 'stopping because nans1'
      stop ! stop if nans
    end if
    end if
    end do
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(1,1:nof_variables,:) - dt* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(1,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta1_2d
  subroutine runge_kutta3_2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        u_c(i)%valdg(2, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
          u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(2,1:nof_variables,:) - dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        u_c(i)%valdg(3, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
          u_c(i)%valdg(1,1:nof_variables,:)=to4*u_c(i)%valdg(2,1:nof_variables,:) + oo4*u_c(i)%valdg(3,1:nof_variables,:) - oo4*dt* transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    u_c(i)%val(1, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
      end if
    end do

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(to4*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo4*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((oo4))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:)=oo3*u_c(i)%valdg(2,1:nof_variables,:) + to3*u_c(i)%valdg(1,1:nof_variables,:) - to3*dt*transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(1, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=((oo3)*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+((to3)*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(((to3))*&
                                                           ((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta3_2d
  subroutine runge_kutta3_2d_mood(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, inds
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    if (mood .eq. 1) then
      inds = 4
    else
      inds = 1
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      if (mood .eq. 1) then
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      end if
      u_c(i)%val(inds, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (mood .eq. 1) then
      call mood_operator_2(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
          u_c(i)%val(inds, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*(rhs(i)%val(1:nof_variables)*oovolume))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
 u_c(i)%val(inds, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(to4*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(oo4*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((oo4))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if

    if (mood .eq. 1) then

      call mood_operator_2(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
 u_c(i)%val(inds, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
    u_c(i)%val(1, 1:nof_variables) = (to4*u_c(i)%val(2, 1:nof_variables)) + (oo4*u_c(i)%val(3, 1:nof_variables)) - (((oo4))*((dt)* &
                                                                                        ((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
        call vortexcalc2d(n)
      end select
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (mood .eq. 1) then
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      end if
   u_c(i)%val(inds, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
    end do
    if (mood .eq. 1) then
      call mood_operator_2(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
   u_c(i)%val(inds, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 2
        end if
      end do
      call mood_operator_1(n)
      do i = 1, kmaxe
        if (ielem(n, i)%recalc .eq. 1) then
          oovolume = 1.0d0/ielem(n, i)%totvolume
      u_c(i)%val(1, 1:nof_variables) = ((oo3)*u_c(i)%val(2, 1:nof_variables)) + ((to3)*u_c(i)%val(1, 1:nof_variables)) - (((to3))* &
                                                                                  ((dt)*((rhs(i)%val(1:nof_variables))*(oovolume))))
          ielem(n, i)%mood_o = 1
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(4, 1:nof_variables)
        end if
      end do
    end if
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
u_ct(i)%val(1,1:turbulenceequations+passivescalar)=((oo3)*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+((to3)*u_ct(i)%val(1,1:turbulenceequations+passivescalar))-(((to3))*&
                                                           ((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta3_2d_mood
  subroutine runge_kutta4(n)
!> @brief
!> ssp runge kutta 4th-order scheme
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe
    real::avrgs, oovolume, to4, oo4, to3, oo3
    real::dumpracein, dumpraceout, flops_count
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(2, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(2,1:nof_variables,:) - dt * 0.391752226571890 * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))!*oovolume
      else
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
     u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*0.391752226571890*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*0.391752226571890*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    if (statistics .eq. 1) then
      !$omp barrier
      !$omp master
      pr_t8 = mpi_wtime()
      prace_t7 = pr_t8 - pr_t7
      dumpracein = prace_t1
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t1 = dumpraceout
      dumpracein = prace_t2
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t2 = dumpraceout
      dumpracein = prace_t3
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t3 = dumpraceout
      dumpracein = prace_t4
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t4 = dumpraceout
      dumpracein = prace_t5
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t5 = dumpraceout
      dumpracein = prace_t6
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t6 = dumpraceout
      dumpracein = prace_t7
      call mpi_allreduce(dumpracein, dumpraceout, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      prace_t7 = dumpraceout
      prace_tx1 = prace_t2 + prace_t4
      prace_tx2 = prace_t1 + prace_t3 + prace_t5 + prace_t6 + prace_t7
      prace_tx3 = prace_tx1 + prace_tx2
      if (n .eq. 0) then
        open (133, file=statfile, form='formatted', status='old', action='write', position='append')
    write(133,'(i6,1x,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4)')it,prace_tx3,prace_tx1,prace_tx2,prace_t1,prace_t2,prace_t3,prace_t4,prace_t5,prace_t6,prace_t7
        close (133)
      end if
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(3, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:) = 0.444370493651235 * u_c(i)%valdg(2,1:nof_variables,:) + 0.555629506348765 * u_c(i)%valdg(3,1:nof_variables,:) - 0.368410593050371 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1,1:nof_variables)=(0.444370493651235 * u_c(i)%val(2,1:nof_variables)) + (0.555629506348765 * u_c(i)%val(3,1:nof_variables)) - 0.368410593050371 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.444370493651235*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.555629506348765*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((0.368410593050371))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(4, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:) = 0.620101851488403 * u_c(i)%valdg(2,1:nof_variables,:) + 0.379898148511597 * u_c(i)%valdg(4,1:nof_variables,:) - 0.251891774271694 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(4, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1,1:nof_variables) = 0.620101851488403 * u_c(i)%val(2,1:nof_variables) + 0.379898148511597 * u_c(i)%val(4,1:nof_variables) - 0.251891774271694 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(4, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.620101851488403*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.379898148511597*u_ct(i)%val(4,1:turbulenceequations+passivescalar))-(((0.251891774271694))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(5, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(6, 1:nof_variables, :) = -dt*transpose(matmul(m_1(i)%val(:, :), rhs(i)%valdg(:, 1:nof_variables)))! * oovolume
        u_c(i)%valdg(1,1:nof_variables,:) = 0.178079954393132 * u_c(i)%valdg(2,1:nof_variables,:) + 0.821920045606868 * u_c(i)%valdg(5,1:nof_variables,:) - 0.544974750228521 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        u_c(i)%val(5, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(6, 1:nof_variables) = -dt*rhs(i)%val(1:nof_variables)*oovolume
        u_c(i)%val(1,1:nof_variables) = 0.178079954393132 * u_c(i)%val(2,1:nof_variables) + 0.821920045606868 * u_c(i)%val(5,1:nof_variables) - 0.544974750228521 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(5, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
   u_ct(i)%val(6, 1:turbulenceequations + passivescalar) = -((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume)))
   u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.178079954393132*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.821920045606868*u_ct(i)%val(5,1:turbulenceequations+passivescalar))-(((0.544974750228521))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_3d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:) = (0.00683325884039 * u_c(i)%valdg(2,1:nof_variables,:)) + (0.517231671970585 * u_c(i)%valdg(4,1:nof_variables,:)) + (0.12759831133288 * u_c(i)%valdg(5,1:nof_variables,:)) + (0.34833675773694 * u_c(i)%valdg(1,1:nof_variables,:)) + (0.08460416338212 * u_c(i)%valdg(6,1:nof_variables,:)) - 0.22600748319395 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(1,1:nof_variables) = (0.00683325884039 * u_c(i)%val(2,1:nof_variables)) + (0.517231671970585 * u_c(i)%val(4,1:nof_variables)) +  (0.12759831133288 * u_c(i)%val(5,1:nof_variables)) + (0.34833675773694 * u_c(i)%val(1,1:nof_variables)) + (0.08460416338212 * u_c(i)%val(6,1:nof_variables)) - (0.22600748319395 * dt * rhs(i)%val(1:nof_variables) * oovolume)
      end if
    end do
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
   u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.00683325884039*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.517231671970585*u_ct(i)%val(4,1:turbulenceequations+passivescalar))+&
                                (0.12759831133288*u_ct(i)%val(5,1:turbulenceequations+passivescalar))+(0.34833675773694*u_ct(i)%val(1,1:turbulenceequations+passivescalar))+&
                                (0.08460416338212*u_ct(i)%val(6,1:turbulenceequations+passivescalar))-(0.22600748319395*(dt)*((rhst(i)%val(1:turbulenceequations+passivescalar))*(oovolume)))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine runge_kutta4
  subroutine call_flux_subroutines_3d
    implicit none
    real::dumpracein, dumpraceout
    integer::kmaxe, i
    kmaxe = xmpielrank(n)
    if (statistics .eq. 1) then
      pr_t1 = mpi_wtime()
    end if
    if (dg .eq. 1) then
      call sol_integ_dg(n) ! calculates cell average of dg solution for fv
    end if
    if (dg .eq. 1) then
      do i = 1, kmaxe
        ielem(n, i)%filtered = 0
      end do
    end if
    if ((dg .eq. 1) .and. (filtering .eq. 1)) then
      call sol_integ_dgx(n)
      call apply_filter_dg(n)
    end if
    if (statistics .eq. 1) then
      !$omp barrier
      !$omp master
      pr_t2 = mpi_wtime()
      prace_t1 = pr_t2 - pr_t1
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
    else
      call exchange_higher(n)
    end if
    if (statistics .eq. 1) then
      pr_t3 = mpi_wtime()
      prace_t2 = pr_t3 - pr_t2
    end if
    if (dg .eq. 1) then
      call reconstruct_dg(n) ! extrapolates solution to faces
      if (multispecies .eq. 1) then
        if (bound_lim .eq. 1) then
          call vfbp_limiter
        end if
      end if
      call trouble_indicator1 ! checks for troubled cells
    end if
    call arbitrary_order(n)
    if (dg .eq. 1) then
      call trouble_indicator2 ! changes dg to fv
    end if
    call exhboundhigher(n)
    if (dg .eq. 1) then
      call exhboundhigher_dg(n)
      if (itestcase .eq. 4) then
        if (br2_yn .eq. 2) then
          call reconstruct_br2_dg
          call exhboundhigher_dg2(n)
        end if
        if (br2_yn .eq. 0) then
          call viscous_dg_ggs(n)
        end if
      end if
    end if
    if (statistics .eq. 1) then
      !$omp barrier
      !$omp master
      pr_t4 = mpi_wtime()
      prace_t3 = pr_t4 - pr_t3
    end if

    if (statistics .eq. 1) then
      pr_t5 = mpi_wtime()
      prace_t4 = pr_t5 - pr_t4
    end if
    if (adda .eq. 1) then
      if (rungekutta .eq. 11) then
        if (iscoun .eq. 1) then
          call fix_dissipation(n)
          call exchange_adda_diss(n)
          call fix_dissipation2(n)
        end if
      else
        call fix_dissipation(n)
        call exchange_adda_diss(n)
        call fix_dissipation2(n)
      end if
    end if
    if (statistics .eq. 1) then
      pr_t6 = mpi_wtime()
      prace_t5 = pr_t6 - pr_t5
    end if
    select case (itestcase)
    case (1, 2)
      call calculate_fluxeshi(n)
    case (3)
      call calculate_fluxeshi_convective(n)
      if ((source_active .eq. 1)) then
        call sources_computation_rot(n)
      end if
    case (4)
      call calculate_fluxeshi_convective(n)
      call calculate_fluxeshi_diffusive(n)
      if ((source_active .eq. 1)) then
        call sources_computation_rot(n)
      end if
      if (turbulence .eq. 1) then
        call sources_computation(n)
      end if
      call vortexcalc(n)
    end select
    if (initcond .eq. 95) then
      call enstrophy_calc(n)
    end if
    if (statistics .eq. 1) then
      pr_t7 = mpi_wtime()
      prace_t6 = pr_t7 - pr_t6
    end if
  end subroutine call_flux_subroutines_3d
  subroutine call_flux_subroutines_2d
    implicit none
    integer::i, iconsidered
    if (dg .eq. 1) then
      call sol_integ_dg(n)
    end if
    if (fastest .eq. 1) then
      call exchange_lower(n)
    else
      call exchange_higher(n)
    end if
    if (dg .eq. 1) then
      call reconstruct_dg(n)
      call trouble_indicator1
    end if
    call arbitrary_order(n)
    if (dg .eq. 1) then
      call trouble_indicator2
    end if
    if (bound_lim .eq. 1) then
      call vfbp_limiter
    end if
    call exhboundhigher(n)
    if (dg .eq. 1) then
      call exhboundhigher_dg(n)
      if (itestcase .eq. 4) then
        if (br2_yn .eq. 2) then
          call reconstruct_br2_dg
          call exhboundhigher_dg2(n)
        end if
        if (br2_yn .eq. 0) then
          call viscous_dg_ggs(n)
        end if
      end if
    end if
    select case (itestcase)
    case (1, 2)
      call calculate_fluxeshi2d(n)
    case (3)
      call calculate_fluxeshi_convective2d(n)
    case (4)
      call calculate_fluxeshi_convective2d(n)
      call calculate_fluxeshi_diffusive2d(n)
      if (turbulence .eq. 1) then
        call sources_computation2d(n)
      end if
    end select
!     if (code_profile.eq.30)then
!         if (probei(n,1).gt.0) then
!           do i=1,xmpielrank(n)
!             if (ielem(n,i)%ihexgl.eq.probei(n,1))then
!             iconsidered=i
!             end if
!           end do
!         call vertex_neighbours_values(n)
!         end if
!     end if
  end subroutine call_flux_subroutines_2d
  subroutine runge_kutta4_2d(n)
!> @brief
!> ssp runge kutta 4th-order scheme in 2d
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, rk_stage
    real::avrgs, oovolume, to4, oo4, to3, oo3
    kmaxe = xmpielrank(n)
    to4 = 3.0d0/4.0d0
    oo4 = 1.0d0/4.0d0
    to3 = 2.0d0/3.0d0
    oo3 = 1.0d0/3.0d0
    rk_stage = 0
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(2, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:)=u_c(i)%valdg(2,1:nof_variables,:) - dt * 0.391752226571890 * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))!*oovolume
      else
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
     u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables) - (dt*0.391752226571890*(rhs(i)%val(1:nof_variables)*oovolume))
      end if
    end do
    rk_stage = rk_stage + 1
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(2, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=u_ct(i)%val(2,1:turbulenceequations+passivescalar)-(dt*0.391752226571890*(rhst(i)%val(1:turbulenceequations+passivescalar)*oovolume))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(3, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:) = 0.444370493651235 * u_c(i)%valdg(2,1:nof_variables,:) + 0.555629506348765 * u_c(i)%valdg(3,1:nof_variables,:) - 0.368410593050371 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1,1:nof_variables)=(0.444370493651235 * u_c(i)%val(2,1:nof_variables)) + (0.555629506348765 * u_c(i)%val(3,1:nof_variables)) - 0.368410593050371 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    rk_stage = rk_stage + 1
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(3, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.444370493651235*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.555629506348765*u_ct(i)%val(3,1:turbulenceequations+passivescalar))-(((0.368410593050371))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(4, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(1,1:nof_variables,:) = 0.620101851488403 * u_c(i)%valdg(2,1:nof_variables,:) + 0.379898148511597 * u_c(i)%valdg(4,1:nof_variables,:) - 0.251891774271694 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(4, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(1,1:nof_variables) = 0.620101851488403 * u_c(i)%val(2,1:nof_variables) + 0.379898148511597 * u_c(i)%val(4,1:nof_variables) - 0.251891774271694 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    rk_stage = rk_stage + 1
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(4, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.620101851488403*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.379898148511597*u_ct(i)%val(4,1:turbulenceequations+passivescalar))-(((0.251891774271694))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_2d
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(5, 1:nof_variables, :) = u_c(i)%valdg(1, 1:nof_variables, :)
        u_c(i)%valdg(6, 1:nof_variables, :) = -dt*transpose(matmul(m_1(i)%val(:, :), rhs(i)%valdg(:, 1:nof_variables)))! * oovolume
        u_c(i)%valdg(1,1:nof_variables,:) = 0.178079954393132 * u_c(i)%valdg(2,1:nof_variables,:) + 0.821920045606868 * u_c(i)%valdg(5,1:nof_variables,:) - 0.544974750228521 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))
      else
        u_c(i)%val(5, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(6, 1:nof_variables) = -dt*rhs(i)%val(1:nof_variables)*oovolume
        u_c(i)%val(1,1:nof_variables) = 0.178079954393132 * u_c(i)%val(2,1:nof_variables) + 0.821920045606868 * u_c(i)%val(5,1:nof_variables) - 0.544974750228521 * dt * rhs(i)%val(1:nof_variables) * oovolume
      end if
    end do
    rk_stage = rk_stage + 1
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
        u_ct(i)%val(5, 1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
   u_ct(i)%val(6, 1:turbulenceequations + passivescalar) = -((dt)*((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume)))
  u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.178079954393132*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.821920045606868*u_ct(i)%val(5,1:turbulenceequations+passivescalar))-(((0.544974750228521))*((dt)*&
                                                                 ((rhst(i)%val(1:turbulenceequations + passivescalar))*(oovolume))))
      end do
    end if
    call call_flux_subroutines_2d
    if (fastest .ne. 1) then
      if (itestcase .eq. 4) then
        call vortexcalc2d(n)
      end if
    end if
    do i = 1, kmaxe
      oovolume = 1.0d0/ielem(n, i)%totvolume
      if (dg .eq. 1) then
        u_c(i)%valdg(1,1:nof_variables,:) = (0.00683325884039 * u_c(i)%valdg(2,1:nof_variables,:)) + (0.517231671970585 * u_c(i)%valdg(4,1:nof_variables,:)) + (0.12759831133288 * u_c(i)%valdg(5,1:nof_variables,:)) + (0.34833675773694 * u_c(i)%valdg(1,1:nof_variables,:)) + (0.08460416338212 * u_c(i)%valdg(6,1:nof_variables,:)) - 0.22600748319395 * dt * transpose(matmul(m_1(i)%val(:,:), rhs(i)%valdg(:,1:nof_variables)))! * oovolume
      else
        u_c(i)%val(1,1:nof_variables) = (0.00683325884039 * u_c(i)%val(2,1:nof_variables)) + (0.517231671970585 * u_c(i)%val(4,1:nof_variables)) +  (0.12759831133288 * u_c(i)%val(5,1:nof_variables)) + (0.34833675773694 * u_c(i)%val(1,1:nof_variables)) + (0.08460416338212 * u_c(i)%val(6,1:nof_variables)) - (0.22600748319395 * dt * rhs(i)%val(1:nof_variables) * oovolume)
      end if
    end do
    rk_stage = rk_stage + 1
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, kmaxe
        oovolume = 1.0d0/ielem(n, i)%totvolume
   u_ct(i)%val(1,1:turbulenceequations+passivescalar)=(0.00683325884039*u_ct(i)%val(2,1:turbulenceequations+passivescalar))+(0.517231671970585*u_ct(i)%val(4,1:turbulenceequations+passivescalar))+&
                                (0.12759831133288*u_ct(i)%val(5,1:turbulenceequations+passivescalar))+(0.34833675773694*u_ct(i)%val(1,1:turbulenceequations+passivescalar))+&
                                (0.08460416338212*u_ct(i)%val(6,1:turbulenceequations+passivescalar))-(0.22600748319395*(dt)*((rhst(i)%val(1:turbulenceequations+passivescalar))*(oovolume)))
      end do
    end if
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
    if (dg .eq. 1) then
    if (all(u_c(1)%valdg(1, :, :) .ne. u_c(1)%valdg(1, :, :))) then
      if (n .eq. 0) print *, 'stopping because nans'
      stop ! stop if nans
    end if
    end if

  end subroutine runge_kutta4_2d
  subroutine implicit_times(n)
!> @brief
!> implicit approximately factored time stepping scheme
    implicit none
    integer::i, k, kmaxe, kill_nan
    integer, intent(in)::n
    real::verysmall
    verysmall = tolsmall
    kmaxe = xmpielrank(n)
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
        if ((source_active .eq. 1)) then
          call sources_computation_rot(n)
        end if
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if ((source_active .eq. 1)) then
          call sources_computation_rot(n)
        end if
        call vortexcalc(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi(n)
      case (3)
        call calculate_fluxeshi_convective(n)
        if ((source_active .eq. 1)) then
          call sources_computation_rot(n)
        end if
      case (4)
        call calculate_fluxeshi_convective(n)
        call calculate_fluxeshi_diffusive(n)
        if ((source_active .eq. 1)) then
          call sources_computation_rot(n)
        end if
        call vortexcalc(n)
        if (turbulence .eq. 1) then
          call sources_computation(n)
        end if
      end select
    end if
    if (relax .eq. 3) then
      call relaxation_lumfree(n)
    else
      if (lowmemory .eq. 0) then
        call relaxation(n)
      else
        call relaxation_lm(n)
      end if
    end if
    kill_nan = 0
    do i = 1, kmaxe
    if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)).or.(impdu(i,5).ne.impdu(i,5)))then
        write (600 + n, *) "nan present", ielem(n, i)%ihexgl, ielem(n, i)%ishape, ielem(n, i)%xxc, ielem(n, i)%yyc, ielem(n, i)%zzc
        write (600 + n, *) ielem(n, i)%dih(:)
        write (500 + n, '(3es14.6)') ielem(n, i)%xxc, ielem(n, i)%yyc, ielem(n, i)%zzc
        if (mrf .eq. 1) then
          write (700 + n, *) 'srf -diagonal', ilocal_recon3(i)%mrf, i
          write (700 + n, '(3es14.6)') ielem(n, i)%xxc, ielem(n, i)%yyc, ielem(n, i)%zzc
          write (700 + n, *) impdu(i, 1), impdu(i, 2), impdu(i, 3), impdu(i, 4), impdu(i, 5)
        end if
        kill_nan = 1
      end if
      u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
    end do
!$omp end do
    if (kill_nan .eq. 1) then
      stop
    end if
    if ((passivescalar .gt. 0) .or. (turbulence .gt. 0)) then
      do i = 1, kmaxe
      do k = 1, turbulenceequations + passivescalar
      if (u_ct(i)%val(1, k) + impdu(i, 5 + k) .ge. zero) then
        u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + 0.4*impdu(i, 5 + k)
      end if
      end do
      end do
    end if
  end subroutine implicit_times
  subroutine implicit_times_2d(n)
!> @brief
!> implicit approximately factored time stepping scheme 2d
    implicit none
    integer::i, k, kmaxe, kill_nan
    integer, intent(in)::n
    real::verysmall
    verysmall = tolsmall
    kmaxe = xmpielrank(n)
    if (fastest .eq. 1) then
      call exchange_lower(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        ! call vortexcalc2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    else
      call exchange_higher(n)
      call arbitrary_order(n)
      call exhboundhigher(n)
      select case (itestcase)
      case (1, 2)
        call calculate_fluxeshi2d(n)
      case (3)
        call calculate_fluxeshi_convective2d(n)
      case (4)
        call calculate_fluxeshi_convective2d(n)
        call calculate_fluxeshi_diffusive2d(n)
        !call vortexcalc2d(n)
        if (turbulence .eq. 1) then
          call sources_computation2d(n)
        end if
      end select
    end if
    if (relax .eq. 3) then
      call relaxation_lumfree(n)
    else
      if (lowmemory .eq. 0) then
        call relaxation2d(n)
      else
        call relaxation_lm2d(n)
      end if
    end if
    kill_nan = 0
!$omp do
    do i = 1, kmaxe
      if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))then
        write (600 + n, *) "nan present", ielem(n, i)%ihexgl, ielem(n, i)%ishape, ielem(n, i)%xxc, ielem(n, i)%yyc
        write (600 + n, *) ielem(n, i)%dih(:)
        kill_nan = 1
      end if
      u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
    end do
    if (kill_nan .eq. 1) then
      stop
    end if
    if ((passivescalar .gt. 0) .or. (turbulence .gt. 0)) then
      do i = 1, kmaxe
      do k = 1, turbulenceequations + passivescalar
      if (ispal .eq. 1) then
      if (u_ct(i)%val(1, k) + impdu(i, 4 + k) .ge. zero) then
        u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + 0.4*impdu(i, 4 + k)
      end if
      else
      u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + 0.4*impdu(i, 4 + k)
      end if
      end do
      end do
    end if
  end subroutine implicit_times_2d
  subroutine dual_time(n)
    implicit none
    integer::i, k, kmaxe, jj, kill_nan
    integer, intent(in)::n
    real::verysmall
    real::firsti, resmaxi, rsumfacei, suml2ri, dummy3i, inner_tol
    verysmall = tolsmall
    inner_tol = reslimit
    kmaxe = xmpielrank(n)
    if (it .eq. restart) then
      do i = 1, kmaxe
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(3, :) = u_ct(i)%val(1, :)
          u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        end if
      end do
    end if
    firsti = 0.0d0
    do jj = 1, upperlimit
      rsumfacei = zero; allresdt = zero; dummy3i = zero; 
      if (jj .eq. 1) then
        iscoun = 1
      else
        iscoun = 2
      end if
      call call_flux_subroutines_3d
      if (relax .eq. 3) then
        call relaxation_lumfree(n)
      else
        if (lowmemory .eq. 0) then
          call relaxation(n)
        else
          call relaxation_lm(n)
        end if
      end if
      kill_nan = 0
      do i = 1, kmaxe
        rsumfacei = sqrt(((impdu(i, 1))**2) + ((impdu(i, 2))**2) + ((impdu(i, 3))**2) + ((impdu(i, 4))**2) + ((impdu(i, 5))**2))
        allresdt = allresdt + (rsumfacei*ielem(n, i)%totvolume)
      if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)).or.(impdu(i,5).ne.impdu(i,5)))then
          kill_nan = 1
        end if
      end do
      if (kill_nan .eq. 1) then
        stop
      end if
      dummy3i = zero
      call mpi_allreduce(allresdt, dummy3i, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      allresdt = dummy3i/totalvolume
      if (allresdt .gt. firsti) then
        firsti = allresdt
      end if
      allresdt = allresdt/firsti
      if (n .eq. 0) then
        open (77, file='res1.txt', form='formatted', action='write', position='append')
        write (77, *) allresdt, jj, it
        close (77)
      end if
      if ((allresdt .le. inner_tol) .or. (jj .eq. upperlimit)) then
        !$omp do
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            if (u_ct(i)%val(1, k) + impdu(i, nof_variables + k) .ge. zero) then
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end if
          end do
          end if
        end do
        exit
      else
        !$omp do
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            if (u_ct(i)%val(1, k) + impdu(i, nof_variables + k) .ge. zero) then
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end if
          end do
          end if
        end do
      end if
    end do
    do i = 1, kmaxe
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables)
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      !u_c(i)%val(1,1:nof_variables)=2.0*u_c(i)%val(2,1:nof_variables)-u_c(i)%val(3,1:nof_variables)

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        u_ct(i)%val(3, :) = u_ct(i)%val(2, :)
        u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        !u_ct(i)%val(1,:)=2.0*u_ct(i)%val(2,:)-u_ct(i)%val(3,:)
      end if
    end do

    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine dual_time

  subroutine dual_time_ex(n)
    implicit none
    integer::i, k, kmaxe, jj
    integer, intent(in)::n
    real::verysmall
    real::firsti, resmaxi, rsumfacei, suml2ri, dummy3i, inner_tol
    verysmall = tolsmall
    inner_tol = reslimit
    kmaxe = xmpielrank(n)
    if (it .eq. restart) then
      do i = 1, kmaxe
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(3, :) = u_ct(i)%val(1, :)
          u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        end if
      end do
    end if
    firsti = 0.0d0
    do jj = 1, upperlimit
      rsumfacei = zero; allresdt = zero; dummy3i = zero; 
      if (fastest .eq. 1) then
        call exchange_lower(n)
        call arbitrary_order(n)
        call exhboundhigher(n)
        select case (itestcase)
        case (1, 2)
          call calculate_fluxeshi(n)
        case (3)
          call calculate_fluxeshi_convective(n)
        case (4)
          call calculate_fluxeshi_convective(n)
          call calculate_fluxeshi_diffusive(n)
          call vortexcalc(n)
          if (turbulence .eq. 1) then
            call sources_computation(n)
          end if
        end select
      else
        call exchange_higher(n)
        call arbitrary_order(n)
        call exhboundhigher(n)
        select case (itestcase)
        case (1, 2)
          call calculate_fluxeshi(n)
        case (3)
          call calculate_fluxeshi_convective(n)
        case (4)
          call calculate_fluxeshi_convective(n)
          call calculate_fluxeshi_diffusive(n)
          call vortexcalc(n)
          if (turbulence .eq. 1) then
            call sources_computation(n)
          end if
        end select
      end if
      call relaxation_ex(n)
      do i = 1, kmaxe
        rsumfacei = sqrt(((impdu(i, 1))**2) + ((impdu(i, 2))**2) + ((impdu(i, 3))**2) + ((impdu(i, 4))**2) + ((impdu(i, 5))**2))
        allresdt = allresdt + (rsumfacei*ielem(n, i)%totvolume)
      end do
      dummy3i = zero
      call mpi_allreduce(allresdt, dummy3i, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      allresdt = dummy3i/totalvolume
      if (allresdt .gt. firsti) then
        firsti = allresdt
      end if
      allresdt = allresdt/firsti
      if (n .eq. 0) then
        write (777, *) allresdt, jj, it
      end if
      if ((allresdt .le. inner_tol) .or. (jj .eq. upperlimit)) then
        !$omp do
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
          end do
          end if
        end do
        exit
      else
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do k = 1, turbulenceequations + passivescalar
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end do
          end if
        end do
      end if
    end do
    do i = 1, kmaxe
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables)
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      u_c(i)%val(1, 1:nof_variables) = 2.0*u_c(i)%val(2, 1:nof_variables) - u_c(i)%val(3, 1:nof_variables)

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        u_ct(i)%val(3, :) = u_ct(i)%val(2, :)
        u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        u_ct(i)%val(1, :) = 2.0*u_ct(i)%val(2, :) - u_ct(i)%val(3, :)
      end if
    end do
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine dual_time_ex
  subroutine dual_time_ex_2d(n)
    implicit none
    integer::i, k, kmaxe, jj
    integer, intent(in)::n
    real::verysmall
    real::firsti, resmaxi, rsumfacei, suml2ri, dummy3i, inner_tol
    verysmall = tolsmall
    inner_tol = reslimit
    kmaxe = xmpielrank(n)
    if (it .eq. restart) then
      do i = 1, kmaxe
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(3, :) = u_ct(i)%val(1, :)
          u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        end if
      end do
    end if
    firsti = 0.0d0
    do jj = 1, upperlimit
      rsumfacei = zero; allresdt = zero; dummy3i = zero; 
      if (fastest .eq. 1) then
        call exchange_lower(n)
        call arbitrary_order(n)
        call exhboundhigher(n)

        select case (itestcase)
        case (1, 2)
          call calculate_fluxeshi2d(n)
        case (3)
          call calculate_fluxeshi_convective2d(n)
        case (4)
          call calculate_fluxeshi_convective2d(n)
          call calculate_fluxeshi_diffusive2d(n)
          call vortexcalc2d(n)
          if (turbulence .eq. 1) then
            call sources_computation2d(n)
          end if
        end select
      else
        call exchange_higher(n)
        call arbitrary_order(n)
        call exhboundhigher(n)
        select case (itestcase)
        case (1, 2)
          call calculate_fluxeshi2d(n)
        case (3)
          call calculate_fluxeshi_convective2d(n)
        case (4)
          call calculate_fluxeshi_convective2d(n)
          call calculate_fluxeshi_diffusive2d(n)
          call vortexcalc2d(n)
          if (turbulence .eq. 1) then
            call sources_computation2d(n)
          end if
        end select
      end if
      call relaxation_ex(n)
      do i = 1, kmaxe
        rsumfacei = sqrt(((impdu(i, 1))**2) + ((impdu(i, 2))**2) + ((impdu(i, 3))**2) + ((impdu(i, 4))**2))
        allresdt = allresdt + (rsumfacei*ielem(n, i)%totvolume)
      end do
      dummy3i = zero
      call mpi_allreduce(allresdt, dummy3i, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      allresdt = dummy3i/totalvolume
      if (allresdt .gt. firsti) then
        firsti = allresdt
      end if
      allresdt = allresdt/firsti
      if (n .eq. 0) then
        write (777, *) allresdt, jj, it
      end if
      if ((allresdt .le. inner_tol) .or. (jj .eq. upperlimit)) then
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
          end do
          end if
        end do
        exit
      else
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do k = 1, turbulenceequations + passivescalar
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end do
          end if
        end do
      end if
    end do
    do i = 1, kmaxe
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables)
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      u_c(i)%val(1, 1:nof_variables) = 2.0*u_c(i)%val(2, 1:nof_variables) - u_c(i)%val(3, 1:nof_variables)
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        u_ct(i)%val(3, :) = u_ct(i)%val(2, :)
        u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        u_ct(i)%val(1, :) = 2.0*u_ct(i)%val(2, :) - u_ct(i)%val(3, :)
      end if
    end do
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine dual_time_ex_2d
  subroutine dual_time_2d(n)
    implicit none
    integer::i, k, kmaxe, nvar, jj, kill_nan
    integer, intent(in)::n
    real::verysmall
    real::firsti, resmaxi, rsumfacei, suml2ri, dummy3i, inner_tol
    verysmall = tolsmall
    inner_tol = reslimit
    kmaxe = xmpielrank(n)
    if (it .eq. restart) then
      do i = 1, kmaxe
        u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(3, :) = u_ct(i)%val(1, :)
          u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        end if
      end do
    end if
    firsti = 0.0d0
    do jj = 1, upperlimit
      rsumfacei = zero; allresdt = zero; dummy3i = zero; 
      if (jj .eq. 1) then
        iscoun = 1
      else
        iscoun = 2
      end if
      call call_flux_subroutines_2d
      if (relax .eq. 3) then
        call relaxation_lumfree(n)
      else
        if (lowmemory .eq. 0) then
          call relaxation2d(n)
        else
          call relaxation_lm2d(n)
        end if
      end if

      do i = 1, kmaxe
        rsumfacei = sqrt(((impdu(i, 1))**2) + ((impdu(i, 2))**2) + ((impdu(i, 3))**2) + ((impdu(i, 4))**2))
        allresdt = allresdt + (rsumfacei*ielem(n, i)%totvolume)
       if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))then
          kill_nan = 1
        end if
      end do
      if (kill_nan .eq. 1) then
        stop
      end if
      dummy3i = zero
      call mpi_allreduce(allresdt, dummy3i, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      allresdt = dummy3i/totalvolume
      if (allresdt .gt. firsti) then
        firsti = allresdt
      end if
      allresdt = allresdt/firsti
      if (n .eq. 0) then
        open (77, file='res1.txt', form='formatted', action='write', position='append')
        write (77, *) allresdt, jj, it
        close (77)
      end if
      if ((allresdt .le. inner_tol) .or. (jj .eq. upperlimit)) then
        !$omp do
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            if (u_ct(i)%val(1, k) + impdu(i, nof_variables + k) .ge. zero) then
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end if
          end do
          end if
        end do
        exit
      else
        do i = 1, kmaxe
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + impdu(i, 1:nof_variables)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
          do k = 1, turbulenceequations + passivescalar
            if (u_ct(i)%val(1, k) + impdu(i, nof_variables + k) .ge. zero) then
              u_ct(i)%val(1, k) = u_ct(i)%val(1, k) + impdu(i, nof_variables + k)
            end if
          end do
          end if
        end do
      end if
    end do
    do i = 1, kmaxe
      u_c(i)%val(3, 1:nof_variables) = u_c(i)%val(2, 1:nof_variables)
      u_c(i)%val(2, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      !u_c(i)%val(1,1:nof_variables)=(2.0*u_c(i)%val(2,1:nof_variables))-u_c(i)%val(3,1:nof_variables)
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        u_ct(i)%val(3, :) = u_ct(i)%val(2, :)
        u_ct(i)%val(2, :) = u_ct(i)%val(1, :)
        !u_ct(i)%val(1,:)=2.0*u_ct(i)%val(2,:)-u_ct(i)%val(3,:)
      end if
    end do
    if (averaging .eq. 1) then
      call averaging_t(n)
    end if
  end subroutine dual_time_2d
  subroutine relaxation_ex(n)
    implicit none
!> @brief
!> this subroutine solves the linear system for implicit time stepping either through matrix free lu-sgs low memory footprint
    integer, intent(in)::n
    integer::i, l, k, ii, sweeps, kmaxe, nvar, igoflux, icaseb, indt1, indt2, indt3, ijk
    real::dt1, dtau
    kmaxe = xmpielrank(n)
    impdu(:, :) = zero
    indt1 = nof_variables + 1
    indt2 = nof_variables + turbulenceequations + passivescalar
    indt3 = turbulenceequations + passivescalar
    do i = 1, kmaxe
      dt1 = ielem(n, i)%totvolume/dt
impdu(i,1:nof_variables)=dt1*(1.5d0*u_c(i)%val(1,1:nof_variables)-2.0d0*u_c(i)%val(2,1:nof_variables)+0.5d0*u_c(i)%val(3,1:nof_variables))+rhs(i)%val(1:nof_variables)
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
impdu(i,indt1:indt2)=dt1*(1.5d0*u_ct(i)%val(1,1:indt3)-2.0d0*u_ct(i)%val(2,1:indt3)+0.5d0*u_ct(i)%val(3,1:indt3))+rhst(i)%val(1:indt3)
      end if
    end do
    do i = 1, kmaxe
      dtau = (ielem(n, i)%dtl/ielem(n, i)%totvolume)*(1.0d0/(1.0d0 + 1.5d0*(ielem(n, i)%dtl/dt)))
      impdu(i, 1:nof_variables) = -impdu(i, 1:nof_variables)*dtau
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        impdu(i, indt1:indt2) = -impdu(i, indt1:indt2)*dtau
      end if
    end do
  end subroutine relaxation_ex
  subroutine averaging_t(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, nvar
    integer::ind1
    kmaxe = xmpielrank(n)
    if (dimensiona .eq. 3) then
    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if
    if (tz1 .gt. zero) then
      !$omp do
      do i = 1, kmaxe
        u_c(i)%val(ind1, :) = (((tz1 - dt)/(tz1))*u_c(i)%val(ind1, :)) + ((dt*u_c(i)%val(1, :))/tz1)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          do nvar = 1, turbulenceequations + passivescalar
      u_ct(i)%val(ind1, nvar) = (((tz1 - dt)/(tz1))*u_ct(i)%val(ind1, nvar)) + ((dt*u_ct(i)%val(1, nvar))/(tz1*u_c(i)%val(ind1, 1)))
          end do
        end if
        !u,v,w,uv,uw,wv,ps
    u_c(i)%rms(1)=sqrt(abs(((u_c(i)%rms(1)**2)*((tz1-dt)/(tz1)))+(((u_c(i)%val(1,2)/u_c(i)%val(1,1)-u_c(i)%val(ind1,2)/u_c(i)%val(ind1,1))**2)*dt/tz1)))
    u_c(i)%rms(2)=sqrt(abs(((u_c(i)%rms(2)**2)*((tz1-dt)/(tz1)))+(((u_c(i)%val(1,3)/u_c(i)%val(1,1)-u_c(i)%val(ind1,3)/u_c(i)%val(ind1,1))**2)*dt/tz1)))
    u_c(i)%rms(3)=sqrt(abs(((u_c(i)%rms(3)**2)*((tz1-dt)/(tz1)))+(((u_c(i)%val(1,4)/u_c(i)%val(1,1)-u_c(i)%val(ind1,4)/u_c(i)%val(ind1,1))**2)*dt/tz1)))
        u_c(i)%rms(4) = (((u_c(i)%rms(4))*((tz1 - dt)/(tz1))) + &
(((((u_c(i)%val(1,2)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,2)/u_c(i)%val(ind1,1)))*((u_c(i)%val(1,3)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,3)/u_c(i)%val(ind1,1)))))*dt/tz1))
        u_c(i)%rms(5) = (((u_c(i)%rms(5))*((tz1 - dt)/(tz1))) + &
(((((u_c(i)%val(1,2)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,2)/u_c(i)%val(ind1,1)))*((u_c(i)%val(1,4)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,4)/u_c(i)%val(ind1,1)))))*dt/tz1))
        u_c(i)%rms(6) = (((u_c(i)%rms(6))*((tz1 - dt)/(tz1))) + &
(((((u_c(i)%val(1,3)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,3)/u_c(i)%val(ind1,1)))*((u_c(i)%val(1,4)/u_c(i)%val(1,1))-(u_c(i)%val(ind1,4)/u_c(i)%val(ind1,1)))))*dt/tz1))
        if ((passivescalar .gt. 0)) then
          u_c(i)%rms(7) = sqrt(abs(((u_c(i)%rms(7)**2)*((tz1 - dt)/(tz1))) + (((u_ct(i)%val(1, turbulenceequations + 1) &
                                                                         - u_ct(i)%val(ind1, turbulenceequations + 1))**2)*dt/tz1)))
        end if
      end do
    else
      !$omp do
      do i = 1, kmaxe
        u_c(i)%val(ind1, :) = zero; u_c(i)%rms = zero
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(ind1, :) = zero
        end if
      end do
    end if
    else
    if (t .gt. 0.0) then
      !$omp do
      do i = 1, kmaxe
        u_c(i)%val(ind1, :) = (((tz1 - dt)/(tz1))*u_c(i)%val(ind1, :)) + ((dt*u_c(i)%val(1, :))/tz1)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          do nvar = 1, turbulenceequations + passivescalar
      u_ct(i)%val(ind1, nvar) = (((tz1 - dt)/(tz1))*u_ct(i)%val(ind1, nvar)) + ((dt*u_ct(i)%val(1, nvar))/(tz1*u_c(i)%val(ind1, 1)))
          end do
        end if
        !u,v,uv,ps
u_c(i)%rms(1)=sqrt(abs(((u_c(i)%rms(1)**2)*((tz1-dt)/(tz1)))+(((u_c(i)%val(1,2)-u_c(i)%val(ind1,2))**2)*dt/tz1)))/u_c(i)%val(ind1,1)
u_c(i)%rms(2)=sqrt(abs(((u_c(i)%rms(2)**2)*((tz1-dt)/(tz1)))+(((u_c(i)%val(1,3)-u_c(i)%val(ind1,3))**2)*dt/tz1)))/u_c(i)%val(ind1,1)

        u_c(i)%rms(3) = (((u_c(i)%rms(4))*((tz1 - dt)/(tz1))) + &
                 ((((u_c(i)%val(1, 2) - u_c(i)%val(ind1, 2))*(u_c(i)%val(1, 3) - u_c(i)%val(ind1, 3))))*dt/tz1))/u_c(i)%val(ind1, 1)

        if ((passivescalar .gt. 0)) then
          u_c(i)%rms(4) = sqrt(abs(((u_c(i)%rms(4)**2)*((tz1 - dt)/(tz1))) + (((u_ct(i)%val(1, turbulenceequations + 1) &
                                                     - u_ct(i)%val(ind1, turbulenceequations + 1))**2)*dt/tz1)))/u_c(i)%val(ind1, 1)
        end if
      end do
    else
      do i = 1, kmaxe
        u_c(i)%val(ind1, :) = zero
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          u_ct(i)%val(ind1, :) = zero
        end if
      end do
    end if
    end if
    if (outsurf .eq. 1) then
      call exchange_higher_av(n)
      call average_stresses(n)
    end if
  end subroutine averaging_t
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!subroutine called to advance solution by one time step size!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine time_marching(n)
    implicit none
    integer, intent(in)::n
    real, dimension(1:5)::dummyout, dummyin
    integer::i, kmaxe, ttime
    real::dtiv
real::cput1,cput2,cput3,cput4,cput5,cput6,cput8,timec3,timec1,timec4,timec8,totv1,totv2,dumetg1,dumetg2,totk,tzx1,tzx2,resolx,totens,totens1,totens2,totensx,totensx1,totensx2
    kill = 0
    t = res_time
    resolx = 0.01
    iscoun = 1
    kmaxe = xmpielrank(n)
    every_time = ((idnint(t/output_freq))*output_freq) + output_freq
    totv1 = 0.0
    if (initcond .eq. 95) then
      call checkpointv3(n)
    end if
    cput1 = cpux1(1)
    cput4 = cpux1(1)
    cput5 = cpux1(1)
    cput8 = cpux1(1)
    it = restart
    if (dg .eq. 1) call sol_integ_dg_init(n)
    if (tecplot .lt. 5) then
      call grid_write
      if (outsurf .eq. 1) then
        call surf_write
      end if
    end if
    if ((average_restart .eq. 0) .and. (averaging .eq. 1)) then
      tz1 = 0.0
    else
      tz1 = t
    end if
    call volume_solution_write
    if (outsurf .eq. 1) then
      call surface_solution_write
    end if
    if ((it .eq. 0) .and. (initcond .eq. 95)) then
      call exchange_higher(n)
      call arbitrary_order(n)
      call enstrophy_calc(n)
    end if
    do
      call calculate_cfl(n)
      if (rungekutta .ge. 5) call calculate_cfll(n)
      if (dg .eq. 1) then
        do i = 1, kmaxe
          ielem(n, i)%condition = 0
          ielem(n, i)%troubled = 0
        end do
      end if
      !$omp barrier
      !$omp master
      dummyout(1) = dt
      cput2 = mpi_wtime()
      timec8 = cput2 - cput8
      timec1 = cput2 - cput1
      dummyout(2) = timec1
      dummyin = 0.0
      timec3 = cput2 - cput4
      dummyout(3) = timec3
      timec4 = cput2 - cput5
      dummyout(4) = timec4
      dummyout(5) = timec8
      call mpi_allreduce(dummyout, dummyin, 5, mpi_double_precision, mpi_min, mpi_comm_world, ierror)
      dtiv = dummyin(1)
      dt = dummyin(1)
      timec1 = dummyin(2)
      timec3 = dummyin(3)
      timec4 = dummyin(4)
      timec8 = dummyin(5)
      if (n .eq. 0) then
        open (63, file='history.txt', form='formatted', status='old', action='write', position='append')
        write (63, *) dt, it, "time step size", t
        close (63)
      end if
      if (initcond .eq. 95) then
        totk = 0; totens = 0; totensx = 0.0d0
        do i = 1, xmpielrank(n)
          totk = totk + ielem(n, i)%totvolume*u_c(i)%val(1, 1)*(1.0/2.0)* &
    (((u_c(i)%val(1, 2)/u_c(i)%val(1, 1))**2) + ((u_c(i)%val(1, 3)/u_c(i)%val(1, 1))**2) + ((u_c(i)%val(1, 4)/u_c(i)%val(1, 1))**2))
          if (boundtype .eq. 1) then
            totens = totens + (ielem(n, i)%totvolume*u_c(i)%val(1, 1)*(1.0/2.0)* &
                               ielem(n, i)%vortex(2))
          else
            totens = totens + (ielem(n, i)%vortex(2))
            totensx = totensx + (ielem(n, i)%vortex(3))
          end if
        end do
        dumetg1 = totk
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totk = dumetg2
        dumetg1 = totens
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totens = dumetg2
        dumetg1 = totensx
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totensx = dumetg2
        if (n .eq. 0) then
          totv1 = totk/((2.0*pi)**3)
          totens1 = totens/(((2.0*pi)**3))
          totensx1 = totensx/((2.0*pi)**3)
          if (it .eq. 0) then
            taylor = totk
            taylor_ens = totens
            taylor_ensx = totensx
          end if
        end if
      end if
      if (rungekutta .ge. 11) then
        dt = timestep
        if (initcond .eq. 95) then
          dt = min(dt, out_time - t, every_time - t)
        else
          dt = min(dt, out_time - t, every_time - t)
        end if
      else
        if (initcond .eq. 95) then
          dt = min(dt, out_time - t, every_time - t)
        else
          dt = min(dt, out_time - t, every_time - t)
        end if
      end if
      if (dg .eq. 1) then
      if (filtering .gt. 0) then
        call filter(n)
      end if
      end if
      select case (rungekutta)
      case (1)
        call runge_kutta1(n)
      case (2)
        call runge_kutta2(n)
      case (3)
        if (mood .eq. 1) then
          call runge_kutta3_mood(n)
        else
          call runge_kutta3(n)
        end if
      case (4)
        call runge_kutta4(n)
      case (5)
        call runge_kutta5(n)
      case (10)
        call implicit_times(n)
      case (11)
        call dual_time(n)
      case (12)
        call dual_time_ex(n)
      end select
      if (dg .eq. 1) call sol_integ_dg(n)
      if (rungekutta .ge. 11) then
        t = t + (dt)
        tz1 = tz1 + (dt)
      else
        t = t + dt
        tz1 = tz1 + dt
      end if
      if (dg .eq. 1) then
        if (code_profile .ne. 102) then
        if (mod(it, 100) .eq. 0) then
          call troubled_history
        end if
        end if
        if (filtering .eq. 1) then
          call filtered_history
        end if
      end if
      if (initcond .eq. 95) then
        totk = 0; totens = 0.0; totensx = 0.0d0
        do i = 1, xmpielrank(n)
          totk = totk + ielem(n, i)%totvolume*u_c(i)%val(1, 1)*(1.0/2.0)* &
    (((u_c(i)%val(1, 2)/u_c(i)%val(1, 1))**2) + ((u_c(i)%val(1, 3)/u_c(i)%val(1, 1))**2) + ((u_c(i)%val(1, 4)/u_c(i)%val(1, 1))**2))
          if (boundtype .eq. 1) then
            totens = totens + (ielem(n, i)%totvolume*u_c(i)%val(1, 1)*(1.0/2.0)* &
                               ielem(n, i)%vortex(2))
          else
            totens = totens + (ielem(n, i)%vortex(2))
            totensx = totensx + (ielem(n, i)%vortex(3))
          end if
        end do
        dumetg1 = totk
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totk = dumetg2
        dumetg1 = totens
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totens = dumetg2
        dumetg1 = totensx
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totensx = dumetg2
        if (n .eq. 0) then
          totv2 = totk/((2.0*pi)**3)
          totens2 = totens/((2.0*pi)**3)
          totensx2 = totensx/((2.0*pi)**3)
          if (it .eq. 0) then
            taylor = totk
            taylor_ens = totens
            taylor_ensx = totensx
          end if
          if (it .eq. 0) then
            open (73, file='energy.dat', form='formatted', status='new', action='write', position='append')
          else
            open (73, file='energy.dat', form='formatted', status='old', action='write', position='append')
          end if
          if (dg .eq. 1) then
            write (73, '(e14.7,1x,e14.7,1x,e14.7)') t, totk/taylor, -(totv2 - totv1)/dt
          else
            if (boundtype .eq. 1) then
              write (73, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, totk/taylor, -(totv2 - totv1)/dt, totens/taylor_ens
            else
              write (73, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, totk/taylor, -(totv2 - totv1)/dt, totens, totensx
            end if
          end if
          close (73)
        end if
        call mpi_barrier(mpi_comm_world, ierror)
        if (adda .eq. 1) then
          totk = 0
          do i = 1, xmpielrank(n)
            totk = totk + ielem(n, i)%er
          end do
          dumetg1 = totk
          dumetg2 = 0.0
          call mpi_barrier(mpi_comm_world, ierror)
          call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
          if (n .eq. 0) then
            totk = dumetg2/imaxe
            if (it .eq. 0) then
              open (123, file='er.dat', form='formatted', status='new', action='write', position='append')
            else
              open (123, file='er.dat', form='formatted', status='old', action='write', position='append')
            end if
            write (123, *) t, totk
            close (123)
          end if
        end if
        call mpi_barrier(mpi_comm_world, ierror)
        end if
      end if
      if (mod(it, iforce) .eq. 0) then
      if (outsurf .eq. 1) then
        call forces
      end if
      end if
      if ((rungekutta .ge. 5) .and. (rungekutta .lt. 11)) then
      if (mod(it, residualfreq) .eq. 0) then
        call residual_compute
      end if
      end if
      if (nprobes .gt. 0) call probing
      if (timec1 .ge. ievery) then
        call volume_solution_write
        if (outsurf .eq. 1) then
          call surface_solution_write
        end if
        cput1 = mpi_wtime()
      end if
      if (initcond .eq. 95) then
      if (abs(t - ((idnint(t/output_freq))*output_freq)) .le. tolsmall) then
        call volume_solution_write
        if (outsurf .eq. 1) then
          call surface_solution_write
        end if
        if (initcond .eq. 95) then
          call checkpointv4(n)
        end if
        every_time = every_time + output_freq
      end if
      else
      if (code_profile .eq. -1) then
        if (abs(t - ((idnint(t/output_freq))*output_freq)) .le. tolsmall) then
          call volume_solution_write
          if (outsurf .eq. 1) then
            call surface_solution_write
          end if
          every_time = every_time + output_freq
        end if
      end if
      end if
      if (timec8 .ge. ieveryav) then
      if (averaging .eq. 1) then
        call volume_solution_write_av
        if (outsurf .eq. 1) then
          call surface_solution_write_av
        end if
      end if
      cput8 = mpi_wtime()
      end if
      if (timec4 .ge. ievery2) then
        call checkpointing
        if (averaging .eq. 1) then
          call checkpointing_av
        end if
        cput5 = mpi_wtime()
      end if
      it = it + 1
      if ((it .eq. ntmax) .or. (timec3 .ge. wallc) .or. (dtiv .gt. out_time)) then
        kill = 1
      end if
      if ((rungekutta .lt. 5) .or. (rungekutta .ge. 11)) then
      if ((t .ge. out_time) .or. (dtiv .gt. out_time)) then
        kill = 1
      end if
      end if
      if (kill .eq. 1) then
        call volume_solution_write
        if (outsurf .eq. 1) then
          call surface_solution_write
        end if
        call checkpointing
        if (averaging .eq. 1) then
          call volume_solution_write_av
          if (outsurf .eq. 1) then
            call surface_solution_write_av
          end if
          call checkpointing_av
        end if
      end if
      if (kill .eq. 1) then
      if (itestcase .le. 3) then
        call calculate_error(n)
      end if
      return
      end if
    end do
  end subroutine time_marching
  subroutine time_marching2(n)
    implicit none
    integer, intent(in)::n
    real, dimension(1:5)::dummyout, dummyin
    integer::i, kmaxe
    real::cput1, cput2, cput3, cput4, cput5, cput6, cput8, timec3, timec1, timec4, timec8, totv1, totv2, dumetg1, dumetg2, totk
    real::dtiv, flort
    kmaxe = xmpielrank(n)
    kill = 0
    t = res_time
    iscoun = 1
    every_time = ((idnint(t/output_freq))*output_freq) + output_freq
!$omp master
    cput1 = cpux1(1)
    cput4 = cpux1(1)
    cput5 = cpux1(1)
    cput8 = cpux1(1)
!$omp end master
!$omp barrier
    it = restart
    if (dg .eq. 1) call sol_integ_dg_init(n)
!$omp barrier
!$omp master
    if (tecplot .lt. 5) then
      call grid_write
    end if
    call volume_solution_write
    if (outsurf .eq. 1) then
      call surf_write
    end if
    if ((average_restart .eq. 0) .and. (averaging .eq. 1)) then
      tz1 = 0.0
    else
      tz1 = t
    end if
    do
      call calculate_cfl2d(n)
      if (rungekutta .ge. 5) call calculate_cfll2d(n)
      if (dg .eq. 1) then
        do i = 1, kmaxe
          ielem(n, i)%condition = 0
          ielem(n, i)%troubled = 0
        end do
      end if
      dummyout(1) = dt
      cput2 = mpi_wtime()
      timec8 = cput2 - cput8
      timec1 = cput2 - cput1
      dummyout(2) = timec1
      dummyin = 0.0d0
      timec3 = cput2 - cput4
      dummyout(3) = timec3
      timec4 = cput2 - cput5
      dummyout(4) = timec4
      dummyout(5) = timec8
      call mpi_allreduce(dummyout, dummyin, 5, mpi_double_precision, mpi_min, mpi_comm_world, ierror)
      dtiv = dummyin(1)
      dt = dummyin(1)
      timec1 = dummyin(2)
      timec3 = dummyin(3)
      timec4 = dummyin(4)
      timec8 = dummyin(5)
      if (n .eq. 0) then
        open (63, file='history.txt', form='formatted', status='old', action='write', position='append')
        write (63, *) dt, it, "time step size", t
        close (63)
      end if
      if (initcond .eq. 95) then
        totk = 0
        do i = 1, kmaxe
          totk = totk + ielem(n, i)%totvolume*(1.0/2.0)* &
                 (((u_c(i)%val(1, 2)/u_c(i)%val(1, 1))**2) + ((u_c(i)%val(1, 3)/u_c(i)%val(1, 1))**2))
        end do
        dumetg1 = totk
        dumetg2 = 0.0
        call mpi_barrier(mpi_comm_world, ierror)
        call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        totk = dumetg2
        if (n .eq. 0) then
          if (it .eq. 0) then
            open (73, file='energy.dat', form='formatted', status='new', action='write', position='append')
          else
            open (73, file='energy.dat', form='formatted', status='old', action='write', position='append')
          end if
          write (73, *) t, totk
          close (73)
        end if
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rungekutta .ge. 11) then
        dt = timestep
        dt = min(dt, out_time - t, every_time - t)
      else
        dt = min(dt, out_time - t, every_time - t)
      end if
      select case (rungekutta)
      case (1)
        call runge_kutta1_2d(n)
      case (2)
        call runge_kutta2_2d(n)
      case (3)
        if (mood .eq. 1) then
          call runge_kutta3_2d_mood(n)
        else
          call runge_kutta3_2d(n)
        end if
      case (4)
        call runge_kutta4_2d(n)
      case (5)
        call runge_kutta5_2d(n)
      case (10)
        call implicit_times_2d(n)
      case (11)
        call dual_time_2d(n)
      case (12)
        call dual_time_ex_2d(n)
      end select
      if (dg .eq. 1) call sol_integ_dg(n)
      if (rungekutta .ge. 11) then
        t = t + (dt)
        tz1 = tz1 + (dt)
      else
        t = t + dt
        tz1 = tz1 + dt
      end if
      if (dg .eq. 1) then
      if (code_profile .ne. 102) then
      if (mod(it, 100) .eq. 0) then
        call troubled_history
      end if
      end if
      end if
      if (mood .gt. 0) then
        call troubled_history
      end if
      if (mod(it, iforce) .eq. 0) then
        if (outsurf .eq. 1) then
          call forces
        end if
      end if
      if ((rungekutta .ge. 5) .and. (rungekutta .lt. 11)) then
        if (mod(it, residualfreq) .eq. 0) then
          call residual_compute
        end if
      end if
      if (nprobes .gt. 0) call probing2d
      if (timec1 .ge. ievery) then
        call volume_solution_write
        if (outsurf .eq. 1) then
          call surface_solution_write
        end if
        cput1 = mpi_wtime()
      end if
      if (timec8 .ge. ieveryav) then
        if (averaging .eq. 1) then
          call volume_solution_write_av
          if (outsurf .eq. 1) then
            call surface_solution_write_av
          end if
        end if
        cput8 = mpi_wtime()
      end if
      if (code_profile .eq. -1) then
        if (abs(t - ((idnint(t/output_freq))*output_freq)) .le. tolsmall) then
          call volume_solution_write
          if (outsurf .eq. 1) then
            call surface_solution_write
          end if
          every_time = every_time + output_freq
        end if
      end if
      it = it + 1
      if ((it .eq. ntmax) .or. (timec3 .ge. wallc) .or. (dtiv .gt. out_time)) then
        kill = 1
      end if
      if ((rungekutta .lt. 5) .or. (rungekutta .ge. 11)) then
      if ((t .ge. out_time) .or. (dtiv .gt. out_time)) then
        kill = 1
      end if
      end if
      if (kill .eq. 1) then
        call volume_solution_write
        if (outsurf .eq. 1) then
          call surface_solution_write
        end if
        call checkpointing
        if (averaging .eq. 1) then
          call volume_solution_write_av
          if (outsurf .eq. 1) then
            call surface_solution_write_av
          end if
          call checkpointing_av
        end if
      end if
      if (kill .eq. 1) then
      if (itestcase .le. 3) then
        call calculate_error(n)
      end if
      return
      end if
    end do
  end subroutine time_marching2
!---------------------------------------------------------------------------------------------!
end module advance
