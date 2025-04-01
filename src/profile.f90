module profile
  use declaration
  use library
  use basis

  implicit none

contains

  real function linear_init3d(n, pox, poy, poz)
    implicit none
    integer, intent(in)::n
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    if (initcond .eq. 0) then
      if (((pox(1) .ge. 0.25d0) .and. (pox(1) .le. 0.75d0)) .and. ((poz(1) .ge. 0.25d0) .and. (poz(1) .le. 0.75d0))) then
        linear_init3d = 1.0d0
      else
        linear_init3d = 1.0d0
      end if
    end if
    if (initcond .eq. 2) then
      linear_init3d = (sin((2.0d0*pi)*(pox(1))))* &
                      (sin((2.0d0*pi)*(poy(1))))*(sin((2.0d0*pi)*(poz(1))))
    end if
  end function linear_init3d

  real function linear_init2d(n, pox, poy, poz)
    implicit none
    integer, intent(in)::n
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real::aadx, aady, sumf, rd
    integer::ixg
    sumf = zero
    if (initcond .eq. 1) then
      if (((pox(1) .ge. 0.25d0) .and. (pox(1) .le. 0.75d0)) .and. ((poy(1) .ge. 0.25d0) .and. (poy(1) .le. 0.75d0))) then
        linear_init2d = 1.0d0
      else
        linear_init2d = 0.0d0
      end if
    end if

    if (initcond .eq. 3) then
      linear_init2d = 0.0d0
      if (sqrt(((pox(1) - 0.25d0)**2) + ((poy(1) - 0.5d0)**2)) .le. 0.15) then
        rd = (1.0d0/0.15d0)*sqrt(((pox(1) - 0.25d0)**2) + ((poy(1) - 0.5d0)**2))
        linear_init2d = 0.25d0*(1.0d0 + cos(pi*min(rd, 1.0d0)))
      end if
      if (sqrt(((pox(1) - 0.5d0)**2) + ((poy(1) - 0.25d0)**2)) .le. 0.15) then
        rd = (1.0d0/0.15d0)*sqrt(((pox(1) - 0.5d0)**2) + ((poy(1) - 0.25d0)**2))
        linear_init2d = 1.0d0 - rd
      end if
      if (sqrt(((pox(1) - 0.5d0)**2) + ((poy(1) - 0.75d0)**2)) .le. 0.15) then
        rd = (1.0d0/0.15d0)*sqrt(((pox(1) - 0.5d0)**2) + ((poy(1) - 0.75d0)**2))
        if ((abs(pox(1) - 0.5) .ge. 0.025d0) .or. (poy(1) .gt. 0.85)) then
          linear_init2d = 1.0d0
        else
          linear_init2d = 0.0d0
        end if
      end if
    end if

    if (initcond .eq. 5) then
      linear_init2d = 0.0d0
      if ((sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.5d0)**2)).gt.0.25).and.(sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.5d0)**2)).lt.0.35))then
        linear_init2d = 1.0
      else
        linear_init2d = 0.0
      end if
    end if

    if (initcond .eq. 2) then
      linear_init2d = (sin((2.0d0*pi)*(pox(1))))*(sin((2.0d0*pi)*(poy(1))))
    end if
  end function linear_init2d
  subroutine initialise_euler3d(n, veccos, pox, poy, poz)
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::veccos
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real, dimension(1:nof_species)::mp_r, mp_a, mp_ie
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx, vhx, amp, dvel
    integer::u_cond1, u_cond2, u_cond3, u_cond4
    veccos(:) = zero
    r1 = rres
    p1 = pres
    s1 = sqrt((gamma*p1)/(r1))
    u1 = uvel
    v1 = vvel
    w1 = wvel
    skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = r1*(skin1 + ie1)
    if (mrf .eq. 1) then
      veccos(1) = r1
      veccos(2) = r1*u1 + 1.0e-15
      veccos(3) = r1*v1 + 1.0e-15
      veccos(4) = r1*w1 + 1.0e-15
      veccos(5) = e1
    else
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
    end if
    if (turbulence .eq. 1) then
      if (turbulencemodel .eq. 1) then
        veccos(6) = visc*turbinit
      end if
      if (turbulencemodel .eq. 2) then
        if (zero_turb_init .eq. 0) then
          if (rframe .eq. 0) then
            veccos(6) = (1.5d0*i_turb_inlet*(ufreestream**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          else
            veccos(6) = (1.5d0*i_turb_inlet*(v_ref**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          end if
        end if
        if (zero_turb_init .eq. 1) then
          if (rframe .eq. 0) then
            veccos(6) = (1.5d0*i_turb_inlet*(ufreestream**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          else
            veccos(6) = (1.5d0*i_turb_inlet*(v_ref**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          end if
        end if
      end if
    end if
    if (passivescalar .gt. 0) then
      veccos(5 + turbulenceequations + 1:5 + turbulenceequations + passivescalar) = zero
    end if
    if (initcond .eq. 10000) then        !shock density interaction
      r1 = 0.5d0
      p1 = 0.4127
      u1 = 0.0
      v1 = 0.0
      w1 = 0.0
      skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
    end if

    if (initcond .eq. 95) then        !taylor green initial profile
      if (boundtype .eq. 1) then
        r1 = 1.0d0
        w1 = 0.0d0
        p1 = 100.0d0 + ((r1/16.0d0)*((cos(2.0d0*poz(1))) + 2.0d0)*((cos(2.0d0*pox(1))) + (cos(2.0d0*poy(1)))))
        u1 = sin(pox(1))*cos(poy(1))*cos(poz(1))
        v1 = -cos(pox(1))*sin(poy(1))*cos(poz(1))
      else
        w1 = 0.0d0
        p1 = (1.0d0/(gamma*1.25*1.25)) + ((1.0d0/16.0d0)*((cos(2.0d0*poz(1))) + 2.0d0)*((cos(2.0d0*pox(1))) + (cos(2.0d0*poy(1)))))
        r1 = (p1*(gamma*1.25*1.25))
        u1 = sin(pox(1))*cos(poy(1))*cos(poz(1))
        v1 = -cos(pox(1))*sin(poy(1))*cos(poz(1))
      end if
      skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
    end if

    if (initcond .eq. 101) then        !shock density interaction
      if (pox(1) .lt. -4.0d0) then
        r1 = 3.8571d0
        u1 = 2.6294d0
        v1 = zero
        w1 = zero
        p1 = 10.333d0
      else
        r1 = (1.0d0 + 0.2d0*sin(5.0d0*pox(1)))
        u1 = zero
        v1 = zero
        w1 = zero
        p1 = 1
      end if
      skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
    end if

    if (initcond .eq. 405) then
      if (pox(1) .lt. -0.1d0) then
        mp_r(1) = 0.166315789d0
        mp_r(2) = 1.658d0
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 114.49d0
        v1 = 0.0d0
        w1 = 0.0d0
        p1 = 159060.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(((pox(1) + 0.05d0)**2) + ((poy(1) - 0.0d0)**2) + ((poz(1) - 0.0d0)**2)) .le. 0.025d0) then
          mp_r(1) = 0.166315789d0
          mp_r(2) = 1.204d0
          mp_a(1) = 0.95d0
          mp_a(2) = 0.05d0
          u1 = 0.0d0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 0.166315789d0
          mp_r(2) = 1.204d0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      veccos(6) = mp_r(1)*mp_a(1)
      veccos(7) = mp_r(2)*mp_a(2)
      veccos(8) = mp_a(1)
    end if

    if (initcond .eq. 470) then
      if (pox(1) .le. 1.0) then   !post shock concidions
        mp_r(2) = 1.0d0             ! water density
        mp_r(1) = 1.0d0                 ! air density
        mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 1.0d0                 ! air volume fraction
        u1 = 0.0                            ! m/s
        v1 = 0.0
        w1 = 0.0
        p1 = 1.0                      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(poy(1)**2 + poz(1)**2) .ge. 1.5d0) then
          mp_r(2) = 1.0d0         ! water density
          mp_r(1) = 0.125d0                 ! air density
          mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 1.0d0                 ! air volume fraction
          u1 = 0.0                            ! m/s
          v1 = 0.0
          w1 = 0.0
          p1 = 0.1                      ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(2) = 1.0d0         ! water density
          mp_r(1) = 0.125d0                 ! air density
          mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 0.0d0                 ! air volume fraction
          u1 = 0.0                            ! m/s
          v1 = 0.0
          w1 = 0.0
          p1 = 0.1                      ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      veccos(6) = mp_r(1)*mp_a(1)
      veccos(7) = mp_r(2)*mp_a(2)
      veccos(8) = mp_a(1)
    end if

    if (initcond .eq. 157) then
      if (sqrt(((pox(1) - 200.0e-6)**2) + ((poy(1) - 150.0e-6)**2) + ((poz(1) - 150.0e-6)**2)) .le. 50.0e-6) then
        mp_r(1) = 1.225
        mp_r(2) = 1000.0
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.0
        v1 = 0.0d0
        w1 = 0.0
        p1 = 100000.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      else
        if (pox(1) .le. 100.0e-6) then
          p1 = 35e6
          mp_r(1) = 1.225
          mp_r(2) = 1000.0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 1647.0
          v1 = 0.0d0
          w1 = 0.0d0
        else
          mp_r(1) = 1.225
          mp_r(2) = 1000.0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 100000.0d0
        end if
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      veccos(6) = mp_r(1)*mp_a(1)
      veccos(7) = mp_r(2)*mp_a(2)
      veccos(8) = mp_a(1)
    end if

    if (initcond .eq. 411) then
      if (pox(1) .le. 0.0066d0) then
        mp_r(2) = 1323.65d0         ! water density
        mp_r(1) = 1d0                 ! air density
        mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 0.0d0                 ! air volume fraction
        u1 = 681.058d0                  ! m/s
        v1 = 0.0d0
        w1 = 0.0d0
        p1 = 1.9e9                      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      else

        if (sqrt(((pox(1) - 0.012)**2) + ((poy(1) - 0.012)**2) + ((poz(1) - 0.012)**2)) .le. 0.003d0) then
          mp_r(2) = 1000.00d0         ! water density
          mp_r(1) = 1d0                 ! air density
          mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 1.0d0                 ! air volume fraction
          u1 = 0.0d0                          ! m/s
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 100000                            ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(2) = 1000.0d0         ! water density
          mp_r(1) = 1d0                 ! air density
          mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 0.0d0                 ! air volume fraction
          u1 = 0.0d0                          ! m/s
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 100000                            ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      veccos(6) = mp_r(1)*mp_a(1)
      veccos(7) = mp_r(2)*mp_a(2)
      veccos(8) = mp_a(1)
    end if

    if (initcond .eq. 408) then
      if (pox(1) .gt. 0.10d0) then
        mp_r(1) = 6.03
        mp_r(2) = 1.658d0
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = -114.49d0
        v1 = 0.0d0
        w1 = 0.0d0
        p1 = 159060.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(((pox(1) - 0.079d0)**2) + ((poy(1) - 0.035d0)**2) + ((poz(1) - 0.035d0)**2)) .le. (0.0325d0/2.0d0)) then
          mp_r(1) = 6.03
          mp_r(2) = 1.204d0
          mp_a(1) = 1.0d0
          mp_a(2) = 0.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 6.03
          mp_r(2) = 1.204d0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      veccos(6) = mp_r(1)*mp_a(1)
      veccos(7) = mp_r(2)*mp_a(2)
      veccos(8) = mp_a(1)
    end if

    if (initcond .eq. 103) then
      r1 = rres
      p1 = pres
      s1 = sqrt((gamma*p1)/(r1))
      v1 = vvel
      w1 = wvel
      if (poy(1) .gt. 0.0d0) then
        u1 = uvel
      else
        u1 = 0.0d0
        v1 = -3.0
        p1 = press_outlet
      end if
      skin1 = (oo2)*((u1**2) + (v1**2) + (w1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = r1*(skin1 + ie1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = r1*w1
      veccos(5) = e1
      if (turbulence .eq. 1) then
        if (turbulencemodel .eq. 1) then
          veccos(6) = visc*turbinit
        end if
        if (turbulencemodel .eq. 2) then
          if (zero_turb_init .eq. 0) then
            veccos(6) = (1.5d0*i_turb_inlet*(ufreestream**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          end if
          if (zero_turb_init .eq. 1) then
            veccos(6) = (1.5d0*i_turb_inlet*(ufreestream**2))*r1
            veccos(7) = r1*veccos(6)/(10.0e-5*visc)
          end if
        end if
      end if
      if (passivescalar .gt. 0) then
        veccos(5 + turbulenceequations + 1:5 + turbulenceequations + passivescalar) = zero
      end if
    end if
  end subroutine initialise_euler3d

  subroutine initialise_euler2d(n, veccos, pox, poy, poz)
    implicit none
    integer, intent(in)::n
    real::acp, mscp, mvcp, vmcp, bcp, rcp, tcp, vfr, theta1
    real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx,vhx,amp,dvel,rgg,tt1,khi_slope,khi_b,theeta,reeta
    real::pr_radius,pr_beta,pr_machnumberfree,pr_pressurefree,pr_temperaturefree,pr_gammafree,pr_rgasfree,pr_xcenter,pr_ylength,pr_xlength,pr_ycenter,pr_densityfree,pr_cpconstant,pr_radiusvar,pr_velocityfree,pr_temperaturevar,drad
    integer::u_cond1, u_cond2, u_cond3, u_cond4, ix
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::veccos
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real, dimension(1:nof_species)::mp_r, mp_a, mp_ie
    veccos(:) = zero

    r1 = rres
    p1 = pres
    s1 = sqrt((gamma*p1)/(r1))
    u1 = uvel
    v1 = vvel
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = r1*(skin1 + ie1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
    if (turbulence .eq. 1) then
      if (turbulencemodel .eq. 1) then
        veccos(5) = visc*turbinit/r1
      end if
      if (turbulencemodel .eq. 2) then
        if (zero_turb_init .eq. 0) then
          veccos(5) = (1.5d0*(i_turb_inlet*ufreestream)**2)*r1
          veccos(6) = ufreestream/l_turb_inlet
          veccos(6) = (c_mu_inlet**(-0.25d0))*sqrt(veccos(5)) &
                      /l_turb_inlet*r1
        end if
        if (zero_turb_init .eq. 1) then
          veccos(5) = zero
          veccos(6) = ufreestream/l_turb_inlet
        end if
      end if
    end if
    if (passivescalar .gt. 0) then
      veccos(4 + turbulenceequations + 1:nof_variables + turbulenceequations + passivescalar) = zero
    end if
    if (initcond .eq. 95) then        !taylor green initial profile
      if ((poy(1) .ge. 0.25d0) .and. (poy(1) .le. 0.75d0)) then
      r1 = 2.0d0
      u1 = -0.5d0
      v1 = 0.01d0*sin(2.0d0*pi*(pox(1) - 0.5))
      p1 = 2.5
    end if
    if ((poy(1) .lt. 0.25d0) .or. (poy(1) .gt. 0.75d0)) then
      r1 = 1.0d0
      u1 = 0.5d0
      v1 = 0.01d0*sin(2.0d0*pi*(pox(1) - 0.5))
      p1 = 2.5
    end if
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = (p1/(gamma - 1)) + (r1*skin1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
    end if

    if (initcond .eq. 75) then        !taylor green initial profile
      if (((pox(1) .ge. (2.0d0)) .and. (pox(1) .le. (4.0d0))) .and. ((poy(1) .ge. (2.0d0)) .and. (poy(1) .le. (4.0d0)))) then
        r1 = 1.0d0
        v1 = 0.0
        u1 = 0.0
        p1 = 10.0d0
      else
        r1 = 0.2d0
        v1 = 0.0
        u1 = 0.0
        p1 = 1.0d0
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 25) then        !taylor green initial profile
      if (sqrt(((pox(1) - 1.0d0)**2) + ((poy(1) - 1.0d0)**2)) .le. 0.4d0) then
        r1 = 1.0d0
        v1 = 0.0
        u1 = 0.0
        p1 = 1.0d0
      else
        r1 = 0.125d0
        v1 = 0.0
        u1 = 0.0
        p1 = 0.1d0
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if
    if (initcond .eq. 31) then        !taylor green initial profile
      if (sqrt(((pox(1) - 1.0d0)**2) + ((poy(1) - 1.0d0)**2)) .le. 0.4d0) then
        r1 = 1.0d0
        v1 = 0.0
        u1 = 0.0
        p1 = 1.0d0
      else
        r1 = 0.125d0
        v1 = 0.0
        u1 = 0.0
        p1 = 0.1d0
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if
    if (initcond .eq. 65) then        !taylor green initial profile
      khi_slope = 15.0d0
      khi_b = tanh(khi_slope*(poy(1) - 1) + 7.5d0) - tanh(khi_slope*(poy(1) - 1) - 7.5d0)
      r1 = 0.5d0 + 0.75d0*khi_b
      u1 = 0.5*(khi_b - 1.d0)
      v1 = 0.1*sin(2.0d0*pi*(pox(1) - 1.0d0))
      p1 = 1.0d0
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 100) then        !taylor green initial profile
      acp = 0.075d0
      bcp = 0.175d0
      mscp = 1.5d0
      mvcp = 0.7d0
      vmcp = sqrt(gamma)*mvcp
      rcp = sqrt((pox(1) - 0.25d0)**2 + (poy(1) - 0.5d0)**2)
      if (rcp .le. acp) then
        vfr = vmcp*rcp/acp
        tcp = ((rcp - 0.175)*(gamma - 1.0d0)*(vfr**2)/(rcp*gamma)) + 1.0d0
        p1 = tcp**(gamma/(gamma - 1.0d0))
        r1 = tcp**(1.0d0/(gamma - 1.0d0))
        theta1 = atan((poy(1) - 0.5d0)/(pox(1) - 0.25d0))
        v1 = vfr*cos(theta1)
        u1 = 1.5d0*sqrt(gamma) - vfr*sin(theta1)
      else

        if ((rcp .ge. acp) .and. (rcp .le. bcp)) then
          vfr = vmcp*(acp/(acp**2 - bcp**2))*(rcp - ((bcp**2)/rcp))
          vfr = vmcp*rcp/acp
          tcp = ((rcp - 0.175)*(gamma - 1.0d0)*(vfr**2)/(rcp*gamma)) + 1.0d0
          p1 = tcp**(gamma/(gamma - 1.0d0))
          r1 = tcp**(1.0d0/(gamma - 1.0d0))
          theta1 = atan((poy(1) - 0.5d0)/(pox(1) - 0.25d0))
          v1 = vfr*cos(theta1)
          u1 = 1.5d0*sqrt(gamma) - vfr*sin(theta1)
        else
          vfr = 0.0d0
          r1 = 1.0d0
          u1 = 1.5d0*sqrt(gamma)
          v1 = 0.0
          p1 = 1.0d0
        end if
      end if

      if (pox(1) .gt. 0.5) then
        r1 = (9.0d0*gamma + 9.0d0)/(9.0d0*gamma - 1.0d0)
        u1 = sqrt(gamma)*((9.0d0*gamma - 1.0d0)/(6.0d0*(gamma + 1.0d0)))
        p1 = (7.0d0*gamma + 2.0d0)/(2.0d0*gamma + 2.0d0)
        v1 = 0.0d0
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 101) then        !shock density interaction
    if (pox(1) .lt. -4.0d0) then
      r1 = 3.8571d0
      u1 = 2.6294d0
      v1 = zero
      p1 = 10.333d0
    else
      r1 = (1.0d0 + 0.2d0*sin(5.0d0*pox(1)))
      u1 = zero
      v1 = zero
      p1 = 1
    end if
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = (p1/(gamma - 1)) + (r1*skin1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
    end if

    if (initcond .eq. 102) then        !shock density interaction
      if (pox(1) .lt. ((1.0d0/6.0d0) + (poy(1)/(sqrt(3.0d0))))) then
        r1 = 8.0d0
        u1 = 8.25*cos(pi/6.0d0)
        v1 = -8.25*sin(pi/6.0d0)
        p1 = 116.5
      else
        r1 = 1.4d0
        u1 = zero
        v1 = zero
        p1 = 1.0d0
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 104) then        !dmr_domain1
      if (pox(1) .lt. (1.0d0/6.0d0)) then
        r1 = 8.0d0
        u1 = 8.25d0
        v1 = 0.0d0
        p1 = 116.5
    else
        r1 = 1.4d0
        u1 = zero
        v1 = zero
        p1 = 1.0d0
    end if
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = (p1/(gamma - 1)) + (r1*skin1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
    end if

    if (initcond .eq. 30) then        !shock density interaction
      if (pox(1) .le. zero) then
        if (poy(1) .le. zero) then
          r1 = 0.138
          u1 = 1.206
          v1 = 1.206
          p1 = 0.029
       end if
       if (poy(1) .gt. zero) then
        r1 = 0.5323
        u1 = 1.206
        v1 = 0.0
        p1 = 0.3
      end if
    end if
    if (pox(1) .gt. zero) then
      if (poy(1) .le. zero) then
        r1 = 0.5323
        u1 = 0.0
        v1 = 1.206
        p1 = 0.3
      end if
    if (poy(1) .gt. zero) then
      r1 = 1.5
      u1 = 0.0
      v1 = 0.0
      p1 = 1.5
    end if
  end if
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = (p1/(gamma - 1)) + (r1*skin1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
 end if

    if (initcond .eq. 3) then        !shock density interaction
      if (pox(1) .ge. zero) then
        r1 = 1.1175d0
        u1 = 0.0d0
        v1 = 0.0d0
        p1 = 95000.0d0
      else
        r1 = 1.7522
        u1 = 166.34345
        v1 = 0.0d0
        p1 = 180219.75d0
    end if
    skin1 = (oo2)*((u1**2) + (v1**2))
    ie1 = ((p1)/((gamma - 1.0d0)*r1))
    e1 = (p1/(gamma - 1)) + (r1*skin1)
    veccos(1) = r1
    veccos(2) = r1*u1
    veccos(3) = r1*v1
    veccos(4) = e1
  end if

    if (initcond .eq. 123) then        !fast inviscid vortex evolution
      pr_radius = 0.005
      pr_beta = 0.2
      pr_machnumberfree = 0.05
      pr_pressurefree = 100000.0
      pr_temperaturefree = 300.0
      pr_gammafree = 1.4
      pr_rgasfree = 287.15
      pr_xcenter = 0.05; pr_ycenter = 0.05
      pr_xlength = 0.1; pr_ylength = 0.1
      pr_densityfree = pr_pressurefree/(pr_rgasfree*pr_temperaturefree)
      pr_cpconstant = pr_rgasfree*(pr_gammafree/(pr_gammafree - 1.0))
      pr_radiusvar = sqrt((pox(1) - pr_xcenter)**2 + (poy(1) - pr_ycenter)**2)/pr_radius
      pr_velocityfree = pr_machnumberfree*sqrt(pr_gammafree*pr_rgasfree*pr_temperaturefree)
      u1 = pr_velocityfree*(1.0 - (pr_beta*((poy(1) - pr_ycenter)/pr_radius)*exp((-pr_radiusvar**2)/2)))
      v1 = pr_velocityfree*((pr_beta*((poy(1) - pr_ycenter)/pr_radius)*exp((-pr_radiusvar**2)/2)))
      pr_temperaturevar = pr_temperaturefree - (((pr_velocityfree**2*pr_beta**2.0)/(2.0*pr_cpconstant))*exp(-pr_radiusvar**2))
      r1 = pr_densityfree*(pr_temperaturevar/pr_temperaturefree)**(1.0/(pr_gammafree - 1.0))
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((pr_gammafree - 1.0d0)*r1))
      e1 = (p1/(pr_gammafree - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 124) then        !fast inviscid vortex evolution
      pr_radius = 0.005
      pr_beta = 0.2
      pr_machnumberfree = 0.2
      pr_pressurefree = 100000.0
      pr_temperaturefree = 300.0
      pr_gammafree = 1.4
      pr_rgasfree = 287.15
      pr_xcenter = 0.05; pr_ycenter = 0.05
      pr_xlength = 0.1; pr_ylength = 0.1
      pr_densityfree = pr_pressurefree/(pr_rgasfree*pr_temperaturefree)
      pr_cpconstant = pr_rgasfree*(pr_gammafree/(pr_gammafree - 1.0))
      pr_radiusvar = sqrt((pox(1) - pr_xcenter)**2 + (poy(1) - pr_ycenter)**2)/pr_radius
      pr_velocityfree = pr_machnumberfree*sqrt(pr_gammafree*pr_rgasfree*pr_temperaturefree)
      u1 = pr_velocityfree*(1.0 - (pr_beta*((poy(1) - pr_ycenter)/pr_radius)*exp((-pr_radiusvar**2)/2)))
      v1 = pr_velocityfree*((pr_beta*((poy(1) - pr_ycenter)/pr_radius)*exp((-pr_radiusvar**2)/2)))
      pr_temperaturevar = pr_temperaturefree - (((pr_velocityfree**2*pr_beta**2.0)/(2.0*pr_cpconstant))*exp(-pr_radiusvar**2))
      r1 = pr_densityfree*(pr_temperaturevar/pr_temperaturefree)**(1.0/(pr_gammafree - 1.0))
      p1 = r1*pr_rgasfree*pr_temperaturevar
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((pr_gammafree - 1.0d0)*r1))
      e1 = (p1/(pr_gammafree - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 401) then
      if (sqrt(((pox(1) - 1.0d0)**2) + ((poy(1) - 1.0d0)**2)) .le. 0.4d0) then
        mp_r(1) = 1.0d0
        mp_r(2) = 0.125
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = zero
        v1 = zero
        p1 = 1.0d0
      else
        mp_r(1) = 1.0d0
        mp_r(2) = 0.125
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = zero
        v1 = zero
        p1 = 0.1d0
      end if
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if
    if (initcond .eq. 402) then
      if ((pox(1) .ge. 0.25) .and. (pox(1) .lt. 0.75)) then
        mp_r(1) = 10.0d0
        mp_r(2) = 1.0d0
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.5d0
        v1 = 0.0d0
        p1 = 1.0d0/1.4d0
      else
        mp_r(1) = 10.0d0
        mp_r(2) = 1.0d0
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 0.5d0
        v1 = 0.0d0
        p1 = 1.0d0/1.4d0
      end if
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 403) then
      mp_r(1) = 7.0d0
      mp_r(2) = 1.0d0
      mp_a(1) = 0.5d0 + 0.25d0*(sin(pi*((pox(1) - 1.0d0) + 0.5d0)))
      mp_a(2) = 1.0d0 - mp_a(1)
      u1 = 1.0d0
      v1 = 0.0d0
      p1 = 1.0d0/1.4d0
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 157) then
      if (sqrt(((pox(1) - 200.0e-6)**2) + ((poy(1) - 150.0e-6)**2)) .le. 50.0e-6) then
        mp_r(1) = 1.225
        mp_r(2) = 1000.0
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.0
        v1 = 0.0d0
        p1 = 100000.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        if (pox(1) .le. 100.0e-6) then
          p1 = 35e6
          mp_r(1) = 1.225
          mp_r(2) = 1000.0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 1647
          v1 = 0.0d0
        else
          mp_r(1) = 1.225
          mp_r(2) = 1000.0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0
          v1 = 0.0d0
          p1 = 100000.0d0
        end if
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 405) then
      drad = sqrt(((pox(1) + 0.05d0)**2) + ((poy(1) - 0.05d0)**2))
      if (pox(1) .lt. -0.1d0) then
        mp_r(1) = 0.166315789
        mp_r(2) = 1.658
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 114.49d0
        v1 = 0.0d0
        p1 = 159060.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        if (drad .le. 0.025d0) then
          mp_r(1) = 0.166315789d0
          mp_r(2) = 1.204d0
          mp_a(1) = 0.95d0
          mp_a(2) = 0.05d0
          u1 = 0.0d0
          v1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 0.166315789
          mp_r(2) = 1.204
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 422) then
      if (pox(1) .le. 0.05) then
        mp_r(1) = 1000.0d0         ! water density
        mp_r(2) = 3.85d0                 ! gas density
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 567.3d0                  ! m/s
        v1 = 0.0d0
        w1 = 0.0d0
        p1 = 664000.0d0                    ! pa
      else
        if ((sqrt(((pox(1) - 0.0576)**2) + ((poy(1) - 0.0576)**2)) .le. 0.0048d0)) then
          mp_r(1) = 1000.0d0         ! water density
          mp_r(2) = 1.2                 ! gas density
          mp_a(1) = 1.0d0
          mp_a(2) = 0.0d0
          u1 = 0.0d0                  ! m/s
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101000                   ! pa
        else
          mp_r(1) = 1000.0d0         ! water density
          mp_r(2) = 1.20                 ! gas density
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0                  ! m/s
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101000                   ! pa
        end if
      end if

      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 406) then
      if (pox(1) .gt. 0.1d0) then
        mp_r(1) = 0.166315789d0
        mp_r(2) = 1.658d0
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = -114.49d0
        v1 = 0.0d0
        p1 = 159060.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        u_cond1 = 0; u_cond2 = 0; u_cond3 = 0; u_cond4 = 0
        if (((pox(1) .ge. 0.01d0) .and. (pox(1) .le. 0.02)) .or. ((pox(1) .ge. 0.03d0) .and. (pox(1) .le. 0.04))) then
          u_cond1 = 1
        end if
        if ((poy(1) .ge. 0.05d0) .and. (poy(1) .le. 0.075)) then
          u_cond2 = 1
        end if
        if ((poy(1) .le. 0.05d0) .and. (poy(1) .gt. 0.02)) then
          u_cond3 = 1
        end if
        if (u_cond3 .eq. 1) then
          if ((sqrt(((pox(1)-0.025d0)**2)+((poy(1)-0.05d0)**2)).le.0.015d0).and.(sqrt(((pox(1)-0.025d0)**2)+((poy(1)-0.05d0)**2)).ge.0.005d0)) then
            u_cond4 = 1
          end if
        end if
        if (((u_cond1 .eq. 1) .and. (u_cond2 .eq. 1)) .or. (u_cond4 .eq. 1)) then
          mp_r(1) = 0.166315789d0
          mp_r(2) = 1.204d0
          mp_a(1) = 0.95d0
          mp_a(2) = 0.05d0
          u1 = 0.0d0
          v1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 0.166315789d0
          mp_r(2) = 1.204d0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          p1 = 101325
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 408) then
      if (pox(1) .gt. 0.10d0) then
        mp_r(1) = 6.03
        mp_r(2) = 1.658d0
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = -114.49d0
        v1 = 0.0d0
        w1 = 0.0d0
        p1 = 159060.0d0
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(((pox(1) - 0.079d0)**2) + ((poy(1) - 0.035d0)**2)) .le. (0.0325d0/2.0d0)) then
          mp_r(1) = 6.03
          mp_r(2) = 1.204d0
          mp_a(1) = 1.0d0
          mp_a(2) = 0.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          w1 = 0.0d0
          p1 = 101325.0d0
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 6.03
          mp_r(2) = 1.204d0
          mp_a(1) = 0.0d0
          mp_a(2) = 1.0d0
          u1 = 0.0d0
          v1 = 0.0d0
          p1 = 101325.0d0
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 410) then
      if (pox(1) .le. -0.007d0) then
        mp_r(1) = 1225.6d0 ! water density
        mp_r(2) = 1.2d0 ! air density
        mp_a(1) = 1.0d0 ! water volume fraction (everything is water here)
        mp_a(2) = 0.0d0 ! air volume fraction
        u1 = 542.76d0          ! m/s
        v1 = 0.0d0
        p1 = 1.6e9      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(((pox(1))**2) + ((poy(1))**2)) .le. 0.003d0) then
          mp_r(1) = 1225.6d0 ! water density
          mp_r(2) = 1.2d0 ! air density
          mp_a(1) = 0.0d0 ! water volume fraction (everything is water here)
          mp_a(2) = 1.0d0 ! air volume fraction
          u1 = 0.0d0          ! m/s
          v1 = 0.0d0
          p1 = 101325    ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(1) = 1000.0d0 ! water density
          mp_r(2) = 1.2d0 ! air density
          mp_a(1) = 1.0d0 ! water volume fraction (everything is water here)
          mp_a(2) = 0.0d0 ! air volume fraction
          u1 = 0.0d0     ! m/s
          v1 = 0.0d0
          p1 = 101325     ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 411) then
      if (pox(1) .le. 0.0066d0) then
        mp_r(2) = 1323.65d0         ! water density
        mp_r(1) = 1d0                 ! air density
        mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 0.0d0                 ! air volume fraction
        u1 = 681.058d0                  ! m/s
        v1 = 0.0d0
        p1 = 1.9e9                      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        if (sqrt(((pox(1) - 0.012)**2) + ((poy(1) - 0.012)**2)) .le. 0.003d0) then
          mp_r(2) = 1000.00d0         ! water density
          mp_r(1) = 1d0                 ! air density
          mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 1.0d0                 ! air volume fraction
          u1 = 0.0d0                          ! m/s
          v1 = 0.0d0
          p1 = 100000                            ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(2) = 1000.0d0         ! water density
          mp_r(1) = 1d0                 ! air density
          mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 0.0d0                 ! air volume fraction
          u1 = 0.0d0                          ! m/s
          v1 = 0.0d0
          p1 = 100000                            ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 412) then
      if (pox(1) .le. 0.00) then
        mp_r(1) = 1.241         ! air density
        mp_r(2) = 0.991                 ! water density
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.0                  ! m/s
        v1 = 0.0d0
        p1 = 2.753                     ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        mp_r(1) = 1.241         ! air density
        mp_r(2) = 0.991                 ! water density
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 0.0                  ! m/s
        v1 = 0.0d0
        p1 = 3.059*10e-4                     ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 420) then
      mp_r(1) = 7.0         ! air density
      mp_r(2) = 1.0                 ! water density
      mp_a(1) = 0.5d0 + 0.25*sin(pi*(((pox(1) - 0.5d0)*2.0d0)))
      mp_a(2) = 1.0d0 - mp_a(1)
      u1 = 1.0d0                  ! m/s
      v1 = 0.0d0
      p1 = 1.0d0/1.4d0                    ! pa
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 421) then
      if (pox(1) .lt. 0.0d0) then
        mp_r(1) = 1.241         ! air density
        mp_r(2) = 0.991                 ! water density
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.0d0                  ! m/s
        v1 = 0.0d0
        p1 = 2.753d0                    ! pa
      else
        mp_r(1) = 1.241         ! air density
        mp_r(2) = 0.991                 ! water density
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 0.0d0                  ! m/s
        v1 = 0.0d0
        p1 = 3.059*(10.0e-4)                    ! pa
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 430) then
      if (poy(1) .lt. 0.386d0) then
        mp_r(1) = 3.483         ! air density
        mp_r(2) = 867                 ! water density
        mp_a(1) = 0.0d0
        mp_a(2) = 1.0d0
        u1 = 0.0d0                  ! m/s
        v1 = 2
        p1 = 300000                    ! pa
      else
        mp_r(1) = 3.483         ! air density
        mp_r(2) = 867                ! water density
        mp_a(1) = 1.0d0
        mp_a(2) = 0.0d0
        u1 = 0.0d0                  ! m/s
        v1 = 0.0d0
        p1 = 300000                    ! pa
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
      mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
      skin1 = (oo2)*((u1**2) + (v1**2))
      e1 = (r1*skin1) + ie1
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 133) then
      if (pox(1) .lt. 0.0d0) then
        p1 = 195557.25
        r1 = p1/(350.5d0*287.058d0)
        u1 = 168.62
        v1 = 0.0d0
        rhc1 = r1
        rhc2 = u1
        rhc3 = v1
        rhc4 = p1
      else
        u1 = 0.0d0
        v1 = 0.0d0
        p1 = 101325
        r1 = p1/(288.15d0*287.058d0)
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 266) then
      if (poy(1) .gt. 1.0d0) then
        p1 = 20000
        r1 = 0.41
        u1 = 850
        v1 = 0.0d0
      else
        u1 = 0.0d0
        v1 = 0.0d0
        p1 = 100000
        r1 = 1.225
      end if
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 444) then
      if (poy(1) .le. 0.00025) then   !post shock concidions
        mp_r(2) = 1323.65d0         ! water density
        mp_r(1) = 1d0                 ! air density
        mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 0.0d0                 ! air volume fraction
        u1 = 0.0                            ! m/s
        v1 = 681.058d0
        p1 = 1.9e9                      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else
        mp_r(2) = 1000.00d0         ! water density
        mp_r(1) = 1d0                 ! air density
        mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 0.0d0                 ! air volume fraction
        u1 = 0.0d0                          ! m/s
        v1 = 0.0d0
        p1 = 100000                            ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
        do ix = 1, nof_bubbles
          if (sqrt(((pox(1) - bubble_centre(ix, 1))**2) + ((poy(1) - bubble_centre(ix, 2))**2)) .le. bubble_radius(ix)) then
            mp_r(2) = 1000.00d0         ! water density
            mp_r(1) = 1d0                 ! air density
            mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
            mp_a(1) = 1.0d0                 ! air volume fraction
            u1 = 0.0d0                          ! m/s
            v1 = 0.0d0
            p1 = 100000                            ! pa
            r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
            mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
            mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
            ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
            skin1 = (oo2)*((u1**2) + (v1**2))
            e1 = (r1*skin1) + ie1
          end if
        end do
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if

    if (initcond .eq. 222) then        !shock density interaction
      r1 = 1.0d0
      p1 = 1.0e-6
      reeta = -1.0d0
      theeta = atan(poy(1)/pox(1))
      u1 = reeta*cos(theeta)
      v1 = reeta*sin(theeta)
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 10000) then        !shock density interaction
      r1 = 0.5d0
      p1 = 0.4127
      u1 = 0.0
      v1 = 0.0
      skin1 = (oo2)*((u1**2) + (v1**2))
      ie1 = ((p1)/((gamma - 1.0d0)*r1))
      e1 = (p1/(gamma - 1)) + (r1*skin1)
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
    end if

    if (initcond .eq. 470) then
      if (pox(1) .lt. 1.0) then   !post shock concidions
        mp_r(2) = 1.0d0         ! water density
        mp_r(1) = 1.0d0                 ! air density
        mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
        mp_a(1) = 1.0d0                 ! air volume fraction
        u1 = 0.0                            ! m/s
        v1 = 0.0
        p1 = 1.0                      ! pa
        r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
        mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
        mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
        ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
        skin1 = (oo2)*((u1**2) + (v1**2))
        e1 = (r1*skin1) + ie1
      else

        if (poy(1) .gt. 1.5) then
          mp_r(2) = 1.0d0         ! water density
          mp_r(1) = 0.125d0                 ! air density
          mp_a(2) = 0.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 1.0d0                 ! air volume fraction
          u1 = 0.0                            ! m/s
          v1 = 0.0
          p1 = 0.1                      ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        else
          mp_r(2) = 1.0d0         ! water density
          mp_r(1) = 0.125d0                 ! air density
          mp_a(2) = 1.0d0                 ! water volume fraction (everything is water here)
          mp_a(1) = 0.0d0                 ! air volume fraction
          u1 = 0.0                            ! m/s
          v1 = 0.0
          p1 = 0.1                      ! pa
          r1 = (mp_r(1)*mp_a(1)) + (mp_r(2)*mp_a(2))
          mp_ie(1) = ((p1 + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
          mp_ie(2) = ((p1 + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
          ie1 = (mp_ie(1)*mp_a(1)) + (mp_ie(2)*mp_a(2))
          skin1 = (oo2)*((u1**2) + (v1**2))
          e1 = (r1*skin1) + ie1
        end if
      end if
      veccos(1) = r1
      veccos(2) = r1*u1
      veccos(3) = r1*v1
      veccos(4) = e1
      veccos(5) = mp_r(1)*mp_a(1)
      veccos(6) = mp_r(2)*mp_a(2)
      veccos(7) = mp_a(1)
    end if
  end subroutine initialise_euler2d
end module profile
