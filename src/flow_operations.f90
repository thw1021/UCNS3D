module flow_operations
  use declaration
  use mpiinfo
  use transform
  implicit none

contains
  subroutine mrfswitch(n, iconsidered, facex, pointx, pox, poy)
    ! @brief
    ! this subroutine  check if the element is on the rotational/stationary reference frame and update the mrf_origin and srf_velocity srf accordingly
    implicit none
    integer, intent(in)::n, iconsidered, facex, pointx
    ! type(local_recon3),allocatable,dimension(:),intent(inout)::ilocal_recon3
    ! output:
    ! ilocal_recon3%mrf_origin; ilocal_recon3%mrf_velocity; ilocal_recon3%rotvel, ilocal_recon3%mrf
    real, dimension(3) ::mrf_origin, mrf_velocity, rotvel
    integer:: rotframe_on
    integer::ninv
    real, dimension(1:dimensiona), intent(inout)::pox, poy
    !internal variables
    real, dimension(3) :: p1p2, pc, popc, po, pgp, poz !element coordinates, roation_axys, cylinder_center_coordinates, vector_element_center, rotational velocity at gaussian points, gausian points coordinates
    real :: d1, d2, r1, theta, dpopc
    po = (pox(1:3))
    pgp = (poy(1:3))
    !body
    do ninv = 1, nrotors
      pc(1:3) = (point1_gl(ninv, 1:3) + point2_gl(ninv, 1:3))/2  !center of cylinder
      p1p2(1:3) = point2_gl(ninv, 1:3) - point1_gl(ninv, 1:3)          !axysvector
      popc(1:3) = po(1:3) - pc(1:3)              ! vector elelement-centre
      dpopc = ((po(1) - pc(1))**2 + (po(2) - pc(2))**2 + (po(3) - pc(3))**2)**0.5 !distance between element and center
      theta = acos((dot_product(popc, p1p2))/(sqrt(popc(1)**2 + popc(2)**2 + popc(3)**2)*sqrt(p1p2(1)**2 + p1p2(2)**2 + p1p2(3)**2))) !angle between element vector and axys
      d2 = dpopc*abs(cos(theta))
      r1 = dpopc*abs(sin(theta))
      d1 = ((point1_gl(ninv, 1) - pc(1))**2 + (point1_gl(ninv, 2) - pc(2))**2 + (point1_gl(ninv, 3) - pc(3))**2)**0.5
      if ((d1 .ge. d2) .and. (r1 .le. radius_gl(ninv))) then
        rotframe_on = 1
        mrf_origin(1:3) = pc(1:3)
        pox(1:3) = pgp(1:3) - mrf_origin(1:3)
        mrf_velocity(1:3) = mrf_rot_gl(ninv)*(p1p2)/(p1p2(1)**2 + p1p2(2)**2 + p1p2(3)**2)**0.5
        poy(1:3) = mrf_velocity(1:3)
        rotvel(1:3) = vect_function(pox, poy)
        go to 606
      else
        mrf_origin(1:3) = 0.0
        rotframe_on = 0
        mrf_velocity(1:3) = 0.0
        rotvel(1:3) = 0.0
      end if
    end do
606 continue
    ilocal_recon3(iconsidered)%mrf_origin = mrf_origin
    ilocal_recon3(iconsidered)%mrf_velocity = mrf_velocity
    ilocal_recon3(iconsidered)%rotvel(facex, pointx, 1:3) = rotvel
    ilocal_recon3(iconsidered)%mrf = rotframe_on
  end subroutine mrfswitch

  function fluxeval2d(leftv)
    implicit none
    real, dimension(1:nof_variables)::fluxeval2d
    real, dimension(1:nof_variables), intent(in)::leftv
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    r = leftv(1)
    u = leftv(2)
    v = leftv(3)
    p = leftv(4)
    gm = gamma
    !kinetic energy first!
    skin = (oo2)*((u**2) + (v**2))
    !internal energy
    ien = ((p)/((gm - 1.0d0)*r))
    !total energy
    e = r*(skin + ien)
    fluxeval2d(1) = r*u
    fluxeval2d(2) = (r*(u**2)) + p
    fluxeval2d(3) = r*u*v
    fluxeval2d(4) = u*(e + p)
  end function

  function fluxeval3d(leftv)
    implicit none
    real, dimension(1:nof_variables)::fluxeval3d
    real, dimension(1:nof_variables), intent(in)::leftv
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    r = leftv(1)
    u = leftv(2)
    v = leftv(3)
    w = leftv(4)
    p = leftv(5)
    gm = gamma
    !kinetic energy first!
    skin = (oo2)*((u**2) + (v**2) + (w**2))
    !internal energy
    ien = ((p)/((gm - 1.0d0)*r))
    !total energy
    e = r*(skin + ien)
    fluxeval3d(1) = r*u
    fluxeval3d(2) = (r*(u**2)) + p
    fluxeval3d(3) = r*u*v
    fluxeval3d(4) = r*u*w
    fluxeval3d(5) = u*(e + p)
  end function

  subroutine cons2prim(n, leftv, mp_pinfl, gammal)
    ! @brief
    ! this subroutine transforms one vector of conservative variables to primitive variables
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables)::temps
    real, dimension(1:nof_variables), intent(inout)::leftv
    real, intent(inout)::mp_pinfl, gammal
    real::oodensity, mp_density, mp_stiff
    real::p_sat, p_tol, rho_g, rho_l, ss_g, ss_l, pp, p_gl, void_frac, p_temp
    real, dimension(nof_species)::mp_ar, mp_ie
    if (nof_variables .gt. 1) then
      if (dimensiona .eq. 3) then
        p_sat = 2000
        p_tol = 10e-5
        if (governingequations .eq. -1) then
          mp_density = (leftv(6) + leftv(7)) !total density of mixture
          mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = leftv(4)*oodensity
          mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(8)*mp_pinf(1)) + ((1.0d0 - leftv(8))*mp_pinf(2))
          temps(5) = (((gammal - 1.0d0))*((leftv(5)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))) - mp_stiff
          temps(6) = leftv(6)
          temps(7) = leftv(7)
          temps(8) = leftv(8)
          if (cavitation .eq. 1) then
            rho_g = leftv(6)/leftv(8)
            rho_l = leftv(7)/leftv(8)
            ss_g = sqrt(gamma_in(1)*(temps(5) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(5) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(leftv(8)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-leftv(8)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(5)
            if ((temps(5) .gt. p_tol) .and. (temps(5) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(5) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(5) = p_temp
          end if
          leftv(1:nof_variables) = temps(1:nof_variables)
        else
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = leftv(4)*oodensity
          temps(5) = ((gamma - 1.0d0))*((leftv(5)) - oo2*leftv(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))
          leftv(1:nof_variables) = temps(1:nof_variables)
        end if
      else
        p_sat = 2000
        p_tol = 10e-5
        if ((governingequations .eq. -1) .and. (viscous_s .ne. 1)) then
          mp_density = (leftv(5) + leftv(6)) !total density of mixture
          mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          ! mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((leftv(7)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammal-1.0d0)
          ! mp_pinfl=(leftv(7)*mp_pinf(1))+((leftv(7)-1.0d0)*mp_pinf(2))
          mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(7)*mp_pinf(1)) + ((1.0d0 - leftv(7))*mp_pinf(2))
          temps(4) = (((gammal - 1.0d0))*((leftv(4)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2)))) - mp_stiff
          temps(5) = leftv(5)
          temps(6) = leftv(6)
          temps(7) = leftv(7)
          if (cavitation .eq. 1) then
            rho_g = leftv(5)/leftv(7)
            rho_l = leftv(6)/leftv(7)
            ss_g = sqrt(gamma_in(1)*(temps(4) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(4) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(leftv(7)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-leftv(7)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(4)
            if ((temps(4) .gt. p_tol) .and. (temps(4) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(4) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(4) = p_temp
          end if
          leftv(1:nof_variables) = temps(1:nof_variables)
        else
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = ((gamma - 1.0d0))*((leftv(4)) - oo2*leftv(1)*(((temps(2))**2) + ((temps(3))**2)))
          leftv(1:nof_variables) = temps(1:nof_variables)
        end if
      end if
    end if
  end subroutine cons2prim

  subroutine cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    ! @brief
    ! this subroutine transforms two vector of conservative variables to primitive variables
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables)::temps
    real::oodensity, mp_density, mp_stiff
    real::p_sat, p_tol, rho_g, rho_l, ss_g, ss_l, pp, p_gl, void_frac, p_temp
    real, dimension(1:nof_variables), intent(inout)::leftv
    real, intent(inout)::mp_pinfl, gammal
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, intent(inout)::mp_pinfr, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    if (nof_variables .gt. 1) then
      if (dimensiona .eq. 3) then
        p_sat = 2000
        p_tol = 10e-5
        if (governingequations .eq. -1) then
          mp_density = (leftv(6) + leftv(7)) !total density of mixture
          mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = leftv(4)*oodensity
          !mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((leftv(8)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammal-1.0d0)
          mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(8)*mp_pinf(1)) + ((1.0d0 - leftv(8))*mp_pinf(2))
          temps(5) = (((gammal - 1.0d0))*((leftv(5)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))) - mp_stiff
          temps(6) = leftv(6)
          temps(7) = leftv(7)
          temps(8) = leftv(8)
          if (cavitation .eq. 1) then
            rho_g = leftv(6)/leftv(8)
            rho_l = leftv(7)/leftv(8)
            ss_g = sqrt(gamma_in(1)*(temps(5) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(5) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(leftv(8)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-leftv(8)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(5)
            if ((temps(5) .gt. p_tol) .and. (temps(5) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(5) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(5) = p_temp
          end if
          leftv(1:nof_variables) = temps(1:nof_variables)
          mp_density = (rightv(6) + rightv(7)) !total density of mixture
          mp_ar(1) = rightv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - rightv(8))/(gamma_in(2) - 1.0d0)
          gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = rightv(2)*oodensity
          temps(3) = rightv(3)*oodensity
          temps(4) = rightv(4)*oodensity
          ! mp_stiff=((rightv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((rightv(8)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammar-1.0d0)
          ! mp_pinfl=(rightv(8)*mp_pinf(1))+((rightv(8)-1.0d0)*mp_pinf(2))
          mp_stiff=((rightv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-rightv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammar-1.0d0)
          mp_pinfr = (rightv(8)*mp_pinf(1)) + ((1.0d0 - rightv(8))*mp_pinf(2))
          temps(5) = (((gammar - 1.0d0))*((rightv(5)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))) - mp_stiff
          temps(6) = rightv(6)
          temps(7) = rightv(7)
          temps(8) = rightv(8)
          if (cavitation .eq. 1) then
            rho_g = rightv(6)/rightv(8)
            rho_l = rightv(7)/rightv(8)
            ss_g = sqrt(gamma_in(1)*(temps(5) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(5) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(rightv(8)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-rightv(8)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(5)
            if ((temps(5) .gt. p_tol) .and. (temps(5) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(5) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(5) = p_temp
          end if
          rightv(1:nof_variables) = temps(1:nof_variables)
        else
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = leftv(4)*oodensity
          temps(5) = ((gamma - 1.0d0))*((leftv(5)) - oo2*leftv(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))
          leftv(1:nof_variables) = temps(1:nof_variables)
          oodensity = 1.0d0/rightv(1)
          temps(1) = rightv(1)
          temps(2) = rightv(2)*oodensity
          temps(3) = rightv(3)*oodensity
          temps(4) = rightv(4)*oodensity
          temps(5) = ((gamma - 1.0d0))*((rightv(5)) - oo2*rightv(1)*(((temps(2))**2) + ((temps(3))**2) + ((temps(4))**2)))
          rightv(1:nof_variables) = temps(1:nof_variables)
        end if
      else    !2d
        p_sat = 2000
        p_tol = 10e-5
        if ((governingequations .eq. -1) .and. (viscous_s .ne. 1)) then
          mp_density = (leftv(5) + leftv(6)) !total density of mixture
          mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumption
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          !mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((leftv(7)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammal-1.0d0)
          !mp_pinfl=(leftv(7)*mp_pinf(1))+((leftv(7)-1.0d0)*mp_pinf(2))
          mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(7)*mp_pinf(1)) + ((1.0d0 - leftv(7))*mp_pinf(2))
          temps(4) = (((gammal - 1.0d0))*((leftv(4)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2)))) - mp_stiff
          temps(5) = leftv(5)
          temps(6) = leftv(6)
          temps(7) = leftv(7)
          if (cavitation .eq. 1) then
            rho_g = leftv(5)/leftv(7)
            rho_l = leftv(6)/leftv(7)
            ss_g = sqrt(gamma_in(1)*(temps(4) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(4) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(leftv(7)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-leftv(7)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(4)
            if ((temps(4) .gt. p_tol) .and. (temps(4) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(4) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(4) = p_temp
          end if
          leftv(1:nof_variables) = temps(1:nof_variables)
          mp_density = (rightv(5) + rightv(6)) !total density of mixture
          mp_ar(1) = rightv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - rightv(7))/(gamma_in(2) - 1.0d0)
          gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          oodensity = 1.0d0/mp_density
          temps(1) = mp_density
          temps(2) = rightv(2)*oodensity
          temps(3) = rightv(3)*oodensity
          !mp_stiff=((rightv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((rightv(7)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammar-1.0d0)
          mp_stiff=((rightv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-rightv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammar-1.0d0)
          mp_pinfr = (rightv(7)*mp_pinf(1)) + ((1.0d0 - rightv(7))*mp_pinf(2))
          temps(4) = (((gammar - 1.0d0))*((rightv(4)) - oo2*temps(1)*(((temps(2))**2) + ((temps(3))**2)))) - mp_stiff
          temps(5) = rightv(5)
          temps(6) = rightv(6)
          temps(7) = rightv(7)
          if (cavitation .eq. 1) then
            rho_g = rightv(5)/rightv(7)
            rho_l = rightv(6)/rightv(7)
            ss_g = sqrt(gamma_in(1)*(temps(4) + mp_pinf(1))/rho_g)
            ss_l = sqrt(gamma_in(2)*(temps(4) + mp_pinf(2))/rho_l)
            p_gl = rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g - rho_l)/((rho_g*rho_g*ss_g*ss_g) - (rho_l*rho_l*ss_l*ss_l))
            void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(rightv(7)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-rightv(7)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
            p_temp = temps(4)
            if ((temps(4) .gt. p_tol) .and. (temps(4) .lt. p_sat)) then
              p_temp = p_sat + p_gl*log(void_frac)
            end if
            if (temps(4) .lt. p_tol) then
              p_temp = p_tol
            end if
            temps(4) = p_temp
          end if
          rightv(1:nof_variables) = temps(1:nof_variables)
        else
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*oodensity
          temps(3) = leftv(3)*oodensity
          temps(4) = ((gamma - 1.0d0))*((leftv(4)) - oo2*leftv(1)*(((temps(2))**2) + ((temps(3))**2)))
          leftv(1:nof_variables) = temps(1:nof_variables)
          oodensity = 1.0d0/rightv(1)
          temps(1) = rightv(1)
          temps(2) = rightv(2)*oodensity
          temps(3) = rightv(3)*oodensity
          temps(4) = ((gamma - 1.0d0))*((rightv(4)) - oo2*rightv(1)*(((temps(2))**2) + ((temps(3))**2)))
          rightv(1:nof_variables) = temps(1:nof_variables)
        end if
      end if
    end if
  end subroutine cons2prim2

  subroutine lmacht(n, leftv, rightv)
    ! @brief
    ! this subroutine applies the low-mach number correction to two vectors of conserved variables
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables), intent(inout)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables), intent(inout)::rightv
    real::mp_pinfr, gammar
    real::q2l, q2r, uul, uur, vvl, vvr, wwr, wwl, rhol, rhor, etal, etar, duu, dvv, dww
    real::mach2, mach, cma, dus, dvs, dws, diff, c1o2, ssl, ssr, ppl, ppr, eel, eer, tole, mlm
    tole = zero
    eel = leftv(5)
    eer = rightv(5)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    c1o2 = 0.50d0
    rhol = leftv(1)
    uul = leftv(2)
    vvl = leftv(3)
    wwl = leftv(4)
    ppl = leftv(5)
    rhor = rightv(1)
    uur = rightv(2)
    vvr = rightv(3)
    wwr = rightv(4)
    ppr = rightv(5)
    q2l = (uul*uul) + (vvl*vvl) + (wwl*wwl)
    q2r = (uur*uur) + (vvr*vvr) + (wwr*wwr)
    if (multispecies .eq. 1) then
      ssl = sqrt((leftv(5) + mp_pinfl)*gammal/leftv(1))
      ssr = sqrt((rightv(5) + mp_pinfr)*gammar/rightv(1))
    else
      ssl = ((gamma*ppl)/(rhol))
      ssr = ((gamma*ppr)/(rhor))
    end if
    cma = 1.0d0
    duu = uur - uul
    dvv = vvr - vvl
    dww = wwr - wwl
    mach2 = max(q2l/ssl, q2r/ssr)
    mach = sqrt(mach2)
    mach = min(cma*mach, 1.0d0)
    dus = uur + uul
    dvs = vvr + vvl
    dws = wwr + wwl
    duu = mach*duu
    dvv = mach*dvv
    dww = mach*dww
    diff = c1o2*duu
    uul = (dus*c1o2 - diff)
    uur = (dus*c1o2 + diff)
    if (lmach_style .eq. 1) then
      diff = c1o2*dvv
      vvl = (dvs*c1o2 - diff)
      vvr = (dvs*c1o2 + diff)
      diff = c1o2*dww
      wwl = (dws*c1o2 - diff)
      wwr = (dws*c1o2 + diff)
    end if
    leftv(1) = rhol
    leftv(2) = uul*rhol
    leftv(3) = vvl*rhol
    leftv(4) = wwl*rhol
    leftv(5) = eel
    rightv(1) = rhor
    rightv(2) = uur*rhor
    rightv(3) = vvr*rhor
    rightv(4) = wwr*rhor
    rightv(5) = eer
  end subroutine lmacht

  subroutine lmacht2d(n, leftv, rightv)
    ! @brief
    ! this subroutine applies the low-mach number correction to two vectors of conserved variables 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables), intent(inout)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables), intent(inout)::rightv
    real::mp_pinfr, gammar
    real::q2l, q2r, uul, uur, vvl, vvr, wwr, wwl, rhol, rhor, etal, etar, duu, dvv, dww
    real::mach2, mach, cma, dus, dvs, dws, diff, c1o2, ssl, ssr, ppl, ppr, eel, eer, tole, mlm
    tole = tolsmall
    eel = leftv(4)
    eer = rightv(4)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    c1o2 = 0.5d0
    rhol = leftv(1)
    uul = leftv(2)
    vvl = leftv(3)
    ppl = leftv(4)
    rhor = rightv(1)
    uur = rightv(2)
    vvr = rightv(3)
    ppr = rightv(4)
    q2l = (uul*uul) + (vvl*vvl)
    q2r = (uur*uur) + (vvr*vvr)
    if (multispecies .eq. 1) then
      ssl = sqrt((leftv(4) + mp_pinfl)*gammal/leftv(1))
      ssr = sqrt((rightv(4) + mp_pinfr)*gammar/rightv(1))
    else
      ssl = ((gamma*ppl)/(rhol))
      ssr = ((gamma*ppr)/(rhor))
    end if
    cma = 1.0d0
    duu = uur - uul
    dvv = vvr - vvl
    mach2 = max(q2l/ssl, q2r/ssr)
    mach = sqrt(mach2)
    mach = min(cma*mach, 1.0d0)
    dus = uur + uul
    dvs = vvr + vvl
    duu = mach*duu
    dvv = mach*dvv
    diff = c1o2*duu
    uul = (dus*c1o2) - diff
    uur = (dus*c1o2) + diff
    if (lmach_style .eq. 1) then
      diff = c1o2*dvv
      vvl = (dvs*c1o2 - diff)
      vvr = (dvs*c1o2 + diff)
    end if
    leftv(1) = rhol
    leftv(2) = uul*rhol
    leftv(3) = vvl*rhol
    leftv(4) = eel
    rightv(1) = rhor
    rightv(2) = uur*rhor
    rightv(3) = vvr*rhor
    rightv(4) = eer
  end subroutine lmacht2d

  subroutine prim2cons(n, leftv)
    ! @brief
    ! ! this subroutine transforms one vector of primitive variables to conservative variables
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables)::temps
    real::oodensity, skin1, ie1, mp_density, mp_stiff
    real, dimension(1:nof_variables), intent(inout)::leftv
    real::mp_pinfl, gammal
    real, dimension(nof_species)::mp_ar, mp_ie
    if (nof_variables .gt. 1) then
      if (dimensiona .eq. 3) then
        if (governingequations .eq. -1) then
          mp_density = (leftv(6) + leftv(7)) !total density of mixture
          mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = leftv(2)*temps(1)
          temps(3) = leftv(3)*temps(1)
          temps(4) = leftv(4)*temps(1)
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2) + (leftv(4)**2))
          !mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((leftv(8)-1.0d0)*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))/(gammal-1.0d0)
          mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          ie1 = ((leftv(5) + mp_stiff)/((gammal - 1.0d0)*temps(1)))
          temps(5) = temps(1)*(ie1 + skin1)
          temps(6:8) = leftv(6:8)
          leftv(1:nof_variables) = temps(1:nof_variables)
        else
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2) + (leftv(4)**2))
          ie1 = ((leftv(5))/((gamma - 1.0d0)*leftv(1)))
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*leftv(1)
          temps(3) = leftv(3)*leftv(1)
          temps(4) = leftv(4)*leftv(1)
          temps(5) = leftv(1)*(ie1 + skin1)
          leftv(1:nof_variables) = temps(1:nof_variables)
        end if
      else
        if ((governingequations .eq. -1) .and. (viscous_s .ne. 1)) then
          mp_density = (leftv(5) + leftv(6)) !total density of mixture
          mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = leftv(2)*temps(1)
          temps(3) = leftv(3)*temps(1)
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2))
          mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(7)*mp_pinf(1)) + ((1.0d0 - leftv(7))*mp_pinf(2))
          ie1 = ((leftv(4) + mp_stiff)/((gammal - 1.0d0)*temps(1)))
          temps(4) = temps(1)*(ie1 + skin1)
          temps(5:7) = leftv(5:7)
          leftv(1:nof_variables) = temps(1:nof_variables)
        else
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2))
          ie1 = ((leftv(4))/((gamma - 1.0d0)*leftv(1)))
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*leftv(1)
          temps(3) = leftv(3)*leftv(1)
          temps(4) = leftv(1)*(ie1 + skin1)
          leftv(1:nof_variables) = temps(1:nof_variables)
        end if
      end if
    end if
  end subroutine prim2cons

  subroutine prim2cons2(n, leftv, rightv)
    ! @brief
    ! this subroutine transforms two vectors of primitive variables to conservative variables
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables)::temps
    real::oodensity, skin1, ie1, mp_density, mp_stiff
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(inout)::leftv, rightv
    real::mp_pinfl, gammal, mp_pinfr, gammar
    if (nof_variables .gt. 1) then
      if (dimensiona .eq. 3) then
        if (governingequations .eq. -1) then
          mp_density = (leftv(6) + leftv(7)) !total density of mixture
          mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = leftv(2)*temps(1)
          temps(3) = leftv(3)*temps(1)
          temps(4) = leftv(4)*temps(1)
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2) + (leftv(4)**2))
          mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.  0d0))*mp_pinf(2)))*(gammal-1.0d0)
          ie1 = ((leftv(5) + mp_stiff)/((gammal - 1.0d0)*temps(1)))
          temps(5) = temps(1)*(ie1 + skin1)
          temps(6:8) = leftv(6:8)
          leftv(1:nof_variables) = temps(1:nof_variables)
          mp_density = (rightv(6) + rightv(7)) !total density of mixture
          mp_ar(1) = rightv(8)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
          gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = rightv(2)*temps(1)
          temps(3) = rightv(3)*temps(1)
          temps(4) = rightv(4)*temps(1)
          skin1 = (oo2)*((rightv(2)**2) + (rightv(3)**2) + (rightv(4)**2))
          mp_stiff=((rightv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-rightv(8))*(gamma_in(2)/(gamma_in(2)  -1.0d0))*mp_pinf(2)))*(gammar-1.0d0)
          ie1 = ((rightv(5) + mp_stiff)/((gammar - 1.0d0)*temps(1)))
          temps(5) = temps(1)*(ie1 + skin1)
          temps(6:8) = rightv(6:8)
          rightv(1:nof_variables) = temps(1:nof_variables)
        else
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2) + (leftv(4)**2))
          ie1 = ((leftv(5))/((gamma - 1.0d0)*leftv(1)))
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*leftv(1)
          temps(3) = leftv(3)*leftv(1)
          temps(4) = leftv(4)*leftv(1)
          temps(5) = leftv(1)*(ie1 + skin1)
          leftv(1:nof_variables) = temps(1:nof_variables)
          skin1 = (oo2)*((rightv(2)**2) + (rightv(3)**2) + (rightv(4)**2))
          ie1 = ((rightv(5))/((gamma - 1.0d0)*rightv(1)))
          oodensity = 1.0d0/rightv(1)
          temps(1) = rightv(1)
          temps(2) = rightv(2)*rightv(1)
          temps(3) = rightv(3)*rightv(1)
          temps(4) = rightv(4)*rightv(1)
          temps(5) = rightv(1)*(ie1 + skin1)
          rightv(1:nof_variables) = temps(1:nof_variables)
        end if
      else
        if ((governingequations .eq. -1) .and. (viscous_s .ne. 1)) then
          mp_density = (leftv(5) + leftv(6)) !total density of mixture
          mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
          gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = leftv(2)*temps(1)
          temps(3) = leftv(3)*temps(1)
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2))
          mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.  0d0))*mp_pinf(2)))*(gammal-1.0d0)
          mp_pinfl = (leftv(7)*mp_pinf(1)) + ((1.0d0 - leftv(7))*mp_pinf(2))
          ie1 = ((leftv(4) + mp_stiff)/((gammal - 1.0d0)*temps(1)))
          temps(4) = temps(1)*(ie1 + skin1)
          temps(5:7) = leftv(5:7)
          leftv(1:nof_variables) = temps(1:nof_variables)
          mp_density = (rightv(5) + rightv(6)) !total density of mixture
          mp_ar(1) = rightv(7)/(gamma_in(1) - 1.0d0)
          mp_ar(2) = (1.0d0 - rightv(7))/(gamma_in(2) - 1.0d0)
          gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
          temps(1) = mp_density
          temps(2) = rightv(2)*temps(1)
          temps(3) = rightv(3)*temps(1)
          skin1 = (oo2)*((rightv(2)**2) + (rightv(3)**2))
          mp_stiff=((rightv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-rightv(7))*(gamma_in(2)/(gamma_in(2)  -1.0d0))*mp_pinf(2)))*(gammar-1.0d0)
          mp_pinfr = (rightv(7)*mp_pinf(1)) + ((1.0d0 - rightv(7))*mp_pinf(2))
          ie1 = ((rightv(4) + mp_stiff)/((gammar - 1.0d0)*temps(1)))
          temps(4) = temps(1)*(ie1 + skin1)
          temps(5:7) = rightv(5:7)
          rightv(1:nof_variables) = temps(1:nof_variables)
        else
          skin1 = (oo2)*((leftv(2)**2) + (leftv(3)**2))
          ie1 = ((leftv(4))/((gamma - 1.0d0)*leftv(1)))
          oodensity = 1.0d0/leftv(1)
          temps(1) = leftv(1)
          temps(2) = leftv(2)*leftv(1)
          temps(3) = leftv(3)*leftv(1)
          temps(4) = leftv(1)*(ie1 + skin1)
          leftv(1:nof_variables) = temps(1:nof_variables)
          skin1 = (oo2)*((rightv(2)**2) + (rightv(3)**2))
          ie1 = ((rightv(4))/((gamma - 1.0d0)*rightv(1)))
          oodensity = 1.0d0/rightv(1)
          temps(1) = rightv(1)
          temps(2) = rightv(2)*rightv(1)
          temps(3) = rightv(3)*rightv(1)
          temps(4) = rightv(1)*(ie1 + skin1)
          rightv(1:nof_variables) = temps(1:nof_variables)
        end if
      end if
    end if
  end subroutine prim2cons2

  function inflow(initcond, pox, poy, poz)
    ! @brief
    ! this function applies a prescribed boundary condition to  the inflow in 3d
    implicit none
    real, dimension(1:nof_variables)::inflow
    integer, intent(in)::initcond
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf
    real:: theta_0, vtang, vradial, gammar
    real::mp_density, mp_stiff
    real, dimension(nof_species)::mp_ar, mp_ie
    if (governingequations .eq. -1) then
      p = pres
      u = uvel
      v = vvel
      w = wvel
      mp_ar(1) = mp_a_in(1)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = mp_a_in(2)/(gamma_in(2) - 1.0d0)
      gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumption
      gm = gammar
      r = (mp_r_in(1)*mp_a_in(1)) + (mp_r_in(2)*mp_a_in(2))
      mp_ie(1) = ((p + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ien = (mp_ie(1)*mp_a_in(1)) + (mp_ie(2)*mp_a_in(2))
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      e = (r*skin) + ien
      inflow(1) = r
      inflow(2) = r*u
      inflow(3) = r*v
      inflow(4) = r*w
      inflow(5) = e
      inflow(6) = mp_r_in(1)*mp_a_in(1)
      inflow(7) = mp_r_in(2)*mp_a_in(2)
      inflow(8) = mp_a_in(1)
    else
      r = rres
      gm = gamma
      p = pres
      u = uvel
      v = vvel
      w = wvel
      if (initcond .eq. 10000) then
        if (sqrt(((poy(1) - 0.0)**2) + ((poz(1) - 0.5)**2)) .le. 0.05) then
          p = 0.4127
          r = 5
          u = 30.0
          v = 0.0d0
          w = 0.0d0
        end if
      end if
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      ien = ((p)/((gm - 1.0d0)*r))
      e = r*(skin + ien)
      inflow(1) = r
      inflow(2) = r*u
      inflow(3) = r*v
      inflow(4) = r*w
      inflow(5) = e
    end if
    if (swirl .eq. 1) then
      if (pox(1) .lt. -0.03) then
        xf = pox(1)
        yf = poy(1)
        zf = poz(1)
        theta_0 = atan2(zf, yf)
        vtang = 18.0375d0
        vradial = -12.63d0
        u = 0.0d0
        v = -vtang*sin(theta_0) + vradial*cos(theta_0)
        w = vtang*cos(theta_0) + vradial*sin(theta_0)
      else
        u = 70.06d0
        v = 0.0d0
        w = 0.0d0
      end if
      r = rres
      gm = gamma
      p = pres
      s = sqrt((gm*p)/(r))
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      ien = ((p)/((gm - 1.0d0)*r))
      e = r*(skin + ien)
      inflow(1) = r
      inflow(2) = r*u
      inflow(3) = r*v
      inflow(4) = r*w
      inflow(5) = e
    end if
  end function inflow

  function vect_function(pox, poy)
    ! @brief
    ! this makes a multipliciation between two vectors
    implicit none
    real, dimension(3)::vect_function
    real, dimension(1:dimensiona), intent(in)::pox, poy
    vect_function(1) = (poy(2)*pox(3)) - (poy(3)*pox(2))
    vect_function(2) = (poy(3)*pox(1)) - (poy(1)*pox(3))
    vect_function(3) = (poy(1)*pox(2)) - (poy(2)*pox(1))
  end function vect_function

  function inflow2d(initcond, pox, poy)
    ! @brief
    ! this function applies a prescribed boundary condition to  the inflow in 2d
    implicit none
    real, dimension(1:nof_variables)::inflow2d
    integer, intent(in)::initcond
    real, dimension(1:2), intent(in)::pox, poy
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ps
    real::xf, yf, zf, lit_a, lit_o
    real:: theta_0, vtang, vradial, gammar
    real::mp_density, mp_stiff
    real, dimension(nof_species)::mp_ar, mp_ie
    if (governingequations .eq. -1) then
      p = pres
      u = uvel
      v = vvel
      if (initcond .eq. 430) then
        v=(179299.375638680*(t**5)) - (82455.0868677361*(t**4)) + (14299.8472891299*(t**3)) - (1281.65548492021*(t**2)) + (62.2260666329356*t) - (0.0419033282181554)
      end if
      mp_ar(1) = mp_a_in(1)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = mp_a_in(2)/(gamma_in(2) - 1.0d0)
      gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumption
      gm = gammar
      r = (mp_r_in(1)*mp_a_in(1)) + (mp_r_in(2)*mp_a_in(2))
      mp_ie(1) = ((p + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ien = (mp_ie(1)*mp_a_in(1)) + (mp_ie(2)*mp_a_in(2))
      skin = (oo2)*((u**2) + (v**2))
      e = (r*skin) + ien
      inflow2d(1) = r
      inflow2d(2) = r*u
      inflow2d(3) = r*v
      inflow2d(4) = e
      inflow2d(5) = mp_r_in(1)*mp_a_in(1)
      inflow2d(6) = mp_r_in(2)*mp_a_in(2)
      inflow2d(7) = mp_a_in(1)
    else
      r = rres
      gm = gamma
      p = pres
      u = uvel
      v = vvel
      if (initcond .eq. 133) then
        p = 195557.25
        r = p/(350.5d0*287.058d0)
        u = 168.62
        v = 0.0d0
      end if
      if (initcond .eq. 10000) then
        if ((poy(1) .ge. -0.05) .and. (poy(1) .le. 0.05)) then
          p = 0.4127
          r = 5
          u = 30.0
          v = 0.0d0
        end if
      end if
      if (initcond .eq. 790) then
        r = (2.4d0*6**2)/((0.4*6**2) + 2)
        u = (6*sqrt(1.4))*(70/(2.4*36))
        v = 0.0d0
        p = (2.8*36 - 0.4)/(2.4)
      end if
      skin = (oo2)*((u**2) + (v**2))
      ien = ((p)/((gm - 1.0d0)*r))
      e = r*(skin + ien)
      inflow2d(1) = r
      inflow2d(2) = r*u
      inflow2d(3) = r*v
      inflow2d(4) = e
    end if
  end function inflow2d

  function outflow2d(initcond, pox, poy)
    ! @brief
    ! this function applies a prescribed boundary condition to  the outflow in 2d
    implicit none
    real, dimension(1:nof_variables)::outflow2d
    integer, intent(in)::initcond
    real, dimension(1:2), intent(in)::pox, poy
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    if (governingequations .eq. -1) then
      p = pres
      u = uvel
      v = vvel
      mp_ar(1) = mp_a_in(1)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = mp_a_in(2)/(gamma_in(2) - 1.0d0)
      gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumption
      gm = gammar
      r = (mp_r_in(1)*mp_a_in(1)) + (mp_r_in(2)*mp_a_in(2))
      mp_ie(1) = ((p + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ien = (mp_ie(1)*mp_a_in(1)) + (mp_ie(2)*mp_a_in(2))
      skin = (oo2)*((u**2) + (v**2))
      e = (r*skin) + ien
      outflow2d(1) = r
      outflow2d(2) = r*u
      outflow2d(3) = r*v
      outflow2d(4) = e
      outflow2d(5) = mp_r_in(1)*mp_a_in(1)
      outflow2d(6) = mp_r_in(2)*mp_a_in(2)
      outflow2d(7) = mp_a_in(1)
    else
      r = rres
      gm = gamma
      p = pres
      u = uvel
      v = vvel
      skin = (oo2)*((u**2) + (v**2))
      ien = ((p)/((gm - 1.0d0)*r))
      e = r*(skin + ien)
      outflow2d(1) = r
      outflow2d(2) = r*u
      outflow2d(3) = r*v
      outflow2d(4) = e
    end if
  end function outflow2d

  function outflow(initcond, pox, poy, poz)
    ! @brief
    ! this function applies a prescribed boundary condition to  the outflow in 3d
    implicit none
    real, dimension(1:nof_variables)::outflow
    integer, intent(in)::initcond
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    if (governingequations .eq. -1) then
      p = pres
      u = uvel
      v = vvel
      w = wvel
      mp_ar(1) = mp_a_in(1)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = mp_a_in(2)/(gamma_in(2) - 1.0d0)
      gammar = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumption
      gm = gammar
      r = (mp_r_in(1)*mp_a_in(1)) + (mp_r_in(2)*mp_a_in(2))
      mp_ie(1) = ((p + (gamma_in(1)*mp_pinf(1)))/((gamma_in(1) - 1.0d0)))
      mp_ie(2) = ((p + (gamma_in(2)*mp_pinf(2)))/((gamma_in(2) - 1.0d0)))
      ien = (mp_ie(1)*mp_a_in(1)) + (mp_ie(2)*mp_a_in(2))
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      e = (r*skin) + ien
      outflow(1) = r
      outflow(2) = r*u
      outflow(3) = r*v
      outflow(4) = r*w
      outflow(5) = e
      outflow(6) = mp_r_in(1)*mp_a_in(1)
      outflow(7) = mp_r_in(2)*mp_a_in(2)
      outflow(8) = mp_a_in(1)
    else
      r = rres
      gm = gamma
      p = pres
      if (initcond .eq. 977) then
        p = 101325
      end if
      u = uvel
      v = vvel
      w = wvel
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      ien = ((p)/((gm - 1.0d0)*r))
      e = r*(skin + ien)
      outflow(1) = r
      outflow(2) = r*u
      outflow(3) = r*v
      outflow(4) = r*w
      outflow(5) = e
    end if
  end function outflow

  function outflow2(initcond, pox, poy, poz)
    ! @brief
    ! this function applies a prescribed boundary condition to  the outflow in 3d
    implicit none
    real, dimension(1:nof_variables)::outflow2
    integer, intent(in)::initcond
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf
    real, dimension(nof_species)::mp_ar, mp_ie
    r = rres
    gm = gamma
    p = press_outlet
    u = uvel
    v = vvel
    w = wvel
    skin = (oo2)*((u**2) + (v**2) + (w**2))
    ien = ((p)/((gm - 1.0d0)*r))
    e = r*(skin + ien)
    outflow2(1) = r
    outflow2(2) = r*u
    outflow2(3) = r*v
    outflow2(4) = r*w
    outflow2(5) = e
  end function outflow2

  function bleed2d(iconsidered, facex, pox, poy)
    ! @brief
    ! this function applies a prescribed boundary condition to  the outflow in 3d
    implicit none
    real, dimension(1:nof_variables)::bleed2d
    integer, intent(in)::iconsidered, facex
    real, dimension(1:dimensiona), intent(in)::pox, poy
    integer::ibleedn
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf
    real::cell_area
    !------coefficients for bleed-----!
    real::bleed_qsonic_s, bleed_mdotsonic_s, bleed_area, bleed_mdotsonic, bleed_region
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    cell_area = ielem(n, iconsidered)%surf(facex)
    !leftv is the input in conservative variables that we transform to primitive variables
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    r = leftv(1)
    u = leftv(2)
    v = leftv(3)
    p = leftv(4)
    !now find the bleed number conditions for this face
    ibleedn = ielem(n, iconsidered)%bleedn(facex)
    bleed_region = sqrt(((bleed_start(ibleedn, 1) - bleed_end(ibleedn, 1))**2) + ((bleed_start(ibleedn, 2) - bleed_end(ibleedn, 2))**2))
    !now compute the bleed area
    bleed_area = bleed_porosity(ibleedn)*bleed_region
    !now compute the bleed m dot sonic-s mass flow rate eq.17 , https://doi.org/10.2514/1.b37474
    bleed_mdotsonic_s = bleed_area*p*(sqrt(gamma*r/p))*(((gamma + 1.0d0)/(2))**((gamma + 1)/(2*(1 - gamma))))
    !now compute qsonic eq. 22
    bleed_qsonic_s = 0.598 + 0.0307*(bleed_plenum(ibleedn)/p) - 0.5936*((bleed_plenum(ibleedn)/p)**2)
    !equation 16
    bleed_mdotsonic = bleed_qsonic_s*bleed_mdotsonic_s
    !now that i have computed the bleed mass flow rate i have to distribute to the surface area of this tagged cell in the normal direction
    call prim2cons2(n, leftv, rightv)
    rightv(1:nof_variables) = leftv(1:nof_variables)
    call rotatef2d(n, cright_rot, rightv, angle1, angle2)
    cright_rot(2) = bleed_mdotsonic/cell_area
    call rotateb2d(n, rightv, cright_rot, angle1, angle2)
    bleed2d(1:nof_variables) = rightv(1:nof_variables)
  end function bleed2d

  function bleed3d(iconsidered, facex, pox, poy, poz)
    ! @brief
    ! this function applies a prescribed boundary condition to  the outflow in 3d
    implicit none
    real, dimension(1:nof_variables)::bleed3d
    integer, intent(in)::iconsidered, facex
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    integer::ibleedn
    real::p, u, v, w, e, r, s, gm, skin, ien, pi
    real::xf, yf, zf
    real::cell_area
    !------coefficients for bleed-----!
    real::bleed_qsonic_s, bleed_mdotsonic_s, bleed_area, bleed_mdotsonic, bleed_region
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    print *, "not ready yet"
    stop
    bleed3d(1:nof_variables) = rightv(1:nof_variables)
  end function bleed3d

  function pass_inlet(initcond, pox, poy, poz)
    ! @brief
    ! this function applies a prescribed boundary condition to  the inlet for a passive scalar
    implicit none
    real, dimension(1:passivescalar)::pass_inlet
    integer, intent(in)::initcond
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    pass_inlet(1:passivescalar) = 1.0d0*rres
  end function pass_inlet

  function pass_inlet2d(initcond, pox, poy)
    ! @brief
    ! this function applies a prescribed boundary condition to  the inlet for a passive scalar in 2d
    implicit none
    real, dimension(1:passivescalar)::pass_inlet2d
    integer, intent(in)::initcond
    real, dimension(1:dimensiona), intent(in)::pox, poy
    pass_inlet2d(1:passivescalar) = 1.0d0
  end function pass_inlet2d

  subroutine shear_x(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the shear stresses in x-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%grads(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    if (rframe .eq. 0) then
      shear_temp = -ssx/(0.5*rres*ufreestream*ufreestream)
    else
      shear_temp = -ssx/(0.5*rres*v_ref*v_ref)
    end if
  end subroutine shear_x

  subroutine shear_y(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the shear stresses in y-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%grads(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    if (rframe .eq. 0) then
      shear_temp = -ssy/(0.5*rres*ufreestream*ufreestream)
    else
      shear_temp = -ssy/(0.5*rres*v_ref*v_ref)
    end if
  end subroutine shear_y

  subroutine shear_z(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the shear stresses in z-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%grads(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    if (rframe .eq. 0) then
      shear_temp = -ssz/(0.5*rres*ufreestream*ufreestream)
    else
      shear_temp = -ssz/(0.5*rres*v_ref*v_ref)
    end if
  end subroutine shear_z

  subroutine shear_x_av(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the average shear stresses in x-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    integer::ind1
    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%gradsav(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    shear_temp = -ssx/(0.5*rres*ufreestream*ufreestream)
  end subroutine shear_x_av

  subroutine shear_y_av(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the average shear stresses in y-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    integer::ind1
    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%gradsav(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    shear_temp = -ssy/(0.5*rres*ufreestream*ufreestream)
  end subroutine shear_y_av

  subroutine shear_z_av(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the avergage shear stresses in z-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    integer::ind1
    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if
    vortet1(1:3, 1:3) = ilocal_recon3(iconsidered)%gradsav(1:3, 1:3)
    ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
    vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
    wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = (cos(angle1)*sin(angle2))
    ny = (sin(angle1)*sin(angle2))
    nz = (cos(angle2))
    leftv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(ind1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
    tauyy = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
    tauzz = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
    tauyx = (uy + vx)
    tauzx = (wx + uz)
    tauzy = (vz + wy)
    ssx = (viscl(1)*((nx*tauxx) + (ny*tauyx) + (nz*tauzx)))
    ssy = (viscl(1)*((nx*tauyx) + (ny*tauyy) + (nz*tauzy)))
    ssz = (viscl(1)*((nx*tauzx) + (ny*tauzy) + (nz*tauzz)))
    shear_temp = -ssz/(0.5*rres*ufreestream*ufreestream)
  end subroutine shear_z_av

  subroutine shear_x2d(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the shear stresses in x-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    integer::ind1
    integer::gqi_points, im
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    ssx = zero; ssp = zero; ssy = zero;
    gqi_points = qp_line_n
    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
    do im = 1, gqi_points
    if (ielem(n, iconsidered)%ggs .eq. 1) then
      vortet1(1:2, 1:2) = ilocal_recon3(iconsidered)%grads(1:2, 1:2)
    else
      vortet1(1, 1:2) = ilocal_recon3(iconsidered)%uleftv(1:2, 2, facex, im)
      vortet1(2, 1:2) = ilocal_recon3(iconsidered)%uleftv(1:2, 3, facex, im)
    end if
    ux = vortet1(1, 1); uy = vortet1(1, 2)
    vx = vortet1(2, 1); vy = vortet1(2, 2)
    angle1 = ielem(n, iconsidered)%faceanglex(facex)
    angle2 = ielem(n, iconsidered)%faceangley(facex)
    nx = angle1
    ny = angle2
    leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    rightv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland2d(n, leftv, rightv, viscl, laml)
    ssx = zero; ssp = zero; ssy = zero; ssz = zero
    tauxx = 2.0d0*ux
    tauyy = 2.0d0*vy
    tauyx = (uy + vx)
    ssx = ssx + ((viscl(1)*((ny*tauyx)))*wequa2d(im))
    ssy = ssy + ((viscl(1)*((nx*tauyx)))*wequa2d(im))
    end do
    shear_temp = -ssx/(0.5*rres*ufreestream*ufreestream)
  end subroutine shear_x2d

  subroutine shear_y2d(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the shear stresses in y-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
    real::ssx, ssy, ssz, ssp
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    integer::ind1
    integer::gqi_points, im
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    ssx = zero; ssp = zero; ssy = zero;
    gqi_points = qp_line_n
    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
    do im = 1, gqi_points
      if (ielem(n, iconsidered)%ggs .eq. 1) then
        vortet1(1:2, 1:2) = ilocal_recon3(iconsidered)%grads(1:2, 1:2)
      else
        vortet1(1, 1:2) = ilocal_recon3(iconsidered)%uleftv(1:2, 2, facex, im)
        vortet1(2, 1:2) = ilocal_recon3(iconsidered)%uleftv(1:2, 3, facex, im)
      end if
      ux = vortet1(1, 1); uy = vortet1(1, 2)
      vx = vortet1(2, 1); vy = vortet1(2, 2)
      angle1 = ielem(n, iconsidered)%faceanglex(facex)
      angle2 = ielem(n, iconsidered)%faceangley(facex)
      nx = angle1
      ny = angle2
      leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      rightv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
      call sutherland2d(n, leftv, rightv, viscl, laml)
      ssx = zero; ssp = zero; ssy = zero; ssz = zero
      tauxx = 2.0d0*ux
      tauyy = 2.0d0*vy
      tauyx = (uy + vx)
      ssx = ssx + ((viscl(1)*((ny*tauyx)))*wequa2d(im))
      ssy = ssy + ((viscl(1)*((nx*tauyx)))*wequa2d(im))
    end do
    shear_temp = -ssx/(0.5*rres*ufreestream*ufreestream)
  end subroutine shear_y2d

  subroutine shear_x2d_av(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the average shear stresses in x-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    shear_temp = 0.0d0
  end subroutine shear_x2d_av

  subroutine shear_y2d_av(iconsidered, facex, shear_temp)
    ! @brief
    ! this subroutine computes the average shear stresses in y-axis
    implicit none
    integer, intent(in)::iconsidered, facex
    real, intent(inout)::shear_temp
    shear_temp = 0.0d0
  end subroutine shear_y2d_av

  subroutine sutherland(n, leftv, rightv, viscl, laml)
    ! @brief
    ! this subroutine computes the viscosity according to sutherland's law
    implicit none
    real, dimension(1:nof_variables), intent(in)::leftv, rightv
    real, dimension(1:4), intent(inout)::viscl, laml
    integer, intent(in)::n
    real::kinetic, u, v, w, t0l, t1l, t0r, t1r
    t1l = leftv(5)/leftv(1)
    t0l = pres/rres
    t1r = rightv(5)/rightv(1)
    t0r = pres/rres
    viscl(1) = visc*((t1l/t0l)**betaas)*((t0l + (suther*t0l))/(t1l + (suther*t0l)))
    viscl(2) = visc*((t1r/t0r)**betaas)*((t0r + (suther*t0r))/(t1r + (suther*t0r)))
    laml(1) = viscl(1)*gamma/(prandtl*(gamma - 1.d0))
    laml(2) = viscl(2)*gamma/(prandtl*(gamma - 1.d0))
  end subroutine sutherland

  subroutine sutherland2d(n, leftv, rightv, viscl, laml)
    ! @brief
    ! this subroutine computes the viscosity according to sutherland's law
    implicit none
    real, dimension(1:nof_variables), intent(in)::leftv, rightv
    real, dimension(1:4), intent(inout)::viscl, laml
    integer, intent(in)::n
    real::kinetic, u, v, w, t0l, t1l, t0r, t1r
    t1l = leftv(4)/leftv(1)
    t0l = pres/rres
    t1r = rightv(4)/rightv(1)
    t0r = pres/rres
    viscl(1) = visc*((t1l/t0l)**betaas)*((t0l + (suther*t0l))/(t1l + (suther*t0l)))
    viscl(2) = visc*((t1r/t0r)**betaas)*((t0r + (suther*t0r))/(t1r + (suther*t0r)))
    laml(1) = viscl(1)*gamma/(prandtl*(gamma - 1.d0))
    laml(2) = viscl(2)*gamma/(prandtl*(gamma - 1.d0))
  end subroutine sutherland2d

  subroutine vortexcalc(n)
    ! @brief
    ! this subroutine computes the q-criterion
    implicit none
    integer, intent(in)::n
    integer::kmaxe, i, ihgt, ihgj
    real::snorm, onorm
    real, dimension(3, 3)::tvort, svort, ovort
    real, dimension(1:dims, 1:dims)::vortet1
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      vortet1(1:3, 1:3) = ilocal_recon3(i)%grads(1:3, 1:3)
      do ihgt = 1, 3; do ihgj = 1, 3
          tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
        end do; end do
      svort = 0.5d0*(vortet1 + tvort)
      ovort = 0.5d0*(vortet1 - tvort)
      snorm = sqrt((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                   (svort(1, 3)*svort(1, 3)) + (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)) + (svort(2, 3)*svort(2, 3)) &
                   + (svort(3, 1)*svort(3, 1)) + (svort(3, 2)*svort(3, 2)) + (svort(3, 3)*svort(3, 3)))
      onorm = sqrt((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + (ovort(1, 3)*ovort(1, 3)) + &
                   (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2)) + (ovort(2, 3)*ovort(2, 3)) + (ovort(3, 1)*ovort(3, 1)) + &
                   (ovort(3, 2)*ovort(3, 2)) + (ovort(3, 3)*ovort(3, 3)))

      ielem(n, i)%vortex(1) = (0.5d0*((onorm**2) - (snorm**2)))
    end do
  end subroutine vortexcalc

  subroutine enstrophy_calc(n)
    ! @brief
    ! this subroutine computes the q-criterion
    implicit none
    integer, intent(in)::n
    integer::kmaxe, i, ihgt, ihgj
    real::snorm, onorm
    real, dimension(3, 3)::tvort, svort, ovort
    real, dimension(3, 3)::taul, taur, tau
    real, dimension(3)::q, nnn, nall
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp, tempxx
    real, dimension(1:dims, 1:dims)::vortet1
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:4)::viscl, laml
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      vortet1(1:3, 1:3) = ilocal_recon3(i)%grads(1:3, 1:3)
      do ihgt = 1, 3
        do ihgj = 1, 3
          tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
        end do
      end do
      ovort = (vortet1 - tvort)
      onorm = ((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + (ovort(1, 3)*ovort(1, 3)) + &
               (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2)) + (ovort(2, 3)*ovort(2, 3)) + (ovort(3, 1)*ovort(3, 1)) + (ovort(3, 2)*ovort(3, 2)) + (ovort(3, 3)*ovort(3, 3)))
      if (boundtype .eq. 1) then
        ielem(n, i)%vortex(2) = (0.5d0*(onorm*u_c(i)%val(1, 1)))
      else
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        rightv(1:nof_variables) = leftv(1:nof_variables)
        call sutherland(n, leftv, rightv, viscl, laml)
        ux = ilocal_recon3(i)%grads(1, 1); uy = ilocal_recon3(i)%grads(1, 2); uz = ilocal_recon3(i)%grads(1, 3);
        vx = ilocal_recon3(i)%grads(2, 1); vy = ilocal_recon3(i)%grads(2, 2); vz = ilocal_recon3(i)%grads(2, 3);
        wx = ilocal_recon3(i)%grads(3, 1); wy = ilocal_recon3(i)%grads(3, 2); wz = ilocal_recon3(i)%grads(3, 3);
        ! tau_xx
        taul(1, 1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
        ! tau_yy
        taul(2, 2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
        ! tau_zz
        taul(3, 3) = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
        ! tau_xy
        taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)
        ! tau_xz
        taul(1, 3) = (wx + uz); taul(3, 1) = taul(1, 3)
        ! tau_yz
        taul(2, 3) = (vz + wy); taul(3, 2) = taul(2, 3)
        snorm = ((wy - vz)**2) + ((uz - wx)**2) + ((vx - uy)**2)
        ielem(n, i)%vortex(2) = ielem(n, i)%totvolume*snorm*viscl(1)
        ielem(n, i)%vortex(3) = (4.0/3.0)*viscl(1)*((ux + vy + wz)**2)*ielem(n, i)%totvolume
      end if
    end do
  end subroutine enstrophy_calc

  subroutine vortexcalc2d(n)
    ! @brief
    ! this subroutine computes the q criterion for 2d
    implicit none
    integer, intent(in)::n
    integer::kmaxe, i, ihgt, ihgj
    real::snorm, onorm
    real, dimension(2, 2)::tvort, svort, ovort
    real, dimension(1:dims, 1:dims)::vortet1
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      vortet1(1:2, 1:2) = ilocal_recon3(i)%grads(1:2, 1:2)
      do ihgt = 1, 2
        do ihgj = 1, 2
          tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
        end do
      end do
      svort = 0.5d0*(vortet1 + tvort)
      ovort = 0.5d0*(vortet1 - tvort)
      snorm = sqrt((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                   (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)))
      onorm = sqrt((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + &
                   (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2)))
      ielem(n, i)%vortex(1) = (0.5d0*((onorm**2) - (snorm**2)))
    end do
  end subroutine vortexcalc2d

  subroutine boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
    ! @brief
    ! this subroutine applies the boundary condition to each bounded cell
    implicit none
    integer, intent(in)::n, b_code, iconsidered, facex
    real, dimension(1:nof_variables), intent(inout)::leftv, rightv
    integer, intent(inout)::ibfc
    real, dimension(1:nof_variables), intent(in)::srf_speedrot, srf_speed
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real, intent(in)::angle1, angle2, nx, ny, nz
    real, dimension(turbulenceequations), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::cright_rot, cleft_rot
    real, dimension(1:nof_variables)::subson1, subson2, subson3
    real::sps, skins, ikins, vel, vnb
    real::mp_pinfl, mp_pinfr, gammal, gammar
    select case(b_code)
      case (1)!inflow subsonic or supersonic will be chosen based on mach number
        if (boundtype .eq. 0) then        !supersonic
          rightv(1:nof_variables) = inflow(initcond, pox, poy, poz)
        else                !subsonic
          rightv(1:nof_variables) = inflow(initcond, pox, poy, poz)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          subson1(1:nof_variables) = rightv(1:nof_variables)
          subson2(1:nof_variables) = leftv(1:nof_variables)
          sps = sqrt((gamma*subson2(5))/(subson2(1)))
          vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
          call prim2cons2(n, leftv, rightv)
          if (vel/(sps + tolsmall) .gt. 1.0d0) then        !supersonic
            rightv(1:nof_variables) = inflow(initcond, pox, poy, poz)
          else                !subsonic
            subson3(5) = 0.5*((subson1(5)) + (subson2(5)) - (subson2(1)*sps*((nx*(subson1(2) - subson2(2))) + (ny*(subson1(3) - subson2(3))) &+ (nz*(subson1(4) - subson2(4))))))
            subson3(1) = subson1(1) + (subson3(5) - subson1(5))/(sps**2)
            subson3(2) = subson1(2) - (nx*(subson1(5) - subson3(5)))/(sps*subson2(1))
            subson3(3) = subson1(3) - (ny*(subson1(5) - subson3(5)))/(sps*subson2(1))
            subson3(4) = subson1(4) - (nz*(subson1(5) - subson3(5)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            rightv(4) = subson3(4)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
            ikins = subson3(5)/((gamma - 1.0d0)*(subson3(1)))
            rightv(5) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulencemodel .eq. 1) then
              cturbr(1) = visc*turbinit
            end if
            if (turbulencemodel .eq. 2) then
              cturbr(1) = (1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
              cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
              if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
                cturbr(1) = (1.5d0*i_turb_inlet*(kinit_srf**2))*rightv(1)!k initialization
                cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
              end if
            end if
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = pass_inlet(initcond, pox, poy, poz)*rightv(1)
            end if
          end if
        end if

      case (2)!outflow subsonic or supersonic will be chosen based on mach number
        if (boundtype .eq. 0) then
          rightv(1:nof_variables) = leftv(1:nof_variables)
        else
          rightv(1:nof_variables) = outflow(initcond, pox, poy, poz)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          subson1(1:nof_variables) = rightv(1:nof_variables)
          subson2(1:nof_variables) = leftv(1:nof_variables)
          sps = sqrt((gamma*subson2(5))/(subson2(1)))
          vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
          sps = sqrt((gamma*subson2(5))/(subson2(1)))
          call prim2cons2(n, leftv, rightv)
          if (vel/(sps + tolsmall) .gt. 1.0d0) then        !supersonic
            call prim2cons2(n, leftv, rightv)
            rightv(1:nof_variables) = leftv(1:nof_variables)
          else
            subson3(5) = subson1(5)
            subson3(1) = subson2(1) + (subson3(5) - subson2(5))/(sps**2)
            subson3(2) = subson2(2) + (nx*(subson2(5) - subson3(5)))/(sps*subson2(1))
            subson3(3) = subson2(3) + (ny*(subson2(5) - subson3(5)))/(sps*subson2(1))
            subson3(4) = subson2(4) + (nz*(subson2(5) - subson3(5)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            rightv(4) = subson3(4)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
            ikins = subson3(5)/((gamma - 1.0d0)*(subson3(1)))
            rightv(5) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
        end if

      case (9)
        !outlets subsonic or supersonic will be chosen based on mach number
        ! if (boundtype.eq.0)then
        ! rightv(1:nof_variables)=leftv(1:nof_variables)
        !else
        rightv(1:nof_variables) = outflow2(initcond, pox, poy, poz)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
        subson1(1:nof_variables) = rightv(1:nof_variables)
        subson2(1:nof_variables) = leftv(1:nof_variables)
        sps = sqrt((gamma*subson2(5))/(subson2(1)))
        vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
        sps = sqrt((gamma*subson2(5))/(subson2(1)))
        call prim2cons2(n, leftv, rightv)
        if (vel/(sps + tolsmall) .gt. 1.0d0) then        !supersonic
          call prim2cons2(n, leftv, rightv)
          rightv(1:nof_variables) = leftv(1:nof_variables)
        else
          subson3(5) = subson1(5)
          subson3(1) = subson2(1) + (subson3(5) - subson2(5))/(sps**2)
          subson3(2) = subson2(2) + (nx*(subson2(5) - subson3(5)))/(sps*subson2(1))
          subson3(3) = subson2(3) + (ny*(subson2(5) - subson3(5)))/(sps*subson2(1))
          subson3(4) = subson2(4) + (nz*(subson2(5) - subson3(5)))/(sps*subson2(1))
          rightv(1) = subson3(1)
          rightv(2) = subson3(2)*subson3(1)
          rightv(3) = subson3(3)*subson3(1)
          rightv(4) = subson3(4)*subson3(1)
          skins = oo2*((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
          ikins = subson3(5)/((gamma - 1.0d0)*(subson3(1)))
          rightv(5) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          !end if
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
        end if

      case (99)    !bleed boundary
        rightv(1:nof_variables) = bleed3d(iconsidered, facex, pox, poy, poz)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
        end if

      case (3)!symmetry
        call rotatef(n, cleft_rot, leftv, angle1, angle2)
        if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
          cright_rot(1) = cleft_rot(1)
          cright_rot(2) = -(cleft_rot(2)) + 2.0d0*cleft_rot(1)*srf_speedrot(2)
          cright_rot(3) = cleft_rot(3)
          cright_rot(4) = cleft_rot(4)
          cright_rot(5) = cleft_rot(5) + 2.0d0*cleft_rot(1)*(srf_speedrot(2)**2) - 2.0d0*cleft_rot(2)*srf_speedrot(2)
        else
          cright_rot(1) = cleft_rot(1)
          cright_rot(2) = -cleft_rot(2)
          cright_rot(3) = cleft_rot(3)
          cright_rot(4) = cleft_rot(4)
          cright_rot(5) = cleft_rot(5)
          if (multispecies .eq. 1) then
            cright_rot(6) = cleft_rot(6)
            cright_rot(7) = cleft_rot(7)
            cright_rot(8) = cleft_rot(8)
          end if
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
          if (passivescalar .gt. 0) then
            cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
            cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)
          end if
        end if
        call rotateb(n, rightv, cright_rot, angle1, angle2)

      case (4)  ! wall
        if (itestcase .eq. 3) then
          call rotatef(n, cleft_rot, leftv, angle1, angle2)
          if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
            cright_rot(1) = cleft_rot(1)
            cright_rot(2) = -(cleft_rot(2)) + 2.0d0*cleft_rot(1)*srf_speedrot(2)
            cright_rot(3) = cleft_rot(3)
            cright_rot(4) = cleft_rot(4)
            cright_rot(5) = cleft_rot(5) + cleft_rot(1)*(srf_speedrot(2)**2)*2.0d0 - 2.0d0*cleft_rot(2)*srf_speedrot(2)
          else
            cright_rot(:) = cleft_rot(:)
            cright_rot(2) = -cleft_rot(2)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cturbr(:) = cturbl(:)
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
              cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)
            end if
          end if
          call rotateb(n, rightv, cright_rot, angle1, angle2)
        else
          if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
            rightv(1) = leftv(1)
            rightv(2) = -leftv(2) + 2.0d0*leftv(1)*srf_speed(2)
            rightv(3) = -leftv(3) + 2.0d0*leftv(1)*srf_speed(3)
            rightv(4) = -leftv(4) + 2.0d0*leftv(1)*srf_speed(4)
            rightv(5) = leftv(5) + 2.0d0*leftv(1)*(srf_speed(2)**2 + srf_speed(3)**2 + srf_speed(4)**2) &- 2.0d0*(leftv(2)  *srf_speed(2) + leftv(3)*srf_speed(3) + leftv(4)*srf_speed(4))
          else
            rightv(1) = leftv(1)
            rightv(2) = -leftv(2)
            rightv(3) = -leftv(3)
            rightv(4) = -leftv(4)
            rightv(5) = leftv(5)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulencemodel .ne. 2) then
              cturbr(:) = -cturbl(:)
              if (passivescalar .gt. 0) then
                cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &-cturbl(turbulenceequations +   1:turbulenceequations + passivescalar)
              end if
            else
              cturbr(1) = -cturbl(1)
              cturbr(2) = 60.0d0*visc/(beta_i1*(ielem(n, iconsidered)%walldist**2))
              if (passivescalar .gt. 0) then
                cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &-cturbl(turbulenceequations +   1:turbulenceequations + passivescalar)
              end if
            end if
          end if
        end if

      case (6)!farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
        call rotatef(n, cleft_rot, leftv, angle1, angle2)
        vnb = cleft_rot(2)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
        subson1(1:nof_variables) = rightv(1:nof_variables)
        subson2(1:nof_variables) = leftv(1:nof_variables)
        sps = sqrt((gamma*subson2(5))/(subson2(1)))
        vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
        call prim2cons2(n, leftv, rightv)
        if (vnb .le. 0.0d0) then                !inflow
          ibfc = -1
          if ((abs(vnb)) .ge. sps) then
            !supersonic
            rightv = inflow(initcond, pox, poy, poz)
          else
            !subsonic
            rightv = inflow(initcond, pox, poy, poz)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
            subson1(1:nof_variables) = rightv(1:nof_variables)
            subson2(1:nof_variables) = leftv(1:nof_variables)
            sps = sqrt((gamma*subson2(5))/(subson2(1)))
            vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
            call prim2cons2(n, leftv, rightv)
            subson3(5) = 0.5*((subson1(5)) + (subson2(5)) - (subson2(1)*sps*((nx*(subson1(2) - subson2(2))) + (ny*(subson1  (3) - subson2(3))) & + (nz*(subson1(4) - subson2(4))))))
            subson3(1) = subson1(1) + (subson3(5) - subson1(5))/(sps**2)
            subson3(2) = subson1(2) - (nx*(subson1(5) - subson3(5)))/(sps*subson2(1))
            subson3(3) = subson1(3) - (ny*(subson1(5) - subson3(5)))/(sps*subson2(1))
            subson3(4) = subson1(4) - (nz*(subson1(5) - subson3(5)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            rightv(4) = subson3(4)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
            ikins = subson3(5)/((gamma - 1.0d0)*(subson3(1)))
            rightv(5) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulencemodel .eq. 1) then
              cturbr(1) = visc*turbinit
            end if
            if (turbulencemodel .eq. 2) then
              cturbr(1) = (1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
              cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
              if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
                cturbr(1) = (1.5d0*i_turb_inlet*(kinit_srf**2))*rightv(1)!k initialization
                cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
              end if
            end if
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = pass_inlet(initcond, pox, poy, poz)  * rightv(1)
            end if
          end if
        else
          ibfc = -2
          if ((abs(vnb)) .ge. sps) then
            rightv(1:nof_variables) = leftv(1:nof_variables)
          else
            rightv(1:nof_variables) = outflow(initcond, pox, poy, poz)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
            subson1(1:nof_variables) = rightv(1:nof_variables)
            subson2(1:nof_variables) = leftv(1:nof_variables)
            call prim2cons2(n, leftv, rightv)
            subson3(5) = subson1(5)
            subson3(1) = subson2(1) + (subson3(5) - subson2(5))/(sps**2)
            subson3(2) = subson2(2) + (nx*(subson2(5) - subson3(5)))/(sps*subson2(1))
            subson3(3) = subson2(3) + (ny*(subson2(5) - subson3(5)))/(sps*subson2(1))
            subson3(4) = subson2(4) + (nz*(subson2(5) - subson3(5)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            rightv(4) = subson3(4)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
            ikins = subson3(5)/((gamma - 1.0d0)*(subson3(1)))
            rightv(5) = (subson3(1)*(ikins)) + (subson3(1)*skins)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              cturbr(:) = cturbl(:)
            end if
          end if
        end if
    end select
  end subroutine boundarys

  subroutine boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,  cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
    ! @brief
    ! this subroutine applies the boundary condition to each bounded cell
    implicit none
    integer, intent(in)::n, b_code, iconsidered, facex
    integer, intent(inout)::ibfc
    real, dimension(1:nof_variables), intent(inout)::leftv, rightv
    real, dimension(1:nof_variables), intent(in)::srf_speedrot, srf_speed
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real, intent(in)::angle1, angle2, nx, ny, nz
    real, dimension(turbulenceequations), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::cright_rot, cleft_rot
    real, dimension(1:nof_variables)::subson1, subson2, subson3
    real::mp_pinfl, mp_pinfr, gammal, gammar
    real::sps, skins, ikins, vel, vnb, theeta, reeta
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx, vhx, amp, dvel, rgg, tt1
    select case (b_code)
      case (1)!inflow subsonic or supersonic will be chosen based on mach number
        if (boundtype .eq. 0) then        !supersonic
          rightv(1:nof_variables) = inflow2d(initcond, pox, poy)
        else                !subsonic
          rightv(1:nof_variables) = inflow2d(initcond, pox, poy)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          subson1(1:nof_variables) = rightv(1:nof_variables)
          subson2(1:nof_variables) = leftv(1:nof_variables)
          sps = sqrt((gamma*subson2(4))/(subson2(1)))
          vel = sqrt(subson2(2)**2 + subson2(3)**2)
          call prim2cons2(n, leftv, rightv)
          if (vel/(sps + tolsmall) .gt. 1.0d0) then        !supersonic
            rightv(1:nof_variables) = inflow2d(initcond, pox, poy)
          else                !subsonic
            subson3(4) = 0.5*((subson1(4)) + (subson2(4)) - (subson2(1)*sps*((nx*(subson1(2) - subson2(2))) + (ny*(subson1(3) - subson2(3))) &)))
            subson3(1) = subson1(1) + (subson3(4) - subson1(4))/(sps**2)
            subson3(2) = subson1(2) - (nx*(subson1(4) - subson3(4)))/(sps*subson2(1))
            subson3(3) = subson1(3) - (ny*(subson1(4) - subson3(4)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2))
            ikins = subson3(4)/((gamma - 1.0d0)*(subson3(1)))
            rightv(4) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (turbulencemodel .eq. 1) then
            cturbr(1) = visc*turbinit
          end if
          if (turbulencemodel .eq. 2) then
            cturbr(1) = (1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
            cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
          end if
          if (passivescalar .gt. 0) then
            cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = pass_inlet2d(initcond, pox, poy)* rightv(1)
          end if
        end if

      case (2)!outflow subsonic or supersonic will be chosen based on mach number
        if (boundtype .eq. 0) then
          rightv(1:nof_variables) = leftv(1:nof_variables)
        else
          rightv(1:nof_variables) = outflow2d(initcond, pox, poy)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          subson1(1:nof_variables) = rightv(1:nof_variables)
          subson2(1:nof_variables) = leftv(1:nof_variables)
          sps = sqrt((gamma*subson2(4))/(subson2(1)))
          vel = sqrt(subson2(2)**2 + subson2(3)**2)
          sps = sqrt((gamma*subson2(4))/(subson2(1)))
          call prim2cons2(n, leftv, rightv)
          if (vel/(sps + tolsmall) .gt. 1.0d0) then        !supersonic
            rightv(1:nof_variables) = leftv(1:nof_variables)
          else
            subson3(4) = subson1(4)
            subson3(1) = subson2(1) + (subson3(4) - subson2(4))/(sps**2)
            subson3(2) = subson2(2) + (nx*(subson2(4) - subson3(4)))/(sps*subson2(1))
            subson3(3) = subson2(3) + (ny*(subson2(4) - subson3(4)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2))
            ikins = subson3(4)/((gamma - 1.0d0)*(subson3(1)))
            rightv(4) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
        end if

      case (99)  !bleed
        rightv(1:nof_variables) = bleed2d(iconsidered, facex, pox, poy)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          cturbr(:) = cturbl(:)
        end if

      case (3)!symmetry
        if ((initcond .eq. 102) .or. (initcond .eq. 30) .or. (initcond .eq. 222)) then        !shock density interaction
          if ((initcond .eq. 102)) then        !shock density interaction
            if (pox(1) .lt. ((1.0d0/6.0d0) + ((1.0d0 + 20.0d0*t)/(sqrt(3.0d0))))) then
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
            !internal energy
            ie1 = ((p1)/((gamma - 1.0d0)*r1))
            !total energy
            e1 = (p1/(gamma - 1)) + (r1*skin1)
            !vector of conserved variables now
            rightv(1) = r1
            rightv(2) = r1*u1
            rightv(3) = r1*v1
            rightv(4) = e1
          end if
          if ((initcond .eq. 222)) then
            if (sqrt((pox(1)**2) + (poy(1)**2)) .lt. (t/3.0d0)) then
              r1 = 16.0d0
              u1 = 0.0
              v1 = 0.0
              p1 = 16.0d0/3.0d0
            else
              r1 = 1.0d0 + (t/sqrt((pox(1)**2) + (poy(1)**2)))
              reeta = -1
              theeta = atan(poy(1)/pox(1))
              u1 = reeta*cos(theeta)
              v1 = reeta*sin(theeta)
              p1 = 1.0e-6
            end if
            skin1 = (oo2)*((u1**2) + (v1**2))
            !internal energy
            ie1 = ((p1)/((gamma - 1.0d0)*r1))
            !total energy
            e1 = (p1/(gamma - 1)) + (r1*skin1)
            !vector of conserved variables now
            rightv(1) = r1
            rightv(2) = r1*u1
            rightv(3) = r1*v1
            rightv(4) = e1
          end if
          if ((initcond .eq. 30)) then
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
            !internal energy
            ie1 = ((p1)/((gamma - 1.0d0)*r1))
            !total energy
            e1 = (p1/(gamma - 1)) + (r1*skin1)
            !vector of conserved variables now
            ! rightv(1)=leftv(1)
            !rightv(2)=0.0
            !rightv(3)=0.0
            !rightv(4)=leftv(4)
            rightv(1) = r1
            rightv(2) = r1*u1
            rightv(3) = r1*v1
            rightv(4) = e1
          end if
        else
          call rotatef2d(n, cleft_rot, leftv, angle1, angle2)
          cright_rot(1) = cleft_rot(1)
          cright_rot(2) = -cleft_rot(2)
          cright_rot(3) = cleft_rot(3)
          cright_rot(4) = cleft_rot(4)
          if (multispecies .eq. 1) then
            cright_rot(5) = cleft_rot(5)
            cright_rot(6) = cleft_rot(6)
            cright_rot(7) = cleft_rot(7)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cturbr(:) = cturbl(:)
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
              cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)
            end if
          end if
          call rotateb2d(n, rightv, cright_rot, angle1, angle2)
        end if

      case (4)!wall
        if (itestcase .eq. 3) then
          call rotatef2d(n, cleft_rot, leftv, angle1, angle2)
          if (governingequations .eq. -1) then
            cright_rot(:) = cleft_rot(:)
            cright_rot(2) = -cleft_rot(2)
          else
            cright_rot(1) = cleft_rot(1)
            cright_rot(2) = -cleft_rot(2)
            cright_rot(3) = cleft_rot(3)
            cright_rot(4) = cleft_rot(4)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cturbr(:) = cturbl(:)
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
              cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)
            end if
          end if
          call rotateb2d(n, rightv, cright_rot, angle1, angle2)
        else
          rightv(1) = leftv(1)
          rightv(2) = -leftv(2)
          rightv(3) = -leftv(3)
          rightv(4) = leftv(4)
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulencemodel .ne. 2) then
              cturbr(:) = -cturbl(:)
              if (passivescalar .gt. 0) then
                cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = & -cturbl(turbulenceequations +   1:turbulenceequations + passivescalar)
              end if
            else
              cturbr(1) = -cturbl(1)
              cturbr(2) = 60.0d0*visc/(beta_i1*(ielem(n, iconsidered)%walldist**2))
              if (passivescalar .gt. 0) then
                cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                  -cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)
              end if
            end if
          end if
        end if

      case (6)!farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
        call rotatef2d(n, cleft_rot, leftv, angle1, angle2)
        vnb = cleft_rot(2)/cleft_rot(1)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
        subson1(1:nof_variables) = rightv(1:nof_variables)
        subson2(1:nof_variables) = leftv(1:nof_variables)
        sps = sqrt((gamma*subson2(4))/(subson2(1)))
        vel = sqrt(subson2(2)**2 + subson2(3)**2)
        call prim2cons2(n, leftv, rightv)
        if (vnb .le. 0.0d0) then                !inflow
          ibfc = -1
          if ((abs(vnb)) .ge. sps) then
            !supersonic
            rightv(1:nof_variables) = inflow2d(initcond, pox, poy)
          else
            !subsonic
            rightv(1:nof_variables) = inflow2d(initcond, pox, poy)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
            subson1(1:nof_variables) = rightv(1:nof_variables)
            subson2(1:nof_variables) = leftv(1:nof_variables)
            sps = sqrt((gamma*subson2(4))/(subson2(1)))
            vel = sqrt(subson2(2)**2 + subson2(3)**2)
            call prim2cons2(n, leftv, rightv)
            subson3(4)=0.5d0*((subson1(4))+(subson2(4))-(subson2(1)*sps*((nx*(subson1(2)-subson2(2)))+(ny*(subson1(3)  -subson2(3))))))
            subson3(1) = subson1(1) + (subson3(4) - subson1(4))/(sps**2)
            subson3(2) = subson1(2) - (nx*(subson1(4) - subson3(4)))/(sps*subson2(1))
            subson3(3) = subson1(3) - (ny*(subson1(4) - subson3(4)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2))
            ikins = subson3(4)/((gamma - 1.0d0)*(subson3(1)))
            rightv(4) = (subson3(1)*(ikins)) + (subson3(1)*skins)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulencemodel .eq. 1) then
              cturbr(1) = visc*turbinit
            end if
            if (turbulencemodel .eq. 2) then
              cturbr(1) = (1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
              cturbr(2) = rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
            end if
            if (passivescalar .gt. 0) then
              cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = pass_inlet2d(initcond, pox, poy)  *rightv(1)
            end if
          end if
        else
          ibfc = -2
          if ((abs(vnb)) .ge. sps) then
            rightv(1:nof_variables) = leftv(1:nof_variables)
          else
            rightv(1:nof_variables) = outflow2d(initcond, pox, poy)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
            subson1(1:nof_variables) = rightv(1:nof_variables)
            subson2(1:nof_variables) = leftv(1:nof_variables)
            call prim2cons2(n, leftv, rightv)
            subson3(4) = subson1(4)
            subson3(1) = subson2(1) + (subson3(4) - subson2(4))/(sps**2)
            subson3(2) = subson2(2) + (nx*(subson2(4) - subson3(4)))/(sps*subson2(1))
            subson3(3) = subson2(3) + (ny*(subson2(4) - subson3(4)))/(sps*subson2(1))
            rightv(1) = subson3(1)
            rightv(2) = subson3(2)*subson3(1)
            rightv(3) = subson3(3)*subson3(1)
            skins = oo2*((subson3(2)**2) + (subson3(3)**2))
            ikins = subson3(4)/((gamma - 1.0d0)*(subson3(1)))
            rightv(4) = (subson3(1)*(ikins)) + (subson3(1)*skins)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              cturbr(:) = cturbl(:)
            end if
          end if
        end if
      end select
  end subroutine boundarys2d

  subroutine compute_eigenvectors(n, rveigl, rveigr, eigvl, eigvr, gamma)
    ! @brief
    ! this subroutine computes the left and right eigenvectors
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables), intent(in)::rveigl, rveigr
    real, intent(in)::gamma
    real, dimension(1:nof_variables, 1:nof_variables), intent(inout)::eigvl, eigvr
    real::rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, oor1, oor2
    integer::ivgt
    eigvr = zero
    gammam1 = gamma - 1.0d0
    oor1 = 1.0d0/rveigl(1)
    oor2 = 1.0d0/rveigr(1)
    rs = oo2*(rveigl(1) + rveigr(1))
    us = oo2*((rveigl(2)*oor1) + (rveigr(2)*oor2))
    vs = oo2*((rveigl(3)*oor1) + (rveigr(3)*oor2))
    ws = oo2*((rveigl(4)*oor1) + (rveigr(4)*oor2))
    es = oo2*(rveigl(5) + rveigr(5))
    vvs = (us**2) + (vs**2) + (ws**2)
    vsd = oo2*vvs
    ps = (gamma - 1.0d0)*(es - oo2*rs*vvs)
    as = sqrt(gamma*ps/rs)
    hs = (oo2*vvs) + ((as**2)/gammam1)
    eigvr(1, 1) = 1.0d0; eigvr(1, 2) = 1.0d0; eigvr(1, 3) = 0.0d0; eigvr(1, 4) = 0.0d0; eigvr(1, 5) = 1.0d0
    eigvr(2, 1) = us - as; eigvr(2, 2) = us; eigvr(2, 3) = 0.0d0; eigvr(2, 4) = 0.0d0; eigvr(2, 5) = us + as
    eigvr(3, 1) = vs; eigvr(3, 2) = vs; eigvr(3, 3) = 1.0d0; eigvr(3, 4) = 0.0d0; eigvr(3, 5) = vs
    eigvr(4, 1) = ws; eigvr(4, 2) = ws; eigvr(4, 3) = 0.0d0; eigvr(4, 4) = 1.0d0; eigvr(4, 5) = ws
    eigvr(5, 1) = hs - (us*as); eigvr(5, 2) = oo2*vvs; eigvr(5, 3) = vs; eigvr(5, 4) = ws; eigvr(5, 5) = hs + (us*as)
    eigvl(1, 1) = (hs*us + as*us**2 + as*vs**2 - as*vsd - us*vsd + as*ws**2)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(1, 2) = (-hs - as*us + vsd)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(1, 3) = -((as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd))
    eigvl(1, 4) = -((as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd))
    eigvl(1, 5) = as/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(2, 1) = (2.0*as*hs - 2.0*as*us**2 - 2.0d0*as*vs**2 - 2.0d0*as*ws**2)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(2, 2) = (2.0d0*as*us)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(2, 3) = (2.0d0*as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(2, 4) = (2.0d0*as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(2, 5) = (-2.0d0*as)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(3, 1) = (-2.0*as*hs*vs + 2.0*as*vs*vsd)/(2.0*as*hs - 2.0*as*vsd)
    eigvl(3, 2) = 0.0d0
    eigvl(3, 3) = 1.0d0
    eigvl(3, 4) = 0.0d0
    eigvl(3, 5) = 0.0d0
    eigvl(4, 1) = (-2.0d0*as*hs*ws + 2.0d0*as*vsd*ws)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(4, 2) = 0.0d0
    eigvl(4, 3) = 0.0d0
    eigvl(4, 4) = 1.0d0
    eigvl(4, 5) = 0.0d0
    eigvl(5, 1) = (-(hs*us) + as*us**2 + as*vs**2 - as*vsd + us*vsd + as*ws**2)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(5, 2) = (hs - as*us - vsd)/(2.0d0*as*hs - 2.0d0*as*vsd)
    eigvl(5, 3) = -((as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd))
    eigvl(5, 4) = -((as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd))
    eigvl(5, 5) = as/(2.0d0*as*hs - 2.0d0*as*vsd)
  end subroutine compute_eigenvectors

  subroutine compute_jacobianse(n, iconsidered, eigvl, rveigl, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
    ! @brief
    ! this subroutine computes the jacobians for the implicit time stepping
    implicit none
    integer, intent(in)::n, iconsidered
    real, intent(in)::angle1, angle2
    real, dimension(1:nof_variables), intent(in)::rveigl, srf_speedrot
    real, intent(in)::gamma
    real, dimension(1:nof_variables, 1:nof_variables), intent(inout)::eigvl
    real::rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, phi, a1, a2, a3, oors, nx, ny, nz
    integer::ivgt
    a2 = gamma - 1.0d0
    a3 = gamma - 2.0d0
    oors = 1.0d0/rveigl(1)
    rs = (rveigl(1))
    us = (rveigl(2)*oors)
    vs = (rveigl(3)*oors)
    ws = (rveigl(4)*oors)
    es = (rveigl(5)*oors)
    phi = oo2*(a2)*((us*us) + (vs*vs) + (ws*ws))
    a1 = gamma*es - phi
    vvs = nx*us + ny*vs + nz*ws
    if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
      eigvl(1, 1) = 0.0d0 - srf_speedrot(2); eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = nz; eigvl(1, 5) = 0.0d0
      eigvl(2, 1) = nx*phi - us*vvs; eigvl(2, 2) = vvs - a3*nx*us - srf_speedrot(2); eigvl(2, 3) = ny*us - a2*nx*vs; eigvl(2, 4) = nz*us - a2*nx*ws; eigvl(2, 5) = a2*nx
      eigvl(3, 1) = ny*phi - vs*vvs; eigvl(3, 2) = nx*vs - a2*ny*us; eigvl(3, 3) = vvs - a3*ny*vs - srf_speedrot(2); eigvl(3, 4) = nz*vs - a2*ny*ws; eigvl(3, 5) = a2*ny
      eigvl(4, 1) = nz*phi - ws*vvs; eigvl(4, 2) = nx*ws - a2*nz*us; eigvl(4, 3) = ny*ws - a2*nz*vs; eigvl(4, 4) = vvs - a3*nz*ws - srf_speedrot(2); eigvl(4, 5) = a2*nz
      eigvl(5, 1) = vvs*(phi - a1); eigvl(5, 2) = nx*a1 - a2*us*vvs; eigvl(5, 3) = ny*a1 - a3*vs*vvs; eigvl(5, 4) = nz*a1 - a2*ws*vvs; eigvl(5, 5) = gamma*vvs - srf_speedrot(2)
    else
      eigvl(1, 1) = 0.0d0; eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = nz; eigvl(1, 5) = 0.0d0
      eigvl(2, 1) = nx*phi - us*vvs; eigvl(2, 2) = vvs - a3*nx*us; eigvl(2, 3) = ny*us - a2*nx*vs; eigvl(2, 4) = nz*us - a2*nx*ws; eigvl(2, 5) = a2*nx
      eigvl(3, 1) = ny*phi - vs*vvs; eigvl(3, 2) = nx*vs - a2*ny*us; eigvl(3, 3) = vvs - a3*ny*vs; eigvl(3, 4) = nz*vs - a2*ny*ws; eigvl(3, 5) = a2*ny
      eigvl(4, 1) = nz*phi - ws*vvs; eigvl(4, 2) = nx*ws - a2*nz*us; eigvl(4, 3) = ny*ws - a2*nz*vs; eigvl(4, 4) = vvs - a3*nz*ws; eigvl(4, 5) = a2*nz
      eigvl(5, 1) = vvs*(phi - a1); eigvl(5, 2) = nx*a1 - a2*us*vvs; eigvl(5, 3) = ny*a1 - a3*vs*vvs; eigvl(5, 4) = nz*a1 - a2*ws*vvs; eigvl(5, 5) = gamma*vvs
    end if
  end subroutine compute_jacobianse

  subroutine compute_eigenvectors2d(n, rveigl, rveigr, eigvl, eigvr, gamma)
    ! @brief
    ! this subroutine computes the left and right eigenvectors  in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables), intent(in)::rveigl, rveigr
    real, intent(in)::gamma
    real, dimension(1:nof_variables, 1:nof_variables), intent(inout)::eigvl, eigvr
    real::rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, g8, s1, s2, vtots, oor1, oor2
    integer::ivgt, j, k
    eigvr = zero
    gammam1 = gamma - 1.0d0
    oor1 = 1.0d0/rveigl(1)
    oor2 = 1.0d0/rveigr(1)
    eigvr = 0.0d0
    gammam1 = gamma - 1.0d0
    rs = 0.5*(rveigl(1) + rveigr(1))
    us = 0.5*((rveigl(2)*oor1) + (rveigr(2)*oor2))
    vs = 0.5*((rveigl(3)*oor1) + (rveigr(3)*oor2))
    es = 0.5*(rveigl(4) + rveigr(4))
    g8 = gamma - 1.0d0
    vtots = us**2 + vs**2
    ps = (gamma - 1)*(es - oo2*rs*vtots)
    as = sqrt(gamma*ps/rs)
    hs = oo2*vtots + (as**2)/(g8)
    eigvr(1, 1) = 1.d0; eigvr(1, 2) = 1.d0; eigvr(1, 3) = 0.d0; eigvr(1, 4) = 1.d0
    eigvr(2, 1) = us - as; eigvr(2, 2) = us; eigvr(2, 3) = 0.d0; eigvr(2, 4) = us + as
    eigvr(3, 1) = vs; eigvr(3, 2) = vs; eigvr(3, 3) = 1.d0; eigvr(3, 4) = vs
    eigvr(4, 1) = hs - us*as; eigvr(4, 2) = oo2*vtots; eigvr(4, 3) = vs; eigvr(4, 4) = hs + us*as
    s1 = as/(gamma - 1d0)
    s2 = as**2/(gamma - 1d0)
    eigvl(1, 1) = hs + s1*(us - as); eigvl(1, 2) = -(us + s1); eigvl(1, 3) = -vs; eigvl(1, 4) = 1.d0
    eigvl(2, 1) = -2d0*hs + 4d0*s2; eigvl(2, 2) = 2d0*us; eigvl(2, 3) = 2d0*vs; eigvl(2, 4) = -2.d0
    eigvl(3, 1) = -2d0*vs*s2; eigvl(3, 2) = 0.d0; eigvl(3, 3) = 2d0*s2; eigvl(3, 4) = 0.d0
    eigvl(4, 1) = hs - s1*(us + as); eigvl(4, 2) = -us + s1; eigvl(4, 3) = -vs; eigvl(4, 4) = 1.d0
    do j = 1, 4
      do k = 1, 4
        eigvl(j, k) = eigvl(j, k)/(2d0*s2)
      end do
    end do
  end subroutine compute_eigenvectors2d

  subroutine compute_jacobianse2d(n, eigvl, rveigl, gamma, angle1, angle2, nx, ny, nz)
    ! @brief
    ! this subroutine computes the jacobians for the implicit time stepping in 2d
    implicit none
    integer, intent(in)::n
    real, intent(in)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables), intent(in)::rveigl
    real, intent(in)::gamma
    real, dimension(1:nof_variables, 1:nof_variables), intent(inout)::eigvl
    real::rs, us, vs, es, ps, vvs, as, hs, gammam1, vsd, phi, a1, a2, a3, oors
    integer::ivgt
    a2 = gamma - 1.0d0
    a3 = gamma - 2.0d0
    oors = 1.0d0/rveigl(1)
    rs = (rveigl(1))
    us = (rveigl(2)*oors)
    vs = (rveigl(3)*oors)
    es = (rveigl(4)*oors)
    phi = oo2*(a2)*((us*us) + (vs*vs))
    a1 = gamma*es - phi
    vvs = nx*us + ny*vs
    eigvl(1, 1) = 0.0d0; eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = 0.0d0
    eigvl(2, 1) = nx*phi - us*vvs; eigvl(2, 2) = vvs - a3*nx*us; eigvl(2, 3) = ny*us - a2*nx*vs; eigvl(2, 4) = a2*nx
    eigvl(3, 1) = ny*phi - vs*vvs; eigvl(3, 2) = nx*vs - a2*ny*us; eigvl(3, 3) = vvs - a3*ny*vs; eigvl(3, 4) = a2*ny
    eigvl(4, 1) = vvs*(phi - a1); eigvl(4, 2) = nx*a1 - a2*us*vvs; eigvl(4, 3) = ny*a1 - a2*vs*vvs; eigvl(4, 4) = gamma*vvs
  end subroutine compute_jacobianse2d

  subroutine eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
    ! @brief
    ! this subroutine computes the tubulent eddy viscosity for turbulence models
    implicit none
    real, dimension(1:2), intent(inout)::turbmv
    real, dimension(1:nof_variables), intent(in)::leftv, rightv
    real, dimension(1), intent(inout)::etvm
    real, dimension(1:4), intent(inout)::viscl, laml
    real, dimension(1:20), intent(inout)::eddyfl, eddyfr
    integer, intent(in)::n
    real::ml, smmm, mtt, chi, tolepsma, chipow3, fv1
    integer:: ihgt, ihgj
    real, dimension(3, 3)::vortet, tvort, svort
    real:: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real:: wally, d_omplus, phi_2, phi_1, f_1, f_2, k_0, om_0, alpha_inf, alpha_star, mu_turb
    real:: sigma_k_l, sigma_om_l, sigma_k_r, sigma_om_r
    real::snorm, dervk_dervom, re_t_sst, rho_0, beta_i
    tolepsma = 10e-16
    if (turbulence .eq. 0) then
      viscl(4) = zero; viscl(3) = zero
      laml(4) = zero; laml(3) = zero
    else
      select case (turbulencemodel)
        case (1)
          turbmv(1) = eddyfl(2)
          turbmv(2) = eddyfr(2)
          chi = abs(max(turbmv(1), tolepsma)/max(viscl(1), tolepsma))
          chipow3 = chi*chi*chi
          fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
          viscl(3) = turbmv(1)*fv1
          chi = abs(max(turbmv(2), tolepsma)/max(viscl(2), tolepsma))
          chipow3 = chi*chi*chi
          fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
          viscl(4) = turbmv(2)*fv1
        case (2)
          vortet(1, 1:3) = eddyfl(4:6)
          vortet(2, 1:3) = eddyfl(7:9)
          vortet(3, 1:3) = eddyfl(10:12)
          ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
          vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
          wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)
          do ihgt = 1, 3
            do ihgj = 1, 3
              tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
            end do
          end do
          svort = 0.5*(vortet + tvort)
          snorm = sqrt(2.0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + (svort(1, 3)*svort(1, 3)) + &
                            (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)) + (svort(2, 3)*svort(2, 3)) + &
                            (svort(3, 1)*svort(3, 1)) + (svort(3, 2)*svort(3, 2)) + (svort(3, 3)*svort(3, 3))))
          wally = eddyfl(1)
          rho_0 = leftv(1)
          k_0 = max(tolepsma, eddyfl(2)/leftv(1))
          om_0 = max(eddyfl(3)/leftv(1), ufreestream/charlength/10.0)
          dervk_dervom = (eddyfl(13)*eddyfl(16)) + (eddyfl(14)*eddyfl(17)) + (eddyfl(15)*eddyfl(18))
          d_omplus = max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)
          phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(1)/(rho_0*wally*wally*om_0))
          phi_1 = min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally))
          f_1 = tanh(phi_1**4)
          f_2 = tanh(phi_2**2)
          re_t_sst = rho_0*k_0/(viscl(1)*om_0)
          beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
          alpha_star0 = beta_i/3.0
          alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
          alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
          viscl(3) = rho_0*k_0/om_0/max(1.0/alpha_star, snorm*f_2/(aa_1*om_0))
          sigma_k_l = sigma_k1/f_1 + sigma_k2/f_2
          sigma_om_l = sigma_om1/f_1 + sigma_om2/f_2
          if (eddyfr(1) .gt. 0.0) then
            vortet(1, 1:3) = eddyfr(4:6)
            vortet(2, 1:3) = eddyfr(7:9)
            vortet(3, 1:3) = eddyfr(10:12)
            ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
            vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
            wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)
            do ihgt = 1, 3
              do ihgj = 1, 3
                tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
              end do
            end do
            svort = 0.5*(vortet + tvort)
            snorm = sqrt(2.0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + (svort(1, 3)*svort(1, 3)) + &
                              (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)) + (svort(2, 3)*svort(2, 3)) + &
                              (svort(3, 1)*svort(3, 1)) + (svort(3, 2)*svort(3, 2)) + (svort(3, 3)*svort(3, 3))))
            wally = eddyfr(1)
            rho_0 = rightv(1)
            k_0 = eddyfr(2)/rightv(1)
            om_0 = max(eddyfr(3)/rightv(1), 1.0e-6)
            dervk_dervom = (eddyfr(13)*eddyfr(16)) + (eddyfr(14)*eddyfr(17)) + (eddyfr(15)*eddyfr(18))
            d_omplus = max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !i need derivative of k
            phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(2)/(rho_0*wally*wally*om_0))
            phi_1 = min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally))
            f_1 = tanh(phi_1**4)
            f_2 = tanh(phi_2**2)
            re_t_sst = rho_0*k_0/(viscl(2)*om_0)
            beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
            alpha_star0 = beta_i/3.0
            alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
            alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
            viscl(4) = rho_0*k_0/om_0/max(1.0/alpha_star, snorm*f_2/(aa_1*om_0))
            sigma_k_r = sigma_k1/f_1 + sigma_k2/f_2
            sigma_om_r = sigma_om1/f_1 + sigma_om2/f_2
          else
            viscl(4) = -viscl(3)
            sigma_k_r = sigma_k_l
            sigma_om_r = sigma_om_l
          end if
      end select
      viscl(3) = min(10000000*visc, viscl(3))
      viscl(4) = min(10000000*visc, viscl(4))
      laml(3) = (viscl(3)*gamma/(prtu*(gamma - 1))) + (viscl(1)*gamma/(prandtl*(gamma - 1)))
      laml(4) = (viscl(4)*gamma/(prtu*(gamma - 1))) + (viscl(2)*gamma/(prandtl*(gamma - 1)))
      viscl(3) = max(0.0d0, viscl(3))
      viscl(4) = max(0.0d0, viscl(4))
      if ((turbmv(1) .lt. zero) .or. (turbmv(2) .lt. zero)) then
        viscl(3) = 0.0d0
        viscl(4) = 0.0d0
      end if
      etvm(1) = (0.5*(viscl(1) + viscl(2))) + (0.5*(viscl(3) + viscl(4)))
      if (turbulencemodel .eq. 2) then
        eddyfl(19) = viscl(1) + viscl(3)/sigma_k_l
        eddyfr(19) = viscl(2) + viscl(4)/sigma_k_r
        eddyfl(20) = viscl(1) + viscl(3)/sigma_om_l
        eddyfr(20) = viscl(2) + viscl(4)/sigma_om_r
      end if
    end if
  end subroutine eddyvisco

  subroutine eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
    ! @brief
    ! this subroutine computes the tubulent eddy viscosity for turbulence models
    implicit none
    real, dimension(1:2), intent(inout)::turbmv
    real, dimension(1:nof_variables), intent(in)::leftv, rightv
    real, dimension(1), intent(inout)::etvm
    real, dimension(1:4), intent(inout)::viscl, laml
    real, dimension(1:20), intent(inout)::eddyfl, eddyfr
    integer, intent(in)::n
    real::ml, smmm, mtt, chi, tolepsma, chipow3, fv1
    integer:: ihgt, ihgj
    real, dimension(2, 2)::vortet, tvort, svort
    real:: ux, uy, uz, vx, vy, vz, wx, wy, wz
    real:: wally, d_omplus, phi_2, phi_1, f_1, f_2, k_0, om_0, alpha_inf, alpha_star, mu_turb
    real:: sigma_k_l, sigma_om_l, sigma_k_r, sigma_om_r
    real::snorm, dervk_dervom, re_t_sst, rho_0, beta_i
    tolepsma = tolsmall
    if (turbulence .eq. 0) then
      viscl(4) = zero; viscl(3) = zero
      laml(4) = zero; laml(3) = zero
    else
      !modified on 19/6/2013
      select case (turbulencemodel)
        case (1)
          if (ispal .eq. 2) then
            turbmv(1) = eddyfl(2)
            turbmv(2) = eddyfr(2)
            chi = abs((turbmv(1))/(viscl(1)))
            chipow3 = chi*chi*chi
            fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
            viscl(3) = turbmv(1)*fv1
            chi = abs((turbmv(2))/(viscl(2)))
            chipow3 = chi*chi*chi
            fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
            viscl(4) = turbmv(2)*fv1
          else
            turbmv(1) = eddyfl(2)
            turbmv(2) = eddyfr(2)
            chi = abs(max(turbmv(1), tolepsma)/max(viscl(1), tolepsma))
            chipow3 = chi*chi*chi
            fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
            viscl(3) = turbmv(1)*fv1
            chi = abs(max(turbmv(2), tolepsma)/max(viscl(2), tolepsma))
            chipow3 = chi*chi*chi
            fv1 = chipow3/(chipow3 + (cv1*cv1*cv1))
            viscl(4) = turbmv(2)*fv1
          end if

        case (2)
          vortet(1, 1:2) = eddyfl(4:5)
          vortet(2, 1:2) = eddyfl(6:7)
          ux = vortet(1, 1); uy = vortet(1, 2)
          vx = vortet(2, 1); vy = vortet(2, 2)
          do ihgt = 1, 2
            do ihgj = 1, 2
              tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
            end do
          end do
          svort = 0.5*(vortet + tvort)
          snorm = sqrt(2.0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                            (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2))))
          wally = eddyfl(1)
          rho_0 = leftv(1)
          k_0 = max(tolepsma, eddyfl(2)/leftv(1))
          om_0 = max(eddyfl(3)/leftv(1), ufreestream/charlength/10.0)
          dervk_dervom = (eddyfl(8)*eddyfl(10)) + (eddyfl(9)*eddyfl(11))
          d_omplus = max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)
          phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(1)/(rho_0*wally*wally*om_0))
          phi_1 = min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally))
          f_1 = tanh(phi_1**4)
          f_2 = tanh(phi_2**2)
          re_t_sst = rho_0*k_0/(viscl(1)*om_0)
          beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
          alpha_star0 = beta_i/3.0
          alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
          alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
          viscl(3) = rho_0*k_0/om_0/max(1.0/alpha_star, snorm*f_2/(aa_1*om_0))
          sigma_k_l = sigma_k1/f_1 + sigma_k2/f_2
          sigma_om_l = sigma_om1/f_1 + sigma_om2/f_2
          if (eddyfr(1) .gt. 0.0) then
            vortet(1, 1:2) = eddyfl(4:5)
            vortet(2, 1:2) = eddyfl(6:7)
            ux = vortet(1, 1); uy = vortet(1, 2)
            vx = vortet(2, 1); vy = vortet(2, 2)
            do ihgt = 1, 2
              do ihgj = 1, 2
                tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
              end do
            end do
            svort = 0.5*(vortet + tvort)
            snorm = sqrt(2.0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                              (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2))))
            wally = eddyfr(1)
            rho_0 = rightv(1)
            k_0 = eddyfr(2)/rightv(1)
            om_0 = max(eddyfr(3)/rightv(1), 1.0e-6)
            dervk_dervom = (eddyfl(8)*eddyfl(10)) + (eddyfl(9)*eddyfl(11))
            d_omplus = max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !i need derivative of k
            phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(2)/(rho_0*wally*wally*om_0))
            phi_1 = min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally))
            f_1 = tanh(phi_1**4)
            f_2 = tanh(phi_2**2)
            re_t_sst = rho_0*k_0/(viscl(2)*om_0)
            beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
            alpha_star0 = beta_i/3.0
            alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
            alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
            viscl(4) = rho_0*k_0/om_0/max(1.0/alpha_star, snorm*f_2/(aa_1*om_0))
            sigma_k_r = sigma_k1/f_1 + sigma_k2/f_2
            sigma_om_r = sigma_om1/f_1 + sigma_om2/f_2
          else
            viscl(4) = -viscl(3)
            sigma_k_r = sigma_k_l
            sigma_om_r = sigma_om_l
          end if
      end select
      viscl(3) = min(10000000*visc, viscl(3))
      viscl(4) = min(10000000*visc, viscl(4))
      laml(3) = (viscl(3)*gamma/(prtu*(gamma - 1))) + (viscl(1)*gamma/(prandtl*(gamma - 1)))
      laml(4) = (viscl(4)*gamma/(prtu*(gamma - 1))) + (viscl(2)*gamma/(prandtl*(gamma - 1)))
      viscl(3) = max(0.0d0, viscl(3))
      viscl(4) = max(0.0d0, viscl(4))
      if ((turbmv(1) .lt. zero) .or. (turbmv(2) .lt. zero)) then
        viscl(3) = 0.0d0
        viscl(4) = 0.0d0
      end if
      etvm(1) = (0.5*(viscl(1) + viscl(2))) + (0.5*(viscl(3) + viscl(4)))
      if (turbulencemodel .eq. 2) then
        eddyfl(12) = viscl(1) + viscl(3)/sigma_k_l
        eddyfr(13) = viscl(2) + viscl(4)/sigma_k_r
        eddyfl(12) = viscl(1) + viscl(3)/sigma_om_l
        eddyfr(13) = viscl(2) + viscl(4)/sigma_om_r
      end if
    end if
  end subroutine eddyvisco2d

  subroutine trajectories
    implicit none
    integer::i, j, k, traj1, traj2, traj3, traj4, kmaxe, writeid, writeconf, num_vg
    real::win1, win2, win3, win4, post, post1, post2, post3, post4
    real, dimension(1:4)::pos_l, pos_g
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    kmaxe = xmpielrank(n)
    post1 = tolbig
    post2 = tolbig
    post3 = -tolbig
    traj1 = 0
    traj2 = 0
    traj3 = 0
    if (dimensiona .eq. 3) then
      num_vg = 5
    else
      num_vg = 4
    end if
    if (initcond .eq. 157) then
      pos_l(1:2) = zero
      pos_g(1:2) = zero
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if (dimensiona .eq. 3) then
          call cons2prim(n, leftv, mp_pinfl, gammal)
        else
          call cons2prim(n, leftv, mp_pinfl, gammal)
        end if
        if ((leftv(nof_variables)) .gt. 0.0d0) then  !total volume of gas evolution
          pos_l(2) = pos_l(2) + (leftv(nof_variables))*ielem(n, i)%totvolume
        end if
        if (leftv(num_vg) .gt. pos_l(1)) then  !total volume of gas evolution
          pos_l(1) = max(pos_l(1), leftv(num_vg))
        end if
      end do
      call mpi_allreduce(pos_l(1), pos_g(1), 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(2), pos_g(2), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      if (n .eq. 0) then
        open(70, file='volumex.dat', form='formatted', action='write', position='append')
        write(70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2)
        close(70)
      end if
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if (initcond .eq. 430) then
      pos_l(1:2) = zero
      pos_g(1:2) = zero
      do i = 1, kmaxe
        if (u_c(i)%val(1, 7) .gt. 0.1d0) then   !volume fraction of gas to be used for lowest location tracking
          if (ielem(n, i)%yyc .le. post3) then
            post3 = ielem(n, i)%yyc
            traj1 = i
          end if
        end if
        if (u_c(i)%val(1, 7) .gt. 0.0d0) then  !total volume of gas evolution
          pos_l(2) = pos_l(2) + u_c(i)%val(1, 7)*ielem(n, i)%totvolume
        end if
      end do
      pos_l(1) = post3  !lowest position in my local  cpu
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(1), pos_g(1), 1, mpi_double_precision, mpi_min, mpi_comm_world, ierror)
      writeid = 100000000
      if (abs(pos_g(1) - pos_l(1)) .le. tolsmall/1000) then
        if (traj1 .gt. 0) then
          writeid = n
        end if
      end if
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
      if (writeid .eq. writeconf) then
        leftv(1:nof_variables) = u_c(traj1)%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        open (70, file='position.dat', form='formatted', action='write', position='append')
        write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
        close (70)
      end if
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(2:2), pos_g(2:2), 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      if (n .eq. 0) then
        open (70, file='volume.dat', form='formatted', action='write', position='append')
        write (70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(2)
        close (70)
      end if
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if (initcond .eq. 405) then
      do i = 1, kmaxe
        if (dimensiona .eq. 3) then
          if (u_c(i)%val(1, 8) .gt. 0.05d0) then
            if (ielem(n, i)%xxc .le. post1) then
              post1 = ielem(n, i)%xxc
              traj1 = i
            end if
            if (((ielem(n,i)%yyc.le.0.052).and.(ielem(n,i)%yyc.ge.0.048)).and.((ielem(n,i)%zzc.le.0.052).and.(ielem(n,i)%zzc.ge.0.048)))then
              if (ielem(n, i)%xxc .le. post2) then
                post2 = ielem(n, i)%xxc
                traj2 = i
              end if
              if (ielem(n, i)%xxc .ge. post3) then
                post3 = ielem(n, i)%xxc
                traj3 = i
              end if
            end if
          end if
        else
          if (u_c(i)%val(1, 7) .gt. 0.4d0) then
            if (ielem(n, i)%xxc .le. post1) then
              post1 = ielem(n, i)%xxc
              traj1 = i
            end if
            if (((ielem(n, i)%yyc .le. 0.055) .and. (ielem(n, i)%yyc .ge. 0.045))) then
              if (ielem(n, i)%xxc .le. post2) then
                post2 = ielem(n, i)%xxc
                traj2 = i
              end if
              if (ielem(n, i)%xxc .ge. post3) then
                post3 = ielem(n, i)%xxc
                traj3 = i
              end if
            end if
          end if
        end if
      end do
      pos_l(1) = post1
      pos_l(2) = post2
      pos_l(3) = post3
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_min, mpi_comm_world, ierror)
      ! now select which cpu has the smallest
      writeid = 100000000
      writeconf = -1
      if (abs(pos_g(1) - pos_l(1)) .le. tolsmall/1000) then
        if (traj1 .gt. 0) then
          writeid = n
        end if
      end if
      ! and if more than one, select the smallest id one
      ! if (traj1.gt.0)then
      call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
      if (writeid .eq. writeconf) then
        if (traj1 .gt. 0) then
          leftv(1:nof_variables) = u_c(traj1)%val(1, 1:nof_variables)
          if (dimensiona .eq. 3) then
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(70, file='pos1.dat', form='formatted', action='write', position='append')
            write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
            close(70)
          else
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(70, file='pos1.dat', form='formatted', action='write', position='append')
            write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
            close(70)
          end if
        end if
      end if
      writeid = 100000000
      writeconf = -1
      if (abs(pos_g(2) - pos_l(2)) .le. tolsmall/1000) then
        if (traj2 .gt. 0) then
          writeid = n
        end if
      end if
      call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
      if (writeid .eq. writeconf) then
        if (traj2 .gt. 0) then
          leftv(1:nof_variables) = u_c(traj2)%val(1, 1:nof_variables)
          if (dimensiona .eq. 3) then
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(71, file='pos2.dat', form='formatted', action='write', position='append')
            write(71,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(2),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
            close(71)
          else
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(70, file='pos2.dat', form='formatted', action='write', position='append')
            write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(2),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
            close(70)
          end if
        end if
      end if
      call mpi_allreduce(pos_l(3), pos_g(3), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      writeid = 100000000
      writeconf = -1
      if (abs(pos_g(3) - pos_l(3)) .le. tolsmall/1000) then
        if (traj3 .gt. 0) then
          writeid = n
        end if
      end if
      call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
      if (writeid .eq. writeconf) then
        if (traj3 .gt. 0) then
          leftv(1:nof_variables) = u_c(traj3)%val(1, 1:nof_variables)
          if (dimensiona .eq. 3) then
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(71, file='pos3.dat', form='formatted', action='write', position='append')
            write(71,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(3),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
            close(71)
          else
            call cons2prim(n, leftv, mp_pinfl, gammal)
            open(70, file='pos3.dat', form='formatted', action='write', position='append')
            write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(3),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
            close(70)
          end if
        end if
      end if
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if ((initcond .eq. 408) .or. (initcond .eq. 422)) then
      pos_l(1:2) = zero
      pos_g(1:2) = zero
      do i = 1, kmaxe
        if (u_c(i)%val(1, 8) .gt. 0.0d0) then
          pos_l(1) = pos_l(1) + ielem(n, i)%totvolume
          pos_l(2) = pos_l(2) + u_c(i)%val(1, 8)*ielem(n, i)%totvolume
        end if
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      if (n .eq. 0) then
        open(70, file='volume.dat', form='formatted', action='write', position='append')
        write(70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2)
        close(70)
      end if
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if ((initcond .eq. 411) .or. (initcond .eq. 444)) then
      pos_l(1:3) = zero
      pos_g(1:3) = zero
      do i = 1, kmaxe
        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if (dimensiona .eq. 2) then
          call cons2prim(n, leftv, mp_pinfl, gammal)
          pos_l(3) = max(abs(leftv(4)), pos_l(3))
          if (u_c(i)%val(1, 7) .gt. 0.0d0) then
            pos_l(1) = pos_l(1) + ielem(n, i)%totvolume
            pos_l(2) = pos_l(2) + u_c(i)%val(1, 7)*ielem(n, i)%totvolume
          end if
        else
          call cons2prim(n, leftv, mp_pinfl, gammal)
          pos_l(3) = max(abs(leftv(5)), pos_l(3))
          if (u_c(i)%val(1, 8) .gt. 0.0d0) then
            pos_l(1) = pos_l(1) + ielem(n, i)%totvolume
            pos_l(2) = pos_l(2) + u_c(i)%val(1, 7)*ielem(n, i)%totvolume
          end if
        end if
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(pos_l(3), pos_g(3), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
      if (n .eq. 0) then
        open(70, file='volume.dat', form='formatted', action='write', position='append')
        write(70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2), pos_g(3)
        close(70)
      end if
      call mpi_barrier(mpi_comm_world, ierror)
    end if
  end subroutine

  subroutine flux2dx(flux_term_x, leftv)
    implicit none
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ie1, mp_stiff, mp_density, gammal, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_x
    if (multispecies .eq. 1) then
      mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
      gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
      mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
      mp_density = leftv(5) + leftv(6)
      r = mp_density
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2))
      !internal energy
      ie1 = ((p + mp_stiff)/((gammal - 1.0d0)*r))
      !total energy
      e = r*(skin + ie1)
      flux_term_x(1) = r*u
      flux_term_x(2) = (r*(u**2)) + p
      flux_term_x(3) = r*u*v
      flux_term_x(4) = u*(e + p)
      flux_term_x(5:7) = leftv(5:7)*u
    else
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)
      gm = gamma
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2))
      !internal energy
      ien = ((p)/((gm - 1.0d0)*r))
      !total energy
      e = r*(skin + ien)
      flux_term_x(1) = r*u
      flux_term_x(2) = (r*(u**2)) + p
      flux_term_x(3) = r*u*v
      flux_term_x(4) = u*(e + p)
    end if
  end subroutine

  subroutine flux2dy(flux_term_y, leftv)
    implicit none
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ie1, mp_stiff, mp_density, gammal, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_y
    if (multispecies .eq. 1) then
      mp_ar(1) = leftv(7)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = (1.0d0 - leftv(7))/(gamma_in(2) - 1.0d0)
      gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
      mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
      mp_density = leftv(5) + leftv(6)
      r = mp_density
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2))
      !internal energy
      ie1 = ((p + mp_stiff)/((gammal - 1.0d0)*r))
      !total energy
      e = r*(skin + ie1)
      flux_term_y(1) = r*v
      flux_term_y(2) = r*u*v
      flux_term_y(3) = (r*(v**2)) + p
      flux_term_y(4) = v*(e + p)
      flux_term_y(5:7) = leftv(5:7)*v
    else
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)
      gm = gamma
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2))
      !internal energy
      ien = ((p)/((gm - 1.0d0)*r))
      !total energy
      e = r*(skin + ien)
      flux_term_y(1) = r*v
      flux_term_y(2) = r*u*v
      flux_term_y(3) = (r*v*v) + p
      flux_term_y(4) = v*(e + p)
    end if
  end subroutine

  subroutine flux3dx(flux_term_x, leftv)
    implicit none
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ie1, mp_stiff, mp_density, gammal, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_x
    if (multispecies .eq. 1) then
      mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
      gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
      mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
      mp_density = leftv(6) + leftv(7)
      r = mp_density
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ie1 = ((p + mp_stiff)/((gammal - 1.0d0)*r))
      !total energy
      e = r*(skin + ie1)
      flux_term_x(1) = r*u
      flux_term_x(2) = (r*(u**2)) + p
      flux_term_x(3) = r*u*v
      flux_term_x(4) = r*u*w
      flux_term_x(5) = u*(e + p)
      flux_term_x(6:8) = leftv(6:8)*u
    else
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      gm = gamma
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ien = ((p)/((gm - 1.0d0)*r))
      !total energy
      e = r*(skin + ien)
      flux_term_x(1) = r*u
      flux_term_x(2) = (r*(u**2)) + p
      flux_term_x(3) = r*u*v
      flux_term_x(4) = r*u*w
      flux_term_x(5) = u*(e + p)
    end if
  end subroutine

  subroutine flux3dy(flux_term_y, leftv)
    implicit none
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ie1, mp_stiff, mp_density, gammal, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_y
    if (multispecies .eq. 1) then
      mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
      gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
      mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
      mp_density = leftv(6) + leftv(7)
      r = mp_density
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ie1 = ((p + mp_stiff)/((gammal - 1.0d0)*r))
      !total energy
      e = r*(skin + ie1)
      flux_term_y(1) = r*v
      flux_term_y(2) = r*u*v
      flux_term_y(3) = (r*(v**2)) + p
      flux_term_y(4) = r*v*w
      flux_term_y(5) = v*(e + p)
      flux_term_y(6:8) = leftv(6:8)*v
    else
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      gm = gamma
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ien = ((p)/((gm - 1.0d0)*r))
      !total energy
      e = r*(skin + ien)
      flux_term_y(1) = r*v
      flux_term_y(2) = r*u*v
      flux_term_y(3) = (r*(v*v)) + p
      flux_term_y(4) = r*v*w
      flux_term_y(5) = v*(e + p)
    end if
  end subroutine

  subroutine flux3dz(flux_term_z, leftv)
    implicit none
    real::p, u, v, w, e, r, s, gm, skin, ien, pi, ie1, mp_stiff, mp_density, gammal, gammar
    real, dimension(nof_species)::mp_ar, mp_ie
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_z
    if (multispecies .eq. 1) then
      mp_ar(1) = leftv(8)/(gamma_in(1) - 1.0d0)
      mp_ar(2) = (1.0d0 - leftv(8))/(gamma_in(2) - 1.0d0)
      gammal = (1.0d0/(mp_ar(1) + mp_ar(2))) + 1.0d0    !mixture gamma isobaric assumptio
      mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
      mp_density = leftv(6) + leftv(7)
      r = mp_density
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ie1 = ((p + mp_stiff)/((gammal - 1.0d0)*r))
      !total energy
      e = r*(skin + ie1)
      flux_term_z(1) = r*w
      flux_term_z(2) = r*u*w
      flux_term_z(3) = r*v*w
      flux_term_z(4) = (r*(w*w)) + p
      flux_term_z(5) = w*(e + p)
      flux_term_z(6:8) = leftv(6:8)*w
    else
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      gm = gamma
      !kinetic energy first!
      skin = (oo2)*((u**2) + (v**2) + (w**2))
      !internal energy
      ien = ((p)/((gm - 1.0d0)*r))
      !total energy
      e = r*(skin + ien)
      flux_term_z(1) = r*w
      flux_term_z(2) = r*u*w
      flux_term_z(3) = r*v*w
      flux_term_z(4) = (r*(w**2)) + p
      flux_term_z(5) = w*(e + p)
    end if
  end subroutine

  subroutine flux_visc2d(flux_term_x, flux_term_y, leftv, leftv_der)
    implicit none
    real:: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_x, flux_term_y
    real, dimension(1:4)::viscl, laml
    real, dimension(1:nof_variables, 1:dimensiona), intent(in)::leftv_der
    call sutherland2d(n, leftv, leftv, viscl, laml)
    u = leftv(2)
    v = leftv(3)
    ux = leftv_der(2, 1)
    uy = leftv_der(2, 2)
    vx = leftv_der(3, 1)
    vy = leftv_der(3, 2)
    tx = leftv_der(4, 1)
    ty = leftv_der(4, 2)
    tauxx = 2.0d0/3.0d0*viscl(1)*(2*ux - vy)
    tauyy = 2.0d0/3.0d0*viscl(1)*(2*vy - ux)
    tauxy = viscl(1)*(uy + vx)
    flux_term_x(2) = flux_term_x(2) - tauxx
    flux_term_x(3) = flux_term_x(3) - tauxy
    flux_term_x(4) = flux_term_x(4) - (u*tauxx + v*tauxy + laml(1)*tx)
    flux_term_y(2) = flux_term_y(2) - tauxy
    flux_term_y(3) = flux_term_y(3) - tauyy
    flux_term_y(4) = flux_term_y(4) - (u*tauxy + v*tauyy + laml(1)*ty)
  end subroutine

  subroutine flux_visc3d(flux_term_x, flux_term_y, flux_term_z, leftv, leftv_der)
    implicit none
    real:: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy, w, uz, vz, wx, wy, wz, tz, tauxz, tauyz, tauzz
    real, dimension(1:nof_variables), intent(in)::leftv
    real, dimension(1:nof_variables), intent(inout)::flux_term_x, flux_term_y, flux_term_z
    real, dimension(1:nof_variables, 1:dimensiona), intent(in)::leftv_der
    real, dimension(1:4)::viscl, laml
    call sutherland(n, leftv, leftv, viscl, laml)
    ! variables extrapolated at boundary
    u = leftv(2)
    v = leftv(3)
    w = leftv(4)
    ux = leftv_der(2, 1)
    uy = leftv_der(2, 2)
    uz = leftv_der(2, 3)
    vx = leftv_der(3, 1)
    vy = leftv_der(3, 2)
    vz = leftv_der(3, 3)
    wx = leftv_der(4, 1)
    wy = leftv_der(4, 2)
    wz = leftv_der(4, 3)
    tx = leftv_der(5, 1)
    ty = leftv_der(5, 2)
    tz = leftv_der(5, 3)
    tauxx = 2.0d0/3.0d0*viscl(1)*(2*ux - vy - wz)
    tauyy = 2.0d0/3.0d0*viscl(1)*(2*vy - ux - wz)
    tauxy = viscl(1)*(uy + vx)
    tauxz = viscl(1)*(uz + wx)
    tauyz = viscl(1)*(wy + vz)
    tauzz = 2.0d0/3.0d0*viscl(1)*(2*wz - ux - vy)
    flux_term_x(2) = flux_term_x(2) - tauxx
    flux_term_x(3) = flux_term_x(3) - tauxy
    flux_term_x(4) = flux_term_x(4) - tauxz
    flux_term_x(5) = flux_term_x(5) - (u*tauxx + v*tauxy + w*tauxz + laml(1)*tx)
    flux_term_y(2) = flux_term_y(2) - tauxy
    flux_term_y(3) = flux_term_y(3) - tauyy
    flux_term_y(4) = flux_term_y(4) - tauyz
    flux_term_y(5) = flux_term_y(5) - (u*tauxy + v*tauyy + w*tauyz + laml(1)*ty)
    flux_term_z(2) = flux_term_z(2) - tauxz
    flux_term_z(3) = flux_term_z(3) - tauyz
    flux_term_z(4) = flux_term_z(4) - tauzz
    flux_term_z(5) = flux_term_z(5) - (u*tauxz + v*tauyz + w*tauzz + laml(1)*tz)
  end subroutine

  subroutine dcons2dprim(leftv_der, leftv)
    implicit none
    real, dimension(1:nof_variables, 1:dimensiona), intent(inout)::leftv_der
    real, dimension(1:nof_variables), intent(in)::leftv
    integer:: i_dim, i_var
    do i_dim = 1, dimensiona
      do i_var = 2, nof_variables - 1
        leftv_der(i_var, i_dim) = (leftv_der(i_var, i_dim) - leftv_der(1, i_dim)*leftv(i_var))/leftv(1) ! ux = (rhoux - rhox * u) / rho
      end do
        leftv_der(nof_variables,i_dim) = (gamma - 1.0d0) * ((leftv(1) * leftv_der(nof_variables,i_dim) - leftv(nof_variables) * leftv_der(1,i_dim)) / leftv(1) ** 2 - dot_product(leftv(2:nof_variables-1), leftv_der(2:nof_variables-1,i_dim)))
    end do
  end subroutine dcons2dprim
end module flow_operations
