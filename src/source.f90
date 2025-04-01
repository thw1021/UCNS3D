module source
  use library
  use transform
  use local
  use riemann
  use flow_operations
  use declaration
  implicit none

contains
  subroutine sources_computation(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, kmaxe, iconsidered
    real, dimension(turbulenceequations) :: source_t
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call sources(n, iconsidered, source_t)
      rhst(i)%val(1:turbulenceequations) = rhst(i)%val(1:turbulenceequations) - (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
    end do
  end subroutine sources_computation

  subroutine sources_computation_rot(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, srf
    real::oodensity
    real, dimension(5)::source_t2
    real, dimension(1:dimensiona)::pox, poy, poz
    kmaxe = xmpielrank(n)
    if (srfg .eq. 1) then
      do i = 1, kmaxe
        oodensity = 1.0d0/u_c(i)%val(1, 1)
        source_t2(1) = 0.0
        source_t2(2) = u_c(i)%val(1, 2)*oodensity - uvel
        source_t2(3) = u_c(i)%val(1, 3)*oodensity - vvel
        source_t2(4) = u_c(i)%val(1, 4)*oodensity - wvel
        source_t2(5) = 0.0
        pox(1:3) = source_t2(2:4)
        poy(1:3) = srf_velocity(1:3)
        source_t2(2:4) = u_c(i)%val(1, 1)*vect_function(pox, poy)
        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + (source_t2(1:nof_variables)*ielem(n, i)%totvolume)
      end do
    end if
    if (mrf .eq. 1) then
      do i = 1, kmaxe
        srf = ilocal_recon3(i)%mrf
        if (ilocal_recon3(i)%mrf .eq. 1) then
          oodensity = 1.0d0/u_c(i)%val(1, 1)
          source_t2(1) = 0.0
          source_t2(2) = u_c(i)%val(1, 2)*oodensity
          source_t2(3) = u_c(i)%val(1, 3)*oodensity
          source_t2(4) = u_c(i)%val(1, 4)*oodensity
          source_t2(5) = 0.0
          pox(1:3) = source_t2(2:4)
          poy(1:3) = ilocal_recon3(i)%mrf_velocity(1:3)
          source_t2(2:4) = u_c(i)%val(1, 1)*vect_function(pox, poy)
          rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + (source_t2(1:nof_variables)*ielem(n, i)%totvolume)
        end if
      end do
    end if
  end subroutine sources_computation_rot

  subroutine sources_derivatives_computation(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, iconsidered
    real, dimension(turbulenceequations)::source_t
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call sources_derivatives(n, iconsidered, source_t)
      sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
    end do
  end subroutine sources_derivatives_computation

  subroutine sources(n, iconsidered, source_t)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(turbulenceequations), intent(inout)::source_t
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
    real::vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
    real::tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
    integer::i, k, j, l, ihgt, ihgj, iex, lowre
    real::snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
    real::gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
    real:: r_des, f_des, ddw, f_des_sst, l_t_des
    real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
    real:: uz, vz, wx, wy, wz
    real:: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz  !for sas only
    real, dimension(3, 3)::vortet, tvort, svort, ovort
    real, dimension(3)::vortem, vertf, vertex
    real, dimension(turbulenceequations, 1:3)::derivturb
    real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
    real:: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
    real:: alpha_raw, alpha_star, re_t_sst, alpha_inf
    real:: beta_stari, beta_i, beta_raw, beta_star
    real:: k_0, om_0, wally
    real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
    real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
    real:: kx, ky, kz, omx, omy, omz
    real, dimension(1:nof_variables)::leftv, rightv
    real::mp_pinfl, gammal
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    i = iconsidered
    verysmall = 10e-16
    vortet(1:3, 1:3) = ilocal_recon3(i)%grads(1:3, 1:3)
    ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
    vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
    wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)
    do ihgt = 1, 3
      do ihgj = 1, 3
        tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
      end do
    end do
    svort = 0.5*(vortet + tvort)
    ovort = 0.5*(vortet - tvort)
    snorm = sqrt(2.0d0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + (svort(1, 3)*svort(1, 3)) + &
                        (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)) + (svort(2, 3)*svort(2, 3)) + &
                        (svort(3, 1)*svort(3, 1)) + (svort(3, 2)*svort(3, 2)) + (svort(3, 3)*svort(3, 3))))
    onorm = sqrt(2.0d0*((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + (ovort(1, 3)*ovort(1, 3)) + &
                        (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2)) + (ovort(2, 3)*ovort(2, 3)) + &
                        (ovort(3, 1)*ovort(3, 1)) + (ovort(3, 2)*ovort(3, 2)) + (ovort(3, 3)*ovort(3, 3))))
    omega = onorm
    divnorm = ux + uy + uz
    usss = sqrt((2.0*((ux*ux) + (vy*vy) + (wz*wz))) &
                + ((uy + vx)*(uy + vx) + (uz + wx)*(uz + wx) + (wy + vz)*(wy + vz)) &
                - (2.0/3.0*(ux + vy + wz)*(ux + vy + wz)))
    derivturb(1, 1:3) = ilocal_recon3(i)%grads(5, 1:3)
    squaret = (sqrt((derivturb(1, 1)**2) + (derivturb(1, 2)**2) + (derivturb(1, 3)**2)))**2
    leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    rightv(1:nof_variables) = leftv(1:nof_variables)
    call sutherland(n, leftv, rightv, viscl, laml)
    turbmv(1) = u_ct(i)%val(1, 1)
    turbmv(2) = turbmv(1)
    select case (turbulencemodel)

    case (1)
      if (rot_corr .eq. 1) then
        omega = omega + 2.0d0*min(0.0d0, snorm - onorm)
      end if
      if (ispal .eq. 1) then
        eddyfl(2) = turbmv(1)
        eddyfr(2) = turbmv(2)
        onesix = 1.0d0/6.0d0
        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        cw1 = cb1/(kappa*kappa)
        cw1 = cw1 + (1.0 + cb2)/sigma
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .lt. verysmall) then
          source_t(1) = 0.0
        else
          tch_x3 = (tch_x)*(tch_x)*(tch_x)
          tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
          tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
          ddw = ielem(n, i)%walldist
          if (des_model .eq. 1) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            ddw = min(ddw, c_des_sa*delta_cell)
          end if
          if (des_model .eq. 2) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            r_des = min(10.0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
            f_des = 1 - tanh((8.0*r_des)**3)
            ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
          end if
          prodterm1 = (turbmv(1))/(leftv(1)*kappa*kappa*ddw*ddw)
          stild = max(omega + (tch_fv2*prodterm1), 0.3*omega)
          prodtermfinal = stild*turbmv(1)*cb1
          rr = min((turbmv(1)/(((leftv(1)*stild*kappa*kappa*(ddw)*(ddw))) + 0.000000001)), 10.0)
          gg = rr + (cw2*(rr**6 - rr))
          fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
          destterm = cw1*fw*(((turbmv(1))/(leftv(1)*(ddw)))**2)
          fodt = cb2*squaret/sigma
          if (d_corr .eq. 0) then
            source_t(1) = prodtermfinal + fodt - destterm
          else
            source_t(1) = prodtermfinal + (fodt - destterm)*leftv(1)
          end if
        end if
      else
        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .gt. 10.0d0) then
          tch_x = tch_x
        else
          tch_x = 0.05d0*log(1.0d0 + exp(20.0d0*(tch_x)))
        end if
        ddw = ielem(n, i)%walldist
        if (des_model .eq. 1) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          ddw = min(ddw, c_des_sa*delta_cell)
        end if
        if (des_model .eq. 2) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          r_des = min(10.0d0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
          f_des = 1 - tanh((8.0d0*r_des)**3)
          ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
        end if
        tch_x3 = (tch_x)*(tch_x)*(tch_x)
        tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
        tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
        prodterm1 = (tch_fv2*tch_x*(turbmv(1)))/(leftv(1)*kappa*kappa*(ddw)*(ddw))
        if (prodterm1 .ge. (-0.7d0*omega)) then
          stild = omega + prodterm1
        end if
        if (prodterm1 .lt. (-0.7d0*omega)) then
          stild = omega + ((omega*(((0.7d0*0.7d0)*(omega)) + (0.9*prodterm1)))/(((0.9 - 1.4)*omega) - prodterm1))
        end if
        prodtermfinal = stild*turbmv(1)*cb1*tch_x
        rr = ((turbmv(1)*tch_x/(leftv(1)*stild*kappa*kappa*(ddw)*(ddw))))
        gg = rr + (cw2*(rr**6 - rr))
        fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
        destterm = cw1*fw*leftv(1)*(((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
        fodt = leftv(1)*cb2*squaret/sigma
        source_t(1) = prodtermfinal + fodt - destterm
      end if

    case (2)                !k omega sst
      eddyfl(1) = ielem(n, i)%walldist
      eddyfl(2) = u_ct(i)%val(1, 1)
      eddyfl(3) = u_ct(i)%val(1, 2)
      eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3)
      eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
      eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3)
      eddyfl(13:15) = ilocal_recon3(i)%grads(4, 1:3)
      eddyfl(16:18) = ilocal_recon3(i)%grads(5, 1:3)
      eddyfr = eddyfl
      call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      k_0 = (max(verysmall, u_ct(i)%val(1, 1)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0 = max(1.0e-1*ufreestream/charlength, u_ct(i)%val(1, 2)/leftv(1))
      wally = ielem(n, i)%walldist
      omx = ilocal_recon3(i)%grads(6, 1)
      omy = ilocal_recon3(i)%grads(6, 2)
      omz = ilocal_recon3(i)%grads(6, 3)
      kx = ilocal_recon3(i)%grads(5, 1)
      ky = ilocal_recon3(i)%grads(5, 2)
      kz = ilocal_recon3(i)%grads(5, 3)
      dervk_dervom = kx*omx + ky*omy + kz*omz
      d_omplus = max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
      phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
      phi_1 = min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally)
      f_1 = tanh(phi_1**4)
      f_2 = tanh(phi_2**2)
      alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
      alpha_star = alpha_starinf
      alpha_raw = alpha_inf
      beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
      alpha_star0 = beta_i/3.0
      beta_stari = beta_starinf
      lowre = 0
      if (lowre .eq. 1) then
        re_t_sst = leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???
        alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
        alpha_raw = alpha_inf/alpha_star*(alpha_0 + re_t_sst/r_om_sst)/(1.0 + re_t_sst/r_om_sst)
        beta_stari = beta_starinf*(4.0/15.0 + (re_t_sst/r_beta)**4)/(1.0 + (re_t_sst/r_beta)**4)
      end if
      beta_star = beta_stari
      beta_raw = beta_i
      if (vort_model .eq. 1) then
        prod_k = viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
        prod_om = alpha_raw*leftv(1)*snorm*onorm
      else
        prod_k = viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
        prod_k = min(prod_k, 10.0*leftv(1)*beta_star*k_0*om_0)
      end if
      ydest_k = leftv(1)*beta_star*k_0*om_0
      ydest_om = leftv(1)*beta_raw*om_0*om_0
      diff_om = 2.0*(1 - f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
      if (qsas_model .eq. 1) then 
        do iex = 1, 3
          vortet(iex, 1:3) = ilocal_recon3(i)%grads(3 + turbulenceequations + iex, 1:3)
        end do
        uxx = vortet(1, 1); uyy = vortet(1, 2); uzz = vortet(1, 3)
        vxx = vortet(2, 1); vyy = vortet(2, 2); vzz = vortet(2, 3)
        wxx = vortet(3, 1); wyy = vortet(3, 2); wzz = vortet(3, 3)
        cell_volume = ielem(n, i)%totvolume
        u_lapl = sqrt((uxx + uyy + uzz)**2 + (vxx + vyy + vzz)**2 + (wxx + wyy + wzz)**2)
        dervk2 = kx*kx + ky*ky + kz*kz
        dervom2 = omx*omx + omy*omy + omz*omz
        delta_cell = cell_volume**0.333333333333333
        l_sas = sqrt(k_0)/(beta_star**0.25*om_0)
        l_vk = max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
                   c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star - alpha_raw)))
        q_sas1 = leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
        q_sas2 = -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2, dervom2/om_0**2)
        q_sas = max(q_sas1 + q_sas2, 0.0)
      else
        q_sas = zero
      end if
      !final source terms
      if (qsas_model .eq. 2) then
        cell_volume = ielem(n, i)%totvolume
        delta_cell = cell_volume**0.333333333333333
        l_t_des = sqrt(k_0)/(beta_star*om_0)
        f_des_sst = max(1.0, l_t_des/(c_des_sst*delta_cell)*(1 - f_2))
        ydest_k = f_des_sst*ydest_k
      end if
      srcfull_k = prod_k - ydest_k
      srcfull_om = prod_om - ydest_om + diff_om + q_sas
      source_t(1) = srcfull_k
      source_t(2) = srcfull_om
    end select
  end subroutine sources

  subroutine sources_derivatives(n, iconsidered, source_t)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(turbulenceequations), intent(inout)::source_t
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
    real::vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
    real::tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
    integer::i, k, j, l, ihgt, ihgj, iex, lowre
    real::snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
    real::gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
    real:: r_des, f_des, ddw, f_des_sst, l_t_des
    real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
    real:: uz, vz, wx, wy, wz
    real:: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz  !for sas only
    real, dimension(3, 3)::vortet, tvort, svort, ovort
    real, dimension(3)::vortem, vertf, vertex
    real, dimension(turbulenceequations, 1:3)::derivturb
    real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
    real:: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
    real:: alpha_raw, alpha_star, re_t_sst, alpha_inf
    real:: beta_stari, beta_i, beta_raw, beta_star
    real:: k_0, om_0, wally
    real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
    real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
    real:: kx, ky, kz, omx, omy, omz
    real, dimension(1:nof_variables)::leftv, rightv
    real::mp_pinfl, gammal
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    i = iconsidered
    verysmall = 10e-16
    vortet(1:3, 1:3) = ilocal_recon3(i)%grads(1:3, 1:3)
    ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
    vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
    wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)
    do ihgt = 1, 3
      do ihgj = 1, 3
        tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
      end do
    end do
    svort = 0.5*(vortet + tvort)
    ovort = 0.5*(vortet - tvort)
    snorm = sqrt(2.0d0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + (svort(1, 3)*svort(1, 3)) + &
                        (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2)) + (svort(2, 3)*svort(2, 3)) + &
                        (svort(3, 1)*svort(3, 1)) + (svort(3, 2)*svort(3, 2)) + (svort(3, 3)*svort(3, 3))))
    onorm = sqrt(2.0d0*((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + (ovort(1, 3)*ovort(1, 3)) + &
                        (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2)) + (ovort(2, 3)*ovort(2, 3)) + &
                        (ovort(3, 1)*ovort(3, 1)) + (ovort(3, 2)*ovort(3, 2)) + (ovort(3, 3)*ovort(3, 3))))
    omega = onorm
    divnorm = ux + uy + uz  !careful with the sign. if it becomes very big, it can produce negative production
    usss = sqrt((2.0*((ux*ux) + (vy*vy) + (wz*wz))) &
                + ((uy + vx)*(uy + vx) + (uz + wx)*(uz + wx) + (wy + vz)*(wy + vz)) &
                - (2.0/3.0*(ux + vy + wz)*(ux + vy + wz)))
    derivturb(1, 1:3) = ilocal_recon3(i)%grads(5, 1:3)
    squaret = (sqrt((derivturb(1, 1)**2) + (derivturb(1, 2)**2) + (derivturb(1, 3)**2)))**2
    leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    rightv(1:nof_variables) = leftv(1:nof_variables)
    call sutherland(n, leftv, rightv, viscl, laml)
    turbmv(1) = u_ct(i)%val(1, 1)
    turbmv(2) = turbmv(1)
    select case (turbulencemodel)

    case (1) !!spalart almaras model
      if (ispal .eq. 1) then
        eddyfl(2) = turbmv(1)
        eddyfr(2) = turbmv(2)
        onesix = 1.0d0/6.0d0
        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        cw1 = cb1/(kappa*kappa)
        cw1 = cw1 + (1.0 + cb2)/sigma
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .lt. verysmall) then
          source_t(1) = zero
        else
          tch_x3 = (tch_x)*(tch_x)*(tch_x)
          tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
          tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
          ddw = ielem(n, i)%walldist
          if (des_model .eq. 1) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            ddw = min(ddw, c_des_sa*delta_cell)
          end if
          if (des_model .eq. 2) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            r_des = min(10.0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
            !previous limiter is just for numerical reasons regarding tanh
            f_des = 1 - tanh((8.0*r_des)**3)
            ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
          end if
          prodterm1 = (turbmv(1))/(leftv(1)*kappa*kappa*ddw*ddw)
          stild = max(omega + (tch_fv2*prodterm1), 0.3*omega)
          prodtermfinal = stild*turbmv(1)*cb1/leftv(1)
          rr = min((turbmv(1)/(((leftv(1)*stild*kappa*kappa*(ddw)*(ddw))) + 0.000000001)), 10.0)
          gg = rr + (cw2*(rr**6 - rr))
          fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
          destterm = cw1*fw*(((turbmv(1))/(leftv(1)*(ddw)))**2)
          fodt = cb2*squaret/sigma
          prodtermfinal = stild*cb1
          destterm = 2.0*cw1*fw*(turbmv(1)/(leftv(1)*ddw**2))
          fodt = 2.0*cb2*(sqrt(squaret))/sigma
          source_t(1) = min(prodtermfinal + fodt - destterm, zero)
        end if
      else
        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .gt. 10.0d0) then
          tch_x = tch_x
        else
          tch_x = 0.05d0*log(1.0d0 + exp(20.0d0*(tch_x)))
        end if
        ddw = ielem(n, i)%walldist
        if (des_model .eq. 1) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          ddw = min(ddw, c_des_sa*delta_cell)
        end if
        if (des_model .eq. 2) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          r_des = min(10.0d0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
          !previous limiter is just for numerical reasons regarding tanh
          f_des = 1 - tanh((8.0d0*r_des)**3)
          ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
        end if
        tch_x3 = (tch_x)*(tch_x)*(tch_x)
        tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
        tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
        prodterm1 = (tch_fv2*tch_x*(turbmv(1)))/(leftv(1)*kappa*kappa*(ddw)*(ddw))
        if (prodterm1 .ge. (-0.7d0*omega)) then
          stild = omega + prodterm1
        end if
        if (prodterm1 .lt. (-0.7d0*omega)) then
          stild = omega + ((omega*(((0.7d0*0.7d0)*(omega)) + (0.9*prodterm1)))/(((0.9 - 1.4)*omega) - prodterm1))
        end if
        prodtermfinal = stild*turbmv(1)*cb1*tch_x
        rr = ((turbmv(1)*tch_x/(leftv(1)*stild*kappa*kappa*(ddw)*(ddw))))
        gg = rr + (cw2*(rr**6 - rr))
        fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
        destterm = cw1*fw*leftv(1)*(((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
        fodt = leftv(1)*cb2*squaret/sigma
        prodtermfinal = stild*cb1*tch_x
        destterm = 2.0*cw1*fw*(turbmv(1)/(leftv(1)*ddw**2))
        fodt = 2.0*cb2*(sqrt(squaret))/sigma
        source_t(1) = min(prodtermfinal + fodt - destterm, zero)
      end if

    case (2)                !k omega sst
      eddyfl(1) = ielem(n, i)%walldist
      eddyfl(2) = u_ct(i)%val(1, 1)
      eddyfl(3) = u_ct(i)%val(1, 2)
      eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3)
      eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
      eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3)
      eddyfl(13:15) = ilocal_recon3(i)%grads(4, 1:3)
      eddyfl(16:18) = ilocal_recon3(i)%grads(5, 1:3)
      eddyfr = eddyfl
      call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      k_0 = (max(verysmall, u_ct(i)%val(1, 1)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0 = max(1.0e-1*ufreestream/charlength, u_ct(i)%val(1, 2)/leftv(1))
      source_t(1) = beta_starinf*om_0
      source_t(2) = beta_starinf*k_0
    end select
  end subroutine sources_derivatives

  subroutine sources_computation2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, iconsidered
    real, dimension(turbulenceequations)::source_t
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call sources2d(n, iconsidered, source_t)
   rhst(i)%val(1:turbulenceequations) = rhst(i)%val(1:turbulenceequations) - (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
    end do
  end subroutine sources_computation2d

  subroutine sources_derivatives_computation2d(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, iconsidered
    real, dimension(turbulenceequations)::source_t
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call sources_derivatives2d(n, iconsidered, source_t)
      sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
    end do
  end subroutine sources_derivatives_computation2d

  subroutine sources2d(n, iconsidered, source_t)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(turbulenceequations), intent(inout)::source_t
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
    real::vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
    real::tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
    integer::i, k, j, l, ihgt, ihgj, iex, lowre
    real::snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
    real::gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
    real:: r_des, f_des, ddw, f_des_sst, l_t_des
    real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
    real:: uz, vz, wx, wy, wz
    real:: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz  !for sas only
    real, dimension(2, 2)::vortet, tvort, svort, ovort
    real, dimension(2)::vortem, vertf, vertex
    real, dimension(turbulenceequations, 1:2)::derivturb
    real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
    real:: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
    real:: alpha_raw, alpha_star, re_t_sst, alpha_inf
    real:: beta_stari, beta_i, beta_raw, beta_star
    real:: k_0, om_0, wally
    real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
    real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
    real:: kx, ky, kz, omx, omy, omz
    real, dimension(1:nof_variables)::leftv, rightv
    real::mp_pinfl, gammal, mp_pinfr, gammar
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    i = iconsidered
    verysmall = 10e-16
    vortet(1:2, 1:2) = ilocal_recon3(i)%grads(1:2, 1:2)
    ux = vortet(1, 1); uy = vortet(1, 2)
    vx = vortet(2, 1); vy = vortet(2, 2)
    do ihgt = 1, 2
      do ihgj = 1, 2
        tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
      end do
    end do
    svort = 0.5*(vortet + tvort)
    ovort = 0.5*(vortet - tvort)
    snorm = sqrt(2.0d0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                        (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2))))
    onorm = sqrt(2.0d0*((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + &
                        (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2))))
    omega = onorm
    divnorm = ux + uy + uz  !careful with the sign. if it becomes very big, it can produce negative production
    usss = sqrt((2.0*((ux*ux) + (vy*vy))) &
                + ((uy + vx)*(uy + vx)) &
                - (2.0/3.0*(ux + vy)*(ux + vy)))
    derivturb(1, 1:2) = ilocal_recon3(i)%grads(4, 1:2)
    squaret = (sqrt((derivturb(1, 1)**2) + (derivturb(1, 2)**2)))**2
    leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    rightv(1:nof_variables) = leftv(1:nof_variables)
    call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
    call sutherland2d(n, leftv, rightv, viscl, laml)
    turbmv(1) = u_ct(i)%val(1, 1)
    turbmv(2) = turbmv(1)
    select case (turbulencemodel)
    case (1) !!spalart almaras model
      if (ispal .eq. 1) then
        eddyfl(2) = turbmv(1)
        eddyfr(2) = turbmv(2)
        onesix = 1.0d0/6.0d0
        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        cw1 = cb1/(kappa*kappa)
        cw1 = cw1 + (1.0 + cb2)/sigma
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .lt. verysmall) then
          source_t(1) = zero
        else
          tch_x3 = (tch_x)*(tch_x)*(tch_x)
          tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
          tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
          ddw = ielem(n, i)%walldist
          if (des_model .eq. 1) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            ddw = min(ddw, c_des_sa*delta_cell)
          end if
          if (des_model .eq. 2) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            r_des = min(10.0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
            !previous limiter is just for numerical reasons regarding tanh
            f_des = 1 - tanh((8.0*r_des)**3)
            ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
          end if
          prodterm1 = (turbmv(1))/(leftv(1)*kappa*kappa*ddw*ddw)
          stild = max(omega + (tch_fv2*prodterm1), 0.3*omega)
          prodtermfinal = stild*turbmv(1)*cb1
          rr = min((turbmv(1)/(((leftv(1)*stild*kappa*kappa*(ddw)*(ddw))) + 0.000000001)), 10.0)
          gg = rr + (cw2*(rr**6 - rr))
          fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
          !  destruction term
          destterm = cw1*fw*(((turbmv(1))/(leftv(1)*(ddw)))**2)
          ! !         !  first order diffusion term
          fodt = cb2*squaret/sigma
          source_t(1) = prodtermfinal + fodt - destterm
        end if
      else
        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .gt. 10.0d0) then
          tch_x = tch_x
        else
          tch_x = 0.05d0*log(1.0d0 + exp(20.0d0*(tch_x)))
        end if
        ddw = ielem(n, i)%walldist
        if (des_model .eq. 1) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          ddw = min(ddw, c_des_sa*delta_cell)
        end if
        if (des_model .eq. 2) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          r_des = min(10.0d0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
          !previous limiter is just for numerical reasons regarding tanh
          f_des = 1 - tanh((8.0d0*r_des)**3)
          ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
        end if
        tch_x3 = (tch_x)*(tch_x)*(tch_x)
        tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
        tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
        prodterm1 = (tch_fv2*tch_x*(turbmv(1)))/(leftv(1)*kappa*kappa*(ddw)*(ddw))
        if (prodterm1 .ge. (-0.7d0*omega)) then
          stild = omega + prodterm1
        end if
        if (prodterm1 .lt. (-0.7d0*omega)) then
          stild = omega + ((omega*(((0.7d0*0.7d0)*(omega)) + (0.9*prodterm1)))/(((0.9 - 1.4)*omega) - prodterm1))
        end if
        prodtermfinal = stild*turbmv(1)*cb1*tch_x
        rr = ((turbmv(1)*tch_x/(leftv(1)*stild*kappa*kappa*(ddw)*(ddw))))
        gg = rr + (cw2*(rr**6 - rr))
        fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
        ! !          destruction term
        destterm = cw1*fw*leftv(1)*(((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
        ! !          first order diffusion term
        fodt = leftv(1)*cb2*squaret/sigma
        source_t(1) = prodtermfinal + fodt - destterm
      end if

    case (2)                !k omega sst
      eddyfl(1) = ielem(n, i)%walldist
      eddyfl(2) = u_ct(i)%val(1, 1)
      eddyfl(3) = u_ct(i)%val(1, 2)
      eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2)
      eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
      eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
      eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)
      eddyfr = eddyfl
      call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      k_0 = (max(verysmall, u_ct(i)%val(1, 1)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0 = max(1.0e-1*ufreestream/charlength, u_ct(i)%val(1, 2)/leftv(1))
      wally = ielem(n, i)%walldist
      !calculate here k and omega gradients
      omx = ilocal_recon3(i)%grads(5, 1)
      omy = ilocal_recon3(i)%grads(5, 2)
      kx = ilocal_recon3(i)%grads(4, 1)
      ky = ilocal_recon3(i)%grads(4, 2)
      dervk_dervom = kx*omx + ky*omy
      !parameters
      !-----needed in diffusion--------
      d_omplus = max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
      phi_2 = max(sqrt(k_0)/(0.09*om_0*wally), 500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
      phi_1 = min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally))
      f_1 = tanh(phi_1**4)
      f_2 = tanh(phi_2**2)
      alpha_inf = f_1*alpha_inf1 + (1.0 - f_1)*alpha_inf2
      alpha_star = alpha_starinf
      alpha_raw = alpha_inf
      beta_i = f_1*beta_i1 + (1.0 - f_1)*beta_i2
      alpha_star0 = beta_i/3.0
      beta_stari = beta_starinf
      lowre = 0
      !low-re correction------------------------------------------------------------------
      if (lowre .eq. 1) then
        re_t_sst = leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???
        alpha_star = alpha_starinf*(alpha_star0 + re_t_sst/r_k_sst)/(1.0 + re_t_sst/r_k_sst)
        alpha_raw = alpha_inf/alpha_star*(alpha_0 + re_t_sst/r_om_sst)/(1.0 + re_t_sst/r_om_sst)
        beta_stari = beta_starinf*(4.0/15.0 + (re_t_sst/r_beta)**4)/(1.0 + (re_t_sst/r_beta)**4)
      end if
      !no mach number corrections for the beta
      beta_star = beta_stari
      beta_raw = beta_i
      !production terms
      !production of k
      if (vort_model .eq. 1) then
        prod_k = viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
        prod_om = alpha_raw*leftv(1)*snorm*onorm
      else
        prod_k = viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
        prod_k = min(prod_k, 10.0*leftv(1)*beta_star*k_0*om_0)
      end if
      !destruction terms
      !destruction of k
      ydest_k = leftv(1)*beta_star*k_0*om_0
      !destruction of omega
      ydest_om = leftv(1)*beta_raw*om_0*om_0
      !crossed-diffusion term
      !crossed diffusion of omega
      diff_om = 2.0*(1 - f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
      !qsas term: scale adaptive
      if (qsas_model .eq. 1) then  
        !calculate here second derivative of u
        !declare all variables
        do iex = 1, 2
          vortet(iex, 1:2) = ilocal_recon3(i)%grads(3 + turbulenceequations + iex, 1:2)
        end do
        uxx = vortet(1, 1); uyy = vortet(1, 2); 
        vxx = vortet(2, 1); vyy = vortet(2, 2); 
        !compute here cell volume
        cell_volume = ielem(n, i)%totvolume
        u_lapl = sqrt((uxx + uyy + uzz)**2 + (vxx + vyy + vzz)**2 + (wxx + wyy + wzz)**2)
        dervk2 = kx*kx + ky*ky + kz*kz
        dervom2 = omx*omx + omy*omy + omz*omz
        delta_cell = cell_volume**0.333333333333333
        l_sas = sqrt(k_0)/(beta_star**0.25*om_0)
        l_vk = max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
                   c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star - alpha_raw)))
        q_sas1 = leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
        q_sas2 = -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2, dervom2/om_0**2)
        q_sas = max(q_sas1 + q_sas2, 0.0)
      else
        q_sas = zero
      end if
      !final source terms
      !des-sst model (if qsas_model=2)
      if (qsas_model .eq. 2) then
        cell_volume = ielem(n, i)%totvolume
        delta_cell = cell_volume**0.333333333333333
        l_t_des = sqrt(k_0)/(beta_star*om_0)
        f_des_sst = max(1.0, l_t_des/(c_des_sst*delta_cell)*(1 - f_2))
        !the (1-f_2) is meant to protect the boundary layer. will result in same
        !separation point that standard s-a
        ydest_k = f_des_sst*ydest_k
      end if
      srcfull_k = prod_k - ydest_k
      srcfull_om = prod_om - ydest_om + diff_om + q_sas
      !filling the output vector
      source_t(1) = srcfull_k
      source_t(2) = srcfull_om
    end select
  end subroutine sources2d

  subroutine sources_derivatives2d(n, iconsidered, source_t)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(turbulenceequations), intent(inout)::source_t
    real::intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
    real::vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
    real::tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
    integer::i, k, j, l, ihgt, ihgj, iex, lowre
    real::snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
    real::gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
    real:: r_des, f_des, ddw, f_des_sst, l_t_des
    real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
    real:: uz, vz, wx, wy, wz
    real:: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz  !for sas only
    real, dimension(2, 2)::vortet, tvort, svort, ovort
    real, dimension(2)::vortem, vertf, vertex
    real, dimension(turbulenceequations, 1:2)::derivturb
    real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
    real:: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
    real:: alpha_raw, alpha_star, re_t_sst, alpha_inf
    real:: beta_stari, beta_i, beta_raw, beta_star
    real:: k_0, om_0, wally
    real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
    real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
    real:: kx, ky, kz, omx, omy, omz
    real, dimension(1:nof_variables)::leftv, rightv
    real::mp_pinfl, gammal
    real, dimension(1:4)::viscl, laml
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    i = iconsidered
    verysmall = 10e-16
    vortet(1:2, 1:2) = ilocal_recon3(i)%grads(1:2, 1:2)
    ux = vortet(1, 1); uy = vortet(1, 2)
    vx = vortet(2, 1); vy = vortet(2, 2)
    do ihgt = 1, 2
      do ihgj = 1, 2
        tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
      end do
    end do
    svort = 0.5*(vortet + tvort)
    ovort = 0.5*(vortet - tvort)
    snorm = sqrt(2.0d0*((svort(1, 1)*svort(1, 1)) + (svort(1, 2)*svort(1, 2)) + &
                        (svort(2, 1)*svort(2, 1)) + (svort(2, 2)*svort(2, 2))))
    onorm = sqrt(2.0d0*((ovort(1, 1)*ovort(1, 1)) + (ovort(1, 2)*ovort(1, 2)) + &
                        (ovort(2, 1)*ovort(2, 1)) + (ovort(2, 2)*ovort(2, 2))))
    omega = onorm
    divnorm = ux + uy !careful with the sign. if it becomes very big, it can produce negative production
    usss = sqrt((2.0*((ux*ux) + (vy*vy))) &
                + ((uy + vx)*(uy + vx)) &
                - (2.0/3.0*(ux + vy)*(ux + vy)))
    derivturb(1, 1:2) = ilocal_recon3(i)%grads(4, 1:2)
    squaret = (sqrt((derivturb(1, 1)**2) + (derivturb(1, 2)**2)))**2
    leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    rightv(1:nof_variables) = leftv(1:nof_variables)
    call sutherland2d(n, leftv, rightv, viscl, laml)
    turbmv(1) = u_ct(i)%val(1, 1)
    turbmv(2) = turbmv(1)
    select case (turbulencemodel)

    case (1) !!spalart almaras model
      if (ispal .eq. 1) then
        eddyfl(2) = turbmv(1)
        eddyfr(2) = turbmv(2)
        onesix = 1.0d0/6.0d0
        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        cw1 = cb1/(kappa*kappa)
        cw1 = cw1 + (1.0 + cb2)/sigma
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .lt. verysmall) then
          source_t(1) = zero
        else
          tch_x3 = (tch_x)*(tch_x)*(tch_x)
          tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
          tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
          ddw = ielem(n, i)%walldist
          if (des_model .eq. 1) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            ddw = min(ddw, c_des_sa*delta_cell)
          end if
          if (des_model .eq. 2) then
            cell_volume = ielem(n, i)%totvolume
            delta_cell = cell_volume**0.333333333333333
            r_des = min(10.0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
            !previous limiter is just for numerical reasons regarding tanh
            f_des = 1 - tanh((8.0*r_des)**3)
            ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
          end if
          prodterm1 = (turbmv(1))/(leftv(1)*kappa*kappa*ddw*ddw)
          stild = max(omega + (tch_fv2*prodterm1), 0.3*omega)
          prodtermfinal = stild*turbmv(1)*cb1/leftv(1)
          rr = min((turbmv(1)/(((leftv(1)*stild*kappa*kappa*(ddw)*(ddw))) + 0.000000001)), 10.0)
          gg = rr + (cw2*(rr**6 - rr))
          fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
          !  destruction term
          destterm = cw1*fw*(((turbmv(1))/(leftv(1)*(ddw)))**2)
          !  first order diffusion term
          fodt = cb2*squaret/sigma
          prodtermfinal = stild*cb1
          destterm = 2.0*cw1*fw*(turbmv(1)/(leftv(1)*(kappa**2)*(ddw**2)))
          fodt = -2.0*cb2*(sqrt(squaret))
        end if
      else
        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
        tch_x = turbmv(1)/viscl(1)
        if (tch_x .gt. 10.0d0) then
          tch_x = tch_x
        else
          tch_x = 0.05d0*log(1.0d0 + exp(20.0d0*(tch_x)))
        end if
        ddw = ielem(n, i)%walldist
        if (des_model .eq. 1) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          ddw = min(ddw, c_des_sa*delta_cell)
        end if
        if (des_model .eq. 2) then
          cell_volume = ielem(n, i)%totvolume
          delta_cell = cell_volume**0.333333333333333
          r_des = min(10.0d0, (viscl(1) + viscl(3))/(snorm*(kappa*ddw)**2 + 1e-16))
          !previous limiter is just for numerical reasons regarding tanh
          f_des = 1 - tanh((8.0d0*r_des)**3)
          ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1e-16), 10.0e-16)
        end if
        tch_x3 = (tch_x)*(tch_x)*(tch_x)
        tch_fv1 = tch_x3/(tch_x3 + (cv1*cv1*cv1))
        tch_fv2 = 1.0d0 - (tch_x/(1.0d0 + tch_x*tch_fv1))
        prodterm1 = (tch_fv2*tch_x*(turbmv(1)))/(leftv(1)*kappa*kappa*(ddw)*(ddw))
        if (prodterm1 .ge. (-0.7d0*omega)) then
          stild = omega + prodterm1
        end if
        if (prodterm1 .lt. (-0.7d0*omega)) then
          stild = omega + ((omega*(((0.7d0*0.7d0)*(omega)) + (0.9*prodterm1)))/(((0.9 - 1.4)*omega) - prodterm1))
        end if
        prodtermfinal = stild*turbmv(1)*cb1*tch_x
        rr = ((turbmv(1)*tch_x/(leftv(1)*stild*kappa*kappa*(ddw)*(ddw))))
        gg = rr + (cw2*(rr**6 - rr))
        fw = gg*(((1.0 + cw3**6)/(gg**6 + cw3**6))**onesix)
        ! destruction term
        destterm = cw1*fw*leftv(1)*(((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
        ! first order diffusion term
        fodt = leftv(1)*cb2*squaret/sigma
        prodtermfinal = stild*cb1*tch_x
        destterm = 2.0*cw1*fw*(turbmv(1)/(leftv(1)*ddw**2))
        fodt = 2.0*cb2*(sqrt(squaret))/sigma
        source_t(1) = min(prodtermfinal + fodt - destterm, zero)
      end if

    case (2)                !k omega sst
      eddyfl(1) = ielem(n, i)%walldist
      eddyfl(2) = u_ct(i)%val(1, 1)
      eddyfl(3) = u_ct(i)%val(1, 2)
      eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2)
      eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
      eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
      eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)
      eddyfr = eddyfl
      call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      k_0 = (max(verysmall, u_ct(i)%val(1, 1)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0 = max(1.0e-1*ufreestream/charlength, u_ct(i)%val(1, 2)/leftv(1))
      !filling the output vector
      source_t(1) = beta_starinf*om_0
      source_t(2) = beta_starinf*k_0
    end select
  end subroutine sources_derivatives2d
end module source