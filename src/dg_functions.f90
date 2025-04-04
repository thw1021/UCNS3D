module dg_functions
  use basis
  use declaration
  use derivatives
  use lapck
  implicit none
  contains
  function dg_sol(n, iconsidered, x1, y1, z1)
    implicit none
    !> @brief
    !> this function returns the dg solution at a given point (x_in, y_in)\n
    !> requires: x_in, y_in: coordinates of the point where the solution is requested, num_variables: number of solution variables, num_dofs: number of basis terms
    real, allocatable, dimension(:)::basis_temp
    integer, intent(in)::n, iconsidered
    integer::i_dof, i_var, icompwrt, number_of_dog, number
    real, dimension(1:nof_variables)::dg_sol
    real, intent(in)::x1, y1, z1
    icompwrt = -2
    number = ielem(n, iconsidered)%iorder
    number_of_dog = ielem(n, iconsidered)%idegfree
    allocate(basis_temp(1:number_of_dog))
    if (dimensiona .eq. 2) then
      basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
    else
      basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
    end if
    do i_var = 1, nof_variables
      dg_sol(i_var) = u_c(iconsidered)%valdg(1,i_var,1)+ dot_product(basis_temp(1:number_of_dog),u_c(iconsidered)%valdg(1,i_var,2:number_of_dog+1))
    end do
    icompwrt = 0
    deallocate(basis_temp)
  end function dg_sol

  function dg_solface(n, facex, pointx, iconsidered, number_of_dog)
    implicit none
    !> @brief
    !> this function returns the dg solution at a given surface point (x_in, y_in)\n
    !> requires: x_in, y_in: coordinates of the point where the solution is requested, num_variables: number of solution variables, num_dofs: number of basis terms
    real, allocatable, dimension(:)::basis_temp
    integer::i_dof, i_var, icompwrt
    integer, intent(in)::n, facex, pointx, iconsidered, number_of_dog
    real, dimension(1:nof_variables)::dg_solface
    real::x1, y1, z1
    integer::number
    icompwrt = -2
    x1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 1)
    y1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 2)
    if (dimensiona .eq. 3) then
      z1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 3)
    end if
    allocate(basis_temp(1:number_of_dog))
    number = ielem(n, iconsidered)%iorder
    if (dimensiona .eq. 2) then
      basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
    else
      basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
    end if
    do i_var = 1, nof_variables
      dg_solface(i_var) = u_c(iconsidered)%valdg(1,i_var,1) + dot_product(basis_temp(1:number_of_dog),u_c(iconsidered)%valdg(1,i_var,2:number_of_dog+1))
    end do
    deallocate(basis_temp)
    icompwrt = 0
  end function dg_solface

  function dg_sol_der(x1, y1, z1, number_of_dog, iconsidered)
    implicit none
    !> @brief
    !> this function returns the derivative of the dg solution at a given point (x1, y1, [z1])
    integer::i_dof, i_var, i_dim, icompwrt, number
    integer, intent(in)::number_of_dog, iconsidered
    real, intent(in)::x1, y1, z1
    real, allocatable, dimension(:, :)::basis_temp
    real, dimension(1:nof_variables, 1:dimensiona)::dg_sol_der
    number = ielem(n, iconsidered)%iorder
    allocate(basis_temp(number_of_dog, dimensiona))
    if (br2_yn .ne. 1) then
      do i_dof = 1, number_of_dog
        if (dimensiona .eq. 2) then
          if (poly .eq. 1) then
            basis_temp(i_dof, 1) = df2dx(x1, y1, i_dof, iconsidered)
            basis_temp(i_dof, 2) = df2dy(x1, y1, i_dof, iconsidered)
          else if (poly .eq. 4) then
            basis_temp(i_dof, 1) = tl2dx(x1, y1, i_dof, iconsidered)
            basis_temp(i_dof, 2) = tl2dy(x1, y1, i_dof, iconsidered)
          end if
        else
          basis_temp(i_dof, 1) = tl3dx(x1, y1, z1, i_dof, iconsidered)
          basis_temp(i_dof, 2) = tl3dy(x1, y1, z1, i_dof, iconsidered)
          basis_temp(i_dof, 3) = tl3dz(x1, y1, z1, i_dof, iconsidered)
        end if
      end do
    end if

    icompwrt = -2
    do i_var = 1, nof_variables
      do i_dim = 1, dimensiona
        if (br2_yn .eq. 1) then
          if (dimensiona .eq. 2) then
            dg_sol_der(i_var, i_dim) = u_c(iconsidered)%br2_aux_var(1,i_var,i_dim) + dot_product(basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt), u_c(iconsidered)%br2_aux_var(2:num_dg_dofs,i_var,i_dim))
          else
            dg_sol_der(i_var, i_dim) = u_c(iconsidered)%br2_aux_var(1,i_var,i_dim) + dot_product(basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog,icompwrt), u_c(iconsidered)%br2_aux_var(2:num_dg_dofs,i_var,i_dim))
          end if
        else
          dg_sol_der(i_var, i_dim) = dot_product(basis_temp(1:number_of_dog, i_dim), u_c(iconsidered)%valdg(1, i_var, 2:num_dg_dofs))
        end if
      end do
    end do
    icompwrt = 0
    deallocate(basis_temp)
  end function dg_sol_der

  function br2_local_lift(n, n_qp, facex, iconsidered, wequa2d)
    implicit none
    !> requires: iconsidered, facex, and weights_dg or wequa2d
    integer, intent(in)::n, n_qp, facex, iconsidered
    real, dimension(1:numberofpoints2), intent(in)::wequa2d
    integer::i, l, ngp, iqp, i_dim, icompwrt, number_of_dog
    real::angle1, angle2
    real, dimension(dimensiona)::nnn
    real, dimension(numberofpoints2)::weights_q, weights_t, weights_temp
    real, dimension(nof_variables, dimensiona)::br2_local_lift
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv, srf_speedrot
    real::mp_pinfr, gammar
    real, allocatable, dimension(:)::basis_temp
    real::x1, y1, z1
    integer::b_code, pointx, number
    real::nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz

    i = iconsidered
    l = facex
    number_of_dog = ielem(n, i)%idegfree
    number = ielem(n, iconsidered)%iorder
    allocate(basis_temp(1:number_of_dog))
    if (dimensiona .eq. 3) then
      angle1 = ielem(n, i)%faceanglex(l)
      angle2 = ielem(n, i)%faceangley(l)
      nnn(1) = (cos(angle1)*sin(angle2))
      nnn(2) = (sin(angle1)*sin(angle2))
      nnn(3) = (cos(angle2))
    else
      nnn(1) = ielem(n, i)%faceanglex(l)
      nnn(2) = ielem(n, i)%faceangley(l)
    end if
    br2_local_lift = 0.0d0
    do ngp = 1, n_qp
      pointx = ngp
      icompwrt = -2
      x1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 1)
      y1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 2)
      if (dimensiona .eq. 3) z1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 3)
        if (dimensiona .eq. 2) then
          basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
        else
          basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
        end if
        call get_left_right_states(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
        leftv = cleft
        rightv = cright
      if (dimensiona .eq. 3) then
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
      else
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
      end if
      ! change pressure to temperature
      leftv(nof_variables) = leftv(nof_variables)/leftv(1)
      rightv(nof_variables) = rightv(nof_variables)/rightv(1)

      do i_dim = 1, dimensiona
        br2_local_lift(1:nof_variables, i_dim) = br2_local_lift(1:nof_variables, i_dim) + nnn(i_dim)*(rightv - leftv)*weights_temp(pointx)
      end do
    end do
    br2_local_lift = oo2*br2_local_lift*ielem(n, iconsidered)%surf(facex)/ielem(n, iconsidered)%totvolume
    deallocate(basis_temp)
  end function br2_local_lift

  function dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
    !> @brief
    !> calculates the rhs flux term to be integrated in the dg formulation
    implicit none
    real, dimension(idegfree + 1, nof_variables)::dg_surf_flux
    integer::i, icompwrt
    integer, intent(in)::n, iconsidered, facex, pointx
    real, dimension(1:nof_variables), intent(in)::rhllcflux
    real, dimension(1:numberofpoints2), intent(in)::weights_temp
    real::x1, y1, z1
    integer::number, number_of_dog
    number_of_dog = ielem(n, iconsidered)%idegfree
    number = ielem(n, iconsidered)%iorder
    icompwrt = -2
    x1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 1)
    y1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 2)
    if (dimensiona .eq. 3) then
      z1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 3)
    end if
      number = ielem(n, iconsidered)%iorder
    if (dimensiona .eq. 2) then
      do i = 1, nof_variables
        dg_surf_flux(1, i) = rhllcflux(i)*weights_temp(pointx)*ielem(n, iconsidered)%surf(facex)
        dg_surf_flux(2:number_of_dog+1,i) = rhllcflux(i) * weights_temp(pointx)*ielem(n,iconsidered)%surf(facex)*basis_rec2d(n,x1,y1,number,iconsidered,number_of_dog,icompwrt)
      end do
    else
      do i = 1, nof_variables
        dg_surf_flux(1, i) = rhllcflux(i)*weights_temp(pointx)*ielem(n, iconsidered)%surf(facex)
        dg_surf_flux(2:number_of_dog+1,i) = rhllcflux(i) * weights_temp(pointx)*ielem(n,iconsidered)%surf(facex)*basis_rec(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt)
      end do
    end if
    icompwrt = 0
  end function dg_surf_flux

  function dg_surf_fluxv(n, iconsidered, facex, pointx, weights_temp, hllcflux)
    !> @brief
    !> calculates the rhs flux term to be integrated in the dg formulation
    implicit none
    real, dimension(idegfree + 1, 1:nof_variables)::dg_surf_fluxv
    integer::i, icompwrt
    integer, intent(in)::n, iconsidered, facex, pointx
    real, dimension(1:nof_variables), intent(in)::hllcflux
    real, dimension(1:nof_variables)::rhllcflux
    real, dimension(1:numberofpoints2), intent(in)::weights_temp
    real::x1, y1, z1
    integer::number, number_of_dog
    icompwrt = -2
    x1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 1)
    y1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 2)
    if (dimensiona .eq. 3) then
      z1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 3)
    end if
    number = ielem(n, iconsidered)%iorder
    if (dimensiona .eq. 2) then
      do i = 1, nof_variables
        dg_surf_fluxv(1, i) = hllcflux(i)*weights_temp(pointx)*ielem(n, iconsidered)%surf(facex)
        dg_surf_fluxv(2:number_of_dog+1,i) = hllcflux(i) * weights_temp(pointx)*ielem(n,iconsidered)%surf(facex)*basis_rec2d(n,x1,y1,number,iconsidered,number_of_dog,icompwrt)
      end do
    else
      do i = 1, nof_variables
        dg_surf_fluxv(1, i) = hllcflux(i)*weights_temp(pointx)*ielem(n, iconsidered)%surf(facex)
        dg_surf_fluxv(2:number_of_dog+1,i) = hllcflux(i) * weights_temp(pointx)*ielem(n,iconsidered)%surf(facex)*basis_rec(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt)
      end do
    end if
      icompwrt = 0
  end function dg_surf_fluxv
  function dg_vol_integral(n, iconsidered)
    !> @brief
    !> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i, nqp, i_qp, icompwrt, number, number_of_dog
    real, dimension(idegfree + 1, nof_variables)::dg_vol_integral
    real, dimension(1:nof_variables)::leftv
    real, dimension(1:nof_variables)::flux_term_x, flux_term_y, flux_term_z
    real, dimension(1:nof_variables, 1:dimensiona)::leftv_der
    real::x1, y1, z1
    real::mp_pinfl, gammal
    icompwrt = -2
    nqp = ielem(n, iconsidered)%itotalpoints
    number = ielem(n, iconsidered)%iorder
    number_of_dog = ielem(n, iconsidered)%idegfree
    dg_vol_integral(:, :) = 0.0d0
    do i = 1, number_of_dog
      do i_qp = 1, nqp
        x1 = qp_array(iconsidered)%x(i_qp)
        y1 = qp_array(iconsidered)%y(i_qp)
        if (dimensiona .eq. 3) then
          z1 = qp_array(iconsidered)%z(i_qp)
        end if
        if (itestcase .lt. 3) then ! linear advection
          flux_term_x = dg_sol(n, iconsidered, x1, y1, z1)*lamx !flux in x-axis of linear advection in 2d (lamx*sol), where lamx is the wave speed for x axis, and sol is the solution
          flux_term_y = dg_sol(n, iconsidered, x1, y1, z1)*lamy !flux in y-axis of linear advection in 2d (lamy*sol), where lamy is the wave speed for y axis, and sol is the solution
          if (dimensiona .eq. 3) then
            flux_term_z = dg_sol(n, iconsidered, x1, y1, z1)*lamz
          end if
        end if
        if (itestcase .eq. 3) then ! euler
          leftv = dg_sol(n, iconsidered, x1, y1, z1)
          if (dimensiona .eq. 2) then
            call cons2prim(n, leftv, mp_pinfl, gammal)
            call flux2dx(flux_term_x, leftv)
            call flux2dy(flux_term_y, leftv)
          else
            call cons2prim(n, leftv, mp_pinfl, gammal)
            call flux3dx(flux_term_x, leftv)
            call flux3dy(flux_term_y, leftv)
            call flux3dz(flux_term_z, leftv)
          end if
        end if
        if (itestcase .eq. 4) then ! ns
          leftv = dg_sol(n, iconsidered, x1, y1, z1)
          leftv_der = dg_sol_der(x1, y1, z1, number_of_dog, iconsidered)
          if (br2_yn .eq. 2) then
            call dcons2dprim(leftv_der, leftv)
            leftv_der = leftv_der + sum(ilocal_recon3(iconsidered)%br2_local_lift, dim=3)
          end if
          if (dimensiona .eq. 2) then
            call cons2prim(n, leftv, mp_pinfl, gammal)
            call flux2dx(flux_term_x, leftv)
            call flux2dy(flux_term_y, leftv)
            call flux_visc2d(flux_term_x, flux_term_y, leftv, leftv_der)
          else
            call cons2prim(n, leftv, mp_pinfl, gammal)
            call flux3dx(flux_term_x, leftv)
            call flux3dy(flux_term_y, leftv)
            call flux3dz(flux_term_z, leftv)
            call flux_visc3d(flux_term_x, flux_term_y, flux_term_z, leftv, leftv_der)
          end if
        end if
        if (dimensiona .eq. 2) then
          if (poly .eq. 1) then
            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array(iconsidered)%qp_weight(i_qp) *(flux_term_x(:)*df2dx(x1,y1,i,iconsidered)+flux_term_y(:)*df2dy(x1,y1,i,iconsidered))
          end if
          if (poly .eq. 4) then
            dg_vol_integral(i+1,1:nof_variables) = dg_vol_integral(i+1,1:nof_variables)+ qp_array(iconsidered)%qp_weight(i_qp)*(flux_term_x(1:nof_variables)*tl2dx(x1,y1,i,iconsidered)+flux_term_y(1:nof_variables)*tl2dy(x1,y1,i,iconsidered))
          end if
        else
          if (poly .eq. 1) then
            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array(iconsidered)%qp_weight(i_qp) *((flux_term_x(:)*dfx(x1,y1,z1,i,iconsidered))&+ (flux_term_y(:)*dfy(x1, y1, z1, i, iconsidered)) + (flux_term_z(:)*dfz(x1, y1, z1, i, iconsidered)))
          end if
          if (poly .eq. 2) then
            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array(iconsidered)%qp_weight(i_qp) *((flux_term_x(:)*dlx(x1,y1,z1,i,iconsidered))&+ (flux_term_y(:)*dly(x1, y1, z1, i, iconsidered)) + (flux_term_z(:)*dlz(x1, y1, z1, i, iconsidered)))
          end if
          if (poly .eq. 4) then
            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array(iconsidered)%qp_weight(i_qp)*((flux_term_x(:)*tl3dx(x1,y1,z1,i,iconsidered))&+ (flux_term_y(:)*tl3dy(x1, y1, z1, i, iconsidered)) + (flux_term_z(:)*tl3dz(x1, y1, z1, i, iconsidered)))
          end if
        end if
      end do
    end do
  end function dg_vol_integral
  function dg_vol_integral2(n, i)
    !> @brief
    !> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
    implicit none
    integer, intent(in)::n, i
    real, dimension(1:nof_variables)::dg_vol_integral2
    integer::j, k, nqp, i_qp, i_var, icompwrt, iconsidered
    real::ph, integ
    real::x1, y1, z1
    integer::number, number_of_dog
    real, dimension(1:nof_variables)::dg_sol2
    real, allocatable, dimension(:)::basis_temp
    allocate(basis_temp(1:idegfree))
    iconsidered = i
    icompwrt = -2
    nqp = ielem(n, iconsidered)%itotalpoints
    number = ielem(n, iconsidered)%iorder
    number_of_dog = ielem(n, iconsidered)%idegfree
    dg_vol_integral2(:) = 0.0d0
    do i_qp = 1, nqp
      x1 = qp_array(iconsidered)%x(i_qp)
      y1 = qp_array(iconsidered)%y(i_qp)
      if (dimensiona .eq. 3) then
        z1 = qp_array(iconsidered)%z(i_qp)
      end if
      if (dimensiona .eq. 2) then
        basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
      else
        basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
      end if
      do i_var = 1, nof_variables
        dg_sol2(i_var) = dot_product(basis_temp(1:number_of_dog), u_c(iconsidered)%valdg(1, i_var, 2:number_of_dog + 1))
      end do
        dg_vol_integral2 = dg_vol_integral2 + (qp_array(iconsidered)%qp_weight(i_qp)*dg_sol2/ielem(n, iconsidered)%totvolume)
    end do
    dg_vol_integral2(:) = u_c(iconsidered)%valdg(1, :, 1) + dg_vol_integral2(:)
    icompwrt = 0
    deallocate(basis_temp)
  end function dg_vol_integral2

  function dg_vol_integral_strong(n, i)
  !> @brief
  !> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
    implicit none
    integer, intent(in)::n, i
    real, dimension(1:nof_variables)::dg_vol_integral_strong
    integer::j, k, nqp, i_qp, i_var, icompwrt, iconsidered
    real::ph, integ
    real::x1, y1, z1
      integer::number, number_of_dog
      real, dimension(1:nof_variables)::dg_sol2
      real, allocatable, dimension(:)::basis_temp
      allocate(basis_temp(1:idegfree))
      iconsidered = i
      icompwrt = -2
      do j = 1, nof_variables
        do k = 1, idegfree
          u_cs(iconsidered)%valdg(1, j, k + 1) = u_c(iconsidered)%valdg(1, j, k + 1)*modal_filter_strong(k)
        end do
      end do
      nqp = ielem(n, iconsidered)%itotalpoints
      number = ielem(n, iconsidered)%iorder
      number_of_dog = ielem(n, iconsidered)%idegfree
      dg_vol_integral_strong(:) = 0.0d0
      do i_qp = 1, nqp
        x1 = qp_array(iconsidered)%x(i_qp)
        y1 = qp_array(iconsidered)%y(i_qp)
        if (dimensiona .eq. 3) then
          z1 = qp_array(iconsidered)%z(i_qp)
        end if
        if (dimensiona .eq. 2) then
          basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
        else
          basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
        end if
        do i_var = 1, nof_variables
          dg_sol2(i_var) = dot_product(basis_temp(1:number_of_dog), u_cs(iconsidered)%valdg(1, i_var, 2:number_of_dog + 1))
        end do
      dg_vol_integral_strong = dg_vol_integral_strong + (qp_array(iconsidered)%qp_weight(i_qp)*dg_sol2/ielem(n, iconsidered)%totvolume)
      end do
      dg_vol_integral_strong(:) = u_c(iconsidered)%valdg(1, :, 1) + dg_vol_integral_strong(:)
      icompwrt = 0
      deallocate(basis_temp)
  end function dg_vol_integral_strong

  function dg_vol_integral_weak(n, i)
  !> @brief
  !> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
    implicit none
    integer, intent(in)::n, i
    real, dimension(1:nof_variables)::dg_vol_integral_weak
    integer::j, k, nqp, i_qp, i_var, icompwrt, iconsidered
    real::ph, integ
    real::x1, y1, z1
    integer::number, number_of_dog
    real, dimension(1:nof_variables)::dg_sol2
    real, allocatable, dimension(:)::basis_temp
    allocate(basis_temp(1:idegfree))
    iconsidered = i
    icompwrt = -2
    do j = 1, nof_variables
      do k = 1, idegfree
        u_cw(iconsidered)%valdg(1, j, k + 1) = u_c(iconsidered)%valdg(1, j, k + 1)*modal_filter_weak(k)
      end do
    end do
    nqp = ielem(n, iconsidered)%itotalpoints
    number = ielem(n, iconsidered)%iorder
    number_of_dog = ielem(n, iconsidered)%idegfree
    dg_vol_integral_weak(:) = 0.0d0
    do i_qp = 1, nqp
      x1 = qp_array(iconsidered)%x(i_qp)
      y1 = qp_array(iconsidered)%y(i_qp)
      if (dimensiona .eq. 3) then
        z1 = qp_array(iconsidered)%z(i_qp)
      end if
      if (dimensiona .eq. 2) then
        basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog, icompwrt)
      else
        basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog, icompwrt)
      end if
      do i_var = 1, nof_variables
        dg_sol2(i_var) = dot_product(basis_temp(1:number_of_dog), u_cw(iconsidered)%valdg(1, i_var, 2:number_of_dog + 1))
      end do
      dg_vol_integral_weak = dg_vol_integral_weak + (qp_array(iconsidered)%qp_weight(i_qp)*dg_sol2/ielem(n, iconsidered)%totvolume)
    end do
    dg_vol_integral_weak(:) = u_c(iconsidered)%valdg(1, :, 1) + dg_vol_integral_weak(:)
    icompwrt = 0

  end function dg_vol_integral_weak
  subroutine reconstruct_dg(n)
    implicit none
    integer, intent(in)::n
    integer::i_face, i_elem, i_qp, iqp, icompwrt
    integer::facex, pointx, iconsidered, number_of_dog
    !$omp do
    do i_elem = 1, xmpielrank(n)
      icompwrt = -2
      do i_face = 1, ielem(n, i_elem)%ifca
        !somewhere prestored the guassian quadrature points for your sides in another subroutine and you build only the cleft states and volume integral
        if (dimensiona .eq. 2) then
          iqp = qp_line_n
        else
          if (ielem(n, i_elem)%types_faces(i_face) .eq. 5) then
            iqp = qp_quad
          else
            iqp = qp_triangle
          end if
        end if
        do i_qp = 1, iqp! qp_line_n
          facex = i_face
          pointx = i_qp
          iconsidered = i_elem
          number_of_dog = ielem(n, i_elem)%idegfree
          ilocal_recon3(iconsidered)%uleft_dg(:, facex, pointx) = dg_solface(n, facex, pointx, iconsidered, number_of_dog)
        end do
      end do
      icompwrt = 0
    end do
  !$omp end do
  end subroutine reconstruct_dg
  subroutine reconstruct_br2_dg
    implicit none
    integer::i_face, i_elem, i_qp, iqp, j, k, iconsidered, facex, pointx, number_of_dog, number
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d, weights_q, weights_t, weights_l, weights_temp
    real, dimension(1:nof_variables, 1:dimensiona)::leftv_der
    real, dimension(1:nof_variables)::leftv
    real::x1, y1, z1
    ! needed for dg_surf_flux, used in br2_local_lift
    if (dimensiona .eq. 3) then
      call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
      weights_q(1:qp_quad) = wequa2d(1:qp_quad)
      call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
      weights_t(1:qp_triangle) = wequa2d(1:qp_triangle)
    else
      call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
      weights_l(1:qp_line_n) = wequa2d(1:qp_line_n)
    end if
    !$omp do
    do i_elem = 1, xmpielrank(n)
      iconsidered = i_elem
      do i_face = 1, ielem(n, i_elem)%ifca
        ! needed for br2_local_lift
        if (dimensiona .eq. 2) then
          iqp = qp_line_n
          weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
        else
          if (ielem(n, i_elem)%types_faces(i_face) .eq. 5) then
            iqp = qp_quad
            weights_temp(1:qp_quad) = weights_t(1:qp_quad)
          else
            iqp = qp_triangle
            weights_temp(1:qp_triangle) = weights_t(1:qp_triangle)
          end if
        end if
        number_of_dog = idegfree
        number = ielem(n, iconsidered)%iorder
        facex = i_face
        ilocal_recon3(iconsidered)%br2_local_lift(:, :, i_face) = br2_local_lift(n, iqp, facex, iconsidered, weights_temp)
        do i_qp = 1, iqp
          pointx = i_qp
          if (br2_yn .ne. 1 .and. itestcase .eq. 4) then
            x1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 1)
            y1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 2)
            if (dimensiona .eq. 3) z1 = ilocal_recon3(iconsidered)%qpoints(facex, pointx, 3)
              leftv_der = dg_sol_der(x1, y1, z1, number_of_dog, iconsidered)
              call dcons2dprim(leftv_der, leftv)
                ilocal_recon3(iconsidered)%br2_aux_var(:,:,i_face,i_qp) = leftv_der + ilocal_recon3(iconsidered)%br2_local_lift(:,:,i_face) * br2_damping
              !now copy to other array
                 ilocal_recon3(iconsidered)%uleftv(1:dims,1,i_face,i_qp)=ilocal_recon3(iconsidered)%br2_aux_var(nof_variables,1:dims,i_face,i_qp)

              do j = 2, nof_variables - 1
                ilocal_recon3(iconsidered)%uleftv(1:dims, j, i_face, i_qp) = ilocal_recon3(iconsidered)%br2_aux_var(j, 1:dims, i_face, i_qp)
              end do
              !end copying
            end if
        end do
      end do
    end do
!$omp end do
  end subroutine reconstruct_br2_dg

subroutine get_left_right_states(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    implicit none
    integer, intent(in)::iconsidered, facex, pointx, n
    integer::i, l, ngp
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    i = iconsidered
    if (dimensiona .eq. 3) then
      if (ielem(n, i)%interior .eq. 0) then
            call get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
      else
            call get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
      end if
    else
      if (ielem(n, i)%interior .eq. 0) then
            call get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
      else
            call get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
      end if
    end if
  end subroutine get_left_right_states

  subroutine  get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    implicit none
    integer, intent(in)::iconsidered, facex, pointx, n
    integer::i, l, ngp
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
    else !fv
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp) !right mean flow state
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
        cturbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
      end if
    end if
    if (ilocal_recon3(i)%mrf .eq. 1) then
      srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, ngp, 1:3)
      call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
    end if
  end subroutine get_states_interior

  subroutine  get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    implicit none
    integer, intent(in)::iconsidered, facex, pointx, n
    integer::i, l, ngp
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
    else !fv
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp) !right mean flow state
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
        cturbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
      end if
    end if
  end subroutine get_states_interior2d

  subroutine calculate_interior_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    integer, intent(in)::iconsidered, facex, pointx, n
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::i, l, ngp, ittt, nvar, iex
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
    else !fv
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp) !right mean flow state
    end if
    lcvgrad(1, 1:3) = ilocal_recon3(i)%uleftv(1:3, 2, l, ngp); lcvgrad(2, 1:3) = ilocal_recon3(i)%uleftv(1:3, 3, l, ngp);
    lcvgrad(3, 1:3) = ilocal_recon3(i)%uleftv(1:3, 4, l, ngp); lcvgrad(4, 1:3) = ilocal_recon3(i)%uleftv(1:3, 1, l, ngp);
    rcvgrad(1, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 2, ielem(n, i)%ineighn(l), ngp); rcvgrad(2, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 3, ielem(n, i)%ineighn(l), ngp);
    rcvgrad(3, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 4, ielem(n, i)%ineighn(l), ngp); rcvgrad(4, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 1, ielem(n, i)%ineighn(l), ngp);
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
        cturbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
      end if
      cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
      cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
      do nvar = 1, turbulenceequations + passivescalar
        rcvgrad_t(nvar, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:3, nvar, ielem(n, i)%ineighn(l), ngp)
        lcvgrad_t(nvar, 1:3) = ilocal_recon3(i)%uleftturbv(1:3, nvar, l, ngp)
      end do
    end if
    if (ilocal_recon3(i)%mrf .eq. 1) then
      srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, ngp, 1:3)
      call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
    end if
  end subroutine calculate_interior_viscous

  subroutine calculate_interior_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    integer, intent(in)::iconsidered, facex, pointx, n
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::i, l, ngp, ittt, nvar, iex
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
      lcvgrad(1, 1:2) = ilocal_recon3(i)%uleftv(1:2, 2, l, ngp); lcvgrad(2, 1:2) = ilocal_recon3(i)%uleftv(1:2, 3, l, ngp);
      lcvgrad(3, 1:2) = ilocal_recon3(i)%uleftv(1:2, 1, l, ngp)
      rcvgrad(1, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 2, ielem(n, i)%ineighn(l), ngp); rcvgrad(2, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 3, ielem(n, i)%ineighn(l), ngp);
      rcvgrad(3, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 1, ielem(n, i)%ineighn(l), ngp)
    else !fv
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
      lcvgrad(1, 1:2) = ilocal_recon3(i)%uleftv(1:2, 2, l, ngp); lcvgrad(2, 1:2) = ilocal_recon3(i)%uleftv(1:2, 3, l, ngp);
      lcvgrad(3, 1:2) = ilocal_recon3(i)%uleftv(1:2, 1, l, ngp)
      cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp) !right mean flow state
      rcvgrad(1, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 2, ielem(n, i)%ineighn(l), ngp); rcvgrad(2, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 3, ielem(n, i)%ineighn(l), ngp);
      rcvgrad(3, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 1, ielem(n, i)%ineighn(l), ngp)
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
        cturbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
      end if
      cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
      cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
      do nvar = 1, turbulenceequations + passivescalar
        rcvgrad_t(nvar, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:2, nvar, ielem(n, i)%ineighn(l), ngp)
        lcvgrad_t(nvar, 1:2) = ilocal_recon3(i)%uleftturbv(1:2, nvar, l, ngp)
      end do
    end if
  end subroutine calculate_interior_viscous2d

subroutine calculate_bounded_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    integer, intent(in)::iconsidered, facex, pointx, n
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::i, l, ngp, ittt, nvar, iex, kk, n_node
    integer::ibfc
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
    else
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
    end if
    lcvgrad(1, 1:3) = ilocal_recon3(i)%uleftv(1:3, 2, l, ngp); lcvgrad(2, 1:3) = ilocal_recon3(i)%uleftv(1:3, 3, l, ngp);
    lcvgrad(3, 1:3) = ilocal_recon3(i)%uleftv(1:3, 4, l, ngp); lcvgrad(4, 1:3) = ilocal_recon3(i)%uleftv(1:3, 1, l, ngp);
    if (ilocal_recon3(i)%mrf .eq. 1) then
      srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, ngp, 1:3)
      call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
      end if
      do nvar = 1, turbulenceequations + passivescalar
        lcvgrad_t(nvar, 1:3) = ilocal_recon3(i)%uleftturbv(1:3, nvar, l, ngp)
      end do
    end if
    if (ielem(n, i)%ineighb(l) .eq. n) then
      if (ielem(n, i)%ibounds(l) .gt. 0) then
        if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then
          if (dg .eq. 1) then
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          else
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          end if
          rcvgrad(1, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 2, l, ngp); rcvgrad(2, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 3, l, ngp);
          rcvgrad(3, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 4, l, ngp); rcvgrad(4, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 1, l, ngp);
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
            end if
            do nvar = 1, turbulenceequations + passivescalar
              rcvgrad_t(nvar, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:3, nvar, ielem(n, i)%ineighn(l), ngp)
            end do
          end if
          if (per_rot .eq. 1) then
            cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
            do kk = 1, 4
              rcvgrad(kk, 1:3) = rotate_per_1(rcvgrad(kk, 1:3), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
            end do
            rcvgrad_t(1, 1:3) = rotate_per_1(rcvgrad_t(1, 1:3), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
          end if
        else
          call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)
          cords(1:3) = zero
          if (ielem(n, iconsidered)%types_faces(facex) .eq. 5) then
            n_node = 4
          else
            n_node = 3
          end if
          cords(1:3) = cordinates3(n, nodes_list, n_node)
          poy(1) = cords(2)
          pox(1) = cords(1)
          poz(1) = cords(3)
          leftv(1:nof_variables) = cleft(1:nof_variables)
          b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
          call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
          cright(1:nof_variables) = rightv(1:nof_variables)
          rcvgrad(:, :) = lcvgrad(:, :)
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rcvgrad_t(:, :) = lcvgrad_t(:, :)
          end if
          if (b_code .eq. 4) then
            rightv = zero
            rightv(2:4) = lcvgrad(4, 1:3)
            leftv = zero
            call rotatef(n, leftv, rightv, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            leftv(2) = -leftv(2)
            call rotateb(n, rightv, leftv, angle1, angle2)
            rcvgrad(4, 1:3) = rightv(2:4)
          end if
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        else
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        end if
        rcvgrad(1, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 2, ielem(n, i)%ineighn(l), ngp); rcvgrad(2, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 3, ielem(n, i)%ineighn(l), ngp);
        rcvgrad(3, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 4, ielem(n, i)%ineighn(l), ngp); rcvgrad(4, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:3, 1, ielem(n, i)%ineighn(l), ngp);
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
          else
           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
          end if
          do nvar = 1, turbulenceequations + passivescalar
            rcvgrad_t(nvar, 1:3) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:3, nvar, ielem(n, i)%ineighn(l), ngp)
          end do
        end if
      end if
    else        !in other cpus they can only be periodic or mpi neighbours
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other
          if (dg .eq. 1) then
        cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          else
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          end if
          ittt = 0
          do iex = 1, nof_variables - 1
            do nvar = 1, dims
              ittt = ittt + 1
              if (dg .eq. 1) then
                rcvgrad(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol_dg(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              else
                if (iex .eq. 1) then
                  rcvgrad(4,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
                else
                  rcvgrad(iex-1,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
                end if
              end if
            end do
          end do
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            end if
          end if
          do iex = 1, turbulenceequations + passivescalar
            do nvar = 1, dims
              ittt = ittt + 1
                rcvgrad_t(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
            end do
          end do
          if (per_rot .eq. 1) then
            cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
            do kk = 1, 4
              rcvgrad(kk, 1:3) = rotate_per_1(rcvgrad(kk, 1:3), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
            end do
            rcvgrad_t(1, 1:3) = rotate_per_1(rcvgrad_t(1, 1:3), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
          end if
        end if
      else
        if (dg .eq. 1) then
        cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        else
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        end if
        ittt = 0
        do iex = 1, nof_variables - 1
          do nvar = 1, dims
            ittt = ittt + 1
            if (dg .eq. 1) then
              rcvgrad(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol_dg(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
            else
              if (iex .eq. 1) then
                rcvgrad(4,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              else
                rcvgrad(iex-1,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              end if
            end if
          end do
        end do
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          end if
          do iex = 1, turbulenceequations + passivescalar
            do nvar = 1, dims
              ittt = ittt + 1
              rcvgrad_t(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
            end do
          end do
        end if
      end if
    end if
  end subroutine calculate_bounded_viscous

  subroutine calculate_bounded_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    integer, intent(in)::iconsidered, facex, pointx, n
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::i, l, ngp, ittt, nvar, iex, n_node
    integer::ibfc
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
    else
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
    end if
    lcvgrad(1, 1:2) = ilocal_recon3(i)%uleftv(1:2, 2, l, ngp); lcvgrad(2, 1:2) = ilocal_recon3(i)%uleftv(1:2, 3, l, ngp);
    lcvgrad(3, 1:2) = ilocal_recon3(i)%uleftv(1:2, 1, l, ngp)
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
      end if
      do nvar = 1, turbulenceequations + passivescalar
        lcvgrad_t(nvar, 1:2) = ilocal_recon3(i)%uleftturbv(1:2, nvar, l, ngp)
      end do
    end if
    if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
          if (dg .eq. 1) then
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          else
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          end if
          rcvgrad(1, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 2, l, ngp); rcvgrad(2, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 3, l, ngp);
          rcvgrad(3, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 1, l, ngp);
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
            end if
            do nvar = 1, turbulenceequations + passivescalar
              rcvgrad_t(nvar, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:2, nvar, ielem(n, i)%ineighn(l), ngp)
            end do
          end if
        else
          !not periodic ones in my cpu
          call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
          n_node = 2
          cords(1:2) = zero
          cords(1:2) = cordinates2(n, nodes_list, n_node)
          poy(1) = cords(2)
          pox(1) = cords(1)
          leftv(1:nof_variables) = cleft(1:nof_variables)
          b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
          call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
          cright(1:nof_variables) = rightv(1:nof_variables)
          rcvgrad(:, :) = lcvgrad(:, :)
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rcvgrad_t(:, :) = lcvgrad_t(:, :)
          end if
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        else
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        end if
        rcvgrad(1, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 2, ielem(n, i)%ineighn(l), ngp); rcvgrad(2, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 3, ielem(n, i)%ineighn(l), ngp);
        rcvgrad(3, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftv(1:2, 1, ielem(n, i)%ineighn(l), ngp);
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
          else
           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
          end if
          do nvar = 1, turbulenceequations + passivescalar
            rcvgrad_t(nvar, 1:2) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturbv(1:2, nvar, ielem(n, i)%ineighn(l), ngp)
          end do
        end if
      end if
    else        !in other cpus they can only be periodic or mpi neighbours
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
          if (dg .eq. 1) then
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          else
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          end if
          ittt = 0
          do iex = 1, nof_variables - 1
            do nvar = 1, dims
              ittt = ittt + 1
              if (dg .eq. 1) then
                rcvgrad(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol_dg(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              else
                if (iex .eq. 1) then
                  rcvgrad(3,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
                else
                  rcvgrad(iex-1,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
                end if
              end if
            end do
          end do
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            end if
            do iex = 1, turbulenceequations + passivescalar
              do nvar = 1, dims
                ittt = ittt + 1
                rcvgrad_t(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              end do
            end do
          end if
        end if
      else
        if (dg .eq. 1) then
        cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        else
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        end if
        ittt = 0
        do iex = 1, nof_variables - 1
          do nvar = 1, dims
            ittt = ittt + 1
            if (dg .eq. 1) then
              rcvgrad(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol_dg(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
            else
              if (iex .eq. 1) then
                rcvgrad(3,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              else
                rcvgrad(iex-1,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
              end if
            end if
          end do
        end do
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          end if
          do iex = 1, turbulenceequations + passivescalar
            do nvar = 1, dims
              ittt = ittt + 1
              rcvgrad_t(iex,nvar)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(ngp),nof_variables+turbulenceequations+passivescalar+ittt)
            end do
          end do
        end if
      end if
    end if
  end subroutine calculate_bounded_viscous2d

  subroutine  get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    implicit none
    integer, intent(in)::iconsidered, facex, pointx, n
    integer::i, l, ngp, n_node
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv, srf_speedrot
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::ibfc
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
    else
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
    end if

    if (ilocal_recon3(i)%mrf .eq. 1) then
      srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, ngp, 1:3)
      call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
      end if
    end if
    if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then
          if (dg .eq. 1) then
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          else
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          end if
          if (per_rot .eq. 1) then
            cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
            end if
          end if
        else
          call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)
          if (ielem(n, iconsidered)%types_faces(facex) .eq. 5) then
            n_node = 4
          else
            n_node = 3
          end if
          cords(1:3) = zero
          cords(1:3) = cordinates3(n, nodes_list, n_node)
          poy(1) = cords(2)
          pox(1) = cords(1)
          poz(1) = cords(3)
          leftv(1:nof_variables) = cleft(1:nof_variables)
          b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
          call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
          cright(1:nof_variables) = rightv(1:nof_variables)
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        else
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
          end if
        end if
      end if
    else        !in other cpus they can only be periodic or mpi neighbours
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then
          if (dg .eq. 1) then
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          else
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          end if
          if (per_rot .eq. 1) then
            cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            end if
          end if
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        else
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          end if
        end if
      end if
    end if
  end subroutine get_states_bounds

  subroutine  get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    implicit none
    integer, intent(in)::iconsidered, facex, pointx, n
    integer::i, l, ngp, n_node
    integer, intent(inout)::b_code
    real, intent(inout)::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables), intent(inout)::cleft, cright, cright_rot, cleft_rot, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar), intent(inout)::cturbl, cturbr
    real, dimension(1:nof_variables), intent(inout)::leftv
    real, dimension(1:nof_variables), intent(inout)::rightv
    real, dimension(1:dimensiona), intent(inout)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    real, dimension(1:6, 1:4, 1:dimensiona)::elem_listd
    integer::ikas
    integer::ibfc
    l = facex
    ngp = pointx
    i = iconsidered
    if (dg .eq. 1) then
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
    else
      cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
    end if
    if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
      if (icoupleturb .eq. 1) then
        cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
      else
        cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
      end if
    end if
    if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
          if (dg .eq. 1) then
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          else !fv
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
            end if
          end if
          ikas = 1
        else
          !not periodic ones in my cpu
          call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
          n_node = 2
          cords(1:2) = zero
          cords(1:2) = cordinates2(n, nodes_list, n_node)
          poy(1) = cords(2)
          pox(1) = cords(1)
          leftv(1:nof_variables) = cleft(1:nof_variables)
          b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
          call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
          cright(1:nof_variables) = rightv(1:nof_variables)
          ikas = 2
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        else !fv
          cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &(1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
          end if
        end if
        ikas = 3
      end if
    else        !in other cpus they can only be periodic or mpi neighbours
      if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
        if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
          if (dg .eq. 1) then
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          else
            cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 1) then
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            else
              cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
            end if
          end if
        end if
      else
        if (dg .eq. 1) then
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        else
          cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          if (icoupleturb .eq. 1) then
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          else
            cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &(ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
          end if
        end if
        ikas = 4
      end if
    end if
  end subroutine get_states_bounds2d

  subroutine allocate_dg
    ! @brief
    ! allocates the gaussian quadrature volume points
    implicit none
    integer::i, k, i_qp, n_qp, i_face, nnd, iqp, idummy, loopc
    real, dimension(1:idegfree + 1)::tempint
    real::tempf, tempx
    allocate(qp_array(xmpielrank(n))); !allocates for 2d
    do i = 1, xmpielrank(n)
      select case (ielem(n, i)%ishape)
      case (1) !hexa
        allocate(qp_array(i)%x(qp_tetra*6))
        allocate(qp_array(i)%y(qp_tetra*6))
        allocate(qp_array(i)%z(qp_tetra*6))
        allocate(qp_array(i)%qp_weight(qp_tetra*6))
      case (2) !tetra
        allocate(qp_array(i)%x(qp_tetra))
        allocate(qp_array(i)%y(qp_tetra))
        allocate(qp_array(i)%z(qp_tetra))
        allocate(qp_array(i)%qp_weight(qp_tetra))
      case (3) !pyramid
        allocate(qp_array(i)%x(qp_tetra*2))
        allocate(qp_array(i)%y(qp_tetra*2))
        allocate(qp_array(i)%z(qp_tetra*2))
        allocate(qp_array(i)%qp_weight(qp_tetra*2))
      case (4) !prism
        allocate(qp_array(i)%x(qp_tetra*3))
        allocate(qp_array(i)%y(qp_tetra*3))
        allocate(qp_array(i)%z(qp_tetra*3))
        allocate(qp_array(i)%qp_weight(qp_tetra*3))
      case (5) !quad
        allocate(qp_array(i)%x(qp_triangle*2))
        allocate(qp_array(i)%y(qp_triangle*2))
        allocate(qp_array(i)%qp_weight(qp_triangle*2))
      case (6) !tetra
        allocate(qp_array(i)%x(qp_triangle))
        allocate(qp_array(i)%y(qp_triangle))
        allocate(qp_array(i)%qp_weight(qp_triangle))
      end select
      if (dg .eq. 1) then
        allocate(ilocal_recon3(i)%uleft_dg(nof_variables, ielem(n, i)%ifca, numberofpoints2))
        if ((itestcase .eq. 4) .or. ((governingequations .eq. -1) .and. (br2_yn .eq. 1))) then !ns
          allocate(ilocal_recon3(i)%br2_aux_var(nof_variables, dimensiona, ielem(n, i)%ifca, numberofpoints2))
          allocate(ilocal_recon3(i)%br2_local_lift(nof_variables, dimensiona, ielem(n, i)%ifca))
          ilocal_recon3(i)%br2_aux_var = zero
          ilocal_recon3(i)%br2_local_lift = zero
        end if
      end if
    end do
  end subroutine allocate_dg

  subroutine prestore_dg1(iconsidered)
    ! @brief
    ! prestores ielem(n,i)%delta_xyz, qp_array, surf_qpoints, mass matrix
    implicit none
    integer, intent(in)::iconsidered
    integer::i, k, i_qp, n_qp, i_face, nnd, iqp, idummy, loopc, icompwrt, eltype, elem_dec, ixx, number
    integer::count_1
    real, dimension(1:idegfree + 1)::tempint
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona)::cords
    real, dimension(1:dimensiona, 1:numberofpoints)::qpoints
    real, dimension(1:numberofpoints)::wequa3d
    real, dimension(1:6, 1:4, 1:dimensiona)::elem_listd
    real::voltemp, x1, y1, z1
    real::tempf, tempx
    number = ielem(n, iconsidered)%iorder
    if (dg .eq. 1) then
      icompwrt = -2
    end if
    i = iconsidered
    !store volume quadrature points
    eltype = ielem(n, i)%ishape
    elem_dec = ielem(n, i)%vdec
    do k = 1, ielem(n, i)%nonodes
      nodes_list(k, 1:dims) = inoder(ielem(n, i)%nodes(k))%cord(1:dims)
      vext(k, 1:dims) = nodes_list(k, 1:dims)
    end do
    if (dimensiona .eq. 2) then
      call decompose2(n, eltype, nodes_list, elem_listd)
    else
      call decompose3(n, eltype, nodes_list, elem_listd)
    end if
    select case (ielem(n, i)%ishape)
    case (1, 2, 3, 4) !hexa
      count_1 = 0
      do k = 1, elem_dec
        vext(1:4, 1:3) = elem_listd(k, 1:4, 1:3)
        call quadraturetetra(n, igqrules, vext, qpoints, wequa3d)
        voltemp = tetravolume(n, vext)
        do i_qp = 1, qp_tetra
          count_1 = count_1 + 1
          qp_array(i)%x(count_1) = (qpoints(1, i_qp) - ielem(n, i)%xxc)
          qp_array(i)%y(count_1) = (qpoints(2, i_qp) - ielem(n, i)%yyc)
          qp_array(i)%z(count_1) = (qpoints(3, i_qp) - ielem(n, i)%zzc)
          qp_array(i)%qp_weight(count_1) = wequa3d(i_qp)*voltemp
        end do
      end do
    case (5)
      count_1 = 0
      do k = 1, elem_dec
        vext(1:3, 1:2) = elem_listd(k, 1:3, 1:2)
        call quadraturetriangle(n, igqrules, vext, qpoints, wequa3d)
        voltemp = trianglevolume(n, vext)
        do i_qp = 1, qp_triangle
          count_1 = count_1 + 1
          qp_array(i)%x(count_1) = (qpoints(1, i_qp) - ielem(n, i)%xxc)
          qp_array(i)%y(count_1) = (qpoints(2, i_qp) - ielem(n, i)%yyc)
          qp_array(i)%qp_weight(count_1) = wequa3d(i_qp)*voltemp
        end do
      end do
    case (6)
      call quadraturetriangle(n, igqrules, vext, qpoints, wequa3d)
      n_qp = qp_triangle
      voltemp = trianglevolume(n, vext)
      do i_qp = 1, n_qp
        qp_array(i)%x(i_qp) = (qpoints(1, i_qp) - ielem(n, i)%xxc)
        qp_array(i)%y(i_qp) = (qpoints(2, i_qp) - ielem(n, i)%yyc)
        qp_array(i)%qp_weight(i_qp) = wequa3d(i_qp)*voltemp
      end do
    end select
    if (dg .eq. 1) then
      !store volume quadrature points
      integ_basis_dg(iconsidered)%value(1:idegfree) = zero
    end if
    tempint = zero
    n_qp = ielem(n, iconsidered)%itotalpoints
    do i_qp = 1, n_qp
      x1 = qp_array(i)%x(i_qp); y1 = qp_array(i)%y(i_qp)
      if (dimensiona .eq. 3) then
        z1 = qp_array(i)%z(i_qp)
      end if
      ixx = iconsidered
      if (dg .eq. 1) then
        if (dimensiona .eq. 2) then
          number = ielem(n, iconsidered)%iorder
          tempint(1:idegfree) = tempint(1:idegfree) + (basis_rec2d(n, x1, y1, number, ixx, idegfree, icompwrt)*qp_array(i)%qp_weight(i_qp))
        else
          number = ielem(n, iconsidered)%iorder
          tempint(1:idegfree) = tempint(1:idegfree) + (basis_rec(n, x1, y1, z1, number, ixx, idegfree, icompwrt)*qp_array(i)%qp_weight(i_qp))
        end if
      end if
    end do
    if (dg .eq. 1) then
      integ_basis_dg(iconsidered)%value(1:idegfree) = tempint(1:idegfree)
    end if
    icompwrt = 0
  end subroutine

  subroutine build_mass_matrix(n)
    ! @brief
    ! assembles the mass matrix
    ! requires: globals: ielem, qp_quad, qp_triangle, mass_matrix
    implicit none
    integer, intent(in)::n
    integer::i_elem, i_qp, n_qp, i_dof, j_dof, kmaxe, icompwrt, number_of_dog, iconsidered, ixx, number
    real::integ_test, integ_sm1, integ_sm2, phx1, phx2, kron, maxs, mins, higher, phx, integ_mm, x1, y1, z1
    real, allocatable, dimension(:, :)::totalmm, invmm
    real, allocatable, dimension(:)::basis_vector
    kmaxe = xmpielrank(n)
    allocate(totalmm(1:num_dg_dofs, 1:num_dg_dofs), invmm(1:num_dg_dofs, 1:num_dg_dofs), basis_vector(1:idegfree))

    do i_elem = 1, kmaxe
      totalmm(:, :) = zero
      invmm(:, :) = zero
      icompwrt = -2
      n_qp = ielem(n, i_elem)%itotalpoints
      iconsidered = i_elem
      number_of_dog = ielem(n, i_elem)%idegfree
      do i_dof = 1, num_dg_dofs
        do j_dof = 1, num_dg_dofs
          integ_mm = zero
          do i_qp = 1, n_qp
            ixx = i_elem;
            x1 = qp_array(i_elem)%x(i_qp);
            y1 = qp_array(i_elem)%y(i_qp)
            if (dimensiona .eq. 3) then
              z1 = qp_array(i_elem)%z(i_qp);
            end if
            kron = 1.0d0
            if (i_dof .ne. j_dof) then
              kron = 1.0d0
            end if
            if (dimensiona .eq. 2) then
              number = ielem(n, iconsidered)%iorder
              basis_vector(1:idegfree) = basis_rec2d(n, x1, y1, number, ixx, number_of_dog, icompwrt)
            else
              number = ielem(n, iconsidered)%iorder
              basis_vector(1:idegfree) = basis_rec(n, x1, y1, z1, number, ixx, number_of_dog, icompwrt)
            end if
            if (i_dof .eq. 1 .and. j_dof .eq. 1) then
              phx = 1.0d0
            else if (i_dof .eq. 1 .and. j_dof .ne. 1) then
              phx = basis_vector(j_dof - 1)
            else if (i_dof .ne. 1 .and. j_dof .eq. 1) then
              phx = basis_vector(i_dof - 1)
            else if (i_dof .ne. 1 .and. j_dof .ne. 1) then
              phx = basis_vector(i_dof - 1)*basis_vector(j_dof - 1)
            end if
            integ_mm = integ_mm + phx*qp_array(i_elem)%qp_weight(i_qp)*kron
          end do
          totalmm(i_dof, j_dof) = totalmm(i_dof, j_dof) + integ_mm
        end do
      end do
      call compmassinv(totalmm, invmm)
      m_1(i_elem)%val(:, :) = invmm(:, :)
      icompwrt = 0
    end do
    deallocate(totalmm, invmm, basis_vector)
  end subroutine build_mass_matrix

  subroutine compmassinv(totalmm, invmm)
    ! calculate the inverse of the input matrix with gauss-jordan elimination
    implicit none
    integer :: i, j, k, l, m, irow, num_dofs
    real:: big, dum
    real, allocatable, dimension(:, :)::a, b
    real, allocatable, dimension(:, :), intent(in)::totalmm
    real, allocatable, dimension(:, :), intent(inout)::invmm
    allocate(a(1:num_dg_dofs, 1:num_dg_dofs), b(1:num_dg_dofs, 1:num_dg_dofs))
    num_dofs = num_dg_dofs
    a(:, :) = totalmm(:, :)
    b(:, :) = zero
    do i = 1, num_dofs
      do j = 1, num_dofs
        b(i, j) = 0.0d0
      end do
      b(i, i) = 1.0d0
    end do
    do i = 1, num_dofs
      big = a(i, i)
      do j = i, num_dofs
        if (a(j, i) .gt. big) then
          big = a(j, i)
          irow = j
        end if
      end do
      ! interchange lines i with irow for both a() and b() matrices
      if (big .gt. a(i, i)) then
        do k = 1, num_dofs
          dum = a(i, k)                      ! matrix a()
          a(i, k) = a(irow, k)
          a(irow, k) = dum
          dum = b(i, k)                 ! matrix b()
          b(i, k) = b(irow, k)
          b(irow, k) = dum
        end do
      end if
      ! divide all entries in line i from a(i,j) by the value a(i,i);
      ! same operation for the identity matrix
      dum = a(i, i)
      do j = 1, num_dofs
        a(i, j) = a(i, j)/dum
        b(i, j) = b(i, j)/dum
      end do
      ! make zero all entries in the column a(j,i); same operation for indent()
      do j = i + 1, num_dofs
        dum = a(j, i)
        do k = 1, num_dofs
          a(j, k) = a(j, k) - dum*a(i, k)
          b(j, k) = b(j, k) - dum*b(i, k)
        end do
      end do
    end do
    do i = 1, num_dofs - 1
      do j = i + 1, num_dofs
        dum = a(i, j)
        do l = 1, num_dofs
          a(i, l) = a(i, l) - dum*a(j, l)
          b(i, l) = b(i, l) - dum*b(j, l)
        end do
      end do
    end do
    invmm(:, :) = b(:, :)
    deallocate(a, b)
  end subroutine compmassinv
end module dg_functions
