module recon
  use declaration
  use derivatives
  use library
  use transform
  use local
  use lapck
  use gradients
  use basis
  implicit none
contains
  subroutine average_stresses(n)
    implicit none
    integer, intent(in)::n
    integer::ii, i, iconsidered
    do ii = 1, nof_interior; i = el_int(ii); iconsidered = i
      call allgrads_inner_av(n, i)
    end do

    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      call allgrads_mix_av(n, i)
    end do
  end subroutine average_stresses

  subroutine memory_fast(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, kmaxe, idummy, l, nnd, iqp, ngp, iex
    integer::iconsidered, facex, pointx
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:dimensiona)::pox, poy, poz
    kmaxe = xmpielrank(n)
    if (dimensiona .eq. 3) then
      do i = 1, kmaxe
        if (ielem(n, i)%ishape .eq. 2) then
          allocate(ilocal_recon3(i)%qpoints(ielem(n, i)%ifca, qp_triangle, 3))
          if (srfg .eq. 1) then
            allocate(ilocal_recon3(i)%rpoints(ielem(n, i)%ifca, qp_triangle, 3))
            allocate(ilocal_recon3(i)%rotvel(ielem(n, i)%ifca, qp_triangle, 3))
          end if
          if (mrf .eq. 1) then
            allocate(ilocal_recon3(i)%rpoints(ielem(n, i)%ifca, qp_triangle, 3))
            allocate(ilocal_recon3(i)%rotvel(ielem(n, i)%ifca, qp_triangle, 3))
            allocate(ilocal_recon3(i)%mrf_origin(1:3))
            allocate(ilocal_recon3(i)%mrf_velocity(1:3))
          end if
        else
          allocate(ilocal_recon3(i)%qpoints(ielem(n, i)%ifca, qp_quad, 3))
          if (srfg .eq. 1) then
            allocate(ilocal_recon3(i)%rpoints(ielem(n, i)%ifca, qp_quad, 3))
            allocate(ilocal_recon3(i)%rotvel(ielem(n, i)%ifca, qp_quad, 3))
          end if
          if (mrf .eq. 1) then
            allocate(ilocal_recon3(i)%rpoints(ielem(n, i)%ifca, qp_quad, 3))
            allocate(ilocal_recon3(i)%rotvel(ielem(n, i)%ifca, qp_quad, 3))
            allocate(ilocal_recon3(i)%mrf_origin(1:3))
            allocate(ilocal_recon3(i)%mrf_velocity(1:3))
          end if
        end if
        iconsidered = i
        do l = 1, ielem(n, i)%ifca
          idummy = 0
          if ((iperiodicity .eq. 1) .and. (ielem(n, i)%interior .eq. 1)) then
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
    if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other cpu
      idummy = 1
    end if
    end if
    if (ielem(n, i)%types_faces(l) .eq. 5) then
      iqp = qp_quad
      nnd = 4
      if (idummy .eq. 0) then
      do k = 1, nnd
      vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
!     vext(k,1:3)=matmul(ilocal_recon3(i)%invccjac(:,:),vext(k,1:3)-ilocal_recon3(i)%vext_ref(1:3))
      end do
      else
      facex = l;
      call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
      end if
      call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
      else
      iqp = qp_triangle
      nnd = 3
      if (idummy .eq. 0) then
      do k = 1, nnd
      vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
      end do
      else
      facex = l;
      call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
    end if
    call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
    end if
    else
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
        nnd = 4
        do k = 1, nnd
          vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
        end do
        call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
       else
          iqp = qp_triangle
          nnd = 3
          do k = 1, nnd
            vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
          end do
          call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
         end if
      end if

    do ngp = 1, iqp                        !for gqp

    if (srfg .eq. 1) then
      ilocal_recon3(i)%rpoints(l, ngp, 1:3) = qpoints2d(1:3, ngp)
        pox(1:3) = ilocal_recon3(i)%rpoints(l, ngp, 1:3) - srf_origin(1:3)
        poy(1:3) = srf_velocity(1:3)
        ilocal_recon3(i)%rotvel(l, ngp, 1:3) = vect_function(pox, poy)
    end if

    if (mrf .eq. 1) then
      ilocal_recon3(i)%rpoints(l, ngp, 1:3) = qpoints2d(1:3, ngp)
      pox(1) = ielem(n, i)%xxc; pox(2) = ielem(n, i)%yyc; pox(3) = ielem(n, i)%zzc
      poy(1:3) = ilocal_recon3(i)%rpoints(l, ngp, 1:3)
      iconsidered = i
      facex = l
      pointx = ngp
      call mrfswitch(n, iconsidered, facex, pointx, pox, poy)
    end if

    end do        !ngp
    end do
    do l = 1, ielem(n, i)%ifca
      idummy = 0
      if ((iperiodicity .eq. 1) .and. (ielem(n, i)%interior .eq. 1)) then
        if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
          if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other cpu
            idummy = 1
          end if
        end if
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
        nnd = 4
        if (idummy .eq. 0) then
          do k = 1, nnd
            vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
            vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
          end do        !ngp
        else
          facex = l;
          call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
          do k = 1, nnd
            vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
          end do
        end if
        call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
        else
          iqp = qp_triangle
          nnd = 3
          if (idummy .eq. 0) then
            do k = 1, nnd
              vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
              vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
            end do

          else ! 2 dimensions
            facex = l;
            call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
            do k = 1, nnd
              vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
            end do
          end if
          call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
          end if
          else
            if (ielem(n, i)%types_faces(l) .eq. 5) then
              iqp = qp_quad
              nnd = 4
              do k = 1, nnd
                vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
                vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
              end do
              call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
            else
              iqp = qp_triangle
              nnd = 3
              do k = 1, nnd
                vext(k, 1:3) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
                vext(k, 1:3) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:3) - ilocal_recon3(i)%vext_ref(1:3))
              end do
              call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
            end if
          end if
          do ngp = 1, iqp                        !for gqp
            ilocal_recon3(i)%qpoints(l, ngp, 1:3) = qpoints2d(1:3, ngp)
          end do        !ngp
        end do
      end do
    else
      do i = 1, kmaxe
        allocate(ilocal_recon3(i)%qpoints(ielem(n, i)%ifca, qp_line, 2))
        iconsidered = i
        do l = 1, ielem(n, i)%ifca
          idummy = 0
          if ((iperiodicity .eq. 1) .and. (ielem(n, i)%interior .eq. 1)) then
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
    if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other cpu
                idummy = 1
              end if
            end if
            iqp = qp_line
            nnd = 2
            if (idummy .eq. 0) then
              do k = 1, nnd
                vext(k, 1:2) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
                !if (dg /= 1) then ! only transforming to reference space if not dg
                vext(k, 1:2) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:2) - ilocal_recon3(i)%vext_ref(1:2))
                !end if
              end do
            else
              facex = l;
              call coordinates_face_period2d1(n, iconsidered, facex, vext, nodes_list)
              do k = 1, nnd
                !if (dg /= 1) then ! only transforming to reference space if not dg
                vext(k, 1:2) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:2) - ilocal_recon3(i)%vext_ref(1:2))
                !end if
              end do
            end if
            call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
          else
            iqp = qp_line
            nnd = 2
            do k = 1, nnd
              vext(k, 1:2) = inoder(ielem(n, i)%nodes_faces(l, k))%cord(1:dims)
              !if (dg /= 1) then ! only transforming to reference space if not dg
              vext(k, 1:2) = matmul(ilocal_recon3(i)%invccjac(:, :), vext(k, 1:2) - ilocal_recon3(i)%vext_ref(1:2))
              !end if
            end do
            call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
          end if

          do ngp = 1, iqp !for gqp
            ilocal_recon3(i)%qpoints(l, ngp, 1:2) = qpoints2d(1:2, ngp) ! storing surface quadrature points
          end do !ngp
        end do
      end do
    end if
  end subroutine memory_fast

  subroutine extrapolate_bound_linear(usol, varcons, facex, pointx, iconsidered)
    implicit none
    integer, intent(in)::varcons, facex, pointx, iconsidered
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:, :, :), intent(in)::usol
    real::mp_pinfl, gammal
    if (wenwrt .eq. 3) then
      leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      leftv(1:nof_variables) = leftv(1:nof_variables) + usol(1:nof_variables, facex, pointx)
      call prim2cons(n, leftv)
    ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)+leftv(1:nof_variables)
    else
      ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)&
                                               + (u_c(iconsidered)%val(1, 1:nof_variables) + (usol(1:nof_variables, facex, pointx)))
    end if
    if (turbulenceequations .ge. 1) then
    ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)=ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)&
    +(u_ct(iconsidered)%val(1,1:turbulenceequations+passivescalar)+(usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,facex,pointx)))
    end if
  end subroutine extrapolate_bound_linear

  subroutine extrapolate_bound_muscl(usol, varcons, facex, pointx, iconsidered, slope)
    implicit none
    integer, intent(in)::varcons, facex, pointx, iconsidered
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:), intent(in)::slope
    real, allocatable, dimension(:, :, :), intent(in)::usol
    real::mp_pinfl, gammal
    if (wenwrt .eq. 3) then
      leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      leftv(1:nof_variables) = leftv(1:nof_variables) + usol(1:nof_variables, facex, pointx)*slope(1:nof_variables)
      call prim2cons(n, leftv)
    ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)+leftv(1:nof_variables)
    else
      ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)&
                        + (u_c(iconsidered)%val(1, 1:nof_variables) + (usol(1:nof_variables, facex, pointx)*slope(1:nof_variables)))
    end if
    if (turbulenceequations .ge. 1) then
    ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)=ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)&
    +(u_ct(iconsidered)%val(1,1:turbulenceequations+passivescalar)+(usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,facex,pointx)*slope(nof_variables+1:nof_variables+turbulenceequations+passivescalar)))
    end if
  end subroutine extrapolate_bound_muscl
  subroutine extrapolate_bound_musclx(usol, varcons, facex, pointx, iconsidered, slope)
    implicit none
    integer, intent(in)::varcons, facex, pointx, iconsidered
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:), intent(in)::slope
    real, allocatable, dimension(:, :, :), intent(in)::usol
    real::mp_pinfl, gammal
    if (wenwrt .eq. 3) then
      leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      leftv(1:nof_variables) = leftv(1:nof_variables) + usol(1:nof_variables, facex, pointx)*slope(1:nof_variables)
      call prim2cons(n, leftv)
    ilocal_recon3(iconsidered)%uleftx(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleftx(1:nof_variables,facex,pointx)+leftv(1:nof_variables)
    else
    ilocal_recon3(iconsidered)%uleftx(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleftx(1:nof_variables,facex,pointx)&
                        + (u_c(iconsidered)%val(1, 1:nof_variables) + (usol(1:nof_variables, facex, pointx)*slope(1:nof_variables)))
    end if
  end subroutine extrapolate_bound_musclx
  subroutine wenoweights(n)
    implicit none
    integer, intent(in)::n
    real::divisionbyzero
    integer::i, j, k, l, m, o, ll, iex, ieul, facx, ieleme, kkd, kmaxe, jf, ngp, iqp, nnd, ii, icd
    integer::idummy, power, itarget, iconsidered
    real::sumomegaatildel
    real::divbyzero, compf, checkf, tau_weno, tempxx
    real, dimension(1:numberofpoints2)::weights_q, weights_t

    kmaxe = xmpielrank(n)
    do ii = 1, nof_interior;
      i = el_int(ii)
      iconsidered = i
      ielem(n, i)%linc = lwci1
      ilocal_recon3(iconsidered)%uleft(:, :, :) = zero
      if (poly .eq. 4) then
        divbyzero = ielem(n, iconsidered)%totvolume**2
      else
        divbyzero = 10e-12
      end if
      power = 4
      if (adda .eq. 1) then
        call adda_filter(n, iconsidered)
      end if
      if (wenwrt .eq. 2) then
        call characteristic_reconstruction(iconsidered, idummy, divbyzero, power)
      else
        call cp_reconstruction(iconsidered, idummy, divbyzero, power)
      end if
      if (((turbulence .eq. 1) .or. (passivescalar .gt. 0)) .and. (icoupleturb .eq. 1)) then
        call cp_reconstruction_turb(iconsidered, idummy, divbyzero, power)
      end if
    end do

    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      ilocal_recon3(iconsidered)%uleft(:, :, :) = zero
      ielem(n, i)%linc = lwci1
      if (poly .eq. 4) then
        divbyzero = ielem(n, iconsidered)%totvolume**2
      else
        divbyzero = 10e-12
      end if
      power = 4
      if (adda .eq. 1) then
        call adda_filter(n, iconsidered)
      end if
      if (wenwrt .eq. 2) then
        call characteristic_reconstruction(iconsidered, idummy, divbyzero, power)
      else
        call cp_reconstruction(iconsidered, idummy, divbyzero, power)
      end if
      if (((turbulence .eq. 1) .or. (passivescalar .gt. 0)) .and. (icoupleturb .eq. 1)) then
        call cp_reconstruction_turb(iconsidered, idummy, divbyzero, power)
      end if
    end do
  end subroutine wenoweights
  subroutine characteristic_reconstruction(iconsidered, idummy, divbyzero, power)
    implicit none
    integer, intent(in)::iconsidered, power
    integer, intent(inout)::idummy
    real, intent(in)::divbyzero
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::veigl, veigr, rveigl, rveigr, leftv, rightv
    real, allocatable, dimension(:, :)::eigvl, eigvr
    integer::facex, kkd, l, i, itarget, iqp, ngp, icompwrt, icd, iex, ll, iadmisx, n_faces, iadmis, k
    real::lwcx1, tau_weno, ax, ay, az
    real::sumomegatilde(1:nof_variables)
    real, allocatable, dimension(:)::lamc
    real, allocatable, dimension(:, :)::limiteddw, consmatrix, consmatrixc, ressolution
    real, allocatable, dimension(:, :, :)::limiteddw_char, gradcharv
    real, allocatable, dimension(:, :, :, :)::lambda, smoothindicator, omegatilde, omega, wenoos, findw
    real, allocatable, dimension(:, :, :, :, :)::findw_char
    iadmis = ielem(n, iconsidered)%admis
    n_faces = ielem(n, iconsidered)%ifca
    allocate(lamc(1:iadmis))
    allocate(lambda(1:nof_variables, 1:iadmis, 1:n_faces, 1:2))
    allocate(smoothindicator(1:nof_variables, 1:iadmis, 1:n_faces, 1:2))
    allocate(omegatilde(1:nof_variables, 1:iadmis, 1:n_faces, 1:2))
    allocate(omega(1:nof_variables, 1:iadmis, 1:n_faces, 1:2))
    allocate(wenoos(1:nof_variables, 1:iadmis, 1:n_faces, 1:2))
    allocate(limiteddw(1:nof_variables, 0:idegfree))
    allocate(limiteddw_char(1:nof_variables, 0:idegfree, 1:iadmis))
    allocate(gradcharv(1:nof_variables, 1:iadmis, 0:idegfree))
    allocate(findw(1:nof_variables, 0:idegfree, 1:n_faces, 1:2))
    allocate(findw_char(1:nof_variables, 0:idegfree, 1:iadmis, 1:n_faces, 1:2))
    allocate(consmatrix(1:numberofpoints2*n_faces, 1:idegfree))
    allocate(consmatrixc(1:numberofpoints2*n_faces, 1:idegfree))
    allocate(ressolution(1:numberofpoints2*n_faces, 1:nof_variables))
    allocate(eigvl(1:nof_variables, 1:nof_variables), eigvr(1:nof_variables, 1:nof_variables))
    i = iconsidered
    lwcx1 = ielem(n, i)%linc
    do l = 1, ielem(n, i)%ifca  !loop faces
      !define
      angle1 = ielem(n, i)%faceanglex(l);
      angle2 = ielem(n, i)%faceangley(l)
      facex = l
      if (dimensiona .eq. 3) then
        nx = (cos(angle1)*sin(angle2));
        ny = (sin(angle1)*sin(angle2));
        nz = (cos(angle2))
      else
        nx = angle1
        ny = angle2
      end if
      veigl(1:nof_variables) = u_c(i)%val(1, 1:nof_variables);
      if (dimensiona .eq. 3) then
        call rotatef(n, rveigl, veigl, angle1, angle2)
      else
        call rotatef2d(n, rveigl, veigl, angle1, angle2)
      end if
      call weno_neighbour(iconsidered, facex, veigl, veigr, nx, ny, nz, angle1, angle2, idummy)
      if (dimensiona .eq. 3) then
        call rotatef(n, rveigr, veigr, angle1, angle2)
      else
        call rotatef2d(n, rveigr, veigr, angle1, angle2)
      end if
      if (dimensiona .eq. 3) then
        call compute_eigenvectors(n, rveigl, rveigr, eigvl, eigvr, gamma)
      else
        call compute_eigenvectors2d(n, rveigl, rveigr, eigvl, eigvr, gamma)
      end if
      lambda(:, :, l, 1) = zero;
      smoothindicator(:, :, l, 1) = zero;
      omegatilde(:, :, l, 1) = zero;
      omega(:, :, l, 1) = zero; facex = l
      call compute_gradcharv_smoothindicator(iconsidered, facex, eigvl, gradcharv, smoothindicator)
      lambda(1:nof_variables, :, l, 1) = 1.0d0; lambda(1:nof_variables, 1, l, 1) = lwcx1
      if (ees .eq. 5) then
      do kkd = 1, nof_variables
        lamc(1) = (1.0d0 - (1.0d0/lwcx1))
        lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
        lambda(kkd, 1:ielem(n, i)%admis, l, 1) = lamc(1:ielem(n, i)%admis)
      end do
      end if
      do kkd = 1, nof_variables
        sumomegatilde(kkd) = zero
        if (ees .eq. 5) then
          tau_weno = zero
          if (wenoz .eq. 1) then
            do ll = 1, ielem(n, i)%admis
              tau_weno = tau_weno + (abs(smoothindicator(kkd, 1, l, 1) - smoothindicator(kkd, ll, l, 1)))
            end do
            tau_weno = (tau_weno/(ielem(n, i)%admis - 1))!**power
            do ll = 1, ielem(n, i)%admis
        omegatilde(kkd, ll, l, 1) = (lambda(kkd, ll, l, 1))*(1.0d0 + (tau_weno/(divbyzero + smoothindicator(kkd, ll, l, 1)))**power)
            end do
          else
            do ll = 1, ielem(n, i)%admis
              omegatilde(kkd, ll, l, 1) = (lambda(kkd, ll, l, 1))/((divbyzero + smoothindicator(kkd, ll, l, 1))**power)
            end do
          end if
        else
          do ll = 1, ielem(n, i)%admis
            omegatilde(kkd, ll, l, 1) = (lambda(kkd, ll, l, 1))/((divbyzero + smoothindicator(kkd, ll, l, 1))**power)
          end do
        end if

        do ll = 1, ielem(n, i)%admis
          sumomegatilde(kkd) = sumomegatilde(kkd) + omegatilde(kkd, ll, l, 1)
        end do
        do ll = 1, ielem(n, i)%admis
          omega(kkd, ll, l, 1) = (omegatilde(kkd, ll, l, 1))/sumomegatilde(kkd)
        end do
        do ll = 1, ielem(n, i)%admis
          wenoos(kkd, ll, l, 1) = omega(kkd, ll, l, 1)
        end do

      end do       !finished the loop for all the variables
      limiteddw(:, :) = zero
      if (ees .eq. 5) then
        limiteddw_char(:, :, :) = zero
        do ll = 1, ielem(n, i)%admis; if (ll .eq. 1) then
            itarget = ielem(n, i)%idegfree
          else
            itarget = idegfree2
          end if
          do k = 0, itarget
                        limiteddw_char(1:nof_variables,k,1)=limiteddw_char(1:nof_variables,k,1)+gradcharv(1:nof_variables,ll,k)*wenoos(1:nof_variables,ll,l,1)
          end do; end do
        findw_char(:, :, l, 1, :) = zero
        do k = 0, ielem(n, i)%idegfree
    findw_char(1:nof_variables, k, l, 1, 1) = matmul(eigvr(1:nof_variables, 1:nof_variables), limiteddw_char(1:nof_variables, k, 1))
        end do
      else
        limiteddw(:, :) = zero
        do k = 0, ielem(n, i)%idegfree; do ll = 1, ielem(n, i)%admis
 limiteddw(1:nof_variables, k) = limiteddw(1:nof_variables, k) + gradcharv(1:nof_variables, ll, k)*wenoos(1:nof_variables, ll, l, 1)
          end do; end do
        findw(:, :, l, 1) = zero
        do k = 0, ielem(n, i)%idegfree
          findw(1:nof_variables, k, l, 1) = matmul(eigvr(1:nof_variables, 1:nof_variables), limiteddw(1:nof_variables, k))
        end do
      end if
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      icd = 0
      do ngp = 1, iqp
        ax = ilocal_recon3(i)%qpoints(l, ngp, 1);
        ay = ilocal_recon3(i)%qpoints(l, ngp, 2);
        if (dimensiona .eq. 3) then
          az = ilocal_recon3(i)%qpoints(l, ngp, 3)
        end if
        icd = icd + 1; icompwrt = 0
        if (dimensiona .eq. 3) then
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        else
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec2d(n, ax, ay, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        end if
        icompwrt = 0
        if (ees .eq. 5) then;
          icompwrt = 1
          if (dimensiona .eq. 3) then
            consmatrixc(icd, 1:idegfree2) = basis_rec(n, ax, ay, az, iorder2, i, idegfree2, icompwrt)
          else
            consmatrixc(icd, 1:idegfree2) = basis_rec2d(n, ax, ay, iorder2, i, idegfree2, icompwrt)
          end if
          icompwrt = 0; end if
      end do
      if (ees .eq. 5) then;
        ilocal_recon3(i)%uleft(1:nof_variables, l, :) = zero
        do ngp = 1, iqp
          ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) &
                                                            + findw_char(1:nof_variables, 0, l, 1, 1)
        end do
        ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),transpose(findw_char(1:nof_variables,1:ielem(n,i)%idegfree,l,1,1)))
        icd = 0;
        do ngp = 1, iqp; icd = icd + 1
          ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) &
                                                            + ressolution(icd, 1:nof_variables)
        end do
      else
        do ngp = 1, iqp
          ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) = findw(1:nof_variables, 0, l, 1)
        end do

      ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),transpose(findw(1:nof_variables,1:ielem(n,i)%idegfree,l,1)))
        icd = 0;
        do ngp = 1, iqp;
          icd = icd + 1
          ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) &
                                                            + ressolution(icd, 1:nof_variables)
        end do
      end if
    end do                        !faces
    deallocate(lamc, lambda, smoothindicator, omegatilde, &
                omega, wenoos, limiteddw, limiteddw_char, gradcharv, findw, &
                findw_char, ressolution, consmatrix, consmatrixc, eigvl, eigvr)

  end subroutine characteristic_reconstruction

  subroutine weno_neighbour(iconsidered, facex, veigl, veigr, nx, ny, nz, angle1, angle2, idummy)
    implicit none
    real, dimension(1:nof_variables), intent(inout)::veigr
    real, dimension(1:nof_variables), intent(inout)::veigl
    real, intent(in)::nx, ny, nz, angle1, angle2
    integer, intent(in)::iconsidered, facex
    integer, intent(inout)::idummy
    real::mp_pinfl, gammal
    integer::i, j, k, l, var2, b_code, n_node
    real, dimension(1:nof_variables)::leftv, srf_speed, srf_speedrot, rightv
    real, dimension(1:dimensiona)::pox, poy, poz, cords
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(turbulenceequations)::cturbl, cturbr
    real, dimension(1:nof_variables)::cright_rot, cleft_rot
    integer::ibfc
    if (dimensiona .eq. 3) then
      l = facex
      i = iconsidered
      if (ielem(n, i)%interior .eq. 0) then
        veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables);
      else
        if (ilocal_recon3(i)%mrf .eq. 1) then
          srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, 1, 1:3)
          call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
        end if
        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in my cpu
              veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
              idummy = 1
              if (per_rot .eq. 1) then
                veigr(2:4) = rotate_per_1(veigr(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
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
              leftv(1:nof_variables) = veigl(1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
                                  call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              veigr(1:nof_variables) = rightv(1:nof_variables)
            end if
          else
            veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
          end if
        else
          !other my cpu
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other cpu
              veigr(1:nof_variables) = (iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables))
              idummy = 1
              if (per_rot .eq. 1) then
                veigr(2:4) = rotate_per_1(veigr(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
              end if
            end if
          else
            veigr(1:nof_variables) = (iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables))
          end if
        end if
      end if
    else
      l = facex
      i = iconsidered
      if (ielem(n, i)%interior .eq. 0) then
        veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables);
      else
        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
              idummy = 1
            else
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              n_node = 2
              cords(1:2) = zero
              cords(1:2) = cordinates2(n, nodes_list, n_node)
              pox(1) = cords(1)
              poy(1) = cords(2)
              leftv(1:nof_variables) = veigl(1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
                                   call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              veigr(1:nof_variables) = rightv(1:nof_variables)
            end if
          else
            !fluid neighbour
            veigr(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
          end if
        else
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
              veigr(1:nof_variables) = (iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables))
              idummy = 1
            end if
          else
            veigr(1:nof_variables) = (iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables))
          end if
        end if
      end if
    end if
  end subroutine weno_neighbour
  subroutine cp_reconstruction(iconsidered, idummy, divbyzero, power)
    implicit none
    integer, intent(in)::iconsidered, power
    integer, intent(inout)::idummy
    real, intent(in)::divbyzero
    integer::facex, kkd, l, i, itarget, iqp, ngp, icompwrt, iex, ll, iadmis, n_faces
    real::lwcx1, ax, ay, az, tau_weno, sumomegaatildel
    integer::icd
    real, dimension(1:nof_variables)::leftv, rightv
    real, allocatable, dimension(:)::grad1al, indicatematrixal, grad3al
    real, allocatable, dimension(:)::lambdaal, omegaatildel, smoothindicatoral, lamc, omegaal
    real, allocatable, dimension(:, :)::consmatrix, consmatrixc, grad5alc, gradssl, weno, ressolution
    iadmis = ielem(n, iconsidered)%admis
    n_faces = ielem(n, iconsidered)%ifca
    allocate(grad1al(1:idegfree), indicatematrixal(1:idegfree))
    allocate(grad3al(idegfree), lambdaal(1:iadmis), omegaatildel(1:iadmis), smoothindicatoral(1:iadmis))
    allocate(lamc(1:iadmis), omegaal(1:iadmis))
    allocate(consmatrix(1:numberofpoints2*n_faces, 1:idegfree), consmatrixc(1:numberofpoints2*n_faces, 1:idegfree))
    allocate(grad5alc(1:idegfree, 1:nof_variables), gradssl(1:idegfree, 1:nof_variables))
    allocate(weno(1:nof_variables + turbulenceequations + passivescalar, 1:iadmis))
    allocate(ressolution(1:numberofpoints2*n_faces, 1:nof_variables))
    i = iconsidered
    lwcx1 = ielem(n, i)%linc
    do iex = 1, nof_variables
      lambdaal = zero; smoothindicatoral = zero; omegaatildel = zero; omegaal = zero
      if (ees .eq. 5) then
        lamc(:) = zero; grad3al(:) = zero; lamc(1) = (1.0d0 - (1.0d0/lwcx1)); lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
        lambdaal(1:ielem(n, i)%admis) = lamc(1:ielem(n, i)%admis)
        !sum the low degree polynomials first
        do ll = 2, ielem(n, i)%admis
          grad3al(1:idegfree2) = grad3al(1:idegfree2) + (lamc(ll)*ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, iex))
        end do
        !this is the zero polynomial
                            grad1al(1:ielem(n,i)%idegfree)=(1.0d0/lamc(1))*(ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,iex)-grad3al(1:ielem(n,i)%idegfree))
        grad5alc(1:ielem(n, i)%idegfree, iex) = grad1al(1:ielem(n, i)%idegfree)
        do ll = 1, ielem(n, i)%admis
          if (ll .eq. 1) then
            indicatematrixal(1:ielem(n,i)%idegfree)=matmul(ilocal_recon3(i)%indicator(1:ielem(n,i)%idegfree,1:ielem(n,i)%idegfree),grad1al(1:ielem(n,i)%idegfree))
            smoothindicatoral(ll) = dot_product(grad1al(1:ielem(n, i)%idegfree), indicatematrixal(1:ielem(n, i)%idegfree))
          else
            grad1al(1:idegfree2) = ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, iex)
            indicatematrixal(1:idegfree2) = matmul(ilocal_recon3(i)%indicatorc(1:idegfree2, 1:idegfree2), grad1al(1:idegfree2))
            smoothindicatoral(ll) = dot_product(grad1al(1:idegfree2), indicatematrixal(1:idegfree2))
          end if
        end do
      else
        do ll = 1, ielem(n, i)%admis
          grad1al(:) = zero
          indicatematrixal(:) = zero
          grad1al(1:ielem(n, i)%idegfree) = ilocal_recon5(iconsidered)%gradients(ll, 1:ielem(n, i)%idegfree, iex)
          indicatematrixal(1:ielem(n,i)%idegfree)=matmul(ilocal_recon3(i)%indicator(1:ielem(n,i)%idegfree,1:ielem(n,i)%idegfree),grad1al(1:ielem(n,i)%idegfree))
          smoothindicatoral(ll) = dot_product(grad1al(1:ielem(n, i)%idegfree), indicatematrixal(1:ielem(n, i)%idegfree))
        end do
      end if
      lambdaal(:) = 1.0d0
      lambdaal(1) = lwcx1
      if (ees .eq. 5) then
        lamc(1) = (1.0d0 - (1.0d0/lwcx1))
        lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
        lambdaal(1:ielem(n, i)%admis) = lamc(1:ielem(n, i)%admis)
      end if
      if (ees .eq. 5) then
        if (wenoz .eq. 1) then
          tau_weno = zero
          do ll = 1, ielem(n, i)%admis
            tau_weno = tau_weno + (abs(smoothindicatoral(1) - smoothindicatoral(ll)))
          end do
          tau_weno = (tau_weno/(ielem(n, i)%admis - 1))
          do ll = 1, ielem(n, i)%admis
            omegaatildel(ll) = (lambdaal(ll))*(1.0d0 + (tau_weno/(divbyzero + smoothindicatoral(ll)))**power)
          end do
        else
          do ll = 1, ielem(n, i)%admis
            omegaatildel(ll) = (lambdaal(ll))/((divbyzero + smoothindicatoral(ll))**power)
          end do

        end if
      else
        do ll = 1, ielem(n, i)%admis
          omegaatildel(ll) = (lambdaal(ll))/((divbyzero + smoothindicatoral(ll))**power)
        end do
      end if
      sumomegaatildel = zero
      do ll = 1, ielem(n, i)%admis
        sumomegaatildel = sumomegaatildel + omegaatildel(ll)
      end do
      do ll = 1, ielem(n, i)%admis
        omegaal(ll) = (omegaatildel(ll))/sumomegaatildel
      end do
      do ll = 1, ielem(n, i)%admis
        weno(iex, ll) = omegaal(ll)

        if (iex .eq. 1) then
          ielem(n, i)%wcx(1) = weno(iex, 1)
        end if
      end do
    end do
    icd = 0
    do l = 1, ielem(n, i)%ifca        !faces
      if (dimensiona .eq. 3) then
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if
      do ngp = 1, iqp                        !for gqp
        icd = icd + 1
        ax = ilocal_recon3(i)%qpoints(l, ngp, 1)
        ay = ilocal_recon3(i)%qpoints(l, ngp, 2)
        if (dimensiona .eq. 3) then
          az = ilocal_recon3(i)%qpoints(l, ngp, 3)
        end if
        icompwrt = 0
        if (dimensiona .eq. 3) then
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        else
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec2d(n, ax, ay, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        end if
        if (ees .eq. 5) then
          icompwrt = 1
          if (dimensiona .eq. 3) then
            consmatrixc(icd, 1:idegfree2) = basis_rec(n, ax, ay, az, iorder2, i, idegfree2, icompwrt)
          else
            consmatrixc(icd, 1:idegfree2) = basis_rec2d(n, ax, ay, iorder2, i, idegfree2, icompwrt)
          end if
          icompwrt = 0
        end if
      end do
    end do        !faces
    ilocal_recon3(i)%uleft(:, :, :) = zero
    if (dg .eq. 1) then
      ilocal_recon6(i)%dg2fv(1:ielem(n, i)%idegfree, :) = zero
    end if
    do ll = 1, ielem(n, i)%admis        !stencils
      if (ees .eq. 5) then
        if (ll .eq. 1) then
          gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = grad5alc(1:ielem(n, i)%idegfree, 1:nof_variables)
   ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
        else
          gradssl(1:idegfree2, 1:nof_variables) = ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, 1:nof_variables)
          ressolution(1:icd, 1:nof_variables) = matmul(consmatrixc(1:icd, 1:idegfree2), gradssl(1:idegfree2, 1:nof_variables))
        end if
      else
gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = ilocal_recon5(iconsidered)%gradients(ll, 1:ielem(n, i)%idegfree, 1:nof_variables)
   ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
      end if
      icd = 0
      do l = 1, ielem(n, i)%ifca
        if (dimensiona .eq. 3) then
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad;
          else
            iqp = qp_triangle;
          end if
        else
          iqp = qp_line;
        end if
        do ngp = 1, iqp
          icd = icd + 1
          call extrapolate_bound(ressolution, iex, l, ngp, i, icd, ll, weno)
        end do
      end do
      do iex = 1, nof_variables        !components
        if (dg .eq. 1) then
          if (ees .eq. 5) then
            if (ll .eq. 1) then
              ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)=ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)+(gradssl(1:idegfree,iex)*weno(iex,ll))
            else
     ilocal_recon6(i)%dg2fv(1:idegfree2, iex) = ilocal_recon6(i)%dg2fv(1:idegfree2, iex) + (gradssl(1:idegfree2, iex)*weno(iex, ll))
            end if
          else
            ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)=ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)+(gradssl(1:idegfree,iex)*weno(iex,ll))
          end if
        end if
      end do
    end do  !stencils finished
    if (wenwrt .eq. 3) then
      do l = 1, ielem(n, i)%ifca
        if (dimensiona .eq. 3) then
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad;
          else
            iqp = qp_triangle;
          end if
        else
          iqp = qp_line
        end if
        do ngp = 1, iqp
          leftv(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
          call prim2cons(n, leftv)
          ilocal_recon3(i)%uleft(1:nof_variables, l, ngp) = leftv(1:nof_variables)
        end do
      end do
    end if
    deallocate(grad1al, indicatematrixal, grad3al, lambdaal, omegaatildel, smoothindicatoral)
    deallocate(lamc, omegaal)
    deallocate(consmatrix, consmatrixc, grad5alc, gradssl)
    deallocate(weno)
    deallocate(ressolution)
  end subroutine cp_reconstruction
  subroutine cp_reconstruction_turb(iconsidered, idummy, divbyzero, power)
    implicit none
    integer, intent(in)::iconsidered, power
    integer, intent(inout)::idummy
    real, intent(in)::divbyzero
    integer::facex, kkd, l, i, itarget, iqp, ngp, icompwrt, icd, iex, ll, iadmis, n_faces
    real::lwcx1, tau_weno, ax, ay, az, sumomegaatildel
    real, dimension(1:nof_variables)::leftv, rightv
    real, allocatable, dimension(:)::grad1al, indicatematrixal, grad3al
    real, allocatable, dimension(:)::lambdaal, omegaatildel, smoothindicatoral, lamc, omegaal
    real, allocatable, dimension(:, :)::consmatrix, consmatrixc, grad5alc, gradssl, weno, ressolution
    iadmis = ielem(n, iconsidered)%admis
    n_faces = ielem(n, iconsidered)%ifca
    allocate(grad1al(1:idegfree), indicatematrixal(1:idegfree))
    allocate(grad3al(idegfree), lambdaal(1:iadmis), omegaatildel(1:iadmis), smoothindicatoral(1:iadmis))
    allocate(lamc(1:iadmis), omegaal(1:iadmis))
    allocate(consmatrix(1:numberofpoints2*n_faces, 1:idegfree), consmatrixc(1:numberofpoints2*n_faces, 1:idegfree))
  allocate(grad5alc(1:idegfree, 1:turbulenceequations + passivescalar), gradssl(1:idegfree, 1:turbulenceequations + passivescalar))
    allocate(weno(1:nof_variables + turbulenceequations + passivescalar, 1:iadmis))
    allocate(ressolution(1:numberofpoints2*n_faces, 1:turbulenceequations + passivescalar))
    i = iconsidered
    lwcx1 = ielem(n, i)%linc
    do iex = 1, turbulenceequations + passivescalar
      lambdaal = zero; smoothindicatoral = zero; omegaatildel = zero; omegaal = zero
      if (ees .eq. 5) then
        lamc(:) = zero; grad3al(:) = zero; lamc(1) = (1.0d0 - (1.0d0/lwcx1)); lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
        lambdaal(1:ielem(n, i)%admis) = lamc(1:ielem(n, i)%admis)
        !sum the low degree polynomials first
        do ll = 2, ielem(n, i)%admis
          grad3al(1:idegfree2) = grad3al(1:idegfree2) + (lamc(ll)*ilocal_recon5(iconsidered)%gradientsc2(ll, 1:idegfree2, iex))
        end do
        !this is the zero polynomial
                            grad1al(1:ielem(n,i)%idegfree)=(1.0d0/lamc(1))*(ilocal_recon5(iconsidered)%gradients2(1,1:ielem(n,i)%idegfree,iex)-grad3al(1:ielem(n,i)%idegfree))
        grad5alc(1:ielem(n, i)%idegfree, iex) = grad1al(1:ielem(n, i)%idegfree)
        do ll = 1, ielem(n, i)%admis
          if (ll .eq. 1) then
            indicatematrixal(1:ielem(n,i)%idegfree)=matmul(ilocal_recon3(i)%indicator(1:ielem(n,i)%idegfree,1:ielem(n,i)%idegfree),grad1al(1:ielem(n,i)%idegfree))
            smoothindicatoral(ll) = dot_product(grad1al(1:ielem(n, i)%idegfree), indicatematrixal(1:ielem(n, i)%idegfree))
          else
            grad1al(1:idegfree2) = ilocal_recon5(iconsidered)%gradientsc2(ll, 1:idegfree2, iex)
            smoothindicatoral(ll) = dot_product(grad1al(1:idegfree2), indicatematrixal(1:idegfree2))
          end if
        end do
      else
        do ll = 1, ielem(n, i)%admis
          grad1al(:) = zero
          indicatematrixal(:) = zero
          grad1al(1:ielem(n, i)%idegfree) = ilocal_recon5(iconsidered)%gradients2(ll, 1:ielem(n, i)%idegfree, iex)
          indicatematrixal(1:ielem(n,i)%idegfree)=matmul(ilocal_recon3(i)%indicator(1:ielem(n,i)%idegfree,1:ielem(n,i)%idegfree),grad1al(1:ielem(n,i)%idegfree))
          smoothindicatoral(ll) = dot_product(grad1al(1:ielem(n, i)%idegfree), indicatematrixal(1:ielem(n, i)%idegfree))
        end do
      end if
      lambdaal(:) = 1.0d0
      lambdaal(1) = lwcx1

      if (ees .eq. 5) then
        lamc(1) = (1.0d0 - (1.0d0/lwcx1))
        lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
        lambdaal(1:ielem(n, i)%admis) = lamc(1:ielem(n, i)%admis)
      end if
      if (ees .eq. 5) then
        if (wenoz .eq. 1) then
          tau_weno = zero
          do ll = 1, ielem(n, i)%admis
            tau_weno = tau_weno + (abs(smoothindicatoral(1) - smoothindicatoral(ll)))
          end do
          tau_weno = (tau_weno/(ielem(n, i)%admis - 1))
          do ll = 1, ielem(n, i)%admis
            omegaatildel(ll) = (lambdaal(ll))*(1.0d0 + (tau_weno/(divbyzero + smoothindicatoral(ll)))**power)
          end do
        else
          do ll = 1, ielem(n, i)%admis
            omegaatildel(ll) = (lambdaal(ll))/((divbyzero + smoothindicatoral(ll))**power)
          end do
        end if
      else
        do ll = 1, ielem(n, i)%admis
          omegaatildel(ll) = (lambdaal(ll))/((divbyzero + smoothindicatoral(ll))**power)
        end do
      end if
      sumomegaatildel = zero
      do ll = 1, ielem(n, i)%admis
        sumomegaatildel = sumomegaatildel + omegaatildel(ll)
      end do
      do ll = 1, ielem(n, i)%admis
        omegaal(ll) = (omegaatildel(ll))/sumomegaatildel
      end do
      do ll = 1, ielem(n, i)%admis
        weno(iex, ll) = omegaal(ll)
      end do
    end do
    icd = 0
    do l = 1, ielem(n, i)%ifca        !faces
      if (dimensiona .eq. 3) then
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if
      do ngp = 1, iqp                        !for gqp
        icd = icd + 1
        ax = ilocal_recon3(i)%qpoints(l, ngp, 1)
        ay = ilocal_recon3(i)%qpoints(l, ngp, 2)
        if (dimensiona .eq. 3) then
          az = ilocal_recon3(i)%qpoints(l, ngp, 3)
        end if
        icompwrt = 0
        if (dimensiona .eq. 3) then
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        else
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec2d(n, ax, ay, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        end if
        if (ees .eq. 5) then
          icompwrt = 1
          if (dimensiona .eq. 3) then
            consmatrixc(icd, 1:idegfree2) = basis_rec(n, ax, ay, az, iorder2, i, idegfree2, icompwrt)
          else
            consmatrixc(icd, 1:idegfree2) = basis_rec2d(n, ax, ay, iorder2, i, idegfree2, icompwrt)
          end if
          icompwrt = 0
        end if
      end do
    end do        !faces
    ilocal_recon3(i)%uleftturb(:, :, :) = zero
    if (dg .eq. 1) then
      ilocal_recon6(i)%dg2fv(1:ielem(n, i)%idegfree, :) = zero
    end if
    do ll = 1, ielem(n, i)%admis        !stencils

      if (ees .eq. 5) then
      if (ll .eq. 1) then
        gradssl(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=grad5alc(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)

        ressolution(1:icd,1:turbulenceequations+passivescalar)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar))
      else
        gradssl(1:idegfree2,1:turbulenceequations+passivescalar)=ilocal_recon5(iconsidered)%gradientsc2(ll,1:idegfree2,1:turbulenceequations+passivescalar)
        ressolution(1:icd,1:turbulenceequations+passivescalar)=matmul(consmatrixc(1:icd,1:idegfree2),gradssl(1:idegfree2,1:turbulenceequations+passivescalar))
      end if
      else
        gradssl(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=ilocal_recon5(iconsidered)%gradients2(ll,1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)
        ressolution(1:icd,1:turbulenceequations+passivescalar)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar))
      end if
      icd = 0
      do l = 1, ielem(n, i)%ifca
        if (dimensiona .eq. 3) then
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad;
          else
            iqp = qp_triangle;
          end if
        else
          iqp = qp_line;
        end if

        do ngp = 1, iqp
          icd = icd + 1
          call extrapolate_boundt(ressolution, iex, l, ngp, i, icd, ll, weno)
        end do
      end do
      do iex = 1, turbulenceequations + passivescalar        !components
      if (dg .eq. 1) then
        if (ees .eq. 5) then
          if (ll .eq. 1) then
            ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)=ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)+(gradssl(1:idegfree,iex)*weno(iex,ll))
          else
        ilocal_recon6(i)%dg2fv(1:idegfree2, iex) = ilocal_recon6(i)%dg2fv(1:idegfree2, iex) + (gradssl(1:idegfree2, iex)*weno(iex, ll))
          end if
        else
          ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)=ilocal_recon6(i)%dg2fv(1:ielem(n,i)%idegfree,iex)+(gradssl(1:idegfree,iex)*weno(iex,ll))
        end if
      end if
      end do
    end do
    deallocate(grad1al, indicatematrixal, grad3al, lambdaal, omegaatildel, smoothindicatoral)
    deallocate(lamc, omegaal)
    deallocate(consmatrix, consmatrixc, grad5alc, gradssl)
    deallocate(weno)
    deallocate(ressolution)
  end subroutine cp_reconstruction_turb
  subroutine extrapolate_bound(ressolution, varcons, facex, pointx, iconsidered, insten, llx, weno)
    implicit none
    integer, intent(in)::varcons, facex, pointx, iconsidered, insten, llx
    real, allocatable, dimension(:, :), intent(in)::weno
    real, allocatable, dimension(:, :), intent(in)::ressolution
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    if (wenwrt .eq. 3) then        !primitive
      leftv(1:nof_variables) = u_c(iconsidered)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)&
                                      + ((leftv(1:nof_variables) + ressolution(insten, 1:nof_variables))*weno(1:nof_variables, llx))
    else
      ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)=ilocal_recon3(iconsidered)%uleft(1:nof_variables,facex,pointx)&
                      + (u_c(iconsidered)%val(1, 1:nof_variables) + ressolution(insten, 1:nof_variables))*weno(1:nof_variables, llx)
    end if
  end subroutine extrapolate_bound

  subroutine extrapolate_boundt(ressolution, varcons, facex, pointx, iconsidered, insten, llx, weno)
    implicit none
    integer, intent(in)::varcons, facex, pointx, iconsidered, insten, llx
    real, allocatable, dimension(:, :), intent(in)::weno
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:, :), intent(in)::ressolution
    real::mp_pinfl, gammal
    ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)=ilocal_recon3(iconsidered)%uleftturb(1:turbulenceequations+passivescalar,facex,pointx)&
    +(u_ct(iconsidered)%val(1,1:turbulenceequations+passivescalar)+ressolution(insten,1:turbulenceequations+passivescalar))*weno(1:turbulenceequations+passivescalar,llx)
  end subroutine extrapolate_boundt
  subroutine diag_at_b_a(iconsidered, a_char, b_char, x_char)
    implicit none
    integer, intent(in):: iconsidered
    real, dimension(:, :, :), allocatable:: ba_char
    real, allocatable, dimension(:, :), intent(inout)::b_char, x_char
    real, allocatable, dimension(:, :, :), intent(inout)::a_char
    integer:: nn, mm
    integer:: i, ll, ics
    nn = size(a_char, 1) ! = size(b,1) = size(b,2)
    mm = size(a_char, 2)
    allocate(ba_char(nn, mm, ielem(n, iconsidered)%admis))
    x_char = zero
    if (ees .eq. 5) then
      do ll = 1, 1

        ba_char(1:ielem(n,iconsidered)%idegfree,1:nof_variables,ll)=matmul(b_char(1:ielem(n,iconsidered)%idegfree,1:ielem(n,iconsidered)%idegfree),a_char(1:ielem(n,iconsidered)%idegfree,1:nof_variables,ll))

      end do
      b_char(1:idegfree2, 1:idegfree2) = ilocal_recon3(iconsidered)%indicatorc(1:idegfree2, 1:idegfree2)
      do ll = 2, ielem(n, iconsidered)%admis
      ba_char(1:idegfree2, 1:nof_variables, ll) = matmul(b_char(1:idegfree2, 1:idegfree2), a_char(1:idegfree2, 1:nof_variables, ll))
      end do
      do ll = 1, 1; do i = 1, mm
          !x_char(i,ll) = dot(a_char(:,i,ll), ba_char(:,i,ll))
      x_char(i, ll) = dot_product(a_char(1:ielem(n, iconsidered)%idegfree, i, ll), ba_char(1:ielem(n, iconsidered)%idegfree, i, ll))
        end do; end do
      do ll = 2, ielem(n, iconsidered)%admis; do i = 1, mm
          !x_char(i,ll) = dot(a_char(1:idegfree2,i,ll), ba_char(1:idegfree2,i,ll))
          x_char(i, ll) = dot_product(a_char(1:idegfree2, i, ll), ba_char(1:idegfree2, i, ll))
        end do; end do
    else
      do ll = 1, ielem(n, iconsidered)%admis
            ba_char(1:ielem(n,iconsidered)%idegfree,1:nof_variables,ll)=matmul(b_char(1:ielem(n,iconsidered)%idegfree,1:ielem(n,iconsidered)%idegfree),a_char(1:ielem(n,iconsidered)%idegfree,1:nof_variables,ll))
      end do
      do ll = 1, ielem(n, iconsidered)%admis; do i = 1, mm
          !x_char(i,ll) = dot(a_char(:,i,ll), ba_char(:,i,ll))
      x_char(i, ll) = dot_product(a_char(1:ielem(n, iconsidered)%idegfree, i, ll), ba_char(1:ielem(n, iconsidered)%idegfree, i, ll))
        end do; end do
    end if

    deallocate(ba_char)
  end subroutine diag_at_b_a

  subroutine compute_gradcharv_smoothindicator(iconsidered, facex, eigvl, gradcharv, smoothindicator)
    implicit none
    integer, intent(in)::iconsidered, facex
    integer:: ll, k, i, l, ifds
    real::lwcx1
    real, dimension(1:nof_variables, 1:nof_variables), intent(in)::eigvl
    real, allocatable, dimension(:, :, :), intent(inout)::gradcharv
    real, allocatable, dimension(:, :, :, :), intent(inout)::smoothindicator
    real, allocatable, dimension(:)::lamc
    real, allocatable, dimension(:, :)::grad5alc, b_char, x_char
    real, allocatable, dimension(:, :, :)::a_char, gradients, gradients_eigvlt

    allocate(lamc(1:typesten))
    allocate(grad5alc(1:idegfree, 1:nof_variables))
    allocate(a_char(1:idegfree, 1:nof_variables, 1:typesten), b_char(1:idegfree, 1:idegfree))
    allocate(x_char(1:nof_variables, 1:typesten))
    allocate(gradients(0:idegfree, 1:nof_variables, 1:typesten))
    allocate(gradients_eigvlt(0:idegfree, 1:nof_variables, 1:typesten))

    i = iconsidered
    l = facex
    gradcharv = zero
    lwcx1 = ielem(n, i)%linc

    gradients(:, :, :) = zero; gradients_eigvlt(:, :, :) = zero
    do ll = 1, ielem(n, i)%admis
      gradients(0, :, ll) = u_c(i)%val(1, 1:nof_variables)
    end do
    if (ees .eq. 5) then
      ; lamc(:) = zero; grad5alc = zero; lamc(1) = (1.0d0 - (1.0d0/lwcx1)); lamc(2:ielem(n, i)%admis) = (1.0d0 - lamc(1))/(ielem(n, i)%admis - 1)
      do ll = 2, ielem(n, i)%admis
        grad5alc(1:idegfree2, 1:nof_variables) = grad5alc(1:idegfree2, 1:nof_variables) &
                                                + (lamc(ll)*ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, 1:nof_variables))
        gradients(1:idegfree2, 1:nof_variables, ll) = ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, 1:nof_variables)
      end do
      do ll = 1, 1
        gradients(1:ielem(n, i)%idegfree, 1:nof_variables, ll) = (1.0d0/lamc(1))* &
     (ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,1:nof_variables)-grad5alc(1:ielem(n,i)%idegfree,1:nof_variables))
      end do
      do ll = 1, 1
        gradients_eigvlt(0:ielem(n,i)%idegfree,1:nof_variables,ll)=matmul(gradients(0:ielem(n,i)%idegfree,1:nof_variables,ll),transpose(eigvl(1:nof_variables,1:nof_variables)))

      end do
      do ll = 2, ielem(n, i)%admis
        gradients_eigvlt(0:idegfree2,1:nof_variables,ll)=matmul(gradients(0:idegfree2,1:nof_variables,ll),transpose(eigvl(1:nof_variables,1:nof_variables)))
      end do
    else
      do ll = 1, ielem(n, i)%admis
        gradients(1:ielem(n, i)%idegfree, :, ll) = ilocal_recon5(iconsidered)%gradients(ll, 1:ielem(n, i)%idegfree, 1:nof_variables)
      end do
      do ll = 1, ielem(n, i)%admis
        gradients_eigvlt(0:ielem(n,i)%idegfree,1:nof_variables,ll)=matmul(gradients(0:ielem(n,i)%idegfree,1:nof_variables,ll),transpose(eigvl(1:nof_variables,1:nof_variables)))
      end do
    end if
    if (ees .ne. 5) then
      do ll = 1, ielem(n, i)%admis; do k = 0, ielem(n, i)%idegfree
          gradcharv(1:nof_variables, ll, k) = gradients_eigvlt(k, 1:nof_variables, ll)
      end do; end do
      a_char(1:ielem(n,i)%idegfree,1:nof_variables,1:ielem(n,i)%admis)=gradients_eigvlt(1:ielem(n,i)%idegfree,1:nof_variables,1:ielem(n,i)%admis)
 b_char(1:ielem(n, i)%idegfree, 1:ielem(n, i)%idegfree) = ilocal_recon3(i)%indicator(1:ielem(n, i)%idegfree, 1:ielem(n, i)%idegfree)
      call diag_at_b_a(iconsidered, a_char, b_char, x_char)
      smoothindicator(1:nof_variables, 1:ielem(n, i)%admis, l, 1) = x_char(1:nof_variables, 1:ielem(n, i)%admis)
    else
      do ll = 1, 1; do k = 0, ielem(n, i)%idegfree
          gradcharv(1:nof_variables, ll, k) = gradients_eigvlt(k, 1:nof_variables, ll)
        end do; end do
      do ll = 2, ielem(n, i)%admis; do k = 0, idegfree2
          gradcharv(1:nof_variables, ll, k) = gradients_eigvlt(k, 1:nof_variables, ll)
        end do; end do
      a_char(1:ielem(n,i)%idegfree,1:nof_variables,1:ielem(n,i)%admis)=gradients_eigvlt(1:ielem(n,i)%idegfree,1:nof_variables,1:ielem(n,i)%admis)
  b_char(1:ielem(n, i)%idegfree, 1:ielem(n, i)%idegfree) = ilocal_recon3(i)%indicator(1:ielem(n, i)%idegfree, 1:ielem(n, i)%idegfree)
      call diag_at_b_a(iconsidered, a_char, b_char, x_char)!,                                     &
      smoothindicator(1:nof_variables, 1:ielem(n, i)%admis, l, 1) = x_char(1:nof_variables, 1:ielem(n, i)%admis)
    end if
    deallocate(lamc, grad5alc, a_char, b_char, x_char, gradients, gradients_eigvlt)
  end subroutine
  subroutine find_bounds(iconsidered, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex, ik, iq, jk
    integer, intent(in)::iconsidered
    real::leftv(1:nof_variables), mp_pinfl, gammal
    real, allocatable, dimension(:, :), intent(inout)::utemp
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout)::maxvars, aver_vars, sumvars, utmin, utmax
    utemp = zero
    aver_vars = zero; sumvars = zero; maxvars = zero
    i = iconsidered
    if (extended_bounds .eq. 0) then        !strict bounds
      utemp(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      if (turbulenceequations .ge. 1) then
        utemp(1, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
          u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
      end if
      k = 1
      if (ielem(n, i)%interior .eq. 0) then
        do l = 1, ielem(n, i)%ifca
          k = k + 1
          utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
          if (turbulenceequations .ge. 1) then
            utemp(k, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
          end if
        end do
      end if
      if (ielem(n, i)%interior .eq. 1) then
        do l = 1, ielem(n, i)%ifca
          if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                k = k + 1
                utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
                if (turbulenceequations .ge. 1) then
                  utemp(k, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                    u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                end if
              else
                !not periodic ones in my cpu
              end if
            else
              k = k + 1
              utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
              if (turbulenceequations .ge. 1) then
                utemp(k, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                  u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
              end if
            end if
          else        !in other cpus they can only be periodic or mpi neighbours
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                k = k + 1
                utemp(k, 1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                            (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                if (turbulenceequations .ge. 1) then
                  utemp(k, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                    iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
           (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
                end if
              end if
            else
              k = k + 1
              utemp(k, 1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              if (turbulenceequations .ge. 1) then
                utemp(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexn(1,ielem(n,i)%indexi(l)))%sol&
           (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
              end if
            end if
          end if
        end do
      end if
      if (wenwrt .eq. 3) then
        do jk = 1, k
          leftv(1:nof_variables) = utemp(jk, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          utemp(jk, 1:nof_variables) = leftv(1:nof_variables)
        end do
      end if
    end if
    if (extended_bounds .eq. 1) then   !extended bounds
      k = 0
      if (ilocal_recon3(i)%local .eq. 1) then
        do iq = 1, ielem(n, i)%inumneighbours
          k = k + 1
          utemp(k, 1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq))%val(1, 1:nof_variables)
          if (turbulenceequations .ge. 1) then
            utemp(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(1,iq))%val(1,1:turbulenceequations+passivescalar)
          end if
        end do
      else
        do iq = 1, ielem(n, i)%inumneighbours
          k = k + 1
          if (ilocal_recon3(i)%ihexb(1, iq) .eq. n) then
            utemp(k, 1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq))%val(1, 1:nof_variables)
            if (turbulenceequations .ge. 1) then
              utemp(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_c(ilocal_recon3(i)%ihexl(1,iq))%val(1,1:turbulenceequations+passivescalar)
            end if
          else
            utemp(k, 1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, iq))%sol(ilocal_recon3(i)%ihexl(1, iq), 1:nof_variables)

            if (turbulenceequations .ge. 1) then
              utemp(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexn(1,iq))%sol(ilocal_recon3(i)%ihexl(1,iq),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if
          end if
        end do
      end if
      if (wenwrt .eq. 3) then
        do jk = 1, k
          leftv(1:nof_variables) = utemp(jk, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          utemp(jk, 1:nof_variables) = leftv(1:nof_variables)
        end do
      end if
    end if
    utmin = zero; utmax = zero
    do iex = 1, nof_variables + turbulenceequations + passivescalar
      utmin(iex) = minval(utemp(1:k, iex))
      utmax(iex) = maxval(utemp(1:k, iex))
      do ik = 2, k
        sumvars(iex) = sumvars(iex) + abs(utemp(ik, iex) - utemp(1, iex))
      end do
      do ik = 1, k
        aver_vars(iex) = aver_vars(iex) + utemp(ik, iex)
      end do
      aver_vars(iex) = aver_vars(iex)/k
      do ik = 1, k
        maxvars(iex) = max(maxvars(iex), abs(utemp(ik, iex)))
      end do
    end do
  end subroutine find_bounds
  subroutine compute_muscl_reconstruction(iconsidered, utmin, utmax, utemp)
    real, allocatable, dimension(:, :, :)::usol, psi
    real::ax, ay, az, mp_pinfl, gammal, limvbg
    integer, intent(in)::iconsidered
    integer::icd, icompwrt, i, ngp, l, iex, iqp
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:, :), intent(in)::utemp
    real, dimension(1:nof_variables), intent(in)::utmin, utmax
    real, allocatable, dimension(:, :)::consmatrix, gradssl, ressolution, gradssl2, ressolution2
    real, allocatable, dimension(:)::slope
    i = iconsidered
    ilocal_recon3(iconsidered)%uleft(:, :, :) = zero
    allocate(slope(1:nof_variables + turbulenceequations + passivescalar))
    allocate(usol(1:nof_variables + turbulenceequations + passivescalar, 1:6, 1:numberofpoints2))
    allocate(psi(nof_variables + turbulenceequations + passivescalar, 1:6, 1:numberofpoints2))
    allocate(consmatrix(1:6*numberofpoints2, 1:idegfree))
    allocate(gradssl(1:idegfree, 1:nof_variables))
    allocate(ressolution(1:6*numberofpoints2, 1:nof_variables))
    if (turbulenceequations .ge. 1) then
      allocate(gradssl2(1:idegfree, 1:turbulenceequations + passivescalar))
      allocate(ressolution2(1:6*numberofpoints2, 1:turbulenceequations + passivescalar))
    end if
    icompwrt = 0
    usol(:, :, :) = zero
    psi = zero
    icd = 0
    do l = 1, ielem(n, i)%ifca        !faces2
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      do ngp = 1, iqp                        !for gqp
        ax = ilocal_recon3(i)%qpoints(l, ngp, 1)
        ay = ilocal_recon3(i)%qpoints(l, ngp, 2)
        if (dimensiona .eq. 3) then
          az = ilocal_recon3(i)%qpoints(l, ngp, 3)
        end if
        icd = icd + 1
        if (dimensiona .eq. 3) then
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        else
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec2d(n, ax, ay, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        end if
      end do
    end do

  gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = ilocal_recon5(iconsidered)%gradients(1, 1:ielem(n, i)%idegfree, 1:nof_variables)
  ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
    if (turbulenceequations .ge. 1) then
    gradssl2(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=ilocal_recon5(iconsidered)%gradients2(1,1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)
    ressolution2(1:icd,1:turbulenceequations+passivescalar)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl2(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar))
    end if
    icd = 0              !initialise counter
    do l = 1, ielem(n, i)%ifca     !loop all faces
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      do ngp = 1, iqp        !all gaussian quadrature points
        icd = icd + 1
        if (wenwrt .eq. 3) then
          leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          usol(1:nof_variables, l, ngp) = ((leftv(1:nof_variables) + ressolution(icd, 1:nof_variables)))
        else
          usol(1:nof_variables, l, ngp) = ((u_c(i)%val(1, 1:nof_variables) + ressolution(icd, 1:nof_variables)))
        end if
        if (turbulenceequations .ge. 1) then
          usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,l,ngp)=((u_ct(i)%val(1,1:turbulenceequations+passivescalar)+ressolution2(icd,1:turbulenceequations+passivescalar)))
        end if
      end do
    end do
    do l = 1, ielem(n, i)%ifca        !faces2
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      do ngp = 1, iqp
        do iex = 1, nof_variables + turbulenceequations + passivescalar
          call slope_limiters(n, i, iex, l, ngp, utmin, utmax, utemp, usol, psi)
        end do
      end do
    end do
    do iex = 1, nof_variables + turbulenceequations + passivescalar
      limvbg = tolbig
      do l = 1, ielem(n, i)%ifca        !faces2
        if (dimensiona .eq. 3) then
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
        else
        iqp = qp_line
        end if
        do ngp = 1, iqp
          limvbg = min(limvbg, psi(iex, l, ngp))
        end do
      end do
      slope(iex) = limvbg
      if (iex .eq. 1) then
        ielem(n, i)%wcx(1) = limvbg
      end if
      if (dg .eq. 1) then
    ilocal_recon6(i)%dg2fv(1:ielem(n, i)%idegfree, iex) = ilocal_recon5(iconsidered)%gradients(1, ielem(n, i)%idegfree, iex)*slope(iex)
      end if
    end do
    ilocal_recon3(i)%uleft(:, :, :) = zero
    if (turbulenceequations .ge. 1) then
      ilocal_recon3(i)%uleftturb(:, :, :) = zero
    end if
    do l = 1, ielem(n, i)%ifca        !faces2
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      do ngp = 1, iqp
        if (wenwrt .eq. 3) then
          leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          usol(1:nof_variables, l, ngp) = usol(1:nof_variables, l, ngp) - leftv(1:nof_variables)
        else
          usol(1:nof_variables, l, ngp) = usol(1:nof_variables, l, ngp) - u_c(i)%val(1, 1:nof_variables)
        end if
        if (turbulenceequations .ge. 1) then
          usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,l,ngp)=usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,l,ngp)&
          - u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
        end if
        call extrapolate_bound_muscl(usol, iex, l, ngp, i, slope)
      end do
    end do

    deallocate(slope)
    deallocate(usol)
    deallocate(psi)
    deallocate(consmatrix)
    deallocate(gradssl)
    deallocate(ressolution)

    if (turbulenceequations .ge. 1) then
      deallocate(gradssl2)
      deallocate(ressolution2)
    end if

  end subroutine compute_muscl_reconstruction

  subroutine muscl(n)
!> @brief
!> subroutine for muscl type reconstruction
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, iex, ieul, jx, lx, kmaxe, iq, ll, ngp, nnd, iqp, idummy, ii, icd, iconsidered
    real::umin, umax, psitot, addc, divg0, limvbg, tempxx
    real, allocatable, dimension(:, :)::utemp
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::maxvars, aver_vars, sumvars, utmin, utmax
    allocate(utemp(imaxdegfree + 1, 1:nof_variables + turbulenceequations + passivescalar))
    kmaxe = xmpielrank(n)
    do ii = 1, nof_interior
      i = el_int(ii)
      iconsidered = i
      if (((ielem(n,i)%troubled.eq.1).and.(ielem(n,i)%reduce.eq.1)).or.((ielem(n,i)%full.eq.0).and.(ielem(n,i)%troubled.eq.1)))then
        if (ielem(n, i)%recalc .gt. 0) then
          if (adda .eq. 1) then
            call adda_filter(n, iconsidered)
          end if
          call find_bounds(iconsidered, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
          call compute_muscl_reconstruction(iconsidered, utmin, utmax, utemp)
        end if
      end if
    end do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      if (adda .eq. 1) then
        call adda_filter(n, iconsidered)
      end if
      if (((ielem(n,i)%troubled.eq.1).and.(ielem(n,i)%reduce.eq.1)).or.((ielem(n,i)%full.eq.0).and.(ielem(n,i)%troubled.eq.1)))then
        if (ielem(n, i)%recalc .gt. 0) then
          if (adda .eq. 1) then
            call adda_filter(n, iconsidered)
          end if
          call find_bounds(iconsidered, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
          call compute_muscl_reconstruction(iconsidered, utmin, utmax, utemp)
        end if
      end if
    end do
    deallocate(utemp)
  end subroutine muscl
  subroutine solutiontriav2(n)
    implicit none
    integer, intent(in)::n
    integer::i,j,k,l,m,ppp,ieul,iex,ihgt,ihgj,kmaxe,decomf,icnn,iqdr,nvar,idummy,iqp,nnd,ngp,icd,icompwrt,iconsidered
    real::raa1, raa2, paa1, paa2, ax, ay, az
    real::solx
    real, dimension(1:dimensiona)::ugradloc
    real, dimension(1:dimensiona, 1:dimensiona)::ainvjt
    real, allocatable, dimension(:)::gradtem
    real, allocatable, dimension(:, :)::xxder, yyder, zzder
    kmaxe = xmpielrank(n);
    allocate(xxder(1:idegfree, 1:numberofpoints2))
    allocate(yyder(1:idegfree, 1:numberofpoints2))
    allocate(zzder(1:idegfree, 1:numberofpoints2))
    allocate(gradtem(1:idegfree))
    do i = 1, kmaxe
      iconsidered = i
      ilocal_recon3(i)%uleftv(:, :, :, :) = zero;
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        ilocal_recon3(i)%uleftturbv(:, :, :, :) = zero; ilocal_recon3(i)%uleftturb(:, :, :) = zero;
      end if
      do ihgt = 1, dimensiona; do ihgj = 1, dimensiona
          ainvjt(ihgt, ihgj) = ilocal_recon3(i)%invccjac(ihgj, ihgt)
        end do; end do
      do l = 1, ielem(n, i)%ifca; idummy = 0
        if (dimensiona .eq. 3) then
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad;
        else
          iqp = qp_triangle;
        end if
        else
        iqp = qp_line;
        end if
        icd = 0
        do ngp = 1, iqp                        !for gqp
          ax = ilocal_recon3(i)%qpoints(l, ngp, 1);
          ay = ilocal_recon3(i)%qpoints(l, ngp, 2);
          if (dimensiona .eq. 3) then
            az = ilocal_recon3(i)%qpoints(l, ngp, 3)
          end if
          icd = icd + 1
          if (dimensiona .eq. 3) then
          do k = 1, ielem(n, i)%idegfree
            if (poly .eq. 1) then
              xxder(k, icd) = dfx(ax, ay, az, k, i); yyder(k, icd) = dfy(ax, ay, az, k, i); zzder(k, icd) = dfz(ax, ay, az, k, i)
            end if
            if (poly .eq. 2) then
              xxder(k, icd) = dlx(ax, ay, az, k, i); yyder(k, icd) = dly(ax, ay, az, k, i); zzder(k, icd) = dlz(ax, ay, az, k, i)
            end if
            if (poly .eq. 4) then
              xxder(k, icd) = tl3dx(ax, ay, az, k, i); yyder(k, icd) = tl3dy(ax, ay, az, k, i); zzder(k, icd) = tl3dz(ax, ay, az, k, i)
            end if
          end do
          else
          do k = 1, ielem(n, i)%idegfree
            if (poly .eq. 4) then
              xxder(k, icd) = tl2dx(ax, ay, k, i); yyder(k, icd) = tl2dy(ax, ay, k, i);
            else
              xxder(k, icd) = df2dx(ax, ay, k, i); yyder(k, icd) = df2dy(ax, ay, k, i);
            end if
          end do
          end if
        end do
        icd = 0
        do ngp = 1, iqp
          icd = icd + 1
          select case (ielem(n, i)%ggs)
          case (0)
            !turbulence first
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then
                do nvar = 1, turbulenceequations + passivescalar
                  ilocal_recon3(i)%uleftturb(nvar, l, ngp) = u_ct(i)%val(1, nvar)
                end do
              end if
              do nvar = 1, turbulenceequations + passivescalar
                gradtem(1:ielem(n, i)%idegfree) = ilocal_recon5(iconsidered)%gradientsturb(1, 1:ielem(n, i)%idegfree, nvar)
                ugradloc = zero
                ugradloc(1) = dot_product(gradtem(1:ielem(n, i)%idegfree), xxder(1:ielem(n, i)%idegfree, icd))
                ugradloc(2) = dot_product(gradtem(1:ielem(n, i)%idegfree), yyder(1:ielem(n, i)%idegfree, icd))
                if (dimensiona .eq. 3) then
                  ugradloc(3) = dot_product(gradtem(1:ielem(n, i)%idegfree), zzder(1:ielem(n, i)%idegfree, icd))
                end if
        ilocal_recon3(i)%uleftturbv(1:dimensiona, nvar, l, ngp) = matmul(ainvjt(1:dimensiona, 1:dimensiona), ugradloc(1:dimensiona))
              end do
            end if
            !now temperature
            gradtem(1:ielem(n, i)%idegfree) = ilocal_recon5(iconsidered)%gradientstemp(1:ielem(n, i)%idegfree)
            ugradloc = zero
            ugradloc(1) = dot_product(gradtem(1:ielem(n, i)%idegfree), xxder(1:ielem(n, i)%idegfree, icd))
            ugradloc(2) = dot_product(gradtem(1:ielem(n, i)%idegfree), yyder(1:ielem(n, i)%idegfree, icd))
            if (dimensiona .eq. 3) then
              ugradloc(3) = dot_product(gradtem(1:ielem(n, i)%idegfree), zzder(1:ielem(n, i)%idegfree, icd))
            end if
            ilocal_recon3(i)%uleftv(1:dimensiona, 1, l, ngp) = matmul(ainvjt(1:dimensiona, 1:dimensiona), ugradloc(1:dimensiona))
            !now velocities
            do iex = 1, dimensiona
              gradtem(1:ielem(n, i)%idegfree) = ilocal_recon5(iconsidered)%velocitydof(iex, 1:ielem(n, i)%idegfree)
              ugradloc = zero
              ugradloc(1) = dot_product(gradtem(1:ielem(n, i)%idegfree), xxder(1:ielem(n, i)%idegfree, icd))
              ugradloc(2) = dot_product(gradtem(1:ielem(n, i)%idegfree), yyder(1:ielem(n, i)%idegfree, icd))
              if (dimensiona .eq. 3) then
                ugradloc(3) = dot_product(gradtem(1:ielem(n, i)%idegfree), zzder(1:ielem(n, i)%idegfree, icd))
              end if
         ilocal_recon3(i)%uleftv(1:dimensiona, iex + 1, l, ngp) = matmul(ainvjt(1:dimensiona, 1:dimensiona), ugradloc(1:dimensiona))
            end do
          case (1)
            !turbulence first
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then
                do nvar = 1, turbulenceequations + passivescalar
                  ilocal_recon3(i)%uleftturb(nvar, l, ngp) = u_ct(i)%val(1, nvar)
                end do
              end if
              do nvar = 1, turbulenceequations + passivescalar
               ilocal_recon3(i)%uleftturbv(1:dimensiona, nvar, l, ngp) = ilocal_recon3(i)%grads(dimensiona + 1 + nvar, 1:dimensiona)
              end do
            end if
            ilocal_recon3(i)%uleftv(1:dimensiona, 1, l, ngp) = ilocal_recon3(i)%grads(dimensiona + 1, 1:dimensiona)
            do iex = 1, dimensiona
              ilocal_recon3(i)%uleftv(1:dimensiona, iex + 1, l, ngp) = ilocal_recon3(i)%grads(iex, 1:dimensiona)
            end do
          end select
        end do
      end do
      if (ielem(n, iconsidered)%ggs .eq. 0) then
        call compute_gradients_center(n, iconsidered)
      end if
    end do
    deallocate(xxder)
    deallocate(yyder)
    deallocate(zzder)
    deallocate(gradtem)
  end subroutine solutiontriav2

  subroutine least_squares(n)
    implicit none
    integer, intent(in)::n
    integer::iconsidered, ii, i

!$omp do
    do ii = 1, nof_interior;
      i = el_int(ii)
      iconsidered = i
      call allgrads_inner(n, i)
    end do
!$omp end do

!$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      call allgrads_mix(n, i)
    end do
!$omp end do
  end subroutine least_squares
  subroutine piecewise_constant(n)
    implicit none
    integer, intent(in)::n
    integer::i, k, kmaxe, iex, l
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
    do iex = 1, nof_variables
      ilocal_recon3(i)%uleft(iex, :, :) = u_c(i)%val(1, iex)
    end do

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
    do iex = 1, turbulenceequations + passivescalar
      ilocal_recon3(i)%uleftturb(iex, :, :) = u_ct(i)%val(1, iex)
    end do
    end if
    end do
  end subroutine piecewise_constant

  subroutine linear_scheme(n)
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, iex, ieul, jx, lx, kmaxe, iq, ll, ngp, nnd, iqp, idummy, ii, icd, iconsidered
    real::umin, umax, psitot, addc, divg0, limvbg, tempxx
    kmaxe = xmpielrank(n)
    do ii = 1, nof_interior
      i = el_int(ii)
      iconsidered = i
      call compute_linear_reconstruction(iconsidered)
    end do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      call compute_linear_reconstruction(iconsidered)
    end do
  end subroutine linear_scheme
  subroutine compute_linear_reconstruction(iconsidered)
    real, allocatable, dimension(:, :, :)::usol, psi
    real::ax, ay, az, mp_pinfl, gammal, limvbg
    integer, intent(in)::iconsidered
    integer::icd, icompwrt, i, ngp, l, iex, iqp
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:, :)::consmatrix, gradssl, ressolution, gradssl2, ressolution2
    i = iconsidered
    ilocal_recon3(iconsidered)%uleft(:, :, :) = zero

    allocate(usol(1:nof_variables + turbulenceequations + passivescalar, 1:6, 1:numberofpoints2))
    allocate(consmatrix(1:6*numberofpoints2, 1:idegfree))
    allocate(gradssl(1:idegfree, 1:nof_variables))
    allocate(ressolution(1:6*numberofpoints2, 1:nof_variables))
    allocate(psi(1:nof_variables + turbulenceequations + passivescalar, 1:6, 1:numberofpoints2))

    if (turbulenceequations .ge. 1) then
      allocate(gradssl2(1:idegfree, 1:turbulenceequations + passivescalar))
      allocate(ressolution2(1:6*numberofpoints2, 1:turbulenceequations + passivescalar))
    end if

    icompwrt = 0

    usol(:, :, :) = zero
    psi = zero
    icd = 0

    do l = 1, ielem(n, i)%ifca        !faces2

      if (dimensiona .eq. 3) then
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if
      do ngp = 1, iqp                        !for gqp
        ax = ilocal_recon3(i)%qpoints(l, ngp, 1)
        ay = ilocal_recon3(i)%qpoints(l, ngp, 2)
        if (dimensiona .eq. 3) then
          az = ilocal_recon3(i)%qpoints(l, ngp, 3)
        end if
        icd = icd + 1
        if (dimensiona .eq. 3) then
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        else
          consmatrix(icd, 1:ielem(n, i)%idegfree) = basis_rec2d(n, ax, ay, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
        end if
      end do
    end do
  gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = ilocal_recon5(iconsidered)%gradients(1, 1:ielem(n, i)%idegfree, 1:nof_variables)


   ressolution(1:icd,1:nof_variables)=matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))

    if (turbulenceequations .ge. 1) then
  gradssl2(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=ilocal_recon5(iconsidered)%gradients2(1,1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)
    ressolution2(1:icd,1:turbulenceequations+passivescalar)= matmul(consmatrix(1:icd,1:ielem(n,i)%idegfree),gradssl2(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar))
    end if
    icd = 0              !initialise counter
    icd = 0              !initialise counter
    do l = 1, ielem(n, i)%ifca     !loop all faces
      if (dimensiona .eq. 3) then
      if (ielem(n, i)%types_faces(l) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
      else
      iqp = qp_line
      end if
      do ngp = 1, iqp        !all gaussian quadrature points
        icd = icd + 1
        usol(1:nof_variables, l, ngp) = ressolution(icd, 1:nof_variables)
        if (turbulenceequations .ge. 1) then
   usol(nof_variables+1:nof_variables+turbulenceequations+passivescalar,l,ngp)=ressolution2(icd,1:turbulenceequations+passivescalar)
        end if
        call extrapolate_bound_linear(usol, iex, l, ngp, i)
      end do
    end do
    deallocate(usol)
    deallocate(psi)
    deallocate(consmatrix)
    deallocate(gradssl)
    deallocate(ressolution)
    if (turbulenceequations .ge. 1) then
      deallocate(gradssl2)
      deallocate(ressolution2)
    end if

  end subroutine compute_linear_reconstruction

  subroutine arbitrary_order(n)
!> @brief
!> subroutine controlling the reconstruction in 3d
    implicit none
    integer, intent(in)::n
    integer::kmaxe, i
    kmaxe = xmpielrank(n)
    ielem(n, 1:kmaxe)%reduce = 0
    call least_squares(n)
    select case (iweno)
    case (1)
      call wenoweights(n)
      call checksol(n)
      call muscl(n)
      call checksolx(n)
    case (-1)
      call muscl(n)
      call checksolx(n)
    case (0)
      if (firstorder .eq. 1) then
        call piecewise_constant(n)
      else
        call linear_scheme(n)
        call checksolx(n)
      end if
    end select
    if (itestcase .eq. 4) then
      call solutiontriav2(n)
    end if
  end subroutine arbitrary_order
  subroutine viscous_dg_ggs(n)
!> @brief
!> subroutine controlling the reconstruction in 3d
    implicit none
    integer, intent(in)::n
    integer::kmaxe, i, iex, l, ngp, iqp, iconsidered
    kmaxe = xmpielrank(n)
    if (dimensiona .eq. 3) then
!$omp do
      do i = 1, kmaxe
        iconsidered = i
        do l = 1, ielem(n, i)%ifca
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad
          else
            iqp = qp_triangle
          end if
          do ngp = 1, iqp
            ilocal_recon3(i)%uleftv(1:3, 1, l, ngp) = ilocal_recon3(i)%grads(4, 1:3)
            do iex = 1, 3
              ilocal_recon3(i)%uleftv(1:3, iex + 1, l, ngp) = ilocal_recon3(i)%grads(iex, 1:3)
            end do
          end do
        end do
      end do
!$omp end do
    else
!$omp do
      do i = 1, kmaxe
        iconsidered = i
        do l = 1, ielem(n, i)%ifca
          iqp = qp_line_n
          do ngp = 1, iqp
            ilocal_recon3(i)%uleftv(1:2, 1, l, ngp) = ilocal_recon3(i)%grads(3, 1:2)
            do iex = 1, 2
              ilocal_recon3(i)%uleftv(1:2, iex + 1, l, ngp) = ilocal_recon3(i)%grads(iex, 1:2)
            end do
          end do
        end do
      end do
!$omp end do
    end if
  end subroutine viscous_dg_ggs
  subroutine checksol(n)
    implicit none
    integer, intent(in)::n
    integer::i, l, ngp, iqp, iex
    integer::reduce1, kmaxe, indx
    real::jump_cond
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    kmaxe = xmpielrank(n)
    jump_cond = 0.85
    if (itestcase .ge. 3) then
!$omp do
      do i = 1, kmaxe        !all elements
        ielem(n, i)%reduce = 0; reduce1 = 0

        if (ielem(n, i)%troubled .eq. 1) then

          if (ielem(n, i)%full .eq. 0) then
            ielem(n, i)%reduce = 1
          end if

          if (ielem(n, i)%full .eq. 1) then

            do l = 1, ielem(n, i)%ifca        !faces2
              if (dimensiona .eq. 3) then
                indx = 5
                if (ielem(n, i)%types_faces(l) .eq. 5) then
                  iqp = qp_quad
                else
                  iqp = qp_triangle
                end if
              else
                indx = 4
                iqp = qp_line
              end if

              do ngp = 1, iqp

                leftv(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
                rightv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
                call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

                if (((abs(leftv(1) - rightv(1))) .ge. (jump_cond*rightv(1)))) then
                  reduce1 = 1
                  ielem(n, i)%reduce = 1
                end if

                if (((abs(leftv(indx) - rightv(indx))) .ge. (jump_cond*rightv(indx)))) then
                  reduce1 = 1
                  ielem(n, i)%reduce = 1
                end if

                if ((leftv(indx) .lt. 0.0) .or. (leftv(1) .lt. 0.0)) then
                  reduce1 = 1
                  ielem(n, i)%reduce = 1
                end if

              end do
            end do

            if (ielem(n, i)%hybrid .eq. 1) then
              reduce1 = 1
              ielem(n, i)%reduce = 1
            end if
          end if
        end if
      end do
    end if
  end subroutine checksol
  subroutine checksolx(n)
    implicit none
    integer, intent(in)::n
    integer::i, l, ngp, iqp, iex
    integer::reduce1, kmaxe, indx
    real::jump_cond
    real, dimension(1:nof_variables)::leftv
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::rightv
    real::mp_pinfr, gammar
    kmaxe = xmpielrank(n)
    jump_cond = 0.95
    if (itestcase .ge. 3) then
      do i = 1, kmaxe        !all elements
        reduce1 = 0
        if (ielem(n, i)%troubled .eq. 1) then
          do l = 1, ielem(n, i)%ifca        !faces2
            if (dimensiona .eq. 3) then
              indx = 5
              if (ielem(n, i)%types_faces(l) .eq. 5) then
                iqp = qp_quad
              else
                iqp = qp_triangle
              end if
            else
              indx = 4
              iqp = qp_line
            end if
            do ngp = 1, iqp
              leftv(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
              rightv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
              call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
              if (((abs(leftv(1) - rightv(1))) .ge. (jump_cond*rightv(1)))) then
                reduce1 = 1
                ielem(n, i)%reduce = 1
              end if
              if (((abs(leftv(indx) - rightv(indx))) .ge. (jump_cond*rightv(indx)))) then
                reduce1 = 1
                ielem(n, i)%reduce = 1
              end if
              if ((leftv(indx) .lt. 0.0) .or. (leftv(1) .lt. 0.0)) then
                reduce1 = 1
                ielem(n, i)%reduce = 1
              end if
            end do
          end do
          if (ielem(n, i)%hybrid .eq. 1) then
            reduce1 = 1
            ielem(n, i)%reduce = 1
          end if
          if (reduce1 .eq. 1) then
            do iex = 1, nof_variables
              ilocal_recon3(i)%uleft(iex, :, :) = u_c(i)%val(1, iex)
            end do
            if (dg .eq. 1) then
              ilocal_recon6(i)%dg2fv(1:ielem(n, i)%idegfree, :) = zero
            end if
          end if
          if (turbulence .eq. 1) then
            if (icoupleturb .eq. 1) then
              reduce1 = 0
              do l = 1, ielem(n, i)%ifca        !faces2
                if (dimensiona .eq. 3) then
                  indx = 5
                  if (ielem(n, i)%types_faces(l) .eq. 5) then
                    iqp = qp_quad
                  else
                    iqp = qp_triangle
                  end if
                else
                  indx = 4
                  iqp = qp_line
                end if
                do ngp = 1, iqp
                  leftv(1) = ilocal_recon3(i)%uleftturb(1, l, ngp)
                  rightv(1) = u_ct(i)%val(1, 1)
                  if (((abs(leftv(1) - rightv(1))) .ge. (0.6*rightv(1))) .or. (leftv(1) .le. zero)) then
                    reduce1 = 1
                  end if
                end do
              end do
              if (ielem(n, i)%hybrid .eq. 1) then
                reduce1 = 1
                ielem(n, i)%reduce = 1
              end if
              if (reduce1 .eq. 1) then
              do iex = 1, 1
                ilocal_recon3(i)%uleftturb(1, :, :) = u_ct(i)%val(1, 1)
              end do
              end if
            end if
          end if

        end if
      end do
    end if
  end subroutine checksolx
  subroutine slope_limiters(n, iconsidered, icons_e, facex, icons_s, utmin, utmax, utemp, usol, psi)
!> @brief
!> subroutine muscl slope limiters
    implicit none
    integer, intent(in)::n, icons_e, facex, icons_s, iconsidered
    real, allocatable, dimension(:, :), intent(in)::utemp
    real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(in)::utmin, utmax
    real, allocatable, dimension(:, :, :), intent(in)::usol
    real, allocatable, dimension(:, :, :), intent(inout)::psi
    integer::i, l, iex, ngp
    real::d2, sfd, epsi2, dmin, dplus, kappa_ven, psi2, sig_1, delu, y_fun, s_y, pol_mog
    iex = icons_e
    l = facex
    i = iconsidered
    ngp = icons_s
    kappa_ven = 10.0d0
    psi2 = zero
    d2 = usol(iex, l, ngp) - utemp(1, iex)
    if (abs(d2) .le. zero) then
      psi(iex, l, ngp) = 1.0d0
    else if (d2 .gt. zero) then
      sfd = (utmax(iex) - utemp(1, iex))/d2
      select case (limiter)
      case (1)
        psi(iex, l, ngp) = min(1.0d0, sfd)                !barth and jespersen
      case (2)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0
        end if
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
        if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
          y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
          s_y = (2.0d0*y_fun**3) - (3.0d0*y_fun**2) + 1.0d0
          sig_1 = s_y
        else
          sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (9)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0
        end if
        !barth and jespersen michalak
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
      if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
        y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
        s_y = (2.0d0*y_fun**3) - (3.0d0*y_fun**2) + 1.0d0
        sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (3)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0

        end if        !barth and jespersen michalak extended
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
      if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
        y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
        s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (4)                  !vkm
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmax(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
      psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
        if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
            sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (5)
        psi(iex, l, ngp) = (sfd**2 + sfd)/(sfd**2 + 1.0d0)                           ! van albada
      case (6)
        psi(iex, l, ngp) = 2.0d0*sfd/(sfd + 1.0d0)                                 ! van leer
      case (7)                                                                !venkatakrishnan
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmax(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
        psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
      case (8)                                                                !venkatakrishnan
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmax(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
      psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
       if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
            sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      end select
    else
      sfd = ((utmin(iex) - utemp(1, iex)))/(d2)
      select case (limiter)
      case (1)
        psi(iex, l, ngp) = min(1.0d0, sfd)                                !minmod limiter
      case (2)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0
        end if                                        !barth and jespersen
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
       if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
            sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (9)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd

        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0
        end if                                        !barth and jespersen
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
       if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then

            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)

            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0

            sig_1 = s_y

          else

            sig_1 = 0.0d0

          end if
        end if

        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (3)
        pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
        if (sfd .lt. 1.5d0) then
          psi(iex, l, ngp) = pol_mog
        else
          psi(iex, l, ngp) = 1.0d0
        end if                                        !barth and jespersen
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
       if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
            sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      case (4)
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmin(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
        psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
       if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
            y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
            s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
            sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2                                 !barth and jespersen (extended stencils)
      case (5)
        psi(iex, l, ngp) = (sfd**2 + sfd)/(sfd**2 + 1.0d0)                           ! van albada
      case (6)
        psi(iex, l, ngp) = 2.0d0*sfd/(sfd + 1.0d0)                                         ! van leer
      case (7)                                                                !venkatakrishnan
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmin(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
        psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
      case (8)                                                                !venkatakrishnan
        dmin = usol(iex, l, ngp) - utemp(1, iex)
        dmin = sign(1.0d0, dmin)*(abs(dmin) + tolsmall)
        dplus = utmin(iex) - utemp(1, iex)
        epsi2 = (kappa_ven*ielem(n, i)%minedge)**3
        psi(iex,l,ngp)=(1.0d0/(dmin))*(((((dplus**2)+epsi2)*dmin)+(2.0d0*(dmin**2)*dplus))/((dplus**2)+(2.0d0*dmin**2)+(dmin*dplus)+epsi2))
        delu = utmax(iex) - utmin(iex)
        if (delu**2 .le. ((kappa_ven*ielem(n, i)%minedge)**3)) then
          sig_1 = 1.0d0
        else
        if ((((kappa_ven*ielem(n, i)%minedge)**3) .lt. delu**2) .and. (delu**2 .lt. 2.0d0*((kappa_ven*ielem(n, i)%minedge)**3))) then
          y_fun = ((delu**2) - ((kappa_ven*ielem(n, i)%minedge)**3))/((kappa_ven*ielem(n, i)%minedge)**3)
          s_y = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
          sig_1 = s_y
          else
            sig_1 = 0.0d0
          end if
        end if
        psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex, l, ngp)
        psi(iex, l, ngp) = psi2
      end select
    end if
  end subroutine slope_limiters
  subroutine trouble_indicator1
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex
    integer::trouble
    integer::iconsidered, facex, pointx
    real, dimension(1:nof_variables)::leftv, rightv
    real, dimension(1:nof_variables)::maxvars, aver_vars, sumvars, utmin, utmax
    real, allocatable, dimension(:, :, :)::usol
    real, allocatable, dimension(:, :)::utemp
    allocate(utemp(imaxdegfree + 1, 1:nof_variables + turbulenceequations + passivescalar))
    allocate(usol(1:nof_variables + turbulenceequations + passivescalar, 1:6, 1:numberofpoints2))
    kmaxe = xmpielrank(n)
    if (code_profile .ne. 102) then

!$omp do
      do i = 1, kmaxe
        iconsidered = i
        call find_bounds(iconsidered, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
        do l = 1, ielem(n, i)%ifca
          if (dimensiona .eq. 2) then
            iqp = qp_line_n
          else
            if (ielem(n, i)%types_faces(l) .eq. 5) then
              iqp = qp_quad
            else
              iqp = qp_triangle
            end if
          end if
          do ngp = 1, iqp!
            facex = l
            pointx = ngp
            usol(:, facex, pointx) = ilocal_recon3(iconsidered)%uleft_dg(:, facex, pointx)
            leftv(:) = usol(:, facex, pointx)
            call pad_dg(iconsidered, leftv)
            call nad_dg(iconsidered, facex, pointx, leftv, rightv, usol, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
          end do
        end do
      end do
!$omp end do
    end if
    deallocate(usol)
    deallocate(utemp)
  end subroutine
  subroutine trouble_indicator2
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex, ndof
    integer::trouble, ifree, i_deg
    integer::iconsidered, facex, pointx
    kmaxe = xmpielrank(n)
    if (code_profile .ne. 102) then
!$omp do
      do i = 1, kmaxe
        iconsidered = i
        if (ielem(n, i)%troubled .eq. 1) then
          do l = 1, ielem(n, i)%ifca
            if (dimensiona .eq. 2) then
              iqp = qp_line_n
            else
              if (ielem(n, i)%types_faces(l) .eq. 5) then
                iqp = qp_quad
              else
                iqp = qp_triangle
              end if
            end if

            do ngp = 1, iqp!
              facex = l
              pointx = ngp
              ilocal_recon3(iconsidered)%uleft_dg(:, facex, pointx) = ilocal_recon3(iconsidered)%uleft(:, facex, pointx)
            end do
          end do
          do iex = 1, nof_variables
            u_c(iconsidered)%valdg(1, iex, 2:idegfree + 1) = ilocal_recon6(iconsidered)%dg2fv(1:idegfree, iex)
          end do
        end if
      end do
!$omp end do
    end if
  end subroutine
  subroutine pad_dg(iconsidered, leftv)
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex
    integer::trouble
    integer, intent(in)::iconsidered
    real, dimension(1:nof_variables), intent(inout)::leftv
    real::mp_pinfl, gammal
    i = iconsidered
    if (itestcase .ge. 3) then
      if (dimensiona .eq. 3) then
        call cons2prim(n, leftv, mp_pinfl, gammal)
        if (multispecies .eq. 1) then
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
          if ((leftv(5) .le. zero) .or. (leftv(5) .ne. leftv(5))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
        else
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1

          end if
          if ((leftv(5) .le. zero) .or. (leftv(5) .ne. leftv(5))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
        end if
      else
        call cons2prim(n, leftv, mp_pinfl, gammal)
        if (multispecies .eq. 1) then

          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
          if ((leftv(4) .le. zero)) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
          if ((leftv(4) .ne. leftv(4))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
        end if
          if ((leftv(nof_variables) .lt. zero) .or. leftv(nof_variables) .gt. 1.0d0) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
        else
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
          if ((leftv(4) .le. zero) .or. (leftv(4) .ne. leftv(4))) then
            ielem(n, i)%troubled = 1; ielem(n, i)%condition = 1
          end if
        end if
      end if
    end if
  end subroutine
  subroutine nad_dg(iconsidered, facex, pointx, leftv, rightv, usol, maxvars, aver_vars, sumvars, utmin, utmax, utemp)
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex
    integer::trouble, img
    integer, intent(in)::iconsidered, facex, pointx
    real::par1, par2, d2, minb, maxb
    real, dimension(1:nof_variables)::nad_dg_el
    real, dimension(1:nof_variables), intent(inout)::leftv, rightv
    real, dimension(1:nof_variables), intent(in)::maxvars, aver_vars, sumvars, utmin, utmax
    real, allocatable, dimension(:, :), intent(in)::utemp
    real, allocatable, dimension(:, :, :), intent(in)::usol
    real::mp_pinfl, gammal

    if (dimensiona .eq. 3) then

      img = 5
    else

      img = 4
    end if

    select case (indicator_type)

    case (1)       !mood indicator

      do iex = 1, nof_variables
        nad_dg_el(iex) = max(indicator_par1, (indicator_par2)*(utmax(iex) - utmin(iex)))
      end do
      leftv(1:nof_variables) = usol(1:nof_variables, facex, pointx)

      if (dimensiona .eq. 2) then
        call cons2prim(n, leftv, mp_pinfl, gammal)
      else
        call cons2prim(n, leftv, mp_pinfl, gammal)
      end if

      do iex = 1, nof_variables
        if ((leftv(iex) .lt. (utmin(iex) - nad_dg_el(iex))) .or. (leftv(iex) .gt. (utmax(iex) + nad_dg_el(iex)))) then
          ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
        end if
      end do

    case (11)       !mood indicator

      do iex = 1, nof_variables
        nad_dg_el(iex) = max(indicator_par1, (indicator_par2)*(utmax(iex) - utmin(iex)))
      end do

      do iex = 1, nof_variables
        if ((usol(iex,facex,pointx).lt.(utmin(iex)-nad_dg_el(iex))).or.(usol(iex,facex,pointx).gt.(utmax(iex)+nad_dg_el(iex))))then
          ielem(n, iconsidered)%condx = 1
        end if
      end do

    case (2)         !shu indicator

      do iex = 1, nof_variables

        if ((sumvars(iex)/maxvars(iex)) .gt. indicator_par1) then
          ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
        end if

      end do

    case (22)         !shock detector, indicator (only density and energy)

      do iex = 1, nof_variables

        if ((iex .eq. 1) .or. (iex .eq. img)) then
          if ((sumvars(iex)/maxvars(iex)) .gt. indicator_par1) then
            ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
          end if
        end if

      end do

    case (3)         !dmp

      do iex = 1, nof_variables

        if ((usol(iex, facex, pointx) .gt. (utmax(iex))) .or. (usol(iex, facex, pointx) .lt. (utmin(iex)))) then
          ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
        end if

      end do

    case (4)    !minmod
      do iex = 1, nof_variables

        if ((iex .eq. 1) .or. (iex .eq. img)) then
          if (abs(usol(iex, facex, pointx) - utemp(1, iex)) .gt. (indicator_par1*utemp(1, iex))) then
            ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
          end if
        end if
      end do

    case (5)    !all troubled

      ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1

    case (6)       !mood indicator only density & energy

      do iex = 1, nof_variables
      if ((iex .eq. 1) .or. (iex .eq. img)) then
        nad_dg_el(iex) = max(indicator_par1, (indicator_par2)*(utmax(iex) - utmin(iex)))
      end if
      end do
      leftv(1:nof_variables) = usol(1:nof_variables, facex, pointx)

      if (dimensiona .eq. 2) then
        call cons2prim(n, leftv, mp_pinfl, gammal)
      else
        call cons2prim(n, leftv, mp_pinfl, gammal)
      end if

      do iex = 1, nof_variables
      if ((iex .eq. 1) .or. (iex .eq. img)) then

        if ((leftv(iex) .lt. (utmin(iex) - nad_dg_el(iex))) .or. (leftv(iex) .gt. (utmax(iex) + nad_dg_el(iex)))) then
          ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
        end if
      end if
      end do

    case (9)       !mood indicator only density & energy

      if (ielem(n, iconsidered)%filtered .eq. 1) then
        ielem(n, iconsidered)%troubled = 1; ielem(n, iconsidered)%condition = 1
      end if

    end select

  end subroutine nad_dg
  subroutine apply_filter(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, j, k
    kmaxe = xmpielrank(n)
!$omp do
    do i = 1, kmaxe
    if (ielem(n, i)%filtered .eq. 1) then

      do j = 1, nof_variables
      do k = 1, idegfree
        rhs(i)%valdg(k + 1, j) = rhs(i)%valdg(k + 1, j)*modal_filter_weak(k)
      end do
      end do
    end if
    end do
!$omp end do
  end subroutine apply_filter
  subroutine apply_filter_dg(n)
    implicit none
    integer, intent(in)::n
    integer::i, j, k, kmaxe
    real::filtered_low
    real::filtered_high
    real::unfiltered, energy_ratio
    real, dimension(1:nof_variables)::en_f_strong, en_f_weak, en_unf, en_average
    integer::fil_i
    real::filx, xorder, ex1, ex2
    kmaxe = xmpielrank(n)
!$omp do
    do i = 1, kmaxe
      ielem(n, i)%filtered = 0
      en_unf(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      en_f_strong(1:nof_variables) = u_cs(i)%val(1, 1:nof_variables)
      en_f_weak(1:nof_variables) = u_cw(i)%val(1, 1:nof_variables)
      en_average(1:nof_variables) = u_c(i)%valdg(1, 1:nof_variables, 1)
      ex2 = (((en_unf(2) - en_f_weak(2))**2) + ((en_unf(3) - en_f_weak(3))**2) + ((en_unf(4) - en_f_weak(4))**2))
      ex1 = (((en_unf(2) - en_f_strong(2))**2) + ((en_unf(3) - en_f_strong(3))**2) + ((en_unf(4) - en_f_strong(4))**2))
      energy_ratio = (ex2 + 10e-32)/(ex1 + 10e-32)
      ielem(n, i)%er1dt = ((ielem(n, i)%er1 - ex1))/dt
      ielem(n, i)%er2dt = ((ielem(n, i)%er2 - ex2))/dt
      ielem(n, i)%er = energy_ratio
      ielem(n, i)%er1 = ex1
      ielem(n, i)%er2 = ex2
      if (ielem(n, i)%er .gt. (1.15)) then
        ielem(n, i)%filtered = 1
      end if
    end do
!$omp end do
  end subroutine apply_filter_dg
  subroutine apply_filter1(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, j, k
    kmaxe = xmpielrank(n)
!$omp do
    do i = 1, kmaxe
    do j = 1, nof_variables
    do k = 1, idegfree
      u_c(i)%valdg(1, j, k + 1) = u_c(i)%valdg(1, j, k + 1)*modal_filter(k)
    end do
    end do
    end do
!$omp end do
  end subroutine apply_filter1
  subroutine apply_filter2(n)
    implicit none
    integer, intent(in)::n
    integer::i, kmaxe, j, k
    real, allocatable, dimension(:)::grad1al
    kmaxe = xmpielrank(n)
    allocate(grad1al(1:idegfree))
!$omp do
    do i = 1, kmaxe
    do j = 1, nof_variables
      grad1al(1:idegfree) = 0.0d0
      do k = 1, idegfree
        grad1al(k) = rhs(i)%valdg(k + 1, j)*modal_filter(k)
      end do
      rhs(i)%valdg(2:idegfree + 1, j) = grad1al(1:idegfree)
    end do
    end do
!$omp end do
    deallocate(grad1al)
  end subroutine apply_filter2
  subroutine filter(n)
    implicit none
    integer, intent(in)::n
    integer::fil_i, i, j
    real::filx, xorder
    real::rfil_alpha, rfil_nc, rfil_s, rfil_i
    real, dimension(1:9)::filter2
    real, dimension(1:200)::dgfr
    integer, dimension(0:9)::filt2
    if (filter_type .eq. 1) then

      do fil_i = 1, iorder
        if (fil_i .le. fil_nc) then
          filter2(fil_i) = 1.0d0
        end if
        if (fil_i .gt. fil_nc) then
          rfil_alpha = fil_alpha
          rfil_nc = fil_nc
          rfil_s = fil_s
          xorder = iorder
          rfil_i = fil_i
          filx = -rfil_alpha*(((rfil_i - rfil_nc)/(xorder - rfil_nc))**rfil_s)
          filter2(fil_i) = exp(filx)
        end if

        filt2(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
      end do
      filt2(0) = 0

      do i = 1, iorder
        do j = filt2(i - 1) + 1, filt2(i)
          dgfr(j) = filter2(i)
        end do
      end do

      do fil_i = 1, idegfree
        modal_filter(fil_i) = dgfr(fil_i)**(1.0d0/(1.0/dt))!
        if ((it .le. 2) .and. (n .eq. 0)) then
          write (200 + n, *) fil_i, modal_filter(fil_i)
        end if

      end do

    end if

    if (filter_type .eq. 2) then

      do fil_i = 1, iorder
        if (fil_i .lt. iorder) then
          filter2(fil_i) = 1.0d0
        end if
        if (fil_i .eq. iorder) then
          filter2(fil_i) = 0.0d0
        end if

        filt2(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
      end do
      filt2(0) = 0

      do i = 1, iorder
        do j = filt2(i - 1) + 1, filt2(i)
          dgfr(j) = filter2(i)
        end do
      end do

      do fil_i = 1, idegfree
        modal_filter(fil_i) = dgfr(fil_i)**(1.0d0/(1.0/dt))!
        if ((it .le. 2) .and. (n .eq. 0)) then
          write (200 + n, *) fil_i, modal_filter(fil_i)
        end if

      end do

    end if

    if (filter_type .eq. 3) then

      do fil_i = 1, iorder
        if (fil_i .lt. iorder) then
          filter2(fil_i) = 1.0d0
        else
          filter2(fil_i) = 0.0d0
        end if
        filt2(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
      end do
      filt2(0) = 0

      do i = 1, iorder
        do j = filt2(i - 1) + 1, filt2(i)
          dgfr(j) = filter2(i)
        end do
      end do

      do fil_i = 1, idegfree
        modal_filter_weak(fil_i) = dgfr(fil_i)!**(1.0d0/(1.0/dt))!

      end do

      do fil_i = 1, iorder
        if (fil_i .lt. iorder - 1) then
          filter2(fil_i) = 1.0d0
        else
          filter2(fil_i) = 0.0d0
        end if

        filt2(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
      end do
      filt2(0) = 0

      do i = 1, iorder
        do j = filt2(i - 1) + 1, filt2(i)
          dgfr(j) = filter2(i)
        end do
      end do

      do fil_i = 1, idegfree
        modal_filter_strong(fil_i) = dgfr(fil_i)!**(1.0d0/(1.0/dt))!
      end do

      do fil_i = 1, idegfree
        if ((it .le. 2) .and. (n .eq. 0)) then
          write (200 + n, *) fil_i, modal_filter_weak(fil_i), modal_filter_strong(fil_i)
        end if

      end do

    end if

  end subroutine filter

  subroutine adda_filter(n, iconsidered)
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i, j, k
    real::filtered_low
    real::filtered_high
    real::unfiltered, energy_ratio
    real, dimension(1:nof_variables)::en_f_strong, en_f_weak, en_unf
    integer::fil_i, countdof, icompwrt, ngp
    real::filx, xorder, ex1, ex2
    real::rfil_alpha, rfil_nc, rfil_s, rfil_i
    real, allocatable, dimension(:)::filter2, filter3
    real, allocatable, dimension(:)::dgfr, dgfr3
    integer, allocatable, dimension(:)::filt2, filt3
    real::ax, ay, az, mp_pinfl, gammal
    real, dimension(1:nof_variables)::leftv
    real, allocatable, dimension(:, :)::consmatrix, gradssl, ressolution

    countdof = ((iorder + 1)*(iorder + 2)*(iorder + 3))/6

    allocate(filter2(1:9), filter3(1:9), dgfr(1:countdof), dgfr3(1:countdof), filt2(0:9), filt3(0:9))
    allocate(consmatrix(1, 1:idegfree))
    allocate(gradssl(1:idegfree, 1:nof_variables))
    allocate(ressolution(1:6*numberofpoints2, 1:nof_variables))

    do fil_i = 1, iorder
      if (iorder .eq. 2) then
        filter2(fil_i) = 0.0d0
      else

        if (fil_i .lt. 2) then  !if (fil_i.le.adda_1)then
          filter2(fil_i) = 1.0d0
        else
          filter2(fil_i) = 0.0d0

        end if
      end if

      filt2(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
    end do
    filt2(0) = 0

    do i = 1, iorder
      do j = filt2(i - 1) + 1, filt2(i)
        dgfr(j) = filter2(i)
      end do
    end do

    do fil_i = 1, iorder
      if (fil_i .le. 2) then!                                        if (fil_i.le.adda_2)then
        filter3(fil_i) = 1.0d0
      else
        filter3(fil_i) = 0.0d0
      end if

      filt3(fil_i) = (((fil_i + 1)*(fil_i + 2)*(fil_i + 3))/6) - 1
    end do
    filt3(0) = 0

    do i = 1, iorder
      do j = filt3(i - 1) + 1, filt3(i)
        dgfr3(j) = filter3(i)
      end do
    end do

    do fil_i = 1, idegfree
      adda_filter_strong(fil_i) = dgfr(fil_i)
      adda_filter_weak(fil_i) = dgfr3(fil_i)
      if ((it .eq. 0) .and. (iconsidered .eq. 1) .and. (n .eq. 0)) then
        write (300 + n, *) fil_i, adda_filter_weak(fil_i), adda_filter_strong(fil_i)

      end if
    end do

    i = iconsidered

    if (dg .ne. 1) then

      if (adda_type .eq. 1) then

        ax = 0.0d0; ay = 0.0d0; az = 0.0d0

        icompwrt = 0

        consmatrix(1, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)

 gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = ilocal_recon5(iconsidered)%gradients(1, 1:ielem(n, i)%idegfree, 1:nof_variables)

       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))

        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)

        en_unf(1:nof_variables) = leftv(1:nof_variables)

        do k = 1, nof_variables
                                gradssl(1:ielem(n,i)%idegfree,k)=ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,k)*adda_filter_strong(1:ielem(n,i)%idegfree)
        end do

       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))

        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)

        en_f_strong(1:nof_variables) = leftv(1:nof_variables)

        do k = 1, nof_variables
                                gradssl(1:ielem(n,i)%idegfree,k)=ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,k)*adda_filter_weak(1:ielem(n,i)%idegfree)
        end do

       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))

        leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)

        en_f_weak(1:nof_variables) = leftv(1:nof_variables)

        ex2 = (((en_unf(2) - en_f_weak(2))**2) + ((en_unf(3) - en_f_weak(3))**2) + ((en_unf(4) - en_f_weak(4))**2))
        ex1 = (((en_unf(2) - en_f_strong(2))**2) + ((en_unf(3) - en_f_strong(3))**2) + ((en_unf(4) - en_f_strong(4))**2))

        energy_ratio = (ex2 + 10e-32)/(ex1 + 10e-32)

      end if

      if (adda_type .eq. 2) then
        en_unf(1:nof_variables) = zero
        en_f_strong(1:nof_variables) = zero
        en_f_weak(1:nof_variables) = zero

        do ngp = 1, ielem(n, iconsidered)%itotalpoints

          ax = qp_array(iconsidered)%x(ngp)
          ay = qp_array(iconsidered)%y(ngp)
          az = qp_array(iconsidered)%z(ngp)

          icompwrt = 0

          consmatrix(1, 1:ielem(n, i)%idegfree) = basis_rec(n, ax, ay, az, ielem(n, i)%iorder, i, ielem(n, i)%idegfree, icompwrt)
 gradssl(1:ielem(n, i)%idegfree, 1:nof_variables) = ilocal_recon5(iconsidered)%gradients(1, 1:ielem(n, i)%idegfree, 1:nof_variables)
       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
          leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)
          en_unf(1:nof_variables) = en_unf(1:nof_variables) + leftv(1:nof_variables)*qp_array(iconsidered)%qp_weight(ngp)
          do k = 1, nof_variables
            gradssl(1:ielem(n,i)%idegfree,k)=ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,k)*adda_filter_strong(1:ielem(n,i)%idegfree)
          end do
       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
          leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)
          en_f_strong(1:nof_variables) = en_f_strong(1:nof_variables) + leftv(1:nof_variables)*qp_array(iconsidered)%qp_weight(ngp)
          do k = 1, nof_variables
            gradssl(1:ielem(n,i)%idegfree,k)=ilocal_recon5(iconsidered)%gradients(1,1:ielem(n,i)%idegfree,k)*adda_filter_weak(1:ielem(n,i)%idegfree)
          end do
       ressolution(1:1,1:nof_variables)=matmul(consmatrix(1:1,1:ielem(n,i)%idegfree),gradssl(1:ielem(n,i)%idegfree,1:nof_variables))
          leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + ressolution(1, 1:nof_variables)
          en_f_weak(1:nof_variables) = en_f_weak(1:nof_variables) + leftv(1:nof_variables)*qp_array(iconsidered)%qp_weight(ngp)
        end do
        ex2 = (((en_unf(2) - en_f_weak(2))**2) + ((en_unf(3) - en_f_weak(3))**2) + ((en_unf(4) - en_f_weak(4))**2))
        ex1 = (((en_unf(2) - en_f_strong(2))**2) + ((en_unf(3) - en_f_strong(3))**2) + ((en_unf(4) - en_f_strong(4))**2))
        energy_ratio = (ex2 + 10e-32)/(ex1 + 10e-32)
      end if
      !energy_ratio=abs(energy_ratio-1.0d0)/(ielem(n,i)%totvolume)
      ielem(n, i)%er1dt = (ielem(n, i)%er1 - ex1)/dt
      ielem(n, i)%er2dt = (ielem(n, i)%er2 - ex2)/dt
      !ielem(n,i)%er=(ielem(n,i)%erx-energy_ratio)/dt
      ielem(n, i)%er = energy_ratio
      ielem(n, i)%er1 = ex1
      ielem(n, i)%er2 = ex2
      ielem(n, i)%er1er2 = abs(ielem(n, i)%er1dt)/ielem(n, i)%er2dt
      call apply_adda_filter(n, iconsidered)
    end if
    deallocate(filter2, filter3, dgfr, dgfr3, filt2, filt3)
    deallocate(consmatrix)
    deallocate(gradssl)
    deallocate(ressolution)
  end subroutine adda_filter
  subroutine apply_adda_filter(n, iconsidered)
    implicit none
    integer, intent(in)::n
    integer::i, iconsidered
    real::lwcx1
    i = iconsidered
    lwcx1 = lwci1
    if (rungekutta .eq. 11) then
    if (iscoun .eq. 1) then
      lwcx1 = lwci1
      if (ielem(n, i)%er .gt. 1.2) then
        lwcx1 = 10!increase dissipation
      end if
      if (ielem(n, i)%er .le. 0.95) then
        lwcx1 = 1000
      end if
      ielem(n, i)%lwcx2 = lwcx1
    else
      lwcx1 = ielem(n, i)%lwcx2
    end if
    else
    if (ielem(n, i)%er .gt. 1.2) then
      lwcx1 = 10!increase dissipation
    end if
    if (ielem(n, i)%er .le. 0.95) then
      lwcx1 = 1000
      !lwcx1=100**(12-(4*ielem(n,i)%er**0.8))
    end if
    ielem(n, i)%lwcx2 = lwcx1
    end if
    if (ielem(n, i)%full .eq. 0) then
      ielem(n, i)%lwcx2 = -10
    end if
    ielem(n, i)%linc = lwcx1
  end subroutine apply_adda_filter

  subroutine fix_dissipation(n)
    implicit none
    real::check1
    real::check
    integer::i, j, k, l, kmaxe
    integer, intent(in)::n

    kmaxe = xmpielrank(n)

!$omp do
    do i = 1, kmaxe

      if (ielem(n, i)%full .eq. 1) then
        !1)reduce dissipation
        if (ielem(n, i)%lwcx2 .gt. 10) then
          if (ielem(n, i)%wcx(1) .ge. 0.999) then
            ielem(n, i)%diss = max(ielem(n, i)%diss - 0.1, 0.5d0)        !reduce dissipation even more
          else
            ielem(n, i)%diss = 1.0d0                                                        !increase dissipation if shock
          end if
        else
          !2) increase dissipation
          if (ielem(n, i)%wcx(1) .ge. 0.999) then
            ielem(n, i)%diss = min(ielem(n, i)%diss + 0.1, 1.0d0)        !increase dissipation even more
          else
            ielem(n, i)%diss = 1.0d0                                                        !increase dissipation if shock
          end if
        end if
      end if

    end do
  end subroutine fix_dissipation
  subroutine fix_dissipation2(n)
    implicit none
    real::check1
    real::check
    integer::i, j, k, l, kmaxe, iconsidered
    integer, intent(in)::n
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      iconsidered = i
      call find_bounds_diss(iconsidered)
    end do
  end subroutine fix_dissipation2
  subroutine find_bounds_diss(iconsidered)
    implicit none
    integer::i, l, j, k, kmaxe, iqp, ngp, iex, ik, iq
    integer, intent(in)::iconsidered
    real, dimension(1:6, 1)::utemp
    i = iconsidered
    utemp(:, :) = 1.0d0
    ielem(n, i)%facediss(:) = 1.0d0
    if (ielem(n, i)%interior .eq. 0) then
      do l = 1, ielem(n, i)%ifca
        utemp(l, 1) = ielem(n, ielem(n, i)%ineigh(l))%diss
      end do
    end if
    if (ielem(n, i)%interior .eq. 1) then
      do l = 1, ielem(n, i)%ifca
        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              utemp(l, 1) = ielem(n, ielem(n, i)%ineigh(l))%diss
            else
              !not periodic ones in my cpu
            end if
          else
            utemp(l, 1) = ielem(n, ielem(n, i)%ineigh(l))%diss
          end if
        else        !in other cpus they can only be periodic or mpi neighbours
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
              utemp(l, 1) = iexsolhird(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                            (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1)
            end if
          else
            utemp(l, 1) = iexsolhird(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1)
          end if
        end if
      end do
    end if
    if (ielem(n, i)%interior .eq. 0) then
      do l = 1, ielem(n, i)%ifca
        ielem(n, i)%facediss(l) = max(utemp(l, 1), ielem(n, i)%diss)
      end do
    end if
    if (ielem(n, i)%interior .eq. 1) then
      do l = 1, ielem(n, i)%ifca
        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              ielem(n, i)%facediss(l) = max(utemp(l, 1), ielem(n, i)%diss)
            else
              !not periodic ones in my cpu
            end if
          else
            ielem(n, i)%facediss(l) = max(utemp(l, 1), ielem(n, i)%diss)
          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              ielem(n, i)%facediss(l) = max(utemp(l, 1), ielem(n, i)%diss)
            end if
          else

            ielem(n, i)%facediss(l) = max(utemp(l, 1), ielem(n, i)%diss)
          end if

        end if

      end do
    end if

  end subroutine find_bounds_diss

  subroutine vfbp_limiter

    implicit none
    real::vf_qpsol, vf_avsol, lthresh, hthresh, scaling, scaling1, scaling2, pd1_qpsol, pd1_avsol, pd2_qpsol, pd2_avsol
    integer::i, ii, icd, k, number_of_nei, idummy, l, nnd, ngp, iqp, i_elem, i_face, img1, img2
    real, dimension(1:numberofpoints2)::scaling_f, scaling_f1, scaling_f2
    integer::facex, pointx, iconsidered, number_of_dog

    !iqp=qp_line_n

    lthresh = 1.0e-16
    hthresh = 1.0d0 - lthresh

!$omp do
    do i = 1, xmpielrank(n)

      do l = 1, ielem(n, i)%ifca

        if (dimensiona .eq. 2) then

          iqp = qp_line_n
          img1 = 5
          img2 = 6
        else
          img1 = 6
          img2 = 7
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad
          else
            iqp = qp_triangle

          end if
        end if

        !compute scaling factors at each face and each qp

        do ngp = 1, iqp! qp_line_n

          facex = l
          pointx = ngp
          iconsidered = i
          number_of_dog = ielem(n, i)%idegfree

          vf_qpsol = ilocal_recon3(iconsidered)%uleft_dg(nof_variables, facex, pointx)
          vf_avsol = u_c(i)%valdg(1, nof_variables, 1)

          pd1_qpsol = ilocal_recon3(iconsidered)%uleft_dg(img1, facex, pointx)
          pd1_avsol = u_c(i)%valdg(1, img1, 1)

          pd2_qpsol = ilocal_recon3(iconsidered)%uleft_dg(img2, facex, pointx)
          pd2_avsol = u_c(i)%valdg(1, img2, 1)

          !volume fraction
          if (vf_qpsol .lt. lthresh) then

            scaling_f(ngp) = (vf_avsol - lthresh)/(vf_avsol - vf_qpsol)

          else if ((vf_qpsol .gt. lthresh) .and. (vf_qpsol .lt. hthresh)) then

            scaling_f(ngp) = 1.0d0

          else if (vf_qpsol .gt. hthresh) then

            scaling_f(ngp) = (hthresh - vf_avsol)/(vf_qpsol - vf_avsol)

          end if

          !partial densities
          if (pd1_qpsol .lt. lthresh) then

            scaling_f1(ngp) = (pd1_avsol - lthresh)/(pd1_avsol - pd1_qpsol)
          else

            scaling_f1(ngp) = 1.0d0

          end if

          if (pd2_qpsol .lt. lthresh) then

            scaling_f2(ngp) = (pd2_avsol - lthresh)/(pd2_avsol - pd2_qpsol)

          else

            scaling_f2(ngp) = 1.0d0

          end if

        end do

        scaling = minval(scaling_f(:))
        scaling1 = minval(scaling_f1(:))
        scaling2 = minval(scaling_f2(:))

        do ngp = 1, iqp! qp_line_n

          facex = l
          pointx = ngp
          iconsidered = i
          number_of_dog = ielem(n, i)%idegfree

          vf_qpsol = ilocal_recon3(iconsidered)%uleft_dg(nof_variables, facex, pointx)
          vf_avsol = u_c(i)%valdg(1, nof_variables, 1)

          pd1_qpsol = ilocal_recon3(iconsidered)%uleft_dg(img1, facex, pointx)
          pd1_avsol = u_c(i)%valdg(1, img1, 1)

          pd2_qpsol = ilocal_recon3(iconsidered)%uleft_dg(img2, facex, pointx)
          pd2_avsol = u_c(i)%valdg(1, img2, 1)

          !volume fraction scaling at qps
          if ((scaling .gt. 0.0d0) .and. (scaling .lt. 1.0d0)) then
            ilocal_recon3(iconsidered)%uleft_dg(nof_variables, facex, pointx) = vf_avsol + scaling*(vf_qpsol - vf_avsol)

          end if

          !partial densities scaling at qps

          if ((scaling1 .gt. 0.0d0) .and. (scaling1 .lt. 1.0d0)) then
            ilocal_recon3(iconsidered)%uleft_dg(img1, facex, pointx) = pd1_avsol + scaling1*(pd1_qpsol - pd1_avsol)

          end if

          if ((scaling2 .gt. 0.0d0) .and. (scaling2 .lt. 1.0d0)) then
            ilocal_recon3(iconsidered)%uleft_dg(img2, facex, pointx) = pd2_avsol + scaling2*(pd2_qpsol - pd2_avsol)

          end if

        end do

      end do

    end do
!$omp end do

  end subroutine vfbp_limiter
end module recon
