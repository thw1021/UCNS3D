module moodr
  use mpiinfo
  use translate
  use declaration
  use memory
  use communications
  use io
  use partition
  use library
  use transform
  use fluxes
  use initialisation
  use boundary
  use recon
  use local
  use profile
  use flow_operations
  use gradients
  use basis
  use prestore
  use riemann
  use source
  use implicit_time
  use implicit_fluxes
  use omp_lib
  implicit none
contains
  subroutine pad_nad(n)
    implicit none
    integer, intent(in) :: n
    integer :: i, l, ngp, iqp, iex, kmaxe, k, ii, iconsidered
    integer :: reduce1
    real, dimension(nof_variables) :: nad_delta1
    integer :: pad_true, nad_true
    real :: relax_mood1, relax_mood2
    real, dimension(1:nof_variables) :: leftv
    real :: mp_pinfl, gammal
    real, dimension(1:nof_variables) :: rightv
    real :: mp_pinfr, gammar
    real, allocatable, dimension(:, :) :: utemp
    real, dimension(1:nof_variables) :: utmin, utmax

    allocate(utemp(imaxdegfree + 1, 1:nof_variables + turbulenceequations + passivescalar))
    kmaxe = xmpielrank(n)

    if (cascade .eq. 1) then
      relax_mood1 = mood_var1
      relax_mood2 = mood_var2
    end if

    if (cascade .eq. 2) then
      relax_mood1 = mood_var3
      relax_mood2 = mood_var4
    end if

    if (itestcase .ge. 3) then
      do ii = 1, nof_interior
        i = el_int(ii)
        iconsidered = i
        ielem(n, i)%mood = 0
        reduce1 = 0
        pad_true = 0
        nad_true = 0
        leftv(1:nof_variables) = u_c(i)%val(4, 1:nof_variables)

        if (dimensiona .eq. 3) then
          call cons2prim(n, leftv, mp_pinfl, gammal)
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            pad_true = 1
          end if
          if ((leftv(5) .le. zero) .or. (leftv(5) .ne. leftv(5))) then
            pad_true = 1
          end if
        else
          call cons2prim(n, leftv, mp_pinfl, gammal)
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            pad_true = 1
          end if
          if ((leftv(4) .le. zero) .or. (leftv(4) .ne. leftv(4))) then
            pad_true = 1
          end if
        end if

        if (pad_true .eq. 0) then
          utemp = zero
          k = 0
          utemp(1, 1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
          k = 1
          do l = 1, ielem(n, i)%ifca
            k = k + 1
            utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
          end do

          do iex = 1, nof_variables
            utmin(iex) = minval(utemp(1:k, iex))
            utmax(iex) = maxval(utemp(1:k, iex))
          end do

          do iex = 1, nof_variables
            nad_delta1(iex) = max(relax_mood1, (relax_mood2)*(utmax(iex) - utmin(iex)))
          end do

          nad_true = 0
          if (mood_mode .gt. 0) then
            do iex = 1, nof_variables
              if ((iex .eq. 1) .or. (iex .eq. nof_variables)) then
                if ((u_c(i)%val(4, iex) .lt. (utmin(iex) - nad_delta1(iex))) .or. &
                    (u_c(i)%val(4, iex) .gt. (utmax(iex) + nad_delta1(iex)))) then
                  nad_true = 1
                end if
              end if
            end do
          else
            do iex = 1, nof_variables
              if ((iex .eq. 1) .or. (iex .eq. nof_variables)) then
                if ((u_c(i)%val(4, iex) .lt. utmin(iex)) .or. &
                    (u_c(i)%val(4, iex) .gt. utmax(iex))) then
                  nad_true = 1
                end if
              end if
            end do
          end if
        end if

        ielem(n, i)%mood = nad_true + pad_true
      end do

      do ii = 1, nof_bounded
        i = el_bnd(ii)
        iconsidered = i
        reduce1 = 0
        pad_true = 0
        nad_true = 0
        ielem(n, i)%mood = 0
        leftv(1:nof_variables) = u_c(i)%val(4, 1:nof_variables)

        if (dimensiona .eq. 3) then
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            pad_true = 1
          end if
          if ((leftv(5) .le. zero) .or. (leftv(5) .ne. leftv(5))) then
            pad_true = 1
          end if
        else
          call cons2prim(n, leftv, mp_pinfl, gammal)
          if ((leftv(1) .le. zero) .or. (leftv(1) .ne. leftv(1))) then
            pad_true = 1
          end if
          if ((leftv(4) .le. zero) .or. (leftv(4) .ne. leftv(4))) then
            pad_true = 1
          end if
        end if

        if (pad_true .eq. 0) then
          utemp = zero
          k = 0
          utemp(1, 1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
          k = 1
          do l = 1, ielem(n, i)%ifca
            if (ielem(n, i)%ineighb(l) .eq. n) then
              if (ielem(n, i)%ibounds(l) .gt. 0) then
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then
                  k = k + 1
                  utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
                end if
              else
                k = k + 1
                utemp(k, 1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
              end if
            else
              if (ielem(n, i)%ibounds(l) .gt. 0) then
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then
                  k = k + 1
                  utemp(k, 1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                end if
              else
                k = k + 1
                utemp(k, 1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                            (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              end if
            end if
          end do

          do iex = 1, nof_variables
            utmin(iex) = minval(utemp(1:k, iex))
            utmax(iex) = maxval(utemp(1:k, iex))
          end do

          do iex = 1, nof_variables
            nad_delta1(iex) = max(relax_mood1, (relax_mood2)*(utmax(iex) - utmin(iex)))
          end do

          nad_true = 0
          if (mood_mode .gt. 0) then
            do iex = 1, nof_variables
              if ((iex .eq. 1) .or. (iex .eq. nof_variables)) then
                if ((u_c(i)%val(4, iex) .lt. (utmin(iex) - nad_delta1(iex))) .or. &
                    (u_c(i)%val(4, iex) .gt. (utmax(iex) + nad_delta1(iex)))) then
                  nad_true = 1
                end if
              end if
            end do
          else
            do iex = 1, nof_variables
              if ((iex .eq. 1) .or. (iex .eq. nof_variables)) then
                if ((u_c(i)%val(4, iex) .lt. utmin(iex)) .or. &
                    (u_c(i)%val(4, iex) .gt. utmax(iex))) then
                  nad_true = 1
                end if
              end if
            end do
          end if
        end if

        ielem(n, i)%mood = nad_true + pad_true
      end do
    end if

    deallocate(utemp)
  end subroutine pad_nad

  subroutine mood_operator_2(n)
    implicit none
    integer :: i, kmaxe
    integer, intent(in) :: n
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      ielem(n, i)%recalc = 0
      ielem(n, i)%mood = 0
    end do

    cascade = 1
    call pad_nad(n)
    call exhboundhigher_mood(n)
    call fix_list(n)
    call muscl(n)
    call exhboundhigher(n)

    if (dimensiona .eq. 2) then
      call calculate_fluxeshi_convective2d_mood(n)
    else
      call calculate_fluxeshi_convective_mood(n)
    end if

    do i = 1, kmaxe
      if (ielem(n, i)%mood .eq. 0) then
        ielem(n, i)%mood_o = iorder + 1  ! Target polynomial/scheme admissible
      else
        ielem(n, i)%mood_o = 2         ! Second order MUSCL admissible
      end if
    end do
  end subroutine mood_operator_2

  subroutine mood_operator_1(n)
    implicit none
    integer :: i, kmaxe
    integer, intent(in) :: n
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      ielem(n, i)%recalc = 0
      ielem(n, i)%mood = 0
    end do

    cascade = 2
    call pad_nad(n)
    call exhboundhigher_mood(n)
    call fix_list(n)

    if (dimensiona .eq. 2) then
      call calculate_fluxeshi_convective2d_mood(n)
    else
      call calculate_fluxeshi_convective_mood(n)
    end if

    do i = 1, kmaxe
      if (ielem(n, i)%mood .eq. 1) then
        ielem(n, i)%mood_o = 1     ! Only first order solution admissible
      end if
    end do
  end subroutine mood_operator_1

  subroutine fix_list(n)
    implicit none
    integer, intent(in) :: n
    real :: godflux2, sum_detect
    integer :: i, l, ngp, kmaxe, iqp
    real, dimension(numberofpoints2) :: weights_temp
    real, dimension(1) :: cright
    kmaxe = xmpielrank(n)

    do i = 1, kmaxe
      ielem(n, i)%recalc = 0
      if (ielem(n, i)%mood .eq. 1) then
        ielem(n, i)%recalc = 1
      else
        if (ielem(n, i)%interior .eq. 0) then
          cright(1) = zero
          do l = 1, ielem(n, i)%ifca
            cright(1) = ielem(n, (ielem(n, i)%ineigh(l)))%mood
            if (cright(1) .gt. 0.5) then
              ielem(n, i)%recalc = 1
            end if
          end do
        end if

        if (ielem(n, i)%interior .eq. 1) then
          cright(1) = zero
          do l = 1, ielem(n, i)%ifca
            if (ielem(n, i)%ineighb(l) .eq. n) then
              if (ielem(n, i)%ibounds(l) .gt. 0) then
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then
                  cright(1) = ielem(n, (ielem(n, i)%ineigh(l)))%mood
                end if
              else
                cright(1) = ielem(n, (ielem(n, i)%ineigh(l)))%mood
              end if
            else
              if (ielem(n, i)%ibounds(l) .gt. 0) then
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then
                  cright(1) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(1), 1)
                end if
              else
                cright(1) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(1), 1)
              end if
            end if

            if (cright(1) .gt. 0.5) then
              ielem(n, i)%recalc = 1
            end if
          end do
        end if
      end if
    end do
  end subroutine fix_list
end module moodr