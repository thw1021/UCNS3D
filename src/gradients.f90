module gradients
  use library
  use flow_operations
  implicit none

contains

  subroutine allgrads_inner(n, iconsidered)
!> @brief
!> this subroutine calls the gradient approximation subroutines for every interior cell
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i
    integer::number_of_dog, number_of_nei, imax
    i = iconsidered
    number_of_dog = ielem(n, i)%idegfree
    number_of_nei = ielem(n, i)%inumneighbours
! imax=number_of_nei-1

    if (dg .eq. 1) then
      call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
    else
      select case (ielem(n, i)%ggs)
      case (0)        !least squares everything
        select case (itestcase)
        case (1, 2, 3)
        call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
          if (initcond .eq. 95) then
            call compute_gradients_inner_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
          end if

        case (4)
          call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
          call compute_gradients_inner_mean_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
          if (turbulence .eq. 1) then
            call compute_gradients_inner_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
            call compute_gradients_turb_lsq(n, iconsidered, number_of_dog, number_of_nei)
            call compute_gradients_turb_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
          end if

        end select

      case (1) !green gauss everything
        select case (itestcase)
        case (1, 2, 3)
          call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
          if (initcond .eq. 95) then
            call compute_gradients_inner_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
          end if

        case (4)
          call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
          call compute_gradients_inner_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)

          if (turbulence .eq. 1) then
            call compute_gradients_turb_lsq(n, iconsidered, number_of_dog, number_of_nei)
            call compute_gradients_inner_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
          end if

        end select

      end select

    end if

  end subroutine allgrads_inner

  subroutine allgrads_mix(n, iconsidered)
!> @brief
!> this subroutine calls the gradient approximation subroutines for every non-interior cell
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i
    integer::number_of_dog, number_of_nei, imax
    i = iconsidered
    number_of_dog = ielem(n, i)%idegfree
    number_of_nei = ielem(n, i)%inumneighbours

    if (dg .eq. 1) then
      call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
    else
      if (fastest .ne. 1) then !least squares
        select case (ielem(n, i)%ggs)

        case (0)        !least squares everything
          select case (itestcase)

          case (1, 2, 3)

            call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
            if (initcond .eq. 95) then
              call compute_gradients_mix_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
            end if

          case (4)
            call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
            if (ielem(n, i)%walls .ne. 1) then
              call compute_gradients_inner_mean_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
            else
              call compute_gradients_wall_mean_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
            end if
            if (turbulence .eq. 1) then
              call compute_gradients_mix_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
              call compute_gradients_turb_lsq(n, iconsidered, number_of_dog, number_of_nei)
              if (ielem(n, i)%walls .ne. 1) then
                call compute_gradients_turb_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
              else
                call compute_gradients_wall_turb_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)
              end if
            end if

          end select

        case (1)
          select case (itestcase)

          case (1, 2, 3)

            call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
            if (initcond .eq. 95) then
              call compute_gradients_mix_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
            end if

          case (4)
            call compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)
            call compute_gradients_mix_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
            if (turbulence .eq. 1) then
              call compute_gradients_turb_lsq(n, iconsidered, number_of_dog, number_of_nei)
              call compute_gradients_mix_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)
            end if

          end select
        end select
      end if
    end if

  end subroutine allgrads_mix

  subroutine allgrads_mix_av(n, iconsidered)
!> @brief
!> this subroutine calls the average gradient approximation subroutines for every non interior cell
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i
    integer::number_of_dog, number_of_nei, imax
    i = iconsidered

    number_of_dog = ielem(n, i)%idegfree
    number_of_nei = ielem(n, i)%inumneighbours

    call compute_gradients_mix_mean_ggs_viscous_av(n, iconsidered, number_of_dog, number_of_nei)

  end subroutine allgrads_mix_av

  subroutine allgrads_inner_av(n, iconsidered)
!> @brief
!> this subroutine calls the average gradient approximation subroutines for every interior cell
    implicit none
    integer, intent(in)::n, iconsidered
    integer::i
    integer::number_of_dog, number_of_nei, imax
    i = iconsidered
    number_of_dog = ielem(n, i)%idegfree
    number_of_nei = ielem(n, i)%inumneighbours

    call compute_gradients_inner_mean_ggs_viscous_av(n, iconsidered, number_of_dog, number_of_nei)

  end subroutine allgrads_inner_av

  subroutine compute_gradients_mean_lsq(n, iconsidered, number_of_dog, number_of_nei)!check all
!> @brief
!> this subroutine computes the gradients of the conserved variables of each cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(nof_variables)::sols1
    real, dimension(nof_variables, typesten)::sols2
    real, allocatable, dimension(:, :, :)::matrix_1, matrix_2
    real, allocatable, dimension(:, :, :)::sol_m
    real, dimension(1:nof_variables)::leftv, rightv
    real::mp_pinfl, gammal
    integer::i, var2, iq, ll, imax
    real::tempxx
    i = iconsidered

    allocate(matrix_1(number_of_nei - 1, nof_variables, ielem(n, iconsidered)%admis))
    allocate(matrix_2(number_of_dog, nof_variables, ielem(n, iconsidered)%admis))
    allocate(sol_m(number_of_dog, nof_variables, ielem(n, iconsidered)%admis))

    imax = number_of_nei - 1

    i = iconsidered
    sols1 = zero; sols2 = zero; matrix_2 = zero; sol_m = zero; matrix_1 = zero; matrix_2 = zero

    sols1(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
    if (wenwrt .eq. 3) then
      leftv(1:nof_variables) = sols1(1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols1(1:nof_variables) = leftv(1:nof_variables)
    end if

    if (ilocal_recon3(i)%local .eq. 1) then

      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          do iq = 1, imax
            sols2(1:nof_variables, ll) = u_c(ilocal_recon3(i)%ihexl(ll, iq + 1))%val(1, 1:nof_variables)

            if (wenwrt .eq. 3) then
              leftv(1:nof_variables) = sols2(1:nof_variables, ll)
              call cons2prim(n, leftv, mp_pinfl, gammal)
              sols2(1:nof_variables, ll) = leftv(1:nof_variables)
            end if

            if (per_rot .eq. 1) then
              if (ilocal_recon3(i)%periodicflag(ll, iq + 1) .eq. 2) then
                tempxx = sols2(2, ll)
                sols2(2, ll) = tempxx*cos(angle_per) - sols2(3, ll)*sin(angle_per)
                sols2(3, ll) = tempxx*sin(angle_per) + sols2(3, ll)*cos(angle_per)
              end if
              if (ilocal_recon3(i)%periodicflag(ll, iq + 1) .eq. 1) then
                tempxx = sols2(2, ll)
                sols2(2, ll) = tempxx*cos(-angle_per) - sols2(3, ll)*sin(-angle_per)
                sols2(3, ll) = tempxx*sin(-angle_per) + sols2(3, ll)*cos(-angle_per)
              end if
            end if
            matrix_1(iq, 1:nof_variables, ll) = (sols2(1:nof_variables, ll) - sols1(1:nof_variables))
          end do

        else

          do iq = 1, numneighbours2 - 1
            sols2(1:nof_variables, ll) = u_c(ilocal_recon3(i)%ihexlc(ll, iq + 1))%val(1, 1:nof_variables)
            if (wenwrt .eq. 3) then
              leftv(1:nof_variables) = sols2(1:nof_variables, ll)
              call cons2prim(n, leftv, mp_pinfl, gammal)
              sols2(1:nof_variables, ll) = leftv(1:nof_variables)
            end if
            matrix_1(iq, 1:nof_variables, ll) = (sols2(1:nof_variables, ll) - sols1(1:nof_variables))
          end do
        end if
      end do
      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then

!          call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables,ll),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables,ll),ielem(n,i)%idegfree)
          sol_m(1:ielem(n,i)%idegfree,1:nof_variables,ll)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables,ll))
        else
!          call dgemm('n','n',idegfree2,nof_variables,numneighbours2-1,&
!          alpha,ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:nof_variables,ll),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:nof_variables,ll),idegfree2)
          sol_m(1:idegfree2,1:nof_variables,ll)=matmul(ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),matrix_1(1:numneighbours2-1,1:nof_variables,ll))
        end if
      end do

      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          ilocal_recon5(iconsidered)%gradients(ll, 1:number_of_dog, 1:nof_variables) = sol_m(1:number_of_dog, 1:nof_variables, ll)
        else
          ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, 1:nof_variables) = sol_m(1:idegfree2, 1:nof_variables, ll)
        end if
      end do

    else

      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          do iq = 1, imax
            if (ilocal_recon3(i)%ihexb(ll, iq + 1) .eq. n) then
              sols2(1:nof_variables, ll) = u_c(ilocal_recon3(i)%ihexl(ll, iq + 1))%val(1, 1:nof_variables)
            else
              sols2(1:nof_variables, ll) = iexsolhir(ilocal_recon3(i)%ihexn(ll, iq + 1))%sol(ilocal_recon3(i)%ihexl(ll, iq + 1), 1:nof_variables)
            end if

            if (wenwrt .eq. 3) then
              leftv(1:nof_variables) = sols2(1:nof_variables, ll)
              call cons2prim(n, leftv, mp_pinfl, gammal)
              sols2(1:nof_variables, ll) = leftv(1:nof_variables)
            end if
            if (per_rot .eq. 1) then
              if (ilocal_recon3(i)%periodicflag(ll, iq + 1) .eq. 1) then
!             write(2900+n,*),ielem(n,i)%xxc,ielem(n,i)%yyc,ielem(n,i)%zzc
                if (ielem(n, i)%xxc .gt. 0.0d0) then
                  tempxx = sols2(2, ll)
                  sols2(2, ll) = sols2(3, ll)
                  sols2(3, ll) = -tempxx
                end if
                if (ielem(n, i)%xxc .lt. 0.0d0) then
                  tempxx = sols2(2, ll)
                  sols2(2, ll) = -sols2(3, ll)
                  sols2(3, ll) = tempxx
                end if
              end if
            end if
            matrix_1(iq, 1:nof_variables, ll) = (sols2(1:nof_variables, ll) - sols1(1:nof_variables))
          end do

        else
          do iq = 1, numneighbours2 - 1
            if (ilocal_recon3(i)%ihexbc(ll, iq + 1) .eq. n) then
              sols2(1:nof_variables, ll) = u_c(ilocal_recon3(i)%ihexlc(ll, iq + 1))%val(1, 1:nof_variables)
            else
              sols2(1:nof_variables,ll)=iexsolhir(ilocal_recon3(i)%ihexnc(ll,iq+1))%sol(ilocal_recon3(i)%ihexlc(ll,iq+1),1:nof_variables)
            end if
            if (wenwrt .eq. 3) then
              leftv(1:nof_variables) = sols2(1:nof_variables, ll)
              call cons2prim(n, leftv, mp_pinfl, gammal)
              sols2(1:nof_variables, ll) = leftv(1:nof_variables)
            end if
              matrix_1(iq, 1:nof_variables, ll) = (sols2(1:nof_variables, ll) - sols1(1:nof_variables))
          end do
        end if
      end do

      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
!          call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables,ll),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables,ll),ielem(n,i)%idegfree)
          sol_m(1:ielem(n,i)%idegfree,1:nof_variables,ll)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables,ll))
        else

!          call dgemm('n','n',idegfree2,nof_variables,numneighbours2-1,&
!          alpha,ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:nof_variables,ll),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:nof_variables,ll),idegfree2)
          sol_m(1:idegfree2,1:nof_variables,ll)=matmul(ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),matrix_1(1:numneighbours2-1,1:nof_variables,ll))
        end if
      end do

      do ll = 1, ielem(n, i)%admis;
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          ilocal_recon5(iconsidered)%gradients(ll, 1:number_of_dog, 1:nof_variables) = sol_m(1:number_of_dog, 1:nof_variables, ll)
        else
          ilocal_recon5(iconsidered)%gradientsc(ll, 1:idegfree2, 1:nof_variables) = sol_m(1:idegfree2, 1:nof_variables, ll)
        end if
      end do

    end if
    deallocate(matrix_1, matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_mean_lsq

  subroutine compute_gradients_inner_mean_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitve variables of each interior cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(nof_variables)::sols1, sols2
    real, allocatable, dimension(:, :)::matrix_1
    real, allocatable, dimension(:, :)::matrix_2
    real, allocatable, dimension(:, :)::sol_m
    integer::i, var2, iq, lq, ll, imax
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::leftv

    allocate(matrix_1(number_of_nei - 1, nof_variables))
    allocate(matrix_2(nof_variables, number_of_dog))
    allocate(sol_m(number_of_dog, nof_variables))

    imax = number_of_nei - 1

    i = iconsidered
    sols1 = zero;
    sols2 = zero

    if (dimensiona .eq. 3) then

      ll = 1

      if (ilocal_recon3(i)%local .eq. 1) then
        matrix_1 = zero; matrix_2 = zero
        leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols1(2:4) = leftv(2:4)
        sols1(1) = leftv(5)/leftv(1)
        do iq = 1, imax
          leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          sols2(2:4) = leftv(2:4)
          sols2(1) = leftv(5)/leftv(1)
          matrix_1(iq, 1:nof_variables) = ((sols2(1:nof_variables) - sols1(1:nof_variables)))
        end do

!             call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables),ielem(n,i)%idegfree)

           sol_m(1:ielem(n,i)%idegfree,1:nof_variables)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables))

        do var2 = 2, 4
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:number_of_dog) = sol_m(1:number_of_dog, var2)
        end do
        ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = sol_m(1:number_of_dog, 1)
      else
        matrix_1 = zero; matrix_2 = zero
        leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols1(2:4) = leftv(2:4)
        sols1(1) = leftv(5)/leftv(1)
        do iq = 1, imax
          if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
            leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
          else
       leftv(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1:nof_variables)
          end if
          call cons2prim(n, leftv, mp_pinfl, gammal)
          sols2(2:4) = leftv(2:4)
          sols2(1) = leftv(5)/leftv(1)
          matrix_1(iq, 1:nof_variables) = ((sols2(1:nof_variables) - sols1(1:nof_variables)))
        end do

!          call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables),ielem(n,i)%idegfree)

          sol_m(1:ielem(n,i)%idegfree,1:nof_variables)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables))
        do var2 = 2, 4
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:number_of_dog) = sol_m(1:number_of_dog, var2)
        end do
        ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = sol_m(1:number_of_dog, 1)
      end if
    else

      ll = 1

      if (ilocal_recon3(i)%local .eq. 1) then
        matrix_1 = zero; matrix_2 = zero
        leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols1(2:3) = leftv(2:3)
        sols1(1) = leftv(4)/leftv(1)
        do iq = 1, imax
          leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
          call cons2prim(n, leftv, mp_pinfl, gammal)
          sols2(2:3) = leftv(2:3)
          sols2(1) = leftv(4)/leftv(1)
          matrix_1(iq, 1:3) = ((sols2(1:3) - sols1(1:3)))
        end do

!             call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables),ielem(n,i)%idegfree)

          sol_m(1:ielem(n,i)%idegfree,1:nof_variables)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables))
        do var2 = 2, 3
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:number_of_dog) = sol_m(1:number_of_dog, var2)
        end do
        ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = sol_m(1:number_of_dog, 1)
      else
        matrix_1 = zero; matrix_2 = zero
        leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols1(2:3) = leftv(2:3)
        sols1(1) = leftv(4)/leftv(1)
        do iq = 1, imax
          if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
            leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
          else
            leftv(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1:nof_variables)
          end if
          call cons2prim(n, leftv, mp_pinfl, gammal)
          sols2(2:3) = leftv(2:3)
          sols2(1) = leftv(4)/leftv(1)
          matrix_1(iq, 1:3) = ((sols2(1:3) - sols1(1:3)))
        end do

!             call dgemm('n','n',ielem(n,i)%idegfree,nof_variables,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:nof_variables),ielem(n,i)%idegfree)

          sol_m(1:ielem(n,i)%idegfree,1:nof_variables)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:nof_variables))
        do var2 = 2, 3
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:number_of_dog) = sol_m(1:number_of_dog, var2)
        end do
        ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = sol_m(1:number_of_dog, 1)
      end if

    end if

    deallocate(matrix_1, matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_inner_mean_lsq_viscous

  subroutine compute_gradients_inner_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(turbulenceequations + passivescalar)::sols1, sols2
    real, dimension(turbulenceequations + passivescalar, nof_variables)::sols_f
    real, dimension(dimensiona)::normal_all
    real::oov2, angle1, angle2
    integer::i, j, k, l, var2
    real, dimension(1:nof_variables)::leftv

    i = iconsidered
    sols_f = zero
    oov2 = 1.0d0/ielem(n, i)%totvolume

    if (dimensiona .eq. 3) then

      sols1(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)/u_c(i)%val(1, 1)

      do j = 1, ielem(n, i)%ifca
        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = (cos(angle1)*sin(angle2))
        normal_all(2) = (sin(angle1)*sin(angle2))
        normal_all(3) = (cos(angle2))
        sols2(1:turbulenceequations+passivescalar)=u_ct(ielem(n,i)%ineigh(j))%val(1,1:turbulenceequations+passivescalar)/u_c(ielem(n,i)%ineigh(j))%val(1,1)
        do k = 1, 3
        sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:3, var2) = sols_f(var2, 1:3)
        ilocal_recon3(i)%grads(4 + var2, 1:3) = sols_f(var2, 1:3)
      end do

    else
      sols1(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)/u_c(i)%val(1, 1)

      do j = 1, ielem(n, i)%ifca
        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = angle1
        normal_all(2) = angle2

        sols2(1:turbulenceequations+passivescalar)=u_ct(ielem(n,i)%ineigh(j))%val(1,1:turbulenceequations+passivescalar)/u_c(ielem(n,i)%ineigh(j))%val(1,1)

        do k = 1, 2
        sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)

        end do
      end do

      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:2, var2) = sols_f(var2, 1:2)
        ilocal_recon3(i)%grads(3 + var2, 1:2) = sols_f(var2, 1:2)
      end do

    end if

  end subroutine compute_gradients_inner_turb_ggs_viscous

  subroutine compute_gradients_turb_lsq(n, iconsidered, number_of_dog, number_of_nei) !check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(turbulenceequations + passivescalar)::sols1, sols2
    real, allocatable, dimension(:, :)::matrix_1
    real, allocatable, dimension(:, :)::matrix_2
    real, allocatable, dimension(:, :)::sol_m
    integer::i, var2, ll, iq, il, ih
    integer::imax

    allocate(matrix_1(number_of_nei - 1, turbulenceequations + passivescalar))
    allocate(matrix_2(turbulenceequations + passivescalar, number_of_dog))
    allocate(sol_m(number_of_dog, turbulenceequations + passivescalar))

    imax = number_of_nei - 1
    i = iconsidered
    sols1 = zero;
    sols2 = zero

    sols1(1:turbulenceequations + passivescalar) = u_ct(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:turbulenceequations + passivescalar)
    if (ilocal_recon3(i)%local .eq. 1) then

      do ll = 1, ielem(n, i)%admis;
        matrix_1 = zero; matrix_2 = zero
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
        do iq = 1, imax
         sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(ll,iq+1))%val(1,1:turbulenceequations+passivescalar)
         matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
        end do
        else
        do iq = 1, numneighbours2 - 1
        sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexlc(ll,iq+1))%val(1,1:turbulenceequations+passivescalar)
        matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
        end do
        end if
        if ((ees .ne. 5) .or. (ll .eq. 1)) then

! c          call dgemm('n','n',ielem(n,i)%idegfree,turbulenceequations+passivescalar,imax,&
! c          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
! c          ielem(n,i)%idegfree,matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! c imax,beta,sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar),ielem(n,i)%idegfree)
          sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:turbulenceequations+passivescalar))
        else

!         call dgemm('n','n',idegfree2,turbulenceequations+passivescalar,numneighbours2-1,&
!          alpha,ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:turbulenceequations+passivescalar),idegfree2)
           sol_m(1:idegfree2,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar))
        end if

        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          ilocal_recon5(iconsidered)%gradients2(ll,1:number_of_dog,1:turbulenceequations+passivescalar)=sol_m(1:number_of_dog,1:turbulenceequations+passivescalar)
        else
          ilocal_recon5(iconsidered)%gradientsc2(ll,1:idegfree2,1:turbulenceequations+passivescalar)=sol_m(1:idegfree2,1:turbulenceequations+passivescalar)
        end if
      end do

    else
      do ll = 1, ielem(n, i)%admis;
        matrix_1 = zero; matrix_2 = zero
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          do iq = 1, imax
            if (ilocal_recon3(i)%ihexb(ll, iq + 1) .eq. n) then
         sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(ll,iq+1))%val(1,1:turbulenceequations+passivescalar)
            else
              sols2(1:turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexn(ll,iq+1))%sol(ilocal_recon3(i)%ihexl(ll,iq+1),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if
              matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
          end do
        else
          do iq = 1, numneighbours2 - 1
            if (ilocal_recon3(i)%ihexbc(ll, iq + 1) .eq. n) then
              sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexlc(ll,iq+1))%val(1,1:turbulenceequations+passivescalar)
            else
              sols2(1:turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexnc(ll,iq+1))%sol(ilocal_recon3(i)%ihexlc(ll,iq+1),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if
              matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
          end do
        end if
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
!          call dgemm('n','n',ielem(n,i)%idegfree,turbulenceequations+passivescalar,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar),ielem(n,i)%idegfree)
          sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:turbulenceequations+passivescalar))
        else

!             call dgemm('n','n',idegfree2,turbulenceequations+passivescalar,numneighbours2-1,&
!          alpha,ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:turbulenceequations+passivescalar),idegfree2)
          sol_m(1:idegfree2,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll),matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar))
        end if
        if ((ees .ne. 5) .or. (ll .eq. 1)) then
          ilocal_recon5(iconsidered)%gradients2(ll,1:number_of_dog,1:turbulenceequations+passivescalar)=sol_m(1:number_of_dog,1:turbulenceequations+passivescalar)
        else
          ilocal_recon5(iconsidered)%gradientsc2(ll,1:idegfree2,1:turbulenceequations+passivescalar)=sol_m(1:idegfree2,1:turbulenceequations+passivescalar)
        end if
      end do
    end if

    deallocate(matrix_1)
    deallocate(matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_turb_lsq

  subroutine compute_gradients_turb_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each interior cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(turbulenceequations + passivescalar)::sols1, sols2
    real, allocatable, dimension(:, :)::matrix_1
    real, allocatable, dimension(:, :)::matrix_2
    real, allocatable, dimension(:, :)::sol_m
    integer::i, var2, iq, lq, ll, imax, ideg

    allocate(matrix_1(number_of_nei - 1, turbulenceequations + passivescalar))
    allocate(matrix_2(turbulenceequations + passivescalar, number_of_dog))
    allocate(sol_m(number_of_dog, turbulenceequations + passivescalar))

    imax = number_of_nei

    i = iconsidered
    sols1 = zero;
    sols2 = zero
    ideg = ielem(n, i)%idegfree
    ll = 1

    if (ilocal_recon3(i)%local .eq. 1) then
      matrix_1 = zero; matrix_2 = zero
  sols1(1:turbulenceequations + passivescalar) = u_ct(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                 u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1)

      do iq = 1, imax
        sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(1,iq+1))%val(1,1:turbulenceequations+passivescalar)/&
                                                   u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1)
        matrix_1(iq,1:turbulenceequations+passivescalar)=((sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
      end do

!         call dgemm('n','n',ielem(n,i)%idegfree,turbulenceequations+passivescalar,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar),ielem(n,i)%idegfree)

        ol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:turbulenceequations+passivescalar))

      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:number_of_dog, var2) = sol_m(1:number_of_dog, var2)
      end do

    else
      matrix_1 = zero; matrix_2 = zero
   sols1(1:turbulenceequations + passivescalar) = u_ct(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:turbulenceequations + passivescalar) &
                                                  /u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1)

      do iq = 1, imax
        if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
         sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(1,iq+1))%val(1,1:turbulenceequations+passivescalar)&
                                                  /u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1)

        else
          sols2(1:turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexn(1,iq+1))%sol(ilocal_recon3(i)%ihexl(1,iq+1),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
          iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1)
        end if
          matrix_1(iq,1:turbulenceequations+passivescalar)=((sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
      end do

!           call dgemm('n','n',ielem(n,i)%idegfree,turbulenceequations+passivescalar,imax,&
!          alpha,ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),&
!          ielem(n,i)%idegfree,matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar),ielem(n,i)%idegfree)
        sol_m(1:ielem(n,i)%idegfree,1:turbulenceequations+passivescalar)=matmul(ilocal_recon3(i)%invmat_stencilt(1:ielem(n,i)%idegfree,1:imax,ll),matrix_1(1:imax,1:turbulenceequations+passivescalar))
      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:number_of_dog, var2) = sol_m(1:number_of_dog, var2)
      end do
    end if

    deallocate(matrix_1, matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_turb_lsq_viscous

  subroutine compute_gradients_inner_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(1:nof_variables)::sols1, sols2, dudl, aver1
    real, dimension(1:nof_variables, dimensiona)::sols_f
    real, dimension(3)::normal_all
    real::oov2, titj, mp_pinfl, gammal, angle1, angle2
    integer::i, j, k, l
    real, dimension(1:nof_variables)::leftv

    i = iconsidered
    sols_f = zero
    oov2 = 1.0d0/ielem(n, i)%totvolume

    if (dimensiona .eq. 3) then

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols1(1:nof_variables) = leftv(1:nof_variables)
      sols1(5) = leftv(5)/leftv(1)

      do j = 1, ielem(n, i)%ifca
        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = (cos(angle1)*sin(angle2))
        normal_all(2) = (sin(angle1)*sin(angle2))
        normal_all(3) = (cos(angle2))

        leftv(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols2(1:nof_variables) = leftv(1:nof_variables)
        sols2(5) = leftv(5)/leftv(1)
        do k = 1, 3
          sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do k = 1, 3
        ilocal_recon3(i)%grads(1:3, k) = sols_f(2:4, k)
        ilocal_recon3(i)%grads(4, k) = sols_f(5, k)
      end do

    else

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols1(1:nof_variables) = leftv(1:nof_variables)
      sols1(4) = leftv(4)/leftv(1)

      do j = 1, ielem(n, i)%ifca
        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = angle1
        normal_all(2) = angle2

        leftv(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols2(1:nof_variables) = leftv(1:nof_variables)
        sols2(4) = leftv(4)/leftv(1)
        do k = 1, 2
          sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do k = 1, 2
        ilocal_recon3(i)%grads(1:2, k) = sols_f(2:3, k)
        ilocal_recon3(i)%grads(3, k) = sols_f(4, k)
      end do

    end if

  end subroutine compute_gradients_inner_mean_ggs_viscous

  subroutine compute_gradients_wall_mean_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei)!check all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each non-interior cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(nof_variables - 1)::sols1, sols2
    real, allocatable, dimension(:, :)::matrix_1
    real, allocatable, dimension(:, :)::matrix_2
    real, allocatable, dimension(:, :)::sol_m
    real, dimension(1:nof_variables)::matrix_3
    integer::i, var2, ii, k0, g0, ttk, ivvm, iq, lq
    real::attt
    integer::ll, imax
    real::mp_pinfl, gammal
    real, dimension(1:nof_variables)::leftv

    allocate(matrix_1(nof_variables, number_of_nei - 1))
    allocate(matrix_2(nof_variables, number_of_dog))
    allocate(sol_m(number_of_dog, nof_variables))

    imax = number_of_nei - 1

    if (dimensiona .eq. 3) then

      ll = 1
      i = iconsidered
      sols1 = zero;
      sols2 = zero

      ll = 1
      k0 = ilocal_recon3(i)%k0
      g0 = ilocal_recon3(i)%g0

      matrix_1 = zero; matrix_2 = zero; sol_m = zero;
      leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)

      sols1(2:4) = leftv(2:4)
      sols1(1) = leftv(5)/leftv(1)

      do iq = 1, imax
      if (ilocal_recon3(i)%local .eq. 1) then
        leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
      else
        if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
          leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
        else
          leftv(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1:nof_variables)
        end if
      end if

      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols2(2:4) = leftv(2:4)
      sols2(1) = leftv(5)/leftv(1)
      matrix_1(1:nof_variables-1,iq)=(ilocal_recon3(i)%volume(1,iq+1)*ilocal_recon3(i)%weightl(1,iq)*(sols2(1:nof_variables-1)-sols1(1:nof_variables-1)))

      matrix_1(2:4, iq) = matrix_1(2:4, iq) + ((sols1(2:4)*ilocal_recon3(i)%stencils(ll, iq, k0))/ilocal_recon3(i)%wallcoeff(k0))

      if (thermal .eq. 1) then
        matrix_1(1,iq)=matrix_1(1,iq)+((sols1(1)*ilocal_recon3(i)%stencils(ll,iq,g0))/ilocal_recon3(i)%wallcoefg(g0))-((wall_temp*ilocal_recon3(i)%stencils(ll,iq,g0))/ilocal_recon3(i)%wallcoefg(g0))
      end if

      end do
      matrix_3(1:nof_variables - 1) = -sols1(1:nof_variables - 1)
      matrix_3(1) = zero

      do var2 = 1, nof_variables - 1
        matrix_2 = zero
        if (var2 .gt. 1) then
        do iq = 1, imax

          do lq = 1, number_of_dog - 1
            matrix_2(var2, lq) = matrix_2(var2, lq) + matrix_1(var2, iq)*ilocal_recon3(i)%vellsq(iq, lq)
          end do
        end do
        else
        do iq = 1, imax
          do lq = 1, number_of_dog - 1
            matrix_2(var2, lq) = matrix_2(var2, lq) + matrix_1(var2, iq)*ilocal_recon3(i)%tempsq(iq, lq)
          end do
        end do
        end if
        if (var2 .eq. 1) then
          sol_m(1:number_of_dog-1,var2)=matmul(ilocal_recon3(i)%tempsqmat(1:number_of_dog-1,1:number_of_dog-1),matrix_2(var2,1:number_of_dog-1))
        else
          sol_m(1:number_of_dog-1,var2)=matmul(ilocal_recon3(i)%velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1),matrix_2(var2,1:number_of_dog-1))
        end if
      end do
      do var2 = 2, 4
        ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:idegfree) = -tolbig
        ivvm = 0
        do ttk = 1, number_of_dog
          if (ttk .eq. k0) cycle
          ivvm = ivvm + 1
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, ttk) = sol_m(ivvm, var2)
        end do
        attt = zero
        attt = -sols1(var2)
        do ttk = 1, number_of_dog
          if (ttk .ne. k0) &
            attt = attt - ilocal_recon5(iconsidered)%velocitydof(var2 - 1, ttk)* &
                  ilocal_recon3(i)%wallcoeff(ttk)
        end do
        attt = attt/ilocal_recon3(i)%wallcoeff(k0)
        ilocal_recon5(iconsidered)%velocitydof(var2 - 1, k0) = attt

      end do

      ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = -tolbig
      ivvm = 0
      do ttk = 1, number_of_dog
        if (ttk .eq. g0) cycle
        ivvm = ivvm + 1
        ilocal_recon5(iconsidered)%gradientstemp(ttk) = sol_m(ivvm, 1)
      end do
      attt = zero
      if (thermal .eq. 1) then
        attt = wall_temp - sols1(1)
      end if

      do ttk = 1, number_of_dog
        if (ttk .ne. g0) &
          attt = attt - ilocal_recon5(iconsidered)%gradientstemp(ttk)* &
                 ilocal_recon3(i)%wallcoefg(ttk)
      end do
      attt = attt/ilocal_recon3(i)%wallcoefg(g0)
      ilocal_recon5(iconsidered)%gradientstemp(g0) = attt

    else        !2d

      i = iconsidered
      sols1 = zero;
      sols2 = zero
      ll = 1

      k0 = ilocal_recon3(i)%k0
      g0 = ilocal_recon3(i)%g0

      matrix_1 = zero; matrix_2 = zero; sol_m = zero;
      leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)

      sols1(2:3) = leftv(2:3)
      sols1(1) = leftv(4)/leftv(1)

      do iq = 1, imax
        if (ilocal_recon3(i)%local .eq. 1) then
          leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
        else
          if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
            leftv(1:nof_variables) = u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1:nof_variables)
          else
            leftv(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1:nof_variables)
          end if
        end if

        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols2(2:3) = leftv(2:3)
        sols2(1) = leftv(4)/leftv(1)
        matrix_1(1:3, iq) = (ilocal_recon3(i)%volume(1, iq + 1)*ilocal_recon3(i)%weightl(1, iq)*(sols2(1:3) - sols1(1:3)))
        matrix_1(2:3, iq) = matrix_1(2:3, iq) + ((sols1(2:3)*ilocal_recon3(i)%stencils(ll, iq, k0))/ilocal_recon3(i)%wallcoeff(k0))

        if (thermal .eq. 1) then
          matrix_1(1,iq)=matrix_1(1,iq)+((sols1(1)*ilocal_recon3(i)%stencils(ll,iq,g0))/ilocal_recon3(i)%wallcoefg(g0))-((wall_temp*ilocal_recon3(i)%stencils(ll,iq,g0))/ilocal_recon3(i)%wallcoefg(g0))
        end if

      end do
      matrix_3(1:3) = -sols1(1:3)
      matrix_3(1) = zero

      do var2 = 1, nof_variables - 1
        matrix_2 = zero
        if (var2 .gt. 1) then
        do iq = 1, imax
          do lq = 1, number_of_dog - 1
            matrix_2(var2, lq) = matrix_2(var2, lq) + matrix_1(var2, iq)*ilocal_recon3(i)%vellsq(iq, lq)
          end do
        end do
        else
        do iq = 1, imax
          do lq = 1, number_of_dog - 1
            matrix_2(var2, lq) = matrix_2(var2, lq) + matrix_1(var2, iq)*ilocal_recon3(i)%tempsq(iq, lq)
          end do
        end do
        end if
        if (var2 .eq. 1) then
          sol_m(1:number_of_dog-1,var2)=matmul(ilocal_recon3(i)%tempsqmat(1:number_of_dog-1,1:number_of_dog-1),matrix_2(var2,1:number_of_dog-1))
        else
          sol_m(1:number_of_dog-1,var2)=matmul(ilocal_recon3(i)%velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1),matrix_2(var2,1:number_of_dog-1))
        end if
      end do
      do var2 = 2, 3

        ilocal_recon5(iconsidered)%velocitydof(var2 - 1, 1:idegfree) = -tolbig
        ivvm = 0
        do ttk = 1, number_of_dog
          if (ttk .eq. k0) cycle
          ivvm = ivvm + 1
          ilocal_recon5(iconsidered)%velocitydof(var2 - 1, ttk) = sol_m(ivvm, var2)
        end do
        attt = zero
        attt = -sols1(var2)
        do ttk = 1, number_of_dog
          if (ttk .ne. k0) &
            attt = attt - ilocal_recon5(iconsidered)%velocitydof(var2 - 1, ttk)* &
                  ilocal_recon3(i)%wallcoeff(ttk)
        end do
        attt = attt/ilocal_recon3(i)%wallcoeff(k0)
        ilocal_recon5(iconsidered)%velocitydof(var2 - 1, k0) = attt
      end do

      ilocal_recon5(iconsidered)%gradientstemp(1:number_of_dog) = -tolbig

      ivvm = 0
      do ttk = 1, number_of_dog
        if (ttk .eq. g0) cycle
        ivvm = ivvm + 1
        ilocal_recon5(iconsidered)%gradientstemp(ttk) = sol_m(ivvm, 1)
      end do
      attt = zero
      if (thermal .eq. 1) then
        attt = wall_temp - sols1(1)
      end if
      do ttk = 1, number_of_dog
        if (ttk .ne. g0) &
          attt = attt - ilocal_recon5(iconsidered)%gradientstemp(ttk)* &
                 ilocal_recon3(i)%wallcoefg(ttk)
      end do
      attt = attt/ilocal_recon3(i)%wallcoefg(g0)
      ilocal_recon5(iconsidered)%gradientstemp(g0) = attt

    end if

    deallocate(matrix_1, matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_wall_mean_lsq_viscous

  subroutine compute_gradients_mix_turb_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulnece variables of each non-interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(turbulenceequations + passivescalar)::sols1, sols2
    real, dimension(turbulenceequations + passivescalar, 3)::sols_f
    real, dimension(3)::normal_all, temp_vert
    real::oov2, titj, mp_pinfl, gammal, angle1, angle2, nx, ny, nz
    integer::i, j, k, l, var2, b_code, facex, n_node, imax
    real, dimension(1:nof_variables)::leftv, srf_speed, srf_speedrot, rightv
    real, dimension(1:dimensiona)::pox, poy, poz, cords
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(turbulenceequations)::cturbl, cturbr
    real, dimension(1:nof_variables)::cright_rot, cleft_rot
    integer::ibfc

    if (dimensiona .eq. 3) then

      i = iconsidered
      sols_f = zero
      oov2 = 1.0d0/ielem(n, i)%totvolume

      sols1(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)/u_c(i)%val(1, 1)

      do j = 1, ielem(n, i)%ifca
        facex = j

        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = (cos(angle1)*sin(angle2))
        normal_all(2) = (sin(angle1)*sin(angle2))
        normal_all(3) = (cos(angle2))
        nx = normal_all(1); ny = normal_all(2); nz = normal_all(3)

        if (ielem(n, i)%ineighb(j) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in my cpu
              sols2(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(j))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                             u_c(ielem(n, i)%ineigh(j))%val(1, 1)
            else
              !not periodic ones in my cpu
              call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)
              if (ielem(n, iconsidered)%types_faces(facex) .eq. 5) then
                n_node = 4
              else
                n_node = 3
              end if

              cords = cordinates3(n, nodes_list, n_node)
              pox(1) = cords(1); poy(1) = cords(2); poz(1) = cords(3)

              leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
              cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
              b_code = ibound(n, ielem(n, i)%ibounds(j))%icode
              call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

              sols2(1:turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)/rightv(1)

            end if
          else
         sols2(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(j))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                        u_c(ielem(n, i)%ineigh(j))%val(1, 1)

          end if
        else        !in other cpus they can only be periodic or mpi neighbours
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in other cpu
              if (fastest .eq. 1) then
                sols2(1:turbulenceequations+passivescalar)=solchanger(ielem(n,i)%ineighn(j))%sol(ielem(n,i)%q_face(j)%q_mapl(1),6:5+turbulenceequations+passivescalar)/&
                                                           solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1)
              else
                sols2(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 6:5 + turbulenceequations + passivescalar)/ &
                iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1)
              end if
            end if
          else

            if (fastest .eq. 1) then
              sols2(1:turbulenceequations+passivescalar)=solchanger(ielem(n,i)%ineighn(j))%sol(ielem(n,i)%q_face(j)%q_mapl(1),6:5+turbulenceequations+passivescalar)/&
                                                         solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1)
            else
              sols2(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 6:5 + turbulenceequations + passivescalar)/ &
              iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1)
            end if
          end if
        end if

        do k = 1, 3
          sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:3, var2) = sols_f(var2, 1:3)
        ilocal_recon3(i)%grads(4 + var2, 1:3) = sols_f(var2, 1:3)
      end do

    else        !2d

      i = iconsidered
      sols_f = zero
      oov2 = 1.0d0/ielem(n, i)%totvolume

      sols1(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)/u_c(i)%val(1, 1)

      do j = 1, ielem(n, i)%ifca
        facex = j
        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = angle1
        normal_all(2) = angle2
        nx = normal_all(1); ny = normal_all(2)

        if (ielem(n, i)%ineighb(j) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in my cpu
              sols2(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(j))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                             u_c(ielem(n, i)%ineigh(j))%val(1, 1)
            else
              !not periodic ones in my cpu
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              cords = cordinates2(n, nodes_list, n_node)
              pox(1) = cords(1); poy(1) = cords(2);
              leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
              cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
              b_code = ibound(n, ielem(n, i)%ibounds(j))%icode
              call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              sols2(1:turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)/rightv(1)
            end if
          else
         sols2(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(j))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                        u_c(ielem(n, i)%ineigh(j))%val(1, 1)

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in other cpu
              if (fastest .eq. 1) then
                sols2(1:turbulenceequations+passivescalar)=solchanger(ielem(n,i)%ineighn(j))%sol(ielem(n,i)%q_face(j)%q_mapl(1),5:4+turbulenceequations+passivescalar)&
                                                          /solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1)
              else
                sols2(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 5:4 + turbulenceequations + passivescalar)/ &
                iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1)
              end if
            end if
          else

            if (fastest .eq. 1) then
              sols2(1:turbulenceequations+passivescalar)=solchanger(ielem(n,i)%ineighn(j))%sol(ielem(n,i)%q_face(j)%q_mapl(1),5:4+turbulenceequations+passivescalar)&
                                                        /solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1)
            else
              sols2(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 5:4 + turbulenceequations + passivescalar)/ &
              iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1)
            end if
          end if
        end if

        do k = 1, 2
          sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do var2 = 1, turbulenceequations + passivescalar
        ilocal_recon5(iconsidered)%gradientsturb(1, 1:2, var2) = sols_f(var2, 1:2)
        ilocal_recon3(i)%grads(3 + var2, 1:2) = sols_f(var2, 1:2)
      end do

    end if

  end subroutine compute_gradients_mix_turb_ggs_viscous

  subroutine compute_gradients_wall_turb_lsq_viscous(n, iconsidered, number_of_dog, number_of_nei) !check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each non-interior cell using the least-squares
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(1:turbulenceequations + passivescalar)::sols1, sols2
    real, allocatable, dimension(:, :)::matrix_1
    real, allocatable, dimension(:, :)::matrix_2
    real, allocatable, dimension(:, :)::sol_m
    real, dimension(1:turbulenceequations + passivescalar)::matrix_3
    integer::i, var2, ii, k0, g0, ttk, ivvm, iq, lq, imax
    real::attt
    integer::ll

    allocate(matrix_1(1:turbulenceequations + passivescalar, number_of_nei - 1))
    allocate(matrix_2(number_of_dog, 1:turbulenceequations + passivescalar))
    allocate(sol_m(number_of_dog, 1:turbulenceequations + passivescalar))

    ll = 1

    i = iconsidered
    sols1 = zero;
    sols2 = zero

    k0 = ilocal_recon3(i)%k0

    matrix_1 = zero; matrix_2 = zero; sol_m = zero;
  sols1(1:turbulenceequations + passivescalar) = u_ct(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1:turbulenceequations + passivescalar)/ &
                                                 u_c(ilocal_recon3(i)%ihexl(1, 1))%val(1, 1)

    do iq = 1, imax
    if (ilocal_recon3(i)%local .eq. 1) then
      sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(1,iq+1))%val(1,1:turbulenceequations+passivescalar)/&
                                                 u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1)
    else
      if (ilocal_recon3(i)%ihexb(1, iq + 1) .eq. n) then
        sols2(1:turbulenceequations+passivescalar)=u_ct(ilocal_recon3(i)%ihexl(1,iq+1))%val(1,1:turbulenceequations+passivescalar)/&
                                                   u_c(ilocal_recon3(i)%ihexl(1, iq + 1))%val(1, 1)
      else
        sols2(1:turbulenceequations+passivescalar)=iexsolhir(ilocal_recon3(i)%ihexn(1,iq+1))%sol(ilocal_recon3(i)%ihexl(1,iq+1),6:5+turbulenceequations+passivescalar)/&
                                                   iexsolhir(ilocal_recon3(i)%ihexn(1, iq + 1))%sol(ilocal_recon3(i)%ihexl(1, iq + 1), 1)
      end if
    end if
      matrix_1(1:turbulenceequations+passivescalar,iq)=(ilocal_recon3(i)%volume(1,iq+1)*ilocal_recon3(i)%weightl(1,iq)*(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
      matrix_1(1:turbulenceequations+passivescalar,iq)=matrix_1(1:turbulenceequations+passivescalar,iq)+((sols1(1:turbulenceequations+passivescalar)*ilocal_recon3(i)%stencils(ll,iq,k0))/ilocal_recon3(i)%wallcoeff(k0))
    end do
    matrix_3(1:turbulenceequations + passivescalar) = -sols1(1:turbulenceequations + passivescalar)
    if (turbulencemodel .eq. 2) then
      matrix_3(2) = 60.0d0*visc/(beta_i1*(ielem(n, iconsidered)%walldist**2))
    end if
    do var2 = 1, turbulenceequations + passivescalar
      matrix_2 = zero
      do iq = 1, imax
        do lq = 1, number_of_dog - 1
          matrix_2(var2, lq) = matrix_2(var2, lq) + matrix_1(var2, iq)*ilocal_recon3(i)%vellsq(iq, lq)
        end do
      end do
        sol_m(1:number_of_dog-1,var2)=matmul(ilocal_recon3(i)%velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1),matrix_2(var2,1:number_of_dog-1))
    end do

    do var2 = 1, turbulenceequations + passivescalar

      ilocal_recon5(iconsidered)%gradientsturb(1, 1:number_of_dog, var2) = -tolbig
      ivvm = 0
      do ttk = 1, number_of_dog
        if (ttk .eq. k0) cycle
        ivvm = ivvm + 1
        ilocal_recon5(iconsidered)%gradientsturb(1, ttk, var2) = sol_m(ivvm, var2)
      end do
      attt = zero
      attt = -sols1(var2)
      do ttk = 1, number_of_dog
        if (ttk .ne. k0) &
          attt = attt - ilocal_recon5(iconsidered)%gradientsturb(1, ttk, var2)* &
          ilocal_recon3(i)%wallcoeff(ttk)
      end do
      attt = attt/ilocal_recon3(i)%wallcoeff(k0)
      ilocal_recon5(iconsidered)%gradientsturb(1, k0, var2) = attt

    end do

    deallocate(matrix_1, matrix_2)
    deallocate(sol_m)

  end subroutine compute_gradients_wall_turb_lsq_viscous

  subroutine compute_gradients_mix_mean_ggs_viscous(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each non-interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(nof_variables)::sols1, sols2, dudl, aver1
    real, dimension(nof_variables, 3)::sols_f
    real, dimension(3)::normal_all
    real::oov2, titj, mp_pinfl, gammal, angle1, angle2, nx, ny, nz
    integer::i, j, k, l, b_code, facex, n_node, imax
    real, dimension(1:nof_variables)::leftv, srf_speed, srf_speedrot, rightv
    real, dimension(1:dimensiona)::pox, poy, poz, cords
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(turbulenceequations)::cturbl, cturbr
    real, dimension(1:nof_variables)::cright_rot, cleft_rot
    integer::ibfc

    if (dimensiona .eq. 3) then

      i = iconsidered
      sols_f = zero
      oov2 = 1.0d0/ielem(n, i)%totvolume

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols1(1:nof_variables) = leftv(1:nof_variables)
      sols1(5) = leftv(5)/leftv(1)

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)

      do j = 1, ielem(n, i)%ifca
        facex = j
        b_code = 0

        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = (cos(angle1)*sin(angle2))
        normal_all(2) = (sin(angle1)*sin(angle2))
        normal_all(3) = (cos(angle2))
        nx = normal_all(1); ny = normal_all(2); nz = normal_all(3)

        if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
          !retrieve rotational velocity in case of rotating reference frame to calculate
          !the correct value of the boundary condition
          srf_speed(2:4) = ilocal_recon3(i)%rotvel(j, 1, 1:3)
          call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
        end if

        if (ielem(n, i)%ineighb(j) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 50)) then        !periodic in my cpu
              sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)
              if (per_rot .eq. 1) then
                sols2(2:4) = rotate_per_1(sols2(2:4), ibound(n, ielem(n, i)%ibounds(j))%icode, angle_per)
              end if
            else
              !not periodic ones in my cpu

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

              leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(j))%icode
              call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              sols2(1:nof_variables) = rightv(1:nof_variables)

            end if
          else
            sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 50)) then        !periodic in other cpu
              if (fastest .eq. 1) then
                sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
              else
                sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
              end if
              if (per_rot .eq. 1) then
                sols2(2:4) = rotate_per_1(sols2(2:4), ibound(n, ielem(n, i)%ibounds(j))%icode, angle_per)
              end if
            end if
          else

            if (fastest .eq. 1) then
              sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
            else
              sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
            end if
          end if
        end if

        leftv(1:nof_variables) = sols2(1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols2(1:nof_variables) = leftv(1:nof_variables)
        sols2(5) = leftv(5)/leftv(1)

        do k = 1, 3
          sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do k = 1, 3
        ilocal_recon3(i)%grads(1:3, k) = sols_f(2:4, k)
        ilocal_recon3(i)%grads(4, k) = sols_f(5, k)
      end do

    else

      i = iconsidered
      sols_f = zero
      oov2 = 1.0d0/ielem(n, i)%totvolume

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols1(1:nof_variables) = leftv(1:nof_variables)
      sols1(4) = leftv(4)/leftv(1)

      leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)

      do j = 1, ielem(n, i)%ifca
        facex = j

        b_code = 0

        angle1 = ielem(n, i)%faceanglex(j)
        angle2 = ielem(n, i)%faceangley(j)
        normal_all(1) = angle1
        normal_all(2) = angle2
        nx = normal_all(1); ny = normal_all(2)

        if (ielem(n, i)%ineighb(j) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in my cpu
              sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)
            else
              !not periodic ones in my cpu
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              cords = cordinates2(n, nodes_list, n_node)
              pox(1) = cords(1); poy(1) = cords(2)
              leftv(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(j))%icode
              call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              sols2(1:nof_variables) = rightv(1:nof_variables)
            end if
          else
            sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(1, 1:nof_variables)
          end if
        else        !in other cpus they can only be periodic or mpi neighbours
          if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in other cpu
              if (fastest .eq. 1) then
                sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
              else
                sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
              end if
            end if
          else

            if (fastest .eq. 1) then
              sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
            else
              sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
            end if
          end if
        end if

        leftv(1:nof_variables) = sols2(1:nof_variables)
        call cons2prim(n, leftv, mp_pinfl, gammal)
        sols2(1:nof_variables) = leftv(1:nof_variables)
        if ((b_code .eq. 4) .and. (thermal .eq. 1)) then
          sols2(4) = wall_temp
        else
          sols2(4) = leftv(4)/leftv(1)
        end if
        do k = 1, 2
          sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
        end do
      end do

      do k = 1, 2
        ilocal_recon3(i)%grads(1:2, k) = sols_f(2:3, k)
        ilocal_recon3(i)%grads(3, k) = sols_f(4, k)
      end do

    end if

  end subroutine compute_gradients_mix_mean_ggs_viscous

  subroutine compute_gradients_center(n, iconsidered)
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(1:nof_variables)::sols1, sols2, dudl, aver1
    real::oov2, titj
    integer::i, j, k, l, iex
    i = iconsidered

    if (dimensiona .eq. 3) then
    do iex = 1, 3
      ilocal_recon3(i)%grads(iex, 1:3) = ilocal_recon3(i)%uleftv(1:3, iex + 1, 1, 1)
    end do
    ilocal_recon3(i)%grads(4, 1:3) = ilocal_recon3(i)%uleftv(1:3, 1, 1, 1)
    else
    do iex = 1, 2
      ilocal_recon3(i)%grads(iex, 1:2) = ilocal_recon3(i)%uleftv(1:2, iex + 1, 1, 1)
    end do
    ilocal_recon3(i)%grads(3, 1:2) = ilocal_recon3(i)%uleftv(1:2, 1, 1, 1)
    end if

  end subroutine compute_gradients_center

  subroutine compute_gradients_mix_mean_ggs_viscous_av(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the averaged primitive variables of each non-interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(1:nof_variables)::sols1, sols2
    real, dimension(1:nof_variables, 3)::sols_f
    real, dimension(3)::normal_all, temp_vert
    real::oov2, titj, mp_pinfl, gammal, angle1, angle2, nx, ny, nz
    integer::i, j, k, l, var2, b_code, facex, n_node, ind1
    real, dimension(1:nof_variables)::leftv, srf_speed, srf_speedrot, rightv
    real, dimension(1:dimensiona)::pox, poy, poz, cords
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(turbulenceequations)::cturbl, cturbr
    real, dimension(1:nof_variables)::cright_rot, cleft_rot
    integer::ibfc

    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if

    i = iconsidered
    sols_f = zero
    oov2 = 1.0d0/ielem(n, i)%totvolume

    leftv(1:nof_variables) = u_c(i)%val(ind1, 1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    sols1(1:nof_variables) = leftv(1:nof_variables)
    sols1(5) = leftv(5)/leftv(1)

    leftv(1:nof_variables) = u_c(i)%val(ind1, 1:nof_variables)

    do j = 1, ielem(n, i)%ifca
      facex = j
      b_code = 0

      angle1 = ielem(n, i)%faceanglex(j)
      angle2 = ielem(n, i)%faceangley(j)
      normal_all(1) = (cos(angle1)*sin(angle2))
      normal_all(2) = (sin(angle1)*sin(angle2))
      normal_all(3) = (cos(angle2))
      nx = normal_all(1); ny = normal_all(2); nz = normal_all(3)

      if (ielem(n, i)%ineighb(j) .eq. n) then        !my cpu only
        if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
          if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in my cpu
            sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(ind1, 1:nof_variables)
          else
            !not periodic ones in my cpu
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
            leftv(1:nof_variables) = u_c(i)%val(ind1, 1:nof_variables)
            b_code = ibound(n, ielem(n, i)%ibounds(j))%icode
            call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
            sols2(1:nof_variables) = rightv(1:nof_variables)
          end if
        else
          sols2(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(ind1, 1:nof_variables)
        end if
      else        !in other cpus they can only be periodic or mpi neighbours
        if (ielem(n, i)%ibounds(j) .gt. 0) then        !check for boundaries
          if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        !periodic in other cpu
            if (fastest .eq. 1) then
              sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
            else
              sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
            end if
          end if
        else
          if (fastest .eq. 1) then
            sols2(1:nof_variables) = solchanger(ielem(n, i)%ineighn(j))%sol(ielem(n, i)%q_face(j)%q_mapl(1), 1:nof_variables)
          else
            sols2(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(j)))%sol &
                                    (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(j)), 1:nof_variables)
          end if
        end if
      end if

      leftv(1:nof_variables) = sols2(1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols2(1:nof_variables) = leftv(1:nof_variables)

      if ((b_code .eq. 4) .and. (thermal .eq. 1)) then
        sols2(5) = wall_temp
      else
        sols2(5) = leftv(5)/leftv(1)
      end if
      do k = 1, 3
        sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
      end do
    end do

    do k = 1, 3
      ilocal_recon3(i)%gradsav(1:3, k) = sols_f(2:4, k)
      ilocal_recon3(i)%gradsav(4, k) = sols_f(5, k)
    end do

  end subroutine compute_gradients_mix_mean_ggs_viscous_av

  subroutine compute_gradients_inner_mean_ggs_viscous_av(n, iconsidered, number_of_dog, number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the averaged primitive variables of each interior cell using the green-gauss algorithm
    implicit none
    integer, intent(in)::n, iconsidered, number_of_dog, number_of_nei
    real, dimension(nof_variables)::sols1, sols2, leftv
    real, dimension(nof_variables, 3)::sols_f
    real, dimension(3)::normal_all
    real::oov2, mp_pinfl, gammal, angle1, angle2
    integer::i, j, k, l, ind1

    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if

    i = iconsidered
    sols_f = zero
    oov2 = 1.0d0/ielem(n, i)%totvolume

    leftv(1:nof_variables) = u_c(i)%val(ind1, 1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    sols1(1:nof_variables) = leftv(1:nof_variables)
    sols1(5) = leftv(5)/leftv(1)

    do j = 1, ielem(n, i)%ifca
      angle1 = ielem(n, i)%faceanglex(j)
      angle2 = ielem(n, i)%faceangley(j)
      normal_all(1) = (cos(angle1)*sin(angle2))
      normal_all(2) = (sin(angle1)*sin(angle2))
      normal_all(3) = (cos(angle2))

      leftv(1:nof_variables) = u_c(ielem(n, i)%ineigh(j))%val(ind1, 1:nof_variables)
      call cons2prim(n, leftv, mp_pinfl, gammal)
      sols2(1:nof_variables) = leftv(1:nof_variables)
      sols2(5) = leftv(5)/leftv(1)

      do k = 1, 3
       sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem(n,i)%surf(j)*oov2)
      end do
    end do

    do k = 1, 3
      ilocal_recon3(i)%gradsav(1:3, k) = sols_f(2:4, k)
      ilocal_recon3(i)%gradsav(4, k) = sols_f(5, k)
    end do

  end subroutine compute_gradients_inner_mean_ggs_viscous_av

end module gradients
