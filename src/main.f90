program ucns3d
  ! @author
  ! Panagiotis Tsoutsanis & Antonis Foivos Antoniadis
  ! copyright: Panagiotis Tsoutsanis & Antonis Foivos Antoniadis
  ! version: 3.0
  ! Main Driver of UCNS3D code
  use mpi_info
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
  use advance
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
  use moodr
  use omp_lib
  use parameters

  implicit none

  external metis_partmeshdual
  external parmetis_v3_partmeshkway

  ! call mpi_init(ierror)
  call mpi_init_thread(mpi_thread_funneled, provided, ierror)
  call mpi_comm_size(mpi_comm_world, isize, ierror)
  call mpi_comm_rank(mpi_comm_world, n, ierror)

  call open_input1(n, itt)    ! open the input files
  call tolerances    ! setup the tolerances values
  call read_ucns3d   ! Read all the parameter files
  call close_input1(n, itt)    ! close the input files

  if (n .eq. 0)
    call translate_mesh    ! translate the mesh from fluent msh format to native format

  call mpi_barrier(mpi_comm_world, ierror)
  call timing(n, cpux1, cpux2, cpux3, cpux4, cpux5, cpux6, timex1, timex2, timex3, timex4, timex5, timex6)    ! start the timers
  cpux1(1) = mpi_wtime()

  call open_arbitrary(n, imaxe, imaxn, imaxb)    ! open the grid files
  call shallocation(ieshape, imaxe)    ! allocatearrays for shape of each element
  call mpi_barrier(mpi_comm_world, ierror)
  call find_shape(n, imaxe, ieshape)    ! find the shape each element
  call checkres    ! check the existence of RESTART/CHECKPOINT files
  call xmpiallocate(xmpie, xmpil, xmpin, xmpielrank, xmpinrank, imaxe, imaxn, nproc)    ! allocatememory for local and global numbering of elements and nodes
  call mpi_barrier(mpi_comm_world, ierror)

  if (emetis .lt. 6) then    ! choose a grid partitioning property if emetis<6 use serial METIS
    if (n .eq. 0) then
      if (emetis .eq. 1) call partitioner1(n, imaxe, imaxn, xmpie, ieshape)
      if (emetis .eq. 2) call partitioner2(n, imaxe, imaxn, xmpie, ieshape)
      if (emetis .eq. 3) call partitioner3(n, imaxe, imaxn, xmpie, ieshape)
      if (emetis .eq. 4) call partitioner4(n, imaxe, imaxn, xmpie, ieshape)
      if (emetis .eq. 5) call partitioner5(n, imaxe, imaxn, xmpie, ieshape)
    end if
    call mpi_bcast(xmpie, imaxe, mpi_integer, 0, mpi_comm_world, ierror)
  else
    if (emetis .eq. 6) call partitioner6(n, imaxe, imaxn, xmpie, ieshape)
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  call xmpifind(xmpie, xmpin, xmpielrank, xmpinrank, imaxe, imaxn, nproc)    ! determine the number of elements in this process
  call elallocation(n, xmpie, xmpielrank, ielem, imaxe, ieshape, itestcase, imaxb, ibound, xmin, xmax, ymin, ymax, zmin, zmax)    ! allocatethe appropriate memory for each elements
  call read_input(n, xmpielrank, xmpinrank, xmpie, xmpin, ielem, inode, imaxn, imaxe, ibound, imaxb, xmpinnumber, scaler, inoder)    ! read the grid files and populate the allocated memory values for vertex coordinates and numbering
  call fix_offsets_local(n)
  call determine_size(n, iorder, iselem, iselemt, ioverst, ioverto, ilx, numneighbours, idegfree, imaxdegfree, iextend)    ! determing the stencil sizes, number of polynomial coefficients etc
  call gaussian_points(igqrules, numberofpoints, numberofpoints2)    ! establish the number of Gausian quadrature points for each element
  call quad_alloc(numberofpoints, numberofpoints2)    ! allocatethe memory required for the Gaussian integration rules
  call shde_allocation(ieshape, imaxe)    ! deallocatememory for the shape allocation
  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_1", cpux3(1) - cpux1(1)    ! write in file the total wall clock time taken so far
  end if

  call allocate2
  call neighbourss(n, ielem, imaxe, imaxn, xmpie, xmpin, xmpielrank, restart, inoder)
  call allocate3
  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_2", cpux3(1) - cpux1(1)
  end if

  call geometry_calc
  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_3", cpux3(1) - cpux1(1)
  end if

  call read_bound(n, imaxb, ibound, xmpielrank)

  if (dg .eq. 1) then
    if (filtering .eq. 1) then
      allocate(modal_filter(1:idegfree), modal_filter_strong(1:idegfree), modal_filter_weak(1:idegfree))
    end if
  end if

  if (adda .eq. 1) then
    allocate(adda_filter_strong(1:idegfree))
    allocate(adda_filter_weak(1:idegfree))
  end if

  call apply_boundary(n, xper, yper, zper, iperiodicity, xmpielrank)
  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finished applying boundary conditions"
    cpux3(1) = mpi_wtime()
    write(63, *) "time1=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  call xmpilocal
  call count_walls
  call mpi_barrier(mpi_comm_world, ierror)

  if (lowmem .eq. 0) call globalistx(n, xmpie, xmpil, xmpielrank, imaxe, isize, centerr, glneigh, ielem)
  if (lowmem .eq. 1) call globalist(n, xmpie, xmpil, xmpielrank, imaxe, isize, centerr, glneigh, glneighper, ielem)

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finished obtaining neighbours within my cpu"
    cpux3(1) = mpi_wtime()
    write(63, *) "time2=", cpux3(1) - cpux1(1)
    close(63)
  end if

  if (dimensiona .eq. 3) then
    call cons(n, iconr, icons, iperiodicity, xmpielrank, isize, iconrpa, iconrpm, iconspo, xper, yper, zper, iconrpf, numneighbours, typesten)
  else
    call cons2d(n, iconr, icons, iperiodicity, xmpielrank, isize, iconrpa, iconrpm, iconspo, xper, yper, zper, iconrpf, numneighbours, typesten)
  end if

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finished obtaining neighbours within my cpu"
    cpux3(1) = mpi_wtime()
    write(63, *) "time22=", cpux3(1) - cpux1(1)
    close(63)
  end if

  if (lowmem .eq. 0)
    call globalistx2(n, xmpie, xmpil, xmpielrank, imaxe, isize, centerr, glneigh, ielem)
  if (lowmem .eq. 1)
    call globalist2(n, xmpie, xmpil, xmpielrank, imaxe, isize, centerr, glneigh, glneighper, ielem)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finished obtaining neighbours across all cpu"
    cpux3(1) = mpi_wtime()
    close(63)
  end if

  if (ischeme .gt. 1) then
    call mpi_barrier(mpi_comm_world, ierror)
    cpux3(1) = mpi_wtime()
    call allocate5
    if (lowmem .eq. 0) call detstenx(n)
    if (lowmem .eq. 1) call detsten(n)
    if (n .eq. 0) then
      open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
      write(63, *) "finished obtaining the central stencils across all cpu"
      cpux3(1) = mpi_wtime()
      write(63, *) "time4=", cpux3(1) - cpux1(1)
      close(63)
    end if
    call localstallocation(n, xmpielrank, ilocalstencil, ilocalstencilper, typesten, numneighbours)
    call mpi_barrier(mpi_comm_world, ierror)
    if (n .eq. 0) then
      cpux3(1) = mpi_wtime()
      write(100 + n, *) "time_i4", cpux3(1) - cpux1(1)
    end if
    if ((ees .eq. 0) .or. (ees .ge. 4)) then
      if (lowmem .eq. 0) call stenciilsx(n)
      if (lowmem .eq. 1) call stenciils(n)
    end if
    if ((ees .gt. 0) .and. (ees .le. 2)) then
      if (lowmem .eq. 0) call stenciils_eesx(n)
      if (lowmem .eq. 1) call stenciils_ees(n)
    end if
    if (ees .eq. 3) call stencils3(n)
    deallocate(ilocalallelg)
    deallocate(ilocalallelgper)
    call globaldea
    call stencils(n, ielem, imaxe, xmpie, xmpielrank, ilocalstencil, typesten, numneighbours, restart)
    if (iadapt .eq. 1) call adapt_criterion
  end if

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finished obtaining the directional stencils across all cpu"
    cpux3(1) = mpi_wtime()
    write(63, *) "time5=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call estabexhange(n, ielem, imaxe, xmpie, xmpin, xmpielrank, ilocalstencil, iexchanger, iexchanges, irecexr, irecexs, numneighbours, ischeme, isize, iperiodicity, typesten, xmpil)
  call renumber_neighbours(n, ielem, xmpie, xmpielrank, iexchanger, iexchanges)
  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_5", cpux3(1) - cpux1(1)
  end if

  call deallocatempi1(n)

  if (nprobes .gt. 0) call probepos(n, probei)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "started prestoring reconstruction matrices"
    cpux3(1) = mpi_wtime()
    write(63, *) "time6=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_6", cpux3(1) - cpux1(1)
  end if

  call local_reconallocation3(n, ilocal_recon3)
  call mpi_barrier(mpi_comm_world, ierror)
  call exch_cords(n)
  call mpi_barrier(mpi_comm_world, ierror)

  if ((rungekutta .ge. 10) .and. (rungekutta .lt. 12))
    call exch_cords2(n, isize, iexboundhiri, iexboundhisi, itestcase, numberofpoints2, iexchanger, iexchanges)

  if ((fastest .ne. 1) .and. (ischeme .ge. 2))
    call allocate_basis_function(n, integ_basis, xmpielrank, idegfree)

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write (100 + n, *) "timei_7", cpux3(1) - cpux1(1)
  end if

  if (stencil_io .eq. 1) call stenprint(n)

  call solex_alloc(n)

  if (rungekutta .ge. 2) then
    if (dimensiona .eq. 3) then
      call direct_side(n)
    else
      call direct_side2d(n)
    end if
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  call exch_cord3(n)

  if (iperiodicity .eq. 1)
    call read_input_period(n, xmpielrank, xmpinrank, xmpie, xmpin, ielem, inode, imaxn, imaxe, ibound, imaxb, xmpinnumber, scaler)

  deallocate(inoder2)
  call mpi_barrier(mpi_comm_world, ierror)

  if ((fastest .ne. 1) .and. (ischeme .ge. 2))
    call walls_higher(n)

  call mpi_barrier(mpi_comm_world, ierror)
  cpux2(1) = mpi_wtime()

  if ((dg .eq. 1) .or. (adda_type .eq. 2))
    call allocate_dg

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write (100 + n, *) "timei_8", cpux3(1) - cpux1(1)
  end if

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "started prestoring reconstruction matrices"
    cpux3(1) = mpi_wtime()
    write(63, *) "time7=", cpux3(1) - cpux1(1)
    close(63)
  end if

  if ((fastest .ne. 1) .and. (ischeme .ge. 2))
    call prestore_1(n)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "finishded prestoring reconstruction matrices"
    cpux3(1) = mpi_wtime()
    write(63, *) "time8=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call dealcordinates2
  call localsdeallocation(n, xmpielrank, ilocalstencil, ilocalstencilper, typesten, numneighbours)
  call deallocatempi2(n)
  call dealcordinates1(n, iexcordr, iexcords)

  call mpi_barrier(mpi_comm_world, ierror)

  if (turbulence .eq. 1) then
    if (dimensiona .eq. 3) then
      call walldistance(n, ielem, imaxe, xmpielrank)
    else
      call walldistance2d(n, ielem, imaxe, xmpielrank)
    end if
  end if

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "started prestoring geometry information"
    cpux3(1) = mpi_wtime()
    write(63, *) "time9=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call mpi_barrier(mpi_comm_world, ierror)

  if (n .eq. 0) then
    cpux3(1) = mpi_wtime()
    write(100 + n, *) "timei_9", cpux3(1) - cpux1(1)
  end if

  call grads_assign(n)
  call find_angles(n)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "allocating solution  and flux variables"
    cpux3(1) = mpi_wtime()
    write(63, *) "time10=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call u_c_allocation(n, xmpielrank, u_c, u_e, itestcase, u_ct)

  if (dg .eq. 1) call build_mass_matrix(n)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "initialising"
    cpux3(1) = mpi_wtime()
    write(63, *) "time10=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  call initialise(n)

  if (restart .gt. 0) call rest_read(n)

  call globaldea2(xmpil, xmpie)

  if (n .eq. 0) then
    open(63, file='history.txt', form='formatted', status='old', action='write', position='append')
    write(63, *) "flux allocation"
    cpux3(1) = mpi_wtime()
    write(63, *) "time10=", cpux3(1) - cpux1(1)
    close(63)
  end if

  call sumflux_allocation(n)

  if (rungekutta .ge. 10) call impallocate(n)

  call mpi_barrier(mpi_comm_world, ierror)
  cpux6(1) = mpi_wtime()
  call timers(n, cpux1, cpux2, cpux3, cpux4, cpux5, cpux6, timex1, timex2, timex3, timex4, timex5, timex6)

  if (statistics .eq. 1) then
    if (n .eq. 0) then
      write(st_n_cpu, fmt='(i10)') isize
      thread_n = omp_get_max_threads()
      write(st_n_threads, fmt='(i10)') thread_n
      statfile = "stats_mpi_"//trim(adjustl(st_n_cpu))//"_threads_"//trim(adjustl(st_n_threads))//".txt"
      open(133, file=statfile, form='formatted', status='replace', action='write')
      write(133,'(5X,A2,1X,A11,A11,A11,A11,A11,A11,A11,A11,A11,A11)')"it","t_time","t_comm","t_comp","t_dgint","t_halo","t_recon","t_bound","t_adda","t_flux","t_update"
      close (133)
    end if
  end if

  if (fastest_q .eq. 1) call memory_fast(n)

  call new_arrays(n)
  call mpi_barrier(mpi_comm_world, ierror)
  call exch_cords_opt(n)

  if (dimensiona .eq. 3) then
    call local_reconallocation4(n)
  else
    call local_reconallocation42d(n)
  end if

  if (fastest .ne. 1) call exchange_higher_pre(n)

  call mpi_barrier(mpi_comm_world, ierror)
  call fix_nodes_local
  call mpi_barrier(mpi_comm_world, ierror)
  call specify_write_variables(n)
  call mpi_barrier(mpi_comm_world, ierror)

  select case (tecplot)
    case (5)
      call partition_preparation(n)
      if (outsurf .eq. 1) then
        call prepare_surfaces_v(n)
        call partition_preparation_wallv(n)
      end if
    case (6)
      call partition_preparation_p(n)
      if (outsurf .eq. 1)
        call partition_preparation_p_wall(n)
  end select

  call mpi_barrier(mpi_comm_world, ierror)
  cpux3(1) = mpi_wtime()
  if (n .eq. 0) write(100 + n, *) cpux3(1) - cpux2(1)

  cpux2(1) = mpi_wtime()

  if (n .eq. 0) print *, "ucns3d running"

  call local_reconallocation5(n)

  if (dimensiona .eq. 3) then
    call time_marching(n)
  else
    call time_marching2(n)
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  cpux3(1) = mpi_wtime()

  if (n .eq. 0)
    write(100 + n, *) "total time taken=", cpux3(1) - cpux2(1), "seconds"

  call mpi_barrier(mpi_comm_world, ierror)
  call mpi_finalize(ierror)

  if (n .eq. 0)
    print *, "ucns3d running running"

end program ucns3d