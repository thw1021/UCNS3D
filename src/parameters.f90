module parameters
  use mpi_info
  use declaration
  implicit none

contains

  subroutine read_ucns3d
    implicit none

    integer :: inv, ix, ibleed
    integer :: inv1
    real    :: angledum
    character(48) :: stamp1, frame
    logical :: here1, here2, here3, here5, here, here4, here7, here8, here9, bleedio

    call mpi_barrier(mpi_comm_world, ierror)
    inquire (file='restartav.dat', exist=here1)
    if (here1) then
      average_restart = 1
    else
      average_restart = 0
    end if
    source_active = 0

    movement = 0

    frame = 'rotframe.dat'
    inquire (file=frame, exist=here)
    if (here) then
      open (16, file=frame, form='formatted', status='old', action='read')
      read (16, *)
      read (16, *)
      read (16, *)
      read (16, *)
      read (16, *) rframe
      read (16, *)
      read (16, *) srf_origin(1), srf_origin(2), srf_origin(3)
      read (16, *)
      read (16, *) srf_velocity(1), srf_velocity(2), srf_velocity(3)
      read (16, *)
      read (16, *) per_rot, angle_per, v_ref
      read (16, *)
      read (16, *) nrotors
      allocate(point1_gl(nrotors, 3), point2_gl(nrotors, 3), radius_gl(nrotors), mrf_rot_gl(nrotors))
      do inv1 = 1, nrotors
        read (16, *) point1_gl(inv1, 1), point1_gl(inv1, 2), point1_gl(inv1, 3)
        read (16, *) point2_gl(inv1, 1), point2_gl(inv1, 2), point2_gl(inv1, 3)
        read (16, *) radius_gl(inv1), mrf_rot_gl(inv1)
      end do
      close (16)
      if (per_rot .eq. 1) then
        tol_per = 1.0e-8
        lowmem = 1
        iperiodicity = 1
      end if
      if (rframe .eq. 2) then
        mrf = 1
        srfg = 0
      end if
      if (rframe .eq. 1) then
        srfg = 1
        mrf = 0
      end if
    else
      rframe = 0
      srfg = 0
      mrf = 0
    end if

    if ((mrf .eq. 1) .or. (srfg .eq. 1)) then
      source_active = 1
      kinit_srf = 0.00001
      if (mrf .eq. 1) then
        rot_corr = 1
        d_corr = 1
      end if
      if (srfg .eq. 1) then
        rot_corr = 0
        d_corr = 1
      end if
    end if

    inquire (file='multispecies.dat', exist=here2)
    if (here2) then
      multispecies = 1
      open (14, file='multispecies.dat', form='formatted', status='old', action='read')
      read (14, *)
      read (14, *)
      read (14, *) nof_species
      allocate(gamma_in(1:nof_species), mp_a_in(1:nof_species), mp_r_in(1:nof_species), mp_pinf(1:nof_species))
      read (14, *) gamma_in(1:nof_species)
      read (14, *) mp_a_in(1:nof_species)
      read (14, *) mp_r_in(1:nof_species)
      read (14, *) mp_pinf(1:nof_species)
      close (14)
    else
      multispecies = 0
    end if

    inquire (file='thermal.dat', exist=here9)
    if (here9) then
      thermal = 1
      open (31, file='thermal.dat', form='formatted', status='old', action='read')
      read (31, *) temp_model
      read (31, *) wall_temp
      wall_temp = wall_temp/287
      close (31)
    else
      thermal = 0
    end if

    inquire (file='mood.dat', exist=here3)
    if (here3) then
      mood = 1
      open (17, file='mood.dat', form='formatted', status='old', action='read')
      read (17, *)
      read (17, *)
      read (17, *) mood_mode
      read (17, *) mood_var1, mood_var2
      read (17, *) mood_var3, mood_var4
      if (n .eq. 0) then
        print *, "mood active"
      end if
      close (17)
    else
      mood = 0
    end if

    inquire (file='indicator.dat', exist=here4)
    if (here4) then
      open (18, file='indicator.dat', form='formatted', status='old', action='read')
      read (18, *)
      read (18, *)
      read (18, *) indicator_type
      read (18, *) indicator_par1
      read (18, *) indicator_par2
      read (18, *) indicator_par3
      close (18)
    else
      indicator_type = 11
      indicator_par1 = 10e-12
      indicator_par2 = 10e-10
      indicator_par3 = zero
    end if

    inquire (file='outlets.dat', exist=here9)
    if (here9) then
      open (29, file='outlets.dat', form='formatted', status='old', action='read')
      read (29, *)
      read (29, *) press_outlet
      print *, "i am reading the pressure for the outlet vents from the file"
      close (29)
    end if

    inquire (file='filter.dat', exist=here7)
    if (here7) then
      filtering = 1
      open (19, file='filter.dat', form='formatted', status='old', action='read')
      read (19, *)
      read (19, *)
      read (19, *) filter_type
      read (19, *) fil_alpha
      read (19, *) fil_s
      read (19, *) fil_nc
      close (19)
    else
      filtering = 0
    end if

    inquire (file='adda.dat', exist=here8)
    if (here8) then
      adda = 1
      open (19, file='adda.dat', form='formatted', status='old', action='read')
      read (19, *)
      read (19, *)
      read (19, *) adda_type
      read (19, *) adda_alpha_1, adda_alpha_2
      read (19, *) adda_1_s, adda_2_s
      read (19, *) adda_1, adda_2
      close (19)
    else
      adda = 0
      adda_type = 0
    end if

    call mpi_barrier(mpi_comm_world, ierror)

    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *) dimensiona, statistics, code_profile
    read (15, *)
    read (15, *)
    read (15, *) governingequations, initcond
    read (15, *)
    read (15, *)
    read (15, *) turbulence, icoupleturb, passivescalar
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *) rres, ufreestream, vvel, wvel, pres
    read (15, *)
    read (15, *)
    read (15, *) aoa, vectorx, vectory, vectorz
    read (15, *)
    read (15, *)
    read (15, *) gamma, prandtl, reynolds, charlength
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *) spatiladiscret, iriemann, spatialorder, limiter, poly
    read (15, *)
    read (15, *)
    read (15, *) wenocnschar, ees, wenoz, wenocentralweight
    read (15, *)
    read (15, *)
    read (15, *) temporder, cfl, timestep, upperlimit, reslimit
    read (15, *)
    read (15, *)
    read (15, *) iboundary, boundtype, scaler
    read (15, *)
    read (15, *)
    read (15, *) greengo, lmach
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *) out_time, ntmax, wallc
    read (15, *)
    read (15, *)
    if (code_profile .lt. 0) then
      read (15, *) tecplot, output_freq, ievery2, ieveryav, stencil_io
      ievery = 10e15
    else
      read (15, *) tecplot, ievery, ievery2, ieveryav, stencil_io
      output_freq = 10e15
    end if
    read (15, *)
    read (15, *)
    read (15, *) averaging
    read (15, *)
    read (15, *)
    read (15, *) outsurf, iforce, surfshear
    read (15, *)
    read (15, *)
    read (15, *) ires_turb, ires_unsteady, lamps, prev_turbmodel
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *)
    read (15, *) nprobes
    read (15, *)

    fastmovie = 0
    dg = 0

    ispal = 1
    cb1 = 0.1355
    cb2 = 0.622
    sigma = 0.666666667
    kappa = 0.41
    cw1 = 3.23886781677
    cw2 = 0.3
    cw3 = 2.0
    cv1 = 7.1
    ct1 = 1.0
    ct2 = 2.0
    ct3 = 1.1
    ct4 = 2.0
    prtu = 0.9
    twall = 0
    turbinit = 3.0
    upturblimit = 1000000
    residualfreq = 10
    irs = 0
    c_des_sa = 0.61

    vort_model = 0
    qsas_model = 0
    zero_turb_init = 0
    sigma_k1 = 1.176470588
    sigma_k2 = 1.0
    sigma_om1 = 2.0
    sigma_om2 = 1.168
    aa_1 = 0.31
    beta_i1 = 0.075
    beta_i2 = 0.0828
    alpha_starinf = 1.0
    alpha_0 = 0.111111111
    beta_starinf = 0.09
    r_beta = 8.0
    r_k_sst = 6.0
    beta_t = 0.072
    kappa_sst = 0.41
    r_om_sst = 2.95
    zeta_star = 1.5
    m_t0 = 0.25
    c_mu_inlet = 0.09
    c_smg = 0.11
    eta2_sas = 3.51
    sigma_phi = 0.666666667
    c_sas = 2.0
    l_turb_inlet = 0.0026
    i_turb_inlet = 0.1
    init_mu_ratio = 0.01
    c_des_sst = 0.61

    schmidt_lam = 10.0
    schmidt_turb = 0.7

    cavitation = 0
    extended_bounds = 0

    select case (code_profile)
    case (0)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 30
      gridar2 = 30
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (501)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 20.0
      gridar2 = 50.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (30)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 20.0
      gridar2 = 50.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (8)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 100.0
      gridar2 = 200.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (88)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 40.0
      gridar2 = 40.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (18)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 1
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 40000.0
      gridar2 = 40000.0
      fastest = 0
      lmach_style = 0
      lamx = 0.0d0
      lamy = 0.0d0
      lamz = 0.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (98)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 1
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 20.0
      gridar2 = 40.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (9)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 7.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (91)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 1
      guassianquadra = 0
      fastest_q = 1
      relax = 3
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 7.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (1)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 7.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      turbinit = 0.1
      des_model = 2

    case (11)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 10.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 2

    case (100)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 3
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 10.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0
      dg = 1
      bound_lim = 0
      cavitation = 0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 2
      br2_damping = 3.0
      br2_yn = 2

    case (101)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 0
      guassianquadra = 6
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 10.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0
      dg = 1
      bound_lim = 0
      cavitation = 0

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 2
      br2_damping = 3.0
      br2_yn = 2

    case (102)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 10.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0
      dg = 1

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 2
      br2_damping = 3.0
      br2_yn = 2

    case default
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      icoupleturb = 0
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 0
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 10.0
      gridar2 = 10.0
      fastest = 0
      lmach_style = 0
      lamx = 1.0d0
      lamy = 1.0d0
      lamz = 1.0d0
      ispal = 1

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    case (17)
      lowmemory = 0
      binio = 1
      lowmem = 0
      reduce_comp = 0
      turbulencemodel = 1
      ihybrid = 0
      hybridist = 0.0d0
      swirl = 0
      iadapt = 0
      icompact = 0
      extf = 2
      weight_lsqr = 1
      guassianquadra = 0
      fastest_q = 1
      relax = 1
      cflmax = 30
      cflramp = 0
      emetis = 6
      itold = 10000
      gridar1 = 200000.0
      gridar2 = 500000.0
      fastest = 0
      lmach_style = 0
      lamx = 0.0d0
      lamy = 0.0d0
      lamz = 0.0d0
      iweno = 5

      if (iboundary .eq. 1) then
        lowmem = 1
      end if

      des_model = 0

    end select

    if (n .eq. 0) then
      open (63, file='history.txt', form='formatted', action='write', position='append')
      call date_and_time(stamp1)
      write (63, *) "ucns3d parallel computation start", stamp1
    end if

    nproc = isize - 1
    vorder = iorder

    !$omp master
    if (nprobes .gt. 0) then
      if (dimensiona .eq. 3) then
        allocate(probec(1:nprobes, 1:3))
        do inv = 1, nprobes
          read (15, *) probec(inv, 1), probec(inv, 2), probec(inv, 3)
        end do
      else
        allocate(probec(1:nprobes, 1:2))
        do inv = 1, nprobes
          read (15, *) probec(inv, 1), probec(inv, 2)
        end do
      end if
    end if
    !$omp end master

    call mpi_barrier(mpi_comm_world, ierror)

    if (governingequations .le. 2) then
      if (dimensiona .eq. 3) then
        nof_variables = 5
        dims = 3
      else
        nof_variables = 4
        dims = 2
      end if
    else
      nof_variables = 1
      if (dimensiona .eq. 3) then
        dims = 3
      else
        dims = 2
      end if
    end if

    if (governingequations .eq. -1) then
      if (dimensiona .eq. 3) then
        nof_variables = 5 + nof_species + (nof_species - 1)
        dims = 3
      else
        nof_variables = 4 + nof_species + (nof_species - 1)
        dims = 2
      end if
    end if

    if (turbulence .eq. 1) then
      if (turbulencemodel .eq. 1) then
        turbulenceequations = 1
      else if (turbulencemodel .eq. 2) then
        turbulenceequations = 2
      end if
    else
      turbulenceequations = 0
      turbulencemodel = 0
    end if

    select case (governingequations)
    case (-1)
      itestcase = 3
      ivortex = 0

    case (1)
      itestcase = 4
      ivortex = 1

    case (2)
      itestcase = 3
      ivortex = 0

    case (3)
      itestcase = 1
      ivortex = 0

    case (4)
      itestcase = 2
      ivortex = 0

    end select

    unwou = 3
    betaas = 1.5d0
    suther = 0.412158681d0
    uvel = ufreestream

    if (initcond .eq. 977) then
      uvel = 0.0d0
      vvel = 0.0d0
      ufreestream = wvel
    end if

    if (pres .lt. 0) pres = rres/gamma

    if (rframe .eq. 0) then
      visc = (rres*ufreestream*charlength)/reynolds
    else
      visc = (rres*v_ref*charlength)/reynolds
    end if

    if (swirl .eq. 1) uvel = zero
    if (aoa .ne. 0.0d0) then
      angledum = (aoa*pi)/180.0d0
      uvel = cos(angledum)*ufreestream*vectorx
      wvel = sin(angledum)*ufreestream*vectorz
      vvel = sin(angledum)*ufreestream*vectory
    end if

    iorder = max(1, spatialorder - 1)
    firstorder = 0
    iextend = 2
    if (icompact .eq. 1) iextend = 10
    iweno = 0
    wenwrt = 1
    ischeme = 3
    typesten = 1
    lwci1 = wenocentralweight

    if (fastest .eq. 1) ischeme = 1

    select case (spatiladiscret)
    case (1)
      if (spatialorder .eq. 1) then
        firstorder = 1
        iorder = 1
      else
        firstorder = 0
      end if

    case (2)
      iweno = -1

    case (3)
      iweno = 1
      iextend = 12
      if (ees .eq. 5) then
        if (icompact .eq. 1) then
          iextend = 10
        else
          iextend = 5
        end if
      end if
      wenwrt = wenocnschar

      if (dimensiona .eq. 2) then
        if ((ees .eq. 0)) typesten = 5
        if ((ees .eq. 5)) typesten = 5
        if (ees .eq. 1) typesten = 9
        if (ees .eq. 2) typesten = 9
        if (ees .eq. 3) typesten = 5
      else
        if ((ees .eq. 0)) typesten = 7
        if (ees .ge. 4) typesten = 7
        if (ees .eq. 1) typesten = 15
        if (ees .eq. 2) typesten = 15
        if (ees .eq. 3) typesten = 7
      end if

    end select

    if (code_profile .eq. 17) iweno = 5

    rungekutta = temporder

    if (iboundary .eq. 0) then
      iperiodicity = -3
    else if (iboundary .eq. 1) then
      iperiodicity = 1
    end if

    if (dg .eq. 1) then
      igqrules = min(iorder + 1, 6)
    else
      igqrules = min(iorder, 6)
    end if

    alls = igqrules**(dimensiona)

    call mpi_barrier(mpi_comm_world, ierror)

    if (n .eq. 0) then
      open (63, file='history.txt', form='formatted', action='write', position='append')
      write (63, *) 'total number of processes:', isize, igqrules
      close (63)
    end if

    if (srfg .eq. 1) then
      if (n .eq. 0) then
        open (63, file='history.txt', form='formatted', action='write', position='append')
        write (63, *) 'single reference frame engaged:'
        close (63)
      end if
    end if

    if (mrf .eq. 1) then
      if (n .eq. 0) then
        open (63, file='history.txt', form='formatted', action='write', position='append')
        write (63, *) 'number of rotating frames engaged:', nrotors
        close (63)
      end if
    end if

    if (per_rot .eq. 1) then
      iboundary = 1
      lowmem = 1
      if (n .eq. 0) then
        open (63, file='history.txt', form='formatted', action='write', position='append')
        write (63, *) 'rotational  periodicity engaged'
        close (63)
      end if
    end if

    if (initcond .eq. 444) then
      inquire (file='bubbles.dat', exist=here5)
      if (here5) then
        open (18, file='bubbles.dat', form='formatted', status='old', action='read')
        read (18, *) nof_bubbles
        do ix = 1, nof_bubbles
          read (18, *) bubble_centre(ix, 1), bubble_centre(ix, 2), bubble_radius(ix)
        end do
        close (18)
      end if
    end if

    if (initcond .eq. 445) then
      inquire (file='bubbles.dat', exist=here5)
      if (here5) then
        open (18, file='bubbles.dat', form='formatted', status='old', action='read')
        read (18, *) nof_bubbles
        do ix = 1, nof_bubbles
          read (18, *) bubble_centre(ix, 1), bubble_centre(ix, 2), bubble_centre(ix, 3), bubble_radius(ix)
        end do
        close (18)
      end if
    end if

    inquire (file='bleed.dat', exist=bleedio)
    if (bleedio) then
      bleed = 1
      open (27, file='bleed.dat', form='formatted', status='old', action='read')
      read (27, *)
      read (27, *)
      read (27, *) bleed_number
      read (27, *) bleed_type
      allocate(bleed_start(1:bleed_number,1:dims),bleed_end(1:bleed_number,1:dims),bleed_plenum(1:bleed_number), bleed_porosity(1:bleed_number))
      do ibleed = 1, bleed_number
        read (27, *) bleed_start(ibleed, 1:dims), bleed_end(ibleed, 1:dims)
      end do
      do ibleed = 1, bleed_number
        read (27, *) bleed_plenum(ibleed), bleed_porosity(ibleed)
      end do
      close (27)
    else
      bleed = 0
    end if

  end subroutine read_ucns3d

end module parameters