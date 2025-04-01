module declaration
  implicit none

  integer::cascade,mood,mood_mode,kmaxn,multispecies,nof_species,dimensiona,lowmem,binio,nof_variables,chunk_n,dims,ires_turb,code_profile,ires_unsteady,lamps,itotalb,totiw,icompact,ees,iscoun,itold,lmach_style,weight_lsqr        !dimensions of problem
  integer::governingequations,filtering,guassianquadra,temporder,iboundary,wenocnschar,required,nodes_i,swirl,iadapt,tecplot,stencil_io,surfshear,issf,n_boundaries,fastest_q,statistics,adda,alls
  integer, allocatable, dimension(:, :)::jtot1, jtot2, jtot3, jtot, el_connect
  integer::jtotal, jtotal1, jtotal2, jtotal3, fastmovie, movement, typ_countn_global, typ_countn, typ_countn_global_w, typ_countn_w
  integer::filter_type, fil_nc, fil_s, fil_alpha                !filter values
  integer::adda_type, adda_alpha_1, adda_alpha_2, adda_1_s, adda_2_s, adda_1, adda_2                !filter values
  integer::kill                         ! flag for killing a simulation
  integer::nderivative                  ! index of the numbering of the component of the polynomials
  real::output_freq                     !output frequency in simulation time
  integer::extended_bounds              !bounds strict or relaxed
  integer::iconsr                       ! index for identifying the considered element
  integer::inwhichel                    ! index for identifying the considered element within a stencil
  integer::subdiv                       ! index for the number of subdivisions of a face
  integer::iconsgvq                     ! index for considered cell
  integer::igianagraps                  !index for considered cell
  integer::it                           ! number of iterations
  integer::iconimp                      !index for considered cell
  integer::ibcode                       ! index for boundary condition code
  integer::ibside                       ! index for which side of cell is bounded
  integer::cfw                          !index for determining from the which section the boundary subroutine is called
  integer::dg, br2_yn                   ! flag for dg discretisation
  real:: br2_damping
  integer::lowmemory                    ! memory usage flag
  integer::fastest                      ! fastest mode flag for the code
  integer::emetis                       ! type of metis partitioning
  integer::allnodesgloball              ! total number of nodes in the domain
  integer::spatiladiscret               ! spatial discretisation
  integer::spatialorder                 ! spatial order
  integer::jk                           !index
  integer::numneighbours, numneighbours2          ! number of neighbours in the stencils
  integer::restart                                ! flag for determining if a restart file is present
  integer::ihybrid                                !hybrid mode where some part of the domain is solved with lower order and another with higher order schemes
  integer::upperlimit !
  integer::residualfreq                                !frequency for writing the residuals
  integer::irs                                        !implicit residual smoothing version
  integer::greengo                                !gradient approximation techniques
  integer::lmach                                        !low mach preconditioning
  integer::ispal                                        !type of spalart allmaras modification
  integer::relax                                        !relaxation schemes gauss seidel, jacobian etc
  integer::icong                                        !index for considered cell
  integer::iloop                                        !index for considered cell
  integer::cflramp                                !cfl ramping activated for implicit time stepping
  integer::averaging                                !flag for activating averaging of transient data
  integer::average_restart                        !flag for having a restart file with the averaged solution
  integer::prev_turbmodel                                !flag of previous turbulence model in case of a restart file
  integer::iforce                                        !index on how often to compute the forces
  integer::boundtype                                !type of boundary condition supersonic,  subsonic
  integer::nproc                                        !number of cpus
  integer::initcond                                !initial condition type
  integer::wenwrt                                        !weno wrt to which variable conserved, characteristics etc
  integer::poly                                        !interpolating polynomials 1 generic, 2 legendre
  integer::ivortex                                !q criterion com putation
  integer::icarlos1                                !index for considered element passive scalars, complicated inflow condition
  integer::icarlos2                                !index for considered face passive scalars, complicated inflow condition
  integer::imaxe                                        !total number of elements of the grid
  integer::ntmax                                        !maximum number of iterations
  integer::ihax1                                        !temporary integer used as pointer for operations regarding hexahedral elements
  integer::limiter                                !slope limiter
  integer::iextend                                !extend parameter for central big stencil only for weno cases
  integer::idegfree, idegfree2, inum2, idegfree3                !degrees of freedom of polynomial it will become loca        l
  integer::outsurf                                !outsurf stands for computing forces
  integer::unwou                                        !type of output files to write all partitions to one file etc
  integer::igqrules                                !gaussian quadrature rule
  integer::imaxdegfree, imaxdegfree2                                !maximum number of neighbours
  integer::qp_hexa, qp_tetra, qp_pyra, qp_prism, qp_quad, qp_triangle, qp_line, qp_line_n, reduce_comp, qp_quad_n, qp_triangle_n !volume gaussian quadrature points
!> number of points allocated to qpoints, wequa (coords and weights for quadrature points)
  integer::numberofpoints, numberofpoints2
  integer::rescounter                                !counter for residuals
  integer::rescountert                                !counter for turbulence residuals
  integer::ilx                                        !degrees of freedom for required order of accuracy
  integer::imaxb                                        !total number of boundary conditions
  integer::imaxn                                        !total number of nodes present in the grid
  integer::in                                        !index for considered cell
  integer::iorder, iorder2                                        !order of accuracy (spatial)
  integer::ioverst                                ! number of elements for the stencil overdetermined matrix required for avoiding ill conditions
  integer::ioverto                                ! number of elements for the estab stencils overdetermined matrix required for avoiding ill conditions
  integer::qrde                                        !type of solution for the least squares problem
  integer::rungekutta                             !time stepping scheme type
  integer::ischeme                                 !type of scheme to be used
  integer::isplit                                        !unsplit finit volume scheme
  integer::iselem                                        !maximum number of elements required for the stencil dependent on the order of accuracy
  integer::istn                                        !variable that takes each time the value of the considered element  to find neighbours in the stencil
  integer::itt, wenoz                                        !integer variable used for opening files
  integer::itestcase              !types of equations to be solved
  integer::iweightlsqr                            !weighted least squres option of inverted distance
  integer::typesten                                !number of stencils for weno
  integer::iperiodicity                                !periodicity detected in domain
  integer::iweno                                        !type of schemes weno, muscl, unlimited
  integer::turbulence                                !turbulence equations
  integer::turbulencemodel                        !which model
  integer::turbulenceequations                        !how many equations to solve
  integer::icoupleturb                                !coupled turbulence model
  integer::iriemann                                !riemann solver
  integer::firstorder                                !first order scheme active
  integer::passivescalar                                !additional transport equations
  integer::qsas_model                                !qsas model
  integer::vort_model                                !vorticity model
  integer::zero_turb_init                                !flag for type of initialization for k-omega
  integer::des_model                                 !integer switches for sst
  integer:: nprobes, totwalls, nof_interior, nof_bounded, mrf                        !number of probes for transient data
  integer:: rot_corr, d_corr   !integer for turbulence corrections
  !--------------------- variables for parallel partitioned output-------!
  integer, allocatable, dimension(:)::dispart1, dispart2, dispart3, dispart4, dispart5, typ_nodesn, typ_nodesn_w
  integer, allocatable, dimension(:)::iarray_part1, iarray_part2, iarray_part3, iarray_part4, iarray_part5, i_array_part2x
  real, allocatable, dimension(:)::rarray_part4, rarray_part2, rarray_part3, rarray_part5
  integer, allocatable, dimension(:)::wdispart1, wdispart2, wdispart3, wdispart4, wdispart5, wallcount_cpu_l, wallcount_cpu_g
  integer, allocatable, dimension(:)::wiarray_part1, wiarray_part2, wiarray_part3, wiarray_part4, wiarray_part5
  integer, allocatable, dimension(:)::offset_vtu, connect_vtu, type_vtu, offset_vtu_w, connect_vtu_w, type_vtu_w
  real, allocatable, dimension(:, :)::sol_vtu, sol_vtu_w
  real, allocatable, dimension(:)::nodes_vtu, nodes_vtu_w
  real, allocatable, dimension(:)::wrarray_part4, wrarray_part2, wrarray_part3, wrarray_part5
  real, allocatable, dimension(:, :)::rarray_part1, wrarray_part1
  integer::part1_end, part2_end, part3_end, part4_end, part5_end
  integer::wpart1_end, wpart2_end, wpart3_end, wpart4_end, wpart5_end
  integer::kdum1, kdum2, write_variables, write_variables_av, nodes_part
  integer::wkdum1, wkdum2, write_variables_w, write_variables_av_w, wnodes_part
  integer::datatypex, datatypey, datatypez, datatypexx, datatypeyy, datatypeint
  integer, dimension(1)::kdum3
  character(len=25)::variable_names(15), variable_names_av(15)
  integer::wdatatypex, wdatatypey, wdatatypez, wdatatypexx, wdatatypeyy, wdatatypeint
  integer, dimension(1)::wkdum3
  character(len=25)::variable_names_w(15), variable_names_av_w(15)
  integer::kloopx, iloopx, totwallsc, iwmaxe
  integer,allocatable,dimension(:)::wallit,offsetwall,wall_nodes,wallcx,offsetwc,offsetwc_g,wallcx_g,wallshape,wallshape_g,wallshape_g2
  integer, allocatable, dimension(:, :)::wall_l
  !--------------------- end of variables for parallel partitioned output-------!
  real, dimension(3)::srf_origin, srf_velocity
  real::press_outlet
  integer::rframe, source_active
  real::per_rot, angle_per, v_ref, kinit_srf, srfg, tol_per
  real, allocatable, dimension(:, :)::point1_gl, point2_gl
  real, allocatable, dimension(:)::radius_gl, mrf_rot_gl
  integer:: nrotors
  integer::num_dg_dofs !number of degrees of freedom for dg polynomial approximation
  integer::num_dg_reconstruct_dofs !number of degrees of freedom for dg order + 1, not including num_dg_dofs
  integer::indicator_type                !troubled indicator type
  integer::viscous_s, jump_cond1, jump_cond2, jump_cond3
  integer::cavitation
  real::indicator_par1, indicator_par2, indicator_par3, bound_lim  !troubled indicator parameters
  real::rhc1, rhc2, rhc3, rhc4
  real::prace_t1,prace_t2,prace_t3,prace_t4,prace_t5,prace_t6,prace_t7,prace_t8,prace_t9,pr_t1,pr_t2,pr_t3,pr_t4,pr_t5,pr_t6,pr_t7,pr_t8,prace_tx1,prace_tx2,prace_tx3
  !------------------start bleed parameters-------------------!
  integer::bleed_number, bleed, bleed_type
  real, allocatable, dimension(:, :)::bleed_start, bleed_end
  real, allocatable, dimension(:)::bleed_plenum, bleed_porosity
  !------------------end bleed parameters-------------------!

  integer, allocatable, dimension(:, :)::probei        !probe positions
  integer, allocatable, dimension(:)::nodes_offset_local, nodes_offset_local2, nodes_offset_localw, nodes_offset_local2w        !offsets for the nodes for each cell in my cpu from the global list
  integer, allocatable, dimension(:)::nodes_offset, nodes_offset2, nodes_offsetw, nodes_offsetw2        !offsets for the nodes for each cell  globally
  integer, allocatable, dimension(:, :)::glneigh        !globalist of all the neighbours in the grid
  integer, allocatable, dimension(:, :)::glneighper        !globalist of all the rotational periodic neighbours in the grid
  integer, allocatable, dimension(:)::xsize, xgo, xgo_v, chunk_size        !allocatable index for describing which stencil is considered within the directional stencil construction routine
  integer, allocatable, dimension(:)::my_nodesl, my_nodesg        !allocatable arrays for local and global indexing of nodes belonging to each cpu (no duplications)

  integer, allocatable, dimension(:)::ibound_t, ibound_t2        !allocatable index for describing the shyape of the considered cell
  integer, allocatable, dimension(:)::iselemt        !total number of elements in each stencil (extended one)
  integer, allocatable, dimension(:)::xmpiall, offset, xmpiwall, woffset, xmpi_re, xmpi_wre, offset_v, xmpiall_v                !allocatable index for stencil routine

  integer, allocatable, dimension(:)::el_int, el_bnd
  integer, allocatable, dimension(:)::pare, dose, pareel, pares
  integer, allocatable, dimension(:, :)::doseel, soseel
  integer, allocatable, dimension(:, :, :, :)::ilocalstencil        !4-d array for stencils
  integer, allocatable, dimension(:, :, :, :)::ilocalstencilper        !4-d array for rotational periodic stencils
  integer, allocatable, dimension(:, :, :, :)::ilocalallelg        !4-d array for stencils
  integer, allocatable, dimension(:, :, :, :)::ilocalallelgper            !4-d array for rotational periodic stencils
  integer, allocatable, dimension(:)::xmpie                        !global list of number of elements that belong to each cpu
  integer, allocatable, dimension(:)::xmpil                        !local list of number of elements that belong to each cpu
  integer, allocatable, dimension(:)::xmpin                        !not used can get rid of
  integer, allocatable, dimension(:)::xmpielrank                !number of elements in this cpu
  integer, allocatable, dimension(:)::xmpinrank                !not used can get rid of
  integer, allocatable, dimension(:, :, :)::xmpinnumber        !not used can get rid of
  integer, allocatable, dimension(:)::stcon                        !dummy variable for recursive subroutine of stencils
  integer, allocatable, dimension(:)::stconc                !dummy variable for recursive subroutine of stencils
  integer, allocatable, dimension(:)::stcons                !dummy variable for recursive subroutine of stencils
  integer, allocatable, dimension(:)::list, ineb, iperb, nodelist                !dummy variable for recursive subroutine of stencils
  integer, allocatable, dimension(:, :, :, :)::ilocalalls      !dummy variable for recursive subroutine of stencils
  integer::thermal, temp_model

  real::wenocentralweight, timestep, oo2, zero, forcex, forcey, forcez, extf, vorder, mood_var1, mood_var2, mood_var3, mood_var4
  real, allocatable, dimension(:)::sumvars, maxvars, aver_vars            !variables for bounds of troubled cell indicator
  real::pi                                        !pi trigonometri
  real::taylor, taylor_ens, taylor_ensx                                        !only to be used for taylor green vortex
  real::voll                                        !total volume of the domain
  real::wall_temp                                        !wall temperature model
  real::upturblimit                                !upper turbulence viscosity ratio
  real::hybridist                                        !upper turbulence viscosity ratio
  real::reynolds                                        !reynolds number
  real::cflmax                                        !maximum number of allowable cfl to be used only with implicit
  real::prevres                                        !previous residual in order to determine ramping strategy
  real::vectorx                                        !setting up with respect to which plane the aoa is defined x
  real::vectory                                        !setting up with respect to which plane the aoa is defined y
  real::vectorz                                        !setting up with respect to which plane the aoa is defined z
  real::gridar1                                        !aspect ratio of grid
  real::gridar2                                        !maximum volume aspect ratio of stencils
  real::t                                                !time
  real::totalvolume                                !total volume of domain
  real::reslimit                                       !limit to stop the simulation
  real::res_time                                        !restart_time
  real::ievery2                                        !how often to write output
  real::wallc                                        !wall clock limit in seconds
  real::betaas                                        !betaas of sutherland exponent
  real::suther                                        !sutherland non dimensional constant of temperature
  real::prandtl                                        !prandtl number
  real::scaler                                        !for scaling the mesh
  real::lwci1                                        !linear weno weight of central stencils                                !
  real::resmax                                        !residual dimensional mean flow
  real::resmaxt                                        !residual dimensional turbulence
  real::cfl                                        !cfl number
  real::ievery, modeio                                        !how often to write output
  real::ieveryav                                        !how often to write output of average solution
  real:: firstresu                                !residual dimensional mean flow
  real::firstresv                                        !residual dimensional mean flow
  real::firstresw                                        !residual dimensional mean flow
  real::firstrese                                        !residual dimensional mean flow
  real::firstresr                                        !residual dimensional mean flow
  real::firstrest                                        !residual dimensional mean flow
  real::firstresk                                        !residual dimensional mean flow
  real::firstresomega                                !residual dimensional mean flow
  real::firstrespass                                !residual dimensional mean flow
  real::gamma
  real::uvel                                        !u velocity
  real::ufreestream                                !u velocity
  real::vvel                                        !v-velocity
  real::wvel                                        !w-velocity
  real::pres                                        !pressure
  real::rres                                        !density
  real::spos                                        !speed of sound
  real::etot                                        !total energy
  real::spkin                                        !specific kinetic energy
  real::lam                                        !heat conductivity (watt/(mk))
  real::visc                                        !viscosity
  real::lamx                                        !flags to be used for restarting
  real::lamy                                        !flags to be used for restarting
  real::lamz                                        !flags to be used for restarting
  real::alpha, beta                !flags to be used for restarting
  real::out_time                                        !final time to write output for unsteady simulations
  real::every_time, ek_time                                        !final time to write output for unsteady simulations
  real::xper                                        !periodicity in x axis
  real::yper                                        !periodicity in y axis
  real::zper                                        !periodicity in z axis
  real::aoa                                        !angle of attack
  real::charlength                                !characteristic length if undefined, and reynolds undefined the free stream will be used for air
  !turbulence model constants
  real::turbinit                                !turbulence initial value vt/v
  real::cb1
  real::cb2
  real::sigma
  real::kappa
  real::cw1
  real::cw2
  real::cw3
  real::cv1
  real::ct1
  real::ct2
  real::ct3
  real::ct4
  real::prtu
  real::twall
  character(len=30)::statfile, st_n_cpu, st_n_threads
  integer::thread_n
  real::sigma_k1
  real::sigma_k2
  real::sigma_om1
  real::sigma_om2
  real::aa_1
  real::beta_i1
  real::beta_i2
  real::alpha_starinf
  real::alpha_0
  real::beta_starinf
  real::r_beta
  real::r_k_sst
  real::beta_t
  real::kappa_sst
  real::r_om_sst
  real::zeta_star
  real::m_t0
  real::c_smg
  real::eta2_sas
  real::sigma_phi
  real::c_sas
  real::alpha_inf1
  real::alpha_inf2
  real::alpha_star0
  real::c_mu_inlet                                                   !functions for sst
  real::l_turb_inlet
  real::i_turb_inlet
  real::init_mu_ratio
  real::c_des_sa
  real::c_des_sst
  real:: schmidt_lam
  real::schmidt_turb
  real::tolsmall, tolbig, allresdt, tz1
  real, dimension(7)::allres, initialres
  real, dimension(100, 3)::bubble_centre
  real, dimension(100)::bubble_radius
  integer::nof_bubbles
  real:: momentx, momenty, momentz

  real, allocatable, dimension(:)::gamma_in, mp_r_in, mp_a_in, mp_pinf !multiphase components
  real, allocatable, dimension(:)::modal_filter, adda_filter_weak, adda_filter_strong, modal_filter_strong, modal_filter_weak
  real::l1norm                !l1 norm of solution for grid convergence studies of euler and linear advection equations
  real::l2norm                !l2 norm of solution for grid convergence studies of euler and linear advection equations
  real::l0norm, stennorm                !l0 norm of solution for grid convergence studies of euler and linear advection equations
  real::dt                                         !real time step size
  real, allocatable, dimension(:)::avrg                !temporary solution averages
  real, allocatable, dimension(:, :)::jac                !jacobians
  real, allocatable, dimension(:, :)::inversejac        !inverse jacobians
  real, allocatable, dimension(:)::deterjac                !determinant jacobians
  real, allocatable, dimension(:, :)::probec                !probe positions
  real, allocatable, dimension(:, :)::impdu                !implicit only change of solution
  real, allocatable, dimension(:, :, :)::impdiag        !implicit only diagonal matrix d
  real, allocatable, dimension(:, :, :, :)::impoff        !implicit only off diagonal matrix d
  real, allocatable, dimension(:)::impdiag_mf        !implicit only diagonal matrix d
  real, allocatable, dimension(:, :)::impoff_mf        !implicit only off diagonal matrix d
  real, allocatable, dimension(:, :, :)::impofft        !implicit only off diagonal matrix d for turbulence
  real, allocatable, dimension(:, :)::impdiagt        !implicit only  diagonal matrix d for turbulence
  real, allocatable, dimension(:, :)::sht                !implicit only  soruce term jacobian
  real, allocatable, dimension(:)::maxdiff                !se
  real, allocatable, dimension(:)::mindiff                !se
  real, allocatable, dimension(:)::tmaxdiff                !se
  real, allocatable, dimension(:)::tmindiff                !se
  real, allocatable, dimension(:)::xmin                !se
  real, allocatable, dimension(:)::xmax                !se
  real, allocatable, dimension(:)::ymin                !se
  real, allocatable, dimension(:)::ymax                !se
  real, allocatable, dimension(:)::zmin                !se
  real, allocatable, dimension(:)::zmax                !se
  real, allocatable, dimension(:)::cpux1                !timer
  real, allocatable, dimension(:)::cpux2                !timer
  real, allocatable, dimension(:)::cpux3                !timer
  real, allocatable, dimension(:)::cpux4                !timer
  real, allocatable, dimension(:)::cpux5                !timer
  real, allocatable, dimension(:)::cpux6                !timer
  real, allocatable, dimension(:)::cpux7                !timer
  real, allocatable, dimension(:)::timex1                !timer
  real, allocatable, dimension(:)::timex2                !timer
  real, allocatable, dimension(:)::timex3                !timer
  real, allocatable, dimension(:)::timex4                !timer
  real, allocatable, dimension(:)::timex5                !timer
  real, allocatable, dimension(:)::timex6                !timer
  real, allocatable, dimension(:, :)::ifin                !temporary pointer for sync output
  real, allocatable, dimension(:, :)::tfin                !temporary pointer for sync output
  real, allocatable, dimension(:, :)::centerr        !for directional stencils
  integer, allocatable, dimension(:)::cand, cands, candr
  integer, allocatable, dimension(:, :)::candxr, candxs, cand2s, cand2rt
  integer, allocatable, dimension(:, :, :)::cand2r
  real, allocatable, dimension(:, :)::xand2s, xand2rt
  real, allocatable, dimension(:, :, :)::xand2r
  real, allocatable, dimension(:)::flux_term_left_z, flux_term_left_x, flux_term_left_y
  real, allocatable, dimension(:)::flux_term_right_z, flux_term_right_x, flux_term_right_y
  real, allocatable, dimension(:, :)::sind1, sind2, sind3, sind4, sind5, sind6

  type::vol_gqp
    real, allocatable, dimension(:)::x, y, z, qp_weight
  end type vol_gqp

  type(vol_gqp), allocatable, dimension(:)::qp_array !cell num, qp num

  type::u_centre                                        !data type that holds the cell center values for the equations to be solved
    real, allocatable, dimension(:, :)::val        !actual values (first index corresponds to runge-kutta stage,the second one corresponds to number of variables for mean flow equations)
    real, allocatable, dimension(:, :, :)::valdg !actual values (rk stage, variable number, cell-averaged solution variables and expansion coefficients)
    real, allocatable, dimension(:)::rms        !rms values of the conserved vector in case of transient simulations
    real, allocatable, dimension(:, :, :)::br2_aux_var ! (num_dg_dofs, nof_variables, dimensiona)
  end type u_centre
  type::mass_matrix                                        !data type that holds the dg_mass_matrix
    real, allocatable, dimension(:, :)::val
  end type mass_matrix
  type::face_d
    integer, allocatable, dimension(:)::q_mapl, qmapd
  end type face_d
  type(u_centre), allocatable, dimension(:)::u_c         !1-d array for type for solution of mean flow equations
  type(u_centre), allocatable, dimension(:)::u_cx, u_cs, u_cw         !1-d array for type for solution of mean flow equations
  type(u_centre), allocatable, dimension(:)::u_ct         !1-d array for type for solution of turbulent equations
  type(mass_matrix), allocatable, dimension(:)::m_1         !1-d array for type for solution of turbulent equations
  type::u_exact                                        !data type that holds the cell center values of the exact solution for linear advection of euler equations if an exact solution exists
    real, allocatable, dimension(:, :)::val        !values of the exact solution
    real, allocatable, dimension(:, :, :)::valdg
  end type u_exact
  type(u_exact), allocatable, dimension(:)::u_e        !1-d array for type for exact solution of mean flow equations
  type::sumfluxes                                        !data type for the fluxes contribution
    real, allocatable, dimension(:)::val        !number of values for the fluxes
    real, allocatable, dimension(:, :)::valdg        !number of values for the fluxes (dofs, variable)
    real, allocatable, dimension(:, :)::sol_mm_dg!number of values for the fluxes (dofs, variable)
  end type sumfluxes
  type(sumfluxes), allocatable, dimension(:)::sumconflux        !convective fluxes mean flow equations
  type(sumfluxes), allocatable, dimension(:)::sumviscflux        !viscous fluxes mean flow equations
  type(sumfluxes), allocatable, dimension(:)::sourceterm        !sourceterm
  type(sumfluxes), allocatable, dimension(:)::sumconfluxt        !convective fluxes additional equations
  type(sumfluxes), allocatable, dimension(:)::sumviscfluxt  !viscous fluxes additional equations
  type(sumfluxes), allocatable, dimension(:)::rhs                !right hand side of the mean flow equations (sum of all fluxes contributions)
  type(sumfluxes), allocatable, dimension(:)::rhst                !right hand side of the additional flow equations (sum of all fluxes contributions)

  type integralbasis !integral basis functions
    real, allocatable, dimension(:)::value, valuec !values for the high and lower-order polynomials respectively
  end type integralbasis
  type(integralbasis), allocatable, dimension(:)::integ_basis, integ_basis_dg

  type local_recon3
    integer::local, mrf
    real, allocatable, dimension(:)::cond  !dummy variable used for gradient approximation estimation
    integer, allocatable, dimension(:, :)::ihexg, ihexgc !global index of cells
    integer, allocatable, dimension(:, :)::ihexl, ihexlc !local index of cells
    integer, allocatable, dimension(:, :)::ihexb, ihexbc !cpu that that each cell belongs to
    integer, allocatable, dimension(:, :)::ihexn, ihexnc !internal index from where to take the values from communicated messages
    real, allocatable, dimension(:, :, :)::gradients, gradientsc !reconstructed gradients for main variables
    real, allocatable, dimension(:, :, :)::gradientsturb   ! unlimited gradients for turbulent variables for
    real, allocatable, dimension(:, :, :)::gradientsturb_wall ! unlimited gradients for turbulent variables for wall cells
    real, allocatable, dimension(:, :, :)::gradients2, gradientsc2 ! reconstructed gradients for turbulent variables for wall cells
    real, allocatable, dimension(:)::gradientstemp !unlimited gradients for temperature
    real, allocatable, dimension(:, :)::gradientstemp_wall !unlimited gradients for temperature for wall cells
    real, allocatable, dimension(:, :, :)::cgradientstemp !unlimited gradients for temperature for wall cells
    real, allocatable, dimension(:, :)::velocitydof !unlimited gradients for velocity
    real, allocatable, dimension(:, :, :)::velocitydof_wall !unlimited gradients for velocity for wall cells
    real, allocatable, dimension(:, :)::invccjac  !inverse jacobian
    real, allocatable, dimension(:, :)::invctjac  !inverse jacobian transposed
    real, allocatable, dimension(:, :, :)::stencils, stencilsc  !stencils entries for matrix a (usually stored only for wall bounded cells)
    real, allocatable, dimension(:, :, :)::invmat_stencilt, invmat_stenciltc !pseudo inverse matrix for least squares reconstruction
    real, allocatable, dimension(:, :)::volume, volumec            !volume of elements in the stencil
    real, allocatable, dimension(:, :)::grads, gradsav !gradients obtained from green gauss approximation
    real, allocatable, dimension(:, :, :)::uleft, uleftx !boundary extrapolated value for main equations variables for considered cell
    real, allocatable, dimension(:, :, :)::qpoints !quadrature points
    real, allocatable, dimension(:, :, :)::uleftturb !boundary extrapolated value for turbulent equations variables for considered cell
    real, allocatable, dimension(:, :, :, :)::uleftv !boundary extrapolated values for main equations gradients for considered cell
    real, allocatable, dimension(:, :, :, :)::uleftturbv !boundary extrapolated values for turbulent equations gradients for considered cell
    real, allocatable, dimension(:, :)::weno  !weno weights for main equations variables
    real, allocatable, dimension(:, :)::weno2 !weno weights for turbulent equations variables
    real, allocatable, dimension(:, :, :, :)::wenos !weno weights for main equations with respect to characteristic variables
    real, allocatable, dimension(:, :, :, :)::findw !weno weights for main equations with respect to characteristic variables
    real, allocatable, dimension(:, :)::indicator, indicatorc !precomputed smoothness indicators
    real, allocatable, dimension(:, :)::dg2fv, weightl
    real, allocatable, dimension(:, :)::tempsq !constrained least squares reconstruction matrix for temperature gradient
    real, allocatable, dimension(:, :)::tempsqmat !constrained least squares reconstruction matrix for temperature gradient
    real, allocatable, dimension(:, :)::vellsq !constrained least squares reconstruction matrix for velocity gradient
    real, allocatable, dimension(:, :)::velinvlsqmat !constrained least squares reconstruction matrix for velocity gradient
    real, allocatable, dimension(:)::wallcoeff !constrained least squares gaussian elimination component
    real, allocatable, dimension(:)::wallcoefg !constrained least squares gaussian elimination component
    integer::g0 !constrained least squares gaussian elimination component
    integer::k0 !constrained least squares gaussian elimination component
    real, allocatable, dimension(:)::vext_ref !reference coordinates of the vertex by which the transformation has been based upon
    real, allocatable, dimension(:, :, :)::surf_qpoints ! physical space surface quadrature points (i_face, i_qp, xy) relative to cell center
    real, allocatable, dimension(:, :, :)::uleft_dg ! boundary extrapolated solution value (var, i_face, i_qp)
    real, allocatable, dimension(:, :, :)::rpoints, rotvel !radius of qpoints, rotational velocity
    real, allocatable, dimension(:)::mrf_origin, mrf_velocity
    integer, allocatable, dimension(:, :)::periodicflag
    real, allocatable, dimension(:, :, :, :)::br2_aux_var ! (var, dim, i_face, i_qp)
    real, allocatable, dimension(:, :, :)::br2_local_lift ! (var, dim, i_face)
  end type local_recon3

  type::neixx
    integer::cpu                              !cpu index
    integer, allocatable, dimension(:, :)::elem1 !elements indexing
    real, allocatable, dimension(:)::centers    !centers of cells for weno directionals stencils and compact stencils
  end type neixx

  type(local_recon3), allocatable, dimension(:)::ilocal_recon3
  type(neixx), allocatable, dimension(:)::neix1, neix2, neix3, neix4, neix5
  type(local_recon3), allocatable, dimension(:)::ilocal_recon4, ilocal_recon5, ilocal_recon6

  type exchange_cord
    integer::procid        !id of processor
    integer::howmany       ! how many elements are needed
    real, allocatable, dimension(:, :, :)::nodecord    !the coordinates of the element
  end type exchange_cord

  type(exchange_cord), allocatable, dimension(:)::iexcordr  !receiving data type
  type(exchange_cord), allocatable, dimension(:)::iexcords  !sending data type

  type exchange_solhi
    integer::procid    !id of processor
    integer::iavt      !number of processors that we need to receive/or send data
    integer::iavc      !number of processors that we need to receive/or send data
    integer::fast
    integer::howmany   !how many elements are needed
    real, allocatable, dimension(:, :)::sol   !array to hold the values to be send/received
  end type exchange_solhi

  type(exchange_solhi), allocatable, dimension(:)::iexsolhir, iexsolhird    !receiving data type for halo cells of stencils
  type(exchange_solhi), allocatable, dimension(:)::iexsolhis, iexsolhisd    !sending data type for halo cells of stencils

  type exchange_boundhi
    integer::procid   !id of processor
    integer::iavt     !number of processors that we need to receive/or send data
    integer::fast
    integer::howmany  !how many elements are needed
    real, allocatable, dimension(:, :)::facesol, facesol_m, facesol_dg  !array holding the number of variables to be sent/received
    integer, allocatable, dimension(:, :)::vertpp !vertex mapping between different processes
  end type exchange_boundhi

  type(exchange_boundhi), allocatable, dimension(:)::iexboundhir !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries
  type(exchange_boundhi), allocatable, dimension(:)::iexboundhis !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries
  type(exchange_boundhi), allocatable, dimension(:)::iexboundhiri !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
  type(exchange_boundhi), allocatable, dimension(:)::iexboundhisi !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
  type(exchange_boundhi), allocatable, dimension(:)::iexboundhirr  !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
  type(exchange_boundhi), allocatable, dimension(:)::iexboundhiss  !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping

  type::node_ne        !name of type for the set of nodes
    integer::itor  !index declaring if the coordinates of this node are needed for each cpu
    real, allocatable, dimension(:)::cord    !coordinates of the node
    integer, allocatable, dimension(:)::bct
    integer::itorm        !index refering to unique list index in cells
  end type node_ne
  type::node_lit        !name of type for the set of nodes
    integer::numberofneib      !number of elements that share this node
    integer, allocatable, dimension(:)::neibids, xne, xneib !ids of the elements, and counters
  end type node_lit

  type(node_ne), allocatable, dimension(:)::inoden, inoder, inoder4
  type(node_lit), allocatable, dimension(:)::inoder2

  type::elem_shape        !name of type for the set of nodes
    integer::ishape        !integer indicating the shape of each element
  end type elem_shape

  integer, allocatable, dimension(:)::ieshape        !1-d array for shape of element

  type::exchange
    integer::procid                                        !processor id
    integer::tot                    ! total number of elements
    integer, allocatable, dimension(:)::whattheyneed   !element numbers
    integer, allocatable, dimension(:)::muchtheyneed  !number of elements they need
    integer, allocatable, dimension(:)::muchineed        !number of elements i need
    integer, allocatable, dimension(:)::whatineed        !element number global
    integer, allocatable, dimension(:)::sideineed !side that i need
    integer, allocatable, dimension(:)::sideineedn, qineed    !from which cpu, which gaussian quadrature point order
    integer, allocatable, dimension(:)::sidetheyneed, qtheyneed !side they need, which gaussian quadrature point order
    integer, allocatable, dimension(:)::sidetheyneedn
    integer, allocatable, dimension(:)::localref !local reference
    integer, allocatable, dimension(:, :)::nodex
  end type exchange

  type(exchange), allocatable, dimension(:)::iexchanger !receive
  type(exchange), allocatable, dimension(:)::iexchanges !send
  type(exchange), allocatable, dimension(:)::iexchanger1
  type(exchange), allocatable, dimension(:)::iexchanges1

  type::solexchange
    integer::procid        !processor id
    real, allocatable, dimension(:, :)::centres   !cell centres
    integer, allocatable, dimension(:, :)::nodes  !nodes
    real, allocatable, dimension(:, :)::sol, sol_dg       !solution
  end type solexchange

  type(solexchange), allocatable, dimension(:)::solchanger  !receives

  type(solexchange), allocatable, dimension(:)::solchanges  !sends

  type::recex
    integer::procid                                        !processor id
    integer::tot
    integer, allocatable, dimension(:)::whattheyneed   !element numbers
    integer, allocatable, dimension(:)::muchtheyneed  !number of elements they need
    integer, allocatable, dimension(:)::muchineed        !number of elements i need
    integer, allocatable, dimension(:)::whatineed        !element number global
    integer, allocatable, dimension(:)::ishape   !shape of elements
    integer, allocatable, dimension(:)::localref !local referencing
    real, allocatable, dimension(:, :)::centers   !barycentres
  end type recex

  type(recex), allocatable, dimension(:)::irecexr                !receive elements due to stencils
  type(recex), allocatable, dimension(:)::irecexr1                !receive elements due to stencils
  type(recex), allocatable, dimension(:)::irecexrg                !receive elements due to stencils
  type(recex), allocatable, dimension(:)::irecexsg                !send elements due to stencils
  type(recex), allocatable, dimension(:)::irecexs1                !send elements due to stencils
  type(recex), allocatable, dimension(:)::irecexs                !send elements due to stencil

  ! name of type for the set of elements
  type::element_number
    integer::ihex        !local index of each cell
    integer::mood, mood_o        !mood flag for every element
    integer::recalc !flag for recalculating the solution post-eriori
    integer::ihexgl !global index of each cell
    integer::full   !specifies if sufficient number of stencils are found to proceed with weno for this cell
    integer::interior !specifies if this cell has any side bounded (interior cells get a value of 0, non interior ones get a value of 1)
    integer::itotalpoints, troubled
    integer::vdec, halo      !number of volume decompositions for each element
    integer::mode      !specifies if decomposition canbe avoided for straight sided element (not implemented yet)
    integer::ggs      ! specifies with what algorithm to compute the gradients (green gauss, least squares or blend of them)
    integer::idegfree  !degrees of freedom for polynomial selected
    integer::inumneighbours    !number of neighbours in each stencil
    integer::nofbc         !number of bounded faces
    integer::iorder    !order of polynomials
    integer::ishape    !> shape of element 1: hex 2: tet 3: pyramid 4: prism 5: quad 6: tri
    integer::ifca      !number of sides
    integer::admis     !number of admissible stencils
    integer::hybrid    !flag for switching to lower order discretisation as a function of wall distance
    integer::nonodes   !number of nodes
    integer::walls     !flag to declare if this is cell bounded by a wall
    integer::lump
    integer::filtered, reduce
    integer, allocatable, dimension(:)::bleedn
    integer, allocatable, dimension(:)::nojecount
    integer, allocatable, dimension(:)::nodes    !nodes index
    integer, allocatable, dimension(:, :)::nodes_neighbours    !nodes index
    integer, allocatable, dimension(:)::types_faces !type of each face (quadrilateral, triangle)
    integer, allocatable, dimension(:)::reorient     !consistency across interface in terms of ordering of quadrature points
    integer, allocatable, dimension(:)::num_of_wall_gqp !not used yet
    type(face_d), allocatable, dimension(:)::q_face  !quadrature face
    integer, allocatable, dimension(:, :)::nodes_faces, nodes_faces_v !counterclockwise numbering of the nodes for each face
    integer, allocatable, dimension(:)::indexi !indexing for the cells
    integer, allocatable, dimension(:)::ibounds, iboundw ! bounded codes for each bounded face
    integer, allocatable, dimension(:)::ineigh !neighbours local numbering
    integer, allocatable, dimension(:)::ineighg !neighbours global numbering
    integer, allocatable, dimension(:)::ineighb !neighbours cpu index
    integer, allocatable, dimension(:)::ineighn !neighbours numbering in other cpus
    integer, allocatable, dimension(:)::nodes_v !nodes_rearranged for output
    real::totvolume                !volume of element
    real::dtl, viscx           !local time step size
    real::minedge       !inscribed sphere radius
    real::stencil_dist    !stencil distance factor
    real::walldist        !wall distance
    real::xxc          !cell centre coordinates in x
    real::yyc           !cell centre coordinates in y
    real::zzc           !cell centre coordinates in z
    real::condition
    real::er, er2, er1, er2dt, er1dt, er1er2, lwcx2, diss, erx        !adda components
    real::linc
    real, dimension(1)::wcx
    real, allocatable, dimension(:)::faceanglex        !faceangle  x,y
    real, allocatable, dimension(:)::facediss
    real, allocatable, dimension(:)::faceangley        !faceangles x,y
    real, allocatable, dimension(:)::dih !distance across cell centres at each face
    real, allocatable, dimension(:, :)::dih2
    real, allocatable, dimension(:)::vortex, avars  !q criterion
    real, allocatable, dimension(:)::surf    !surface area
    real, allocatable, dimension(:)::delta_xyz ! 0.5(x_max - x_min),  0.5(y_max - y_min),  0.5(z_max - z_min)
    integer::condx
  end type element_number

  type::faces
    integer::ishape   !shape of each face
    integer, allocatable, dimension(:)::neighg, neighl   !global and local indexing of cells
    integer, allocatable, dimension(:)::nodes           !number of nodes
  end type faces

  type(element_number), allocatable, dimension(:, :)::ielem  !element number
  type(faces), allocatable, dimension(:, :)::iface

  type::bound_number        !name of type for the set of elements
    integer::ishape, icode, inum, which, face, ibid        !indices for the shape, bc code, number, and which face is bounded
    integer, allocatable, dimension(:)::ibl   !bc flag
    integer, allocatable, dimension(:)::localn, cpun   !local number and cpu for each bounded surface
    character(len=12)::description
  end type bound_number

  type(bound_number), allocatable, dimension(:, :)::ibound          !1-d array for boundary condition type

  type::connx        !name of type for the set of elements
    integer::procid        !cpu id
    integer, allocatable, dimension(:)::howmanyi        !how many element i need
    integer, allocatable, dimension(:)::retm     !how many entries to receive
    integer, allocatable, dimension(:)::ret      !how many entries to send
    integer, allocatable, dimension(:)::howmanythey !how many element they need
    integer, allocatable, dimension(:, :)::whichi !which element i need
    integer, allocatable, dimension(:, :)::whichthey  !which element they need
    real, allocatable, dimension(:, :)::facx      !fac
  end type connx

  type(connx), allocatable, dimension(:)::iconr
  type(connx), allocatable, dimension(:)::icons
  type(connx), allocatable, dimension(:)::iconrpa
  type(connx), allocatable, dimension(:)::iconrpm
  type(connx), allocatable, dimension(:)::iconspo
  type(connx), allocatable, dimension(:)::iconrpf

  type::node_number        !name of type for the set of nodes
    integer::noden        !identification number that can be used as a pointer inside an array
    integer::blockn !block number or cpu number of this element
    integer::nodegl
    real::x                !coordinates in x axis
    real::y                !coordinates in y axis
    real::z                !coordinates in z axis
  end type node_number

  type(node_number), allocatable, dimension(:, :)::inode          !1-d array for pointer type for nodes

  type(node_number)::itemp                                  !temporary pointer node iteration

end module declaration
