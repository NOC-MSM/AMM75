!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  : namelist for BDY generation tool
!!            
!!             User inputs for generating open boundary conditions
!!             employed by the BDY module in NEMO. Boundary data
!!             can be set up for v3.2 NEMO and above.
!!            
!!             More info here.....
!!            
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!------------------------------------------------------------------------------
!   vertical coordinate
!------------------------------------------------------------------------------
   ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)  
   ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
   ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
   rn_hmin     =   -10     !  min depth of the ocean (>0) or 
                           !  min number of ocean level (<0)

!------------------------------------------------------------------------------
!   s-coordinate or hybrid z-s-coordinate
!------------------------------------------------------------------------------
   rn_sbot_min =   10.     !  minimum depth of s-bottom surface (>0) (m)
   rn_sbot_max = 7000.     !  maximum depth of s-bottom surface 
                           !  (= ocean depth) (>0) (m)
   ln_s_sigma  = .false.   !  hybrid s-sigma coordinates
   rn_hc       =  150.0    !  critical depth with s-sigma

!------------------------------------------------------------------------------
!  grid information 
!------------------------------------------------------------------------------
   sn_src_hgr = './mesh_mask_eORCA025-GO6_subset.nc' 
   sn_src_zgr = './inputs_src_zgr.ncml'
   sn_dst_hgr = './domain_cfg.nc'
   sn_dst_zgr = './inputs_dst.ncml'
   sn_src_msk = './mesh_mask_eORCA025-GO6_subset.nc'
   sn_bathy   = './bathy_meter.nc'

!------------------------------------------------------------------------------
!  I/O 
!------------------------------------------------------------------------------
   sn_src_dir = 'GLOSEA6.ncml' ! src_files
   sn_dst_dir = './OUTPUT/V1/ATLANTIC/TIDE'
   sn_fn      = 'AMM75'                ! prefix for output files
   nn_fv      = -1e20                 ! set fill value for output files
   nn_src_time_adj = 0                ! src time adjustment
   sn_dst_metainfo = 'GLOSEA6-AMM75'

!------------------------------------------------------------------------------
!  unstructured open boundaries                         
!------------------------------------------------------------------------------
    ln_coords_file = .false.               !  =T : produce bdy coordinates files
    cn_coords_file = 'coords_open.bdy.nc' !  name of bdy coordinates files 
                                          !  (if ln_coords_file=.TRUE.)
    ln_mask_file   = .true.               !  =T : read mask from file
    cn_mask_file   = 'bdy_open_mask.nc'   !  name of mask file 
                                          !  (if ln_mask_file=.TRUE.)
    ln_dyn2d       = .false.               !  boundary conditions for 
                                          !  barotropic fields
    ln_dyn3d       = .false.               !  boundary conditions for 
                                          !  baroclinic velocities
    ln_tra         = .false.               !  boundary conditions for T and S
    ln_ice         = .false.              !  ice boundary condition   
    nn_rimwidth    = 6                   !  width of the relaxation zone

!------------------------------------------------------------------------------
!  unstructured open boundaries tidal parameters
!------------------------------------------------------------------------------
    ln_tide        = .true.              !  =T : produce bdy tidal conditions
    sn_tide_model  = 'FES2014'            !  Name of tidal model (FES2014|TPXO7p2)

    clname(1)      = 'SA'
    clname(2)      = 'SSA'
    clname(3)      = 'MM'
    clname(4)      = 'MF'
    clname(5)      = 'MTM'
    clname(6)      = 'MSF'
    clname(7)      = 'MSQM'
    clname(8)      = 'K1'
    clname(9)      = 'O1'
    clname(10)     = 'Q1'
    clname(11)     = 'P1'
    clname(12)     = 'S1'
    clname(13)     = 'J1'
    clname(14)     = 'M2'
    clname(15)     = 'N2'
    clname(16)     = 'S2'
    clname(17)     = 'K2'
    clname(18)     = 'L2'
    clname(19)     = 'T2'
    clname(20)     = 'R2'
    clname(21)     = 'MU2'
    clname(22)     = 'NU2'
    clname(23)     = '2N2'
    clname(24)     = 'MKS2'
    clname(25)     = 'LA2'
    clname(26)     = 'EPS2'
    clname(27)     = 'M3'
    clname(28)     = 'M4'
    clname(29)     = 'M6'
    clname(30)     = 'M8'
    clname(31)     = 'N4'
    clname(32)     = 'S4'
    clname(33)     = 'MN4'
    clname(34)     = 'MS4'

    ln_trans       = .true.                !  interpolate transport rather than velocities

!------------------------------------------------------------------------------
!  Time information
!------------------------------------------------------------------------------
    nn_year_000     = 1993           !  year start
    nn_year_end     = 1993           !  year end
    nn_month_000    = 1              !  month start (default = 1 is years>1)
    nn_month_end    = 1             !  month end (default = 12 is years>1)
    sn_dst_calendar = 'gregorian'    !  output calendar format
    nn_base_year    = 1900           !  base year for time counter
    ln_time_interpolation = .true.  !  temporal interpolation parent to child
    ! TPXO
	sn_tide_grid   = ''
	sn_tide_h      = ''
	sn_tide_u      = ''
    ! FES
    sn_tide_fes    = '/gws/nopw/j04/jmmp/jmmp_collab/AMM7/FES2014/'
 
!------------------------------------------------------------------------------
!  Additional parameters
!------------------------------------------------------------------------------
    nn_wei  = 1                   !  smoothing filter weights 
    rn_r0   = 0.041666666         !  decorrelation distance use in gauss
                                  !  smoothing onto dst points. Need to 
                                  !  make this a funct. of dlon
    sn_history  = 'test case'
                                  !  history for netcdf file
    ln_nemo3p4  = .true.          !  else presume v3.2 or v3.3
    nn_alpha    = 0               !  Euler rotation angle
    nn_beta     = 0               !  Euler rotation angle
    nn_gamma    = 0               !  Euler rotation angle
	rn_mask_max_depth = 100.0     !  Maximum depth to be ignored for the mask
	rn_mask_shelfbreak_dist = 20000.0 !  Distance from the shelf break
