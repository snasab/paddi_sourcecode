SUBROUTINE read_parameter(restarted)
!+----------------------------------------------------------------------------------------------+
!| The code solves the following equations:                              						|
!|                                                                       						|
!| Du                                                                    						|
!| -- = - \nabla P + D_visc \nabla^2 u + B_therm T e_z - B_comp C e_z + R_part ((u_p - u)/T_part)*r   	|
!| Dt         																					|
!|																								|
!| Du_p                                                                    						|
!| ---- = (u_p - u)/T_part - G_part e_z   + Dv_part \nabla^2 up														|
!| Dt                                                                    						|
!|                                                                       						|
!| DT                                                                    						|
!| -- = D_therm \nabla^2 T - S_therm w                                   						|
!| Dt                                                                    						|
!|                                                                       						|
!| DC                                                                    						|
!| -- = D_comp \nabla^2 C - S_comp w                                     						|
!| Dt      																						|
!|                                                                       						|
!| dr                                                                    						|
!| -- = - div (u_p r) +  D_part \nabla^2 r - S_part (u_p e_z)                                   |					
!| dt                                                                     						|
!|                                                                       						|
!| \nabla * u = 0                                                        +----------------------+
!|                                                                       |
!| in a triply periodic cube.                                            |
!|                                                                       |
!| Coefficients:                                                         |
!|                                                                       |
!|   B_therm, B_comp represent the coefficients in fromt of the buoyancy |
!|                   terms;                                              |
!|   D_visc,D_therm,D_comp,D_part are the viscous,thermal,compositional, |
!|                   particulate diffusion parameters;                   |
!|   S_therm, S_comp, S_part are the coefficients resulting from the     |
!|                   background stratification;  		         |
!|   T_part represents the stopping time parameter;			 |
!|	 G_part represents the gravity parameter.			 |
!|   R_part scales the drag term accordingly.                            |
!|   Dv_part is the particle viscous diffusion parameter in particle     |
!|       momentum equation.                                              |
!|                                                                       |
!| Implementing the algorithm for these general equations has the        |
!| advantage that switching between different non-dimensionalizations    |
!| becomes a trivial matter.                                             |
!+-----------------------------------------------------------------------+  
!| This subroutine reads all the input parameters for the run.           |
!|                                                                       |
!| The namelist feature of FORTRAN 90 is used. The following             |
!| namelist entries can be specified by the user:                        |
!|                                                                       |
!| Restart_from_dumped_data (logical):                                   |
!|            True if current run should be started from saved           |
!|            restart data. False otherwise.                             |
!| Restart_from_netCDF_output_file(logical):                             |
!|            True if current run should be started from                 |
!|            a netcdf output file instead of a dedicated restart file.  |
!| max_degree_of_x_fourier_modes (integer):                              |
!|            Maximum degree of Fourier modes used for the               |
!|            expansion in x-direction.                                  |
!| max_degree_of_y_fourier_modes (integer):                              |
!|            Maximum degree of Fourier modes used for the               |
!|            expansion in y-direction.                                  |
!| max_degree_of_z_fourier_modes (integer):                              |
!|            Maximum degree of Fourier modes used for the               |
!|            expansion in z-direction.                                  |
!| dealiasing(logical): Apply the 3/2-rule to avoid aliasing errors      |
!| Thermal_buoyancy_parameter (real):                                    |
!|            The coefficient in front of the thermal buoyancy term.     |
!| Compositional_buoyancy_param (real):                                  |
!|            The coefficient in front of the compos. buoyancy term.     |
!| Viscous_diffusion_coeff (real):                                       |
!|            The coefficient in front of the viscous diffusion term.    |
!| Thermal_diffusion_coeff (real):                                       |
!|            The coefficient in front of the thermal diffusion term.    |
!| Compositional_diffusion_coeff (real):                                 |
!|            The coefficient in front of the compos. diffusion term.    |
!| Particle_viscous_diffusion_coeff (real):	        		 |
!|	      The coefficient in front of the particle diffusion term    |
!|            in momentum equation.                                      |
!| Particle_diffusion_coeff (real):                                      |
!|            The coefficient in front of the particle diffusion term.   |
!| Thermal_stratif_param (real):                                         |
!|            The coefficient in front of the thermal stratification     |
!|            advection term.                                            |
!| Compositional_stratif_param (real):                                   |
!|            The coefficient in front of the thermal stratification     |
!|            advection term.                                            |
!| Particle_stratif_param (real):                                        |
!|            The coefficient in front of the particle stratification    |
!|            advection term.                                            |
!| Stopping_time_param (real):                                           |
!|            The coefficient in front of the particle stopping time     |
!|            coefficient in front of the drag term.                     |
!| Gravity_param (real):                                                 |
!|            The coefficient in front of the particle gravity term.     |
!| Particle_scaling_coeff (real):    					 |
!|            The coefficient in front of the drag term.                 |
!| x_extent_of_the_box,y_extent_of_the_box,z_extent_of_the_box (real):   |
!|            The governing equations are solved in the cube             |
!|            [0,x_extent_of_the_box] x [0,y_extent_of_the_box] x        |
!|            [0,z_extent_of_the_box]                                    |
!| initial_time_step_length (real):                                      |
!|            Time step size used for the first time step                |
!| maximum_time_step_length (real):                                      |
!|            Maximum allowed time step size                             |
!| CFL_safety_factor (real):                                             |
!|            Safety factor for Courant-Friedrichs-Levi-Criterion:       |
!|            It dt_CFL is the predicted stability limit, the time step  |
!|            dt is always chosen so that                                |
!|              dt <= CFL_safety_factor * dt_CFL                         |
!| number_of_time_steps (integer):                                       |
!|            The total number of time steps to be performed.            |
!| save_state_every_nth_timestep (integer):                              |
!|            The computed fields are saved in compressed form every     |
!|            save_state_every_nth_timestep-th timestep. The jc library  |
!|            written by Joerg Schmalzl is used, which uses JPEG         |
!|            compression algorithms.                                    |
!|            CAUTION: The compression is LOSSY and SERIAL. The fields   |
!|            are collected on a master node before they are written     |
!|            to the file. We assume that the master node has enough     |
!|            memory for this...                                         |
!| save_state_netCDF_ev_nth_step (integer):                              |
!|            The computed fields are saved in netCDF format every       |
!|            save_state_netcdf_ev_nth_step-th timestep. The pnectCDF    |
!|            library is used for this.                                  |
!| restart_info_every_nth_timestep (integer):                            |
!|            Restart information is written to disk every               |
!|            restart_info_every_nth_timestep-th timestep.               |
!| comp_diagno_every_nth_timestep (integer):                             |
!|            Diagnostical parameters such as Nusselt number,            |
!|            RMS value of velocity, ... are computed every              |
!|            comp_diagno_every_nth_timestep-th timestep.                |
!| write_spec_every_nth_timestep (integer):                              |
!|            Information concerning the spectra of the solution are     |
!|            written to disc every  write_spec_every_nth_timestep-th    |
!|            time step.                                                 |
!| write_prof_every_nth_timestep (integer):                              |
!|            Profiles are written to disk every                         |
!|            write_prof_every_nth_timestep-th time step                 |
!| Name_of_input_restart_file (character (len=100)):                     |
!|            The name of the file from which the restart information    |
!|            can be read. This can either be an explicit restart file   |
!|            or it can be a netCDF output file. In the latter case      |
!|            the last saved state is used to restart the run.           |
!| Name_of_output_restart_file (character (len=100)):                    |
!|            The name of the file to which the restart information is   |
!|            written.                                                   |
!| write_state_to_compressed_file (logical):                             |
!|            Determindes if the computed fields should be written       |
!|            to a compressed data file.                                 |
!|            CAUTION: The compression is LOSSY and SERIAL. The fields   |
!|            are collected on a master node before they are written     |
!|            to the file. We assume that the master node has enough     |
!|            memory for this...                                         |
!| write_state_to_netCDF_file (logical):                                 |
!|            Determines if the computed fields should be written       |
!|            to a netCDF file. The parallel netCDF library is used.     |
!| Name_of_compressed_data_file (character (len=100)):                   |
!|            The name of the file to which the computed fields are      |
!|            written in compressed form.                                |
!| Name_of_netCDF_data_file (character (len=100)):                       |
!|            The name of the file to which the computed fields are      |
!|            written in netCDF format.                                  |
!| Name_of_diagnostics_data_file (character (len=100)):                  |
!|            The name of the file to which global diagnostical values   |
!|            (such as various energies, mean temperature, etc) are      |
!|            written.                                                   |
!| Name_of_horizontal_spectra_file (character (len=100)):                |
!|            The name of the file to which vertically avaraged Fourier  |
!|            spectra of various energies are written.                   |
!| Name_of_vertical_spectra_file (character (len=100)):                  |
!|            The name of the file to which the vertical spectra of      |
!|            the unknowns are written.                                  |
!| Name_of_z_profile_file (character (len=100)):                         |
!|            The name of the file to which vertical profiles are        |
!|            written.                                                   |
!| Use_an_FFTW_wisdom_file (logical):                                    |
!|            If TRUE, the program tries to import information about     |
!|            the fastest FFT algorithm for the given problem size       |
!|            from disk. See FFTW manual for details.                    |
!| Name_of_FFTW_wisdom_file (logical):                                   |
!|            The name of the file from which FFTW tries to import       |
!|            information about the fastest FFT algorithm. If no         |
!|            valuable information is contained in that file, FFTW       |
!|            tries a huge number of algorithms and measures their       |
!|            performance, which could take quite a long time. When the  |
!|            fastest available algorithm has been determined, the file  |
!|            Name_of_FFTW_wisdom_file is updated.                       |
!| number_of_tasks_1st_transpose(integer):                               |
!|            The number of processes used for the first transpose.      |
!| number_of_tasks_2nd_transpose(integer):                               |
!|            The number of processes used for the second transpose.     |
!|                                                                       |
!| Author: Stephan Stellmach                                             |
!+-----------------------------------------------------------------------+
  USE parameter_module
  USE IO_module, ONLY    : outfile,jc_out,nout,write_compressed_fields,  &
  &                        n_comp_diag,n_wrt_jc,n_wrt_netCDF,n_wrt_dump, &
  &                        n_wrt_spec,n_wrt_prof
  USE pnetCDF_IO_module, ONLY : write_pnetCDF_sim_dat,              &
  &                             netCDF_simdat_file_name,            &
  &                             netCDF_in_simdat_file_name,         &
  &                             netCDF_in_dump_file_name,           &
  &                             netCDF_out_dump_file_name,          &
  &                             pn_read_size_and_pa_from_dump,      &
  &                             pn_read_size_and_pa_from_simdat
  USE message_passing_module, ONLY : nprocs1,nprocs2,  &
  &                                  myid,numtasks,stop_mpi
  USE defprecision_module
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  CHARACTER          :: restarted
!
  INTEGER (kind=ki)  :: ierr
!
  LOGICAL            :: Restart_from_dumped_data        = .FALSE.
  LOGICAL            :: Restart_from_netCDF_output_file = .FALSE.
  INTEGER (kind=ki)  :: max_degree_of_x_fourier_modes
  INTEGER (kind=ki)  :: max_degree_of_y_fourier_modes
  INTEGER (kind=ki)  :: max_degree_of_z_fourier_modes
  LOGICAL            :: dealiasing = .TRUE.
  LOGICAL            :: RandNum_using_myid = .TRUE. 
  REAL (kind=kr)     :: Thermal_buoyancy_param 
  REAL (kind=kr)     :: Compositional_buoyancy_param
  REAL (kind=kr)     :: Viscous_diffusion_coeff 
  REAL (kind=kr)     :: Thermal_diffusion_coeff 
  REAL (kind=kr)     :: Compositional_diffusion_coeff
  REAL (kind=kr)     :: Thermal_stratif_param
  REAL (kind=kr)     :: Compositional_stratif_param
  REAL (kind=kr)     :: Stopping_time_param
  REAL (kind=kr)     :: Gravity_param
  REAL (kind=kr)     :: Particle_visc_diffusion_coeff
  REAL (kind=kr)     :: Particle_diffusion_coeff
  REAL (kind=kr)     :: Particle_stratif_param
  REAL (kind=kr)     :: Particle_scaling_coeff
  REAL (kind=kr)     :: x_extent_of_the_box             = 2._kr*pi
  REAL (kind=kr)     :: y_extent_of_the_box             = 2._kr*pi
  REAL (kind=kr)     :: z_extent_of_the_box             = 2._kr*pi
  REAL (kind=kr)     :: initial_time_step_length        = 1.E-5_kr
  REAL (kind=kr)     :: maximum_time_step_length        = 1.E-5_kr
  REAL (kind=kr)     :: CFL_safety_factor               = 0.4_kr
  INTEGER (kind=ki)  :: number_of_time_steps            = 100
  INTEGER (kind=ki)  :: save_state_every_nth_timestep   = 10 
  INTEGER (kind=ki)  :: save_state_netcdf_ev_nth_step   = 100
  INTEGER (kind=ki)  :: comp_diagno_every_nth_timestep  = 1 
  INTEGER (kind=ki)  :: restart_info_every_nth_timestep = 100
  INTEGER (kind=ki)  :: write_spec_every_nth_timestep   = 100
  INTEGER (kind=ki)  :: write_prof_every_nth_timestep   = 100 
  INTEGER (kind=ki)  :: number_of_tasks_1st_transpose   
  INTEGER (kind=ki)  :: number_of_tasks_2nd_transpose 
  CHARACTER (LEN=100) :: Name_of_input_restart_file      = "DUMP_IN"
  CHARACTER (LEN=100) :: Name_of_output_restart_file     = "DUMP_OUT"
  CHARACTER (LEN=100) :: Name_of_compressed_data_file    = "j__data"
  CHARACTER (LEN=100) :: Name_of_netCDF_data_file        = "Sim_data_netCDF"
  CHARACTER (LEN=100) :: Name_of_diagnostics_data_file   = "OUT" 
  CHARACTER (LEN=100) :: Name_of_horizontal_spectra_file = "XYSPEC"
  CHARACTER (LEN=100) :: Name_of_vertical_spectra_file   = "XYSPEC"
  CHARACTER (LEN=100) :: Name_of_z_profile_file          = "ZPROF"
  LOGICAL             :: write_state_to_compressed_file  = .TRUE.
  LOGICAL             :: write_state_to_netCDF_file      = .TRUE.
  LOGICAL             :: Use_an_FFTW_wisdom_file         = .TRUE.
  CHARACTER (LEN=100) :: Name_of_FFTW_wisdom_file        = "FFTW_WISDOM" 

  CHARACTER (len=100) :: in_dumpfile,out_dumpfile
  LOGICAL             :: restarted_dump,restarted_simdat
! 
  NAMELIST /input_values/ Restart_from_dumped_data
  NAMELIST /input_values/ Restart_from_netCDF_output_file
  NAMELIST /input_values/ max_degree_of_x_fourier_modes
  NAMELIST /input_values/ max_degree_of_y_fourier_modes
  NAMELIST /input_values/ max_degree_of_z_fourier_modes
  NAMELIST /input_values/ dealiasing
  NAMELIST /input_values/ RandNum_using_myid
  NAMELIST /input_values/ Thermal_buoyancy_param 
  NAMELIST /input_values/ Compositional_buoyancy_param
  NAMELIST /input_values/ Viscous_diffusion_coeff 
  NAMELIST /input_values/ Thermal_diffusion_coeff 
  NAMELIST /input_values/ Compositional_diffusion_coeff
  NAMELIST /input_values/ Thermal_stratif_param
  NAMELIST /input_values/ Compositional_stratif_param
  NAMELIST /input_values/ Stopping_time_param
  NAMELIST /input_values/ Gravity_param
  NAMELIST /input_values/ Particle_visc_diffusion_coeff
  NAMELIST /input_values/ Particle_diffusion_coeff
  NAMELIST /input_values/ Particle_stratif_param
  NAMELIST /input_values/ Particle_scaling_coeff
  NAMELIST /input_values/ x_extent_of_the_box
  NAMELIST /input_values/ y_extent_of_the_box
  NAMELIST /input_values/ z_extent_of_the_box
  NAMELIST /input_values/ initial_time_step_length
  NAMELIST /input_values/ maximum_time_step_length
  NAMELIST /input_values/ CFL_safety_factor
  NAMELIST /input_values/ number_of_time_steps
  NAMELIST /input_values/ save_state_every_nth_timestep
  NAMELIST /input_values/ save_state_netcdf_ev_nth_step
  NAMELIST /input_values/ comp_diagno_every_nth_timestep
  NAMELIST /input_values/ restart_info_every_nth_timestep
  NAMELIST /input_values/ write_spec_every_nth_timestep
  NAMELIST /input_values/ write_prof_every_nth_timestep
  NAMELIST /input_values/ number_of_tasks_1st_transpose
  NAMELIST /input_values/ number_of_tasks_2nd_transpose
  NAMELIST /input_values/ write_state_to_compressed_file
  NAMELIST /input_values/ write_state_to_netCDF_file
  NAMELIST /input_values/ Name_of_input_restart_file
  NAMELIST /input_values/ Name_of_output_restart_file
  NAMELIST /input_values/ Name_of_compressed_data_file
  NAMELIST /input_values/ Name_of_netCDF_data_file
  NAMELIST /input_values/ Name_of_diagnostics_data_file
  NAMELIST /input_values/ Name_of_horizontal_spectra_file
  NAMELIST /input_values/ Name_of_vertical_spectra_file
  NAMELIST /input_values/ Name_of_z_profile_file
  NAMELIST /input_values/ Use_an_FFTW_wisdom_file
  NAMELIST /input_values/ Name_of_FFTW_wisdom_file
!
  IF (myid.EQ.0) THEN 
     ! 
     ! Read data from given namelist
     ! -----------------------------
     READ(*,NML=input_values,IOSTAT=ierr)  
     IF (ierr .NE. 0 ) THEN 
        PRINT*,myid,'READ ERROR IN INPUT NAMELIST'
        STOP
     ENDIF

     !
     restarted_dump       =  Restart_from_dumped_data
     restarted_simdat     =  Restart_from_netCDF_output_file
     nsteps               =  number_of_time_steps
     n_wrt_jc             =  save_state_every_nth_timestep
     n_wrt_netCDF         =  save_state_netcdf_ev_nth_step
     n_wrt_dump           =  restart_info_every_nth_timestep
     n_comp_diag          =  comp_diagno_every_nth_timestep
     n_wrt_spec           =  write_spec_every_nth_timestep
     n_wrt_prof           =  write_prof_every_nth_timestep
     outfile(1)           =  Name_of_diagnostics_data_file
     outfile(2)           =  Name_of_horizontal_spectra_file
     outfile(3)           =  Name_of_vertical_spectra_file
     outfile(4)           =  Name_of_z_profile_file
     in_dumpfile          =  Name_of_input_restart_file
     out_dumpfile         =  Name_of_output_restart_file
     jc_out               =  Name_of_compressed_data_file
     rn_using_myid        =  RandNum_using_myid
     netCDF_simdat_file_name = Name_of_netCDF_data_file
     use_FFTW_wisdom_file =  Use_an_FFTW_wisdom_file
     FFTW_wisdom_file     =  Name_of_FFTW_wisdom_file
     dt_max               =  maximum_time_step_length
     CFL_safety_fac       =  CFL_safety_factor
     nprocs1              =  number_of_tasks_1st_transpose
     nprocs2              =  number_of_tasks_2nd_transpose
     Lmax                 = max_degree_of_x_fourier_modes
     Mmax                 = max_degree_of_y_fourier_modes
     Nmax                 = max_degree_of_z_fourier_modes
     B_therm              = Thermal_buoyancy_param
     B_comp               = Compositional_buoyancy_param
     D_visc               = Viscous_diffusion_coeff
     D_therm              = Thermal_diffusion_coeff
     D_comp               = Compositional_diffusion_coeff
     S_therm              = Thermal_stratif_param
     S_comp               = Compositional_stratif_param
     T_part               = Stopping_time_param
     G_part               = Gravity_param
     Dv_part              = Particle_visc_diffusion_coeff
     D_part               = Particle_diffusion_coeff
     S_part               = Particle_stratif_param
     R_part               = Particle_scaling_coeff
     Gammax               = x_extent_of_the_box
     Gammay               = y_extent_of_the_box
     Gammaz               = z_extent_of_the_box
     dt_initial           = initial_time_step_length     
     write_compressed_fields   = write_state_to_compressed_file
     write_pnetCDF_sim_dat     = write_state_to_netCDF_file

  ENDIF
! scatter input parameters to all the other processes
  CALL MPI_BCAST(restarted_dump        ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(restarted_simdat      ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nsteps                ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_wrt_jc              ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_wrt_netCDF          ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_wrt_dump            ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_comp_diag           ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_wrt_spec            ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(n_wrt_prof            ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(outfile               ,nout*100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(in_dumpfile           ,100, MPI_CHARACTER    ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(out_dumpfile          ,100, MPI_CHARACTER    ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(jc_out                ,100, MPI_CHARACTER    ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(rn_using_myid         ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(use_FFTW_wisdom_file  ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(FFTW_wisdom_file      ,100, MPI_CHARACTER    ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt_max                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(CFL_safety_fac        ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nprocs1               ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nprocs2               ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Lmax                  ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Mmax                  ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Nmax                  ,1,   MPI_INTEGER      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dealiasing            ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(B_therm               ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(B_comp                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(D_visc                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(D_therm               ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(D_comp                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(S_therm               ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(S_comp                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(T_part                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(G_part                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Dv_part               ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(D_part                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(S_part                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(R_part                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Gammax                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Gammay                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(Gammaz                ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt_initial            ,1,PM_MPI_FLOAT_TYPE,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(write_compressed_fields  ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(write_pnetCDF_sim_dat ,1,   MPI_LOGICAL      ,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(netCDF_simdat_file_name,100, MPI_CHARACTER    ,0,MPI_COMM_WORLD,ierr)

! Compute derived quantities
! --------------------------
!

  IF (dealiasing) THEN 
     Nx = MAX( 3*(2*Lmax) / 2, 1) 
     Ny = MAX( 3*(2*Mmax) / 2, 1) 
     Nz = MAX( 3*(2*Nmax) / 2, 1) 
  ELSE
     Nx = MAX( 2*Lmax, 1)
     Ny = MAX( 2*Mmax, 1)
     Nz = MAX( 2*Nmax, 1)
  ENDIF

#ifdef TWO_DIMENSIONAL
! In the 2d case, set Ny = 1 and Mmax = 0.  
! Gammay is set to one. This provides additional safety in case 
! I forget to remove it from a normalization factor.
! beta is set to zero.   
  Ny = 1
  Mmax = 0
  Gammay = 1._kr
  beta = 0._kr
  nprocs1 = 1
  nprocs2 = numtasks
#else 
  beta =2._kr*pi/Gammay
#endif

  alpha=2._kr*pi/Gammax
  gamma=2._kr*pi/Gammaz


! Open restart files if required

  netCDF_out_dump_file_name = TRIM(out_dumpfile) // '.cdf'
  netCDF_simdat_file_name = TRIM(netCDF_simdat_file_name) // '.cdf'
  restarted="N"
  IF (restarted_dump) THEN 
     ! Read parameter values from DUMP file and check if they agree with 
     ! those given by the user. 
     restarted="D"
     netCDF_in_dump_file_name  = TRIM(in_dumpfile) // '.cdf'
     CALL pn_read_size_and_pa_from_dump
  ELSE 
     IF (restarted_simdat) THEN 
        ! Read parameter values from NetCDF simdat file and check if they agree with 
        ! those given by the user.
        restarted="S" 
        netCDF_in_simdat_file_name  = TRIM(in_dumpfile) // '.cdf'
        CALL pn_read_size_and_pa_from_simdat
     ENDIF
  ENDIF


  IF (myid.EQ.0) THEN 
     !
     ! Write input to standart output for safety
     ! -----------------------------------------
     IF (restarted_dump) PRINT*, '*** Restarted form DUMP File ***'
     IF (restarted_simdat) PRINT*, '*** Restarted form SIMDAT File ***'
     PRINT*,'Retain x-modes with harmonic degrees up to Lmax =',Lmax
     PRINT*,'Retain y-modes with harmonic degrees up to Mmax =',Mmax
     PRINT*,'Retain z-modes with harmonic degrees up to Nmax =',Nmax
     PRINT*,'Dealiasing of Fourier modes ?',dealiasing
     PRINT*,'Length (x-direction) of the box =',Gammax
     PRINT*,'Width (y-direction) of the box =',Gammay
     PRINT*,'Height (z-direction) of the box =',Gammaz
     PRINT*,'Initial time step size =',dt_initial
     PRINT*,'Maximum time step size =',dt_max
     PRINT*,'CFL-Safety factor =', CFL_safety_fac
     PRINT*,'Number of time steps =',nsteps
     PRINT*,'Write out every nth time step, n=',n_wrt_jc
#ifdef TWO_DIMENSIONAL
     PRINT*,"********************** CAUTION: 2D Version *************************"
     PRINT*,"** Parameter max_degree_of_y_fourier_modes automatically set to 0 **"
     PRINT*,"** Parameters Gammay, number_of_tasks_1st_transpose and           **"
     PRINT*,"** number_of_tasks_2nd_transpose have no influence on the code    **"
     PRINT*,"** behavior.                                                      **" 
     PRINT*,"********************************************************************"
#endif
  ENDIF

END SUBROUTINE read_parameter
