!+------------------------------------------------------------+
!| This subroutine performs all the necessary initializations |
!| like setting initial:q condition, setting up the FFT, ...    |
!+------------------------------------------------------------+
SUBROUTINE init(u,Temp,Chem,up,Part,t,dt,dt1,dt2,startstep,restarted)
  USE defprecision_module
  USE MPI
  USE parameter_module, ONLY: Lmax,Mmax,Nmax,Nx,Ny,Nz,Gammax,Gammay,Gammaz,dt_initial, &
     &                        FFTW_wisdom_file,uwisdom,              &
     &                        allocate_fft_storage_scheme,           &
     &                        init_fft_storage_scheme,use_FFTW_wisdom_file
  USE mpi_transf_module, ONLY: init_transforms
  USE message_passing_module, ONLY: nprocs1,nprocs2
  USE state_module, ONLY: velocity,buoyancy,allocate_uTC,                     &
      &                   init_u_phys,init_Temp_phys,init_Chem_phys,          &
      &                   init_up_phys,init_Part_phys,set_initial_condition
  USE pnetCDF_IO_module, ONLY: pn_read_state_from_dump,pn_read_state_from_simdat
  USE minit_module, ONLY: M_init, compute_particle_mass
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  CHARACTER      :: restarted
  REAL(kind=kr)  :: t,dt,dt1,dt2
  INTEGER(kind=ki) :: startstep,ierr

! initialize the parallel FFTs  
  PRINT*,Nx,Ny,Nz,restarted
  IF (use_FFTW_wisdom_file) THEN 
     CALL init_transforms(Lmax,Mmax,Nmax,Nx,Ny,Nz,               &
          &               nprocs1,nprocs2,MPI_COMM_WORLD,        &
          &               FFTW_wisdom_file,uwisdom)
  ELSE
     CALL init_transforms(Lmax,Mmax,Nmax,Nx,Ny,Nz,               &
          &               nprocs1,nprocs2,MPI_COMM_WORLD)
  ENDIF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! initialize the FFT storage scheme
  CALL allocate_fft_storage_scheme
  CALL init_fft_storage_scheme

! allocate memory to store the system state
  CALL allocate_uTC(u,Temp,Chem,up,Part)

! set initial state
  IF (restarted.EQ."D" .OR. restarted.EQ."d") THEN 
     CALL pn_read_state_from_dump(u,Temp,Chem,up,Part,startstep,t,dt)
     CALL compute_phys_space_vars(u,Temp,Chem,up,Part)
  ELSE 
     ! initialize the fields in physical space
     IF (restarted.EQ."S".OR.restarted.EQ."s") THEN
        CALL pn_read_state_from_simdat(u,Temp,Chem,up,Part,startstep,t,dt)
        PRINT*,"DT=",dt
     ELSE
        CALL init_u_phys(u)
#ifdef TEMPERATURE_FIELD
        CALL init_Temp_phys(Temp)
#endif
#ifdef CHEMICAL_FIELD
        CALL init_Chem_phys(Chem)
#endif
#ifdef PARTICLE_FIELD
        CALL init_Part_phys(Part)
        CALL init_up_phys(up)
#endif
        startstep=0
        t = 0._kr
        dt =  dt_initial
     ENDIF

     CALL set_initial_condition(u,Temp,Chem,up,Part)
  ENDIF
!#ifdef PARTICLE_FIELD
!  call compute_particle_mass(M_init,Part)
!  M_init = Gammax*Gammay*Gammaz
!#endif
  dt1 = 0._kr  ! use RK2/CN as starting scheme in multistep AB/BDF3 version
  dt2 = 0._kr  ! there is no meaningful size of previous time steps

END SUBROUTINE init
