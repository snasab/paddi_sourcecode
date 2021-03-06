! This subroutine advances the solution on step in time 
! Currently, a semi-implicit BDF3/AB3 time stepping scheme is implemented.
! (third order accuracy, implicit part with backward-differencing, 
!  explicit part with third order Adams-Bashforth.)
SUBROUTINE timestep_AB_BDF3(u,Temp,Chem,up,Part,drag,istep,t,dt,dt1,dt2)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,shift_time_pointers,ltime0,ltime1
  USE parameter_module, ONLY : Nx,D_therm,D_comp,S_therm,S_comp,T_part,G_part,D_part,S_part
  USE message_passing_module, ONLY: myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  IMPLICIT NONE
  TYPE(velocity) :: u, up
  TYPE(buoyancy) :: Temp,Chem,Part
  INTEGER(kind=ki) :: istep
  REAL(kind=kr) :: t                ! dimensionless time
  REAL(kind=kr) :: dt,dt1,dt2       ! last three time step sizes
  REAL(kind=kr),POINTER :: drag(:,:,:,:) ! RHS 
!  ALLOCATE(drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec))
  

  ! shift time pointers (t_{n-1} -> t_{n-2}, t_n -> t_{n-1})
  CALL shift_time_pointers
  dt2=dt1
  dt1=dt

  ! time step size control
  CALL adapt_dt_AB_BDF3(u,up,dt,dt1,Part)
  t=t+dt

  ! Compute coefficients in the time stepping scheme 
  CALL comp_coeff_AB_BDF3(dt,dt1,dt2)
  
  ! compute new state in spectral space 
#ifdef TEMPERATURE_FIELD
  CALL tmstp_buoyancy_AB_BDF3(Temp,u,D_therm,S_therm) ! Compute FCs of Temp at new time level
#endif
#ifdef CHEMICAL_FIELD
  CALL tmstp_buoyancy_AB_BDF3(Chem,u,D_comp,S_comp) ! Compute FCs of Chem at new time level
#endif
#ifdef PARTICLE_FIELD
  CALL tmstp_particle_AB_BDF3(Part,up) ! Compute FCs of Part at new time level
  CALL tmstp_particle_vel_AB_BDF3(u,up,Part,drag)
#endif
  CALL tmstp_velocity_AB_BDF3(u,Temp,Chem,up,Part,drag)     ! Compute FCs of u at new time level


  ! compute "new" state in physical space 
  CALL compute_phys_space_vars(u,Temp,Chem,up,Part)
  
END SUBROUTINE timestep_AB_BDF3
