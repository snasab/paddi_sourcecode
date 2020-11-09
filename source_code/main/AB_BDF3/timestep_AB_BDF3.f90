! This subroutine advances the solution on step in time 
! Currently, a semi-implicit BDF3/AB3 time stepping scheme is implemented.
! (third order accuracy, implicit part with backward-differencing, 
!  explicit part with third order Adams-Bashforth.)
SUBROUTINE timestep_AB_BDF3(u,Temp,Chem,istep,t,dt,dt1,dt2)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,shift_time_pointers,ltime0,ltime1
  USE parameter_module, ONLY : D_therm,D_comp,S_therm,S_comp
  USE message_passing_module, ONLY: myid
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(buoyancy) :: Temp,Chem
  INTEGER(kind=ki) :: istep
  REAL(kind=kr) :: t                ! dimensionless time
  REAL(kind=kr) :: dt,dt1,dt2       ! last three time step sizes

  ! shift time pointers (t_{n-1} -> t_{n-2}, t_n -> t_{n-1})
  CALL shift_time_pointers
  dt2=dt1
  dt1=dt

  ! time step size control
  CALL adapt_dt_AB_BDF3(u,dt,dt1)
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
  CALL tmstp_velocity_AB_BDF3(u,Temp,Chem)     ! Compute FCs of u at new time level

  ! compute "new" state in physical space 
  CALL compute_phys_space_vars(u,Temp,Chem)
  
END SUBROUTINE timestep_AB_BDF3
