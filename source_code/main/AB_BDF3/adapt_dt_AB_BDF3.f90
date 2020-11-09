SUBROUTINE adapt_dt_AB_BDF3(u,up,dt,dt1,Part)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy
  USE parameter_module, ONLY: dt_max,CFL_safety_fac,T_part
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Part
  REAL(kind=kr)  :: dt,dt1
  REAL(kind=kr)  :: dt_cfl
  
  dt_cfl = CFL_safety_fac*CFL_AB_BDF3(u,up,Part)

  IF (dt_cfl .GT. 1.05*dt1) THEN !allow to increase the time step by 5% maximum
     dt = MIN(1.05*dt1,dt_max)
  ELSE
     dt = MIN(dt_cfl,dt_max)
  ENDIF
 
END SUBROUTINE adapt_dt_AB_BDF3

