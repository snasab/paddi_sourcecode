!+--------------------------------------------------------------------+
!| The following subroutine computes various diagnostical parameters  |
!| which are then written to an external output file.                 |
!| Caution: It is assumed that this routine is not called at every    |
!|          time step so that it does not need to be optimized.       |
!|          In fact it is likely to be highly inefficient...          |
!|          To lazy to optimize now, any volunteers welcome!          |
!+--------------------------------------------------------------------+
!| Author: Stephan Stellmach                                          |
!+--------------------------------------------------------------------+
SUBROUTINE write_diagnostics_file(u,Temp,Chem,up,Part,istep,t,dt)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity,ltime0
  USE diagnostics_module, ONLY: compute_uTC_rms,compute_uTC_minmax,compute_average_flux, &
      &                         dissipation_buo,compute_upP_minmax,compute_upP_rms
  USE message_passing_module, ONLY : myid
  USE testing_module, ONLY: error_Temp_adv_diff,peak_div_u
  IMPLICIT NONE
  TYPE(velocity)   :: u,up
  TYPE(buoyancy)   :: Temp,Chem,Part
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: t,dt
  REAL(kind=kr)    :: urms,VORTrms,TEMPrms,CHEMrms,flux_Temp,flux_Chem,err,err_part
  REAL(kind=kr)    :: div_u,div_up
  REAL(kind=kr)    :: uprms,VORTprms,PARTrms
  REAL(kind=kr)    :: uxrms,uyrms,uzrms,VORTXrms,VORTYrms,VORTZrms
  REAL(kind=kr)    :: upxrms,upyrms,upzrms,VORTpXrms,VORTpYrms,VORTpZrms
  REAL(kind=kr)    :: Temp_min,Temp_max,Chem_min,Chem_max,           &
     &                u_min(3),u_max(3),VORT_min(3),VORT_max(3),     &
     &                u_max_abs,VORT_max_abs,                        &
	 &                flux_Part,Part_min,Part_max,                   &
     &                up_min(3),up_max(3),VORTp_min(3),VORTp_max(3), &
     &                up_max_abs,VORTp_max_abs
  REAL(kind=kr)    :: diss_Temp,diss_Chem,diss_Part,diss_mom
! compute rms-values
  CALL compute_uTC_rms(urms,VORTrms,TEMPrms,CHEMrms,uxrms,uyrms,uzrms, &
  &                    VORTXrms,VORTYrms,VORTZrms,u,Temp,Chem)
  
  CALL compute_upP_rms(uprms,VORTprms,PARTrms,upxrms,upyrms,upzrms,    &
  &                        VORTpXrms,VORTpYrms,VORTpZrms,up,Part)

! compute minimum and maximum values 
  CALL compute_uTC_minmax(Temp_min,Temp_max,Chem_min,Chem_max,   &
       &                  u_min,u_max,VORT_min,VORT_max,         &
       &                  u_max_abs,VORT_max_abs,                &
       &                  u,Temp,Chem)

  CALL   compute_upP_minmax(Part_min,Part_max,up_min,up_max,     &
	   &                  VORTp_min,VORTp_max,                   &
	   &                  up_max_abs,VORTp_max_abs,              &
	   &                  up,Part                                )

! compute fluxes 
#ifdef TEMPERATURE_FIELD
  CALL compute_average_flux(flux_Temp,Temp,u)
#else 
  flux_Temp = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  CALL compute_average_flux(flux_Chem,Chem,u)
#else
  flux_Chem = 0._kr
#endif
#ifdef PARTICLE_FIELD
  CALL compute_average_flux(flux_Part,Part,up)
#else 
  flux_Part = 0._kr
#endif


! compute momentum, thermal, and chemical dissipation rates
#ifdef TWO_DIMENSIONAL
diss_mom = dissipation_buo(u%spec(:,:,:,vec_x,ltime0)) + &
dissipation_buo(u%spec(:,:,:,vec_z,ltime0))
#else 
diss_mom = dissipation_buo(u%spec(:,:,:,vec_x,ltime0)) + &
dissipation_buo(u%spec(:,:,:,vec_y,ltime0)) &
         + dissipation_buo(u%spec(:,:,:,vec_z,ltime0))
#endif
#ifdef TEMPERATURE_FIELD
  diss_Temp = dissipation_buo(Temp%spec(:,:,:,ltime0))
#else
  diss_Temp = 0._kr
#endif
#ifdef CHEMICAL_FIELD
  diss_Chem = dissipation_buo(Chem%spec(:,:,:,ltime0))
#else
  diss_Chem = 0._kr
#endif
#ifdef PARTICLE_FIELD
  diss_Part = dissipation_buo(Part%spec(:,:,:,ltime0))
#else
  diss_Part = 0._kr
#endif

! check for peak in pectrum of div(u)

  !call error_Temp_adv_diff(err,Temp,t)
  CALL peak_div_u(u,err)
#ifdef PARTICLE_FIELD
  CALL peak_div_u(up,err_part)
#else
    err_part = 0._kr
#endif
  

  IF (myid.EQ.0) THEN
     !PRINT*,"Peak spectrum div(u) = ",err
     WRITE(uout(1),'(i7,65E20.7)')                                                         &
          !         1     2 3  4    5       6       7       8       9
          &         istep,t,dt,urms,VORTrms,TEMPrms,CHEMrms,flux_Temp,flux_Chem,           &
          !         10       11       12       13 
          &         Temp_min,Temp_max,Chem_min,Chem_max,                                   &
          !         14       15       16       17       18       19 
          &         u_min(1),u_max(1),u_min(2),u_max(2),u_min(3),u_max(3),                 &
          !         20           21           22           23
          &         VORT_min(1),VORT_max(1),VORT_min(2),VORT_max(2),                       &
          !         24           25           26        27 
          &         VORT_min(3),VORT_max(3),u_max_abs,VORT_max_abs,                        &
          !         28    29    30    31       32       33   
          &         uxrms,uyrms,uzrms,VORTXrms,VORTYrms,VORTZrms,                          &
          !         34        35
          &         diss_Temp,diss_Chem,                                                   &          
          !         36        37       38        39          40        41
          &         flux_Part,Part_min,Part_max,up_min(1),up_max(1),up_min(2),             &
		  !            42        43        44         45           46
		  &         up_max(2),up_min(3),up_max(3),VORTp_min(1),VORTp_max(1),               &
		  !            47            48          49           50
		  &         VORTp_min(2),VORTp_max(2),VORTp_min(3),VORTp_max(3),                   &
		  !            51         52           53
		  &         up_max_abs,VORTp_max_abs,diss_Part,                                    &
		  !            54    55      56      57     58     59     60        61
		  &         uprms,VORTprms,PARTrms,upxrms,upyrms,upzrms,VORTpXrms,VORTpYrms,       &
		  !            62      63   64         65
		  &         VORTpZrms,err,err_part, diss_mom                                        
  ENDIF

  IF (myid.EQ.0) THEN                  ! Use for testing routines which compute power spectra
     PRINT*,"urms**2=",urms**2
     PRINT*,"Temp_rms**2=",TEMPrms**2
     PRINT*,"Chem_rms**2_",CHEMrms**2
	 PRINT*,"Part_rms**2_",PARTrms**2
  ENDIF

END SUBROUTINE write_diagnostics_file
