!+-----------------------------------------------------------------------------+
!| The following function computes the L2-norm of the particle velocity,       |
!| part-vorticity, and particle field.                                         |
!+-----------------------------------------------------------------------------+
SUBROUTINE compute_upP_rms(uprms,VORTprms,PARTrms,upxrms,upyrms,upzrms, &
  &                        VORTpXrms,VORTpYrms,VORTpZrms,up,Part)
  USE defprecision_module
  USE state_module, ONLY: velocity, buoyancy,ltime0
  USE message_passing_module, ONLY : myid
  IMPLICIT NONE
  TYPE(velocity)   :: up
  TYPE(buoyancy)   :: Part
  REAL(kind=kr)    :: uprms,VORTprms,PARTrms
  REAL(kind=kr)    :: upxrms,upyrms,upzrms,VORTpXrms,VORTpYrms,VORTpZrms

#ifdef PARTICLE_FIELD 
  ! compute L2-Norm of up
  uprms    = rms(up%spec(:,:,:,:,ltime0)) !urms = rms(up%phys) more expensive 
  IF (myid.EQ.0) PRINT*,"Part. kinetic Energy: = ",uprms**2
#else
  uprms    = 0._kr
  ! compute L2-Norm of Particle
#endif
#ifdef PARTICLE_FIELD
  PARTrms = rms(Part%spec(:,:,:,ltime0)) !  Partrms = rms(Part%phys) more epensive
#else
  PARTrms = 0._kr
#endif

  ! compute L2-Norm of Particle VORTicity
#ifdef PARTICLE_FIELD
#ifdef TWO_DIMENSIONAL
  VORTprms = rms_scalar(up%curl(:,:,:,curl_y)) ! No spectral representation stored. Use physical one. 
#else
  VORTprms = rms_vector(up%curl) ! No spectral representation stored. Use physical one. 
#endif
  ! compute L2-Norm of the three components of velocity
  upxrms = rms(up%spec(:,:,:,vec_x,ltime0)) 
#ifndef TWO_DIMENSIONAL
  upyrms = rms(up%spec(:,:,:,vec_y,ltime0)) 
#else
  upyrms = 0._kr
#endif
  upzrms = rms(up%spec(:,:,:,vec_z,ltime0)) 
  ! compute L2-Norm of the three components of vorticity
#ifndef TWO_DIMENSIONAL
  VORTpXrms = rms(up%curl(:,:,:,curl_x)) 
  VORTpZrms = rms(up%curl(:,:,:,curl_z)) 
#else
  VORTpZrms = 0._kr
  VORTpXrms = 0._kr
#endif
  VORTpYrms = rms(up%curl(:,:,:,curl_y)) 
#else
  VORTprms  = 0._kr
  VORTpXrms = 0._kr
  VORTpYrms = 0._kr
  VORTpZrms = 0._kr
  upxrms    = 0._kr
  upyrms    = 0._kr
  upzrms    = 0._kr
#endif

  


END SUBROUTINE compute_upP_rms
