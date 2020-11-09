SUBROUTINE set_initial_condition(u,Temp,Chem,up,Part)
  USE defprecision_module
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec,mysz_spec,FFT_r2c,FFT_c2r
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY: Nmax 
  USE diff_op_module, ONLY: curl
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  COMPLEX(kind=kr), POINTER :: work(:,:,:,:)
  real :: start, finish
  
! transform to spectral space 

#ifdef TEMPERATURE_FIELD
  Temp%spec=(0._kr,0._kr)
  CALL FFT_r2c(Temp%phys,Temp%spec(:,:,:,ltime0))
#endif

#ifdef CHEMICAL_FIELD
  Chem%spec=(0._kr,0._kr)
  CALL FFT_r2c(Chem%phys,Chem%spec(:,:,:,ltime0))
#endif

#ifdef PARTICLE_FIELD
  Part%spec=(0._kr,0._kr)
  CALL FFT_r2c(Part%phys,Part%spec(:,:,:,ltime0))
#endif

  u%spec=(0._kr,0._kr)
  up%spec=(0._kr,0._kr)

#ifdef TWO_DIMENSIONAL
if (myid .EQ. 7) print*, "START FFT r2c uphys"
if (myid .EQ. 7) call cpu_time(start) 
  CALL FFT_r2c(u%phys(:,:,:,vec_x),u%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_z),u%spec(:,:,:,vec_z,ltime0))
#else 
  CALL FFT_r2c(u%phys(:,:,:,vec_x),u%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_y),u%spec(:,:,:,vec_y,ltime0))
  CALL FFT_r2c(u%phys(:,:,:,vec_z),u%spec(:,:,:,vec_z,ltime0))
#endif
if (myid .EQ. 7) print*, "END FFT r2c uphys"
if (myid .EQ. 7) call cpu_time(finish)
if (myid .EQ. 7) print*, "FFT TIME: ", finish-start 

#ifdef PARTICLE_FIELD
#ifdef TWO_DIMENSIONAL
  CALL FFT_r2c(up%phys(:,:,:,vec_x),up%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(up%phys(:,:,:,vec_z),up%spec(:,:,:,vec_z,ltime0))
#else 
  CALL FFT_r2c(up%phys(:,:,:,vec_x),up%spec(:,:,:,vec_x,ltime0))
  CALL FFT_r2c(up%phys(:,:,:,vec_y),up%spec(:,:,:,vec_y,ltime0))
  CALL FFT_r2c(up%phys(:,:,:,vec_z),up%spec(:,:,:,vec_z,ltime0))
#endif
#endif  

! Make initial velocity solenoidal and compute phsical space values of corrected velocity
 CALL make_solenoidal(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(u%spec(:,:,:,vec_x,ltime0),u%phys(:,:,:,vec_x))
  CALL FFT_c2r(u%spec(:,:,:,vec_z,ltime0),u%phys(:,:,:,vec_z))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(u%spec(:,:,:,vec_y,ltime0),u%phys(:,:,:,vec_y))
#endif

  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_curl))
  work = curl(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),u%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(work(:,:,:,curl_x),u%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,curl_z),u%curl(:,:,:,curl_z))
#endif
!for particle field: compute curl(up) in spectral space and transform to physical space
#ifdef PARTICLE_FIELD
  work = curl(up%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),up%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(work(:,:,:,curl_x),up%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,curl_z),up%curl(:,:,:,curl_z))
#endif
#endif
  DEALLOCATE(work)
  

! set to zero mean modes
  IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0 .AND. mysz_spec.EQ.0) THEN 
     u%spec(0,0,0,:,:) = (0._kr,0._kr)
#ifdef TEMPERATURE_FIELD
     Temp%spec(0,0,0,:) = (0._kr,0._kr)
#endif
#ifdef CHEMICAL_FIELD
     Chem%spec(0,0,0,:) = (0._kr,0._kr)
#endif
!#ifdef PARTICLE_FIELD
!     Part%spec(0,0,0,:) = (0._kr,0._kr)
! This is the mean particle #, and we don't want this to be zero. 
!#endif
  ENDIF

#ifdef AB_BDF3
! set right hand sides to zero
  u%rhs = (0._kr,0._kr)
#ifdef TEMPERATURE_FIELD
  Temp%rhs = (0._kr,0._kr)
#endif
#ifdef CHEMICAL_FIELD
  Chem%rhs = (0._kr,0._kr)
#endif
#ifdef PARTICLE_FIELD
  up%rhs = (0._kr,0._kr)
  Part%rhs = (0._kr,0._kr)
#endif
#endif

END SUBROUTINE set_initial_condition
