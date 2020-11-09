SUBROUTINE compute_phys_space_vars(u,Temp,Chem,up,Part)
  USE defprecision_module
  USE parameter_module, ONLY: Nmax,Nx,Ny,Nz,Gammax,Gammay,Gammaz
  USE state_module, ONLY: velocity,buoyancy,ltime0,shift_time_pointers
  USE minit_module, ONLY: M_init,compute_particle_mass
  USE mpi_transf_module, ONLY : FFT_c2r,FFT_r2c,mysx_spec,myex_spec,mysy_spec,myey_spec &
      & ,mysy_phys,myey_phys,mysz_phys,myez_phys
  USE diff_op_module, ONLY: curl,d_by_dx,d_by_dz
  USE MPI
  IMPLICIT NONE
  INTEGER(kind=ki) :: ierr
  
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  COMPLEX(kind=kr), POINTER :: work(:,:,:,:)
  INTEGER(kind=ki) :: i,j,k
  REAL(kind=kr) :: M_curr

! Transform Temp,Chem,Part,u to physical space
#ifdef TEMPERATURE_FIELD
  CALL FFT_c2r(Temp%spec(:,:,:,ltime0),Temp%phys)
#endif
#ifdef CHEMICAL_FIELD
  CALL FFT_c2r(Chem%spec(:,:,:,ltime0),Chem%phys)
#endif

  CALL FFT_c2r(u%spec(:,:,:,vec_x,ltime0),u%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(u%spec(:,:,:,vec_y,ltime0),u%phys(:,:,:,vec_y))
#endif
  CALL FFT_c2r(u%spec(:,:,:,vec_z,ltime0),u%phys(:,:,:,vec_z))
! compute curl(u) in spectral space and transform to physical space
#ifdef TWO_DIMENSIONAL
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,curl_y:curl_y))
  work(:,:,:,curl_y) = d_by_dz(u%spec(:,:,:,vec_x,ltime0)) - d_by_dx(u%spec(:,:,:,vec_z,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),u%curl(:,:,:,curl_y))
#else
  ALLOCATE(work(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,curl_x:curl_z))
  work = curl(u%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,vec_x),u%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,vec_y),u%curl(:,:,:,curl_y))
  CALL FFT_c2r(work(:,:,:,vec_z),u%curl(:,:,:,curl_z))
#endif

#ifdef PARTICLE_FIELD
  CALL FFT_c2r(Part%spec(:,:,:,ltime0),Part%phys)
  do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
         do i=0,Nx-1
           Part%phys(i,j,k) = max(0._kr, Part%phys(i,j,k))
         enddo
     enddo
  enddo
CALL compute_particle_mass(M_curr,Part)
M_init = Gammax*Gammay*Gammaz
Part%phys = Part%phys * M_init/M_curr
!!CHECK CALLING SEQUENCE FOR R2C
call FFT_r2c(Part%phys,Part%spec(:,:,:,ltime0))

  
  CALL FFT_c2r(up%spec(:,:,:,vec_x,ltime0),up%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  CALL FFT_c2r(up%spec(:,:,:,vec_y,ltime0),up%phys(:,:,:,vec_y))
#endif
  CALL FFT_c2r(up%spec(:,:,:,vec_z,ltime0),up%phys(:,:,:,vec_z))
! compute curl(up) in spectral space and transform to physical space
#ifdef TWO_DIMENSIONAL
  work(:,:,:,curl_y) = d_by_dz(up%spec(:,:,:,vec_x,ltime0)) - d_by_dx(up%spec(:,:,:,vec_z,ltime0))
  CALL FFT_c2r(work(:,:,:,curl_y),up%curl(:,:,:,curl_y))
#else
  work = curl(up%spec(:,:,:,:,ltime0))
  CALL FFT_c2r(work(:,:,:,vec_x),up%curl(:,:,:,curl_x))
  CALL FFT_c2r(work(:,:,:,vec_y),up%curl(:,:,:,curl_y))
  CALL FFT_c2r(work(:,:,:,vec_z),up%curl(:,:,:,curl_z))
#endif
#endif

  DEALLOCATE(work)
END SUBROUTINE compute_phys_space_vars
