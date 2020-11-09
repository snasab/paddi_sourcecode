! Compute the RHS of the velocity equation
SUBROUTINE crhs_velocity(rhs,u,Temp,Chem)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity
  USE parameter_module, ONLY: Nx,Nmax,kx,ky,kz,B_therm,B_comp
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  TYPE(buoyancy) :: Temp,Chem
  TYPE(velocity) :: u
  COMPLEX(kind=kr) ::  rhs(0:,mysx_spec:,mysy_spec:,1:) ! RHS 
  REAL (kind=kr),POINTER     :: work_phys(:,:,:)
  REAL (kind=kr) :: hkx,hky,hkz
  REAL (kind=kr) :: ksquare
  COMPLEX (kind=kr) :: kN,fac
  INTEGER (kind=ki) :: i,j,k

! compute Fourier coeff. of -curl(u) \times u and add buoyancy force if present
  ALLOCATE(work_phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))

#ifdef TWO_DIMENSIONAL
  work_phys = - u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_z) 
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))
  ! z-component 
  work_phys =   u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_x) 
#else
  ! x-component
  work_phys = - u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_z) &
            & + u%curl(:,:,:,curl_z) * u%phys(:,:,:,vec_y) 
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))
  ! y-component 
  work_phys = - u%curl(:,:,:,curl_z) * u%phys(:,:,:,vec_x) &
            & + u%curl(:,:,:,curl_x) * u%phys(:,:,:,vec_z) 
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_y))
  ! z-component 
  work_phys = - u%curl(:,:,:,curl_x) * u%phys(:,:,:,vec_y) &
            & + u%curl(:,:,:,curl_y) * u%phys(:,:,:,vec_x) 
#endif

#ifdef TEMPERATURE_FIELD
  work_phys = work_phys + B_therm * Temp%phys 
#endif
#ifdef CHEMICAL_FIELD
  work_phys = work_phys - B_comp * Chem%phys
#endif

  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_z))
  DEALLOCATE(work_phys)

! add contribution from pressure gradient 
  ! This does not look too cache efficient...
#ifdef TWO_DIMENSIONAL
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     hkx = kx(i)
     DO k=0,2*Nmax-1
        hkz = kz(k)
        kN =   hkx * rhs(k,i,j,vec_x) &
             + hkz * rhs(k,i,j,vec_z)
        ksquare = hkx**2 + hkz**2
        ksquare = MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later...
        fac = kN / ksquare
        rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - hkx * fac
        rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - hkz * fac
     ENDDO
  ENDDO
#else
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           hkz = kz(k)
           kN =   hkx * rhs(k,i,j,vec_x) &
           &    + hky * rhs(k,i,j,vec_y) &
           &    + hkz * rhs(k,i,j,vec_z)
           ksquare = hkx**2 + hky**2 + hkz**2
           ksquare = MAX(ksquare,EPSILON(1._kr)) ! Avoid floating exception in 1/ksquare later...
           fac = kN / ksquare
           rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - hkx * fac
           rhs(k,i,j,vec_y) =   rhs(k,i,j,vec_y) - hky * fac
           rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - hkz * fac
        ENDDO
     ENDDO
  ENDDO
#endif
  ! Set to zero constant part. 
  IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0) rhs(0,0,0,:) = (0._kr,0._kr)

END SUBROUTINE crhs_velocity
