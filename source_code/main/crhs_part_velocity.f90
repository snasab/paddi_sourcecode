! Compute the RHS of the velocity equation
SUBROUTINE crhs_part_velocity(rhs,u,up,Part,drag)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity
  USE parameter_module, ONLY: Nx,Nmax,kx,ky,kz,T_part,G_part
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        FFT_r2c
  IMPLICIT NONE
  TYPE(buoyancy) :: Part
  TYPE(velocity) :: u,up
  REAL(kind=kr) :: drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec) 
  COMPLEX(kind=kr) ::  rhs(0:,mysx_spec:,mysy_spec:,1:) ! RHS 
  COMPLEX (kind=kr), POINTER :: work_spec(:,:,:)
  REAL (kind=kr),POINTER     :: work_phys(:,:,:)
  REAL (kind=kr) :: hkx,hky,hkz
  REAL (kind=kr) :: ksquare
  COMPLEX (kind=kr) :: kN,fac
  INTEGER (kind=ki) :: i,j,k
  COMPLEX(kind=kr), PARAMETER :: iu = (0._kr,1._kr)

! compute Fourier coeff. of -curl(u) \times u and add buoyancy force if present
  ALLOCATE(work_phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  ALLOCATE(work_spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  
#ifdef TWO_DIMENSIONAL
  work_phys = - up%curl(:,:,:,curl_y) * up%phys(:,:,:,vec_z) &
  			  & + drag(:,:,:,vec_x) 
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))
  ! z-component 
  work_phys =   up%curl(:,:,:,curl_y) * up%phys(:,:,:,vec_x) &
                & + drag(:,:,:,vec_z) - G_part 
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_z))
#else
  ! x-component
  work_phys = - up%curl(:,:,:,curl_y) * up%phys(:,:,:,vec_z) &
            & + up%curl(:,:,:,curl_z) * up%phys(:,:,:,vec_y) &
			& + drag(:,:,:,vec_x)
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_x))
  ! y-component 
  work_phys = - up%curl(:,:,:,curl_z) * up%phys(:,:,:,vec_x) &
            & + up%curl(:,:,:,curl_x) * up%phys(:,:,:,vec_z) &
			& + drag(:,:,:,vec_y)
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_y))
  ! z-component 
  work_phys = - up%curl(:,:,:,curl_x) * up%phys(:,:,:,vec_y) &
            & + up%curl(:,:,:,curl_y) * up%phys(:,:,:,vec_x) &
			& + drag(:,:,:,vec_z) - G_part
  CALL FFT_r2c(work_phys,rhs(:,:,:,vec_z))
#endif

!u_x^2 + u_z^2 (+ u_y^2)
work_phys = (up%phys(:,:,:,vec_x)*up%phys(:,:,:,vec_x)) + (up%phys(:,:,:,vec_z)*up%phys(:,:,:,vec_z))
#ifndef TWO_DIMENSIONAL
work_phys = work_phys + (up%phys(:,:,:,vec_y)*up%phys(:,:,:,vec_y))
#endif 
CALL FFT_r2c(work_phys,work_spec)

#ifdef TWO_DIMENSIONAL
  j=mysy_spec !=myey_spec=0
  DO i=mysx_spec,myex_spec
     hkx = kx(i)
     DO k=0,2*Nmax-1
        hkz = kz(k)
        rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - 0.5 * iu * hkx * work_spec(k,i,j)
        rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - 0.5 * iu * hkz * work_spec(k,i,j)
     ENDDO
  ENDDO
#else
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           hkz = kz(k)
           rhs(k,i,j,vec_x) =   rhs(k,i,j,vec_x) - 0.5 * iu * hkx * work_spec(k,i,j)
           rhs(k,i,j,vec_y) =   rhs(k,i,j,vec_y) - 0.5 * iu * hky * work_spec(k,i,j)
           rhs(k,i,j,vec_z) =   rhs(k,i,j,vec_z) - 0.5 * iu * hkz * work_spec(k,i,j)
        ENDDO
     ENDDO
  ENDDO
#endif
  
  DEALLOCATE(work_phys,work_spec)

  ! Set to zero constant part. 
  !IF (mysx_spec.EQ.0 .AND. mysy_spec.EQ.0) rhs(0,0,0,:) = (0._kr,0._kr)

END SUBROUTINE crhs_part_velocity
