! Given the velocity field u, this subroutine computes an estimate for the CFL-Limit
! for the third order Adams-Bashforth/Backward-Differencing (AB/BDF3) scheme.
! The estimation is based on results for the advection-diffusion equation given in 
! the book of Peyret(2002). 
FUNCTION CFL_AB_BDF3(u,up,Part)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy
  USE parameter_module, ONLY: Nx,Nmax,D_visc,D_therm,D_comp,D_part,T_part,R_part,kx,ky,kz
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec
  USE MPI
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Part
  REAL(kind=kr) :: CFL_AB_BDF3
  REAL(kind=kr) :: lumax(dim_vec),gumax(dim_vec),ku,ksquared,la(4),ga(4), &
                   & lupmax(dim_vec),gupmax(dim_vec),lpartmax, &
		           & gpartmax,kup
  REAL(kind=kr) :: hkx,hky,hkz
  REAL(kind=kr), PARAMETER :: a=5.61_kr, b=0.69
  REAL(kind=kr), PARAMETER :: safety = 1.0
  INTEGER(kind=ki) :: l,i,j,k,ierr
  ! Find the maximm of ux,uy and uz
  lumax  = 0._kr
  DO l=vec_x,vec_z
     DO k=mysz_phys,myez_phys
        DO j=mysy_phys,myey_phys
           DO i=0,Nx-1
              lumax(l)  = MAX(ABS(u%phys(i,j,k,l)),lumax(l))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(lumax,gumax,dim_vec,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
  
#ifdef PARTICLE_FIELD
  lupmax = 0._kr
  lpartmax = 0._kr
  
  DO l=vec_x,vec_z
     DO k=mysz_phys,myez_phys
        DO j=mysy_phys,myey_phys
           DO i=0,Nx-1
	      lupmax(l) = MAX(ABS(up%phys(i,j,k,l)),lupmax(l))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(lupmax,gupmax,dim_vec,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)

  DO k=mysz_phys,myez_phys
        DO j=mysy_phys,myey_phys
           DO i=0,Nx-1
              lpartmax = MAX(ABS(Part%phys(i,j,k)),lpartmax)
           ENDDO
        ENDDO
     ENDDO
  CALL MPI_ALLREDUCE(lpartmax,gpartmax,1,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
  gpartmax = R_part*gpartmax
  if(gpartmax.lt.1) gpartmax = 1.

#endif 
  
  ! The stability region boundary is nearly rectilinear and may be parameterized according
  ! to Peyret (2002). 
  la = 0._kr
  DO j=mysy_spec,myey_spec
     hky = ky(j)
     DO i=mysx_spec,myex_spec
        hkx = kx(i)
        DO k=0,2*Nmax-1
           hkz = kz(k)
#ifdef TWO_DIMENSIONAL
           ku  = ABS(hkx * gumax(vec_x) + hkz*gumax(vec_z))
	   kup = ABS(hkx * gupmax(vec_x) + hkz*gupmax(vec_z))
           ksquared = hkx**2 + hkz**2
#else
           ku  = ABS(hkx * gumax(vec_x) + hky*gumax(vec_y) + hkz*gumax(vec_z))
	   kup = ABS(hkx * gupmax(vec_x) + hky*gupmax(vec_y) + hkz*gupmax(vec_z))
           ksquared = hkx**2 + hky**2 + hkz**2
#endif
           la(1) = MAX(la(1), a * ku - D_visc *ksquared      ) ! velocity
#ifdef TEMPERATURE_FIELD
           la(2) = MAX(la(2), a * ku - D_therm*ksquared      ) ! temperature
#else 
           la(2) = 0._kr
#endif
#ifdef CHEMICAL_FIELD
           la(3) = MAX(la(3), a * ku - D_comp*ksquared       ) ! concentration
#else
           la(3) = 0._kr
#endif
#ifdef PARTICLE_FIELD
           la(4) = MAX(la(4), a * kup - D_part*ksquared       ) ! particle concentration
#else
           la(4) = 0._kr
#endif
        ENDDO
     ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(la,ga,4,PM_MPI_FLOAT_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)

#ifdef PARTICLE_FIELD
  CFL_AB_BDF3 = MIN( a * b / MAX(MAXVAL(ga),EPSILON(1._kr)), safety*T_part/gpartmax)
#else
  CFL_AB_BDF3 = a * b / MAX(MAXVAL(ga),EPSILON(1._kr))
#endif 


  ! The scheme may be unconditionally stable. In this case, return at huge value
  IF (CFL_AB_BDF3.LE.0._kr) CFL_AB_BDF3=1.E20_kr

END FUNCTION CFL_AB_BDF3
