! This subroutine computes the new Fourier coefficients of the velocity field with the 
! AB3/BDF3 method. 
SUBROUTINE tmstp_particle_vel_AB_BDF3(u,up,Part,drag)
  USE defprecision_module
  USE state_module, ONLY: buoyancy,velocity,ltime0,ltime1,ltime2,rtime1,rtime2
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module, ONLY : Lmax,Mmax,Nmax,Nx,kx,ky,kz,Dv_part
  IMPLICIT NONE
  TYPE(buoyancy)   :: Part ! Temperature and Chemical field
  TYPE(velocity)   :: u,up ! velocity field
  REAL(kind=kr)    :: drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec) 
  REAL(kind=kr)    :: fac,kxsquared,kysquared,kzsquared
  INTEGER(kind=ki) :: i,j,k,l,nqfac,nqfacx,nqfacy,nqfacz


! compute term from BDF part   
  up%spec(:,:,:,:,ltime0) = - ta3 * up%spec(:,:,:,:,ltime0) & ! still contains time level n-3
  &                        - ta1 * up%spec(:,:,:,:,ltime1) &
  &                        - ta2 * up%spec(:,:,:,:,ltime2)                         

! add AB3 part 
  up%spec(:,:,:,:,ltime0) =   up%spec(:,:,:,:,ltime0)        &
  &                        + tb3 * up%rhs(:,:,:,:,rtime1)   & ! still contains rhs from level n-3
  &                        + tb2 * up%rhs(:,:,:,:,rtime2)   

  CALL crhs_part_velocity(up%rhs(:,:,:,:,rtime1),u,up,Part,drag) ! compute terms from time level n-1 ...
  up%spec(:,:,:,:,ltime0) =   up%spec(:,:,:,:,ltime0) + tb1 * up%rhs(:,:,:,:,rtime1) !...and add

! compute new velocity field at timelevel n 
  ! Set Nynquist mode to zero. This is done by multiplying with nqfac
  ! which is 0 for Nynquist modes and 1 otherwise. This is a little crytic, but 
  ! presumably more efficient than if statements
  DO l=vec_x,vec_z
     DO j=mysy_spec,myey_spec
#ifdef TWO_DIMENSIONAL
        nqfacy = 1
        kysquared = 0._kr
#else
        nqfacy = ( 1 - (j-2*Mmax*((j-1)/(Mmax)))/Mmax )
        kysquared = ky(j)**2
#endif
        DO i=mysx_spec,myex_spec
           nqfacx = ( 1 - i/Lmax)
           kxsquared = kx(i)**2
           DO k=0,2*Nmax-1
              nqfacz = ( 1 - (k-2*Nmax*((k-1)/(Nmax)))/Nmax )
              kzsquared = kz(k)**2
              nqfac = nqfacx * nqfacy * nqfacz
              !This is to check that the nqfac stuff works...
              !if (nqfac.eq.0. .and. (i.ne.Lmax .and. j.ne.Mmax .and. k.ne.Nmax)) stop
              !if (nqfac.ne.0. .and. (i.eq.Lmax .or. j.eq.Mmax .or. k.eq.Nmax .or. nqfac.ne.1)) stop
              fac = nqfac * 1._kr / ( ta0 + Dv_part * ( kxsquared + kysquared + kzsquared ) )
              up%spec(k,i,j,l,ltime0) = fac * up%spec(k,i,j,l,ltime0) 
           ENDDO
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE tmstp_particle_vel_AB_BDF3
