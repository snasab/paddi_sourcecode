SUBROUTINE init_Temp_phys(Temp)
  USE defprecision_module
  USE parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,alpha,beta,gamma,rn_using_myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY: myid
  IMPLICIT NONE
  TYPE(buoyancy) :: Temp
  INTEGER(kind=ki) :: i,j,k,idum
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz,rn

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

IF (rn_using_myid) THEN
  idum = myid
ELSE
  idum = -7
ENDIF

  DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
           !rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           Temp%phys(i,j,k) = 0._kr
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE init_Temp_phys
