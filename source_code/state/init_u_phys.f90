SUBROUTINE init_u_phys(u)
  USE defprecision_module
  USE parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,pi,alpha,beta,gamma,rn_using_myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY: myid
  IMPLICIT NONE
  TYPE(velocity) :: u
  INTEGER(kind=ki) :: i,j,k,idum,rand,vec
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz,rn
  REAL(kind=kr) :: start, finish
  
  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

IF (rn_using_myid) THEN 
  idum = myid
ELSE
  idum = -5
ENDIF

DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
#ifdef TWO_DIMENSIONAL
           ! Initial velocity for 2D case
           rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           u%phys(i,j,k,vec_x) = (1.E-2)*rn
           u%phys(i,j,k,vec_z) = (1.E-2)*rn !(1.E-7)*rn*sin(xc+zc)
#else
           ! Initial velocity for 3D case 
           rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           u%phys(i,j,k,vec_x) =  (1.E-2)*rn
           u%phys(i,j,k,vec_y) =  (1.E-2)*rn
           u%phys(i,j,k,vec_z) =  (1.E-2)*rn
#endif
        ENDDO
     ENDDO
ENDDO


END SUBROUTINE init_u_phys
