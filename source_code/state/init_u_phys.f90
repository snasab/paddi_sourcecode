SUBROUTINE init_u_phys(u)
  USE defprecision_module
  USE parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,pi,alpha,beta,gamma
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  IMPLICIT NONE
  TYPE(velocity) :: u
  INTEGER(kind=ki) :: i,j,k
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

  DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
#ifdef TWO_DIMENSIONAL
           ! Initial velocity for 2D case 
           u%phys(i,j,k,vec_x) =  0._kr 
           u%phys(i,j,k,vec_z) =  0._kr 
#else
           ! Initial velocity for 3D case 
           u%phys(i,j,k,vec_x) =  0._kr 
           u%phys(i,j,k,vec_y) =  0._kr 
           u%phys(i,j,k,vec_z) =  0._kr 
#endif
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE init_u_phys
