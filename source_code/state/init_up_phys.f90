SUBROUTINE init_up_phys(up)
  USE defprecision_module
  USE parameter_module, ONLY: Nx,Ny,Nz,Gammax,Gammay,Gammaz,pi,G_part,T_part,rn_using_myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  USe message_passing_module, ONLY: myid
  IMPLICIT NONE
  TYPE(velocity) :: up
  INTEGER(kind=ki) :: i,j,k,idum
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

IF (rn_using_myid) THEN
  idum = myid
ELSE
  idum = -9
ENDIF

  DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           xc = i*dx
#ifdef TWO_DIMENSIONAL
           ! Initial velocity for 2D case 
           up%phys(i,j,k,vec_x) =  0._kr 
           up%phys(i,j,k,vec_z) =  0._kr
#else
           ! Initial velocity for 3D case 
           up%phys(i,j,k,vec_x) =  0._kr 
           up%phys(i,j,k,vec_y) =  0._kr 
           up%phys(i,j,k,vec_z) =  0._kr
#endif
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE init_up_phys
