SUBROUTINE init_Part_phys(Part)
  USE defprecision_module
  USE parameter_module,ONLY:Nx,Ny,Nz,Gammax,Gammay,Gammaz,alpha,beta,gamma,&
                            &G_part,D_visc,B_therm,R_part,rn_using_myid
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys
  USE message_passing_module, ONLY: myid
 
  IMPLICIT NONE
  INTEGER(kind=ki) :: ierr
  
  TYPE(buoyancy) :: Part
  INTEGER(kind=ki) :: i,j,k,idum
  REAL(kind=kr) :: xc,yc,zc,dx,dy,dz,rn,sigma,amplitude,dens_grad,temp_grad
  REAL(kind=kr),POINTER :: r(:,:,:)
  REAL(kind=kr) :: start,finish,rn_myid

  dx = Gammax / Nx
  dy = Gammay / Ny
  dz = Gammaz / Nz 

IF (rn_using_myid) THEN
  idum = myid 
ELSE
  idum = -8
ENDIF

Part%phys = 1.0
go to 50
DO k=mysz_phys,myez_phys
     zc = k*dz
     DO j=mysy_phys,myey_phys
        yc = j*dy
        DO i=0,Nx-1
           !xc = i*dx
           !rn = decomp_independent_random(i,j,k,idum) - 0.5_kr
           !Part%phys(i,j,k) = (1+ 1.E-2*rn)*exp(-(zc-gammaz/2)**2/(2*sigma**2))
           !Part%phys(i,j,k) = 50.*(D_visc/G_part)*exp(-(zc-gammaz/2)**2/(2*30.**2))
           Part%phys(i,j,k) = 1.0
        ENDDO
     ENDDO
ENDDO
50 continue

 
END SUBROUTINE init_Part_phys
