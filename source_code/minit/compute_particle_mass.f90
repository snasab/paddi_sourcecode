subroutine compute_particle_mass(mass,Part)
  use defprecision_module
  use parameter_module, ONLY : Nx,Ny,Nz,Gammax,Gammay,Gammaz
  USE mpi_transf_module, ONLY:  mysy_phys,myey_phys,mysz_phys,myez_phys
  USE state_module, ONLY: buoyancy
  use MPI
  implicit none
  real(kind=kr)  :: mass
  type(buoyancy) :: Part
  real (kind=kr) :: lmass,dx,dy,dz,dv!,volume
  integer(kind=ki) :: i,j,k,ierr
  dx = Gammax / Nx
  dz = Gammaz / Nz
#ifdef TWO_DIMENSIONAL
  dv = dx*dz
  !volume = Gammax * Gammaz
#else
  dy = Gammay / Ny
  dv = dx*dy*dz
  !volume = Gammax * Gammay * Gammaz
#endif

  lmass = 0._kr
  do k=mysz_phys,myez_phys
     do j=mysy_phys,myey_phys
        do i=0,Nx-1
           lmass = lmass + dv*Part%phys(i,j,k)
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(lmass,mass,1,PM_MPI_FLOAT_TYPE,MPI_SUM,MPI_COMM_WORLD,ierr)
  
end subroutine compute_particle_mass
