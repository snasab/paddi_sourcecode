subroutine allocate_uTC(u,Temp,Chem,up,Part)
! +------------------------------------------------------------+
! | Allocate memory needed to store the state variables and     |
! | related quantities.                                        |
! +------------------------------------------------------------+
  use defprecision_module
  use parameter_module, ONLY:Nx,Ny,Nz,Lmax,Mmax,Nmax
  use mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  implicit none
  type(velocity) :: u, up
  type(buoyancy) :: Temp,Chem,Part

#ifdef AB_BDF3 
!Adams-Bashforth-Backward_Differencing-3rd-order ImEx method
  allocate(u%phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec))
  allocate(u%curl(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_curl))
  allocate(u%spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_vec,time_levels)) 
  allocate( u%rhs(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_vec,time_levels-1))

#ifdef TEMPERATURE_FIELD
  allocate(Temp%phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  allocate(Temp%spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels)) 
  allocate( Temp%rhs(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels-1))
#endif

#ifdef CHEMICAL_FIELD
  allocate(Chem%phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  allocate(Chem%spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels)) 
  allocate( Chem%rhs(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels-1))
#endif

#ifdef PARTICLE_FIELD
  allocate(Part%phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys))
  allocate(Part%spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels)) 
  allocate(Part%rhs(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,time_levels-1))
  allocate(up%phys(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec))
  allocate(up%curl(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_curl))
  allocate(up%spec(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_vec,time_levels)) 
  allocate(up%rhs(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,dim_vec,time_levels-1))
#endif

#endif

end subroutine allocate_uTC
