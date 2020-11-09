subroutine deallocate_uTC(u,Temp,Chem,up,Part)
  implicit none
  type(velocity) :: u,up
  type(buoyancy) :: Temp, Chem, Part

  deallocate(u%phys,u%spec,u%curl)
#ifdef TEMPERATURE_FIELD
  deallocate(Temp%phys,Temp%spec)
#endif
#ifdef CHEMICAL_FIELD
  deallocate(Chem%phys,Chem%spec)
#endif
#ifdef PARTICLE_FIELD
  deallocate(Part%phys,Part%spec)
  deallocate(up%phys,up%spec,up%curl)
#endif

#ifdef AB_BDF3
  deallocate(u%rhs,Temp%rhs,Chem%rhs,up%rhs,Part%rhs)
#endif

end subroutine deallocate_uTC
