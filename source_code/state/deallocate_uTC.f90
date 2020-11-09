subroutine deallocate_uTC(u,Temp,Chem)
  implicit none
  type(velocity) :: u
  type(buoyancy) :: Temp, Chem

  deallocate(u%phys,u%spec,u%curl)
#ifdef TEMPERATURE_FIELD
  deallocate(Temp%phys,Temp%spec)
#endif
#ifdef CHEMICAL_FIELD
  deallocate(Chem%phys,Chem%spec)
#endif

#ifdef AB_BDF3
  deallocate(u%rhs,Temp%rhs,Chem%rhs)
#endif

end subroutine deallocate_uTC
