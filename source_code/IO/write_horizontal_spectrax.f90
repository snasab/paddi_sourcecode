subroutine write_horizontal_spectra(u,Temp,Chem,up,Part,istep,t)
  use defprecision_module
  use state_module, ONLY: buoyancy, velocity, ltime0
  use message_passing_module, ONLY : myid
  use parameter_module, ONLY: kx,ky,Lmax,Mmax
  USE mpi_transf_module,  ONLY:  mysx_spec,myex_spec,mysy_spec,myey_spec
  implicit none
  type(velocity)   :: u,up
  type(buoyancy)   :: Temp,Chem,Part
  integer(kind=ki) :: istep
  real(kind=kr)    :: t
  REAL(kind=kr), POINTER   :: Ener_spec_u(:,:,:),Ener_spec_up(:,:,:)
  REAL(Kind=kr), POINTER   :: Ener_spec_Temp(:,:),Ener_spec_Chem(:,:),Ener_spec_Part(:,:)
  REAL(kind=kr), POINTER   :: spec_local(:,:)
  REAL(kind=kr)    :: Ekin,ETemp,EChem,EPart,Ekinp
  integer(kind=ki) :: l,m
  ! allocate arrays for the horizontal energy spectra
  ALLOCATE(Ener_spec_u(0:Lmax,0:max(2*Mmax-1,0),3))
  ALLOCATE(Ener_spec_up(0:Lmax,0:max(2*Mmax-1,0),3))
  ALLOCATE(Ener_spec_Temp(0:Lmax,0:max(2*Mmax-1,0)))
  ALLOCATE(Ener_spec_Chem(0:Lmax,0:max(2*Mmax-1,0)))
  ALLOCATE(Ener_spec_Part(0:Lmax,0:max(2*Mmax-1,0)))
  ALLOCATE(spec_local(mysx_spec:myex_spec,mysy_spec:myey_spec))
  
  ! compute spectra
  spec_local = horizontal_power_spectrum(u%spec(:,:,:,vec_x,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_u(:,:,1))
#ifdef TWO_DIMENSIONAL
  Ener_spec_u(:,:,2)=0._kr
#else
  spec_local = horizontal_power_spectrum(u%spec(:,:,:,vec_y,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_u(:,:,2))
#endif
  spec_local = horizontal_power_spectrum(u%spec(:,:,:,vec_z,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_u(:,:,3))

#ifdef TEMPERATURE_FIELD
  spec_local = horizontal_power_spectrum(Temp%spec(:,:,:,ltime0))
  call gather_E_spec(spec_local,Ener_spec_Temp)
#else
  Ener_spec_Temp = 0._kr
#endif

#ifdef CHEMICAL_FIELD
  spec_local = horizontal_power_spectrum(Chem%spec(:,:,:,ltime0))
  call gather_E_spec(spec_local,Ener_spec_Chem)
#else
  Ener_spec_Chem = 0._kr
#endif  

#ifdef PARTICLE_FIELD
  ! compute spectra
  spec_local = horizontal_power_spectrum(up%spec(:,:,:,vec_x,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_up(:,:,1))
#ifdef TWO_DIMENSIONAL
  Ener_spec_up(:,:,2)=0._kr
#else
  spec_local = horizontal_power_spectrum(up%spec(:,:,:,vec_y,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_up(:,:,2))
#endif
  spec_local = horizontal_power_spectrum(up%spec(:,:,:,vec_z,ltime0)) 
  call gather_E_spec(spec_local,Ener_spec_up(:,:,3))
  
  spec_local = horizontal_power_spectrum(Part%spec(:,:,:,ltime0))
  call gather_E_spec(spec_local,Ener_spec_Part)
#else
    Ener_spec_Part = 0._kr
#endif

  ! write spectrum
  IF (myid.EQ.0) THEN 
     WRITE (uout(2),'(a,I8,a,E20.7)') '# Timstep =',istep,' time =',t
     Ekin=0._kr ! used for testing: see below
	 Ekinp=0._kr
     Etemp=0._kr
     EChem=0._kr
	 EPart=0._kr
     DO l=0,Lmax
        DO m=0,max(2*Mmax-1,0)
           WRITE(uout(2),'(18E20.7)') kx(l),ky(m),                      &
                & Ener_spec_u(l,m,1) , Ener_spec_u(l,m,2) , Ener_spec_u(l,m,3), &
                & Ener_spec_u(l,m,1) + Ener_spec_u(l,m,2) + Ener_spec_u(l,m,3), &
                & Ener_spec_Temp(l,m),Ener_spec_Chem(l,m),Ener_spec_Part(l,m),  &
                & Ener_spec_up(l,m,1) , Ener_spec_up(l,m,2) , Ener_spec_up(l,m,3), &
                & Ener_spec_up(l,m,1) + Ener_spec_up(l,m,2) + Ener_spec_up(l,m,3)
               Ekin=Ekin+(Ener_spec_u(l,m,1) + Ener_spec_u(l,m,2) + Ener_spec_u(l,m,3))
			   Ekinp=Ekinp+(Ener_spec_up(l,m,1) + Ener_spec_up(l,m,2) + Ener_spec_up(l,m,3))
               Etemp=Etemp+Ener_spec_Temp(l,m)
               Echem=Echem+Ener_spec_Chem(l,m)
			   EPart=EPart+Ener_spec_Part(l,m)
        ENDDO
     ENDDO
     WRITE(*,'(a,E30.16)') "Ekin(horiz.Spec) =",Ekin   ! this might be used for testing:
	 WRITE(*,'(a,E30.16)') "Ekinp(horiz.Spec) =",Ekinp   
     WRITE(*,'(a,E30.16)') "Etemp(horiz.Spec) =",Etemp ! volume avaraged energy? 
     WRITE(*,'(a,E30.16)') "Echem(horiz.Spec) =",Echem ! Compare with the the squared rms values 
	 WRITE(*,'(a,E30.16)') "EPart(horiz.Spec) =",EPart
     !                                                ! computed in routine wrtite_diagnostics_file
     WRITE(uout(2),*) 
     WRITE(uout(2),*)
  ENDIF

  !deallocate memory
  !DEALLOCATE(Ener_spec_u,Ener_spec_Temp,Ener_spec_Chem,spec_local)
  DEALLOCATE(Ener_spec_u,Ener_spec_up,Ener_spec_Temp,Ener_spec_Chem,Ener_spec_Part,spec_local)
  
end subroutine write_horizontal_spectra


! Subroutine to compute local power spectrum 
FUNCTION horizontal_power_spectrum(x)
  use defprecision_module
  use parameter_module, ONLY: Nmax
  use mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec
  implicit none
  complex(kind=kr) :: x(0:,mysx_spec:,mysy_spec:)
  real(kind=kr)    :: horizontal_power_spectrum(mysx_spec:myex_spec,mysy_spec:myey_spec)
  integer(kind=ki) :: i,j,k
  real(kind=kr)    :: zsum,fac

  do j=mysy_spec,myey_spec
     do i=mysx_spec,myex_spec
        zsum = 0._kr  
        do k=0,2*Nmax-1
           fac = REAL((2 - DIM(1,i)),kr) 
           zsum = zsum + x(k,i,j) * conjg(x(k,i,j)) * fac
        enddo
        horizontal_power_spectrum(i,j) = zsum
     enddo
  enddo

end FUNCTION horizontal_power_spectrum
  
  
! Subroutine that collects a distributed 2d array in_local 
! at node with id zero.  
! not optimized for efficiency...
SUBROUTINE gather_E_spec(in_local,out_global)
  USE defprecision_module
  USE message_passing_module, ONLY : myid,numtasks
  USE parameter_module , ONLY : Lmax,Mmax
  USE mpi_transf_module, ONLY : mysx_spec,myex_spec,mysy_spec,myey_spec 
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif

  REAL(kind=kr)    :: in_local(mysx_spec:,mysy_spec:),out_global(0:,0:)
  REAL(kind=kr)    :: work((Lmax+1)*(max(2*Mmax,1)))
  INTEGER(kind=ki) :: sendcount
  INTEGER(kind=ki) :: sx_spec(0:numtasks-1),ex_spec(0:numtasks-1)
  INTEGER(kind=ki) :: sy_spec(0:numtasks-1),ey_spec(0:numtasks-1)
  INTEGER(kind=ki) :: recvcounts(0:numtasks-1)
  INTEGER(kind=ki) :: rdispls(0:numtasks-1)
  INTEGER(kind=ki) :: count,process,i,j,ierr

! node 0 finds out which part of the array the other processes hold
  CALL MPI_GATHER(mysx_spec,1,MPI_INTEGER,sx_spec,1,MPI_INTEGER,0, &
  &               MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(myex_spec,1,MPI_INTEGER,ex_spec,1,MPI_INTEGER,0, &
  &               MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(mysy_spec,1,MPI_INTEGER,sy_spec,1,MPI_INTEGER,0, &
  &               MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(myey_spec,1,MPI_INTEGER,ey_spec,1,MPI_INTEGER,0, &
  &               MPI_COMM_WORLD,ierr)

! now gather vertically integrated spectrum at node 0 in array work
  sendcount  = (myey_spec-mysy_spec+1) * (myex_spec-mysx_spec+1)
  IF (myid.EQ.0) THEN 
     recvcounts(:) = (ey_spec(:)-sy_spec(:)+1) * (ex_spec(:)-sx_spec(:)+1) 
     count=1
     DO process = 0,numtasks-1
        rdispls(process) = count - 1 
        count = count +  (ey_spec(process)-sy_spec(process)+1) & 
        &              * (ex_spec(process)-sx_spec(process)+1)
     ENDDO
  ENDIF
     
  CALL mpi_gatherv(in_local,sendcount,PM_MPI_FLOAT_TYPE,      &
  &                work,recvcounts,rdispls,PM_MPI_FLOAT_TYPE, &
  &                0,MPI_COMM_WORLD,ierr)
  
! sort stuff into array out_global
  IF (myid.EQ.0) THEN 
     count=1
     DO process=0,numtasks-1
        DO j=sy_spec(process),ey_spec(process)
           DO i=sx_spec(process),ex_spec(process)
              out_global(i,j) = work(count)
              count=count+1
           ENDDO
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE gather_E_spec
