MODULE parameter_module
!+----------------------------------------------------------------------------------------------+
!| This module contains all the parameters necessary for the calculation. 						|
!| The code solves the following equations:                               						|
!|                                                                        						|
!| Du                                                                    						|
!| -- = - \nabla P + D_visc \nabla^2 u + B_therm T e_z - B_comp C e_z + R_part((u_p - u)/T_part)*r   	|
!| Dt         																					|
!|																								|
!| Du_p                                                                    						|
!| ---- = (u_p - u)/T_part - G_part e_z  + Dv_part \nabla^2 up 														|
!| Dt                                                                     +---------------------+
!|                                                                        |
!| DT                                                                     |
!| -- = D_therm \nabla^2 T + S_therm w                                    |
!| Dt                                                                     |
!|                                                                        |
!| DC                                                                     |
!| -- = D_comp \nabla^2 C + S_comp w                                      |
!| Dt                                                                     |
!| 																		  |
!| dr                                                                     |
!| -- = - div (u_p r) + D_part \nabla^2 r - S_part (u_p e_z)                |
!| dt                                                                     |
!|                                                                        |
!| \nabla * u = 0                                                         |
!|                                                                        |
!| in a triply periodic cube.                                             |
!|                                                                        |
!| Coefficients:                                                          |
!|                                                                        |
!|   B_therm, B_comp represent the coefficients in fromt of the buoyancy  |
!|                   terms;                                               |
!|   D_visc,D_therm,D_comp,D_part are the viscous,thermal,compositional,  |
!|                   particulate diffusion parameters;                    |
!|   S_therm, S_comp, S_part are the coefficients resulting from the      |
!|                   background stratification;                           |
!|   T_part represents the stopping time parameter;			  |
!|	 G_part represents the gravity parameter.			  |
!|   R_part scales the drag term accordingly.                             |
!|   Dv_part is the particle viscous diffusion parameter in particle      |
!|       momentum equation.                                               |
!|                                                                        |
!| Implementing the algorithm for these general equations has the         |
!| advantage that switching between different non-dimensionalizations     |
!| becomes a trivial matter.                                              |  
!+------------------------------------------------------------------------+
  USE defprecision_module
  IMPLICIT NONE 
  SAVE
!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************
!
! Parameters which result from user settings
! ------------------------------------------ 
  INTEGER(kind=ki) :: Lmax       ! highest order in x-direction
  INTEGER(kind=ki) :: Mmax       ! highest order in y-direction
  INTEGER(kind=ki) :: Nmax       ! highest irder in Chebyshev expansion
  INTEGER(kind=ki) :: Nx         ! = Number of grid points in x-direction
  INTEGER(kind=ki) :: Ny         ! = Number of grid points in y-direction
  INTEGER(kind=ki) :: Nz         ! = Number of grid points in z-direction
  INTEGER(kind=ki) :: Nsteps     ! Number of time steps 
  REAL(kind=kr)    :: B_therm    ! Thermal buoyancy coefficient
  REAL(kind=kr)    :: B_comp     ! Compositional buoyancy coefficient
  REAL(kind=kr)    :: D_visc     ! Viscous diffusion parameter
  REAL(kind=kr)    :: D_therm    ! Thermal diffusion parameter
  REAL(kind=kr)    :: D_comp     ! Compositional diffusion parameter
  REAL(kind=kr)    :: S_therm    ! Thermal background stratification parameter
  REAL(kind=kr)    :: S_comp     ! Compositional background stratification parameter
  REAL(kind=kr)    :: T_part     ! Particle stopping time parameter
  REAL(kind=kr)    :: G_part     ! Gravity parameter
  REAL(kind=kr)    :: Dv_part     ! Particle viscous diffusion parameter
  REAL(kind=kr)    :: D_part     ! Particle diffusion parameter
  REAL(kind=kr)    :: S_part     ! Particle background stratification parameter
  REAL(kind=kr)    :: R_part     ! Particle scaling coefficient  
  REAL(kind=kr)    :: Gammax     ! length (x-direction) of the box
  REAL(kind=kr)    :: Gammay     ! width (y-direction) of the box
  REAL(kind=kr)    :: Gammaz     ! height (z-direction) of the box
  REAL(kind=kr)    :: dt_initial ! Initial time step
  REAL(kind=kr)    :: dt_max     ! maximum time step
  REAL(kind=kr)    :: CFL_safety_fac ! Safety Factor for CFL-Condition
  REAL(kind=kr)    :: alpha      ! =2 * pi / Gammax
  REAL(kind=kr)    :: beta       ! =2 * pi / Gammay
  REAL(kind=kr)    :: gamma      ! =2 * pi / Gammaz
  ! FFTW wisdom stuff
  CHARACTER(LEN=100):: FFTW_wisdom_file     ! File to store FFTW wisdom
  LOGICAL           :: use_FFTW_wisdom_file ! Should wisdom be imported from this file?
  LOGICAL           :: rn_using_myid        ! Generating random numbers using myid
  INTEGER(kind=ki),PARAMETER :: uwisdom=60  ! File unit of FFTW widsom file
  ! FFT storage scheme
  REAL (kind=kr) , POINTER   :: kx(:),ky(:),kz(:) !Mapping tables - array index to wavenumber

!
! Constants which are frequently needed by the program
! ----------------------------------------------------

  REAL   (kind=kr),PARAMETER :: pi = 4._kr*ATAN(1._kr) ! Pi

CONTAINS
!*********************************************************************
!***                       Subroutines                             ***
!*********************************************************************
  
#include "parameter/allocate_fft_storage_scheme.f90"
#include "parameter/init_fft_storage_scheme.f90"
#include "parameter/deallocate_fft_storage_scheme.f90"

END MODULE parameter_module

