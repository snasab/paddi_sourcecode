!+--------------------------------------------------------------------+
!| This module contains global data M_init                            |
!+--------------------------------------------------------------------+
#include "defs_MPI.h"
MODULE minit_module
  USE defprecision_module

  IMPLICIT NONE

!*********************************************************************
!***                        GLOBAL DATA                            ***
!*********************************************************************
  
  REAL (kind=kr) :: M_init !Calculated in init_Part
CONTAINS

#include "minit/compute_particle_mass.f90"
END MODULE minit_module
