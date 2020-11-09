SUBROUTINE init_random_seed()
            USe message_passing_module, ONLY: myid
            INTEGER :: i, n, clock, idum
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            idum = -5

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            !CALL SYSTEM_CLOCK(COUNT=clock)

            seed = myid !clock + 37 * (/ (i - 1, i = 1, n) /)

            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
END SUBROUTINE

