MODULE partition_mod
   USE par_oce, ONLY: jpni, jpnj, jpnij, jpi, jpj, jpim1, jpjm1, jpij, &
                      jpreci, jprecj, jpk, jpkm1, jperio, jpiglo, jpjglo
   USE dom_oce, ONLY: ln_zco, nbondi, nbondj, nidom, npolj, &
                      nlci, nlcj,  &  ! i- & j-dimss of the local subdomain
                      nldi, nlei,  &  ! first and last indoor i- and j-indexes
                      nldj, nlej,  &  !
                      nlcit, nlcjt,&  ! 
                      nldit, nldjt,&  ! (Unity-indexed) arrays storing above 
                      nleit, nlejt,&  ! values for each domain.
                      nimpp,njmpp, &  ! i- & j-indices for mpp-subdom. left bottom
                      nimppt,njmppt,& ! Unity-indexed arrays storing above 
                                      ! values for each domain.
                      nperio,      &  ! Local periodicity 
                      nwidthmax,   &  ! Width of widest northern domain
                      narea,       &  ! ID of local area (= rank + 1)
                      mbkmax          ! Deepest level above ocean floor
#if defined key_mpp_mpi
   USE lib_mpp,        ONLY: mppsize, mppsync, mpi_comm_opa,                &
                             MAX_FACTORS, xfactors, yfactors, nn_pttrim,    &
                             nn_cpnode, nn_readpart
#endif
   USE lib_mpp,        ONLY: ctl_stop, ctl_warn
   USE in_out_manager, ONLY: numout, lwp
   USE mapcomm_mod,    ONLY: ielb, ieub, mapcomms, pielb, pjelb, pieub, pjeub,&
                             iesub, jesub, jeub, ilbext, iubext, jubext,      &
                             jlbext, pnactive, piesub, pjesub, jelb, pilbext, &
                             piubext, pjlbext, pjubext, nextra,              &
                             nprocp   ! No. of PEs to partition over
   USE iom,            ONLY: wp, jpdom_unknown, iom_open, iom_get, iom_close
   IMPLICIT NONE
   PRIVATE

   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: imask ! Mask used for partitioning 
                                                 ! (1 for ocean, 0 for land)
                                                 ! set in nemogcm.F90
   ! Holds the bottom level of the ocean at each grid point - used for 
   ! trimming halos in z direction
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:), TARGET :: ibotlevel 

   ! Parameters for the cost function used when evaluating different 
   ! possible domain partitions

   !     Mnemonics:
   !      p = process (i.e. core)
   !      n = node
   !      l = length (number of points)
   !      c = count  (number of messages)
   !      i = internal (intra-node, or on-node)
   !      x = external (inter-node, or off-node)

   INTEGER, PARAMETER :: pv_index_overall = 1
   INTEGER, PARAMETER :: pv_index_wet     = 2
   INTEGER, PARAMETER :: pv_index_dry     = 3
   INTEGER, PARAMETER :: pv_index_pli     = 4
   INTEGER, PARAMETER :: pv_index_plx     = 5
   INTEGER, PARAMETER :: pv_index_pci     = 6
   INTEGER, PARAMETER :: pv_index_pcx     = 7
   INTEGER, PARAMETER :: pv_index_nli     = 8
   INTEGER, PARAMETER :: pv_index_nlx     = 9
   INTEGER, PARAMETER :: pv_index_nci     = 10
   INTEGER, PARAMETER :: pv_index_ncx     = 11
   INTEGER, PARAMETER :: pv_index_tli     = 12
   INTEGER, PARAMETER :: pv_index_tlx     = 13
   INTEGER, PARAMETER :: pv_index_tci     = 14
   INTEGER, PARAMETER :: pv_index_tcx     = 15

   INTEGER, PARAMETER :: pv_index_pcomms  = 16
   INTEGER, PARAMETER :: pv_index_ncomms  = 17
   INTEGER, PARAMETER :: pv_index_tcomms  = 18
     
   INTEGER, PARAMETER :: pv_num_scores    = 18
   
   REAL(wp),PARAMETER :: pv_awful = 1.0e20

#define PARTIT_DEBUG

   PUBLIC imask, ibotlevel, smooth_global_bathy, global_bot_level, partition_mask_alloc
   PUBLIC mpp_init3, partition_rk, partition_mca_rk, write_partition_map
   PUBLIC read_partition, write_partition

CONTAINS

   SUBROUTINE partition_mask_alloc(xsize, ysize, ierr)
      !!------------------------------------------------------------------
      !!                  ***  ROUTINE partition_mask_alloc  ***
      !!
      !! Called from nemogcm to allocate the masks that are members of 
      !! this module
      !!
      !!------------------------------------------------------------------
      INTEGER, INTENT(in) :: xsize, ysize
      INTEGER, INTENT(out):: ierr

      ALLOCATE(imask(xsize,ysize), ibotlevel(xsize,ysize), Stat=ierr)

   END SUBROUTINE partition_mask_alloc


   SUBROUTINE mpp_init3()
      !!------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init3  ***
      !!
      !! History:
      !!       Code adapted from POLCOMS code    17/01/2008  A. Porter
      !!
      !!------------------------------------------------------------------
#if defined key_mpp_mpi
!$AGRIF_DO_NOT_TREAT
      USE mpi ! For better interface checking
#endif
      USE exchtestmod, ONLY : mpp_test_comms
      IMPLICIT NONE
      ! Local vars
      INTEGER :: inum                          ! temporary logical unit
      INTEGER :: iproc, jproc                  ! Loop index
      INTEGER :: ierr                          ! Error flag
      INTEGER :: ji, jj
      CHARACTER(LEN=4) :: intStr

      ! Also set original NEMO arrays as they store internal limits of
      ! each domain in local coordinates
      nldit(narea)  = nldi
      nleit(narea)  = nlei
      nlcit(narea)  = nlci
      nimppt(narea) = nimpp
      nldjt(narea)  = nldj
      nlejt(narea)  = nlej
      nlcjt(narea)  = nlcj
      njmppt(narea) = njmpp
      ! Make sure all PEs have all these values
      ! ARPDBG - wrap this MPI in lib_mpp ??
#if defined key_mpp_mpi
      CALL MPI_ALLGATHER(nldi,1,MPI_INTEGER, &
                         nldit,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nlei,1,MPI_INTEGER, &
                         nleit,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nlci,1,MPI_INTEGER, &
                         nlcit,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nimpp,1,MPI_INTEGER, &
                         nimppt,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nldj,1,MPI_INTEGER, &
                         nldjt,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nlej,1,MPI_INTEGER, &
                         nlejt,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(nlcj,1,MPI_INTEGER, &
                         nlcjt,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(njmpp,1,MPI_INTEGER, &
                         njmppt,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(iesub,1,MPI_INTEGER, &
                         piesub,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(ielb,1,MPI_INTEGER, &
                         pielb,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(ieub,1,MPI_INTEGER, &
                         pieub,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(jesub,1,MPI_INTEGER, &
                         pjesub,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(jelb,1,MPI_INTEGER, &
                         pjelb,1,MPI_INTEGER,mpi_comm_opa,ierr)
      CALL MPI_ALLGATHER(jeub,1,MPI_INTEGER, &
                         pjeub,1,MPI_INTEGER,mpi_comm_opa,ierr)
#endif

      IF(lwp)THEN
         ! Write out domains in postscript

         OPEN(UNIT=997, FILE="domain_decomp.ps", &
              STATUS='REPLACE', ACTION='WRITE', IOSTAT=iproc)

         IF(iproc .EQ. 0)THEN ! Check file opened OK

            ! Header of ps file
            WRITE (997,FMT='("%!PS-Adobe-1.0")')
            WRITE (997,FMT='("/Helvetica findfont %load the font dictionary")')
            WRITE (997,FMT='("12 scalefont        %scale to 12pt")')
            WRITE (997,FMT='("setfont             %make this the current font")')
            WRITE (997,FMT='("/u { 2 mul } def    %set up a scaling factor")')

            ! Put green circles at dry points
            WRITE (997,FMT='("% Filled green circles at dry points")')
            WRITE (997,FMT='("0.1 setlinewidth")') ! Thin green line
            WRITE (997,FMT='("0 1 0 setrgbcolor")')
            WRITE (997,FMT='("newpath")')
            DO jj = 1,jpjglo,1
               DO ji = 1,jpiglo,1
                  IF(imask(ji,jj) /= 1)THEN
                     WRITE (997,FMT='(I3," u ",I3," u 1 u 0 360 arc closepath")') ji, jj
                     ! Use 'fill' instead of 'stroke' to get filled circles
                     WRITE (997,FMT='("fill")')
                  END IF
               END DO
            END DO

            ! Draw the outline of the global domain in red
            WRITE (997,FMT='("% Draw the outline of the global domain in red")')
            WRITE (997,FMT='("3.0 setlinewidth")') ! Fat line of 3mm for global dom.
            WRITE (997,FMT='("1 0 0 setrgbcolor")')
            WRITE (997,FMT='("newpath")')
            WRITE (997,FMT='("% Cursor initialization")')
            WRITE (997,FMT='(I3," u ",1x,I3," u moveto")') 1,1
            WRITE (997,FMT='("% Drawing the rectangle")')
            WRITE (997,FMT='(I3," u ",1x,I3," u lineto")') jpiglo,1
            WRITE (997,FMT='(I3," u ",1x,I3," u lineto")') jpiglo,jpjglo
            WRITE (997,FMT='(I3," u ",1x,I3," u lineto")') 1,jpjglo
            WRITE (997,FMT='("closepath stroke")')

            ! Now draw the outline of each individual domain, nprocp should have been
            ! set at the very beginning of the recursive partioning process.
            WRITE (997,FMT='("0.7 setlinewidth")') ! Back to default width
            WRITE (997,FMT='("0 0 0 setrgbcolor")')! and colour
            DO iproc=1,nprocp
               WRITE (997,FMT='("newpath")')
               WRITE (997,FMT='("% Cursor initialization")')
               WRITE (997,FMT='(I3," u ",1x,I3," u moveto %BL Corner")') pielb(iproc),pjelb(iproc)
               WRITE (997,FMT='("% Drawing the rectangle")')
               WRITE (997,FMT='(I3," u ",1x,I3," u lineto %BR Corner")') pieub(iproc),pjelb(iproc)
               WRITE (997,FMT='(I3," u ",1x,I3," u lineto %TR Corner")') pieub(iproc),pjeub(iproc)
               WRITE (997,FMT='(I3," u ",1x,I3," u lineto %TL Corner")') pielb(iproc),pjeub(iproc)
               WRITE (997,FMT='("closepath stroke")')
               WRITE (997,FMT='(I3," u ",1x,I3," u moveto")') &
                    INT(0.5*(pieub(iproc)+pielb(iproc))), &
                    INT(0.5*(pjeub(iproc)+pjelb(iproc)-1))
               ! Write processor label, left justified
               WRITE (intStr, FMT='(I4)') iproc-1
               ! Use postscipt operator 'stringwidth' to get the approximate width of the
               ! string containing the PE id and then adjust position so it is centred
               ! about point we've just moveto'd. This doesn't attempt to deal with
               ! vertical centering.
               WRITE (997,FMT='("(",(A),") dup stringwidth pop 2 div neg 0 rmoveto show")') TRIM(ADJUSTL(intStr))
            END DO
            WRITE (997,FMT='("showpage")')
            WRITE (997,FMT='("%%EOF")')
            CLOSE(997)

          END IF ! File opened OK
       END IF ! lwp

       ! Test for overlaps of domains (there shouldn't be any!)
       DO iproc=1, nprocp,1
          DO jproc=iproc+1, nprocp, 1
             IF( ( ( (pielb(iproc) .LE. pieub(jproc)) .AND.  &
                    (pielb(iproc) .GE. pielb(jproc)) )      &
                    .OR.                                    &
                  ( (pieub(iproc) .LE. pieub(jproc)) .AND.  &
                    (pieub(iproc) .GE. pielb(jproc)) )  )   &
                .AND.                                       &
                ( ( (pjelb(iproc) .LE. pjeub(jproc)) .AND.  &
                    (pjelb(iproc) .GE. pjelb(jproc)) )      &
                    .OR.                                    &
                  ( (pjeub(iproc) .LE. pjeub(jproc)) .AND.  &
                    (pjeub(iproc) .GE. pjelb(jproc)) )      &
               ) )THEN
                WRITE(*,"('ERROR: domains ',I3,' and ',I3,' overlap!')") &
                      iproc, jproc
                WRITE(*,"(I3,': [',I3,':',I3,'][',I3,':',I3,']')") &
                      iproc, pielb(iproc), pieub(iproc), &
                      pjelb(iproc), pjeub(iproc)
                WRITE(*,"(I3,': [',I3,':',I3,'][',I3,':',I3,']')") &
                      jproc, pielb(jproc), pieub(jproc), pjelb(jproc), pjeub(jproc)

            END IF
         END DO
      END DO

      ! Map out the communications for the partitioned domain.
      CALL mapcomms (imask, ibotlevel, jpiglo, jpjglo, jperio, ierr)
      IF ( ierr.NE.0 ) THEN
        IF ( lwp ) WRITE(numout,*) 'Communications mapping failed : ',ierr
        RETURN
      ENDIF

      ! Prepare mpp north fold
#if defined key_mpp_mpi
      ! This invokes the version of the routine contained in this module
      ! and not the original in lib_mpp.F90
      CALL mpp_ini_north()
#endif

! From mppini_2.h90:
! Defined npolj, either 0, 3 , 4 , 5 , 6
! In this case the important thing is that npolj /= 0
! Because if we go through these line it is because jpni >1 and thus
! we must use lbcnorthmpp, which tests only npolj =0 or npolj /= 0
      npolj = 0
      IF( jperio == 3 .OR. jperio == 4 ) THEN

         IF( pjubext(narea) ) npolj = 3 ! Top row?
         IF( pjubext(narea) .AND. piubext(narea) ) npolj = 4 ! Top right? ARPDBG almost certainly wrong!
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN
         IF( pjubext(narea) ) npolj = 5 ! Top left?
         IF( pjubext(narea) .AND. piubext(narea) ) npolj = 6 ! Top right? ARPDBG almost certainly wrong!
      ENDIF

      ! Prepare NetCDF output file (if necessary)
      CALL mpp_init_ioipsl()

      ! Write this partition to file in the format that the code can
      ! read
      CALL write_partition()

      ! Initialise mbkmax because it's needed in any halo swap calls but we didn't have
      ! jpi and jpj until the partitioning was done.
      DO jj=1,jpj,1
         DO ji=1,jpi,1
            mbkmax(ji,jj) = ibotlevel(MIN(jpiglo,ji+nimpp-1), &
                                      MIN(jpjglo,jj+njmpp-1))
         END DO
      END DO

      ! ARPDBG - test comms setup
      CALL mpp_test_comms(imask, ibotlevel)

      ! Free array holding mask used for partitioning
      DEALLOCATE(imask)

    END SUBROUTINE mpp_init3

# if defined key_dimgout
   !!----------------------------------------------------------------------
   !!   'key_dimgout'                  NO use of NetCDF files
   !!----------------------------------------------------------------------
    SUBROUTINE mpp_init_ioipsl       ! Dummy routine
    END SUBROUTINE mpp_init_ioipsl
# else
    SUBROUTINE mpp_init_ioipsl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init_ioipsl  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! History :
      !!   9.0  !  04-03  (G. Madec)  MPP-IOIPSL 
      !!----------------------------------------------------------------------
      USE ioipsl
      IMPLICIT NONE
      !! Local declarations
      INTEGER, DIMENSION(2) :: iglo, iloc, iabsf, iabsl, ihals, ihale, idid
      !!----------------------------------------------------------------------

      ! The domain is decomposed in 2D only along i- or/and j- direction
      ! So we need at the most only 1D arrays with 2 elements
      iglo(1) = jpiglo
      iglo(2) = jpjglo
      iloc(1) = iesub
      iloc(2) = jesub
      iabsf(1) = ielb
      iabsf(2) = jelb
      iabsl(:) = iabsf(:) + iloc(:) - 1
      ihals(1) = jpreci
      ihals(2) = jprecj
      ihale(1) = jpreci
      ihale(2) = jprecj
      idid(1) = 1
      idid(2) = 2

      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'partmod: mpp_init_ioipsl :   iloc  = ', iloc (1), iloc (2)
          WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
          WRITE(numout,*) '                             ihals = ', ihals(1), ihals(2)
          WRITE(numout,*) '                             ihale = ', ihale(1), ihale(2)
      ENDIF

#if defined key_mpp_mpi
      CALL flio_dom_set ( mppsize, narea-1, idid, iglo, iloc, iabsf, iabsl, &
                          ihals, ihale, 'BOX', nidom)
#endif

    END SUBROUTINE mpp_init_ioipsl  
# endif

    SUBROUTINE partition_rk ( mask, istart, istop, jstart, jstop, ierr )
      !!------------------------------------------------------------------
      ! Partition the domain of nx x ny points
      ! according to the land/sea mask array
      ! using a recursive k-section algorithm,
      ! into nprocx x nprocy sub-domains.
      !
      !     mask                    real  input     Land/sea mask.
      !     gnx                     int   input     Size of the full domain.
      !     gny                     int   input
      !     ierr                    int   output    Error flag.
      !!------------------------------------------------------------------

      USE iom,     ONLY: jpiglo, jpjglo, wp
      USE par_oce, ONLY: jpni, jpnj
#if defined key_mpp_mpi
      USE lib_mpp, ONLY: mppsize
#endif
      IMPLICIT NONE

      ! Subroutine arguments.
      INTEGER, INTENT(out)       :: ierr
      INTEGER, INTENT(in)        :: istart, istop, jstart, jstop
      INTEGER, INTENT(in)        :: mask(:,:)
      ! Local variables
#if defined key_mpp_mpi
      INTEGER, DIMENSION(MAX_FACTORS) :: fx,fy
#endif
      INTEGER                    :: f,gnactive &
            ,i,ifax,ifin,ifx,ify,ilb,iproc,ist,isub,isub_old &
            ,isub_new,iub                                    &
            ,j,jfin,jlb,jst,jub,line                         &
            ,nfax,nfx,nfy,ngone,nsub_old,nsub_new,ntarget,ntry
      LOGICAL :: partx

      ! Clear the error flag.
      ierr = 0

#if defined key_mpp_mpi

#if defined PARTIT_DEBUG
      IF(lwp)WRITE(*,*) 'ARPDBG partition_rk: jpn{i,j} = ',jpni,jpnj
#endif
      ! Factorise the nprocx and nprocy parameters.
      CALL factor (fx,MAX_FACTORS,nfx,jpni,ierr)
      IF ( lwp .AND. ierr.NE.0 ) THEN
        WRITE (numout,*) 'partition_rk: factorisation in x failed ',ierr
        RETURN
      ENDIF
      CALL factor (fy,MAX_FACTORS,nfy,jpnj,ierr)
      IF ( lwp .AND. ierr.NE.0 ) THEN
        WRITE (numout,*) 'partition_rk: factorisation in y failed ',ierr
        RETURN
      ENDIF

#if defined PARTIT_DEBUG
      IF(lwp)THEN
         WRITE(*,*) 'ARPDBG partition_rk: nf{x,y} = ',nfx,nfy
         WRITE(*,*) 'ARPDBG partition_rk: fx = ',fx(:nfx)
         WRITE(*,*) 'ARPDBG partition_rk: fy = ',fx(:nfy)
      END IF
#endif

      CALL partition_rk_core(mask, jpiglo, jpjglo, &
                             MAX_FACTORS, fx, nfx, fy, nfy, ierr)

      CALL finish_partition()

#endif

    END SUBROUTINE partition_rk


    SUBROUTINE partition_mca_rk(mask, istart, istop, jstart, jstop, ierr)
#if defined key_mpp_mpi
       USE mpi
       USE lib_mpp, ONLY: mppsize, mpi_comm_opa, &
                          nxfactors, nyfactors, xfactors, yfactors
#endif
       USE lib_mpp, ONLY: ctl_stop
       USE dom_oce, ONLY: narea
       IMPLICIT NONE
       !!------------------------------------------------------------------
       !! Multi-Core Aware recursive partitioning of the domain. As for 
       !! partition_rk but choose the partion 
       !! so as to minimize off-node MPI communication
       !!
       !! Original by Stephen Pickles for POLCOMS code.
       !! Implementation in NEMO by Andrew Porter, 26/01/2012
       !!------------------------------------------------------------------
       ! Subroutine arguments.
       INTEGER, INTENT(in)        :: istart, istop, jstart, jstop
       INTEGER, INTENT(in)        :: mask(:,:)
       INTEGER, INTENT(out)       :: ierr
       ! Local variables
       INTEGER :: ii
#if defined key_mpp_mpi
       INTEGER, DIMENSION(MAX_FACTORS) :: fx, fy, factors
       INTEGER, DIMENSION(MAX_FACTORS) :: df, multiplicity
#endif
       INTEGER :: nfx, nfy, nfactors, ndf, nperms
       INTEGER :: check_nprocx, check_nprocy, check_nprocp
       INTEGER :: iperm
       CHARACTER(LEN=256) :: fstr
       INTEGER :: myinst                     ! MPI process ID of this process
       INTEGER :: nprocx, nprocy

       ! A high score is bad
       REAL(wp), DIMENSION(pv_num_scores)   :: score
       REAL(wp) :: best_score
       INTEGER  :: best_perm
       REAL(wp), DIMENSION(2,pv_num_scores) :: best, gbest, wrst, gwrst

#if defined key_mpp_mpi

       ! NEMO only has narea public and not the actual PE rank so
       ! set that here
       myinst = narea - 1

       ! Factorise the total number of MPI processes that we have
       CALL factor (factors,MAX_FACTORS,nfactors,nprocp,ierr)

       IF ( lwp ) THEN
          IF ( ierr.NE.0 ) THEN
             WRITE (numout,*) 'partition_mca_rk: factorisation failed ',ierr
             RETURN
          ELSE
             WRITE (numout,*) 'partition_mca_rk: factors are: ', factors(:nfactors)
          ENDIF
       ENDIF

       CALL calc_perms( nfactors, factors,     &
                        ndf, df, multiplicity, &
                        nperms )

       DO ii=1,pv_num_scores
          best(1,ii) = pv_awful
          best(2,ii) = -1.0_wp
       END DO
       DO ii=1,pv_num_scores
          wrst(1,ii) = 0.0_wp
          wrst(2,ii) = -1.0_wp
       END DO

       IF (lwp) THEN
          WRITE(numout,"('% Partn',2X,10(A4,2X),4(A9,1X),A7)")        &
                        'Wet', 'Dry',                      &
                        'pli', 'plx', 'pci', 'pcx',                 &
                        'nlx', 'ncx', 'tlx', 'tcx',                 &
                        'P comms', 'N comms', 'T comms', 'Overall', &
                        'Factors'
       END IF

       perm_loop: DO iperm=myinst, nperms-1, nprocp

          CALL get_perm_factors( iperm, nfactors, ndf, df, multiplicity, &
                                 fx, nfx, fy, nfy,                       &
                                 nprocx, nprocy, .FALSE. )

          CALL partition_rk_core(mask, jpiglo, jpjglo,    &
                                 MAX_FACTORS, fx, nfx, fy, nfy, ierr)

          IF (ierr .NE. 0) CYCLE perm_loop
          CALL finish_partition()

          ! Compute the cost function for this partition
          !
          CALL eval_partition(jpiglo, jpjglo, mask, score)
          CALL factor_string(fstr,nfx,fx,nfy,fy)

          WRITE (6,'(''%'',I6,1X,10(I5,1X),3(F9.2,1X),E12.4,1x,(A))') &
                    iperm,                                     &
                    INT(score(pv_index_wet)),                  &
                    INT(score(pv_index_dry)),                  &
                    INT(score(pv_index_pli)),                  &
                    INT(score(pv_index_plx)),                  &
                    INT(score(pv_index_pci)),                  &
                    INT(score(pv_index_pcx)),                  &
                    INT(score(pv_index_nlx)),                  &
                    INT(score(pv_index_ncx)),                  &
                    INT(score(pv_index_tlx)),                  &
                    INT(score(pv_index_tcx)),                  &
                    score(pv_index_pcomms),                    &
                    score(pv_index_ncomms),                    &
                    score(pv_index_tcomms),                    &
                    score(pv_index_overall),                   &
                    TRIM(fstr)

          DO ii=1,pv_num_scores
             IF (score(ii) .LT.  best(1,ii)) THEN
                best(1,ii) = score(ii)
                best(2,ii) = iperm
             END IF
             IF (score(ii) .GT. wrst(1,ii)) THEN
                wrst(1,ii) = score(ii)
                wrst(2,ii) = iperm
             END IF
          END DO
         
      END DO perm_loop

      !  Now choose the "best" partition

#if defined key_mpp_mpi
      CALL MPI_ALLREDUCE(best, gbest, pv_num_scores, &
                         MPI_2DOUBLE_PRECISION,      &
                         MPI_MINLOC, mpi_comm_opa, ierr)
      CALL MPI_ALLREDUCE(wrst, gwrst, pv_num_scores, &
                         MPI_2DOUBLE_PRECISION,      &
                         MPI_MAXLOC, mpi_comm_opa, ierr)
#else
      CALL ctl_stop('STOP', 'partition_mca_rk: this version requires MPI')
#endif
      best_score = gbest(1,pv_index_overall)
      best_perm  = gbest(2,pv_index_overall)

      IF ( lwp ) THEN
         WRITE (numout,'(A32,A20,A20)')  &
                ' ','  best score / perm ','  worst score / perm'
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)') 'overall: ',    &
                gbest(1,pv_index_overall), gbest(2,pv_index_overall), &
                gwrst(1,pv_index_overall), gwrst(2,pv_index_overall)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)') 'wet points: ', &
                gbest(1,pv_index_wet), gbest(2,pv_index_wet),         &
                gwrst(1,pv_index_wet), gwrst(2,pv_index_wet)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)') 'dry points: ', &
                gbest(1,pv_index_dry), gbest(2,pv_index_dry),         &
                gwrst(1,pv_index_dry), gwrst(2,pv_index_dry)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'proc max  on-node wet points: ',                     &
                gbest(1,pv_index_pli), gbest(2,pv_index_pli),         &
                gwrst(1,pv_index_pli), gwrst(2,pv_index_pli)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'proc max off-node wet points: ',                     &
                gbest(1,pv_index_plx), gbest(2,pv_index_plx),         &
                gwrst(1,pv_index_plx), gwrst(2,pv_index_plx)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'proc max  on-node   messages: ',                     &
                gbest(1,pv_index_pci), gbest(2,pv_index_pci),         &
                gwrst(1,pv_index_pci), gwrst(2,pv_index_pci)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'proc max off-node   messages: ',                     &
                gbest(1,pv_index_pcx), gbest(2,pv_index_pcx),         &
                gwrst(1,pv_index_pcx), gwrst(2,pv_index_pcx)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'node max off-node wet points: ',                     &
                gbest(1,pv_index_nlx), gbest(2,pv_index_nlx),         &
                gwrst(1,pv_index_nlx), gwrst(2,pv_index_nlx)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'node max off-node   messages: ',                     &
                gbest(1,pv_index_ncx), gbest(2,pv_index_ncx),         &
                gwrst(1,pv_index_ncx), gwrst(2,pv_index_ncx)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'total    off-node wet points: ',                     &
                gbest(1,pv_index_tlx), gbest(2,pv_index_tlx),         &
                gwrst(1,pv_index_tlx), gwrst(2,pv_index_tlx)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'per core communications cost: ',                     &
                gbest(1,pv_index_pcomms), gbest(2,pv_index_pcomms),   &
                gwrst(1,pv_index_pcomms), gwrst(2,pv_index_pcomms)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'per node communications cost: ',                     &
                gbest(1,pv_index_ncomms), gbest(2,pv_index_ncomms),   &
                gwrst(1,pv_index_ncomms), gwrst(2,pv_index_ncomms)
         WRITE (numout,'(A32,F12.3,F8.0,E12.3,F8.0)')                 &
                'network communications cost: ',                      &
                gbest(1,pv_index_tcomms), gbest(2,pv_index_tcomms),   &
                gwrst(1,pv_index_tcomms), gwrst(2,pv_index_tcomms)

         WRITE (numout,"('partition_mca_rk: overall best perm is ',I6,', score=',F12.3)") &
                best_perm, best_score
      END IF

      ! Use the same partition on all processes

      ! If a particular factorisation has been forced by
      ! the nn_{x,y}factors fields in the nammpp section of the namelist
      ! then use that one instead

      IF ((nxfactors + nyfactors) > 0) THEN
         check_nprocx = 1
         check_nprocy = 1
         DO ii=1,nxfactors
            check_nprocx = check_nprocx*xfactors(ii)
         END DO
         DO ii=1,nyfactors
            check_nprocy = check_nprocy*yfactors(ii)
         END DO
         check_nprocp = check_nprocx*check_nprocy
         IF (check_nprocp .EQ. nprocp) THEN
            nprocx = check_nprocx
            nprocy = check_nprocy
            nfx = nxfactors
            nfy = nyfactors
            fx(1:nfx) = xfactors(1:nfx)
            fy(1:nfy) = yfactors(1:nfy)
         ELSE
            CALL ctl_stop('STOP', 'part_mca_rk: '//   &
                          'requested factorisation does not match nprocp')
         END IF
      ELSE
         ! Use the best factorisation found above
         IF (best_perm < 0.0_wp) THEN
            IF (lwp) THEN
               WRITE (numout,*) 'partition_mca_rk: no feasible partition found'
            END IF
            ierr = 1
            RETURN
         END IF
         CALL get_perm_factors(best_perm, nfactors, ndf, df, multiplicity, &
                               fx, nfx, fy, nfy,                           &
                               nprocx, nprocy, lwp )
      END IF

      ! Set corresponding NEMO variables for PE grid, even though it is now
      ! rather irregular
      jpni = nprocx
      jpnj = nprocy
      jpnij = jpni*jpnj

      IF (lwp) THEN
         WRITE (numout,'(A39)',advance='no') &
                'partition_mca_rk: partitioning with factors '
         CALL print_factors(numout, nfx, fx, nfy, fy)
      END IF


      CALL partition_rk_core ( mask, jpiglo, jpjglo, &
                               MAX_FACTORS,          &
                               fx, nfx, fy, nfy,     &
                               ierr )

      IF (ierr .NE. 0) THEN
         IF (lwp) THEN
            WRITE (numout,*) 'partition_mca_rk: failed to apply selection partition'
         END IF
         RETURN
      END IF
      CALL finish_partition()
 
      IF (lwp) THEN
         CALL eval_partition(jpiglo, jpjglo, mask, score)
         CALL factor_string(fstr,nfx,fx,nfy,fy)
         WRITE(numout,'(10(A7,1X),4(A9,1X),A7)')           &
               'Wet', 'Dry',                               &
               'pli', 'plx', 'pci', 'pcx',                 &
               'nlx', 'ncx', 'tlx', 'tcx',                 &
               'P comms', 'N comms', 'T comms', 'Overall', &
               'Factors'
        WRITE (6,'(10(F7.0,1X),4(F9.2,1X),A50)') &
               score(pv_index_wet),              &
               score(pv_index_dry),              &
               score(pv_index_pli),              &
               score(pv_index_plx),              &
               score(pv_index_pci),              &
               score(pv_index_pcx),              &
               score(pv_index_nlx),              &
               score(pv_index_ncx),              &
               score(pv_index_tlx),              &
               score(pv_index_tcx),              &
               score(pv_index_pcomms),           &
               score(pv_index_ncomms),           &
               score(pv_index_tcomms),           &
               score(pv_index_overall),          &
               fstr
     END IF

#endif

   END SUBROUTINE partition_mca_rk


   SUBROUTINE partition_rk_core( mask, nx, ny, maxfax, fx, nfx, fy, nfy,   &
                                 ierr )
#if defined key_mpp_mpi
       USE lib_mpp, ONLY: mppsize
#endif
       IMPLICIT NONE
       !!------------------------------------------------------------------
       !!------------------------------------------------------------------
       INTEGER, INTENT(in)  :: nx, ny
       INTEGER, INTENT(in)  :: mask(nx,ny)
       INTEGER, INTENT(in)  :: maxfax, nfx, nfy
       INTEGER, INTENT(in)  :: fx(maxfax), fy(maxfax)
       INTEGER, INTENT(out) :: ierr
       ! Local variables
       INTEGER :: istart, jstart, istop, jstop
       INTEGER :: f,gnactive 
       INTEGER :: i,ifax,nfax,ifin,ifx,ify,ilb,iproc
       INTEGER :: ist,isub,isub_old,isub_new,iub    
       INTEGER :: j,jfin,jlb,jst,jub,line
       INTEGER :: ngone,nsub_old,nsub_new,ntarget,ntry
       LOGICAL :: partx

       ! Zero the error flag
       ierr = 0

       ! Count the significant (non-zero) factors.
       nfax = nfx+nfy
#if defined PARTIT_DEBUG
       IF ( lwp ) THEN
          WRITE (numout,"('jpni = ',I3,1x,'nfx = ',I3,/'Factorized into: ',(I3,1x))") &
                 jpni, nfx, fx(:nfx)
          WRITE (numout,"('jpnj = ',I3,1x,'nfy = ',I3,/'Factorized into: ',(I3,1x))") & 
                 jpnj, nfy, fy(:nfy)
          WRITE (numout,"('Partitioning a total of ',I3,' times')") nfax
       ENDIF
#endif

       ! Define the full domain as the one and only sub-domain.
       istart = 1
       istop = nx
       jstart = 1
       jstop = ny

       nsub_old = 1
       pielb(nsub_old) = istart
       pjelb(nsub_old) = jstart
       pieub(nsub_old) = istop
       pjeub(nsub_old) = jstop
       ! Count the number of active points.
       gnactive = 0
       DO j=jstart,jstop,1
          DO i=istart,istop,1
             IF ( mask(i,j) == 1 ) gnactive = gnactive+1
          ENDDO
       ENDDO
#if defined PARTIT_DEBUG
       IF ( lwp )WRITE (numout,*) 'Total wet points ',gnactive
#endif
       pnactive(nsub_old) = gnactive

       ! Partition for each factor.
       DO ifax=1,nfax
          IF ( ifax.EQ.1 ) THEN
             ! Start by partitioning the dimension with more factors.
             partx = nfx.GE.nfy
             ifx = 0
             ify = 0
          ELSE
             ! Try to toggle the partitioning direction.
             partx = .NOT.partx
             IF ( partx ) THEN
                ! If we have run out of factors in x, switch to y.
                partx = .NOT. ifx+1.GT.nfx
             ELSE
                ! If we have run out of factors in y, switch to x.
                partx =       ify+1.GT.nfy
             ENDIF
          ENDIF
#if defined PARTIT_DEBUG
          IF ( lwp ) THEN
             WRITE (numout,*)
             WRITE (numout,*) '###########################################'
             WRITE (numout,*)
             WRITE (numout,*) 'Partition ',ifax,'partx = ',partx
          ENDIF
#endif
          IF ( partx ) THEN
             ! ============================================================
             !         Partition in x.
             ! ============================================================
             ifx = ifx+1
             f = fx(ifx)
             nsub_new = nsub_old*f
#if defined PARTIT_DEBUG
             IF ( lwp ) THEN
                WRITE (numout,*) 'Partitioning in x from ',nsub_old &
                                ,' to ',nsub_new
             ENDIF
#endif
             DO isub_old=1,nprocp,nprocp/nsub_old
                ! Check that there are sufficient points to factorise.
                IF ( pieub(isub_old)-pielb(isub_old)+1.LT.f ) THEN
                   WRITE (numout,*)  &
                        'partition_rk: insufficient points to partition'
                   ierr = 1
                   RETURN
                ENDIF

                ! Set the target number of points in the new sub-domains.
                ntarget = NINT(REAL(pnactive(isub_old))/REAL(f))
                ist = pielb(isub_old)
                iub = pieub(isub_old)
                jlb = pjelb(isub_old)
                jub = pjeub(isub_old)
#if defined PARTIT_DEBUG
                IF ( lwp ) THEN
!                   WRITE (numout,*)
                   WRITE (numout,"(/'=======================================')")
!                   WRITE (numout,*)
                   WRITE (numout,"(/'Partitioning sub-domain ',I3,': (',I3,':',I3,')(',I3,':',I3,')')") &
                                    isub_old,ist,iub,jlb,jub
                   WRITE (numout,"('Old partition has ',I8,' points')") pnactive(isub_old)
                   WRITE (numout,"('Target is ',I8,' points')") ntarget
                ENDIF
#endif
                ! Create the new sub-domains.
                ngone = 0
                DO isub=1,f,1
                   isub_new = isub_old + (isub-1)*nprocp/nsub_new
#if defined PARTIT_DEBUG
                   IF ( lwp ) THEN
                      WRITE (numout,*)
                      WRITE (numout,*) 'Creating new sub-domain ',isub_new
                   ENDIF
#endif
                   ! The last new sub-domain takes the remaining data.
                   IF ( isub.EQ.f ) THEN
                      ifin = iub
                   ELSE
                      ! Scan lines in x counting active grid points until we
                      ! exceed the target.
                      ntry = 0
                      ifin = -1
                      scanrows: DO i=ist,iub,1
                         ! Count up the contribution from the next line.
                         line = 0
                         DO j=jlb,jub,1
                            IF ( mask(i,j) == 1 ) line = line+1
                         ENDDO

                         ! Check against the target.
                         IF ( ntry+line.GT.ntarget ) THEN
                            ! Is it nearer the target to include this line or not ?
                            IF ( ntry+line-ntarget.GT.ntarget-ntry ) THEN
                               ! Nearer start than end. Finish at previous line.
                               ifin = i-1
                            ELSE
                               ! Nearer end than start. Include this line.
                               ifin = i
                               ntry = ntry + line
                            ENDIF
#if defined PARTIT_DEBUG
                            IF ( lwp ) THEN
                               WRITE (numout,*) 'Reached target at ',ifin  &
                                               ,' ntry = ',ntry
                            ENDIF
#endif
                            EXIT scanrows
                         ENDIF
                         ! Add in the line count to the running total.
                         ntry = ntry + line
                      ENDDO scanrows
                      IF ( ifin.LT.0 ) ifin = iub
                   ENDIF

                   ! Set up the parameters of the new sub-domain.
                   pnactive(isub_new) = 0
                   DO j=jlb,jub
                      DO i=ist,ifin
                         IF ( mask(i,j) == 1 ) THEN
                            pnactive(isub_new) = pnactive(isub_new)+1
                         ENDIF
                      ENDDO
                   ENDDO
                   pielb(isub_new) = ist
                   pieub(isub_new) = ifin
                   pjelb(isub_new) = jlb
                   pjeub(isub_new) = jub
#if defined PARTIT_DEBUG
                   IF ( lwp ) THEN
                      WRITE (numout,*) 'New subdomain is ',ist,ifin,jlb,jub &
                           ,' with ',pnactive(isub_new),' points'
                   ENDIF
#endif
!             Reset the target based on the actual number of points
!             remaining, split between the partitions remaining.
                   ngone = ngone+ntry
#if defined PARTIT_DEBUG
                   IF ( lwp ) THEN
                      !               write (numout,*) 'Target reset to ',ntarget
                   ENDIF
#endif
!             Start the next search at the next point.
                   ist = ifin+1
                ENDDO
             ENDDO

          ELSE

!         ============================================================
!         Partition in y.
!         ============================================================
             ify = ify+1
             f = fy(ify)
             nsub_new = nsub_old*f
#if defined PARTIT_DEBUG
             IF ( lwp ) THEN
                WRITE (numout,*) 'Partitioning in y from ',nsub_old &
                                ,' to ',nsub_new
             ENDIF
#endif
             DO isub_old=1,nprocp,nprocp/nsub_old
!           Check that there are sufficient points to factorise.
                IF ( pjeub(isub_old)-pjelb(isub_old)+1.LT.f ) THEN
                   WRITE (numout,*)  &
                        'partition_rk: insufficient points to partition'
                   ierr = 1
                   RETURN
                ENDIF
!           Set the target number of points in the new sub-domains.
                ntarget = NINT(REAL(pnactive(isub_old))/REAL(f))
                ilb = pielb(isub_old)
                iub = pieub(isub_old)
                jst = pjelb(isub_old)
                jub = pjeub(isub_old)
#if defined PARTIT_DEBUG
                IF ( lwp ) THEN
                   WRITE (numout,*)
                   WRITE (numout,*) '======================================='
                   WRITE (numout,*)
                   WRITE (numout,*) 'Partitioning sub-domain ',isub_old,' : ' &
                       ,ilb,iub,jst,jub
                   WRITE (numout,*) 'Old partition has ',pnactive(isub_old)   &
                       ,' points'
                   WRITE (numout,*) 'Target is ',ntarget
                ENDIF
#endif
!           Create the new sub-domains.
                ngone = 0
                DO isub=1,f
                   isub_new = isub_old + (isub-1)*nprocp/nsub_new
#if defined PARTIT_DEBUG
                   IF ( lwp ) THEN
                      WRITE (numout,*)
                      WRITE (numout,*) 'Creating new sub-domain ',isub_new
                   ENDIF
#endif
!             The last new sub-domain takes the remaining data.
                   IF ( isub.EQ.f ) THEN
                      jfin = jub
                   ELSE
!               Scan lines in y counting active grid points until we
!               exceed the target.
                      ntry = 0
                      jfin = -1
                      scancols: DO j=jst,jub
!                 Count up the contribution from the next line.
                         line = 0
                         DO i=ilb,iub
                            IF ( mask(i,j) == 1 ) line = line+1
                         ENDDO
#if defined PARTIT_DEBUG
                         IF ( lwp ) THEN
                            !dbg    write (numout,*) 'Line ',j,' has ',line
                         ENDIF
#endif
                         ! Check against the target.
                         IF ( ntry+line.GT.ntarget ) THEN
                            ! Is it nearer the target to include this line or not ?
                            IF ( ntry+line-ntarget.GT.ntarget-ntry ) THEN
                               ! Nearer start than end. Finish at previous line.
                               jfin = j-1
                            ELSE
                               ! Nearer end than start. Include this line.
                               jfin = j
                               ntry = ntry + line
                            ENDIF
#if defined PARTIT_DEBUG
                            IF ( lwp ) THEN
                               WRITE (numout,*) 'Reached target at ',jfin &
                                    ,' ntry = ',ntry
                            ENDIF
#endif
                            EXIT scancols
                         ENDIF
                         ! Add in the line count to the running total.
                         ntry = ntry + line
                      ENDDO scancols
                      IF ( jfin.LT.0 ) jfin = jub
                   ENDIF
                   ! Set up the parameters of the new sub-domain.
                   pnactive(isub_new) = 0
                   DO j=jst,jfin
                      DO i=ilb,iub
                         IF ( mask(i,j) == 1 ) THEN
                            pnactive(isub_new) = pnactive(isub_new)+1
                         ENDIF
                      ENDDO
                   ENDDO
                   pielb(isub_new) = ilb
                   pieub(isub_new) = iub
                   pjelb(isub_new) = jst
                   pjeub(isub_new) = jfin
#if defined PARTIT_DEBUG
                   IF ( lwp ) THEN
                      WRITE (numout,*) 'New subdomain is ',ilb,iub,jst,jfin &
                         ,' with ',pnactive(isub_new),' points'
                   ENDIF
#endif
!             Reset the target based on the actual number of points
!             remaining, split between the partitions remaining.
                   ngone = ngone+ntry
#if defined PARTIT_DEBUG
                   IF(lwp)WRITE (numout,*) 'Target reset to ',ntarget
#endif
!             Start the next search at the next point.
                   jst = jfin+1
                ENDDO
             ENDDO
          ENDIF
!       Update the number of sub-domains ready for the next loop.
          nsub_old = nsub_new
       ENDDO
#if defined PARTIT_DEBUG
       IF ( lwp ) THEN
          WRITE (numout,*)
          WRITE (numout,*) '========================================='
          WRITE (numout,*)
       ENDIF
#endif
!     Set the size of each sub-domain.
       DO iproc=1,nprocp
          piesub(iproc) = pieub(iproc)-pielb(iproc)+1
          pjesub(iproc) = pjeub(iproc)-pjelb(iproc)+1
#if defined PARTIT_DEBUG
          IF(lwp)THEN
             WRITE (numout,*) 'Domain ',iproc,' has ',pnactive(iproc), ' ocean points'
             WRITE(*,"(I3,': [',I3,':',I3,'][',I3,':',I3,']')") &
                    iproc, pielb(iproc), pieub(iproc), pjelb(iproc), pjeub(iproc)

          END IF
#endif
       ENDDO

     END SUBROUTINE partition_rk_core


     SUBROUTINE factor ( ifax, maxfax, nfax, n, ierr )

        !!------------------------------------------------------------------
        ! Subroutine to return the prime factors of n.
        ! nfax factors are returned in array ifax which is of maximum
        ! dimension maxfax.
        !!------------------------------------------------------------------
        IMPLICIT NONE
        ! Subroutine arguments
        INTEGER ierr, n, nfax, maxfax
        INTEGER ifax(maxfax)
        ! Local variables.
        INTEGER i, ifac, l, nu
        INTEGER lfax(20)
        ! lfax contains the set of allowed factors.
        DATA (lfax(i),i=1,20) / 71,67,59,53,47,43,41,37,31,29 &
                               ,23,19,17,13,11, 7, 5, 3, 2, 1 /
        ! Clear the error flag.
        ierr = 0
        ! Find the factors of n.
        IF ( n.EQ.1 ) THEN
           nfax = 0
           GOTO 20
        ENDIF
        ! nu holds the unfactorised part of the number.
        ! nfax holds the number of factors found.
        ! l points to the allowed factor list.
        ! ifac holds the current factor.
        nu = n
        nfax = 0
        l = 1
        ifac = lfax(l)
        ! Label 10 is the start of the factor search loop.
10      CONTINUE
        ! Test whether the factor will divide.
        IF ( MOD(nu,ifac).EQ.0 ) THEN
           ! Add the factor to the list.
           nfax = nfax+1
           IF ( nfax.GT.maxfax ) THEN
              ierr = 6
              WRITE (*,*)  &
                   'FACTOR: insufficient space in factor array ',nfax
              RETURN
           ENDIF
           ifax(nfax) = ifac
           ! Divide out the factor from the remaining number.
           ! If unity remains, the number has been
           ! successfully factorized.
           nu = nu/ifac
           IF ( nu.EQ.1   ) GO TO 20
           ! Loop to try the factor again.
           GOTO 10
        ENDIF
        ! Move on to the next factor in the list.
        l = l+1
        ifac = lfax(l)
        ! Unless the end of the factor list has been reached, loop.
        IF ( ifac.GT.1 ) go to 10
        ! There is a factor in n which is not in the lfac list.
        ! Add the remaining number onto the end of the list.
        nfax = nfax+1
        IF ( nfax.GT.maxfax ) THEN
           ierr = 6
           WRITE (*,*) 'FACTOR: insufficient space in factor array ',nfax
           RETURN
        ENDIF
        ifax(nfax) = nu
        ! Label 20 is the exit point from the factor search loop.
20      CONTINUE

      END SUBROUTINE factor

!#define TRIM_DEBUG

      SUBROUTINE part_trim ( depth, isTrimmed, ierr )
        !!------------------------------------------------------------------
        !!
        !! Examines all the sub-domains and trims off boundary rows and
        !! columns which are all land.
        !!
        !!     depth                   real  input     Depth map.
        !!     isTrimmed               logical output  Whether N,E,S,W edge
        !!                                             of domain has been 
        !!                                             trimmed
        !!     ierr                    int   output    Error flag.
        !!
        !!        Mike Ashworth, CLRC Daresbury Laboratory, July 1999
        !!------------------------------------------------------------------
        USE par_oce, ONLY: jpreci, jprecj
        USE iom,     ONLY: jpiglo, jpjglo
        IMPLICIT NONE

        ! Subroutine arguments.
        INTEGER, INTENT(in)  :: depth(jpiglo,jpjglo)
        LOGICAL, DIMENSION(:,:), INTENT(out) :: isTrimmed
        INTEGER, INTENT(out) :: ierr
        ! Local variables.       
        INTEGER :: i, iproc, j, newbound

        ! Clear the error code.
        ierr = 0

        ! Clear all flags
        ! N E S W
        ! 1 2 3 4
        isTrimmed(1:4,1:nprocp) = .FALSE.

        ! Examine each sub-domain in turn.
        
        subdomain_loop: DO iproc=1,nprocp

           ! Do not trim if there are no active points at all.
           ! Otherwise we will trim away the whole sub-domain and we
           ! will be in big trouble.

           ! Look at the low-i columns (Western edge of domain)

           newbound = pielb(iproc)
           lowi: DO i=pielb(iproc),pieub(iproc)

              ! Scan this column for wet points

              DO j=MAX(1,pjelb(iproc)-jprecj),MIN(jpjglo,pjeub(iproc)+jprecj)

                 IF ( depth(i,j) == 1 ) THEN
                     newbound = MAX(i - jpreci - nextra, pielb(iproc))
#if defined TRIM_DEBUG
                 IF ( lwp ) THEN
                    WRITE (numout,*) 'Sub-domain ',iproc,': Low-i loop: row ',i &
                                     ,' is land. New bound is ',newbound
                 ENDIF
#endif
                     EXIT lowi
                 ENDIF
               ENDDO
            ENDDO lowi

            ! Update with the new boundary and monitor.

            IF ( newbound.NE.pielb(iproc) ) THEN
#if defined TRIM_DEBUG
               IF ( lwp ) THEN
                  IF ( newbound-pielb(iproc).GT.1 ) THEN
                     WRITE(numout,'(a,i5,3(a,i6))') ' Process ',iproc-1 &
                                    ,' trimmed ',newbound-pielb(iproc) &
                          ,' cols: ',pielb(iproc),' to ',newbound-1
                  ELSE
                     WRITE(numout,'(a,i5,3(a,i6))') ' Process ',iproc-1 &
                                    ,' trimmed ',newbound-pielb(iproc) &
                                    ,' col : ',pielb(iproc)
                  ENDIF
               ENDIF
#endif
               pielb(iproc) = newbound
               isTrimmed(4,iproc) = .TRUE.
            ENDIF

            ! Look at the high-i columns (Eastern edge of domain).

            newbound = pieub(iproc)
            highi: DO i=pieub(iproc),pielb(iproc),-1

               DO j=MAX(1,pjelb(iproc)-jprecj), MIN(jpjglo,pjeub(iproc)+jprecj)

                  IF ( depth(i,j) == 1 ) THEN
                     ! We've found a wet point in this column so this is as far 
                     ! as we can trim.
                     newbound = MIN(i + jpreci + nextra, pieub(iproc))
#if defined TRIM_DEBUG
                     IF ( lwp ) THEN
                        WRITE (numout,"('Sub-domain ',I3,': High-i loop: col ',I3, &
                          &   ' contains water. New bound is ',I3)") iproc,i,newbound
                     ENDIF
#endif
                     EXIT highi

                  ENDIF
               ENDDO
            ENDDO highi

            ! Update with the new boundary and monitor.

          IF ( newbound.NE.pieub(iproc) ) THEN
#if defined TRIM_DEBUG
             IF ( lwp ) THEN
                IF ( (pieub(iproc)-newbound) .GT.1 ) THEN
                   WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                    ,' trimmed ',pieub(iproc)-newbound &
                          ,' cols: ',newbound+1,' to ',pieub(iproc)
                ELSE
                   WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                    ,' trimmed ',pieub(iproc)-newbound &
                                    ,' col : ',pieub(iproc)
              ENDIF
            ENDIF
#endif
            pieub(iproc) = newbound
            isTrimmed(2,iproc) = .TRUE.
          ENDIF

          ! Look at the low-j rows (Southern edge of domain)

          newbound = pjelb(iproc)
          lowj: DO j=pjelb(iproc),pjeub(iproc),1
             
             ! Scan this row for wet points

             DO i=MAX(1,pielb(iproc)-jpreci),MIN(jpiglo,pieub(iproc)+jpreci)
                IF ( depth(i,j) == 1) THEN
                   newbound = MAX(j - jpreci - nextra, pjelb(iproc))
#if defined TRIM_DEBUG
                   IF ( lwp ) THEN
                      WRITE (numout,*) 'Sub-domain ',iproc,': Low-j loop: column ',j &
                            ,' is land. New bound is ',newbound
                   ENDIF
#endif
                   EXIT lowj
                ENDIF
             ENDDO

          ENDDO lowj


          ! Update with the new boundary and monitor.

          IF ( newbound.NE.pjelb(iproc) ) THEN
#if defined TRIM_DEBUG
             IF ( lwp ) THEN
                IF ( (newbound-pjelb(iproc)) .GT.1 ) THEN
                   WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                     ,' trimmed ',newbound-pjelb(iproc) &
                                     ,' rows: ',pjelb(iproc),' to ',newbound-1
                ELSE
                   WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                     ,' trimmed ',newbound-pjelb(iproc) &
                                     ,' row : ',pjelb(iproc)
                ENDIF
             ENDIF
#endif
             pjelb(iproc) = newbound
             isTrimmed(3,iproc) = .TRUE.
         ENDIF

         ! Look at the high-j rows (Northern edge of domain)

         newbound = pjeub(iproc)
         highj: DO j=pjeub(iproc),pjelb(iproc),-1

            ! Scan this row for wet points

            DO i=MAX(1,pielb(iproc)-jpreci),MIN(jpiglo,pieub(iproc)+jpreci)
               IF ( depth(i,j) == 1 ) THEN
                  newbound = MIN(j + jpreci + nextra, pjeub(iproc))
#if defined TRIM_DEBUG
                  IF ( lwp ) then
                     WRITE (numout,*) 'Sub-domain ',iproc,': High-j loop: column ',j &
                                             ,' is land. New bound is ',newbound
                  ENDIF
#endif
                  EXIT highj
               ENDIF
            ENDDO
         ENDDO highj

          ! Update with the new boundary and monitor.

          IF ( newbound.NE.pjeub(iproc) ) THEN
#if defined TRIM_DEBUG
            IF ( lwp ) THEN
              IF ( pjeub(iproc)-newbound.GT.1 ) THEN
                WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                  ,' trimmed ',pjeub(iproc)-newbound &
                                  ,' rows: ',newbound+1,' to ',pjeub(iproc)
              ELSE
                WRITE (numout,'(a,i5,3(a,i6))') ' Sub-domain ',iproc &
                                  ,' trimmed ',pjeub(iproc)-newbound &
                                  ,' row : ',pjeub(iproc)
              ENDIF
            ENDIF
#endif
             pjeub(iproc) = newbound
             isTrimmed(1,iproc) = .TRUE.
          ENDIF

          ! Reset the size of the sub-domain.

          piesub(iproc) = pieub(iproc)-pielb(iproc)+1
          pjesub(iproc) = pjeub(iproc)-pjelb(iproc)+1

          ! endif active_subdomain

       ENDDO subdomain_loop

   END SUBROUTINE part_trim


   SUBROUTINE finish_partition(fromFile)
      USE mapcomm_mod, ONLY: ielb,ieub,pielb,pjelb,pieub,pjeub,          &
                             iesub,jesub,jeub,ilbext,iubext,jubext,      &
                             jlbext,pnactive,piesub,pjesub,jelb,pilbext, &
                             piubext, pjlbext, pjubext,                  &
                             trimmed, nidx,eidx,sidx,widx
      IMPLICIT NONE
      LOGICAL, INTENT(in), OPTIONAL :: fromFile
      ! Locals
      INTEGER :: iproc, ierr
      LOGICAL :: lFromFile

      ! Check to see whether we're dealing with a partion that has been
      ! read from file. If we are then there are some things we don't
      ! calculate here.
      lFromFile = .FALSE.
      IF( PRESENT(fromFile) ) lFromFile = fromFile

      IF(.NOT. lFromFile)THEN
         ! Set the external boundary flags before boundaries are
         ! altered by the trimming process and it becomes more difficult
         ! to recognize which were the external boundaries.
      
         DO iproc=1, nprocp, 1
            pilbext(iproc) = pielb(iproc).EQ.1
            piubext(iproc) = pieub(iproc).EQ.jpiglo
            pjlbext(iproc) = pjelb(iproc).EQ.1
            pjubext(iproc) = pjeub(iproc).EQ.jpjglo
         ENDDO

         ! Trim off redundant rows and columns containing all land.
         IF(.NOT. ALLOCATED(trimmed) )THEN
            ALLOCATE(trimmed(4,nprocp), Stat=ierr)
            IF(ierr /= 0)THEN
               CALL ctl_stop('STOP',    &
                             'Failed to allocate memory for domain trimming')
            END IF
         END IF

#if defined key_mpp_mpi
        IF ( nn_pttrim ) THEN
           nextra = 2
           CALL part_trim ( imask, trimmed, ierr )
        ELSE
           ! Need non-zero nextra because otherwise hit trouble with fields
           ! not being read from disk over land regions
           nextra = 2
           !nextra = 0 ! Don't need to back-off on message trimming
                      ! if we're not trimming the domains
           trimmed(1:4,1:nprocp) = .FALSE.
        ENDIF
#else
         trimmed(1:4,1:nprocp) = .FALSE.
#endif
      END IF ! not read from file

      ! Lower boundary (long.) of sub-domain, GLOBAL coords
      ! before correction for global halos
      nimpp = pielb(narea)

      ! Is the domain on an external LONGITUDE boundary?
      nbondi = 0
      ilbext = pilbext(narea)
      IF(ilbext)THEN
         nbondi = -1
      END IF

      IF( (.NOT. ilbext) .OR. (ilbext .AND. trimmed(widx,narea)) )THEN
         ! It isn't, which means we must shift its lower boundary by
         ! -jpreci to allow for the overlap of this domain with its
         ! westerly neighbour.
         nimpp = nimpp - jpreci
      END IF

      iubext = piubext(narea)
      IF(iubext)THEN
         nbondi = 1
         IF(ilbext)nbondi = 2 ! Both East and West boundaries are at 
                              ! edges of global domain
      END IF

      ! Set local values for limits in global coords of the sub-domain 
      ! owned by this process.
      ielb   = pielb (narea)
      ieub   = pieub (narea)
      iesub  = piesub(narea)

      jpi  = iesub + 2*jpreci ! jpi is the same for all domains - this is
                              ! what original decomposition did
      nlci = jpi

      ! If the domain is at the edge of the model domain and a cyclic 
      ! East-West b.c. is in effect then it already incorporates one 
      ! extra halo column (because of the way the model domain itself is
      ! set-up). If we've trimmed-off dry rows and columns then, even if
      ! a domain is on the model boundary, it may still need a halo so
      ! we add one.
      IF( nbondi == -1 .AND. (.NOT. trimmed(widx,narea))  )THEN
           ! Western boundary
           ! First column of global domain is actually a halo but NEMO
           ! still sets nldi to 1.
           nldi = 1            ! Lower bnd of int. region of sub-domain, LOCAL
           nlei = nldi + iesub - 1
           nlci = nlei + jpreci
           jpi  = nlci

        ELSE IF( nbondi == 1 .AND. (.NOT. trimmed(eidx,narea)) )THEN
           ! Eastern boundary
           nldi = 1 + jpreci
           ! Last column of global domain is actually a halo
           nlei = nldi + iesub - 1
           nlci = nlei
           jpi  = nlei

        ELSE IF( nbondi == 2)THEN
           ! Both boundaries are on the edges of the global domain
           IF(trimmed(widx,narea))THEN
              nldi = 1 + jpreci
           ELSE
              nldi = 1
           END IF
           nlei = nldi + iesub - 1

           IF(trimmed(eidx,narea))THEN
              nlci = nlei + jpreci
           ELSE
              nlci = nlei
           END IF
           jpi  = nlci

        ELSE
           ! We have no external boundaries to worry about
           nldi = 1 + jpreci 
           nlei = nldi + iesub - 1 !
        END IF

        jpim1 = jpi - 1

        ! Lower ext. boundary (lat.) of sub-domain, global coords
        njmpp= pjelb (narea)

        ! Is the domain on an external LATITUDE boundary?
        nbondj = 0
        jlbext = pjlbext(narea)
        IF(jlbext)THEN
           nbondj = -1
        ENDIF

        IF((.NOT. jlbext) .OR. (jlbext .AND. trimmed(sidx,narea)) )THEN
           ! It isn't, which means we must shift its lower boundary by
           ! -jprecj to allow for the overlap of this domain with its
           ! southerly neighbour.
           njmpp = njmpp - jprecj
        END IF
        ! ARPDBG - should we allow for trimming of northern edge of
        ! sub-domains here?
        jubext = pjubext(narea)
        IF(jubext)THEN
           nbondj = 1
           IF(jlbext)nbondj = 2
        END IF

        jelb   = pjelb (narea) ! Lower bound of internal domain
        jeub   = pjeub (narea) ! Upper bound of internal domain
        jesub  = pjesub(narea) ! Extent of internal domain

        jpj  = jesub + 2*jprecj ! jpj is the same for all domains - this is
                                ! what original decomposition did
        nlcj = jpj

      ! Unlike the East-West boundaries, the global domain does not include
      ! halo (rows) at the Northern and Southern edges. In fact, there is no
      ! cyclic boundary condition in the North-South direction so there are no
      ! halos at all at the edges of the global domain.
      IF( nbondj == -1 .AND. (.NOT. trimmed(sidx,narea)) )THEN
         ! Southern edge
         nldj = 1                ! Start of internal region (local coords)
         nlej = nldj + jesub - 1 ! Upper bound of int. region of sub-domain, local
         nlcj = nlej + jprecj
         jpj  = nlcj
      ELSE IF( nbondj ==  1 .AND. (.NOT. trimmed(nidx,narea)) )THEN
         ! Northern edge
         nldj = 1+jprecj       ! Start of internal region (local coords)
         nlej = nldj + jesub - 1
         nlcj = nlej
         jpj  = nlcj
         jpj = jpj + 1 ! Add one extra row for zero BC along N edge as
                       ! happens in orig. code when MOD(jpjglo,2)!=0
                       ! Many loops go up to j=jpjm1 so unless jpj>nlej
                       ! they won't cover the whole domain.
      ELSE IF( nbondj == 2)THEN
         ! Both boundaries are on the edges of the global domain

         IF( trimmed(sidx,narea) )THEN
            nldj = 1+jprecj
         ELSE
            nldj = 1
         END IF
         nlej = nldj + jesub - 1

         IF( trimmed(nidx,narea) )THEN
            nlcj = nlej + jprecj
            jpj  = nlcj
         ELSE
            nlcj = nlej
            jpj  = nlcj
         END IF
         jpj = jpj + 1 ! Add one extra row for zero BC along N edge as
                       ! happens in orig. code when MOD(jpjglo,2)!=0
      ELSE
         ! We have no external boundaries to worry about
         nldj = 1+jprecj    ! Lower bound of int. region of sub-domain, local
         nlej = nldj + jesub - 1
      END IF

      jpjm1 = jpj-1
      jpij  = jpi*jpj 


   END SUBROUTINE finish_partition


   SUBROUTINE mpp_ini_north
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_north  ***
      !!
      !! ** Purpose :   Initialize special communicator for north folding 
      !!      condition together with global variables needed in the mpp folding
      !!
      !! ** Method  : - Look for northern processors
      !!              - Put their number in nrank_north
      !!              - Create groups for the world processors and the north 
      !!                processors
      !!              - Create a communicator for northern processors
      !!
      !! ** output
      !!      njmppmax = njmpp for northern procs
      !!      ndim_rank_north = number of processors in the northern line
      !!      nrank_north (ndim_rank_north) = number  of the northern procs.
      !!      ngrp_world = group ID for the world processors
      !!      ngrp_north = group ID for the northern processors
      !!      ncomm_north = communicator for the northern procs.
      !!      north_root = number (in the world) of proc 0 in the northern comm.
      !!      nwidthmax = width of widest northern domain
      !!
      !! History :
      !!        !  03-09 (J.M. Molines, MPI only )
      !!        !  08-09 (A.R. Porter - for new decomposition)
      !!----------------------------------------------------------------------
      USE par_oce, ONLY: jperio, jpni
      USE exchmod, ONLY: nrank_north, north_root, ndim_rank_north, &
                         ncomm_north, ngrp_world, ngrp_north,      &
                         do_nfold, num_nfold_rows, nfold_npts
      USE dom_oce, ONLY: narea
      IMPLICIT none
#ifdef key_mpp_shmem
      CALL ctl_stop('STOP', ' mpp_ini_north not available in SHMEM' )
# elif key_mpp_mpi
      INTEGER :: ierr
      INTEGER :: jproc
      INTEGER :: ii,ji
      !!----------------------------------------------------------------------

      ! Look for how many procs on the northern boundary
      !
      ndim_rank_north = 0
      nwidthmax       = 0
      do_nfold        = .FALSE.

      IF (.NOT. (jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1) ) THEN
         ! No northern boundary to worry about
         RETURN
      END IF

      DO jproc=1,mppsize,1
         IF ( pjubext(jproc) ) THEN

            ! If trimming of dry land from sub-domains is enabled
            ! then check that this PE does actually have data to
            ! contribute to the N-fold. If trimming is not enabled
            ! then this condition will always be true for northern
            ! PEs.
            IF( pjeub(jproc) > (jpjglo - num_nfold_rows) )THEN

               ndim_rank_north = ndim_rank_north + 1

               ! and for the width of the widest northern domain...
               nwidthmax = MAX(nwidthmax, piesub(jproc))
            ENDIF

         END IF
      END DO
      nwidthmax = nwidthmax + 2*jpreci ! Allow for halos

      ! Allocate the right size to nrank_north
      !
      ALLOCATE(nrank_north(ndim_rank_north), nfold_npts(ndim_rank_north), &
               Stat=ierr)
      IF( ierr /= 0 )THEN
         CALL ctl_stop('STOP','mpp_ini_north: failed to allocate arrays')
      END IF

#if defined PARTIT_DEBUG
      IF(lwp)THEN
         WRITE(*,*) 'mpp_ini_north: no. of northern PEs = ',ndim_rank_north
         WRITE(*,*) 'mpp_ini_north: nwidthmax = ',nwidthmax
      END IF
#endif
      ! Fill the nrank_north array with proc. number of northern procs.
      ! Note : ranks start at 0 in MPI
      !
      ii=0
      DO ji = 1, mppsize, 1
         IF (  pjubext(ji)       .AND.          &
              (pjeub(ji) > (jpjglo - num_nfold_rows)) ) THEN
            ii=ii+1
            nrank_north(ii)=ji-1

            ! Flag that this PE does do North-fold (with trimming, checking
            ! npolj is no longer sufficient)
            IF(ji == narea) do_nfold = .TRUE.

#if defined NO_NFOLD_GATHER
            ! How many data points will this PE have to send for N-fold?

            ! No. of valid rows for n-fold = num_nfold_rows - <no. trimmed rows>
            !                              = num_nfold_rows - jpjglo + pjeub(ji)

            ! ARPDBG - could trim land-only rows/cols from this...
            nfold_npts(ii) = MAX(num_nfold_rows - jpjglo + pjeub(ji), 0) * &
                             ( nleit(ji) - nldit(ji) + 1 )
#endif 
         END IF
      END DO
      ! create the world group
      !
      CALL MPI_COMM_GROUP(mpi_comm_opa,ngrp_world,ierr)
      !
      ! Create the North group from the world group
      CALL MPI_GROUP_INCL(ngrp_world,ndim_rank_north,nrank_north, &
                          ngrp_north,ierr)

      ! Create the North communicator , ie the pool of procs in the north group
      !
      CALL MPI_COMM_CREATE(mpi_comm_opa,ngrp_north,ncomm_north,ierr)


      ! find proc number in the world of proc 0 in the north
      CALL MPI_GROUP_TRANSLATE_RANKS(ngrp_north,1,0,ngrp_world,north_root,ierr)

#endif

   END SUBROUTINE mpp_ini_north


   SUBROUTINE eval_partition( nx, ny, mask, score )

      ! Compute the cost function for the current partition
      !
      ! Assume that the time taken for a run is proportional
      ! to the maximum over processors of:
      !     w_processing * cost_processing
      !   + w_communications * cost_communications
      ! Assume further that cost_processing goes as
      !   (number of wet points) + f_proc * (number of dry points)
      ! (with f_proc << 1)
      ! and that cost_communications goes as
      !   (cost of intra-node communications) +
      !   f_comm * (cost of inter-node communications)
      ! (with f_comm << 1)
      !
      ! However, because of the possiblity of network contention,
      ! other factors may also matter, especially:
      !   total over sub-domains of halo points with off-node neighbours
      !   max over nodes of total off-node halo points and message counts
      !
      ! With this in mind, we construct the ansatz
      !  maximum over processors of {
      !     w_1 * (number of wet points)
      !   + w_2 * (number of dry points)
      !   + w_3 * (halo points with off-node neighbours)
      !   + w_4 * (halo points with on-node neighbours)
      !   + ...
      ! }
#if defined key_mpp_mpi
      USE lib_mpp,     ONLY: mppsize
#endif
      USE mapcomm_mod, ONLY: iprocmap, land
      IMPLICIT NONE
      !     Arguments
        INTEGER, INTENT(in) :: nx, ny
        INTEGER, INTENT(in) :: mask(nx,ny)
        REAL(wp), DIMENSION(pv_num_scores), INTENT(out) :: score
        !     Local variables
        INTEGER :: iproc, inode, i, j

        REAL(wp) :: score_wet, score_dry
        REAL(wp) :: score_pli, score_plx
        REAL(wp) :: score_pci, score_pcx
        REAL(wp) :: score_nli, score_nlx
        REAL(wp) :: score_nci, score_ncx
        REAL(wp) :: score_tli, score_tlx
        REAL(wp) :: score_tci, score_tcx

        REAL(wp) :: score_too_narrow

        REAL(wp) :: proc_overall_score
        REAL(wp) :: proc_comm_score
        REAL(wp) :: node_comm_score

        REAL(wp), PARAMETER :: weight_wet  =  1.00D0
        REAL(wp), PARAMETER :: weight_dry  =  0.9D0
 
        REAL(wp), PARAMETER :: weight_pli  =  0.01D0
        REAL(wp), PARAMETER :: weight_plx  =  0.20D0
        REAL(wp), PARAMETER :: weight_pci  =  0.50D0
        REAL(wp), PARAMETER :: weight_pcx  = 10.00D0

        REAL(wp), PARAMETER :: weight_nli  =  0.00D0
        REAL(wp), PARAMETER :: weight_nlx  =  0.20D0
        REAL(wp), PARAMETER :: weight_nci  =  0.00D0
        REAL(wp), PARAMETER :: weight_ncx  = 10.00D0
        
        REAL(wp), PARAMETER :: weight_tli  =  0.00D0
        REAL(wp), PARAMETER :: weight_tlx  =  0.01D0
        REAL(wp), PARAMETER :: weight_tci  =  0.00D0
        REAL(wp), PARAMETER :: weight_tcx  =  0.50D0

        INTEGER :: peer, last_peer

        ! Which node is each process on?
        ! Assumes that first nn_cpnode ranks are assigned to node 0,
        ! next nn_cpnode ranks are assigned to node 1, etc
        INTEGER, ALLOCATABLE :: node(:)

#if defined key_mpp_mpi

        ALLOCATE(node(nprocp))

        DO iproc=1, nprocp
           node(iproc) = (iproc-1)/nn_cpnode
        END DO

        ! Calculate maximum per processor score

        score(:) = 0.0_wp

        ! Total (over all processors) off node comms 
        score_tli  = 0.0_wp
        score_tlx  = 0.0_wp
        score_tci  = 0.0_wp
        score_tcx  = 0.0_wp

        ! Set this to pv_awful if any sub-domain is too narrow.
        score_too_narrow = 0.0_wp

        ! loop over nodes in job, counting from 0
        node_loop: DO inode=0, (nprocp-1)/nn_cpnode

        score_nli  = 0.0_wp
        score_nlx  = 0.0_wp
        score_nci  = 0.0_wp
        score_ncx  = 0.0_wp

        ! loop over processes in the node
        proc_loop: DO iproc=1+inode*nn_cpnode, &
                            MIN(nprocp,(inode+1)*nn_cpnode) 

           score_wet  = DBLE(pnactive(iproc))
           score_dry  = DBLE(piesub(iproc)*pjesub(iproc)-score_wet)

           score_pli  = 0.0_wp
           score_plx  = 0.0_wp
           score_pci  = 0.0_wp
           score_pcx  = 0.0_wp

           ! Things sometimes go wrong when a sub-domain has very
           ! narrow partitions (2 grid points or less).
           ! Prevent such problematic partitions from being selected.
           IF (      ((pieub(iproc)-pielb(iproc)) .LE. 2) &
                .OR. ((pjeub(iproc)-pjelb(iproc)) .LE. 2) ) THEN
              score_too_narrow = pv_awful
           END IF

           IF (.NOT. pjlbext(iproc)) THEN
              j=pjelb(iproc)
              IF (j .GT. 1) THEN
                 last_peer = -1
                 DO i=pielb(iproc),pieub(iproc)
                    IF (       (mask(i,j) .NE. land) &
                         .AND. (mask(i,j-1) .NE. land)) THEN
                       peer=iprocmap(i,j-1)
                       ! Total points involved in messages
                       IF (node(peer) .EQ. inode) THEN
                          score_pli = score_pli+1.0_wp
                       ELSE
                          score_plx = score_plx+1.0_wp
                       END IF
                       ! Total number of messages
                       IF (peer .NE. last_peer) THEN
                          last_peer = peer
                          IF (node(peer) .EQ. inode) THEN
                             score_pci = score_pci+1.0_wp
                          ELSE
                             score_pcx = score_pcx+1.0_wp
                          END IF
                       END IF
                    END IF
                 END DO
              END IF
           END IF
  
           IF (.NOT. pjubext(iproc)) THEN
            j=pjeub(iproc)
            IF (j .LT. ny) THEN
              last_peer = -1
              DO i=pielb(iproc),pieub(iproc)
                IF (      (mask(i,j)   .NE. land)  &
                    .AND. (mask(i,j+1) .NE. land)) THEN
                  peer=iprocmap(i,j+1)
!                 Total points involved in messages
                  IF (node(peer) .EQ. inode) THEN
                    score_pli = score_pli+1.0_wp
                  ELSE
                    score_plx = score_plx+1.0_wp
                  END IF
!                 Total number of messages
                  IF (peer .NE. last_peer) THEN
                    last_peer = peer
                    IF (node(peer) .EQ. inode) THEN
                      score_pci = score_pci+1.0_wp
                    ELSE
                      score_pcx = score_pcx+1.0_wp
                    END IF
                  END IF
                END IF
              END DO
            END IF
          END IF

          IF (.NOT. pilbext(iproc)) THEN
            i=pielb(iproc)
            IF (i .GT. 1) THEN
              last_peer = -1
              DO j=pjelb(iproc),pjeub(iproc)
                IF (      (mask(i,j)   .NE. land) &
                    .AND. (mask(i-1,j) .NE. land)) THEN
                  peer=iprocmap(i-1,j)
!                 Total points involved in messages
                  IF (node(peer) .EQ. inode) THEN
                    score_pli = score_pli+1.0_wp
                  ELSE
                    score_plx = score_plx+1.0_wp
                  END IF
!                 Total number of messages
                  IF (peer .NE. last_peer) THEN
                    last_peer = peer
                    IF (node(peer) .EQ. inode) THEN
                      score_pci = score_pci+1.0_wp
                    ELSE
                      score_pcx = score_pcx+1.0_wp
                    END IF
                  END IF
                END IF
              END DO
            END IF
          END IF

          IF (.NOT. piubext(iproc)) THEN
            i=pieub(iproc)
            IF (i .LT. nx) THEN
              last_peer = -1
              DO j=pjelb(iproc),pjeub(iproc)
                IF (      (mask(i,j)   .NE. land)  &
                    .AND. (mask(i+1,j) .NE. land)) THEN
                  peer=iprocmap(i+1,j)
                  ! Total points involved in messages
                  IF (node(peer) .EQ. inode) THEN
                    score_pli = score_pli+1.0_wp
                  ELSE
                    score_plx = score_plx+1.0_wp
                  END IF
                  ! Total number of messages
                  IF (peer .NE. last_peer) THEN
                    last_peer = peer
                    IF (node(peer) .EQ. inode) THEN
                      score_pci = score_pci+1.0_wp
                    ELSE
                      score_pcx = score_pcx+1.0_wp
                    END IF
                  END IF
                END IF
              END DO
            END IF
          END IF

          score_nli = score_nli + score_pli
          score_nlx = score_nlx + score_plx
          score_nci = score_nci + score_pci
          score_ncx = score_ncx + score_pcx

          proc_overall_score = weight_wet*score_wet &
                            + weight_dry*score_dry  &
                            + weight_pli*score_pli  &
                            + weight_plx*score_plx  &
                            + weight_pci*score_pci  &
                            + weight_pcx*score_pcx

          proc_comm_score    = weight_pli*score_pli &
                             + weight_plx*score_plx &
                             + weight_pci*score_pci &
                             + weight_pcx*score_pcx

          IF (score(pv_index_overall) < proc_overall_score) &
              score(pv_index_overall) = proc_overall_score
          IF (score(pv_index_pcomms) < proc_comm_score)     &
              score(pv_index_pcomms) = proc_comm_score
          IF (score(pv_index_wet) < score_wet) &
              score(pv_index_wet) = score_wet
          IF (score(pv_index_dry) < score_dry) &
              score(pv_index_dry) = score_dry
          IF (score(pv_index_pli) < score_pli) &
              score(pv_index_pli) = score_pli
          IF (score(pv_index_plx) < score_plx) &
              score(pv_index_plx) = score_plx
          IF (score(pv_index_pci) < score_pci) &
              score(pv_index_pci) = score_pci
          IF (score(pv_index_pcx) < score_pcx) &
              score(pv_index_pcx) = score_pcx
  
        END DO proc_loop

        score_tli = score_tli + score_nli
        score_tlx = score_tlx + score_nlx
        score_tci = score_tci + score_nci
        score_tcx = score_tcx + score_ncx

        node_comm_score    = weight_nli*score_nli &
                           + weight_nlx*score_nlx &
                           + weight_nci*score_nci &
                           + weight_ncx*score_ncx

        IF (score(pv_index_ncomms) < node_comm_score) &
            score(pv_index_ncomms) = node_comm_score
        IF (score(pv_index_nli) < score_nli) &
            score(pv_index_nli) = score_nli
        IF (score(pv_index_nlx) < score_nlx) &
            score(pv_index_nlx) = score_nlx
        IF (score(pv_index_nci) < score_nci) &
            score(pv_index_nci) = score_nci
        IF (score(pv_index_ncx) < score_ncx) &
            score(pv_index_ncx) = score_ncx

      END DO node_loop

      score(pv_index_tcomms) = weight_tli*score_tli &
                             + weight_tlx*score_tlx &
                             + weight_tci*score_tci &
                             + weight_tcx*score_tcx
      
      score(pv_index_tli) = score_tli
      score(pv_index_tlx) = score_tlx
      score(pv_index_tci) = score_tci
      score(pv_index_tcx) = score_tcx

      score(pv_index_overall) = score(pv_index_overall) &
                              + score(pv_index_ncomms)  &
                              + score(pv_index_tcomms)  &
                              + score_too_narrow

      DEALLOCATE(node)

#endif

     END SUBROUTINE eval_partition


     SUBROUTINE calc_perms( nfactors, factors,     &
                            ndf, df, multiplicity, &
                            nperms )
      IMPLICIT NONE

!     Subroutine arguments
!    
!     Number of factors (including repetitions)
      INTEGER, INTENT(in) :: nfactors

!     Factors (in descending order)
      INTEGER, INTENT(in) :: factors(nfactors)

!     Number of distinct factors
      INTEGER, INTENT(out) :: ndf

!     Distinct factors (in ascending order)
      INTEGER :: df(nfactors)

!     Multiplicity of each distinct factor
      INTEGER :: multiplicity(nfactors)

!     Number of distinct permutations
      INTEGER, INTENT(out) :: nperms

!     Local variables

      INTEGER :: product
      INTEGER :: i, j

      product = factors(nfactors)
      ndf = 1
      df(:) = 0
      df(1) = factors(nfactors)
      multiplicity(:) = 0
      multiplicity(1) = 1
      
      DO i=nfactors-1,1,-1
        product = product*factors(i)
        IF (factors(i) .EQ. df(ndf)) THEN
          multiplicity(ndf) = multiplicity(ndf)+1
        ELSE
          ndf = ndf+1
          df(ndf) = factors(i)
          multiplicity(ndf) = 1
        END IF
      END DO
!     write (*,*) 'number: ', product

!     A careful code would try to avoid overflow here
      nperms = 1
      DO i=1, nfactors
        nperms = nperms*i
      END DO
      DO i=1, ndf
        DO j=1, multiplicity(i)
          nperms = nperms/j
        END DO
      END DO

!     At this point, nperms is the number of distinct permutations
!     of the factors provided. But each of these permutations can
!     be split between x and y in (nfactors+1) ways, e.g.
!       x=(2,2,3), y=()
!       x=(2,3),   y=(2)
!       x=(3),     y=(2,2)
!       x=(),      y=(2,2,3)

      nperms = nperms*(nfactors+1)
      IF (lwp) THEN
        WRITE (numout,*) 'distinct factorisations: ', nperms
      END IF
          
      END SUBROUTINE calc_perms


    
      SUBROUTINE get_perm_factors( iperm, nf, ndf, df, m, &
                                   fx, nfx, fy, nfy,      &
                                   prodx, prody, printit )
         USE dom_oce, ONLY: narea
         IMPLICIT NONE
         !!------------------------------------------------------------------
         !     Our goal is to visit a particular permutation, 
         !     number perm ( 0 <= perm <= nperms-1 )
         !     of nf things, of ndf distinct values (actually the prime
         !     factors of number of MPI tasks), each of which can be repeated
         !     with multiplicity m_i
         !     assert nf = sum_{i=1..ndf} m(i)
         !     
         !     We don't care about lexicographic ordering, but we do
         !     need to ensure that we don't visit any permutation twice,
         !     in a loop over 0..nperms-1.
         !     Textbook methods typically assume that all the things being
         !     permuted are distinct.
         !
         !     We use what I call a nested variable radix method.
         !
         !     Stephen Pickles, STFC
         !     Taken from POLCOMS code by Andrew Porter.
         !!------------------------------------------------------------------
         !     Arguments
         INTEGER, INTENT(in)  :: iperm, nf, ndf
         INTEGER, INTENT(in), DIMENSION(ndf) :: df, m
         INTEGER, INTENT(out), DIMENSION(nf) :: fx, fy
         INTEGER, INTENT(out) :: nfx, nfy
         INTEGER, INTENT(out) :: prodx, prody
         LOGICAL, INTENT(in)  :: printit
         !
         !     Local variables
         !
         INTEGER :: perm, split
         INTEGER, DIMENSION(nf) :: bin, a
         INTEGER, DIMENSION(ndf) :: ways
         INTEGER, DIMENSION(0:ndf) :: root, representation
         INTEGER :: i, j, k, v, p, q, r
         INTEGER :: unfilled, pm, rm
         INTEGER :: myinst
         LOGICAL, PARAMETER :: debug=.FALSE.
         !!------------------------------------------------------------------

         ! MPI rank of this process
         myinst = narea - 1

         perm = iperm / (nf+1)
         !     Where to split between x and y 
         split = MOD(iperm, (nf+1))

         !     interpret ways(i) is the number of ways of distributing 
         !     m(i) copies of the i'th distinct factor into the remaining
         !     bins
         k = nf
         DO i=1,ndf
            ways(i) = k
            DO j=2,m(i)
               ways(i) = ways(i)*(k-j+1)/j
            END DO
            k = k-m(i)
         END DO

         !     compute (outer) radices
         !     root is the variable radix basis corresponding to ways
         !     root(ndf) is always 1
         root(ndf) = 1
         root(0:ndf-1) = ways(1:ndf)
         DO i=ndf-1,0,-1
            root(i)=root(i)*root(i+1)
         END DO

         bin(:) = 0
         unfilled = nf

         r = perm
         !     Next line is valid as long as perm < nperms
         representation(0) = 0
         DO i=1, ndf
            p = r/root(i)
            r = r - p*root(i)
            representation(i) = p

            !       At this point, we are considering distinct ways to
            !       distribute m(i) copies of the i'th distinct factor into
            !       the unfilled remaining bins. We want to select the p'th one.

            DO j=1,unfilled
               a(j) = j
            END DO
            q = 0
            find_p: DO
               IF (q .GE. p) EXIT find_p

               k=m(i)
               k_loop: DO
                  IF (a(k) .EQ. (unfilled - m(i) + k)) THEN
                     k=k-1
                  ELSE
                     EXIT k_loop
                  END IF
               END DO k_loop
               a(k) = a(k) + 1
               DO v=k+1,m(i)
                  a(v) = a(k) + v - k
               END DO
               q = q+1
            END DO find_p

            IF (debug) THEN
               WRITE (*,'(A3)',advance='no') 'a=('
               DO j=1, m(i)-1
                  WRITE (*,'(I3,A1)',advance='no') a(j), ','
               END DO
               WRITE (*,'(I3,A1)') a(m(i)), ')'
            END IF

            DO j=1, m(i)
               pm=a(j)-j+1

               ! put this factor in the pm'th empty bin
               v = 1
               find_bin: DO k=1, nf
                  IF (bin(k) .EQ. 0) THEN
                     IF (v .EQ. pm) THEN
                        bin(k) = df(i)
                        unfilled = unfilled-1
                        EXIT find_bin
                     ELSE
                        v=v+1
                     END IF
                  END IF
               END DO find_bin

            END DO
            
         END DO

         !     We have identified the (perm)th distinct permutation,
         !     but we still need to split the factors between x and y.
         fx(:) = 0
         prodx = 1
         DO i=1,split
            fx(i)=bin(i)
            prodx=prodx*fx(i)
         END DO
         nfx=split

         fy(:) = 0 
         prody = 1
         j=1
         DO i=split+1,nf
            fy(j)=bin(i)
            prody=prody*fy(j)
            j=j+1
         END DO
         nfy=nf-nfx

         IF (printit) THEN

            WRITE (6,'(A14,I6,A1,I6,A2)',advance='no') &
                      'factorisation[', myinst, ']', iperm, ' ('
            DO k=1,ndf-1
               WRITE (6,'(I4,A1)',advance='no') representation(k), ','
            END DO
            WRITE (6,'(I4,A1)',advance='no') representation(ndf), ')'
     
            CALL print_factors(6,nfx,fx,nfy,fy)

         END IF

      END SUBROUTINE get_perm_factors


      SUBROUTINE print_factors(lu,nfx,fx,nfy,fy)
         IMPLICIT NONE
         INTEGER, INTENT(in) :: lu
         INTEGER, INTENT(in) :: nfx, nfy
         INTEGER, INTENT(in) :: fx(nfx), fy(nfy)
         INTEGER :: k

         IF (nfx .GT. 0) THEN
            WRITE (lu,'(A1)',advance='no') ' '
            DO k=1,nfx-1
               WRITE (lu,'(I4,A1)',advance='no') fx(k), ','
            END DO
            WRITE (lu,'(I4)',advance='no') fx(nfx)
         END IF
         WRITE (lu,'(A1)',advance='no') ':'
         IF (nfy .GT. 0) THEN
            DO k=1,nfy-1
               WRITE (lu,'(I4,A1)',advance='no') fy(k), ','
            END DO
            WRITE (lu,'(I4)',advance='no') fy(nfy)
         END IF
         WRITE (lu,'(A1)') ' '

      END SUBROUTINE print_factors


      CHARACTER(len=16) FUNCTION num_to_string(number)
         INTEGER, INTENT(in) :: number

         CHARACTER*16 :: outs
      
         WRITE (outs,'(I15)') number
         num_to_string = ADJUSTL(outs)

      END FUNCTION num_to_string


      SUBROUTINE factor_string(fstr,nfx,fx,nfy,fy)
         IMPLICIT NONE
         CHARACTER*256, INTENT(out) :: fstr
         INTEGER, INTENT(in) :: nfx, nfy
         INTEGER, INTENT(in) :: fx(nfx), fy(nfy)
         !!----------------------------------------------------------------------
         !!----------------------------------------------------------------------
         INTEGER :: k

         fstr = ' '
         IF (nfx .GT. 0) THEN
            DO k=1,nfx-1
               fstr = TRIM(fstr)//TRIM(num_to_string(fx(k)))//'x'
            END DO
            fstr = TRIM(fstr)//TRIM(num_to_string(fx(nfx)))
         END IF
         fstr = TRIM(fstr)//'-'
         IF (nfy .GT. 0) THEN
            DO k=1,nfy-1
               fstr = TRIM(fstr)//TRIM(num_to_string(fy(k)))//'x'
            END DO
            fstr = TRIM(fstr)//TRIM(num_to_string(fy(nfy)))
         END IF
      END SUBROUTINE factor_string


    SUBROUTINE write_partition_map(depth)
       IMPLICIT NONE
       INTEGER, DIMENSION(:,:), INTENT(in) :: depth
       !!----------------------------------------------------------------------
       !!     Writes an ASCII and PPM format map of the domain decomposition
       !!     to separate files.
       !!----------------------------------------------------------------------
       ! Locals
       INTEGER, ALLOCATABLE, DIMENSION(:,:) :: map
       INTEGER :: nx, ny
       INTEGER :: iproc, i, j, icol
       INTEGER :: lumapout
       CHARACTER(LEN=6)  :: mode
       CHARACTER(LEN=40) :: mapout
       INTEGER, PARAMETER :: ncol=15
       INTEGER, DIMENSION(3,-2:ncol-1) :: rgbcol
       !!----------------------------------------------------------------------

       nx = jpiglo
       ny = jpjglo

       ! Generate an integer map of the partitioning.

       ALLOCATE (map(nx,ny))
       map = -1
       DO iproc=1,nprocp
          DO j=pjelb(iproc),pjeub(iproc)
             DO i=pielb(iproc),pieub(iproc)
                IF ( depth(i,j) == 0 ) THEN

                   ! Identify the land using -2 which is set to black
                   ! in the colour map below.
                   map(i,j) = -2
                ELSE

                   ! Identify to which process the point is assigned.
                   map(i,j) = iproc-1
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       ! Write the map to a file for plotting.

       IF ( lwp ) THEN

          ! ASCII format map file.

          lumapout = 9
          mode = 'simple'
          IF ( nprocp.LT.10 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i1,''.dat'')') mode,nprocp
          ELSEIF ( nprocp.LT.100 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i2,''.dat'')') mode,nprocp
          ELSEIF ( nprocp.LT.1000 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i3,''.dat'')') mode,nprocp
          ELSE
             WRITE (mapout,'(''Map_'',a6,''_'',i4,''.dat'')') mode,nprocp
          ENDIF
          OPEN (lumapout,file=mapout)
          WRITE (lumapout,*) nx,ny
          DO j=1,ny
          !            write (lumapout,'(15i5)') (map(i,j),i=1,ny)
             DO i=1,nx,1
                WRITE (lumapout,'(3i5)') i ,j, map(i,j)
             END DO
          ENDDO
          CLOSE (lumapout)

          ! PPM format map file.

          lumapout = 10
          mode = 'partrk'

          IF ( nprocp.LT.10 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i1,''.ppm'')') mode,nprocp
          ELSEIF ( nprocp.LT.100 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i2,''.ppm'')') mode,nprocp
          ELSEIF ( nprocp.LT.1000 ) THEN
             WRITE (mapout,'(''Map_'',a6,''_'',i3,''.ppm'')') mode,nprocp
          ELSE
             WRITE (mapout,'(''Map_'',a6,''_'',i4,''.ppm'')') mode,nprocp
          ENDIF
          OPEN (lumapout,file=mapout)

          ! PPM magic number.
          ! Comment line
          ! width and height of image (same as that of the domain)
          ! Maximum colour value.

          WRITE (lumapout,'(a)') 'P3'
          WRITE (lumapout,'(a)') '# '//mapout
          WRITE (lumapout,'(2i6)') nx,ny
          WRITE (lumapout,'(i6)') 255

          ! Define RGB colours. 0 is grey for the land. 1-16 for the sub-domains.
          ! When there are more than 16 sub-domains the colours wrap around.

          rgbcol(:,-2) = (/   0,   0,   0 /)
          rgbcol(:,-1) = (/ 170, 170, 170 /)
          rgbcol(:, 0) = (/   0,   0, 255 /) ! dark blue
          rgbcol(:, 1) = (/   0,  85, 255 /) ! blue
          rgbcol(:, 2) = (/   0, 170, 255 /) ! pale blue
          rgbcol(:, 3) = (/   0, 255, 255 /) ! cyan
          rgbcol(:, 4) = (/   0, 170,   0 /) ! dark green
          rgbcol(:, 5) = (/   0, 255,   0 /) ! green
          rgbcol(:, 6) = (/   0, 255, 170 /) ! blue-green
          rgbcol(:, 7) = (/ 128, 255,   0 /) ! yellow-green
          rgbcol(:, 8) = (/ 128, 170,   0 /) ! khaki
          rgbcol(:, 9) = (/ 255, 255,   0 /) ! yellow
          rgbcol(:,10) = (/ 255,  85,   0 /) ! orange
          rgbcol(:,11) = (/ 255,   0,  85 /) ! pink-ish
          rgbcol(:,12) = (/ 128,   0, 255 /) ! violet
          rgbcol(:,13) = (/ 255,   0, 255 /) ! magenta
          rgbcol(:,14) = (/ 170,   0, 128 /) ! purple
          !ma     rgbcol(:,15) = (/ 255,   0,  85 /) ! red

          ! Write out the RGB pixels, one per point in the domain.

          DO j=ny,1,-1
             DO i=1,nx
                IF ( map(i,j).LT.0 ) THEN
                   icol = map(i,j)
                ELSE
                   icol = MOD(map(i,j),ncol)
                ENDIF
                WRITE (lumapout,'(3i4)')  &
                     rgbcol(1,icol),rgbcol(2,icol),rgbcol(3,icol)
             ENDDO
          ENDDO
          CLOSE (lumapout)
       ENDIF ! (lwp)

       DEALLOCATE (map)

    END SUBROUTINE write_partition_map


    SUBROUTINE smooth_global_bathy(inbathy, imask)
       USE dom_oce
       USE domzgr, ONLY: rn_sbot_min, rn_sbot_max, rn_theta, rn_thetb, &
                         rn_rmax, ln_s_sigma, rn_bb, rn_hc, fssig1, &
                         namzgr_sco
       USE in_out_manager, ONLY: numnam
       IMPLICIT none
       !!----------------------------------------------------------------------
       !!                      Routine smooth_global_bathy
       !!   Replicates the smoothing done on the decomposed domain in zgr_sco()
       !!   in domzgr.F90. However, here the domain is NOT decomposed and
       !!   every processor performs an identical smoothing on the whole model
       !!   bathymetry. This is to ensure that the domain decomposition
       !!   is done using a mask that is the same as that which is eventually
       !!   computed after zgr_sco() has been called. (The smoothing process
       !!   below can (erroneously) change whether grid points are wet or dry.)
       !!----------------------------------------------------------------------
       REAL(wp), INTENT(inout), DIMENSION(:,:) :: inbathy ! The bathymetry to 
                                                          ! be smoothed
       INTEGER, INTENT(inout), DIMENSION(:,:)  :: imask   ! Mask holding index of
                                                          ! bottom level
       ! Locals
       INTEGER  :: ji, jj, jk, jl, ierr
       INTEGER  :: iip1, ijp1, iim1, ijm1   ! temporary integers
       INTEGER  :: x_size, y_size
       REAL(wp) :: zrmax, zri, zrj, zcoeft
       REAL(wp), PARAMETER :: TOL_ZERO = 1.0E-20_wp ! Any value less than 
                                                    ! this assumed zero
       REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: zenv, ztmp, zmsk, zbot, &
                                                  zscosrf, zhbatt
       REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zgsigt3, zgdept
       !
       !!----------------------------------------------------------------------

       ! Do this because we've not decomposed the domain yet and therefore
       ! jpi,jpj,nlc{i,j} etc. are not set.
       x_size = SIZE(inbathy, 1)
       y_size = SIZE(inbathy, 2)

       ALLOCATE(zenv(x_size,y_size), ztmp(x_size,y_size), zmsk(x_size,y_size), &
                zbot(x_size,y_size), zgdept(x_size,y_size,jpkdta), zhbatt(x_size, y_size), &
                zscosrf(x_size,y_size), zgsigt3(x_size,y_size,jpkdta), Stat=ierr)
       IF( ierr /= 0 ) THEN
          CALL ctl_stop('smooth_global_bathy: ERROR - failed to allocate workspace arrays')
          RETURN
       ENDIF

       REWIND( numnam )              ! Read Namelist namzgr_sco : sigma-stretching
                                     ! parameters
       READ  ( numnam, namzgr_sco )

       zscosrf(:,:) = 0._wp            ! ocean surface depth (here zero: no under ice-shelf sea)
       zbot(:,:) = inbathy(:,:)        ! ocean bottom depth
       !                               ! set maximum ocean depth
       inbathy(:,:) = MIN( rn_sbot_max, inbathy(:,:) )

       WHERE( inbathy(:,:) > TOL_ZERO ) inbathy(:,:) = MAX( rn_sbot_min, inbathy(:,:) )

       ! use r-value to create hybrid coordinates
       zenv(:,:) = MAX( inbathy(:,:), rn_sbot_min )
       ! 
       ! Smooth the bathymetry (if required)
       !
       jl = 0
       zrmax = 1._wp
       !                                                     ! ================ !
       DO WHILE( jl <= 10000 .AND. zrmax > rn_rmax )         !  Iterative loop  !
          !                                                  ! ================ !
          jl = jl + 1
          zrmax = 0._wp
          zmsk(:,:) = 0._wp

          DO jj = 1, y_size
             DO ji = 1, x_size
                iip1 = MIN( ji+1, x_size )      ! force zri = 0 on last line (ji=ncli+1 to jpi)
                ijp1 = MIN( jj+1, y_size )      ! force zrj = 0 on last row  (jj=nclj+1 to jpj)
                zri = ABS( zenv(iip1,jj  ) - zenv(ji,jj) ) / ( zenv(iip1,jj  ) + zenv(ji,jj) )
                zrj = ABS( zenv(ji  ,ijp1) - zenv(ji,jj) ) / ( zenv(ji  ,ijp1) + zenv(ji,jj) )
                zrmax = MAX( zrmax, zri, zrj )

                IF( zri > rn_rmax )   zmsk(ji  ,jj  ) = 1._wp
                IF( zri > rn_rmax )   zmsk(iip1,jj  ) = 1._wp
                IF( zrj > rn_rmax )   zmsk(ji  ,jj  ) = 1._wp
                IF( zrj > rn_rmax )   zmsk(ji  ,ijp1) = 1._wp
            END DO
         END DO

         !
         IF(lwp)WRITE(numout,"('smooth_global_bathy : iter=',I5,' rmax=',F8.4,' nb of pt= ',I8)") &
                                                         jl, zrmax, INT( SUM(zmsk(:,:) ) )
         !

         ! Copy current surface before next smoothing iteration 
         ztmp(:,:) = zenv(:,:)

         DO jj = 1, y_size
            DO ji = 1, x_size
               iip1 = MIN( ji+1, x_size )     ! last  line (ji=nlci)
               ijp1 = MIN( jj+1, y_size )     ! last  raw  (jj=nlcj)
               iim1 = MAX( ji-1,  1  )      ! first line (ji=nlci)
               ijm1 = MAX( jj-1,  1  )      ! first raw  (jj=nlcj)
               IF( zmsk(ji,jj) == 1._wp ) THEN
                ztmp(ji,jj) =   (                                                                                  &
             &    zenv(iim1,ijp1)*zmsk(iim1,ijp1) + zenv(ji,ijp1)*zmsk(ji,ijp1) + zenv(iip1,ijp1)*zmsk(iip1,ijp1)  &
             &  + zenv(iim1,jj  )*zmsk(iim1,jj  ) + zenv(ji,jj  )*    2._wp     + zenv(iip1,jj  )*zmsk(iip1,jj  )  &
             &  + zenv(iim1,ijm1)*zmsk(iim1,ijm1) + zenv(ji,ijm1)*zmsk(ji,ijm1) + zenv(iip1,ijm1)*zmsk(iip1,ijm1)  &
             &                  ) / (                                                                              &
             &                    zmsk(iim1,ijp1) +               zmsk(ji,ijp1) +                 zmsk(iip1,ijp1)  &
             &  +                 zmsk(iim1,jj  ) +                   2._wp     +                 zmsk(iip1,jj  )  &
             &  +                 zmsk(iim1,ijm1) +               zmsk(ji,ijm1) +                 zmsk(iip1,ijm1)  &
             &                      )
               ENDIF
            END DO
         END DO
         !
         DO jj = 1,y_size
            DO ji = 1,x_size
               IF( zmsk(ji,jj) >= 1._wp-TOL_ZERO ) zenv(ji,jj) = MAX( ztmp(ji,jj), inbathy(ji,jj) )
            END DO
         END DO
         !
         !                                                ! ================ !
      END DO                                              !     End loop     !
      !                                                   ! ================ !
      !
      !                                        ! envelop bathymetry saved in zhbatt
      zhbatt(:,:) = zenv(:,:) 
      ! gphit calculated in nemo_init->dom_init->dom_hgr and dom_hgr requires that 
      ! partitioning already done. Could repeat its calculation here but since AMM doesn't
      ! require it we leave it out for the moment ARPDBG
      CALL ctl_warn( ' ARPDBG - NOT checking whether s-coordinates are tapered in vicinity of the Equator' )
!!$      IF( MINVAL( gphit(:,:) ) * MAXVAL( gphit(:,:) ) <= 0._wp ) THEN
!!$         CALL ctl_warn( ' s-coordinates are tapered in vicinity of the Equator' )
!!$         DO jj = 1, jpj
!!$            DO ji = 1, jpi
!!$               ztaper = EXP( -(gphit(ji,jj)/8._wp)**2 )
!!$               hbatt(ji,jj) = rn_sbot_max * ztaper + hbatt(ji,jj) * ( 1._wp - ztaper )
!!$            END DO
!!$         END DO
!!$      ENDIF

      ! Subtract off rn_sbot_min so can check for land using zenv = LAND (0)
      inbathy(:,:) = zenv(:,:) - rn_sbot_min


      !                                            ! =======================
      !                                            !   s-ccordinate fields     (gdep., e3.)
      !                                            ! =======================
      !
      ! non-dimensional "sigma" for model level depth at w- and t-levels

      IF( ln_s_sigma ) THEN        ! Song and Haidvogel style stretched sigma for depths
         !                         ! below rn_hc, with uniform sigma in shallower waters
         DO ji = 1, x_size
            DO jj = 1, y_size

               IF( zhbatt(ji,jj) > rn_hc ) THEN    !deep water, stretched sigma
                  DO jk = 1, jpk
                     zgsigt3(ji,jj,jk) = -fssig1( REAL(jk,wp)       , rn_bb )
                  END DO
               ELSE ! shallow water, uniform sigma
                  DO jk = 1, jpk
                     zgsigt3(ji,jj,jk) = ( REAL(jk-1,wp) + 0.5_wp ) / REAL(jpk-1,wp)
                  END DO
               ENDIF
               !
               DO jk = 1, jpk
                  zcoeft = ( REAL(jk,wp) - 0.5_wp ) / REAL(jpkm1,wp)
                  zgdept (ji,jj,jk) = zscosrf(ji,jj) + (zhbatt(ji,jj)-rn_hc)*zgsigt3(ji,jj,jk)+rn_hc*zcoeft
               END DO
               !
            END DO   ! for all jj's
         END DO    ! for all ji's
      ELSE
         CALL ctl_stop('STOP', &
                       'partition_mod::smooth_global_bathy() only supports ln_s_sigma = .TRUE. currently!')
      END IF

      ! HYBRID scheme
      DO jj = 1, y_size
         DO ji = 1, x_size
            DO jk = 1, jpkm1
               IF( zbot(ji,jj) >= zgdept(ji,jj,jk) )  imask(ji,jj) = MAX( 2, jk )
               IF( zbot(ji,jj) == 0._wp           )   imask(ji,jj) = 0
            END DO
         END DO
      END DO

      ! Dump to file for debugging ARPDBG
      IF(lwp)THEN
         OPEN(UNIT=1098, FILE='smoothed_bathy.dat', STATUS='REPLACE', &
              ACTION='WRITE', IOSTAT=jj)
         IF(jj == 0)THEN
            DO jj = 1, y_size
               DO ji = 1, x_size
                  WRITE (1098,"(I4,1x,I4,3(E14.4,1x),I4)") ji, jj, &
                        inbathy(ji,jj),             zbot(ji,jj),   &
                        inbathy(ji,jj)-zbot(ji,jj), imask(ji,jj)
               END DO
               WRITE (1098,*)
            END DO
            CLOSE(1098)
         END IF
      END IF

    END SUBROUTINE smooth_global_bathy


    SUBROUTINE global_bot_level(imask)
      USE par_oce, ONLY: jperio
      IMPLICIT none
      !!----------------------------------------------------------------------
      !! Compute the deepest level for any of the u,v,w or T grids. (Code
      !! taken from zgr_bot_level() and intermediate arrays for U and V
      !! removed.)
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(inout) :: imask
      ! Locals
      INTEGER :: ji, jj
      INTEGER :: x_size, y_size

       ! Do this because we've not decomposed the domain yet and therefore
       ! jpi,jpj,nlc{i,j} etc. are not set.
       x_size = SIZE(imask, 1)
       y_size = SIZE(imask, 2)

      imask(:,:) = MAX( imask(:,:) , 1 )  ! bottom k-index of T-level (=1 over land)

      !
      ! Compute and store the deepest bottom k-index of any grid-type at 
      ! each grid point.
      ! For use in removing data below ocean floor from halo exchanges.
      DO jj = 1, y_size-1
         DO ji = 1, x_size-1
            imask(ji,jj) = MAX(imask(ji,jj)+1,                           & ! W (= T-level + 1)
                               MIN(  imask(ji+1,jj  ) , imask(ji,jj)  ), & ! U
                               MIN(  imask(ji  ,jj+1) , imask(ji,jj)  ) )  ! V
         END DO
         imask(x_size,jj) = imask(x_size-1,jj)
      END DO

      ! Check on jperio because we've not set cyclic_bc in mapcomms yet
      IF(jperio == 1 .OR. jperio == 4 .OR. jperio == 6)THEN
         ! Impose global cyclic boundary conditions on the array holding the
         ! deepest level
         imask(1,:)      = imask(x_size - 1, :)
         imask(x_size,:) = imask(2,:)
      END IF

      ! Dump to file for debugging ARPDBG
      IF(lwp)THEN
         OPEN(UNIT=1098, FILE='bathy_bottom.dat', STATUS='REPLACE', &
              ACTION='WRITE', IOSTAT=jj)
         IF(jj == 0)THEN
            DO jj = 1, y_size
               DO ji = 1, x_size
                  WRITE (1098,"(I4,1x,I4,1x,I4)") ji, jj, imask(ji,jj)
               END DO
               WRITE (1098,*)
            END DO
            CLOSE(1098)
         END IF
      END IF

    END SUBROUTINE global_bot_level


    SUBROUTINE read_partition(ierr)
      USE par_oce, ONLY: jpni, jpnj, jpnij
      USE mapcomm_mod, ONLY: eidx, widx, sidx, nidx, trimmed, &
                             pilbext, piubext, pjlbext, pjubext
      IMPLICIT none
      INTEGER, INTENT(out) :: ierr
      ! Locals
      INTEGER, PARAMETER :: funit = 1099
      INTEGER :: idom, ndom
      CHARACTER(len=200) :: linein
      !======================================================================

      ierr = 0

      OPEN(UNIT=funit, file='partition.dat', status='OLD', &
           ACTION='READ', IOSTAT=ierr)
      IF(ierr /= 0)THEN
         CALL ctl_warn('read_partition: failed to read partitioning from file - will calculate it instead.')
         RETURN
      END IF

      ! Number of procs in x and y
      CALL read_next_line(funit, linein, ierr)
      READ(linein,FMT=*) jpni, jpnj

      ! Store their product
      jpnij = jpni*jpnj

      ! Check that the implied number of PEs matches that
      ! in our MPI world
      ndom = jpni*jpnj
      IF(ndom /= mppsize)THEN
         CALL ctl_stop('STOP', &
                       'read_partition: no. of PEs specified in partition.dat does not match no. of PEs in use by this job.')
      END IF

      ! Read the description of each sub-domain
      domains: DO idom = 1, ndom, 1

         ! Coordinates of bottom-left (SW) corner of domain
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) pielb(idom), pjelb(idom)
         ! Top-right (NE) corner
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) pieub(idom), pjeub(idom)
         ! Whether this domain has external boundaries and has been trimmed
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) pilbext(idom), trimmed(widx,idom)
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) piubext(idom), trimmed(eidx,idom)
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) pjlbext(idom), trimmed(sidx,idom)
         CALL read_next_line(funit, linein, ierr)
         READ(linein,FMT=*) pjubext(idom), trimmed(nidx,idom)

         piesub(idom) = pieub(idom) - pielb(idom) + 1
         pjesub(idom) = pjeub(idom) - pjelb(idom) + 1

      END DO domains

      ! All done - close the file
      CLOSE(UNIT=funit)

      CALL finish_partition(fromFile=.TRUE.)

    END SUBROUTINE read_partition

    !===================================================================

    SUBROUTINE write_partition
      USE par_oce, ONLY: jpni, jpnj
      USE mapcomm_mod, ONLY: eidx, widx, sidx, nidx, trimmed,    &
                             pjubext, pjlbext, piubext, pilbext, &
                             pielb, pieub, pjelb, pjeub
      IMPLICIT none
      INTEGER, PARAMETER :: funit = 1099
      INTEGER :: ierr
      INTEGER :: idom

      ! Only PE 0 (narea==1) writes this file
      IF(narea /= 1) RETURN

      OPEN(UNIT=funit, file='partition.dat.new', status='REPLACE', &
           ACTION='WRITE', IOSTAT=ierr)
      IF(ierr /= 0)THEN
         CALL ctl_warn('write_partition: failed to write partition description to file.')
         RETURN
      END IF
      WRITE(UNIT=funit,FMT="('#  jpni  jpnj')")
      WRITE(UNIT=funit,FMT="(I5,1x,I5)") jpni, jpnj

      DO idom = 1, mppsize, 1
         WRITE(UNIT=funit,FMT="('# Domain: ',I5)") idom
         IF(idom==1)WRITE(UNIT=funit,FMT="('# Lower bounds: x  y')")
         WRITE(UNIT=funit,FMT="(I5,1x,I5)") pielb(idom), pjelb(idom)
         IF(idom==1)WRITE(UNIT=funit,FMT="('# Upper bounds: x  y')")
         WRITE(UNIT=funit,FMT="(I5,1x,I5)") pieub(idom), pjeub(idom)
         IF(idom==1)WRITE(UNIT=funit,FMT="('# x: Lower bound external, trimmed')")
         WRITE(UNIT=funit,FMT="(L5,1x,L5)") pilbext(idom), trimmed(widx,idom)
         IF(idom==1)WRITE(UNIT=funit,FMT="('# x: Upper bound external, trimmed')")
         WRITE(UNIT=funit,FMT="(L5,1x,L5)") piubext(idom), trimmed(eidx,idom)
         IF(idom==1)WRITE(UNIT=funit,FMT="('# y: Lower bound external, trimmed')")
         WRITE(UNIT=funit,FMT="(L5,1x,L5)") pjlbext(idom), trimmed(sidx,idom)
         IF(idom==1)WRITE(UNIT=funit,FMT="('# y: Upper bound external, trimmed')")
         WRITE(UNIT=funit,FMT="(L5,1x,L5)") pjubext(idom), trimmed(nidx,idom)
      END DO

      CLOSE(UNIT=funit)

    END SUBROUTINE write_partition

    SUBROUTINE read_next_line(funit, linein, ierr)
      IMPLICIT none
      !!------------------------------------------------------------------
      INTEGER,            INTENT( in) :: funit  ! Unit no. to read
      CHARACTER(len=200), INTENT(out) :: linein ! String containing next
                                                ! non-comment line in file
      INTEGER,            INTENT(out) :: ierr   ! Error flag (0==OK)
      !!------------------------------------------------------------------

      ierr = 0

      READ(UNIT=funit,FMT="(200A)") linein

      ! Comment lines begin with '#'. Skip those plus any blank
      ! lines...
      DO WHILE( INDEX( TRIM(ADJUSTL(linein)),'#') /= 0 .OR. &
                LEN_TRIM(linein) == 0 )
         READ(UNIT=funit,FMT="(200A)") linein   
      END DO

      WRITE(*,*)'returning linein >>'//linein//'<<'

    END SUBROUTINE read_next_line

END MODULE partition_mod
