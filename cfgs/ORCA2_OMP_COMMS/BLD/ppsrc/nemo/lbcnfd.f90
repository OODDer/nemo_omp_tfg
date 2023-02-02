












MODULE lbcnfd
   !!======================================================================
   !!                       ***  MODULE  lbcnfd  ***
   !! Ocean        : north fold  boundary conditions
   !!======================================================================
   !! History :  3.2  ! 2009-03  (R. Benshila)  Original code 
   !!            3.5  ! 2013-07  (I. Epicoco, S. Mocavero - CMCC) MPP optimization
   !!            4.0  ! 2017-04  (G. Madec) automatique allocation of array argument (use any 3rd dimension)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   lbc_nfd       : generic interface for lbc_nfd_sp and lbc_nfd_dp routines that is doing the north fold in a non-mpi case 
   !!   mpp_nfd       : generic interface for mpp_nfd_sp and mpp_nfd_dp routines that will use lbc_nfd directly or indirectly
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain 
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE MPI

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_nfd            ! called by mpp_nfd, lbc_lnk_pt2pt or lbc_lnk_neicoll
      MODULE PROCEDURE   lbc_nfd_sp, lbc_nfd_ext_sp
      MODULE PROCEDURE   lbc_nfd_dp, lbc_nfd_ext_dp
   END INTERFACE

   INTERFACE mpp_nfd            ! called by lbc_lnk_pt2pt or lbc_lnk_neicoll
      MODULE PROCEDURE   mpp_nfd_sp, mpp_nfd_dp
   END INTERFACE
   
   PUBLIC   mpp_nfd            ! mpi north fold conditions
   PUBLIC   lbc_nfd            ! north fold conditions

   INTEGER, PUBLIC                               :: nfd_nbnei
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION (:  ) :: nfd_rknei
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION (:,:) :: nfd_rksnd
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION (:,:) :: nfd_jisnd
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lbcnfd.F90 15267 2021-09-17 09:04:34Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!                   ***  routine lbc_nfd_[sd]p  ***
   !!                   ***  routine lbc_nfd_ext_[sd]p  ***
   !!----------------------------------------------------------------------
   !!
   !! ** Purpose :   lateral boundary condition 
   !!                North fold treatment without processor exchanges. 
   !!
   !! ** Method  :   
   !!
   !! ** Action  :   ptab with updated values along the north fold
   !!----------------------------------------------------------------------
   !
   !                       !==  SINGLE PRECISION VERSIONS
   !

   SUBROUTINE lbc_nfd_sp( ptab, cd_nat, psgn, khls, kfld )
      TYPE(PTR_4d_sp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(sp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      !
      INTEGER  ::    ji,  jj,  jk,  jl,  jf   ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipk, ipl, ipf   ! dimension of the input array
      INTEGER  ::   ii1, ii2, ij1, ij2
      !!----------------------------------------------------------------------
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( ipi /= Ni0glo+2*khls ) THEN
         WRITE(ctmp1,*) 'lbc_nfd input array does not match khls', ipi, khls, Ni0glo
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      DO jf = 1, ipf                      ! Loop on the number of arrays to be treated
         !
         IF( c_NFtype == 'T' ) THEN            ! *  North fold  T-point pivot
            !
            SELECT CASE ( cd_nat(jf) )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls+1
                        ii1 = khls + ji
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =   2 + khls + ji - 1        ! ends at: 2 + khls + ipi - 2*khls - 1 - 1 = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls + 1
                        ii1 = ipi - khls + ji
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! last khls-1 points
                        ii1 = ipi - khls + 1 + ji        ! ends at: ipi - khls + 1 + khls - 1 = ipi
                        ii2 = ipi - khls + 1 - ji        ! ends at: ipi - khls + 1 - khls + 1 = ipi - 2*khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2-1        ! points from ipi/2+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji + 1             ! ends at: ipi/2 + (ipi/2 - khls - 1) + 1 = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls - 1) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls-1 points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'U' )                               ! U-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2          ! points from ipi/2+1 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls-1 points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'V' )                               ! V-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls+1 lines (from ipj to ipj-khls) : full
               	  DO jj = 1, khls+1
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - ( khls + 1 ) + 1 = ipj - khls
                     ij2 = ipj - 2*khls + jj - 2         ! ends at: ipj - 2*khls + khls + 1 - 2 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls+1
                        ii1 = khls + ji
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =   2 + khls + ji - 1        ! ends at: 2 + khls + ipi - 2*khls - 1 - 1 = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls + 1
                        ii1 = ipi - khls + ji
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! last khls-1 points
                        ii1 = ipi - khls + 1 + ji        ! ends at: ipi - khls + 1 + khls - 1 = ipi
                        ii2 = ipi - khls + 1 - ji        ! ends at: ipi - khls + 1 - khls + 1 = ipi - 2*khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'F' )                               ! F-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls+1 lines (from ipj to ipj-khls) : full
               	  DO jj = 1, khls+1
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - ( khls + 1 ) + 1 = ipj - khls
                     ij2 = ipj - 2*khls + jj - 2         ! ends at: ipj - 2*khls + khls + 1 - 2 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            END SELECT   ! cd_nat(jf)
            !
         ENDIF   ! c_NFtype == 'T'
         !
         IF( c_NFtype == 'F' ) THEN            ! *  North fold  F-point pivot
            !
            SELECT CASE ( cd_nat(jf) )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! first: line number ipj-khls : 3 points
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, 1                 ! points from ipi/2+1
                        ii1 = ipi/2 + ji
                        ii2 = ipi/2 - ji + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                     DO ji = 1, 1                 ! points ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                     DO ji = 1, 1                 ! point khls: redo it just in case (if e-w periodocity already done)
                        !                         ! as we just changed point ipi - khls
                        ii1 = khls + ji - 1
                        ii2 = khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                  END DO
                  !
                  ! Second: last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj + 1      - jj             ! ends at: ipj + 1 - khls
                     ij2 = ipj - 2*khls + jj             ! ends at: ipj - 2*khls + khls = ipj - khls
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'U' )                               ! U-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj + 1      - jj             ! ends at: ipj + 1 - khls
                     ij2 = ipj - 2*khls + jj             ! ends at: ipj - 2*khls + khls = ipj - khls
                     !
                     DO ji = 1, khls-1            ! first khls-1 points
                        ii1 =          ji                ! ends at: khls-1
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls
                        ii1 = khls + ji - 1
                        ii2 = ipi - ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+1 to ipi - khls - 1  (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ( ipi - 2*khls - 1 ) = ipi - khls - 1
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - khls = ipi - 2*khls
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'V' )                               ! V-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji          ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1      ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls            ! last khls points
                        ii1 = ipi - khls + ji          ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1      ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO   
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2          ! points from ipi/2+1 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'F' )                               ! F-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls-1            ! first khls-1 points
                        ii1 =          ji                ! ends at: khls-1
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls
                        ii1 = khls + ji - 1
                        ii2 = ipi - ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+1 to ipi - khls - 1  (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ( ipi - 2*khls - 1 ) = ipi - khls - 1
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - khls = ipi - 2*khls
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO   
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2-1        ! points from ipi/2+1 to ipi - khls-1  (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji                 ! ends at: ipi/2 - (ipi/2 - khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! first khls-1 points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-nn_hl-1  
                        ii1 =          ji                ! ends at: khls
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            END SELECT   ! cd_nat(jf)
            !
         ENDIF   ! c_NFtype == 'F'
         !
      END DO   ! ipf
      !
   END SUBROUTINE lbc_nfd_sp


   SUBROUTINE lbc_nfd_ext_sp( ptab, cd_nat, psgn, kextj )
      !!----------------------------------------------------------------------
      REAL(sp), DIMENSION(:,1-kextj:),INTENT(inout) ::   ptab
      CHARACTER(len=1), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(sp),  INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER,          INTENT(in   ) ::   kextj       ! extra halo width at north fold
      !
      INTEGER  ::    ji,  jj,  jh   ! dummy loop indices
      INTEGER  ::   ipj
      INTEGER  ::   ijt, iju, ipjm1
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( jpni )
      CASE ( 1 )     ;   ipj = jpj        ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   ipj = 4          ! several proc along the i-direction
      END SELECT
      !
      ipjm1 = ipj-1
      !
      IF( c_NFtype == 'T' ) THEN            ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_nat  )
         CASE ( 'T' , 'W' )                         ! T-, W-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-2-jh)
               END DO
               ptab(1,ipj+jh) = psgn * ptab(3,ipj-2-jh)
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+2
               ptab(ji,ipjm1) = psgn * ptab(ijt,ipjm1)
            END DO
         CASE ( 'U' )                               ! U-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo-1
                  iju = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(iju,ipj-2-jh)
               END DO
               ptab(   1  ,ipj+jh) = psgn * ptab(    2   ,ipj-2-jh)
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-1,ipj-2-jh) 
            END DO
            DO ji = jpiglo/2, jpiglo-1
               iju = jpiglo-ji+1
               ptab(ji,ipjm1) = psgn * ptab(iju,ipjm1)
            END DO
         CASE ( 'V' )                               ! V-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  ptab(ji,ipj-1+jh) = psgn * ptab(ijt,ipj-2-jh)
                  ptab(ji,ipj+jh  ) = psgn * ptab(ijt,ipj-3-jh)
               END DO
               ptab(1,ipj+jh) = psgn * ptab(3,ipj-3-jh) 
            END DO
         CASE ( 'F' )                               ! F-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  ptab(ji,ipj-1+jh) = psgn * ptab(iju,ipj-2-jh)
                  ptab(ji,ipj+jh  ) = psgn * ptab(iju,ipj-3-jh)
               END DO
            END DO
            DO jh = 0, kextj
               ptab(   1  ,ipj+jh) = psgn * ptab(    2   ,ipj-3-jh)
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-1,ipj-3-jh)
            END DO
         END SELECT
         !
      ENDIF   ! c_NFtype == 'T'
      !
      IF( c_NFtype == 'F' ) THEN            ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_nat  )
         CASE ( 'T' , 'W' )                         ! T-, W-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-1-jh)
               END DO
            END DO
         CASE ( 'U' )                               ! U-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  ptab(ji,ipj+jh) = psgn * ptab(iju,ipj-1-jh)
               END DO
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-2,ipj-1-jh)
            END DO
         CASE ( 'V' )                               ! V-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-2-jh)
               END DO
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+1
               ptab(ji,ipjm1) = psgn * ptab(ijt,ipjm1)
            END DO
         CASE ( 'F' )                               ! F-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  ptab(ji,ipj+jh  ) = psgn * ptab(iju,ipj-2-jh)
               END DO
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-2,ipj-2-jh)
            END DO
            DO ji = jpiglo/2+1, jpiglo-1
               iju = jpiglo-ji
               ptab(ji,ipjm1) = psgn * ptab(iju,ipjm1)
            END DO
         END SELECT
         !
      ENDIF   ! c_NFtype == 'F'
      !
   END SUBROUTINE lbc_nfd_ext_sp

   !
   !                       !==  DOUBLE PRECISION VERSIONS
   !

   SUBROUTINE lbc_nfd_dp( ptab, cd_nat, psgn, khls, kfld )
      TYPE(PTR_4d_dp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(dp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      !
      INTEGER  ::    ji,  jj,  jk,  jl,  jf   ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipk, ipl, ipf   ! dimension of the input array
      INTEGER  ::   ii1, ii2, ij1, ij2
      !!----------------------------------------------------------------------
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( ipi /= Ni0glo+2*khls ) THEN
         WRITE(ctmp1,*) 'lbc_nfd input array does not match khls', ipi, khls, Ni0glo
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      DO jf = 1, ipf                      ! Loop on the number of arrays to be treated
         !
         IF( c_NFtype == 'T' ) THEN            ! *  North fold  T-point pivot
            !
            SELECT CASE ( cd_nat(jf) )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls+1
                        ii1 = khls + ji
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =   2 + khls + ji - 1        ! ends at: 2 + khls + ipi - 2*khls - 1 - 1 = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls + 1
                        ii1 = ipi - khls + ji
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! last khls-1 points
                        ii1 = ipi - khls + 1 + ji        ! ends at: ipi - khls + 1 + khls - 1 = ipi
                        ii2 = ipi - khls + 1 - ji        ! ends at: ipi - khls + 1 - khls + 1 = ipi - 2*khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2-1        ! points from ipi/2+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji + 1             ! ends at: ipi/2 + (ipi/2 - khls - 1) + 1 = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls - 1) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls-1 points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'U' )                               ! U-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2          ! points from ipi/2+1 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls-1 points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'V' )                               ! V-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls+1 lines (from ipj to ipj-khls) : full
               	  DO jj = 1, khls+1
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - ( khls + 1 ) + 1 = ipj - khls
                     ij2 = ipj - 2*khls + jj - 2         ! ends at: ipj - 2*khls + khls + 1 - 2 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 2 - ji            ! ends at: 2*khls + 2 - khls = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls+1
                        ii1 = khls + ji
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+2 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =   2 + khls + ji - 1        ! ends at: 2 + khls + ipi - 2*khls - 1 - 1 = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) + 1 = khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls + 1
                        ii1 = ipi - khls + ji
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! last khls-1 points
                        ii1 = ipi - khls + 1 + ji        ! ends at: ipi - khls + 1 + khls - 1 = ipi
                        ii2 = ipi - khls + 1 - ji        ! ends at: ipi - khls + 1 - khls + 1 = ipi - 2*khls + 2
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'F' )                               ! F-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls+1 lines (from ipj to ipj-khls) : full
               	  DO jj = 1, khls+1
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - ( khls + 1 ) + 1 = ipj - khls
                     ij2 = ipj - 2*khls + jj - 2         ! ends at: ipj - 2*khls + khls + 1 - 2 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            END SELECT   ! cd_nat(jf)
            !
         ENDIF   ! c_NFtype == 'T'
         !
         IF( c_NFtype == 'F' ) THEN            ! *  North fold  F-point pivot
            !
            SELECT CASE ( cd_nat(jf) )
            CASE ( 'T' , 'W' )                         ! T-, W-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! first: line number ipj-khls : 3 points
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, 1                 ! points from ipi/2+1
                        ii1 = ipi/2 + ji
                        ii2 = ipi/2 - ji + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                     DO ji = 1, 1                 ! points ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 =       khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                     DO ji = 1, 1                 ! point khls: redo it just in case (if e-w periodocity already done)
                        !                         ! as we just changed point ipi - khls
                        ii1 = khls + ji - 1
                        ii2 = khls + ji
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) =            ptab(jf)%pt4d(ii2,ij2,jk,jl)   ! Warning: pb with sign...
                     END DO
                  END DO
                  !
                  ! Second: last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj + 1      - jj             ! ends at: ipj + 1 - khls
                     ij2 = ipj - 2*khls + jj             ! ends at: ipj - 2*khls + khls = ipj - khls
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1        ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'U' )                               ! U-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj + 1      - jj             ! ends at: ipj + 1 - khls
                     ij2 = ipj - 2*khls + jj             ! ends at: ipj - 2*khls + khls = ipj - khls
                     !
                     DO ji = 1, khls-1            ! first khls-1 points
                        ii1 =          ji                ! ends at: khls-1
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls
                        ii1 = khls + ji - 1
                        ii2 = ipi - ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+1 to ipi - khls - 1  (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ( ipi - 2*khls - 1 ) = ipi - khls - 1
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - khls = ipi - 2*khls
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO
                  !
               END DO; END DO
            CASE ( 'V' )                               ! V-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls              ! first khls points
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo            ! points from khls to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji          ! ends at: khls + ipi - 2*khls = ipi - khls
                        ii2 = ipi - khls - ji + 1      ! ends at: ipi - khls - ( ipi - 2*khls ) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls            ! last khls points
                        ii1 = ipi - khls + ji          ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji + 1      ! ends at: ipi - khls + 1 - khls = ipi - 2*khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO   
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2          ! points from ipi/2+1 to ipi - khls   (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji + 1             ! ends at: ipi/2 - (ipi/2 - khls) + 1 = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! first khls points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-khls  
                        ii1 =              ji            ! ends at: khls
                        ii2 = 2*khls + 1 - ji            ! ends at: 2*khls + 1 - khls = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            CASE ( 'F' )                               ! F-point
               DO jl = 1, ipl   ;   DO jk = 1, ipk
                  !
                  ! last khls lines (from ipj to ipj-khls+1) : full
               	  DO jj = 1, khls
               	     ij1 = ipj          - jj + 1         ! ends at: ipj - khls + 1
                     ij2 = ipj - 2*khls + jj - 1         ! ends at: ipj - 2*khls + khls - 1 = ipj - khls - 1
                     !
                     DO ji = 1, khls-1            ! first khls-1 points
                        ii1 =          ji                ! ends at: khls-1
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point khls
                        ii1 = khls + ji - 1
                        ii2 = ipi - ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, Ni0glo - 1        ! points from khls+1 to ipi - khls - 1  (note: Ni0glo = ipi - 2*khls)
                        ii1 =       khls + ji            ! ends at: khls + ( ipi - 2*khls - 1 ) = ipi - khls - 1
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - ( ipi - 2*khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, 1                 ! point ipi - khls
                        ii1 = ipi - khls + ji - 1
                        ii2 = ii1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls              ! last khls points
                        ii1 = ipi - khls + ji            ! ends at: ipi - khls + khls = ipi
                        ii2 = ipi - khls - ji            ! ends at: ipi - khls - khls = ipi - 2*khls
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                  END DO   
                  !
                  ! line number ipj-khls : right half
               	  DO jj = 1, 1
                     ij1 = ipj - khls
                     ij2 = ij1   ! same line
                     !
                     DO ji = 1, Ni0glo/2-1        ! points from ipi/2+1 to ipi - khls-1  (note: Ni0glo = ipi - 2*khls)
                        ii1 = ipi/2 + ji                 ! ends at: ipi/2 + (ipi/2 - khls) = ipi - khls
                        ii2 = ipi/2 - ji                 ! ends at: ipi/2 - (ipi/2 - khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     DO ji = 1, khls-1            ! first khls-1 points: redo them just in case (if e-w periodocity already done)
                        !                         ! as we just changed points ipi-2khls+1 to ipi-nn_hl-1  
                        ii1 =          ji                ! ends at: khls
                        ii2 = 2*khls - ji                ! ends at: 2*khls - ( khls - 1 ) = khls + 1
                        ptab(jf)%pt4d(ii1,ij1,jk,jl) = psgn(jf) * ptab(jf)%pt4d(ii2,ij2,jk,jl)
                     END DO
                     !                            ! last khls points: have been / will done by e-w periodicity 
                  END DO
                  !
               END DO; END DO
            END SELECT   ! cd_nat(jf)
            !
         ENDIF   ! c_NFtype == 'F'
         !
      END DO   ! ipf
      !
   END SUBROUTINE lbc_nfd_dp


   SUBROUTINE lbc_nfd_ext_dp( ptab, cd_nat, psgn, kextj )
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,1-kextj:),INTENT(inout) ::   ptab
      CHARACTER(len=1), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(dp),  INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER,          INTENT(in   ) ::   kextj       ! extra halo width at north fold
      !
      INTEGER  ::    ji,  jj,  jh   ! dummy loop indices
      INTEGER  ::   ipj
      INTEGER  ::   ijt, iju, ipjm1
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( jpni )
      CASE ( 1 )     ;   ipj = jpj        ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   ipj = 4          ! several proc along the i-direction
      END SELECT
      !
      ipjm1 = ipj-1
      !
      IF( c_NFtype == 'T' ) THEN            ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_nat  )
         CASE ( 'T' , 'W' )                         ! T-, W-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-2-jh)
               END DO
               ptab(1,ipj+jh) = psgn * ptab(3,ipj-2-jh)
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+2
               ptab(ji,ipjm1) = psgn * ptab(ijt,ipjm1)
            END DO
         CASE ( 'U' )                               ! U-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo-1
                  iju = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(iju,ipj-2-jh)
               END DO
               ptab(   1  ,ipj+jh) = psgn * ptab(    2   ,ipj-2-jh)
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-1,ipj-2-jh) 
            END DO
            DO ji = jpiglo/2, jpiglo-1
               iju = jpiglo-ji+1
               ptab(ji,ipjm1) = psgn * ptab(iju,ipjm1)
            END DO
         CASE ( 'V' )                               ! V-point
            DO jh = 0, kextj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  ptab(ji,ipj-1+jh) = psgn * ptab(ijt,ipj-2-jh)
                  ptab(ji,ipj+jh  ) = psgn * ptab(ijt,ipj-3-jh)
               END DO
               ptab(1,ipj+jh) = psgn * ptab(3,ipj-3-jh) 
            END DO
         CASE ( 'F' )                               ! F-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  ptab(ji,ipj-1+jh) = psgn * ptab(iju,ipj-2-jh)
                  ptab(ji,ipj+jh  ) = psgn * ptab(iju,ipj-3-jh)
               END DO
            END DO
            DO jh = 0, kextj
               ptab(   1  ,ipj+jh) = psgn * ptab(    2   ,ipj-3-jh)
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-1,ipj-3-jh)
            END DO
         END SELECT
         !
      ENDIF   ! c_NFtype == 'T'
      !
      IF( c_NFtype == 'F' ) THEN            ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_nat  )
         CASE ( 'T' , 'W' )                         ! T-, W-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-1-jh)
               END DO
            END DO
         CASE ( 'U' )                               ! U-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  ptab(ji,ipj+jh) = psgn * ptab(iju,ipj-1-jh)
               END DO
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-2,ipj-1-jh)
            END DO
         CASE ( 'V' )                               ! V-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  ptab(ji,ipj+jh) = psgn * ptab(ijt,ipj-2-jh)
               END DO
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+1
               ptab(ji,ipjm1) = psgn * ptab(ijt,ipjm1)
            END DO
         CASE ( 'F' )                               ! F-point
            DO jh = 0, kextj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  ptab(ji,ipj+jh  ) = psgn * ptab(iju,ipj-2-jh)
               END DO
               ptab(jpiglo,ipj+jh) = psgn * ptab(jpiglo-2,ipj-2-jh)
            END DO
            DO ji = jpiglo/2+1, jpiglo-1
               iju = jpiglo-ji
               ptab(ji,ipjm1) = psgn * ptab(iju,ipjm1)
            END DO
         END SELECT
         !
      ENDIF   ! c_NFtype == 'F'
      !
   END SUBROUTINE lbc_nfd_ext_dp


   !!======================================================================
   !
   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_nfd_[sd]p  ***
   !!
   !!   * Argument : dummy argument use in mpp_nfd_... routines
   !!                ptab      :   pointer of arrays on which the boundary condition is applied
   !!                cd_nat    :   nature of array grid-points
   !!                psgn      :   sign used across the north fold boundary
   !!                kfld      :   optional, number of pt3d arrays
   !!                kfillmode :   optional, method to be use to fill the halos (see jpfill* variables)
   !!                pfillval  :   optional, background value (used with jpfillcopy)
   !!----------------------------------------------------------------------
   !!
   !!   ----   SINGLE PRECISION VERSIONS
   !!

   SUBROUTINE mpp_nfd_sp( ptab, cd_nat, psgn, kfillmode, pfillval, khls, kfld, fTag )
      TYPE(PTR_4d_sp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(sp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfillmode   ! filling method for halo over land 
      REAL(sp)               , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER                       , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER, OPTIONAL, INTENT(in)                 ::   fTag        ! if present, there may be multithreaded messaging
      !
      LOGICAL  ::   ll_add_line
      INTEGER  ::   ji,  jj,  jk,  jl, jf, jr, jg, jn   ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipj2, ipk, ipl, ipf   ! dimension of the input array
      INTEGER  ::   ierr, ibuffsize, iis0, iie0, impp
      INTEGER  ::   ii1, ii2, ij1, ij2, iis, iie, iib, iig, iin
      INTEGER  ::   i0max
      INTEGER  ::   ij, iproc, ipni, ijnr
      INTEGER, DIMENSION (:), ALLOCATABLE ::   ireq_s, ireq_r   ! for mpi_isend when avoiding mpi_allgather
      INTEGER                             ::   ipjtot           ! sum of lines for all multi fields
      INTEGER                             ::   i012             ! 0, 1 or 2
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijsnd  ! j-position of sent lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijbuf  ! j-position of send buffer lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijrcv  ! j-position of recv buffer lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ii1st, iiend
      INTEGER , DIMENSION(:)          , ALLOCATABLE ::   ipjfld ! number of sent lines for each field
      REAL(sp), DIMENSION(:,:,:,:)    , ALLOCATABLE ::   zbufs  ! buffer, receive and work arrays
      REAL(sp), DIMENSION(:,:,:,:,:)  , ALLOCATABLE ::   zbufr  ! buffer, receive and work arrays
      REAL(sp), DIMENSION(:,:,:,:,:)  , ALLOCATABLE ::   znorthloc
      REAL(sp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE ::   znorthglo
      TYPE(PTR_4D_sp), DIMENSION(:), ALLOCATABLE ::   ztabglo        ! array or pointer of arrays on which apply the b.c.
      !!----------------------------------------------------------------------
      !
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( ln_nnogather ) THEN      !==  no allgather exchanges  ==!

         !   ---   define number of exchanged lines   ---
         !
         ! In theory we should exchange only nn_hls lines.
         !
         ! However, some other points are duplicated in the north pole folding:
         !  - c_NFtype='T', grid=T : half of the last line (jpiglo/2+2:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=U : half of the last line (jpiglo/2+1:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=V : all the last line nn_hls+1 and (nn_hls+2:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=F : all the last line (nn_hls+1:jpiglo-nn_hls)
         !  - c_NFtype='F', grid=T : 2 points of the last line (jpiglo/2+1 and jpglo-nn_hls)
         !  - c_NFtype='F', grid=U : no points are duplicated
         !  - c_NFtype='F', grid=V : half of the last line (jpiglo/2+1:jpiglo-nn_hls)
         !  - c_NFtype='F', grid=F : half of the last line (jpiglo/2+1:jpiglo-nn_hls-1)
         ! The order of the calculations may differ for these duplicated points (as, for example jj+1 becomes jj-1)
         ! This explain why these duplicated points may have different values even if they are at the exact same location.
         ! In consequence, we may want to force the folding on these points by setting l_full_nf_update = .TRUE.
         ! This is slightly slower but necessary to avoid different values on identical grid points!!
         !
         !!!!!!!!!           temporary switch off this optimisation ==> force TRUE           !!!!!!!!
         !!!!!!!!!  needed to get the same results without agrif and with agrif and no zoom  !!!!!!!!
         !!!!!!!!!                    I don't know why we must do that...                    !!!!!!!!
         l_full_nf_update = .TRUE.
         ! also force it if not restart during the first 2 steps (leap frog?)
         ll_add_line = l_full_nf_update .OR. ( ncom_stp <= nit000+1 .AND. .NOT. ln_rstart )
         
         ALLOCATE(ipjfld(ipf))                 ! how many lines do we exchange for each field?
         IF( ll_add_line ) THEN
            DO jf = 1, ipf                     ! Loop over the number of arrays to be processed
               ipjfld(jf) = khls + COUNT( (/ c_NFtype == 'T' .OR. cd_nat(jf) == 'V' .OR. cd_nat(jf) == 'F' /) )
            END DO
         ELSE
            ipjfld(:) = khls
         ENDIF
         
         ipj    = MAXVAL(ipjfld(:))            ! Max 2nd dimension of message transfers
         ipjtot = SUM(   ipjfld(:))            ! Total number of lines to be exchanged

         ! Index of modifying lines in input
         ALLOCATE( ijsnd(ipj, ipf), ijbuf(ipj, ipf), ijrcv(ipj, ipf), ii1st(ipj, ipf), iiend(ipj, ipf) )

         ij1 = 0
         DO jf = 1, ipf                        ! Loop over the number of arrays to be processed
            !
            DO jj = 1, khls   ! first khls lines (starting from top) must be fully defined
               ii1st(jj, jf) = 1
               iiend(jj, jf) = jpi
            END DO
            !
            ! what do we do with line khls+1 (starting from top)
            IF( c_NFtype == 'T' ) THEN          ! *  North fold  T-point pivot
               SELECT CASE ( cd_nat(jf) )
               CASE ('T','W')   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+2)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls)
               CASE ('U'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls)
               CASE ('V'    )   ;   i012 = 2   ;   ii1st(khls+1, jf) = 1                 ;   iiend(khls+1, jf) = jpi
               CASE ('F'    )   ;   i012 = 2   ;   ii1st(khls+1, jf) = 1                 ;   iiend(khls+1, jf) = jpi
               END SELECT
            ENDIF
            IF( c_NFtype == 'F' ) THEN          ! *  North fold  F-point pivot
               SELECT CASE ( cd_nat(jf) )
               CASE ('T','W')   ;   i012 = 0   ! we don't touch line khls+1
               CASE ('U'    )   ;   i012 = 0   ! we don't touch line khls+1
               CASE ('V'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls  )
               CASE ('F'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls-1)
               END SELECT
            ENDIF
            !
            DO jj = 1, ipjfld(jf)
               ij1 = ij1 + 1
               ijsnd(jj,jf) = jpj - 2*khls + jj - i012   ! sent lines (from bottom of sent lines)
               ijbuf(jj,jf) = ij1                        ! gather all lines in the snd/rcv buffers
               ijrcv(jj,jf) = jpj - jj + 1               ! recv lines (from the top -> reverse order for jj)
            END DO
            !
         END DO
         !
         i0max = jpimax - 2 * khls                                    ! we are not sending the halos
         ALLOCATE( zbufs(i0max,ipjtot,ipk,ipl), ireq_s(nfd_nbnei) )   ! store all the data to be sent in a buffer array
         ibuffsize = i0max * ipjtot * ipk * ipl
         !
         ! fill the send buffer with all the lines
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
            DO jj = 1, ipjfld(jf)
               ij1 = ijbuf(jj,jf)
               ij2 = ijsnd(jj,jf)
               DO ji = Nis0, Nie0       ! should not use any other value
                  iib = ji - Nis0 + 1
                  zbufs(iib,ij1,jk,jl) = ptab(jf)%pt4d(ji,ij2,jk,jl)
               END DO
               DO ji = Ni_0+1, i0max    ! avoid sending uninitialized values (make sure we don't use it)
                  zbufs(ji,ij1,jk,jl) = HUGE(0._sp)   ! make sure we don't use it...
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         ! start waiting time measurement
         IF( ln_timing ) CALL tic_tac(.TRUE.)
         !
         ! send the same buffer data to all neighbourgs as soon as possible
         DO jn = 1, nfd_nbnei
            iproc = nfd_rknei(jn)
            IF( iproc /= narea-1 .AND. iproc /= -1 ) THEN
               CALL MPI_Isend( zbufs, ibuffsize, MPI_REAL, iproc, 5, mpi_comm_oce, ireq_s(jn), ierr )
            ELSE
               ireq_s(jn) = MPI_REQUEST_NULL
            ENDIF
         END DO
         !
         ALLOCATE( zbufr(i0max,ipjtot,ipk,ipl,nfd_nbnei), ireq_r(nfd_nbnei) ) 
         !
         DO jn = 1, nfd_nbnei
            !
            iproc = nfd_rknei(jn)
            !
            IF(           iproc == -1 ) THEN   ! No neighbour (land proc that was suppressed)
               !
               ireq_r(jn) = MPI_REQUEST_NULL                ! no message to be received
               zbufr(:,:,:,:,jn) = HUGE(0._sp)   ! default: define it and make sure we don't use it...
               SELECT CASE ( kfillmode )
               CASE ( jpfillnothing )                       ! no filling 
               CASE ( jpfillcopy    )                       ! filling with inner domain values
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipjfld(jf)
                        ij1 = ijbuf(jj,jf)
                        ij2 = ijsnd(jj,jf)                                      ! we will use only the first value, see init_nfdcom
                        zbufr(1,ij1,jk,jl,jn) = ptab(jf)%pt4d(Nis0,ij2,jk,jl)   ! chose to take the 1st inner domain point
                     END DO
                  END DO   ;   END DO   ;   END DO
               CASE ( jpfillcst     )                       ! filling with constant value
                  zbufr(1,:,:,:,jn) = pfillval              ! we will use only the first value, see init_nfdcom
               END SELECT
               !
            ELSE IF( iproc == narea-1 ) THEN   ! get data from myself!
               !
               ireq_r(jn) = MPI_REQUEST_NULL                ! no message to be received
               DO jf = 1, ipf   ;   DO jl = 1, ipl  ;   DO jk = 1, ipk
                  DO jj = 1, ipjfld(jf)
                     ij1 = ijbuf(jj,jf)
                     ij2 = ijsnd(jj,jf)
                     DO ji = Nis0, Nie0                     ! should not use any other value
                        iib = ji - Nis0 + 1
                        zbufr(iib,ij1,jk,jl,jn) = ptab(jf)%pt4d(ji,ij2,jk,jl)
                     END DO
                  END DO
               END DO   ;   END DO   ;   END DO
               !
            ELSE                               ! get data from a neighbour trough communication
               CALL MPI_Irecv( zbufr(:,:,:,:,jn), ibuffsize, MPI_REAL, iproc, 5, mpi_comm_oce, ireq_r(jn), ierr )
            ENDIF
            !
         END DO   ! nfd_nbnei
         !
         CALL mpi_waitall(nfd_nbnei, ireq_r, MPI_STATUSES_IGNORE, ierr)   ! wait for all Irecv
         !
         IF( ln_timing ) CALL tic_tac(.FALSE.)
         !
         ! North fold boundary condition
         !
         DO jf = 1, ipf
            !
            SELECT CASE ( cd_nat(jf) )     ! which grid number?
            CASE ('T','W')   ;   iig = 1   ! T-, W-point
            CASE ('U')       ;   iig = 2   ! U-point
            CASE ('V')       ;   iig = 3   ! V-point
            CASE ('F')       ;   iig = 4   ! F-point
            END SELECT
            !
            DO jl = 1, ipl   ;   DO jk = 1, ipk
               !
               ! if T point with F-point pivot : must be done first
               !    --> specific correction of 3 points near the 2 pivots (to be clean, usually masked -> so useless) 
               IF( c_NFtype == 'F' .AND. iig == 1 ) THEN
                  ij1 = jpj - khls     ! j-index in the receiving array
                  ij2 = 1              ! only 1 line in the buffer
                  DO ji = mi0(khls), mi1(khls)               ! change because of EW periodicity as we also change jpiglo-khls
                     iib = nfd_jisnd(mi0(       khls),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0(       khls),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
                  DO ji = mi0(jpiglo/2+1), mi1(jpiglo/2+1)
                     iib = nfd_jisnd(mi0( jpiglo/2+1),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0( jpiglo/2+1),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
                  DO ji = mi0(jpiglo-khls), mi1(jpiglo-khls)
                     iib = nfd_jisnd(mi0(jpiglo-khls),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0(jpiglo-khls),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
               ENDIF
               !
               ! Apply the North pole folding.
               DO jj = 1, ipjfld(jf)   ! for all lines to be exchanged for this field
                  ij1 = ijrcv(jj,jf)   ! j-index in the receiving array
                  ij2 = ijbuf(jj,jf)   ! j-index in the buffer
                  iis = ii1st(jj,jf)   ! stating i-index in the receiving array
                  iie = iiend(jj,jf)   !  ending i-index in the receiving array
                  DO ji = iis, iie 
                     iib = nfd_jisnd(ji,iig)   ! i-index in the buffer
                     iin = nfd_rksnd(ji,iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = psgn(jf) * zbufr(iib,ij2,jk,jl,iin)
                  END DO
               END DO
               !
               ! re-apply periodocity when we modified the eastern side of the inner domain (and not the full line)
               IF( c_NFtype == 'T' ) THEN          ! *  North fold  T-point pivot
                  IF(     iig <= 2 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls)   ! 'T','W','U': update west halo
                  ELSE                      ;   iis = 1        ;   iie = 0           ! 'V','F'    : full line already exchanged
                  ENDIF
               ENDIF
               IF( c_NFtype == 'F' ) THEN          ! *  North fold  F-point pivot
                  IF(     iig <= 2 ) THEN   ;   iis = 1        ;   iie = 0           ! 'T','W','U': nothing to do
                  ELSEIF( iig == 3 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls)   ! 'V'        : update west halo
                  ELSEIF( khls > 1 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls-1) ! 'F' and khls > 1
                  ELSE                      ;   iis = 1        ;   iie = 0           ! 'F' and khls == 1 : nothing to do
                  ENDIF
               ENDIF
               jj  = ipjfld(jf)     ! only for the last line of this field
               ij1 = ijrcv(jj,jf)   ! j-index in the receiving array
               ij2 = ijbuf(jj,jf)   ! j-index in the buffer
               DO ji = iis, iie
                  iib = nfd_jisnd(ji,iig)   ! i-index in the buffer
                  iin = nfd_rksnd(ji,iig)   ! neigbhour-index in the buffer
                  IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                  ptab(jf)%pt4d(ji,ij1,jk,jl) = psgn(jf) * zbufr(iib,ij2,jk,jl,iin)
               END DO
               !               
            END DO   ;   END DO   ! ipl   ; ipk
            !               
         END DO   ! ipf
       
         !
         DEALLOCATE( zbufr, ireq_r, ijsnd, ijbuf, ijrcv, ii1st, iiend, ipjfld )
         !
         CALL mpi_waitall(nfd_nbnei, ireq_s, MPI_STATUSES_IGNORE, ierr)   ! wait for all Isend
         !
         DEALLOCATE( zbufs, ireq_s )
         !
      ELSE                             !==  allgather exchanges  ==!
         !
         ! how many lines do we exchange at max? -> ipj    (no further optimizations in this case...)
         ipj =      khls + 2
         ! how many lines do we     need at max? -> ipj2   (no further optimizations in this case...)
         ipj2 = 2 * khls + 2
         !
         i0max = jpimax - 2 * khls
         ibuffsize = i0max * ipj * ipk * ipl * ipf
         ALLOCATE( znorthloc(i0max,ipj,ipk,ipl,ipf), znorthglo(i0max,ipj,ipk,ipl,ipf,ndim_rank_north) )
         !
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk               ! put in znorthloc ipj j-lines of ptab
            DO jj = 1, ipj
               ij2 = jpj - ipj2 + jj                        ! the first ipj lines of the last ipj2 lines
               DO ji = 1, Ni_0
                  ii2 = Nis0 - 1 + ji                       ! inner domain: Nis0 to Nie0
                  znorthloc(ji,jj,jk,jl,jf) = ptab(jf)%pt4d(ii2,ij2,jk,jl)
               END DO
               DO ji = Ni_0+1, i0max
                  znorthloc(ji,jj,jk,jl,jf) = HUGE(0._sp)   ! avoid sending uninitialized values (make sure we don't use it)
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         ! start waiting time measurement
         IF( ln_timing ) CALL tic_tac(.TRUE.)
         CALL MPI_ALLGATHER( znorthloc, ibuffsize, MPI_REAL, znorthglo, ibuffsize, MPI_REAL, ncomm_north, ierr )
         ! stop waiting time measurement
         IF( ln_timing ) CALL tic_tac(.FALSE.)
         DEALLOCATE( znorthloc )
         ALLOCATE( ztabglo(ipf) )
         DO jf = 1, ipf
            ALLOCATE( ztabglo(jf)%pt4d(jpiglo,ipj2,ipk,ipl) )
         END DO
         !
         ! need to fill only the first ipj lines of ztabglo as lbc_nfd don't use the last khls lines
         ijnr = 0
         DO jr = 1, jpni                                                        ! recover the global north array
            iproc = nfproc(jr)
            impp  = nfimpp(jr)
            ipi   = nfjpi( jr) - 2 * khls                       ! corresponds to Ni_0 but for subdomain iproc
            IF( iproc == -1 ) THEN   ! No neighbour (land proc that was suppressed)
              !
               SELECT CASE ( kfillmode )
               CASE ( jpfillnothing )               ! no filling
                  CALL ctl_stop( 'STOP', 'mpp_nfd_generic : cannot use jpfillnothing with ln_nnogather = F')
               CASE ( jpfillcopy    )               ! filling with inner domain values
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipj
                        ij2 = jpj - ipj2 + jj                    ! the first ipj lines of the last ipj2 lines
                        DO ji = 1, ipi
                           ii1 = impp + khls + ji - 1            ! corresponds to mig(khls + ji) but for subdomain iproc
                           ztabglo(jf)%pt4d(ii1,jj,jk,jl) = ptab(jf)%pt4d(Nis0,ij2,jk,jl) ! chose to take the 1st inner domain point
                        END DO
                     END DO
                  END DO   ;   END DO   ;   END DO
               CASE ( jpfillcst     )               ! filling with constant value
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipj
                        DO ji = 1, ipi
                           ii1 = impp + khls + ji - 1            ! corresponds to mig(khls + ji) but for subdomain iproc
                           ztabglo(jf)%pt4d(ii1,jj,jk,jl) = pfillval
                        END DO
                     END DO
                 END DO   ;   END DO   ;   END DO
               END SELECT
               !
            ELSE
               ijnr = ijnr + 1
               DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                  DO jj = 1, ipj
                     DO ji = 1, ipi
                        ii1 = impp + khls + ji - 1               ! corresponds to mig(khls + ji) but for subdomain iproc
                        ztabglo(jf)%pt4d(ii1,jj,jk,jl) = znorthglo(ji,jj,jk,jl,jf,ijnr)
                     END DO
                  END DO
               END DO   ;   END DO   ;   END DO
            ENDIF
            !
         END DO   ! jpni
         DEALLOCATE( znorthglo )
         !
         DO jf = 1, ipf
            CALL lbc_nfd( ztabglo(jf:jf), cd_nat(jf:jf), psgn(jf:jf), khls, 1 )   ! North fold boundary condition
            DO jl = 1, ipl   ;   DO jk = 1, ipk                  ! e-w periodicity
               DO jj = 1, khls + 1
                  ij1 = ipj2 - (khls + 1) + jj                   ! need only the last khls + 1 lines until ipj2
                  ztabglo(jf)%pt4d(            1:  khls,ij1,jk,jl) = ztabglo(jf)%pt4d(jpiglo-2*khls+1:jpiglo-khls,ij1,jk,jl)
                  ztabglo(jf)%pt4d(jpiglo-khls+1:jpiglo,ij1,jk,jl) = ztabglo(jf)%pt4d(         khls+1:     2*khls,ij1,jk,jl)
               END DO
            END DO   ;   END DO
         END DO     
         !
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk               ! Scatter back to ARRAY_IN
            DO jj = 1, khls + 1
               ij1 = jpj  - (khls + 1) + jj   ! last khls + 1 lines until jpj
               ij2 = ipj2 - (khls + 1) + jj   ! last khls + 1 lines until ipj2
               DO ji= 1, jpi
                  ii2 = mig(ji)
                  ptab(jf)%pt4d(ji,ij1,jk,jl) = ztabglo(jf)%pt4d(ii2,ij2,jk,jl)
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         DO jf = 1, ipf
            DEALLOCATE( ztabglo(jf)%pt4d )
         END DO
         DEALLOCATE( ztabglo )
         !
      ENDIF   ! ln_nnogather
      !
   END SUBROUTINE mpp_nfd_sp

   !!
   !!   ----   DOUBLE PRECISION VERSIONS
   !!

   SUBROUTINE mpp_nfd_dp( ptab, cd_nat, psgn, kfillmode, pfillval, khls, kfld, fTag)
      TYPE(PTR_4d_dp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(dp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfillmode   ! filling method for halo over land 
      REAL(dp)               , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER                       , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER, OPTIONAL, INTENT(in)                 ::   fTag        ! if present, there may be multithreaded messaging
      !
      LOGICAL  ::   ll_add_line
      INTEGER  ::   ji,  jj,  jk,  jl, jf, jr, jg, jn   ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipj2, ipk, ipl, ipf   ! dimension of the input array
      INTEGER  ::   ierr, ibuffsize, iis0, iie0, impp
      INTEGER  ::   ii1, ii2, ij1, ij2, iis, iie, iib, iig, iin
      INTEGER  ::   i0max
      INTEGER  ::   ij, iproc, ipni, ijnr
      INTEGER, DIMENSION (:), ALLOCATABLE ::   ireq_s, ireq_r   ! for mpi_isend when avoiding mpi_allgather
      INTEGER                             ::   ipjtot           ! sum of lines for all multi fields
      INTEGER                             ::   i012             ! 0, 1 or 2
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijsnd  ! j-position of sent lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijbuf  ! j-position of send buffer lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ijrcv  ! j-position of recv buffer lines for each field
      INTEGER , DIMENSION(:,:)        , ALLOCATABLE ::   ii1st, iiend
      INTEGER , DIMENSION(:)          , ALLOCATABLE ::   ipjfld ! number of sent lines for each field
      REAL(dp), DIMENSION(:,:,:,:)    , ALLOCATABLE ::   zbufs  ! buffer, receive and work arrays
      REAL(dp), DIMENSION(:,:,:,:,:)  , ALLOCATABLE ::   zbufr  ! buffer, receive and work arrays
      REAL(dp), DIMENSION(:,:,:,:,:)  , ALLOCATABLE ::   znorthloc
      REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE ::   znorthglo
      TYPE(PTR_4D_dp), DIMENSION(:), ALLOCATABLE ::   ztabglo        ! array or pointer of arrays on which apply the b.c.
      !!----------------------------------------------------------------------
      !
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( ln_nnogather ) THEN      !==  no allgather exchanges  ==!

         !   ---   define number of exchanged lines   ---
         !
         ! In theory we should exchange only nn_hls lines.
         !
         ! However, some other points are duplicated in the north pole folding:
         !  - c_NFtype='T', grid=T : half of the last line (jpiglo/2+2:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=U : half of the last line (jpiglo/2+1:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=V : all the last line nn_hls+1 and (nn_hls+2:jpiglo-nn_hls)
         !  - c_NFtype='T', grid=F : all the last line (nn_hls+1:jpiglo-nn_hls)
         !  - c_NFtype='F', grid=T : 2 points of the last line (jpiglo/2+1 and jpglo-nn_hls)
         !  - c_NFtype='F', grid=U : no points are duplicated
         !  - c_NFtype='F', grid=V : half of the last line (jpiglo/2+1:jpiglo-nn_hls)
         !  - c_NFtype='F', grid=F : half of the last line (jpiglo/2+1:jpiglo-nn_hls-1)
         ! The order of the calculations may differ for these duplicated points (as, for example jj+1 becomes jj-1)
         ! This explain why these duplicated points may have different values even if they are at the exact same location.
         ! In consequence, we may want to force the folding on these points by setting l_full_nf_update = .TRUE.
         ! This is slightly slower but necessary to avoid different values on identical grid points!!
         !
         !!!!!!!!!           temporary switch off this optimisation ==> force TRUE           !!!!!!!!
         !!!!!!!!!  needed to get the same results without agrif and with agrif and no zoom  !!!!!!!!
         !!!!!!!!!                    I don't know why we must do that...                    !!!!!!!!
         l_full_nf_update = .TRUE.
         ! also force it if not restart during the first 2 steps (leap frog?)
         ll_add_line = l_full_nf_update .OR. ( ncom_stp <= nit000+1 .AND. .NOT. ln_rstart )
         
         ALLOCATE(ipjfld(ipf))                 ! how many lines do we exchange for each field?
         IF( ll_add_line ) THEN
            DO jf = 1, ipf                     ! Loop over the number of arrays to be processed
               ipjfld(jf) = khls + COUNT( (/ c_NFtype == 'T' .OR. cd_nat(jf) == 'V' .OR. cd_nat(jf) == 'F' /) )
            END DO
         ELSE
            ipjfld(:) = khls
         ENDIF
         
         ipj    = MAXVAL(ipjfld(:))            ! Max 2nd dimension of message transfers
         ipjtot = SUM(   ipjfld(:))            ! Total number of lines to be exchanged

         ! Index of modifying lines in input
         ALLOCATE( ijsnd(ipj, ipf), ijbuf(ipj, ipf), ijrcv(ipj, ipf), ii1st(ipj, ipf), iiend(ipj, ipf) )

         ij1 = 0
         DO jf = 1, ipf                        ! Loop over the number of arrays to be processed
            !
            DO jj = 1, khls   ! first khls lines (starting from top) must be fully defined
               ii1st(jj, jf) = 1
               iiend(jj, jf) = jpi
            END DO
            !
            ! what do we do with line khls+1 (starting from top)
            IF( c_NFtype == 'T' ) THEN          ! *  North fold  T-point pivot
               SELECT CASE ( cd_nat(jf) )
               CASE ('T','W')   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+2)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls)
               CASE ('U'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls)
               CASE ('V'    )   ;   i012 = 2   ;   ii1st(khls+1, jf) = 1                 ;   iiend(khls+1, jf) = jpi
               CASE ('F'    )   ;   i012 = 2   ;   ii1st(khls+1, jf) = 1                 ;   iiend(khls+1, jf) = jpi
               END SELECT
            ENDIF
            IF( c_NFtype == 'F' ) THEN          ! *  North fold  F-point pivot
               SELECT CASE ( cd_nat(jf) )
               CASE ('T','W')   ;   i012 = 0   ! we don't touch line khls+1
               CASE ('U'    )   ;   i012 = 0   ! we don't touch line khls+1
               CASE ('V'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls  )
               CASE ('F'    )   ;   i012 = 1   ;   ii1st(khls+1, jf) = mi0(jpiglo/2+1)   ;   iiend(khls+1, jf) = mi1(jpiglo-khls-1)
               END SELECT
            ENDIF
            !
            DO jj = 1, ipjfld(jf)
               ij1 = ij1 + 1
               ijsnd(jj,jf) = jpj - 2*khls + jj - i012   ! sent lines (from bottom of sent lines)
               ijbuf(jj,jf) = ij1                        ! gather all lines in the snd/rcv buffers
               ijrcv(jj,jf) = jpj - jj + 1               ! recv lines (from the top -> reverse order for jj)
            END DO
            !
         END DO
         !
         i0max = jpimax - 2 * khls                                    ! we are not sending the halos
         ALLOCATE( zbufs(i0max,ipjtot,ipk,ipl), ireq_s(nfd_nbnei) )   ! store all the data to be sent in a buffer array
         ibuffsize = i0max * ipjtot * ipk * ipl
         !
         ! fill the send buffer with all the lines
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
            DO jj = 1, ipjfld(jf)
               ij1 = ijbuf(jj,jf)
               ij2 = ijsnd(jj,jf)
               DO ji = Nis0, Nie0       ! should not use any other value
                  iib = ji - Nis0 + 1
                  zbufs(iib,ij1,jk,jl) = ptab(jf)%pt4d(ji,ij2,jk,jl)
               END DO
               DO ji = Ni_0+1, i0max    ! avoid sending uninitialized values (make sure we don't use it)
                  zbufs(ji,ij1,jk,jl) = HUGE(0._dp)   ! make sure we don't use it...
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         ! start waiting time measurement
         IF( ln_timing ) CALL tic_tac(.TRUE.)
         !
         ! send the same buffer data to all neighbourgs as soon as possible
         DO jn = 1, nfd_nbnei
            iproc = nfd_rknei(jn)
            IF( iproc /= narea-1 .AND. iproc /= -1 ) THEN
               IF(PRESENT(fTag)) THEN
                  CALL MPI_Isend( zbufs, ibuffsize, MPI_DOUBLE_PRECISION, iproc, fTag, mpi_comm_oce, ireq_s(jn), ierr )
               ELSE
                  CALL MPI_Isend( zbufs, ibuffsize, MPI_DOUBLE_PRECISION, iproc, 5, mpi_comm_oce, ireq_s(jn), ierr )
               ENDIF
            ELSE
               ireq_s(jn) = MPI_REQUEST_NULL
            ENDIF
         END DO
         !
         ALLOCATE( zbufr(i0max,ipjtot,ipk,ipl,nfd_nbnei), ireq_r(nfd_nbnei) ) 
         !
         DO jn = 1, nfd_nbnei
            !
            iproc = nfd_rknei(jn)
            !
            IF(           iproc == -1 ) THEN   ! No neighbour (land proc that was suppressed)
               !
               ireq_r(jn) = MPI_REQUEST_NULL                ! no message to be received
               zbufr(:,:,:,:,jn) = HUGE(0._dp)   ! default: define it and make sure we don't use it...
               SELECT CASE ( kfillmode )
               CASE ( jpfillnothing )                       ! no filling 
               CASE ( jpfillcopy    )                       ! filling with inner domain values
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipjfld(jf)
                        ij1 = ijbuf(jj,jf)
                        ij2 = ijsnd(jj,jf)                                      ! we will use only the first value, see init_nfdcom
                        zbufr(1,ij1,jk,jl,jn) = ptab(jf)%pt4d(Nis0,ij2,jk,jl)   ! chose to take the 1st inner domain point
                     END DO
                  END DO   ;   END DO   ;   END DO
               CASE ( jpfillcst     )                       ! filling with constant value
                  zbufr(1,:,:,:,jn) = pfillval              ! we will use only the first value, see init_nfdcom
               END SELECT
               !
            ELSE IF( iproc == narea-1 ) THEN   ! get data from myself!
               !
               ireq_r(jn) = MPI_REQUEST_NULL                ! no message to be received
               DO jf = 1, ipf   ;   DO jl = 1, ipl  ;   DO jk = 1, ipk
                  DO jj = 1, ipjfld(jf)
                     ij1 = ijbuf(jj,jf)
                     ij2 = ijsnd(jj,jf)
                     DO ji = Nis0, Nie0                     ! should not use any other value
                        iib = ji - Nis0 + 1
                        zbufr(iib,ij1,jk,jl,jn) = ptab(jf)%pt4d(ji,ij2,jk,jl)
                     END DO
                  END DO
               END DO   ;   END DO   ;   END DO
               !
            ELSE                               ! get data from a neighbour trough communication
               IF(PRESENT(fTag)) THEN
                  CALL MPI_Irecv( zbufr(:,:,:,:,jn), ibuffsize, MPI_DOUBLE_PRECISION, iproc, fTag, mpi_comm_oce, ireq_r(jn), ierr )
               ELSE
                  CALL MPI_Irecv( zbufr(:,:,:,:,jn), ibuffsize, MPI_DOUBLE_PRECISION, iproc, 5, mpi_comm_oce, ireq_r(jn), ierr )
               ENDIF
            ENDIF
            !
         END DO   ! nfd_nbnei
         !
         CALL mpi_waitall(nfd_nbnei, ireq_r, MPI_STATUSES_IGNORE, ierr)   ! wait for all Irecv
         !
         IF( ln_timing ) CALL tic_tac(.FALSE.)
         !
         ! North fold boundary condition
         !
         DO jf = 1, ipf
            !
            SELECT CASE ( cd_nat(jf) )     ! which grid number?
            CASE ('T','W')   ;   iig = 1   ! T-, W-point
            CASE ('U')       ;   iig = 2   ! U-point
            CASE ('V')       ;   iig = 3   ! V-point
            CASE ('F')       ;   iig = 4   ! F-point
            END SELECT
            !
            DO jl = 1, ipl   ;   DO jk = 1, ipk
               !
               ! if T point with F-point pivot : must be done first
               !    --> specific correction of 3 points near the 2 pivots (to be clean, usually masked -> so useless) 
               IF( c_NFtype == 'F' .AND. iig == 1 ) THEN
                  ij1 = jpj - khls     ! j-index in the receiving array
                  ij2 = 1              ! only 1 line in the buffer
                  DO ji = mi0(khls), mi1(khls)               ! change because of EW periodicity as we also change jpiglo-khls
                     iib = nfd_jisnd(mi0(       khls),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0(       khls),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
                  DO ji = mi0(jpiglo/2+1), mi1(jpiglo/2+1)
                     iib = nfd_jisnd(mi0( jpiglo/2+1),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0( jpiglo/2+1),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
                  DO ji = mi0(jpiglo-khls), mi1(jpiglo-khls)
                     iib = nfd_jisnd(mi0(jpiglo-khls),iig)   ! i-index in the buffer
                     iin = nfd_rksnd(mi0(jpiglo-khls),iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = zbufr(iib,ij2,jk,jl,iin)   ! no psgn(jf)
                  END DO
               ENDIF
               !
               ! Apply the North pole folding.
               DO jj = 1, ipjfld(jf)   ! for all lines to be exchanged for this field
                  ij1 = ijrcv(jj,jf)   ! j-index in the receiving array
                  ij2 = ijbuf(jj,jf)   ! j-index in the buffer
                  iis = ii1st(jj,jf)   ! stating i-index in the receiving array
                  iie = iiend(jj,jf)   !  ending i-index in the receiving array
                  DO ji = iis, iie 
                     iib = nfd_jisnd(ji,iig)   ! i-index in the buffer
                     iin = nfd_rksnd(ji,iig)   ! neigbhour-index in the buffer
                     IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                     ptab(jf)%pt4d(ji,ij1,jk,jl) = psgn(jf) * zbufr(iib,ij2,jk,jl,iin)
                  END DO
               END DO
               !
               ! re-apply periodocity when we modified the eastern side of the inner domain (and not the full line)
               IF( c_NFtype == 'T' ) THEN          ! *  North fold  T-point pivot
                  IF(     iig <= 2 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls)   ! 'T','W','U': update west halo
                  ELSE                      ;   iis = 1        ;   iie = 0           ! 'V','F'    : full line already exchanged
                  ENDIF
               ENDIF
               IF( c_NFtype == 'F' ) THEN          ! *  North fold  F-point pivot
                  IF(     iig <= 2 ) THEN   ;   iis = 1        ;   iie = 0           ! 'T','W','U': nothing to do
                  ELSEIF( iig == 3 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls)   ! 'V'        : update west halo
                  ELSEIF( khls > 1 ) THEN   ;   iis = mi0(1)   ;   iie = mi1(khls-1) ! 'F' and khls > 1
                  ELSE                      ;   iis = 1        ;   iie = 0           ! 'F' and khls == 1 : nothing to do
                  ENDIF
               ENDIF
               jj  = ipjfld(jf)     ! only for the last line of this field
               ij1 = ijrcv(jj,jf)   ! j-index in the receiving array
               ij2 = ijbuf(jj,jf)   ! j-index in the buffer
               DO ji = iis, iie
                  iib = nfd_jisnd(ji,iig)   ! i-index in the buffer
                  iin = nfd_rksnd(ji,iig)   ! neigbhour-index in the buffer
                  IF( nfd_rknei(iin) == -1 .AND. kfillmode == jpfillnothing )   CYCLE
                  ptab(jf)%pt4d(ji,ij1,jk,jl) = psgn(jf) * zbufr(iib,ij2,jk,jl,iin)
               END DO
               !               
            END DO   ;   END DO   ! ipl   ; ipk
            !               
         END DO   ! ipf
       
         !
         DEALLOCATE( zbufr, ireq_r, ijsnd, ijbuf, ijrcv, ii1st, iiend, ipjfld )
         !
         CALL mpi_waitall(nfd_nbnei, ireq_s, MPI_STATUSES_IGNORE, ierr)   ! wait for all Isend
         !
         DEALLOCATE( zbufs, ireq_s )
         !
      ELSE                             !==  allgather exchanges  ==!
         !
         ! how many lines do we exchange at max? -> ipj    (no further optimizations in this case...)
         ipj =      khls + 2
         ! how many lines do we     need at max? -> ipj2   (no further optimizations in this case...)
         ipj2 = 2 * khls + 2
         !
         i0max = jpimax - 2 * khls
         ibuffsize = i0max * ipj * ipk * ipl * ipf
         ALLOCATE( znorthloc(i0max,ipj,ipk,ipl,ipf), znorthglo(i0max,ipj,ipk,ipl,ipf,ndim_rank_north) )
         !
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk               ! put in znorthloc ipj j-lines of ptab
            DO jj = 1, ipj
               ij2 = jpj - ipj2 + jj                        ! the first ipj lines of the last ipj2 lines
               DO ji = 1, Ni_0
                  ii2 = Nis0 - 1 + ji                       ! inner domain: Nis0 to Nie0
                  znorthloc(ji,jj,jk,jl,jf) = ptab(jf)%pt4d(ii2,ij2,jk,jl)
               END DO
               DO ji = Ni_0+1, i0max
                  znorthloc(ji,jj,jk,jl,jf) = HUGE(0._dp)   ! avoid sending uninitialized values (make sure we don't use it)
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         ! start waiting time measurement
         IF( ln_timing ) CALL tic_tac(.TRUE.)
         CALL MPI_ALLGATHER( znorthloc, ibuffsize, MPI_DOUBLE_PRECISION, znorthglo, ibuffsize, MPI_DOUBLE_PRECISION, ncomm_north, ierr )
         ! stop waiting time measurement
         IF( ln_timing ) CALL tic_tac(.FALSE.)
         DEALLOCATE( znorthloc )
         ALLOCATE( ztabglo(ipf) )
         DO jf = 1, ipf
            ALLOCATE( ztabglo(jf)%pt4d(jpiglo,ipj2,ipk,ipl) )
         END DO
         !
         ! need to fill only the first ipj lines of ztabglo as lbc_nfd don't use the last khls lines
         ijnr = 0
         DO jr = 1, jpni                                                        ! recover the global north array
            iproc = nfproc(jr)
            impp  = nfimpp(jr)
            ipi   = nfjpi( jr) - 2 * khls                       ! corresponds to Ni_0 but for subdomain iproc
            IF( iproc == -1 ) THEN   ! No neighbour (land proc that was suppressed)
              !
               SELECT CASE ( kfillmode )
               CASE ( jpfillnothing )               ! no filling
                  CALL ctl_stop( 'STOP', 'mpp_nfd_generic : cannot use jpfillnothing with ln_nnogather = F')
               CASE ( jpfillcopy    )               ! filling with inner domain values
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipj
                        ij2 = jpj - ipj2 + jj                    ! the first ipj lines of the last ipj2 lines
                        DO ji = 1, ipi
                           ii1 = impp + khls + ji - 1            ! corresponds to mig(khls + ji) but for subdomain iproc
                           ztabglo(jf)%pt4d(ii1,jj,jk,jl) = ptab(jf)%pt4d(Nis0,ij2,jk,jl) ! chose to take the 1st inner domain point
                        END DO
                     END DO
                  END DO   ;   END DO   ;   END DO
               CASE ( jpfillcst     )               ! filling with constant value
                  DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                     DO jj = 1, ipj
                        DO ji = 1, ipi
                           ii1 = impp + khls + ji - 1            ! corresponds to mig(khls + ji) but for subdomain iproc
                           ztabglo(jf)%pt4d(ii1,jj,jk,jl) = pfillval
                        END DO
                     END DO
                 END DO   ;   END DO   ;   END DO
               END SELECT
               !
            ELSE
               ijnr = ijnr + 1
               DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk
                  DO jj = 1, ipj
                     DO ji = 1, ipi
                        ii1 = impp + khls + ji - 1               ! corresponds to mig(khls + ji) but for subdomain iproc
                        ztabglo(jf)%pt4d(ii1,jj,jk,jl) = znorthglo(ji,jj,jk,jl,jf,ijnr)
                     END DO
                  END DO
               END DO   ;   END DO   ;   END DO
            ENDIF
            !
         END DO   ! jpni
         DEALLOCATE( znorthglo )
         !
         DO jf = 1, ipf
            CALL lbc_nfd( ztabglo(jf:jf), cd_nat(jf:jf), psgn(jf:jf), khls, 1 )   ! North fold boundary condition
            DO jl = 1, ipl   ;   DO jk = 1, ipk                  ! e-w periodicity
               DO jj = 1, khls + 1
                  ij1 = ipj2 - (khls + 1) + jj                   ! need only the last khls + 1 lines until ipj2
                  ztabglo(jf)%pt4d(            1:  khls,ij1,jk,jl) = ztabglo(jf)%pt4d(jpiglo-2*khls+1:jpiglo-khls,ij1,jk,jl)
                  ztabglo(jf)%pt4d(jpiglo-khls+1:jpiglo,ij1,jk,jl) = ztabglo(jf)%pt4d(         khls+1:     2*khls,ij1,jk,jl)
               END DO
            END DO   ;   END DO
         END DO     
         !
         DO jf = 1, ipf   ;   DO jl = 1, ipl   ;   DO jk = 1, ipk               ! Scatter back to ARRAY_IN
            DO jj = 1, khls + 1
               ij1 = jpj  - (khls + 1) + jj   ! last khls + 1 lines until jpj
               ij2 = ipj2 - (khls + 1) + jj   ! last khls + 1 lines until ipj2
               DO ji= 1, jpi
                  ii2 = mig(ji)
                  ptab(jf)%pt4d(ji,ij1,jk,jl) = ztabglo(jf)%pt4d(ii2,ij2,jk,jl)
               END DO
            END DO
         END DO   ;   END DO   ;   END DO
         !
         DO jf = 1, ipf
            DEALLOCATE( ztabglo(jf)%pt4d )
         END DO
         DEALLOCATE( ztabglo )
         !
      ENDIF   ! ln_nnogather
      !
   END SUBROUTINE mpp_nfd_dp


   !!======================================================================
END MODULE lbcnfd
