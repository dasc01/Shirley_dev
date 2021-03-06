!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine phq_summary
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output the quantities which have been read
  !    from the punch file, and the quantities computed in the phq_setup
  !    file.
  !
  !    if iverbosity = 0 only a partial summary is done.
  !
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp, atm, tau, ntyp => nsp, amass
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at, bg, ibrav, alat, omega, celldm
  USE klist,         ONLY : lgauss, smearing, degauss, ngauss, nkstot, xk, wk
  USE grid_dimensions,ONLY: nr1, nr2, nr3
  USE gvect,         ONLY : gcutm, ngm
  USE gvecs,       ONLY : doublegrid, dual, gcutms, ngms
  USE smooth_grid_dimensions, ONLY : nr1s, nr2s, nr3s
  USE symm_base,     ONLY : s, sr, ftau, sname, t_rev
  USE constants,     ONLY : amconv
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,      ONLY : lspinorb, domag
  USE funct,         ONLY : write_dft_name
  USE printout_base, ONLY : title
  USE gamma_gamma,   ONLY : with_symmetry, nasr
  USE control_ph,    ONLY : lgamma_gamma, lnoloc, lrpa, zue, epsil, ldisp, &
                            nmix_ph, alpha_mix, tr2_ph, zeu
  USE freq_ph,       ONLY : fpol, nfs, fiu
  USE partial,       ONLY : atomo, nat_todo, all_comp, done_irr, comp_irr
  USE modes,         ONLY : u, npert, irotmq, irgq, minus_q, nsymq, nirr, &
                            name_rap_mode
  USE qpoint,        ONLY : xq
  USE ramanm,        ONLY : lraman, elop
  USE control_flags, ONLY : iverbosity
  USE wvfct,         ONLY : ecutwfc

  implicit none

  integer :: i, mu, nu, ipol, apol, na, isymq, isym, nsymtot, &
       ik, irr, imode0, iu
  ! generic counter
  ! counter on modes
  ! counter on modes
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on symmetries
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points
  ! counter on beta functions
  ! counter on irreducible representation
  ! the first mode

  real(DP) :: ft1, ft2, ft3, xkg (3)
  ! fractionary translations
  ! k point in crystal coordinates

  !
  WRITE( stdout, 100) title, ibrav, alat, omega, nat, ntyp, &
       ecutwfc, ecutwfc * dual, tr2_ph, alpha_mix (1), &
       nmix_ph
100 format (/,5x,a75,/,/,5x, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'charge density cut-off    = ',f12.4,'  Ry',/,5x, &
       &     'convergence threshold     = ',1pe12.1,/,5x, &
       &     'beta                      = ',0pf12.4,/,5x, &
       &     'number of iterations used = ',i12)

  CALL write_dft_name ( )

  !
  !  Here add a message if this is a noncollinear or a spin_orbit calculation
  !
  IF (noncolin) THEN
     IF (lspinorb) THEN
        IF (domag) THEN
           WRITE( stdout, '(5x, "Noncollinear calculation with spin-orbit",/)')
        ELSE
           WRITE( stdout, '(5x, "Non magnetic calculation with spin-orbit",/)')
        ENDIF
     ELSE
        WRITE( stdout, '(5x, "Noncollinear calculation without spin-orbit",/)')
     END IF
  ELSE
     WRITE(stdout,'(/)')
  END IF
  !
  !    and here more detailed information. Description of the unit cell
  !
  WRITE( stdout, '(2(3x,3(2x,"celldm(",i1,")=",f11.5),/))') (i, &
       celldm (i) , i = 1, 6)
  WRITE( stdout, '(5x, &
       &  "crystal axes: (cart. coord. in units of a_0)",/, &
       &         3(15x,"a(",i1,") = (",3f8.4," )  ",/ ) )')  (apol, &
       & (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  WRITE( stdout, '(5x, &
       &"reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &         3(15x,"b(",i1,") = (",3f8.4," )  ",/ ) )')  (apol, &
       & (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  !
  !    description of the atoms inside the unit cell
  !
  WRITE( stdout, '(/, 5x,"Atoms inside the unit cell: ")')
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.  atom      mass ", &
       &                "          positions (a_0 units)")')

  WRITE( stdout, '(7x,i2,5x,a6,f8.4,"   tau(",i2, &
       &                              ") = (",3f11.5,"  )")')  &
       &(na, atm (ityp (na) ) , amass (ityp (na) ), na,  &
       &(tau (ipol, na) , ipol = 1, 3) , na = 1, nat)
  WRITE( stdout, '(/,5x,"Computing dynamical matrix for ")')
  WRITE( stdout, '(20x,"q = (",3f12.7," )")') (xq (ipol) , ipol = 1, 3)
  !
  !   description of symmetries
  !
  WRITE( stdout, * )
  if (nsymq.le.1.and..not.minus_q) then
     WRITE( stdout, '(5x,"No symmetry!")')
  else
     if (minus_q) then
        WRITE( stdout, '(5x,i2," Sym.Ops. (with q -> -q+G )",/)') &
             nsymq + 1
     else
        WRITE( stdout, '(5x,i2," Sym.Ops. (no q -> -q+G )",/)') nsymq
     endif

  endif
  if (iverbosity.eq.1) then

     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     if (minus_q) then
        nsymtot = nsymq + 1
     else
        nsymtot = nsymq
     endif
     do isymq = 1, nsymtot
        if (isymq.gt.nsymq) then
           isym = irotmq
           WRITE( stdout, '(/,5x,"This transformation sends q -> -q+G")')
        else
           isym = irgq (isymq)
        endif
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isymq, sname (isym)
        IF (noncolin.and.domag) &
            WRITE(stdout,'(1x, "Time Reversal",i3)') t_rev(isym)

        if (ftau (1, isym) .ne.0.or.ftau (2, isym) .ne.0.or.ftau (3, &
             isym) .ne.0) then
           ft1 = at (1, 1) * ftau (1, isym) / nr1 + at (1, 2) * ftau ( &
                2, isym) / nr2 + at (1, 3) * ftau (3, isym) / nr3
           ft2 = at (2, 1) * ftau (1, isym) / nr1 + at (2, 2) * ftau ( &
                2, isym) / nr2 + at (2, 3) * ftau (3, isym) / nr3
           ft3 = at (3, 1) * ftau (1, isym) / nr1 + at (3, 2) * ftau ( &
                2, isym) / nr2 + at (3, 3) * ftau (3, isym) / nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &                    " )    f =( ",f10.7," )")') isymq,  (s (1, &
                & ipol, isym) , ipol = 1, 3) , DBLE (ftau (1, isym) )  / DBLE (nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )")')  (s (2, ipol, &
                &isym) , ipol = 1, 3) , DBLE (ftau (2, isym) )  / DBLE (nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), &
                &                    " )       ( ",f10.7," )"/)')  (s (3, ipol, &
                & isym) , ipol = 1, 3) , DBLE (ftau (3, isym) )  / DBLE (nr3)
           WRITE( stdout, '(1x,"cart.",4x,"s(",i2,") = (",3f11.7, &
                &                    " )    f =( ",f10.7," )")') isymq,  &
                &  (sr (1, ipol,isym) , ipol = 1, 3) , ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                   (sr (2, ipol,isym) , ipol = 1, 3) , ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
               & (sr (3, ipol,isym) , ipol = 1, 3) , ft3
        else
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                &                    " )")') isymq,  (s (1, ipol, isym) , ipol = &
                &1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")') (s (2, ipol, isym) &
                , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s (3, ipol, &
                isym) , ipol = 1, 3)
           WRITE( stdout, '(1x,"cart.",4x,"s(",i2,") = (",3f11.7, " )")') &
               isymq,  (sr (1, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")') &
               (sr (2, ipol,isym) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') &
               (sr (3, ipol,isym) , ipol = 1, 3)
        endif
     enddo
  endif
  !
  !     Description of the reciprocal lattice vectors
  !
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngm, nr1, nr2, nr3

  if (doublegrid) WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
       &                      i7," G-vectors)","  smooth grid: (",i3, &
       &                      ",",i3,",",i3,")")') gcutms, ngms, nr1s, nr2s, nr3s
  if (.NOT.lgauss) then
     WRITE( stdout, '(5x,"number of k points=",i6)') nkstot
  else
     WRITE( stdout, '(/5x,"number of k points=", i6, 2x, &
          &             a," smearing, width (Ry)=",f8.4)') &
          &             nkstot, TRIM(smearing), degauss
  endif
  IF (iverbosity==1 .or. (nkstot < 100 .and. .not.ldisp) ) then
     WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
     do ik = 1, nkstot
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
             (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
     enddo
  ENDIF
  if (iverbosity.eq.1) then
     WRITE( stdout, '(/23x,"cryst. coord.")')
     do ik = 1, nkstot
        do ipol = 1, 3
           xkg (ipol) = at (1, ipol) * xk (1, ik) + at (2, ipol) * xk (2, &
                ik) + at (3, ipol) * xk (3, ik)
           ! xkg are the components  of xk in the reciprocal lattice basis
        enddo
        WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     enddo

  endif

  CALL print_ps_info ( )

  IF (lgamma_gamma) &
        WRITE(stdout,'(/5x,"k=gamma and q=gamma tricks are used")')

  IF (epsil) THEN
     WRITE( stdout, '(//5x,"Electric field:")')
     IF (lgamma_gamma) THEN
        WRITE(stdout,'(5x,"Dielectric constant and polarizability")')
     ELSE
        WRITE( stdout, '(5x,"Dielectric constant")')
     END IF
     IF (zue.AND.zeu) THEN
             WRITE( stdout, '(5x,"Born effective charges in two ways ")' )
     ELSEIF (zue) THEN
             WRITE( stdout, '(5x,"Born effective charges as d P / d u")')
     ELSEIF (zeu) THEN
             WRITE( stdout, '(5x,"Born effective charges as d Force / d E")')
     END IF
     IF (lraman) &
          WRITE( stdout, '(5x,"Raman tensor")')
     IF (elop) &
          WRITE( stdout, '(5x,"Electro-optic tensor")')
     IF (fpol)  THEN
        WRITE( stdout, '(5x,"Frequency Dependent Polarizability at (Ry) ")' )
        WRITE( stdout, '(5x,8(f9.4,"i"))') (fiu(iu), iu=nfs,1,-1)
     ENDIF
  ENDIF

  WRITE( stdout, '(//5x,"Atomic displacements:")')
  WRITE( stdout, '(5x,"There are ",i3," irreducible representations")') nirr
  imode0 = 0
  DO irr = 1, nirr
     IF (done_irr (irr) .eq.1) then
        WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &                  " modes -",a,"  Done")') irr, npert (irr),&
                      TRIM( name_rap_mode(irr) )
     ELSEIF (comp_irr (irr) .eq.1) then
        WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &             " modes -",a," To be done")') irr, npert (irr), &
                      TRIM( name_rap_mode(irr) )
     ELSEIF (comp_irr (irr) .eq.0) THEN
        IF (lgamma_gamma) THEN
           IF ((irr-1)/3+1==nasr) THEN
              WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
                 &     " modes - Calculated using asr")') irr, npert (irr)
              done_irr(irr) = 1
           ELSEIF (with_symmetry(irr)==1) THEN
              WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
                 &     " modes - Calculated using symmetry")') irr, npert (irr)
              done_irr(irr) = 1
           ELSE
              WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
                &     " modes - Not done in this run")') irr, npert (irr)
           ENDIF
        ELSE
           WRITE( stdout, '(/, 5x,"Representation ",i5,i7, &
             &     " modes -",a,"  Not done in this run")') irr, npert (irr), &
                      TRIM( name_rap_mode(irr) )
        ENDIF
     ENDIF

     if (iverbosity.eq.1) then
        WRITE( stdout, '(5x,"Phonon polarizations are as follows:",/)')
        if (npert (irr) .eq.1) then
           WRITE( stdout, '(20x," mode # ",i3)') imode0 + 1
           WRITE( stdout, '(20x," (",2f10.5,"   ) ")')  ( (u (mu, nu) ,&
                &nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        elseif (npert (irr) .eq.2) then
           WRITE( stdout, '(2(10x," mode # ",i3,16x))') imode0 + 1, &
                imode0 + 2
           WRITE( stdout, '(2(10x," (",2f10.5,"   ) "))')  ( (u (mu, nu) , nu &
                &= imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        else
           WRITE( stdout, '(4x,3(" mode # ",i3,13x))') imode0 + 1, imode0 &
                + 2, imode0 + 3
           WRITE( stdout, '((5x,3("(",2f10.5," ) ")))') ( (u (mu, nu) , &
                nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
        endif
        imode0 = imode0 + npert (irr)
     endif
  enddo
  if (.not.all_comp) then
     WRITE( stdout, '(/,5x,"Compute atoms: ",8(i5,","))') (atomo (na) &
          , na = 1, nat_todo)
  endif
  write(stdout,'(/)')
  !
  CALL flush_unit( stdout )
  !
  return
end subroutine phq_summary
