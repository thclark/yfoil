INCLUDE 'splprocs.f90'
INCLUDE 'epspsi.f90'
INCLUDE 'nacax.f90'
!+
PROGRAM Naca456Fixture
! ------------------------------------------------------------------------------
! PURPOSE - Write the reference ordinates of one NACA section at full precision, for the
!   yFoil fixtures (tests/fixtures/naca456/). The section is described by the same NAMELIST
!   /NACA/ that PDAS naca456 reads (see input.txt in the naca456 distribution), read from
!   standard input; the thickness form, mean line and their slopes, and the combined upper and
!   lower surfaces, are evaluated at naca456's "very fine" station table (LoadX, denCode = 3)
!   by naca456's own module procedures (NacaAuxilary, EpsilonPsi, SplineProcedures), unchanged.
!   Only the output differs from naca456: ES24.16 (17 significant figures, the minimum that
!   round-trips a double) instead of F10.6, and one file with every column.
!
!   Build (scripts/naca456-fixtures.sh):
!     gfortran -fdefault-real-8 -ffp-contract=off -O0 -I <naca456 sources> driver.f90
!
! AUTHOR - yFoil project; the included modules are public domain (Ralph L. Carmichael, PDAS)
! ------------------------------------------------------------------------------
USE NACAauxilary
IMPLICIT NONE

  REAL:: a = 1.0            ! chordwise extent of uniform loading (6-series mean lines)
  CHARACTER(LEN=2):: camber = '0'   ! '0', '2', '3', '3R', '6', '6A'
  REAL:: cl = 0.0           ! design lift coefficient (3, 3R, 6, 6A mean lines)
  REAL:: cmax = 0.0         ! maximum camber (2-digit mean line)
  INTEGER:: denCode = 3     ! station table density (3 = very fine, 98 stations)
  REAL:: leIndex = 6.0      ! leading-edge radius index (4-digit modified)
  CHARACTER(LEN=80):: name = ' '    ! the designation, e.g. 'NACA 63-415'
  CHARACTER(LEN=3):: profile = '4'  ! '4', '4M', '63' ... '67', '63A', '64A', '65A'
  REAL:: toc = 0.1          ! thickness ratio
  REAL:: xmaxc = 0.4        ! position of maximum camber (2, 3, 3R mean lines)
  REAL:: xmaxt = 0.3        ! position of maximum thickness (4-digit modified)

  INTEGER:: k, nTable
  REAL, DIMENSION(1000):: xTable
  REAL, ALLOCATABLE, DIMENSION(:):: x, yt, ytp, ym, ymp, xu, yu, xl, yl

  NAMELIST /NACA/ a, camber, cl, cmax, denCode, leIndex, name, profile, toc, xmaxc, xmaxt
!-------------------------------------------------------------------------------
  READ(*, NML=NACA)
  profile = AdjustL(profile)
  camber = AdjustL(camber)

  CALL LoadX(denCode, nTable, xTable)
  ALLOCATE(x(nTable), yt(nTable), ytp(nTable), ym(nTable), ymp(nTable))
  ALLOCATE(xu(nTable), yu(nTable), xl(nTable), yl(nTable))
  x = xTable(1:nTable)

  SELECT CASE(profile)
    CASE('4')
      CALL Thickness4(toc, x, yt, ytp)
    CASE('4M', '4m')
      CALL Thickness4M(toc, leIndex, xmaxt, x, yt, ytp)
    CASE('63')
      CALL Thickness6(1, toc, x, yt, ytp)
    CASE('64')
      CALL Thickness6(2, toc, x, yt, ytp)
    CASE('65')
      CALL Thickness6(3, toc, x, yt, ytp)
    CASE('66')
      CALL Thickness6(4, toc, x, yt, ytp)
    CASE('67')
      CALL Thickness6(5, toc, x, yt, ytp)
    CASE('63A', '63a')
      CALL Thickness6(6, toc, x, yt, ytp)
    CASE('64A', '64a')
      CALL Thickness6(7, toc, x, yt, ytp)
    CASE('65A', '65a')
      CALL Thickness6(8, toc, x, yt, ytp)
    CASE DEFAULT
      WRITE(0, *) "Not a valid profile: " // profile
      STOP 1
  END SELECT

  SELECT CASE(camber)
    CASE('0')
      ym = 0.0
      ymp = 0.0
    CASE('2')
      CALL MeanLine2(cmax, xmaxc, x, ym, ymp)
    CASE('3')
      CALL MeanLine3(cl, xmaxc, x, ym, ymp)
    CASE('3R', '3r')
      CALL MeanLine3Reflex(cl, xmaxc, x, ym, ymp)
    CASE('6')
      CALL MeanLine6(a, cl, x, ym, ymp)
    CASE('6A', '6a', '6M', '6m')
      CALL MeanLine6M(cl, x, ym, ymp)
    CASE DEFAULT
      WRITE(0, *) "Not a valid mean line: " // camber
      STOP 1
  END SELECT

  CALL CombineThicknessAndCamber(x, yt, ym, ymp, xu, yu, xl, yl)

  WRITE(*, '(A)') "# naca456 reference ordinates (PDAS naca456 modules, NASA TM-4741 algorithm), " // &
    "driver scripts/naca456-fixtures/driver.f90"
  WRITE(*, '(A)') "# naca456_aux_version " // AUX_VERSION
  WRITE(*, '(A)') "# name " // Trim(name)
  WRITE(*, '(A)') "# profile " // Trim(profile)
  WRITE(*, '(A, ES24.16)') "# toc", toc
  WRITE(*, '(A, ES24.16)') "# leindex", leIndex
  WRITE(*, '(A, ES24.16)') "# xmaxt", xmaxt
  WRITE(*, '(A)') "# camber " // Trim(camber)
  WRITE(*, '(A, ES24.16)') "# cmax", cmax
  WRITE(*, '(A, ES24.16)') "# xmaxc", xmaxc
  WRITE(*, '(A, ES24.16)') "# cl", cl
  WRITE(*, '(A, ES24.16)') "# a", a
  WRITE(*, '(A, I4)') "# stations", nTable
  WRITE(*, '(A)') "# columns x yt ytp ym ymp xu yu xl yl"
  DO k = 1, nTable
    WRITE(*, '(9ES24.16)') x(k), yt(k), ytp(k), ym(k), ymp(k), xu(k), yu(k), xl(k), yl(k)
  END DO
  STOP
END PROGRAM Naca456Fixture
