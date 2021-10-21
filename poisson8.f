
C**********************************************************************
      PROGRAM ADAPTIVE
C**********************************************************************
C
C     This program uses the new diagonal translation FMM
C     to solve a two dimensional Poisson equation on an
C     adaptive grid.  The grid is a quad tree structure 
C     with the restriction that no two touching boxes can
C     be more than one level apart.
C
C     All points lie in the unit box centered
C     at the origin.
C
C     FRIGHT(j,i)  = right-hand side
C        This denotes that value of the right hand
C        side at the jth point in the ith box.
C
C     POT(j,i) =  solution
C        This denotes that value of the solution
C        at the jth point in the ith box.
C
C     XF(j,i) and YF(j,i) denote the corresponding
C        x and y coordinates for the above values.
C
C     The error is checked by comparison with direct calculation.
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  MAXWRK, MAXBOXES, MAXLEVEL, LENMAPS, IWRK
      INTEGER  NCNTR
      PARAMETER(MAXBOXES=30000)
      PARAMETER(MAXWRK=20000000)
      PARAMETER(IWRK=20000000)
      PARAMETER(LENMAPS=30000)
      PARAMETER(MAXLEVEL=15)
      INTEGER  NLEV
      INTEGER  NINBOX, NBOXES, IBOX, JCNTR
      INTEGER  I, J, L
      INTEGER  IPREC, IPERIODTEMP, IPERIOD
      INTEGER  LEVELBOX(MAXBOXES), IPARENTBOX(MAXBOXES)
      INTEGER  ICHILDBOX(4,MAXBOXES), ICOLLEAGBOX(9,MAXBOXES)
      INTEGER  ICOLBOX(MAXBOXES), IROWBOX(MAXBOXES)
      INTEGER  NBLEVEL(0:MAXLEVEL), IBOXLEV(MAXBOXES)
      INTEGER  ITEMPARRAY(MAXBOXES)
      INTEGER  ISTARTLEV(0:MAXLEVEL)
      INTEGER  IFLAG(MAXBOXES)
      INTEGER  IWORK(MAXWRK)
      INTEGER  IFIXFLAG
      INTEGER  LADDER(0:MAXLEVEL+1)
      REAL *8  FRIGHT(64,MAXBOXES), POT(64,MAXBOXES)
      REAL *8  HEXACT, HEXACTX, HEXACTY
      REAL *8  MAP(LENMAPS)
      REAL *8  TEST(64,MAXBOXES)
      REAL *8  POTX(64,MAXBOXES), POTY(64,MAXBOXES)
      REAL *8  LAP(64,MAXBOXES)
      REAL *8  XF(64,MAXBOXES), YF(64,MAXBOXES)
      REAL *8  XFPRACTICAL(81,MAXBOXES), YFPRACTICAL(81,MAXBOXES)
      REAL *8  POTPRACTICAL(81,MAXBOXES)
      REAL *8  POTPRACTICALX(81,MAXBOXES), POTPRACTICALY(81,MAXBOXES)
      REAL *8  WORK(MAXWRK)
      REAL *8  ERRR, ERRSCALE
      REAL *8  XMAX, DIFFMAX, DIFF
      REAL *8  ERRRX, ERRSCALEX
      REAL *8  XMAXX, DIFFMAXX, DIFFX
      REAL *8  ERRRY, ERRSCALEY
      REAL *8  XMAXY, DIFFMAXY, DIFFY
      REAL *8  COEFFSU(0:7,0:7,MAXBOXES)
      REAL *8  HH, EPS
      REAL *8  H2, HLEFT, HRIGHT, HBOTTOM, HTOP
      REAL *8  TIME0,TIME1,SECOND
      EXTERNAL H2, HLEFT, HRIGHT, HBOTTOM, HTOP

C-----Open output file:
      OPEN(25,FILE='outputclassic.m')
      OPEN(26,FILE='outdiff.m')
      OPEN(27,FILE='outputclassicx.m')
      OPEN(28,FILE='outputclassicy.m')
      OPEN(29,FILE='dot.m')
      OPEN(30,FILE='outputlap.m')
      OPEN(32,FILE='outputpractical.m')
      OPEN(33,FILE='outputpracticalx.m')
      OPEN(34,FILE='outputpracticaly.m')

C     Set the number of points in the box to 64.
C     Each box will contain 64 chebyshev points.
      NINBOX = 64

C     Now let's set the flag variable that determines whether
C     we are solving the periodic problem or the free space
C     problem, IPERIOD is set to distinguish between each
C     of the following cases:
C     IPERIOD = 0 for the free space case
C     IPERIOD = 1 for the periodic case
C     IPERIOD = 2 for the homogeneous dirichlet case
C     IPERIOD = 3 for the homogeneous neumann case
C     IPERIOD = 4 for the dirichlet left/right and periodic top/bottom
C     IPERIOD = 5 for the dirichlet left/right and neumann top/bottom
C     IPERIOD = 6 for the neumann left/right and periodic top/bottom
C     IPERIOD = 7 for the purely dirichlet (inhomogeneous)
C     IPERIOD = 8 for the purely neumann (inhomogeneous)
C     IPERIOD = 9 for the dirichlet left/right and periodic top/bottom
C     IPERIOD = 10 for the dirichlet left/right and neumann top/bottom
C     IPERIOD = 11 for the neumann left/right and periodic top/bottom

      IPERIOD = 7
      IPERIOD = 0

C     IPREC = 0  3  digits
C     IPREC = 1  6  digits
C     IPREC = 2  9  digits
C     IPREC = 3  12 digits

      IPREC = 3


C     Call an initialization routine in which
C     the tree is set to be an arbitrary tree.
      WRITE(*,*)'Setting tree'
      CALL SETTREE(LEVELBOX,ICOLBOX,IROWBOX,NBOXES,NLEV)

ccc      WRITE(*,*)'Setting a uniform tree'
ccc      NLEV = 3
ccc      CALL UNITREE(LEVELBOX,ICOLBOX,IROWBOX,NBOXES,NLEV,
ccc     1             ICHILDBOX, IPARENTBOX, LADDER)

      WRITE(*,*)'Generating parents and children'
      CALL MKCHILD(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1                              IROWBOX,NBOXES,NLEV)

C     Call an initialization routine in which
C     the tree is set to be an arbitrary tree.
      WRITE(*,*)'Making tree'
C     Set the error bound (it is scaled later):
C     (This bound is used to determine whether or not
C     further subdivision is needed by looking at the
C     error between the right hand side and an
C     approximating polynomial.)
ccc      EPS = 10.0D0**(-1)
ccc      EPS = 10.0D0**(-2)
ccc      EPS = 10.0D0**(-3)
      EPS = 10.0D0**(-4)
ccc      EPS = 10.0D0**(-5)
ccc      EPS = 10.0D0**(-6)
ccc      EPS = 10.0D0**(-7)
ccc      EPS = 10.0D0**(-8)
ccc      EPS = 10.0D0**(-9)
      EPS = 10.0D0**(-10)
      CALL MKTREE8(LEVELBOX,ICOLBOX,IROWBOX,NBOXES,NLEV,
     1          IPARENTBOX, ICHILDBOX, 
     2          NBLEVEL, IBOXLEV, ISTARTLEV,
     3          MAXBOXES, ITEMPARRAY, MAXLEVEL, EPS, H2)
      WRITE(*,*)'NBOXES = ',NBOXES


      IF(NBOXES .GT. MAXBOXES)THEN
       WRITE(*,*)'There are more boxes needed than are allotted'
       WRITE(*,*)'for in the memory space.'
       WRITE(*,*)'Readjust the memory.'
       WRITE(*,*)'I am stopping'
       STOP
      ENDIF

C     Call a routine to determine whether or not we need to call the
C     fixtree routine.  If the flag IFIXTREE is set to 1, FIXTREE
C     is called, if IFIXTREE is set to 0, FIXTREE does not need to
C     be called.
      IFIXFLAG = 0
      IF (IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN
         CALL RESTRICTION(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1                IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2                NBLEVEL,IBOXLEV,ISTARTLEV,IPERIOD,IFIXFLAG)
      ELSEIF (IPERIOD .EQ. 2 .OR. IPERIOD .EQ. 3
     1   .OR. IPERIOD .EQ. 5 .OR. IPERIOD .EQ. 7
     2   .OR. IPERIOD .EQ. 8 .OR. IPERIOD .EQ. 10)THEN
         IPERIODTEMP = 0
         CALL RESTRICTION(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1                IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2                NBLEVEL,IBOXLEV,ISTARTLEV,IPERIODTEMP,
     3                IFIXFLAG)
      ELSEIF (IPERIOD .EQ. 4 .OR. IPERIOD .EQ. 6 .OR.
     1        IPERIOD .EQ. 9 .OR. IPERIOD .EQ. 11)THEN
        IPERIODTEMP = 3
        CALL RESTRICTION(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1               IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2               NBLEVEL,IBOXLEV,ISTARTLEV,IPERIODTEMP,
     3               IFIXFLAG)
      ENDIF


      IF(IFIXFLAG .EQ. 1)THEN
C       The correction routine is used only in the adaptive case:
        WRITE(*,*)'Correcting tree'
        IF (IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN
          CALL FIXTREE(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1               IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2               NBLEVEL, IBOXLEV, ISTARTLEV,IPERIOD,
     3               IFLAG, MAXBOXES,ITEMPARRAY)
        ELSEIF (IPERIOD .EQ. 2 .OR. IPERIOD .EQ. 3
     1     .OR. IPERIOD .EQ. 5 .OR. IPERIOD .EQ. 7
     2     .OR. IPERIOD .EQ. 8 .OR. IPERIOD .EQ. 10)THEN
          IPERIODTEMP = 0
          CALL FIXTREE(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1               IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2               NBLEVEL, IBOXLEV, ISTARTLEV,IPERIODTEMP,
     3               IFLAG, MAXBOXES,ITEMPARRAY)
        ELSEIF (IPERIOD .EQ. 4 .OR. IPERIOD .EQ. 6 .OR.
     1          IPERIOD .EQ. 9 .OR. IPERIOD .EQ. 11)THEN
          IPERIODTEMP = 3
          CALL FIXTREE(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1               IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2               NBLEVEL, IBOXLEV, ISTARTLEV,IPERIODTEMP,
     3               IFLAG, MAXBOXES,ITEMPARRAY)
        ENDIF
      ENDIF



C     Test to see if the total number of boxes
C     is over the allotted memory space.
      IF(NBOXES .GT. MAXBOXES)THEN
       WRITE(*,*)'There are more boxes needed than are allotted'
       WRITE(*,*)'for in the memory space.'
       WRITE(*,*)'Readjust the memory.'
       WRITE(*,*)'I am stopping'
       STOP
      ENDIF
 
      WRITE(*,*)'NBOXES = ',NBOXES

C     Let's set the right hand side, FRIGHT:
      WRITE(*,*)'Initializing the right hand side'
      CALL SETF8(FRIGHT,
     1          ICOLBOX, IROWBOX, ICHILDBOX,
     2          NLEV, NBLEVEL, IBOXLEV, ISTARTLEV,H2)

      WRITE(*,*)'Precomputing tables'
      CALL COMPUTEMAPS8(MAP, IPREC, LENMAPS)


C     Call a routine that will divide up the necessary workspace
C     and then call the routine adaptive subroutine which 
C     contains the bulk of the algorithm.
      WRITE(*,*)'Doing FMM'
      TIME0 = SECOND()


      CALL FMMSTART8_WRAP(NLEV,LEVELBOX,IPARENTBOX,ICHILDBOX,
     1   ICOLBOX,IROWBOX,NBOXES,NBLEVEL,IBOXLEV,ISTARTLEV,
     2   FRIGHT,POT,IPREC)

c      CALL FMMSTART8(WORK,MAXWRK,IWORK,IWRK,
c     1   NLEV,LEVELBOX,IPARENTBOX,ICHILDBOX,
c     2   ICOLBOX,IROWBOX,NBOXES,NBLEVEL,IBOXLEV,
c     3   ISTARTLEV,IPERIOD,FRIGHT,POT,IPREC,MAP,
c     4   HLEFT,HRIGHT,HBOTTOM,HTOP)
      TIME1 = SECOND()



      CALL GETXYCLASSICAL8(XF, YF,
     1       ICOLBOX, IROWBOX, ICHILDBOX,NLEV,
     2       NBLEVEL, IBOXLEV, ISTARTLEV)

      CALL GETXYPRACTICAL8(XFPRACTICAL, YFPRACTICAL,
     1       ICOLBOX, IROWBOX, ICHILDBOX,NLEV,
     2       NBLEVEL, IBOXLEV, ISTARTLEV)

      WRITE(*,*)'Getting derivatives'
      CALL GETDERIVATIVES8(POT, 
     1      NLEV, ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2      POTX, POTY, LAP, COEFFSU)

      CALL GETPOTPRACTICAL(POT, POTX, POTY,
     1           POTPRACTICAL, POTPRACTICALX, 
     2           POTPRACTICALY,
     3           NLEV, ICHILDBOX, NBLEVEL,
     4           IBOXLEV, ISTARTLEV,COEFFSU)


C     Output all of the data to a matlab script so that the
C     solution can be viewed.  Only the data for childless
C     boxes is outputted.
C     The data for a childless box is written out as a four
C     by four grid with the appropriate x, y, and z coordinates.
      DO I = 0, NLEV
      DO 100 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100
       DO L = 1, 8
       WRITE(25,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(25,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(25,*)' ];'
       WRITE(25,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(25,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(25,*)' ];'
       WRITE(25,*)' pot(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(25,*)POT(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(25,*)' ];' 
      END DO
       WRITE(25,*)'mesh(x,y,pot)'
       WRITE(25,*)'hold on'
       WRITE(25,*)'axis(''square'')'
100   CONTINUE
      END DO


C     Now let's test the accuracy of our output.
C     For each of the childless boxes, use a gaussian
C     quadrature routine to calculate the exact integral
C     for each output point.
      IF(IPERIOD .EQ. 0)THEN
       WRITE(*,*)'Computing test values'
       DO I = 0, NLEV
        DO 200 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
        IBOX = IBOXLEV(J)
        WRITE(*,*)IBOX
        IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 200
         DO L = 1, 64
ccc           CALL MKSUR(XF(L,IBOX),YF(L,IBOX),TEST(L,IBOX))
           TEST(L,IBOX) = HEXACT(XF(L,IBOX),YF(L,IBOX))
         END DO
200     CONTINUE
       END DO
      ENDIF



C     Output the error:
      IF(IPERIOD .EQ. 0)THEN
      ERRR  =  0.0D0
      ERRSCALE  =  0.0D0
      DIFFMAX = 0.0D0
      XMAX = 0.0D0
      NCNTR  =  0
      DO I = 0, NLEV
      HH = 1.0D0 / DBLE(4**I)
      DO 300 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 300
       NCNTR = NCNTR + 1
       DO L = 1, 8
       WRITE(26,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(26,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(26,*)' ];'
       WRITE(26,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(26,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(26,*)' ];'
       WRITE(26,*)' pot(:,',L,') = ['
        DO JCNTR = 1, 8
          WRITE(26,*)DABS(TEST(8*(JCNTR-1)+L,IBOX) 
     1    -  POT(8*(JCNTR-1)+L,IBOX)) * DSQRT(HH)

C      First sum up the L2 error:
          ERRR = ERRR + (DABS(TEST(8*(JCNTR-1)+L,IBOX)
     1    - POT(8*(JCNTR-1)+L,IBOX)))**2 * HH

          ERRSCALE = ERRSCALE
     1    + DABS(TEST(8*(JCNTR-1)+L,IBOX))**2 * HH


C       Now find the L- infinity error:
          DIFF = TEST(8*(JCNTR-1)+L,IBOX)
     1    -  POT(8*(JCNTR-1)+L,IBOX)

          IF(DABS(DIFF) .GE. DABS(DIFFMAX))THEN
              DIFFMAX = DABS(DIFF)
          ENDIF

          IF(DABS(TEST(8*(JCNTR-1)+L,IBOX)) .GE. DABS(XMAX))THEN
             XMAX = DABS(TEST(8*(JCNTR-1)+L,IBOX))
          ENDIF

        END DO
       WRITE(26,*)' ];'
       END DO
       WRITE(26,*)'mesh(x,y,pot)'
       WRITE(26,*)'hold on'
       WRITE(26,*)'axis(''square'')'
300   CONTINUE
      END DO

      WRITE(*,*)'Number of childless boxes = ',NCNTR
      WRITE(*,*)'The L2 error is = ',DSQRT(ERRR)/DSQRT(ERRSCALE)
      WRITE(*,*)'The L infinity error is = ',DIFFMAX/XMAX

      ELSE
      ERRR  =  0.0D0
      ERRRX  =  0.0D0
      ERRRY  =  0.0D0
      ERRSCALE  =  0.0D0
      ERRSCALEX  =  0.0D0
      ERRSCALEY  =  0.0D0
      DIFFMAX = 0.0D0
      DIFFMAXX = 0.0D0
      DIFFMAXY = 0.0D0
      XMAX = 0.0D0
      XMAXX = 0.0D0
      XMAXY = 0.0D0
      NCNTR  =  0
      DO I = 0, NLEV
      HH = 1.0D0 / DBLE(4**I)
      DO 700 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 700
       NCNTR = NCNTR + 1
       DO L = 1, 8
       WRITE(26,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(26,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(26,*)' ];'
       WRITE(26,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(26,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(26,*)' ];'
       WRITE(26,*)' pot(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(26,*)POT(8*(JCNTR-1)+L,IBOX) -
     1    HEXACT(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX)) 

          ERRR = ERRR + 
     1  (DABS(HEXACT(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     2    - POT(8*(JCNTR-1)+L,IBOX)))**2 * HH

          ERRSCALE = ERRSCALE
     1 + DABS(HEXACT(XF(8*(JCNTR-1)+L,IBOX),
     2               YF(8*(JCNTR-1)+L,IBOX))**2) * HH


          ERRRX = ERRRX +
     1  (DABS(HEXACTX(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     2    - POTX(8*(JCNTR-1)+L,IBOX)))**2 * HH

          ERRSCALEX = ERRSCALEX
     1 + DABS(HEXACTX(XF(8*(JCNTR-1)+L,IBOX),
     2                YF(8*(JCNTR-1)+L,IBOX))**2) * HH


          ERRRY = ERRRY +
     1  (DABS(HEXACTY(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     2    - POTY(8*(JCNTR-1)+L,IBOX)))**2 * HH

          ERRSCALEY = ERRSCALEY
     1 + DABS(HEXACTY(XF(8*(JCNTR-1)+L,IBOX),
     2                YF(8*(JCNTR-1)+L,IBOX))**2) * HH


C       Now sum up the L- infinity error:
          DIFF = HEXACT(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     1    -  POT(8*(JCNTR-1)+L,IBOX)

          IF(DABS(DIFF) .GE. DABS(DIFFMAX))THEN
              DIFFMAX = DABS(DIFF)
          ENDIF

          IF(DABS(HEXACT(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX)))
     1              .GE. DABS(XMAX))THEN
             XMAX = DABS(HEXACT(XF(8*(JCNTR-1)+L,IBOX),
     1                         YF(8*(JCNTR-1)+L,IBOX)))
          ENDIF


          DIFFX = HEXACTX(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     1    -  POTX(8*(JCNTR-1)+L,IBOX)

          IF(DABS(DIFFX) .GE. DABS(DIFFMAXX))THEN
              DIFFMAXX = DABS(DIFFX)
          ENDIF

          IF(DABS(HEXACTX(XF(8*(JCNTR-1)+L,IBOX),
     1                    YF(8*(JCNTR-1)+L,IBOX)))
     2              .GE. DABS(XMAXX))THEN
             XMAXX = DABS(HEXACTX(XF(8*(JCNTR-1)+L,IBOX),
     1                         YF(8*(JCNTR-1)+L,IBOX)))
          ENDIF


          DIFFY = HEXACTY(XF(8*(JCNTR-1)+L,IBOX),YF(8*(JCNTR-1)+L,IBOX))
     1    -  POTY(8*(JCNTR-1)+L,IBOX)

          IF(DABS(DIFFY) .GE. DABS(DIFFMAXY))THEN
              DIFFMAXY = DABS(DIFFY)
          ENDIF

          IF(DABS(HEXACTY(XF(8*(JCNTR-1)+L,IBOX),
     1                    YF(8*(JCNTR-1)+L,IBOX)))
     2              .GE. DABS(XMAXY))THEN
             XMAXY = DABS(HEXACTY(XF(8*(JCNTR-1)+L,IBOX),
     1                         YF(8*(JCNTR-1)+L,IBOX)))
          ENDIF


        END DO
       WRITE(26,*)' ];'
      END DO
       WRITE(26,*)'mesh(x,y,pot)'
       WRITE(26,*)'hold on'
       WRITE(26,*)'axis(''square'')'
700   CONTINUE
      END DO
      WRITE(*,*)'Number of childless boxes = ',NCNTR
      WRITE(*,*)'The L2 error is = ',DSQRT(ERRR)/DSQRT(ERRSCALE)
      WRITE(*,*)'The L infinity error is = ',DIFFMAX/XMAX
      WRITE(*,*)'The L2 error in the x derivative is = ',
     1                    DSQRT(ERRRX)/DSQRT(ERRSCALEX)
      WRITE(*,*)'The L infinity error in the x derivative is = ',
     1                    DIFFMAXX/XMAXX
      WRITE(*,*)'The L2 error in the y derivative is = ',
     1                    DSQRT(ERRRY)/DSQRT(ERRSCALEY)
      WRITE(*,*)'The L infinity error in the y derivative is = ',
     1                    DIFFMAXY/XMAXY
      ENDIF

      WRITE(*,*)'The time for doing the fast multipole'
      WRITE(*,*)'calculation is:',TIME1-TIME0
      WRITE(*,*)'The time per grid point is:',
     1                 (TIME1-TIME0)/DBLE(NINBOX*NCNTR)



C     Output the x derivatives:
      DO I = 0, NLEV
      DO 400 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 400
       DO L = 1, 8
       WRITE(27,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(27,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(27,*)' ];'
       WRITE(27,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(27,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(27,*)' ];'
       WRITE(27,*)' potx(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(27,*)POTX(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(27,*)' ];'
       END DO
       WRITE(27,*)'mesh(x,y,potx)'
       WRITE(27,*)'hold on'
       WRITE(27,*)'axis(''square'')'
400   CONTINUE
      END DO

C     Output the y derivatives:
      DO I = 0, NLEV
      DO 500 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 500
       DO L = 1, 8
       WRITE(28,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(28,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(28,*)' ];'
       WRITE(28,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(28,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(28,*)' ];'
       WRITE(28,*)' poty(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(28,*)POTY(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(28,*)' ];'
       END DO
       WRITE(28,*)'mesh(x,y,poty)'
       WRITE(28,*)'hold on'
       WRITE(28,*)'axis(''square'')'
500   CONTINUE
      END DO

C     Output the laplacian:
      DO I = 0, NLEV
      DO 600 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 600
       DO L = 1, 8
       WRITE(30,*)' y(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(30,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(30,*)' ];'
       WRITE(30,*)' x(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(30,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(30,*)' ];'
       WRITE(30,*)' lap(:,',L,') = ['
        DO JCNTR = 1, 8
         WRITE(30,*)LAP(8*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(30,*)' ];'
       END DO
       WRITE(30,*)'mesh(x,y,lap)'
       WRITE(30,*)'hold on'
       WRITE(30,*)'axis(''square'')'
600   CONTINUE
      END DO


      WRITE(29,*)' xx = ['
      DO I = 0, NLEV
      DO 1200 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 1200
       DO L = 1, 8
        DO JCNTR = 1, 8
         WRITE(29,*)XF(8*(JCNTR-1)+L,IBOX)
        END DO
       END DO
1200  CONTINUE
      END DO
      WRITE(29,*)' ];'

      WRITE(29,*)
      WRITE(29,*)

      WRITE(29,*)' yy = ['
      DO I = 0, NLEV
      DO 1300 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 1300
       DO L = 1, 8
        DO JCNTR = 1, 8
         WRITE(29,*)YF(8*(JCNTR-1)+L,IBOX)
        END DO
       END DO
1300  CONTINUE
      END DO
      WRITE(29,*)' ];'
      WRITE(29,*)'plot(xx,yy,''.'')'
      WRITE(29,*)'axis(''square'')'

      DO I = 0, NLEV
      DO 1400 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 1400
       DO L = 1, 9
       WRITE(32,*)' ypractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(32,*)YFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(32,*)' ];'
       WRITE(32,*)' xpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(32,*)XFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(32,*)' ];'
       WRITE(32,*)' potpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(32,*)POTPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(32,*)' ];'
      END DO
      WRITE(32,*)'mesh(xpractical,ypractical,potpractical)'
      WRITE(32,*)'hold on'
      WRITE(32,*)'axis(''square'')'
1400  CONTINUE
      END DO

      DO I = 0, NLEV
      DO 1500 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 1500
       DO L = 1, 9
       WRITE(33,*)' ypractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(33,*)YFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(33,*)' ];'
       WRITE(33,*)' xpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(33,*)XFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(33,*)' ];'
       WRITE(33,*)' potpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(33,*)POTPRACTICALX(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(33,*)' ];'
      END DO
      WRITE(33,*)'mesh(xpractical,ypractical,potpractical)'
      WRITE(33,*)'hold on'
      WRITE(33,*)'axis(''square'')'
1500  CONTINUE
      END DO

      DO I = 0, NLEV
      DO 1600 J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
       IBOX = IBOXLEV(J)
       IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 1600
       DO L = 1, 9
       WRITE(34,*)' ypractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(34,*)YFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(34,*)' ];'
       WRITE(34,*)' xpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(34,*)XFPRACTICAL(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(34,*)' ];'
       WRITE(34,*)' potpractical(:,',L,') = ['
        DO JCNTR = 1, 9
         WRITE(34,*)POTPRACTICALY(9*(JCNTR-1)+L,IBOX)
        END DO
       WRITE(34,*)' ];'
      END DO
      WRITE(34,*)'mesh(xpractical,ypractical,potpractical)'
      WRITE(34,*)'hold on'
      WRITE(34,*)'axis(''square'')'
1600  CONTINUE
      END DO


      STOP
      END



      SUBROUTINE FMMSTART8_WRAP(NLEV,LEVELBOX,IPARENTBOX,ICHILDBOX,
     1   ICOLBOX,IROWBOX,NBOXES,NBLEVEL,IBOXLEV,ISTARTLEV,
     2   FRIGHT,POT,IPREC)
      implicit real *8 (a-h,o-z)
      INTEGER NLEV,NBOXES
      INTEGER LEVELBOX(NBOXES)
      INTEGER IPARENTBOX(NBOXES)
      INTEGER ICHILDBOX(4,NBOXES)
      INTEGER ICOLBOX(NBOXES)
      INTEGER IROWBOX(NBOXES)
      INTEGER NBLEVEL(0:NLEV)
      INTEGER IBOXLEV(NBOXES)
      INTEGER ISTARTLEV(0:NLEV)
      REAL *8 FRIGHT(64,NBOXES)
      REAL *8 POT(64,NBOXES)
      REAL *8, allocatable :: WORK(:)
      INTEGER, allocatable :: IWORK(:)
      REAL *8, allocatable :: MAPS(:)
      EXTERNAL HLEFT,HRIGHT,HBOTTOM,HTOP

      LENMAPS = 30000
      allocate(MAPS(LENMAPS))
      CALL COMPUTEMAPS8(MAPS, IPREC, LENMAPS)
      
      MAXWRK = 1000*NBOXES*64
      IWRK = 1000*NBOXES*64

      allocate(WORK(MAXWRK),IWORK(IWRK))

      IPERIOD = 0
      
      

      CALL FMMSTART8(WORK,MAXWRK,IWORK,IWRK,
     1   NLEV,LEVELBOX,IPARENTBOX,ICHILDBOX,
     2   ICOLBOX,IROWBOX,NBOXES,NBLEVEL,IBOXLEV,
     3   ISTARTLEV,IPERIOD,FRIGHT,POT,IPREC,MAPS,
     4   HLEFT,HRIGHT,HBOTTOM,HTOP)
      RETURN
      END


C**********************************************************************
C     The following subroutine is just set up to initialize the
C     tree structure.  For each box, the level, parent, row, column,
C     and children must be specified.  If any of these things are 
C     not present, then that quantity is set to -1.
C     (This is used as a flag in later routines.)
C     Nothing is set upon input to this routine.  All of the above
C     arrays are just set in the subroutine below.
C
C
C     INPUT:
C     
C     Nothing is defined on input as this is an initialization routine.
C
C
C     OUTPUT:
C
C     LEVELBOX is an array defining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C**********************************************************************
      SUBROUTINE SETTREE(LEVELBOX,ICOLBOX,IROWBOX,NBOXES,NLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER LEVELBOX(1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER NBOXES, NLEV

C-----Local variables
      INTEGER IBOX

C     Check to see if we are in the adaptive
C     or the uniform case:

ccc      GOTO 400
C     Mock tree number one:
      IBOX = 1
      LEVELBOX(IBOX) = 0
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 2
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 3
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1

      IBOX = 4
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 2

      IBOX = 5
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 2

      IBOX = 6
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 7
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1

      IBOX = 8
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 2

      IBOX = 9
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 2

      NBOXES = 9
      NLEV = 2
      RETURN

300   CONTINUE
      IBOX = 1
      LEVELBOX(IBOX) = 0
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1
 
      IBOX = 2
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 3
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1
 
      IBOX = 4
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 2
 
      IBOX = 5
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 2

      IBOX = 6
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 7
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1
 
      IBOX = 8
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 3
      IROWBOX(IBOX) = 1

      IBOX = 9
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 4
      IROWBOX(IBOX) = 1

      NBOXES = 9
      NLEV = 2
      RETURN

400   CONTINUE
      IBOX = 1
      LEVELBOX(IBOX) = 0
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 2
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 3
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1

      IBOX = 4
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 2

      IBOX = 5
      LEVELBOX(IBOX) = 1
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 2

      IBOX = 6
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1

      IBOX = 7
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 1

      IBOX = 8
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 2

      IBOX = 9
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 2

      IBOX = 10
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 3

      IBOX = 11
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 3

      IBOX = 12
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 4

      IBOX = 13
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 2
      IROWBOX(IBOX) = 4

      IBOX = 14
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 3
      IROWBOX(IBOX) = 1

      IBOX = 15
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 4
      IROWBOX(IBOX) = 1

      IBOX = 16
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 3
      IROWBOX(IBOX) = 2

      IBOX = 17
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 4
      IROWBOX(IBOX) = 2

      IBOX = 18
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 3
      IROWBOX(IBOX) = 3

      IBOX = 19
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 4
      IROWBOX(IBOX) = 3

      IBOX = 20
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 3
      IROWBOX(IBOX) = 4

      IBOX = 21
      LEVELBOX(IBOX) = 2
      ICOLBOX(IBOX) = 4
      IROWBOX(IBOX) = 4

      IBOX = 22
      LEVELBOX(IBOX) = 3
      ICOLBOX(IBOX) = 5
      IROWBOX(IBOX) = 7

      IBOX = 23
      LEVELBOX(IBOX) = 3
      ICOLBOX(IBOX) = 5
      IROWBOX(IBOX) = 8

      IBOX = 24
      LEVELBOX(IBOX) = 3
      ICOLBOX(IBOX) = 6
      IROWBOX(IBOX) = 7

      IBOX = 25
      LEVELBOX(IBOX) = 3
      ICOLBOX(IBOX) = 6
      IROWBOX(IBOX) = 8

      IBOX = 26
      LEVELBOX(IBOX) = 4
      ICOLBOX(IBOX) = 9
      IROWBOX(IBOX) = 15

      IBOX = 27
      LEVELBOX(IBOX) = 4
      ICOLBOX(IBOX) = 10
      IROWBOX(IBOX) = 15

      IBOX = 28
      LEVELBOX(IBOX) = 4
      ICOLBOX(IBOX) = 9
      IROWBOX(IBOX) = 16

      IBOX = 29
      LEVELBOX(IBOX) = 4
      ICOLBOX(IBOX) = 10
      IROWBOX(IBOX) = 16


      NBOXES = 29
      NLEV = 4
      RETURN
      END


C**********************************************************************
C     The following subroutine is used to generate a uniform tree.
C     This is only for testing purposes.
C
C
C     INPUT:
C
C     NLEV is the finest level
C
C     OUTPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C
C**********************************************************************
      SUBROUTINE UNITREE(LEVELBOX,ICOLBOX,IROWBOX,NBOXES,NLEV, 
     1                   ICHILDBOX, IPARENTBOX, LADDER)
      IMPLICIT NONE
C-----Global variables
      INTEGER  LEVELBOX(1)
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  NBOXES, NLEV
      INTEGER  ICHILDBOX(4,1), IPARENTBOX(1)
      INTEGER  LADDER(0:1)
C-----Local variables
      INTEGER  I, IBOX, J, K, L, NSIDE
      INTEGER  ICOL, IROW
      INTEGER  IC1, IC2, IC3, IC4
      INTEGER  ILENGTH,  ILEV, ISTART, IDIM
      INTEGER  MYBOX, ISTARTC

C     Just create a uniform tree, given NLEV as the input:
C     (In this loop, let's initialize the parents
C     and children arrays to be -1.)
      IBOX = 1
      NSIDE = 1
      DO I = 0, NLEV
        DO J = 1, NSIDE
          DO K = 1, NSIDE
            LEVELBOX(IBOX) = I
            ICOLBOX(IBOX)  = J
            IROWBOX(IBOX)  = K
            IPARENTBOX(IBOX) = -1
            DO L = 1, 4
              ICHILDBOX(L,IBOX) = -1
            END DO
            IBOX = IBOX + 1
          END DO
        END DO
       NSIDE = 2 * NSIDE
      END DO
      NBOXES = IBOX - 1

C     Now initialize the LADDER array:
      J = 1
      LADDER(0) = 0
      LADDER(1) = 1
      DO I = 2, NLEV + 1
        J = 4*J
        LADDER(I) = LADDER(I-1) +  J
      END DO


C     Now set all of the parents and children:
      DO I = 0, NLEV - 1
        ISTART = LADDER(NLEV-I-1)
        ILENGTH = LADDER(NLEV-I) - ISTART
        ILEV = NLEV-I-1
        ISTARTC = LADDER(NLEV-I)
        IDIM = 2**ILEV
        NSIDE = 2*IDIM

        DO J = 1, ILENGTH
          ICOL = 1 + MOD(J-1,IDIM )
          IROW = 1 + (J-1)/IDIM

          MYBOX = ISTART+J
          IC1 = ISTARTC + (2*ICOL) + (2*IROW-2)*NSIDE
          IC2 = ISTARTC + (2*ICOL) + (2*IROW-1)*NSIDE
          IC3 = ISTARTC + (2*ICOL-1) + (2*IROW-1)*NSIDE
          IC4 = ISTARTC + (2*ICOL-1) + (2*IROW-2)*NSIDE

          ICHILDBOX(1,MYBOX) = IC1
          ICHILDBOX(2,MYBOX) = IC2
          ICHILDBOX(3,MYBOX) = IC3
          ICHILDBOX(4,MYBOX) = IC4

          IPARENTBOX(IC1) = MYBOX
          IPARENTBOX(IC2) = MYBOX
          IPARENTBOX(IC3) = MYBOX
          IPARENTBOX(IC4) = MYBOX

        END DO
      END DO

      RETURN
      END


C***************************************************************************
C     The following subroutine is used to generate a tree given a
C     right hand side, FRIGHT.  
C     The algorithm works by testing the function values vs. an approximating 
C     polynomial and dividing if necessary.  It is important to note that
C     the tree generated by this algorithm may not satisfy the level 
C     restriction.  It may be necessary to call the routine FIXTREE 
C     after this to make sure that the tree can be handled properly by the
C     ADAPFMM6  and BOUNDFMM6 routines.
C
C
C     INPUT:
C
C     MAXBOXES denotes the maximum number of boxes allowed
C
C     MAXLEVEL denotes the deepest level allowed
C
C     ITEMPARRAY is just a dummy array
C
C     EPS denotes the error bound that determines 
C         when refinement take place
C
C     H is the real function that is the right hand side
C       of the poisson equation
C
C     OUTPUT:
C
C     ISTARTLEV is the pointer to where each level
C               begins in the IBOXLEV array
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C***************************************************************************
      SUBROUTINE MKTREE8(LEVELBOX, ICOLBOX, IROWBOX, NBOXES, NLEV,
     1           IPARENTBOX, ICHILDBOX,
     2           NBLEVEL, IBOXLEV, ISTARTLEV,
     3           MAXBOXES, ITEMPARRAY, MAXLEVEL, EPS, H)
C-----Global variables
      IMPLICIT NONE
      INTEGER  LEVELBOX(1), MAXBOXES
      INTEGER  NLEV, NBOXES,  MAXLEVEL
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  ITEMPARRAY(1)
      REAL *8  EPS, H
C-----Local variables
      INTEGER  I, IBOX, IFLAG
      INTEGER  J, ISTART, IEND, JJ
      INTEGER  LEVFLAG, II
      REAL *8  COEFFTEMP(8,8)
      REAL *8  EPSCALED
      REAL *8  ERROR, HH
      REAL *8  XF(8), YF(8)
      REAL *8  FTEMP(8,8)
      REAL *8  WSAVE2(1000)

      DO I = 0, MAXLEVEL
        NBLEVEL(I) = 0
        ISTARTLEV(I) = 0
      END DO
      DO I = 1, MAXBOXES
        IBOXLEV(I) = 0
      END DO

C     First set the big parent box to the 
C     appropriate settings:
C     (Initially, there is just one box and
C     it is the big parent box at level 0)
      IBOX = 1
      LEVELBOX(IBOX) = 0
      ICOLBOX(IBOX) = 1
      IROWBOX(IBOX) = 1
      IPARENTBOX(IBOX) = -1
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1

      NBOXES = 1
      NLEV = 0

C     We also need to initialize the adaptive 'Ladder'
C     structures to the correct initial values:
      NBLEVEL(0) = 1
      ISTARTLEV(0) = 1
      IBOXLEV(1) = 1

      CALL CHXCIN(8,WSAVE2)
      DO I = 0, MAXLEVEL - 1
      IFLAG = 0
      LEVFLAG = 0
      ISTART = ISTARTLEV(I)
      IEND = ISTART + NBLEVEL(I) - 1
      DO J = ISTART, IEND
       IBOX = IBOXLEV(J)

       CALL MKGRID(XF,YF,ICOLBOX(IBOX),IROWBOX(IBOX),LEVELBOX(IBOX))

       DO JJ = 1, 8
         DO II = 1, 8
           FTEMP(JJ,II) = H(XF(JJ),YF(II))
         END DO
       END DO


C      compute Chebyshev transforms
       CALL GETCOEFF(8,8,FTEMP,COEFFTEMP,WSAVE2)

       CALL GETERROR(COEFFTEMP,ERROR)

       HH = DBLE(4**LEVELBOX(IBOX))
       EPSCALED = EPS * HH

       IF(ERROR .GE. EPSCALED)THEN
C        Call subdivide
         CALL SUBDIVIDE1(IBOX,IPARENTBOX,ICHILDBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,MAXLEVEL,
     2         ISTARTLEV, NBLEVEL, IBOXLEV,ITEMPARRAY)
         IF(LEVFLAG .EQ. 0)THEN
           NLEV = NLEV + 1
           LEVFLAG = 1
         ENDIF
         IFLAG = 1
       ENDIF
      END DO
      IF(IFLAG .EQ. 0)THEN
C      Nothing was divided at the
C      last level, so exit the loop.
       RETURN
      ENDIF
      END DO
      RETURN
      END


C***************************************************************************
C     This subroutine is only called within the MKTREE6 subroutine. 
C     This subroutine generates two grids: (XF,YF) and (XF2,YF2). 
C     They represent the 36 point grid that the polynomial 
C     coefficients are obtained from and the 144  point grid (the children's
C     grids) that is used to test the accuracy of the polynomial
C     approximation.  (Both grids are face centered.) 
C
C     INPUT:
C
C     ICOL denotes the column of the given box
C
C     IROW denotes the row of the given box
C
C     LEVEL denotes the level of the given box
C
C     EPS denotes the error bound that determines
C         when refinement take place
C
C     OUTPUT:
C
C     XF denotes the X values of 36 cell centered points in the box
C
C     YF denotes the Y values of 36 cell centered points in the box
C
C     XF2 denotes the X values of 144 cell centered points in the box
C
C     YF2 denotes the Y values of 144 cell centered points in the box
C
C***************************************************************************
      SUBROUTINE MKGRID(XF,YF,ICOL,IROW,LEVEL) 
      IMPLICIT NONE
C-----Global variables
      INTEGER  LEVEL
      INTEGER  ICOL, IROW
      REAL *8  XF(8), YF(8)
C-----Local variables
      INTEGER  J
      REAL *8  XSTART
      REAL *8  XSHIFT, YSHIFT
      REAL *8  TEMP1, PI16
      REAL *8  XX(8), XSCALE(8)

       PI16 = DACOS(-1.0D0) / 16.0D0
       XX(1) = DCOS(15.0D0*PI16) / 2.0D0
       XX(2) = DCOS(13.0D0*PI16) / 2.0D0
       XX(3) = DCOS(11.0D0*PI16) / 2.0D0
       XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
       XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
       XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
       XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
       XX(8) = DCOS( 1.0D0*PI16) / 2.0D0

       TEMP1 = 1.0D0 / DBLE(2**LEVEL)
       XSTART = 0.5D0 * (1.0D0 - TEMP1)

       XSCALE(1) = XX(1) * TEMP1 - XSTART
       XSCALE(2) = XX(2) * TEMP1 - XSTART
       XSCALE(3) = XX(3) * TEMP1 - XSTART
       XSCALE(4) = XX(4) * TEMP1 - XSTART
       XSCALE(5) = XX(5) * TEMP1 - XSTART
       XSCALE(6) = XX(6) * TEMP1 - XSTART
       XSCALE(7) = XX(7) * TEMP1 - XSTART
       XSCALE(8) = XX(8) * TEMP1 - XSTART


       XSHIFT = DBLE(ICOL - 1) * TEMP1
       YSHIFT = DBLE(IROW - 1) * TEMP1

       DO J = 1, 8
        XF(J) = XSCALE(J) + XSHIFT
        YF(J) = XSCALE(J) + YSHIFT
       END DO

      RETURN
      END


C***************************************************************************
C     The following subroutine calculates the L2 error between two vectors
C     that is needed in order to generate the tree structure.  This is used 
C     to determine how accurate the approximating polynomial is.
C
C     INPUT:
C
C     XF denotes the X values of 36 cell centered points in the box
C
C     YF denotes the Y values of 36 cell centered points in the box
C
C     XF2 denotes the X values of 144 cell centered points in the box
C
C     YF2 denotes the Y values of 144 cell centered points in the box
C
C     LEVEL denotes the level of the given box
C
C     A is the matrix that maps 36 function values
C       onto the 21 basis function coefficients
C
C     OUTPUT:
C
C     ERROR is the L2 error between the exact solution evaluated
C     at the 144 points and the approximation given by the polynomial
C     (coefficients found from the 36 points)
C
C***************************************************************************
      SUBROUTINE GETERROR(COEFFTEMP,ERROR)
      IMPLICIT NONE
      REAL *8  ERROR, COEFFTEMP(8,8)
C-----Local variables

C      Just sum over the antidiagonal: (all 7th degree coefficients)
        ERROR = 0.0D0
        ERROR = ERROR + DABS(COEFFTEMP(8,1))
        ERROR = ERROR + DABS(COEFFTEMP(7,2))
        ERROR = ERROR + DABS(COEFFTEMP(6,3))
        ERROR = ERROR + DABS(COEFFTEMP(5,4))
        ERROR = ERROR + DABS(COEFFTEMP(4,5))
        ERROR = ERROR + DABS(COEFFTEMP(3,6))
        ERROR = ERROR + DABS(COEFFTEMP(2,7))
        ERROR = ERROR + DABS(COEFFTEMP(1,8))

C      With a weighted version of the off diagonal terms
        ERROR = ERROR + DABS(COEFFTEMP(7,1))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(6,2))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(5,3))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(4,4))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(3,5))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(2,6))/4.0D0
        ERROR = ERROR + DABS(COEFFTEMP(1,7))/4.0D0

      RETURN
      END


C**********************************************************************
C     The following subroutine is used to generate the colleagues
C     for all of the boxes in the tree structure.  If a colleague
C     doesn't exist it is set to -1.  Each box has nine colleagues
C     and they are ordered as follows:
C
C                        7     8     9
C               
C                        4     5     6
C
C                        1     2     3
C
C     You are your own colleague number 5.
C     The algorithm used here is recursive and takes advantage of
C     the fact that your colleagues can only be the children of 
C     your parents colleagues.  There is no need to scan all of the
C     boxes.  IPERIOD denotes whether or not we are in a periodic
C     or free space case.  The basic algorithm is the same, but in the
C     periodic case we have to look for boxes that are 'outside' of
C     the standard size box.
C
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     IPERIOD denotes what kind of colleagues are to be generated
C             IPERIOD = 0 : free space
C             IPERIOD = 1 or 2 : periodic
C             IPERIOD = 3 : periodic up/down and free space left/right
C
C     OUTPUT:
C
C     ICOLLEAGBOX denotes the colleagues of a given box
C
C**********************************************************************
      SUBROUTINE MKCOLLS(ICOLBOX,
     1      IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2      IPARENTBOX,ICHILDBOX,NBLEVEL,
     3      IBOXLEV,ISTARTLEV,IPERIOD)
      IMPLICIT NONE
C-----Global variables
      INTEGER ICOLLEAGBOX(9,1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER NBOXES, NLEV, IPARENTBOX(1)
      INTEGER ICHILDBOX(4,1)
      INTEGER NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER IPERIOD
C-----Local variables
      INTEGER COLLEAGUE, PARTEMP
      INTEGER JCNTR, IBOX, ITEST
      INTEGER ICNTR, ILEV, J, L, NSIDE
      INTEGER IROWTEMP, ICOLTEMP
      INTEGER IROWTEST, ICOLTEST


C     Initialize colleague number 5 to
C     yourself and all other colleagues to
C     -1.  -1 is the flag for the case when
C     the colleagues don't exist.  It can 
C     be overwritten below. 
      DO IBOX = 1, NBOXES
       ICOLLEAGBOX(5,IBOX) = IBOX
       DO J = 1, 4
         ICOLLEAGBOX(J,IBOX) = -1
       END DO
       DO J = 6, 9
         ICOLLEAGBOX(J,IBOX) = -1
       END DO
      END DO


C     Scan through all of the levels except the coarsest level.
C     The one box at this level cannot have any colleagues.
C     Do the uniform case first:
      IF(IPERIOD .EQ. 0)THEN
      DO ILEV = 1, NLEV
C      Scan through all of the boxes on each level.  For each test
C      box, scan the parent's colleagues and test to see if 
C      their children are in contact with the box being tested.
C      Each colleague is placed in the correct order, so there is
C      no need to 'shuffle' them later on.
       DO L = ISTARTLEV(ILEV), ISTARTLEV(ILEV) + NBLEVEL(ILEV) - 1
         IBOX    = IBOXLEV(L)
         PARTEMP = IPARENTBOX(IBOX)

C        IROWTEMP and ICOLTEMP denote the row and column of
C        the test box.
         IROWTEMP = IROWBOX(IBOX)
         ICOLTEMP = ICOLBOX(IBOX)
      
         DO 100 JCNTR = 1, 9
C          COLLEAGUE denotes the colleague of the parent box.
           COLLEAGUE = ICOLLEAGBOX(JCNTR,PARTEMP)
C          If the colleague doesn't exist
C          or is childless, skip it:
           IF (COLLEAGUE .LT. 0)GOTO 100
           IF (ICHILDBOX(1,COLLEAGUE) .LT. 0)GOTO 100
           DO ICNTR = 1, 4
             J = ICHILDBOX(ICNTR,COLLEAGUE)
C            IROWTEST and ICOLTEST denote the row and column of
C            the box being compared to the test box.
             IROWTEST = IROWBOX(J)
             ICOLTEST = ICOLBOX(J)

             IF(IROWTEMP .EQ. IROWTEST+1)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(1,IBOX) = J
               ELSEIF(ICOLTEMP .EQ. ICOLTEST)THEN
                 ICOLLEAGBOX(2,IBOX) = J
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(3,IBOX) = J
               ENDIF
             ELSEIF(IROWTEMP .EQ. IROWTEST)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(4,IBOX) = J
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(6,IBOX) = J
               ENDIF
             ELSEIF(IROWTEMP .EQ. IROWTEST-1)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(7,IBOX) = J
               ELSEIF(ICOLTEMP .EQ. ICOLTEST)THEN
                 ICOLLEAGBOX(8,IBOX) = J
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(9,IBOX) = J
               ENDIF
             ENDIF
          END DO 
100      CONTINUE
       END DO
      END DO 

C     Now compute the colleagues in
C     the periodic case, if necessary:
      ELSEIF(IPERIOD .EQ. 1 .OR. IPERIOD .EQ. 2)THEN
C     Initialize the first box (level 0) so
C     that it has its own colleagues.
C     This is necessary, because at deeper 
C     levels, the algorithm works by scanning 
C     the parent boxes colleagues.
      IBOX = IBOXLEV(ISTARTLEV(0))
      ICOLLEAGBOX(1,IBOX) = IBOX
      ICOLLEAGBOX(2,IBOX) = IBOX
      ICOLLEAGBOX(3,IBOX) = IBOX
      ICOLLEAGBOX(4,IBOX) = IBOX
      ICOLLEAGBOX(5,IBOX) = IBOX
      ICOLLEAGBOX(6,IBOX) = IBOX
      ICOLLEAGBOX(7,IBOX) = IBOX
      ICOLLEAGBOX(8,IBOX) = IBOX
      ICOLLEAGBOX(9,IBOX) = IBOX

      DO ILEV = 1, NLEV
      NSIDE = 2**ILEV
C      Scan through all of the boxes on each level.  For each test
C      box, scan the parent's colleagues and test to see if
C      their children are in contact with the box being tested.
C      Each colleague is placed in the correct order, so there is
C      no need to 'shuffle' them later on.
       DO L = ISTARTLEV(ILEV), ISTARTLEV(ILEV) + NBLEVEL(ILEV) - 1
        IBOX = IBOXLEV(L)

C       IROWTEMP and ICOLTEMP denote the
C       row and column of the test box.
        IROWTEMP = IROWBOX(IBOX)
        ICOLTEMP = ICOLBOX(IBOX)

C       IROWTEST and ICOLTEST denote the row and column of
C       the box being compared to the test box.

        DO 300 JCNTR = 1, 9
C         First determine the column and row numbers
C         of all of the potential colleagues:
          IF(JCNTR .EQ. 5)GOTO 300
          IF(JCNTR .EQ. 1)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 2)THEN
            ICOLTEST = ICOLTEMP
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 3)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 4)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP
          ELSEIF(JCNTR .EQ. 6)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP
          ELSEIF(JCNTR .EQ. 7)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP + 1
          ELSEIF(JCNTR .EQ. 8)THEN
            ICOLTEST = ICOLTEMP
            IROWTEST = IROWTEMP + 1
          ELSEIF(JCNTR .EQ. 9)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP + 1
          ENDIF

C         Now test to see if the test parameters 
C         lie in the domain:
C         (If they are outside of the domain, just
C         add or subtract the appropriate number so
C         that the boxes 'wrap around.')
          IF(ICOLTEST .LT. 1)THEN
            ICOLTEST = ICOLTEST + NSIDE
          ELSEIF(ICOLTEST .GT. NSIDE)THEN
            ICOLTEST = ICOLTEST - NSIDE
          ENDIF
          IF(IROWTEST .LT. 1)THEN
            IROWTEST = IROWTEST + NSIDE
          ELSEIF(IROWTEST .GT. NSIDE)THEN
            IROWTEST = IROWTEST - NSIDE
          ENDIF


       DO 200 J = 1, 9
        IF(ICOLLEAGBOX(J,IPARENTBOX(IBOX)) .LT. 0)GOTO 200
        IF(ICHILDBOX(1,ICOLLEAGBOX(J,IPARENTBOX(IBOX))) .LT. 0)GOTO 200
          DO ICNTR = 1, 4
            ITEST = ICHILDBOX(ICNTR,ICOLLEAGBOX(J,IPARENTBOX(IBOX)))
            IF(IROWBOX(ITEST) .EQ. IROWTEST
     1         .AND. ICOLBOX(ITEST) .EQ. ICOLTEST)THEN
               ICOLLEAGBOX(JCNTR,IBOX) = ITEST
            ENDIF
          END DO
200    CONTINUE
300    CONTINUE
      END DO
      END DO
      ELSEIF(IPERIOD .EQ. 3)THEN
C     Initialize the first box (level 0) so
C     that it has its own colleagues.
C     This is necessary, because at deeper 
C     levels, the algorithm works by scanning 
C     the parent boxes colleagues.
      IBOX = IBOXLEV(ISTARTLEV(0))
      ICOLLEAGBOX(1,IBOX) = IBOX
      ICOLLEAGBOX(2,IBOX) = IBOX
      ICOLLEAGBOX(3,IBOX) = IBOX
      ICOLLEAGBOX(4,IBOX) = IBOX
      ICOLLEAGBOX(5,IBOX) = IBOX
      ICOLLEAGBOX(6,IBOX) = IBOX
      ICOLLEAGBOX(7,IBOX) = IBOX
      ICOLLEAGBOX(8,IBOX) = IBOX
      ICOLLEAGBOX(9,IBOX) = IBOX

      DO ILEV = 1, NLEV
      NSIDE = 2**ILEV
C      Scan through all of the boxes on each level.  For each test
C      box, scan the parent's colleagues and test to see if
C      their children are in contact with the box being tested.
C      Each colleague is placed in the correct order, so there is
C      no need to 'shuffle' them later on.
       DO L = ISTARTLEV(ILEV), ISTARTLEV(ILEV) + NBLEVEL(ILEV) - 1
        IBOX = IBOXLEV(L)

C       IROWTEMP and ICOLTEMP denote the
C       row and column of the test box.
        IROWTEMP = IROWBOX(IBOX)
        ICOLTEMP = ICOLBOX(IBOX)

C       IROWTEST and ICOLTEST denote the row and column of
C       the box being compared to the test box.

        DO 500 JCNTR = 1, 9
C         First determine the column and row numbers
C         of all of the potential colleagues:
          IF(JCNTR .EQ. 5)GOTO 500
          IF(JCNTR .EQ. 1)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 2)THEN
            ICOLTEST = ICOLTEMP
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 3)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP - 1
          ELSEIF(JCNTR .EQ. 4)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP
          ELSEIF(JCNTR .EQ. 6)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP
          ELSEIF(JCNTR .EQ. 7)THEN
            ICOLTEST = ICOLTEMP - 1
            IROWTEST = IROWTEMP + 1
          ELSEIF(JCNTR .EQ. 8)THEN
            ICOLTEST = ICOLTEMP
            IROWTEST = IROWTEMP + 1
          ELSEIF(JCNTR .EQ. 9)THEN
            ICOLTEST = ICOLTEMP + 1
            IROWTEST = IROWTEMP + 1
          ENDIF

C         Now test to see if the test parameters 
C         lie in the domain:
C         (If they are outside of the domain, just
C         add or subtract the appropriate number so
C         that the boxes 'wrap around.')
          IF(IROWTEST .LT. 1)THEN
            IROWTEST = IROWTEST + NSIDE
          ELSEIF(IROWTEST .GT. NSIDE)THEN
            IROWTEST = IROWTEST - NSIDE
          ENDIF


       DO 400 J = 1, 9
        IF(ICOLLEAGBOX(J,IPARENTBOX(IBOX)) .LT. 0)GOTO 400
        IF(ICHILDBOX(1,ICOLLEAGBOX(J,IPARENTBOX(IBOX))) .LT. 0)GOTO 400
          DO ICNTR = 1, 4
            ITEST = ICHILDBOX(ICNTR,ICOLLEAGBOX(J,IPARENTBOX(IBOX)))
            IF(IROWBOX(ITEST) .EQ. IROWTEST
     1         .AND. ICOLBOX(ITEST) .EQ. ICOLTEST)THEN
               ICOLLEAGBOX(JCNTR,IBOX) = ITEST
            ENDIF
          END DO
400    CONTINUE
500    CONTINUE
      END DO
      END DO
      ENDIF
      RETURN
      END



C*************************************************************************
C
C     This subroutine will determine whether or not a given tree satisfies 
C     the level restriction.  If it doesn't, call FIXTREE to fix the
C     tree.
C
C*************************************************************************
      SUBROUTINE RESTRICTION(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX, 
     1             IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2             NBLEVEL,IBOXLEV,ISTARTLEV,IPERIOD,IFIXFLAG)
      IMPLICIT NONE
C-----Global variables
      INTEGER LEVELBOX(1), ICOLLEAGBOX(9,1)
      INTEGER IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER NBOXES, NLEV
      INTEGER NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER IPERIOD, IFIXFLAG
C-----Local variables
      INTEGER ICHILD(4),ICOLL(9), IBOX
      INTEGER I, IPAR, ITEST, J, NB
      INTEGER ITEMP


C     Let's sort all of the boxes by level.
C     This takes the place of the LADDER structure
C     in the uniform case.  All boxes are sorted
C     into the array and placed in their proper places.
      CALL SORTBOXES(LEVELBOX,NBOXES,NLEV,
     1           NBLEVEL,IBOXLEV,ISTARTLEV)


C     First let's call a subroutine that will
C     generate all of the colleagues for each
C     box.  The colleagues are generated in the
C     correct order so there is no need to 'shuffle'
C     them later on.
      CALL MKCOLLS(ICOLBOX,
     1       IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2       IPARENTBOX,ICHILDBOX,NBLEVEL,
     3       IBOXLEV, ISTARTLEV,IPERIOD)


      DO I = NLEV, 2, -1
        DO J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
          IBOX = IBOXLEV(J)
          IPAR  = IPARENTBOX(IBOX)
          ITEST = IPARENTBOX(IPAR)

            ICOLL(1) = ICOLLEAGBOX(1,ITEST)
            ICOLL(2) = ICOLLEAGBOX(2,ITEST)
            ICOLL(3) = ICOLLEAGBOX(3,ITEST)
            ICOLL(4) = ICOLLEAGBOX(4,ITEST)
            ICOLL(5) = ICOLLEAGBOX(5,ITEST)
            ICOLL(6) = ICOLLEAGBOX(6,ITEST)
            ICOLL(7) = ICOLLEAGBOX(7,ITEST)
            ICOLL(8) = ICOLLEAGBOX(8,ITEST)
            ICOLL(9) = ICOLLEAGBOX(9,ITEST)

            ICHILD(1) = ICHILDBOX(1,ITEST)
            ICHILD(2) = ICHILDBOX(2,ITEST)
            ICHILD(3) = ICHILDBOX(3,ITEST)
            ICHILD(4) = ICHILDBOX(4,ITEST)


          DO NB = 1, 9
            ITEMP = ICOLL(NB)
            IF(ICHILDBOX(1,ITEMP) .LT. 0)THEN
C             The neighboring box is not divided
C             we could have problems.
              IF (NB .EQ. 1)THEN
                IF(IPAR .EQ. ICHILD(4))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 2)THEN
                IF(IPAR .EQ. ICHILD(3) .OR. IPAR .EQ. ICHILD(4))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 3)THEN
                IF(IPAR .EQ. ICHILD(3))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 4)THEN
                IF(IPAR .EQ. ICHILD(4) .OR. IPAR .EQ. ICHILD(1))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 6)THEN
                IF(IPAR .EQ. ICHILD(2) .OR. IPAR .EQ. ICHILD(3))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 7)THEN 
                IF(IPAR .EQ. ICHILD(1))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 8)THEN
                IF(IPAR .EQ. ICHILD(1) .OR. IPAR .EQ. ICHILD(2))THEN
                    IFIXFLAG = 1
                END IF
              ELSEIF (NB .EQ. 9)THEN
                IF(IPAR .EQ. ICHILD(2))THEN
                    IFIXFLAG = 1
                END IF
              ENDIF
            ENDIF
          END DO
        END DO
      END DO

      RETURN
      END


C*************************************************************************
C     The following subroutine is designed to take a correctly defined
C     tree and alter it so that no two boxes that contact each other
C     are more than one level apart.  This is corrected only by adding
C     boxes.  The procedure involves flagging down bigger boxes and
C     dividing them and their children as necessary.
C     This routine also produces an array of colleagues for each box
C     that is numbered in the correct order.  All of the children are set
C     so that they satisfy our ordering convention.
C     The algorithm in the periodic case works the same way, it is just 
C     that upon subdivision the new colleagues must be put down to 
C     account for the periodicity.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICOLLEAGBOX denotes the colleagues of a given box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     IPERIOD denotes what kind of colleagues are to be generated
C             IPERIOD = 0 : free space
C             IPERIOD = 1 or 2 : periodic
C             IPERIOD = 3 : periodic up/down and free space left/right
C
C     IFLAG is just a dummy array
C
C     MAXBOXES is the maximum number of boxes we have storage for
C
C     OUTPUT:
C
C     ICOLBOX, IROWBOX, ICOLLEAGBOX, NBOXES, and all other
C     arrays describing the tree may be change on output
C
C*************************************************************************
      SUBROUTINE FIXTREE(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX, 
     1             IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2             NBLEVEL,IBOXLEV,ISTARTLEV,IPERIOD,
     3             IFLAG, MAXBOXES,ITEMPARRAY)
      IMPLICIT NONE
C-----Global variables
      INTEGER LEVELBOX(1), ICOLLEAGBOX(9,1)
      INTEGER MAXBOXES
      INTEGER IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER NBOXES, NLEV
      INTEGER NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER IPERIOD
      INTEGER ITEMPARRAY(1)
C-----Local variables
      INTEGER IFLAG(MAXBOXES)
      INTEGER ICHILD(4),ICOLL(9), IBOX
      INTEGER I, IPAR, ITEST, J, NB
      INTEGER ITEMP, NTEMP, JCNTR, ICNTR
      INTEGER START, ISTOP


C     Let's sort all of the boxes by level.
C     This takes the place of the LADDER structure
C     in the uniform case.  All boxes are sorted
C     into the array and placed in their proper places.
      CALL SORTBOXES(LEVELBOX,NBOXES,NLEV,
     1           NBLEVEL,IBOXLEV,ISTARTLEV)


C     First let's call a subroutine that will
C     generate all of the colleagues for each
C     box.  The colleagues are generated in the
C     correct order so there is no need to 'shuffle'
C     them later on.
      CALL MKCOLLS(ICOLBOX,
     1       IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2       IPARENTBOX,ICHILDBOX,NBLEVEL,
     3       IBOXLEV, ISTARTLEV,IPERIOD)


C     Let's initialize all of the flags to zero.
      DO I = 1, NBOXES
        IFLAG(I) = 0
      END DO

C     Find all of the boxes that need to be
C     flagged.  A flagged box will be denoted by 
C     setting IFLAG(box) = 1.
C     This refers to any box that is directly touching 
C     a box that is more than one level smaller than
C     it.  It is found by performing an upward pass
C     and looking a box's parents parents and seeing
C     if they are childless and contact the given box.
C     Note that we only need to get up to level two, as
C     we will not find a violation at a coarser level
C     than that.
      DO I = NLEV, 2, -1
        DO J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
          IBOX = IBOXLEV(J)
          IPAR  = IPARENTBOX(IBOX)
          ITEST = IPARENTBOX(IPAR)

            ICOLL(1) = ICOLLEAGBOX(1,ITEST)
            ICOLL(2) = ICOLLEAGBOX(2,ITEST)
            ICOLL(3) = ICOLLEAGBOX(3,ITEST)
            ICOLL(4) = ICOLLEAGBOX(4,ITEST)
            ICOLL(5) = ICOLLEAGBOX(5,ITEST)
            ICOLL(6) = ICOLLEAGBOX(6,ITEST)
            ICOLL(7) = ICOLLEAGBOX(7,ITEST)
            ICOLL(8) = ICOLLEAGBOX(8,ITEST)
            ICOLL(9) = ICOLLEAGBOX(9,ITEST)

            ICHILD(1) = ICHILDBOX(1,ITEST)
            ICHILD(2) = ICHILDBOX(2,ITEST)
            ICHILD(3) = ICHILDBOX(3,ITEST)
            ICHILD(4) = ICHILDBOX(4,ITEST)


          DO NB = 1, 9
            ITEMP = ICOLL(NB)
            IF(ICHILDBOX(1,ITEMP) .LT. 0)THEN
C             The neighboring box is not divided
C             we could have problems.
              IF (NB .EQ. 1)THEN
                IF(IPAR .EQ. ICHILD(4))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 2)THEN
                IF(IPAR .EQ. ICHILD(3) .OR. IPAR .EQ. ICHILD(4))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 3)THEN
                IF(IPAR .EQ. ICHILD(3))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 4)THEN
                IF(IPAR .EQ. ICHILD(4) .OR. IPAR .EQ. ICHILD(1))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 6)THEN
                IF(IPAR .EQ. ICHILD(2) .OR. IPAR .EQ. ICHILD(3))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 7)THEN 
                IF(IPAR .EQ. ICHILD(1))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 8)THEN
                IF(IPAR .EQ. ICHILD(1) .OR. IPAR .EQ. ICHILD(2))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ELSEIF (NB .EQ. 9)THEN
                IF(IPAR .EQ. ICHILD(2))THEN
                    IFLAG(ITEMP) = 1
                END IF
              ENDIF
            ENDIF
          END DO
        END DO
      END DO


C     Find all of the boxes that need to be
C     given a flag+.  A flag+ box will be denoted by 
C     setting IFLAG(box) = 2.
C     This refers to any box that is not already flagged
C     and is bigger than and is contacting a flagged box
C     or another box that has already been given a flag+.
C     It is found by performing an upward pass
C     and looking at a flagged box's parents colleagues
C     and a flag+ box's parents colleagues and seeing if
C     they are childless and present the case where a 
C     bigger box is contacting a flagged or a flag+ box.
      DO I = NLEV, 2, -1
        DO J = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
         IBOX = IBOXLEV(J)
          IF(IFLAG(IBOX) .EQ. 1 .OR. IFLAG(IBOX) .EQ. 2)THEN

          IPAR  = IPARENTBOX(IBOX)
 
            ICOLL(1) = ICOLLEAGBOX(1,IPAR)
            ICOLL(2) = ICOLLEAGBOX(2,IPAR)
            ICOLL(3) = ICOLLEAGBOX(3,IPAR)
            ICOLL(4) = ICOLLEAGBOX(4,IPAR)
            ICOLL(5) = ICOLLEAGBOX(5,IPAR)
            ICOLL(6) = ICOLLEAGBOX(6,IPAR)
            ICOLL(7) = ICOLLEAGBOX(7,IPAR)
            ICOLL(8) = ICOLLEAGBOX(8,IPAR)
            ICOLL(9) = ICOLLEAGBOX(9,IPAR)

            ICHILD(1) = ICHILDBOX(1,IPAR)
            ICHILD(2) = ICHILDBOX(2,IPAR)
            ICHILD(3) = ICHILDBOX(3,IPAR)
            ICHILD(4) = ICHILDBOX(4,IPAR)
          

          DO NB = 1, 9
            ITEMP = ICOLL(NB)
C           Let's check using the same criteria as above, but noting that
C           a flag will take precedence over a flag+.
            IF(ICHILDBOX(1,ITEMP) .LT. 0 
     1          .AND. IFLAG(ITEMP) .NE. 1)THEN
C             The neighboring box is not divided
C             we could have problems.
              IF (NB .EQ. 1)THEN
                IF(IBOX .EQ. ICHILD(4))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 2)THEN
                IF(IBOX .EQ. ICHILD(3) .OR. IBOX .EQ. ICHILD(4))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 3)THEN
                IF(IBOX .EQ. ICHILD(3))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 4)THEN
                IF(IBOX .EQ. ICHILD(4) .OR. IBOX .EQ. ICHILD(1))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 6)THEN
                IF(IBOX .EQ. ICHILD(2) .OR. IBOX .EQ. ICHILD(3))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 7)THEN
                IF(IBOX .EQ. ICHILD(1))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 8)THEN
                IF(IBOX .EQ. ICHILD(1) .OR. IBOX .EQ. ICHILD(2))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ELSEIF (NB .EQ. 9)THEN
                IF(IBOX .EQ. ICHILD(2))THEN
                    IFLAG(ITEMP) = 2
                END IF
              ENDIF
            ENDIF
           END DO
          ENDIF
        END DO
      END DO


C     Now let's divide the boxes that need to be immediately
C     divided up.  All of the flagged and flag+ boxes need to
C     be divided one time.  The distinction lies in the fact
C     that the children of a flag+ box will never need to be
C     divided but the children of a flagged box may need to 
C     be divided further.
C     Below, all flagged and flag+ boxes are divided once.  The
C     children of a flag+ box are left unflagged while those of
C     the flagged boxes are given a flag++ (denoted by setting
C     IFLAG(box) = 3) which will be needed in the downward pass.     
      NTEMP = NBOXES
      DO I = 1, NTEMP
C      Divide flagged boxes:
       IF (IFLAG(I) .EQ. 1)THEN

         IF(ICHILDBOX(1,I) .LT. 0)THEN
         CALL SUBDIVIDE(I,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C        Give flag++ to children of flagged boxes.
         ITEMP = ICHILDBOX(1,I)
         IFLAG(ITEMP) = 3

         ITEMP = ICHILDBOX(2,I)
         IFLAG(ITEMP) = 3

         ITEMP = ICHILDBOX(3,I)
         IFLAG(ITEMP) = 3

         ITEMP = ICHILDBOX(4,I)
         IFLAG(ITEMP) = 3

C      Divide flag+ boxes.
       ELSEIF (IFLAG(I) .EQ. 2)THEN
  
         IF(ICHILDBOX(1,I) .LT. 0)THEN
         CALL SUBDIVIDE(I,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF

       ENDIF
      END DO 


C     Now we need to do a downward pass.
C     We will concern ourselves only with the children of
C     flagged boxes and their children.  At each level,
C     for each flag++ box, test colleagues children and see
C     if they have children that are contacting you.  If so,
C     divide and flag++ all children that are created.     

      DO I = 0, NLEV
      NTEMP = NBOXES
      START = ISTARTLEV(I)
      ISTOP  = ISTARTLEV(I) + NBLEVEL(I) - 1
      DO 500 J = START, ISTOP
       IBOX = IBOXLEV(J)
C      Only be concerned with boxes on this level and
C      boxes that are given a flag++:
       IF(IFLAG(IBOX) .NE. 3)GOTO 500

         ICOLL(1) = ICOLLEAGBOX(1,IBOX)
         ICOLL(2) = ICOLLEAGBOX(2,IBOX)
         ICOLL(3) = ICOLLEAGBOX(3,IBOX)
         ICOLL(4) = ICOLLEAGBOX(4,IBOX)
         ICOLL(5) = ICOLLEAGBOX(5,IBOX)
         ICOLL(6) = ICOLLEAGBOX(6,IBOX)
         ICOLL(7) = ICOLLEAGBOX(7,IBOX)
         ICOLL(8) = ICOLLEAGBOX(8,IBOX)
         ICOLL(9) = ICOLLEAGBOX(9,IBOX)


C       Scan colleagues.
        DO 400 JCNTR = 1, 9
        IF(ICOLL(JCNTR) .LT. 0)GOTO 400
        IF(ICHILDBOX(1,ICOLL(JCNTR)) .LT. 0)GOTO 400

         ICHILD(1) = ICHILDBOX(1,ICOLL(JCNTR))
         ICHILD(2) = ICHILDBOX(2,ICOLL(JCNTR))
         ICHILD(3) = ICHILDBOX(3,ICOLL(JCNTR))
         ICHILD(4) = ICHILDBOX(4,ICOLL(JCNTR))


C          Scan colleague's children.
           DO 300 ICNTR = 1, 4
           IF (ICHILDBOX(1,ICHILD(ICNTR)) .LT. 0)GOTO 300 

           IF(JCNTR .EQ. 1 .AND. ICNTR .EQ. 2)THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF



C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 2 .AND. 
     1        (ICNTR .EQ. 1 .OR. ICNTR .EQ. 2))THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF



C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 3 .AND. ICNTR .EQ. 1)THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF



C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 4 .AND. 
     1        (ICNTR .EQ. 2 .OR. ICNTR .EQ. 3))THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 6 .AND. 
     1        (ICNTR .EQ. 1 .OR. ICNTR .EQ. 4))THEN
             
             
C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 7 .AND. ICNTR .EQ. 3)THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 8 .AND. 
     1        (ICNTR .EQ. 3 .OR. ICNTR .EQ. 4))THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3

           ELSEIF(JCNTR .EQ. 9 .AND. ICNTR .EQ. 4)THEN

C           Call subdivide
         IF(ICHILDBOX(1,IBOX) .LT. 0)THEN
            CALL SUBDIVIDE(IBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,
     3         ITEMPARRAY)
         ENDIF


C           flag++ all children created
            ITEMP = ICHILDBOX(1,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(2,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(3,IBOX)
            IFLAG(ITEMP) = 3
            ITEMP = ICHILDBOX(4,IBOX)
            IFLAG(ITEMP) = 3
           ENDIF
300     CONTINUE   
400     CONTINUE
500     CONTINUE
      END DO
      RETURN
      END


C**********************************************************************
C     The following subroutine is used to set up the arrays of parents
C     and children for each box.  If a box is childless, all of its
C     children are set to -1.
C     Each box that has children has four.  They are numbered
C     as follows:
C
C                        1     2
C
C                        4     3
C
C
C     This is done merely by scanning through all of the boxes and
C     comparing the column and row numbers of boxes that are one level
C     apart.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     OUTPUT:
C
C     ICHILDBOX denotes the four children of each box
C
C**********************************************************************
      SUBROUTINE MKCHILD(LEVELBOX,IPARENTBOX,ICHILDBOX,ICOLBOX,
     1         IROWBOX,NBOXES,NLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER LEVELBOX(1)
      INTEGER IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER NBOXES, NLEV
C-----Local variables
      INTEGER I, J, L
      INTEGER ICOL, IROW
      INTEGER ICOLTEST, IROWTEST

C     First let's initialize all of the ICHILDBOX
C     arrays to zero:
      DO J = 1, NBOXES
        IPARENTBOX(J) = -1
        DO I = 1, 4
          ICHILDBOX(I,J) = -1
        END DO
      END DO

      DO I = 0, NLEV
        DO 200 J = 1, NBOXES
        IF (LEVELBOX(J) .NE. I)GOTO 200
C         Figure out where the first child should be
C         by looking at the parents row and column.
          ICOLTEST = 2*(ICOLBOX(J) - 1) + 1
          IROWTEST = 2*(IROWBOX(J) - 1) + 1
C         Now scan all of the boxes on the next
C         finest level.
          DO 100 L = 1, NBOXES
          IF(LEVELBOX(L) .NE. I+1)GOTO 100
            ICOL = ICOLBOX(L)
            IROW = IROWBOX(L)
            IF(ICOL .EQ. ICOLTEST .AND. IROW .EQ. IROWTEST)THEN
              ICHILDBOX(4,J) = L
              IPARENTBOX(L) = J
            ELSEIF(ICOL .EQ. ICOLTEST + 1 
     1                        .AND. IROW .EQ. IROWTEST + 1)THEN
              ICHILDBOX(2,J) = L
              IPARENTBOX(L) = J
            ELSEIF(ICOL .EQ. ICOLTEST 
     1                        .AND. IROW .EQ. IROWTEST + 1)THEN
              ICHILDBOX(1,J) = L
              IPARENTBOX(L) = J
            ELSEIF(ICOL .EQ. ICOLTEST + 1 
     1                        .AND. IROW .EQ. IROWTEST)THEN
              ICHILDBOX(3,J) = L
              IPARENTBOX(L) = J
            ENDIF
100       CONTINUE
200     CONTINUE
      END DO
      RETURN
      END


C**********************************************************************
C     The following subroutine is designed to divide up a childless
C     box into four children.  
C     The children are placed in correct order (clockwise starting
C     from the upper left corner) so there is no need to 'shuffle' 
C     the child order later on. In the periodic case, the colleagues
C     must be obtained by looking at the potential colleague numbers
C     and their row and column and seeing if they lie outside of 
C     the domain. If they do it must be readjusted to account for the 
C     periodicity.
C
C     INPUT:
C
C     IPARBOX denotes the box being divided
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     ICOLLEAGBOX denotes the colleagues of a given box
C
C     NBOXES is the total number of boxes
C
C     IROWBOX denotes the row of each box
C
C     ICOLBOX denotes the column of each box
C
C     LEVELBOX is an array determining the level of each box
C
C     NLEV is the finest level
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     IPERIOD denotes what kind of colleagues are to be generated
C             IPERIOD = 1 : free space
C             IPERIOD = 2 : periodic
C             IPERIOD = 3 : periodic up/down and free space left/right
C
C     ITEMPARRAY is just a dummy array
C
C     OUTPUT:
C     
C     NBOXES and ICHILDBOX are altered to
C            reflect the addition of new boxes
C
C**********************************************************************
      SUBROUTINE SUBDIVIDE(IPARBOX,IPARENTBOX,ICHILDBOX,ICOLLEAGBOX, 
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV, IPERIOD,ITEMPARRAY)
      IMPLICIT NONE
C-----Global variables
      INTEGER  IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  LEVELBOX(1), NBOXES
      INTEGER  IPARBOX, NLEV, IPERIOD
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  ICOLLEAGBOX(9,1)
C-----Local variables
      INTEGER  LEVEL, IBOX, ITEMPARRAY(1)
      INTEGER  ICOLUMN, IROW, I, J
      INTEGER  ICNTR, JCNTR, ISISTER
      INTEGER  ICOLTEMP, ICOLTEST
      INTEGER  IROWTEMP, IROWTEST
      INTEGER  PARTEMP, COLLEAGUE
      INTEGER  NSIDE, ILEV, ITEST, L

C     Let's initialize the array ITEMPARRAY to zero:
      DO I = 1, NBOXES + 4
        ITEMPARRAY(I) = 0
      END DO

C     Level, icolumn, and irow refer to the level, column,
C     and row of the parent box, respectively.
      LEVEL   = LEVELBOX(IPARBOX)
      ICOLUMN = ICOLBOX(IPARBOX)
      IROW    = IROWBOX(IPARBOX)

C     Here are the new boxes placed in the
C     correct positions.  They are all childless.
C     There columns and rows are determined from
C     the parents columns and rows.  The level is
C     obviously one level finer than the parent.
      IBOX = NBOXES + 1
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 1
      IROWBOX(IBOX) = 2*(IROW-1) + 1

      IBOX = NBOXES + 2
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 2
      IROWBOX(IBOX) = 2*(IROW-1) + 1

      IBOX = NBOXES + 3
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 1
      IROWBOX(IBOX) = 2*(IROW-1) + 2

      IBOX = NBOXES + 4
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 2
      IROWBOX(IBOX) = 2*(IROW-1) + 2

      ICHILDBOX(1,IPARBOX) = NBOXES + 3
      ICHILDBOX(2,IPARBOX) = NBOXES + 4
      ICHILDBOX(3,IPARBOX) = NBOXES + 2
      ICHILDBOX(4,IPARBOX) = NBOXES + 1

C     Set up a temporary array to store the old one in:
      DO I = 1, NBOXES
        ITEMPARRAY(I) = IBOXLEV(I)
      END DO

C     Now let's rearrange the ladder structure:
      NBLEVEL(LEVEL + 1) = NBLEVEL(LEVEL + 1) + 4

      IBOXLEV(ISTARTLEV(LEVEL+2))   = NBOXES + 1
      IBOXLEV(ISTARTLEV(LEVEL+2)+1) = NBOXES + 2
      IBOXLEV(ISTARTLEV(LEVEL+2)+2) = NBOXES + 3
      IBOXLEV(ISTARTLEV(LEVEL+2)+3) = NBOXES + 4

      DO I = ISTARTLEV(LEVEL+2) + 4, NBOXES + 4
        IBOXLEV(I) = ITEMPARRAY(I - 4)
      END DO 

      DO I = LEVEL + 2, NLEV
        ISTARTLEV(I) = ISTARTLEV(I) + 4
      END DO
      NBOXES = NBOXES + 4


C      Now let's go through the process of reforming any 
C      necessary colleagues.  For each of the child boxes 
C      that we just formed, all we need to do is scan through
C      the boxes that are children of the above parent boxes 
C      colleagues and test the column and row numbers.  We can 
C      also take advantage of the fact that for every one of 
C      the newly formed boxes colleagues, that box will list 
C      the newly formed box as one of its colleagues.  
C      The colleague numbers can be found easily if we think 
C      of a 'reflection.'  Colleague 1 and 9 are opposites, 
C      3 and 7 etc.
C      First do the free space case:
       IF(IPERIOD .EQ. 0)THEN
       DO 200 I = 1, 4
         IF(ICHILDBOX(1,IPARBOX) .LT. 0)GOTO 200
         IBOX = ICHILDBOX(I,IPARBOX)
         ICOLLEAGBOX(5,IBOX) = IBOX
         DO J = 1, 4
           ICOLLEAGBOX(J,IBOX) = -1
         END DO
         DO J = 6, 9
           ICOLLEAGBOX(J,IBOX) = -1
         END DO

         PARTEMP = IPARENTBOX(IBOX)

C        IROWTEMP and ICOLTEMP denote the
C        row and column of the test box.
         IROWTEMP = IROWBOX(IBOX)
         ICOLTEMP = ICOLBOX(IBOX)

         DO 100 JCNTR = 1, 9
C          COLLEAGUE denotes the colleague of the parent box.
           COLLEAGUE = ICOLLEAGBOX(JCNTR,PARTEMP)
C          If the colleague doesn't exist
C          or is childless, skip it:
           IF (COLLEAGUE .LT. 0)GOTO 100
           IF (ICHILDBOX(1,COLLEAGUE) .LT. 0)GOTO 100
C          Otherwise scan the four children:
           DO ICNTR = 1, 4
             J = ICHILDBOX(ICNTR,COLLEAGUE)
C            IROWTEST and ICOLTEST denote the row and column of
C            the box being compared to the test box.
             IROWTEST = IROWBOX(J)
             ICOLTEST = ICOLBOX(J)

             IF(IROWTEMP .EQ. IROWTEST+1)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(1,IBOX) = J
                 ICOLLEAGBOX(9,J) = IBOX
               ELSEIF(ICOLTEMP .EQ. ICOLTEST)THEN
                 ICOLLEAGBOX(2,IBOX) = J
                 ICOLLEAGBOX(8,J) = IBOX
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(3,IBOX) = J
                 ICOLLEAGBOX(7,J) = IBOX
               ENDIF
             ELSEIF(IROWTEMP .EQ. IROWTEST)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(4,IBOX) = J
                 ICOLLEAGBOX(6,J) = IBOX
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(6,IBOX) = J
                 ICOLLEAGBOX(4,J) = IBOX
               ENDIF
             ELSEIF(IROWTEMP .EQ. IROWTEST-1)THEN
               IF(ICOLTEMP .EQ. ICOLTEST+1)THEN
                 ICOLLEAGBOX(7,IBOX) = J
                 ICOLLEAGBOX(3,J) = IBOX
               ELSEIF(ICOLTEMP .EQ. ICOLTEST)THEN
                 ICOLLEAGBOX(8,IBOX) = J
                 ICOLLEAGBOX(2,J) = IBOX
               ELSEIF(ICOLTEMP .EQ. ICOLTEST-1)THEN
                 ICOLLEAGBOX(9,IBOX) = J
                 ICOLLEAGBOX(1,J) = IBOX
               ENDIF
             ENDIF
          END DO
100      CONTINUE
200   CONTINUE
C     Next do the periodic case:
      ELSEIF(IPERIOD .EQ. 1 .OR. IPERIOD .EQ. 2)THEN
       ILEV = LEVEL + 1
       NSIDE = 2**ILEV

       DO L = 1, 4
        IBOX = ICHILDBOX(L,IPARBOX)
C       Initialize colleague number 5 to
C       yourself and all other colleagues to
C       -1.  -1 is the flag for the case when
C       the colleagues don't exist.  This is 
C       for all of the newly created boxes.
         ICOLLEAGBOX(5,IBOX) = IBOX
         DO J = 1, 4
           ICOLLEAGBOX(J,IBOX) = -1
         END DO
         DO J = 6, 9
           ICOLLEAGBOX(J,IBOX) = -1
         END DO

C       IROWTEMP and ICOLTEMP denote the row and column of
C       the test box.
        IROWTEMP = IROWBOX(IBOX)
        ICOLTEMP = ICOLBOX(IBOX)

C       IROWTEST and ICOLTEST denote the row and column of
C       the box being compared to the test box.
             DO 400 JCNTR = 1, 9
               IF(JCNTR .EQ. 5)GOTO 400
               IF(JCNTR .EQ. 1)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 9
               ELSEIF(JCNTR .EQ. 2)THEN 
                 ICOLTEST = ICOLTEMP
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 8
               ELSEIF(JCNTR .EQ. 3)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 7
               ELSEIF(JCNTR .EQ. 4)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP
                 ISISTER = 6
               ELSEIF(JCNTR .EQ. 6)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP
                 ISISTER = 4
               ELSEIF(JCNTR .EQ. 7)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 3
               ELSEIF(JCNTR .EQ. 8)THEN
                 ICOLTEST = ICOLTEMP
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 2
               ELSEIF(JCNTR .EQ. 9)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 1
               ENDIF

C         Now test to see if the test parameters
C         lie in the domain.
            IF(ICOLTEST .LT. 1)THEN
               ICOLTEST = ICOLTEST + NSIDE
            ELSEIF(ICOLTEST .GT. NSIDE)THEN
               ICOLTEST = ICOLTEST - NSIDE
            ENDIF
            IF(IROWTEST .LT. 1)THEN
               IROWTEST = IROWTEST + NSIDE
            ELSEIF(IROWTEST .GT. NSIDE)THEN
               IROWTEST = IROWTEST - NSIDE
            ENDIF


       DO 300 J = 1, 9
        IF(ICOLLEAGBOX(J,IPARENTBOX(IBOX)) .LT. 0)GOTO 300
        IF(ICHILDBOX(1,ICOLLEAGBOX(J,IPARENTBOX(IBOX))) .LT. 0)GOTO 300
          DO ICNTR = 1, 4
            ITEST = ICHILDBOX(ICNTR,ICOLLEAGBOX(J,IPARENTBOX(IBOX)))
            IF(IROWBOX(ITEST) .EQ. IROWTEST
     1         .AND. ICOLBOX(ITEST) .EQ. ICOLTEST)THEN
               ICOLLEAGBOX(JCNTR,IBOX) = ITEST
               ICOLLEAGBOX(ISISTER,ITEST) = IBOX
            ENDIF
          END DO
300    CONTINUE
400    CONTINUE
      END DO
      ELSEIF(IPERIOD .EQ. 3)THEN
       ILEV = LEVEL + 1
       NSIDE = 2**ILEV

       DO L = 1, 4
        IBOX = ICHILDBOX(L,IPARBOX)
C       Initialize colleague number 5 to
C       yourself and all other colleagues to
C       -1.  -1 is the flag for the case when
C       the colleagues don't exist.  This is 
C       for all of the newly created boxes.
         ICOLLEAGBOX(5,IBOX) = IBOX
         DO J = 1, 4
           ICOLLEAGBOX(J,IBOX) = -1
         END DO
         DO J = 6, 9
           ICOLLEAGBOX(J,IBOX) = -1
         END DO

C       IROWTEMP and ICOLTEMP denote the row and column of
C       the test box.
        IROWTEMP = IROWBOX(IBOX)
        ICOLTEMP = ICOLBOX(IBOX)

C       IROWTEST and ICOLTEST denote the row and column of
C       the box being compared to the test box.
             DO 600 JCNTR = 1, 9
               IF(JCNTR .EQ. 5)GOTO 600
               IF(JCNTR .EQ. 1)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 9
               ELSEIF(JCNTR .EQ. 2)THEN 
                 ICOLTEST = ICOLTEMP
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 8
               ELSEIF(JCNTR .EQ. 3)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP - 1
                 ISISTER = 7
               ELSEIF(JCNTR .EQ. 4)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP
                 ISISTER = 6
               ELSEIF(JCNTR .EQ. 6)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP
                 ISISTER = 4
               ELSEIF(JCNTR .EQ. 7)THEN
                 ICOLTEST = ICOLTEMP - 1
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 3
               ELSEIF(JCNTR .EQ. 8)THEN
                 ICOLTEST = ICOLTEMP
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 2
               ELSEIF(JCNTR .EQ. 9)THEN
                 ICOLTEST = ICOLTEMP + 1
                 IROWTEST = IROWTEMP + 1
                 ISISTER = 1
               ENDIF

C         Now test to see if the test parameters
C         lie in the domain.
            IF(IROWTEST .LT. 1)THEN
               IROWTEST = IROWTEST + NSIDE
            ELSEIF(IROWTEST .GT. NSIDE)THEN
               IROWTEST = IROWTEST - NSIDE
            ENDIF


       DO 500 J = 1, 9
        IF(ICOLLEAGBOX(J,IPARENTBOX(IBOX)) .LT. 0)GOTO 500
        IF(ICHILDBOX(1,ICOLLEAGBOX(J,IPARENTBOX(IBOX))) .LT. 0)GOTO 500
          DO ICNTR = 1, 4
            ITEST = ICHILDBOX(ICNTR,ICOLLEAGBOX(J,IPARENTBOX(IBOX)))
            IF(IROWBOX(ITEST) .EQ. IROWTEST
     1         .AND. ICOLBOX(ITEST) .EQ. ICOLTEST)THEN
               ICOLLEAGBOX(JCNTR,IBOX) = ITEST
               ICOLLEAGBOX(ISISTER,ITEST) = IBOX
            ENDIF
          END DO
500    CONTINUE
600    CONTINUE
      END DO
      ENDIF
      RETURN
      END


C********************************************************************
C     This routine is identical to the SUBDIVIDE routine except that
C     it does not concern itself at all with generating colleagues.
C     This routine is used only within the MKTREE8 routine.
C
C     INPUT:
C
C     IPARBOX denotes the box being divided
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NBOXES is the total number of boxes
C
C     IROWBOX denotes the row of each box
C
C     ICOLBOX denotes the column of each box
C
C     LEVELBOX is an array determining the level of each box
C
C     NLEV is the finest level
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ITEMPARRAY is just a dummy array
C
C     OUTPUT:
C
C     NBOXES and ICHILDBOX are altered to reflect the addition
C            of new boxes
C
C********************************************************************
      SUBROUTINE SUBDIVIDE1(IPARBOX,IPARENTBOX,ICHILDBOX,
     1         NBOXES,IROWBOX,ICOLBOX,LEVELBOX,NLEV,
     2         ISTARTLEV, NBLEVEL, IBOXLEV,ITEMPARRAY)
      IMPLICIT NONE
C-----Global variables
      INTEGER IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER ICOLBOX(1), IROWBOX(1)
      INTEGER LEVELBOX(1), NBOXES
      INTEGER IPARBOX,NLEV
      INTEGER NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
C-----Local variables
      INTEGER LEVEL, IBOX, ITEMPARRAY(1)
      INTEGER ICOLUMN, IROW, I, J
      INTEGER NCNTR

C     Let's initialize the array ITEMPARRAY to zero:
      DO I = 1, NBOXES + 4
        ITEMPARRAY(I) = 0
      END DO

C     Level, icolumn, and irow refer to the level, column,
C     and row of the parent box, respectively.
      LEVEL   = LEVELBOX(IPARBOX)
      ICOLUMN = ICOLBOX(IPARBOX)
      IROW    = IROWBOX(IPARBOX)

C     Here are the new boxes placed in the
C     correct positions.  They are all childless.
C     There columns and rows are determined from
C     the parents columns and rows.  The level is
C     obviously one level finer than the parent.
      IBOX = NBOXES + 1
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 1
      IROWBOX(IBOX) = 2*(IROW-1) + 1

      IBOX = NBOXES + 2
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 2
      IROWBOX(IBOX) = 2*(IROW-1) + 1

      IBOX = NBOXES + 3
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 1
      IROWBOX(IBOX) = 2*(IROW-1) + 2

      IBOX = NBOXES + 4
      LEVELBOX(IBOX) = LEVEL + 1
      IPARENTBOX(IBOX) = IPARBOX
      ICHILDBOX(1,IBOX) = -1
      ICHILDBOX(2,IBOX) = -1
      ICHILDBOX(3,IBOX) = -1
      ICHILDBOX(4,IBOX) = -1
      ICOLBOX(IBOX) = 2*(ICOLUMN-1) + 2
      IROWBOX(IBOX) = 2*(IROW-1) + 2

      ICHILDBOX(1,IPARBOX) = NBOXES + 3
      ICHILDBOX(2,IPARBOX) = NBOXES + 4
      ICHILDBOX(3,IPARBOX) = NBOXES + 2
      ICHILDBOX(4,IPARBOX) = NBOXES + 1

C     Set up a temporary array to store the old one in:
      DO I = 1, NBOXES
        ITEMPARRAY(I) = IBOXLEV(I)
      END DO

C     Now let's rearrange the ladder structure:
      NBLEVEL(LEVEL + 1) = NBLEVEL(LEVEL + 1) + 4

      IBOXLEV(ISTARTLEV(LEVEL+2))   = NBOXES + 1
      IBOXLEV(ISTARTLEV(LEVEL+2)+1) = NBOXES + 2
      IBOXLEV(ISTARTLEV(LEVEL+2)+2) = NBOXES + 3
      IBOXLEV(ISTARTLEV(LEVEL+2)+3) = NBOXES + 4

      DO I = ISTARTLEV(LEVEL+2) + 4, NBOXES + 4
        IBOXLEV(I) = ITEMPARRAY(I - 4)
      END DO 

      DO I = LEVEL + 2, NLEV
        ISTARTLEV(I) = ISTARTLEV(I) + 4
      END DO
      NBOXES = NBOXES + 4


      NCNTR = 1
      DO I = 0, NLEV
        NBLEVEL(I) = 0
        ISTARTLEV(I) = NCNTR
        DO J = 1, NBOXES
          IF(LEVELBOX(J) .EQ. I)THEN
           IBOXLEV(NCNTR) = J
           NCNTR = NCNTR + 1
           NBLEVEL(I) = NBLEVEL(I) + 1
          ENDIF
        END DO
      END DO
      RETURN
      END


C***********************************************************************
C     The following subroutine sets up a structure that is analogous
C     to the 'LADDER' structure in the nonadaptive case.
C     It is just a way of organizing the boxes by level in one long
C     array and denoting where in the array the levels change.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICHILDBOX denotes the four children of each box
C
C     OUTPUT:
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C***********************************************************************
      SUBROUTINE SORTBOXES(LEVELBOX,NBOXES,NLEV,
     1           NBLEVEL,IBOXLEV,ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER LEVELBOX(1), NBOXES, NLEV
      INTEGER NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
C-----Local variables
      INTEGER NCNTR, I, J

        NCNTR = 1
        DO I = 0, NLEV
          NBLEVEL(I) = 0
          ISTARTLEV(I) = NCNTR
          DO J = 1, NBOXES
            IF(LEVELBOX(J) .EQ. I)THEN
             IBOXLEV(NCNTR) = J
             NCNTR = NCNTR + 1
             NBLEVEL(I) = NBLEVEL(I) + 1
            ENDIF
          END DO
        END DO
      RETURN
      END


C***********************************************************************
C     The following subroutine carves up the workspace array before it
C     is sent to the ADAPFMM6 subroutine, where the bulk of the work in
C     the algorithm is done.
C
C     IPERIOD is set to distinguish between each of the following cases:
C     IPERIOD = 0 for the free space case
C     IPERIOD = 1 for the periodic case
C     IPERIOD = 2 for the homogeneous dirichlet case
C     IPERIOD = 3 for the homogeneous neumann case
C     IPERIOD = 4 for the dirichlet left/right and periodic top/bottom
C     IPERIOD = 5 for the dirichlet left/right and neumann top/bottom
C     IPERIOD = 6 for the neumann left/right and periodic top/bottom
C     IPERIOD = 7 for the purely dirichlet (inhomogeneous)
C     IPERIOD = 8 for the purely neumann (inhomogeneous)
C     IPERIOD = 9 for the dirichlet left/right and periodic top/bottom
C     IPERIOD = 10 for the dirichlet left/right and neumann top/bottom
C     IPERIOD = 11 for the neumann left/right and periodic top/bottom
C
C     ICOLBOX, IROWBOX, NBOXES, LEVELBOX, IPARENTBOX, AND ICHILDBOX
C     define the necessary parts of the tree.
C***********************************************************************
      SUBROUTINE FMMSTART8(WORK, LENW, IWORK, ILEN,
     1     NLEV, LEVELBOX, IPARENTBOX, ICHILDBOX,
     2     ICOLBOX, IROWBOX, NBOXES, NBLEVEL, IBOXLEV,
     3     ISTARTLEV, IPERIOD, FRIGHT, POT, IPREC, MAP,
     4     HL, HR, HB, HT)
      IMPLICIT NONE
C-----Global variables
      INTEGER  LENW, ILEN
      INTEGER  NLEV, IPERIOD
      INTEGER  NNODES, NDEG
      INTEGER  NTERMS, NBOXES
      INTEGER  LEVELBOX(1), IPARENTBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  IROWBOX(1), ICOLBOX(1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      INTEGER  IWORK(ILEN), IPREC
      REAL *8  WORK(LENW)
      REAL *8  POT(64,1)
      REAL *8  MAP(1)
      REAL *8  FRIGHT(64,1)
      REAL *8  HL, HR, HB, HT
C-----Local variables
      INTEGER  NBOXES1, NLEV1
      INTEGER  LLEVELBOX, MLEVELBOX
      INTEGER  LPARENTBOX, MPARENTBOX
      INTEGER  LCOLBOX, MCOLBOX
      INTEGER  LROWBOX, MROWBOX
      INTEGER  LCHILDBOX, MCHILDBOX
      INTEGER  LSTARTLEV, MSTARTLEV
      INTEGER  LNBLEVEL, MNBLEVEL
      INTEGER  LBOXLEV, MBOXLEV
      INTEGER  LMPOLE, LCOEFF, MCOEFF
      INTEGER  IEXPN, IEXPS
      INTEGER  IEXPE, IEXPW, LEXP
      INTEGER  IEXPNBIG, IEXPSBIG
      INTEGER  IEXPEBIG, IEXPWBIG
      INTEGER  ITOT, LLOCXP
      INTEGER  LNODES, MXNODES, MWNODES, MTEMP
      INTEGER  IPERIODTEMP
      INTEGER  MCOLLEAG, LCOLLEAG
      INTEGER  MCOMP, LCOMP
      INTEGER  MCHOOSE, LCHOOSE
      INTEGER  MCOEFFSTOP, MCOEFFSSIDE
      INTEGER  MCOEFFDTOP, MCOEFFDSIDE
      INTEGER  MZS, LZS
      INTEGER  MZSEAST, MZSWEST
      INTEGER  MZSSOUTH, MZSNORTH
      INTEGER  MPOLE, LOCEXP
      INTEGER  MFTOPD, MFTOPS
      INTEGER  MFSIDED, MFSIDES
      INTEGER  LF, LCOEFFSIDE
      INTEGER  LVOLMAPS, LSIDEMAPS, LWINT
      INTEGER  MMAPSOUTH, MMAPNORTH, MMAPEAST, MMAPWEST
      INTEGER  MEDBLETOP, MEDBLESIDE, MESNGLETOP, MESNGLESIDE
      INTEGER  MWDBLETOP, MWDBLESIDE, MWSNGLETOP, MWSNGLESIDE
      INTEGER  MNDBLETOP, MNDBLESIDE, MNSNGLETOP, MNSNGLESIDE
      INTEGER  MSDBLETOP, MSDBLESIDE, MSSNGLETOP, MSSNGLESIDE
      INTEGER  MWINT
      INTEGER  MFLAGEAST, MFLAGWEST, LFLAG
      INTEGER  MFLAGNORTH, MFLAGSOUTH
      INTEGER  MLOCALONOFF, LLOCALONOFF
      INTEGER  MDOUBLETOPONOFF, LMULTONOFF
      INTEGER  MDOUBLESIDEONOFF
      INTEGER  MSINGLETOPONOFF
      INTEGER  MSINGLESIDEONOFF
   
      NDEG = 7

      IF(IPREC .EQ. 0)THEN
        NTERMS = 8
        NNODES = 8
      ELSEIF(IPREC .EQ. 1)THEN
        NTERMS = 16
        NNODES = 16
      ELSEIF(IPREC .EQ. 2)THEN
        NTERMS = 24
        NNODES = 24
      ELSEIF(IPREC .EQ. 3)THEN
        NTERMS = 40
        NNODES = 40
      ENDIF

C     First divide up the workspace appropriately
      IF(IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN

C        First break up the integer array:
C        (In this first case, the only workspace
C        needed is for the colleagues)
         LCOLLEAG   =  9 * NBOXES
         LFLAG    = NBOXES
         LLOCALONOFF   = NBOXES

         MCOLLEAG    =  1
         MLOCALONOFF =  MCOLLEAG + LCOLLEAG
         MFLAGEAST   =  MLOCALONOFF + LLOCALONOFF
         MFLAGWEST   =  MFLAGEAST + LFLAG
         MFLAGNORTH  =  MFLAGWEST + LFLAG
         MFLAGSOUTH  =  MFLAGNORTH + LFLAG


         ITOT = MFLAGSOUTH + LFLAG

         IF ( ITOT .GE. ILEN ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the integer array'
            WRITE(*,*)'to be ',ilen
            WRITE(*,*)'The workspace needed is: ',itot
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF

         LMPOLE   = 2*(NTERMS+1)*NBOXES
         LLOCXP   = 2*(NTERMS+1)*NBOXES
         LEXP     = 2*(NNODES+1)*NBOXES
         LCOEFF   = 64*NBOXES
         LNODES   = NNODES
         LCOMP    = NNODES*(NTERMS+1)
         LCHOOSE  = 2*NTERMS * 2*NTERMS
         LZS      = 2 * 7 * 7 * (NNODES + 1)
         LVOLMAPS = 2*(NNODES + 1) * 36
         LWINT    = 2*(NTERMS + 1) * 8 * 8

C
C        Now carve up the real workspace:
C
         MPOLE    = 1
         LOCEXP   = MPOLE   + LMPOLE
         IEXPN    = LOCEXP  + LLOCXP
         IEXPS    = IEXPN   + LEXP
         IEXPE    = IEXPS   + LEXP
         IEXPW    = IEXPE   + LEXP
         IEXPNBIG = IEXPW   + LEXP
         IEXPSBIG = IEXPNBIG   + LEXP
         IEXPEBIG = IEXPSBIG   + LEXP
         IEXPWBIG = IEXPEBIG   + LEXP
         MCOMP    = IEXPWBIG   + LEXP
         MCHOOSE  = MCOMP   + LCOMP
         MZS      = MCHOOSE + LCHOOSE
         MZSEAST  = MZS + LZS
         MZSWEST  = MZSEAST  + LZS
         MZSNORTH = MZSWEST  + LZS
         MZSSOUTH = MZSNORTH + LZS
         MTEMP    = MZSSOUTH + LZS
         MWNODES  = MTEMP   + LNODES
         MXNODES  = MWNODES + LNODES
         MCOEFF   = MXNODES + LNODES

         ITOT = MCOEFF + LCOEFF

         IF ( ITOT .GE. LENW ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the real array'
            WRITE(*,*)'to be ',lenw
            WRITE(*,*)'The workspace needed is: ',itot
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF

         MMAPSOUTH   =  1
         MMAPNORTH   =  MMAPSOUTH   + LVOLMAPS
         MMAPEAST    =  MMAPNORTH   + LVOLMAPS
         MMAPWEST    =  MMAPEAST    + LVOLMAPS
         MWINT       =  MMAPWEST    + LVOLMAPS

      ELSEIF(IPERIOD .EQ. 2 .OR. IPERIOD .EQ. 3 .OR.
     1       IPERIOD .EQ. 4 .OR. IPERIOD .EQ. 5 .OR.
     2       IPERIOD .EQ. 6)THEN

C        First break up the integer array:
C        (In this case, we need to copy the
C        initial data structure over again so
C        that there is not a problem with
C        altering the structure on input)
         LCOLLEAG    =  9 * (4*NBOXES  + 1)
         LLEVELBOX   =  4 * NBOXES + 1
         LPARENTBOX  =  4 * NBOXES + 1
         LCOLBOX     =  4 * NBOXES + 1
         LROWBOX     =  4 * NBOXES + 1
         LCHILDBOX   =  4 * (4 * NBOXES + 1)
         LSTARTLEV   =  NLEV + 2
         LNBLEVEL    =  NLEV + 2
         LBOXLEV     =  4*(NBOXES + 1)
         LFLAG       =  4*(NBOXES + 1)
         LLOCALONOFF =  4*(NBOXES + 1)

         MCOLLEAG   =  1
         MLEVELBOX  =  MCOLLEAG + LCOLLEAG
         MLOCALONOFF=  MLEVELBOX + LLEVELBOX
         MPARENTBOX =  MLOCALONOFF + LLOCALONOFF
         MROWBOX    =  MPARENTBOX + LPARENTBOX
         MCOLBOX    =  MROWBOX + LROWBOX
         MCHILDBOX  =  MCOLBOX + LCOLBOX
         MSTARTLEV  =  MCHILDBOX + LCHILDBOX
         MNBLEVEL   =  MSTARTLEV + LSTARTLEV
         MBOXLEV    =  MNBLEVEL + LNBLEVEL
         MFLAGEAST  =  MBOXLEV + LBOXLEV
         MFLAGWEST  =  MFLAGEAST + LFLAG
         MFLAGNORTH =  MFLAGWEST + LFLAG
         MFLAGSOUTH =  MFLAGNORTH + LFLAG

         ITOT = MFLAGSOUTH + LFLAG

         IF ( ITOT .GE. ILEN ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the integer array'
            WRITE(*,*)'to be ',ilen
            WRITE(*,*)'The workspace needed is: ',itot 
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF

C
C        Now carve up the real workspace:
C
         LMPOLE     =  2*(NTERMS+1)*(4*NBOXES + 1)
         LLOCXP     =  2*(NTERMS+1)*(4*NBOXES + 1)
         LEXP       =  2*(NNODES+1)*(4*NBOXES + 1)
         LCOEFF     =  64*(4*NBOXES + 1)
         LNODES     =  NNODES
         LCOMP      =  NNODES*(NTERMS+1)
         LCHOOSE    =  2*NTERMS * 2*NTERMS
         LZS        =  2 * 7 * 7 * (NNODES + 1)
         LCOEFFSIDE =  8 * (4*NBOXES + 1)
         LF         =  8 * (4*NBOXES + 1)
         LVOLMAPS   =  2*(NNODES + 1) * 36
         LWINT      =  2*(NTERMS + 1) * 8 * 8

         MMAPSOUTH   =  1
         MMAPNORTH   =  MMAPSOUTH   + LVOLMAPS
         MMAPEAST    =  MMAPNORTH   + LVOLMAPS
         MMAPWEST    =  MMAPEAST    + LVOLMAPS
         MWINT       =  MMAPWEST    + LVOLMAPS
         MPOLE       =  1
         LOCEXP      =  MPOLE + LMPOLE
         IEXPN       =  LOCEXP + LLOCXP
         IEXPS       =  IEXPN + LEXP
         IEXPE       =  IEXPS + LEXP
         IEXPW       =  IEXPE + LEXP
         IEXPNBIG    =  IEXPW + LEXP
         IEXPSBIG    =  IEXPNBIG   + LEXP
         IEXPEBIG    =  IEXPSBIG   + LEXP
         IEXPWBIG    =  IEXPEBIG   + LEXP
         MCOMP       =  IEXPWBIG   + LEXP
         MCHOOSE     =  MCOMP + LCOMP
         MZS         =  MCHOOSE + LCHOOSE
         MZSEAST     =  MZS + LZS
         MZSWEST     =  MZSEAST  + LZS
         MZSNORTH    =  MZSWEST  + LZS
         MZSSOUTH    =  MZSNORTH + LZS
         MTEMP       =  MZSSOUTH + LZS
         MWNODES     =  MTEMP + LNODES
         MXNODES     =  MWNODES + LNODES
         MCOEFF      =  MXNODES + LNODES
         MCOEFFSTOP  =  MCOEFF + LCOEFF
         MCOEFFSSIDE =  MCOEFFSTOP + LCOEFFSIDE
         MCOEFFDTOP  =  MCOEFFSSIDE + LCOEFFSIDE
         MCOEFFDSIDE =  MCOEFFDTOP + LCOEFFSIDE
         MFTOPS      =  MCOEFFDSIDE + LF
         MFTOPD      =  MFTOPS + LF
         MFSIDED     =  MFTOPD + LF
         MFSIDES     =  MFSIDED + LF

         ITOT = MFSIDES + LF

         IF ( ITOT .GE. LENW ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the real array'
            WRITE(*,*)'to be ',lenw
            WRITE(*,*)'The workspace needed is: ',itot
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF


      ELSEIF(IPERIOD .EQ. 7 .OR. IPERIOD .EQ. 8 .OR.
     1       IPERIOD .EQ. 9 .OR. IPERIOD .EQ. 10 .OR.
     2       IPERIOD .EQ. 11)THEN

C        First break up the integer array:
C        (In this case, we need to copy the
C        initial data structure over again so
C        that there is not a problem with
C        altering the structure on input)
         LCOLLEAG     =  9 * (4*NBOXES  + 1)
         LLEVELBOX    =  4 * NBOXES + 1
         LPARENTBOX   =  4 * NBOXES + 1
         LCOLBOX      =  4 * NBOXES + 1
         LROWBOX      =  4 * NBOXES + 1
         LCHILDBOX    =  4 * (4 * NBOXES + 1)
         LSTARTLEV    =  NLEV + 2
         LNBLEVEL     =  NLEV + 2
         LBOXLEV      =  4*(NBOXES + 1)
         LFLAG        =  4*(NBOXES + 1)
         LLOCALONOFF  =  4*(NBOXES + 1)
         LMULTONOFF   =  4*(NBOXES + 1)

         MCOLLEAG   =  1
         MLEVELBOX  =  MCOLLEAG + LCOLLEAG
         MLOCALONOFF =  MLEVELBOX + LLEVELBOX
         MPARENTBOX  =  MLOCALONOFF + LLOCALONOFF
         MROWBOX    =  MPARENTBOX + LPARENTBOX
         MCOLBOX    =  MROWBOX + LROWBOX
         MCHILDBOX  =  MCOLBOX + LCOLBOX
         MSTARTLEV  =  MCHILDBOX + LCHILDBOX
         MNBLEVEL   =  MSTARTLEV + LSTARTLEV
         MBOXLEV    =  MNBLEVEL + LNBLEVEL
         MFLAGEAST  =  MBOXLEV + LBOXLEV
         MFLAGWEST  =  MFLAGEAST + LFLAG
         MFLAGNORTH =  MFLAGWEST + LFLAG
         MFLAGSOUTH =  MFLAGNORTH + LFLAG
         MDOUBLETOPONOFF  =  MFLAGSOUTH + LFLAG
         MDOUBLESIDEONOFF =  MDOUBLETOPONOFF + LMULTONOFF
         MSINGLETOPONOFF  =  MDOUBLESIDEONOFF + LMULTONOFF
         MSINGLESIDEONOFF =  MSINGLETOPONOFF + LMULTONOFF

         ITOT = MSINGLESIDEONOFF + LMULTONOFF

         IF ( ITOT .GE. ILEN ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the integer array'
            WRITE(*,*)'to be ',ilen
            WRITE(*,*)'The workspace needed is: ',itot 
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF

C
C        Now carve up the real workspace:
C
         LMPOLE     =  2*(NTERMS+1)*(4*NBOXES + 1)
         LLOCXP     =  2*(NTERMS+1)*(4*NBOXES + 1)
         LEXP       =  2*(NNODES+1)*(4*NBOXES + 1)
         LCOEFF     =  64*(4*NBOXES + 1)
         LNODES     =  NNODES
         LCOMP      =  NNODES*(NTERMS+1)
         LCHOOSE    =  2*NTERMS * 2*NTERMS
         LZS        =  2 * 7 * 7 * (NNODES + 1)
         LCOEFFSIDE =  8 * (4*NBOXES + 1)
         LF         =  8 * (4*NBOXES + 1)
         LVOLMAPS   =  2*(NNODES + 1) * 36
         LWINT      =  2*(NTERMS + 1) * 8 * 8
         LSIDEMAPS  =  2*(NNODES + 1) * 8

         MMAPSOUTH   =  1
         MMAPNORTH   =  MMAPSOUTH   + LVOLMAPS
         MMAPEAST    =  MMAPNORTH   + LVOLMAPS
         MMAPWEST    =  MMAPEAST    + LVOLMAPS
         MWINT       =  MMAPWEST    + LVOLMAPS
         MEDBLETOP   =  MWINT       + LWINT
         MEDBLESIDE  =  MEDBLETOP   + LSIDEMAPS
         MESNGLETOP  =  MEDBLESIDE  + LSIDEMAPS
         MESNGLESIDE =  MESNGLETOP  + LSIDEMAPS
         MWDBLETOP   =  MESNGLESIDE + LSIDEMAPS
         MWDBLESIDE  =  MWDBLETOP   + LSIDEMAPS
         MWSNGLETOP  =  MWDBLESIDE  + LSIDEMAPS
         MWSNGLESIDE =  MWSNGLETOP  + LSIDEMAPS
         MNDBLETOP   =  MWSNGLESIDE + LSIDEMAPS
         MNDBLESIDE  =  MNDBLETOP   + LSIDEMAPS
         MNSNGLETOP  =  MNDBLESIDE  + LSIDEMAPS
         MNSNGLESIDE =  MNSNGLETOP  + LSIDEMAPS
         MSDBLETOP   =  MNSNGLESIDE + LSIDEMAPS
         MSDBLESIDE  =  MSDBLETOP   + LSIDEMAPS
         MSSNGLETOP  =  MSDBLESIDE  + LSIDEMAPS
         MSSNGLESIDE =  MSSNGLETOP  + LSIDEMAPS

         MPOLE       =  1
         LOCEXP      =  MPOLE   + LMPOLE
         IEXPN       =  LOCEXP  + LLOCXP
         IEXPS       =  IEXPN   + LEXP
         IEXPE       =  IEXPS   + LEXP
         IEXPW       =  IEXPE   + LEXP
         IEXPNBIG    =  IEXPW   + LEXP
         IEXPSBIG    =  IEXPNBIG   + LEXP
         IEXPEBIG    =  IEXPSBIG   + LEXP
         IEXPWBIG    =  IEXPEBIG   + LEXP
         MCOMP       =  IEXPWBIG   + LEXP
         MCHOOSE     =  MCOMP   + LCOMP
         MZS         =  MCHOOSE + LCHOOSE
         MZSEAST     =  MZS + LZS
         MZSWEST     =  MZSEAST  + LZS
         MZSNORTH    =  MZSWEST  + LZS
         MZSSOUTH    =  MZSNORTH + LZS
         MTEMP       =  MZSSOUTH + LZS
         MWNODES     =  MTEMP   + LNODES
         MXNODES     =  MWNODES + LNODES
         MCOEFF      =  MXNODES + LNODES
         MCOEFFSTOP  =  MCOEFF  + LCOEFF
         MCOEFFSSIDE =  MCOEFFSTOP  + LCOEFFSIDE
         MCOEFFDTOP  =  MCOEFFSSIDE + LCOEFFSIDE
         MCOEFFDSIDE =  MCOEFFDTOP  + LCOEFFSIDE
         MFTOPS      =  MCOEFFDSIDE + LCOEFFSIDE
         MFTOPD      =  MFTOPS  + LF
         MFSIDED     =  MFTOPD  + LF
         MFSIDES     =  MFSIDED + LF


         ITOT = MFSIDES + LF


         IF ( ITOT .GE. LENW ) THEN
            WRITE(*,*)'The total workspace needed exceeds the'
            WRITE(*,*)'amount allotted.'
            WRITE(*,*)'You have set the length of the real array'
            WRITE(*,*)'to be ',lenw
            WRITE(*,*)'The workspace needed is: ',itot
            WRITE(*,*)'I am stopping.'
            STOP
         ENDIF

      ENDIF



      IF(IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN

        CALL ADAPFMM8(NLEV,NDEG,WORK(MXNODES),WORK(MWNODES),
     1     WORK(MTEMP),WORK(MZS),WORK(MZSEAST),WORK(MZSWEST),
     1     WORK(MZSNORTH),WORK(MZSSOUTH),
     2     WORK(MCOMP),NTERMS,NNODES,POT,
     2     WORK(MPOLE),WORK(LOCEXP),
     3     WORK(IEXPN),WORK(IEXPS),WORK(IEXPE),WORK(IEXPW),
     3     WORK(IEXPNBIG),WORK(IEXPSBIG),WORK(IEXPEBIG),WORK(IEXPWBIG),
     4     WORK(MCHOOSE), WORK(MCOEFF), 
     5     LEVELBOX,IPARENTBOX,ICHILDBOX,
     6     ICOLBOX,IROWBOX,IWORK(MCOLLEAG),NBOXES,
     7     NBLEVEL,IBOXLEV,ISTARTLEV,IPERIOD,
     8     FRIGHT,MAP(MMAPNORTH), MAP(MMAPSOUTH), MAP(MMAPEAST),
     9     MAP(MMAPWEST), MAP(MWINT),IWORK(MFLAGEAST),
     1     IWORK(MFLAGNORTH),IWORK(MFLAGWEST),IWORK(MFLAGSOUTH),
     2     IWORK(MLOCALONOFF))



      ELSEIF(IPERIOD .EQ. 2 .OR. IPERIOD .EQ. 3 .OR.
     1       IPERIOD .EQ. 4 .OR. IPERIOD .EQ. 5 .OR.
     2       IPERIOD .EQ. 6 .OR. IPERIOD .EQ. 7 .OR.
     3       IPERIOD .EQ. 8 .OR. IPERIOD .EQ. 9 .OR.
     4       IPERIOD .EQ. 10 .OR. IPERIOD .EQ. 11)THEN

         CALL COPY(LEVELBOX,IWORK(MLEVELBOX),ICOLBOX,IWORK(MCOLBOX),
     1     IROWBOX,IWORK(MROWBOX),IPARENTBOX, IWORK(MPARENTBOX),
     2     ICHILDBOX,IWORK(MCHILDBOX),NBOXES,ISTARTLEV,IWORK(MSTARTLEV),
     3     IBOXLEV, IWORK(MBOXLEV), NBLEVEL, IWORK(MNBLEVEL),
     4     NLEV, NLEV1, NBOXES1)

         CALL FOLDOVER(IWORK(MLEVELBOX),IWORK(MCOLBOX),
     1     IWORK(MROWBOX),IWORK(MPARENTBOX),
     2     IWORK(MCHILDBOX),NBOXES1, IPERIOD, NLEV1, FRIGHT,
     3     IWORK(MSTARTLEV),IWORK(MBOXLEV),IWORK(MNBLEVEL))

         CALL SORTBOXES(IWORK(MLEVELBOX),NBOXES1,
     1     NLEV1, IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV))

         CALL MERGE(IWORK(MLEVELBOX),NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX), IWORK(MPARENTBOX),
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV))


C       Set IPERIODTEMP to the appropriate
C       call for the ADAPFMM6 subroutine.
        IPERIODTEMP = 2

        CALL ADAPFMM8(NLEV1,NDEG,WORK(MXNODES),WORK(MWNODES), 
     1     WORK(MTEMP),WORK(MZS),WORK(MZSEAST),WORK(MZSWEST),
     1     WORK(MZSNORTH),WORK(MZSSOUTH),
     2     WORK(MCOMP), 
     2     NTERMS,NNODES,POT, WORK(MPOLE),WORK(LOCEXP),
     3     WORK(IEXPN),WORK(IEXPS),WORK(IEXPE),WORK(IEXPW),
     3     WORK(IEXPNBIG),WORK(IEXPSBIG),WORK(IEXPEBIG),WORK(IEXPWBIG),
     4     WORK(MCHOOSE), WORK(MCOEFF),
     5     IWORK(MLEVELBOX),IWORK(MPARENTBOX),IWORK(MCHILDBOX),
     6     IWORK(MCOLBOX),IWORK(MROWBOX),IWORK(MCOLLEAG),NBOXES1,
     7     IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     8     IPERIODTEMP,FRIGHT,MAP(MMAPNORTH), MAP(MMAPSOUTH),
     9     MAP(MMAPEAST),  MAP(MMAPWEST), MAP(MWINT),
     1     IWORK(MFLAGEAST), IWORK(MFLAGNORTH),
     2     IWORK(MFLAGWEST), IWORK(MFLAGSOUTH),
     3     IWORK(MLOCALONOFF))


        IF(IPERIOD .EQ. 7 .OR. IPERIOD .EQ. 8 .OR.
     1     IPERIOD .EQ. 9 .OR. IPERIOD .EQ. 10 .OR.
     2     IPERIOD .EQ. 11)THEN

        IF(IPERIOD .EQ. 7)THEN
         CALL SETBCS8DIR(NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX),
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     4           WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     5           IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     6           IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF),
     7           HL, HR, HB, HT)

        ELSEIF(IPERIOD .EQ. 8)THEN
         CALL SETBCS8NEU(NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX), 
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     4           WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     5           IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     6           IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF),
     7           HL, HR, HB, HT)

        ELSEIF(IPERIOD .EQ. 9)THEN
         CALL SETBCS8DIRPER(NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX),
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     4           WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     5           IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     6           IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF),
     7           HL, HR)

        ELSEIF(IPERIOD .EQ. 10)THEN
         CALL SETBCS8DIRNEU(NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX),
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     4           WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     5           IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     6           IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF),
     7           HL, HR, HB, HT)

        ELSEIF(IPERIOD .EQ. 11)THEN
         CALL SETBCS8NEUPER(NBOXES1,NLEV1,
     1           IWORK(MCOLBOX),IWORK(MROWBOX),
     2           IWORK(MCHILDBOX),
     3           IWORK(MNBLEVEL),IWORK(MBOXLEV),IWORK(MSTARTLEV),
     4           WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     5           IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     6           IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF),
     7           HL, HR)
        ENDIF


        CALL BOUNDFMM8(NLEV1,NDEG,WORK(MXNODES),WORK(MWNODES),
     1     WORK(MTEMP),WORK(MZS),WORK(MZSEAST),WORK(MZSWEST),
     1     WORK(MZSNORTH),WORK(MZSSOUTH), WORK(MCOMP),NTERMS,NNODES,POT,
     2     WORK(MPOLE),WORK(LOCEXP),WORK(IEXPN),WORK(IEXPS),
     3     WORK(IEXPE),WORK(IEXPW),
     2     WORK(IEXPNBIG),WORK(IEXPSBIG),
     3     WORK(IEXPEBIG),WORK(IEXPWBIG),
     3     WORK(MCHOOSE),
     4     IWORK(MLEVELBOX),IWORK(MPARENTBOX),
     5     IWORK(MCHILDBOX),IWORK(MCOLBOX),IWORK(MROWBOX),
     6     IWORK(MCOLLEAG),NBOXES1,IWORK(MNBLEVEL),IWORK(MBOXLEV),
     7     IWORK(MSTARTLEV),WORK(MCOEFFSTOP),
     8     WORK(MCOEFFSSIDE),WORK(MCOEFFDTOP),WORK(MCOEFFDSIDE),
     9     WORK(MFTOPS),WORK(MFSIDES),WORK(MFTOPD),WORK(MFSIDED),
     1     MAP(MEDBLETOP), MAP(MEDBLESIDE),  MAP(MESNGLETOP),
     2     MAP(MESNGLESIDE), MAP(MWDBLETOP), MAP(MWDBLESIDE),
     3     MAP(MWSNGLETOP), MAP(MWSNGLESIDE), MAP(MNDBLETOP),
     4     MAP(MNDBLESIDE), MAP(MNSNGLETOP),  MAP(MNSNGLESIDE),
     5     MAP(MSDBLETOP), MAP(MSDBLESIDE), MAP(MSSNGLETOP),
     6     MAP(MSSNGLESIDE),IWORK(MFLAGEAST),IWORK(MFLAGNORTH),
     7     IWORK(MFLAGWEST),IWORK(MFLAGSOUTH),IWORK(MLOCALONOFF),
     8     IWORK(MDOUBLETOPONOFF), IWORK(MDOUBLESIDEONOFF),
     9     IWORK(MSINGLETOPONOFF), IWORK(MSINGLESIDEONOFF))

        ENDIF

       ENDIF
      RETURN
      END


C***********************************************************************
C     The following subroutine defines the array that represents the 
C     right hand side of the Poisson equation. 
C     In each childless box, there are 16 cell centered points 
C     where the right hand side values (and later
C     the output values) are defined.  
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NLEV is the finest level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     NBLEVEL is the total number of boxes per level
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     H is the real function that is the right hand side
C       of the poisson equation
C
C     OUTPUT:
C
C     FRIGHT is the right hand side defined on the tree
C
C***********************************************************************
      SUBROUTINE SETF8(FRIGHT,
     1       ICOLBOX, IROWBOX, ICHILDBOX,NLEV,
     2       NBLEVEL, IBOXLEV, ISTARTLEV, H)
      IMPLICIT NONE
C-----Global variables
      INTEGER  ICOLBOX(1), IROWBOX(1), NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  FRIGHT(64,1), XF(64), YF(64)
C-----Local variables
      INTEGER  I, IBOX, J, K, L
      REAL *8  TEMP1, PI16
      REAL *8  XSTART
      REAL *8  XSHIFT, XX(8)
      REAL *8  XSCALE(8)
      REAL *8  YSHIFT
C-----External functions
      REAL *8  H

      PI16 = DACOS(-1.0D0) / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0

      TEMP1 = 1.0D0
      DO K = 0, NLEV
      XSTART = (1.0D0 - TEMP1) / 2.0D0

      XSCALE(1) = XX(1) * TEMP1 - XSTART
      XSCALE(2) = XX(2) * TEMP1 - XSTART
      XSCALE(3) = XX(3) * TEMP1 - XSTART
      XSCALE(4) = XX(4) * TEMP1 - XSTART
      XSCALE(5) = XX(5) * TEMP1 - XSTART
      XSCALE(6) = XX(6) * TEMP1 - XSTART
      XSCALE(7) = XX(7) * TEMP1 - XSTART
      XSCALE(8) = XX(8) * TEMP1 - XSTART


      DO 100 I = ISTARTLEV(K), ISTARTLEV(K) + NBLEVEL(K) - 1
        IBOX = IBOXLEV(I)
        IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100

        XSHIFT  =  DBLE(ICOLBOX(IBOX) - 1) * TEMP1
        YSHIFT  =  DBLE(IROWBOX(IBOX) - 1) * TEMP1

        DO J = 1, 8
          DO L = 1, 8
            XF(8*(L-1)+J) = XSCALE(J) + XSHIFT
            YF(8*(J-1)+L) = XSCALE(J) + YSHIFT
          END DO
        END DO
       
        DO J = 1, 64
          FRIGHT(J,IBOX) = H(XF(J),YF(J))
        END DO
100   CONTINUE
      TEMP1 = TEMP1 / 2.0D0
      END DO
      RETURN
      END



C***********************************************************************
C     The following subroutine defines the array COEFFS that contains
C     the polynomial coefficients for the polynomial that approximates
C     the right hand side of the Poisson equation.
C
C     INPUT:
C
C     FRIGHT is the right hand side defined on the old grid
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICHILDBOX denotes the four children of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins
C               in the IBOXLEV array
C
C     IPERIOD denotes which of the boundary cases we are in
C
C     A is the matrix that maps from 36 function
C          values to 21 polynomial coefficients
C
C     OUTPUT:
C  
C     COEFFS is the array of coefficients for the basis functions
C
C***********************************************************************
      SUBROUTINE MKCOEFFS8(COEFFS,FRIGHT,
     1         NLEV, ICHILDBOX, NBLEVEL, 
     2         IBOXLEV, ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  COEFFS(0:7,0:7,1), FRIGHT(64,1)
C-----Local variables
      INTEGER  I, IBOX, J, L, II
      REAL *8  FTEMP(8,8)
      REAL *8  COEFFTEMP(8,8)
      REAL *8  WSAVE2(1000)

      DO L = 0, NLEV
       DO I = ISTARTLEV(L), ISTARTLEV(L) + NBLEVEL(L) - 1
       IBOX = IBOXLEV(I)
        IF (ICHILDBOX(1,IBOX) .LT.0)THEN


C       initialize WSAVE2 array (needed by GETCOEFF).
        CALL CHXCIN(8,WSAVE2)

        DO J = 1, 8
          DO II = 1, 8
            FTEMP(J,II) = FRIGHT(8*(II-1)+J,IBOX)
          END DO
        END DO


C       compute Chebyshev transforms
        CALL GETCOEFF(8,8,FTEMP,COEFFTEMP,WSAVE2)

C       Now set these values to out coefficients
        DO J = 0, 7
          DO II = 0, 7
            COEFFS(J,II,IBOX) = COEFFTEMP(J+1,II+1)
          END DO
        END DO
        ENDIF
       END DO
      END DO
      RETURN
      END


      SUBROUTINE GETCOEFF(N,M,FDAT,COEFF,WSAVE2)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 FDAT(N,M),COEFF(N,M)
      REAL *8 WSAVE2(1000)
      REAL *8 WORK(1000)
      REAL *8 F(1000),TEXP(1000)
C
C     transform rows
C
      DO J = 1,M
         DO I = 1,N
            F(I) = FDAT(I,J)
         ENDDO
         CALL CHEXFC(F,N,TEXP,WSAVE2,WORK)
         DO I = 1,N
            COEFF(I,J) = TEXP(I)
         ENDDO
      ENDDO
C
C     transform columns
C
      DO I = 1,N
         DO J = 1,M
            F(J) = COEFF(I,J)
         ENDDO
         CALL CHEXFC(F,N,TEXP,WSAVE2,WORK)
         DO J = 1,M
            COEFF(I,J) = TEXP(J)
         ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE GETCOEFF1(FDAT,COEFF,WSAVE2)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 FDAT(8,1),COEFF(8,1)
      REAL *8 WSAVE2(1000)
      REAL *8 WORK(1000)
      REAL *8 F(1000),TEXP(1000)
C
C     transform rows
C
         DO I = 1,8
            F(I) = FDAT(I,1)
         ENDDO
         CALL CHEXFC(F,8,TEXP,WSAVE2,WORK)
         DO I = 1,8
            COEFF(I,1) = TEXP(I)
         ENDDO
      RETURN
      END


C**********************************************************************
      SUBROUTINE ADDEXP(B,A,NTERMS)
C**********************************************************************
C     INPUT :     Multipole Expansions  A and B of length NTERMS
C     OUTPUT:     A is over written by (A+B).
C***********************************************************************
      IMPLICIT NONE
      INTEGER  NTERMS
      COMPLEX *16 B(0:1),A(0:1)
      INTEGER I
C --------------------------------------------------------------------
C
      DO I = 0,NTERMS
         A(I) = A(I) + B(I)
      END DO
      RETURN
      END

 
C********************************************************************
C     The main subroutine of multipole algorithm. Two passes are
C     executed. In the first (upward) pass, multipole expansions for
C     all boxes at all levels are computed. In the second (downward) 
C     pass, interactions are computed at successively finer levels.
C
C     INPUT:
C
C     NLEV is the finest level
C
C     NDEG is the degree of the approximating polynomial
C
C     XNODES is a blank array that is set to 
C            the nodes in the plane wave expansions
C
C     WNODES is a blank array that is set to 
C            the weights in the plane wave expansions
C
C     TEMP  is a blank array that is set to 
C            the ratio of the weights and nodes
C
C     ZS represents the shifts for the exponential expansions
C
C     COMP is a blank array
C
C     NTERMS is the order of the multipole expansions
C
C     NNODES is the order of the plane wave expansions
C
C     MPOLE is the array that the multipole expansions are stored in
C
C     LOCEXP is the array that the local expansions are stored in
C
C     EXPN is the array that the north plane wave expansions are stored in
C
C     EXPS is the array that the south plane wave expansions are stored in
C
C     EXPE is the array that the east plane wave expansions are stored in
C
C     EXPW is the array that the west plane wave expansions are stored in
C
C     C is the array that the binomial coefficients are stored in
C
C     COEFFS is the array that the basis function coefficients are stored in
C
C     LEVELBOX is an array determining the level of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICOLLEAGBOX denotes the colleagues of a given box
C
C     NBOXES is the total number of boxes
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     IPERIOD denotes which of the boundary cases we are in
C
C     FRIGHT is the right hand side defined on the old grid
C
C     MAPNORTH is the array containing the north BTOSFAR map
C
C     MAPSOUTH is the array containing the south BTOSFAR map
C
C     MAPEAST is the array containing the east BTOSFAR map
C
C     MAPWEST is the array containing the west BTOSFAR map
C
C     WINT is the array containing the polynomial to multipole weights
C
C     OUTPUT:
C  
C     POT represents the solution (defined on the same tree as the
C         right hand side on input)
C
C********************************************************************
      SUBROUTINE ADAPFMM8(NLEV,NDEG,XNODES,WNODES,TEMP, 
     1         ZS,ZSEAST,ZSWEST,ZSNORTH,ZSSOUTH,
     2         COMP,NTERMS,NNODES,POT,
     3         MPOLE,LOCEXP,EXPN,EXPS,EXPE,EXPW,
     4         EXPNBIG,EXPSBIG,EXPEBIG,EXPWBIG,
     5         C,COEFFS,LEVELBOX,IPARENTBOX,ICHILDBOX,
     6         ICOLBOX,IROWBOX,ICOLLEAGBOX,NBOXES,
     7         NBLEVEL,IBOXLEV,ISTARTLEV,IPERIOD,
     8         FRIGHT,MAPNORTH, MAPSOUTH, MAPEAST,
     9         MAPWEST, WINT,IFLAGEAST,IFLAGNORTH,
     1         IFLAGWEST,IFLAGSOUTH,LOCALONOFF)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NLEV,NTERMS
      INTEGER  NDEG
      INTEGER  NNODES, NBOXES
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      INTEGER  LEVELBOX(1), IPARENTBOX(1)
      INTEGER  ICHILDBOX(4,1), ICOLLEAGBOX(9,1)
      INTEGER  IFLAGEAST(1), IFLAGNORTH(1)
      INTEGER  IFLAGWEST(1), IFLAGSOUTH(1)
      INTEGER  LOCALONOFF(1)
      REAL *8  FRIGHT(64,1)
      REAL *8  XNODES(1), WNODES(1)
      REAL *8  C(2*NTERMS,2*NTERMS)
      REAL *8  POT(64,1)
      REAL *8  COEFFS(0:NDEG,0:NDEG,1)
      REAL *8  COMP(NNODES,0:NTERMS)
      REAL *8  TEMP(1)
      COMPLEX *16  MPOLE(0:NTERMS,1),LOCEXP(0:NTERMS,1)
      COMPLEX *16  EXPN(0:NNODES,1), EXPS(0:NNODES,1)
      COMPLEX *16  EXPE(0:NNODES,1), EXPW(0:NNODES,1)
      COMPLEX *16  EXPNBIG(0:NNODES,1), EXPSBIG(0:NNODES,1)
      COMPLEX *16  EXPEBIG(0:NNODES,1), EXPWBIG(0:NNODES,1)
      COMPLEX *16  ZS(-3:3,-3:3,0:1)
      COMPLEX *16  ZSEAST(-3:3,-3:3,0:1)
      COMPLEX *16  ZSWEST(-3:3,-3:3,0:1)
      COMPLEX *16  ZSNORTH(-3:3,-3:3,0:1)
      COMPLEX *16  ZSSOUTH(-3:3,-3:3,0:1)
      COMPLEX *16  WINT(0:NTERMS,0:7,0:7)
C-----Local variables
      INTEGER  IFAR1, IFAR2, IFAR3 
      INTEGER  ICLOSE1, ICLOSE2, IPERIOD
      INTEGER  I, J, IP
      INTEGER  IOUT, II, JJ, IER
      INTEGER  IC1,IC2,IC3,IC4
      INTEGER  INALL(6),IYNALL(6),IN12(4),IY12(4)
      INTEGER  ISALL(6),IYSALL(6),IS34(4),IY34(4)
      INTEGER  NNALL,NN12,NSALL,NS34
      INTEGER  IEALL(4),IYEALL(4),IE13(2),IY13(2)
      INTEGER  IE1(1),IY1(1),IE3(1),IY3(1)
      INTEGER  IWALL(4),IYWALL(4),IW24(2),IY24(2)
      INTEGER  IW2(1),IY2(1),IW4(1),IY4(1)
      INTEGER  NEALL,NE13,NE1,NE3,NWALL,NW24,NW2,NW4
      INTEGER  NSIDE, NSIDEMARK
      INTEGER  INBIG12(3), ISBIG34(3)
      INTEGER  IEBIG13(1), IWBIG24(1)
      INTEGER  IEBIG1(1), IWBIG2(1)
      INTEGER  IEBIG3(1), IWBIG4(1)
      INTEGER  NB, ISTART, IEND
      REAL *8  XP(64),YP(64)
      REAL *8  XLENGTH
      REAL *8  T(3), SCALE(0:100), TEMP1, TLOG(36), TLOGS(36)
      REAL *8  SUM, ZERO
      REAL *8  SCALETEMP
      REAL *8  WBTOS(12,64,36)
      REAL *8  WSTOB(12,64,36)
      REAL *8  W(9,64,36), X
      REAL *8  PI
      COMPLEX *16  ZSHIFT, ZPOT
      COMPLEX *16  FTARGET1, FTARGET2, FTARGET3, IMAG
      COMPLEX *16  B(0:50)
      COMPLEX *16  MEXNALL(0:50)
      COMPLEX *16  MEXN12(0:50)
      COMPLEX *16  MEXSALL(0:50)
      COMPLEX *16  MEXS34(0:50)
      COMPLEX *16  MEXEALL(0:50)
      COMPLEX *16  MEXE13(0:50)
      COMPLEX *16  MEXE1(0:50)
      COMPLEX *16  MEXE3(0:50)
      COMPLEX *16  MEXWALL(0:50)
      COMPLEX *16  MEXW24(0:50)
      COMPLEX *16  MEXW2(0:50)
      COMPLEX *16  MEXW4(0:50)
      COMPLEX *16  SPIN
      COMPLEX *16  EXPNALL(0:60), EXPSALL(0:60)
      COMPLEX *16  EXPEALL(0:60), EXPWALL(0:60)
      COMPLEX *16  MAPNORTH(0:NNODES,36), MAPWEST(0:NNODES,36)
      COMPLEX *16  MAPSOUTH(0:NNODES,36), MAPEAST(0:NNODES,36)
      COMPLEX *16  TEMPEAST(64,40), TEMPNORTH(64,40)
      COMPLEX *16  TEMPWEST(64,40), TEMPSOUTH(64,40)
      COMPLEX *16  TEMPSHIFTWEST, TEMPSHIFTNORTH
      COMPLEX *16  TEMPSHIFTEAST, TEMPSHIFTSOUTH
      DATA ZERO/0.0D0/
      DATA IMAG/(0.0D0,1.0D0)/
      PI = DACOS(-1.0D0)

      IF(IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN
        SCALETEMP = 1.0D0
      ELSEIF(IPERIOD .EQ. 2)THEN 
        SCALETEMP = 2.0D0
      ENDIF

C     Let's initially set the potential to zero:
      DO II = 1, NBOXES
       DO JJ = 1, 64
         POT(JJ,II) = 0.0D0
       END DO
      END DO

      DO II = 1, NBOXES
        IFLAGEAST(II) = 0
        IFLAGWEST(II) = 0
        IFLAGNORTH(II) = 0
        IFLAGSOUTH(II) = 0
        LOCALONOFF(II) = 0
      END DO

C     Initialize the LOCALONOFF switch to the correct values
      IF(IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN
       DO II = 1, NBOXES
          LOCALONOFF(II) = 1
       END DO
      ELSEIF(IPERIOD .EQ. 2)THEN
       NSIDE = 1
       DO JJ = 0, NLEV
         ISTART = ISTARTLEV(JJ)
         IEND = ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
         IF(NSIDE .EQ. 1)THEN
           NSIDEMARK = 1
         ELSE
           NSIDEMARK = NSIDE / 2
         ENDIF
         DO II = ISTART, IEND
           I = IBOXLEV(II)
           IF(IROWBOX(I) .LE. NSIDEMARK .AND.
     1        ICOLBOX(I) .LE. NSIDEMARK)THEN
             LOCALONOFF(I) = 1
           ENDIF
         END DO
        NSIDE = 2 * NSIDE
       END DO
      ENDIF


C     Call a set of routines that will precompute various 
C     things needed later in the adapfmm, mkcoeffs, and
C     mkshifts routines.
      CALL PWTS4(XNODES, WNODES, NNODES)

      CALL PRECOMPUTE(COMP,NTERMS,NNODES,
     1  XNODES,WNODES, TEMP,C, SUM)

      DO I = 1, 8
        XP((I-1)*8+1) = DCOS(15D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+2) = DCOS(13D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+3) = DCOS(11D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+4) = DCOS(9D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+5) = DCOS(7D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+6) = DCOS(5D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+7) = DCOS(3D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+8) = DCOS(1D0*PI/16.0D0) / 2.0D0

        YP(I) =   DCOS(15D0*PI/16.0D0) / 2.0D0
        YP(8+I) = DCOS(13D0*PI/16.0D0) / 2.0D0
        YP(16+I) = DCOS(11D0*PI/16.0D0)/ 2.0D0
        YP(24+I) = DCOS(9D0*PI/16.0D0) / 2.0D0
        YP(32+I) = DCOS(7D0*PI/16.0D0) / 2.0D0
        YP(40+I) = DCOS(5D0*PI/16.0D0) / 2.0D0
        YP(48+I) = DCOS(3D0*PI/16.0D0) / 2.0D0
        YP(56+I) = DCOS(1D0*PI/16.0D0) / 2.0D0
      ENDDO

      DO JJ = 1, 64
        TEMPSHIFTWEST = 2.0D0*(XP(JJ) + IMAG*YP(JJ))
        TEMPSHIFTEAST = -TEMPSHIFTWEST
        TEMPSHIFTNORTH = IMAG*TEMPSHIFTWEST
        TEMPSHIFTSOUTH = -TEMPSHIFTNORTH
        DO II = 1, NNODES
          TEMPEAST(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTEAST)
          TEMPWEST(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTWEST)
          TEMPNORTH(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTNORTH)
          TEMPSOUTH(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTSOUTH)
        END DO
      END DO


C     Call a routine that will generate the array ZS
C     which will be needed later on when shifting the
C     exponential expansions.
      CALL MKSHIFTS2D(XNODES,NNODES,ZS,
     1    ZSEAST,ZSWEST,ZSNORTH,ZSSOUTH)


      CALL MKCOEFFS8(COEFFS,FRIGHT,NLEV,
     1    ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV)

      CALL MKCOLLS(ICOLBOX,
     1    IROWBOX,ICOLLEAGBOX,NBOXES,NLEV,
     2    IPARENTBOX,ICHILDBOX,NBLEVEL,
     3    IBOXLEV, ISTARTLEV,IPERIOD)

C
C     create SCALE array
C
      X = 1.0D0 / SCALETEMP
      DO I = 0,NLEV
       SCALE(I) = X
       X = X*2.0D0
      ENDDO

C     Initialize the arrays of weights for all of the local coefficients.  W is 
C     the array for local interactions between boxes of the same size, 
C     W36ADAPBTOS is
C     the array for local interactions going from the big to small case, and 
C     W36ADAPSTOB is the array for local interactions going from the small to
C     big case.
      CALL W36ADAPBTOS(WBTOS)
      CALL W36ADAPSTOB(WSTOB)
      CALL WEIGHTS36(W)

C**********************************************************************
C
C     =============================================
C     UPWARD PASS
C     =============================================
C
C-----initialize multipole and local expansions to zero.
C
      DO I = 1, NBOXES
       DO J = 0, NTERMS
         MPOLE(J,I) = ZERO 
         LOCEXP(J,I) = ZERO 
       END DO
       DO J = 0, NNODES
         EXPE(J,I) = ZERO
         EXPW(J,I) = ZERO
         EXPN(J,I) = ZERO
         EXPS(J,I) = ZERO
         EXPEBIG(J,I) = ZERO
         EXPWBIG(J,I) = ZERO
         EXPNBIG(J,I) = ZERO
         EXPSBIG(J,I) = ZERO
       END DO
      END DO


      
      DO I = NLEV, 0, -1
         XLENGTH = 1.0D0/SCALE(I)
         ISTART = ISTARTLEV(I)
         IEND = ISTART + NBLEVEL(I) - 1
         DO II = ISTART,IEND
            J = IBOXLEV(II)
            IF(ICHILDBOX(1,J) .LT. 0)THEN

             CALL MULTIPOLE(COEFFS(0,0,J),NDEG,NTERMS,XLENGTH,
     1                                      MPOLE(0,J), WINT)
            ELSEIF(ICHILDBOX(1,J) .GT. 0)THEN

             IC1 = ICHILDBOX(1,J)
             IC2 = ICHILDBOX(2,J)
             IC3 = ICHILDBOX(3,J)
             IC4 = ICHILDBOX(4,J)

             CALL CHILDPAR(MPOLE(0,J),
     1         MPOLE(0,IC1),MPOLE(0,IC2),
     2         MPOLE(0,IC3),MPOLE(0,IC4),
     3         NTERMS, C)
            ENDIF
         END DO
      END DO

C     For the periodic case, let's call a routine that 
C     will set the initial local expansion to the correct 
C     value to account for the periodicity in the far field terms:
      IF (IPERIOD .EQ. 1 .OR. IPERIOD .EQ. 2)THEN
C      Set IFLAG2 to zero to indicate that we are 
C      working with the potential and not the 
C      electric field.
        ISTART = ISTARTLEV(0)
        J = IBOXLEV(ISTART)
        CALL PSPPIN(IER)
        CALL PSPPTA(MPOLE(0,J),LOCEXP(0,J), NTERMS)
      ENDIF


C     =============================================
C     DOWNWARD PASS
C     =============================================
      XLENGTH = SCALETEMP
      DO I = 0,NLEV
C        First let's set up an array of scaling scale
C        factors needed in the local interaction:
         IF(IPERIOD .EQ. 0 .OR. IPERIOD .EQ. 1)THEN
          T(1)  =  4.0D0**(I+1)
          T(2)  =  1.0D0/2.0D0**(I+1)
          T(3)  =  4.0D0 * T(1)

          TEMP1 =  2.0D0 * DLOG(T(2))/PI
          TLOG(1)  =  TEMP1 
          TLOG(2)  =  0.0D0
          TLOG(3)  =  0.0D0
          TLOG(4)  = -TEMP1 / 3.0D0
          TLOG(5)  =  0.0D0
          TLOG(6)  = -TEMP1 / 3.0D0
          TLOG(7)  =  0.0D0
          TLOG(8)  =  0.0D0
          TLOG(9)  =  0.0D0
          TLOG(10) =  0.0D0
          TLOG(11) = -TEMP1 /15.0D0
          TLOG(12) =  0.0D0
          TLOG(13) =  TEMP1 /9.0D0
          TLOG(14) =  0.0D0
          TLOG(15) = -TEMP1 /15.0D0
          TLOG(16) =  0.0D0
          TLOG(17) =  0.0D0
          TLOG(18) =  0.0D0
          TLOG(19) =  0.0D0
          TLOG(20) =  0.0D0
          TLOG(21) =  0.0D0
          TLOG(22) = -TEMP1 /35.0D0
          TLOG(23) =  0.0D0
          TLOG(24) =  TEMP1 /45.0D0
          TLOG(25) =  0.0D0
          TLOG(26) =  TEMP1 /45.0D0
          TLOG(27) =  0.0D0
          TLOG(28) = -TEMP1 /35.0D0
          TLOG(29) =  0.0D0
          TLOG(30) =  0.0D0
          TLOG(31) =  0.0D0
          TLOG(32) =  0.0D0
          TLOG(33) =  0.0D0
          TLOG(34) =  0.0D0
          TLOG(35) =  0.0D0
          TLOG(36) =  0.0D0

          TEMP1 =  2.0D0 * DLOG(.50D0*T(2))/PI
          TLOGS(1)  =  TEMP1
          TLOGS(2)  =  0.0D0
          TLOGS(3)  =  0.0D0
          TLOGS(4)  = -TEMP1 / 3.0D0
          TLOGS(5)  =  0.0D0
          TLOGS(6)  = -TEMP1 / 3.0D0
          TLOGS(7)  =  0.0D0
          TLOGS(8)  =  0.0D0
          TLOGS(9)  =  0.0D0
          TLOGS(10) =  0.0D0
          TLOGS(11) = -TEMP1 /15.0D0
          TLOGS(12) =  0.0D0
          TLOGS(13) =  TEMP1 /9.0D0
          TLOGS(14) =  0.0D0
          TLOGS(15) = -TEMP1 /15.0D0
          TLOGS(16) =  0.0D0
          TLOGS(17) =  0.0D0
          TLOGS(18) =  0.0D0
          TLOGS(19) =  0.0D0
          TLOGS(20) =  0.0D0
          TLOGS(21) =  0.0D0
          TLOGS(22) = -TEMP1 /35.0D0
          TLOGS(23) =  0.0D0
          TLOGS(24) =  TEMP1 /45.0D0
          TLOGS(25) =  0.0D0
          TLOGS(26) =  TEMP1 /45.0D0
          TLOGS(27) =  0.0D0
          TLOGS(28) = -TEMP1 /35.0D0
          TLOGS(29) =  0.0D0
          TLOGS(30) =  0.0D0
          TLOGS(31) =  0.0D0
          TLOGS(32) =  0.0D0
          TLOGS(33) =  0.0D0
          TLOGS(34) =  0.0D0
          TLOGS(35) =  0.0D0
          TLOGS(36) =  0.0D0

         ELSEIF(IPERIOD .EQ. 2)THEN
          T(1)  =  4.0D0**I
          T(2)  =  1.0D0/2.0D0**I
          T(3)  =  4.0D0 * T(1)

          TEMP1 =  2.0D0 * DLOG(T(2))/PI
          TLOG(1)  =  TEMP1
          TLOG(2)  =  0.0D0
          TLOG(3)  =  0.0D0
          TLOG(4)  = -TEMP1 / 3.0D0
          TLOG(5)  =  0.0D0
          TLOG(6)  = -TEMP1 / 3.0D0
          TLOG(7)  =  0.0D0
          TLOG(8)  =  0.0D0
          TLOG(9)  =  0.0D0
          TLOG(10) =  0.0D0
          TLOG(11) = -TEMP1 / 15.0D0
          TLOG(12) =  0.0D0
          TLOG(13) =  TEMP1 / 9.0D0
          TLOG(14) =  0.0D0
          TLOG(15) = -TEMP1 / 15.0D0
          TLOG(16) =  0.0D0
          TLOG(17) =  0.0D0
          TLOG(18) =  0.0D0
          TLOG(19) =  0.0D0
          TLOG(20) =  0.0D0
          TLOG(21) =  0.0D0
          TLOG(22) = -TEMP1 / 35.0D0
          TLOG(23) =  0.0D0
          TLOG(24) =  TEMP1 / 45.0D0
          TLOG(25) =  0.0D0
          TLOG(26) =  TEMP1 / 45.0D0
          TLOG(27) =  0.0D0
          TLOG(28) = -TEMP1 / 35.0D0
          TLOG(29) =  0.0D0
          TLOG(30) =  0.0D0
          TLOG(31) =  0.0D0
          TLOG(32) =  0.0D0
          TLOG(33) =  0.0D0
          TLOG(34) =  0.0D0
          TLOG(35) =  0.0D0
          TLOG(36) =  0.0D0

          TEMP1 =  2.0D0 * DLOG(.50D0*T(2))/PI
          TLOGS(1)  =  TEMP1
          TLOGS(2)  =  0.0D0
          TLOGS(3)  =  0.0D0
          TLOGS(4)  = -TEMP1 / 3.0D0
          TLOGS(5)  =  0.0D0
          TLOGS(6)  = -TEMP1 / 3.0D0
          TLOGS(7)  =  0.0D0
          TLOGS(8)  =  0.0D0
          TLOGS(9)  =  0.0D0
          TLOGS(10) =  0.0D0
          TLOGS(11) = -TEMP1 / 15.0D0
          TLOGS(12) =  0.0D0
          TLOGS(13) =  TEMP1 / 9.0D0
          TLOGS(14) =  0.0D0
          TLOGS(15) = -TEMP1 / 15.0D0
          TLOGS(16) =  0.0D0
          TLOGS(17) =  0.0D0
          TLOGS(18) =  0.0D0
          TLOGS(19) =  0.0D0
          TLOGS(20) =  0.0D0
          TLOGS(21) =  0.0D0
          TLOGS(22) = -TEMP1 / 35.0D0
          TLOGS(23) =  0.0D0
          TLOGS(24) =  TEMP1 / 45.0D0
          TLOGS(25) =  0.0D0
          TLOGS(26) =  TEMP1 / 45.0D0
          TLOGS(27) =  0.0D0
          TLOGS(28) = -TEMP1 / 35.0D0
          TLOGS(29) =  0.0D0
          TLOGS(30) =  0.0D0
          TLOGS(31) =  0.0D0
          TLOGS(32) =  0.0D0
          TLOGS(33) =  0.0D0
          TLOGS(34) =  0.0D0
          TLOGS(35) =  0.0D0
          TLOGS(36) =  0.0D0

         ENDIF

         ISTART = ISTARTLEV(I)
         IEND = ISTART + NBLEVEL(I) - 1
         DO II = ISTART,IEND

           J = IBOXLEV(II)
           IF(ICHILDBOX(1,J) .GT. 0)THEN
C          If the box has children, do all of the work involving local 
C          expansions here.

            IC1 = ICHILDBOX(1,J)
            IC2 = ICHILDBOX(2,J)
            IC3 = ICHILDBOX(3,J)
            IC4 = ICHILDBOX(4,J)

            CALL PARENTCHILD(LOCEXP(0,J),LOCEXP(0,IC1),
     1       LOCEXP(0,IC2),LOCEXP(0,IC3),
     2       LOCEXP(0,IC4),NTERMS,C,LOCALONOFF(J))

            CALL MKEXP2D(J,NTERMS,MPOLE,
     1       NNODES,MEXNALL,MEXN12,MEXSALL,MEXS34,MEXEALL,
     2       MEXE13,MEXE1,MEXE3,MEXWALL,MEXW24,MEXW2,MEXW4,
     3       ZS,COMP,WNODES,SCALE(I+1),SUM,TEMP,
     4       ICHILDBOX)
 

            CALL MKLISTS(J,INALL,NNALL,IYNALL,
     1       IN12,NN12,IY12,
     2       ISALL,NSALL,IYSALL,IS34,NS34,IY34,
     3       IEALL,NEALL,IYEALL,IE13,NE13,IY13,
     4       IWALL,NWALL,IYWALL,IW24,NW24,IY24,
     5       IW2,IY2,NW2,IW4,IY4,NW4,
     6       IE1,IY1,NE1,IE3,IY3,NE3,
     7       INBIG12,ISBIG34,IEBIG13,IWBIG24,
     8       IEBIG1, IWBIG2, IEBIG3, IWBIG4,
     9       ICOLLEAGBOX,ICHILDBOX,
     1       ICOLBOX, IROWBOX, IPERIOD,
     2       IFLAGEAST, IFLAGWEST, IFLAGNORTH,
     3       IFLAGSOUTH, LOCALONOFF)


            CALL PROCESSNO(EXPN,INALL,NNALL,IYNALL,
     1       IN12,NN12,IY12,MEXNALL,MEXN12,ZS,NNODES,
     2       INBIG12,EXPNBIG,ZSNORTH,LOCALONOFF)
            CALL PROCESSSO(EXPS,ISALL,NSALL,IYSALL,
     1       IS34,NS34,IY34,MEXSALL,MEXS34,ZS,NNODES,
     2       ISBIG34,EXPSBIG,ZSSOUTH, LOCALONOFF)
            CALL PROCESSEA(EXPE,IEALL,NEALL,IYEALL,
     1       IE13,NE13,IY13,IE1,NE1,IY1,IE3,NE3,IY3,
     2       MEXEALL,MEXE13,MEXE1,MEXE3,ZS,NNODES,
     3       IEBIG13,EXPEBIG,IEBIG1,IEBIG3,ZSEAST,
     4       LOCALONOFF)
            CALL PROCESSWE(EXPW,IWALL,NWALL,IYWALL,
     1       IW24,NW24,IY24,IW2,NW2,IY2,IW4,NW4,IY4,
     2       MEXWALL,MEXW24,MEXW2,MEXW4,ZS,NNODES,
     3       IWBIG24,EXPWBIG,IWBIG2,IWBIG4,ZSWEST,
     4       LOCALONOFF)


           ELSEIF (ICHILDBOX(1,J) .LT. 0)THEN
C           NOW LET'S SCAN THE COLLEAGUES
 
            DO 250 NB = 1, 9
              IOUT = ICOLLEAGBOX(NB,J)
              IF(IOUT .LT. 0)GOTO 250
              IF(ICHILDBOX(1,IOUT) .LT. 0)THEN
C               The colleague is childless, so just do the
C               local interaction as in the uniform case
C               (just outgoing).

                CALL COLLOC8(POT(1,IOUT),COEFFS(0,0,J),NDEG,
     1                       NB,T,TLOG,W,LOCALONOFF(IOUT))


              ELSEIF(ICHILDBOX(1,IOUT) .GT. 0)THEN
C               Colleague has children (have to go to big to small and
C               small to big stuff)

                IC1 = ICHILDBOX(2,IOUT)
                IC2 = ICHILDBOX(3,IOUT)
                IC3 = ICHILDBOX(1,IOUT)
                IC4 = ICHILDBOX(4,IOUT)

C              Form the four expansions needed in the
C              big to small and small to big process.


                CALL MKEXPBTOS8(COEFFS(0,0,J), EXPNALL,
     1           EXPSALL,EXPEALL,EXPWALL, MAPSOUTH,MAPNORTH,
     2           MAPEAST,MAPWEST, SUM, XLENGTH, NNODES)


                IF(NB .EQ. 1)THEN
C                 Colleague with small boxes is in the lower left corner,
C                 one box is not well separated and 3 are.
 
                  IFAR1 = IC4
                  IFAR2 = IC2
                  IFAR3 = IC3
                  ICLOSE1 = IC1

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                       NDEG,1,T,TLOGS,WSTOB,LOCALONOFF(J))
 
 
C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                     NDEG,1,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))


C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (-2.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1), 
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))


                  FTARGET2 = (-1.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))


 
                  SPIN = 1.0D0
                  FTARGET3 = (-2.0D0,-1.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR3))

 

                ELSEIF(NB .EQ. 2)THEN
C                 Colleague with small boxes is below this box
C                 two boxes are not well separated and two are.

                  IFAR1 = IC4
                  IFAR2 = IC2
                  ICLOSE1 = IC3
                  ICLOSE2 = IC1

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                          NDEG,2,T,TLOGS,WSTOB,LOCALONOFF(J))
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE2),
     1                          NDEG,3,T,TLOGS,WSTOB,LOCALONOFF(J))
 

C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                        NDEG,2,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
                  CALL BTOSLOC8(POT(1,ICLOSE2),COEFFS(0,0,J),
     1                        NDEG,3,T,TLOG,WBTOS,LOCALONOFF(ICLOSE2))
 

C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (0.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))


                  FTARGET2 = (1.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))

 
                ELSEIF(NB .EQ. 3)THEN
C                 Colleague with small boxes is in the lower right corner,
C                 one box is not well separated and 3 are.

                  ICLOSE1 = IC3
                  IFAR1 = IC4  
                  IFAR2 = IC2 
                  IFAR3 = IC1

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                          NDEG,4,T,TLOGS,WSTOB,LOCALONOFF(J))


C                 First do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                      NDEG,4,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))


C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (2.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1), 
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2), 
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))

 
                  SPIN = -1.0D0
                  FTARGET3 = (3.0D0,-1.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR3))

            
                ELSEIF(NB .EQ. 4)THEN
C                 Colleague with small boxes is left of this box
C                 two boxes are not well separated and two are.
 
 
                  ICLOSE1 = IC2
                  ICLOSE2 = IC1
                  IFAR1 = IC4 
                  IFAR2 = IC3
 
C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                          NDEG,5,T,TLOGS,WSTOB,LOCALONOFF(J))
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE2),
     1                          NDEG,7,T,TLOGS,WSTOB,LOCALONOFF(J))
       

C                 Next do the local work, big to small 
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                          NDEG,5,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
                  CALL BTOSLOC8(POT(1,ICLOSE2),COEFFS(0,0,J),
     1                          NDEG,7,T,TLOG,WBTOS,LOCALONOFF(ICLOSE2))
 

C                 Finally do the far work, big to small
                  SPIN = 1.0D0
                  FTARGET1 = (-2.0D0,0.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (-2.0D0,1.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 6)THEN
C                 Colleague with small boxes is right this box
C                 two boxes are not well separated and two are.
            
                  IFAR1 = IC1
                  IFAR2 = IC2
                  ICLOSE1 = IC3
                  ICLOSE2 = IC4

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                         NDEG,8,T,TLOGS,WSTOB,LOCALONOFF(J))
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE2),
     1                         NDEG,6,T,TLOGS,WSTOB,LOCALONOFF(J))
 

C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                         NDEG,8,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
                  CALL BTOSLOC8(POT(1,ICLOSE2),COEFFS(0,0,J),
     1                         NDEG,6,T,TLOG,WBTOS,LOCALONOFF(ICLOSE2))
 

C                 Finally do the far work, big to small
                  SPIN = -1.0D0
                  FTARGET1 = (3.0D0,1.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,0.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 7)THEN
C                 Colleague with small boxes is in the upper left corner,
C                 one box is not well separated and 3 are.
 
                  ICLOSE1 = IC2
                  IFAR1 = IC4
                  IFAR2 = IC3
                  IFAR3 = IC1

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                         NDEG,9,T,TLOGS,WSTOB,LOCALONOFF(J))


C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                         NDEG,9,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
 

C                 Finally do the far work, big to small
                  SPIN = 1.0D0
                  FTARGET1 = (-2.0D0,2.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR1))

 
                  SPIN = IMAG
                  FTARGET2 = (-2.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))

 
                  FTARGET3 = (-1.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR3))


                ELSEIF(NB .EQ. 8)THEN
C                 Colleague with small boxes is above this box
C                 two boxes are not well separated and two are.

                  ICLOSE1 = IC4
                  ICLOSE2 = IC2
                  IFAR1 = IC3
                  IFAR2 = IC1

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                         NDEG,10,T,TLOGS,WSTOB,LOCALONOFF(J))
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE2),
     1                         NDEG,11,T,TLOGS,WSTOB,LOCALONOFF(J))


C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J),
     1                         NDEG,10,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
                  CALL BTOSLOC8(POT(1,ICLOSE2),COEFFS(0,0,J),
     1                         NDEG,11,T,TLOG,WBTOS,LOCALONOFF(ICLOSE2))


C                 Finally do the far work, big to small
                  SPIN = IMAG
                  FTARGET1 = (0.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (1.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 9)THEN
C                 Colleague with small boxes is in the upper right corner,
C                 one box is not well separated and 3 are.

                  ICLOSE1 = IC4
                  IFAR1 = IC3  
                  IFAR2 = IC1
                  IFAR3 = IC2

C                 First do the local work, small to big
                  CALL STOBLOC8(POT(1,J),COEFFS(0,0,ICLOSE1),
     1                        NDEG,12,T,TLOGS,WSTOB,LOCALONOFF(J))


C                 Next do the local work, big to small
                  CALL BTOSLOC8(POT(1,ICLOSE1),COEFFS(0,0,J), 
     1                      NDEG,12,T,TLOG,WBTOS,LOCALONOFF(ICLOSE1))
 

C                 Finally do the far work, big to small
                  SPIN = IMAG
                  FTARGET1 = (2.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))
 
                  SPIN = -1.0D0
                  FTARGET3 = (3.0D0,2.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR3))
                ENDIF
              ENDIF
250         CONTINUE
           ENDIF
         END DO

C        For boxes with children, convert the exponential expansions
C        to one local expansion:
         DO JJ = ISTART,IEND
          J = IBOXLEV(JJ)
          IF(ICHILDBOX(1,J) .GT. 0 .AND. LOCALONOFF(J) .EQ. 1)THEN

              DO II = 1, 4
                IC1 = ICHILDBOX(II,J)
                  CALL EXP4LOCAL(B,NTERMS,NNODES,EXPW(0,IC1),
     1               EXPS(0,IC1), EXPE(0,IC1), EXPN(0,IC1),COMP)
                  CALL ADDEXP(B,LOCEXP(0,IC1),NTERMS)
              END DO

          ELSEIF((ICHILDBOX(1,J) .LT. 0) .AND. LOCALONOFF(J) .EQ. 1)THEN
C          Evaluate the exponential expansions at the target points

               CALL EXPBIG4EVAL(POT(1,J),NNODES,EXPWBIG(0,J),
     1           EXPSBIG(0,J), EXPEBIG(0,J), EXPNBIG(0,J),
     2           IFLAGEAST(J),IFLAGWEST(J),
     3           IFLAGNORTH(J),IFLAGSOUTH(J),TEMPEAST,
     4           TEMPNORTH, TEMPWEST, TEMPSOUTH)
          ENDIF
         END DO
      XLENGTH = XLENGTH / 2.0D0
      ENDDO

      DO JJ = 0, NLEV
       ISTART = ISTARTLEV(JJ)
       IEND = ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
        DO II = ISTART, IEND
          I = IBOXLEV(II)
          IF(ICHILDBOX(1,I) .LT. 0 .AND. LOCALONOFF(I) .EQ. 1)THEN
           DO IP = 1,64
             ZSHIFT = DCMPLX(XP(IP),YP(IP))
             ZPOT = ZERO
              DO J = 1,NTERMS
                ZPOT = ZSHIFT*ZPOT + LOCEXP(NTERMS-J,I)
              END DO
             POT(IP,I) = POT(IP,I) + DREAL(ZPOT)
           END DO
          ENDIF
        END DO
      END DO
      RETURN
      END


C**********************************************************************
C     This subroutine is designed to merge the multipole expansions
C     of four child boxes together to form the multipole expansion
C     of the parent.  
C     The algorithm works by passing through the loop mod 4 and
C     rescaling the new coefficients so that they are on the parent's
C     level.
C
C     INPUT:
C 
C     MPOLE1 is the multipole coefficients for the child box in the
C            upper left corner
C 
C     MPOLE2 is the multipole coefficients for the child box in the
C            upper right corner
C 
C     MPOLE3 is the multipole coefficients for the child box in the
C            lower right corner
C 
C     MPOLE4 is the multipole coefficients for the child box in the
C            lower left corner
C 
C     NTERMS is the order of the multipole expansions
C 
C     C is the array that the binomial coefficients are stored in
C
C     OUTPUT:
C 
C     MPOLEPAR is the multipole coefficients of the parent box
C
C**********************************************************************
      SUBROUTINE CHILDPAR(MPOLEPAR, MPOLE1, MPOLE2, MPOLE3,
     1                    MPOLE4, NTERMS, C)
      IMPLICIT NONE
C-----Global variables
      INTEGER NTERMS
      REAL *8 C(2*NTERMS, 1)
      COMPLEX *16 MPOLE1(0:NTERMS), MPOLE2(0:NTERMS)
      COMPLEX *16 MPOLE3(0:NTERMS), MPOLE4(0:NTERMS)
      COMPLEX *16 MPOLEPAR(0:NTERMS)
C-----Local variables
      INTEGER  I, M, K, M1, M2, M3, K1, K2, K3
      REAL *8  XNTEMP, XNTEMP1
      COMPLEX *16  ATEMP(0:60), ATEMP2(0:60), ATEMP3(0:60), ATEMP4(0:60)
      COMPLEX *16 B1(60)
      COMPLEX *16 Z00, Z0P(0:60), Z0P2(0:60)
      COMPLEX *16 CD, CDD, Z0, IMAG
      COMPLEX *16 TEMP1, TEMP2, TEMP3, TEMP4
      DATA IMAG/(0.0D0,1.0D0)/

      Z0 = -.5D0 + .5D0*IMAG
      Z00=Z0
      CD=1.0D0 / Z0
      CDD=CD
      Z0P(0) =1.0d0
      Z0P2(0)=1.0d0
      DO I=1,NTERMS
         Z0P(I)  = Z00
         Z0P2(I) = CDD
         CDD     = CDD*CD
         Z00     = Z00*Z0
      END DO

      DO I=1,NTERMS
         B1(I)    = 0.0D0
      END DO

      MPOLEPAR(0) = MPOLEPAR(0) + MPOLE1(0)
     1      + MPOLE2(0) + MPOLE3(0) + MPOLE4(0)

      DO M = 1, NTERMS
       TEMP1 = MPOLE1(M) + MPOLE3(M)
       TEMP2 = MPOLE1(M) - MPOLE3(M)
       TEMP3 = MPOLE2(M) + MPOLE4(M)
       TEMP4 = IMAG*(MPOLE2(M) - MPOLE4(M))

       ATEMP(M) =  (TEMP1 + TEMP3)*Z0P2(M)
       ATEMP2(M) = (TEMP1 - TEMP3)*Z0P2(M)
       ATEMP3(M) = (TEMP2 - TEMP4)*Z0P2(M)
       ATEMP4(M) = (TEMP2 + TEMP4)*Z0P2(M)
      END DO

      XNTEMP = 2.0D0
      DO M=1,NTERMS-3,4
      M1 = M + 1
      M2 = M + 2
      M3 = M + 3
      DO K=1,NTERMS-3,4
        K1 = K + 1
        K2 = K + 2
        K3 = K + 3

        B1(M)=B1(M) + ATEMP(K)*C(M,K)
     1   + ATEMP4(K1)*C(M,K1)
     2   + ATEMP2(K2)*C(M,K2)
     3   + ATEMP3(K3)*C(M,K3)

        B1(M1)=B1(M1) + ATEMP3(K)*C(M1,K)
     1   + ATEMP(K1)*C(M1,K1)
     2   + ATEMP4(K2)*C(M1,K2)
     3   + ATEMP2(K3)*C(M1,K3)

        B1(M2)=B1(M2) + ATEMP2(K)*C(M2,K)
     1   + ATEMP3(K1)*C(M2,K1)
     2   + ATEMP(K2)*C(M2,K2)
     3   + ATEMP4(K3)*C(M2,K3)

        B1(M3)=B1(M3) + ATEMP4(K)*C(M3,K)
     1   + ATEMP2(K1)*C(M3,K1)
     2   + ATEMP3(K2)*C(M3,K2)
     3   + ATEMP(K3)*C(M3,K3)
      END DO

       TEMP1 = MPOLE1(0) - MPOLE3(0)
       TEMP2 = IMAG*(MPOLE2(0) - MPOLE4(0))
       TEMP3 = MPOLE1(0) + MPOLE3(0)
       TEMP4 = MPOLE2(0) + MPOLE4(0)

       XNTEMP1 = XNTEMP

       MPOLEPAR(M) = MPOLEPAR(M) + Z0P(M)*(B1(M)
     1     +(TEMP2 - TEMP1)/DBLE(M))/ XNTEMP1

       XNTEMP1 = 2.0D0*XNTEMP1
       MPOLEPAR(M1) = MPOLEPAR(M1) + Z0P(M1)*(B1(M1)
     1    + (TEMP4 - TEMP3)/DBLE(M1))/ XNTEMP1

       XNTEMP1 = 2.0D0*XNTEMP1
       MPOLEPAR(M2) = MPOLEPAR(M2) + Z0P(M2)*(B1(M2)
     1 -(TEMP1 + TEMP2)/DBLE(M2))/ XNTEMP1

       XNTEMP1 = 2.0D0*XNTEMP1
       MPOLEPAR(M3) = MPOLEPAR(M3) + Z0P(M3)*(B1(M3)
     1  -(TEMP3 + TEMP4)/DBLE(M3))/ XNTEMP1
    
      XNTEMP = 16.0D0*XNTEMP

      END DO
      RETURN
      END

C**********************************************************************
C     This subroutine is designed to shift the local expansions of
C     a parent box to its four children.  This is used in the downward
C     pass. 
C
C     INPUT:
C 
C     BETAHATPAR is the multipole coefficients of the parent box
C
C     NTERMS is the order of the multipole expansions
C
C     C is the array that the binomial coefficients are stored in
C
C     OUTPUT:
C
C     BETAHAT1 is the multipole coefficients for the child box in the
C            upper left corner
C
C     BETAHAT2 is the multipole coefficients for the child box in the
C            upper right corner
C
C     BETAHAT3 is the multipole coefficients for the child box in the
C            lower right corner
C
C     BETAHAT4 is the multipole coefficients for the child box in the
C            lower left corner
C
C**********************************************************************
      SUBROUTINE PARENTCHILD(BETAHATPAR, 
     1      BETA1HAT, BETA2HAT, BETA3HAT,
     2      BETA4HAT, NTERMS, C, ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NTERMS, ISWITCH
      REAL *8 C(2*NTERMS, 1)
      COMPLEX *16 BETA1HAT(0:NTERMS), BETA2HAT(0:NTERMS)
      COMPLEX *16 BETA3HAT(0:NTERMS), BETA4HAT(0:NTERMS)
      COMPLEX *16 BETAHATPAR(0:NTERMS)
C-----Local variables
      INTEGER I, M, K
      INTEGER K1,K2,K3,K4
      INTEGER M1,M2,M3
      COMPLEX *16 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5
      COMPLEX *16 TEMP6, TEMP7, TEMP8, TEMP9, TEMP10
      COMPLEX *16 Z0,Z00,Z0P(0:60),Z0P2(0:60), ANEW(0:60)
      COMPLEX *16 CD,CDD,IMAG

      DATA IMAG/(0.0D0, 1.0D0)/

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF

C     Now let's initialize an array containing
C     powers of the shift vector (-.25 + i*.25)
C     Set Z0P(i) equal to Z0^i and
C     Z0P2(i) equal to 1 / Z0^i:
      Z0 = -.25D0 + IMAG*.25D0
      Z00= Z0
      Z0P(0)=1.0D0
      CD = 1.0D0 /(Z0*2.0D0)
      CDD = CD
      Z0P2(0) = 1.0D0

      DO  I=1,NTERMS
         Z0P(I)=Z00
         Z00= Z00*Z0
         Z0P2(I) = CDD
         CDD = CDD*CD
      END DO

C     Initialize all local expansions to 
C     zero and set up the array ANEW:
      DO I=0, NTERMS
          BETA1HAT(I) = 0.0D0
          BETA2HAT(I) = 0.0D0
          BETA3HAT(I) = 0.0D0
          BETA4HAT(I) = 0.0D0
          ANEW(I) = BETAHATPAR(I)*Z0P(I)
      END DO

C     Now go through a mod 4 loop in which the
C     expansions are shifted:
      DO K=0, NTERMS-3, 4
      K1 = K + 1
      K2 = K + 2
      K3 = K + 3
      K4 = K + 4
       DO M=0, K
       M1 = M + 1
        TEMP5 = ANEW(K)*C(K1,M1)
        TEMP6 = ANEW(K1)*C(K2,M1)
        TEMP7 = ANEW(K2)*C(K3,M1)
        TEMP8 = ANEW(K3)*C(K4,M1)

        TEMP1 = TEMP5 + TEMP7
        TEMP2 = TEMP5 - TEMP7
        TEMP3 = TEMP6 + TEMP8
        TEMP4 = IMAG*(TEMP6 - TEMP8)

        BETA1HAT(M) = BETA1HAT(M) + TEMP1 + TEMP3
        BETA2HAT(M) = BETA2HAT(M) + TEMP2 - TEMP4
        BETA3HAT(M) = BETA3HAT(M) + TEMP1 - TEMP3
        BETA4HAT(M) = BETA4HAT(M) + TEMP2 + TEMP4
       END DO


        TEMP7 = ANEW(K2)*C(K3,K2)
        TEMP8 = ANEW(K3)*C(K4,K2)
        TEMP9 = IMAG*(ANEW(K1) - TEMP8)
        TEMP10 = ANEW(K1) + TEMP8
        BETA1HAT(K1) = BETA1HAT(K1) + TEMP7 + TEMP10
        BETA2HAT(K1) = BETA2HAT(K1) - TEMP9 - TEMP7
        BETA3HAT(K1) = BETA3HAT(K1) - TEMP10 + TEMP7
        BETA4HAT(K1) = BETA4HAT(K1) + TEMP9  - TEMP7


        TEMP8 = ANEW(K3)*C(K4,K3)
        TEMP9 = IMAG*TEMP8
        BETA1HAT(K2) = BETA1HAT(K2) + ANEW(K2) + TEMP8
        BETA2HAT(K2) = BETA2HAT(K2) - ANEW(K2) + TEMP9
        BETA3HAT(K2) = BETA3HAT(K2) + ANEW(K2) - TEMP8
        BETA4HAT(K2) = BETA4HAT(K2) - ANEW(K2) - TEMP9


        TEMP8 = IMAG*ANEW(K3)
        BETA1HAT(K3) = BETA1HAT(K3) + ANEW(K3)
        BETA2HAT(K3) = BETA2HAT(K3) + TEMP8
        BETA3HAT(K3) = BETA3HAT(K3) - ANEW(K3)
        BETA4HAT(K3) = BETA4HAT(K3) - TEMP8
      END DO

      DO M = 1, NTERMS-3, 4
      M1 = M + 1
      M2 = M + 2
      M3 = M + 3
         TEMP1 = Z0P2(M)*IMAG
         TEMP2 = Z0P2(M2)*IMAG

         BETA1HAT(M)  = BETA1HAT(M)*Z0P2(M)
         BETA1HAT(M1) = BETA1HAT(M1)*Z0P2(M1)
         BETA1HAT(M2) = BETA1HAT(M2)*Z0P2(M2)
         BETA1HAT(M3) = BETA1HAT(M3)*Z0P2(M3)

         BETA2HAT(M)  =  BETA2HAT(M)*TEMP1
         BETA2HAT(M1) = -BETA2HAT(M1)*Z0P2(M1)
         BETA2HAT(M2) = -BETA2HAT(M2)*TEMP2
         BETA2HAT(M3) =  BETA2HAT(M3)*Z0P2(M3)

         BETA3HAT(M)  = -BETA3HAT(M)*Z0P2(M)
         BETA3HAT(M1) =  BETA3HAT(M1)*Z0P2(M1)
         BETA3HAT(M2) = -BETA3HAT(M2)*Z0P2(M2)
         BETA3HAT(M3) =  BETA3HAT(M3)*Z0P2(M3)

         BETA4HAT(M)  = -BETA4HAT(M)*TEMP1
         BETA4HAT(M1) = -BETA4HAT(M1)*Z0P2(M1)
         BETA4HAT(M2) =  BETA4HAT(M2)*TEMP2
         BETA4HAT(M3) =  BETA4HAT(M3)*Z0P2(M3)
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine computes several things and arrays that
C     are needed for later use in the code.  
C
C     INPUT:
C 
C     Everything is blank on input
C
C     OUTPUT:
C 
C     COMP is an array that is a combination of factorial and nodal terms
C 
C     NTERMS is the order of the multipole expansions
C 
C     NNODES is the order of the plane wave expansions
C 
C     XNODES is a blank array that is set to
C            the nodes in the plane wave expansions
C 
C     WNODES is a blank array that is set to
C            the weights in the plane wave expansions
C 
C     TEMP  is a blank array that is set to
C            the ratio of the weights and nodes
C 
C     CHOOSE is the array that the binomial coefficients are stored in
C 
C     SUM is a precomputed terms needed in the exponential expansions
C
C**************************************************************************
      SUBROUTINE PRECOMPUTE(COMP, NTERMS, 
     1         NNODES, XNODES, WNODES, TEMP,
     2         CHOOSE, SUM)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES, NTERMS
      REAL *8  CHOOSE(2*NTERMS, 2*NTERMS)
      REAL *8  WNODES(1), XNODES(1), TEMP(1)
      REAL *8  COMP(NNODES,0:NTERMS), SUM
C-----Local variables
      INTEGER  I, J
      REAL *8  FACT(0:100), X

C     Set up the factorial array:
      FACT(0) = 1.0D0
      DO  I = 1, 2*NTERMS
       FACT(I) = DBLE(I)
        DO  J = I-1, 2, -1
          FACT(I) = FACT(I) * DBLE(J)
       END DO
      END DO
C     FACT(i) is now equal to i!

C     Now that we have set up the factorial
C     array, let's set up the choose array.
      DO I = 1, 2*NTERMS
       DO J = 1, I
         CHOOSE(I,J) = FACT(I-1)/(FACT(I-J)*FACT(J-1))
       END DO
      END DO
C     CHOOSE(i,j) is now equal to (i-1) choose (j-1).

C     Call a subroutine that returns precomputed nodes that will
C     accurately approximate the integral representing 1/z, note
C     that these nodes will only give accurate results for certain
C     specific values of z (see comments in the pwts4 subroutine):

      DO I = 1, NNODES
       TEMP(I) = WNODES(I) / XNODES(I)
      END DO


C     Now let's precompute the array of terms COMP(i,j).
      DO I = 1, NNODES
        X = 1.0D0
        DO J = 0, NTERMS
          COMP(I,J) = X / FACT(J)
          X = X * XNODES(I)
        END DO
      END DO
C     COMP(i,j) is now set equal to xnodes(i)^j / j!
C     This is used later on in the parent to child and
C     child to parent subroutines.

      SUM = 0.0D0
      DO I = 1, NNODES
       SUM = SUM + DEXP(-XNODES(I))*WNODES(I)/XNODES(I)
      END DO
      RETURN
      END


C*************************************************************************
C     The following subroutine takes four multipole expansions (from four
C     different directions) as input after they have been shifted to a 
C     common target point.  The four expansions are joined into one 
C     common local expansion.  
C     None of the input parameters are altered.
C     The expansions are merged by counting through the loop mod 4.
C
C     INPUT:
C
C     P is the order of the multipole expansions
C
C     NNODES is the order of the plane wave expansions
C
C     BETAW are the incoming west exponential coefficients
C
C     BETAS are the incoming south exponential coefficients
C
C     BETAE are the incoming east exponential coefficients
C
C     BETAN are the incoming north exponential coefficients
C
C     COMP is a precomputed term involving the ratio of
C          the weights and factorial terms
C
C     OUTPUT:
C
C     BETAHAT is the local coefficients of the parent box
C
C*************************************************************************
      SUBROUTINE EXP4LOCAL(BETAHAT, P, NNODES, BETAW,
     1                     BETAS, BETAE, BETAN, COMP)
      IMPLICIT NONE
C-----Global variables
      INTEGER NNODES, P
      REAL *8 COMP(NNODES,0:P)
      COMPLEX *16 BETAHAT(0:P)
      COMPLEX *16 BETAN(0:NNODES), BETAS(0:NNODES)
      COMPLEX *16 BETAE(0:NNODES), BETAW(0:NNODES)
C-----Local variables
      INTEGER I, J
      COMPLEX *16 BSUM1, BSUM2, BSUM3, BSUM4
      COMPLEX *16 BS14, BS3M2, BS1M4, BS32
      COMPLEX *16 IMAG

      DATA IMAG/(0.0D0, 1.0D0)/

C     The local coefficients are computed below.
C     They are obtained just by taking the exponential
C     expansions and performing a Taylor expansion of
C     each term.  These Taylor expansions are then combined
C     to generate one local expansion.
C
C     Upon input, the exponential expansions are in the 
C     form sum(i 0 to nnodes)beta(i)*exp(spin*xnodes(i)*z).
C     Upon output, the local expansion is of the form
C     sum(i 0 to p) betahat(i)*z^i.

      BETAHAT(0) = 0.0D0
      DO I = 1, NNODES
        BETAHAT(0) = BETAHAT(0) - (BETAW(I) + BETAE(I) +
     1          BETAS(I) + BETAN(I))
      END DO

      DO I = 1, NNODES
         BSUM1 = IMAG*(BETAN(I) - BETAS(I))
         BSUM2 = BETAE(I) + BETAW(I)
         BSUM3 = BETAN(I) + BETAS(I)
         BSUM4 = BETAW(I) - BETAE(I)
	 BS14  = BSUM4 + BSUM1
	 BS3M2 = BSUM3 - BSUM2
	 BS1M4 = BSUM1 - BSUM4
	 BS32  = BSUM3 + BSUM2

C        Initialize all terms to zero the
C        first time through the loop.
         IF(I .EQ. 1)THEN
            DO J = 1, NNODES
              BETAHAT(J) = 0.0D0
            END DO
         ENDIF

         DO J = 1, NNODES-3, 4
           BETAHAT(J)   = BETAHAT(J)   - COMP(I,J)*BS14
           BETAHAT(J+1) = BETAHAT(J+1) + COMP(I,J+1)*BS3M2
           BETAHAT(J+2) = BETAHAT(J+2) + COMP(I,J+2)*BS1M4
           BETAHAT(J+3) = BETAHAT(J+3) - COMP(I,J+3)*BS32
         END DO
      END DO
      BETAHAT(0) = BETAHAT(0) + BETAW(0) + BETAE(0)
     1           + BETAS(0) + BETAN(0)
      RETURN
      END


C*************************************************************************
C
C     The following subroutine is used to evaluate the 'Big'
C     expansions.  These are the expansions that result from the
C     Stobfar interactions.
C
C     IFLAGEAST, IFLAGWEST, IFLAGNORTH, and IFLAGSOUTH are flags that
C     indicate whether or not a given expansion needs to be evaluated.
C
C     TEMPNORTH, TEMPSOUTH, TEMPWEST, and TEMPEAST are precomputed arrays 
C     that may be needed on the evaluation.
C
C     EXPW, EXPS, EXPE, and EXPN are the exponential expansions that
C     represent the incoming Stobfar interactions.  POT is the potential
C     of the box in question and NNODES is the number of terms in the 
C     exponential expansion.
C
C*************************************************************************
      SUBROUTINE EXPBIG4EVAL(POT,NNODES,
     1      EXPW, EXPS, EXPE, EXPN,
     2      IFLAGEAST,IFLAGWEST,
     3      IFLAGNORTH,IFLAGSOUTH, TEMPEAST, 
     4      TEMPNORTH, TEMPWEST, TEMPSOUTH)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES
      INTEGER  IFLAGEAST, IFLAGNORTH
      INTEGER  IFLAGWEST, IFLAGSOUTH
      REAL *8  POT(64)
      COMPLEX *16  EXPN(0:NNODES), EXPS(0:NNODES)
      COMPLEX *16  EXPE(0:NNODES), EXPW(0:NNODES)
      COMPLEX *16  TEMPEAST(64,1), TEMPNORTH(64,1)
      COMPLEX *16  TEMPWEST(64,1), TEMPSOUTH(64,1)
C-----Local variables
      INTEGER  I, J

        IF(IFLAGEAST  .EQ. 0 .AND.  IFLAGWEST .EQ. 0 .AND.
     1     IFLAGNORTH .EQ. 0 .AND. IFLAGSOUTH .EQ. 0)THEN
              RETURN
        ENDIF


        IF(IFLAGWEST .EQ. 1)THEN
         DO J = 1, 64
          POT(J) = POT(J) + DREAL(EXPW(0))
          DO I = 1, NNODES
           POT(J) = POT(J) + DREAL(EXPW(I)*TEMPWEST(J,I))
          END DO
         END DO
        ENDIF

        IF(IFLAGEAST .EQ. 1)THEN
         DO J = 1, 64
          POT(J) = POT(J) + DREAL(EXPE(0))
          DO I = 1, NNODES
           POT(J) = POT(J) + DREAL(EXPE(I)*TEMPEAST(J,I))
          END DO
         END DO
        ENDIF

        IF(IFLAGNORTH .EQ. 1)THEN
         DO J = 1, 64
          POT(J) = POT(J) + DREAL(EXPN(0))
          DO I = 1, NNODES
           POT(J) = POT(J) + DREAL(EXPN(I)*TEMPNORTH(J,I))
          END DO
         END DO
        ENDIF

        IF(IFLAGSOUTH .EQ. 1)THEN
         DO J = 1, 64
          POT(J) = POT(J) + DREAL(EXPS(0))
          DO I = 1, NNODES
           POT(J) = POT(J) + DREAL(EXPS(I)*TEMPSOUTH(J,I))
          END DO
         END DO
        ENDIF

      RETURN
      END


C***********************************************************************
C     This file contains all of the expansion creation routines for 
C     a parent box from its four children.
C
C     MKEXPS creates the table of shifting coefficients in the physical
C     domain e^{lambda z} for a given lambda discretization.
C
C     MKNSEXP creates all NORTH and SOUTH expansions centered at child 1
C
C     MKEWEXP creates all EAST and WEST expansions centered at child 1
C
C***********************************************************************
      SUBROUTINE MKSHIFTS2D(XNODES,NNODES,ZS,
     1    ZSEAST,ZSWEST,ZSNORTH,ZSSOUTH)
      IMPLICIT NONE
      COMPLEX *16 ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16 ZSEAST(-3:3,-3:3,0:NNODES)
      COMPLEX *16 ZSWEST(-3:3,-3:3,0:NNODES)
      COMPLEX *16 ZSNORTH(-3:3,-3:3,0:NNODES)
      COMPLEX *16 ZSSOUTH(-3:3,-3:3,0:NNODES)
      COMPLEX *16 TEMPEXP
      REAL *8     XNODES(NNODES)
      INTEGER  NNODES, K, M, N
C     Loop over each lambda value 
      DO K = 1,NNODES
         DO N = -3,3
            DO M = -3,3
               ZS(N,M,K)=CDEXP(-XNODES(K)*DCMPLX(DBLE(N),DBLE(M)))
            ENDDO
         ENDDO
      ENDDO
      DO N = -3,3
         DO M = -3,3
            ZS(N,M,0)=1.0D0
         ENDDO
      ENDDO

      DO K = 1,NNODES
         TEMPEXP = CDEXP(-XNODES(K)*DCMPLX(0.5D0,0.5D0))
         DO N = -3,3
            DO M = -3,3
               ZSEAST(N,M,K)=ZS(N,M,K) * TEMPEXP
            ENDDO
         ENDDO
      ENDDO
      DO N = -3,3
         DO M = -3,3
            ZSEAST(N,M,0)=1.0D0
         ENDDO
      ENDDO

      DO K = 1,NNODES
         TEMPEXP = CDEXP(-XNODES(K)*DCMPLX(-0.5D0,-0.5D0))
         DO N = -3,3
            DO M = -3,3
               ZSWEST(N,M,K)=ZS(N,M,K) * TEMPEXP
            ENDDO
         ENDDO
      ENDDO
      DO N = -3,3
         DO M = -3,3
            ZSWEST(N,M,0)=1.0D0
         ENDDO
      ENDDO

      DO K = 1,NNODES
         TEMPEXP = CDEXP(-XNODES(K)*DCMPLX(0.5D0,-0.5D0))
         DO N = -3,3
            DO M = -3,3
               ZSNORTH(N,M,K)=ZS(N,M,K) * TEMPEXP
            ENDDO
         ENDDO
      ENDDO
      DO N = -3,3
         DO M = -3,3
            ZSNORTH(N,M,0)=1.0D0
         ENDDO
      ENDDO

      DO K = 1,NNODES
         TEMPEXP = CDEXP(-XNODES(K)*DCMPLX(-0.5D0,0.5D0))
         DO N = -3,3
            DO M = -3,3
               ZSSOUTH(N,M,K)=ZS(N,M,K) * TEMPEXP
            ENDDO
         ENDDO
      ENDDO
      DO N = -3,3
         DO M = -3,3
            ZSSOUTH(N,M,0)=1.0D0
         ENDDO
      ENDDO

      RETURN
      END

C***********************************************************************
      SUBROUTINE MKEXP2D(IBOX,NTERMS,MPOLE,
     1           NNODES,MEXNALL,MEXN12,MEXSALL,MEXS34,
     2           MEXEALL,MEXE13,MEXE1,MEXE3,MEXWALL,MEXW24,MEXW2,MEXW4,
     3           ZS,COMP,WNODES,SCALE, SUM, TEMP,
     4           ICHILDBOX)
C***********************************************************************
      IMPLICIT NONE
      INTEGER  IBOX
      INTEGER  NNODES,IC
      INTEGER  ICHILD(4)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NTERMS, JJ
      REAL *8     XLENGTH
      COMPLEX *16  MPOLE(0:NTERMS,*)
      COMPLEX *16  ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16  EXPN(0:50),EXPS(0:50),EXPE(0:50),EXPW(0:50)
      COMPLEX *16  MEXNALL(0:NNODES),MEXN12(0:NNODES)
      COMPLEX *16  MEXSALL(0:NNODES),MEXS34(0:NNODES)
      COMPLEX *16  MEXEALL(0:NNODES),MEXE13(0:NNODES)
      COMPLEX *16  MEXE1(0:NNODES),MEXE3(0:NNODES)
      COMPLEX *16  MEXWALL(0:NNODES),MEXW24(0:NNODES)
      COMPLEX *16  MEXW2(0:NNODES),MEXW4(0:NNODES)
      REAL *8  COMP(NNODES,0:NTERMS)
      REAL *8  WNODES(NNODES),SCALE, SUM, SUM1
      REAL *8  TEMP(NNODES)
C***********************************************************************
C
C     This subroutine creates the north (+y)  and south (-y) exponential 
C     expansions for a parent box due to all four children. 
C
C     Some intelligence is used in the order of summation. 
C
C***********************************************************************

      XLENGTH = 1.0D0/SCALE
      SUM1 = SUM + DLOG(XLENGTH)

      ICHILD(1) = ICHILDBOX(1,IBOX)
      ICHILD(2) = ICHILDBOX(2,IBOX)
      ICHILD(3) = ICHILDBOX(3,IBOX)
      ICHILD(4) = ICHILDBOX(4,IBOX)

C
C     include contributions from child 1
C

      IC = ICHILD(4)

      CALL EXPCOEFF(EXPE,EXPW,EXPN,EXPS,
     1                    NNODES, MPOLE(0,IC), NTERMS, COMP,
     2                    WNODES, SUM1, TEMP)

      DO JJ = 0,NNODES
         MEXN12(JJ) = EXPN(JJ)
         MEXSALL(JJ) = EXPS(JJ)
      ENDDO
      DO JJ = 0,NNODES
         MEXE1(JJ) = EXPE(JJ)
         MEXWALL(JJ) = EXPW(JJ)
      ENDDO

C
C     include contributions from child 2
C
      IC = ICHILD(3)


      CALL EXPCOEFF(EXPE,EXPW,EXPN,EXPS,
     1                    NNODES, MPOLE(0,IC), NTERMS, COMP,
     2                    WNODES, SUM1, TEMP)


      DO JJ = 0,NNODES
         EXPN(JJ) = EXPN(JJ)*ZS(0,1,JJ)
         MEXN12(JJ) = MEXN12(JJ) + EXPN(JJ)
         EXPS(JJ) = EXPS(JJ)*ZS(0,-1,JJ)
         MEXSALL(JJ) = MEXSALL(JJ) + EXPS(JJ)
      ENDDO

      DO JJ = 0,NNODES
         MEXEALL(JJ) = EXPE(JJ)*ZS(-1,0,JJ)
         MEXW2(JJ) = EXPW(JJ)*ZS(1,0,JJ)
      ENDDO
C
C     include contributions from child 3
C
      IC = ICHILD(1)


      CALL EXPCOEFF(EXPE,EXPW,EXPN,EXPS,
     1                    NNODES, MPOLE(0,IC), NTERMS, COMP,
     2                    WNODES, SUM1, TEMP)


      DO JJ = 0,NNODES
         EXPN(JJ) = EXPN(JJ)*ZS(-1,0,JJ)
         MEXNALL(JJ) = EXPN(JJ) + MEXN12(JJ)
         EXPS(JJ) = EXPS(JJ)*ZS(1,0,JJ)
         MEXS34(JJ) = EXPS(JJ)
      ENDDO

      DO JJ = 0,NNODES
         MEXE3(JJ) = EXPE(JJ)*ZS(0,-1,JJ)
         MEXE13(JJ) = MEXE1(JJ) + MEXE3(JJ)
         MEXEALL(JJ) = MEXEALL(JJ) + MEXE13(JJ)
         EXPW(JJ) = EXPW(JJ)*ZS(0,1,JJ)
         MEXWALL(JJ) = MEXWALL(JJ) + EXPW(JJ)
      ENDDO
C
C     include contributions from child 4
C
      IC = ICHILD(2)

      CALL EXPCOEFF(EXPE,EXPW,EXPN,EXPS,
     1                    NNODES, MPOLE(0,IC), NTERMS, COMP,
     2                    WNODES, SUM1, TEMP)



      DO JJ = 0,NNODES
         EXPN(JJ) = EXPN(JJ)*ZS(-1,1,JJ)
         MEXNALL(JJ) = MEXNALL(JJ) + EXPN(JJ)
         EXPS(JJ) = EXPS(JJ)*ZS(1,-1,JJ)
         MEXS34(JJ) = MEXS34(JJ) + EXPS(JJ)
         MEXSALL(JJ) = MEXSALL(JJ) + MEXS34(JJ)
      ENDDO

      DO JJ = 0,NNODES
         EXPE(JJ) = EXPE(JJ)*ZS(-1,-1,JJ)
         MEXEALL(JJ) = MEXEALL(JJ) + EXPE(JJ)
         MEXW4(JJ) = EXPW(JJ)*ZS(1,1,JJ)
         MEXW24(JJ) = MEXW2(JJ) + MEXW4(JJ)
         MEXWALL(JJ) = MEXWALL(JJ) + MEXW24(JJ)
      ENDDO
      RETURN
      END


C***********************************************************************
C     This subroutine is set up to compute all of the north, south,
C     east, and west interaction lists.  Because this is for the
C     adaptive case, all of these lists have to be computed in one
C     loop.  North and south take precedence over the east and
C     west and this is reflected  in the else statements within
C     the routine.  Processing is done at the parent level.
C     In the periodic case, the procedure is basically 
C     the same, but if a column or row number lies outside the center
C     box, it is readjusted to account for the periodicity.
C     When the lists are actually processed, there is no difference
C     between free space and periodic case.
C    
C     INPUT:
C
C     IBOX denotes the box being considered
C
C     LEVEL is the level of IBOX
C
C     ICOLLEAGBOX, IROWBOX, AND ICOLBOX define the tree
C
C     IPERIOD denotes whether or not the case is periodic or free space
C
C
C     OUTPUT:
C
C     The naming convention for the lists is is as follows:
C
C     INALL is an array that denotes the boxes in the
C     north all list.  
C
C     NNALL is the number of boxes in the north all
C     list.  
C
C     IYNALL represents the corresponding offsets of the boxes
C     in the north all list.
C
C     The same convention is used for the south all, east all, west all,
C     north12, south34, east13, west24, west2, west4, east1,
C     and east3 lists.
C
C***********************************************************************
      SUBROUTINE MKLISTS(IBOX,INALL,NNALL,IYNALL, IN12,NN12,IY12,                        
     1    ISALL,NSALL,IYSALL,IS34,NS34,IY34,IEALL,NEALL,
     2    IYEALL,IE13,NE13,IY13,IWALL,NWALL,IYWALL,
     3    IW24,NW24,IY24,IW2,IY2,NW2,IW4,IY4,NW4,
     4    IE1,IY1,NE1,IE3,IY3,NE3,
     5    INBIG12,ISBIG34,IEBIG13,
     6    IWBIG24,IEBIG1, IWBIG2, 
     7    IEBIG3, IWBIG4, ICOLLEAGBOX,ICHILDBOX,
     8    ICOLBOX, IROWBOX, IPERIOD,IFLAGEAST, IFLAGWEST, 
     9    IFLAGNORTH,IFLAGSOUTH, LOCALONOFF)
      IMPLICIT NONE
C-----Global variables
      INTEGER  INALL(6),NNALL,IYNALL(6)
      INTEGER  ISALL(6),NSALL,IYSALL(6)
      INTEGER  IEALL(4),NEALL,IYEALL(4)
      INTEGER  IWALL(4),NWALL,IYWALL(4)
      INTEGER  IN12(2),NN12,IY12(2)
      INTEGER  IS34(2),NS34,IY34(2)
      INTEGER  IW24(2),NW24,IY24(2)
      INTEGER  IE13(2),NE13,IY13(2)
      INTEGER  IW4(1),IY4(1)
      INTEGER  NW4
      INTEGER  IW2(1),IY2(1)
      INTEGER  NW2
      INTEGER  IE1(1),IY1(1)
      INTEGER  NE1
      INTEGER  IE3(1),IY3(1)
      INTEGER  NE3
      INTEGER  LOCALONOFF(1)
      INTEGER  IBOX
      INTEGER  INBIG12(3), ISBIG34(3)
      INTEGER  IEBIG13(1), IWBIG24(1)
      INTEGER  IEBIG1(1), IWBIG2(1)
      INTEGER  IEBIG3(1), IWBIG4(1)
      INTEGER  ICOLLEAGBOX(9,1), ICHILDBOX(4,1)
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  IPERIOD
      INTEGER  IFLAGEAST(1), IFLAGNORTH(1)
      INTEGER  IFLAGWEST(1), IFLAGSOUTH(1)
C-----Local variables
      INTEGER  I, J, IOUT
      INTEGER  ICHILD, NCNTR, SCNTR, ECNTR, WCNTR
      INTEGER  N12CNTR, S34CNTR
      INTEGER  W24CNTR, E13CNTR
      INTEGER  W4CNTR, W2CNTR
      INTEGER  E1CNTR, E3CNTR
      INTEGER  ICOLTEST, IROWTEST
      INTEGER  ICOL, IROW

C     Initially, set all list entries
C     and offsets to zero.
      DO J = 1, 6
        INALL(J)  = 0
        IYNALL(J) = 0
        ISALL(J)  = 0
        IYSALL(J) = 0
      END DO

      DO J = 1, 4
        IEALL(J)  = 0
        IYEALL(J) = 0
        IWALL(J)  = 0
        IYWALL(J) = 0
      END DO

      DO J = 1, 4
        ISBIG34(J) = -1
        INBIG12(J) = -1
      END DO

      DO J = 1, 2
        IN12(J) = 0
        IY12(J) = 0
        IE13(J) = 0
        IY13(J) = 0
        IW24(J) = 0
        IY24(J) = 0
        IS34(J) = 0
        IY34(J) = 0
      END DO

      IE1(1) = 0
      IY1(1) = 0
      IW2(1) = 0
      IY2(1) = 0
      IE3(1) = 0
      IY3(1) = 0
      IW4(1) = 0
      IY4(1) = 0

      IEBIG13(1) = -1
      IWBIG24(1) = -1
      IEBIG1(1)  = -1
      IEBIG3(1)  = -1
      IWBIG2(1)  = -1
      IWBIG4(1)  = -1


C     All of the offsets are set based from the
C     box in the lower left corner (child 4 in
C     out ordering convention).
C     ICOL and IROW are the rows and columns of
C     the box whose list we are trying to generate.

      ICOL = ICOLBOX(ICHILDBOX(4,IBOX))
      IROW = IROWBOX(ICHILDBOX(4,IBOX))


C     First do the free space case:
      IF(IPERIOD .EQ. 0)THEN
C     Set all of the counters to 1
      NCNTR   =  1
      SCNTR   =  1
      ECNTR   =  1
      WCNTR   =  1
      N12CNTR =  1
      S34CNTR =  1
      W24CNTR =  1
      E13CNTR =  1
      W4CNTR  =  1
      W2CNTR  =  1
      E1CNTR  =  1
      E3CNTR  =  1
C     First scan through all nine of the boxes colleagues
      DO 100 I = 1, 9
      IOUT = ICOLLEAGBOX(I,IBOX)
C     Test to see if this colleague doesn't exist or is
C     childless, if so, skip it.
      IF(IOUT .LT. 0)GOTO 100
      IF(ICHILDBOX(1,IOUT) .GT. 0)THEN
C       Scan all four of the colleagues children.
        DO J = 1, 4
C       ICOLTEST and IROWTEST represent the row and column
C       of the box being checked.
        ICHILD = ICHILDBOX(J,IOUT)
        ICOLTEST = ICOLBOX(ICHILD)
        IROWTEST = IROWBOX(ICHILD)


        IF(IROWTEST .EQ. IROW+3)THEN
           INALL(NCNTR) = ICHILD
           IYNALL(NCNTR) = ICOLTEST - ICOL
           NCNTR = NCNTR + 1
        ELSEIF(IROWTEST .EQ. IROW-2)THEN
           ISALL(SCNTR) = ICHILD
           IYSALL(SCNTR) = ICOLTEST - ICOL
           SCNTR = SCNTR + 1
        ELSEIF(ICOLTEST .EQ. ICOL-2)THEN
           IWALL(WCNTR) = ICHILD
           IYWALL(WCNTR) = IROWTEST - IROW
           WCNTR = WCNTR + 1
        ELSEIF(ICOLTEST .EQ. ICOL+3)THEN
           IEALL(ECNTR) = ICHILD
           IYEALL(ECNTR) = IROWTEST - IROW
           ECNTR = ECNTR + 1
        ELSEIF(IROWTEST .EQ. IROW+2)THEN
           IN12(N12CNTR) = ICHILD
           IY12(N12CNTR) = ICOLTEST - ICOL
           N12CNTR = N12CNTR + 1
            IF(ICOLTEST .EQ. ICOL-1)THEN
                IW4(W4CNTR) = ICHILD
                IY4(W4CNTR) = IROWTEST - IROW
                W4CNTR = W4CNTR + 1
            ENDIF
            IF(ICOLTEST .EQ. ICOL+2)THEN
                IE3(E3CNTR) = ICHILD
                IY3(E3CNTR) = IROWTEST - IROW
                E3CNTR = E3CNTR + 1
            ENDIF
        ELSEIF(IROWTEST .EQ. IROW-1)THEN
           IS34(S34CNTR) = ICHILD
           IY34(S34CNTR) = ICOLTEST - ICOL
           S34CNTR = S34CNTR + 1
            IF(ICOLTEST .EQ. ICOL-1)THEN
                IW2(W2CNTR) = ICHILD
                IY2(W2CNTR) = IROWTEST - IROW
                W2CNTR = W2CNTR + 1
            ENDIF
            IF(ICOLTEST .EQ. ICOL+2)THEN
                IE1(E1CNTR) = ICHILD
                IY1(E1CNTR) = IROWTEST - IROW
                E1CNTR = E1CNTR + 1
            ENDIF
        ELSEIF(ICOLTEST .EQ. ICOL-1)THEN
           IW24(W24CNTR) = ICHILD
           IY24(W24CNTR) = IROWTEST - IROW
           W24CNTR = W24CNTR + 1
        ELSEIF(ICOLTEST .EQ. ICOL+2)THEN
           IE13(E13CNTR) = ICHILD
           IY13(E13CNTR) = IROWTEST - IROW
           E13CNTR = E13CNTR + 1
        ENDIF
        ENDDO
      ELSEIF(ICHILDBOX(1,IOUT) .LT. 0)THEN
        IF(I .EQ. 1)THEN
          ISBIG34(1) = IOUT
          IWBIG2(1) = IOUT

          IFLAGSOUTH(IOUT) = 1
          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 2)THEN
          ISBIG34(2) = IOUT

          IFLAGSOUTH(IOUT) = 1
        ELSEIF(I .EQ. 3)THEN
          ISBIG34(3) = IOUT
          IEBIG1(1) = IOUT

          IFLAGEAST(IOUT) = 1
          IFLAGSOUTH(IOUT) = 1
        ELSEIF(I .EQ. 4)THEN
          IWBIG24(1) = IOUT

          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 6)THEN
          IEBIG13(1) = IOUT

          IFLAGEAST(IOUT) = 1
        ELSEIF(I .EQ. 7)THEN
          INBIG12(1) = IOUT
          IWBIG4(1) = IOUT

          IFLAGNORTH(IOUT) = 1
          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 8)THEN
          INBIG12(2) = IOUT

          IFLAGNORTH(IOUT) = 1
        ELSEIF(I .EQ. 9)THEN
          INBIG12(3) = IOUT
          IEBIG3(1) = IOUT

          IFLAGNORTH(IOUT) = 1
          IFLAGEAST(IOUT) = 1
        ENDIF
      ENDIF
100   CONTINUE
      NNALL = NCNTR   -  1
      NSALL = SCNTR   -  1
      NEALL = ECNTR   -  1
      NWALL = WCNTR   -  1
      NN12  = N12CNTR -  1
      NS34  = S34CNTR -  1
      NW24  = W24CNTR -  1
      NE13  = E13CNTR -  1
      NW4   = W4CNTR  -  1
      NW2   = W2CNTR  -  1
      NE1   = E1CNTR  -  1
      NE3   = E3CNTR  -  1

C     Now let's do the periodic case:
      ELSEIF(IPERIOD .EQ. 1 .OR. IPERIOD .EQ. 2)THEN
C     Set all of the counters to 1
      NCNTR   =  1
      SCNTR   =  1
      ECNTR   =  1
      WCNTR   =  1
      N12CNTR =  1
      S34CNTR =  1
      W24CNTR =  1
      E13CNTR =  1
      W4CNTR  =  1
      W2CNTR  =  1
      E1CNTR  =  1
      E3CNTR  =  1
C     First scan through all nine of the boxes colleagues
      DO 200 I = 1, 9
      IOUT = ICOLLEAGBOX(I,IBOX)
C     Test to see if this colleague doesn't exist or is
C     childless, if so, skip it.
      IF(IOUT .LT. 0)GOTO 200
      IF(ICHILDBOX(1,IOUT) .GT. 0)THEN
C       Scan all four of the colleagues children.
        DO J = 1, 4
         IF(LOCALONOFF(ICHILDBOX(J,IOUT)) .EQ. 1)THEN
          IF(I .EQ. 7)THEN
             IF(J .EQ. 1)THEN
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = -2
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 2)THEN 
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = -1
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 4)THEN
               IWALL(WCNTR) = ICHILDBOX(J,IOUT)
               IYWALL(WCNTR) = 2
               WCNTR = WCNTR + 1
             ELSEIF(J .EQ. 3)THEN 
               IN12(N12CNTR) = ICHILDBOX(J,IOUT)
               IY12(N12CNTR) = -1
               N12CNTR = N12CNTR + 1
               IW4(W4CNTR) = ICHILDBOX(J,IOUT)
               IY4(W4CNTR) = 2
               W4CNTR = W4CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 8)THEN
             IF(J .EQ. 1)THEN
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = 0
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = 1
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 3)THEN 
               IN12(N12CNTR) = ICHILDBOX(J,IOUT)
               IY12(N12CNTR) = 1
               N12CNTR = N12CNTR + 1
             ELSEIF(J .EQ. 4)THEN 
               IN12(N12CNTR) = ICHILDBOX(J,IOUT)
               IY12(N12CNTR) = 0
               N12CNTR = N12CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 9)THEN
             IF(J .EQ. 1)THEN
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = 2
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               INALL(NCNTR)  = ICHILDBOX(J,IOUT)
               IYNALL(NCNTR) = 3
               NCNTR = NCNTR + 1
             ELSEIF(J .EQ. 3)THEN
               IEALL(ECNTR) = ICHILDBOX(J,IOUT)
               IYEALL(ECNTR) = 2
               ECNTR = ECNTR + 1
             ELSEIF(J .EQ. 4)THEN 
               IN12(N12CNTR) = ICHILDBOX(J,IOUT)
               IY12(N12CNTR) = 2
               N12CNTR = N12CNTR + 1
               IE3(E3CNTR) = ICHILDBOX(J,IOUT)
               IY3(E3CNTR) = 2
               E3CNTR = E3CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 1)THEN
             IF(J .EQ. 4)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = -2
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 3)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = -1
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 1)THEN
               IWALL(WCNTR) = ICHILDBOX(J,IOUT)
               IYWALL(WCNTR) = -1
               WCNTR = WCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               IS34(S34CNTR) = ICHILDBOX(J,IOUT)
               IY34(S34CNTR) = -1
               S34CNTR = S34CNTR + 1
               IW2(W2CNTR) = ICHILDBOX(J,IOUT)
               IY2(W2CNTR) = -1
               W2CNTR = W2CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 2)THEN
             IF(J .EQ. 4)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = 0
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 3)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = 1
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               IS34(S34CNTR) = ICHILDBOX(J,IOUT)
               IY34(S34CNTR) = 1
               S34CNTR = S34CNTR + 1
             ELSEIF(J .EQ. 1)THEN
               IS34(S34CNTR) = ICHILDBOX(J,IOUT)
               IY34(S34CNTR) = 0
               S34CNTR = S34CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 3)THEN
             IF(J .EQ. 4)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = 2
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 3)THEN
               ISALL(SCNTR) = ICHILDBOX(J,IOUT)
               IYSALL(SCNTR) = 3
               SCNTR = SCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               IEALL(ECNTR) = ICHILDBOX(J,IOUT)
               IYEALL(ECNTR) = -1
               ECNTR = ECNTR + 1
             ELSEIF(J .EQ. 1)THEN
               IS34(S34CNTR) = ICHILDBOX(J,IOUT)
               IY34(S34CNTR) = 2
               S34CNTR = S34CNTR + 1
               IE1(E1CNTR) = ICHILDBOX(J,IOUT)
               IY1(E1CNTR) = -1
               E1CNTR = E1CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 4)THEN
             IF(J .EQ. 1)THEN
               IWALL(WCNTR) = ICHILDBOX(J,IOUT)
               IYWALL(WCNTR) = 1
               WCNTR = WCNTR + 1
             ELSEIF(J .EQ. 4)THEN
               IWALL(WCNTR) = ICHILDBOX(J,IOUT)
               IYWALL(WCNTR) = 0
               WCNTR = WCNTR + 1
             ELSEIF(J .EQ. 2)THEN
               IW24(W24CNTR) = ICHILDBOX(J,IOUT)
               IY24(W24CNTR) = 1
               W24CNTR = W24CNTR + 1
             ELSEIF(J .EQ. 3)THEN
               IW24(W24CNTR) = ICHILDBOX(J,IOUT)
               IY24(W24CNTR) = 0
               W24CNTR = W24CNTR + 1
             ENDIF
          ELSEIF(I .EQ. 6)THEN
             IF(J .EQ. 3)THEN
               IEALL(ECNTR) = ICHILDBOX(J,IOUT)
               IYEALL(ECNTR) = 0
               ECNTR = ECNTR + 1
             ELSEIF(J .EQ. 2)THEN
               IEALL(ECNTR) = ICHILDBOX(J,IOUT)
               IYEALL(ECNTR) = 1
               ECNTR = ECNTR + 1
             ELSEIF(J .EQ. 4)THEN
               IE13(E13CNTR) = ICHILDBOX(J,IOUT)
               IY13(E13CNTR) = 0
               E13CNTR = E13CNTR + 1
             ELSEIF(J .EQ. 1)THEN
               IE13(E13CNTR) = ICHILDBOX(J,IOUT)
               IY13(E13CNTR) = 1
               E13CNTR = E13CNTR + 1
             ENDIF
           ENDIF
          ENDIF
        ENDDO
      ELSEIF(ICHILDBOX(1,IOUT) .LT. 0 .AND.
     1       LOCALONOFF(IOUT) .EQ. 1)THEN
        IF(I .EQ. 1)THEN
          ISBIG34(1) = IOUT
          IWBIG2(1) = IOUT

          IFLAGSOUTH(IOUT) = 1
          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 2)THEN
          ISBIG34(2) = IOUT

          IFLAGSOUTH(IOUT) = 1
        ELSEIF(I .EQ. 3)THEN
          ISBIG34(3) = IOUT
          IEBIG1(1) = IOUT

          IFLAGEAST(IOUT) = 1
          IFLAGSOUTH(IOUT) = 1
        ELSEIF(I .EQ. 4)THEN
          IWBIG24(1) = IOUT

          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 6)THEN
          IEBIG13(1) = IOUT

          IFLAGEAST(IOUT) = 1
        ELSEIF(I .EQ. 7)THEN
          INBIG12(1) = IOUT
          IWBIG4(1) = IOUT

          IFLAGNORTH(IOUT) = 1
          IFLAGWEST(IOUT) = 1
        ELSEIF(I .EQ. 8)THEN
          INBIG12(2) = IOUT

          IFLAGNORTH(IOUT) = 1
        ELSEIF(I .EQ. 9)THEN
          INBIG12(3) = IOUT
          IEBIG3(1) = IOUT

          IFLAGNORTH(IOUT) = 1
          IFLAGEAST(IOUT) = 1
        ENDIF
      ENDIF
200   CONTINUE
      NNALL = NCNTR   -  1
      NSALL = SCNTR   -  1
      NEALL = ECNTR   -  1
      NWALL = WCNTR   -  1
      NN12  = N12CNTR -  1
      NS34  = S34CNTR -  1
      NW24  = W24CNTR -  1
      NE13  = E13CNTR -  1
      NW4   = W4CNTR  -  1
      NW2   = W2CNTR  -  1
      NE1   = E1CNTR  -  1
      NE3   = E3CNTR  -  1
      ENDIF
      RETURN
      END 


C***********************************************************************
C
C     This subroutine processes the NORTH interaction lists.
C
C     INPUT:
C
C     INALL(NNALL), IYNALL(NNALL) are the boxes 
C          receiving all child box data and the X and Y offsets from
C          child 1, respectively.
C     The other north lists are similarly defined (see MKNOLIST).
C
C     MEXNALL is the exponential expansion for all eight children, etc.
C     ZS are the shift coefficients computed by subroutine
C          MKEXPS.
C
C     OUTPUT:
C
C     LEXP1, which contains the local NORTH expansion information for
C            all boxes, is incremented for each box in the 
C            interaction lists.
C
C
C***********************************************************************
      SUBROUTINE PROCESSNO(LEXP1,INALL,NNALL,IYNALL,
     1           IN12,NN12,IY12,MEXNALL,MEXN12,ZS,NNODES,
     2           INBIG12,EXPNBIG,ZSNORTH,LOCALONOFF)
      IMPLICIT NONE
      INTEGER  INALL(NNALL),NNALL,IYNALL(NNALL)
      INTEGER  IN12(NN12),NN12,IY12(NN12)
      INTEGER  I, J, LOCALONOFF(1)
      INTEGER  NNODES
      INTEGER  INBIG12(3)
      COMPLEX *16  LEXP1(0:NNODES,1)
      COMPLEX *16  EXPNBIG(0:NNODES,1)
      COMPLEX *16  MEXNALL(0:NNODES)
      COMPLEX *16  MEXN12(0:NNODES)
      COMPLEX *16  ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZSNORTH(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZMUL
C
      DO I = 1,NNALL
         IF(LOCALONOFF(INALL(I)) .EQ. 1)THEN
          DO J = 0,NNODES
            ZMUL = ZS(3,-IYNALL(I),J)
            LEXP1(J,INALL(I)) = LEXP1(J,INALL(I)) +
     1      MEXNALL(J)*ZMUL
          ENDDO
         ENDIF
      ENDDO
C
      DO I = 1,NN12
         IF(LOCALONOFF(IN12(I)) .EQ. 1)THEN
          DO J = 0,NNODES
             ZMUL = ZS(2,-IY12(I),J)
             LEXP1(J,IN12(I)) = LEXP1(J,IN12(I)) +
     1            MEXN12(J)*ZMUL
          ENDDO
         ENDIF
      ENDDO

      DO I = 1, 3
       IF(LOCALONOFF(INBIG12(I)) .EQ. 1 .AND.
     1    INBIG12(I) .GT. 0)THEN
         J = 0
         ZMUL = ZSNORTH(2,-2*(I-2),J) 
         EXPNBIG(J,INBIG12(I)) = EXPNBIG(J,INBIG12(I)) +
     1          MEXN12(J)*ZMUL

         DO J = 1, NNODES
           ZMUL = ZSNORTH(2,-2*(I-2),J)
           EXPNBIG(J,INBIG12(I)) = EXPNBIG(J,INBIG12(I))
     1                         + MEXN12(J)*ZMUL
         ENDDO
       ENDIF
      ENDDO
      RETURN
      END


C***********************************************************************
C
C     This subroutine processes the SOUTH interaction lists.
C
C     INPUT:
C
C     ISALL(NSALL), IYSALL(NSALL) are the boxes 
C          receiving all child box data and the X and Y offsets from
C          child 1, respectively.
C     The other south lists are similarly defined (see MKSOLIST).
C
C     MEXSALL is the exponential expansion for all eight children, etc.
C     ZS are the shift coefficients computed by subroutine
C          MKEXPS.
C
C     OUTPUT:
C
C     LEXP2, which contains the local NORTH expansion information for
C            all boxes, is incremented for each box in the 
C            interaction lists.
C
C
C***********************************************************************
      SUBROUTINE PROCESSSO(LEXP2,ISALL,NSALL,IYSALL,
     1           IS34,NS34,IY34,MEXSALL,MEXS34,ZS,NNODES,
     2           ISBIG34,EXPSBIG,ZSSOUTH,LOCALONOFF)
      IMPLICIT NONE
      INTEGER  ISALL(NSALL),NSALL,IYSALL(NSALL)
      INTEGER  IS34(NS34),NS34,IY34(NS34)
      INTEGER  I, J,LOCALONOFF(1)
      INTEGER  NNODES
      INTEGER  ISBIG34(3)
      COMPLEX *16  LEXP2(0:NNODES,1)
      COMPLEX *16  MEXSALL(0:NNODES)
      COMPLEX *16  MEXS34(0:NNODES)
      COMPLEX *16  ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZSSOUTH(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZMUL
      COMPLEX *16  EXPSBIG(0:NNODES,1)

C
      DO I = 1,NSALL
        IF(LOCALONOFF(ISALL(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(2,IYSALL(I),J)
            LEXP2(J,ISALL(I)) = LEXP2(J,ISALL(I)) +
     1            MEXSALL(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NS34
        IF(LOCALONOFF(IS34(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(1,IY34(I),J)
            LEXP2(J,IS34(I)) = LEXP2(J,IS34(I)) +
     1            MEXS34(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO

      DO I = 1, 3
       IF(ISBIG34(I) .GT. 0 .AND.
     1    LOCALONOFF(ISBIG34(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSSOUTH(2,2*(I-2),J)
          EXPSBIG(J,ISBIG34(I)) = EXPSBIG(J,ISBIG34(I)) +
     1          MEXS34(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSSOUTH(2,2*(I-2),J)
          EXPSBIG(J,ISBIG34(I)) = EXPSBIG(J,ISBIG34(I))
     1                          + MEXS34(J)*ZMUL
        ENDDO
      ENDIF
      ENDDO
      RETURN
      END

C***********************************************************************
C
C     This subroutine processes the EAST interaction lists.
C
C     INPUT:
C
C     IEALL(NEALL), IYEALL(NEALL) are the boxes 
C          receiving all child box data and the X and Y offsets from
C          child 1, respectively.
C     The other east lists are similarly defined (see MKEALIST).
C
C     MEXEALL is the exponential expansion for all eight children, etc.
C     ZS are the shift coefficients computed by subroutine
C          MKEXPS.
C
C     OUTPUT:
C
C     LEXP1, which contains the local EAST expansion information for
C            all boxes, is incremented for each box in the 
C            interaction lists.
C
C
C***********************************************************************
      SUBROUTINE PROCESSEA(LEXP1,IEALL,NEALL,IYEALL,
     1           IE13,NE13,IY13,IE1,NE1,IY1,IE3,NE3,IY3,
     2           MEXEALL,MEXE13,MEXE1,MEXE3,ZS,NNODES,
     3           IEBIG13,EXPEBIG,IEBIG1,IEBIG3,ZSEAST,
     4           LOCALONOFF)
      IMPLICIT NONE
      INTEGER  IEALL(NEALL),NEALL,IYEALL(NEALL)
      INTEGER  IE13(NE13),NE13,IY13(NE13)
      INTEGER  IE1(NE1),NE1,IY1(NE1)
      INTEGER  IE3(NE3),NE3,IY3(NE3)
      INTEGER  IEBIG13(1)
      INTEGER  IEBIG1(1)
      INTEGER  IEBIG3(1)
      INTEGER  I, J, LOCALONOFF(1)
      INTEGER  NNODES
      COMPLEX *16  LEXP1(0:NNODES,1)
      COMPLEX *16  MEXEALL(0:NNODES)
      COMPLEX *16  MEXE13(0:NNODES)
      COMPLEX *16  MEXE1(0:NNODES)
      COMPLEX *16  MEXE3(0:NNODES)
      COMPLEX *16  ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZSEAST(-3:3,-3:3,0:NNODES)
      COMPLEX *16  EXPEBIG(0:NNODES,1)
      COMPLEX *16  ZMUL
C
      DO I = 1,NEALL
        IF(LOCALONOFF(IEALL(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(3,IYEALL(I),J)
            LEXP1(J,IEALL(I)) = LEXP1(J,IEALL(I)) +
     1            MEXEALL(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NE13
        IF(LOCALONOFF(IE13(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(2,IY13(I),J)
            LEXP1(J,IE13(I)) = LEXP1(J,IE13(I)) +
     1            MEXE13(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NE1
        IF(LOCALONOFF(IE1(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(2,IY1(I),J)
            LEXP1(J,IE1(I)) = LEXP1(J,IE1(I)) +
     1            MEXE1(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NE3
        IF(LOCALONOFF(IE3(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(2,IY3(I),J)
            LEXP1(J,IE3(I)) = LEXP1(J,IE3(I)) +
     1            MEXE3(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO

      I = 1
      IF(IEBIG13(I) .GT. 0 .AND.
     1   LOCALONOFF(IEBIG13(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSEAST(2,0,J)
          EXPEBIG(J,IEBIG13(I)) = EXPEBIG(J,IEBIG13(I)) +
     1          MEXE13(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSEAST(2,0,J)
          EXPEBIG(J,IEBIG13(I)) = EXPEBIG(J,IEBIG13(I)) +
     1          MEXE13(J)*ZMUL
        ENDDO
      ENDIF

      I = 1
      IF(IEBIG3(I) .GT. 0 .AND.
     1   LOCALONOFF(IEBIG3(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSEAST(2,2,J)
          EXPEBIG(J,IEBIG3(I)) = EXPEBIG(J,IEBIG3(I)) +
     1          MEXE3(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSEAST(2,2,J)
          EXPEBIG(J,IEBIG3(I)) = EXPEBIG(J,IEBIG3(I)) +
     1          MEXE3(J)*ZMUL
        ENDDO
      ENDIF

      I = 1
      IF(IEBIG1(I) .GT. 0 .AND.
     1   LOCALONOFF(IEBIG1(I)).EQ. 1)THEN
        J = 0
         ZMUL = ZSEAST(2,-2,J)
          EXPEBIG(J,IEBIG1(I)) = EXPEBIG(J,IEBIG1(I)) +
     1          MEXE1(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSEAST(2,-2,J)
         EXPEBIG(J,IEBIG1(I)) = EXPEBIG(J,IEBIG1(I)) +
     1          MEXE1(J)*ZMUL
        ENDDO
      ENDIF
      RETURN
      END


C***********************************************************************
C
C     This subroutine processes the WEST interaction lists.
C
C     INPUT:
C
C     IWALL(NWALL), IYWALL(NWALL) are the boxes 
C          receiving all child box data and the X and Y offsets from
C          child 1, respectively.
C     The other west lists are similarly defined (see MKWELIST).
C
C     MEXEALL is the exponential expansion for all eight children, etc.
C     ZS are the shift coefficients computed by subroutine
C          MKEXPS.
C
C     OUTPUT:
C
C     LEXP1, which contains the local WEST expansion information for
C            all boxes, is incremented for each box in the 
C            interaction lists.
C
C
C***********************************************************************
      SUBROUTINE PROCESSWE(LEXP2,IWALL,NWALL,IYWALL,
     1           IW24,NW24,IY24,IW2,NW2,IY2,IW4,NW4,IY4,
     2           MEXWALL,MEXW24,MEXW2,MEXW4,ZS,NNODES,
     3           IWBIG24,EXPWBIG,IWBIG2,IWBIG4,ZSWEST,
     4           LOCALONOFF)
      IMPLICIT NONE
      INTEGER  IWALL(NWALL),NWALL,IYWALL(NWALL)
      INTEGER  IW24(NW24),NW24,IY24(NW24)
      INTEGER  IW2(NW2),NW2,IY2(NW2)
      INTEGER  IW4(NW4),NW4,IY4(NW4)
      INTEGER  I, J, LOCALONOFF(1)
      INTEGER  NNODES
      INTEGER  IWBIG24(1)
      INTEGER  IWBIG2(1)
      INTEGER  IWBIG4(1)
      COMPLEX *16  LEXP2(0:NNODES,1)
      COMPLEX *16  MEXWALL(0:NNODES)
      COMPLEX *16  MEXW24(0:NNODES)
      COMPLEX *16  MEXW2(0:NNODES)
      COMPLEX *16  MEXW4(0:NNODES)
      COMPLEX *16  ZS(-3:3,-3:3,0:NNODES)
      COMPLEX *16  ZSWEST(-3:3,-3:3,0:NNODES)
      COMPLEX *16  EXPWBIG(0:NNODES,1)
      COMPLEX *16  ZMUL
C
      DO I = 1,NWALL
        IF(LOCALONOFF(IWALL(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(2,-IYWALL(I),J)
            LEXP2(J,IWALL(I)) = LEXP2(J,IWALL(I)) +
     1            MEXWALL(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NW24
        IF(LOCALONOFF(IW24(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(1,-IY24(I),J)
            LEXP2(J,IW24(I)) = LEXP2(J,IW24(I)) +
     1            MEXW24(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NW2
        IF(LOCALONOFF(IW2(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(1,-IY2(I),J)
            LEXP2(J,IW2(I)) = LEXP2(J,IW2(I)) +
     1            MEXW2(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO
C
      DO I = 1,NW4
        IF(LOCALONOFF(IW4(I)) .EQ. 1)THEN
         DO J = 0,NNODES
            ZMUL = ZS(1,-IY4(I),J)
            LEXP2(J,IW4(I)) = LEXP2(J,IW4(I)) +
     1            MEXW4(J)*ZMUL
         ENDDO
        ENDIF
      ENDDO

      I = 1
      IF(IWBIG24(I) .GT. 0 .AND. 
     1   LOCALONOFF(IWBIG24(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSWEST(2,0,J)
          EXPWBIG(J,IWBIG24(I)) = EXPWBIG(J,IWBIG24(I)) +
     1          MEXW24(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSWEST(2,0,J)
          EXPWBIG(J,IWBIG24(I)) = EXPWBIG(J,IWBIG24(I)) +
     1          MEXW24(J)*ZMUL
        ENDDO
      ENDIF

      I = 1
      IF(IWBIG2(I) .GT. 0 .AND.
     1   LOCALONOFF(IWBIG2(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSWEST(2,2,J)
          EXPWBIG(J,IWBIG2(I)) = EXPWBIG(J,IWBIG2(I)) +
     1          MEXW2(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSWEST(2,2,J)
          EXPWBIG(J,IWBIG2(I)) = EXPWBIG(J,IWBIG2(I)) +
     1          MEXW2(J)*ZMUL
        ENDDO
      ENDIF


      I = 1
      IF(IWBIG4(I) .GT. 0 .AND.
     1   LOCALONOFF(IWBIG4(I)) .EQ. 1)THEN
        J = 0
         ZMUL = ZSWEST(2,-2,J)
          EXPWBIG(J,IWBIG4(I)) = EXPWBIG(J,IWBIG4(I)) +
     1          MEXW4(J)*ZMUL

        DO J = 1, NNODES
         ZMUL = ZSWEST(2,-2,J)
          EXPWBIG(J,IWBIG4(I)) = EXPWBIG(J,IWBIG4(I)) +
     1          MEXW4(J)*ZMUL
        ENDDO
      ENDIF
      RETURN
      END




C************************************************************************
C     The following subroutine calculates the multipole
C     coefficients for a given set of polynomial coefficients.
C
C     INPUT:
C
C     COEFFS is the array of coefficients for the basis functions
C
C     NDEG is the degree of the approximating polynomial
C
C     NTERMS is the number order of the multipole expansion
C
C     XLENGTH is the length of a side of the box
C
C     W are the basis function to multipole weights.
C
C     OUTPUT:
C
C     MPOLE is the array of multipole coefficients for the box
C
C************************************************************************
      SUBROUTINE MULTIPOLE(COEFF, NDEG, NTERMS, XLENGTH, MPOLE, W)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, NTERMS
      REAL *8 COEFF(0:NDEG,0:NDEG), XLENGTH
      COMPLEX *16 MPOLE(0:NTERMS), W(0:NTERMS,0:NDEG,0:NDEG)
C-----Local variables
      INTEGER I, J, K
      REAL *8 COEFFSC(0:7,0:7), Y
      REAL *8 PI2
      COMPLEX *16 SUM

C     W contains a set of precomputed integrals for each of the
C     basis functions being considered.  W(i,j,k) represents the
C     integral of X^i * Y^j /(X + imag*Y)^k over the square 
C     -1 < X < 1 and -1 < Y < 1.
C     These integrals can be rescaled as needed.
C     The multipole coefficients are obtained by summing over
C     all of the basis functions and multiplying the above
C     rescaled integral by the corresponding polynomial
C     coefficient.

 
      PI2 = 2.0D0*DACOS(-1.0D0)

      DO I = 0, NDEG
        DO J = 0, NDEG - I
          COEFFSC(I,J) = COEFF(I,J) * XLENGTH**2
        END DO
      END DO

      SUM = 0.0D0
      DO I = 0, NDEG
        DO J = 0, NDEG-I
         SUM = SUM + COEFFSC(I,J) * W(0,I,J)
        END DO
      END DO
      MPOLE(0) = SUM/PI2


      Y =  -1.0D0 / (PI2*2.0D0)
      DO K = 1, NTERMS
       SUM = 0.0D0
       DO I = 0, NDEG
         DO J = 0, NDEG-I
           SUM = SUM  + COEFFSC(I,J) * W(K,I,J)
         END DO
       END DO
       MPOLE(K) = Y * SUM / DBLE(K)
       Y = Y / 2.0D0
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine does the local work (just outgoing)
C     between two boxes that are of the same size.  This routine does
C     exactly the same type of work as is doe for the local work in the
C     uniform case.
C
C     INPUT:
C
C     COEFFS is the array of  basis functions coefficients for the box
C                 being considered
C
C     NDEG is the degree of the approximating polynomial
C
C     NB denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the ADAPFMM6 routine)
C
C     W is the table of weights needed to determine the local contributions
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C**************************************************************************
      SUBROUTINE COLLOC8(POT,COEFFS,NDEG,NB,T,TLOG,W,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, NB, ISWITCH
      REAL *8 POT(64), COEFFS(0:NDEG,0:NDEG)
      REAL *8 T(3), W(9,64,36)
      REAL *8 TLOG(36)
C-----Local variables
      INTEGER I, J
      REAL *8 C(36), CLOG(36)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF
  
       C(1)  = COEFFS(0,0) / T(1)
       C(2)  = COEFFS(1,0) / T(1)
       C(3)  = COEFFS(0,1) / T(1)
       C(4)  = COEFFS(2,0) / T(1)
       C(5)  = COEFFS(1,1) / T(1)
       C(6)  = COEFFS(0,2) / T(1)
       C(7)  = COEFFS(3,0) / T(1)
       C(8)  = COEFFS(2,1) / T(1)
       C(9)  = COEFFS(1,2) / T(1)
       C(10) = COEFFS(0,3) / T(1)
       C(11) = COEFFS(4,0) / T(1)
       C(12) = COEFFS(3,1) / T(1)
       C(13) = COEFFS(2,2) / T(1)
       C(14) = COEFFS(1,3) / T(1)
       C(15) = COEFFS(0,4) / T(1)
       C(16) = COEFFS(5,0) / T(1)
       C(17) = COEFFS(4,1) / T(1)
       C(18) = COEFFS(3,2) / T(1)
       C(19) = COEFFS(2,3) / T(1)
       C(20) = COEFFS(1,4) / T(1)
       C(21) = COEFFS(0,5) / T(1)
       C(22) = COEFFS(6,0) / T(1)
       C(23) = COEFFS(5,1) / T(1)
       C(24) = COEFFS(4,2) / T(1)
       C(25) = COEFFS(3,3) / T(1)
       C(26) = COEFFS(2,4) / T(1)
       C(27) = COEFFS(1,5) / T(1)
       C(28) = COEFFS(0,6) / T(1)
       C(29) = COEFFS(7,0) / T(1)
       C(30) = COEFFS(6,1) / T(1)
       C(31) = COEFFS(5,2) / T(1)
       C(32) = COEFFS(4,3) / T(1)
       C(33) = COEFFS(3,4) / T(1)
       C(34) = COEFFS(2,5) / T(1)
       C(35) = COEFFS(1,6) / T(1)
       C(36) = COEFFS(0,7) / T(1)

 
       CLOG(1)  = C(1) * TLOG(1)
       CLOG(2)  = 0.0D0
       CLOG(3)  = 0.0D0
       CLOG(4)  = C(4) * TLOG(4)
       CLOG(5)  = 0.0D0
       CLOG(6)  = C(6) * TLOG(6)
       CLOG(7)  = 0.0D0
       CLOG(8)  = 0.0D0
       CLOG(9)  = 0.0D0
       CLOG(10) = 0.0D0
       CLOG(11) = C(11) * TLOG(11)
       CLOG(12) = 0.0D0
       CLOG(13) = C(13) * TLOG(13)
       CLOG(14) = 0.0D0
       CLOG(15) = C(15) * TLOG(15) 
       CLOG(16) = 0.0D0
       CLOG(17) = 0.0D0
       CLOG(18) = 0.0D0
       CLOG(19) = 0.0D0
       CLOG(20) = 0.0D0
       CLOG(21) = 0.0D0
       CLOG(22) = C(22) * TLOG(22)
       CLOG(23) = 0.0D0
       CLOG(24) = C(24) * TLOG(24)
       CLOG(25) = 0.0D0
       CLOG(26) = C(26) * TLOG(26)
       CLOG(27) = 0.0D0
       CLOG(28) = C(28) * TLOG(28)
       CLOG(29) = 0.0D0
       CLOG(30) = 0.0D0
       CLOG(31) = 0.0D0
       CLOG(32) = 0.0D0
       CLOG(33) = 0.0D0
       CLOG(34) = 0.0D0
       CLOG(35) = 0.0D0
       CLOG(36) = 0.0D0
 

       DO I = 1,64
        DO J = 1,36
          POT(I) = POT(I) + C(J)*W(NB,I,J) + CLOG(J)
        ENDDO
       ENDDO
       RETURN
       END


C***************************************************************************
C     The following subroutine does the local work that figures the
C     potential in a bigger box due to a small box that is touching it.
C     The boxes are only allowed to be one level apart.
C
C     INPUT:
C
C     COEFFS is the array of  basis functions coefficients for the box
C                 being considered
C
C     NDEG is the degree of the approximating polynomial
C
C     IFLAG denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the ADAPFMM6 routine)
C
C     W is the table of weights needed to determine the local contributions
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C***************************************************************************
      SUBROUTINE STOBLOC8(POT,COEFFS,NDEG,IFLAG,T,TLOG,W,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, IFLAG, ISWITCH
      REAL *8 POT(36), COEFFS(0:NDEG,0:NDEG)
      REAL *8 T(3), W(12,64,36)
      REAL *8 TLOG(36)
C-----Local variables
      INTEGER I, J
      REAL *8 C(36), CLOG(36)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF

       C(1)  = COEFFS(0,0) / T(3)
       C(2)  = COEFFS(1,0) / T(3)
       C(3)  = COEFFS(0,1) / T(3)
       C(4)  = COEFFS(2,0) / T(3)
       C(5)  = COEFFS(1,1) / T(3)
       C(6)  = COEFFS(0,2) / T(3)
       C(7)  = COEFFS(3,0) / T(3)
       C(8)  = COEFFS(2,1) / T(3)
       C(9)  = COEFFS(1,2) / T(3)
       C(10) = COEFFS(0,3) / T(3)
       C(11) = COEFFS(4,0) / T(3)
       C(12) = COEFFS(3,1) / T(3)
       C(13) = COEFFS(2,2) / T(3)
       C(14) = COEFFS(1,3) / T(3)
       C(15) = COEFFS(0,4) / T(3)
       C(16) = COEFFS(5,0) / T(3)
       C(17) = COEFFS(4,1) / T(3)
       C(18) = COEFFS(3,2) / T(3)
       C(19) = COEFFS(2,3) / T(3)
       C(20) = COEFFS(1,4) / T(3)
       C(21) = COEFFS(0,5) / T(3)
       C(22) = COEFFS(6,0) / T(3)
       C(23) = COEFFS(5,1) / T(3)
       C(24) = COEFFS(4,2) / T(3)
       C(25) = COEFFS(3,3) / T(3)
       C(26) = COEFFS(2,4) / T(3)
       C(27) = COEFFS(1,5) / T(3)
       C(28) = COEFFS(0,6) / T(3)
       C(29) = COEFFS(7,0) / T(3)
       C(30) = COEFFS(6,1) / T(3)
       C(31) = COEFFS(5,2) / T(3)
       C(32) = COEFFS(4,3) / T(3)
       C(33) = COEFFS(3,4) / T(3)
       C(34) = COEFFS(2,5) / T(3)
       C(35) = COEFFS(1,6) / T(3)
       C(36) = COEFFS(0,7) / T(3)

       CLOG(1)  = C(1) * TLOG(1)
       CLOG(2)  = 0.0D0
       CLOG(3)  = 0.0D0
       CLOG(4)  = C(4) * TLOG(4)
       CLOG(5)  = 0.0D0
       CLOG(6)  = C(6) * TLOG(6)
       CLOG(7)  = 0.0D0
       CLOG(8)  = 0.0D0
       CLOG(9)  = 0.0D0
       CLOG(10) = 0.0D0
       CLOG(11) = C(11) * TLOG(11)
       CLOG(12) = 0.0D0
       CLOG(13) = C(13) * TLOG(13)
       CLOG(14) = 0.0D0
       CLOG(15) = C(15) * TLOG(15)
       CLOG(16) = 0.0D0
       CLOG(17) = 0.0D0
       CLOG(18) = 0.0D0
       CLOG(19) = 0.0D0
       CLOG(20) = 0.0D0
       CLOG(21) = 0.0D0
       CLOG(22) = C(22) * TLOG(22)
       CLOG(23) = 0.0D0
       CLOG(24) = C(24) * TLOG(24)
       CLOG(25) = 0.0D0
       CLOG(26) = C(26) * TLOG(26)
       CLOG(27) = 0.0D0
       CLOG(28) = C(28) * TLOG(28)
       CLOG(29) = 0.0D0
       CLOG(30) = 0.0D0
       CLOG(31) = 0.0D0
       CLOG(32) = 0.0D0
       CLOG(33) = 0.0D0
       CLOG(34) = 0.0D0
       CLOG(35) = 0.0D0
       CLOG(36) = 0.0D0

       DO I = 1, 64
         DO J = 1, 36
           POT(I) = POT(I) + C(J)*W(IFLAG,I,J) + CLOG(J)
         END DO
       END DO
       RETURN
       END


C***************************************************************************
C     The following subroutine does the local work that figures the
C     potential in a smaller box due to a big box that is touching it.
C     The boxes are only allowed to be one level apart.
C
C     INPUT:
C
C     COEFFS is the array of  basis functions coefficients for the box
C                 being considered
C
C     NDEG is the degree of the approximating polynomial
C
C     IFLAG denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the ADAPFMM6 routine)
C
C     W is the table of weights needed to determine the local contributions
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C***************************************************************************
      SUBROUTINE BTOSLOC8(POT,COEFFS,NDEG,IFLAG,T,TLOG,W,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, IFLAG, ISWITCH
      REAL *8 POT(36), COEFFS(0:NDEG,0:NDEG), T(3), TLOG(36)
      REAL *8 W(12,64,36)
C-----Local variables
      INTEGER I, J
      REAL *8 C(36), CLOG(36) 

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF
       
       C(1)  = COEFFS(0,0)/T(1)
       C(2)  = COEFFS(1,0)/T(1)
       C(3)  = COEFFS(0,1)/T(1)
       C(4)  = COEFFS(2,0)/T(1)
       C(5)  = COEFFS(1,1)/T(1)
       C(6)  = COEFFS(0,2)/T(1)
       C(7)  = COEFFS(3,0)/T(1)
       C(8)  = COEFFS(2,1)/T(1)
       C(9)  = COEFFS(1,2)/T(1)
       C(10) = COEFFS(0,3)/T(1)
       C(11) = COEFFS(4,0)/T(1)
       C(12) = COEFFS(3,1)/T(1)
       C(13) = COEFFS(2,2)/T(1)
       C(14) = COEFFS(1,3)/T(1)
       C(15) = COEFFS(0,4)/T(1)
       C(16) = COEFFS(5,0)/T(1)
       C(17) = COEFFS(4,1)/T(1)
       C(18) = COEFFS(3,2)/T(1)
       C(19) = COEFFS(2,3)/T(1)
       C(20) = COEFFS(1,4)/T(1)
       C(21) = COEFFS(0,5)/T(1)
       C(22) = COEFFS(6,0)/T(1)
       C(23) = COEFFS(5,1)/T(1)
       C(24) = COEFFS(4,2)/T(1)
       C(25) = COEFFS(3,3)/T(1)
       C(26) = COEFFS(2,4)/T(1)
       C(27) = COEFFS(1,5)/T(1)
       C(28) = COEFFS(0,6)/T(1)
       C(29) = COEFFS(7,0)/T(1)
       C(30) = COEFFS(6,1)/T(1)
       C(31) = COEFFS(5,2)/T(1)
       C(32) = COEFFS(4,3)/T(1)
       C(33) = COEFFS(3,4)/T(1)
       C(34) = COEFFS(2,5)/T(1)
       C(35) = COEFFS(1,6)/T(1)
       C(36) = COEFFS(0,7)/T(1)

       CLOG(1)  = C(1) * TLOG(1)
       CLOG(2)  = 0.0D0
       CLOG(3)  = 0.0D0
       CLOG(4)  = C(4) * TLOG(4)
       CLOG(5)  = 0.0D0
       CLOG(6)  = C(6) * TLOG(6)
       CLOG(7)  = 0.0D0
       CLOG(8)  = 0.0D0
       CLOG(9)  = 0.0D0
       CLOG(10) = 0.0D0
       CLOG(11) = C(11) * TLOG(11)
       CLOG(12) = 0.0D0
       CLOG(13) = C(13) * TLOG(13)
       CLOG(14) = 0.0D0
       CLOG(15) = C(15) * TLOG(15)
       CLOG(16) = 0.0D0
       CLOG(17) = 0.0D0
       CLOG(18) = 0.0D0
       CLOG(19) = 0.0D0
       CLOG(20) = 0.0D0
       CLOG(21) = 0.0D0
       CLOG(22) = C(22) * TLOG(22)
       CLOG(23) = 0.0D0
       CLOG(24) = C(24) * TLOG(24)
       CLOG(25) = 0.0D0
       CLOG(26) = C(26) * TLOG(26)
       CLOG(27) = 0.0D0
       CLOG(28) = C(28) * TLOG(28)
       CLOG(29) = 0.0D0
       CLOG(30) = 0.0D0
       CLOG(31) = 0.0D0
       CLOG(32) = 0.0D0
       CLOG(33) = 0.0D0
       CLOG(34) = 0.0D0
       CLOG(35) = 0.0D0
       CLOG(36) = 0.0D0
 
      DO I = 1, 64
        DO J = 1, 36
           POT(I) = POT(I) + C(J)*W(IFLAG,I,J) + CLOG(J)
        END DO
      END DO
      RETURN
      END


C************************************************************************
C     The following subroutine takes a specified version of the multipole
C     expansion for box B and converts it to a local expansion about
C     one of the * boxes.  The specified multipole expansion is the
C     one obtained by interpolating down on B to four 'ghost' boxes.
C     This is explained in detail in MKEXPBTOS.
C           -------------------------
C           | * | * | * | * | * | * |
C           |___|___|___|___|___|___|
C           |   |   |   |   |   |   |
C           | * |   |   |   |   | * |
C           -------------------------
C           | * |   |       |   | * |
C           |___|___|   B   |___|___|
C           |   |   |       |   | * |
C           | * |   |       |   |   |
C           -------------------------
C           | * |   |   |   |   | * |
C           |___|___|___|___|___|___|
C           |   |   |   |   |   |   |
C           | * | * | * | * | * | * |
C           -------------------------
C
C     INPUT:
C
C     EXPON is the smaller boxes exponential expansion
C     
C     Z is the distance from the center of box B to the center
C     * target box,  where the scaling units are that the side of
C     a * box is of length one.
C     
C     NTERMS is the order of the multipole expansion
C
C     XNODES is a blank array that is set to
C            the nodes in the plane wave expansions
C
C     COMP is a precomputed term involving the ratio of
C          the weights and factorial terms
C
C     BETALL is the larger boxes exponential expansion that
C           has been precomputed
C     
C     SPIN is a parameter that determines what direction the exponential
C          expansions decay in
C
C     OUTPUT:
C
C     EXPON is altered to account for contribution of the larger box
C
C************************************************************************
      SUBROUTINE BTOSFAR(EXPON, Z, NTERMS, XNODES,
     1                 COMP, BETALL, SPIN,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS, ISWITCH
      REAL *8  XNODES(1), COMP(NTERMS,0:NTERMS)
      COMPLEX *16  EXPON(0:NTERMS), BETALL(0:NTERMS)
      COMPLEX *16  SPIN, Z
C-----Local variables
      INTEGER  I
      COMPLEX *16  EXPALL(0:60)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF

        DO I = 0, NTERMS
          EXPALL(I) = BETALL(I)
        END DO

        CALL EXPSHIFT(EXPALL, NTERMS, XNODES, SPIN, Z)
        CALL ADDEXP(EXPALL, EXPON, NTERMS)

      RETURN
      END


C*******************************************************************
C     The following subroutine precomputes some stuff that is needed
C     in the 'big to small far' routine later on.  Let's suppose
C     that we have the following picture:
C           -------------------------
C           | * | * | * | * | * | * |
C           |___|___|___|___|___|___|
C           |   |   |   |   |   |   |
C           | * |   |   |   |   | * |
C           -------------------------
C           | * |   |       |   | * |
C           |___|___|   B   |___|___|
C           |   |   |       |   | * |
C           | * |   |       |   |   |
C           -------------------------
C           | * |   |   |   |   | * |
C           |___|___|___|___|___|___|
C           |   |   |   |   |   |   |
C           | * | * | * | * | * | * |
C           -------------------------
C     and we need to get compute the interaction from box B to
C     any one of the * boxes.  Note that from the perspective
C     of B, the * boxes are not well separated, and also note
C     that there is no reason that one of the * boxes could
C     not be divided.  The expansions are obtained using a matrix
C     that was precomputed earlier.
C
C     INPUT:
C
C     COEFFS are the basis function coefficients for the larger box
C
C     MAPSOUTH maps from the basis function coefficients to the 
C              south exponential coefficients
C
C     MAPNORTH maps from the basis function coefficients to the 
C              north exponential coefficients
C
C     MAPEAST maps from the basis function coefficients to the 
C              east exponential coefficients
C
C     MAPWEST maps from the basis function coefficients to the 
C              west exponential coefficients
C
C     SUM2 is a precomputed terms needed in the exponential expansions
C
C     XLENGTH is the length of the larger box
C
C     NNODES is the number of nodes in the exponential expansion
C
C     OUTPUT:
C
C     EXPNALL is north expansion for the larger box
C
C     EXPSALL is south expansion for the larger box
C
C     EXPEALL is east expansion for the larger box
C
C     EXPWALL is west expansion for the larger box
C
C*******************************************************************
      SUBROUTINE MKEXPBTOS8(COEFFS, EXPNALL,EXPSALL,EXPEALL,EXPWALL,
     1      MAPSOUTH,MAPNORTH,MAPEAST,MAPWEST,SUM2, XLENGTH, NNODES)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES
      REAL *8  COEFFS(0:7,0:7), SUM2,  XLENGTH
      COMPLEX *16  EXPNALL(0:NNODES), EXPSALL(0:NNODES)
      COMPLEX *16  EXPEALL(0:NNODES), EXPWALL(0:NNODES)
      COMPLEX *16  MAPWEST(0:NNODES,36), MAPSOUTH(0:NNODES,36)
      COMPLEX *16  MAPEAST(0:NNODES,36), MAPNORTH(0:NNODES,36)
C-----Local variables
      INTEGER  I, K
      REAL *8  TCOEFF(36)
      REAL *8  TEMP1, TEMP2
      COMPLEX *16  SUM

C     Now initialize all of the output
C     arrays of coefficients to zero:
      DO I = 0, NNODES
        EXPEALL(I) = 0.0D0
        EXPWALL(I) = 0.0D0
        EXPNALL(I) = 0.0D0
        EXPSALL(I) = 0.0D0
      END DO

C     Now let's place the coefficients in a
C     one dimensional array that is indexed 
C     correctly:
      TEMP1 = XLENGTH**2

      TCOEFF(1) = COEFFS(0,0) * TEMP1
      TCOEFF(2) = COEFFS(1,0) * TEMP1
      TCOEFF(3) = COEFFS(0,1) * TEMP1
      TCOEFF(4) = COEFFS(2,0) * TEMP1
      TCOEFF(5) = COEFFS(1,1) * TEMP1
      TCOEFF(6) = COEFFS(0,2) * TEMP1
      TCOEFF(7) = COEFFS(3,0) * TEMP1
      TCOEFF(8) = COEFFS(2,1) * TEMP1
      TCOEFF(9) = COEFFS(1,2) * TEMP1
      TCOEFF(10)= COEFFS(0,3) * TEMP1
      TCOEFF(11)= COEFFS(4,0) * TEMP1
      TCOEFF(12)= COEFFS(3,1) * TEMP1
      TCOEFF(13)= COEFFS(2,2) * TEMP1
      TCOEFF(14)= COEFFS(1,3) * TEMP1
      TCOEFF(15)= COEFFS(0,4) * TEMP1
      TCOEFF(16)= COEFFS(5,0) * TEMP1
      TCOEFF(17)= COEFFS(4,1) * TEMP1
      TCOEFF(18)= COEFFS(3,2) * TEMP1
      TCOEFF(19)= COEFFS(2,3) * TEMP1
      TCOEFF(20)= COEFFS(1,4) * TEMP1
      TCOEFF(21)= COEFFS(0,5) * TEMP1
      TCOEFF(22)= COEFFS(6,0) * TEMP1
      TCOEFF(23)= COEFFS(5,1) * TEMP1
      TCOEFF(24)= COEFFS(4,2) * TEMP1
      TCOEFF(25)= COEFFS(3,3) * TEMP1
      TCOEFF(26)= COEFFS(2,4) * TEMP1
      TCOEFF(27)= COEFFS(1,5) * TEMP1
      TCOEFF(28)= COEFFS(0,6) * TEMP1
      TCOEFF(29)= COEFFS(7,0) * TEMP1
      TCOEFF(30)= COEFFS(6,1) * TEMP1
      TCOEFF(31)= COEFFS(5,2) * TEMP1
      TCOEFF(32)= COEFFS(4,3) * TEMP1
      TCOEFF(33)= COEFFS(3,4) * TEMP1
      TCOEFF(34)= COEFFS(2,5) * TEMP1
      TCOEFF(35)= COEFFS(1,6) * TEMP1
      TCOEFF(36)= COEFFS(0,7) * TEMP1

C     Now just multiply by a precomputed matrix that
C     will map from the ten polynomial coefficients
C     to the four expansions (they are centered 
C     around the 'child' box in the lower left
C     hand corner of the parent box)
      DO K = 0, NNODES
       SUM = 0.0D0
       DO I = 1, 36
            SUM = SUM + MAPNORTH(K,I)*TCOEFF(I)
       END DO
       EXPNALL(K) = SUM
      END DO

      DO K = 0, NNODES
       SUM = 0.0D0
       DO I = 1, 36
            SUM = SUM + MAPWEST(K,I)*TCOEFF(I)
       END DO
       EXPWALL(K) = SUM
      END DO

      DO K = 0, NNODES
       SUM = 0.0D0
       DO I = 1, 36
            SUM = SUM + MAPEAST(K,I)*TCOEFF(I)
       END DO
       EXPEALL(K) = SUM
      END DO

      DO K = 0, NNODES
       SUM = 0.0D0
       DO I = 1, 36
            SUM = SUM + MAPSOUTH(K,I)*TCOEFF(I)
       END DO
       EXPSALL(K) = SUM
      END DO

      TEMP2 = SUM2 + DLOG(XLENGTH/2.0D0)
      EXPEALL(0) = EXPEALL(0) * TEMP2
      EXPWALL(0) = EXPWALL(0) * TEMP2
      EXPNALL(0) = EXPNALL(0) * TEMP2
      EXPSALL(0) = EXPSALL(0) * TEMP2
      RETURN
      END


C**********************************************************************
C     The following subroutine is designed to evaluate the
C     multipole expansion given an array of multipole coefficients.
C
C     INPUT:
C
C     MPOLE is the array of multipole coefficients being evaluated
C
C     NTERMS is the order of the multipole expansion
C
C     Z is the target point (scaled from the standpoint of the side of
C                       the above box being one)
C
C     XLENGTH is the length of a side of the box being considered
C
C     OUTPUT:
C
C     PHI is value of the multipole expansion evaluated at Z
C
C**********************************************************************
      SUBROUTINE MULTEVAL(MPOLE, NTERMS, Z, PHI, XLENGTH)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS
      REAL *8  XLENGTH 
      COMPLEX *16  MPOLE(0:NTERMS), PHI, Z
C-----Local variables
      INTEGER  K
      COMPLEX *16  ZK

C      The basic formula for PHI is:
C      PHI = MPOLE(0)*log(Z) - sum(K=1 to NTERMS)(MPOLE(K)*(1/Z^K)
C      The scaling factor is included below.
 
C      ZK is just set to be Z^K within the loop below.
        PHI = MPOLE(0) * (CDLOG(Z) + DLOG(XLENGTH))
        ZK = Z
        DO K = 1, NTERMS
          PHI = PHI  +  MPOLE(K) / ZK
          ZK = ZK * Z
        END DO

      RETURN
      END


C***********************************************************************
C     The following subroutine is designed to shift an exponential
C     expansion to a new target point, where it can then be converted
C     to a local expansion.
C
C     INPUT:
C
C     BETA is the array of plane-wave coefficients being shifted
C
C     NNODES is the order of the plane-wave expansion
C
C     XNODES are the weights associated with the plane wave expansion
C
C     SPIN is a parameter that determines what direction the exponential
C          expansions decay in
C
C     ZSHIFT is the distance by which the expansion is to be shifted 
C                                      (scaled)
C
C     OUTPUT:
C
C     BETA is altered to account for the shift
C
C***********************************************************************
      SUBROUTINE EXPSHIFT(BETA, NNODES, XNODES, SPIN, ZSHIFT)
      IMPLICIT NONE
C-----Global variables
      INTEGER NNODES
      REAL *8 XNODES(NNODES)
      COMPLEX *16 BETA(0:NNODES), SPIN, ZSHIFT
C-----Local variables
      INTEGER I
      COMPLEX *16 ZSPIN
 
C      Now adjust the coefficients to account for the shift.
C      Leave the first one unchanged.
 
        ZSPIN = SPIN*ZSHIFT
        DO I = 1, NNODES
          BETA(I) = BETA(I) * CDEXP(ZSPIN*XNODES(I))
        END DO

      RETURN
      END


C**********************************************************************
      SUBROUTINE PSPPTA(A,B,NTERMS)
C**********************************************************************
      IMPLICIT NONE
      COMPLEX *16 A(0:1),B(0:1)
      REAL *8  C(100,100)
      REAL *8  ZERO,PI,ZPOW(25)
      REAL *8  CONPOT(0:50,0:50)
      INTEGER  NTERMS
      INTEGER  I,J,IPOW,IER
      DATA  ZERO/0.0D0/
      SAVE
C
C        This entry converts the given multipole expansion (A) into
C        a power series expansion (B) about the origin which describes
C        the potential due to all well-separated image boxes.
C
C     A      - The original multipole decomposition
C
C     NOTE:    We ignore A(0) = net charge since the resulting
C              potential is not meaningful.
C
C     B      - The evaluated power series coefficients
C
C     NTERMS - The order of the decompositions A and B
C     ____________________________________________________________
C
C
      PI = DACOS(-1.0D0)
      DO I=0,NTERMS-1
         B(I) = ZERO
         DO J=1,NTERMS-1
            B(I) = B(I) + CONPOT(I,J)*A(J)
         END DO
      END DO

         B(0) = B(0) + PI*CONJG(A(2))
         B(1) = B(1) - PI*CONJG(A(1))
      RETURN
C
C*****************************************************************
C
      ENTRY PSPPIN(IER)
C
C     This is the initialization entry point.
C     It precomputes the binomial coefficients and the potential
C     conversion matrices in full storage format.
C
C     Create the binomial coefficients
C                C(m,k) = (m-1) choose (k-1)
C
      C(1,1) = 1.0D0
      DO I = 2,100
         C(1,I) = 0.0D0
         C(I,1) = 1.0D0
      END DO

      DO I = 2,100
        DO J = 2,100
         C(I,J)=C(I-1,J-1)+C(I-1,J)
        END DO
      END DO
C
C----- create ZPOW, the scaled sums over S of 1/z^k. (The periodic
C      structure has width 64 instead of 1.)
C

      ZPOW(1) = 0.151212002153897D+00
      ZPOW(2) = 0.577303536518952D-02
      ZPOW(3) = 0.134901282797037D-02
      ZPOW(4) = 0.700330250248559D-04
      ZPOW(5) = 0.300317628955958D-05
      ZPOW(6) = 0.242803838628101D-06
      ZPOW(7) = 0.160994099657098D-07
      ZPOW(8) = 0.897580666668933D-09
      ZPOW(9) = 0.570451512710591D-10
      ZPOW(10)= 0.371803105861787D-11
      ZPOW(11)= 0.227440251455375D-12
      ZPOW(12)= 0.140812857755456D-13
      ZPOW(13)= 0.890974169335439D-15
      ZPOW(14)= 0.556558380974353D-16
      ZPOW(15)= 0.346173274816875D-17
      ZPOW(16)= 0.216781734001496D-18
      ZPOW(17)= 0.135661847686001D-19
      ZPOW(18)= 0.846820937604322D-21
      ZPOW(19)= 0.529224560383524D-22
      ZPOW(20)= 0.330944477656778D-23
      ZPOW(21)= 0.206806338091673D-24
      ZPOW(22)= 0.129232908060200D-25
      ZPOW(23)= 0.807807171359492D-27
      ZPOW(24)= 0.504890432194844D-28
      ZPOW(25)= 0.315537827943134D-29

C
C-----Generate the conversion matrices
C
      DO I=0,50
       DO J=1,50
        CONPOT(I,J) = 0.0D0
         DO IPOW = 1,25
           IF (I+J .EQ. 4*IPOW) THEN
            CONPOT(I,J) = ZPOW(IPOW)*C(I+J,J)*(-1)**J
           ENDIF
         END DO
       END DO
      END DO
      RETURN
      END



C***********************************************************************
C     The following subroutine is designed to copy the initial tree
C     structure.  In the cases where we need to alter the tree structure
C     and work on a different tree (I am referring here to cases where the
C     boundary conditions are being imposed)  the original structure is 
C     copied and the work is done on the copy.  The original structure is
C     returned to the used, unchanged, on output.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     IROWBOX denotes the row of each box
C
C     ICOLBOX denotes the column of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NLEV is the finest level
C
C     ICHILDBOX denotes the four children of each box
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     NBLEVEL is the total number of boxes per level
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C
C     OUTPUT: (all things in the output refer to the copy)
C
C     LEVELBOX1 is an array determining the level of each box
C
C     NBOXES1 is the total number of boxes
C
C     IROWBOX1 denotes the row of each box
C
C     ICOLBOX1 denotes the column of each box
C
C     IPARENTBOX1 denotes the parent of each box
C
C     NLEV1 is the finest level
C
C     ICHILDBOX1 denotes the four children of each box
C
C     IBOXLEV1 is the array in which the boxes are arranged
C
C     NBLEVEL1 is the total number of boxes per level
C
C     ISTARTLEV1 is the pointer to where each level begins in the
C               IBOXLEV array
C
C***********************************************************************
      SUBROUTINE COPY(LEVELBOX,LEVELBOX1,ICOLBOX,ICOLBOX1,
     1    IROWBOX,IROWBOX1,IPARENTBOX, IPARENTBOX1,
     2    ICHILDBOX,ICHILDBOX1,NBOXES,ISTARTLEV,ISTARTLEV1,
     3    IBOXLEV, IBOXLEV1, NBLEVEL, NBLEVEL1, NLEV, NLEV1,
     4    NBOXES1)
      IMPLICIT NONE
C-----Global variables
      INTEGER  LEVELBOX(1), LEVELBOX1(1)
      INTEGER  ICOLBOX(1), ICOLBOX1(1)
      INTEGER  IROWBOX(1), IROWBOX1(1)
      INTEGER  IPARENTBOX(1), IPARENTBOX1(1)
      INTEGER  ICHILDBOX(4,1), ICHILDBOX1(4,1)
      INTEGER  ISTARTLEV(0:1), NBLEVEL(0:1)
      INTEGER  IBOXLEV(1)
      INTEGER  ISTARTLEV1(0:1), NBLEVEL1(0:1)
      INTEGER  IBOXLEV1(1)
      INTEGER  NBOXES, NBOXES1, NLEV, NLEV1
C-----Local variables
      INTEGER  I, J

         NBOXES1 =  NBOXES
         NLEV1   =  NLEV
C     For the special cases, let's copy the workspace over
C     to a new one:
      DO I = 1, NBOXES
        LEVELBOX1(I)   = LEVELBOX(I)
        IPARENTBOX1(I) = IPARENTBOX(I)
        ICOLBOX1(I) = ICOLBOX(I)
        IROWBOX1(I) = IROWBOX(I)
        IBOXLEV1(I) = IBOXLEV(I)
        DO J = 1, 4
         ICHILDBOX1(J,I) = ICHILDBOX(J,I)
        END DO
      END DO

       DO I = 0, NLEV
         ISTARTLEV1(I) =   ISTARTLEV(I)
         NBLEVEL1(I)   =   NBLEVEL(I)
       END DO

      RETURN
      END

C***********************************************************************
C     The following subroutine takes one of the boundary conditions 
C     options as input and makes the necessary changes to the initial
C     box so that one of the boundary condition cases can be solved.
C     Also inputed are the tree structure and the right hand side,
C     FRIGHT.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NBOXES is the total number of boxes
C
C     IPERIOD denotes which of the boundary cases we are in
C
C     NLEV is the finest level
C
C     FRIGHT is the right hand side defined on the old grid
C
C     OUTPUT:
C
C     All of the above is altered to reflect the new,
C     four fold tree
C
C***********************************************************************
      SUBROUTINE FOLDOVER(LEVELBOX,ICOLBOX,IROWBOX,IPARENTBOX,
     1            ICHILDBOX,NBOXES, IPERIOD, NLEV, FRIGHT, ISTARTLEV,
     2            IBOXLEV, NBLEVEL)
      IMPLICIT NONE
C-----Global variables
      INTEGER  LEVELBOX(1)
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  IPARENTBOX(1), ICHILDBOX(4,1)
      INTEGER  NBOXES, NLEV, IPERIOD
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      REAL *8  FRIGHT(64,1)
C-----Local variables
      INTEGER  IPLUS

      IF(IPERIOD .EQ. 2 .OR. IPERIOD .EQ. 7)THEN

        IPLUS = -1
        CALL FLIPUP8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV,IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS, ISTARTLEV, IBOXLEV, NBLEVEL)

        IPLUS = -1
        CALL FLIPRIGHT8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS)

      ELSEIF(IPERIOD .EQ. 3 .OR. IPERIOD .EQ. 8)THEN

        IPLUS = 1
        CALL FLIPUP8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV,IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS, ISTARTLEV, IBOXLEV, NBLEVEL)

        IPLUS = 1
        CALL FLIPRIGHT8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS)

      ELSEIF(IPERIOD .EQ. 4 .OR. IPERIOD .EQ. 9)THEN

        CALL SHIFTUP8(LEVELBOX, NBOXES,
     1   ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2   ICHILDBOX, FRIGHT, ISTARTLEV, IBOXLEV, NBLEVEL)

        IPLUS = -1
        CALL FLIPRIGHT8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS)

      ELSEIF(IPERIOD .EQ. 5 .OR. IPERIOD .EQ. 10)THEN

        IPLUS = 1
        CALL FLIPUP8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV,IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS, ISTARTLEV, IBOXLEV, NBLEVEL)

        IPLUS = -1
        CALL FLIPRIGHT8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS)


      ELSEIF(IPERIOD .EQ. 6 .OR. IPERIOD .EQ. 11)THEN

        CALL SHIFTUP8(LEVELBOX, NBOXES,
     1   ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2   ICHILDBOX, FRIGHT, ISTARTLEV, IBOXLEV, NBLEVEL)

        IPLUS = 1
        CALL FLIPRIGHT8(LEVELBOX, NBOXES,
     1    ICOLBOX, IROWBOX, NLEV, IPARENTBOX,
     2    ICHILDBOX, FRIGHT, IPLUS)

      ENDIF
      RETURN 
      END




C**********************************************************************
C     The following subroutine takes a tree structure and right hand
C     side as input and reflects it over the top, doubling
C     the size of the total number of boxes. 
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NLEV is the finest level
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     FRIGHT represents the right hand side defined on the leaf nodes
C
C     IPLUS is a parameter that determines whether or not the sign 
C     of the right hand side is reversed or stays the
C     same in its reflection
C
C     OUTPUT:
C
C     FRIGHT is altered to account for the flipping
C
C**********************************************************************
      SUBROUTINE FLIPUP8(LEVELBOX, NBOXES, ICOLBOX, IROWBOX, NLEV,
     1      IPARENTBOX, ICHILDBOX, FRIGHT, IPLUS,
     2      ISTARTLEV, IBOXLEV, NBLEVEL)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1), IPARENTBOX(1)
      INTEGER  LEVELBOX(1), IPLUS
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      REAL *8  FRIGHT(64,1)
C-----Local variables
      INTEGER  I, IBOX, II, J, NSIDE

        NSIDE = 1
        DO I = 0, NLEV
          DO II = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
           J = IBOXLEV(II)
           IBOX = J + NBOXES
           LEVELBOX(IBOX) = LEVELBOX(J)
           ICOLBOX(IBOX)  = ICOLBOX(J)

           IF (I  .EQ. 0)THEN
            IROWBOX(IBOX) = IROWBOX(J)
           ELSEIF (I  .NE. 0)THEN
            IROWBOX(IBOX) = IROWBOX(J) + 2*(NSIDE/2 - IROWBOX(J)) + 1
           ENDIF

           IROWBOX(IBOX)  = IROWBOX(IBOX) + NSIDE

           IPARENTBOX(IBOX) = IPARENTBOX(J) + NBOXES
           IF(ICHILDBOX(1,J) .EQ. -1)THEN
              ICHILDBOX(1,IBOX) = -1
              ICHILDBOX(2,IBOX) = -1
              ICHILDBOX(3,IBOX) = -1
              ICHILDBOX(4,IBOX) = -1
           ELSE
              ICHILDBOX(4,IBOX) = ICHILDBOX(1,J) + NBOXES
              ICHILDBOX(3,IBOX) = ICHILDBOX(2,J) + NBOXES
              ICHILDBOX(2,IBOX) = ICHILDBOX(3,J) + NBOXES
              ICHILDBOX(1,IBOX) = ICHILDBOX(4,J) + NBOXES
           ENDIF


           FRIGHT(57,IBOX)  = DBLE(IPLUS)*FRIGHT(1,J)
           FRIGHT(58,IBOX)  = DBLE(IPLUS)*FRIGHT(2,J)
           FRIGHT(59,IBOX)  = DBLE(IPLUS)*FRIGHT(3,J)
           FRIGHT(60,IBOX)  = DBLE(IPLUS)*FRIGHT(4,J)
           FRIGHT(61,IBOX)  = DBLE(IPLUS)*FRIGHT(5,J)
           FRIGHT(62,IBOX)  = DBLE(IPLUS)*FRIGHT(6,J)
           FRIGHT(63,IBOX)  = DBLE(IPLUS)*FRIGHT(7,J)
           FRIGHT(64,IBOX)  = DBLE(IPLUS)*FRIGHT(8,J)

           FRIGHT(49,IBOX)  = DBLE(IPLUS)*FRIGHT(9,J)
           FRIGHT(50,IBOX)  = DBLE(IPLUS)*FRIGHT(10,J)
           FRIGHT(51,IBOX)  = DBLE(IPLUS)*FRIGHT(11,J)
           FRIGHT(52,IBOX)  = DBLE(IPLUS)*FRIGHT(12,J)
           FRIGHT(53,IBOX)  = DBLE(IPLUS)*FRIGHT(13,J)
           FRIGHT(54,IBOX)  = DBLE(IPLUS)*FRIGHT(14,J)
           FRIGHT(55,IBOX)  = DBLE(IPLUS)*FRIGHT(15,J)
           FRIGHT(56,IBOX)  = DBLE(IPLUS)*FRIGHT(16,J)

           FRIGHT(41,IBOX)  = DBLE(IPLUS)*FRIGHT(17,J)
           FRIGHT(42,IBOX)  = DBLE(IPLUS)*FRIGHT(18,J)
           FRIGHT(43,IBOX)  = DBLE(IPLUS)*FRIGHT(19,J)
           FRIGHT(44,IBOX)  = DBLE(IPLUS)*FRIGHT(20,J)
           FRIGHT(45,IBOX)  = DBLE(IPLUS)*FRIGHT(21,J)
           FRIGHT(46,IBOX)  = DBLE(IPLUS)*FRIGHT(22,J)
           FRIGHT(47,IBOX)  = DBLE(IPLUS)*FRIGHT(23,J)
           FRIGHT(48,IBOX)  = DBLE(IPLUS)*FRIGHT(24,J)

           FRIGHT(33,IBOX)  = DBLE(IPLUS)*FRIGHT(25,J)
           FRIGHT(34,IBOX)  = DBLE(IPLUS)*FRIGHT(26,J)
           FRIGHT(35,IBOX)  = DBLE(IPLUS)*FRIGHT(27,J)
           FRIGHT(36,IBOX)  = DBLE(IPLUS)*FRIGHT(28,J)
           FRIGHT(37,IBOX)  = DBLE(IPLUS)*FRIGHT(29,J)
           FRIGHT(38,IBOX)  = DBLE(IPLUS)*FRIGHT(30,J)
           FRIGHT(39,IBOX)  = DBLE(IPLUS)*FRIGHT(31,J)
           FRIGHT(40,IBOX)  = DBLE(IPLUS)*FRIGHT(32,J)

           FRIGHT(25,IBOX)  = DBLE(IPLUS)*FRIGHT(33,J)
           FRIGHT(26,IBOX)  = DBLE(IPLUS)*FRIGHT(34,J)
           FRIGHT(27,IBOX)  = DBLE(IPLUS)*FRIGHT(35,J)
           FRIGHT(28,IBOX)  = DBLE(IPLUS)*FRIGHT(36,J)
           FRIGHT(29,IBOX)  = DBLE(IPLUS)*FRIGHT(37,J)
           FRIGHT(30,IBOX)  = DBLE(IPLUS)*FRIGHT(38,J)
           FRIGHT(31,IBOX)  = DBLE(IPLUS)*FRIGHT(39,J)
           FRIGHT(32,IBOX)  = DBLE(IPLUS)*FRIGHT(40,J)

           FRIGHT(17,IBOX)  = DBLE(IPLUS)*FRIGHT(41,J)
           FRIGHT(18,IBOX)  = DBLE(IPLUS)*FRIGHT(42,J)
           FRIGHT(19,IBOX)  = DBLE(IPLUS)*FRIGHT(43,J)
           FRIGHT(20,IBOX)  = DBLE(IPLUS)*FRIGHT(44,J)
           FRIGHT(21,IBOX)  = DBLE(IPLUS)*FRIGHT(45,J)
           FRIGHT(22,IBOX)  = DBLE(IPLUS)*FRIGHT(46,J)
           FRIGHT(23,IBOX)  = DBLE(IPLUS)*FRIGHT(47,J)
           FRIGHT(24,IBOX)  = DBLE(IPLUS)*FRIGHT(48,J)

           FRIGHT(9,IBOX)   = DBLE(IPLUS)*FRIGHT(49,J)
           FRIGHT(10,IBOX)  = DBLE(IPLUS)*FRIGHT(50,J)
           FRIGHT(11,IBOX)  = DBLE(IPLUS)*FRIGHT(51,J)
           FRIGHT(12,IBOX)  = DBLE(IPLUS)*FRIGHT(52,J)
           FRIGHT(13,IBOX)  = DBLE(IPLUS)*FRIGHT(53,J)
           FRIGHT(14,IBOX)  = DBLE(IPLUS)*FRIGHT(54,J)
           FRIGHT(15,IBOX)  = DBLE(IPLUS)*FRIGHT(55,J)
           FRIGHT(16,IBOX)  = DBLE(IPLUS)*FRIGHT(56,J)

           FRIGHT(1,IBOX)  = DBLE(IPLUS)*FRIGHT(57,J)
           FRIGHT(2,IBOX)  = DBLE(IPLUS)*FRIGHT(58,J)
           FRIGHT(3,IBOX)  = DBLE(IPLUS)*FRIGHT(59,J)
           FRIGHT(4,IBOX)  = DBLE(IPLUS)*FRIGHT(60,J)
           FRIGHT(5,IBOX)  = DBLE(IPLUS)*FRIGHT(61,J)
           FRIGHT(6,IBOX)  = DBLE(IPLUS)*FRIGHT(62,J)
           FRIGHT(7,IBOX)  = DBLE(IPLUS)*FRIGHT(63,J)
           FRIGHT(8,IBOX)  = DBLE(IPLUS)*FRIGHT(64,J)

         END DO
         NSIDE = 2*NSIDE
        END DO
        NBOXES = NBOXES + NBOXES

      RETURN
      END


C**************************************************************************
C     The following subroutine takes a tree structure and right hand
C     side as input and shifts it up, doubling the size of the total 
C     number of boxes.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NLEV is the finest level
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     FRIGHT represents the right hand side defined on the leaf nodes
C
C     IPLUS is a parameter that determines whether or not the sign 
C     of the right hand side is reversed or stays the
C     same upon being shifted
C
C     OUTPUT:
C
C     FRIGHT is altered to account for the flipping
C
C**************************************************************************
      SUBROUTINE SHIFTUP8(LEVELBOX, NBOXES, ICOLBOX, IROWBOX, NLEV,
     1                     IPARENTBOX, ICHILDBOX,FRIGHT,
     2                     ISTARTLEV, IBOXLEV, NBLEVEL)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1), IPARENTBOX(1)
      INTEGER  LEVELBOX(1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      REAL *8  FRIGHT(64,1)
C-----Local variables
      INTEGER  I, IBOX, II, J, NSIDE, JJ

        NSIDE = 1
        DO I = 0, NLEV
         DO JJ = ISTARTLEV(I), ISTARTLEV(I) + NBLEVEL(I) - 1
           J = IBOXLEV(JJ)
           IBOX = J + NBOXES
           LEVELBOX(IBOX) = LEVELBOX(J)
           ICOLBOX(IBOX)  = ICOLBOX(J)
           IROWBOX(IBOX)  = IROWBOX(J) + NSIDE

           IPARENTBOX(IBOX) = IPARENTBOX(J) + NBOXES
           IF(ICHILDBOX(1,J) .EQ. -1)THEN
              ICHILDBOX(1,IBOX) = -1
              ICHILDBOX(2,IBOX) = -1
              ICHILDBOX(3,IBOX) = -1
              ICHILDBOX(4,IBOX) = -1
           ELSE
              ICHILDBOX(1,IBOX) = ICHILDBOX(1,J) + NBOXES
              ICHILDBOX(2,IBOX) = ICHILDBOX(2,J) + NBOXES
              ICHILDBOX(3,IBOX) = ICHILDBOX(3,J) + NBOXES
              ICHILDBOX(4,IBOX) = ICHILDBOX(4,J) + NBOXES
           ENDIF


           DO II = 1, 64
             FRIGHT(II,IBOX) = FRIGHT(II,J)
           END DO

         END DO
        NSIDE = 2*NSIDE
        END DO
        NBOXES = NBOXES + NBOXES

      RETURN
      END


C**************************************************************************
C     The following subroutine takes a tree structure and right hand
C     side as input and reflects it over the right side, doubling
C     the size of the total number of boxes.  IPLUS determines whether
C     or not the sign of the right hand side is reversed or stays the
C     same in its reflection.
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NLEV is the finest level
C
C     IPARENTBOX denotes the parent of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     FRIGHT represents the right hand side defined on the leaf nodes
C
C     IPLUS is a parameter that determines whether or not the sign 
C     of the right hand side is reversed or stays the same in its
C     reflection
C
C     OUTPUT:
C
C     FRIGHT is altered to account for the flipping
C
C**************************************************************************
      SUBROUTINE FLIPRIGHT8(LEVELBOX, NBOXES, ICOLBOX, IROWBOX, NLEV,
     1                    IPARENTBOX, ICHILDBOX, FRIGHT, IPLUS)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  LEVELBOX(1), IPARENTBOX(1)
      INTEGER  ICHILDBOX(4,1), IPLUS
      REAL *8  FRIGHT(64,1)
C-----Local variables
      INTEGER  I, IBOX, J, NSIDE

        DO I = 0, NLEV
         NSIDE = 2**I
         DO 100 J = 1, NBOXES
           IF(LEVELBOX(J) .NE. I)GOTO 100
           IBOX = J + NBOXES
           LEVELBOX(IBOX) = LEVELBOX(J)
           IROWBOX(IBOX)  = IROWBOX(J)

           IF (I  .EQ. 0)THEN
            ICOLBOX(IBOX) = ICOLBOX(J)
           ELSEIF (I  .NE. 0)THEN
            ICOLBOX(IBOX) = ICOLBOX(J) + 2*(NSIDE/2 - ICOLBOX(J)) + 1
           ENDIF

           ICOLBOX(IBOX)  = ICOLBOX(IBOX) + NSIDE

           IPARENTBOX(IBOX) = IPARENTBOX(J) + NBOXES
           IF(ICHILDBOX(1,J) .EQ. -1)THEN
              ICHILDBOX(1,IBOX) = -1
              ICHILDBOX(2,IBOX) = -1
              ICHILDBOX(3,IBOX) = -1
              ICHILDBOX(4,IBOX) = -1
           ELSE
              ICHILDBOX(2,IBOX) = ICHILDBOX(1,J) + NBOXES
              ICHILDBOX(1,IBOX) = ICHILDBOX(2,J) + NBOXES
              ICHILDBOX(4,IBOX) = ICHILDBOX(3,J) + NBOXES
              ICHILDBOX(3,IBOX) = ICHILDBOX(4,J) + NBOXES
           ENDIF


           FRIGHT(8,IBOX)  = DBLE(IPLUS)* FRIGHT(1,J)
           FRIGHT(7,IBOX)  = DBLE(IPLUS)* FRIGHT(2,J)
           FRIGHT(6,IBOX)  = DBLE(IPLUS)* FRIGHT(3,J)
           FRIGHT(5,IBOX)  = DBLE(IPLUS)* FRIGHT(4,J)
           FRIGHT(4,IBOX)  = DBLE(IPLUS)* FRIGHT(5,J)
           FRIGHT(3,IBOX)  = DBLE(IPLUS)* FRIGHT(6,J)
           FRIGHT(2,IBOX)  = DBLE(IPLUS)* FRIGHT(7,J)
           FRIGHT(1,IBOX)  = DBLE(IPLUS)* FRIGHT(8,J)

           FRIGHT(16,IBOX)  = DBLE(IPLUS)*FRIGHT(9,J)
           FRIGHT(15,IBOX)  = DBLE(IPLUS)*FRIGHT(10,J)
           FRIGHT(14,IBOX)  = DBLE(IPLUS)*FRIGHT(11,J)
           FRIGHT(13,IBOX)  = DBLE(IPLUS)*FRIGHT(12,J)
           FRIGHT(12,IBOX)  = DBLE(IPLUS)*FRIGHT(13,J)
           FRIGHT(11,IBOX)  = DBLE(IPLUS)*FRIGHT(14,J)
           FRIGHT(10,IBOX)  = DBLE(IPLUS)*FRIGHT(15,J)
           FRIGHT(9,IBOX)   = DBLE(IPLUS)*FRIGHT(16,J)

           FRIGHT(24,IBOX)  = DBLE(IPLUS)*FRIGHT(17,J)
           FRIGHT(23,IBOX)  = DBLE(IPLUS)*FRIGHT(18,J)
           FRIGHT(22,IBOX)  = DBLE(IPLUS)*FRIGHT(19,J)
           FRIGHT(21,IBOX)  = DBLE(IPLUS)*FRIGHT(20,J)
           FRIGHT(20,IBOX)  = DBLE(IPLUS)*FRIGHT(21,J)
           FRIGHT(19,IBOX)  = DBLE(IPLUS)*FRIGHT(22,J)
           FRIGHT(18,IBOX)  = DBLE(IPLUS)*FRIGHT(23,J)
           FRIGHT(17,IBOX)  = DBLE(IPLUS)*FRIGHT(24,J)

           FRIGHT(32,IBOX)  = DBLE(IPLUS)*FRIGHT(25,J)
           FRIGHT(31,IBOX)  = DBLE(IPLUS)*FRIGHT(26,J)
           FRIGHT(30,IBOX)  = DBLE(IPLUS)*FRIGHT(27,J)
           FRIGHT(29,IBOX)  = DBLE(IPLUS)*FRIGHT(28,J)
           FRIGHT(28,IBOX)  = DBLE(IPLUS)*FRIGHT(29,J)
           FRIGHT(27,IBOX)  = DBLE(IPLUS)*FRIGHT(30,J)
           FRIGHT(26,IBOX)  = DBLE(IPLUS)*FRIGHT(31,J)
           FRIGHT(25,IBOX)  = DBLE(IPLUS)*FRIGHT(32,J)

           FRIGHT(40,IBOX)  = DBLE(IPLUS)*FRIGHT(33,J)
           FRIGHT(39,IBOX)  = DBLE(IPLUS)*FRIGHT(34,J)
           FRIGHT(38,IBOX)  = DBLE(IPLUS)*FRIGHT(35,J)
           FRIGHT(37,IBOX)  = DBLE(IPLUS)*FRIGHT(36,J)
           FRIGHT(36,IBOX)  = DBLE(IPLUS)*FRIGHT(37,J)
           FRIGHT(35,IBOX)  = DBLE(IPLUS)*FRIGHT(38,J)
           FRIGHT(34,IBOX)  = DBLE(IPLUS)*FRIGHT(39,J)
           FRIGHT(33,IBOX)  = DBLE(IPLUS)*FRIGHT(40,J)

           FRIGHT(48,IBOX)  = DBLE(IPLUS)*FRIGHT(41,J)
           FRIGHT(47,IBOX)  = DBLE(IPLUS)*FRIGHT(42,J)
           FRIGHT(46,IBOX)  = DBLE(IPLUS)*FRIGHT(43,J)
           FRIGHT(45,IBOX)  = DBLE(IPLUS)*FRIGHT(44,J)
           FRIGHT(44,IBOX)  = DBLE(IPLUS)*FRIGHT(45,J)
           FRIGHT(43,IBOX)  = DBLE(IPLUS)*FRIGHT(46,J)
           FRIGHT(42,IBOX)  = DBLE(IPLUS)*FRIGHT(47,J)
           FRIGHT(41,IBOX)  = DBLE(IPLUS)*FRIGHT(48,J)

           FRIGHT(56,IBOX)  = DBLE(IPLUS)*FRIGHT(49,J)
           FRIGHT(55,IBOX)  = DBLE(IPLUS)*FRIGHT(50,J)
           FRIGHT(54,IBOX)  = DBLE(IPLUS)*FRIGHT(51,J)
           FRIGHT(53,IBOX)  = DBLE(IPLUS)*FRIGHT(52,J)
           FRIGHT(52,IBOX)  = DBLE(IPLUS)*FRIGHT(53,J)
           FRIGHT(51,IBOX)  = DBLE(IPLUS)*FRIGHT(54,J)
           FRIGHT(50,IBOX)  = DBLE(IPLUS)*FRIGHT(55,J)
           FRIGHT(49,IBOX)  = DBLE(IPLUS)*FRIGHT(56,J)

           FRIGHT(64,IBOX)  = DBLE(IPLUS)*FRIGHT(57,J)
           FRIGHT(63,IBOX)  = DBLE(IPLUS)*FRIGHT(58,J)
           FRIGHT(62,IBOX)  = DBLE(IPLUS)*FRIGHT(59,J)
           FRIGHT(61,IBOX)  = DBLE(IPLUS)*FRIGHT(60,J)
           FRIGHT(60,IBOX)  = DBLE(IPLUS)*FRIGHT(61,J)
           FRIGHT(59,IBOX)  = DBLE(IPLUS)*FRIGHT(62,J)
           FRIGHT(58,IBOX)  = DBLE(IPLUS)*FRIGHT(63,J)
           FRIGHT(57,IBOX)  = DBLE(IPLUS)*FRIGHT(64,J)

100      CONTINUE
        END DO
        NBOXES = NBOXES + NBOXES

      RETURN
      END
 
C**************************************************************************
C     The following subroutine is set up to merge all of the boxes (after
C     the data structure has been copied over) into one new tree.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     OUTPUT:
C
C     All inputs may be altered to reflect the new tree
C
C**************************************************************************
      SUBROUTINE MERGE(LEVELBOX,NBOXES,NLEV,ICOLBOX,IROWBOX,
     1           ICHILDBOX, IPARENTBOX,NBLEVEL,IBOXLEV,ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  LEVELBOX(1), NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1), IPARENTBOX(1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
C-----Local variables
      INTEGER  IBOX, NCNTR, I, J
      INTEGER  ICOL, IROW, L, LL
      INTEGER  ICOLTEST, IROWTEST
      INTEGER  ISTART, II

C       Because the reflection routines consist of shifting the
C       parent box (at level 0) everything must be shifted down
C       by one level.
        DO IBOX = 1, NBOXES
          LEVELBOX(IBOX) = LEVELBOX(IBOX) + 1
        END DO
        NLEV = NLEV + 1

C       Now everybody has been moved down by a level, so
C       let's place a box at level zero that is the parent
C       box of the boxes that are now at level one.  By making
C       the boxes fit inside a square of side one, we can
C       just feed this into the old periodic subroutine and
C       then recover the correct function values at the end
C       by inverting the process performed in this subroutine:
        IBOX = NBOXES + 1
        LEVELBOX(IBOX) = 0
        ICOLBOX(IBOX) = 1
        IROWBOX(IBOX) = 1
        IPARENTBOX(IBOX) = -1

        NBOXES = NBOXES + 1

C       Now resort all of the boxes
        DO I = NLEV, 1, -1
          NBLEVEL(I) = NBLEVEL(I-1)
          ISTARTLEV(I) = ISTARTLEV(I-1) + 1
        END DO
        NBLEVEL(0) = 1
        ISTARTLEV(0) = 1

        NCNTR = 1
        DO I = 0, NLEV
          DO J = 1, NBOXES
            IF(LEVELBOX(J) .EQ. I)THEN
             IBOXLEV(NCNTR) = J
             NCNTR = NCNTR + 1
            ENDIF
          END DO
        END DO


         
          ISTART = ISTARTLEV(0)
          DO II = ISTART, ISTART + NBLEVEL(0) - 1
           J = IBOXLEV(II)
           ICHILDBOX(1,J) = -1
           ICHILDBOX(2,J) = -1
           ICHILDBOX(3,J) = -1
           ICHILDBOX(4,J) = -1

            ICOLTEST = 2*(ICOLBOX(J) - 1) + 1
            IROWTEST = 2*(IROWBOX(J) - 1) + 1

            DO LL = ISTARTLEV(1),ISTARTLEV(1) + NBLEVEL(1) - 1
              L = IBOXLEV(LL)

              ICOL = ICOLBOX(L)
              IROW = IROWBOX(L)
              IF(ICOL .EQ. ICOLTEST .AND. IROW .EQ. IROWTEST)THEN
                ICHILDBOX(4,J) = L
              ELSEIF(ICOL .EQ. ICOLTEST + 1
     1                        .AND. IROW .EQ. IROWTEST + 1)THEN
                ICHILDBOX(2,J) = L
              ELSEIF(ICOL .EQ. ICOLTEST
     1                        .AND. IROW .EQ. IROWTEST + 1)THEN
                ICHILDBOX(1,J) = L
              ELSEIF(ICOL .EQ. ICOLTEST + 1
     1                        .AND. IROW .EQ. IROWTEST)THEN
                ICHILDBOX(3,J) = L
             ENDIF
           END DO
          END DO


        ISTART = ISTARTLEV(1)
        DO II = ISTART, ISTART + NBLEVEL(1) - 1
          J = IBOXLEV(II)
          IPARENTBOX(J) = IBOXLEV(ISTARTLEV(0))
        END DO

      RETURN
      END


C**************************************************************************
C     The following subroutine generates the x and y values on which a 
C     function is defined, given the tree structure as input.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     NLEV is the finest level
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     OUTPUT:
C
C     XF denotes the x values on the leaf nodes
C
C     YF denotes the y values on the leaf nodes
C
C     Note the arrays XF and YF are defined on childless boxes
C     exactly the same way as the right hand side.  The 36 points
C     on each node are numbered the same way.
C
C**************************************************************************
      SUBROUTINE GETXYCLASSICAL8(XF, YF,
     1       ICOLBOX, IROWBOX, ICHILDBOX,NLEV,
     2       NBLEVEL, IBOXLEV, ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  ICOLBOX(1), IROWBOX(1), NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  XF(64,1), YF(64,1)
C-----Local variables
      INTEGER  I, IBOX, J, K, L
      REAL *8  TEMP1
      REAL *8  XSTART, PI16
      REAL *8  XX(8), XSHIFT, YSHIFT
      REAL *8  XSCALE(8)

      PI16 = DACOS(-1.0D0) / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      TEMP1 = 1.0D0
      DO K = 0, NLEV
      XSTART = (1.0D0 - TEMP1) / 2.0D0

      XSCALE(1) = XX(1) * TEMP1 - XSTART
      XSCALE(2) = XX(2) * TEMP1 - XSTART
      XSCALE(3) = XX(3) * TEMP1 - XSTART
      XSCALE(4) = XX(4) * TEMP1 - XSTART
      XSCALE(5) = XX(5) * TEMP1 - XSTART
      XSCALE(6) = XX(6) * TEMP1 - XSTART
      XSCALE(7) = XX(7) * TEMP1 - XSTART
      XSCALE(8) = XX(8) * TEMP1 - XSTART

      DO 100 I = ISTARTLEV(K), ISTARTLEV(K) + NBLEVEL(K) - 1
        IBOX = IBOXLEV(I)
        IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100

        XSHIFT  =  DBLE(ICOLBOX(IBOX) - 1) * TEMP1
        YSHIFT  =  DBLE(IROWBOX(IBOX) - 1) * TEMP1

        DO J = 1, 8
          DO L = 1, 8
            XF(8*(L-1)+J,IBOX) = XSCALE(J) + XSHIFT
            YF(8*(J-1)+L,IBOX) = XSCALE(J) + YSHIFT
          END DO
        END DO
100   CONTINUE
      TEMP1 = TEMP1/2.0D0
      END DO
      RETURN
      END


      SUBROUTINE GETXYPRACTICAL8(XFP, YFP,
     1       ICOLBOX, IROWBOX, ICHILDBOX,NLEV,
     2       NBLEVEL, IBOXLEV, ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  ICOLBOX(1), IROWBOX(1), NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  XFP(81,1), YFP(81,1)
C-----Local variables
      INTEGER  I, IBOX, J, K, L
      REAL *8  TEMP1
      REAL *8  XSTART, PI16
      REAL *8  XX(9), XSHIFT, YSHIFT
      REAL *8  XSCALE(9)

      PI16 = DACOS(-1.0D0) / 16.0D0
      XX(1) = -1.0D0 / 2.0D0
      XX(2) = DCOS(14.0D0*PI16) / 2.0D0
      XX(3) = DCOS(12.0D0*PI16) / 2.0D0
      XX(4) = DCOS(10.0D0*PI16) / 2.0D0
      XX(5) = 0.0D0
      XX(6) = DCOS( 6.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 4.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 2.0D0*PI16) / 2.0D0
      XX(9) = 1.0D0 / 2.0D0


      TEMP1 = 1.0D0
      DO K = 0, NLEV
      XSTART = (1.0D0 - TEMP1) / 2.0D0

      XSCALE(1) = XX(1) * TEMP1 - XSTART
      XSCALE(2) = XX(2) * TEMP1 - XSTART
      XSCALE(3) = XX(3) * TEMP1 - XSTART
      XSCALE(4) = XX(4) * TEMP1 - XSTART
      XSCALE(5) = XX(5) * TEMP1 - XSTART
      XSCALE(6) = XX(6) * TEMP1 - XSTART
      XSCALE(7) = XX(7) * TEMP1 - XSTART
      XSCALE(8) = XX(8) * TEMP1 - XSTART
      XSCALE(9) = XX(9) * TEMP1 - XSTART

      DO 100 I = ISTARTLEV(K), ISTARTLEV(K) + NBLEVEL(K) - 1
        IBOX = IBOXLEV(I)
        IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100

        XSHIFT  =  DBLE(ICOLBOX(IBOX) - 1) * TEMP1
        YSHIFT  =  DBLE(IROWBOX(IBOX) - 1) * TEMP1

        DO J = 1, 9
          DO L = 1, 9
            XFP(9*(L-1)+J,IBOX) = XSCALE(J) + XSHIFT
            YFP(9*(J-1)+L,IBOX) = XSCALE(J) + YSHIFT
          END DO
        END DO
100   CONTINUE
      TEMP1 = TEMP1/2.0D0
      END DO
      RETURN
      END



C********************************************************************
C     The main subroutine of multipole algorithm. Two passes are
C     executed. In the first (upward) pass, multipole expansions for
C     all boxes at all levels are computed. In the second (downward) 
C     pass, interactions are computed at successively finer levels.
C
C********************************************************************
      SUBROUTINE BOUNDFMM8(NLEV,NDEG,XNODES,WNODES,TEMP,
     1     ZS,ZSEAST,ZSWEST, ZSNORTH,ZSSOUTH,
     1     COMP,NTERMS,NNODES,POT,
     2     MPOLE,LOCEXP,EXPN,EXPS,EXPE,EXPW,
     2     EXPNBIG,EXPSBIG,EXPEBIG,EXPWBIG,
     3     C,LEVELBOX,IPARENTBOX,ICHILDBOX,
     4     ICOLBOX,IROWBOX,ICOLLEAGBOX,NBOXES,
     5     NBLEVEL,IBOXLEV,ISTARTLEV,
     6     COEFFSTOP,COEFFSSIDE,COEFFDTOP,
     7     COEFFDSIDE, FTOP,FSIDE,FTOPD,FSIDED,
     8     EDBLETOP, EDBLESIDE,  ESNGLETOP,
     9     ESNGLESIDE,  WDBLETOP, WDBLESIDE,
     1     WSNGLETOP,  WSNGLESIDE, NDBLETOP,
     2     NDBLESIDE, NSNGLETOP,  NSNGLESIDE,
     3     SDBLETOP, SDBLESIDE, SSNGLETOP,
     4     SSNGLESIDE,IFLAGEAST,IFLAGNORTH,
     5     IFLAGWEST,IFLAGSOUTH,LOCALONOFF,
     8     DOUBLETOPONOFF, DOUBLESIDEONOFF,
     9     SINGLETOPONOFF, SINGLESIDEONOFF)
      IMPLICIT NONE 
C-----Global variables 
      INTEGER  NLEV,NTERMS 
      INTEGER  NDEG
      INTEGER  NNODES, NBOXES
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1)
      INTEGER  ISTARTLEV(0:1)
      INTEGER  LEVELBOX(1), IPARENTBOX(1)
      INTEGER  ICHILDBOX(4,1), ICOLLEAGBOX(9,1)
      INTEGER  IFLAGEAST(1), IFLAGNORTH(1)
      INTEGER  IFLAGWEST(1), IFLAGSOUTH(1)
      INTEGER  LOCALONOFF(1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XNODES(1), WNODES(1)
      REAL *8  C(2*NTERMS,2*NTERMS)
      REAL *8  POT(64,1)
      REAL *8  COMP(NNODES,0:NTERMS), TEMP(1)
      COMPLEX *16  MPOLE(0:NTERMS,1),LOCEXP(0:NTERMS,1)
      COMPLEX *16  DUMMY(0:60)
      COMPLEX *16  EXPN(0:NNODES,1), EXPS(0:NNODES,1)
      COMPLEX *16  EXPE(0:NNODES,1), EXPW(0:NNODES,1)
      COMPLEX *16  EXPNBIG(0:NNODES,1), EXPSBIG(0:NNODES,1)
      COMPLEX *16  EXPEBIG(0:NNODES,1), EXPWBIG(0:NNODES,1)
      COMPLEX *16  ZS(-3:3,-3:3,0:1)
      COMPLEX *16  ZSEAST(-3:3,-3:3,0:1)
      COMPLEX *16  ZSWEST(-3:3,-3:3,0:1)
      COMPLEX *16  ZSNORTH(-3:3,-3:3,0:1)
      COMPLEX *16  ZSSOUTH(-3:3,-3:3,0:1)
      COMPLEX *16  EDBLETOP(0:NNODES,0:7)
      COMPLEX *16  EDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  ESNGLETOP(0:NNODES,0:7)
      COMPLEX *16  ESNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  WDBLETOP(0:NNODES,0:7)
      COMPLEX *16  WDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  WSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  WSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  NDBLETOP(0:NNODES,0:7)
      COMPLEX *16  NDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  NSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  NSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  SDBLETOP(0:NNODES,0:7)
      COMPLEX *16  SDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  SSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  SSNGLESIDE(0:NNODES,0:7)
C-----Local variables
      INTEGER  IPERIOD
      INTEGER  I, J, IP
      INTEGER  IC1,IC2,IC3,IC4
      INTEGER  INALL(6),IYNALL(6),IN12(4),IY12(4)
      INTEGER  ISALL(6),IYSALL(6),IS34(4),IY34(4)
      INTEGER  NNALL,NN12,NSALL,NS34
      INTEGER  IEALL(4),IYEALL(4),IE13(2),IY13(2)
      INTEGER  IE1(1),IY1(1),IE3(1),IY3(1)
      INTEGER  IWALL(4),IYWALL(4),IW24(2),IY24(2)
      INTEGER  IW2(1),IY2(1),IW4(1),IY4(1)
      INTEGER  NEALL,NE13,NE1,NE3,NWALL,NW24,NW2,NW4
      INTEGER  NSIDE, NSIDEMARK
      INTEGER  INBIG12(3), ISBIG34(3)
      INTEGER  IEBIG13(1), IWBIG24(1)
      INTEGER  IEBIG1(1), IWBIG2(1)
      INTEGER  IEBIG3(1), IWBIG4(1)
      INTEGER  NB, ISTART, IEND
      INTEGER  IOUT, II, JJ
      INTEGER  IFAR1, IFAR2, IFAR3
      INTEGER  ICLOSE1, ICLOSE2, IER
      REAL *8  COEFFDTOP(0:7,1)
      REAL *8  COEFFDSIDE(0:7,1)
      REAL *8  COEFFSTOP(0:7,1)
      REAL *8  COEFFSSIDE(0:7,1)
      REAL *8  XP(64),YP(64)
      REAL *8  XLENGTH
      REAL *8  T(10), SCALE(0:100), TEMP1, TEMP2
      REAL *8  SUM, ZERO
      REAL *8  WDSIDE(9,64,8), WDTOP(9,64,8)
      REAL *8  WSSIDE(9,64,8), WSTOP(9,64,8)
      REAL *8  WDTOPBTOS(12,64,8), WDSIDEBTOS(12,64,8)
      REAL *8  WSTOPBTOS(12,64,8), WSSIDEBTOS(12,64,8)
      REAL *8  WDTOPSTOB(12,64,8), WDSIDESTOB(12,64,8)
      REAL *8  WSTOPSTOB(12,64,8), WSSIDESTOB(12,64,8)
      REAL *8  X
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  PI
      COMPLEX *16  FTARGET1, FTARGET2, FTARGET3, IMAG
      COMPLEX *16  ZSHIFT, ZPOT
      COMPLEX *16  B(0:50)
      COMPLEX *16  MEXNALL(0:50)
      COMPLEX *16  MEXN12(0:50)
      COMPLEX *16  MEXSALL(0:50)
      COMPLEX *16  MEXS34(0:50)
      COMPLEX *16  MEXEALL(0:50)
      COMPLEX *16  MEXE13(0:50)
      COMPLEX *16  MEXE1(0:50)
      COMPLEX *16  MEXE3(0:50)
      COMPLEX *16  MEXWALL(0:50)
      COMPLEX *16  MEXW24(0:50)
      COMPLEX *16  MEXW2(0:50)
      COMPLEX *16  MEXW4(0:50)
      COMPLEX *16  SPIN
      COMPLEX *16  EXPNALL(0:60), EXPSALL(0:60)
      COMPLEX *16  EXPEALL(0:60), EXPWALL(0:60)
      COMPLEX *16  TEMPEAST(64,40), TEMPNORTH(64,40)
      COMPLEX *16  TEMPWEST(64,40), TEMPSOUTH(64,40)
      COMPLEX *16  TEMPSHIFTWEST, TEMPSHIFTNORTH
      COMPLEX *16  TEMPSHIFTEAST, TEMPSHIFTSOUTH
      DATA IMAG/(0.0D0,1.0D0)/
      DATA ZERO/0.0D0/
      PI = DACOS(-1.0D0)
      IPERIOD = 2

C     Let's initially set the potential to zero:
      DO I = 1, NBOXES
       DO J = 0, 7
         COEFFDTOP(J,I)  = ZERO
         COEFFSTOP(J,I)  = ZERO
         COEFFDSIDE(J,I) = ZERO
         COEFFSSIDE(J,I) = ZERO
       END DO
      END DO

      DO I = 1, NBOXES
        IFLAGEAST(I) = 0
        IFLAGWEST(I) = 0
        IFLAGNORTH(I) = 0
        IFLAGSOUTH(I) = 0
        LOCALONOFF(I) = 0
      END DO

C     Initialize the LOCALONOFF switch to the correct values
      NSIDE = 1
      DO JJ = 0, NLEV
        ISTART = ISTARTLEV(JJ)
        IEND = ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
        IF(NSIDE .EQ. 1)THEN
          NSIDEMARK = 1
        ELSE
          NSIDEMARK = NSIDE / 2
        ENDIF
        DO II = ISTART, IEND
          I = IBOXLEV(II)
          IF(IROWBOX(I) .LE. NSIDEMARK .AND.
     1       ICOLBOX(I) .LE. NSIDEMARK)THEN
            LOCALONOFF(I) = 1
          ENDIF
        END DO
       NSIDE = 2 * NSIDE
      END DO


C
C     create SCALE array
C
      X = 0.50D0
      DO I = 0,NLEV
       SCALE(I) = X
       X = X*2.0D0
      ENDDO


C     Call a set of routines that will precompute various 
C     things needed later in the adapfmm, mkcoeffs, and
C     mkshifts routines.
      CALL PWTS4(XNODES, WNODES, NNODES)

      CALL PRECOMPUTE(COMP,NTERMS,NNODES,
     1  XNODES,WNODES, TEMP,C, SUM)

      DO I = 1, 8
        XP((I-1)*8+1) = DCOS(15D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+2) = DCOS(13D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+3) = DCOS(11D0*PI/16.0D0) / 2.0D0
        XP((I-1)*8+4) = DCOS(9D0*PI/16.0D0)  / 2.0D0
        XP((I-1)*8+5) = DCOS(7D0*PI/16.0D0)  / 2.0D0
        XP((I-1)*8+6) = DCOS(5D0*PI/16.0D0)  / 2.0D0
        XP((I-1)*8+7) = DCOS(3D0*PI/16.0D0)  / 2.0D0
        XP((I-1)*8+8) = DCOS(1D0*PI/16.0D0)  / 2.0D0

        YP(I)    = DCOS(15D0*PI/16.0D0) / 2.0D0
        YP(8+I)  = DCOS(13D0*PI/16.0D0) / 2.0D0
        YP(16+I) = DCOS(11D0*PI/16.0D0) / 2.0D0
        YP(24+I) = DCOS(9D0*PI/16.0D0)  / 2.0D0
        YP(32+I) = DCOS(7D0*PI/16.0D0)  / 2.0D0
        YP(40+I) = DCOS(5D0*PI/16.0D0)  / 2.0D0
        YP(48+I) = DCOS(3D0*PI/16.0D0)  / 2.0D0
        YP(56+I) = DCOS(1D0*PI/16.0D0)  / 2.0D0
      ENDDO

      DO JJ = 1, 64
        TEMPSHIFTWEST = 2.0D0*(XP(JJ) + IMAG*YP(JJ))
        TEMPSHIFTEAST = -TEMPSHIFTWEST
        TEMPSHIFTNORTH = IMAG*TEMPSHIFTWEST
        TEMPSHIFTSOUTH = -TEMPSHIFTNORTH
        DO II = 1, NNODES
          TEMPEAST(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTEAST)
          TEMPWEST(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTWEST)
          TEMPNORTH(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTNORTH)
          TEMPSOUTH(JJ,II) = -CDEXP(XNODES(II) * TEMPSHIFTSOUTH)
        END DO
      END DO


      DO I = 0, NLEV
       ISTART = ISTARTLEV(I)
       IEND = ISTART + NBLEVEL(I) - 1
       DO II = ISTART,IEND
        J = IBOXLEV(II)
C       We only need to concern ourselves 
C       with childless boxes.
        IF(ICHILDBOX(1,J) .LT. 0)THEN

           IF(SINGLETOPONOFF(J) .EQ. 1)THEN
             CALL INTERPOLATE8SIDE(FTOP(1,J),  COEFFSTOP(0,J))
           ENDIF

           IF(DOUBLETOPONOFF(J) .EQ. 1)THEN
             CALL INTERPOLATE8SIDE(FTOPD(1,J), COEFFDTOP(0,J))
           ENDIF

           IF(SINGLESIDEONOFF(J) .EQ. 1)THEN
             CALL INTERPOLATE8SIDE(FSIDE(1,J), COEFFSSIDE(0,J))
           ENDIF

           IF(DOUBLESIDEONOFF(J) .EQ. 1)THEN
             CALL INTERPOLATE8SIDE(FSIDED(1,J),COEFFDSIDE(0,J))
           ENDIF

        ENDIF
       END DO
      END DO


      CALL DLAYSIDE(WDSIDE)
      CALL DLAYTOP(WDTOP)
      CALL SLAYSIDE(WSSIDE)
      CALL SLAYTOP(WSTOP)

      CALL DLAYTOPBTOS(WDTOPBTOS)
      CALL DLAYSIDEBTOS(WDSIDEBTOS)
      CALL SLAYTOPBTOS(WSTOPBTOS)
      CALL SLAYSIDEBTOS(WSSIDEBTOS)

      CALL DLAYTOPSTOB(WDTOPSTOB)
      CALL DLAYSIDESTOB(WDSIDESTOB)
      CALL SLAYTOPSTOB(WSTOPSTOB)
      CALL SLAYSIDESTOB(WSSIDESTOB)


C     Call a routine that will generate the array ZS
C     which will be needed later on when shifting the
C     exponential expansions.
      CALL MKSHIFTS2D(XNODES,NNODES,ZS,
     1    ZSEAST,ZSWEST,ZSNORTH,ZSSOUTH)

      CALL MKCOLLS(ICOLBOX,
     1    IROWBOX, ICOLLEAGBOX, NBOXES, NLEV,
     2    IPARENTBOX, ICHILDBOX, NBLEVEL,
     3    IBOXLEV, ISTARTLEV, IPERIOD)

 
C**********************************************************************
C     =============================================
C     UPWARD PASS
C     =============================================
C
C-----initialize multipole and local expansions to zero.
C
      DO I = 1, NBOXES
       DO J = 0, NTERMS
         MPOLE(J,I)  = ZERO 
         LOCEXP(J,I) = ZERO 
       END DO
       DO J = 0, NNODES
         EXPE(J,I) = ZERO
         EXPW(J,I) = ZERO
         EXPN(J,I) = ZERO
         EXPS(J,I) = ZERO
         EXPEBIG(J,I) = ZERO
         EXPWBIG(J,I) = ZERO
         EXPNBIG(J,I) = ZERO
         EXPSBIG(J,I) = ZERO
       END DO
      END DO


      DO I = NLEV, 0, -1
         XLENGTH = 1.0D0/SCALE(I)
         ISTART = ISTARTLEV(I)
         IEND = ISTART + NBLEVEL(I) - 1
         DO II = ISTART,IEND
            J = IBOXLEV(II)
            IF(ICHILDBOX(1,J) .LT. 0)THEN

              IF(SINGLETOPONOFF(J) .EQ. 1)THEN
                CALL SINGLETOPTOMULT8(COEFFSTOP(0,J),
     1                        NTERMS,C,XLENGTH,DUMMY)
    
                DO JJ = 0, NTERMS
                  MPOLE(JJ,J) = MPOLE(JJ,J) + DUMMY(JJ)
                END DO
              ENDIF

              IF(SINGLESIDEONOFF(J) .EQ. 1)THEN
                CALL SINGLESIDETOMULT8(COEFFSSIDE(0,J),
     1                        NTERMS,C,XLENGTH, DUMMY)

                DO JJ = 0, NTERMS
                  MPOLE(JJ,J) = MPOLE(JJ,J) + DUMMY(JJ)
                END DO
              ENDIF


              IF(DOUBLETOPONOFF(J) .EQ. 1)THEN
                CALL DOUBLETOPTOMULT8(COEFFDTOP(0,J),
     1                        NTERMS,C,DUMMY)
   
                DO JJ = 0, NTERMS
                  MPOLE(JJ,J) = MPOLE(JJ,J) + DUMMY(JJ)
                END DO
              ENDIF

              IF(DOUBLESIDEONOFF(J) .EQ. 1)THEN
                CALL DOUBLESIDETOMULT8(COEFFDSIDE(0,J),
     1                        NTERMS,C,DUMMY)

                DO JJ = 0, NTERMS
                  MPOLE(JJ,J) = MPOLE(JJ,J) + DUMMY(JJ)
                END DO
              ENDIF


            ELSEIF(ICHILDBOX(1,J) .GT. 0)THEN

              IC1 = ICHILDBOX(1,J)
              IC2 = ICHILDBOX(2,J)
              IC3 = ICHILDBOX(3,J)
              IC4 = ICHILDBOX(4,J)

            CALL CHILDPAR(MPOLE(0,J),
     1        MPOLE(0,IC1),MPOLE(0,IC2),
     2        MPOLE(0,IC3),MPOLE(0,IC4),
     3        NTERMS, C)
            ENDIF
         END DO
      END DO


C     Let's call a routine that will set the initial local 
C     expansion to the correct value to account for the 
C     periodicity in the far field terms:

C      Set IFLAG2 to zero to indicate that we are
C      working with the potential and not the
C      electric field.
        ISTART = ISTARTLEV(0)
        J = IBOXLEV(ISTART)
        CALL PSPPIN(IER)
        CALL PSPPTA(MPOLE(0,J),LOCEXP(0,J), NTERMS)


C     =============================================
C     DOWNWARD PASS
C     =============================================
      XLENGTH = 1.0D0
      DO I = 0,NLEV
C        First let's set up an array of scaling scale
C        factors needed in the local interaction:
          T(1)  =  2.0D0**I 
          T(2)  =  2.0D0 * T(1) 
          TEMP1 =  -DLOG(T(1))
          TEMP2 =  -DLOG(T(2))
          T(3)  =  TEMP1 * 2.0D0
          T(4)  = -TEMP1 * 2.0D0/3.0D0
          T(5)  = -TEMP1 * 2.0D0/15.0D0
          T(6)  = -TEMP1 * 2.0D0/35.0D0
          T(7)  =  TEMP2 * 2.0D0
          T(8)  = -TEMP2 * 2.0D0/3.0D0
          T(9)  = -TEMP2 * 2.0D0/15.0D0
          T(10) = -TEMP2 * 2.0D0/35.0D0


         ISTART = ISTARTLEV(I)
         IEND = ISTART + NBLEVEL(I) - 1
         DO II = ISTART,IEND

           J = IBOXLEV(II)
           IF(ICHILDBOX(1,J) .GT. 0)THEN
C          If the box has children, do all of the work involving local 
C          expansions here.

            IC1 = ICHILDBOX(1,J)
            IC2 = ICHILDBOX(2,J)
            IC3 = ICHILDBOX(3,J)
            IC4 = ICHILDBOX(4,J)

            CALL PARENTCHILD(LOCEXP(0,J),LOCEXP(0,IC1),
     1       LOCEXP(0,IC2),LOCEXP(0,IC3),
     2       LOCEXP(0,IC4),NTERMS,C,LOCALONOFF(J))

            CALL MKEXP2D(J,NTERMS,MPOLE,
     1       NNODES,MEXNALL,MEXN12,MEXSALL,MEXS34,MEXEALL,
     2       MEXE13,MEXE1,MEXE3,MEXWALL,MEXW24,MEXW2,MEXW4,
     3       ZS,COMP,WNODES,SCALE(I+1),SUM,TEMP,
     4       ICHILDBOX)
 

            CALL MKLISTS(J,INALL,NNALL,IYNALL,
     1       IN12,NN12,IY12,
     2       ISALL,NSALL,IYSALL,IS34,NS34,IY34,
     3       IEALL,NEALL,IYEALL,IE13,NE13,IY13,
     4       IWALL,NWALL,IYWALL,IW24,NW24,IY24,
     5       IW2,IY2,NW2,IW4,IY4,NW4,
     6       IE1,IY1,NE1,IE3,IY3,NE3,
     7       INBIG12,ISBIG34,IEBIG13,IWBIG24,
     8       IEBIG1, IWBIG2, IEBIG3, IWBIG4,
     9       ICOLLEAGBOX,ICHILDBOX,
     1       ICOLBOX, IROWBOX, IPERIOD,
     2       IFLAGEAST, IFLAGWEST, IFLAGNORTH,
     3       IFLAGSOUTH, LOCALONOFF)


            CALL PROCESSNO(EXPN,INALL,NNALL,IYNALL,
     1       IN12,NN12,IY12,MEXNALL,MEXN12,ZS,NNODES,
     2       INBIG12,EXPNBIG,ZSNORTH,LOCALONOFF)
            CALL PROCESSSO(EXPS,ISALL,NSALL,IYSALL,
     1       IS34,NS34,IY34,MEXSALL,MEXS34,ZS,NNODES,
     2       ISBIG34,EXPSBIG,ZSSOUTH, LOCALONOFF)
            CALL PROCESSEA(EXPE,IEALL,NEALL,IYEALL,
     1       IE13,NE13,IY13,IE1,NE1,IY1,IE3,NE3,IY3,
     2       MEXEALL,MEXE13,MEXE1,MEXE3,ZS,NNODES,
     3       IEBIG13,EXPEBIG,IEBIG1,IEBIG3,ZSEAST,
     4       LOCALONOFF)
            CALL PROCESSWE(EXPW,IWALL,NWALL,IYWALL,
     1       IW24,NW24,IY24,IW2,NW2,IY2,IW4,NW4,IY4,
     2       MEXWALL,MEXW24,MEXW2,MEXW4,ZS,NNODES,
     3       IWBIG24,EXPWBIG,IWBIG2,IWBIG4,ZSWEST,
     4       LOCALONOFF)


           ELSEIF (ICHILDBOX(1,J) .LT. 0)THEN
C           NOW LET'S SCAN THE COLLEAGUES
 
            DO 250 NB = 1, 9
              IOUT = ICOLLEAGBOX(NB,J)
              IF(IOUT .LT. 0)GOTO 250
              IF(ICHILDBOX(1,IOUT) .LT. 0)THEN


           CALL COLLOCLAY8(POT(1,IOUT),COEFFSTOP(0,J),
     1          COEFFDTOP(0,J),COEFFSSIDE(0,J),COEFFDSIDE(0,J),
     2          NDEG,NB,T,WDTOP,WDSIDE,WSSIDE,WSTOP,LOCALONOFF(IOUT),
     3          DOUBLESIDEONOFF(J),DOUBLETOPONOFF(J),
     4          SINGLESIDEONOFF(J),SINGLETOPONOFF(J))

              ELSEIF(ICHILDBOX(1,IOUT) .GT. 0)THEN
C               Colleague has children (have to go to big to small and
C               small to big stuff)


                IC1 = ICHILDBOX(2,IOUT)
                IC2 = ICHILDBOX(3,IOUT)
                IC3 = ICHILDBOX(1,IOUT)
                IC4 = ICHILDBOX(4,IOUT)

C              Form the four expansions needed in the
C              big to small and small to big process.

              CALL MKEXPBTOSLAY8(COEFFDTOP(0,J),COEFFSTOP(0,J),
     1         COEFFDSIDE(0,J), COEFFSSIDE(0,J),XLENGTH,
     2         SUM, EXPNALL,EXPSALL,EXPEALL,EXPWALL,
     3         EDBLETOP, WDBLETOP, NDBLETOP, SDBLETOP,
     4         EDBLESIDE, WDBLESIDE, NDBLESIDE, SDBLESIDE,
     5         ESNGLETOP, WSNGLETOP, SSNGLETOP, NSNGLETOP,
     6         ESNGLESIDE, WSNGLESIDE, SSNGLESIDE, NSNGLESIDE,NNODES,
     7         DOUBLESIDEONOFF(J),DOUBLETOPONOFF(J),
     8         SINGLESIDEONOFF(J),SINGLETOPONOFF(J))


                IF(NB .EQ. 1)THEN
C                 Colleague with small boxes is in the lower left corner,
C                 one box is not well separated and 3 are.
 
                  IFAR1 = IC4
                  IFAR2 = IC2
                  IFAR3 = IC3
                  ICLOSE1 = IC1

C                 First do the local work, small to big

                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,1,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))
 
C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,1,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))


C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (-2.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))

                  FTARGET2 = (-1.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))


 
                  SPIN = 1.0D0
                  FTARGET3 = (-2.0D0,-1.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR3))

 

                ELSEIF(NB .EQ. 2)THEN
C                 Colleague with small boxes is below this box
C                 two boxes are not well separated and two are.

                  IFAR1 = IC4
                  IFAR2 = IC2
                  ICLOSE1 = IC3
                  ICLOSE2 = IC1

C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,2,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE2),
     1               COEFFSSIDE(0,ICLOSE2),COEFFDTOP(0,ICLOSE2),
     2               COEFFDSIDE(0,ICLOSE2),NDEG,3,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

 

C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,2,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))


                  CALL BTOSLOCLAY8(POT(1,ICLOSE2),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,3,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE2))

 
C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (0.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))


                  FTARGET2 = (1.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))

 
                ELSEIF(NB .EQ. 3)THEN
C                 Colleague with small boxes is in the lower right corner,
C                 one box is not well separated and 3 are.

                  ICLOSE1 = IC3
                  IFAR1 = IC4  
                  IFAR2 = IC2 
                  IFAR3 = IC1

C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,4,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))


C                 First do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,4,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))
            

C                 Finally do the far work, big to small
                  SPIN = -IMAG
                  FTARGET1 = (2.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR1), 
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,-2.0D0)
                  CALL BTOSFAR(EXPS(0,IFAR2), 
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPSALL,SPIN,LOCALONOFF(IFAR2))

 
                  SPIN = -1.0D0
                  FTARGET3 = (3.0D0,-1.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR3))

            
                ELSEIF(NB .EQ. 4)THEN
C                 Colleague with small boxes is left of this box
C                 two boxes are not well separated and two are.
 
 
                  ICLOSE1 = IC2
                  ICLOSE2 = IC1
                  IFAR1 = IC4 
                  IFAR2 = IC3
 
C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,5,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE2),
     1               COEFFSSIDE(0,ICLOSE2),COEFFDTOP(0,ICLOSE2),
     2               COEFFDSIDE(0,ICLOSE2),NDEG,7,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

       

C                 Next do the local work, big to small 
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,5,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))

                  CALL BTOSLOCLAY8(POT(1,ICLOSE2),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,7,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE2))

 
C                 Finally do the far work, big to small
                  SPIN = 1.0D0
                  FTARGET1 = (-2.0D0,0.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR1))

                  FTARGET2 = (-2.0D0,1.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 6)THEN
C                 Colleague with small boxes is right this box
C                 two boxes are not well separated and two are.
            
                  IFAR1 = IC1
                  IFAR2 = IC2
                  ICLOSE1 = IC3
                  ICLOSE2 = IC4


C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,8,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE2),
     1               COEFFSSIDE(0,ICLOSE2),COEFFDTOP(0,ICLOSE2),
     2               COEFFDSIDE(0,ICLOSE2),NDEG,6,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

 

C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,8,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))

                  CALL BTOSLOCLAY8(POT(1,ICLOSE2),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,6,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE2))

 

C                 Finally do the far work, big to small
                  SPIN = -1.0D0
                  FTARGET1 = (3.0D0,1.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,0.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 7)THEN
C                 Colleague with small boxes is in the upper left corner,
C                 one box is not well separated and 3 are.
 
                  ICLOSE1 = IC2
                  IFAR1 = IC4
                  IFAR2 = IC3
                  IFAR3 = IC1

C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,9,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))



C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,9,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))

 

C                 Finally do the far work, big to small
                  SPIN = 1.0D0
                  FTARGET1 = (-2.0D0,2.0D0)
                  CALL BTOSFAR(EXPW(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPWALL,SPIN,LOCALONOFF(IFAR1))

 
                  SPIN = IMAG
                  FTARGET2 = (-2.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))

 
                  FTARGET3 = (-1.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR3))


                ELSEIF(NB .EQ. 8)THEN
C                 Colleague with small boxes is above this box
C                 two boxes are not well separated and two are.

                  ICLOSE1 = IC4
                  ICLOSE2 = IC2
                  IFAR1 = IC3
                  IFAR2 = IC1

C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,10,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))

                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE2),
     1               COEFFSSIDE(0,ICLOSE2),COEFFDTOP(0,ICLOSE2),
     2               COEFFDSIDE(0,ICLOSE2),NDEG,11,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))



C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,10,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))

                  CALL BTOSLOCLAY8(POT(1,ICLOSE2),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,11,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE2))



C                 Finally do the far work, big to small
                  SPIN = IMAG
                  FTARGET1 = (0.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (1.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2), 
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))


                ELSEIF(NB .EQ. 9)THEN
C                 Colleague with small boxes is in the upper right corner,
C                 one box is not well separated and 3 are.

                  ICLOSE1 = IC4
                  IFAR1 = IC3  
                  IFAR2 = IC1
                  IFAR3 = IC2

C                 First do the local work, small to big
                  CALL STOBLOCLAY8(POT(1,J),COEFFSTOP(0,ICLOSE1),
     1               COEFFSSIDE(0,ICLOSE1),COEFFDTOP(0,ICLOSE1),
     2               COEFFDSIDE(0,ICLOSE1),NDEG,12,T,WDTOPSTOB,
     3               WDSIDESTOB,WSTOPSTOB,WSSIDESTOB,LOCALONOFF(J))



C                 Next do the local work, big to small
                  CALL BTOSLOCLAY8(POT(1,ICLOSE1),
     1                COEFFSTOP(0,J), COEFFSSIDE(0,J),
     2                COEFFDTOP(0,J), COEFFDSIDE(0,J),
     3                NDEG,12,T,WDTOPBTOS,WDSIDEBTOS,
     4                WSTOPBTOS,WSSIDEBTOS,LOCALONOFF(ICLOSE1))

 

C                 Finally do the far work, big to small
                  SPIN = IMAG
                  FTARGET1 = (2.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR1),
     1             FTARGET1,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR1))

 
                  FTARGET2 = (3.0D0,3.0D0)
                  CALL BTOSFAR(EXPN(0,IFAR2),
     1             FTARGET2,NTERMS,
     2             XNODES,COMP,
     3             EXPNALL,SPIN,LOCALONOFF(IFAR2))
 
                  SPIN = -1.0D0
                  FTARGET3 = (3.0D0,2.0D0)
                  CALL BTOSFAR(EXPE(0,IFAR3),
     1             FTARGET3,NTERMS,
     2             XNODES,COMP,
     3             EXPEALL,SPIN,LOCALONOFF(IFAR3))
                ENDIF
              ENDIF
250         CONTINUE
           ENDIF
         END DO

C        For boxes with children, convert the exponential expansions
C        to one local expansion:
         DO JJ = ISTART,IEND
          J = IBOXLEV(JJ)
          IF(ICHILDBOX(1,J) .GT. 0 .AND. LOCALONOFF(J) .EQ. 1)THEN
              DO II = 1, 4
                IC1 = ICHILDBOX(II,J)
                  CALL EXP4LOCAL(B,NTERMS,NNODES,EXPW(0,IC1),
     1               EXPS(0,IC1), EXPE(0,IC1), EXPN(0,IC1),COMP)
                  CALL ADDEXP(B,LOCEXP(0,IC1),NTERMS)
              END DO


          ELSEIF(ICHILDBOX(1,J) .LT. 0 .AND. LOCALONOFF(J) .EQ. 1)THEN
C          Evaluate the exponential expansions at the target points

               CALL EXPBIG4EVAL(POT(1,J),NNODES,EXPWBIG(0,J),
     1           EXPSBIG(0,J), EXPEBIG(0,J), EXPNBIG(0,J),
     2           IFLAGEAST(J),IFLAGWEST(J),
     3           IFLAGNORTH(J),IFLAGSOUTH(J),TEMPEAST,
     4           TEMPNORTH, TEMPWEST, TEMPSOUTH)
          ENDIF
         END DO
      XLENGTH = XLENGTH / 2.0D0
      ENDDO

      DO JJ = 0, NLEV
       ISTART = ISTARTLEV(JJ)
       IEND = ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
        DO II = ISTART, IEND
          I = IBOXLEV(II)
          IF(ICHILDBOX(1,I) .LT. 0 .AND. LOCALONOFF(I) .EQ. 1)THEN
           DO IP = 1,64
             ZSHIFT = DCMPLX(XP(IP),YP(IP))
             ZPOT = ZERO
              DO J = 1,NTERMS
                ZPOT = ZSHIFT*ZPOT + LOCEXP(NTERMS-J,I)
              END DO
             POT(IP,I) = POT(IP,I) + DREAL(ZPOT)
           END DO
          ENDIF
        END DO
      END DO

      RETURN
      END


C**************************************************************************
C     The following subroutine maps from four function values on the side
C     of the box (representing charge densities) to four polynomial 
C     coefficients approximating the polynomial.
C
C     INPUT:
C
C     F is the array of 6 function values arranged from left to right
C     
C     XLENGTH is the length of a side of the box being considered
C
C     B is a precomputed matrix that maps from the 6 function values to
C       the 6 coefficients values
C
C     OUTPUT:
C
C     COEFFS contains the coefficients approximating the polynomial
C
C**************************************************************************
      SUBROUTINE INTERPOLATE8SIDE(F ,COEFF)
      IMPLICIT NONE
C-----Global variables
      REAL *8  F(8), COEFF(0:7)
C-----Local variables
      INTEGER  J
      REAL *8  FTEMP(8,1)
      REAL *8  COEFFTEMP(8,1)
      REAL *8  WSAVE2(1000)

C       initialize WSAVE2 array (needed by GETCOEFF).
        CALL CHXCIN(8,WSAVE2)

        DO J = 1, 8
            FTEMP(J,1) = F(J)
        END DO

C       compute Chebyshev transforms
        CALL GETCOEFF1(FTEMP,COEFFTEMP,WSAVE2)

C       Now set these values to out coefficients
        DO J = 0, 7
            COEFF(J) = COEFFTEMP(J+1,1)
          END DO
      RETURN
      END



C**************************************************************************
C     The following subroutine maps from four polynomial coefficients to
C     the multipole expansion caused by a single layer potential on the
C     top of the box.  The multipole expansion is shifted within this 
C     routine so that it lies at the center of the box.
C
C     INPUT:
C
C     SCOEFFS is are the coefficients approximating the charge density
C     
C     NTERMS is the order of the multipole expansion
C
C     XLENGTH is the length of a side of the box being considered
C
C     C is an array of binomial coefficients
C
C     OUTPUT:
C
C     MPOLE contains the approximating multipole expansion
C
C**************************************************************************
      SUBROUTINE SINGLETOPTOMULT8(SCOEFFS,NTERMS,C,XLENGTH,MPOLE) 
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS
      REAL *8  SCOEFFS(0:7)
      REAL *8  C(2*NTERMS,2*NTERMS)
      REAL *8  XLENGTH
      COMPLEX *16  MPOLE(0:NTERMS)
C-----Local variables
      INTEGER  I, K, M
      REAL *8  SCOEFFSCALED(0:7)
      COMPLEX *16  Z0, Z00, Z0P(0:60), Z0P2(0:60)
      COMPLEX *16  CD, CDD
      COMPLEX *16  ANEW(0:60), IMAG
      COMPLEX *16  SUM
      DATA IMAG/(0.0,1.0)/

      SCOEFFSCALED(0) = SCOEFFS(0) * XLENGTH
      SCOEFFSCALED(1) = SCOEFFS(1) * XLENGTH
      SCOEFFSCALED(2) = SCOEFFS(2) * XLENGTH
      SCOEFFSCALED(3) = SCOEFFS(3) * XLENGTH
      SCOEFFSCALED(4) = SCOEFFS(4) * XLENGTH
      SCOEFFSCALED(5) = SCOEFFS(5) * XLENGTH
      SCOEFFSCALED(6) = SCOEFFS(6) * XLENGTH
      SCOEFFSCALED(7) = SCOEFFS(7) * XLENGTH

      DO I = 0, NTERMS
         MPOLE(I) = 0.0D0
      END DO

      Z0 = IMAG * .5D0
      Z0P(0) = 1.0D0
      Z00 = Z0
      CD = 1.0D0 / Z0
      CDD = CD
      Z0P2(0) = 1.0D0
      DO I = 1, NTERMS
         Z0P(I)  = Z00
         Z0P2(I) = CDD
         CDD     = CDD*CD
         Z00     = Z00*Z0
      END DO

C     First, get the zero term in the multipole expansion
C     (Note that even though we are summing over all of the
C     coefficients, that we only need be concerned with 0 and
C     2, because the integrals in the other two cases are zero.)
      MPOLE(0) = SCOEFFSCALED(0) 
     1         + SCOEFFSCALED(2) * (2.0D0/3.0D0 - 1.0D0)
     2         + SCOEFFSCALED(4) * (8.0D0/5.0D0 - 8.0D0/3.0D0 + 1.0D0)
     3         + SCOEFFSCALED(6) * (32.0D0/7.0D0 - 48.0D0/5.0D0 
     4                           +  18.0D0/3.0D0 - 1.0D0)


      DO K = 1, NTERMS
        SUM = 0.0D0
C         We only need to concern ourselves with the cases
C         where i+k is even:
          IF(MOD(K,2) .EQ. 0)THEN

           SUM = SUM + SCOEFFSCALED(0) / DBLE(K+1)

           SUM = SUM + SCOEFFSCALED(2) * 
     1         (2.0D0/DBLE(K+3) - 1.0D0/DBLE(K+1))

           SUM = SUM + SCOEFFSCALED(4) *
     1       (8.0D0/DBLE(K+5) - 8.0D0/DBLE(K+3) + 1.0D0/DBLE(K+1))

           SUM = SUM + SCOEFFSCALED(6) *
     1       (32.0D0/DBLE(K+7) - 48.0D0/DBLE(K+5) 
     2      + 18.0D0/DBLE(K+3) -  1.0D0/DBLE(K+1))

          ELSE

           SUM = SUM + SCOEFFSCALED(1) / DBLE(K+2)

           SUM = SUM + SCOEFFSCALED(3) *
     1         (4.0D0/DBLE(K+4) - 3.0D0/DBLE(K+2))

           SUM = SUM + SCOEFFSCALED(5) *
     1      (16.0D0/DBLE(K+6) - 20.0D0/DBLE(K+4) + 5.0D0/DBLE(K+2))

           SUM = SUM + SCOEFFSCALED(7) *
     1      (64.0D0/DBLE(K+8) - 112.0D0/DBLE(K+6) 
     2     + 56.0D0/DBLE(K+4) -   7.0D0/DBLE(K+2))

          ENDIF
        ANEW(K) = -SUM * Z0P2(K) / (DBLE(K)*2.0D0**K)
      END DO

      DO M=1, NTERMS
         DO K=1,M
            MPOLE(M) = MPOLE(M) + ANEW(K)*C(M,K)
         END DO
      END DO

      DO M = 1, NTERMS
        MPOLE(M) = MPOLE(M)*Z0P(M) - MPOLE(0)*Z0P(M)/DBLE(M)
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine maps from four polynomial coefficients to
C     the multipole expansion caused by a single layer potential on the
C     side of the box.  The multipole expansion is shifted within this 
C     routine so that it lies at the center of the box.
C
C     INPUT:
C
C     SCOEFFS is are the coefficients approximating the charge density
C     
C     NTERMS is the order of the multipole expansion 
C
C     XLENGTH is the length of a side of the box being considered
C
C     C is an array of binomial coefficients
C
C     OUTPUT:
C
C     MPOLE contains the approximating multipole expansion
C
C**************************************************************************
      SUBROUTINE SINGLESIDETOMULT8(SCOEFFS,NTERMS,C,XLENGTH,MPOLE) 
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS
      REAL *8  SCOEFFS(0:7)
      REAL *8  C(2*NTERMS,2*NTERMS)
      REAL *8  XLENGTH
      COMPLEX *16  MPOLE(0:NTERMS)
C-----Local variables
      INTEGER  I, K, M
      REAL *8  SCOEFFSCALED(0:7)
      COMPLEX *16  Z0, Z00, Z0P(0:60), Z0P2(0:60)
      COMPLEX *16  CD, CDD
      COMPLEX *16  ANEW(0:60), IMAG
      COMPLEX *16  SUM
      DATA IMAG/(0.0,1.0)/

      SCOEFFSCALED(0) = SCOEFFS(0) * XLENGTH
      SCOEFFSCALED(1) = SCOEFFS(1) * XLENGTH
      SCOEFFSCALED(2) = SCOEFFS(2) * XLENGTH
      SCOEFFSCALED(3) = SCOEFFS(3) * XLENGTH
      SCOEFFSCALED(4) = SCOEFFS(4) * XLENGTH
      SCOEFFSCALED(5) = SCOEFFS(5) * XLENGTH
      SCOEFFSCALED(6) = SCOEFFS(6) * XLENGTH
      SCOEFFSCALED(7) = SCOEFFS(7) * XLENGTH

      DO I = 0, NTERMS
         MPOLE(I) = 0.0D0
      END DO

      Z0 = .5D0
      Z0P(0) = 1.0D0
      Z00 = Z0
      CD = 1.0D0 / Z0
      CDD = CD
      Z0P2(0) = 1.0D0
      DO I = 1, NTERMS
         Z0P(I)  = Z00
         Z0P2(I) = CDD
         CDD     = CDD*CD
         Z00     = Z00*Z0
      END DO

C     First, get the zero term in the multipole expansion
C     (Note that even though we are summing over all of the
C     coefficients, that we only need be concerned with 0 and
C     2, because the integrals in the other two cases are zero.)
      MPOLE(0) = SCOEFFSCALED(0) 
     1         + SCOEFFSCALED(2) * (2.0D0/3.0D0 - 1.0D0)
     2         + SCOEFFSCALED(4) * (8.0D0/5.0D0 - 8.0D0/3.0D0 + 1.0D0)
     3         + SCOEFFSCALED(6) * (32.0D0/7.0D0 - 48.0D0/5.0D0 
     4                           +  18.0D0/3.0D0 - 1.0D0)


      DO K = 1, NTERMS
        SUM = 0.0D0
C         We only need to concern ourselves with the cases
C         where i+k is even:
          IF(MOD(K,2) .EQ. 0)THEN

           SUM = SUM + SCOEFFSCALED(0) / DBLE(K+1)

           SUM = SUM + SCOEFFSCALED(2) * 
     1         (2.0D0/DBLE(K+3) - 1.0D0/DBLE(K+1))

           SUM = SUM + SCOEFFSCALED(4) *
     1       (8.0D0/DBLE(K+5) - 8.0D0/DBLE(K+3) + 1.0D0/DBLE(K+1))

           SUM = SUM + SCOEFFSCALED(6) *
     1       (32.0D0/DBLE(K+7) - 48.0D0/DBLE(K+5) 
     2      + 18.0D0/DBLE(K+3) -  1.0D0/DBLE(K+1))

          ELSE

           SUM = SUM + SCOEFFSCALED(1) / DBLE(K+2)

           SUM = SUM + SCOEFFSCALED(3) *
     1         (4.0D0/DBLE(K+4) - 3.0D0/DBLE(K+2))

           SUM = SUM + SCOEFFSCALED(5) *
     1      (16.0D0/DBLE(K+6) - 20.0D0/DBLE(K+4) + 5.0D0/DBLE(K+2))

           SUM = SUM + SCOEFFSCALED(7) *
     1      (64.0D0/DBLE(K+8) - 112.0D0/DBLE(K+6) 
     2     + 56.0D0/DBLE(K+4) -   7.0D0/DBLE(K+2))

          ENDIF
        ANEW(K) = -SUM * Z0P2(K) * IMAG**K / (DBLE(K)*2.0D0**K)
      END DO

      DO M=1, NTERMS
         DO K=1,M
            MPOLE(M) = MPOLE(M) + ANEW(K)*C(M,K)
         END DO
      END DO

      DO M = 1, NTERMS
        MPOLE(M) = MPOLE(M)*Z0P(M) - MPOLE(0)*Z0P(M)/DBLE(M)
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine maps from four polynomial coefficients to
C     the multipole expansion caused by a double potential on the side of
C     the box.  The multipole expansion is shifted within this routine so
C     that it lies at the center of the box.
C
C     INPUT:
C
C     DCOEFFS is are the coefficients approximating the charge density
C     
C     NTERMS is the order of the multipole expansion 
C
C     XLENGTH is the length of a side of the box being considered
C
C     C is an array of binomial coefficients
C
C     OUTPUT:
C
C     MPOLE contains the approximating multipole expansion
C
C**************************************************************************
      SUBROUTINE DOUBLESIDETOMULT8(DCOEFFS,NTERMS,C,MPOLE)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS
      REAL *8  DCOEFFS(0:7)
      REAL *8  C(2*NTERMS,2*NTERMS)
      COMPLEX *16  MPOLE(0:NTERMS)
C-----Local variables
      INTEGER  I, K, M
      COMPLEX *16  Z0, Z00, Z0P(0:60), Z0P2(0:60)
      COMPLEX *16  CD, CDD
      COMPLEX *16  ANEW(0:60), IMAG
      COMPLEX *16  SUM
      DATA IMAG/(0.0,1.0)/

      DO I = 0, NTERMS
         MPOLE(I) = 0.0D0
      END DO

      Z0 = .5D0
      Z0P(0) = 1.0D0
      Z00 = Z0
      CD = 1.0D0 / Z0
      CDD = CD
      Z0P2(0) = 1.0D0
      DO I = 1, NTERMS
         Z0P(I)  = Z00
         Z0P2(I) = CDD
         CDD     = CDD*CD
         Z00     = Z00*Z0
      END DO

C     First, get the zero term in the multipole expansion
C     (Note that even though we are summing over all of the
C     coefficients, that we only need be concerned with 0 and
C     2, because the integrals in the other two cases are zero.)
      MPOLE(0) = DCOEFFS(0) 
     1         + DCOEFFS(2) * (2.0D0/3.0D0 - 1.0D0)
     2         + DCOEFFS(4) * (8.0D0/5.0D0 - 8.0D0/3.0D0 + 1.0D0)
     3         + DCOEFFS(6) * (32.0D0/7.0D0 - 48.0D0/5.0D0 
     4                           +  18.0D0/3.0D0 - 1.0D0)


      DO K = 1, NTERMS
        SUM = 0.0D0
C         We only need to concern ourselves with the cases
C         where i+k is even:
          IF(MOD(K,2) .EQ. 0)THEN

           SUM = SUM + DCOEFFS(0) / DBLE(K+1)

           SUM = SUM + DCOEFFS(2) * 
     1         (2.0D0/DBLE(K+3) - 1.0D0/DBLE(K+1))

           SUM = SUM + DCOEFFS(4) *
     1       (8.0D0/DBLE(K+5) - 8.0D0/DBLE(K+3) + 1.0D0/DBLE(K+1))

           SUM = SUM + DCOEFFS(6) *
     1       (32.0D0/DBLE(K+7) - 48.0D0/DBLE(K+5) 
     2      + 18.0D0/DBLE(K+3) -  1.0D0/DBLE(K+1))

          ELSE

           SUM = SUM + DCOEFFS(1) / DBLE(K+2)

           SUM = SUM + DCOEFFS(3) *
     1         (4.0D0/DBLE(K+4) - 3.0D0/DBLE(K+2))

           SUM = SUM + DCOEFFS(5) *
     1      (16.0D0/DBLE(K+6) - 20.0D0/DBLE(K+4) + 5.0D0/DBLE(K+2))

           SUM = SUM + DCOEFFS(7) *
     1      (64.0D0/DBLE(K+8) - 112.0D0/DBLE(K+6) 
     2     + 56.0D0/DBLE(K+4) -   7.0D0/DBLE(K+2))

          ENDIF
        ANEW(K) = -SUM * Z0P2(K) * IMAG**K / (DBLE(K)*2.0D0**K)
      END DO

      DO M=1, NTERMS
         DO K=1,M
            MPOLE(M) = MPOLE(M) + ANEW(K)*C(M,K)
         END DO
      END DO

      DO M = 1, NTERMS
        MPOLE(M) = MPOLE(M)*Z0P(M) - MPOLE(0)*Z0P(M)/DBLE(M)
      END DO

      DO I = NTERMS, 2, -1
        MPOLE(I) = -DBLE(I-1)*MPOLE(I-1)
      END DO
      MPOLE(1) = MPOLE(0)
      MPOLE(0) = 0.0D0
      RETURN
      END


C**************************************************************************
C     The following subroutine maps from four polynomial coefficients to
C     the multipole expansion caused by a double potential on the top of
C     the box.  The multipole expansion is shifted within this routine so
C     that it lies at the center of the box.
C
C     INPUT:
C
C     DCOEFFS is are the coefficients approximating the charge density
C
C     NTERMS is the order of the multipole expansion
C
C     XLENGTH is the length of a side of the box being considered
C
C     C is an array of binomial coefficients
C
C     OUTPUT:
C
C     MPOLE contains the approximating multipole expansion
C
C**************************************************************************
      SUBROUTINE DOUBLETOPTOMULT8(DCOEFFS,NTERMS,C,MPOLE)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NTERMS
      REAL *8  DCOEFFS(0:7)
      REAL *8  C(2*NTERMS,2*NTERMS)
      COMPLEX *16  MPOLE(0:NTERMS)
C-----Local variables
      INTEGER  I, K, M
      COMPLEX *16  Z0, Z00, Z0P(0:60), Z0P2(0:60)
      COMPLEX *16  CD, CDD
      COMPLEX *16  ANEW(0:60), IMAG
      COMPLEX *16  SUM
      DATA IMAG/(0.0,1.0)/

      DO I = 0, NTERMS
         MPOLE(I) = 0.0D0
      END DO

      Z0 = IMAG * .5D0
      Z0P(0) = 1.0D0
      Z00=Z0
      CD = 1.0D0 / Z0
      CDD = CD
      Z0P2(0) = 1.0D0
      DO I = 1, NTERMS
         Z0P(I)  = Z00
         Z0P2(I) = CDD
         CDD     = CDD*CD
         Z00     = Z00*Z0
      END DO

C     First, get the zero term in the multipole expansion
C     (Note that even though we are summing over all of the
C     coefficients, that we only need be concerned with 0 and
C     2, because the integrals in the other two cases are zero.)
      MPOLE(0) = DCOEFFS(0) 
     1         + DCOEFFS(2) * (2.0D0/3.0D0 - 1.0D0)
     2         + DCOEFFS(4) * (8.0D0/5.0D0 - 8.0D0/3.0D0 + 1.0D0)
     3         + DCOEFFS(6) * (32.0D0/7.0D0 - 48.0D0/5.0D0 
     4                           +  18.0D0/3.0D0 - 1.0D0)


      DO K = 1, NTERMS
        SUM = 0.0D0
C         We only need to concern ourselves with the cases
C         where i+k is even:
          IF(MOD(K,2) .EQ. 0)THEN

           SUM = SUM + DCOEFFS(0) / DBLE(K+1)

           SUM = SUM + DCOEFFS(2) * 
     1         (2.0D0/DBLE(K+3) - 1.0D0/DBLE(K+1))

           SUM = SUM + DCOEFFS(4) *
     1       (8.0D0/DBLE(K+5) - 8.0D0/DBLE(K+3) + 1.0D0/DBLE(K+1))

           SUM = SUM + DCOEFFS(6) *
     1       (32.0D0/DBLE(K+7) - 48.0D0/DBLE(K+5) 
     2      + 18.0D0/DBLE(K+3) -  1.0D0/DBLE(K+1))

          ELSE

           SUM = SUM + DCOEFFS(1) / DBLE(K+2)

           SUM = SUM + DCOEFFS(3) *
     1         (4.0D0/DBLE(K+4) - 3.0D0/DBLE(K+2))

           SUM = SUM + DCOEFFS(5) *
     1      (16.0D0/DBLE(K+6) - 20.0D0/DBLE(K+4) + 5.0D0/DBLE(K+2))

           SUM = SUM + DCOEFFS(7) *
     1      (64.0D0/DBLE(K+8) - 112.0D0/DBLE(K+6) 
     2     + 56.0D0/DBLE(K+4) -   7.0D0/DBLE(K+2))

          ENDIF
        ANEW(K) = -SUM * Z0P2(K) / (DBLE(K)*2.0D0**K)
      END DO

      DO M=1, NTERMS
         DO K=1,M
            MPOLE(M) = MPOLE(M) + ANEW(K)*C(M,K)
         END DO
      END DO

      DO M = 1, NTERMS
        MPOLE(M) = MPOLE(M)*Z0P(M) - MPOLE(0)*Z0P(M)/DBLE(M)
      END DO

      DO I = NTERMS, 2, -1
        MPOLE(I) = -DBLE(I-1)*IMAG*MPOLE(I-1)
      END DO
      MPOLE(1) = IMAG*MPOLE(0)
      MPOLE(0) = 0.0D0
      RETURN
      END


C*****************************************************************************
C     The following subroutine does the local work (just outgoing) involving
C     layer potentials between two boxes that are of the same size.  
C
C     INPUT:
C
C     COEFFSTOP is the array of  basis functions coefficients for the single
C               layer potential on the top
C
C     COEFFSSIDE is the array of  basis functions coefficients for the single
C               layer potential on the side
C
C     COEFFDTOP is the array of  basis functions coefficients for the double
C               layer potential on the top
C
C     COEFFDSIDE is the array of  basis functions coefficients for the double
C               layer potential on the side
C
C     NDEG is the degree of the approximating polynomial
C
C     NB denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the ADAPFMM6 routine)
C
C     WDTOP is the table of weights needed to determine the local 
C           contributions for the double layer on the top
C
C     WDSIDE is the table of weights needed to determine the local 
C           contributions for the double layer on the side
C
C     WSSIDE is the table of weights needed to determine the local 
C           contributions for the single layer on the side
C
C     WSTOP is the table of weights needed to determine the local 
C           contributions for the single layer on the top
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C*****************************************************************************
      SUBROUTINE COLLOCLAY8(POT,COEFFSTOP, 
     1           COEFFDTOP,COEFFSSIDE,COEFFDSIDE,
     2           NDEG,NB,T,WDTOP,WDSIDE,
     3           WSSIDE,WSTOP,ISWITCH,
     3           DSIDESWITCH,DTOPSWITCH,
     4           SSIDESWITCH,STOPSWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NDEG, NB, ISWITCH
      INTEGER  DSIDESWITCH, DTOPSWITCH
      INTEGER  SSIDESWITCH, STOPSWITCH
      REAL *8  POT(64)
      REAL *8  COEFFSTOP(0:NDEG)
      REAL *8  COEFFDTOP(0:NDEG)
      REAL *8  COEFFSSIDE(0:NDEG)
      REAL *8  COEFFDSIDE(0:NDEG)
      REAL *8  T(10)
      REAL *8  WDTOP(9,64,8)
      REAL *8  WDSIDE(9,64,8)
      REAL *8  WSSIDE(9,64,8)
      REAL *8  WSTOP(9,64,8)
C-----Local variables
      INTEGER  I, J
      REAL *8  C(8), CLOG(8)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF


        IF(STOPSWITCH .EQ. 1)THEN
         C(1) = COEFFSTOP(0) / T(1)
         C(2) = COEFFSTOP(1) / T(1)
         C(3) = COEFFSTOP(2) / T(1)
         C(4) = COEFFSTOP(3) / T(1)
         C(5) = COEFFSTOP(4) / T(1)
         C(6) = COEFFSTOP(5) / T(1)
         C(7) = COEFFSTOP(6) / T(1)
         C(8) = COEFFSTOP(7) / T(1)

         CLOG(1)  = C(1) * T(3)
         CLOG(2)  = 0.0D0
         CLOG(3)  = C(3) * T(4)
         CLOG(4)  = 0.0D0
         CLOG(5)  = C(5) * T(5)
         CLOG(6)  = 0.0D0
         CLOG(7)  = C(7) * T(6)
         CLOG(8)  = 0.0D0
       
         DO I = 1,64
          DO J = 1,8
            POT(I) = POT(I) + C(J)*WSTOP(NB,I,J) + CLOG(J)
          ENDDO
         ENDDO
        ENDIF


        IF(SSIDESWITCH .EQ. 1)THEN
         C(1) = COEFFSSIDE(0) / T(1)
         C(2) = COEFFSSIDE(1) / T(1)
         C(3) = COEFFSSIDE(2) / T(1)
         C(4) = COEFFSSIDE(3) / T(1)
         C(5) = COEFFSSIDE(4) / T(1)
         C(6) = COEFFSSIDE(5) / T(1)
         C(7) = COEFFSSIDE(6) / T(1)
         C(8) = COEFFSSIDE(7) / T(1)

         CLOG(1)  = C(1) * T(3)
         CLOG(2)  = 0.0D0
         CLOG(3)  = C(3) * T(4)
         CLOG(4)  = 0.0D0
         CLOG(5)  = C(5) * T(5)
         CLOG(6)  = 0.0D0
         CLOG(7)  = C(7) * T(6)
         CLOG(8)  = 0.0D0

         DO I = 1,64
          DO J = 1,8
            POT(I) = POT(I) + C(J)*WSSIDE(NB,I,J) + CLOG(J)
          ENDDO
         ENDDO
        ENDIF


        IF(DTOPSWITCH .EQ. 1)THEN
         C(1) = COEFFDTOP(0)
         C(2) = COEFFDTOP(1)
         C(3) = COEFFDTOP(2)
         C(4) = COEFFDTOP(3)
         C(5) = COEFFDTOP(4)
         C(6) = COEFFDTOP(5)
         C(7) = COEFFDTOP(6)
         C(8) = COEFFDTOP(7)

         DO I = 1,64
          DO J = 1,8
            POT(I) = POT(I) + C(J)*WDTOP(NB,I,J)
          ENDDO
         ENDDO
        ENDIF


        IF(DSIDESWITCH .EQ. 1)THEN
         C(1) = COEFFDSIDE(0)
         C(2) = COEFFDSIDE(1)
         C(3) = COEFFDSIDE(2)
         C(4) = COEFFDSIDE(3)
         C(5) = COEFFDSIDE(4)
         C(6) = COEFFDSIDE(5)
         C(7) = COEFFDSIDE(6)
         C(8) = COEFFDSIDE(7)

         DO I = 1,64
          DO J = 1,8
            POT(I) = POT(I) + C(J)*WDSIDE(NB,I,J)
          ENDDO
         ENDDO
        ENDIF
       RETURN
       END



C*****************************************************************************
C     The following subroutine does the local work involving layer potentials
C     that figures the potential in a smaller box due to a big box that is
C     touching it.
C     The boxes are only allowed to be one level apart.
C
C     INPUT:
C
C     COEFFSTOP is the array of  basis functions coefficients for the single
C                 layer potential on the top
C
C     COEFFSSIDE is the array of  basis functions coefficients for the single
C                 layer potential on the side
C
C     COEFFDTOP is the array of  basis functions coefficients for the double
C                 layer potential on the top
C
C     COEFFDSIDE is the array of  basis functions coefficients for the double
C                 layer potential on the side
C
C     NDEG is the degree of the approximating polynomial
C
C     IFLAG denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the BOUNDFMM6 routine)
C
C     WDTOP is the table of weights needed to determine the local 
C           contributions from the double layer on top
C
C     WDSIDE is the table of weights needed to determine the local 
C           contributions from the double layer on side
C
C     WSTOP is the table of weights needed to determine the local 
C           contributions from the single layer on top
C
C     WSSIDE is the table of weights needed to determine the local 
C           contributions from the single layer on side
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C*****************************************************************************
      SUBROUTINE BTOSLOCLAY8(POT,COEFFSTOP,
     1            COEFFSSIDE,COEFFDTOP,COEFFDSIDE,
     2            NDEG,IFLAG,T,WDTOP,WDSIDE,
     3            WSTOP,WSSIDE,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, IFLAG, ISWITCH
      REAL *8 POT(64), T(10)
      REAL *8 COEFFSTOP(0:NDEG), COEFFSSIDE(0:NDEG)
      REAL *8 COEFFDTOP(0:NDEG), COEFFDSIDE(0:NDEG)
      REAL *8 WDTOP(12,64,8), WDSIDE(12,64,8)
      REAL *8 WSTOP(12,64,8), WSSIDE(12,64,8)
C-----Local variables
      INTEGER I, J
      REAL *8 C(8), CLOG(8)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF

      C(1)  = COEFFDTOP(0)
      C(2)  = COEFFDTOP(1)
      C(3)  = COEFFDTOP(2)
      C(4)  = COEFFDTOP(3)
      C(5)  = COEFFDTOP(4)
      C(6)  = COEFFDTOP(5)
      C(7)  = COEFFDTOP(6)
      C(8)  = COEFFDTOP(7)

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WDTOP(IFLAG,I,J)
        END DO
      END DO


      C(1)  = COEFFDSIDE(0)
      C(2)  = COEFFDSIDE(1)
      C(3)  = COEFFDSIDE(2)
      C(4)  = COEFFDSIDE(3)
      C(5)  = COEFFDSIDE(4)
      C(6)  = COEFFDSIDE(5)
      C(7)  = COEFFDSIDE(6)
      C(8)  = COEFFDSIDE(7)

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WDSIDE(IFLAG,I,J)
        END DO
      END DO


      C(1) = COEFFSTOP(0)  / T(1)
      C(2) = COEFFSTOP(1)  / T(1)
      C(3) = COEFFSTOP(2)  / T(1)
      C(4) = COEFFSTOP(3)  / T(1)
      C(5) = COEFFSTOP(4)  / T(1)
      C(6) = COEFFSTOP(5)  / T(1)
      C(7) = COEFFSTOP(6)  / T(1)
      C(8) = COEFFSTOP(7)  / T(1)

      CLOG(1)  = C(1) * T(3)
      CLOG(2)  = 0.0D0
      CLOG(3)  = C(3) * T(4)
      CLOG(4)  = 0.0D0
      CLOG(5)  = C(5) * T(5)
      CLOG(6)  = 0.0D0
      CLOG(7)  = C(7) * T(6)
      CLOG(8)  = 0.0D0

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)* WSTOP(IFLAG,I,J) + CLOG(J)
        END DO
      END DO


      C(1) = COEFFSSIDE(0) / T(1)
      C(2) = COEFFSSIDE(1) / T(1)
      C(3) = COEFFSSIDE(2) / T(1)
      C(4) = COEFFSSIDE(3) / T(1)
      C(5) = COEFFSSIDE(4) / T(1)
      C(6) = COEFFSSIDE(5) / T(1)
      C(7) = COEFFSSIDE(6) / T(1)
      C(8) = COEFFSSIDE(7) / T(1)

      CLOG(1)  = C(1) * T(3)
      CLOG(2)  = 0.0D0
      CLOG(3)  = C(3) * T(4)
      CLOG(4)  = 0.0D0
      CLOG(5)  = C(5) * T(5)
      CLOG(6)  = 0.0D0
      CLOG(7)  = C(7) * T(6)
      CLOG(8)  = 0.0D0

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J) * WSSIDE(IFLAG,I,J) + CLOG(J)
        END DO
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine does the local work involving layer potentials
C     that figures the potential in a bigger box due to a small box that is
C     touching it.
C     The boxes are only allowed to be one level apart.
C
C     INPUT:
C
C     COEFFSTOP is the array of  basis functions coefficients for the single
C                 layer potential on the top
C
C     COEFFSSIDE is the array of  basis functions coefficients for the single
C                 layer potential on the side
C
C     COEFFDTOP is the array of  basis functions coefficients for the double
C                 layer potential on the top
C
C     COEFFDSIDE is the array of  basis functions coefficients for the double
C                 layer potential on the side
C
C     NDEG is the degree of the approximating polynomial
C
C     IFLAG denotes the colleague number being considered
C
C     T is a set of precomputed scalings (set in the BOUNDFMM6 routine)
C
C     WDTOP is the table of weights needed to determine the local 
C           contributions from the double layer on top
C
C     WDSIDE is the table of weights needed to determine the local 
C           contributions from the double layer on side
C
C     WSTOP is the table of weights needed to determine the local
C           contributions from the single layer on top
C
C     WSSIDE is the table of weights needed to determine the local
C           contributions from the single layer on side
C
C     OUTPUT:
C
C     POT is altered to account for the local contribution
C
C**************************************************************************
      SUBROUTINE STOBLOCLAY8(POT,COEFFSTOP, 
     1            COEFFSSIDE,COEFFDTOP,COEFFDSIDE, 
     2            NDEG,IFLAG,T,WDTOP,WDSIDE,
     3            WSTOP,WSSIDE,ISWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER NDEG, IFLAG, ISWITCH
      REAL *8 POT(64), T(10)
      REAL *8 COEFFSTOP(0:NDEG), COEFFSSIDE(0:NDEG)
      REAL *8 COEFFDTOP(0:NDEG), COEFFDSIDE(0:NDEG)
      REAL *8 WDTOP(12,64,8), WDSIDE(12,64,8)
      REAL *8 WSTOP(12,64,8), WSSIDE(12,64,8)
C-----Local variables
      INTEGER I, J
      REAL *8 C(8), CLOG(8)

      IF(ISWITCH .EQ. 0)THEN
        RETURN
      ENDIF
    
      C(1)  = COEFFDTOP(0)
      C(2)  = COEFFDTOP(1)
      C(3)  = COEFFDTOP(2)
      C(4)  = COEFFDTOP(3)
      C(5)  = COEFFDTOP(4)
      C(6)  = COEFFDTOP(5)
      C(7)  = COEFFDTOP(6)
      C(8)  = COEFFDTOP(7)

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WDTOP(IFLAG,I,J)
        END DO
      END DO


      C(1)  = COEFFDSIDE(0)
      C(2)  = COEFFDSIDE(1)
      C(3)  = COEFFDSIDE(2)
      C(4)  = COEFFDSIDE(3)
      C(5)  = COEFFDSIDE(4)
      C(6)  = COEFFDSIDE(5)
      C(7)  = COEFFDSIDE(6)
      C(8)  = COEFFDSIDE(7)

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WDSIDE(IFLAG,I,J)
        END DO
      END DO


      C(1) = COEFFSTOP(0)  / T(2)
      C(2) = COEFFSTOP(1)  / T(2)
      C(3) = COEFFSTOP(2)  / T(2)
      C(4) = COEFFSTOP(3)  / T(2)
      C(5) = COEFFSTOP(4)  / T(2)
      C(6) = COEFFSTOP(5)  / T(2)
      C(7) = COEFFSTOP(6)  / T(2)
      C(8) = COEFFSTOP(7)  / T(2)

      CLOG(1)  = C(1) * T(7)
      CLOG(2)  = 0.0D0
      CLOG(3)  = C(3) * T(8)
      CLOG(4)  = 0.0D0
      CLOG(5)  = C(5) * T(9)
      CLOG(6)  = 0.0D0
      CLOG(7)  = C(7) * T(10)
      CLOG(8)  = 0.0D0

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WSTOP(IFLAG,I,J) + CLOG(J)
        END DO
      END DO


      C(1) = COEFFSSIDE(0) / T(2)
      C(2) = COEFFSSIDE(1) / T(2)
      C(3) = COEFFSSIDE(2) / T(2)
      C(4) = COEFFSSIDE(3) / T(2)
      C(5) = COEFFSSIDE(4) / T(2)
      C(6) = COEFFSSIDE(5) / T(2)
      C(7) = COEFFSSIDE(6) / T(2)
      C(8) = COEFFSSIDE(7) / T(2)

      CLOG(1)  = C(1) * T(7)
      CLOG(2)  = 0.0D0
      CLOG(3)  = C(3) * T(8)
      CLOG(4)  = 0.0D0
      CLOG(5)  = C(5) * T(9)
      CLOG(6)  = 0.0D0
      CLOG(7)  = C(7) * T(10)
      CLOG(8)  = 0.0D0

      DO I = 1, 64
        DO J = 1, 8
           POT(I) = POT(I) + C(J)*WSSIDE(IFLAG,I,J) + CLOG(J)
        END DO
      END DO
      RETURN
      END


C*****************************************************************************
C     The following subroutine generates four exponential expansions, north,
C     south, east, and west respectively, that are use when passing from a 
C     larger box to a smaller box.  These expansions represent the layer
C     potentials of the given box.  The exponential expansions are obtained
C     by using a set of precomputed matrices that map from the polynomial
C     coefficients (that approximate the charge density) to the exponential
C     coefficients.
C     
C     INPUT:
C
C     COEFFDTOP is the array of  basis functions coefficients for the double
C               layer potential on the top
C
C     COEFFSTOP is the array of  basis functions coefficients for the single
C               layer potential on the top
C
C     COEFFDSIDE is the array of  basis functions coefficients for the double
C               layer potential on the side
C
C     COEFFSSIDE is the array of  basis functions coefficients for the single
C               layer potential on the side
C
C     LEVEL is the level of the larger box
C
C     NNODES is the order of the exponential expansion
C
C     SUM1 is a precomputed terms needed in the exponential expansions
C
C     The other input entries are just the precomputed maps from the
C     polynomial coefficients to the exponential coefficients.  The 
C     are written as direction, double or single, and top or side.
C     So EASTDBLETOP would be the map from the four polynomial
C     coefficients for the double layer on top to the expansion going in the
C     east direction.  The other matrices are defined similarly.
C
C     OUTPUT:
C
C     EXPNALL is the outgoing north expansion for the big box
C
C     EXPSALL is the outgoing south expansion for the big box
C
C     EXPEALL is the outgoing east expansion for the big box
C
C     EXPWALL is the outgoing west expansion for the big box
C
C*****************************************************************************
      SUBROUTINE MKEXPBTOSLAY8(COEFFDTOP,COEFFSTOP, 
     1      COEFFDSIDE, COEFFSSIDE,XLENGTH,SUM1,
     2      EXPNALL,EXPSALL,EXPEALL,EXPWALL,
     4      EASTDBLETOPMAT,   WESTDBLETOPMAT,
     5      NORTHDBLETOPMAT,  SOUTHDBLETOPMAT,
     5      EASTDBLESIDEMAT,  WESTDBLESIDEMAT,
     5      NORTHDBLESIDEMAT, SOUTHDBLESIDEMAT,
     6      EASTSNGLETOPMAT,  WESTSNGLETOPMAT,
     5      SOUTHSNGLETOPMAT, NORTHSNGLETOPMAT,
     7      EASTSNGLESIDEMAT, WESTSNGLESIDEMAT, 
     8      SOUTHSNGLESIDEMAT,NORTHSNGLESIDEMAT,NNODES,
     9      DSIDESWITCH,DTOPSWITCH,
     1      SSIDESWITCH,STOPSWITCH)
      IMPLICIT NONE
C-----Global variables
      INTEGER  I, II, NNODES
      INTEGER  DSIDESWITCH, DTOPSWITCH
      INTEGER  SSIDESWITCH, STOPSWITCH
      REAL *8  COEFFDTOP(0:7), COEFFSTOP(0:7)
      REAL *8  COEFFDSIDE(0:7), COEFFSSIDE(0:7)
      REAL *8  SUM1
      COMPLEX *16  EXPNALL(0:60), EXPSALL(0:60)
      COMPLEX *16  EXPEALL(0:60), EXPWALL(0:60)
C-----Local variables
      REAL *8  XLENGTH, SUM2
      REAL *8  COEFFSTOPCALED(0:7)
      REAL *8  COEFFSSIDECALED(0:7)
      COMPLEX *16  EXPNTEST(0:60), EXPSTEST(0:60)
      COMPLEX *16  EXPETEST(0:60), EXPWTEST(0:60)
      COMPLEX *16  EASTDBLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  EASTDBLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  EASTSNGLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  EASTSNGLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  WESTDBLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  WESTDBLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  WESTSNGLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  WESTSNGLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  NORTHDBLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  NORTHDBLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  NORTHSNGLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  NORTHSNGLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  SOUTHDBLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  SOUTHDBLESIDEMAT(0:NNODES,0:7)
      COMPLEX *16  SOUTHSNGLETOPMAT(0:NNODES,0:7)
      COMPLEX *16  SOUTHSNGLESIDEMAT(0:NNODES,0:7)

      SUM2 = SUM1 + DLOG(XLENGTH)

      DO I = 0, NNODES
        EXPEALL(I) = 0.0D0
        EXPWALL(I) = 0.0D0
        EXPNALL(I) = 0.0D0
        EXPSALL(I) = 0.0D0
      END DO

      IF(DTOPSWITCH .EQ. 1)THEN
        DO I = 0, NNODES
        EXPETEST(I) = 0.0D0
        DO II = 0, 7
         EXPETEST(I) = EXPETEST(I)
     1        + EASTDBLETOPMAT(I,II) * COEFFDTOP(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPWTEST(I) = 0.0D0
        DO II = 0, 7
         EXPWTEST(I) = EXPWTEST(I)
     1        + WESTDBLETOPMAT(I,II) * COEFFDTOP(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPNTEST(I) = 0.0D0
        DO II = 0, 7
         EXPNTEST(I) = EXPNTEST(I)
     1        + NORTHDBLETOPMAT(I,II) * COEFFDTOP(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPSTEST(I) = 0.0D0
        DO II = 0, 7
         EXPSTEST(I) = EXPSTEST(I)
     1        + SOUTHDBLETOPMAT(I,II) * COEFFDTOP(II)
        END DO
        END DO

        DO I = 0, NNODES
          EXPEALL(I) = EXPEALL(I) + EXPETEST(I)
          EXPWALL(I) = EXPWALL(I) + EXPWTEST(I)
          EXPNALL(I) = EXPNALL(I) + EXPNTEST(I)
          EXPSALL(I) = EXPSALL(I) + EXPSTEST(I)
        END DO
      ENDIF


      IF(DSIDESWITCH .EQ. 1)THEN
        DO I = 0, NNODES
        EXPETEST(I) = 0.0D0
        DO II = 0, 7
         EXPETEST(I) = EXPETEST(I)
     1        + EASTDBLESIDEMAT(I,II) * COEFFDSIDE(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPWTEST(I) = 0.0D0
        DO II = 0, 7
         EXPWTEST(I) = EXPWTEST(I)
     1        + WESTDBLESIDEMAT(I,II) * COEFFDSIDE(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPNTEST(I) = 0.0D0
        DO II = 0, 7
         EXPNTEST(I) = EXPNTEST(I)
     1        + NORTHDBLESIDEMAT(I,II) * COEFFDSIDE(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPSTEST(I) = 0.0D0
        DO II = 0, 7
         EXPSTEST(I) = EXPSTEST(I)
     1        + SOUTHDBLESIDEMAT(I,II) * COEFFDSIDE(II)
        END DO
        END DO

        DO I = 0, NNODES
          EXPEALL(I) = EXPEALL(I) + EXPETEST(I)
          EXPWALL(I) = EXPWALL(I) + EXPWTEST(I)
          EXPNALL(I) = EXPNALL(I) + EXPNTEST(I)
          EXPSALL(I) = EXPSALL(I) + EXPSTEST(I)
        END DO
      ENDIF


      IF(STOPSWITCH .EQ. 1)THEN
        COEFFSTOPCALED(0) = COEFFSTOP(0) * XLENGTH
        COEFFSTOPCALED(1) = COEFFSTOP(1) * XLENGTH
        COEFFSTOPCALED(2) = COEFFSTOP(2) * XLENGTH
        COEFFSTOPCALED(3) = COEFFSTOP(3) * XLENGTH
        COEFFSTOPCALED(4) = COEFFSTOP(4) * XLENGTH
        COEFFSTOPCALED(5) = COEFFSTOP(5) * XLENGTH
        COEFFSTOPCALED(6) = COEFFSTOP(6) * XLENGTH
        COEFFSTOPCALED(7) = COEFFSTOP(7) * XLENGTH

        DO I = 0, NNODES
        EXPETEST(I) = 0.0D0
        DO II = 0, 7
         EXPETEST(I) = EXPETEST(I)
     1        + EASTSNGLETOPMAT(I,II) * COEFFSTOPCALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPWTEST(I) = 0.0D0
        DO II = 0, 7
         EXPWTEST(I) = EXPWTEST(I)
     1        + WESTSNGLETOPMAT(I,II) * COEFFSTOPCALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPSTEST(I) = 0.0D0
        DO II = 0, 7
         EXPSTEST(I) = EXPSTEST(I)
     1        + SOUTHSNGLETOPMAT(I,II) * COEFFSTOPCALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPNTEST(I) = 0.0D0
        DO II = 0, 7
         EXPNTEST(I) = EXPNTEST(I)
     1        + NORTHSNGLETOPMAT(I,II) * COEFFSTOPCALED(II)
        END DO
        END DO

        DO I = 0, NNODES
          EXPEALL(I) = EXPEALL(I) + EXPETEST(I)
          EXPWALL(I) = EXPWALL(I) + EXPWTEST(I)
          EXPSALL(I) = EXPSALL(I) + EXPSTEST(I)
          EXPNALL(I) = EXPNALL(I) + EXPNTEST(I)
        END DO
      ENDIF


      IF(SSIDESWITCH .EQ. 1)THEN
        COEFFSSIDECALED(0) = COEFFSSIDE(0) * XLENGTH
        COEFFSSIDECALED(1) = COEFFSSIDE(1) * XLENGTH
        COEFFSSIDECALED(2) = COEFFSSIDE(2) * XLENGTH
        COEFFSSIDECALED(3) = COEFFSSIDE(3) * XLENGTH
        COEFFSSIDECALED(4) = COEFFSSIDE(4) * XLENGTH
        COEFFSSIDECALED(5) = COEFFSSIDE(5) * XLENGTH
        COEFFSSIDECALED(6) = COEFFSSIDE(6) * XLENGTH
        COEFFSSIDECALED(7) = COEFFSSIDE(7) * XLENGTH

        DO I = 0, NNODES
        EXPETEST(I) = 0.0D0
        DO II = 0, 7
         EXPETEST(I) = EXPETEST(I)
     1        + EASTSNGLESIDEMAT(I,II) * COEFFSSIDECALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPWTEST(I) = 0.0D0
        DO II = 0, 7
         EXPWTEST(I) = EXPWTEST(I)
     1        + WESTSNGLESIDEMAT(I,II) * COEFFSSIDECALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPSTEST(I) = 0.0D0
        DO II = 0, 7
         EXPSTEST(I) = EXPSTEST(I)
     1        + SOUTHSNGLESIDEMAT(I,II) * COEFFSSIDECALED(II)
        END DO
        END DO

        DO I = 0, NNODES
        EXPNTEST(I) = 0.0D0
        DO II = 0, 7
         EXPNTEST(I) = EXPNTEST(I)
     1        + NORTHSNGLESIDEMAT(I,II) * COEFFSSIDECALED(II)
        END DO
        END DO

        DO I = 0, NNODES
          EXPEALL(I) = EXPEALL(I) + EXPETEST(I)
          EXPWALL(I) = EXPWALL(I) + EXPWTEST(I)
          EXPNALL(I) = EXPNALL(I) + EXPNTEST(I)
          EXPSALL(I) = EXPSALL(I) + EXPSTEST(I)
        END DO
      ENDIF

      EXPEALL(0) = EXPEALL(0) * SUM2
      EXPWALL(0) = EXPWALL(0) * SUM2
      EXPNALL(0) = EXPNALL(0) * SUM2
      EXPSALL(0) = EXPSALL(0) * SUM2
      RETURN
      END



C**************************************************************************
C     The following subroutine sets up boundary conditions for the
C     case with purely dirichlet conditions. 
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C    
C     HL is the function representing the boundary condition imposed
C        on the left
C    
C     HR is the function representing the boundary condition imposed
C        on the right
C    
C     HB is the function representing the boundary condition imposed
C        on the bottom
C    
C     HT is the function representing the boundary condition imposed
C        on the top
C
C     OUTPUT:
C
C     FTOP is the array of function values representing the single layer
C          densities on the top of each box
C
C     FSIDE is the array of function values representing the single layer
C          densities on the side of each box
C
C     FTOPD is the array of function values representing the double layer
C          densities on the top of each box
C
C     FSIDED is the array of function values representing the double layer
C          densities on the side of each box
C
C**************************************************************************
      SUBROUTINE SETBCS8DIR(NBOXES, NLEV, ICOLBOX, IROWBOX,
     1           ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2           FTOP, FSIDE, FTOPD, FSIDED, 
     5           DOUBLETOPONOFF, DOUBLESIDEONOFF,
     6           SINGLETOPONOFF, SINGLESIDEONOFF,
     7           HL, HR, HB, HT)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XSHIFT
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  HL, HR, HB, HT
C-----Local variables
      INTEGER  I, JJ, II, ILEV
      INTEGER  NSIDE
      INTEGER  ISTART, IEND
      REAL *8  TEMP3
      REAL *8  XUSE, PI, PI16, X, XX(8)

      PI = DACOS(-1.0D0)

      PI16 = PI / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      DO II = 1, NBOXES
        DOUBLETOPONOFF(II) = 0
        DOUBLESIDEONOFF(II) = 0
        SINGLETOPONOFF(II) = 0
        SINGLESIDEONOFF(II) = 0
        DO JJ = 1, 8
           FTOP(JJ,II) = 0.0D0
           FTOPD(JJ,II) = 0.0D0
           FSIDE(JJ,II) = 0.0D0
           FSIDED(JJ,II) = 0.0D0
        END DO
      END DO

      TEMP3 = .5D0
      NSIDE = 1
      DO ILEV = 0, NLEV
      ISTART = ISTARTLEV(ILEV)
      IEND = ISTART + NBLEVEL(ILEV) - 1
      DO I = ISTART, IEND
      II = IBOXLEV(I)
      IF(ICHILDBOX(1,II) .LT. 0)THEN

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ) / TEMP3
           FSIDED(JJ,II) = HR(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDED(JJ,II) = -HR(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FSIDED(JJ,II) = -HL(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDED(JJ,II) = HL(XUSE)
         END DO
       ENDIF
       ENDIF


C     Now set the top BCs
       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         DOUBLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3  
           FTOPD(JJ,II) = HT(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         DOUBLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOPD(JJ,II) = -HT(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         DOUBLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FTOPD(JJ,II) = -HB(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         DOUBLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOPD(JJ,II) = HB(XUSE)
         END DO
       ENDIF
       ENDIF

      ENDIF
      END DO
      TEMP3 = 2.0D0*TEMP3
      NSIDE = 2*NSIDE
      END DO

C     Now scale everything appropriately
C     and account for the sign changes:
      DO II = 1, NBOXES
        DO JJ = 1, 8
           FTOPD(JJ,II)  = -FTOPD(JJ,II) / PI
           FSIDED(JJ,II) = -FSIDED(JJ,II) / PI
        END DO
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine sets up boundary conditions for the
C     case with purely neumann conditions.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     HL is the function representing the boundary condition imposed
C        on the left
C
C     HR is the function representing the boundary condition imposed
C        on the right
C
C     HB is the function representing the boundary condition imposed
C        on the bottom
C
C     HT is the function representing the boundary condition imposed
C        on the top
C
C     OUTPUT:
C
C     FTOP is the array of function values representing the single layer
C          densities on the top of each box
C
C     FSIDE is the array of function values representing the single layer
C          densities on the side of each box
C
C     FTOPD is the array of function values representing the double layer
C          densities on the top of each box
C
C     FSIDED is the array of function values representing the double layer
C          densities on the side of each box
C
C**************************************************************************
      SUBROUTINE SETBCS8NEU(NBOXES, NLEV, ICOLBOX, IROWBOX,
     1           ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2           FTOP, FSIDE, FTOPD, FSIDED, 
     5           DOUBLETOPONOFF, DOUBLESIDEONOFF,
     6           SINGLETOPONOFF, SINGLESIDEONOFF,
     7           HL, HR, HB, HT)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XSHIFT
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  HL, HR, HB, HT
C-----Local variables
      INTEGER  I, JJ, II, ILEV
      INTEGER  NSIDE
      INTEGER  ISTART, IEND
      REAL *8  TEMP3
      REAL *8  XUSE, PI, PI16, X, XX(8)

      PI = DACOS(-1.0D0)

      PI16 = PI / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      DO II = 1, NBOXES
        DOUBLETOPONOFF(II) = 0
        DOUBLESIDEONOFF(II) = 0
        SINGLETOPONOFF(II) = 0
        SINGLESIDEONOFF(II) = 0
        DO JJ = 1, 8
           FTOP(JJ,II) = 0.0D0
           FTOPD(JJ,II) = 0.0D0
           FSIDE(JJ,II) = 0.0D0
           FSIDED(JJ,II) = 0.0D0
        END DO
      END DO

      TEMP3 = .5D0
      NSIDE = 1
      DO ILEV = 0, NLEV
      ISTART = ISTARTLEV(ILEV)
      IEND = ISTART + NBLEVEL(ILEV) - 1
      DO I = ISTART, IEND
      II = IBOXLEV(I)
      IF(ICHILDBOX(1,II) .LT. 0)THEN

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ) / TEMP3
           FSIDE(JJ,II) = HR(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDE(JJ,II) = HR(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FSIDE(JJ,II) = HL(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDE(JJ,II) = HL(XUSE)
         END DO
       ENDIF
       ENDIF


C     Now set the top BCs
       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3  
           FTOP(JJ,II) = HT(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOP(JJ,II) = HT(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FTOP(JJ,II) = HB(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOP(JJ,II) = HB(XUSE)
         END DO
       ENDIF
       ENDIF

      ENDIF
      END DO
      TEMP3 = 2.0D0*TEMP3
      NSIDE = 2*NSIDE
      END DO

C     Now scale everything appropriately
C     and account for the sign changes:
      DO II = 1, NBOXES
        DO JJ = 1, 8
           FTOP(JJ,II)  = -FTOP(JJ,II) / PI
           FSIDE(JJ,II) = -FSIDE(JJ,II) / PI
        END DO
      END DO
      RETURN
      END

C**************************************************************************
C     The following subroutine sets up boundary conditions for the
C     case with neumann conditions on the top and bottom and dirichlet
C     conditions on the sides.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     HL is the function representing the boundary condition imposed
C        on the left
C
C     HR is the function representing the boundary condition imposed
C        on the right
C
C     HB is the function representing the boundary condition imposed
C        on the bottom
C
C     HT is the function representing the boundary condition imposed
C        on the top
C
C     OUTPUT:
C
C     FTOP is the array of function values representing the single layer
C          densities on the top of each box
C
C     FSIDE is the array of function values representing the single layer
C          densities on the side of each box
C
C     FTOPD is the array of function values representing the double layer
C          densities on the top of each box
C
C     FSIDED is the array of function values representing the double layer
C          densities on the side of each box
C
C**************************************************************************
      SUBROUTINE SETBCS8DIRNEU(NBOXES, NLEV, ICOLBOX, IROWBOX,
     1           ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2           FTOP, FSIDE, FTOPD, FSIDED,
     5           DOUBLETOPONOFF, DOUBLESIDEONOFF,
     6           SINGLETOPONOFF, SINGLESIDEONOFF,
     7           HL, HR, HB, HT)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XSHIFT
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  HL, HR, HB, HT
C-----Local variables
      INTEGER  I, JJ, II, ILEV
      INTEGER  NSIDE
      INTEGER  ISTART, IEND
      REAL *8  TEMP3
      REAL *8  XUSE, PI, PI16, X, XX(8)

      PI = DACOS(-1.0D0)

      PI16 = PI / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      DO II = 1, NBOXES
        DOUBLETOPONOFF(II) = 0
        DOUBLESIDEONOFF(II) = 0
        SINGLETOPONOFF(II) = 0
        SINGLESIDEONOFF(II) = 0
        DO JJ = 1, 8
           FTOP(JJ,II) = 0.0D0
           FTOPD(JJ,II) = 0.0D0
           FSIDE(JJ,II) = 0.0D0
           FSIDED(JJ,II) = 0.0D0
        END DO
      END DO

      TEMP3 = .5D0
      NSIDE = 1
      DO ILEV = 0, NLEV
      ISTART = ISTARTLEV(ILEV)
      IEND = ISTART + NBLEVEL(ILEV) - 1
      DO I = ISTART, IEND
      II = IBOXLEV(I)
      IF(ICHILDBOX(1,II) .LT. 0)THEN

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ) / TEMP3
           FSIDED(JJ,II) = HR(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDED(JJ,II) = HR(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FSIDED(JJ,II) = -HL(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FSIDED(JJ,II) = -HL(XUSE)
         END DO
       ENDIF
       ENDIF


C     Now set the top BCs
       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3  
           FTOP(JJ,II) = HT(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE/2)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOP(JJ,II) = -HT(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .LE. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FTOP(JJ,II) = HB(X)
         END DO
       ENDIF
       ENDIF

       IF(IROWBOX(II) .EQ. NSIDE)THEN
       IF(ICOLBOX(II) .GT. NSIDE/2)THEN
         SINGLETOPONOFF(II) = 1
         XSHIFT  = DBLE(ICOLBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           XUSE = -XUSE
           FTOP(JJ,II) = -HB(XUSE)
         END DO
       ENDIF
       ENDIF

      ENDIF
      END DO
      TEMP3 = 2.0D0*TEMP3
      NSIDE = 2*NSIDE
      END DO

C     Now scale everything appropriately
C     and account for the sign changes:
      DO II = 1, NBOXES
        DO JJ = 1, 8
           FTOP(JJ,II)  = -FTOP(JJ,II) / PI
           FSIDED(JJ,II) = -FSIDED(JJ,II) / PI
        END DO
      END DO
      RETURN
      END

C**************************************************************************
C     The following subroutine sets up boundary conditions for the
C     case with dirichlet conditions on the side and periodic top bottom.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     HL is the function representing the boundary condition imposed
C        on the left
C
C     HR is the function representing the boundary condition imposed
C        on the right
C
C     OUTPUT:
C
C     FTOP is the array of function values representing the single layer
C          densities on the top of each box
C
C     FSIDE is the array of function values representing the single layer
C          densities on the side of each box
C
C     FTOPD is the array of function values representing the double layer
C          densities on the top of each box
C
C     FSIDED is the array of function values representing the double layer
C          densities on the side of each box
C
C**************************************************************************
      SUBROUTINE SETBCS8DIRPER(NBOXES, NLEV, ICOLBOX, IROWBOX,
     1           ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2           FTOP, FSIDE, FTOPD, FSIDED,
     5           DOUBLETOPONOFF, DOUBLESIDEONOFF,
     6           SINGLETOPONOFF, SINGLESIDEONOFF,
     7           HL, HR)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XSHIFT
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  HL, HR
C-----Local variables
      INTEGER  I, JJ, II, ILEV
      INTEGER  NSIDE
      INTEGER  ISTART, IEND
      REAL *8  TEMP3
      REAL *8  XUSE, PI, PI16, X, XX(8)

      PI = DACOS(-1.0D0)

      PI16 = PI / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      DO II = 1, NBOXES
        DOUBLETOPONOFF(II) = 0
        DOUBLESIDEONOFF(II) = 0
        SINGLETOPONOFF(II) = 0
        SINGLESIDEONOFF(II) = 0
        DO JJ = 1, 8
           FTOP(JJ,II) = 0.0D0
           FTOPD(JJ,II) = 0.0D0
           FSIDE(JJ,II) = 0.0D0
           FSIDED(JJ,II) = 0.0D0
        END DO
      END DO

      TEMP3 = .5D0
      NSIDE = 1
      DO ILEV = 0, NLEV
      ISTART = ISTARTLEV(ILEV)
      IEND = ISTART + NBLEVEL(ILEV) - 1
      DO I = ISTART, IEND
      II = IBOXLEV(I)
      IF(ICHILDBOX(1,II) .LT. 0)THEN

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ) / TEMP3
           FSIDED(JJ,II) = HR(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           FSIDED(JJ,II) = HR(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FSIDED(JJ,II) = -HL(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         DOUBLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           FSIDED(JJ,II) = -HL(XUSE)
         END DO
       ENDIF
       ENDIF

      ENDIF
      END DO
      TEMP3 = 2.0D0*TEMP3
      NSIDE = 2*NSIDE
      END DO

C     Now scale everything appropriately
C     and account for the sign changes:
      DO II = 1, NBOXES
        DO JJ = 1, 8
           FSIDED(JJ,II) = -FSIDED(JJ,II) / PI
        END DO
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine sets up boundary conditions for the
C     case with neumann conditions on the side and periodic top bottom.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     ICHILDBOX denotes the four children of each box
C
C     IPARENTBOX denotes the parent of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     HL is the function representing the boundary condition imposed
C        on the left
C
C     HR is the function representing the boundary condition imposed
C        on the right
C
C     OUTPUT:
C
C     FTOP is the array of function values representing the single layer
C          densities on the top of each box
C
C     FSIDE is the array of function values representing the single layer
C          densities on the side of each box
C
C     FTOPD is the array of function values representing the double layer
C          densities on the top of each box
C
C     FSIDED is the array of function values representing the double layer
C          densities on the side of each box
C
C**************************************************************************
      SUBROUTINE SETBCS8NEUPER(NBOXES, NLEV, ICOLBOX, IROWBOX,
     1           ICHILDBOX,  NBLEVEL, IBOXLEV, ISTARTLEV,
     2           FTOP, FSIDE, FTOPD, FSIDED, 
     5           DOUBLETOPONOFF, DOUBLESIDEONOFF,
     6           SINGLETOPONOFF, SINGLESIDEONOFF,
     7           HL, HR)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NBOXES, NLEV
      INTEGER  ICOLBOX(1), IROWBOX(1)
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      INTEGER  DOUBLETOPONOFF(1), DOUBLESIDEONOFF(1)
      INTEGER  SINGLETOPONOFF(1), SINGLESIDEONOFF(1)
      REAL *8  XSHIFT
      REAL *8  FTOP(8,1), FSIDE(8,1)
      REAL *8  FTOPD(8,1), FSIDED(8,1)
      REAL *8  HL, HR
C-----Local variables
      INTEGER  I, JJ, II, ILEV
      INTEGER  NSIDE
      INTEGER  ISTART, IEND
      REAL *8  TEMP3
      REAL *8  XUSE, PI, PI16, X, XX(8)

      PI = DACOS(-1.0D0)

      PI16 = PI / 16.0D0
      XX(1) = DCOS(15.0D0*PI16) / 2.0D0
      XX(2) = DCOS(13.0D0*PI16) / 2.0D0
      XX(3) = DCOS(11.0D0*PI16) / 2.0D0
      XX(4) = DCOS( 9.0D0*PI16) / 2.0D0
      XX(5) = DCOS( 7.0D0*PI16) / 2.0D0
      XX(6) = DCOS( 5.0D0*PI16) / 2.0D0
      XX(7) = DCOS( 3.0D0*PI16) / 2.0D0
      XX(8) = DCOS( 1.0D0*PI16) / 2.0D0


      DO II = 1, NBOXES
        DOUBLETOPONOFF(II) = 0
        DOUBLESIDEONOFF(II) = 0
        SINGLETOPONOFF(II) = 0
        SINGLESIDEONOFF(II) = 0
        DO JJ = 1, 8
           FTOP(JJ,II) = 0.0D0
           FTOPD(JJ,II) = 0.0D0
           FSIDE(JJ,II) = 0.0D0
           FSIDED(JJ,II) = 0.0D0
        END DO
      END DO

      TEMP3 = .5D0
      NSIDE = 1
      DO ILEV = 0, NLEV
      ISTART = ISTARTLEV(ILEV)
      IEND = ISTART + NBLEVEL(ILEV) - 1
      DO I = ISTART, IEND
      II = IBOXLEV(I)
      IF(ICHILDBOX(1,II) .LT. 0)THEN

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ) / TEMP3
           FSIDE(JJ,II) = HR(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE/2)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           FSIDE(JJ,II) = HR(XUSE)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .LE. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           FSIDE(JJ,II) = HL(X)
         END DO
       ENDIF
       ENDIF

       IF(ICOLBOX(II) .EQ. NSIDE)THEN
       IF(IROWBOX(II) .GT. NSIDE/2)THEN
         SINGLESIDEONOFF(II) = 1
         XSHIFT  = DBLE(IROWBOX(II) - .50D0) / TEMP3 - .50D0
         DO JJ = 1, 8
           X = XSHIFT + XX(JJ)/TEMP3
           XUSE = X - 1.0D0
           FSIDE(JJ,II) = HL(XUSE)
         END DO
       ENDIF
       ENDIF

      ENDIF
      END DO
      TEMP3 = 2.0D0*TEMP3
      NSIDE = 2*NSIDE
      END DO

C     Now scale everything appropriately
C     and account for the sign changes:
      DO II = 1, NBOXES
        DO JJ = 1, 8
           FSIDE(JJ,II) = -FSIDE(JJ,II) / PI
        END DO
      END DO
      RETURN
      END



C**************************************************************************
C--------------------------------------------------------------------------
C
C     Function Module
C
C--------------------------------------------------------------------------
C**************************************************************************

C************************************************************************
C     The following function is just the right right hand side of
C     the Poisson equation.
C************************************************************************
      FUNCTION H2(X,Y) 
      IMPLICIT NONE
C-----Global variables
      REAL *8  H2, X, Y
C-----Local variables
      REAL *8  ALPHA, PI, R, R0, R1, R2

      PI = DACOS(-1.0D0)
      R = DSQRT(X**2 + Y**2)
      R0 = DSQRT((X-.25D0)**2 + (Y-.25D0)**2)
      R1 = DSQRT(X**2 + Y**2)
      R2 = DSQRT((X+.25D0)**2 + (Y-.10D0)**2)
      ALPHA = 200.0D0

ccc      H2 = DEXP(-ALPHA*R**2)

C      Basis functions for level zero below:
ccc      H2 = 1.0D0
ccc      H2 = 2.0D0*X
ccc      H2 = 8.0D0*X**2 - 1.0D0
ccc      H2 = 32.0D0*X**3 - 6.0D0*X
ccc      H2 = 128.0D0*X**4 - 32.0D0*X**2 + 1.0D0
ccc      H2 = 512.0D0*X**5 - 160.0D0*X**3 + 10.0D0*X
ccc      H2 = 2048.0D0*X**6 - 768.0D0*X**4 + 72.0D0*X**2 - 1.0D0
ccc      H2 = 8192.0D0*X**7 - 3584.0D0*X**5 + 448.0D0* X**3 - 14.0D0*X

ccc         H2 =  1.0D0
ccc     1   + 2.0D0* X - Y
ccc     2   +  X**2 + X*Y + Y**2
ccc     3   +  X**3 + 3.0D0* X**2*Y + X*Y**2 + Y**3
ccc     4   +  X**4 + X**3*Y + X**2*Y**2 + X*Y**3 + Y**4
ccc     5   +  X**5 + X**4*Y - X**3*Y**2 + X**2*Y**3 + X*Y**4 + Y**5
ccc     6   +  X**6 + X**5*Y + X**4*Y**2 + X**3*Y**3
ccc     7                         + X**2*Y**4 + X*Y**5 + Y**6
ccc     8   +  X**7 + X**6*Y - X**5*Y**2 + X**4*Y**3 + X**3*Y**4
ccc     9                 - 5.0D0* X**2*Y**5 + X*Y**6 + Y**7

ccc      H2 = 2.0D0
ccc      H2 = X
ccc      H2 = X**2
ccc      H2 = X**3
ccc      H2 = -4.0D0 * PI**2 * X * DCOS(2.0D0*PI*Y)
ccc      H2 = 4.0D0
      H2 = 42.0D0 * (X**5 + Y**5)
ccc      H2 = X**7 + Y**7 + X**5 + Y**5 + X - Y + 1.0D0
      H2 = -8.0D0 * PI**2 * DSIN(2.0D0*PI*X) * DSIN(2.0D0*PI*Y)
      H2 = (4*X*X*ALPHA**2 -2*ALPHA)*DEXP(-ALPHA*R**2) 
      H2 = H2 + (4*Y*Y*ALPHA**2 -2*ALPHA)*DEXP(-ALPHA*R**2) 


      RETURN
      END


      FUNCTION HEXACT(X,Y) 
      IMPLICIT NONE
C-----Global variables
      REAL *8  HEXACT, X, Y
C-----Local variables
      REAL *8  PI
      REAL *8  ALPHA, R, R0, R1, R2

      PI = DACOS(-1.0D0)
      R  = DSQRT(X**2 + Y**2)
      R0 = DSQRT((X-.25D0)**2 + (Y-.25D0)**2)
      R1 = DSQRT(X**2 + Y**2)
      R2 = DSQRT((X+.25D0)**2 + (Y-.10D0)**2)
      ALPHA = 200.0D0

ccc      HEXACT = 8.0D0*X**2 - 1.0D0
ccc      HEXACT = X * DSIN(2.0D0*PI*Y)
ccc      HEXACT = X**7 + Y**7
ccc      HEXACT = X*Y
ccc      HEXACT = X * DCOS(2.0D0*PI*Y)
ccc      HEXACT = (X**2 - .25D0)
ccc      HEXACT = X**5 + Y**5
ccc      HEXACT = X**3 + Y**3
ccc      HEXACT = X**5 * Y**2
ccc      HEXACT = X * Y**3
ccc      HEXACT = X**5 + Y**5
ccc      HEXACT = X**5 * Y
ccc      HEXACT = X**7
ccc      HEXACT = X**2*Y
ccc      HEXACT = X**2 + Y**2
      HEXACT = X**7 + Y**7
ccc      HEXACT = X**2 * Y**2
ccc      HEXACT = X * Y
ccc      HEXACT = DSIN(8.0D0*PI*X) * DSIN(8.0D0*PI*Y)
      HEXACT = DSIN(2.0D0*PI*X) * DSIN(2.0D0*PI*Y)
      HEXACT = DEXP(-ALPHA*R**2) 

      RETURN
      END


      FUNCTION HEXACTX(X,Y)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HEXACTX, X, Y
C-----Local variables
      REAL *8  PI

      PI = DACOS(-1.0D0)

ccc      HEXACTX = 2.0D0 * PI * DCOS(2.0D0*PI*X) * DSIN(2.0D0*PI*Y)
      HEXACTX = 7.0D0*X**6

      RETURN
      END


      FUNCTION HEXACTY(X,Y)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HEXACTY, X, Y
C-----Local variables
      REAL *8  PI

      PI = DACOS(-1.0D0)

ccc      HEXACTY = 2.0D0 * PI * DSIN(2.0D0*PI*X) * DCOS(2.0D0*PI*Y)
      HEXACTY = 7.0D0*Y**6

      RETURN
      END


      FUNCTION HTOP(X)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HTOP, X, PI

      PI = DACOS(-1.0D0)

ccc        HTOP = .25d0 * X**2
ccc        HTOP = DEXP(2.0D0*X)
ccc        HTOP = 7.0D0 * X * .50D0**6
ccc        HTOP = 8.0D0 * X**2 - 1.0D0
ccc        HTOP = 2*X**3 - X**2 + X  - 3.0D0
ccc        HTOP = .125D0 + X**3
ccc        HTOP = 8.0D0 * X**4 - 8.0D0*X**2 + 1.0D0
ccc        HTOP = 0.0D0
ccc        HTOP = 7.0D0 * .50D0**6 * X
ccc        HTOP = 0.5D0
ccc        HTOP = X**2 + .25D0
        HTOP = X**7 + .50D0**7
ccc        HTOP = X
ccc        HTOP = .50D0*X
ccc        HTOP = 0.5D0
      RETURN
      END

      FUNCTION HBOTTOM(X)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HBOTTOM, X, PI

      PI = DACOS(-1.0D0)

ccc        HBOTTOM = .25d0 * X**2
ccc        HBOTTOM = DEXP(2.0D0 * X)
ccc        HBOTTOM = -7.0D0 * X * .50D0**6
ccc        HBOTTOM = 2*X**3 - X**2 + X  - 3.0D0
ccc        HBOTTOM = 2*X**3 - X
ccc        HBOTTOM = X**3 - .75D0 * X
ccc        HBOTTOM = 8.0D0 * X**2 - 1.0D0
ccc        HBOTTOM = -X
ccc        HBOTTOM = -.75D0 * X
ccc        HBOTTOM = -1.0D0
        HBOTTOM = X**7 - .50D0**7
ccc        HBOTTOM = 8.0D0 * X**4 - 8.0D0*X**2 + 1.0D0
ccc        HBOTTOM = -.50D0
ccc        HBOTTOM = -.50D0*X
ccc        HBOTTOM = -0.5D0
      RETURN
      END


      FUNCTION HRIGHT(Y)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HRIGHT, Y, PI

      PI = DACOS(-1.0D0)

ccc        HRIGHT = 0.0D0
        HRIGHT = .50D0**7 + Y**7
ccc        HRIGHT = DEXP(1.0D0)
ccc        HRIGHT = 1.0D0 / 2.0D0**7
ccc        HRIGHT = Y**7
ccc        HRIGHT = .50D0
ccc        HRIGHT = -2.0D0 * PI * DSINH(2.0D0*PI*Y)
ccc        HRIGHT = 2.0D0 * PI * DSIN(2.0D0*PI*Y) * DCOSH(PI)
ccc        HRIGHT = 4d0* (.5d0)**3
ccc        HRIGHT = .0625D0 + Y**4
ccc        HRIGHT = .50D0 * Y**3
ccc        HRIGHT = .015625D0
ccc        HRIGHT = Y**2 + .25D0
ccc        HRIGHT = .50D0 * Y
ccc        HRIGHT = DCOS(2.0D0*PI*Y)
ccc        HRIGHT = .50D0*Y
      RETURN
      END


      FUNCTION HLEFT(Y)
      IMPLICIT NONE
C-----Global variables
      REAL *8  HLEFT, Y, PI

      PI = DACOS(-1.0D0)

ccc        HLEFT = 0.0D0
        HLEFT = -.50D0**7 + Y**7
ccc        HLEFT = DEXP(-1.0D0)
ccc        HLEFT = 1.0D0
ccc        HLEFT = Y**7
ccc        HLEFT = -SIN(2.0D0*PI*Y)
ccc        HLEFT = 0.0D0
ccc        HLEFT = Y**5 - .03125D0
ccc        HLEFT = -1.0D0 / 2.0D0**7
ccc        HLEFT = .015625D0
ccc        HLEFT = -.03125D0 * Y
ccc        HLEFT = -.50D0 * Y**3
ccc        HLEFT = Y**2 + .25D0
ccc        HLEFT = -.50D0 * Y
ccc        HLEFT = -DCOS(2.0D0*PI*Y)
ccc        HLEFT = -.50D0*Y
      RETURN
      END


      FUNCTION CHEB(X,I)
      IMPLICIT NONE  
      INTEGER  I
      REAL *8  CHEB, X

        IF(I .EQ. 0)THEN
          CHEB = 1.0D0
        ELSEIF(I .EQ. 1)THEN
           CHEB = X
        ELSEIF(I .EQ. 2)THEN
           CHEB = 2.0D0*X**2 - 1.0D0
        ELSEIF(I .EQ. 3)THEN
           CHEB = 4.0D0*X**3 - 3.0D0*X
        ELSEIF(I .EQ. 4)THEN
           CHEB = 8.0D0*X**4 - 8.0D0*X**2 + 1.0D0
        ELSEIF(I .EQ. 5)THEN
           CHEB = 16.0D0*X**5 - 20.0D0*X**3 + 5.0D0*X
        ELSEIF(I .EQ. 6)THEN
           CHEB = 32.0D0*X**6 - 48.0D0*X**4 + 18.0D0*X**2 - 1.0D0
        ELSEIF(I .EQ. 7)THEN
           CHEB = 64.0D0*X**7 - 112.0D0*X**5 + 56.0D0*X**3 - 7.0D0*X
        ENDIF

      RETURN
      END


C**************************************************************************
C--------------------------------------------------------------------------
C
C     Initialization Module
C
C--------------------------------------------------------------------------
C**************************************************************************

C**************************************************************************
C     The following subroutine precomputes all of the maps that depend on 
C     the specified precision.
C
C     INPUT:
C
C     IPREC is a parameter that determines the precision (and in turn the
C           number of terms in the multipole and exponential expansions)
C
C     MAP is a real array that is blank
C
C     LENMAPS is the length of the array MAP
C
C     OUTPUT:
C   
C     MAP contains all of the maps that depend on precision.
C         When MAP is passed to the FMMSTART6 routine, it
C         is divided up in a similar manner.
C
C**************************************************************************
      SUBROUTINE COMPUTEMAPS8(MAP, IPREC, LENMAPS)
      IMPLICIT NONE
C-----Global variables
      INTEGER  IPREC, LENMAPS
      REAL *8  MAP(1)
C-----Local variables
      INTEGER  NTERMS, NNODES, ITOT
      INTEGER  LVOLMAPS, LSIDEMAPS, LWINT
      INTEGER  MMAPSOUTH, MMAPNORTH, MMAPEAST, MMAPWEST
      INTEGER  MEDBLETOP, MEDBLESIDE, MESNGLETOP, MESNGLESIDE
      INTEGER  MWDBLETOP, MWDBLESIDE, MWSNGLETOP, MWSNGLESIDE
      INTEGER  MNDBLETOP, MNDBLESIDE, MNSNGLETOP, MNSNGLESIDE
      INTEGER  MSDBLETOP, MSDBLESIDE, MSSNGLETOP, MSSNGLESIDE
      INTEGER  MWINT

      IF(IPREC .EQ. 0)THEN
        NTERMS = 8
        NNODES = 8
      ELSEIF(IPREC .EQ. 1)THEN
        NTERMS = 16
        NNODES = 16
      ELSEIF(IPREC .EQ. 2)THEN
        NTERMS = 24
        NNODES = 24
      ELSEIF(IPREC .EQ. 3)THEN
        NTERMS = 40
        NNODES = 40
      ENDIF

      LVOLMAPS = 2*(NNODES + 1) * 36
      LSIDEMAPS = 2*(NNODES + 1) * 8
      LWINT = 2*(NTERMS + 1) * 8 * 8

      MMAPSOUTH   =  1
      MMAPNORTH   =  MMAPSOUTH   + LVOLMAPS
      MMAPEAST    =  MMAPNORTH   + LVOLMAPS
      MMAPWEST    =  MMAPEAST    + LVOLMAPS
      MWINT       =  MMAPWEST    + LVOLMAPS
      MEDBLETOP   =  MWINT       + LWINT
      MEDBLESIDE  =  MEDBLETOP   + LSIDEMAPS
      MESNGLETOP  =  MEDBLESIDE  + LSIDEMAPS
      MESNGLESIDE =  MESNGLETOP  + LSIDEMAPS
      MWDBLETOP   =  MESNGLESIDE + LSIDEMAPS
      MWDBLESIDE  =  MWDBLETOP   + LSIDEMAPS
      MWSNGLETOP  =  MWDBLESIDE  + LSIDEMAPS
      MWSNGLESIDE =  MWSNGLETOP  + LSIDEMAPS
      MNDBLETOP   =  MWSNGLESIDE + LSIDEMAPS
      MNDBLESIDE  =  MNDBLETOP   + LSIDEMAPS
      MNSNGLETOP  =  MNDBLESIDE  + LSIDEMAPS
      MNSNGLESIDE =  MNSNGLETOP  + LSIDEMAPS
      MSDBLETOP   =  MNSNGLESIDE + LSIDEMAPS
      MSDBLESIDE  =  MSDBLETOP   + LSIDEMAPS
      MSSNGLETOP  =  MSDBLESIDE  + LSIDEMAPS
      MSSNGLESIDE =  MSSNGLETOP  + LSIDEMAPS

      ITOT = MSSNGLESIDE + LSIDEMAPS

      IF(ITOT .GT. LENMAPS)THEN
        WRITE(*,*)'There is not enough workspace needed to'
        WRITE(*,*)'form the maps.'
        WRITE(*,*)'You have set the length of the map'
        WRITE(*,*)'array to be:',LENMAPS
        WRITE(*,*)'The workspace needed to form maps is:',ITOT
        WRITE(*,*)'I am stopping.'
        STOP
      ENDIF


      CALL GETMAPS8(MAP(MMAPSOUTH), MAP(MMAPNORTH),
     1     MAP(MMAPEAST), MAP(MMAPWEST), MAP(MEDBLETOP),
     2     MAP(MEDBLESIDE), MAP(MESNGLETOP), MAP(MESNGLESIDE),
     3     MAP(MWDBLETOP), MAP(MWDBLESIDE), MAP(MWSNGLETOP),
     4     MAP(MWSNGLESIDE), MAP(MNDBLETOP), MAP(MNDBLESIDE),
     5     MAP(MNSNGLETOP), MAP(MNSNGLESIDE), MAP(MSDBLETOP),
     6     MAP(MSDBLESIDE), MAP(MSSNGLETOP),
     7     MAP(MSSNGLESIDE), MAP(MWINT), NNODES, NTERMS)

      RETURN
      END


C**************************************************************************
C     The following subroutine precomputes the maps that are needed in the
C     big to small far case (layer and volume) as well as the weights that
C     map from the polynomial basis coefficients to the multipole terms.
C
C     INPUT:
C     
C     NTERMS is the number of terms in the multipole expansions
C     
C     NNODES is the number of terms in the exponential expansions
C
C     OUTPUT:
C
C     MAPSOUTH is the map for getting the big to small south exponential
C
C     MAPNORTH is the map for getting the big to small north exponential
C
C     MAPEAST  is the map for getting the big to small east  exponential
C
C     MAPWEST  is the map for getting the big to small west  exponential
C
C     WINT represents the weights that pass from the coefficients to the
C          multipole terms
C
C     EDBLETOP is the east going layer map for the double layer on the top
C              (The other matrices are defined using the same convention)
C
C**************************************************************************
      SUBROUTINE GETMAPS8(MAPSOUTH, MAPNORTH, MAPEAST, MAPWEST,
     1     EDBLETOP, EDBLESIDE, ESNGLETOP, ESNGLESIDE,
     2     WDBLETOP, WDBLESIDE, WSNGLETOP, WSNGLESIDE,
     3     NDBLETOP, NDBLESIDE, NSNGLETOP, NSNGLESIDE,
     5     SDBLETOP, SDBLESIDE, SSNGLETOP, SSNGLESIDE,
     6     WINT, NNODES, NTERMS)
      IMPLICIT NONE
C-----Global variables
      INTEGER   NTERMS, NNODES
      COMPLEX *16  MAPNORTH(0:NNODES,36), MAPSOUTH(0:NNODES,36)
      COMPLEX *16  MAPEAST(0:NNODES,36),  MAPWEST(0:NNODES,36)
      COMPLEX *16  EDBLETOP(0:NNODES,0:7),  EDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  ESNGLETOP(0:NNODES,0:7), ESNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  WDBLETOP(0:NNODES,0:7),  WDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  WSNGLETOP(0:NNODES,0:7), WSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  NDBLETOP(0:NNODES,0:7),  NDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  NSNGLETOP(0:NNODES,0:7), NSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  SDBLETOP(0:NNODES,0:7),  SDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  SSNGLETOP(0:NNODES,0:7), SSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  WINT(0:NTERMS,0:7,0:7)

      CALL GENWEIGHTSMULT8(WINT,NTERMS)

      CALL GENBTOSVOLMAT8(MAPEAST,MAPWEST,MAPSOUTH,MAPNORTH,NNODES)

      CALL GENBTOSLAYMAT8(ESNGLETOP, WSNGLETOP, NSNGLETOP, 
     1       SSNGLETOP, ESNGLESIDE, WSNGLESIDE, NSNGLESIDE, 
     2       SSNGLESIDE, EDBLETOP, WDBLETOP, NDBLETOP,
     3       SDBLETOP, EDBLESIDE, WDBLESIDE, NDBLESIDE,
     4       SDBLESIDE, NNODES)

      RETURN
      END


C**************************************************************************
C     The following subroutine precomputes the weights that pass from 
C     the polynomial coefficients to the multipole terms.
C
C     INPUT:
C     
C     P is the number of terms in the multipole expansions
C
C     OUTPUT:
C
C     W represents the weights that pass from the coefficients to the
C          multipole terms
C
C**************************************************************************
      SUBROUTINE GENWEIGHTSMULT8(W, P)
      IMPLICIT NONE
C-----Global variables
      INTEGER  P
      COMPLEX *16  W(0:P,0:7,0:7)
C-----Local variables
      INTEGER  I, J, K
      INTEGER  KPTS, KIND, M, L
      REAL *8  ALPHA, BETA, ENDPTS(2), B(30)
      REAL *8  T(30), WNODE(30)
      REAL *8  CHEB
      COMPLEX *16  SUM, IMAG
      DATA IMAG/(0.0,1.0)/

      ALPHA = 0.0D0
      BETA = 0.0D0
      KPTS = 0
      ENDPTS(1) = 0.0D0
      ENDPTS(2) = 0.0D0
      KIND = 1

      CALL GAUSSQ(KIND, 30, ALPHA, BETA, KPTS, ENDPTS, B, T, WNODE)

      DO K = 0,P
       DO I = 0,7
         DO J = 0,7-I


          SUM = 0.0D0
          DO M = 1,30
            DO L = 1,30
             SUM = SUM + WNODE(L)*WNODE(M)*CHEB(T(M),J)*CHEB(T(L),I)
     1                  * (T(L) + IMAG*T(M))**K
            END DO
          END DO

         W(K,I,J) = SUM / 4.0D0

        END DO
       END DO
      END DO

      RETURN
      END


C**************************************************************************
C     The following subroutine precomputes the matrices that map from the
C     coefficients values (that approximating the charge densities) to
C     the representative exponential expansions needed in the big to small
C     far case.
C
C     INPUT:
C     
C     NNODES is the number of terms in the exponential expansions
C
C     OUTPUT:
C
C     The matrices are defined as direction, single or double, top or side.
C     So EASTSNGLETOP, would represent the map from the single layer
C     coefficients on the top of the box to the expansion going in the 
C     east direction.
C
C**************************************************************************
      SUBROUTINE GENBTOSLAYMAT8(EASTSNGLETOP, WESTSNGLETOP,
     1    NORTHSNGLETOP, SOUTHSNGLETOP, EASTSNGLESIDE,
     2    WESTSNGLESIDE, NORTHSNGLESIDE, SOUTHSNGLESIDE, 
     3    EASTDBLETOP, WESTDBLETOP, NORTHDBLETOP,
     4    SOUTHDBLETOP, EASTDBLESIDE, WESTDBLESIDE, 
     5    NORTHDBLESIDE, SOUTHDBLESIDE, NNODES)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES
      COMPLEX *16  EASTSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  WESTSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  NORTHSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  SOUTHSNGLETOP(0:NNODES,0:7)
      COMPLEX *16  EASTSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  WESTSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  NORTHSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  SOUTHSNGLESIDE(0:NNODES,0:7)
      COMPLEX *16  EASTDBLETOP(0:NNODES,0:7)
      COMPLEX *16  WESTDBLETOP(0:NNODES,0:7)
      COMPLEX *16  NORTHDBLETOP(0:NNODES,0:7)
      COMPLEX *16  SOUTHDBLETOP(0:NNODES,0:7)
      COMPLEX *16  EASTDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  WESTDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  NORTHDBLESIDE(0:NNODES,0:7)
      COMPLEX *16  SOUTHDBLESIDE(0:NNODES,0:7)
C-----Local variables
      INTEGER  I, J
      REAL *8  WNODES(40), XNODES(40)
      REAL *8  TEMP1, TEMP2, TEMP5
      REAL *8  TX1D
      COMPLEX *16  TEMP3, TEMP4
      COMPLEX *16  TY1D
      COMPLEX *16  TEMPNORTHTOP, TEMPSOUTHTOP
      COMPLEX *16  TEMPEASTTOP, TEMPWESTTOP
      COMPLEX *16  TEMPNORTHSIDE, TEMPSOUTHSIDE
      COMPLEX *16  TEMPEASTSIDE, TEMPWESTSIDE
      COMPLEX *16  ZSOUTHTOP, ZNORTHTOP
      COMPLEX *16  ZWESTTOP, ZEASTTOP
      COMPLEX *16  ZSOUTHSIDE, ZNORTHSIDE
      COMPLEX *16  ZWESTSIDE, ZEASTSIDE
      COMPLEX *16  ZSHIFTBIG
      COMPLEX *16  IMAG
      DATA IMAG/(0.0D0,1.0D0)/

      CALL PWTS4(XNODES,WNODES,NNODES)

      TEMP1 = -2.0D0 / 3.0D0
      TEMP2 = (16.0D0/5.0D0 - 16.0D0/3.0D0 + 2.0D0)
      TEMP5 = (64.0D0/7.0D0 - 96.0D0/5.0D0 + 10.0D0)
      EASTSNGLETOP(0,0) = 2.0D0
      EASTSNGLETOP(0,1) = 0.0D0
      EASTSNGLETOP(0,2) = TEMP1
      EASTSNGLETOP(0,3) = 0.0D0
      EASTSNGLETOP(0,4) = TEMP2
      EASTSNGLETOP(0,5) = 0.0D0
      EASTSNGLETOP(0,6) = TEMP5
      EASTSNGLETOP(0,7) = 0.0D0

      WESTSNGLETOP(0,0) = 2.0D0
      WESTSNGLETOP(0,1) = 0.0D0
      WESTSNGLETOP(0,2) = TEMP1
      WESTSNGLETOP(0,3) = 0.0D0
      WESTSNGLETOP(0,4) = TEMP2
      WESTSNGLETOP(0,5) = 0.0D0
      WESTSNGLETOP(0,6) = TEMP5
      WESTSNGLETOP(0,7) = 0.0D0

      NORTHSNGLETOP(0,0) = 2.0D0
      NORTHSNGLETOP(0,1) = 0.0D0
      NORTHSNGLETOP(0,2) = TEMP1
      NORTHSNGLETOP(0,3) = 0.0D0
      NORTHSNGLETOP(0,4) = TEMP2
      NORTHSNGLETOP(0,5) = 0.0D0
      NORTHSNGLETOP(0,6) = TEMP5
      NORTHSNGLETOP(0,7) = 0.0D0

      SOUTHSNGLETOP(0,0) = 2.0D0
      SOUTHSNGLETOP(0,1) = 0.0D0
      SOUTHSNGLETOP(0,2) = TEMP1
      SOUTHSNGLETOP(0,3) = 0.0D0
      SOUTHSNGLETOP(0,4) = TEMP2
      SOUTHSNGLETOP(0,5) = 0.0D0
      SOUTHSNGLETOP(0,6) = TEMP5
      SOUTHSNGLETOP(0,7) = 0.0D0

      EASTSNGLESIDE(0,0) = 2.0D0
      EASTSNGLESIDE(0,1) = 0.0D0
      EASTSNGLESIDE(0,2) = TEMP1
      EASTSNGLESIDE(0,3) = 0.0D0
      EASTSNGLESIDE(0,4) = TEMP2
      EASTSNGLESIDE(0,5) = 0.0D0
      EASTSNGLESIDE(0,6) = TEMP5
      EASTSNGLESIDE(0,7) = 0.0D0

      WESTSNGLESIDE(0,0) = 2.0D0
      WESTSNGLESIDE(0,1) = 0.0D0
      WESTSNGLESIDE(0,2) = TEMP1
      WESTSNGLESIDE(0,3) = 0.0D0
      WESTSNGLESIDE(0,4) = TEMP2
      WESTSNGLESIDE(0,5) = 0.0D0
      WESTSNGLESIDE(0,6) = TEMP5
      WESTSNGLESIDE(0,7) = 0.0D0

      NORTHSNGLESIDE(0,0) = 2.0D0
      NORTHSNGLESIDE(0,1) = 0.0D0
      NORTHSNGLESIDE(0,2) = TEMP1
      NORTHSNGLESIDE(0,3) = 0.0D0
      NORTHSNGLESIDE(0,4) = TEMP2
      NORTHSNGLESIDE(0,5) = 0.0D0
      NORTHSNGLESIDE(0,6) = TEMP5
      NORTHSNGLESIDE(0,7) = 0.0D0

      SOUTHSNGLESIDE(0,0) = 2.0D0
      SOUTHSNGLESIDE(0,1) = 0.0D0
      SOUTHSNGLESIDE(0,2) = TEMP1
      SOUTHSNGLESIDE(0,3) = 0.0D0
      SOUTHSNGLESIDE(0,4) = TEMP2
      SOUTHSNGLESIDE(0,5) = 0.0D0
      SOUTHSNGLESIDE(0,6) = TEMP5
      SOUTHSNGLESIDE(0,7) = 0.0D0

      EASTDBLETOP(0,0) = 0.0D0
      EASTDBLETOP(0,1) = 0.0D0
      EASTDBLETOP(0,2) = 0.0D0 
      EASTDBLETOP(0,3) = 0.0D0
      EASTDBLETOP(0,4) = 0.0D0
      EASTDBLETOP(0,5) = 0.0D0
      EASTDBLETOP(0,6) = 0.0D0
      EASTDBLETOP(0,7) = 0.0D0

      WESTDBLETOP(0,0) = 0.0D0
      WESTDBLETOP(0,1) = 0.0D0
      WESTDBLETOP(0,2) = 0.0D0
      WESTDBLETOP(0,3) = 0.0D0
      WESTDBLETOP(0,4) = 0.0D0
      WESTDBLETOP(0,5) = 0.0D0
      WESTDBLETOP(0,6) = 0.0D0
      WESTDBLETOP(0,7) = 0.0D0

      NORTHDBLETOP(0,0) = 0.0D0
      NORTHDBLETOP(0,1) = 0.0D0
      NORTHDBLETOP(0,2) = 0.0D0 
      NORTHDBLETOP(0,3) = 0.0D0
      NORTHDBLETOP(0,4) = 0.0D0
      NORTHDBLETOP(0,5) = 0.0D0
      NORTHDBLETOP(0,6) = 0.0D0
      NORTHDBLETOP(0,7) = 0.0D0

      SOUTHDBLETOP(0,0) = 0.0D0
      SOUTHDBLETOP(0,1) = 0.0D0
      SOUTHDBLETOP(0,2) = 0.0D0
      SOUTHDBLETOP(0,3) = 0.0D0
      SOUTHDBLETOP(0,4) = 0.0D0
      SOUTHDBLETOP(0,5) = 0.0D0
      SOUTHDBLETOP(0,6) = 0.0D0
      SOUTHDBLETOP(0,7) = 0.0D0

      EASTDBLESIDE(0,0) = 0.0D0
      EASTDBLESIDE(0,1) = 0.0D0
      EASTDBLESIDE(0,2) = 0.0D0
      EASTDBLESIDE(0,3) = 0.0D0
      EASTDBLESIDE(0,4) = 0.0D0
      EASTDBLESIDE(0,5) = 0.0D0
      EASTDBLESIDE(0,6) = 0.0D0
      EASTDBLESIDE(0,7) = 0.0D0

      WESTDBLESIDE(0,0) = 0.0D0
      WESTDBLESIDE(0,1) = 0.0D0
      WESTDBLESIDE(0,2) = 0.0D0
      WESTDBLESIDE(0,3) = 0.0D0
      WESTDBLESIDE(0,4) = 0.0D0
      WESTDBLESIDE(0,5) = 0.0D0
      WESTDBLESIDE(0,6) = 0.0D0
      WESTDBLESIDE(0,7) = 0.0D0

      NORTHDBLESIDE(0,0) = 0.0D0
      NORTHDBLESIDE(0,1) = 0.0D0
      NORTHDBLESIDE(0,2) = 0.0D0
      NORTHDBLESIDE(0,3) = 0.0D0
      NORTHDBLESIDE(0,4) = 0.0D0
      NORTHDBLESIDE(0,5) = 0.0D0
      NORTHDBLESIDE(0,6) = 0.0D0
      NORTHDBLESIDE(0,7) = 0.0D0

      SOUTHDBLESIDE(0,0) = 0.0D0
      SOUTHDBLESIDE(0,1) = 0.0D0
      SOUTHDBLESIDE(0,2) = 0.0D0
      SOUTHDBLESIDE(0,3) = 0.0D0
      SOUTHDBLESIDE(0,4) = 0.0D0
      SOUTHDBLESIDE(0,5) = 0.0D0
      SOUTHDBLESIDE(0,6) = 0.0D0
      SOUTHDBLESIDE(0,7) = 0.0D0

      ZSHIFTBIG =  -0.5D0 - IMAG * 1.5D0
      ZEASTTOP  =      -ZSHIFTBIG
      ZWESTTOP  =       ZSHIFTBIG
      ZNORTHTOP =  IMAG*ZSHIFTBIG
      ZSOUTHTOP = -IMAG*ZSHIFTBIG

      ZSHIFTBIG =  -1.5D0 - IMAG * 0.5D0
      ZEASTSIDE  =      -ZSHIFTBIG
      ZWESTSIDE  =       ZSHIFTBIG
      ZNORTHSIDE =  IMAG*ZSHIFTBIG
      ZSOUTHSIDE = -IMAG*ZSHIFTBIG

      DO I = 1, NNODES
      TEMP1 =  WNODES(I) / XNODES(I)
      TEMP2 =  WNODES(I)
      TEMP3 =  TEMP2 * IMAG

      TEMPEASTTOP  = CDEXP(ZEASTTOP*XNODES(I))
      TEMPWESTTOP  = CDEXP(ZWESTTOP*XNODES(I))
      TEMPNORTHTOP = CDEXP(ZNORTHTOP*XNODES(I))
      TEMPSOUTHTOP = CDEXP(ZSOUTHTOP*XNODES(I))

      TEMPEASTSIDE  = CDEXP(ZEASTSIDE*XNODES(I))
      TEMPWESTSIDE  = CDEXP(ZWESTSIDE*XNODES(I))
      TEMPNORTHSIDE = CDEXP(ZNORTHSIDE*XNODES(I))
      TEMPSOUTHSIDE = CDEXP(ZSOUTHSIDE*XNODES(I))

       DO J = 0, 7 

         CALL GETINTSX(XNODES(I), TX1D, J)

         TEMP4 = TX1D * TEMPEASTTOP
         EASTSNGLETOP(I,J) =  TEMP4 * TEMP1
         EASTDBLETOP(I,J)  = -TEMP4 * TEMP3
             

         TEMP4 = TX1D * TEMPNORTHSIDE
         NORTHSNGLESIDE(I,J) = TEMP4 * TEMP1
         NORTHDBLESIDE(I,J)  = TEMP4 * TEMP3



         CALL GETINTSX(-XNODES(I), TX1D, J)

         TEMP4 = TX1D * TEMPWESTTOP
         WESTSNGLETOP(I,J) = TEMP4 * TEMP1
         WESTDBLETOP(I,J)  = TEMP4 * TEMP3

         TEMP4 = TX1D * TEMPSOUTHSIDE
         SOUTHSNGLESIDE(I,J) =  TEMP4 * TEMP1
         SOUTHDBLESIDE(I,J)  = -TEMP4 * TEMP3



         CALL GETINTSY(XNODES(I), TY1D, J)

         TEMP4 = TY1D * TEMPSOUTHTOP
         SOUTHSNGLETOP(I,J) =  TEMP4 * TEMP1
         SOUTHDBLETOP(I,J)  =  TEMP4 * TEMP2

         TEMP4 = TY1D * TEMPEASTSIDE
         EASTSNGLESIDE(I,J) =  TEMP4 * TEMP1
         EASTDBLESIDE(I,J)  = -TEMP4 * TEMP2



         CALL GETINTSY(-XNODES(I), TY1D, J)

         TEMP4 = TY1D * TEMPNORTHTOP
         NORTHSNGLETOP(I,J) =  TEMP4 * TEMP1
         NORTHDBLETOP(I,J)  = -TEMP4 * TEMP2

         TEMP4 = TY1D * TEMPWESTSIDE
         WESTSNGLESIDE(I,J) =  TEMP4 * TEMP1
         WESTDBLESIDE(I,J)  =  TEMP4 * TEMP2

        END DO
      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine precomputes the matrices that map from the
C     coefficients values, representing the right hand side of the poisson
C     equation to the representative exponential expansions needed in the
C     big to small far case.
C
C     INPUT:
C
C     NNODES is the number of terms in the exponential expansions
C
C     OUTPUT:
C
C     MAPEAST  is the map for east  going expansions
C
C     MAPWEST  is the map for west  going expansions
C
C     MAPSOUTH is the map for south going expansions
C
C     MAPNORTH is the map for north going expansions
C
C**************************************************************************
      SUBROUTINE GENBTOSVOLMAT8(MAPEAST, MAPWEST, 
     1                MAPSOUTH, MAPNORTH, NNODES)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES
      COMPLEX *16  MAPEAST(0:NNODES,36)
      COMPLEX *16  MAPWEST(0:NNODES,36)
      COMPLEX *16  MAPNORTH(0:NNODES,36)
      COMPLEX *16  MAPSOUTH(0:NNODES,36)
C-----Local variables
      INTEGER  I
      REAL *8  PI2
      REAL *8  WNODES(40), XNODES(40)
      REAL *8  TXINT0, TXINT1
      REAL *8  TXINT2, TXINT3
      REAL *8  TXINT4, TXINT5
      REAL *8  TXINT6, TXINT7
      COMPLEX *16  TYINT0, TYINT1
      COMPLEX *16  TYINT2, TYINT3
      COMPLEX *16  TYINT4, TYINT5
      COMPLEX *16  TYINT6, TYINT7
      COMPLEX *16  ZEAST, ZWEST, ZNORTH, ZSOUTH
      COMPLEX *16  TEMP1, TEMP2
      COMPLEX *16  TEMP3, TEMP4
      COMPLEX *16  TEMP5, TEMP6
      COMPLEX *16  IMAG 
      DATA IMAG/(0.0D0,1.0D0)/

      PI2 = 2.0D0*DACOS(-1.0D0)

      CALL PWTS4(XNODES,WNODES,NNODES)
  
      ZEAST  =  0.5D0 + IMAG * 0.5D0
      ZWEST  = -0.5D0 - IMAG * 0.5D0
      ZNORTH =  0.5D0 - IMAG * 0.5D0
      ZSOUTH = -0.5D0 + IMAG * 0.5D0

      TEMP1 =  1.0D0 / PI2
      TEMP2 = -1.0D0 / (3.0D0*PI2)
      TEMP3 = -1.0D0 / (15.0D0*PI2)
      TEMP4 =  1.0D0 / (9.0D0*PI2)
      TEMP5 = -1.0D0 / (35.0D0*PI2)
      TEMP6 =  1.0D0 / (45.0D0*PI2)
      MAPEAST(0,1) = TEMP1
      MAPEAST(0,2) = 0.0D0
      MAPEAST(0,3) = 0.0D0
      MAPEAST(0,4) = TEMP2
      MAPEAST(0,5) = 0.0D0
      MAPEAST(0,6) = TEMP2
      MAPEAST(0,7) = 0.0D0
      MAPEAST(0,8) = 0.0D0
      MAPEAST(0,9) = 0.0D0
      MAPEAST(0,10)= 0.0D0
      MAPEAST(0,11)= TEMP3
      MAPEAST(0,12)= 0.0D0
      MAPEAST(0,13)= TEMP4
      MAPEAST(0,14)= 0.0D0
      MAPEAST(0,15)= TEMP3
      MAPEAST(0,16)= 0.0D0
      MAPEAST(0,17)= 0.0D0
      MAPEAST(0,18)= 0.0D0
      MAPEAST(0,19)= 0.0D0
      MAPEAST(0,20)= 0.0D0
      MAPEAST(0,21)= 0.0D0
      MAPEAST(0,22)= TEMP5
      MAPEAST(0,23)= 0.0D0
      MAPEAST(0,24)= TEMP6
      MAPEAST(0,25)= 0.0D0
      MAPEAST(0,26)= TEMP6
      MAPEAST(0,27)= 0.0D0
      MAPEAST(0,28)= TEMP5
      MAPEAST(0,29)= 0.0D0
      MAPEAST(0,30)= 0.0D0
      MAPEAST(0,31)= 0.0D0
      MAPEAST(0,32)= 0.0D0
      MAPEAST(0,33)= 0.0D0
      MAPEAST(0,34)= 0.0D0
      MAPEAST(0,35)= 0.0D0
      MAPEAST(0,36)= 0.0D0

      MAPWEST(0,1) = TEMP1
      MAPWEST(0,2) = 0.0D0
      MAPWEST(0,3) = 0.0D0
      MAPWEST(0,4) = TEMP2
      MAPWEST(0,5) = 0.0D0
      MAPWEST(0,6) = TEMP2
      MAPWEST(0,7) = 0.0D0
      MAPWEST(0,8) = 0.0D0
      MAPWEST(0,9) = 0.0D0
      MAPWEST(0,10)= 0.0D0
      MAPWEST(0,11)= TEMP3
      MAPWEST(0,12)= 0.0D0
      MAPWEST(0,13)= TEMP4
      MAPWEST(0,14)= 0.0D0
      MAPWEST(0,15)= TEMP3
      MAPWEST(0,16)= 0.0D0
      MAPWEST(0,17)= 0.0D0
      MAPWEST(0,18)= 0.0D0
      MAPWEST(0,19)= 0.0D0
      MAPWEST(0,20)= 0.0D0
      MAPWEST(0,21)= 0.0D0
      MAPWEST(0,22)= TEMP5
      MAPWEST(0,23)= 0.0D0
      MAPWEST(0,24)= TEMP6
      MAPWEST(0,25)= 0.0D0
      MAPWEST(0,26)= TEMP6
      MAPWEST(0,27)= 0.0D0
      MAPWEST(0,28)= TEMP5
      MAPWEST(0,29)= 0.0D0
      MAPWEST(0,30)= 0.0D0
      MAPWEST(0,31)= 0.0D0
      MAPWEST(0,32)= 0.0D0
      MAPWEST(0,33)= 0.0D0
      MAPWEST(0,34)= 0.0D0
      MAPWEST(0,35)= 0.0D0
      MAPWEST(0,36)= 0.0D0

      MAPNORTH(0,1) = TEMP1
      MAPNORTH(0,2) = 0.0D0
      MAPNORTH(0,3) = 0.0D0
      MAPNORTH(0,4) = TEMP2
      MAPNORTH(0,5) = 0.0D0
      MAPNORTH(0,6) = TEMP2
      MAPNORTH(0,7) = 0.0D0
      MAPNORTH(0,8) = 0.0D0
      MAPNORTH(0,9) = 0.0D0
      MAPNORTH(0,10)= 0.0D0
      MAPNORTH(0,11)= TEMP3
      MAPNORTH(0,12)= 0.0D0
      MAPNORTH(0,13)= TEMP4
      MAPNORTH(0,14)= 0.0D0
      MAPNORTH(0,15)= TEMP3
      MAPNORTH(0,16)= 0.0D0
      MAPNORTH(0,17)= 0.0D0
      MAPNORTH(0,18)= 0.0D0
      MAPNORTH(0,19)= 0.0D0
      MAPNORTH(0,20)= 0.0D0
      MAPNORTH(0,21)= 0.0D0
      MAPNORTH(0,22)= TEMP5
      MAPNORTH(0,23)= 0.0D0
      MAPNORTH(0,24)= TEMP6
      MAPNORTH(0,25)= 0.0D0
      MAPNORTH(0,26)= TEMP6
      MAPNORTH(0,27)= 0.0D0
      MAPNORTH(0,28)= TEMP5
      MAPNORTH(0,29)= 0.0D0
      MAPNORTH(0,30)= 0.0D0
      MAPNORTH(0,31)= 0.0D0
      MAPNORTH(0,32)= 0.0D0
      MAPNORTH(0,33)= 0.0D0
      MAPNORTH(0,34)= 0.0D0
      MAPNORTH(0,35)= 0.0D0
      MAPNORTH(0,36)= 0.0D0


      MAPSOUTH(0,1) = TEMP1
      MAPSOUTH(0,2) = 0.0D0
      MAPSOUTH(0,3) = 0.0D0
      MAPSOUTH(0,4) = TEMP2
      MAPSOUTH(0,5) = 0.0D0
      MAPSOUTH(0,6) = TEMP2
      MAPSOUTH(0,7) = 0.0D0
      MAPSOUTH(0,8) = 0.0D0
      MAPSOUTH(0,9) = 0.0D0
      MAPSOUTH(0,10)= 0.0D0
      MAPSOUTH(0,11)= TEMP3
      MAPSOUTH(0,12)= 0.0D0
      MAPSOUTH(0,13)= TEMP4
      MAPSOUTH(0,14)= 0.0D0
      MAPSOUTH(0,15)= TEMP3
      MAPSOUTH(0,16)= 0.0D0
      MAPSOUTH(0,17)= 0.0D0
      MAPSOUTH(0,18)= 0.0D0
      MAPSOUTH(0,19)= 0.0D0
      MAPSOUTH(0,20)= 0.0D0
      MAPSOUTH(0,21)= 0.0D0
      MAPSOUTH(0,22)= TEMP5
      MAPSOUTH(0,23)= 0.0D0
      MAPSOUTH(0,24)= TEMP6
      MAPSOUTH(0,25)= 0.0D0
      MAPSOUTH(0,26)= TEMP6
      MAPSOUTH(0,27)= 0.0D0
      MAPSOUTH(0,28)= TEMP5
      MAPSOUTH(0,29)= 0.0D0
      MAPSOUTH(0,30)= 0.0D0
      MAPSOUTH(0,31)= 0.0D0
      MAPSOUTH(0,32)= 0.0D0
      MAPSOUTH(0,33)= 0.0D0
      MAPSOUTH(0,34)= 0.0D0
      MAPSOUTH(0,35)= 0.0D0
      MAPSOUTH(0,36)= 0.0D0



      DO I = 1, NNODES
        TEMP1 = WNODES(I) / (4.0D0*PI2 * XNODES(I)) 

         CALL GETINTSX(XNODES(I), TXINT0, 0)
         CALL GETINTSX(XNODES(I), TXINT1, 1)
         CALL GETINTSX(XNODES(I), TXINT2, 2)
         CALL GETINTSX(XNODES(I), TXINT3, 3)
         CALL GETINTSX(XNODES(I), TXINT4, 4)
         CALL GETINTSX(XNODES(I), TXINT5, 5)
         CALL GETINTSX(XNODES(I), TXINT6, 6)
         CALL GETINTSX(XNODES(I), TXINT7, 7)

         CALL GETINTSY(XNODES(I), TYINT0, 0)
         CALL GETINTSY(XNODES(I), TYINT1, 1)
         CALL GETINTSY(XNODES(I), TYINT2, 2)
         CALL GETINTSY(XNODES(I), TYINT3, 3)
         CALL GETINTSY(XNODES(I), TYINT4, 4)
         CALL GETINTSY(XNODES(I), TYINT5, 5)
         CALL GETINTSY(XNODES(I), TYINT6, 6)
         CALL GETINTSY(XNODES(I), TYINT7, 7)

        TEMP2 = TEMP1 * CDEXP(ZEAST*XNODES(I)) 

         MAPEAST(I,1) = TXINT0 * TYINT0 * TEMP2
         MAPEAST(I,2) = TXINT1 * TYINT0 * TEMP2
         MAPEAST(I,3) = TXINT0 * TYINT1 * TEMP2
         MAPEAST(I,4) = TXINT2 * TYINT0 * TEMP2
         MAPEAST(I,5) = TXINT1 * TYINT1 * TEMP2
         MAPEAST(I,6) = TXINT0 * TYINT2 * TEMP2
         MAPEAST(I,7) = TXINT3 * TYINT0 * TEMP2
         MAPEAST(I,8) = TXINT2 * TYINT1 * TEMP2
         MAPEAST(I,9) = TXINT1 * TYINT2 * TEMP2
         MAPEAST(I,10)= TXINT0 * TYINT3 * TEMP2
         MAPEAST(I,11)= TXINT4 * TYINT0 * TEMP2
         MAPEAST(I,12)= TXINT3 * TYINT1 * TEMP2
         MAPEAST(I,13)= TXINT2 * TYINT2 * TEMP2
         MAPEAST(I,14)= TXINT1 * TYINT3 * TEMP2
         MAPEAST(I,15)= TXINT0 * TYINT4 * TEMP2
         MAPEAST(I,16)= TXINT5 * TYINT0 * TEMP2
         MAPEAST(I,17)= TXINT4 * TYINT1 * TEMP2
         MAPEAST(I,18)= TXINT3 * TYINT2 * TEMP2
         MAPEAST(I,19)= TXINT2 * TYINT3 * TEMP2
         MAPEAST(I,20)= TXINT1 * TYINT4 * TEMP2
         MAPEAST(I,21)= TXINT0 * TYINT5 * TEMP2
         MAPEAST(I,22)= TXINT6 * TYINT0 * TEMP2
         MAPEAST(I,23)= TXINT5 * TYINT1 * TEMP2
         MAPEAST(I,24)= TXINT4 * TYINT2 * TEMP2
         MAPEAST(I,25)= TXINT3 * TYINT3 * TEMP2
         MAPEAST(I,26)= TXINT2 * TYINT4 * TEMP2
         MAPEAST(I,27)= TXINT1 * TYINT5 * TEMP2
         MAPEAST(I,28)= TXINT0 * TYINT6 * TEMP2
         MAPEAST(I,29)= TXINT7 * TYINT0 * TEMP2
         MAPEAST(I,30)= TXINT6 * TYINT1 * TEMP2
         MAPEAST(I,31)= TXINT5 * TYINT2 * TEMP2
         MAPEAST(I,32)= TXINT4 * TYINT3 * TEMP2
         MAPEAST(I,33)= TXINT3 * TYINT4 * TEMP2
         MAPEAST(I,34)= TXINT2 * TYINT5 * TEMP2
         MAPEAST(I,35)= TXINT1 * TYINT6 * TEMP2
         MAPEAST(I,36)= TXINT0 * TYINT7 * TEMP2


        TEMP2 = TEMP1 * CDEXP(ZWEST*XNODES(I))

         MAPWEST(I,1) =  TXINT0 * TYINT0 * TEMP2
         MAPWEST(I,2) = -TXINT1 * TYINT0 * TEMP2
         MAPWEST(I,3) = -TXINT0 * TYINT1 * TEMP2
         MAPWEST(I,4) =  TXINT2 * TYINT0 * TEMP2
         MAPWEST(I,5) =  TXINT1 * TYINT1 * TEMP2
         MAPWEST(I,6) =  TXINT0 * TYINT2 * TEMP2
         MAPWEST(I,7) = -TXINT3 * TYINT0 * TEMP2
         MAPWEST(I,8) = -TXINT2 * TYINT1 * TEMP2
         MAPWEST(I,9) = -TXINT1 * TYINT2 * TEMP2
         MAPWEST(I,10)= -TXINT0 * TYINT3 * TEMP2
         MAPWEST(I,11)=  TXINT4 * TYINT0 * TEMP2
         MAPWEST(I,12)=  TXINT3 * TYINT1 * TEMP2
         MAPWEST(I,13)=  TXINT2 * TYINT2 * TEMP2
         MAPWEST(I,14)=  TXINT1 * TYINT3 * TEMP2
         MAPWEST(I,15)=  TXINT0 * TYINT4 * TEMP2
         MAPWEST(I,16)= -TXINT5 * TYINT0 * TEMP2
         MAPWEST(I,17)= -TXINT4 * TYINT1 * TEMP2
         MAPWEST(I,18)= -TXINT3 * TYINT2 * TEMP2
         MAPWEST(I,19)= -TXINT2 * TYINT3 * TEMP2
         MAPWEST(I,20)= -TXINT1 * TYINT4 * TEMP2
         MAPWEST(I,21)= -TXINT0 * TYINT5 * TEMP2
         MAPWEST(I,22)=  TXINT6 * TYINT0 * TEMP2
         MAPWEST(I,23)=  TXINT5 * TYINT1 * TEMP2
         MAPWEST(I,24)=  TXINT4 * TYINT2 * TEMP2
         MAPWEST(I,25)=  TXINT3 * TYINT3 * TEMP2
         MAPWEST(I,26)=  TXINT2 * TYINT4 * TEMP2
         MAPWEST(I,27)=  TXINT1 * TYINT5 * TEMP2
         MAPWEST(I,28)=  TXINT0 * TYINT6 * TEMP2
         MAPWEST(I,29)= -TXINT7 * TYINT0 * TEMP2
         MAPWEST(I,30)= -TXINT6 * TYINT1 * TEMP2
         MAPWEST(I,31)= -TXINT5 * TYINT2 * TEMP2
         MAPWEST(I,32)= -TXINT4 * TYINT3 * TEMP2
         MAPWEST(I,33)= -TXINT3 * TYINT4 * TEMP2
         MAPWEST(I,34)= -TXINT2 * TYINT5 * TEMP2
         MAPWEST(I,35)= -TXINT1 * TYINT6 * TEMP2
         MAPWEST(I,36)= -TXINT0 * TYINT7 * TEMP2


        TEMP2 = TEMP1 * CDEXP(ZNORTH*XNODES(I)) 

         MAPNORTH(I,1) =  TXINT0 * TYINT0 * TEMP2
         MAPNORTH(I,2) = -TXINT0 * TYINT1 * TEMP2
         MAPNORTH(I,3) =  TXINT1 * TYINT0 * TEMP2
         MAPNORTH(I,4) =  TXINT0 * TYINT2 * TEMP2
         MAPNORTH(I,5) = -TXINT1 * TYINT1 * TEMP2
         MAPNORTH(I,6) =  TXINT2 * TYINT0 * TEMP2
         MAPNORTH(I,7) = -TXINT0 * TYINT3 * TEMP2
         MAPNORTH(I,8) =  TXINT1 * TYINT2 * TEMP2
         MAPNORTH(I,9) = -TXINT2 * TYINT1 * TEMP2
         MAPNORTH(I,10)=  TXINT3 * TYINT0 * TEMP2
         MAPNORTH(I,11)=  TXINT0 * TYINT4 * TEMP2
         MAPNORTH(I,12)= -TXINT1 * TYINT3 * TEMP2
         MAPNORTH(I,13)=  TXINT2 * TYINT2 * TEMP2
         MAPNORTH(I,14)= -TXINT3 * TYINT1 * TEMP2
         MAPNORTH(I,15)=  TXINT4 * TYINT0 * TEMP2
         MAPNORTH(I,16)= -TXINT0 * TYINT5 * TEMP2
         MAPNORTH(I,17)=  TXINT1 * TYINT4 * TEMP2
         MAPNORTH(I,18)= -TXINT2 * TYINT3 * TEMP2
         MAPNORTH(I,19)=  TXINT3 * TYINT2 * TEMP2
         MAPNORTH(I,20)= -TXINT4 * TYINT1 * TEMP2
         MAPNORTH(I,21)=  TXINT5 * TYINT0 * TEMP2
         MAPNORTH(I,22)=  TXINT0 * TYINT6 * TEMP2
         MAPNORTH(I,23)= -TXINT1 * TYINT5 * TEMP2
         MAPNORTH(I,24)=  TXINT2 * TYINT4 * TEMP2
         MAPNORTH(I,25)= -TXINT3 * TYINT3 * TEMP2
         MAPNORTH(I,26)=  TXINT4 * TYINT2 * TEMP2
         MAPNORTH(I,27)= -TXINT5 * TYINT1 * TEMP2
         MAPNORTH(I,28)=  TXINT6 * TYINT0 * TEMP2
         MAPNORTH(I,29)= -TXINT0 * TYINT7 * TEMP2
         MAPNORTH(I,30)=  TXINT1 * TYINT6 * TEMP2
         MAPNORTH(I,31)= -TXINT2 * TYINT5 * TEMP2
         MAPNORTH(I,32)=  TXINT3 * TYINT4 * TEMP2
         MAPNORTH(I,33)= -TXINT4 * TYINT3 * TEMP2
         MAPNORTH(I,34)=  TXINT5 * TYINT2 * TEMP2
         MAPNORTH(I,35)= -TXINT6 * TYINT1 * TEMP2
         MAPNORTH(I,36)=  TXINT7 * TYINT0 * TEMP2


        TEMP2 = TEMP1 * CDEXP(ZSOUTH*XNODES(I))

         MAPSOUTH(I,1) =  TXINT0 * TYINT0 * TEMP2
         MAPSOUTH(I,2) =  TXINT0 * TYINT1 * TEMP2
         MAPSOUTH(I,3) = -TXINT1 * TYINT0 * TEMP2
         MAPSOUTH(I,4) =  TXINT0 * TYINT2 * TEMP2
         MAPSOUTH(I,5) = -TXINT1 * TYINT1 * TEMP2
         MAPSOUTH(I,6) =  TXINT2 * TYINT0 * TEMP2
         MAPSOUTH(I,7) =  TXINT0 * TYINT3 * TEMP2
         MAPSOUTH(I,8) = -TXINT1 * TYINT2 * TEMP2
         MAPSOUTH(I,9) =  TXINT2 * TYINT1 * TEMP2
         MAPSOUTH(I,10)= -TXINT3 * TYINT0 * TEMP2
         MAPSOUTH(I,11)=  TXINT0 * TYINT4 * TEMP2
         MAPSOUTH(I,12)= -TXINT1 * TYINT3 * TEMP2
         MAPSOUTH(I,13)=  TXINT2 * TYINT2 * TEMP2
         MAPSOUTH(I,14)= -TXINT3 * TYINT1 * TEMP2
         MAPSOUTH(I,15)=  TXINT4 * TYINT0 * TEMP2
         MAPSOUTH(I,16)=  TXINT0 * TYINT5 * TEMP2
         MAPSOUTH(I,17)= -TXINT1 * TYINT4 * TEMP2
         MAPSOUTH(I,18)=  TXINT2 * TYINT3 * TEMP2
         MAPSOUTH(I,19)= -TXINT3 * TYINT2 * TEMP2
         MAPSOUTH(I,20)=  TXINT4 * TYINT1 * TEMP2
         MAPSOUTH(I,21)= -TXINT5 * TYINT0 * TEMP2
         MAPSOUTH(I,22)=  TXINT0 * TYINT6 * TEMP2
         MAPSOUTH(I,23)= -TXINT1 * TYINT5 * TEMP2
         MAPSOUTH(I,24)=  TXINT2 * TYINT4 * TEMP2
         MAPSOUTH(I,25)= -TXINT3 * TYINT3 * TEMP2
         MAPSOUTH(I,26)=  TXINT4 * TYINT2 * TEMP2
         MAPSOUTH(I,27)= -TXINT5 * TYINT1 * TEMP2
         MAPSOUTH(I,28)=  TXINT6 * TYINT0 * TEMP2
         MAPSOUTH(I,29)=  TXINT0 * TYINT7 * TEMP2
         MAPSOUTH(I,30)= -TXINT1 * TYINT6 * TEMP2
         MAPSOUTH(I,31)=  TXINT2 * TYINT5 * TEMP2
         MAPSOUTH(I,32)= -TXINT3 * TYINT4 * TEMP2
         MAPSOUTH(I,33)=  TXINT4 * TYINT3 * TEMP2
         MAPSOUTH(I,34)= -TXINT5 * TYINT2 * TEMP2
         MAPSOUTH(I,35)=  TXINT6 * TYINT1 * TEMP2
         MAPSOUTH(I,36)= -TXINT7 * TYINT0 * TEMP2

      END DO
      RETURN
      END


C**************************************************************************
C     The following subroutine computes an integral needed when
C     forming the big to small maps for the volume part
C
C     INPUT:
C
C     XNODE is the value of xnodes(i) here
C
C     XPOW is as power that X/2 is raised to
C
C     OUTPUT:
C
C     XINT is set equal to int -1 < x < 1 (x/2)^xpow * exp(xnode*x) dx
C     Note that XINT must be real here
C
C**************************************************************************
      SUBROUTINE GETINTSX(XNODE, XINT, XPOW)
      IMPLICIT NONE
C-----Global variables
      INTEGER  XPOW
      REAL *8  XNODE
      REAL *8  XINT
C-----Local variables
      INTEGER  KPTS, KIND, M
      REAL *8  ALPHA, BETA, ENDPTS(2), B(60), CHEB
      REAL *8  T(60), WNODE(60)

      ALPHA = 0.0D0
      BETA = 0.0D0
      KPTS = 0
      ENDPTS(1) = 0.0D0
      ENDPTS(2) = 0.0D0
      KIND = 1

      CALL GAUSSQ(KIND, 60, ALPHA, BETA, KPTS, ENDPTS, B, T, WNODE)

        XINT = 0.0D0
        DO M = 1,60
         XINT = XINT + WNODE(M) * CHEB(T(M),XPOW) * DEXP(XNODE*T(M))
       END DO
      RETURN
      END



C**************************************************************************
C     The following subroutine computes an integral needed when
C     forming the big to small maps for the volume part
C
C     INPUT:
C
C     XNODE is the value of xnodes(i) here
C
C     YPOW is as power that X/2 is raised to
C
C     OUTPUT:
C
C     XINT is set equal to int -1 < x < 1 (x/2)^ypow * exp(i*xnode*x) dx
C     Note that XINT can be complex here
C
C**************************************************************************
      SUBROUTINE GETINTSY(XNODE, XINT, YPOW)
      IMPLICIT NONE
C-----Global variables
      INTEGER  YPOW
      REAL *8  XNODE
      COMPLEX *16  XINT
C-----Local variables
      INTEGER  KPTS, KIND, L
      REAL *8  ALPHA, BETA, ENDPTS(2), B(30), CHEB
      REAL *8  T(30), WNODE(30)
      COMPLEX *16  IMAG
      DATA IMAG/(0.0,1.0)/

      ALPHA = 0.0D0
      BETA = 0.0D0
      KPTS = 0
      ENDPTS(1) = 0.0D0
      ENDPTS(2) = 0.0D0
      KIND = 1

      CALL GAUSSQ(KIND, 30, ALPHA, BETA, KPTS, ENDPTS, B, T, WNODE)

        XINT = 0.0D0
        DO L = 1,30
           XINT = XINT + WNODE(L) * CHEB(T(L),YPOW) 
     1                            * CDEXP(XNODE*IMAG*T(L))
        END DO
      RETURN
      END


C**************************************************************************
C     This subroutine generates four exponential expansions given
C     one multipole expansion as input.
C
C     INPUT:
C
C     NNODES is the order of the plane wave expansions
C
C     MPOLE is the array or multipole coefficients
C
C     NTERMS is the order of the multipole expansions
C
C     COMP is a precomputed term involving the ratio of
C          the weights and factorial terms
C
C     WNODES is the array of exponential weights
C
C     XNODES is the array of exponential nodes
C
C     SUM is a precomputed real term
C
C     TEMP(I) is set equal to WNODES(I) / XNODES(I)
C
C     OUTPUT:
C
C     BETAE is set to the east exponential coefficients
C
C     BETAW is set to the west exponential coefficients
C
C     BETAN is set to the north exponential coefficients
C
C     BETAS is set to the south exponential coefficients
C
C**************************************************************************
      SUBROUTINE EXPCOEFF(BETAE, BETAW, BETAN, BETAS,
     1                    NNODES, MPOLE, NTERMS, COMP,
     2                    WNODES, SUM, TEMP)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NNODES, NTERMS
      REAL *8  COMP(NNODES,0:NTERMS)
      REAL *8  WNODES(NNODES)
      REAL *8  TEMP(NNODES)
      REAL *8  SUM
      COMPLEX *16  BETAE(0:NNODES), BETAW(0:NNODES)
      COMPLEX *16  BETAN(0:NNODES), BETAS(0:NNODES)
      COMPLEX *16  MPOLE(0:NTERMS)
C-----Local variables
      INTEGER  I, J
      COMPLEX *16  IMAG

      COMPLEX *16 SUM1, SUM2, SUM5, SUM6
      COMPLEX *16 SUM7, SUM8
      COMPLEX *16 TEMP1

      DATA IMAG/(0.0D0, 1.0D0)/

      TEMP1 = MPOLE(0)*SUM

      BETAE(0) = TEMP1
      BETAW(0) = TEMP1
      BETAS(0) = TEMP1
      BETAN(0) = TEMP1

      DO I = 1, NNODES

         SUM1 = 0.0D0
         SUM2 = 0.0D0
         SUM5 = 0.0D0
         SUM6 = 0.0D0

         DO J = 1, NTERMS-3, 4
            SUM1 = SUM1 - COMP(I,J)   * MPOLE(J+1)
            SUM2 = SUM2 - COMP(I,J-1) * MPOLE(J)
            SUM5 = SUM5 - COMP(I,J+1) * MPOLE(J+2)
            SUM6 = SUM6 - COMP(I,J+2) * MPOLE(J+3)
         END DO

         SUM8 = IMAG*(SUM2 - SUM5)
         SUM7 = SUM6 - SUM1
         SUM1 = SUM1 + SUM6
         SUM2 = SUM2 + SUM5
         BETAE(I) = WNODES(I) * (SUM1 + SUM2)
         BETAW(I) = WNODES(I) * (SUM1 - SUM2)
         BETAS(I) = WNODES(I) * (SUM7 + SUM8)
         BETAN(I) = WNODES(I) * (SUM7 - SUM8)

      END DO

C     NOW LET'S READJUST THE EAST COEFFICIENTS TO ACCOUNT FOR THE
C     DIFFERENCE IN THE LOG TERM
      DO I = 1, NNODES
       TEMP1 =  MPOLE(0)*TEMP(I)
       BETAE(I) = BETAE(I) + TEMP1
       BETAW(I) = BETAW(I) + TEMP1
       BETAS(I) = BETAS(I) + TEMP1
       BETAN(I) = BETAN(I) + TEMP1
      END DO
      RETURN
      END


      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
C
C           THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS
C        W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
C
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
C
C                              N
C        BY                   SUM W  F(T )
C                             J=1  J    J
C
C        (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
C
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
C        ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.
C
C     INPUT PARAMETERS (ALL REAL NUMBERS ARE IN REAL*8)
C
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
C                 QUADRATURE RULE:
C
C        KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
C        KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
C        KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)
C        KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
C                   (-INFINITY, +INFINITY)
C        KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.
C                   NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
C        KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1
C
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.D0).
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
C                 (OTHERWISE USE 0.D0)
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
C                 ENDPOINTS (1 OR 2).
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
C        B        REAL SCRATCH ARRAY OF LENGTH N
C
C     OUTPUT PARAMETERS (BOTH REAL*8 ARRAYS OF LENGTH N)
C
C        T        WILL CONTAIN THE DESIRED NODES.
C        W        WILL CONTAIN THE DESIRED WEIGHTS W(J).
C
C     SUBROUTINES REQUIRED
C
C        SOLVE, CLASS, AND IMTQL2 ARE PROVIDED.  UNDERFLOW MAY SOMETIMES
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
C        TURNED OFF.  TO DO THIS, THE FIRST CALL OF THE MAIN PROGRAM
C        SHOULD BE
C                  CALL TRAPS (0, 0, 10000, 0, 0)    IN WATFIV
C        OR
C                  CALL INIT                         IN FORTRAN G OR H.
C
C     ACCURACY
C
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
C        COMPLETELY HARMLESS.
C
C     METHOD
C
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
C
C     REFERENCES
C
C        1.  GOLUB, G. H., AND WELSCH, J. H., "CALCULATION OF GAUSSIAN
C            QUADRATURE RULES," MATHEMATICS OF COMPUTATION 23 (APRIL,
C            1969), PP. 221-230.
C        2.  GOLUB, G. H., "SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,"
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.
C
C     ..................................................................
C
      REAL *8 B(N), T(N), W(N), ENDPTS(2), MUZERO, T1,
     X GAM, SOLVE, DSQRT, ALPHA, BETA
C
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
C
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
C           B THE OFF-DIAGONAL ELEMENTS.
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
C           SUBMATRIX.
C
      IF (KPTS.EQ.0)  GO TO 100
      IF (KPTS.EQ.2)  GO TO  50
C
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED
C
      T(N) = SOLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
      GO TO 100
C
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
C
   50 GAM = SOLVE(ENDPTS(1), N, T, B)
      T1 = ((ENDPTS(1) - ENDPTS(2))/(SOLVE(ENDPTS(2), N, T, B) - GAM))
      B(N-1) = DSQRT(T1)
      T(N) = ENDPTS(1) + GAM*T1
C
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
C
  100 W(1) = 1.0D0
      DO 105 I = 2, N
  105    W(I) = 0.0D0
C
      CALL IMTQL2 (N, T, B, W, IERR)
      DO 110 I = 1, N
  110    W(I) = MUZERO * W(I) * W(I)
C
      RETURN
      END
C
C
C
      REAL*8 FUNCTION SOLVE(SHIFT, N, A, B)
C
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
C
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,
C
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
C       THE N-TH POSITION.
C
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
C
C
      REAL*8 SHIFT, A(N), B(N), ALPHA
C
      ALPHA = A(1) - SHIFT
      NM1 = N - 1
      DO 10 I = 2, NM1
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
      SOLVE = 1.0D0/ALPHA
      RETURN
      END
C
C
C
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
C
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
C        RECURRENCE RELATION
C
C             B P (X) = (X - A ) P   (X) - B   P   (X)
C              J J            J   J-1       J-1 J-2
C
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
C        AND THE ZERO-TH MOMENT
C
C             MUZERO = INTEGRAL W(X) DX
C
C        OF THE GIVEN POLYNOMIAL'S WEIGHT FUNCTION W(X).  SINCE THE
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
C        GUARANTEED TO BE SYMMETRIC.
C
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
C        REQUIRE THE GAMMA FUNCTION.
C
C     ..................................................................
C
      REAL*8 A(N), B(N), MUZERO, ALPHA, BETA
      REAL*8 ABI, A2B2, PI, DSQRT, AB
      DATA PI / 3.141592653589793D0/
C
      NM1 = N - 1
      GO TO (10, 20, 30, 40, 50, 60), KIND
C
C              KIND = 1:  LEGENDRE POLYNOMIALS P(X)
C              ON (-1, +1), W(X) = 1.
C
   10 MUZERO = 2.0D0
      DO 11 I = 1, NM1
         A(I) = 0.0D0
         ABI = I
   11    B(I) = ABI/DSQRT(4*ABI*ABI - 1.0D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
C
   20 MUZERO = PI
      DO 21 I = 1, NM1
         A(I) = 0.0D0
   21    B(I) = 0.5D0
      B(1) = DSQRT(0.5D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
C              ON (-1, +1), W(X) = SQRT(1 - X*X)
C
   30 MUZERO = PI/2.0D0
      DO 31 I = 1, NM1
         A(I) = 0.0D0
   31    B(I) = 0.5D0
      A(N) = 0.0D0
      RETURN
C
C              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
C              +INFINITY), W(X) = EXP(-X**2)
C
   40 MUZERO = DSQRT(PI)
      DO 41 I = 1, NM1
         A(I) = 0.0D0
   41    B(I) = DSQRT(I/2.0D0)
      A(N) = 0.0D0
      RETURN
C
C              KIND = 5:  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
C              BETA GREATER THAN -1
C
   50 AB = ALPHA + BETA
      ABI = 2.0D0 + AB
      MUZERO = 2.0D0 ** (AB + 1.0D0)
      A(1) = (BETA - ALPHA)/ABI
      B(1) = DSQRT(4.0D0*(1.0D0 + ALPHA)*(1.0D0 + BETA)/((ABI + 1.0D0)*
     1  ABI*ABI))
      A2B2 = BETA*BETA - ALPHA*ALPHA
      DO 51 I = 2, NM1
         ABI = 2.0D0*I + AB
         A(I) = A2B2/((ABI - 2.0D0)*ABI)
   51    B(I) = DSQRT (4.0D0*I*(I + ALPHA)*(I + BETA)*(I + AB)/
     1   ((ABI*ABI - 1)*ABI*ABI))
      ABI = 2.0D0*N + AB
      A(N) = A2B2/((ABI - 2.0D0)*ABI)
      RETURN
C
C              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
C              THAN -1.
C
   60 MUZERO = 1.0D0
      DO 61 I = 1, NM1
         A(I) = 2.0D0*I - 1.0D0 + ALPHA
   61    B(I) = DSQRT(I*(I + ALPHA))
      A(N) = 2.0D0*N - 1 + ALPHA
      RETURN
      END
 
 
      SUBROUTINE IMTQL2(N, D, E, Z, IERR)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C     THIS IS A MODIFIED VERSION OF THE 'EISPACK' ROUTINE IMTQL2.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
C     METHOD.
C
C     ON INPUT:
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
C
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
C
C      ON OUTPUT:
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
C
C        E HAS BEEN DESTROYED;
C
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ------------------------------------------------------------------
C
      INTEGER I, J, K, L, M, N, II, MML, IERR
      REAL*8 D(N), E(N), Z(N), B, C, F, G, P, R, S, MACHEP
      REAL*8 DSQRT, DABS, DSIGN
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
CCCC  DATA MACHEP/Z3410000000000000/
      DATA MACHEP/1.0D-14/
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (DABS(E(M)) .LE. MACHEP * (DABS(D(M)) + DABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = DSQRT(G*G+1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R, G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (DABS(F) .LT. DABS(G)) GO TO 150
            C = G / F
            R = DSQRT(C*C+1.0D0)
            E(I+1) = F * R
            S = 1.0D0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = DSQRT(S*S+1.0D0)
            E(I+1) = G * R
            C = 1.0D0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     :::::::::: FORM FIRST COMPONENT OF VECTOR ::::::::::
            F = Z(I+1)
            Z(I+1) = S * Z(I) + C * F
  200       Z(I) = C * Z(I) - S * F
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
  300 CONTINUE
C
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
      END


C**************************************************************************
C--------------------------------------------------------------------------
C
C     Differentiate Module
C
C--------------------------------------------------------------------------
C**************************************************************************

C**********************************************************************
C     The following subroutine is set up to calculate both the x and y
C     derivatives of a given function. 
C
C     INPUT:
C
C     POT(16,NBOXES) is the quantity being differentiated 
C         defined on the leaf nodes
C
C     LEVELBOX is an array determining the level of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C     ICHILDBOX denotes the four children of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     COEFFS is a dummy array needed to store the polynomial coefficients
C            that approximate POT
C
C     OUTPUT:
C
C     POTX(16,NBOXES) is the x derivative defined on the leaf nodes
C
C     POTY(16,NBOXES) is the y derivative defined on the leaf nodes
C
C**********************************************************************
      SUBROUTINE GETDERIVATIVES8(POT,
     1      NLEV, ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV,
     2      POTX, POTY, LAP, COEFFS)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  COEFFS(0:7,0:7,1)
      REAL *8  POT(64,1)
      REAL *8  POTX(64,1), POTY(64,1)
      REAL *8  LAP(64,1) 

C       Call a routine to get the matrix A, which will be 
C       needed to calculate the interpolating polynomial
C       for POT.

C       Now call a routine that will approximate POT with
C       a polynomial. 
        CALL MKCOEFFS8(COEFFS,POT,NLEV,
     1    ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV)

C       Now call a routine that will calculate the x and y
C       derivatives of POT, given the approximating polynomial
C       COEFFS as input.
        CALL DIFFERENTIATE8(COEFFS,POTX,POTY,LAP,
     1           NLEV, ICHILDBOX,
     2           NBLEVEL, IBOXLEV, ISTARTLEV)

       RETURN
       END


      SUBROUTINE GETPOTPRACTICAL(POT, POTX, POTY,
     1                POTPRACTICAL, POTPRACTICALX, 
     2                POTPRACTICALY,
     3                NLEV, ICHILDBOX, NBLEVEL,
     4                IBOXLEV, ISTARTLEV, COEFFS)
      IMPLICIT NONE
      INTEGER  NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  POT(64,1), POTX(64,1), POTY(64,1)
      REAL *8  POTPRACTICAL(81,1)
      REAL *8  POTPRACTICALX(81,1), POTPRACTICALY(81,1)
      REAL *8  COEFFS(0:7,0:7,1)

C       Now call a routine that will approximate POT with
C       a polynomial.
        CALL MKCOEFFS8(COEFFS,POT,NLEV,
     1    ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV)


        CALL PRACTICAL(COEFFS, POTPRACTICAL,
     1              NLEV, ICHILDBOX,
     2              NBLEVEL, IBOXLEV, ISTARTLEV)

C       Now do the same for the x derivative
        CALL MKCOEFFS8(COEFFS,POTX,NLEV,
     1    ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV)


        CALL PRACTICAL(COEFFS, POTPRACTICALX,
     1              NLEV, ICHILDBOX,
     2              NBLEVEL, IBOXLEV, ISTARTLEV)

C       Finally, do the same for the y derivative
        CALL MKCOEFFS8(COEFFS,POTY,NLEV,
     1    ICHILDBOX, NBLEVEL, IBOXLEV, ISTARTLEV)


        CALL PRACTICAL(COEFFS, POTPRACTICALY,
     1              NLEV, ICHILDBOX,
     2              NBLEVEL, IBOXLEV, ISTARTLEV)


      RETURN
      END


      SUBROUTINE PRACTICAL(COEFFS,POTPRACTICAL,
     1               NLEV, ICHILDBOX,
     2               NBLEVEL, IBOXLEV, ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  COEFFS(0:7,0:7,1)
      REAL *8  POTPRACTICAL(81,1)
      REAL *8  CHEB
C-----Local variables
      INTEGER  I
      INTEGER  J, II, JJ, IBOX, III, JJJ
      REAL *8  PI16
      REAL *8  XF(81), YF(81), X(9), TX(0:7), TY(0:7)

      PI16 = DACOS(-1.0D0) / 16.0D0
      X(1) = -1.0D0
      X(2) = DCOS(14.0D0*PI16)
      X(3) = DCOS(12.0D0*PI16)
      X(4) = DCOS(10.0D0*PI16)
      X(5) = 0.0D0
      X(6) = DCOS( 6.0D0*PI16)
      X(7) = DCOS( 4.0D0*PI16)
      X(8) = DCOS( 2.0D0*PI16)
      X(9) = 1.0D0


      DO JJ = 0,NLEV
       DO I = 1, 9
       DO J = 1, 9
         XF(I+9*(J-1)) = X(I)
         YF(J+9*(I-1)) = X(I)
       END DO
       END DO


        DO 100 II = ISTARTLEV(JJ), ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
         IBOX = IBOXLEV(II)
C        Only find the derivative for childless boxes:
         IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100

         DO I = 1, 81

         DO III = 0, 7
         DO JJJ = 0, 7
           TX(III) = CHEB(XF(I),III)
           TY(JJJ) = CHEB(YF(I),JJJ)
         END DO
         END DO

         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,0,IBOX)*TX(0)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,0,IBOX)*TX(1)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,1,IBOX)*TX(0)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,0,IBOX)*TX(2)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,1,IBOX)*TX(1)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,2,IBOX)*TX(0)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(3,0,IBOX)*TX(3)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,1,IBOX)*TX(2)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,2,IBOX)*TX(1)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,3,IBOX)*TX(0)* TY(3)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(4,0,IBOX)*TX(4)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(3,1,IBOX)*TX(3)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,2,IBOX)*TX(2)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,3,IBOX)*TX(1)* TY(3)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,4,IBOX)*TX(0)* TY(4)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(5,0,IBOX)*TX(5)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(4,1,IBOX)*TX(4)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(3,1,IBOX)*TX(3)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,3,IBOX)*TX(2)* TY(3)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,4,IBOX)*TX(1)* TY(4)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,5,IBOX)*TX(0)* TY(5)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(6,0,IBOX)*TX(6)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(5,1,IBOX)*TX(5)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(4,2,IBOX)*TX(4)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(3,3,IBOX)*TX(3)* TY(3)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,4,IBOX)*TX(2)* TY(4)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,5,IBOX)*TX(1)* TY(5)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,6,IBOX)*TX(0)* TY(6)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(7,0,IBOX)*TX(7)* TY(0)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(6,1,IBOX)*TX(6)* TY(1)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(5,2,IBOX)*TX(5)* TY(2)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(4,3,IBOX)*TX(4)* TY(3)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(3,4,IBOX)*TX(3)* TY(4)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(2,5,IBOX)*TX(2)* TY(5)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(1,6,IBOX)*TX(1)* TY(6)
         POTPRACTICAL(I,IBOX) = POTPRACTICAL(I,IBOX) 
     1                        + COEFFS(0,7,IBOX)*TX(0)* TY(7)

        END DO
100     CONTINUE
      END DO
      RETURN
      END


C**********************************************************************
C     The following subroutine computes the x and y derivatives of the 
C     solution, given the approximating polynomial and tree as input.
C     The algorithm works just by differentiating the approximating
C     polynomial.
C
C
C     INPUT:
C
C     COEFFS is a the array containing polynomial coefficients
C            that approximate POT
C
C     NLEV is the finest level
C
C     NBOXES is the total number of boxes
C
C     ICHILDBOX denotes the four children of each box
C
C     NBLEVEL is the total number of boxes per level
C
C     IBOXLEV is the array in which the boxes are arranged
C
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     OUTPUT:
C
C     POTX(64,NBOXES) is the x derivative defined on the leaf nodes
C
C     POTY(64,NBOXES) is the y derivative defined on the leaf nodes
C
C**********************************************************************
      SUBROUTINE DIFFERENTIATE8(COEFFS,POTX,POTY,LAP, 
     1               NLEV, ICHILDBOX,
     2               NBLEVEL, IBOXLEV, ISTARTLEV)
      IMPLICIT NONE
C-----Global variables
      INTEGER  NLEV
      INTEGER  ICHILDBOX(4,1)
      INTEGER  NBLEVEL(0:1), IBOXLEV(1), ISTARTLEV(0:1)
      REAL *8  COEFFS(0:7,0:7,1)
      REAL *8  POTX(64,1), POTY(64,1)
      REAL *8  LAP(64,1)
      REAL *8  CHEB
C-----Local variables
      INTEGER  I
      INTEGER  J, II, JJ, IBOX, III, JJJ
      REAL *8  PI16, XLENGTH
      REAL *8  XF(64), YF(64), X(8), TX(0:7), TY(0:7)
      REAL *8  TXDIFF(7), TYDIFF(7)
      REAL *8  TXLAP(2:7), TYLAP(2:7)

      PI16 = DACOS(-1.0D0) / 16.0D0
      X(1) = DCOS(15.0D0*PI16)
      X(2) = DCOS(13.0D0*PI16)
      X(3) = DCOS(11.0D0*PI16)
      X(4) = DCOS( 9.0D0*PI16)
      X(5) = DCOS( 7.0D0*PI16)
      X(6) = DCOS( 5.0D0*PI16)
      X(7) = DCOS( 3.0D0*PI16)
      X(8) = DCOS( 1.0D0*PI16)


      XLENGTH = 0.5D0
      DO JJ = 0,NLEV
       DO I = 1, 8
       DO J = 1, 8
         XF(I+8*(J-1)) = X(I)
         YF(J+8*(I-1)) = X(I)
       END DO
       END DO


        DO 100 II = ISTARTLEV(JJ), ISTARTLEV(JJ) + NBLEVEL(JJ) - 1
         IBOX = IBOXLEV(II)
C        Only find the derivative for childless boxes:
         IF(ICHILDBOX(1,IBOX) .GT. 0)GOTO 100

         DO I = 1, 64

         DO III = 0, 7
         DO JJJ = 0, 7
           TX(III) = CHEB(XF(I),III)
           TY(JJJ) = CHEB(YF(I),JJJ)
         END DO
         END DO

         TXDIFF(1) =        TX(0)
         TXDIFF(2) =  4.0D0*TX(1)
         TXDIFF(3) =  6.0D0*TX(2) +  3.0D0*TX(0)
         TXDIFF(4) =  8.0D0*TX(3) +  8.0D0*TX(1)
         TXDIFF(5) = 10.0D0*TX(4) + 10.0D0*TX(2) +  5.0D0*TX(0)
         TXDIFF(6) = 12.0D0*TX(5) + 12.0D0*TX(3) + 12.0D0*TX(1)
         TXDIFF(7) = 14.0D0*TX(6) + 14.0D0*TX(4) 
     1             + 14.0D0*TX(2) +  7.0D0*TX(0)

         TYDIFF(1) =        TY(0)
         TYDIFF(2) =  4.0D0*TY(1)
         TYDIFF(3) =  6.0D0*TY(2) +  3.0D0*TY(0)
         TYDIFF(4) =  8.0D0*TY(3) +  8.0D0*TY(1)
         TYDIFF(5) = 10.0D0*TY(4) + 10.0D0*TY(2) +  5.0D0*TY(0)
         TYDIFF(6) = 12.0D0*TY(5) + 12.0D0*TY(3) + 12.0D0*TY(1)
         TYDIFF(7) = 14.0D0*TY(6) + 14.0D0*TY(4) 
     1             + 14.0D0*TY(2) +  7.0D0*TY(0)

         DO III = 1, 7
           TXDIFF(III) = TXDIFF(III) / XLENGTH
           TYDIFF(III) = TYDIFF(III) / XLENGTH
         END DO

         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,0,IBOX)*TXDIFF(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,0,IBOX)*TXDIFF(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(3,0,IBOX)*TXDIFF(3)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(4,0,IBOX)*TXDIFF(4)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(5,0,IBOX)*TXDIFF(5)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(6,0,IBOX)*TXDIFF(6)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(7,0,IBOX)*TXDIFF(7)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,1,IBOX)*TXDIFF(1)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,1,IBOX)*TXDIFF(2)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(3,1,IBOX)*TXDIFF(3)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(4,1,IBOX)*TXDIFF(4)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(5,1,IBOX)*TXDIFF(5)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(6,1,IBOX)*TXDIFF(6)*TY(1)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,2,IBOX)*TXDIFF(1)*TY(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,2,IBOX)*TXDIFF(2)*TY(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(3,2,IBOX)*TXDIFF(3)*TY(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(4,2,IBOX)*TXDIFF(4)*TY(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(5,2,IBOX)*TXDIFF(5)*TY(2)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,3,IBOX)*TXDIFF(1)*TY(3)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,3,IBOX)*TXDIFF(2)*TY(3)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(3,3,IBOX)*TXDIFF(3)*TY(3)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(4,3,IBOX)*TXDIFF(4)*TY(3)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,4,IBOX)*TXDIFF(1)*TY(4)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,4,IBOX)*TXDIFF(2)*TY(4)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(3,4,IBOX)*TXDIFF(3)*TY(4)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,5,IBOX)*TXDIFF(1)*TY(5)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(2,5,IBOX)*TXDIFF(2)*TY(5)
         POTX(I,IBOX) = POTX(I,IBOX) + COEFFS(1,6,IBOX)*TXDIFF(1)*TY(6)

         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,1,IBOX)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,2,IBOX)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,3,IBOX)*TYDIFF(3)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,4,IBOX)*TYDIFF(4)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,5,IBOX)*TYDIFF(5)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,6,IBOX)*TYDIFF(6)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(0,7,IBOX)*TYDIFF(7)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,1,IBOX)*TX(1)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,2,IBOX)*TX(1)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,3,IBOX)*TX(1)*TYDIFF(3)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,4,IBOX)*TX(1)*TYDIFF(4)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,5,IBOX)*TX(1)*TYDIFF(5)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(1,6,IBOX)*TX(1)*TYDIFF(6)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(2,1,IBOX)*TX(2)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(2,2,IBOX)*TX(2)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(2,3,IBOX)*TX(2)*TYDIFF(3)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(2,4,IBOX)*TX(2)*TYDIFF(4)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(2,5,IBOX)*TX(2)*TYDIFF(5)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(3,1,IBOX)*TX(3)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(3,2,IBOX)*TX(3)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(3,3,IBOX)*TX(3)*TYDIFF(3)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(3,4,IBOX)*TX(3)*TYDIFF(4)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(4,1,IBOX)*TX(4)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(4,2,IBOX)*TX(4)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(4,3,IBOX)*TX(4)*TYDIFF(3)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(5,1,IBOX)*TX(5)*TYDIFF(1)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(5,2,IBOX)*TX(5)*TYDIFF(2)
         POTY(I,IBOX) = POTY(I,IBOX) + COEFFS(6,1,IBOX)*TX(6)*TYDIFF(1)


         TXLAP(2) =   4.0D0*TX(0)
         TXLAP(3) =  24.0D0*TX(1)
         TXLAP(4) =  48.0D0*TX(2) +  32.0D0*TX(0)
         TXLAP(5) =  80.0D0*TX(3) + 120.0D0*TX(1)
         TXLAP(6) = 120.0D0*TX(4) + 192.0D0*TX(2) + 108.0D0*TX(0)
         TXLAP(7) = 168.0D0*TX(5) + 280.0D0*TX(3) + 336.0D0*TX(1)


         TYLAP(2) =   4.0D0*TY(0)
         TYLAP(3) =  24.0D0*TY(1)
         TYLAP(4) =  48.0D0*TY(2) +  32.0D0*TY(0)
         TYLAP(5) =  80.0D0*TY(3) + 120.0D0*TY(1)
         TYLAP(6) = 120.0D0*TY(4) + 192.0D0*TY(2) + 108.0D0*TY(0)
         TYLAP(7) = 168.0D0*TY(5) + 280.0D0*TY(3) + 336.0D0*TY(1)

         DO III = 2, 7
           TXLAP(III) = TXLAP(III) / XLENGTH**2
           TYLAP(III) = TYLAP(III) / XLENGTH**2
         END DO


         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,0,IBOX)*TXLAP(2)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,2,IBOX)*TYLAP(2)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(3,0,IBOX)*TXLAP(3) 
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,1,IBOX)*TXLAP(2)*TY(1) 
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(1,2,IBOX)*TX(1)*TYLAP(2)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,3,IBOX)*TYLAP(3)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(4,0,IBOX)*TXLAP(4)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(3,1,IBOX)*TXLAP(3)*TY(1)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,2,IBOX)*
     1                             (TX(2)*TYLAP(2) + TXLAP(2)*TY(2))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(1,3,IBOX)*TX(1)*TYLAP(3)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,4,IBOX)*TYLAP(4)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(5,0,IBOX)*TXLAP(5)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(4,1,IBOX)*TXLAP(4)*TY(1)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(3,2,IBOX)*
     1                             (TX(3)*TYLAP(2) + TXLAP(3)*TY(2))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,3,IBOX)*
     1                             (TX(2)*TYLAP(3) + TXLAP(2)*TY(3))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(1,4,IBOX)*TX(1)*TYLAP(4)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,5,IBOX)*TYLAP(5)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(6,0,IBOX)*TXLAP(6)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(5,1,IBOX)*TXLAP(5)*TY(1)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(4,2,IBOX)*
     1                             (TXLAP(4)*TY(2) + TX(4)*TYLAP(2))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(3,3,IBOX)*
     1                             (TXLAP(3)*TY(3) + TX(3)*TYLAP(3))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,4,IBOX)*
     1                             (TXLAP(2)*TY(4) + TX(2)*TYLAP(4))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(1,5,IBOX)*TX(1)*TYLAP(5)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,6,IBOX)*TYLAP(6)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(7,0,IBOX)*TXLAP(7)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(6,1,IBOX)*TXLAP(6)*TY(1)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(5,2,IBOX)*
     1                             (TXLAP(5)*TY(2) + TX(5)*TYLAP(2))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(4,3,IBOX)*
     1                             (TXLAP(4)*TY(3) + TX(4)*TYLAP(3))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(3,4,IBOX)*
     1                             (TXLAP(3)*TY(4) + TX(3)*TYLAP(4))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(2,5,IBOX)*
     1                             (TXLAP(2)*TY(5) + TX(2)*TYLAP(5))
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(1,6,IBOX)*TX(1)*TYLAP(6)
         LAP(I,IBOX) = LAP(I,IBOX) + COEFFS(0,7,IBOX)*TYLAP(7)


        END DO
100     CONTINUE
       XLENGTH = XLENGTH / 2.0D0
      END DO
      RETURN
      END

