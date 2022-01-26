
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
c     exact solution to the Poisson equation for the above RHS function
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


