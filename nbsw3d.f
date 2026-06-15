      SUBROUTINE NBSW3D
     1(   ICOOR      ,LACCOC     ,LACCOP     ,LACCOT     ,
     2    MCNTS      ,NDIMN       ,NDIM2     ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	LOGICAL
     1    UPDMTA
      CHARACTER
     1    NAME*6
      DIMENSION
     1    ICOOR(NDIM2,N)    ,LACCOC(N)          ,
     2    LACCOP(N)         ,LACCOT(MCNTS) 
      DIMENSION
     1    ICELSZ(3),MARK(N),NEXTX(N),NEXTY(N),NEXTZ(N),IPAIR(2,MCNTS),
     2    IXSLIC(N),IYSLIC(N),LSTAPP(N),LSTARC(N),LSTARP(N),LSTARY(N)  
      integer, allocatable :: IHEADX(:),IHEADY(:),IHEADZ(:),MARKX(:),
     1                        MARKY(:),MARKZ(:),ITAILZ(:)                    
      DATA
     1    R2P5       ,R0         /
     2    2.5D0      ,0.0D0      /
      DATA NAME/'NBSW3D'/
C***********************************************************************
C*ACRONYM
C NBS contact dtection 3D
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICOOR  - Bounding box coordinates
C=LACCOC - List of number of contacted bodies
C=LACCOP - List of pointers to contacted bodies
C=LACCOT - List of targetor numbers
C LNODS  - List of element nodes
C*Variables 
C ICTYPE - Dem contact type
C MCNTS  - Max number of contacting couples
C MFNOD  - Maximum number of segment nodes
C N      - Number of bounding boxes
C NCPOIN - Number of contact points
C NCSEG  - Number of contact segments
C NDIMN  - Problem dimension
C NDIM2  - Search space dimension
C NDISKS - Number of grains
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      UPDMTA=.FALSE.
C Comput average cell size
      XMEAN=R0
      YMEAN=R0
      ZMEAN=R0
      DO 10 I=1,N
        XMEAN=XMEAN+ICOOR(4,I)-ICOOR(1,I)
        YMEAN=YMEAN+ICOOR(5,I)-ICOOR(2,I)
        ZMEAN=ZMEAN+ICOOR(6,I)-ICOOR(3,I)
   10 CONTINUE
      ICELSZ(1)=XMEAN/N
      ICELSZ(2)=YMEAN/N
      ICELSZ(3)=ZMEAN/N
C Evaluate minumum and maximum coordinates of region and divide into cells
      CALL NDVCL3
     1( ICELSZ     ,ICOOR      ,
     2  IXMIN      ,IYMIN      ,IZMIN      ,N          ,NDIM2       ,
     3  NX         ,NY         ,NZ         )
C Allocate spaces
C      allocate(IHEADX(NX))
      allocate(IHEADX(NX),IHEADY(NY),IHEADZ(NZ),MARKX(NX),MARKY(NY),
     1         MARKZ(NZ),ITAILZ(NZ))
C Buld iz list
      CALL NBIZL3
	1( IHEADZ     ,ICELSZ     ,ICOOR      ,ITAILZ      ,IXSLIC     ,
     2  IYSLIC     ,MARK       ,MARKZ      ,NEXTZ       ,
     3  IXMIN      ,IYMIN      ,IZMIN      ,N           ,NDIM2      ,
     4  NZ         )
C     Contact detection
      MP=0
      CALL NZCNT3
     1( IHEADX     ,IHEADY     ,IHEADZ     ,ICELSZ     ,ICOOR       ,
     2  IPAIR      ,ITAILZ     ,IXSLIC     ,IYSLIC     ,LSTARC      ,
     3  LSTARP     ,LSTARY     ,MARK       ,MARKX      ,MARKY       ,
     4  MARKZ      ,NEXTX      ,NEXTY      ,NEXTZ      ,
     5  IXMIN      ,IYMIN      ,IZMIN      ,MCNTS      ,MP          ,
     6  N          ,NDIM2      ,NX         ,NY         ,NZ          ,
     7  UPDMTA     ) 
C Convert to ELFEN contact list format
      CALL NOUTCV
     1( ICOOR      ,IPAIR      ,LACCOC     ,LACCOP     ,LACCOT     ,
     2  MP         ,NDIMN      ,NDIM2      ,N          )
	deallocate(IHEADX,IHEADY,IHEADZ,MARKX,MARKY,MARKZ,ITAILZ)
      RETURN
      END


      SUBROUTINE NDVCL3
     1( ICELSZ     ,ICOOR      ,
     2  IXMIN      ,IYMIN      ,IZMIN      ,N          ,NDIM2      ,
     3  NX         ,NY         ,NZ         )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  ICELSZ(3)          ,ICOOR(NDIM2,N) 
      DATA NAME/'NDVCL3'/
C***********************************************************************
C Evaluate minumum and maximum coordinates of region and divide into cells
C*ACRONYM
C Divide search space into cells for 3D NBS contact dtection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICELSZ - Cell sizes
C ICOOR  - Bounding box coordinates
C*Variables
C IXMIN  - Minimum x-coordinate of search space 
C IYMIN  - Minimum y-coordinate of search space 
C IZMIN  - Minimum z-coordinate of search space 
C N      - Number of bounding boxes
C NDIM2  - Search space dimension
C NX     - Number of cells in x-direction
C NY     - Number of cells in y-direction
C NZ     - Number of cells in z-direction
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      IXMIN=ICOOR(1,1)
      IYMIN=ICOOR(2,1)
      IZMIN=ICOOR(3,1)
      IXMAX=IXMIN
      IYMAX=IYMIN
      IZMAX=IZMIN
C Evaluate minumum and maximum coordinates of region
      DO 10 I=2,N
        IF(ICOOR(1,I).LT.IXMIN)IXMIN=ICOOR(1,I)
        IF(ICOOR(2,I).LT.IYMIN)IYMIN=ICOOR(2,I)
        IF(ICOOR(3,I).LT.IZMIN)IZMIN=ICOOR(3,I)
        IF(ICOOR(1,I).GT.IXMAX)IXMAX=ICOOR(1,I)
        IF(ICOOR(2,I).GT.IYMAX)IYMAX=ICOOR(2,I)
        IF(ICOOR(3,I).GT.IZMAX)IZMAX=ICOOR(3,I)
10    CONTINUE
C Compute number of cells in each direction
      NX=(IXMAX-IXMIN)/ICELSZ(1)+1
      NY=(IYMAX-IYMIN)/ICELSZ(2)+1
      NZ=(IZMAX-IZMIN)/ICELSZ(3)+1
C D     CALL SEXIT(MODEDB)
      RETURN
      END



      SUBROUTINE NBIZL3
     1( IHEADZ     ,ICELSZ     ,ICOOR      ,ITAILZ      ,IXSLIC     ,
     2  IYSLIC     ,MARK       ,MARKZ      ,NEXTZ       ,
     3  IXMIN      ,IYMIN      ,IZMIN      ,N           ,NDIM2      ,
     4  NZ         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IHEADZ(NZ)         ,ICELSZ(3)        ,ICOOR(NDIM2,N)     ,
     2  ITAILZ(NZ)         ,IXSLIC(N)        ,IYSLIC(N)          ,
     3  MARK(N)            ,MARKZ(NZ)        ,NEXTZ(N)             
      DATA NAME/'NBIZL3'/
C***********************************************************************
C*ACRONYM
C Build iz list for 3D NBS contact dtection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IHEADZ - Head of iz list
C ICELSZ - Cell sizes
C ICOOR  - Bounding box coordinates
C ITAILZ - Tail of iz list
C IXSLIC - Integerised x coordinates
C IYSLIC - Integerised y coordinates
C*Variables
C IXMIN  - Minimum x-coordinate of search space 
C IYMIN  - Minimum y-coordinate of search space 
C IZMIN  - Minimum z-coordinate of search space 
C N      - Number of bounding boxes
C NDIM2  - Search space dimension
C NZ     - Number of cells in z-direction
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
C Initialisation
      DO 10 I=1,NZ
              IHEADZ(I)=0
              ITAILZ(I)=0
              MARKZ(I)=0
10    CONTINUE
C Compute integerised coordinates 
      DO 20 I=1,N
        IX=(ICOOR(1,I)-IXMIN)/ICELSZ(1)+1
        IY=(ICOOR(2,I)-IYMIN)/ICELSZ(2)+1
        IZ=(ICOOR(3,I)-IZMIN)/ICELSZ(3)+1
        IXSLIC(I)=IX
        IYSLIC(I)=IY
c     Build iz list
        IF(IHEADZ(IZ).NE.0)THEN
          NEXTZ(I)=IHEADZ(IZ)
        ELSE
          NEXTZ(I)=0
          ITAILZ(IZ)=I
        ENDIF
        IHEADZ(IZ)=I
        MARKZ(IZ)=MARKZ(IZ)+1
        MARK(I)=0
20    CONTINUE
      RETURN
      END



      SUBROUTINE NZCNT3
     1( IHEADX     ,IHEADY     ,IHEADZ     ,ICELSZ     ,ICOOR       ,
     2  IPAIR      ,ITAILZ     ,IXSLIC     ,IYSLIC     ,LSTARC      ,
     3  LSTARP     ,LSTARY     ,MARK       ,MARKX      ,MARKY       ,
     4  MARKZ      ,NEXTX      ,NEXTY      ,NEXTZ      ,
     5  IXMIN      ,IYMIN      ,IZMIN      ,MCNTS      ,MP          ,
     6  N          ,NDIM2      ,NX         ,NY         ,NZ          ,
     7  UPDMTA     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL
     1  UPDMTA
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IHEADX(NX)         ,IHEADY(NY)          ,IHEADZ(NZ)         ,
     2  ICELSZ(3)          ,ICOOR(NDIM2,N)      ,IPAIR(2,MCNTS)     ,
     3  IXSLIC(N)          ,IYSLIC(N)           ,LSTARC(N)          ,
     4  LSTARP(N)          ,LSTARY(N)           ,MARK(N)            ,
     5  ITAILZ(NZ)         ,MARKX(NX)           ,MARKY(NY)          ,
     6  MARKZ(NZ)          ,NEXTX(N)            ,NEXTY(N)           ,
     7  NEXTZ(N)                 
      DATA NAME/'NZCNT3'/
C***********************************************************************
C*ACRONYM
C Contact dtection for 3D NBS algorithm
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IHEADX - Head of ix list
C IHEADY - Head of iy list
C IHEADZ - Head of iz list
C ICELSZ - Cell sizes
C ICOOR  - Bounding box coordinates
C IPAIR  - Contact pair list
C IXSLIC - Integerised x coordinates
C IYSLIC - Integerised y coordinates
C*Variables
C IXMIN  - Minimum x-coordinate of search space 
C IYMIN  - Minimum y-coordinate of search space
C IZMIN  - Minimum z-coordinate of search space  
C MCNTS  - Max number of contacting couples
C MP     - Number of contact pairs
C N      - Number of bounding boxes
C NDIM2  - Search space dimension
C NX     - Number of cells in x-direction
C NY     - Number of cells in y-direction 
C NZ     - Number of cells in z-direction 
C=UPDMTA - Logical flag for extending MCNTS
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
C Sort out non-empty layers in z-direction
      CALL NCOMPS
     1( IHEADZ     ,MARKZ      , 
     2  NZ         ,NNZ        ) 
C     Loop over non-empty z-layers    
      DO 10 II=1,NNZ
C Initialise number of potential contact in previous list
        NP=0
C Get layer number
        IZ=MARKZ(II)
C March boxes with potential contact from previous layers
        IF(II.NE.1)THEN
          MAXIZ=IZMIN+(IZ-1)*ICELSZ(3)
          CALL NMARL3
     1( ICOOR      ,MARK       ,NEXTZ      ,  
     2  3          ,IHEADZ0    ,ITAILZ(IZ) ,MAXIZ      ,N           ,
     3  NDIM2      )
        ENDIF
C Build iy list
        CALL NBIYL3
     1( IHEADX     ,IHEADY     ,IYSLIC     ,MARKX      ,MARKY       ,
     2  NEXTY      ,NEXTZ      ,
     3  IHEADZ(IZ) ,N          ,NX         ,NY         )
C Contact detection in the layer
        CALL NYCNT3
     1( IHEADX     ,IHEADY     ,ICELSZ     ,ICOOR      ,IPAIR       ,
     2  IXSLIC     ,LSTARC     ,LSTARP     ,LSTARY     ,MARK        ,
     3  MARKX      ,MARKY      ,NEXTX      ,NEXTY      ,     
     4  IXMIN      ,IYMIN      ,MCNTS      ,MP         ,N           ,
     5  NDIM2      ,NX         ,NY         ,UPDMTA     )
        IF(UPDMTA)GOTO 100
        IHEADZ0=IHEADZ(IZ)
   10 CONTINUE
  100 CONTINUE
C D     CALL SEXIT(MODEDB)
      RETURN
      END


      SUBROUTINE NCOMPS
     1( IA         ,IB        ,  
     2  N          ,NN         ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IA(N)          ,IB(N)     
      DATA NAME/'NCOMPS'/
C***********************************************************************
C*ACRONYM
C Sort non-zero elements in array ia 
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IA     - Array
C IB     - Array
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      NN=0
      DO 10 I=1,N
        IF(IA(I).NE.0)THEN
          NN=NN+1
          IB(NN)=I
        ENDIF
10    CONTINUE
C D     CALL SEXIT(MODEDB)
      RETURN
      END



      SUBROUTINE NMARL3
     1( ICOOR      ,MARK       ,NEXT       ,  
     2  IDIM       ,IHEAD      ,JTAIL      ,MAXCOO     ,N           ,
     3  NDIM2      )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  ICOOR(NDIM2,N)     ,MARK(N)            ,NEXT(N)             
      DATA NAME/'NMARL3'/
C***********************************************************************
C*ACRONYM
C March boxes with potential contact from previous layers for 3D NBS algorithm
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICOOR  - Bounding box coordinates
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      MDIM=NDIM2/2+IDIM
      I=IHEAD
10    IF(I.NE.0)THEN
        IF(ICOOR(MDIM,I).GE.MAXCOO)THEN
          INEXT=NEXT(I)
            NEXT(JTAIL)=I
            NEXT(I)=0
            JTAIL=I
            MARK(I)=4
            I=INEXT
        ELSE
          I=NEXT(I)
        ENDIF
        GOTO 10
      ENDIF
      RETURN
      END


      SUBROUTINE NBIYL3
     1( IHEADX     ,IHEADY     ,IYSLIC     ,MARKX      ,MARKY       ,
     2  NEXTY      ,NEXTZ      ,
     3  IHEADZ     ,N          ,NX         ,NY         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IHEADX(NX)         ,IHEADY(NY)         ,IYSLIC(N)        ,
     2  MARKX(NX)          ,MARKY(NY)          ,NEXTY(N)         ,
     3  NEXTZ(N)    
      DATA NAME/'NBIYL3'/
C***********************************************************************
C*ACRONYM
C Build iy list for 3D NBS contact dtection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IHEADX - Head of ix list
C IHEADY - Head of iy list
C IYSLIC - Integerised y coordinates
C*Variables
C N      - Number of bounding boxes
C NX     - Number of cells in x-direction
C NY     - Number of cells in y-direction 
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
C Initialisation
      DO 10 I=1,NX
        IHEADX(I)=0
        MARKX(I)=0
10    CONTINUE
      DO 20 I=1,NY
        IHEADY(I)=0
        MARKY(I)=0
20    CONTINUE
C Get box number
      I=IHEADZ
100   IF(I.NE.0)THEN 
C Get row number
        IY=IYSLIC(I)
C Build iy list
        IF(IHEADY(IY).NE.0)THEN
          NEXTY(I)=IHEADY(IY)
        ELSE
          NEXTY(I)=0
        ENDIF
        IHEADY(IY)=I
        MARKY(IY)=1
        I=NEXTZ(I)
        GOTO 100
      ENDIF
      RETURN
      END


      SUBROUTINE NYCNT3
     1( IHEADX     ,IHEADY     ,ICELSZ     ,ICOOR      ,IPAIR       ,
     2  IXSLIC     ,LSTARC     ,LSTARP     ,LSTARY     ,MARK        ,
     3  MARKX      ,MARKY      ,NEXTX      ,NEXTY      ,
     4  IXMIN      ,IYMIN      ,MCNTS      ,MP         ,N           ,
     5  NDIM2      ,NX         ,NY         ,UPDMTA     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL
     1  UPDMTA
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IHEADX(NX)         ,IHEADY(NY)          ,ICELSZ(3)          ,
     2  ICOOR(NDIM2,N)     ,IPAIR(2,MCNTS)      ,IXSLIC(N)          ,
     3  LSTARC(N)          ,LSTARP(N)           ,LSTARY(N)          ,
     4  MARK(N)            ,MARKX(NX)           ,MARKY(NY)          ,
     5  NEXTX(N)           ,NEXTY(N)
      DATA NAME/'NYCNT3'/
C***********************************************************************
C*ACRONYM
C Contact dtection in the layer for 3D NBS algorithm
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IHEADX - Head of ix list
C IHEADY - Head of iy list
C ICELSZ - Cell sizes
C ICOOR  - Bounding box coordinates
C IPAIR  - Contact pair list
C IXSLIC - Integerised x coordinates
C*Variables
C IXMIN  - Minimum x-coordinate of search space 
C IYMIN  - Minimum y-coordinate of search space 
C MCNTS  - Max number of contacting couples
C MP     - Number of contact pairs
C N      - Number of bounding boxes
C NDIM2  - Search space dimension
C NX     - Number of cells in x-direction
C NY     - Number of cells in y-direction 
C=UPDMTA - Logical flag for extending MCNTS
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
C Sort out non-empty rows
      CALL NCOMPS
     1( IHEADY     ,MARKY      ,
     2  NY         ,NNY        ) 
C Loop over non-empty rows
      DO 20 II=1,NNY
C Initialise number of potential contact in previous list
        NP=0
C Get row number
        IY=MARKY(II)
C March boxes with potential contact from previous rows
        IF(II.NE.1)THEN
          MAXIY=IYMIN+(IY-1)*ICELSZ(2)
          CALL NMARR3
     1( LSTARY     ,ICOOR      ,MARK       ,NEXTY      ,      
     2  2          ,IHEADY(IY) ,JK         ,MAXIY      ,N           ,
     3  NDIM2      )
        ENDIF
C Create a temporary array to store boxes in the list
        CALL NLSTAY
     1( LSTARY     ,NEXTY     ,  
     2  IHEADY(IY) ,JK        ,N          )
C Build ix list
        CALL NBDIXL
     1( IHEADX     ,IXSLIC    ,LSTARY    ,NEXTX     ,  
     2  JK         ,N         ,NX        ,NY        ) 
C Sort out non-empty cells
        CALL NCOMPS
     1( IHEADX     ,MARKX      , 
     2  NX         ,NNX        )      
C Loop over non-empty cells                           
        DO 10 JJ=1,NNX
C Get cell number
          IX=MARKX(JJ)
C Create a temporary array to store boxes in the current cell
C     (m:  total number of boxes in the cell
C      mm: number of boxes marched from previous rows)
          CALL NLSTAC
     1( LSTARC     ,MARK      ,NEXTX      ,  
     2  IHEADX(IX) ,M         ,MM         ,N          ) 
C Check overlap within the current cell
          IF(M.GT.1)CALL NOVPC3
     1( ICOOR      ,IPAIR      ,LSTARC     ,MARK       ,
     2  M          ,MCNTS      ,MM         ,MP         ,N          ,
     3  NDIM2      ,UPDMTA     )
          IF(UPDMTA)GOTO 100
C Check overlap with previous cells
          IF(NP.GT.0)CALl NOVPP3
     1( ICOOR      ,IPAIR      ,LSTARC     ,LSTARP     ,MARK        ,
     3  MCNTS      ,M          ,NP         ,MP         ,N           ,
     4  NDIM2      ,UPDMTA     )
          IF(UPDMTA)GOTO 100
C Update potential contact from previous cells
          IF(JJ.NE.NNX)THEN
            MAXIX=IXMIN+(MARKX(JJ+1)-1)*ICELSZ(1)
            CALL NUPDA3
     1( ICOOR      ,LSTARC     ,LSTARP     ,MARK      ,
     2  1          ,M          ,MAXIX      ,N         ,NDIM2     ,
     3  NP         )
          ENDIF
          IHEADX(IX)=0
   10   CONTINUE
   20 CONTINUE
  100 CONTINUE
      RETURN
      END


      SUBROUTINE NMARR3
     1( IARRAY     ,ICOOR      ,MARK       ,NEXT       ,      
     2  IDIM       ,IHEAD      ,M          ,MAXCOO     ,N          ,
     3  NDIM2      )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IARRAY(M)          ,ICOOR(NDIM2,N)     ,MARK(N)            ,
     2  NEXT(N)             
      DATA NAME/'NMARR3'/
C***********************************************************************
C*ACRONYM
C March boxes with potential contact from previous rows for 3D NBS algorithm
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICOOR  - Bounding box coordinates
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
      MDIM=NDIM2/2+IDIM
      DO 10 I=1,M
        J=IARRAY(I)
        IF(ICOOR(MDIM,J).GE.MAXCOO)THEN
          NEXT(J)=IHEAD
          IHEAD=J
          MARK(J)=IBSET(MARK(J), 1)
          MARK(J)=IBCLR(MARK(J), 0)
        ENDIF
10    CONTINUE
      RETURN
      END


      SUBROUTINE NLSTAY
     1( IARRAY     ,NEXT      ,  
     2  IHEAD      ,M         ,N          ) 
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IARRAY(N)          ,NEXT(N)     
      DATA NAME/'NLSTAY'/
C***********************************************************************
C*ACRONYM
C Store elements of list in an arary
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      I=IHEAD
      M=0
10    IF(I.NE.0)THEN 
        M=M+1
        IARRAY(M)=I
        I=NEXT(I)
        GOTO 10
      ENDIF
      RETURN
      END


      SUBROUTINE NBDIXL
     1( IHEADX     ,IXSLIC    ,LSTARY      ,NEXTX      ,
     2  M          ,N         ,NX          ,NY         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IHEADX(NX)     ,IXSLIC(N)     ,LSTARY(M)      ,
     2  NEXTX(N)
      DATA NAME/'NBDIXL'/
C***********************************************************************
C*ACRONYM
C Build ix list for NBS contact dtection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IHEADX - Head of ix list
C IXSLIC - Integerised x coordinates
C LSTARY - Array of iy list
C NEXTX  - Next element of ix list
C*Variables
C M      - Number of boxes in the list
C N      - Number of bounding boxes
C NX     - Number of cells in x-direction
C NY     - Number of cells in y-direction 
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
C Loop over boxes in lstary
      DO 10 IK=1,M
C Get box number
        J=LSTARY(IK)
C Get cell (column) number
        IX=IXSLIC(J)
C Build ix list
        IF(IHEADX(IX).NE.0)THEN
          NEXTX(J)=IHEADX(IX)
        ELSE
          NEXTX(J)=0
        ENDIF
        IHEADX(IX)=J
10    CONTINUE
      RETURN
      END


      SUBROUTINE NLSTAC
     1( IARRAY     ,MARK      ,NEXT       ,  
     2  IHEAD      ,M         ,MM         ,N          ) 
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  IARRAY(N)          ,MARK(N)       ,NEXT(N)     
      DATA NAME/'NLSTAC'/
C***********************************************************************
C*ACRONYM
C Store elements of list in an arary and count number of marked elements
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      I=IHEAD
      M=0
      MM=0
10    IF(I.NE.0)THEN 
        M=M+1
        IARRAY(M)=I
        IF(MARK(I).NE.0)MM=MM+1
        I=NEXT(I)
        GOTO 10
      ENDIF
      RETURN
      END


      SUBROUTINE NOVPC3
     1( ICOOR      ,IPAIR      ,LSTARC     ,MARK       ,
     2  M          ,MCNTS      ,MM         ,MP         ,N          ,
     3  NDIM2      ,UPDMTA     )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      LOGICAL
     1  UPDMTA
      DIMENSION
     1  ICOOR(NDIM2,N)     ,IPAIR(2,MCNTS)    ,LSTARC(M)     ,
     2  MARK(N)          
      DATA NAME/'NOVPC3'/
C***********************************************************************
C*ACRONYM
C Check overlap within the cell for 3D NBS contact detection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IPAIR  - Contact pair list
C ICOOR  - Bounding box coordinates
C LSTARC - Current cell list array
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
C Loop over boxes excluding those marched from previous rows
      DO 20 I=1,M-MM
        II=LSTARC(I)
        DO 10 J=I+1,M
          JJ=LSTARC(J)
          IF(NOVCK3(ICOOR(1,II),ICOOR(1,JJ)).EQ.1)THEN
            MP=MP+1
C Check if MCNTS is big enough
            IF(MP.GT.MCNTS)THEN
              UPDMTA=.TRUE.
              GOTO 100
            ENDIF
            IPAIR(1,MP)=II
            IPAIR(2,MP)=JJ
          ENDIF
10      CONTINUE
20    CONTINUE
      DO 40 I=M-MM+1,M
        II=LSTARC(I)
        MARKI=MARK(II)
        DO 30 J=I+1,M
          JJ=LSTARC(J)
          MARKJ=MARK(JJ)
         IF(MARKI.EQ.0.OR.MARKJ.EQ.0)
     1         PRINT *,'ERROR!', MARKI, MARKJ
            IF(IAND(MARKI,MARKJ).EQ.0)THEN
              IF(NOVCK3(ICOOR(1,II),ICOOR(1,JJ)).EQ.1)THEN
                MP=MP+1
C Check if MCNTS is big enough
                IF(MP.GT.MCNTS)THEN
                  UPDMTA=.TRUE.
                  GOTO 100
                ENDIF
                IPAIR(1,MP)=II
                IPAIR(2,MP)=JJ
              ENDIF
            ENDIF
30      CONTINUE
40    CONTINUE
100   CONTINUE
      RETURN
      END


      FUNCTION NOVCK3
     1( ICOOR      ,JCOOR      )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  ICOOR(6)     ,JCOOR(6)
      DATA NAME/'NOVCK3'/
C***********************************************************************
C*ACRONYM
C Function of overlap check between two bounding boxes
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
      NOVCK3=0
      IF(ICOOR(1).GT.JCOOR(4)) GOTO 900
      IF(ICOOR(4).LT.JCOOR(1)) GOTO 900
      IF(ICOOR(2).GT.JCOOR(5)) GOTO 900
      IF(ICOOR(5).LT.JCOOR(2)) GOTO 900
      IF(ICOOR(3).GT.JCOOR(6)) GOTO 900
      IF(ICOOR(6).LT.JCOOR(3)) GOTO 900
      NOVCK3=1
900   CONTINUE
C D     CALL SEXIT(MODEDB)
      RETURN
      END


	SUBROUTINE NOVPP3
     1( ICOOR      ,IPAIR      ,LSTARC     ,LSTARP     ,MARK        ,
     3  MCNTS      ,MI         ,MJ         ,MP         ,N           ,
     4  NDIM2      ,UPDMTA     )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
	LOGICAL
     1  UPDMTA
      DIMENSION
     1  ICOOR(NDIM2,N)     ,IPAIR(2,MCNTS)    ,LSTARC(MI)           ,
     2  LSTARP(MJ)         ,MARK(N)      
      DATA NAME/'NOVPP3'/
C***********************************************************************
C*ACRONYM
C Check overlap with previous cells for 3D NBS contact detection
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C IPAIR  - Contact pair list
C ICOOR  - Bounding box coordinates
C LSTARC - Current cell list array
C LSTARP - Previous cells list array
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)
	DO 20 I=1,MI
	  II=LSTARC(I)
	  MARKI=MARK(II)
	  DO 10 J=1,MJ
	    JJ=LSTARP(J)
	    MARKJ=MARK(JJ)
	    IF(IAND(MARKI,MARKJ).EQ.0)THEN
	      IF(NOVCK3(ICOOR(1,II),ICOOR(1,JJ)).EQ.1)THEN
	        MP=MP+1
C Check if MCNTS is big enough
	        IF(MP.GT.MCNTS)THEN
                UPDMTA=.TRUE.
                GOTO 100
              ENDIF
	        IPAIR(1,MP)=II
	        IPAIR(2,MP)=JJ
	      ENDIF
	    ENDIF				
10	  CONTINUE
20	CONTINUE
100	CONTINUE
	RETURN
	END

      SUBROUTINE NUPDA3
     1( ICOOR      ,LSTARC     ,LSTARP     ,MARK      ,
     2  IDIM       ,M          ,MAXCOO     ,N          ,NDIM2     ,
     3  NP         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  ICOOR(NDIM2,N)     ,LSTARC(M)      ,LSTARP(N)          ,
     2  MARK(N)
      DATA NAME/'NOVPP3'/
C***********************************************************************
C*ACRONYM
C Update potential contact in current and previous lists for 3D NBS algorithm
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICOOR  - Bounding box coordinates
C LSTARC - Current cell list array
C LSTARP - Previous cells list array
C*Variables
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)                
      MDIM=NDIM2/2+IDIM
C Check boxes in previous list
      MP=0
      DO 10 I=1,NP
        IP=LSTARP(I)
        IF(ICOOR(MDIM,IP).GE.MAXCOO)THEN
          MP=MP+1
          LSTARP(MP)=IP
          MARK(IP)=IBSET(MARK(IP),0)
        ENDIF
10    CONTINUE
C Check boxes in current list
      DO 20 I=1,M
        IC=LSTARC(I)
        IF(ICOOR(MDIM,IC).GE.MAXCOO)THEN
          MP=MP+1
          LSTARP(MP)=IC
          MARK(IC)=IBSET(MARK(IC),0)
        ENDIF
20    CONTINUE
      NP=MP
      RETURN
      END


      SUBROUTINE NOUTCV
     1( ICOOR      ,IPAIR      ,LACCOC     ,LACCOP     ,LACCOT     ,
     2  MP         ,NDIMN      ,NDIM2      ,N          )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1  ICOOR(NDIM2,N)        ,IPAIR(2,MP)      ,LACCOC(N) ,
     2  LACCOP(N)             ,LACCOT(MP)       
      DATA NAME/'NOUTCV'/
C***********************************************************************
C*ACRONYM
C Convert output of NBS contact dtection to ELFEN contact list format
C*HYSTORY
C*Name            Date            Comment
C K. Han          Feb 2003        Initial coding
C*EXTERNAL
C*Arrays
C ICOOR  - Bounding box coordinates
C IPAIR  - Contact pair list
C=LACCOC - List of number of contacted bodies
C=LACCOP - List of pointers to contacted bodies
C=LACCOT - List of targetor numbers
C LNODS  - List of element nodes
C*Variables 
C ICTYPE - Dem contact type
C MFNOD  - Maximum number of segment nodes
C MP     - Number of contact pairs
C NCPOIN - Number of contact points
C NCSEG  - Number of contact segments
C NDIMN  - Problem dimension
C NDIM2  - Search space dimension
C NDISKS - Number of grains
C NN     - NCPOIN+NDISKS
C (c) Copyright 2003, Rockfield Software Limited, Swansea, UK
C***********************************************************************
C D     CALL SENTRY(NAME,MODEDB)        
C     Initialisation
      DO 10 I=1,N
        LACCOC(I)=0
        LACCOP(I)=0
10    CONTINUE
      DO 20 I=1,MP
        LACCOT(I)=0
20    CONTINUE
C Calculate number of targets for each contactor
      DO 30 I=1,MP
        ICON=IPAIR(1,I)
        ITAR=IPAIR(2,I)
        IF(ICON.GT.ITAR)THEN
          LACCOC(ICON)=LACCOC(ICON)+1
        ELSE
          LACCOC(ITAR)=LACCOC(ITAR)+1
        ENDIF
30    CONTINUE
C Pointers to target list
      LACCOP(1)=1
      DO 40 I=2,N
        LACCOP(I)=LACCOP(I-1)+LACCOC(I-1)
40    CONTINUE
C Build up target list
      DO 50 I=1,MP
        ICON=IPAIR(1,I)
        ITAR=IPAIR(2,I)
        IF(ICON.GT.ITAR)THEN
          LACCOT(LACCOP(ICON))=ITAR
          LACCOP(ICON)=LACCOP(ICON)+1
        ELSE
          LACCOT(LACCOP(ITAR))=ICON
          LACCOP(ITAR)=LACCOP(ITAR)+1
        ENDIF
50    CONTINUE
C     Rewind the pointer
      LACCOP(1)=1
      DO 60 I=2,N
        LACCOP(I)=LACCOP(I-1)+LACCOC(I-1)
60    CONTINUE
      RETURN
      END

      