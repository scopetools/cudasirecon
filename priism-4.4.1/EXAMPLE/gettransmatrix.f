c     Computes the general 2D transformation matrix (rotation,
c     translation, magnification, skewing) for interpo.c.  The
c     rotation from the old axis to NEW axis is positive for an
c     anti-clockwise rotation.
      subroutine gettransmatrix(dx, dy, rotx, angle, dxp, dyp, amat2)
      DIMENSION AMAT(2,2), AMAT2(2,2)
      DIMENSION ROT(2,2), T1(2,2), T2(2,2)
      SAVE CNV
      DATA CNV/ .017453291/

c
c  SET UP CONVERION MATRICES
c
      AXANG = angle
      AROTX = rotx
      AROTX = CNV*AROTX
      AXANG = AXANG*CNV
      GAMMA = angle*CNV
      CR = COS(AROTX)
      SR = SIN(AROTX)
      CA = COS(AXANG)
      SA = SIN(AXANG)
      CG = COS(GAMMA)
      SG = SIN(GAMMA)
c
      ROT(1,1) = CR
      ROT(1,2) = SR
      ROT(2,1) = -SR
      ROT(2,2) = CR
      T1(1,1) = DX
      T1(1,2) = DY*CG
      T1(2,1) = 0.0
      T1(2,2) = DY*SG
c
      T2(1,1) = 1./DXP
      T2(1,2) = CA/(DXP*SA)
      if( T2(1,2) .lt. .0000001) T2(1,2) = 0.0
      T2(2,1) = 0.0
      T2(2,2) = 1./(DYP*SA)
c     
      CALL MM(AMAT,ROT,T1)
      CALL MM(AMAT2,T2,AMAT)
      return
      end

c
c*MM  MULTIPLY 2X2 MATRICES
c
      SUBROUTINE MM(CMAT,AMAT,BMAT)
      DIMENSION AMAT(2,2),BMAT(2,2),CMAT(2,2)
c
      CMAT(1,1) = 0
      CMAT(1,2) = 0
      CMAT(2,1) = 0
      CMAT(2,2) = 0
      CMAT(1,1) = AMAT(1,1)*BMAT(1,1) + AMAT(1,2)*BMAT(2,1)
      CMAT(1,2) = AMAT(1,1)*BMAT(1,2) + AMAT(1,2)*BMAT(2,2)
      CMAT(2,1) = AMAT(2,1)*BMAT(1,1) + AMAT(2,2)*BMAT(2,1)
      CMAT(2,2) = AMAT(2,1)*BMAT(1,2) + AMAT(2,2)*BMAT(2,2)
      RETURN
      END
