C****************************************************************************
C       This is a user material subroutine based on an isotropic,
C             rate-independent, gradient plasticity model 
C                to be used with ABAQUS/Explicit
C
C       LAST MODIFIED  Jan 5, 2018  BY YIN ZHANG
C
C****************************************************************************

      subroutine vuhard(
C Read only -
     *     nblock, 
     *     jElem, kIntPt, kLayer, kSecPt, 
     *     lAnneal, stepTime, totalTime, dt, cmname,
     *     nstatev, nfieldv, nprops, 
     *     props, tempOld, tempNew, fieldOld, fieldNew,
     *     stateOld,
     *     eqps, eqpsRate,
C Write only -
     *     yield, dyieldDtemp, dyieldDeqps,
     *     stateNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), tempOld(nblock), tempNew(nblock),
     1   fieldOld(nblock,nfieldv), fieldNew(nblock,nfieldv),
     2   stateOld(nblock,nstatev), eqps(nblock), eqpsRate(nblock),
     3   yield(nblock), dyieldDtemp(nblock), dyieldDeqps(nblock,2),
     4   stateNew(nblock,nstatev), jElem(nblock)
C
      character*80 cmname
      parameter( NOEL=80)
      real*8    XYZ0(NOEL,8,3), pl_strain(NOEL,8)
      real*8    coords(3,8), coords_rh(3,8), p(8), h
      real*8    dpdx(NOEL,3) 
      common    /grad/ XYZ0, pl_strain, dpdx
      real*8    l, ep_c, n, PI, e_grad, wavelength, SY, xs, ys
      PI    = 4.D0*DATAN(1.D0)
      l     = PROPS(1)
      ep_c  = PROPS(2)
      n     = PROPS(3)
      wavelength = PROPS(4)
C
      do 100 km = 1,nblock
        e_grad = dsqrt(dpdx(jElem(km),1)*dpdx(jElem(km),1)*0+
     1                dpdx(jElem(km),2)*dpdx(jElem(km),2)+
     2                dpdx(jElem(km),3)*dpdx(jElem(km),3)*0)
        yield(km) = 425.D6
C        dyieldDeqps(km,1) = 0
C       dyieldDeqps(km,1)=450.D6*dsqrt(1+l*e_grad/(1+(eqps(km)/ep_c)**n))
       dyieldDeqps(km,1)=2000.D6*
     1 (1.0+(l*e_grad)**(0.5)/(1+(eqps(km)/ep_c)**n))/
     2 (1.0+(eqps(km)/0.015)**0.6)

        dyieldDeqps(km,2)=0


        if(stepTime.LE.1.0D-2) then
            xs    = 0.5-DABS(XYZ0(jElem(km),kIntPt,2)/wavelength - 0.5)
C            SY    = 425D6 + (1060*xs*xs*xs-180*xs*xs-600*xs)*1D6
            SY    = 446D6 - 446*xs*1D6

            stateNew(km,1) = SY
        else
            SY    = stateOld(km,1) +
     1                   dyieldDeqps(km,1)*eqpsRate(km)*dt
        end if

        yield(km) = SY
        stateNew(km,1) = SY
        stateNew(km,2) = dyieldDeqps(km,1)
        stateNew(km,3) = e_grad
  100 continue
C
      return
      end


C ------------------------------------------------------------------
C This subroutine is used to read info at the end of an increment
C only works for C3D8 element with 8 int points.
C ------------------------------------------------------------------
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElem, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
c
      include 'vaba_param.inc'
c
      dimension jElem(nblock), coordMp(nblock,*), 
     *          direct(nblock,3,3), T(nblock,3,3), 
     *          charLength(nblock), props(nprops), 
     *          stateOld(nblock,nstatev), 
     *          stateNew(nblock,nstatev),
     *          field(nblock,nfieldv)
      character*80 cmname
c
c     Local arrays from vgetvrm are dimensioned to 
c     maximum block size (maxblk)
c
      parameter( nrData=6 )
      character*3 cData(maxblk*nrData)
      dimension   rData(maxblk*nrData), jData(maxblk*nrData)
      parameter( NOEL=80)
      real*8    XYZ0(NOEL,8,3), XYZ_t(NOEL,8,3), pl_strain(NOEL,8)
      real*8    p(8), dpdx(NOEL,3)
      common    /grad/ XYZ_t, pl_strain, dpdx
      real*8    coords(3,8), temp(3,4), coords_rh(3,8), h
C
      jStatus = 1
      call vgetvrm( 'PEEQ', rData, jData, cData, jStatus )
c
      if( jStatus .ne. 0 ) then
         call xplb_abqerr(-2,'Utility routine VGETVRM '//
     *      'failed to get variable.',0,zero,' ')
         call xplb_exit
      end if


c
      do 100 K=1,nblock
        pl_strain(jElem(K), kIntPt) = rData(K)
        XYZ0(jElem(K),kIntPt,1) = coordMp(K,1)
        XYZ0(jElem(K),kIntPt,2) = coordMp(K,2)
        XYZ0(jElem(K),kIntPt,3) = coordMp(K,3) 
        XYZ_t(jElem(K),kIntPt,1) = coordMp(K,1)
        XYZ_t(jElem(K),kIntPt,2) = coordMp(K,2)
        XYZ_t(jElem(K),kIntPt,3) = coordMp(K,3) 

        if (kIntPt.EQ.8) then
            do 13 K1=1,8
                p(K1) = pl_strain(jElem(K), K1)
                coords(:,K1) = XYZ_t(jElem(K),K1,:)
  13        continue
            temp(:,1)   = coords(:,3) 
            coords(:,3) = coords(:,4)
            coords(:,4) = temp(:,1)

            temp(:,1)   = coords(:,7)          
            coords(:,7) = coords(:,8)
            coords(:,8) = temp(:,1)

            temp(1,1)   = p(3)
            p(3)        = p(4)
            p(4)        = temp(1,1)

            temp(1,1)   = p(7)
            p(7)        = p(8)
            p(8)        = temp(1,1)

            call right_hand(coords_rh, coords, h, jElem(K))
            call strain_grad(p, dpdx(jElem(K),:), coords_rh)

            coords(:,:) = 0.0
            p(:)        = 0.0
        end if

        do 99 K2=1,nstatev
            stateNew(K,K2) = stateOld(K,K2)
  99    continue 

  100 continue

      return
      end

C --------------------------------------------------------------------
C Cross product function
C --------------------------------------------------------------------
        subroutine cross(a,b,m)
        implicit real*8 (a-h,o-z)
        real*8 a(3), b(3)
        real*8 m(3)

        m(1) = a(2)*b(3) - a(3)*b(2)
        m(2) = a(3)*b(1) - a(1)*b(3)
        m(3) = a(1)*b(2) - a(2)*b(1)
        return
        end subroutine
C --------------------------------------------------------------------
C Convert the coordinate system into right handed system.
C --------------------------------------------------------------------
        subroutine right_hand(coords_rh, coords, h, kelem)
        implicit real*8 (a-h,o-z)

        real*8 coords(3,8), coords_rh(3,8)
        integer i,j, kelem
        real*8 x(3), y(3), z(3), m(3), h

        coords_rh = coords

        x = coords(:,2) - coords(:,1)
        y = coords(:,4) - coords(:,1)
        z = coords(:,5) - coords(:,1)

        Call cross(x,y,m) 
        h = dot_product(z,m)

        if (h<=0)  then
            print*, 'This is not a right-handed system.', kelem
            coords_rh(:,1:4) = coords(:,5:8)
            coords_rh(:,5:8) = coords(:,1:4)
        endif

        end subroutine


C --------------------------------------------------------------------
C Evaluate the strain gradient dpdx_int at the centor of the integration pts.
C --------------------------------------------------------------------
        subroutine strain_grad(p,dpdx_int,coords)
        implicit real*8 (a-h,o-z)
        real*8 p(8)
        real*8 N(8)
        real*8 dpdx_int(1,3)
        real*8 coords(3,8)
        real*8 dNdx(3,8)
        real*8 dNdx_parr(3,8)
        real*8 pt_parr(3)
        real*8 jacob(3,3), inv_jacob(3,3)
        real*8 J
        integer i

        pt_parr(:)  = 0.0 
        dpdx_int(:,:)   = 0.0
 
        Call Shape3D8N_parrent(dNdx_parr, N, pt_parr)

        Call GradientShape3D8N(dNdx, dNdx_parr, coords, J,
     $       jacob, inv_jacob)


        do i = 1,8
            dpdx_int(1,1)  = dpdx_int(1,1) + dNdx(1,i)*p(i)
            dpdx_int(1,2)  = dpdx_int(1,2) + dNdx(2,i)*p(i)
            dpdx_int(1,3)  = dpdx_int(1,3) + dNdx(3,i)*p(i)
        end do
        
      end
C --------------------------------------------------------------------
C It returns the gradient of shape functions in current
C coordinates (dN/dx) and the determination of jacobian.
C --------------------------------------------------------------------
       SUBROUTINE GradientShape3D8N(dNdx, dNdx_parr, coords_elem, 
     $ J, jacob, inv_jacob)
        implicit real*8 (a-h,o-z)
    
        real*8 dNdx_parr(3,8)
        real*8 coords_elem(3,8)
        real*8 J                 
        real*8 dNdx(3,8)
        integer i
        real*8 jacob(3,3), inv_jacob(3,3)
        
        jacob(:,:)= 0.0
        jacob = matmul(coords_elem,transpose(dNdx_parr))
    
        CALL InverseMatrix(inv_jacob, jacob, J)
    
        dNdx(:,:) = 0
        dNdx = transpose(matmul(transpose(dNdx_parr),inv_jacob))
    
        return
        END SUBROUTINE
C ---------------------------------------------------------------
C Shape3D8Np_parrent returns the shape functions as well as gradient of shape functions
C of 3D8nodes elements for a given point. The coordinates are all in parrent
C coordinates.
C The dimensions of dNdx_parr(3,8) implies the meaning of the two axises. 3
C means 3 dimensions, 8 means 8 nodes.
C
C      8/------/|7      
C     5|------|6|
C      | 4    | /3
C	   |      |/
C     1-------/2
C            
C From Zhanli:
C  z(xi)/|\ /y(eta)
C        | /
C        |/------>x(epu)   
C
C x(epu)/|\ 
C        | 
C <------|y(eta)   
C       /
C      /
C     /z(xi)
C ---------------------------------------------------------------
        SUBROUTINE Shape3D8N_parrent(dNdx_parr, N, pt_parr)
        implicit real*8 (a-h,o-z)
        real*8 pt_parr(3)
        real*8 N(8)
        real*8 dNdx_parr(3,8)
    
        integer i
        real*8 epu, eta, xi
        real*8 epu0(8), eta0(8), xi0(8)
        
        epu=pt_parr(1)
        eta=pt_parr(2)
        xi =pt_parr(3)
        epu0 = (/-1, 1, 1, -1, -1, 1, 1, -1/)
        eta0 = (/-1, -1, 1, 1, -1, -1, 1, 1/)
        xi0  = (/-1, -1, -1, -1, 1, 1, 1, 1/)
    
        do i = 1,8
        N(i)     = 1/8.0d0*(1+epu0(i)*epu)*(1+eta0(i)*eta)*(1+xi0(i)*xi)
        dNdx_parr(1,i)=1/8.0d0*(epu0(i))*(1+eta0(i)*eta)*(1+xi0(i)*xi)
        dNdx_parr(2,i)=1/8.0d0*(1+epu0(i)*epu)*(eta0(i))*(1+xi0(i)*xi)
        dNdx_parr(3,i)=1/8.0d0*(1+epu0(i)*epu)*(1+eta0(i)*eta)*(xi0(i))
        end do
        return
        END SUBROUTINE
C --------------------------------------------------------------------------
C Inverse a 3 by 3 matrix.
C --------------------------------------------------------------------------
        SUBROUTINE InverseMatrix(inv_M, M, J)
        implicit real*8 (a-h,o-z)
        real*8 M(3,3), inv_M(3,3)
        real*8 J
        integer k,i
        
        inv_M(1:3,1:3) = 0

        inv_M(1,1) =  M(2,2)*M(3,3)-M(2,3)*M(3,2)
        inv_M(2,1) = -M(2,1)*M(3,3)+M(2,3)*M(3,1)
        inv_M(3,1) =  M(2,1)*M(3,2)-M(2,2)*M(3,1)
        
        inv_M(1,2) = -M(1,2)*M(3,3)+M(1,3)*M(3,2)
        inv_M(2,2) =  M(1,1)*M(3,3)-M(1,3)*M(3,1)
        inv_M(3,2) = -M(1,1)*M(3,2)+M(1,2)*M(3,1)
        
        inv_M(1,3) =  M(1,2)*M(2,3)-M(2,2)*M(1,3)
        inv_M(2,3) = -M(1,1)*M(2,3)+M(2,1)*M(1,3)
        inv_M(3,3) =  M(1,1)*M(2,2)-M(2,1)*M(1,2)
        
        J = M(1,1)*inv_M(1,1)+M(1,2)*inv_M(2,1)+M(1,3)*inv_M(3,1)

        inv_M = inv_M/J
        return
        END SUBROUTINE

