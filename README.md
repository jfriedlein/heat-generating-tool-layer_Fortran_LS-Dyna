# heat-generating-tool-layer_Fortran_LS-Dyna
User-defined subroutines to model frictional heat generation by a substitution model that applies user-defined material and thermal subroutines


intro text based on publication

## Changes and extensions to object files

First, the overarching dyn21umat.f must be extended to transmit the index i of the current quadrature point to the desired umat subroutine, here umat41.

[dyn21umat.f]:
```fortran
...
subroutine urmats(...
...
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,eltype,tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject,i)
...
```   

The actual umat subroutine is next extended by calling the HGTL_get_traction_radius subroutine. This call should be placed after updating the stress "sig" to input the current stress. HGTL_get_traction_radius returns the normal traction and x-coordinate of the element centre, which are last saved exemplarily into the first (hsv(1)) and second (hsv(2)) entry of the list of mechanical history variables hsv.

[dyn21umats.f]:
```fortran
      include "heat-generating tool layer_Fortran_LS-Dyna/
     &HGTL_subroutine.f"
c
      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject,inlq)
...
      real*8 :: xCoord,tnormal
...
c Computation/update of stress "sig"
...
          ! Call the subroutine HGTL_get_traction_radius from "HGTL_subroutine.f"
          !  to return the scalar normal traction tnormal and
          !  the x-coordinate xCoord of the element centre
          ! @note Currently only valid for 2D axisymmetric elements
           call HGTL_get_traction_radius(inlq,sig,tnormal,xCoord)

          ! Store the surface normal stress component as history variable to access in thumat
           hsv(1) = tnormal
          ! Store the x-coordinate (radius) as history variable to access in thumat
           hsv(2) = xCoord
...
      return
      end
```      

The thermal subroutine is extended to retrieve the material parameters and mechanical history variables. This is used to compute the possibly temperature-/pressure-/velocity-dependent friction coefficient and eventually compute the frictional heat. Be aware of the needed correction factor "1. / HGTL_height" to compensate for the surface vs volume heat source/integration.

[dyn21tumat.f]:
```fortran
      subroutine thumat11(c1,c2,c3,cvl,dcvdtl,hsrcl,dhsrcdtl,
     1 hsv,hsvm,nmecon,r_matp,crv,nnpcrv,
     2 nel,nep,iep,eltype,dt,atime,ihsrcl,temp,ttconv,ttnext,
     3 tmconv,tmnext,inithis)
c
      include 'nlqparm'
      include 'iounits.inc'
c
      character*(*) eltype
      logical inithis
      dimension hsv(*),hsvm(*),r_matp(*),crv(lq1,2,*)
      integer nnpcrv(*)
      real*8 rotationSpeed, HGTL_height, contactTractionThreshold
      real*8 frictionCoeff_mu0, TzeroF, fricExpM, fricMin,
     &       frictionEvolutionSign, frictionCoeff
      real*8 tnormal, rCoord
      PARAMETER (PI = 3.141592653589793D0)
c
c     Thermal usermaterial number 11.
c
      ! Assign specific heat capacity from thermal material parameters (first parameter)
      cvl = r_matp(1 +8)
      dcvdtl = 0
      ! Assign thermal conductivity from thermal material parameters (second parameter)
      c1 =  r_matp(2 +8)
c
      ! rotation speed n in 1/s (omega=2 pi n)
       rotationSpeed = r_matp(3 +8) 
      ! height/thickness of heat-generating tool layer
       HGTL_height = r_matp(4 +8)
      ! tractioin threshold to eliminate tension and small false compressive tractions
       contactTractionThreshold = r_matp(5 +8)
c
      ! Extract material parameters for friction model
       frictionCoeff_mu0 = r_matp(9 +8) 
       TzeroF = r_matp(10 +8)
       fricExpM = r_matp(11 +8)
       fricMin = r_matp(12 +8)
       frictionEvolutionSign = r_matp(13 +8)
c
      ! Extract mechanical history variables
       ! Scalar contact normal traction
        tnormal = hsvm(1 +7)
       ! Radial coordinate of element centre
        rCoord  = hsvm(2 +7)
c
      ! Set default values for heat source (off)
       ihsrcl = 0
       dhsrcdtl = 0.
      ! Check normal traction for contact (negative/compressive normal traction)
       if ( tnormal < contactTractionThreshold ) then
            ! Compute friction coefficient based on exemplary friction model
             frictionCoeff = frictionCoeff_mu0
     &                       * ( 1. +
     &                           frictionEvolutionSign *
     &                           (
     &                             min(
     &                             max(
     &                                  (temp-293.15) / (TzeroF-293.15)
     &                                   , 0.0
     &                                )
     &                                 , 1.0
     &                                )
     &                           )
     &                            ** fricExpM
     &                         )
             frictionCoeff = max( fricMin, frictionCoeff )
c             
            ! Activate heat source term
             ihsrcl = 1
            ! Compute friction-induced heating as heat source
             hsrcl = frictionCoeff * abs(tnormal) ! frictional traction, "abs(*)": heat should be positive, but for contact tnormal<0
     &               * rotationSpeed * 2. * PI * rCoord ! rotational velocity (mm/s)
     &               / HGTL_height ! remove effect of HGTL height (surface vs volume integral)
                     !* 2*pi*rCoord ! for axisymmetric integration is already done by LS-Dyna internal assembly
c            
            ! Compute linearisation  of heat source term with respect to temperature
             if ( frictionCoeff > fricMin+1e-8
     &            .and. frictionCoeff < frictionCoeff_mu0-1e-8 ) then
               dhsrcdtl = frictionEvolutionSign * abs(tnormal)
     &                 * rotationSpeed * 2. * PI * rCoord
     &                 / HGTL_height
     &                 * 1.
     &                 * frictionCoeff_mu0 * fricExpM
     &                 * ( (temp-293.15) / (TzeroF-293.15) )
     &                   **(fricExpM-1.)
     &                 * 1./(TzeroF-293.15)
             endif
       endif
c
      return
      end
```        


## Keywords

```    
*KEYWORD
$
$ ########################## MATERIAL PARAMETERS ###############################
$
$ Avoid using reduced integration for HGTL layer since e.g. hourglassing can cause strong oscillating normal stress fields, which is regarded as contact in the thermal model
$
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R DENSITY   1.435E-8
$
$ ELASTICITY: ------
$ Young's modulus
R YOUNGMOD    590000
$ Poisson's ratio
R POISSON       0.22
$ Bulk modulus
R BulkMod  351190.48
$ Shear modulus
R ShearMod 241803.28
$
$ THERMAL PROPERTIES OF HGTL:
$ specific heat capacity: (assumed zero)
R heatCap          0
$ heat conductivity: (chosen small)
R heatCond         1
$
$ FRICTION MODEL AND PARAMETERS:
$ friction coefficient between rod and tool: (needs to be identified)
R fricCoef      0.81
$ temperature in K at zero friction: (needs to be identified, typically close to melting temperature)
R TzeroF     1773.15
$ exponent m of friction model: (needs to be identified)
R fricExpM      2.09
$ minimum friction coefficient:
R fricMin        0.0
$ temperature-dependent friction evolution (if changing: -1.0, else 0.0)
R fricEvol      -1.0
$
$ PROCESS PARAMETERS:
$ rotation speed n of rod in 1/s:
R rotSpeed   66.6666
$
$ NUMERICAL PARAMETERS:
$ stress threshold to turn on contact only for a minimum amount of normal traction
$  (choose negative value to trigger for some minium compression)
R strssThr      -1.0
$
R emiWp          0.5
$
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         2         2         1         0         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0STRESS              
$#    dtdt    resplt     neipb   quadsld    cubsld   deleres
         0         0         0         0         0         0
$
$ ###################### USER MATERIAL MODEL CARD ##############################
$
*MAT_USER_DEFINED_MATERIAL_MODELS_TITLE
UMAT - elasticity with history generation for HGTL
$#     mid        ro        mt       lmc       nhv    iortho     ibulk        ig
        11  &density        41         8         2         0         3         4
$#   ivect     ifail    itherm    ihyper      ieos      lmca       ext     epshv
         0         0         1         0         0         0         0         0
$ ELASTICITY:
$ YoungMod  PoissonR   BulkMod  ShearMod    unused    unused    unused    unused
$#      p1        p2        p3        p4        p5        p6        p7        p8
 &youngmod  &poisson  &bulkmod &shearmod         0         0         0         0
$
*MAT_THERMAL_USER_DEFINED_TITLE
THUMAT - frictional heat generation
$#    tmid        ro        mt       lmc       nhv      aopt    iortho      ihve
        11       0.0        11        16         1       0.0         0         1
$  heatCap  heatCond  rotSpeed HGTL_hgth  strssThr
$#      p1        p2        p3        p4        p5        p6        p7        p8
  &heatcap &heatcond &rotspeed   &hgtl_h &strssthr       0.0       0.0       0.0
$ fricCoef    TzeroF  fricExpM   fricMin  fricEvol
$#      p1        p2        p3        p4        p5        p6        p7        p8
 &friccoef   &tzerof &fricexpm  &fricmin &fricevol       0.0       0.0       0.0
$
*INITIAL_TEMPERATURE_SET
$#    nsid      temp       loc
         0    293.15         0
*END
```    


