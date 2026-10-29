!DEC$ IF .NOT. DEFINED (HGTL_get_traction_radius_f)
!DEC$ DEFINE HGTL_get_traction_radius_f
      subroutine HGTL_get_traction_radius(inlq,sig,tnormal,xCoord)
c Return the scalar normal traction tnormal
c  and the x-coordinate xCoord of the element centre
c INPUT:
c - "inlq": index of current quadrature point (between lft and llt)
c - "sig": components of Cauchy stress tensor (ordering: xx,yy,zz,xy,yz,zx)
c OUTPUT:
c - "tnormal": scalar normal traction computed from normal vector of
c               longest element edge and stress tensor
c - "xCoord": x-coordinate (radius) of element centre
c USAGE:
c - "call HGTL_get_traction_radius(inlq,sig,tnormal,xCoord)"
c   integer :: inlq ! loop index "i" from loop over quadrature points "do i=lft,llt"
c   real*8, dimension(6) :: sig
c   real*8 :: tnormal, xCoord
c NOTE:
c - compatible with Linux and Windows
c - tested under:
c       * 9_ls-dyna_smp_d_R14_1_0_x64_centos79_ifort190_sse2
c - implemented for Release R10.2 and newer
c   (Windows: older releases require different accessor method for "dm_rots";
c    Linux: for older releases such as R9 this approach is currently not available)
c - hardcoded for dimension=2D and 4-noded bilinear elements
c
c ToDo:
c - measure computation time and improve performance
c - ensure correct usage based on element type 'etype' e.g. for axisym elements
c - extend to 3D
c
      use userinterface ! additionally necessary for use of getreal8ptr
c
      ! Include "nlqparm" to get the size "nlq"
       include 'nlqparm'
c
c Retrieve the index map for the current element
C_TASKCOMMON (aux33loc)
      common/aux33loc/ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq)
c
      integer, intent(in) :: inlq ! loop index "i" from loop over quadrature points "do i=lft,llt"
      real*8, dimension(6), intent(in) :: sig
      real*8, intent(out) :: tnormal, xCoord
      real*8 :: x1,y1,x2,y2,x3,y3,
     &          edgeLength21,edgeLength32,nx,ny,nLength
c
      ! retrieve the initial nodal coordinates "dm_rots"
      ! from module "mod_nodes"
      ! "dm_rots" is a vector:
      !  shape: number of elements times dimension
      !  entries: for instance: ( x1, y1, z1, x2, y2, z2, ... ) for 3D   
       real*8, CONTIGUOUS, POINTER :: dm_rots(:) ! [POINTER_CONTIGUOUS]
c$omp threadprivate (/aux33loc/)
       call getreal8ptr('dm_rots',dm_rots)
c
      ! Compute the centre x-coordinate (radius) of the present element
       xCoord = 0.25*( dm_rots((ix1(inlq)-1)*3+1)
     &                 + dm_rots((ix2(inlq)-1)*3+1)
     &                 + dm_rots((ix3(inlq)-1)*3+1)
     &                 + dm_rots((ix4(inlq)-1)*3+1) )
c
      ! Compute normal vector of surface (assuming a rectangular HGTL element)
        ! Determine the nodal coordinates
         x1 = dm_rots((ix1(inlq)-1)*3+1)
         y1 = dm_rots((ix1(inlq)-1)*3+1 +1)
         x2 = dm_rots((ix2(inlq)-1)*3+1)
         y2 = dm_rots((ix2(inlq)-1)*3+1 +1)
         x3 = dm_rots((ix3(inlq)-1)*3+1)
         y3 = dm_rots((ix3(inlq)-1)*3+1 +1)
        ! Find the longest edge
         edgeLength21 = sqrt( (x2-x1)**2 + (y2-y1)**2 )
         edgeLength32 = sqrt( (x3-x2)**2 + (y3-y2)**2 )
         !write(*,*) "edgeLengths=",edgeLength21,edgeLength32
         if ( edgeLength32 > edgeLength21 ) then
          ! If the edge 2-3 is the longest edge, we overwrite the first point with third point, so we can subsequently always use the 1-2 edge
           x1 = x3
           y1 = y3
         endif
         !write(*,*) "coords x1,y1,x2,y2=",x1,y1,x2,y2
         ! Compute the unit normal vector to edge 1-2
          ! For the scalar normal traction component the normal vector is sufficient (not necessarily the outward pointing)
          ! @note indices cross (2-1 for y vs 1-2 for x)
           nx = y2-y1
           ny = x1-x2
           nLength = sqrt( nx**2 + ny**2 )
           nx = nx / nLength
           ny = ny / nLength
           !write(*,*) "nx,ny=",nx,ny

         ! Compute the scalar normal traction component tn = n_i * sigma_ij * n_j
         ! @note n_3 = 0, sigma_23 = sigma_31 = 0 (sig=xx,yy,zz,xy,yz,zx)
          tnormal = nx*sig(1)*nx + ny*sig(2)*ny + 2.*nx*sig(4)*ny
          !write(*,*) "tnormal=",tnormal,"for sig=",sig(1),sig(2),sig(4)
c
      end subroutine HGTL_get_traction_radius
!DEC$ ENDIF
