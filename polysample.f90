module polysample

  use file_functions, only: stderr
  use precision

  implicit none

  private

  character(len=*), parameter :: version = "$Id: polysample.f90,v 1.8 2005/06/23 23:23:31 aidan Exp $"

  !! This is a shameless reproduction of the Java class PolySamp developed by 
  !! Tom McGlynn, NASA/GSFC (http://skyview.gsfc.nasa.gov/polysamp.html).
  !!
  !! Here is the description from the original java class:
  !!
  !! ---------------------------------------------------------------------------
  !! The class implements a fast flux conserving resampling based on the 
  !! Sutherland-Hodgman clipping algorithm.
  !!
  !! Consider an original image of data in some projection and a new map 
  !! of an overlapping region with some projection, coordinate system, etc.  
  !! Assume that the flux within a given pixel in the original image is constant
  !! over the area of the pixel. The flux within each pixel in the new map 
  !! should be the integral of the flux in the region occupied by each pixel 
  !! in the map.
  !! 
  !! * Find all of the corners of the pixels in resampled map and project the 
  !!   positions of these corners into the projection of the original image. 
  !!   These should be scaled such that the coordinates run from 0 to mx and my
  !!   where mx and my are the dimensions of the original image. I.e., in these 
  !!   coordinate each pixel in the original image is a unit square.  This step 
  !!   is done prior to call the primary polySamp methods in this class.  Note
  !!   that the corners of the pixels are required rather than the pixel centers.
  !! 
  !! * For each pixel in the resampled map, find a bounding box in the original 
  !!   image's coordinates.  This is used to find a rectangle of candidate pixels 
  !!   in the original image that may contribute flux to this pixel in the 
  !!   resampled map.
  !! 
  !! * For each candidate image pixel clip the resampled pixel to the image pixels 
  !!   boundaries.
  !! 
  !! * Calculate the area of the clipped region.  This is easy since the clipped 
  !!   region is a convex polygon and we have the vertices. Triangulating the 
  !!   polygon allows us to calculate its area.
  !! 
  !! * Add a flux to the resampled pixel equal to the area of the clipped resampled 
  !!   pixel times the flux in the original map pixel
  !! 
  !! * Repeat for all candidate original pixels
  !! 
  !! * Go to the next resampling pixel.
  !!
  !! The instance methods of this class are not thread-safe, however it is possible 
  !! to generate a separate PolySamp object for each thread to resample the same 
  !! input image.
  !!
  !! Developed by Tom McGlynn, NASA/GSFC
  !! October 3, 2002
  !! ---------------------------------------------------------------------------

  !! $Log: polysample.f90,v $
  !! Revision 1.8  2005/06/23 23:23:31  aidan
  !! Fixed a simple error in the expression calculating the covered fraction
  !! of a pixel -- was calculating the fraction covered by *missing* pixels. Ooops.
  !!
  !! Removed the hard limit, of one pixel, on minimum coverage. This really
  !! only makes sense when the output pixel is of the same order in size
  !! as the input. Where output >> input this limit would be silly.
  !!
  !! Revision 1.7  2005/06/15 01:22:15  aidan
  !! Changed to parameterised real kinds from the precision module. Also made
  !! the min_coverage variable single precision.
  !!
  !! Revision 1.6  2005/06/06 04:57:17  aidan
  !! Ooops! Version 1.5 didn't compile. Forgot to include the module which
  !! defines stderr. Fixed.
  !!
  !! Revision 1.5  2005/06/02 00:08:47  aidan
  !! Added support for a "minimum coverage" option. If the coverage option is passed
  !! to polysamp it specifies the minimum area of the final pixel that must have
  !! non-missing data. Below this fraction the entire pixel will be set to the missing
  !! value. There is also a hard limit of one pixel which, if met, means the output
  !! pixel will have a real value, regardless of the specified coverage fraction.
  !!
  !! Fixed a logical flaw when detecting pixels that are entirely off image.
  !!
  !! Revision 1.4  2004/12/07 05:12:44  aidan
  !! Added ability to specify intensive polysampling.
  !!
  !! Revision 1.3  2004/10/21 04:57:29  aidan
  !! Fixed a bug with totally bounded pixels. We were indexing the original data
  !! with the lower bound pixel (min) rather than the upper bound pixel (max). A
  !! few moer debug statements were added and some spurious semicolons left over
  !! from the Java translation were deleted.
  !!
  !! Revision 1.2  2004/08/02 07:28:22  aidan
  !! Added missing value support. If a missing value is not supplied then the
  !! previous behaviour (flux conservation) is used. If a missing value is
  !! specified then pixels which match this value do not contribute to the final
  !! value of any output pixel they intersect. Except, and this is the rub, where
  !! there are *nothing but missing value pixels contributing to the output*. In
  !! this case the output pixel is given the missing value, irrespective of size.
  !! This is important -- missing pixels are simply propagated with the same value
  !! and are not altered. We don't want -9999 becoming -123.678 do we!? The value
  !! of a missing pixel is a flag, nothing more. It is not real data and should
  !! not be confused with real data.
  !!
  !! Made the debug value a parameter. This means it is "optimized out" by the
  !! compiler, so we can keep all the debug flags in there for the time being
  !! until we are sure everything is working ok and not incur and performance
  !! penalties.
  !!
  !! Fixed a bug in polySamp (the result of samplePixel was assigned to the wrong
  !! index -- needed nx*(j-1) rather than nx*j).
  !!
  !! Revision 1.1  2004/07/23 05:08:51  aidan
  !! Initial revision
  !!

  ! This image to be resampled
  real(kind=rd_kind), allocatable, dimension(:,:) :: original
  
  ! X and Y dimensions of image
  integer :: mx, my

  ! Is this an intensive image
  logical :: isintensive = .false. 
    
  ! Do we have missing pixels, i.e. with "no-data"?
  logical :: missing = .false.

  ! What is the value we use to represent missing data
  ! and what is the minimum coverage (proportion) of a pixel of
  ! *not-missing* data before we will accept this as ok?
  real(kind=rd_kind) :: missing_value
  real(kind=rs_kind) :: min_coverage
    
  ! Drizzle offset and area
  real(kind=rd_kind) :: drizzOffset, drizzArea
    
  ! Intermediate storage used by rectClip.
  ! The maximum number of vertices we will get if we start with
  ! a convex quadrilateral is 12, but we use larger
  ! arrays in case this routine is used is some other context.
  ! If we were good we'd be checking this when we add points in
  ! the clipping process.
  real(kind=rd_kind), dimension(100) :: rcX0, rcX1, rcY0, rcY1

  ! Intermediate storage used by polySamp 
  real(kind=rd_kind), dimension(12) :: psX0, psY0, psX1, psY1
    
  ! logical, parameter :: debug =  .TRUE.
  logical, parameter :: debug =  .FALSE.

  ! Public routines
  public :: polysamp

contains
    
  ! Calculate the area of a convex polygon.
  ! This function calculates the area of a convex polygon
  ! by deconvolving the polygon into triangles and summing
  ! the areas of the consituents.  The user provides the
  ! coordinates of the vertices of the polygon in sequence
  ! along the circumference (in either direction and starting
  ! at any point).
  ! 
  ! Only distinct vertices should be given, i.e., the
  ! first vertex should not be repeated at the end of the list.

  real(kind=rd_kind) function convexArea( n, x, y )

    ! Interface variables

    ! Number of vertices
    integer, intent(in) :: n
    ! The x and y coordinates of the vertices
    real(kind=rd_kind), dimension(:), intent(in) :: x, y

    ! Internal variables
    integer :: i

    convexArea = 0_rd_kind
    do i = 2, n-1
       convexArea = convexArea + triangleArea((/x(1),x(i:i+1)/),(/y(1), y(i:i+1)/))
    end do

  end function convexArea
    
  ! Calculate the area of an arbitrary triangle.
  !  Use the vector formula
  !     A = 1/2 sqrt(X^2 Y^2 - (X-Y)^2)
  !  where X and Y are vectors describing two sides
  !  of the triangle.
  ! 
  ! 
  !  @return         Area of the triangle.
    
  real(kind=rd_kind) function triangleArea(x, y)

    ! Interface variables
    ! x,y coordinates of the first, second and third vertices
    real(kind=rd_kind), dimension(3), intent(in) :: x, y

    ! Local variables
    real(kind=rd_kind) :: a, b, e, f

    a = x(1)-x(2)
    b = y(1)-y(2)
    e = x(1)-x(3)
    f = y(1)-y(3)
    
    if (a == 0 .and. b == 0 .and. e == 0 .and. f == 0) then
       triangleArea = 0.0_rd_kind  ! Assume this is due to roundoff errors
    else
       triangleArea = (a*a+b*b)*(e*e+f*f) - (a*e+b*f)*(a*e+b*f)

       if (triangleArea <= 0.0_rd_kind) then
          triangleArea = 0.0_rd_kind  ! Assume this is due to roundoff errors
       else
          triangleArea = sqrt(triangleArea)/2.d0
       end if
    end if

  end function triangleArea

  ! Is the test value on the proper side of a line.
  logical function inPlane(test, divider, direction)

    ! Interface variables
    ! 
    ! test - Value to be tested
    ! divider - Critical value
    ! direction True if values greater than divider are 'in'
    !           False if smaller values are 'in'.
    real(kind=rd_kind), intent(in) :: test, divider
    logical, intent(in)      :: direction

    ! Note that since we always include
    ! points on the dividing line as 'in'.  Not sure
    ! if this is important though...
    
    if (direction) then
       inPlane = (test >= divider)
    else
       inPlane = (test <= divider)
    end if

  end function inPlane
    
  ! Clip a polygon to a half-plane bounded by a vertical line.
  ! Users can flip the input axis order to clip by a horizontal line.
  ! 
  ! This function uses pre-allocated arrays for
  ! output so that no new objects are generated
  ! during a call.
  !  
  !  @param n Number of vertices in the polygon
  !  @param x X coordinates of the vertices
  !  @param y Y coordinates of the vertices
  !  @param nx New X coordinates
  !  @param ny      New Y coordinates
  !  @param val Value at which clipping is to occur.
  !  @param dir     Direction for which data is to be
  !                  clipped.  true-> clip below val, false->clip above val.
  ! 
  !  @return         The number of new vertices.
  ! 
  !    
  integer function lineClip ( n, x, y, nx, ny, value, direction ) result(nout)

    integer, intent(in) :: n

    real(kind=rd_kind), dimension(:), intent(in)  :: x, y
    real(kind=rd_kind), dimension(:), intent(out) :: nx, ny

    real(kind=rd_kind), intent(in) :: value
    logical, intent(in) :: direction

    ! Local variables
    real(kind=rd_kind) :: ycross
    integer :: i
    logical :: last

    nout = 0

    ! Need to handle first segment specially
    ! since we don't want to duplicate vertices.
    last = inPlane(x(n), value, direction)
    ! print *,'last? ',last,' x(n): ',x(n),' value: ',value,' direction: ',direction

    if (debug) print '(A,I4,L,100(F8.2))','LINECLIP:   ',n,last,x(:n),y(:n)

    do i = 1, n
       
       if (last) then

          if (inPlane(x(i), value, direction)) then
             ! Both endpoints in, just add the new point
             nout = nout + 1
             nx(nout) = x(i)
             ny(nout) = y(i)
          else
             ! Moved out of the clip region, add the point we moved out
             if (i == 1) then
                ycross = y(n) + (y(1)-y(n))*(value-x(n))/(x(1)-x(n))
             else
                ycross = y(i-1) + (y(i)-y(i-1))*(value-x(i-1))/(x(i)-x(i-1))
             end if
             nout     = nout + 1
             nx(nout) = value
             ny(nout) = ycross
             last     = .false.
          end if
          
       else
          
          if (inPlane(x(i), value, direction)) then
             ! Moved into the clip region.  Add the point
             ! we moved in, and the end point.
             if (i == 1) then
                ycross = y(n) + (y(1)-y(n))*(value-x(n))/(x(1)-x(n))
             else
                ycross = y(i-1) + (y(i)-y(i-1))*(value-x(i-1))/(x(i)-x(i-1))
             end if
             nout = nout + 1
             nx(nout) = value
             ny(nout) = ycross

             nout = nout + 1
             nx(nout) = x(i)
             ny(nout) = y(i)
             last     = .true.

          else
             ! Segment entirely clipped.
          end if
       end if
    end do
    if (debug) print '(A,I5,100(F8.2))','NEWLINECLIP:',nout,nx(:nout),ny(:nout)

  end function lineClip
    
    
  ! Clip a polygon by a non-rotated rectangle.
  ! 
  !  This uses a simplified version of the Sutherland-Hodgeman polygon
  !  clipping method.  We assume that the region to be clipped is
  !  convex.  This implies that we will not need to worry about
  !  the clipping breaking the input region into multiple
  !  disconnected areas.
  !    [Proof: Suppose the resulting region is not convex.  Then
  !     there is a line between two points in the region that
  !     crosses the boundary of the clipped region.  However the
  !     clipped boundaries are all lines from one of the two
  !     figures we are intersecting.  This would imply that
  !     this line crosses one of the boundaries in the original
  !     image.  Hence either the original polygon or the clipping
  !     region would need to be non-convex.]
  ! 
  !  Private arrays are used for intermediate results to minimize
  !  allocation costs.
  ! 
  !  @param n Number of vertices in the polygon.
  !  @param x X values of vertices
  !  @param y        Y values of vertices
  !  @param nx X values of clipped polygon
  !  @param ny       Y values of clipped polygon
  ! 
  !  @param          minX Minimum X-value
  !  @param  minY Minimum Y-value
  !  @param          maxX MAximum X-value
  !  @param          maxY Maximum Y-value
  ! 
  !  @return  Number of vertices in clipped polygon.
  !/
  integer function rectClip(n, x, y, nx, ny, minX, minY, maxX, maxY) result(nCurr)

    ! Interface variables
    integer, intent(in) :: n
    real(kind=rd_kind), dimension(:), intent(in)  :: x, y
    real(kind=rd_kind), dimension(:), intent(out) :: nx, ny

    real(kind=rd_kind), intent(in) :: minX, minY, maxX, maxY

    ! debug = .FALSE.
    ! if (.NOT. any( (/minX, minY, maxX, maxY/) /= (/ 0, 2, 1, 3 /))) debug = .TRUE.

    ! lineClip is called four times, once for each constraint.
    ! Note the inversion of order of the arguments when
    ! clipping vertically.

    nCurr = lineClip(n, x, y, rcX0, rcY0, minX, .true.)

    if (nCurr > 0) then
       nCurr = lineClip(nCurr, rcX0, rcY0, rcX1, rcY1, maxX, .false.)
       if (nCurr > 0) then
          nCurr = lineClip(nCurr, rcY1, rcX1, rcY0, rcX0, minY, .true.)
          if (nCurr > 0) then
             nCurr = lineClip(nCurr, rcY0, rcX0, ny, nx, maxY, .false.)
          end if
       end if
    end if

    ! We don't need to worry that we might not have set the output arrays.
    ! If nCurr == 0, then it doesn't matter that
    ! we haven't set nx and ny.  And if it is then we've gone
    ! all the way and they are already set.

  end function rectClip
    
  ! Sample a single map pixel.
  ! 
  !  @param x The x values of the corners of the pixel [4]
  !  @param y The y values of the corners of the pixel [4]
  ! 
  !  @return  The total flux in the resampled pixel.
  !
  real(kind=rd_kind) function samplePixel(x, y) result(value)

    ! public double samplePixel(double[] x, double[] y) {

    ! Interface variables
    real(kind=rd_kind), dimension(4), intent(in) :: x, y

    ! Local variables
    real(kind=rd_kind) :: minX, maxX, minY, maxY, vminX, vmaxX, vminY, vmaxY
    real(kind=rd_kind) :: tArea, area, factor, missingArea
    integer :: k, n, m, nv

    value = 0.d0

    ! Find a bounding box for the pixel coordinates.
    minX = minval(x)
    maxX = maxval(x)
    minY = minval(y)
    maxY = maxval(y)

    ! Round the extrema of the pixel coordinates to
    ! integer values.
    minX = floor(minX)
    maxX = ceiling(maxX)

    minY = floor(minY)
    maxY = ceiling(maxY)

    if (debug) print *,'Min/Max x and y: ',minX,maxX,minY,maxY

    ! Check to see if pixel is entirely off original image.
    ! If so we don't need to do anything further.
    if ((maxX <= 0.d0) .OR. (minX >= mx) .OR. (maxY <= 0.d0) .OR. (minY >= my)) then
       if (debug) print *,'mx and my: ',mx, my
       if (debug) print *,'Off image: ',x,y
       if (debug) print *,'Set to misssing? ',missing,missing_value
       if (missing) value = missing_value
       return
    end if

    ! Check if the resampling pixel is entirely enclosed in
    ! the image pixel.  Need to check this before
    ! we 'clip' our bounding box.  This check
    ! should significantly increase the speed
    ! for oversampling images, but it will cause
    ! a tiny slowdown when the sampling pixels are as large
    ! or larger than the original pixels.
    ! 
    ! We're doing equalities with
    ! double values, but they are guaranteed to be
    ! integers.
    if (minX == (maxX-1) .AND. minY == (maxY-1)) then
       if (debug) print *,'Totally bounded pixel man at ',int(maxX),int(maxY)
       if (debug) print *,'Original: ',original(int(maxX),int(maxY))
       if (debug) print *,'Scale: ',convexArea(4,x,y)
       value = original(int(maxX),int(maxY))*convexArea(4,x,y)
       return
    end if

    ! Clip the bounding box to the original image dimensions
    if (maxX >= mx) maxX = mx
    if (minX < 0) minX = 0

    if (maxY >= my) maxY = my
    if (minY < 0) minY = 0

    if (debug) print *,'Clipped Min/Max x and y: ',minX,maxX,minY,maxY

    ! Loop over the potentially overlapping pixels.
    tArea = 0.d0
    missingArea = 0.d0

    ! Loop over the *start* values of all the pixels that are 
    ! potentially inside our sample pixel. Below we make the
    ! range n to n + 1, hence the loop is only over minY to maxY - 1
    do n = int(minY), int(maxY)-1

       ! the vmin/max values are the areas in
       ! which the 'flux' of a given pixel is treated
       ! as being.

       vminY = n + drizzOffset
       vmaxY = n + 1 - drizzOffset

       do m = int(minX), int(maxX)-1

          vminX = m + drizzOffset
          vmaxX = m + 1 - drizzOffset

          ! Clip the quadrilaterel given by the coordinates
          ! of the resampling pixel, to this particular
          ! image pixel.
          nv = rectClip(4, x, y, psX1, psY1, vminX, vminY, vmaxX,vmaxY)
          ! print '(3(4F5.2),)', x, y, vminX, vmaxX, vminY, vmaxY, nv
          ! if (debug) print '(A,F0.3,A,F0.3,F0.3,A,3F0.3,I5)', 'rectclip', vminX,' -> ', vmaxX, vminY,' -> ', vmaxY, x, y, nv
          if (debug) print *, 'rectclip', vminX,' -> ', vmaxX, vminY,' -> ', vmaxY, x, y, nv

          ! If there is no overlap we won't get any
          ! vertices back in the clipped set.
          if (nv > 0) then
             ! Calculate the area of the clipped pixel and compare
             ! it to the area in which the flux of the original
             ! pixel is found.  The returned area should
             ! never be greater than the drizzArea.
             area = convexArea(nv,psX1,psY1)
             factor = area/drizzArea
             
             ! Check if we are processing missing pixels
             if (missing) then
                ! If our original pixel has the "missing" value then we'll
                ! calculate the area of this pixel and then go to the next 
                ! pixel (note that this is an equality between two real numbers, 
                ! so the pixel must match the missing value *exactly*
                if (original(m+1,n+1) == missing_value) then
                   missingArea = missingArea + area
                   cycle
                end if
             end if
             tArea  = tArea + area

             ! Add the appropriate fraction of the original
             ! flux into the output pixel.
             value = value + factor*original(m+1,n+1)
             if (debug) print *,'nv area ...',nv,area,psX1(1:nv),psY1(1:nv),m+1,n+1,original(m+1,n+1),factor,value
          end if
       end do
    end do
    if (missing) then
       if (tArea == 0.d0) then
          ! If we had no signal in our new pixel, then propagate our
          ! missing_value into it
          value = missing_value
       else if (missingArea > 0.d0) then
          if (tArea/(tArea + missingArea) >= min_coverage) then
             ! Scale value up to account for areas of missing data. This
             ! assumes that data would have been roughly homogenous. It 
             ! is a greater sin to *not* do this -- it will make some
             ! areas far too weak otherwise (holes in peaks etc)
             value = ((tArea + missingArea)/tArea)*value
          else
             value = missing_value
          end if
       end if
    end if

  end function samplePixel
    
  ! Perform a flux conserving resampling of one image to a
  ! resampled map.  In this version we are given the
  ! the pixel corners as matrices of X and Y arrays.
  !
  ! Resampling is done by clipping each output pixel to the box
  ! specified by one of the original pixels.  The area of the
  ! overlap between the pixels is computed and an appropriate fraction
  ! of the flux from the original pixel is added to the output pixel.
  !
  ! To minimize object allocation costs private arrays are
  ! used for intermediate storage.
  !
  ! A number of the arrays used here are linearized versions
  ! of two-d arrays.  The 'real' dimensions are indicated
  ! in the comments.
  !
  !  @param nx       X dimension of resampled map
  !  @param ny       Y dimension of resampled map
  !  @param x        X coordinates of corners of pixels in resampled map
  !                    [(nx+1)*(ny+1)]. These coordinates are in terms
  !                    of the pixels in the original image.
  !  @param y        Y coordinates of corners of pixels in resampled map
  !                    [(nx+1)*(ny+1)]. These coordinates are in terms
  !                    of the pixels in the original image.
  !  @return         Array of resampled pixels [nx x ny].
  
  
  subroutine polySamp (input, x, y, output, missingvalue, intensive, coverage)  

    ! subroutine polySamp (input, nx, ny, x, y, output)  
    ! function polySamp (input, nx, ny, x, y) result(output)  

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: input
    real(kind=rd_kind), dimension(:), intent(in)    :: x, y
    real(kind=rd_kind), dimension(:,:), intent(out) :: output
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    logical, intent(in), optional                   :: intensive
    real(kind=rs_kind), intent(in), optional        :: coverage
  
    ! Local variables
    integer :: nx, ny
    real(kind=rd_kind), dimension(size(output)) :: output_1d
    real(kind=rd_kind) :: area
    real(kind=rd_kind), dimension(0:3) :: px, py
    integer :: i, j, p0, p1, p2, p3, dest

    nx = size(output,1)
    ny = size(output,2)
    
    if (debug) print '(A,I0,A,I0)','Output is ',nx,' x ',ny

    call setDrizzle(1.d0)
    if (present(intensive)) then
       call setIntensive(intensive)
    else
       call setIntensive(.FALSE.)
    end if

    if (present(missingvalue)) then
       if (present(coverage)) then
          call setMissing(missingvalue,coverage)
       else
          call setMissing(missingvalue)
       end if
    else
       if (present(coverage)) then
          write(stderr,*) 'POLYSAMPLE :: Specifying coverage without specifying'&
               &' a missing value is meaningless. Ignoring this option'
       end if
       call setMissing()
    end if

    allocate(original(size(input,1),size(input,2)))

    original = input
    mx = size(original,1)
    my = size(original,2)

    ! Loop over the output pixels.
    do j = 1, ny
       if (debug) write(*,'(A,I4)') 'j: ',j
       do i=1, nx
          if (debug) write(*,'(A,I4)',advance='NO') ' i: ',i

          ! The p's need to go around the vertices of
          ! the pixel sequentially as if we were traversing
          ! the pixel circumference.  Either direction is fine.
          ! 
          ! Remember that the X and Y arrays are one larger
          ! than the number of pixels in each direction.
          
          p0 = i  + (j-1)*(nx+1)
          p1 = p0 + 1
          p2 = p1 + (nx+1)
          p3 = p2 - 1
          if (debug) print *,p0,p1,p2,p3
          
          ! Intialize the arrays of pixel vertex coordinates
          px(0) = x(p0)
          px(1) = x(p1)
          px(2) = x(p2)
          px(3) = x(p3)
          
          py(0) = y(p0)
          py(1) = y(p1)
          py(2) = y(p2)
          py(3) = y(p3)
          
          dest = i + nx*(j-1)
          
          if (debug) print *,px
          if (debug) print *,py
          
          if (debug) print *,'index: ',(i+nx*(j-1))
          if (dest > size(output_1d)) print *,'PANIC!',i,j,(i+nx*(j-1))
          output_1d(dest) = samplePixel(px, py)
          if (debug) print *,'output: ',output_1d(dest)
          
          ! For intensive quantities convert from sum to weighted average.
          ! We just need to divide by the total area of the output pixel.
          if (isintensive) then
             if (.NOT.(missing .and. (output_1d(dest) == missing_value))) then
                area = convexArea(4, px, py)
                if (dest > size(output_1d)) then
                   print *,'PANIC!',i,j,p0,dest
                   print *,p0,p1,p2,p3
                end if
                if (area > 0) then
                   output_1d(dest) = output_1d(dest) / area
                else
                   if (missing) then
                      output_1d(dest) = missing_value
                   else
                      output_1d(dest) = 0
                   end if
                end if
             end if
          end if
          
          ! write(*,'(I4)',advance='NO') i
          
       end do ! Y Ends of loops over output pixels.
    end do ! X
    
    if (debug) print *,shape(output),size(output),size(output_1d)
    
    ! We should have filled in each map pixel.
    if (debug) print *,'max: ',maxval(output_1d)
    output = 0.d0
    output = reshape(output_1d,(/nx,ny/))
    if (debug) print *,'Out of loop'
    if (debug) print *,'max: ',maxval(output)
    
    deallocate(original)
    
    return
    
  end subroutine polySamp
  ! end function polySamp

  ! Set the drizzle factor for sampling
  !  @param drizzle The drizzle factor should range from >0 to 1 and
  !                 indicates the length of the side of the pixel
  !                 inside the original image pixel in which the
  !                 flux is assumed to be contained.
  subroutine setDrizzle(drizzle)
    
    real(kind=rd_kind) :: drizzle
    ! public void setDrizzle(double drizzle) {
    
    if (drizzle <= 0) drizzle = .01
    if (drizzle > 1) drizzle = 1
    
    drizzArea = drizzle*drizzle
    drizzOffset = (1.d0-drizzle)/2.d0
    
  end subroutine setDrizzle
  
  ! Indicate whether the image is intensive.
  !  If intensive is set to true, then the drizzle is set to 1.
  !  @param intensive	Indicates whether the image
  !                          has a dimensionality like temperature
  !                          independent of the size of the pixel.
  subroutine setIntensive(intensive)
    ! public void setIntensive(boolean intensive) {
    
    logical, intent(in) :: intensive
    
    isintensive = intensive
    
    if (isintensive) call setDrizzle(1.d0)
    
  end subroutine setIntensive
  
  subroutine setMissing(value, coverage)
    
    ! Used to set a missing value, and minimum coverage
    
    real(kind=rd_kind), intent(in), optional :: value
    real(kind=rs_kind), intent(in), optional :: coverage
        
    missing = .FALSE.
    min_coverage = 0.0
    
    ! Set the missing value if one is specified
    if (present(value)) then
       missing = .TRUE.
       missing_value = value
    end if
    
    ! Set the minimum coverage value if one is specified
    if (present(coverage)) then
       min_coverage = coverage
    end if
    
  end subroutine setMissing
  
end module polysample
