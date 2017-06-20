program test1

  use polysample
!!$  use pnm_class

  implicit none

!!$  type (pnm_object) :: pnm

  integer :: width, height, i
  integer, parameter :: newxsize = 401, newysize = 401
  integer, dimension(:,:), allocatable :: data
  real(kind=8), dimension(:,:), allocatable :: realdata
  real(kind=8), dimension(newxsize,newysize) :: output
  integer, dimension(newxsize,newysize) :: intoutput

  real(kind=8), dimension((newxsize+1)*(newysize+1)) :: xvertices, yvertices
  real(kind=8), dimension((newxsize+1),(newysize+1)) :: tmp

  real(kind=8) :: xinc, yinc
  real(kind=8) :: smalltest(5,5), smallerresult(1,1), smallresult(2,2)

  real :: tolerance = 0.0001

  ! Let us test this thing!

  call announce("Testing simple re-gridding")

  ! We start with a 5x5 grid filled with 5, which makes a total of 125
  smalltest = 5.d0

  ! Now we resample this to a 1x1 grid which is rotated by 45 deg, like
  ! so:
  !
  !      ----------*----------
  !      |   |   *   *   |   |
  !      ------*-------*------
  !      |   *   |   |   *   |
  !      --*---------------*--
  !      *   |   |   |   |   *
  !      --*---------------*--
  !      |   *   |   |   *   |
  !      ------*-------*------
  !      |   |   *   *   |   |
  !      ----------*----------
  !
  ! Trigonometry tells us that this should result in a single cell whose
  ! "flux" is half the total contained in the 5x5 grid.

  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult)

  if (abs(sum(smallerresult) - 62.5) > 0.00001) then
     print *,'This should be 125/2 = 62.5:     ',smallerresult
     call failed()
     stop
  end if

  ! Now we set the central pixel to zero ...

  smalltest(3,3) = 0.d0

  ! now the grid looks like so:
  !
  !      ----------*----------
  !      |   |   *   *   |   |
  !      ------*-------*------
  !      |   *   |   |   *   |
  !      --*---------------*--
  !      *   |   |-0-|   |   *
  !      --*---------------*--
  !      |   *   |   |   *   |
  !      ------*-------*------
  !      |   |   *   *   |   |
  !      ----------*----------
  !

  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult)

  ! And the resulting flux should be 5 less than it was:
  if (sum(smallerresult) /= 57.5) then
     print *,'This should be 125/2 - 5 = 57.5: ',smallerresult
     call failed
     stop
  end if

  ! Now try a different sized output grid like so:
  !
  !      ----------*----------
  !      |   |   *   *   |   |
  !      ------*-------*------
  !      |   * * |   | * *   |
  !      --*-----*---*-----*--
  !      *   |   |-0-|   |   *
  !      --*-----*---*-----*--
  !      |   * * |   | * *   |
  !      ------*-------*------
  !      |   |   *   *   |   |
  !      ----------*----------
  !

  call polySamp(smalltest, &
       real((/0.0,1.25,2.5,1.25,2.5,3.75,2.5,3.75,5.0/),8), &
       real((/2.5,1.25,0.0,3.75,2.5,1.25,5.0,3.75,2.5/),8), &
       smallresult)

  ! And the resulting flux should be 57.5
  if (abs(sum(smallerresult) - 57.5) > tolerance) then
     print *,'This should be 125/2 - 5 = 57.5: ',sum(smallresult)
     call failed
     stop
  end if

  call ok()

  call announce("Testing dealing with missing data")

  ! With the grid as above, we now specify a "missing data" value
  ! of zero -- we should get the same as the original answer, as
  ! the algorithm effectively replaces the missing pixel with the 
  ! mean value of the rest of the pixels.

  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult,missingvalue=0.d0)

  if (abs(sum(smallerresult) - 62.5) > tolerance) then
     print *,'This should be 125/2 = 62.5:     ',smallerresult
     call failed()
     stop
  end if

  ! Now try a different sized output grid
  call polySamp(smalltest, &
       real((/0.0,1.25,2.5,1.25,2.5,3.75,2.5,3.75,5.0/),8), &
       real((/2.5,1.25,0.0,3.75,2.5,1.25,5.0,3.75,2.5/),8), &
       smallresult,missingvalue=0.d0)

  ! And the resulting flux should be 57.5
  if (abs(sum(smallerresult) - 62.5) > tolerance) then
     print *,'This should be 125/2 = 62.5: ',sum(smallresult)
     call failed
     stop
  end if

  ! Re-initialise the grid to all '5's and make a different pixel
  ! zero
  smalltest = 5.d0
  smalltest(2,2) = 0.d0

  ! now the grid looks like so:
  !
  !      ----------*----------
  !      |   |   *   *   |   |
  !      ------*-------*------
  !      |   *   |   |   *   |
  !      --*---------------*--
  !      *   |   |   |   |   *
  !      --*---------------*--
  !      |   *-0-|   |   *   |
  !      ------*-------*------
  !      |   |   *   *   |   |
  !      ----------*----------
  !

  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult)

  ! And the resulting flux should be 4.625 (5 - 1/8*5) less than 62.5
  if (sum(smallerresult) /= 58.125) then
     print *,'This should be 125/2 - 2.5 = 58.125: ',smallerresult
     call failed
     stop
  end if

  ! Now 'filling-in' missing values
  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult,missingvalue=0.d0)

  ! And the resulting flux should be 4.625 (5 - 1/8*5) less than 62.5
  if (abs(sum(smallerresult) - 62.5) > tolerance) then
     print *,'This should be 125/2 = 62.5: ',smallerresult
     call failed
     stop
  end if

  ! And again ..
  smalltest = 5.d0
  smalltest(1,2) = 0.d0

  ! now the grid looks like so:
  !
  !      ----------*----------
  !      |   |   *   *   |   |
  !      ------*-------*------
  !      |   *   |   |   *   |
  !      --*---------------*--
  !      *   |   |   |   |   *
  !      --*---------------*--
  !      |-0-*   |   |   *   |
  !      ------*-------*------
  !      |   |   *   *   |   |
  !      ----------*----------
  !

  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult)

  ! And the resulting flux should be .625 (1/8*5) less than 62.5
  if (abs(sum(smallerresult) - 61.875) > tolerance) then
     print *,'This should be 125/2 - 2.5 = 61.875: ',smallerresult
     call failed
     stop
  end if

  ! Now 'filling-in' missing values
  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult,missingvalue=0.d0)

  ! And the resulting flux should be 4.625 (5 - 1/8*5) less than 62.5
  if (abs(sum(smallerresult) - 62.5) > tolerance) then
     print *,'This should be 125/2 = 62.5: ',smallerresult
     call failed
     stop
  end if

  call ok()

  call announce('Testing minimum coverage')

  ! Now specifying a minimum coverage of 50%
  call polySamp(smalltest, &
       real((/0.5,1.5,0.5,1.5/),8), &
       real((/1.5,1.5,2.0,2.0/),8), &
       smallerresult,missingvalue=0.d0,coverage=0.50)

  ! And the resulting flux should be 2.5
  if (sum(smallerresult) /= 2.5) then
     print *,'This should be 2.5: ',sum(smallerresult)
     call failed
     stop
  end if

  ! Same minimum coverage of 50% but intensive
  call polySamp(smalltest, &
       real((/0.5,1.5,0.5,1.5/),8), &
       real((/1.5,1.5,2.0,2.0/),8), &
       smallerresult,missingvalue=0.d0,intensive=.TRUE.,coverage=0.5)

  ! And the resulting flux should be 5.0
  if (sum(smallerresult) /= 5.0) then
     print *,'This should be 5.0: ',sum(smallerresult)
     call failed
     stop
  end if

  ! Now the minimum coverage is 51%
  call polySamp(smalltest, &
       real((/0.5,1.5,0.5,1.5/),8), &
       real((/1.5,1.5,2.0,2.0/),8), &
       smallerresult,missingvalue=0.d0,coverage=0.51)

  ! And the resulting flux should be 0.0
  if (sum(smallerresult) /= 0.0) then
     print *,'This should be 0.0: ',sum(smallerresult)
     call failed
     stop
  end if

  call ok()

  call announce('Testing different missing values')

  ! Now we'll experiment with different missing values
  smalltest = 5.d0
  smalltest(1,2) = 1.d0

  ! First make 1.d0 the missing value ...
  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult,missingvalue=1.d0)

  ! And the resulting flux should be 4.625 (5 - 1/8*5) less than 62.5
  if (abs(sum(smallerresult) - 62.5) > tolerance) then
     print *,'This should be 125/2 = 62.5: ',smallerresult
     call failed
     stop
  end if

  ! Now make 5.d0 the missing value ...
  call polySamp(smalltest, &
       real((/0.0,2.5,2.5,5.0/),8), &
       real((/2.5,0.0,5.0,2.5/),8), &
       smallerresult,missingvalue=5.d0)

  ! And the resulting flux should be 12.5. Why? Well we are ignoring 
  ! pixels with the value 5.d0, i.e. 24/25 pixels, and the only pixel
  ! we are counting has a value of 1.d0. Then we, effectively, fill
  ! in the rest of the area with the mean value, 1.d0, and get a total
  ! which is the total area of the new pixel.
  if (sum(smallerresult) /= 12.5) then
     print *,'This should be 25/2 = 12.5: ',smallerresult
     call failed
     stop
  end if

  call ok()

  call announce('Testing different sized grids')

  ! Now try a different sized output grid
  call polySamp(smalltest, &
       real((/0.0,1.25,2.5,1.25,2.5,3.75,2.5,3.75,5.0/),8), &
       real((/2.5,1.25,0.0,3.75,2.5,1.25,5.0,3.75,2.5/),8), &
       smallresult,missingvalue=5.d0)

  ! And the resulting flux should be 18.125. WHA!? 
  !
  ! Here goes ... we are now dividing our output into four grids, each
  ! of which is one quarter the size of the single pixel grid, which is
  ! half the size of the original. So each grid is (125/2)/4 = 3.125. 
  ! Three of our grids have "no data". How do we propagate no data? We
  ! bung the no data value into them, irrespective of their relative 
  ! size. The fourth? Well it has a single pixel of value 1.d0, and this
  ! is multiplied by it's relative size (as described above) and so 
  ! becomes 3.125.
  if (abs(sum(smallresult) - 18.125) > tolerance) then
     print *,'This should be 3 * 5 + (25/2)/4 = 18.125: ',sum(smallresult)
     call failed
     stop
  end if

  ! Just to prove this does as it claims, if we go back to the original
  ! example and make the central pixel the missing one ...
  smalltest = 5.d0
  smalltest(3,3) = 1.d0

  call polySamp(smalltest, &
       real((/0.0,1.25,2.5,1.25,2.5,3.75,2.5,3.75,5.0/),8), &
       real((/2.5,1.25,0.0,3.75,2.5,1.25,5.0,3.75,2.5/),8), &
       smallresult,missingvalue=5.d0)

  ! And the resulting flux should be 12.5. As there are no pixels with
  ! exclusively missing pixels, so all four pixels should have some
  ! portion of their pixels with value 1.d0, so 
  if (abs(sum(smallerresult) - 12.5) > tolerance) then
     print *,'This should be 4 * (25/2)/4 = 12.5: ',sum(smallresult)
     call failed
     stop
  end if

  call ok()

  stop

!!$  call announce("Testing scaling of an image")

!!$  call read(pnm,'test.pgm')
!!$
!!$  width = size(pnm,1)
!!$  height = size(pnm,2) 
!!$
!!$  allocate(data(width,height))
!!$  allocate(realdata(width,height))
!!$
!!$  data = pnm
!!$
!!$  realdata = real(data,8)
!!$
!!$  xinc = real(width,8) / real(newxsize,8)
!!$  yinc = real(height,8) / real(newysize,8)
!!$
!!$  ! print *,shape(spread( real((/ (i, i = 0, newxsize) /),8)*xinc, 2, newysize+1))
!!$  ! print *,size(pack(spread( real((/ (i, i = 0, newxsize) /),8)*xinc, 2, newysize+1),.TRUE.))
!!$  ! print *,size(xvertices)
!!$  ! tmp = spread( real((/ (i, i = 0, newxsize) /),8)*xinc, 2, newysize+1)
!!$  xvertices = pack(tmp,.TRUE.)
!!$
!!$  xvertices = pack(spread( real((/ (i, i = 0, newxsize) /),8)*xinc, 2, newysize+1),.TRUE.)
!!$  yvertices = pack(spread( real((/ (i, i = 0, newysize) /),8)*yinc, 1, newxsize+1),.TRUE.)
!!$
!!$  call polySamp(realdata, xvertices, yvertices, output, missingvalue = 0.d0)
!!$
!!$  ! print *,maxval(output)
!!$
!!$  output = output /(maxval(output)/65535)
!!$
!!$  intoutput = int(output)
!!$
!!$  pnm = intoutput
!!$  call write(pnm,'output.pgm')
!!$
!!$  call ok()

  ! Tests completed

  contains 

    subroutine announce(text)
      character(len=*) :: text
      write (*,'(A)',advance='no') text
   end subroutine announce
    
    subroutine ok
      print *,' ........... ok'
    end subroutine ok

    subroutine failed
      print *,' ........... FAILED!'
   end subroutine failed


end program test1
