c     Draws a rectangular box in window 2 (it is assumed that another
c     program has created that window).  The coordinates are based on
c     the window coordinates: (0,0) is the lower lefthand corner
c     and the units are in pixels.

      program drawbox
        
      real*4 plist(10)

      if (imopen(1, "2", "old") .ne. 0) then
         print *, 'Could not open window 2.'
         stop
      endif

      plist(1) = 20
      plist(2) = 20
      plist(3) = 20
      plist(4) = 100
      plist(5) = 100
      plist(6) = 100
      plist(7) = 100
      plist(8) = 20
      plist(9) = 20
      plist(10) = 20

      i1 = IWGrAddLns2D(1, plist, 5, 1, 0, 1)
      do j = 1,10
         plist(j) = plist(j) + 10
      enddo
      i2 = IWGrAddLns2D(1, plist, 5, 1, 0, 1)

      plist(1) = 20
      plist(2) = 20
      plist(3) = 300
      plist(4) = 300
      plist(5) = 20
      plist(6) = 300
      plist(7) = 300
      plist(8) = 20

      i3 = IWGrAddMultLns2D(1, plist, 4, 3, 0, 2)

      plist(1) = 0
      plist(2) = 0
      plist(3) = 100
      plist(4) = 0

      i4 = IWGrAddLns2D(1, plist, 2, 1, 0, 1)
      call IWDisplay(1)
      call sleep(2)
      call IWGrRmGrID(1, i1)
      call IWAlClrBkg(1, 1)
      call IWDisplay(1)
      call sleep(2)
      call IWGrRmGrID(1, i2)
      call IWAlClrBkg(1, 1)
      call IWDisplay(1)
      call sleep(2)
      call IWGrRmGrID(1, i3)
      call IWAlClrBkg(1, 1)
      call IWDisplay(1)
      call sleep(2)
      call IWGrRmGrID(1, i4)
      call IWAlClrBkg(1, 1)
      call IWDisplay(1)

      end
