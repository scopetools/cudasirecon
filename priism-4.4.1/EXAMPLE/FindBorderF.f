c     Given a dataset, returns a dataset where points which are above
c     a threshold but which have one 4-connected neighbor below the
c     threshold are set to 1 and all other points are set to zero.
c
c     Does not handle complex-valued data.
c
c     Identical to FindBorder.c but written in Fortran.


      subroutine IPAppSpecifics()
      include 'ip.inc'

      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      external NullOp
      external CheckParams, ProcData
      external SetMenu
      external ReadCmdLine, DisplayUsage, WriteCmdLine
      external PromptForParameters

      thresh = 0.0

      call IPAddInput(1, 'Open', 1, 0, NullOp, 0)
      call IPAddOutput(2, 'Save', 1, 0, NullOp, 0)
c
c     Since complex values are not handled, do not allow the user to
c     choose complex output.
c
      call IPAllowMode(2, 3, 0)
      call IPAllowMode(2, 4, 0)
      call IPSetCustRoutine(CHECK_PARAMS, CheckParams, 0)
      call IPSetCustRoutine(PROC_FUNC, ProcData, 0)
      call IPSetMenuTitle('Find Border')
      call IPSetMenuLoc(450, 400)
      call IPSetCustRoutine(CUS_MENU, SetMenu, 0)
      call IPSetCustRoutine(CMD_LINE_FUNC, ReadCmdLine, 0)
      call IPSetCustRoutine(USAGE_FUNC, DisplayUsage, 0)
      call IPSetCustRoutine(PROMPT_FUNC, PromptForParameters, 0)
      call IPSetCustRoutine(WRITE_CMD, WriteCmdLine, 0)
      call IPEnableResolutionSelection()

      return
      end


      subroutine CheckParams(idum)
      include 'ip.inc'
      integer mode

      if (IPIsOpen(1) .ne. 0) then
         call irtmod(1, mode)
         if (mode .eq. 3 .or. mode .eq. 4) then
            call IPDisplayMessage(
     &           'Can not handle complex-valued input', IP_LOG_ERROR
     &      )
            call IPDesist()
         endif
      endif

      return
      end


      subroutine ProcData(nxyz, dmmm, idum)
      include 'ip.inc'

      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      integer nxyz(3), idum
      real dmmm(3)
      real src(1), des(1)
      integer i, ix, j, iy, imy, ipy, npix
      pointer (psrc, src)
      pointer (pdes, des)

      psrc = IPGetDataPtr(1)
      pdes = IPGetDataPtr(2)
      
      imy = 0
      iy = nxyz(1)
      ipy = nxyz(1) * 2
      npix = 0
      do i = 1, nxyz(1)
         des(i) = 0.0
      enddo
      do j = 2, nxyz(2) - 1
         des(1 + iy) = 0.0
         do i = 2, nxyz(1) - 1
            ix = i + iy
            des(ix) = 0.0
            if (src(ix) .ge. thresh) then
               if (src(ix - 1) .lt. thresh .or.
     &              src(ix + 1) .lt. thresh .or.
     &              src(i + ipy) .lt. thresh .or.
     &              src(i + imy) .lt. thresh) then
                  des(ix) = 1.0
                  npix = npix + 1
               endif
            endif
         enddo
         iy = iy + nxyz(1)
         des(iy) = 0.0
         imy = imy + nxyz(1)
         ipy = ipy + nxyz(1)
      enddo
      do i = iy + 1, iy + nxyz(1)
         des(i) = 0.0
      enddo

      dmmm(1) = 0.0
      if (npix .eq. 0) then
         dmmm(2) = 0.0
      else
         dmmm(2) = 1.0
      endif
      dmmm(3) = dmmm(3) + npix

      return
      end


      subroutine PromptForParameters(idum)
      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      i = IPPromptReal('Threshold', 1, 1, 1, thresh)
      return
      end


      subroutine ReadCmdLine(idum)
      include 'ip.inc'
      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      integer iarg, istat
      character*80 arg
      character*100 msg

      iarg = 1
 890  continue
      istat = IPGetArg(iarg, arg)
      if (istat .eq. -1) then
         goto 900
      endif
      iarg = iarg + 1
      istat = IPParseRealArg(arg, 'threshold', 1, thresh)
      if (istat .eq. 1) then
         write(string, 690) thresh
 690     format(g)
      else if (istat .eq. 0) then
         write(msg, 700) arg
 700     format('Unrecognized option ', a)
         call IPDisplayMessage(msg, IP_LOG_ERROR)
         call IPDesist()
      else
         call IPDisplayMessage(
     &        '-threshold requires one real-valued argument',
     &        IP_LOG_ERROR
     &   )
         call IPDesist()
      endif
      goto 890
      
 900  continue

      return
      end


      subroutine DisplayUsage(idum)
      include 'ip.inc'

      call IPDisplayMessage(
     &     'Application-specific arguments:', IP_LOG_ERROR
     &)
      call IPDisplayMessage(
     &     '     [-threshold=value]', IP_LOG_ERROR
     &)
      return
      end


      subroutine WriteCmdLine(idum)
      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      character*40 arg

      arg = '-threshold='
      call LeftJustify(thresh, arg(12:40))
      call IPAppendCommand(arg)
      return
      end


      subroutine SetMenu(idum)
      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      character*24 key
      integer WMAddInfoButton, WMAddCharField, WMNewRow
      external GetString
      save key

      call LeftJustify(thresh, string)

c     Some compilers (notably Absoft's for Mac OS X) load a
c     string constant passed to a function into a temporary variable.
c     The lookup key passed to WMAddInfo has to persist as long
c     as the widget created does so use a static character array
c     to ensure that.  The label for the button is copied so it does
c     not matter if it is temporary or not.
      key = 'FindBorder Threshold'
      i = WMAddInfoButton('Threshold:', key)
      i = WMAddCharField(string, 5, 24, GetString, 0, 0, 4)
      call WMAttachRightSide()
      i = WMNewRow()

      return
      end


      integer function GetString(idum)
      character*24 string
      real thresh
      common /FindBorder/ string, thresh

      integer WMUpdateGroup

      GetString = 0
      read(string, 700, err=900, end=900) thresh
 700  format(g23.0)
      return

 900  continue
      call LeftJustify(thresh, string)
      i = WMUpdateGroup(4)

      return
      end
      

      subroutine LeftJustify(value, str)
      real value
      character*(*) str
      character*30 numstr
      integer offset

      write(numstr, 600) value
 600  format(g)
      offset = 0
 100  continue
      offset = offset + 1
      if (numstr(offset:offset) .eq. ' ') goto 100
      str = numstr(offset:30)

      return
      end


      subroutine NullOp(idum)
      return
      end
