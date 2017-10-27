c     Copyright(C) 2007
c     Macromolecular Structure Gropu of the Biochemistry
c     Department at the University of California, San
c     Francisco.

c     Declares constants for the event loop functions.

      integer IVE_IOEVENT_READ
      parameter (IVE_IOEVENT_READ = 1)
      integer IVE_IOEVENT_WRITE
      parameter (IVE_IOEVENT_WRITE = 2)
      integer IVE_IOEVENT_EXCPT
      parameter (IVE_IOEVENT_EXCPT = 4)

      integer IVE_MESSAGE_MIN
      parameter (IVE_MESSAGE_MIN = 1)

      integer IVE_MESSAGE_MAX
      parameter (IVE_MESSAGE_MAX = 255)
