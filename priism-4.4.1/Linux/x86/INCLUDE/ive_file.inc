c     Copyright(C) 2002
c     Macromolecular Structure Group of the Biochemistry Department
c     at the University of California, San Francisco

c     Declares constants for IVECheckFileAttributes.

      integer IVE_FILE_CHECK_EXIST
      parameter (IVE_FILE_CHECK_EXIST = 1)

      integer IVE_FILE_CHECK_DIR
      parameter (IVE_FILE_CHECK_DIR = 2)

      integer IVE_FILE_CHECK_ACCESS
      parameter (IVE_FILE_CHECK_ACCESS = 4)

      integer IVE_FILE_CHECK_ALL
      parameter (IVE_FILE_CHECK_ALL = 7)

      integer IVE_FILE_INDETERM
      parameter (IVE_FILE_INDETERM = 1)

      integer IVE_FILE_EXIST
      parameter (IVE_FILE_EXIST = 2)

      integer IVE_FILE_DIR
      parameter (IVE_FILE_DIR = 4)

      integer IVE_FILE_READ
      parameter (IVE_FILE_READ = 8)

      integer IVE_FILE_WRITE
      parameter (IVE_FILE_WRITE = 16)

      integer IVE_FILE_EXEC
      parameter (IVE_FILE_EXEC = 32)
