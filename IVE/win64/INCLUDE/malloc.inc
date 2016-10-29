c     @(#) $Id: malloc.inc_m4,v 1.1 2001/11/20 20:02:11 eric Exp $
c     $Name:  $
c
c     Define appropriate return type for malloc and size of pointer in
c     bytes.

      integer iptrsz
      parameter (iptrsz = 8)

      integer*8 malloc
