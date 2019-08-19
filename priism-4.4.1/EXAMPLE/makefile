## Uncomment out the macro definitions appropriate for the operating
## system and compiler combination you use.  The optimization and warning flags
## are just suggestions and can be changed to match what you want.  If IVE_BASE
## (the root directory of the Priism distribution) is not set via the
## environment you will need to define it.  You may also need to modify the
## search paths for X and OpenGL.


## IRIX
## MipsPro compilers; also possible to use -mips4 but will not run on R4000
## based systems.
# CC = /usr/bin/cc
# CFLAGS_ARCH = -n32 -mips3
# CFLAGS_PIC = -KPIC
# CFLAGS_FULLWARN = -fullwarn
# CFLAGS_STDOPT = -O

# FC = /usr/bin/f77
# FFLAGS_ARCH = -n32 -mips3
# FFLAGS_PIC = -KPIC
# FFLAGS_FULLWARN = -fullwarn
# FFLAGS_STDOPT = -O2
# FFLAGS_STATIC = -static
# FFLAGS_EXTSRC = -extend_source
# FLINKFLAGS_NOFORMAIN =
# FUNIX_LIB =

## For IRIX, the X and OpenGL libraries and headers are in the default
## search paths; so unless you want another version, there is no need to set
## X_INC_DIR_FLAG, X_LIB_DIR_FLAG, GL_INC_DIR_FLAG, and GL_LIB_DIR_FLAG.
# X_INC_DIR_FLAG =
# X_LIB_DIR_FLAG =
# GL_INC_DIR_FLAG =
# GL_LIB_DIR_FLAG =
# MATH_LIB = -lm
# IVE_INC_DIR = $(IVE_BASE)/IRIX/INCLUDE
# IVE_LIB_DIR = $(IVE_BASE)/IRIX/LIB

## When explicitly linking against -lXm, you made need to specify additional
## libraries (besides libXt).  Using IRIX's libXm with the IRIX linker does
## not require anything extra.
# MOTIF_INDIRECT_DEP =


## Linux (x86 or x86_64)
## gcc
# CC = /usr/bin/cc
# CFLAGS_PIC = -fpic
# CFLAGS_FULLWARN = -Wall -W -Wpointer-arith -Wmissing-declarations
# CFLAGS_STDOPT = -O

## The Portland Group's C compiler
# CC = $(PGI)/linux86/bin/pgcc
# CFLAGS_PIC = -fpic
# CFLAGS_FULLWARN = -Minform,inform
# CFLAGS_STDOPT = -O

## Intel's C compiler
# CC = icc
# CFLAGS_PIC = -KPIC
# CFLAGS_FULLWARN = -w2
# CFLAGS_STDOPT = -O

## g77 (also look for IVE_LIB_DIR_FLAG for the modifications needed there)
# FC = /usr/bin/g77
# FFLAGS_PIC = -fpic
# FFLAGS_FULLWARN = -Wall -W -Wsurprising
# FFLAGS_STDOPT = -O
# FFLAGS_STATIC = -fno-automatic
# FFLAGS_EXTSRC = -ffixed-line-length-132
# FFLAGS_PASS_TO_LINKER = -Wl
# FLINKFLAGS_NOFORMAIN =
# FUNIX_LIB =

## The Portland Group's Fortran compiler
# FC = $(PGI)/linux86/bin/pgf77
# FFLAGS_PIC = -fpic
# FFLAGS_FULLWARN = -Minform,inform
# FFLAGS_STDOPT = -O
# FFLAGS_STATIC = -Msave
# FFLAGS_EXTSRC = -Mextend
# FFLAGS_PASS_TO_LINKER = -Wl
# FLINKFLAGS_NOFORMAIN = -Mnomain
# FUNIX_LIB =

## Intel's Fortran compiler (note that if you are using this on an x86
## rather than x86_64 system, Priism does not include versions of libimcompat,
## libiomenu, or libifft that are compatible with the Intel compiler).
## For version 7 of the compiler FC needs to be ifc rather than ifort;
## do not know what is needed for version 8.
# FC = ifort
# FFLAGS_PIC = -KPIC
# FFLAGS_FULLWARN = -W1
# FFLAGS_STDOPT = -O
# FFLAGS_STATIC = -save -zero
# FFLAGS_EXTSRC = -extend_source
# FFLAGS_PASS_TO_LINKER = -Qoption,link
## -nofor_main is for version 9 of the compiler (not needed for 7; do not
## know about version 8).
# FLINKFLAGS_NOFORMAIN = -nofor_main
## For version 7 of the compiler, FUNIX_LIB needs to be -posixlib -Valib;
## do not know what is needed for version 8.
# FUNIX_LIB = -lifport

## The GNU linker does not look at -L options to set the search path for
## shared libraries that are used indirectly; instead it uses -rpath-link.
# CFLAGS_RLINK = -Wl,-rpath-link,$(RLINK_PATH)
# FFLAGS_RLINK = $(FFLAGS_PASS_TO_LINKER),-rpath-link,$(RLINK_PATH)

## If X is not in /usr (older distributions put those files in /usr/X11R6,
## for instance), you will need to change the definitions of X_INC_DIR_FLAG
## and X_LIB_DIR.  If the OpenGL headers are not in the same place as the X
## headers, define GL_INC_DIR_FLAG.  If the OpenGL libraries are not in the
## same place as the X libraries, change GL_LIB_DIR and define
## GL_LIB_DIR_FLAG.  On some distributions, /usr/lib has the libraries
## appropriate for the current hardware (i.e. 64-bit if it's a 64-bit system
## and 32-bit if it's a 32-bit system), /usr/lib32 has the 32-bit libraries,
## and /usr/lib64 has the 64-bit libraries.  Other distributions put 32-bit
## libraries in /usr/lib and 64-bit libraries in /usr/lib64.  Depending on
## your distribution and whether you are compiling for 32- or 64-bit, you
## may have to change X_LIB_DIR.
# X_INC_DIR_FLAG = -I/usr/include
# X_LIB_DIR = /usr/lib
# X_LIB_DIR_FLAG = -L$(X_LIB_DIR)
# GL_INC_DIR_FLAG =
# GL_LIB_DIR = $(X_LIB_DIR)
# GL_LIB_DIR_FLAG =
# MATH_LIB = -lm
## These two are for the x86 Linux version of Priism.
# IVE_INC_DIR = $(IVE_BASE)/Linux/x86/INCLUDE
# IVE_LIB_DIR = $(IVE_BASE)/Linux/x86/LIB
## These two are for the x86_64 Linux version of Priism.
# IVE_INC_DIR = $(IVE_BASE)/Linux/x86_64/INCLUDE
# IVE_LIB_DIR = $(IVE_BASE)/Linux/x86_64/LIB
# RLINK_PATH = $(IVE_LIB_DIR):$(GL_LIB_DIR):$(X_LIB_DIR)

## When explicitly linking against the Motif/LessTif librar (-lXm), you may
## need to specify additional libraries (besides libXt).  Using OpenMotif or
## LessTif with the GNU linker does not appear to require anything extra.
# MOTIF_INDIRECT_DEP =


## Mac OS X (PowerPC or Intel)
## Apple's developer tools
# CC = /usr/bin/cc
## This value of CFLAGS_ARCH is only for compiling 64-bit applications on Intel
## Macs.  To compile 32-bit applications on OS X 10.6 or later you'd want to
## use -arch i386 .
# CFLAGS_ARCH = -arch x86_64
# CFLAGS_PIC = -dynamic
# CFLAGS_FULLWARN =  -Wall -W -Wpointer-arith -Wmissing-declarations
# CFLAGS_STDOPT = -O

## Absoft Pro Fortran (pass -V so it will accept VAX tab-format style
## files; -f and -N15 force Unix-style mangling of names in the object files).
# FC = $(ABSOFT)/bin/f77 -V -f -N15
# FFLAGS_PIC =
# FFLAGS_FULLWARN =
# FFLAGS_STDOPT = -O
# FFLAGS_STATIC = -s
# FFLAGS_EXTSRC = -W
# FLINKFLAGS_NOFORMAIN =
# FUNIX_LIB = -lU77

## IBM XL Fortran (pass -qextname to force Unix-style mangling of names in
## the object files; also look for IVE_LIB_DIR_FLAG for the modifications
## needed there).
# FC = /opt/ibmcmp/xlf/8.1/bin/xlf -qextname
# FFLAGS_PIC = -qpic
# FFLAGS_FULLWARN = -qflag=i:i
# FFLAGS_STDOPT = -O3
# FFLAGS_STATIC = -qsave=all
# FFLAGS_EXTSRC = -qfixed=132
# FLINKFLAGS_NOFORMAIN =
# FUNIX_LIB =

## Intel's Fortran compiler
# FC = ifort
# FFLAGS_PIC =
# FFLAGS_FULLWARN = -W1
# FFLAGS_STDOPT = -O
# FFLAGS_STATIC = -save -zero
# FFLAGS_EXTSRC = -extend_source
# FFLAGS_PASS_TO_LINKER = -Wl
# FLINKFLAGS_NOFORMAIN = -nofor_main
# FUNIX_LIB = -lifport -lifcore -limf -lirc

## The IVE libraries were linked without hardwired paths to the other IVE
## libraries or to the Motif library; this forces you to specify to the linker
## how to find indirect dependencies.  Below is one way; you could also
## explicitly list the indirect dependencies when linking.
# CFLAGS_RLINK = -Wl,-dylib_file,libGLw.dylib:$(IVE_LIB_DIR)/libGLw.dylib -Wl,-dylib_file,libive.dylib:$(IVE_LIB_DIR)/libive.dylib -Wl,-dylib_file,libIWL.dylib:$(IVE_LIB_DIR)/libIWL.dylib -Wl,-dylib_file,libWM.dylib:$(IVE_LIB_DIR)/libWM.dylib -Wl,-dylib_file,libXm.2.dylib:$(IVE_LIB_DIR)/libXm.2.dylib
## This one is for the Absoft compiler.
# FFLAGS_RLINK = -X -dylib_file -X libGLw.dylib:$(IVE_LIB_DIR)/libGLw.dylib -X -dylib_file -X libive.dylib:$(IVE_LIB_DIR)/libive.dylib -X -dylib_file -X libIWL.dylib:$(IVE_LIB_DIR)/libIWL.dylib -X -dylib_file -X libWM.dylib:$(IVE_LIB_DIR)/libWM.dylib -X -dylib_file -X libXm.2.dylib:$(IVE_LIB_DIR)/libXm.2.dylib
## This one is for the IBM XL Fortran or Intel Fortran compilers.
# FFLAGS_RLINK = $(CFLAGS_RLINK)

## If X is not in /usr/X11R6, you will need to change the definitions
## of X_INC_DIR_FLAG and X_LIB_DIR_FLAG.  If the OpenGL headers are not
## in the same place as the X headers, define GL_INC_DIR_FLAG.  If the
## OpenGLlibraries are not in the same place as the X libraries, define
## GL_LIB_DIR_FLAG.
# X_INC_DIR_FLAG = -I/usr/X11R6/include
# X_LIB_DIR_FLAG = -L/usr/X11R6/lib
# GL_LIB_DIR =
# GL_LIB_DIR_FLAG =
# MATH_LIB = 
## These two are for the 32-bit version (PowerPC or Intel).
# IVE_INC_DIR = $(IVE_BASE)/Darwin/INCLUDE
# IVE_LIB_DIR = $(IVE_BASE)/Darwin/LIB
## These two are for the 64-bit version (Intel only; Mac OS X 10.5 or later).
# IVE_INC_DIR = $(IVE_BASE)/Darwin64/INCLUDE
# IVE_LIB_DIR = $(IVE_BASE)/Darwin64/LIB


## When explicitly linking against -lXm, you may need to specify additional
## libraries (besides libXt).  On Mac OS X, libXp and libXext appear to be
## required with either OpenMotif or LessTif.
# MOTIF_INDIRECT_DEP = -lXp -lXext


IVE_INC_DIR_FLAG = -I$(IVE_INC_DIR)
## If you are using the IBM XL Fortran compiler for Mac OS X or g77,
## uncomment out the one of the definitions for IVE_LIB_DIR_FLAG below (the
## first is for the IBM XL Fortran compiler, and the third is for g77) and
## comment out the normal definition of IVE_LIB_DIR_FLAG.
# IVE_LIB_DIR_FLAG = -L$(IVE_LIB_DIR)/xlf -L$(IVE_LIB_DIR)
# IVE_LIB_DIR_FLAG = -L$(IVE_LIB_DIR)/g77 -L$(IVE_LIB_DIR)
IVE_LIB_DIR_FLAG = -L$(IVE_LIB_DIR)

MODULES = \
	ArrowButton \
	CenterImg \
	Circle \
	createmrc \
	Cross \
	DrawBox \
	FindBorder \
	FindBorderF \
	interpo \
	InvertPseudo \
	Point \
	ProcWinEvent \
	ResetLUT \
	sample_iomenu_c \
	sample_iomenu_f \
	ScrollText \
	Select \
	simple_menu \
	TestWMOverlay \
	Trace

CFLAGS_ALL = $(IVE_INC_DIR_FLAG) $(GL_INC_DIR_FLAG) $(X_INC_DIR_FLAG) \
	$(CFLAGS_ARCH) $(CFLAGS_PIC) $(CFLAGS_FULLWARN) $(CFLAGS_STDOPT) \
	$(CFLAGS)

CLINKFLAGS_ALL = $(CFLAGS_ARCH) $(CFLAGS_PIC) $(CFLAGS_FULLWARN) \
	$(CFLAGS_STDOPT) $(CFLAGS_RLINK) $(CLINKFLAGS_NAMESPACE) $(CLINKFLAGS)

FFLAGS_ALL = $(IVE_INC_DIR_FLAG) $(FFLAGS_ARCH) $(FFLAGS_PIC) \
	$(FFLAGS_FULLWARN) $(FFLAGS_STDOPT) $(FFLAGS_EXTSRC) $(FFLAGS)

FLINKFLAGS_ALL = $(FFLAGS_ARCH) $(FFLAGS_PIC) $(FFLAGS_FULLWARN) \
	$(FFLAGS_STDOPT) $(FFLAGS_RLINK) $(FLINKFLAGS_NAMESPACE) $(FLINKFLAGS)


.c.o:
	$(CC) -c $(CFLAGS_ALL) $<

.f.o:
	$(FC) -c $(FFLAGS_ALL) $<


all: $(MODULES)

clean:
	-rm -f *.o core

clobber: clean
	-rm -f $(MODULES)

ArrowButton: ArrowButton.o
ScrollText: ScrollText.o
simple_menu: simple_menu.o
ArrowButton \
ScrollText \
simple_menu \
:
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) -lWM

CenterImg: CenterImg.o
Circle: Circle.o
ProcWinEvent: ProcWinEvent.o
Trace: Trace.o
CenterImg \
Circle \
ProcWinEvent \
Trace \
:
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) -lWM -lIWL

Cross: Cross.o
Point: Point.o
Cross \
Point \
:
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o \
		$(IVE_LIB_DIR_FLAG) -lWM -lIWL $(MATH_LIB)

createmrc: createmrc.o
InvertPseudo: InvertPseudo.o
ResetLUT: ResetLUT.o
createmrc \
InvertPseudo \
ResetLUT \
:
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) -lIWL

FindBorder: FindBorder.o
Select: Select.o
FindBorder \
Select \
:
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) \
		-lProcFunc -lWM -lIWL

DrawBox: DrawBox.o
	$(FC) $(FLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) -lIWL \
		$(FUNIX_LIB)

# Need FLINKFLAGS_NOFORMAIN since the main program (which is in libFProcFunc)
# is not a Fortran main program.  Use FProcFunc rather than ProcFunc since
# IPAppSpecifics is written in Fortran.
FindBorderF: FindBorderF.o
	$(FC) $(FLINKFLAGS_NOFORMAIN) $(FLINKFLAGS_ALL) \
		-o $@ $@.o $(IVE_LIB_DIR_FLAG) -lFProcFunc -lWM -lIWL

# Use the Fortran compiler to link since that makes it easier to pull in
# any Fortran libraries used by gettransmatrix.f.  Need FLINKFLAGS_NOFORMAIN
# since the main program (which is in libProcFunc) is not a Fortran main
# program.
interpo: interpo.o gettransmatrix.o
	$(FC) $(FLINKFLAGS_NOFORMAIN) $(FLINKFLAGS_ALL) \
		-o $@ $@.o gettransmatrix.o \
		$(IVE_LIB_DIR_FLAG) -lProcFunc -lWM -lIWL

sample_iomenu_c: sample_iomenu_c.o
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) \
		-liomenu -lWM -lIWL -live

sample_iomenu_f: sample_iomenu_f.o
	$(FC) $(FLINKFLAGS_ALL) -o $@ $@.o $(IVE_LIB_DIR_FLAG) \
		-liomenu -lWM -lIWL -live -limcompat $(FUNIX_LIB)

TestWMOverlay: TestWMOverlay.o
	$(CC) $(CLINKFLAGS_ALL) -o $@ $@.o \
		$(IVE_LIB_DIR_FLAG) -lWM -lGLw \
		$(X_LIB_DIR_FLAG) -lXm -lXt $(MOTIF_INDIRECT_DEP) -lX11 -lGL
