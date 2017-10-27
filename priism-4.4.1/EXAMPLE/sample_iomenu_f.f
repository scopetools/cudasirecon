c     Demonstrates the C/C++ interface to the iomenu utility.
c
c     Assume that you had a command-line application which expected
c     the following command-line:
c
c         ffilter in out optional_otf -x=start:end -y=start:end \
c             -z=start:end -t=start:end:step -w=wave1:wave2:... \
c             -bandpass=l1:l2 -enhance=e -smooth=s -fraction=f \
c             -apodize=n -nzpad=n -2D -norm -orig -mode2
c     (this if the command-line used by the ffilter application in
c     Priism), and want to write a graphical user interface as a
c     wrapper.  The iomenu utility will automatically handle the
c     region selection arguments (-x, -y, -z, -t, -w); the code
c     below will handle the others.  If you want to see what the
c     corresponding graphical user interface wrapper in Priism
c     looks like, go to the Processing->Fourier Tools->FFilter
c     entry in the Priism menu bar.


      program ffilteri

c     Pick up the definition of the iomenu_struct structure and
c     the declaration of iom as an instance of that structure.
      include 'iomenu.inc'
      character*32 spclab
      integer i

c     Set up the iomenu_struct structure; use the iom variable
c     declared in iomenu.inc since it is handy.

c     The first entry on the command line is the input image file.
c     I use this file to set the default region and output file names
c     so the type is 1 (other input files are 2).  Set the file filter
c     to '*' (accept any name) and the file name to " " (do not know a
c     value yet).  The label shown in the user interface will be
c     "InFile".
      iom.filelabels(1) = 'InFile'
      iom.ifiletype(1) = 1
      iom.filefilter(1) = '*'
      iom.filenames(1) = ' '

c     The second entry on the command line is the output image file.c
c     I want the default name of this file to be set from the input
c     file so I use a type of 3 (for other output files use a type of
c     4); the default name will be the name of the input file plus
c     the contents of the filefilter element, "_flt" in this case.
c     Set the initial file name to be empty; once the user selects
c     an input file this value will be overwritten.  The label that
c     will appear in the user interface is "OutFile".
      iom.filelabels(2) = 'OutFile'
      iom.ifiletype(2) = 3
      iom.filefilter(2) = '_flt'
      iom.filenames(2) = ' '

c     The third entry on the command line is an optional input file
c     which will have the special value, "none", if is is not
c     specified.  If I could change the command-line application,
c     this optional file would be better handle via a command-line
c     switch, but sicne I can not, I'll declare is as type (output
c     without a special default name) so that iomenu does not check
c     that the file exists (which it would if I used type 2).  The
c     initial value for the file name is set to "none" and the file
c     filter is set to "*" (though the value will not be used).  The
c     label that appears in the user interface will be "OTFfile".
      iom.filelabels(3) = 'OTFfile'
      iom.ifiletype(3) = 4
      iom.filefilter(3) = '*'
      iom.filenames(3) = 'none' // char(0)

c     Terminate the list of required files.
      iom.ifiletype(4) = 0

c     -bandpass= takes 2 values (iom.itypes(1,1) is 2) which
c     are floating-point numbers (iom.itypes(2,1) is 2).  Initialize
c     the two values to zero.
      iom.labels(1) = "-bandpass="
      iom.itypes(1,1) = 2
      iom.itypes(2,1) = 2
      iom.rvals(1) = 0.0
      iom.rvals(2) = 0.0

c     -enhance= takes 1 value (iom.itypes(1,2) is 1) which is
c     a floating-point number (iom.itypes(2,2) is 2).  Initialize
c     the value to 0.9.
      iom.labels(2) = "-enhance="
      iom.itypes(1,2) = 1
      iom.itypes(2,2) = 2
      iom.rvals(3) = 0.9

c     -smooth= takes 1 value (iom.itypes(1,3) is 1) which is
c     a floating-point number (iom.itypes(2,3) is 2).  Initialize
c     the value to 0.15.
      iom.labels(3) = "-smooth="
      iom.itypes(1,3) = 1
      iom.itypes(2,3) = 2
      iom.rvals(4) = 0.15

c     -fraction= takes 1 value (iom.itypes(1,4) is 1) which is
c     a floating-point number (iom.itypes(2,4) is 2).  Initialize
c     the value to 0.15.
      iom.labels(4) = "-fraction="
      iom.itypes(1,4) = 1
      iom.itypes(2,4) = 2
      iom.rvals(5) = 0.0

c     -apodize= takes 1 value (iom.itypes(1,5) is 1) which is
c     an integer (iom.itypes(2,5) is 1).  Initialize the value to
c     zero.
      iom.labels(5) = "-apodize="
      iom.itypes(1,5) = 1
      iom.itypes(2,5) = 1
      iom.ivals(1) = 0

c     -apodize= takes 1 value (iom.itypes(1,6) is 1) which is
c     an integer (iom.itypes(2,6) is 1).  Initialize the value to
c     zero (do no know the z size yet so set the padded size to
c     zero).
      iom.labels(6) = "-nzpad="
      iom.itypes(1,6) = 1
      iom.itypes(2,6) = 1
      iom.ivals(2) = 0

c     Because iom.itypes(2,7) is 3, print out the -2D switch if
c     the corresponding integer is nonzero.  iom.itypes(1,7)
c     could be any positive integer, but use 1 since require
c     one integer from iom.ivals.  Initially, 2D processing
c     is not forced so set iom.ivals(3) to zero.
      iom.labels(7) = "-2D"
      iom.itypes(1,7) = 1
      iom.itypes(2,7) = 3
      iom.ivals(3) = 0

c     Because iom.itypes(2,8) is 3, print out the -norm switch
c     if the corresponding integer is nonzero.  iom.itypes(1,8)
c     could be any positive integer, but use 1 since require
c     one integer from iom.ivals.  Initially, CTF normalization
c     is on so set iom.ivals(4) to a nonzero value.
      iom.labels(8) = "-norm"
      iom.itypes(1,8) = 1
      iom.itypes(2,8) = 3
      iom.ivals(4) = 1

c     Because iom.itypes(2,9) is 3, print out the -orig switch
c     if the corresponding integer is nonzero.  iom.itypes(1,9)
c     could be any positive integer, but use 1 since require
c     one integer from iom.ivals.  Initially, preservation of
c     the mean is on so set iom.ivals(5) to a nonzero value.
      iom.labels(9) = "-orig"
      iom.itypes(1,9) = 1
      iom.itypes(2,9) = 3
      iom.ivals(5) = 1

c     Because iom.itypes(2,10) is 3, print out the -mode2 switch
c     if the corresponding integer is nonzero.  iom.itypes(1,10)
c     could be any positive integer, but use 1 since require
c     one integer from iom.ivals.  Initially, forcing floating-
c     point output is off so set iom.ivals(6) to zero.
      iom.labels(10) = "-mode2"
      iom.itypes(1,10) = 1
      iom.itypes(2,10) = 3
      iom.ivals(6) = 0

c     Terminate the list of switches.
      iom.itypes(1,11) = 0

c     Use elements of iuser for user interface IDs so that
c     it is possible to hide, redisplay, or update items
c     in response to user input.  Until the controls are
c     created, set the values to zero.  iom.iuser(1) is
c     the bandpass help button id.  iom.iuser(2) is the
c     bandpass field id.  iom.iuser(3) is the enhancement
c     help button id.  iom.iuser(4) is the the enhancement
c     field id.  iom.iuser(5) is the smoothing help
c     button id.  iom.iuser(6) is the smoothing field id.
c     iom.iuser(7) is the padded z size help button id.
c     iom.iuser(8) is the padded z size field.
      iom.iuser(1) = 0
      iom.iuser(2) = 0
      iom.iuser(3) = 0
      iom.iuser(4) = 0
      iom.iuser(5) = 0
      iom.iuser(6) = 0
      iom.iuser(7) = 0
      iom.iuser(8) = 0

c     Use iom.iuser(9) to hold the minimum padded size
c     allowed.  Until know the size of the z region, use
c     zero.
      iom.iuser(9) = 0

c     Alternatively, could compile the program in such a way so that
c     the iom structure is initialized to zeros.
      do i = 1, 9
         iom.ixyztw(1,i) = 0
         iom.ixyztw(2,i) = 0
      enddo
      do i = 1, 5
         iom.nxyztw(i) = 0
      enddo

c     The first argument is the title of the dialog.  The second
c     is the label for the button to open the child dialog for
c     special parameters.  The third is the name of the
c     command-line executable.  The fourth is the iomenu_struct
c     structure.  Use a character array with some extra space as
c     the second argument rather than a string literal to be safe
c     (since it will be used as the label for a function button
c     created by WMAddFuncButton and iomenu automatically adds
c     a WMAddSaveButton control, if the label were ever to change,
c     the WM library would write to the label when the user
c     restored the dialog settings from an older save file).
      spclab = "Special Parameters"
      call iomenu('2D & 3D Fourier Filter Setup Menu',
     &            spclab,
     &            'ffilter',
     &            iom)

      end


c     Required routine for iomenu.  Adds controls to the main dialog.
      subroutine personal(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'

      integer i
      character*32 label

c     Store fixed bounds for input parameters.
      integer i1, i2
      real fzero, fone, f1k
      save i1, i2, fzero, fone, f1k

c     Some compilers (namely the Absoft compiler on the Mac) pass
c     literal strings by making a temporary copy.  Since the help
c     lookup values have to exist as long as the help button
c     does, create static character arrays for the help lookup
c     values.
      character bndkey*20, enhkey*20, smtkey*20, frakey*16, frckey*20
      save bndkey, enhkey, smtkey, frakey, frckey

      integer handle_dummy, handle_bandpass, handle_force2d
      external handle_dummy, handle_bandpass, handle_force2d

      integer WMAddFloatField, WMAddInfoButton
      integer WMAddToggleButton, WMNewRow, WMSetOffset
      external WMAddFloatField, WMAddInfoButton
      external WMAddToggleButton, WMNewRow, WMSetOffset

      data i1, i2 / 1, 2 /
      data fzero, fone, f1k / 0.0, 1.0, 1000.0 /

      i = WMNewRow()
      bndkey = 'ffilter bandpass' // char(0)
      label = 'Bandpass (res min,max)' // char(0)
      iom.iuser(1) = WMAddInfoButton(label, bndkey)
      iom.iuser(2) = WMAddFloatField(iom.rvals(1), 2, 10, 20,
     &    fzero, f1k, i1, handle_bandpass, iom, 0, 0)

      i = WMNewRow()
      enhkey = 'ffilter enhancement' // char(0)
      label = 'Enhancement (0-1)' // char(0)
      iom.iuser(3) = WMAddInfoButton(label, enhkey)
      iom.iuser(4) = WMAddFloatField(iom.rvals(3), 1, 10, 10,
     &    fzero, fone, i2, handle_dummy, 0, 0, 0)

      i = WMNewRow()
      smtkey = 'ffilter smoothing' // char(0)
      label = 'Smoothing   (0-1)' // char(0)
      iom.iuser(5) = WMAddInfoButton(label, smtkey)
      iom.iuser(6) = WMAddFloatField(iom.rvals(4), 1, 10, 10,
     &    fzero, fone, i2, handle_dummy, 0, 0, 0)

      i = WMNewRow()
      frakey = 'ffilter fract' // char(0)
      label = 'Fraction original image (0-1)' // char(0)
      i = WMAddInfoButton(label, frakey)
      i = WMAddFloatField(iom.rvals(5), 1, 8, 8,
     &     fzero, fone, i2, handle_dummy, 0, 0, 0)

      i = WMNewRow()
      frckey = 'ffilter twod' // char(0)
      label = 'Force Processing as 2D data' // char(0)
      i = WMAddInfoButton(label, frckey)
      label = char(0)
      i = WMAddToggleButton(label, iom.ivals(3),
     &     handle_force2d, iom, 0, 0)

c     Force the special parameters button to be offset from the
c     parameter controls and the left side of the dialog.
      i = WMSetOffset(100, 0, 20, 0)

      return
      end


c     Required routine for iomenu.  Adds controls to the special
c     parameters dialog.
      subroutine special_parameters(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'

      integer i
      character*32 label

c     Store fixed bounds for input parameters
      integer i0, i1, i2, imx
      save i0, i1, i2, imx

c     Some compilers (namely the Absoft compiler on the Mac) pass
c     literal strings by making a temporary copy.  Since the help
c     lookup values have to exist as long as the help button
c     does, create static character arrays for the help lookup
c     values.
      character rllkey*16, padkey*16, nrmkey*16, orgkey*16, modkey*16
      save rllkey, padkey, nrmkey, orgkey, modkey

      integer handle_dummy
      external handle_dummy

      integer WMAddInfoButton, WMAddIntField, WMAddToggleButton
      integer WMNewRow, WMSetOffset
      external WMAddInfoButton, WMAddIntField, WMAddToggleButton
      external WMNewRow, WMSetOffset

      data i0, i1, i2, imx / 0, 1, 2, 2147483647 /

      i = WMNewRow()
      rllkey = 'ffilter apodize' // char(0)
      label = 'Border rolloff, # pixels' // char(0)
      i = WMAddInfoButton(label, rllkey)
      i = WMAddIntField(iom.ivals(1), 1, 12, 12, i0, imx, i1,
     &    handle_dummy, 0, 0, 0)

      i = WMNewRow()
      padkey = 'ffilter nzpad' // char(0)
      label = 'Size for Z transforms' // char(0)
      iom.iuser(7) = WMAddInfoButton(label, padkey)
      iom.iuser(8) = WMAddIntField(iom.ivals(2), 1, 12, 12,
     &    iom.iuser(9), imx, i2, handle_dummy, 0, 0, 0)

      i = WMNewRow()
      nrmkey = 'ffilter norm' // char(0)
      label = 'Normalize CTF' // char(0)
      i = WMAddInfoButton(label, nrmkey)
      label = char(0)
      i = WMAddToggleButton(label, iom.ivals(4), handle_dummy, 0, 0, 0)
      i = WMSetOffset(20, 0, 0, 0)
      orgkey = 'ffilter orig' // char(0)
      label = 'Preserve mean' // char(0)
      i = WMAddInfoButton(label, orgkey)
      i = WMSetOffset(0, 0, 0, 0)
      label = char(0)
      i = WMAddToggleButton(label, iom.ivals(5), handle_dummy, 0, 0, 0)

      i = WMNewRow()
      modkey = 'ffilter mode2' // char(0)
      label = 'Force floating-point' // char(0)
      i = WMAddInfoButton(label, modkey)
      label = char(0)
      i = WMAddToggleButton(label, iom.ivals(6), handle_dummy, 0, 0, 0)

c     Force the special parameters button to be offset from the
c     parameter controls and the left side of the dialog.
      i = WMSetOffset(70, 0, 20, 0)

      return
      end


c     Required routine for iomenu.  Called when the names of
c     required files change or the region bounds change.
      subroutine custom_ifile(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'
      integer mz, nc, i
c     lenstring is in libimcompat.
c     WMPasteWidget, WMUnpasteWidget, and WMUpdateField in libWM.
      integer lenstring, WMPasteWidget, WMUnpasteWidget, WMUpdateField
      external lenstring, WMPasteWidget, WMUnpasteWidget, WMUpdateField

c     If the first file (main input) changed or the z
c     region changed, update the padded size.
      if (iom.iparm .eq. 1 .or. iom.iparm .eq. -3) then
         mz = iom.ixyztw(2,3) - iom.ixyztw(1,3) + 1
         iom.iuser(9) = mz
         iom.ivals(2) = 2
         do while (iom.ivals(2) .lt. 1.3 * mz)
             iom.ivals(2) = iom.ivals(2) * 2
         enddo
         if (iom.iuser(8) .ne. 0) then
             i = WMUpdateField(iom.iuser(8))
         endif
      else if (iom.iparm .eq. 3) then
c         If the third file (optional ctf) changed, display the
c         bandpass widgets if no CTF has been supplied or hide
c         them if not using a CTF (the bandpass filter is not
c         used in that case).
c
c         It is probably better to grey out the unused controls
c         (i.e  WMDisableField and WMEnableField rather than
c         WMUnpasteWidget and WMPasteWidget) because pasting and
c         unpasting can cause geometry problems if the dialog has
c         not been realized.  Hide and display the controls
c         because that is what ffilter_i has traditionally done.
          nc = lenstring(iom.filenames(3))
          if (nc .eq. 0 .or. (nc .eq. 4 .and.
     &        iom.filenames(3)(1:4) .eq. 'none')) then
              if (iom.iuser(1) .ne. 0) then
                  i = WMPasteWidget(iom.iuser(1))
                  i = WMPasteWidget(iom.iuser(2))
              endif
          else
              if (iom.iuser(1) .ne. 0) then
                  i = WMUnpasteWidget(iom.iuser(1))
                  i = WMUnpasteWidget(iom.iuser(2))
              endif
          endif
      endif

      return
      end


c     Required routine for iomenu.  Called when the user
c     presses the Exit button.
      subroutine call_exit(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'

c     Simply stop.
      stop
      end


c     If the bandpass filter is enabled, hide the enhancement
c     and smoothing controls (since they are not relevant).
c     Otherwise, display the enhancement and smoothing controls.
c
c     It is probably better to grey out the unused controls
c     (i.e  WMDisableField and WMEnableField rather than
c     WMUnpasteWidget and WMPasteWidget) because pasting and
c     unpasting can cause geometry problems if the dialog has not
c     been realized.  Hide and display the controls because
c     that is what ffilter_i has traditionally done.
      integer function handle_bandpass(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'

      integer i

      integer WMPasteWidget, WMUnpasteWidget
      external WMPasteWidget, WMUnpasteWidget

      if (iom.rvals(1) + iom.rvals(2) .ge. 1e-4) then
         i = WMUnpasteWidget(iom.iuser(3))
         i = WMUnpasteWidget(iom.iuser(4))
         i = WMUnpasteWidget(iom.iuser(5))
         i = WMUnpasteWidget(iom.iuser(6))
      else
         i = WMPasteWidget(iom.iuser(3))
         i = WMPasteWidget(iom.iuser(4))
         i = WMPasteWidget(iom.iuser(5))
         i = WMPasteWidget(iom.iuser(6))
      endif

      handle_bandpass = 0
      return
      end


c     Define a dummy callback for WM calls.
      integer function handle_dummy(i)
      integer i

      handle_dummy = 0
      return
      end


c     If forcing 2D processing, hide the padding controls (they
c     are not relevant); otherwise, display the padding controls.
c
c     It is probably better to grey out the unused controls
c     (i.e  WMDisableField and WMEnableField rather than
c     WMUnpasteWidget and WMPasteWidget) because pasting and
c     unpasting can cause geometry problems if the dialog has not
c     been realized.  Hide and display the controls because
c     that is what ffilter_i has traditionally done.
      integer function handle_force2d(iom)

c     Pick up definition of iomenu_struct.  Declare iom.
      include 'iomenu.inc'

      integer i

      integer WMPasteWidget, WMUnpasteWidget
      external WMPasteWidget, WMUnpasteWidget

      if (iom.ivals(3) .ne. 0) then
         if (iom.iuser(7) .ne. 0) then
            i = WMUnpasteWidget(iom.iuser(7))
            i = WMUnpasteWidget(iom.iuser(8))
         endif
      else
         if (iom.iuser(7) .ne. 0) then
            i = WMPasteWidget(iom.iuser(7))
            i = WMPasteWidget(iom.iuser(8))
         endif
      endif

      handle_force2d = 0
      return
      end
