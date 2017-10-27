


      parameter (KeyPressMask = 1)
      parameter (KeyReleaseMask = 2)
      parameter (ButtonPressMask = 4)
      parameter (ButtonReleaseMask = 8)
      parameter (EnterWindowMask = 16)
      parameter (LeaveWindowMask = 32)
      parameter (PointerMotionMask = 64)
      parameter (ButtonMotionMask = 4096)
      parameter (ExposureMask = 32768)

      parameter (WM_UP_ARROW = 0)
      parameter (WM_DOWN_ARROW = 1)
      parameter (WM_LEFT_ARROW = 2)
      parameter (WM_RIGHT_ARROW = 3)
      parameter (WM_OVERLAY_NO = 0)
      parameter (WM_OVERLAY_YES = 1)

      parameter (IW_BLACK = 0)
      parameter (IW_RED = 1)
      parameter (IW_GREEN = 2)
      parameter (IW_YELLOW = 3)
      parameter (IW_BLUE = 4)
      parameter (IW_PINK = 5)
      parameter (IW_CYAN = 6)
      parameter (IW_WHITE = 7)
      parameter (IW_SANDY_BROWN = 8)
      parameter (IW_BROWN = 9)
      parameter (IW_LT_GREY = 10)
      parameter (IW_PURPLE = 11)

	
      integer XK_VoidSymbol
      parameter (XK_VoidSymbol = 'FFFFFF'X)

      integer XK_BackSpace
      parameter (XK_BackSpace =   'FF08'X)
      integer XK_Tab
      parameter (XK_Tab =         'FF09'X)
      integer XK_Linefeed
      parameter (XK_Linefeed =    'FF0A'X)
      integer XK_Clear
      parameter (XK_Clear =       'FF0B'X)
      integer XK_Return
      parameter (XK_Return =      'FF0D'X)
      integer XK_Pause
      parameter (XK_Pause =       'FF13'X)
      integer XK_Scroll_Lock
      parameter (XK_Scroll_Lock = 'FF14'X)
      integer XK_Sys_Req
      parameter (XK_Sys_Req =     'FF15'X)
      integer XK_Escape
      parameter (XK_Escape =      'FF1B'X)
      integer XK_Delete
      parameter (XK_Delete =      'FFFF'X)

      integer XK_Multi_key
      parameter (XK_Multi_key =   'FF20'X)

      integer XK_Kanji
      parameter (XK_Kanji =             'FF21'X)
      integer XK_Muhenkan
      parameter (XK_Muhenkan =          'FF22'X)  
      integer XK_Henkan_Mode
      parameter (XK_Henkan_Mode =       'FF23'X)
      integer XK_Henkan
      parameter (XK_Henkan =            'FF23'X)
      integer XK_Romaji
      parameter (XK_Romaji =            'FF24'X)
      integer XK_Hiragana
      parameter (XK_Hiragana =          'FF25'X)
      integer XK_Katakana
      parameter (XK_Katakana =          'FF26'X)
      integer XK_Hiragana_Katakana
      parameter (XK_Hiragana_Katakana = 'FF27'X) ! Hiragana/Katakana toggle */
      integer XK_Zenkaku
      parameter (XK_Zenkaku =           'FF28'X) ! to Zenkaku */
      integer XK_Hankaku
      parameter (XK_Hankaku =           'FF29'X) ! to Hankaku */
      integer XK_Zenkaku_Hankaku
      parameter (XK_Zenkaku_Hankaku =   'FF2A'X) ! Zenkaku/Hankaku toggle */
      integer XK_Touroku
      parameter (XK_Touroku =           'FF2B'X) ! Add to Dictionary */
      integer XK_Massyo
      parameter (XK_Massyo =            'FF2C'X) ! Delete from Dictionary */
      parameter (XK_Kana_Lock =         'FF2D'X) ! Kana Lock */
      integer XK_Kana_Shift
      parameter (XK_Kana_Shift =        'FF2E'X) ! Kana Shift */
      integer XK_Eisu_Shift
      parameter (XK_Eisu_Shift =        'FF2F'X) ! Alphanumeric Shift */
      integer XK_Eisu_toggle
      parameter (XK_Eisu_toggle =       'FF30'X) ! Alphanumeric toggle */

      integer XK_Home
      parameter (XK_Home =      'FF50'X)
      integer XK_Left
      parameter (XK_Left =      'FF51'X)
      parameter (XK_Up =        'FF52'X) ! Move up, up arrow */
      integer XK_Right
      parameter (XK_Right =     'FF53'X) ! Move right, right arrow */
      integer XK_Down
      parameter (XK_Down =      'FF54'X)
      integer XK_Prior
      parameter (XK_Prior =     'FF55'X) ! Prior, previous */
      integer XK_Page_Up
      parameter (XK_Page_Up =   'FF55'X)
      integer XK_Next
      parameter (XK_Next =      'FF56'X) ! Next */
      integer XK_Page_Down
      parameter (XK_Page_Down = 'FF56'X)
      integer XK_End
      parameter (XK_End =       'FF57'X) ! EOL */
      integer XK_Begin
      parameter (XK_Begin =     'FF58'X) ! BOL */


      integer XK_Select
      parameter (XK_Select =        'FF60'X) ! Select, mark */
      integer XK_Print
      parameter (XK_Print =         'FF61'X)
      integer XK_Execute
      parameter (XK_Execute =       'FF62'X) ! Execute, run, do */
      integer XK_Insert
      parameter (XK_Insert =        'FF63'X) ! Insert, insert here */
      integer XK_Undo
      parameter (XK_Undo =          'FF65'X) ! Undo, oops */
      integer XK_Redo
      parameter (XK_Redo =          'FF66'X) ! redo, again */
      integer XK_Menu
      parameter (XK_Menu =          'FF67'X)
      integer XK_Find
      parameter (XK_Find =          'FF68'X) ! Find, search */
      integer XK_Cancel
      parameter (XK_Cancel =        'FF69'X) ! Cancel, stop, abort, exit */
      integer XK_Help
      parameter (XK_Help =          'FF6A'X) ! Help */
      integer XK_Break
      parameter (XK_Break =         'FF6B'X)
      integer XK_Mode_switch
      parameter (XK_Mode_switch =   'FF7E'X) ! Character set switch */
      integer XK_script_switch
      parameter (XK_script_switch = 'FF7E'X) ! Alias for mode_switch */
      integer XK_Num_Lock
      parameter (XK_Num_Lock =      'FF7F'X)


      integer XK_KP_Space
      parameter (XK_KP_Space =     'FF80'X) ! space */
      integer XK_KP_Tab
      parameter (XK_KP_Tab =       'FF89'X)
      integer XK_KP_Enter
      parameter (XK_KP_Enter =     'FF8D'X) ! enter */
      integer XK_KP_F1
      parameter (XK_KP_F1 =        'FF91'X) ! PF1, KP_A, ... */
      integer XK_KP_F2
      parameter (XK_KP_F2 =        'FF92'X)
      integer XK_KP_F3
      parameter (XK_KP_F3 =        'FF93'X)
      integer XK_KP_F4
      parameter (XK_KP_F4 =        'FF94'X)
      integer XK_KP_Home
      parameter (XK_KP_Home =      'FF95'X)
      integer XK_KP_Left
      parameter (XK_KP_Left =      'FF96'X)
      integer XK_KP_Up
      parameter (XK_KP_Up =        'FF97'X)
      integer XK_KP_Right
      parameter (XK_KP_Right =     'FF98'X)
      integer XK_KP_Down
      parameter (XK_KP_Down =      'FF99'X)
      integer XK_KP_Prior
      parameter (XK_KP_Prior =     'FF9A'X)
      integer XK_KP_Page_Up
      parameter (XK_KP_Page_Up =   'FF9A'X)
      integer XK_KP_Next
      parameter (XK_KP_Next =      'FF9B'X)
      integer XK_KP_Page_Down
      parameter (XK_KP_Page_Down = 'FF9B'X)
      integer XK_KP_End
      parameter (XK_KP_End =       'FF9C'X)
      integer XK_KP_Begin
      parameter (XK_KP_Begin =     'FF9D'X)
      integer XK_KP_Insert
      parameter (XK_KP_Insert =    'FF9E'X)
      integer XK_KP_Delete
      parameter (XK_KP_Delete =    'FF9F'X)
      integer XK_KP_Equal
      parameter (XK_KP_Equal =     'FFBD'X) ! equals */
      integer XK_KP_Multiply
      parameter (XK_KP_Multiply =  'FFAA'X)
      integer XK_KP_Add
      parameter (XK_KP_Add =       'FFAB'X)
      integer XK_KP_Separator
      parameter (XK_KP_Separator = 'FFAC'X) ! separator, often comma */
      integer XK_KP_Subtract
      parameter (XK_KP_Subtract =  'FFAD'X)
      integer XK_KP_Decimal
      parameter (XK_KP_Decimal =   'FFAE'X)
      integer XK_KP_Divide
      parameter (XK_KP_Divide =    'FFAF'X)

      integer XK_KP_0
      parameter (XK_KP_0 = 'FFB0'X)
      integer XK_KP_1
      parameter (XK_KP_1 = 'FFB1'X)
      integer XK_KP_2
      parameter (XK_KP_2 = 'FFB2'X)
      integer XK_KP_3
      parameter (XK_KP_3 = 'FFB3'X)
      integer XK_KP_4
      parameter (XK_KP_4 = 'FFB4'X)
      integer XK_KP_5
      parameter (XK_KP_5 = 'FFB5'X)
      integer XK_KP_6
      parameter (XK_KP_6 = 'FFB6'X)
      integer XK_KP_7
      parameter (XK_KP_7 = 'FFB7'X)
      integer XK_KP_8
      parameter (XK_KP_8 = 'FFB8'X)
      integer XK_KP_9
      parameter (XK_KP_9 = 'FFB9'X)

      integer XK_F1
      parameter (XK_F1 =  'FFBE'X)
      integer XK_F2
      parameter (XK_F2 =  'FFBF'X)
      integer XK_F3
      parameter (XK_F3 =  'FFC0'X)
      integer XK_F4
      parameter (XK_F4 =  'FFC1'X)
      integer XK_F5
      parameter (XK_F5 =  'FFC2'X)
      integer XK_F6
      parameter (XK_F6 =  'FFC3'X)
      integer XK_F7
      parameter (XK_F7 =  'FFC4'X)
      integer XK_F8
      parameter (XK_F8 =  'FFC5'X)
      integer XK_F9
      parameter (XK_F9 =  'FFC6'X)
      integer XK_F10
      parameter (XK_F10 = 'FFC7'X)
      integer XK_F11
      parameter (XK_F11 = 'FFC8'X)
      integer XK_L1
      parameter (XK_L1 =  'FFC8'X)
      integer XK_F12
      parameter (XK_F12 = 'FFC9'X)
      integer XK_L2
      parameter (XK_L2 =  'FFC9'X)
      integer XK_F13
      parameter (XK_F13 = 'FFCA'X)
      integer XK_L3
      parameter (XK_L3 =  'FFCA'X)
      integer XK_F14
      parameter (XK_F14 = 'FFCB'X)
      integer XK_L4
      parameter (XK_L4 =  'FFCB'X)
      integer XK_F15
      parameter (XK_F15 = 'FFCC'X)
      integer XK_L5
      parameter (XK_L5 =  'FFCC'X)
      integer XK_F16
      parameter (XK_F16 = 'FFCD'X)
      integer XK_L6
      parameter (XK_L6 =  'FFCD'X)
      integer XK_F17
      parameter (XK_F17 = 'FFCE'X)
      integer XK_L7
      parameter (XK_L7 =  'FFCE'X)
      integer XK_F18
      parameter (XK_F18 = 'FFCF'X)
      integer XK_L8
      parameter (XK_L8 =  'FFCF'X)
      integer XK_F19
      parameter (XK_F19 = 'FFD0'X)
      integer XK_L9
      parameter (XK_L9 =  'FFD0'X)
      integer XK_F20
      parameter (XK_F20 = 'FFD1'X)
      integer XK_L10
      parameter (XK_L10 = 'FFD1'X)
      integer XK_F21
      parameter (XK_F21 = 'FFD2'X)
      integer XK_R1
      parameter (XK_R1 =  'FFD2'X)
      integer XK_F22
      parameter (XK_F22 = 'FFD3'X)
      integer XK_R2
      parameter (XK_R2 =  'FFD3'X)
      integer XK_F23
      parameter (XK_F23 = 'FFD4'X)
      integer XK_R3
      parameter (XK_R3 =  'FFD4'X)
      integer XK_F24
      parameter (XK_F24 = 'FFD5'X)
      integer XK_R4
      parameter (XK_R4 =  'FFD5'X)
      integer XK_F25
      parameter (XK_F25 = 'FFD6'X)
      integer XK_R5
      parameter (XK_R5 =  'FFD6'X)
      integer XK_F26
      parameter (XK_F26 = 'FFD7'X)
      integer XK_R6
      parameter (XK_R6 =  'FFD7'X)
      integer XK_F27
      parameter (XK_F27 = 'FFD8'X)
      integer XK_R7
      parameter (XK_R7 =  'FFD8'X)
      integer XK_F28
      parameter (XK_F28 = 'FFD9'X)
      integer XK_R8
      parameter (XK_R8 =  'FFD9'X)
      integer XK_F29
      parameter (XK_F29 = 'FFDA'X)
      integer XK_R9
      parameter (XK_R9 =  'FFDA'X)
      integer XK_F30
      parameter (XK_F30 = 'FFDB'X)
      integer XK_R10
      parameter (XK_R10 = 'FFDB'X)
      integer XK_F31
      parameter (XK_F31 = 'FFDC'X)
      integer XK_R11
      parameter (XK_R11 = 'FFDC'X)
      integer XK_F32
      parameter (XK_F32 = 'FFDD'X)
      integer XK_R12
      parameter (XK_R12 = 'FFDD'X)
      integer XK_F33
      parameter (XK_F33 = 'FFDE'X)
      integer XK_R13
      parameter (XK_R13 = 'FFDE'X)
      integer XK_F34
      parameter (XK_F34 = 'FFDF'X)
      integer XK_R14
      parameter (XK_R14 = 'FFDF'X)
      integer XK_F35
      parameter (XK_F35 = 'FFE0'X)
      integer XK_R15
      parameter (XK_R15 = 'FFE0'X)

      integer XK_Shift_L      
      parameter (XK_Shift_L       = 'FFE1'X) ! Left  */
      integer XK_Shift_R      
      parameter (XK_Shift_R       = 'FFE2'X) ! Right  */
      integer XK_Control_L    
      parameter (XK_Control_L     = 'FFE3'X) ! Left control */
      integer XK_Control_R    
      parameter (XK_Control_R     = 'FFE4'X) ! Right control */
      integer XK_Caps_Lock    
      parameter (XK_Caps_Lock     = 'FFE5'X) ! Caps lock */
      integer XK_Shift_Lock   
      parameter (XK_Shift_Lock    = 'FFE6'X) ! Shift lock */

      integer XK_Meta_L       
      parameter (XK_Meta_L        = 'FFE7'X) ! Left meta */
      integer XK_Meta_R       
      parameter (XK_Meta_R        = 'FFE8'X) ! Right meta */
      integer XK_Alt_L        
      parameter (XK_Alt_L         = 'FFE9'X) ! Left alt */
      integer XK_Alt_R        
      parameter (XK_Alt_R         = 'FFEA'X) ! Right alt */
      integer XK_Super_L      
      parameter (XK_Super_L       = 'FFEB'X) ! Left super */
      integer XK_Super_R      
      parameter (XK_Super_R       = 'FFEC'X) ! Right super */
      integer XK_Hyper_L      
      parameter (XK_Hyper_L       = 'FFED'X) ! Left hyper */
      integer XK_Hyper_R      
      parameter (XK_Hyper_R       = 'FFEE'X) ! Right hyper */

      integer XK_space        
      parameter (XK_space         = '020'X)
      integer XK_exclam       
      parameter (XK_exclam        = '021'X)
      integer XK_quotedbl     
      parameter (XK_quotedbl      = '022'X)
      integer XK_numbersign   
      parameter (XK_numbersign    = '023'X)
      integer XK_dollar       
      parameter (XK_dollar        = '024'X)
      integer XK_percent      
      parameter (XK_percent       = '025'X)
      integer XK_ampersand    
      parameter (XK_ampersand     = '026'X)
      integer XK_apostrophe   
      parameter (XK_apostrophe    = '027'X)
      integer XK_quoteright   
      parameter (XK_quoteright    = '027'X) ! deprecated */
      integer XK_parenleft    
      parameter (XK_parenleft     = '028'X)
      integer XK_parenright   
      parameter (XK_parenright    = '029'X)
      integer XK_asterisk     
      parameter (XK_asterisk      = '02a'X)
      integer XK_plus         
      parameter (XK_plus          = '02b'X)
      integer XK_comma        
      parameter (XK_comma         = '02c'X)
      integer XK_minus        
      parameter (XK_minus         = '02d'X)
      integer XK_period       
      parameter (XK_period        = '02e'X)
      integer XK_slash        
      parameter (XK_slash         = '02f'X)
      integer XK_0            
      parameter (XK_0             = '030'X)
      integer XK_1            
      parameter (XK_1             = '031'X)
      integer XK_2            
      parameter (XK_2             = '032'X)
      integer XK_3            
      parameter (XK_3             = '033'X)
      integer XK_4            
      parameter (XK_4             = '034'X)
      integer XK_5            
      parameter (XK_5             = '035'X)
      integer XK_6            
      parameter (XK_6             = '036'X)
      integer XK_7            
      parameter (XK_7             = '037'X)
      integer XK_8            
      parameter (XK_8             = '038'X)
      integer XK_9            
      parameter (XK_9             = '039'X)
      integer XK_colon        
      parameter (XK_colon         = '03a'X)
      integer XK_semicolon    
      parameter (XK_semicolon     = '03b'X)
      integer XK_less         
      parameter (XK_less          = '03c'X)
      integer XK_equal        
      parameter (XK_equal         = '03d'X)
      integer XK_greater      
      parameter (XK_greater       = '03e'X)
      integer XK_question     
      parameter (XK_question      = '03f'X)
      integer XK_at           
      parameter (XK_at            = '040'X)
      integer XK_A            
      parameter (XK_A             = '041'X)
      integer XK_B            
      parameter (XK_B             = '042'X)
      integer XK_C            
      parameter (XK_C             = '043'X)
      integer XK_D            
      parameter (XK_D             = '044'X)
      integer XK_E            
      parameter (XK_E             = '045'X)
      integer XK_F            
      parameter (XK_F             = '046'X)
      integer XK_G            
      parameter (XK_G             = '047'X)
      integer XK_H            
      parameter (XK_H             = '048'X)
      integer XK_I            
      parameter (XK_I             = '049'X)
      integer XK_J            
      parameter (XK_J             = '04a'X)
      integer XK_K            
      parameter (XK_K             = '04b'X)
      integer XK_L            
      parameter (XK_L             = '04c'X)
      integer XK_M            
      parameter (XK_M             = '04d'X)
      integer XK_N            
      parameter (XK_N             = '04e'X)
      integer XK_O            
      parameter (XK_O             = '04f'X)
      integer XK_P            
      parameter (XK_P             = '050'X)
      integer XK_Q            
      parameter (XK_Q             = '051'X)
      integer XK_R            
      parameter (XK_R             = '052'X)
      integer XK_S            
      parameter (XK_S             = '053'X)
      integer XK_T            
      parameter (XK_T             = '054'X)
      integer XK_U            
      parameter (XK_U             = '055'X)
      integer XK_V            
      parameter (XK_V             = '056'X)
      integer XK_W            
      parameter (XK_W             = '057'X)
      integer XK_X            
      parameter (XK_X             = '058'X)
      integer XK_Y            
      parameter (XK_Y             = '059'X)
      integer XK_Z            
      parameter (XK_Z             = '05a'X)
      integer XK_bracketleft  
      parameter (XK_bracketleft   = '05b'X)
      integer XK_backslash    
      parameter (XK_backslash     = '05c'X)
      integer XK_bracketright 
      parameter (XK_bracketright  = '05d'X)
      integer XK_asciicircum  
      parameter (XK_asciicircum   = '05e'X)
      integer XK_underscore   
      parameter (XK_underscore    = '05f'X)
      integer XK_grave        
      parameter (XK_grave         = '060'X)
      integer XK_quoteleft    
      parameter (XK_quoteleft     = '060'X)

      integer XK_la           
      parameter (XK_la            = '061'X)
      integer XK_lb           
      parameter (XK_lb            = '062'X)
      integer XK_lc           
      parameter (XK_lc            = '063'X)
      integer XK_ld           
      parameter (XK_ld            = '064'X)
      integer XK_le           
      parameter (XK_le            = '065'X)
      integer XK_lf           
      parameter (XK_lf            = '066'X)
      integer XK_lg           
      parameter (XK_lg            = '067'X)
      integer XK_lh           
      parameter (XK_lh            = '068'X)
      integer XK_li           
      parameter (XK_li            = '069'X)
      integer XK_lj           
      parameter (XK_lj            = '06a'X)
      integer XK_lk           
      parameter (XK_lk            = '06b'X)
      integer XK_ll           
      parameter (XK_ll            = '06c'X)
      integer XK_lm           
      parameter (XK_lm            = '06d'X)
      integer XK_ln           
      parameter (XK_ln            = '06e'X)
      integer XK_lo           
      parameter (XK_lo            = '06f'X)
      integer XK_lp           
      parameter (XK_lp            = '070'X)
      integer XK_lq           
      parameter (XK_lq            = '071'X)
      integer XK_lr           
      parameter (XK_lr            = '072'X)
      integer XK_ls           
      parameter (XK_ls            = '073'X)
      integer XK_lt           
      parameter (XK_lt            = '074'X)
      integer XK_lu           
      parameter (XK_lu            = '075'X)
      integer XK_lv           
      parameter (XK_lv            = '076'X)
      integer XK_lw           
      parameter (XK_lw            = '077'X)
      integer XK_lx           
      parameter (XK_lx            = '078'X)
      integer XK_ly           
      parameter (XK_ly            = '079'X)
      integer XK_lz           
      parameter (XK_lz            = '07a'X)
      integer XK_braceleft    
      parameter (XK_braceleft     = '07b'X)
      integer XK_bar          
      parameter (XK_bar           = '07c'X)
      integer XK_braceright   
      parameter (XK_braceright    = '07d'X)
      integer XK_asciitilde   
      parameter (XK_asciitilde    = '07e'X)
