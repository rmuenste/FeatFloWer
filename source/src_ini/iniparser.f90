module iniparser

  implicit none

  !private

  !<constants>

  integer, parameter, public :: INIP_STRLEN = 256
  ! flag for appending data to a file
  integer, parameter, public :: INIP_APPEND = 0

  ! flag for replacing a file
  integer, parameter, public :: INIP_REPLACE = 1

  character(len=1), parameter, public :: NEWLINE = achar(10)

  ! logical value 'true'
  integer, parameter, public :: YES = 0

  ! logical value 'false'
  integer, parameter, public :: NO = 1


  !<constantblock>

  ! Maximum length of a section name.
  integer, parameter, public :: INIP_MLSECTION = 64

  ! Maximum length of parameter names: 32 characters
  integer, parameter, public :: INIP_MLNAME = 32

  ! Default length of parameter data: 256 characters
  integer, parameter, public :: INIP_MLDATA = 256

  ! Minimum number of free parameter 'slots' per parameter section.
  ! If there are too many parameters in a parameter section, the
  ! structure is dynamically extended in terms of INIP_NPARSPERBLOCK
  ! entries.
  integer, parameter, public :: INIP_NPARSPERBLOCK = 32

  ! Minimum number of parameter sections.
  ! If there are too many parameter sections in a parameter block, the
  ! structure is dynamically extended in terms of INIP_NSECTIONS
  ! entries.
  integer, parameter, public :: INIP_NSECTIONS = 8

  ! Maximum length of a line in a INI file. Lines longer than this
  ! are truncated.
  integer, parameter, public :: INIP_LENLINEBUF = 1024

  ! Comment character
  character, parameter, public :: INIP_COMMENT = "#"

  !</constantblock>

  !</constants>

  !<types>

  !<typeblock>

  ! This structure realises a value associated to a parameter name.
  ! A value consists of one string or an array of strings.

  type t_parlstValue

    private

    ! Number of strings. If set to 0, the value consists of one
    ! string, to be found in svalue. If > 0, there are nsize
    ! strings to be found in p_SentryList.
    integer :: nsize = 0

    ! Single string; contains the value in case nsize=0
    character, dimension(:), pointer :: p_sentry => null()

    ! Array of strings in case nsize>0
    character, dimension(:,:), pointer :: p_SentryList => null()

  end type

  public :: t_parlstValue

  !</typeblock>

  !<typeblock>

  ! This structure realises a parameter section. It contains an
  ! array with parameter names and an array with parameter values
  ! to these names. The arrays are dynamically allocated.

  type t_parlstSection

    ! The name of the section.
    character(LEN=INIP_MLSECTION) :: ssectionName = ''

    ! Actual number of parameters in this section.
    integer :: iparamCount = 0

    ! A list of parameter names. Each name contains INIP_MLNAME
    ! characters.
    character(LEN=INIP_MLNAME), dimension(:), pointer :: p_Sparameters => null()

    ! A list of t_parlstValue structures corresponding to the parameters
    ! in p_Sparameters.
    type(t_parlstValue), dimension(:), pointer :: p_Rvalues

  end type

  public :: t_parlstSection

  !</typeblock>

  !<typeblock>

  ! This structure realises a parameter list. Parameters can be read into
  ! it from a file. Parameters can be obtained from the structure using
  ! the query/get routines.

  type t_parlist

    private

    ! Actual number of sections in the parameter list. There is at least
    ! one section - the unnamed section. If this value is =0, the parameter
    ! list is not initialised.
    integer :: isectionCount = 0

    ! A list of sections. The first section is always the unnamed section.
    type(t_parlstSection), dimension(:), pointer :: p_Rsections => null()

  end type

  public :: t_parlist

  !</typeblock>

  !</types>

  private :: inip_initsection, INIP_reallocsection, INIP_realloclist
  private :: INIP_reallocSubVariables
  private :: INIP_fetchparameter,INIP_readlinefromfile,INIP_parseline

  interface INIP_queryvalue
    module procedure INIP_queryvalue_direct
    module procedure INIP_queryvalue_indir
  end interface

  interface INIP_querysubstrings
    module procedure INIP_querysubstrings_direct
    module procedure INIP_querysubstrings_indir
  end interface

  interface INIP_addvalue
    module procedure INIP_addvalue_direct
    module procedure INIP_addvalue_indir
  end interface

  interface INIP_setvalue
    module procedure INIP_setvalue_fetch
    module procedure INIP_setvalue_indir
    module procedure INIP_setvalue_direct
  end interface

  interface INIP_getvalue_string
    module procedure INIP_getvalue_string_fetch
    module procedure INIP_getvalue_string_indir
    module procedure INIP_getvalue_string_direct
  end interface

  interface INIP_getvalue_int
    module procedure INIP_getvalue_int_fetch
    module procedure INIP_getvalue_int_direct
    module procedure INIP_getvalue_int_indir
  end interface

  interface INIP_getvalue_single
    module procedure INIP_getvalue_single_fetch
    module procedure INIP_getvalue_single_indir
    module procedure INIP_getvalue_single_direct
  end interface

  interface INIP_getvalue_double
    module procedure INIP_getvalue_double_fetch
    module procedure INIP_getvalue_double_indir
    module procedure INIP_getvalue_double_direct
  end interface INIP_getvalue_double

  interface INIP_getvalue_logical
    module procedure INIP_getvalue_logical_fetch
    module procedure INIP_getvalue_logical_indir
    module procedure INIP_getvalue_logical_direct
  end interface

  interface INIP_findvalue
    module procedure INIP_findvalue_direct
    module procedure INIP_findvalue_indir
  end interface

  integer, parameter :: out_linelength = 80
  integer, public    :: inip_iunit_protfile = -1
  integer, public    :: inip_iunit_terminal = -1
  integer, public    :: inip_showid = 1
  integer, public    :: inip_id = -1

  private

  public :: inip_output_init
  public :: inip_getFreeUnit
  public :: INIP_readfromfile
  public :: INIP_readfromsinglefile
  public :: INIP_info
  public :: INIP_dumpToFile
  public :: INIP_init
  public :: INIP_done
  public :: INIP_clear
  public :: INIP_getvalue_double
  public :: INIP_getvalue_single
  public :: INIP_getvalue_string
  public :: INIP_getvalue_int
  public :: INIP_getvalue_logical
  public :: INIP_querysection
  public :: INIP_addsection
  public :: INIP_addvalue
  public :: INIP_setvalue
  public :: INIP_findvalue
  public :: INIP_getStringRepresentation
  public :: INIP_querysubstrings
  public :: inip_toupper_replace

  interface inip_toupper
    module procedure inip_toupper_replace
    module procedure inip_toupper_copy
  end interface

contains

  ! Some auxiliary routines

  !<function>
  function inip_isDirectory (spath) result (bisDir)

    !<description>
    ! Checks whether a given string is a directory
    !</description>

    !<input>
    ! Path to the file.
    character(len=*), intent(in) :: spath
    !</input>

    !<result>
    ! Is .TRUE. if the given string is an existing directory
    logical :: bisDir
    !</result>

    !</function>

    integer :: ierr
    integer :: iisdir

    ! use wrapper for C system call stat() in isdirectory.c
    external isdirectory
    ierr = 0
    bisDir = .false.
    call isdirectory(trim(spath) // achar(0), iisdir, ierr)
    if (ierr .ge. 0) then
      bisDir = (iisdir .gt. 0)
    endif

  end function inip_isDirectory

  !<subroutine>
  subroutine inip_makeDirectory (spath)

    !<description>
    ! Portably creates directories (even recursively)
    !</description>

    !<input>
    ! Path to the file.
    character(len=*), intent(in) :: spath
    !</input>

    !</subroutine>

    external mkdir_recursive

    ! status variable for (recursive) directory creation
    ! > 0: number of subdirectories created
    ! < 0: -errno
    integer :: istatus


    if (len_trim(spath) .eq. 0) then
      call inip_output_line('Path name <'//trim(adjustl(spath))//'> empty.')
      return
    end if

    ! Ensure to explicitly NULL-terminate the string when calling C function!
    call mkdir_recursive(trim(spath) // achar(0), istatus)
    if (istatus .lt. 0) then
      call inip_output_line('Could not (auto-)create path <'//trim(adjustl(spath))//'>.')
      call inip_sys_halt()
    end if

  end subroutine inip_makeDirectory

  !<subroutine>
  subroutine inip_openFileForReading(sfilename, iunit, bformatted)

    !<description>
    !This routine tries to open a file for reading. If succesful, on can read from it
    !via unit "iunit". Otherwise, iunit is -1.
    !</description>

    !<input>

    !filename
    character(*), intent(in) :: sfilename

    ! OPTIONAL:
    ! TRUE : Open the file formatted, i.e. in human readable form
    ! FALSE: Open the file in unformatted, machine dependent form
    ! If not specified, the default system dependent setting is used.
    logical, intent(in), optional :: bformatted

    !</input>

    !<output>

    !number of unit
    integer, intent(out) :: iunit
    !</output>
    !</subroutine>

    logical :: bexists !true, if a file with name sfilename exists
    integer :: istatus !status variable for opening. 0, if opening succesful


    if (trim(sfilename) .eq. "") then
      call inip_output_line('File name <'//trim(adjustl(sfilename))//'> empty.')
      return
    endif

    iunit = inip_getFreeUnit()
    if (iunit .eq. -1) then
      call inip_output_line('No free unit found. Not able to open the file.')
      !give it up
      return
    endif

    inquire(file=trim(sfilename), exist=bexists)

    if (bexists) then

      if (.not. present(bformatted)) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read")
      else if (bformatted) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
          form="formatted")
      else
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="read",&
          form="unformatted")
      end if
      if (istatus .ne. 0) then
        call inip_output_line('Error while opening file <'//trim(adjustl(sfilename))//'> for reading.')
        iunit = -1
      end if

    else

      call inip_output_line('File <'//trim(adjustl(sfilename))//'> does not exist.')
      call inip_sys_halt()
    endif

  end subroutine inip_openFileForReading

  !<subroutine>
  subroutine inip_openFileForWriting(sfilename, iunit, cflag, bfileExists, bformatted)

    !<description>
    !This routine tries to open a file for writing. If succesful, one can write to it
    !via unit "iunit". Otherwise, iunit is -1. cflag specifies if an already existing file
    !should be replaced or if the output should be appended. bfileExists is an optional
    !parameter, which will be set to true, if the file already existed, otherwise false.
    !</description>

    !<input>

    !filename
    character(*), intent(in) :: sfilename

    !mode: SYS_APPEND or SYS_REPLACE
    integer, intent(in) :: cflag

    ! OPTIONAL:
    ! TRUE : Open the file formatted, i.e. in human readable form
    ! FALSE: Open the file in unformatted, machine dependent form
    ! If not specified, the default system dependent setting is used.
    logical, intent(in), optional :: bformatted
    !</input>

    !<output>

    ! unit of the opened file
    integer, intent(out) :: iunit

    ! optional parameter (see description)
    logical, intent(out), optional :: bfileExists

    !</output>
    !</subroutine>

    logical :: bexists !true, if the file to be written in exists
    integer :: istatus !status variable for opening procedure

    ! the result "dirname(sfilename)" would yield
    character(len=len(sfilename)) :: sfilepath

    if (len_trim(sfilename) .eq. 0) then
      call inip_output_line('File name <'//trim(adjustl(sfilename))//'> empty.')
      return
    end if

    iunit = inip_getFreeUnit()
    if (iunit .eq. -1) then
      call inip_output_line('No free unit found. Not able to open the file.')
      return
    end if

    inquire(file=trim(sfilename), exist=bexists)
    if (.not. present(bformatted)) then
      if (bexists .and. cflag .eq. INIP_REPLACE) then
        open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
          action="write")
      else
        ! Ensure that the path up to the given file name does exist
        call inip_pathExtract(sfilename, sfilepath)
        if (.not. inip_isDirectory(sfilepath)) then
          call inip_makeDirectory(sfilepath)
        end if
        open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
          position="append")
      end if
    else
      if (bexists .and. cflag .eq. INIP_REPLACE) then
        if (bformatted) then
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
            action="write", form="formatted")
        else
          open(unit=iunit, file=trim(sfilename), iostat=istatus, status="replace", &
            action="write", form="unformatted")
        end if
      else
        ! Ensure that the path up to the given file name does exist. If it does not,
        ! create it
        call inip_pathExtract(sfilename, sfilepath)
        if (.not. inip_isDirectory(sfilepath)) then
          call inip_makeDirectory(sfilepath)
        end if
        if (bformatted) then
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
            position="append", form="formatted")
        else
          open(unit=iunit, file=trim(sfilename), iostat=istatus, action="write", &
            position="append", form="unformatted")
        end if
      end if
    end if
    if (present(bfileExists)) then
      bfileExists = bexists
    end if
    if (istatus .ne. 0) then
      call inip_output_line('Error while opening file <'//trim(adjustl(sfilename))//'> for writing.')
      iunit = -1
    end if

  end subroutine inip_openFileForWriting

  ! ***************************************************************************

  integer function inip_getFreeUnit()

    !<description>
    !This routine tries to find a free unit (for file input/output). If a free unit is
    !found, it is returned, otherwise -1 is returned.
    !</description>

    !<result>
    !number of free unit (-1 if no free unit available)
    !</result>

    !</function>

    logical :: bexists, bopened!flags indicating errors
    integer :: itry !free unit candidate

    inip_getFreeUnit = -1
    do itry = 20,10000
      !does unit exist?
      inquire(unit=itry, exist=bexists)
      if (bexists) then
        !is unit already opened?
        inquire(unit=itry, opened=bopened)
        if (.not. bopened) then
          !free unit found
          inip_getFreeUnit = itry
          !exit do-loop
          exit
        endif
      endif
    enddo
    if (inip_getFreeUnit .eq. -1) then
      call inip_output_line("*** WARNING! No free unit between 20 and 10000 found! ***")
    endif

  end function inip_getFreeUnit


  ! init the output by setting the ID, showid etc.
  subroutine inip_output_init(myid,showid,unit_protfile,unit_terminal)
    integer, intent(in) :: myid, showid, unit_protfile, unit_terminal
    inip_iunit_protfile = unit_protfile
    inip_iunit_protfile = unit_terminal
    inip_showid = showid
    inip_id = myid
  end subroutine


  subroutine inip_output_line(smessage)
    character(LEN=*), intent(in) :: smessage

    if (inip_id .eq. inip_showid) then
      if (inip_iunit_terminal .ne. -1) then
        write(inip_iunit_terminal,'(A)') trim(adjustl(smessage))
      end if

      if(inip_iunit_protfile .ne. -1) then
        write(inip_iunit_protfile,'(A)') trim(adjustl(smessage))
      end if
    end if

  end subroutine

  subroutine inip_output_lbrk()
    call inip_output_line('')
  end subroutine

  subroutine inip_sys_halt()
    stop
  end subroutine inip_sys_halt

  !<subroutine>

  subroutine inip_toupper_replace (str)

    !<description>
    ! Convert a string to upper case.
    ! The given string is replaced by its uppercase version.
    !</description>

    !<inputoutput>

    ! The string that is to make uppercase
    character(LEN=*), intent(inout) :: str

    !</inputoutput>

    !</subroutine>

    ! local variables
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c

    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "a") .and. (c .le. "z")) then
        str(i:i) = achar (iachar(c) - up2low)
      end if
    end do

  end subroutine inip_toupper_replace

  !************************************************************************

  !<subroutine>

  subroutine inip_toupper_copy (str,strUpper)

    !<description>
    ! Convert a string to upper case.
    !</description>

    !<input>

    ! The string that is to make uppercase
    character(LEN=*), intent(in) :: str

    !</input>

    !<output>

    ! Uppercase version of the given string
    character(LEN=*), intent(out) :: strUpper

    !</output>

    !</subroutine>

    ! local variables
    integer, parameter :: up2low = iachar("a") - iachar("A")
    integer :: i
    character    :: c

    if (len(str) .gt. len(strUpper)) then
      call inip_output_line("inip_toupper_copy: target string is too short")
      call inip_sys_halt()
    end if

    ! Initialise string
    strUpper = ""

    do i=1,len(str)
      c = str(i:i)
      if ((c .ge. "a") .and. (c .le. "z")) then
        strUpper(i:i) = achar (iachar(c) - up2low)
      else
        strUpper(i:i) = c
      end if
    end do

  end subroutine inip_toupper_copy


  !<subroutine>

  subroutine inip_stringToCharArray (sstring,schararray,slength)

    !<description>
    ! Converts a string to a character array.
    !</description>

    !<input>
    ! String to convert
    character(len=*), intent(in) :: sstring

    ! OPTIONAL: Length of the string.
    ! If not present, the default string length is used.
    integer, intent(in), optional :: slength
    !</input>

    !<output>
    ! Character array that receives the converted string.
    ! Must be at least as long as the string or as slength.
    character, dimension(:), intent(out) :: schararray
    !</output>

    !</subroutine>

    integer :: i,j

    if (present(slength)) then
      j = slength
    else
      j = min(len_trim(sstring),size(schararray))
    end if

    ! Copy all characters.
    do i=1,j
      schararray(i) = sstring(i:i)
    end do

    ! Fill up the rest with spaces. This emulates a string copy.
    schararray(j+1:) = " "

  end subroutine

  !************************************************************************

  !<subroutine>

  subroutine inip_charArrayToString (schararray,sstring,slength)

    !<description>
    ! Converts a character array to a string.
    !</description>

    !<input>
    ! Character array to convert
    character, dimension(:), intent(in) :: schararray

    ! OPTIONAL: Length of the string.
    ! If not present, the default string length is used.
    integer, intent(in), optional :: slength
    !</input>

    !<output>
    ! Character array that receives the converted string.
    ! Must be at least as long as the character array or slength.
    character(len=*), intent(out) :: sstring
    !</output>

    !</subroutine>

    integer :: i,j

    if (present(slength)) then
      j = slength
    else
      j = min(size(schararray),len(sstring))
    end if

    ! Copy all characters.
    sstring = ""
    do i=1,j
      sstring(i:i) = schararray(i)
    end do

  end subroutine


  subroutine inip_dequote (sstring)

    !<description>
    ! Removes possible quotation marks around a string.
    !</description>

    !<inputoutput>
    ! String to de-quote
    character(len=*), intent(inout) :: sstring
    !</inputoutput>

    !</subroutine>

    character(len=len(sstring)+1) :: sstring2

    ! Adjust the string
    sstring2=trim(adjustl(sstring))

    ! Does the string start with a quotation mark?
    if ((sstring2(1:1) .eq. "'") .or. &
        (sstring2(1:1) .eq. """")) then
      ! Re-read the string, remove them.
      read (sstring2,*) sstring
    else
      ! Just transfer the string, it is ok.
      sstring = sstring2
    end if

  end subroutine

  !<subroutine>

  subroutine inip_pathExtract (sfile, sfilepath, sfilename, babsolute)

    !<description>
    ! Extracts the path of a file from a path+filename string.
    !</description>

    !<input>
    ! Filename + path of a specific file (or directory).
    character(len=*), intent(in) :: sfile
    !</input>

    !<output>
    ! OPTIONAL: Receives the directory that contains the specific file,
    ! or "" if no directory was specified in sfile.
    character(len=*), intent(out), optional :: sfilepath

    ! OPTIONAL: Receives the name of the file without a probably preceding
    ! directory string.
    character(len=*), intent(out), optional :: sfilename

    ! OPTIONAL: Returns TRUE if the path specification in sfile points to an
    ! absolute path. Returns FALSE if the path in sfile is relative.
    logical, intent(out), optional :: babsolute
    !</output>

    !</subroutine>

    integer :: i
    character(len=10) :: ssubpath

    ! Find the last "/" or "\" in sfile.                                (!" cpp fix)
    ! Note that we specified "\\" and not "\" because the PGI compiler  (!" cpp fix)
    ! (stupid thing) would otherwise use the backslash to escape the quote
    ! character. So PGI sees "/\" and other compiler see "/\\", but this (!" cpp fix)
    ! does not matter since the string must only contain a couple of
    ! delimiters which may occur more than once in the string.
    i = scan(sfile,"/\\",.true.)
    if (i .ne. 0) then
      ! Directory ends at position i.
      if (present(sfilepath)) sfilepath = sfile(1:i-1)
      if (present(sfilename)) sfilename = sfile(i+1:)
    else
      ! No directory specified.
      if (present(sfilepath)) sfilepath = ""
      if (present(sfilename)) sfilename = sfile
    end if

    if (present(babsolute)) then
      ! Take a look if this is an absolute or relative path.
      i = scan(trim(adjustl(sfile)),"/\\",.false.)
      babsolute = i .eq. 1

      ! In Windows environments, the path is also absolute if
      ! a volume descriptor like "C:" precedes the (back-)slash.
      if (.not. babsolute) then
        if (i .eq. 3) then
          ! Extract the first 10 characters and check
          ssubpath = trim(adjustl(sfile))
          if (ssubpath(2:2) .eq. ":") then
            babsolute = .true.
          end if
        end if
      end if
    end if

  end subroutine


  ! ***************************************************************************

  ! Internal subroutine: Initialise a newly created parameter section.

  subroutine INIP_initsection (rparlstSection,sname)

    type(t_parlstSection), intent(inout) :: rparlstSection
    character(LEN=*), intent(in) :: sname

    ! Simply allocate the pointers with an empty list
    allocate(rparlstSection%p_Sparameters(INIP_NPARSPERBLOCK))
    allocate(rparlstSection%p_Rvalues(INIP_NPARSPERBLOCK))

    ! and set the section name
    rparlstSection%ssectionName = sname

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Reallocate a section.
  ! This increases the size of a parameter section by reallocation of the
  ! arrays.

  subroutine INIP_reallocsection (rparlstSection, inewsize)

    ! The section to reallocate.
    type(t_parlstSection), intent(inout) :: rparlstSection

    ! The new 'size' of the section, i.e. the new number of parameters,
    ! the section should be able to handle.
    integer, intent(in) :: inewsize

    ! local variables

    integer :: sz,oldsize

    ! Pointers to new lists for replacing the old.
    character(LEN=INIP_MLNAME), dimension(:), pointer :: p_Sparameters
    type(t_parlstValue), dimension(:), pointer :: p_Rvalues

    oldsize = size(rparlstSection%p_Sparameters)
    sz = max(oldsize,inewsize)

    if (size(rparlstSection%p_Sparameters) .eq. sz) return ! nothing to do

    ! Allocate the pointers for the new lists
    allocate(p_Sparameters(sz))
    allocate(p_Rvalues(sz))

    ! Copy the content of the old ones
    p_Sparameters(1:oldsize) = rparlstSection%p_Sparameters (1:oldsize)
    p_Rvalues(1:oldsize) = rparlstSection%p_Rvalues (1:oldsize)

    ! Throw away the old arrays, replace by the new ones
    deallocate(rparlstSection%p_Rvalues)
    deallocate(rparlstSection%p_Sparameters)

    rparlstSection%p_Sparameters => p_Sparameters
    rparlstSection%p_Rvalues => p_Rvalues

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Reallocates a sub-parameter list
  ! This increases the size of the sub-parameters of a parameter.
  ! Old strings are copied.

  subroutine INIP_reallocSubVariables (rvalue, inewsize)

    ! THe parameter item where to reallocate sub-parameters.
    type(t_parlstValue), intent(inout) :: rvalue

    ! The new 'length' of the sub-parameters.
    integer, intent(in) :: inewsize

    ! local variables
    integer :: i,iold
    character, dimension(:,:), pointer :: p_Sdata

    if (ubound(rvalue%p_SentryList,2) .eq. 0) return ! nothing to do

    ! Old size
    iold = ubound(rvalue%p_SentryList,1)

    ! Allocate memory for the strings.
    allocate(p_Sdata(inewsize,ubound(rvalue%p_SentryList,2)))
    p_Sdata(:,:) = ' '

    ! Copy the substrings.
    do i=1,ubound(rvalue%p_SentryList,2)
      p_Sdata(1:min(inewsize,iold),i) = rvalue%p_SentryList(1:min(inewsize,iold),i)
    end do

    ! Replace the data.
    deallocate(rvalue%p_SentryList)
    rvalue%p_SentryList => p_Sdata

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Release a section.
  ! Removes all temporary memory that is allocated by a section.

  subroutine INIP_releasesection (rparlstSection)

    ! The section to release.
    type(t_parlstSection), intent(inout) :: rparlstSection

    ! local variables
    integer :: i

    ! Loop through all values in the current section if there is
    ! an array-value. Release them.
    do i=size(rparlstSection%p_Rvalues),1,-1
      if (associated(rparlstSection%p_Rvalues(i)%p_sentry)) then
        deallocate(rparlstSection%p_Rvalues(i)%p_sentry)
      end if
      if (rparlstSection%p_Rvalues(i)%nsize .gt. 0) then
        deallocate(rparlstSection%p_Rvalues(i)%p_SentryList)
      end if
    end do

    ! Remove the content of the section.
    deallocate(rparlstSection%p_Rvalues)
    deallocate(rparlstSection%p_Sparameters)
    rparlstSection%iparamCount = 0

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Reallocate the section list
  ! This increases the size of a section list by reallocation of the
  ! arrays.

  subroutine INIP_realloclist (rparlist, inewsize)

    ! The section list to reallocate.
    type(t_parlist), intent(inout) :: rparlist

    ! The new 'size' of the section, i.e. the new number of parameters,
    ! the section should be able to handle.
    integer, intent(in) :: inewsize

    ! local variables

    integer :: sz

    ! Pointers to new lists for replacing the old.
    type(t_parlstSection), dimension(:), pointer :: p_Rsections

    ! Allocate the pointers for the new lists
    allocate(p_Rsections(inewsize))

    sz = min(size(rparlist%p_Rsections),inewsize)

    ! Copy the content of the old ones
    p_Rsections(1:sz) = rparlist%p_Rsections (1:sz)

    ! Throw away the old arrays, replace by the new ones
    deallocate(rparlist%p_Rsections)

    rparlist%p_Rsections => p_Rsections

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Search in a section for a parameter
  ! and return the index - or 0 if the parameter does not exist.

  subroutine INIP_fetchparameter(rsection, sname, iparamnum)

    ! The section.
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name to look for. Must be uppercase.
    character(LEN=*), intent(in) :: sname

    ! The number of the parameter in the list or 0 if it does not exist.
    integer, intent(out) :: iparamnum

    ! local variables
    integer :: i

    iparamnum = 0

    ! If the parameter list is empty, the section does not exist for sure
    if (rsection%iparamCount .eq. 0) return

    ! Loop through all sections to see if the section exists
    do i=1,rsection%iparamCount
      if (rsection%p_Sparameters(i) .eq. sname) then
        iparamnum = i
        return
      end if
    end do

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_init (rparlist)

    !<description>

    ! This routine initialises a parameter list. It must be applied to a
    ! parameter list structure before doing anything to it, just to initialise.

    !</description>

    !<inputoutput>

    ! The parameter list to initialise.
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !</subroutine>

    ! Set the section-count to 1.
    rparlist%isectionCount = 1

    ! Allocate a first set of sections
    allocate(rparlist%p_Rsections(INIP_NSECTIONS))

    ! Initialise the first section - it is the unnamed one.
    call INIP_initsection (rparlist%p_Rsections(1),'')

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_clear (rparlist)

    !<description>
    ! This routine cleans up a parameter list. All parameters in rparlist are
    ! removed.
    !</description>

    !<inputoutput>
    ! The parameter list to clean up.
    type(t_parlist), intent(inout) :: rparlist
    !</inputoutput>

    !</subroutine>

    ! Clean up = done+reinit. We make that simple here...
    call INIP_done (rparlist)
    call INIP_init (rparlist)

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_done (rparlist)

    !<description>

    ! This routine releases a parameter list. All memory allocated by the
    ! parameter list is released.

    !</description>

    !<inputoutput>

    ! The parameter list to release.
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !</subroutine>

    ! local variables
    integer :: i

    ! Probably nothing to do
    if (rparlist%isectionCount .eq. 0) return

    ! Loop through the parameter lists and release the content
    do i=rparlist%isectionCount,1,-1
      call INIP_releasesection (rparlist%p_Rsections(i))
    end do

    ! Release all sections
    deallocate(rparlist%p_Rsections)

    ! Mark the structure as 'empty', finish
    rparlist%isectionCount = 0

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_querysection(rparlist, sname, p_rsection)

    !<description>

    ! Searches for a section and return a pointer to it -
    ! or NULL() of the section does not exist.

    !</description>

    !<input>

    ! The parameter list to scan for the section.
    type(t_parlist), intent(in) :: rparlist

    ! The section name to look for.
    character(LEN=*), intent(in) :: sname

    !</input>

    !<output>

    ! A pointer to the section.
    type(t_parlstSection), pointer :: p_rsection

    !</output>

    !</subroutine>

    ! local variables
    integer :: i
    character(LEN=INIP_MLSECTION) :: sectionname

    nullify(p_rsection)

    ! If the parameter list is empty, the section does not exist for sure
    if (rparlist%isectionCount .eq. 0) return

    ! If the section name is '', return a pointer to the first section.
    if (sname .eq. '') then
      p_rsection => rparlist%p_Rsections(1)
      return
    end if

    ! Create the upper-case section name
    sectionname = adjustl(sname)
    call inip_toupper (sectionname)

    ! Loop through all sections to see if the section exists
    do i=1,rparlist%isectionCount
      if (rparlist%p_Rsections(i)%ssectionName .eq. sectionname) then
        p_rsection => rparlist%p_Rsections(i)
        return
      end if
    end do

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_addsection (rparlist, sname)

    !<description>

    ! Adds a section with the name sname to the list of sections in the
    ! parameter list rparlist. The name must NOT contain brackets ('[',']')
    ! in front and at the end!

    !</description>

    !<inputoutput>

    ! The parameter list where to add the section.
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !<input>

    ! The section name to add - without brackets in front and at the end!
    character(LEN=*), intent(in) :: sname

    !</input>

    !</subroutine>

    ! local variables
    character(LEN=INIP_MLSECTION) :: sectionname

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Create the upper-case section name
    sectionname = adjustl(sname)
    call inip_toupper (sectionname)

    ! Add a new section - reallocate the section list if necessary
    if (rparlist%isectionCount .eq. size(rparlist%p_Rsections)) then
      call INIP_realloclist (rparlist, size(rparlist%p_Rsections)+INIP_NSECTIONS)
    end if
    rparlist%isectionCount = rparlist%isectionCount + 1

    ! Initialise the new section.
    call INIP_initsection(rparlist%p_Rsections(rparlist%isectionCount),sectionname)

  end subroutine

  ! ***************************************************************************

  !<function>

  integer function INIP_queryvalue_indir (rsection, sparameter) &
      result (exists)

    !<description>
    ! Checks whether a parameter sparameter exists in the section rsection.
    !</description>

    !<result>
    ! The index of the parameter in the section ssection or =0, if the
    ! parameter does not exist within the section.
    !</result>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    !</input>

    !</function>

    ! local variables
    character(LEN=INIP_MLNAME) :: paramname

    exists = 0

    if (sparameter .eq. '') then
      call inip_output_line ('Empty parameter name!')
      call inip_sys_halt()
    end if

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter index into 'exists', finish.
    call INIP_fetchparameter(rsection, paramname, exists)

  end function

  ! ***************************************************************************

  !<function>

  integer function INIP_queryvalue_direct (rparlist, ssectionName, sparameter) &
      result (exists)

    !<description>
    ! Checks whether a parameter sparameter exists in the section ssectionname
    ! in the parameter list rparlist.
    !</description>

    !<result>
    ! The index of the parameter in the section ssectionName or =0, if the
    ! parameter does not exist within the section.
    !</result>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    !</input>

    !</function>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    exists = 0

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call inip_output_line ('Section not found: '//trim(ssectionName))
      return
    end if

    ! Search for the parameter
    exists = INIP_queryvalue_indir (p_rsection, sparameter)

  end function

  ! ***************************************************************************

  !<function>

  integer function INIP_querysubstrings_indir (rsection, sparameter) &
      result (iresult)

    !<description>
    ! Returns the number of substrings of a parameter.
    !</description>

    !<result>
    ! The number of substrings of parameter sparameter in section rsection.
    !</result>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    !</input>

    !</function>

    ! local variables
    integer :: idx
    character(LEN=INIP_MLNAME) :: paramname

    if (sparameter .eq. '') then
      call inip_output_line ('Empty parameter name!')
      call inip_sys_halt()
    end if

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter index into 'idx', finish.
    call INIP_fetchparameter(rsection, paramname, idx)

    ! Return number of substrings
    if (idx .eq. 0) then
      iresult = 0
    else
      iresult = rsection%p_Rvalues(idx)%nsize
    end if

  end function

  ! ***************************************************************************

  !<function>

  integer function INIP_querysubstrings_direct (rparlist, ssectionName, sparameter) &
      result (iresult)

    !<description>
    ! Checks whether a parameter sparameter exists in the section ssectionname
    ! in the parameter list rparlist.
    !</description>

    !<result>
    ! The index of the parameter in the section ssectionName or =0, if the
    ! parameter does not exist within the section.
    !</result>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    !</input>

    !</function>

    ! local variables
    integer :: idx
    type(t_parlstSection), pointer :: p_rsection

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call inip_output_line ('Section not found: '//trim(ssectionName))
      return
    end if

    ! Get the parameter index
    idx = INIP_queryvalue_indir (p_rsection, sparameter)

    ! Return number of substrings
    if (idx .eq. 0) then
      iresult = 0
    else
      iresult = p_rsection%p_Rvalues(idx)%nsize
    end if

  end function

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_string_indir (rsection, sparameter, svalue, &
      sdefault, isubstring, bdequote)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, sdefault is returned.
    ! If sdefault is not given, an error will be thrown.
    !
    ! If the value is an array of strings, the optional parameter isubstring>=0
    ! allows to specify the number of the substring to be returned;
    ! isubstring=0 returns the value directly
    ! behind the '=' sign in the line of the parameter, isubstring>0 returns
    ! the array-entry in the lines below the parameter.
    !
    ! When omitting isubstring, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    character(LEN=*), intent(in), optional :: sdefault

    ! OPTIONAL: The number of the substring to be returned.
    ! =0: returns the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns substring isubstring.
    integer, intent(in), optional :: isubstring

    ! OPTIONAL: De-quote the string.
    ! This provides a save way of removing quotation marks around a string in case
    ! the parameter contains exactly one string.
    ! =false: Return the string as it is (standard)
    ! =true: Re-read the string and remove any leading and trailing quotation marks
    !   (if there are any).
    logical, intent(in), optional :: bdequote
    !</input>

    !<output>

    ! The value of the parameter
    character(LEN=*), intent(out) :: svalue

    !</output>

    !</subroutine>

    ! local variables
    integer :: i,isub
    character(LEN=INIP_MLNAME) :: paramname

    if (sparameter .eq. '') then
      call inip_output_line ('Empty parameter name!')
      call inip_sys_halt()
    end if

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter index into 'exists', finish.
    call INIP_fetchparameter(rsection, paramname, i)

    if (i .eq. 0) then
      if (present(sdefault)) then
        svalue = sdefault
      else
        call inip_output_line ('Parameter not found: '//trim(paramname))
        call inip_sys_halt()
      end if
    else
      ! Depending on isubstring, return either the 'headline' or one
      ! of the substrings.
      isub = 0
      if (present(isubstring)) isub = isubstring

      if ((isub .le. 0) .or. (isub .gt. rsection%p_Rvalues(i)%nsize)) then
        call inip_chararraytostring(rsection%p_Rvalues(i)%p_sentry,svalue)
      else
        call inip_chararraytostring(rsection%p_Rvalues(i)%p_SentryList(:,isub),svalue)
      end if
    end if

    if (present(bdequote)) then
      if (bdequote) then
        call inip_dequote(svalue)
      end if
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_string_fetch (rsection, iparameter, svalue,&
      bexists, isubstring, bdequote)
    !<description>

    ! Returns the value of a parameter in the section rsection.
    ! iparameter specifies the number of the parameter in section rsection.
    ! If bexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    !
    ! If bexists is given, it will be set to TRUE if the parameter number
    ! iparameter exists, otherwise it will be set to FALSE and svalue=''.
    !
    ! If the value is an array of strings, the optional parameter isubstring>=0
    ! allows to specify the number of the substring to be returned;
    ! isubstring=0 returns the value directly
    ! behind the '=' sign in the line of the parameter, isubstring>0 returns
    ! the array-entry in the lines below the parameter.
    !
    ! When omitting isubstring, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The number of the parameter.
    integer, intent(in) :: iparameter

    ! OPTIONAL: The number of the substring to be returned.
    ! =0: returns the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns substring isubstring.
    integer, intent(in), optional :: isubstring

    ! OPTIONAL: De-quote the string.
    ! This provides a save way of removing quotation marks around a string in case
    ! the parameter contains exactly one string.
    ! =false: Return the string as it is (standard)
    ! =true: Re-read the string and remove any leading and trailing quotation marks
    !   (if there are any).
    logical, intent(in), optional :: bdequote
    !</input>

    !<output>

    ! The value of the parameter
    character(LEN=*), intent(out) :: svalue

    ! OPTIONAL: Parameter existance check
    ! Is set to TRUE/FALSE, depending on whether the parameter exists.
    logical, intent(out), optional :: bexists

    !</output>

    !</subroutine>

    integer :: isub

    ! Check if iparameter is out of bounds. If yes, probably
    ! throw an error.

    if ((iparameter .lt. 0) .or. (iparameter .gt. rsection%iparamCount)) then

      if (.not. present(bexists)) then
        call inip_output_line ('Error. Parameter '//trim(inip_siL(iparameter,10))//&
          ' does not exist!')
        call inip_sys_halt()
      else
        svalue = ''
        bexists = .false.
        return
      end if

    end if

    ! Get the parameter value.
    ! Depending on isubstring, return either the 'headline' or one
    ! of the substrings.
    isub = 0
    if (present(isubstring)) isub = isubstring

    if ((isub .le. 0) .or. &
        (isub .gt. rsection%p_Rvalues(iparameter)%nsize)) then
      call inip_charArrayToString(rsection%p_Rvalues(iparameter)%p_sentry,svalue)
    else
      call inip_charArrayToString(rsection%p_Rvalues(iparameter)%p_SentryList(:,isub),svalue)
    end if

    if (present(bexists)) bexists = .true.

    if (present(bdequote)) then
      if (bdequote) then
        call inip_dequote(svalue)
      end if
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_string_direct (rparlist, ssectionName, &
      sparameter, svalue, sdefault,&
      isubstring,bdequote)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, sdefault is returned.
    ! If sdefault is not given, an error will be thrown.
    !
    ! If the value is an array of strings, the optional parameter isubstring>=0
    ! allows to specify the number of the substring to be returned;
    ! isubstring=0 returns the value directly
    ! behind the '=' sign in the line of the parameter, isubstring>0 returns
    ! the array-entry in the lines below the parameter.
    !
    ! When omitting isubstring, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    character(LEN=*), intent(in), optional :: sdefault

    ! OPTIONAL: The number of the substring to be returned.
    ! =0: returns the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns substring isubstring.
    integer, intent(in), optional :: isubstring

    ! OPTIONAL: De-quote the string.
    ! This provides a save way of removing quotation marks around a string in case
    ! the parameter contains exactly one string.
    ! =false: Return the string as it is (standard)
    ! =true: Re-read the string and remove any leading and trailing quotation marks
    !   (if there are any).
    logical, intent(in), optional :: bdequote

    !</input>

    !<output>

    ! The value of the parameter
    character(LEN=*), intent(out) :: svalue

    !</output>

    !</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      if (present(sdefault)) then
        svalue = sdefault
        return
      else
        call inip_output_line ('Section not found: '//trim(ssectionName))
        call inip_sys_halt()
      end if
    end if

    ! Get the value
    call INIP_getvalue_string_indir (p_rsection, sparameter, svalue, sdefault,&
      isubstring,bdequote)

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_single_indir (rsection, sparameter, fvalue, &
      fdefault, iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, idefault is returned.
    ! If idefault is not given, an error will be thrown.
    !
    ! If the value is an array of singles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the single to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    real*4, intent(in), optional :: fdefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*4, intent(out) :: fvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(fdefault)) then
      write (sdefault,'(E17.10E2)') fdefault
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        sdefault, iarrayindex)
    else
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    fvalue = inip_StringToSingle(svalue,'(E17.10E2)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_getvalue_single_fetch (rsection, iparameter, fvalue, &
      bexists, iarrayindex)

    !<description>

    ! Returns the value of a parameter in the section rsection.
    ! iparameter specifies the number of the parameter in section rsection.
    !
    ! If bexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    ! If bexists is given, it will be set to TRUE if the parameter number
    ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
    !
    ! If the value is an array of singles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the single to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The number of the parameter.
    integer, intent(in) :: iparameter

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*4, intent(out) :: fvalue

    ! OPTIONAL: Parameter existance check
    ! Is set to TRUE/FALSE, depending on whether the parameter exists.
    logical, intent(out), optional :: bexists

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: svalue

    svalue = '0.0E0'
    call INIP_getvalue_string_fetch (rsection, iparameter, svalue, &
      bexists, iarrayindex)

    fvalue = inip_StringToSingle(svalue,'(E17.10E2)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_single_direct (rparlist, ssectionName, &
      sparameter, fvalue, fdefault,&
      iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection. If the
    ! value does not exist, ddefault is returned.  If fdefault is not
    ! given, an error will be thrown.
    !
    ! If the value is an array of singles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the single to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.
    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    real*4, intent(in), optional :: fdefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*4, intent(out) :: fvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(fdefault)) then
      write (sdefault,'(E17.10E2)') fdefault
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, sdefault, &
        iarrayindex)
    else
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    fvalue = inip_StringToSingle(svalue,'(E17.10E2)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_double_indir (rsection, sparameter, dvalue, &
      ddefault, iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, ddefault is returned.
    ! If ddefault is not given, an error will be thrown.
    !
    ! If the value is an array of doubles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the double to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    real*8, intent(in), optional :: ddefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*8, intent(out) :: dvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(ddefault)) then
      write (sdefault,'(E27.19E3)') ddefault
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        sdefault, iarrayindex)
    else
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    dvalue = inip_StringToDouble(svalue,'(E27.19E3)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_getvalue_double_fetch (rsection, iparameter, dvalue, &
      bexists, iarrayindex)

    !<description>

    ! Returns the value of a parameter in the section rsection.
    ! iparameter specifies the number of the parameter in section rsection.
    !
    ! If bexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    ! If bexists is given, it will be set to TRUE if the parameter number
    ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
    !
    ! If the value is an array of doubles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the double to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The number of the parameter.
    integer, intent(in) :: iparameter

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*8, intent(out) :: dvalue

    ! OPTIONAL: Parameter existance check
    ! Is set to TRUE/FALSE, depending on whether the parameter exists.
    logical, intent(out), optional :: bexists

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: svalue

    svalue = '0.0E0'
    call INIP_getvalue_string_fetch (rsection, iparameter, svalue, &
      bexists, iarrayindex)

    dvalue = inip_StringToDouble(svalue,'(E27.19E3)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_double_direct (rparlist, ssectionName, &
      sparameter, dvalue, ddefault,&
      iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, ddefault is returned.
    ! If ddefault is not given, an error will be thrown.
    !
    ! If the value is an array of doubles, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the double to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    real*8, intent(in), optional :: ddefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    real*8, intent(out) :: dvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(ddefault)) then
      write (sdefault,'(E27.19E3)') ddefault
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, sdefault, &
        iarrayindex)
    else
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    dvalue = inip_StringToDouble(svalue,'(E27.19E3)')

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_logical_indir (rsection, sparameter, bvalue, &
      bdefault, iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, bdefault is returned.
    ! If bdefault is not given, an error will be thrown.
    !
    ! If the value is an array of logicals, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the logical to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    logical, intent(in), optional :: bdefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    logical, intent(out) :: bvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(bdefault)) then
      write (sdefault,'(L1)') bdefault
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        sdefault, iarrayindex)
    else
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    read (svalue,'(L1)') bvalue

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_getvalue_logical_fetch (rsection, iparameter, bvalue, &
      bexists, iarrayindex)

    !<description>

    ! Returns the value of a parameter in the section rsection.
    ! iparameter specifies the number of the parameter in section rsection.
    !
    ! If bexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    ! If bexists is given, it will be set to TRUE if the parameter number
    ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
    !
    ! If the value is an array of logicals, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the logical to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The number of the parameter.
    integer, intent(in) :: iparameter

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    logical, intent(out) :: bvalue

    ! OPTIONAL: Parameter existance check
    ! Is set to TRUE/FALSE, depending on whether the parameter exists.
    logical, intent(out), optional :: bexists

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: svalue

    write (svalue,'(L1)') .false.
    call INIP_getvalue_string_fetch (rsection, iparameter, svalue, &
      bexists, iarrayindex)

    read (svalue,'(L1)') bvalue

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_logical_direct (rparlist, ssectionName, &
      sparameter, bvalue, bdefault,&
      iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, bdefault is returned.
    ! If bdefault is not given, an error will be thrown.
    !
    ! If the value is an array of logicals, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the logical to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    logical, intent(in), optional :: bdefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    logical, intent(out) :: bvalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(bdefault)) then
      write (sdefault,'(L1)') bdefault
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, sdefault, &
        iarrayindex)
    else
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    read (svalue,'(L1)') bvalue

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_int_fetch (rsection, iparameter, ivalue, &
      bexists, iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section rsection.
    ! iparameter specifies the number of the parameter in section rsection.
    ! If bexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    !
    ! If bexists is given, it will be set to TRUE if the parameter number
    ! iparameter exists, otherwise it will be set to FALSE and ivalue=0.
    !
    ! If the value is an array of integers, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the integer to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The number of the parameter.
    integer, intent(in) :: iparameter

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    integer, intent(out) :: ivalue

    ! OPTIONAL: Parameter existance check
    ! Is set to TRUE/FALSE, depending on whether the parameter exists.
    logical, intent(out), optional :: bexists

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: svalue

    svalue = '0'
    call INIP_getvalue_string_fetch (rsection, iparameter, svalue, &
      bexists, iarrayindex)
    read(svalue,*) ivalue

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_getvalue_int_direct (rparlist, ssectionName, &
      sparameter, ivalue, idefault,&
      iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection. If the
    ! value does not exist, idefault is returned.  If idefault is not
    ! given, an error will be thrown.
    !
    ! If the value is an array of integers, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the integer to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.

    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    integer, intent(in), optional :: idefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    integer, intent(out) :: ivalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(idefault)) then
      write (sdefault,*) idefault
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, sdefault, &
        iarrayindex)
    else
      call INIP_getvalue_string_direct (rparlist, ssectionName, &
        sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    read(svalue,*) ivalue

  end subroutine


  !<subroutine>
  subroutine INIP_getvalue_int_indir (rsection, sparameter, ivalue, &
      idefault, iarrayindex)
    !<description>

    ! Returns the value of a parameter in the section ssection.
    ! If the value does not exist, idefault is returned.
    ! If idefault is not given, an error will be thrown.
    !
    ! If the value is an array of integers, the optional parameter
    ! iarrayindex>=0 allows to specify the number of the integer to be
    ! returned; iarrayindex=0 returns the value directly behind the '='
    ! sign in the line of the parameter, iarrayindex>0 returns the
    ! array-entry in the lines below the parameter.
    !
    ! When omitting iarrayindex, the value directly behind the '=' sign
    ! is returned.
    !</description>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! OPTIONAL: A default value
    integer, intent(in), optional :: idefault

    ! OPTIONAL: The number of the arrayindex to be returned.
    ! =0: returns the integer directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: returns array index  iarrayindex.
    integer, intent(in), optional :: iarrayindex

    !</input>

    !<output>

    ! The value of the parameter
    integer, intent(out) :: ivalue

    !</output>

    !</subroutine>

    ! local variables
    character (LEN=INIP_LENLINEBUF) :: sdefault,svalue

    ! Call the string routine, perform a conversion afterwards.
    if (present(idefault)) then
      write (sdefault,*) idefault
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        sdefault, iarrayindex)
    else
      call INIP_getvalue_string_indir (rsection, sparameter, svalue, &
        isubstring=iarrayindex)
    end if

    read(svalue,*) ivalue

  end subroutine


  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_addvalue_indir (rsection, sparameter, svalue, nsubstrings)

    !<description>
    ! Adds a parameter to a section rsection.
    ! If the parameter exists, it is overwritten.
    !</description>

    !<inputoutput>

    ! The section where to arr the parameter
    type(t_parlstSection), intent(inout) :: rsection

    !</inputoutput>

    !<input>

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! The value of the parameter
    character(LEN=*), intent(in) :: svalue

    ! OPTIONAL: Number of substrings. This allows a parameter to have
    ! multiple substrings, which can be accessed via the 'isubstring'
    ! parameter in the GET-routines.
    integer, intent(in), optional :: nsubstrings

    !</input>

    !</subroutine>

    ! local variables
    character(LEN=INIP_MLNAME) :: paramname
    integer :: i,j

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter index into 'exists', finish.
    call INIP_fetchparameter(rsection, paramname, i)

    if (i .eq. 0) then

      ! Does not exist. Append.
      !
      ! Enough space free? Otherwise reallocate the parameter list
      if (rsection%iparamCount .eq. size(rsection%p_Sparameters)) then
        call INIP_reallocsection (rsection, size(rsection%p_Sparameters)+INIP_NPARSPERBLOCK)
      end if

      ! Add the parameter - without any adjustment of the 'value' string
      rsection%iparamCount = rsection%iparamCount + 1

      ! Set i to the index of the parameter
      i = rsection%iparamCount

    else

      ! Check if there are substrings. If yes, deallocate.
      ! Will be allocated later if necessary.
      if (associated(rsection%p_Rvalues(i)%p_SentryList)) then
        deallocate(rsection%p_Rvalues(i)%p_SentryList)
        rsection%p_Rvalues(i)%nsize = 0
      end if

      ! Deallocate memory for the value, will be reallocated later.
      if (associated(rsection%p_Rvalues(i)%p_sentry)) then
        deallocate(rsection%p_Rvalues(i)%p_sentry)
      end if

    end if

    rsection%p_Sparameters(i) = paramname
    j = len_trim(svalue)
    allocate(rsection%p_Rvalues(i)%p_sentry(max(1,j)))
    rsection%p_Rvalues(i)%p_sentry(:) = ' '
    call inip_stringtochararray(svalue,rsection%p_Rvalues(i)%p_sentry,j)

    ! Add a list for the substrings if the parameter should have substrings.
    if (present(nsubstrings)) then
      if (nsubstrings .gt. 0) then
        allocate(rsection%p_Rvalues(i)%p_SentryList(max(1,j),nsubstrings))
        rsection%p_Rvalues(i)%p_SentryList(:,:) = ' '
        rsection%p_Rvalues(i)%nsize = nsubstrings
      else
        nullify(rsection%p_Rvalues(i)%p_SentryList)
        rsection%p_Rvalues(i)%nsize = 0
      end if
    else
      nullify(rsection%p_Rvalues(i)%p_SentryList)
      rsection%p_Rvalues(i)%nsize = 0
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_addvalue_direct (rparlist, ssectionName, sparameter, svalue,&
      nsubstrings)
    !<description>

    ! Adds a parameter to a section with name ssectionName in the parameter list
    ! rparlist. If ssectionName='', the parameter is added to the unnamed
    ! section.

    !</description>

    !<inputoutput>

    ! The parameter list.
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !<input>

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! The value of the parameter
    character(LEN=*), intent(in) :: svalue

    ! OPTIONAL: Number of substrings. This allows a parameter to have
    ! multiple substrings, which can be accessed via the 'isubstring'
    ! parameter in the GET-routines.
    integer, intent(in), optional :: nsubstrings

    !</input>

    !</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call inip_output_line ('Section not found: '//trim(ssectionName))
      return
    end if

    ! Add the parameter

    call INIP_addvalue_indir (p_rsection, sparameter, svalue, nsubstrings)

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_setvalue_fetch (rsection, iparameter, svalue, iexists,&
      isubstring)

    !<description>

    ! Modifies the value of a parameter in the section rsection.
    ! The value of parameter iparameter in the section rsection is modified.
    ! If iexists does not appear, an error is thrown if a nonexisting
    ! parameter is accessed.
    ! If iexists is given, it will be set to YES if the parameter number
    ! iparameter exists and was modified, otherwise it will be set to NO.
    !
    ! isubstring allows to specify the numer of a substring of the parameter to
    ! change. If omitted or = 0, the 'headline' directly behind the '='
    ! sign of the line 'name=value' is modified. Otherwise, the corresponding
    ! substring is changed.

    !</description>

    !<inputoutput>

    ! The section where to arr the parameter
    type(t_parlstSection), intent(inout) :: rsection

    !</inputoutput>

    !<input>

    ! The parameter name.
    integer, intent(in) :: iparameter

    ! The new value of the parameter
    character(LEN=*), intent(in) :: svalue

    ! OPTIONAL: The number of the substring to be changed.
    ! =0: changes the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: changes substring isubstring.
    integer, intent(in), optional :: isubstring

    !</input>

    !<output>

    ! Optional parameter. Is set to YES/NO, depending on whether
    ! the parameter exists.
    integer, intent(out), optional :: iexists

    !</output>

    !</subroutine>

    integer :: isub,j

    ! Check if iparameter is out of bounds. If yes, probably
    ! throw an error.

    if ((iparameter .lt. 0) .or. (iparameter .gt. rsection%iparamCount)) then

      if (.not. present(iexists)) then
        call inip_output_line ('Error. Parameter '//trim(inip_siL(iparameter,10))//&
          ' does not exist!')
        call inip_sys_halt()
      else
        iexists = NO
        return
      end if

    end if

    ! Depending on isubstring, change either the 'headline' or one
    ! of the substrings.
    isub = 0
    if (present(isubstring)) isub = isubstring

    j = len_trim(svalue)
    if ((isub .le. 0) .or. &
        (isub .gt. rsection%p_Rvalues(iparameter)%nsize)) then
      ! Reallocate memory
      deallocate(rsection%p_Rvalues(iparameter)%p_sentry)
      allocate(rsection%p_Rvalues(iparameter)%p_sentry(max(1,j)))
      rsection%p_Rvalues(iparameter)%p_sentry(:) = ' '
      call inip_stringToCharArray(svalue,rsection%p_Rvalues(iparameter)%p_sentry,j)
    else
      ! Check that there is enough memory to save the string.
      if (ubound(rsection%p_Rvalues(iparameter)%p_SentryList,1) .le. j) then
        call INIP_reallocSubVariables(rsection%p_Rvalues(iparameter),j)
      end if
      call inip_stringToCharArray(svalue,rsection%p_Rvalues(iparameter)%p_SentryList(:,isub),j)
    end if

    if (present(iexists)) iexists = YES

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_setvalue_indir (rsection, sparameter, svalue, isubstring)

    !<description>

    ! Modifies the value of a parameter in the section rsection.
    ! If the parameter does not exist, an error is thrown.
    !
    ! isubstring allows to specify the numer of a substring of the parameter to
    ! change. If omitted or = 0, the 'headline' directly behind the '='
    ! sign of the line 'name=value' is modified. Otherwise, the corresponding
    ! substring is changed.

    !</description>

    !<inputoutput>

    ! The section where to arr the parameter
    type(t_parlstSection), intent(inout) :: rsection

    !</inputoutput>

    !<input>

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! The new value of the parameter
    character(LEN=*), intent(in) :: svalue

    ! OPTIONAL: The number of the substring to be changed.
    ! =0: changes the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: changes substring isubstring.
    integer, intent(in), optional :: isubstring

    !</input>

    !</subroutine>

    ! local variables
    integer :: i,isub,j
    character(LEN=INIP_MLNAME) :: paramname

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter position
    i = INIP_queryvalue_indir (rsection, paramname)

    if (i .eq. 0) then
      call inip_output_line ('Parameter '//trim(paramname)//&
        ' does not exist, cannot be modified!')
      call inip_sys_halt()
    else

      ! Depending on isubstring, change either the 'headline' or one
      ! of the substrings.
      isub = 0
      if (present(isubstring)) isub = isubstring

      j = len_trim(svalue)
      if ((isub .le. 0) .or. (isub .gt. rsection%p_Rvalues(i)%nsize)) then
        ! Reallocate memory
        deallocate(rsection%p_Rvalues(i)%p_sentry)
        allocate(rsection%p_Rvalues(i)%p_sentry(max(1,j)))
        rsection%p_Rvalues(i)%p_sentry(:) = ' '
        call inip_stringToCharArray(svalue,rsection%p_Rvalues(i)%p_sentry,j)
      else
        ! Check that there is enough memory to save the string.
        if (ubound(rsection%p_Rvalues(i)%p_SentryList,1) .le. j) then
          call INIP_reallocSubVariables(rsection%p_Rvalues(i),j)
        end if
        call inip_stringToCharArray(svalue,rsection%p_Rvalues(i)%p_SentryList(:,isub),j)
      end if

    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  subroutine INIP_setvalue_direct (rparlist, ssectionName, sparameter, svalue,&
      isubstring)
    !<description>

    ! Modifies the value of a parameter in the section with name ssectionName
    ! in the parameter list rparlist.
    ! If the parameter does not exist, an error is thrown.
    !
    ! isubstring allows to specify the numer of a substring of the parameter to
    ! change. If omitted or = 0, the 'headline' directly behind the '='
    ! sign of the line 'name=value' is modified. Otherwise, the corresponding
    ! substring is changed.

    !</description>

    !<inputoutput>

    ! The parameter list.
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !<input>

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name.
    character(LEN=*), intent(in) :: sparameter

    ! The new value of the parameter
    character(LEN=*), intent(in) :: svalue

    ! OPTIONAL: The number of the substring to be changed.
    ! =0: changes the string directly behind the '=' sign in the line
    !     'name=value'.
    ! >0: changes substring isubstring.
    integer, intent(in), optional :: isubstring

    !</input>

    !</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call inip_output_line ('Section not found: '//trim(ssectionName))
      return
    end if

    ! Set the parameter

    call INIP_setvalue_indir (p_rsection, sparameter, svalue, isubstring)

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Read a line from a text file.

  subroutine INIP_readlinefromfile (iunit, sdata, ilinelen, ios)

    ! The unit where to read from; must be connected to a file.
    integer, intent(in) :: iunit

    ! The string where to write data to
    character(LEN=*), intent(out) :: sdata

    ! Length of the output
    integer, intent(out) :: ilinelen

    ! Status of the reading process. Set to a value <> 0 if the end
    ! of the file is reached.
    integer, intent(out) :: ios

    ! local variables
    character :: c

    sdata = ''
    ilinelen = 0

    ! Read the data - as long as the line/file does not end.
    do

      ! Read a character.
      ! Unfortunately, Fortran forces me to use this dirty GOTO
      ! to decide processor-independently whether the line or
      ! the record ends.
      read (unit=iunit, fmt='(A1)', iostat=ios, advance='NO',&
        end=10, eor=20) c

      ! Do not do anything in case of an error
      if (ios .eq. 0) then

        ilinelen = ilinelen + 1
        sdata (ilinelen:ilinelen) = c

      end if

      ! Proceed to next character
      cycle

      ! End of file.
      10  ios = -1
      exit

      ! End of record = END OF LINE.
      20  ios = 0
      exit

    end do

  end subroutine

  ! ***************************************************************************

  ! Internal subroutine: Parse a text line.
  ! This parses the text line sdata.
  ! Return values:
  !  ityp = 0 -> The line is a comment
  !  ityp = 1 -> The line is a section. ssecname is the name of the section
  !              without '[]'
  !  ityp = 2 -> The line is a parameter. sparamname is the uppercase
  !              parameter name. svalue is the value of the parameter,
  !              trimmed and left adjusted.
  !  ityp = 3 -> Line is the beginning of a multi-valued parameter.
  !              The next isubstring lines contain additional substrings.
  !  ityp = 4 -> Line is a substring of a multi-valued parameter.

  subroutine INIP_parseline (sdata, ityp, isubstring, ilinenum, &
      ssecname, sparamname, svalue, sfilename)

    ! The line to be parsed
    character(LEN=*), intent(in) :: sdata

    ! The typ of the line
    integer, intent(out) :: ityp

    ! input: =0: parse line as parameter. isubstring is changed to a value > 0
    !            is the parameter has multiple values attached.
    !        >0: parse line as substring of a multi-valued parameter, not
    !            containing a leading 'name='.
    ! output: If the 'headline' of a multi-valued parameter is read, isubstring is
    !         changed to the number of substrings (the k in 'name(k)=...').
    !         Otherwise unchanged.
    integer, intent(inout) :: isubstring

    ! Line number
    integer, intent(in) :: ilinenum

    ! Section name, if it is a section
    character(LEN=*), intent(inout) :: ssecname

    ! Parameter name, if it is a parameter
    character(LEN=*), intent(inout) :: sparamname

    ! Parameter value, if it is a parameter
    character(LEN=*), intent(inout) :: svalue

    ! OPTIONAL: Filename of the file to be parsed.
    ! Will be printed out in error messages if present.
    character(LEN=*), intent(in), optional :: sfilename

    ! local variables
    integer :: i,j1,j2,ltr
    character(LEN=INIP_LENLINEBUF) :: sbuf,slen

    ityp = 0

    ! Do we have data in sdata?
    if (sdata .eq. '') return

    ! Copy the input string - left adjusted - and get the string length
    sbuf = adjustl(sdata)

    ! Should we parse the line as first line of a parameter or as substring
    ! of a multi-valued parameter?
    if (isubstring .eq. 0) then

      ! Standard parameter or section header.
      !
      ! Do we start with '[' and end with ']'?
      if (sbuf(1:1) .eq. "[") then

        ! Find the final ']'.
        do ltr = 1,len(sbuf)
          if (sbuf(ltr:ltr) .eq. "]") exit
        end do

        if (sbuf(ltr:ltr) .ne. ']') then
          if (present(sfilename)) then
            call inip_output_line ('File: '//trim(sfilename))
          end if
          call inip_output_line ('Wrong syntax of section name. Line '//&
            trim(inip_siL(ilinenum,10))//':')
          call inip_output_line (sbuf)
          call inip_sys_halt()
        end if

        ! Get the section name
        ssecname = sbuf(2:ltr-1)
        ityp = 1
        return

      else if (sbuf(1:1) .eq. INIP_COMMENT) then

        ! Comment sign
        return

      else

        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        call linelength(sbuf, ltr)

        ! ltr=0 means: empty line. Ignore that.
        if (ltr .eq. 0) return

        ! Is there a '(..)' that is indicating a multi-valued parameter?
        j1 = index(sbuf(1:ltr),'(')
        j2 = index(sbuf(1:ltr),')')

        ! Is there a '=' sign?
        i = index(sbuf(1:ltr),'=')

        if (i .eq. 0) then
          if (present(sfilename)) then
            call inip_output_line ('File: '//trim(sfilename))
          end if
          call inip_output_line ('Invalid parameter syntax. Line '&
            //trim(inip_siL(ilinenum,10))//':')
          call inip_output_line (trim(sbuf))
          call inip_sys_halt()
        end if

        if ((j1 .eq. 0) .or. (j2 .le. j1) .or. (i .le. j1)) then

          ityp = 2

          ! Get the name of the parameter
          sparamname = adjustl(sbuf(1:i-1))

          ! Get the parameter value
          svalue = adjustl(sbuf(i+1:ltr))

        else

          ! Probably multi-valued parameter with substrings in the
          ! following lines.

          ! Get the name of the parameter
          sparamname = adjustl(sbuf(1:j1-1))

          ! Get the parameter value
          svalue = adjustl(sbuf(i+1:ltr))

          ! Get the length of the parameter list.
          slen = sbuf (j1+1:min(j2-1,len(slen)))

          isubstring = 0
          read(slen,*) isubstring

          if (isubstring .le. 0) then
            ! Oh, only one line. User wants to cheat :-)
            isubstring = 0

            ityp = 2
          else
            ! Real multi-valued parameter.
            ityp = 3
          end if

        end if

      end if

    else

      ! Substring of a multi-valued parameter.
      if (sbuf(1:1) .eq. INIP_COMMENT) then

        ! Comment sign
        return

      else

        ! Must be a parameter. Get the length of the string without comment
        ! at the end.
        call linelength(sbuf, ltr)

        ! ltr=0 means: empty line. Ignore that.
        if (ltr .eq. 0) return

        ityp = 4

        ! Get the parameter value. Do not get a parameter name; there is none.
        svalue = adjustl(sbuf(1:ltr))

      end if

    end if

  contains

    ! Sub-subroutine: find the length of the line, removing comments
    ! at the end.

    subroutine linelength (sdata, l)

      ! The string to parse. Must not be ''!
      character(LEN=*), intent(in) :: sdata

      ! The index of the last character without any comment at the end.
      integer, intent(out) :: l

      ! local variables
      logical :: bflag   ! Set to true if we are in apostroph mode
      integer :: lsdata

      bflag = .false.

      ! Go through all characters
      l = 0
      lsdata = len(sdata)
      do while (l .lt. lsdata)

        ! next character
        l = l+1

        ! A comment character while we are not in apostroph mode? Stop.
        if ((.not. bflag) .and. (sdata(l:l) .eq. INIP_COMMENT)) then
          l = l-1
          exit
        end if

        ! An apostroph?
        if (sdata(l:l) .eq. "'") then

          ! Switch the apostroph mode.
          ! Btw.: Two subsequent apostrophes will switch the mode off and on again.
          bflag = .not. bflag

        end if

      end do

    end subroutine

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_readfromfile (rparlist, sfile, sdirectory, bexpandVars)

    !<description>

    ! This routine parses a text file for data of the INI-file form.
    ! sfile must be the name of a file on the hard disc.
    ! The file may have references to subfiles in the unnamed section:
    ! If there is a parameter "simportdatafiles(n)" present in the
    ! unnamed section, the routine expects a list of filenames in this
    ! parameter. All the files listed there are read in as well.
    ! Parameters at the end of the master file will overwrite
    ! the parameters from the files in simportdatafiles.
    !
    ! Sub-files specified in the main file are searched in the following
    ! directories:
    ! 1.) sdirectory (if specified)
    ! 2.) the directory that contains sfile (if sfile specifies
    !     a directory)
    ! 3.) current directory
    !
    ! The parameters read from the file(s) are added to the parameter list
    ! rparlist, which has to be initialised with INIP_init before
    ! calling the routine.
    !
    ! Remark: When adding parameters/sections to rparlist, the routine
    !   checks whether the parameters/sections already exist.
    !   Adding a parameter/section which exists does not result in an error -
    !   the first instance of the parameter/section will just be overwritten.

    !</description>

    !<inputoutput>

    ! The parameter list which is filled with data from the file
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !<input>

    ! The filename of the file to read.
    character(LEN=*), intent(in) :: sfile

    ! OPTIONAL: Directory containing other data files in case
    ! the main file sfile contains references to subfiles.
    character(LEN=*), intent(in), optional :: sdirectory

    ! OPTIONAL: Expand references to subvariables.
    ! TRUE: Subvariables like "%{section.varname} are expanded to actual values.
    !       This is the standard setting.
    ! FALSE: Subvariables are left as they are.
    logical, intent(in), optional :: bexpandVars
    !</input>

    ! local variables
    integer :: iidxsubfiles,icurrentsubfile
    character(LEN=INIP_LENLINEBUF), dimension(:), pointer :: p_Ssubfiles,p_SsubfilesTemp
    character(LEN=INIP_LENLINEBUF) :: sstring,smainfile
    integer :: nsubfiles,nnewsubfiles,j,ilensubf
    logical :: bexists,bmainpath,babsolute
    character(LEN=INIP_LENLINEBUF) :: smainpath,sfilepath,sfilename

    ! Filename/path of the master dat file.
    ! Search at first in the specified path.
    bexists = .false.
    if (present(sdirectory)) then
      smainfile = trim(sdirectory)//"/"//sfile
      inquire(file=smainfile, exist=bexists)
    end if

    if (.not. bexists) then
      smainfile = sfile
      inquire(file=smainfile, exist=bexists)
    end if

    if (.not. bexists) then
      ! Cancel if the file does not exist.
      return
    end if

    ! Get the main path/name of the file.
    call inip_pathExtract (smainfile, smainpath)
    bmainpath = smainpath .ne. ""

    ! Create a list of files to be read.
    ! They contain filenames including the directory.
    allocate(p_Ssubfiles(1))
    p_Ssubfiles(1) = smainfile

    icurrentsubfile = 0
    nsubfiles = 1

    ! Now read all files in the list. Append new files to the list if necessary.
    do while (icurrentsubfile .lt. nsubfiles)

      ! Read the unnamed section from the next file.
      icurrentsubfile = icurrentsubfile + 1

      ! Get the filename including the path.
      bexists = .false.

      call inip_pathExtract (trim(p_Ssubfiles(icurrentsubfile)), sfilepath, sfilename, babsolute)

      if (babsolute) then
        ! Ok, we have an absolute path given. Test it.
        sstring = trim(p_Ssubfiles(icurrentsubfile))
        inquire(file=sstring, exist=bexists)
      else
        ! Path is relative -- a little bit more complicated.
        ! Is there a directory given?
        if (sfilepath .ne. "") then
          ! Directory specified. We add "sdirectory" if we have it.
          if (present(sdirectory)) then
            sstring = trim(sdirectory)//"/"//trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          if (.not. bexists) then
            ! No, not there. Then take the path directly.
            sstring = trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          if (bmainpath .and. (.not. bexists)) then
            ! No, not there. Add the master directory and test there.
            sstring = trim(smainpath)//"/"//trim(sfilepath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if
        else
          ! No directory given. Then we search in the directory
          ! of the master file...
          if (bmainpath .and. (.not. bexists)) then
            sstring = trim(smainpath)//"/"//trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if

          ! And in the current directory.
          if (.not. bexists) then
            sstring = trim(sfilename)
            inquire(file=sstring, exist=bexists)
          end if
        end if
      end if

      if (bexists) then
        call INIP_readfromsinglefile (rparlist, sstring, .false., .false.)

        ! Replace the filename with the string including the path.
        ! Then the actual read process at the end of the routine can be
        ! handled easier.
        p_Ssubfiles(icurrentsubfile) = trim(sstring)

      else
        call inip_output_line ('Specified data-subfile does not exist: '//sstring)
      end if

      ! Check if there is a parameter "simportdatafiles" available in the
      ! parameter list. It must be present in the unnamed section.
      call INIP_fetchparameter(rparlist%p_Rsections(1), "SIMPORTDATAFILES", iidxsubfiles)

      if (iidxsubfiles .ne. 0) then
        ! Append the new files to the file list.
        !
        ! Get the number of new files. The parameter definitely exists as
        ! it was created when reading the 'master' INI file.
        nnewsubfiles = INIP_querysubstrings (rparlist%p_Rsections(1), &
          "SIMPORTDATAFILES")

        ! if nnewsubfiles=0, there is (hopefully) only one string here.
        if (nnewsubfiles .eq. 0) then
          call INIP_getvalue_string(rparlist%p_Rsections(1), iidxsubfiles, sstring,bdequote=.true.)
          if (trim(sstring) .ne. "") then
            ! Append the data.
            allocate(p_SsubfilesTemp(nsubfiles+1))
            p_SsubfilesTemp(1:nsubfiles) = p_Ssubfiles(:)
            deallocate(p_Ssubfiles)
            p_Ssubfiles => p_SsubfilesTemp
            nullify(p_SsubfilesTemp)

            ! Expand subvariables and environment variables here.
            ! This point is independent of a parameter bexpandVars
            ! as subfiles may otherwise not be found!
            call INIP_expandEnvVariable(sstring)
            call INIP_expandSubvariable(rparlist,sstring)

            p_Ssubfiles(nsubfiles+1) = sstring
            nsubfiles = nsubfiles + 1
          end if
        else
          ! Get all the filenames.
          allocate(p_SsubfilesTemp(nsubfiles+nnewsubfiles))
          p_SsubfilesTemp(1:nsubfiles) = p_Ssubfiles(:)
          deallocate(p_Ssubfiles)
          p_Ssubfiles => p_SsubfilesTemp
          nullify(p_SsubfilesTemp)

          do j=1,nnewsubfiles
            call INIP_getvalue_string(rparlist%p_Rsections(1), iidxsubfiles, &
              p_Ssubfiles(nsubfiles+j),isubstring=j,bdequote=.true.)

            ! Expand subvariables and environment variables here.
            ! This point is independent of a parameter bexpandVars
            ! as subfiles may otherwise not be found!
            call INIP_expandEnvVariable(p_Ssubfiles(nsubfiles+j))
            call INIP_expandSubvariable(rparlist,p_Ssubfiles(nsubfiles+j))
          end do

          nsubfiles = nsubfiles + nnewsubfiles
        end if

        ! Remove all substrings from the simportdatafiles parameter, so the parameter
        ! is filled with new data upon the next read statement.
        call INIP_addvalue(rparlist%p_Rsections(1), "SIMPORTDATAFILES", "")

      end if

    end do

    ! Ok, at that point we know which files to read -- so read them, one after
    ! the other. The 'master' file must be read at last!
    ilensubf = 1
    do icurrentsubfile = 2,nsubfiles
      sstring = trim(p_Ssubfiles(icurrentsubfile))
      ilensubf = max(ilensubf,len_trim(sstring))
      inquire(file=sstring, exist=bexists)

      if (bexists) then
        ! Read the sub-files...
        ! Do not yet expand variables.
        call INIP_readfromsinglefile (rparlist, sstring, .true., .false.)
      end if
    end do

    ! ... and the master
    icurrentsubfile = 1
    sstring = trim(p_Ssubfiles(icurrentsubfile))
    inquire(file=sstring, exist=bexists)

    if (bexists) then
      ! Do not yet expand variables.
      call INIP_readfromsinglefile (rparlist, sstring, .true., .false.)
    end if

    if (nsubfiles .gt. 1) then
      ! There have a couple of subfiles been read from disc.
      ! All subfiles can be found in p_Ssubfiles.

      ! Incorporate the complete list to the parameter "SIMPORTDATAFILES".
      if (associated(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList)) then
        deallocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList)
      end if

      allocate(rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(ilensubf+2,nsubfiles-1))
      rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(:,:) = ' '
      do icurrentsubfile = 1,nsubfiles-1
        call inip_stringToCharArray('"'//trim(p_Ssubfiles(1+icurrentsubfile))//'"',&
          rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%p_SentryList(:,icurrentsubfile));
      end do
      rparlist%p_Rsections(1)%p_Rvalues(iidxsubfiles)%nsize = nsubfiles-1
    end if

    ! Release memory, finish
    deallocate(p_Ssubfiles)

    ! Now expand all subvariables and environment variables to the actual values.
    if (.not. present(bexpandVars)) then
      call INIP_expandEnvVariables(rparlist)
      call INIP_expandSubvars(rparlist)
    else if (bexpandVars) then
      call INIP_expandEnvVariables(rparlist)
      call INIP_expandSubvars(rparlist)
    end if

  end subroutine

  ! ***************************************************************************

  subroutine INIP_readfromsinglefile (rparlist, sfilename, bimportSections, bexpandVars)

    !<description>

    ! This routine parses a text file for data of the INI-file form.
    ! sfilename must be the name of a file on the hard disc.
    ! The parameters read from the file are added to the parameter list
    ! rparlist, which has to be initialised with INIP_init before
    ! calling the routine.
    ! Remark: When adding parameters/sections to rparlist, the routine
    !   checks whether the parameters/sections already exist.
    !   Adding a parameter/section which exists does not result in an error -
    !   the first instance of the parameter/section will just be overwritten.

    !</description>

    !<inputoutput>

    ! The parameter list which is filled with data from the file
    type(t_parlist), intent(inout) :: rparlist

    !</inputoutput>

    !<input>

    ! The filename of the file to read.
    character(LEN=*), intent(in) :: sfilename

    ! TRUE: Import all sections in the DAT file.
    ! FALSE: Import only the main (unnamed) section and ignore all other
    ! sections.#
    logical, intent(in) :: bimportSections

    ! OPTIONAL: Expand references to subvariables.
    ! TRUE: Subvariables like "%{section.varname} are expanded to actual values.
    !       This is the standard setting.
    ! FALSE: Subvariables are left as they are.
    logical, intent(in), optional :: bexpandVars

    !</input>

    ! local variables
    integer :: iunit,ios,isbuflen,ityp,ilinenum,isubstring,nsubstrings,iparpos
    type(t_parlstSection), pointer :: p_currentsection
    character(LEN=INIP_LENLINEBUF) :: sdata
    character(LEN=INIP_MLSECTION) :: ssectionname
    character(LEN=INIP_MLNAME) :: sparname
    character(LEN=INIP_LENLINEBUF) :: svalue

    ! Try to open the file
    call inip_openFileForReading(sfilename, iunit)

    ! Oops...
    if (iunit .eq. -1) then
      call inip_output_line ('Error opening .INI file: '//trim(sfilename))
      call inip_sys_halt()
    end if

    ! Start adding parameters to the unnamed section
    p_currentsection => rparlist%p_Rsections(1)

    ! Read all lines from the file
    ios = 0
    ilinenum = 0
    isubstring = 0
    nsubstrings = 0
    do while (ios .eq. 0)

      ! Read a line from the file into sbuf
      call INIP_readlinefromfile (iunit, sdata, isbuflen, ios)
      ilinenum = ilinenum + 1

      if (isbuflen .ne. 0) then

        ! Parse the line
        call INIP_parseline (sdata, ityp, nsubstrings, ilinenum, ssectionname, &
          sparname, svalue,sfilename)

        select case (ityp)
          case (1)
            ! Stop parsing the file here if bimportSections tells us to do so.
            ! When the first section starts, the unnamed section is finished.
            if (.not. bimportSections) exit

            ! Check if the section exists; if not, create a new one.
            call INIP_querysection(rparlist, ssectionname, p_currentsection)

            if (.not. associated(p_currentsection)) then
              ! A new section name. Add a section, set the current section
              ! to the new one.
              call INIP_addsection (rparlist, ssectionname)
              p_currentsection => rparlist%p_Rsections(rparlist%isectionCount)
            end if

          case (2)
            ! A new parameter. Add it to the current section.
            call INIP_addvalue (p_currentsection, sparname, svalue)

          case (3)
            ! 'Headline' of a multi-valued parameter. Add the parameter with
            ! isubstring subvalues
            call INIP_addvalue (p_currentsection, sparname, svalue, nsubstrings)

            ! Fetch the parameter for later adding of subvalues.
            iparpos = INIP_queryvalue(p_currentsection, sparname)

            ! isubstring counts the current readed substring.
            ! Set it to 0, it will be increased up to nsubstrings in 'case 4'.
            isubstring = 0

          case (4)
            ! Increase number of current substring
            isubstring = isubstring + 1

            ! Sub-parameter of a multi-valued parameter. Add the value to
            ! the last parameter that was added in case 3.
            call INIP_setvalue_fetch (p_currentsection, iparpos, svalue, &
              isubstring=isubstring)

            ! Decrement the substring counter. If we reach 0, INIP_parseline
            ! continues to parse standard parameters.
            nsubstrings = nsubstrings - 1

            ! Other cases: comment.
        end select

      end if

    end do

    ! Close the file.
    close (iunit)

    if (.not. present(bexpandVars)) then
      ! Expand all subvariables and environment variables to the actual values.
      call INIP_expandEnvVariables(rparlist)
      call INIP_expandSubvars(rparlist)
    else if (bexpandVars) then
      ! Expand all subvariables and environment variables to the actual values.
      call INIP_expandEnvVariables(rparlist)
      call INIP_expandSubvars(rparlist)
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_getStringRepresentation (rparlist, p_sconfiguration)

    !<description>
    ! Creates a string representation of the given parameter list rparlist.
    ! p_sconfiguration will be created as "array[1..*] of char" on
    ! the heap containing this representation. The memory must be manually
    ! released by the caller using DEALLOCATE when finished using the string
    ! representation.
    !</description>

    !<input>
    ! The parameter list which is filled with data from the file
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<output>
    ! A pointer to a character array containing all lines of the parameter list.
    ! Points to NULL() if there is no data in the parameter list.
    ! Each line is terminated by NEWLINE.
    ! If there is data, a new pointer is allocated for this on the heap.
    ! The user must manually release the memory when finished using it.
    character, dimension(:), pointer :: p_sconfiguration
    !</output>

    !</subroutine>

    integer :: ilength,isection,ivalue,ientry,icount
    character, dimension(:), pointer :: p_sbuf
    character(len=INIP_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised')
      call inip_sys_halt()
    end if

    nullify(p_sbuf)

    ! Number of characters in the buffer
    ilength = 0

    ! Loop through all sections
    do isection = 1,rparlist%isectionCount

      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      if (isection .gt. 1) then
        ! Empty line before
        if (ilength .gt. 0) call appendString(p_sbuf,ilength,'')
        call appendString(p_sbuf,ilength,&
          '['//trim(rparlist%p_Rsections(isection)%ssectionName)//']')
      end if

      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Write "name=value"
          call inip_charArrayToString(&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call appendString(p_sbuf,ilength,&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="//trim(sbuf))
        else
          ! Write "name(icount)="
          call appendString(p_sbuf,ilength,&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(inip_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call inip_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call appendString(p_sbuf,ilength,trim(sbuf))
          end do
        end if

      end do ! ivalue

    end do ! isection

    ! Allocate a new character array with the correct size, copy p_sbuf to
    ! that ald release the old p_sbuf.
    ! Return NULL() if there is no data.
    nullify(p_sconfiguration)
    if (ilength .gt. 0) then
      allocate(p_sconfiguration(ilength))
      p_sconfiguration = p_sbuf(1:ilength)
    end if

    ! Release our temp buffer
    if (associated(p_sbuf)) deallocate(p_sbuf)

  contains

    ! Makes sure, the character buffer points to a character memory block of
    ! size nsize. If not, the block is reallocated to have that size.
    subroutine assumeBufSize(p_sconfig,nsize)

      character, dimension(:), pointer :: p_sconfig
      integer, intent(in) :: nsize

      character, dimension(:), pointer :: p_sconfignew

      if (.not. associated(p_sconfig)) then
        allocate(p_sconfig(nsize))
      else if (size(p_sconfig) .lt. nsize) then
        allocate(p_sconfignew(nsize))
        p_sconfignew(1:size(p_sconfig)) = p_sconfig
        deallocate(p_sconfig)
        p_sconfig => p_sconfignew
      end if

    end subroutine

    ! Appends sstring to the buffer p_sconfig, followed by a NEWLINE
    ! character. Reallocates memory if necessary.
    subroutine appendString(p_sconfig,iconfigLength,sstring)

      ! Pointer to character data
      character, dimension(:), pointer :: p_sconfig

      ! In: Current length of data stream in p_sconfig.
      ! Out: New length of data stream in p_sconfig
      integer, intent(inout) :: iconfigLength

      ! The string to be added.
      character(LEN=*), intent(in) :: sstring

      integer :: nblocks,nblocksneeded,i

      ! How many memory blocks do we need for the current configuration?
      ! We work block-wise to prevent too often reallocation.
      if (.not. associated(p_sconfig)) then
        nblocks = 0
      else
        nblocks = size(p_sconfig) / INIP_STRLEN
      end if
      nblocksneeded = 1 + (iconfigLength+len(sstring)+1) / INIP_STRLEN
      if (nblocksneeded .gt. nblocks) then
        call assumeBufSize(p_sconfig,nblocksneeded*INIP_STRLEN)
      end if

      ! Append the data
      do i=1,len(sstring)
        iconfigLength = iconfigLength+1
        p_sconfig(iconfigLength) = sstring(i:i)
      end do

      ! Append NEWLINE as line-end character
      iconfigLength = iconfigLength+1
      p_sconfig(iconfigLength) = NEWLINE

    end subroutine

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_info (rparlist)

    !<description>
    ! Prints the parameter list rparlist to the terminal.
    !</description>

    !<input>
    ! The parameter list which is to be printed to the terminal.
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !</subroutine>


    integer :: isection,ivalue,ientry,icount
    character(len=INIP_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised')
      call inip_sys_halt()
    end if

    ! Loop through all sections
    do isection = 1,rparlist%isectionCount

      ! Append the section name. May be empty for the unnamed section,
      ! which is always the first one.
      if (isection .gt. 1) then
        ! Empty line before
        call inip_output_lbrk()
        call inip_output_line('['//trim(rparlist%p_Rsections(isection)%ssectionName)//']')
      end if

      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Write "name=value"
          call inip_charArrayToString(&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call inip_output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"="//trim(sbuf))
        else
          ! Write "name(icount)="
          call inip_output_line(&
            trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
            //"("//trim(inip_siL(icount, 10))//")=")
          ! Write all the entries of that value, one each line.
          do ientry = 1,icount
            call inip_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
              p_SentryList(:,ientry),sbuf)
            call inip_output_line(trim(sbuf))
          end do
        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

  !<function>

  integer function INIP_findvalue_indir (rsection, sparameter, svalue) &
      result (isubstring)

    !<description>
    ! Checks whether the parameter sparameter in the section rsection
    ! has the given value svalue.
    !</description>

    !<result>
    ! The index of the substring in the parameter sparameter which has value
    ! svalue or =-1, if the parameter does not exist within the section.
    !</result>

    !<input>

    ! The section where to search for the parameter
    type(t_parlstSection), intent(in) :: rsection

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    ! The value to search for
    character(LEN=*), intent(in) :: svalue

    !</input>

    !</function>

    ! local variables
    integer :: idx
    character(LEN=INIP_MLNAME) :: paramname
    character(len(svalue)) :: sbuf

    if (sparameter .eq. '') then
      call inip_output_line ('Empty parameter name!')
      call inip_sys_halt()
    end if

    ! Create the upper-case parameter name
    paramname = adjustl(sparameter)
    call inip_toupper (paramname)

    ! Get the parameter index into 'idx', finish.
    call INIP_fetchparameter(rsection, paramname, idx)

    ! Check if value svalue exists in some substring and return its
    ! index; of the value does not exist return -1
    if (idx .eq. 0) then
      call inip_charArrayToString(&
        rsection%p_Rvalues(idx)%p_sentry, sbuf)
      if (trim(sbuf) .eq. trim(svalue)) then
        isubstring = 0
      else
        isubstring = -1
      end if
    else
      do isubstring = 0, rsection%p_Rvalues(idx)%nsize
        call inip_charArrayToString(&
          rsection%p_Rvalues(idx)%p_SentryList(:,isubstring), sbuf)
        if (trim(sbuf) .eq. trim(svalue)) then
          return
        end if
      end do

      ! We did not find the desired value
      isubstring = -1
    end if

  end function

  ! ***************************************************************************

  !<function>

  integer function INIP_findvalue_direct (rparlist, ssectionName, sparameter, svalue) &
      result (isubstring)

    !<description>
    ! Checks whether the parameter sparameter in the section ssectionname
    ! in the parameter list rparlist has the given value.
    !</description>

    !<result>
    ! The index of the substring in the parameter sparameter which has value
    ! svalue or =-1, if the parameter does not exist within the section.
    !</result>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The section name - '' identifies the unnamed section.
    character(LEN=*), intent(in) :: ssectionName

    ! The parameter name to search for.
    character(LEN=*), intent(in) :: sparameter

    ! The value to search for
    character(LEN=*), intent(in) :: svalue

    !</input>

    !</function>

    ! local variables
    integer :: idx
    type(t_parlstSection), pointer :: p_rsection
    character(len(svalue)) :: sbuf

    ! Cancel if the list is not initialised.
    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised!')
      call inip_sys_halt()
    end if

    ! Get the section
    call INIP_querysection(rparlist, ssectionName, p_rsection)
    if (.not. associated(p_rsection)) then
      call inip_output_line ('Section not found: '//trim(ssectionName))
      return
    end if

    ! Get the parameter index
    idx = INIP_queryvalue_indir (p_rsection, sparameter)

    ! Check if value svalue exists in some substring and return its
    ! index; of the value does not exist return -1
    if (idx .eq. 0) then
      call inip_charArrayToString(&
        p_rsection%p_Rvalues(idx)%p_sentry, sbuf)
      if (trim(sbuf) .eq. trim(svalue)) then
        isubstring = 0
      else
        isubstring = -1
      end if
    else
      do isubstring = 0, p_rsection%p_Rvalues(idx)%nsize
        call inip_charArrayToString(&
          p_rsection%p_Rvalues(idx)%p_SentryList(:,isubstring), sbuf)
        if (trim(sbuf) .eq. trim(svalue)) then
          return
        end if
      end do

      ! We did not find the desired value
      isubstring = -1
    end if

  end function

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_expandEnvVariables(rparlist)

    !<description>
    ! This subroutine expands variables when referring to subvariables:
    ! The value of a parameter may refer to another variable in the parameter
    ! list. This can be specified by tokens of the form "%{NAME}" or "{SECTION.NAME}",
    ! depending on whether it is part of the main section or not.
    ! Example: "NLMIN = %{NLMAX}"; in this case, NLMIN is always set to
    ! the same value as NLMAX.
    !
    ! The routine parses all variables in the DAT file to resolve such
    ! references. Note that recursive definitions are not allowed!
    !</description>

    !<inputoutput>
    ! The parameter list which is filled with data from the file
    type(t_parlist), intent(inout) :: rparlist
    !</inputoutput>

    !</subroutine>

    integer :: isection,ivalue,ientry,icount,j
    character(len=INIP_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised')
      call inip_sys_halt()
    end if

    ! Loop through all sections
    do isection = 1,rparlist%isectionCount

      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
        if (icount .eq. 0) then
          ! Expand the value if is refers to subvalues.
          call inip_charArrayToString(&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call INIP_expandEnvVariable(sbuf)

          j = len_trim(sbuf)
          deallocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry)
          allocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(j))
          rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(:) = ' '
          call inip_stringToCharArray(sbuf,&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,j)
        else
          ! Loop through the subvalues.
          do ientry = 1,icount
            ! Expand the value if is refers to subvalues.
            call inip_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call INIP_expandEnvVariable(sbuf)

            ! Reallocate before writing back if necessary
            j = len_trim(sbuf)
            if (j .gt. ubound(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList,1)) then
              call INIP_reallocSubVariables(rparlist%p_Rsections(isection)%p_Rvalues(ivalue),j)
            end if
            call inip_stringToCharArray(sbuf,&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),j)
          end do
        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_expandEnvVariable(sbuffer)

    !<description>
    ! This subroutine recursively expands all environment variables in the given
    ! string sbuffer.
    !</description>

    !<input>
    ! string; all environment variables in here are replaced
    character(len=*), intent(inout) :: sbuffer
    !</input>
    !</subroutine>

    ! flag
    logical :: bfoundInEnv

    ! start and end position of variable
    integer :: istartPos, istopPosRelative

    ! variable to expand environment variable on-the-fly to if found
    character(len=INIP_STRLEN) :: sauxEnv

    ! Buffer for the result
    character(len=INIP_LENLINEBUF) :: sresult

    ! Initialise return value
    sresult = trim(sbuffer)

    ! check for a $ character
    istartPos = index(sresult, "$")
    do while (istartPos .gt. 0)
      ! Detect end of variable: a variable ends at the first character that
      ! is neither in '[A-Z]', '[0-9]' nor '_'.
      istopPosRelative = verify(sresult(istartPos+1:), &
        "abcdefghijklmnopqrstuvwxyz" // &
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ" // &
        "0123456789_")

      bfoundInEnv = .false.
      ! Retrieve value of environment variable
      ! (Do not forget to cut the dollar sign.)
      bfoundInEnv = &
        inip_getenv_string(trim(&
        sresult(istartPos + 1 : istartPos + istopPosRelative - 1)), sauxEnv)
      if (bfoundInEnv) then
        ! Replace environment variable by its content
        sresult = sresult(1:istartPos-1) // &
          trim(sauxEnv) // &
          trim(sresult(istartPos + istopPosRelative:))
      else
        call inip_output_line ('Environment variable <'//&
          trim(sresult(istartPos + 1 : istartPos + istopPosRelative - 1))//&
          '> not found!')
        call inip_sys_halt()
      endif

      ! check for next $ character
      istartPos = index(sresult, "$")
    enddo

    ! Replace by the result
    sbuffer = sresult

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_expandSubvars(rparlist)

    !<description>
    ! This subroutine expands variables when referring to subvariables:
    ! The value of a parameter may refer to another variable in the parameter
    ! list. This can be specified by tokens of the form "%{NAME}" or "{SECTION.NAME}",
    ! depending on whether it is part of the main section or not.
    ! Example: "NLMIN = %{NLMAX}"; in this case, NLMIN is always set to
    ! the same value as NLMAX.
    !
    ! The routine parses all variables in the DAT file to resolve such
    ! references. Note that recursive definitions are not allowed!
    !</description>

    !<inputoutput>
    ! The parameter list which is filled with data from the file
    type(t_parlist), intent(inout) :: rparlist
    !</inputoutput>

    !</subroutine>

    integer :: isection,ivalue,ientry,icount,j
    character(len=INIP_LENLINEBUF) :: sbuf

    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised')
      call inip_sys_halt()
    end if

    ! Loop through all sections
    do isection = 1,rparlist%isectionCount

      ! Loop through the values in the section
      do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

        ! Do we have one or multiple entries to that parameter?
        icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize

        if (icount .eq. 0) then
          ! Expand the value if is refers to subvalues.
          call inip_charArrayToString(&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
          call INIP_expandSubvariable(rparlist,sbuf)

          j = len_trim(sbuf)
          deallocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry)
          allocate(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(j))
          rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry(:) = ' '
          call inip_stringToCharArray(sbuf,&
            rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,j)
        else
          ! Loop through the subvalues.
          do ientry = 1,icount
            ! Expand the value if is refers to subvalues.
            call inip_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),sbuf)
            call INIP_expandSubvariable(rparlist,sbuf)

            ! Reallocate before writing back if necessary
            j = len_trim(sbuf)
            if (j .gt. ubound(rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList,1)) then
              call INIP_reallocSubVariables(rparlist%p_Rsections(isection)%p_Rvalues(ivalue),j)
            end if
            call inip_stringToCharArray(sbuf,&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_SentryList(:,ientry),j)
          end do

        end if

      end do ! ivalue

    end do ! isection

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_findSubvariable(sstring,istart,iend,ssection,sname,ivalue)

    !<description>
    ! Searchs in the string sstring for the first occurance of a subvariable.
    ! A subvariable has the form "%{NAME}" or "{SECTION.NAME}"
    ! or "%{NAME:INDEX}" or "%{SECTION.NAME:INDEX}"
    ! depending on whether it is part of the main section or not.
    !</description>

    !<input>
    ! The string where a subvariable is searched.
    character(len=*), intent(in) :: sstring
    !</input>

    !<output>
    ! Returns the start of the variable in the string or 0 if no subvariable
    ! is found.
    integer, intent(out) :: istart

    ! Returns the end of the variable in the string or 0 if no subvariable
    ! is found.
    integer, intent(out) :: iend

    ! Returns the name of the section or "" if either no subvariable is found
    ! or the unnamed section is referred to.
    character(len=*), intent(out) :: ssection

    ! Returns the name of the subvariable or "" if no subvariable is found.
    character(len=*), intent(out) :: sname

    ! Returns the number/index INDEX of the subvalue or 0, if there is
    ! no index or if no subvariable is found.
    integer, intent(out) :: ivalue
    !</output>

    !</subroutine>

    ! local variables
    integer :: i,j,istrlen,idotpos,icolonpos
    logical :: bstartfound

    ! Ok, this is a parser. We have the following rules:
    ! "%%" means one "%" and is not interpreted as the begin of a
    ! token.
    ! "%{NAME}" is a token referring to a variable in the unnamed section.
    ! "%{NAME:INDEX}" is a token referring to a subvariable
    ! of variable NAME in section SECTION.
    ! "%{SECTION.NAME}" is a token referring to a name in a named section
    ! which must exist.
    ! "%{SECTION.NAME:INDEX}" is a token referring to a subvariable
    ! of variable NAME in section SECTION.
    !
    istart = 0
    iend = 0
    ivalue = 0
    bstartfound = .false.

    ! Lets loop through the characters.
    istrlen = len(sstring)
    do i=1,istrlen
      if (sstring(i:i) .eq. "%") then
        ! Did we already found the "%"?
        if (bstartfound) then
          ! That is our escape sequence. Do not do anything, just
          ! return to 'normal' mode.
          bstartfound = .false.
        else
          ! That is probably the beginning of a token.
          bstartfound = .true.
        end if

        ! Go on.
        cycle
      end if

      ! The next things only execute if we are close to a token...
      if (bstartfound) then
        if (sstring(i:I) .eq. "{") then
          ! Yes, that is a token.
          istart = i-1

          ! Find the end of the token and probably the dot/colon
          idotpos = 0
          icolonpos = 0
          do j=istart+1,istrlen
            if (sstring(j:J) .eq. ".") then
              ! Here is the dot.
              idotpos = j
            end if

            if (sstring(j:J) .eq. ":") then
              ! Here is the dot.
              icolonpos = j
            end if

            if (sstring(j:j) .eq. "}") then
              ! Here, the token ends.
              iend = j

              ! Extract name and probably the section
              if (idotpos .eq. 0) then
                ssection = ""
                if (icolonpos .eq. 0) then
                  sname = sstring(istart+2:iend-1)
                else
                  sname = sstring(istart+2:icolonpos-1)

                  ! Get the value
                  read(sstring(icolonpos+1:iend-1),*) ivalue
                end if
              else
                ssection = sstring(istart+2:idotpos-1)
                if (icolonpos .eq. 0) then
                  sname = sstring(idotpos+1:iend-1)
                else
                  sname = sstring(idotpos+1:icolonpos-1)

                  ! Get the value
                  read(sstring(icolonpos+1:iend-1),*) ivalue
                end if
              end if

              ! That is it.
              return

            end if
          end do
        end if
      end if
    end do

    ! Nothing found.
    ssection = ""
    sname = ""

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_expandSubvariable(rparlist,sstring)

    !<description>
    ! Expands the subvariables in sstring to a fully qualified string.
    !</description>

    !<input>
    ! The parameter list containing the variables that can be used
    ! as subvariables.
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    ! The string to be expanded. Receives the expanded string upon return.
    character(len=*), intent(inout) :: sstring
    !</inputoutput>

    !</subroutine>

    ! local variables
    integer :: istart,iend,ivalue
    character(len=INIP_MLSECTION) :: ssection
    character(len=INIP_MLNAME) :: sname
    character(len=len(sstring)) :: sbuffer
    character(len=INIP_LENLINEBUF) :: sdata

    ! Repeat until we found all subvariables
    istart = 1
    iend = 0
    do while (istart .ne. 0)

      ! Find the first subvariable
      call INIP_findSubvariable(sstring,istart,iend,ssection,sname,ivalue)

      if (istart .ne. 0) then
        ! Copy the string to the buffer
        sbuffer = sstring

        ! Now copy back and replace the variable by the stuff from the
        ! parameter list.
        if (ivalue .eq. 0) then
          call INIP_getvalue_string (rparlist, ssection, sname, sdata, bdequote=.true.)
        else
          call INIP_getvalue_string (rparlist, ssection, sname, sdata, &
            isubstring=ivalue, bdequote=.true.)
        end if
        sstring = sbuffer(1:istart-1)//trim(sdata)//sbuffer(iend+1:)

      end if

    end do

  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine INIP_dumpToFile(rparlist, sfilename, cflag)

    !<description>
    ! This subroutine dumps a given parameter list into a text file in
    ! INI-file form. Note that no additional comments are exported but
    ! just the tuples `parameter = value` in the corresponding sections.
    !</description>

    !<input>

    ! The parameter list.
    type(t_parlist), intent(in) :: rparlist

    ! The name of the output file
    character(*), intent(in) :: sfilename

    ! mode: INIP_APPEND or INIP_REPLACE
    integer, intent(in) :: cflag

    !</input>
    !</subroutine>

    ! local variables
    character(len=INIP_LENLINEBUF) :: sbuf
    integer :: iunit,isection,ivalue,ientry,icount

    if (rparlist%isectionCount .eq. 0) then
      call inip_output_line ('Parameter list not initialised')
      call inip_sys_halt()
    end if

    if (inip_id .eq. inip_showid) then

      ! Open file for output
      call inip_openFileForWriting(sfilename, iunit, cflag, bformatted=.true.)
      if (iunit .eq. -1) then
        call inip_output_line ('Unable to open file for output!')
        call inip_sys_halt()
      end if

      ! Loop through all sections
      do isection = 1,rparlist%isectionCount

        ! Append the section name. May be empty for the unnamed section,
        ! which is always the first one.
        if (isection .gt. 1) then
          ! Empty line before
          write(iunit,*)
          write(iunit,*) '['//trim(rparlist%p_Rsections(isection)%ssectionName)//']'
        end if

        ! Loop through the values in the section
        do ivalue = 1,rparlist%p_Rsections(isection)%iparamCount

          ! Do we have one or multiple entries to that parameter?
          icount = rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%nsize
          if (icount .eq. 0) then
            ! Write "name=value"
            call inip_charArrayToString(&
              rparlist%p_Rsections(isection)%p_Rvalues(ivalue)%p_sentry,sbuf)
            write(iunit,*)&
              trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
              //" = "//trim(sbuf)
          else
            ! Write "name(icount)="
            write(iunit,*)&
              trim(rparlist%p_Rsections(isection)%p_Sparameters(ivalue)) &
              //"("//trim(inip_siL(icount, 10))//") = "
            ! Write all the entries of that value, one each line.
            do ientry = 1,icount
              call inip_charArrayToString(&
                rparlist%p_Rsections(isection)%p_Rvalues(ivalue)% &
                p_SentryList(:,ientry),sbuf)
              write(iunit,*) trim(sbuf)
            end do
          end if

        end do ! ivalue

      end do ! isection

      ! Close file
      close(iunit)

    end if

  end subroutine INIP_dumpToFile

  !<function>

  function inip_siL(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with white spaces.
    !</description>

    !<input>
    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
    !</input>

    !<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=idigits) :: soutput
    !</result>

    !</function>

    soutput = adjustl(inip_si(ivalue, idigits))
  end function inip_siL

  !<function>

  function inip_stringToSingle(svalue,sformat) result(fvalue)

    !<description>
    ! This routine converts a given string that provides a valid
    ! IEEE 745 representation of a real number into a single value.
    !</description>

    !<input>

    ! string containing the real number
    character(LEN=*), intent(in) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(in) :: sformat

    !</input>

    !<result>
    ! single precision value
    real*4 :: fvalue
    !</result>
    !</function>

    ! local variables
    character(LEN=len(svalue)+3) :: svalueTemp
    integer :: ipos

    ! Check if string contains a 'dot'
    if (scan(svalue,".") .ne. 0) then

      ! Read original string
      read(svalue,sformat) fvalue
    else
      ! Check if string is given in scientific notation
      ipos = scan(svalue,"dDeE")

      if (ipos .eq. 0) then
        ! Append '.E0' to convert string into scientific notation
        svalueTemp = trim(svalue)//".E0"

      elseif (ipos .eq. 1) then
        ! Prepend '1.' to convert string into scientific notation
        svalueTemp = "1."//adjustl(svalue)
      else
        ! Insert '.' to convert string into scientific notation
        svalueTemp = svalue(1:ipos-1)//"."//svalue(ipos:)
      end if

      ! Read modified string
      read(svalueTemp,sformat) fvalue
    end if

  end function inip_stringToSingle

  !<function>

  function inip_stringToDouble(svalue,sformat) result(dvalue)

    !<description>
    ! This routine converts a given string that provides a valid
    ! IEEE 745 representation of a real number into a double value.
    !</description>

    !<input>

    ! string containing the real number
    character(LEN=*), intent(in) :: svalue

    ! format description to use for conversion
    character(LEN=*), intent(in) :: sformat

    !</input>

    !<result>
    ! double precision value
    real*8 :: dvalue
    !</result>
    !</function>

    ! local variables
    character(LEN=len(svalue)+3) :: svalueTemp
    integer :: ipos

    ! Check if string contains a 'dot'
    if (scan(svalue,".") .ne. 0) then

      ! Read original string
      read(svalue,sformat) dvalue
    else
      ! Check if string is given in scientific notation
      ipos = scan(svalue,"dDeE")

      if (ipos .eq. 0) then
        ! Append '.E0' to convert string into scientific notation
        svalueTemp = trim(svalue)//".E0"

      elseif (ipos .eq. 1) then
        ! Prepend '1.' to convert string into scientific notation
        svalueTemp = "1."//adjustl(svalue)
      else
        ! Insert '.' to convert string into scientific notation
        svalueTemp = svalue(1:ipos-1)//"."//svalue(ipos:)
      end if

      ! Read modified string
      read(svalueTemp,sformat) dvalue
    end if

  end function inip_stringToDouble

  !<function>

  character (len=32) function inip_si(ivalue, idigits) result(soutput)

    !<description>
    ! This routine converts an integer value to a string of length idigits.
    !</description>

    !<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
    !</result>

    !<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
    !</input>
    !</function>

    character (len=16) :: sformat
    character (len=2)  :: saux
    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   inip_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
      call inip_output_line("*** WARNING! Too many decimal places requested in inip_si! ***")
      write(saux, "(i2)") 16
    else if (idigits .lt. 10) then
      write(saux, "(i1)") idigits
    else
      write(saux, "(i2)") idigits
    endif

    sformat = "(i" // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue

  end function inip_si

  !<function>
  logical function inip_getenv_string(svar, sresult)

    !<description>
    ! This functions returns the string value of a given enviroment variable. The routine
    ! returns .TRUE., if the variable exists, otherwise .FALSE. .
    !</description>

    !<input>
    ! Name of the enviroment variable
    character(len=*), intent(in) :: svar
    !</input>

    !<output>
    ! Value of the enviroment variable
    character(len=*), intent(out) :: sresult
    !</output>

    !<result>
    ! exit status
    !</result>
    !</function>

    character(len=max(INIP_STRLEN,len(sresult))) :: svalueInEnv

    integer :: nstatus

    call get_environment_variable(trim(svar), svalueInEnv, status=nstatus)

    select case (nstatus)
      case (0)
        ! Copy string only up to first whitespace character
        !      read(svalueInEnv, '(A)') sresult

        ! Copy complete string
        sresult = svalueInEnv
        inip_getenv_string = .true.

      case (1)
        ! Environment variable does not exist
        sresult = ""
        inip_getenv_string = .false.

      case default
        !  2: Processor does not support environment variables
        ! >2: Some error occurred
        ! -1: variable svalueInEnv too short to absorb environment variables` content
        sresult = ""
        inip_getenv_string = .false.

    end select

  end function inip_getenv_string

end module
