program parserdemonstration
    use iniparser

    implicit none

    integer :: someInt, someOtherInt
    real*8 :: someReal8
    real*4 :: someReal4
    real*8 :: someOtherReal, defaultReal
    character(len=INIP_STRLEN) somestring1, somestring2
    character(len=INIP_STRLEN) somestring3, somestring4
    logical :: somebool
    integer :: iunit

    ! Usually the application that uses the parser is parallel.
    ! Only one of the threads should write out, the one that
    ! has the id "showid"
    ! Both values are always set in the applications, the just
    ! need to be set for the parser too as it is a standalone module.
    integer :: myid = 1
    integer :: showid = 1
    ! Define two units (for the protocolfile and for the terminal)
    ! Every call to inip_output_line (an internal routine in the
    ! parser) will be written out to both units.
    ! If a unit is set to -1 then nothing is written out to that.
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    call inip_output_init(myid,showid,unitProtfile,unitTerminal)

    ! Init the parameterlist
    call inip_init(parameterlist)

    ! Now compare to the file testcfg.dat
    call inip_readfromfile(parameterlist,'testcfg.dat')
    write(*,*) "Hello. I read in something"

    ! Demonstration how to get values
    call inip_getvalue_int(parameterlist,"SECTION1","VAR",someInt)
    ! This section and variable do not exist. Therefore, take the default value
    ! which is set to 1
    call INIP_getvalue_int(parameterlist,"SECTION2","SomeDemo",someOtherInt,1)
    call INIP_getvalue_single(parameterlist,"SECTION1","VAR2", someReal4)
    call INIP_getvalue_double(parameterlist,"SECTION1","VAR2", someReal8)
    ! We can set a default that is taken if the value is not found
    defaultReal = 20.5
    call INIP_getvalue_double(parameterlist,"SECTION1","VAR100",someOtherReal,defaultReal)
    call INIP_getvalue_string(parameterlist,"SECTION1","VAR4",somestring1,bdequote=.TRUE.)
    call INIP_getvalue_string(parameterlist,"SECTION1","VAR4",somestring2,bdequote=.FALSE.)

    call INIP_getvalue_string(parameterlist,"SECTION1","VAR5",somestring3,bdequote=.TRUE.)
    call INIP_getvalue_string(parameterlist,"SECTION1","VAR5",somestring4,bdequote=.FALSE.)

    call INIP_getvalue_logical(parameterlist,"TRALALA","myfunnyname",somebool)


    ! Control output so you know what the values are

    write(*,*) "Control output:"
    write(*,*) someInt
    write(*,*) someOtherInt
    write(*,*) someReal4
    write(*,*) someReal8
    write(*,*) someOtherReal
    write(*,*) somestring1
    write(*,*) somestring2
    write(*,*) somestring3
    write(*,*) somestring4
    write(*,*) somebool
    write(*,*) "End Control output"

    ! We can modify values in the parameterlist
    call INIP_setvalue(parameterlist,"TRALALA","FOO","BAR")

    ! We can even set new values in the parameterlist
    ! Only condition is: The parameter value is presented as a string
    call INIP_addvalue(parameterlist,"tralala","setFromProgram",'2.5')

    ! And even create new sections
    call INIP_addsection(parameterlist,"SomeWeirdSection")
    call INIP_addvalue(parameterlist,'SomeWeirdSection','something','anything')

    ! Make some output on the terminal and write the content of the
    ! parameterlist to the screen
    ! Only works if the output system is initialised
    call inip_info(parameterlist)

    ! Write it into a file
    call inip_dumpToFile(parameterlist,"demo_dump_parlst.dat",INIP_REPLACE)

    ! We can also indent the sections one level in some sense
    ! (Some folks consider B a subsection of A if the sectionname
    ! is A/B )
    ! All but one section
    call INIP_indentAllSectionsButOne(parameterlist,"SECTION1","FORALLBUTONE")
    ! Only one section:
    call INIP_indentSection(parameterlist,"SECTION1","myIndentionForSectionOne")
    ! OR all sections
    call INIP_indentAllSections(parameterlist,"myIndentionForALL")
    call inip_dumpToFile(parameterlist,"demo_dump_parlst_indented.dat",INIP_REPLACE)

    ! We can also dump to a unit
    ! First we need to open the unit
    call inip_openFileForWriting("someFile.dat", iunit, INIP_REPLACE)
    call inip_dumpToUnit(parameterlist,iunit)

    ! Clean up the parameterlist
    call inip_done(parameterlist)

end program
