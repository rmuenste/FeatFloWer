SUBROUTINE LoadMPIDumpFiles(iOutO,cList)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,master,coarse,MPI_SEEK_SET,MPI_REAL8,MPI_MODE_RDONLY,MPI_MAX
USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL,mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind,mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
! USE PP3D_MPI
USE var_QuadScalar,ONLY:QuadSc,LinSc,myDump,mg_mesh,Screw,MaterialDistribution,temperature,mySegmentIndicator
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:GenLinScalar,bParallel

IMPLICIT NONE
character(len=*),intent(in) :: cLIST
INTEGER iOutO
!----------------------------------------------------------------------
INTEGER  :: mpiFile,ierr
CHARACTER cPOutFile*256,cDumpFolder*256
INTEGER iOut,myCoarseElem,nLengthE,nLengthV
INTEGER iP,lCmnd
INTEGER ifilen,i
DATA ifilen/0/
real*8, allocatable :: ElementOffsets(:)
integer iFld,jFld,nFLD
character*(1),allocatable :: cFLD(:)
INTEGER iGlobalError,jGlobalError
INTEGER :: dblesize=8,intsize=4

IF (myid.ne.0) CALL CreateDumpStructures(1)

nFLD = (LEN(cLIST)+1)/2
allocate(cFLD(nFLD))
read(cLIST,*,IOSTAT=iERR) cFLD

allocate(ElementOffsets(subnodes+1))
ElementOffsets = 0d0

IF (bParallel) THEN
 if (myid.ne.0) ElementOffsets(myid+1) = DBLE(KNEL(NLMIN))
 CALL Comm_SummN(ElementOffsets,subnodes+1)
 do i=2,subnodes+1
  ElementOffsets(i) = ElementOffsets(i) + ElementOffsets(i-1)
 end do
ELSE

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

iOut = iOutO

IF (myid.NE.0) THEN

 NLMAX = NLMAX+1

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-2)
 nLengthV = (2**(NLMAX-1)+1)**3
 
 iGlobalError = 0
 DO jFld = 1,nFLD
  if (cFLD(jFld).eq.'p'.or.cFLD(jFld).eq.'P') CALL LoadMPIFieldP1('pressure',LinSc%valP(NLMAX-1)%x)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','pressure')
  
  if (cFLD(jFld).eq.'m'.or.cFLD(jFld).eq.'M') CALL LoadMPIFieldP0_INT('material',MaterialDistribution(NLMAX-1)%x)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','material')
  
  if (cFLD(jFld).eq.'x'.or.cFLD(jFld).eq.'X') THEN
   if (bParallel) then
    CALL LoadMPIFieldQ2_X3('coordinates',mg_mesh%level(NLMAX)%dcorvg)
   else
    CALL ReLoadMPIFieldQ2_X3('coordinates',mg_mesh%level(NLMAX)%dcorvg)
   end if
  end if
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','coordinates')
  
  if (cFLD(jFld).eq.'v'.or.cFLD(jFld).eq.'V') THEN
   if (bParallel) then
    CALL LoadMPIFieldQ2_NX('velocity',3,QuadSc%valU,QuadSc%valV,QuadSc%valW)
   else
    CALL ReLoadMPIFieldQ2_NX('velocity',3,QuadSc%valU,QuadSc%valV,QuadSc%valW)
   end if
  end if
  
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','velocity')

  if (cFLD(jFld).eq.'d'.or.cFLD(jFld).eq.'D') CALL LoadMPIFieldQ2_NX('distance',1,Screw)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','distance')
  
  if (cFLD(jFld).eq.'s'.or.cFLD(jFld).eq.'S') CALL LoadMPIFieldQ2_NX('segment',1,mySegmentIndicator(2,:))
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','segment')

  if (cFLD(jFld).eq.'t'.or.cFLD(jFld).eq.'T') CALL LoadMPIFieldQ2_NX('temperature',1,temperature)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','temperature')
  
  if (cFLD(jFld).eq.'q'.or.cFLD(jFld).eq.'Q') THEN
   if (allocated(GenLinScalar%Fld)) then
    DO iFld=1,GenLinScalar%nOfFields
     CALL LoadMPIFieldQ2_NX(adjustl(trim(GenLinScalar%prm%cField(iFld))),1,GenLinScalar%fld(iFld)%Val)
     CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
     if (jGlobalError.ne.0) CALL ProcessError('R',adjustl(trim(GenLinScalar%prm%cField(iFld))))
    END DO     
   END IF
  END IF
 END DO
 
!  DO jFld = 1,nFLD
!   if (cFLD(jFld).eq.'p'.or.cFLD(jFld).eq.'P') CALL LoadMPIFieldP1('pressure',LinSc%valP(NLMAX-1)%x)
!   if (cFLD(jFld).eq.'m'.or.cFLD(jFld).eq.'M') CALL LoadMPIFieldP0_INT('material',MaterialDistribution(NLMAX-1)%x)
!   if (cFLD(jFld).eq.'x'.or.cFLD(jFld).eq.'X') CALL LoadMPIFieldQ2_X3('coordinates',mg_mesh%level(NLMAX)%dcorvg)
!   if (cFLD(jFld).eq.'v'.or.cFLD(jFld).eq.'V') CALL LoadMPIFieldQ2_NX('velocity',3,QuadSc%valU,QuadSc%valV,QuadSc%valW)
!   if (cFLD(jFld).eq.'d'.or.cFLD(jFld).eq.'D') CALL LoadMPIFieldQ2_NX('distance',1,Screw)
!   if (cFLD(jFld).eq.'s'.or.cFLD(jFld).eq.'S') CALL LoadMPIFieldQ2_NX('segment',1,mySegmentIndicator(2,:))
!   if (cFLD(jFld).eq.'t'.or.cFLD(jFld).eq.'T') CALL LoadMPIFieldQ2_NX('temperature',1,temperature)
!   if (cFLD(jFld).eq.'q'.or.cFLD(jFld).eq.'Q') THEN
!    if (allocated(GenLinScalar%Fld)) then
!     DO iFld=1,GenLinScalar%nOfFields
!      CALL LoadMPIFieldQ2_NX(adjustl(trim(GenLinScalar%prm%cField(iFld))),1,GenLinScalar%fld(iFld)%Val)
!     END DO     
!    END IF
!   END IF
!  END DO
 
 NLMAX = NLMAX-1
  
END IF
 
!!! Exchange trhe coordinates with the master !!!!
IF (bParallel) THEN
 CALL CommCoordinatesWithMaster()
end if
IF (myid.ne.0) CALL CreateDumpStructures(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(ElementOffsets)

 CONTAINS
 
 SUBROUTINE ReloadMPIFieldP1(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=4
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nF*nLengthE*knel(nlmin)
  allocate(daux(ndof))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)
  
  myFieldOffset = offset + intsize*INT(ElementOffsets(subnodes+1)) + dblesize*INT(ElementOffsets(myid))*nLengthE*nF
  
  CALL MPI_File_read_at(mpiFile,myFieldOffset,daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  CALL mpi_file_close(mpiFile,ierr)
  
  iEntry = 0
  DO iel=1,knel(nlmin)
   DO jel=1,nLengthE
    kel = myDump%Elements(IEL,jel)
    iP = nF*(kel-1)
    Field(iP+1:iP+nF) = daux(iEntry+1:iEntry+nF) 
    iEntry = iEntry + nF
   END DO
  end do
 
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux)

 END SUBROUTINE ReloadMPIFieldP1

 SUBROUTINE ReLoadMPIFieldQ2_NX(cF,nF,Field1,Field2,Field3,Field4,Field5,Field6,Field7,Field8,Field9,Field10)
 character cF*(*)
 integer :: nF
 REAL*8 Field1(*)
 REAL*8, OPTIONAL :: Field2(*)
 REAL*8, OPTIONAL :: Field3(*)
 REAL*8, OPTIONAL :: Field4(*)
 REAL*8, OPTIONAL :: Field5(*)
 REAL*8, OPTIONAL :: Field6(*)
 REAL*8, OPTIONAL :: Field7(*)
 REAL*8, OPTIONAL :: Field8(*)
 REAL*8, OPTIONAL :: Field9(*)
 REAL*8, OPTIONAL :: Field10(*)
 
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nF*nLengthV*knel(nlmin)
  allocate(iaux(knel(nlmin)))
  allocate(daux(ndof))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  myFieldOffset = offset
  CALL MPI_File_read_at(mpiFile, myFieldOffset,iaux, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)

  myFieldOffset = offset + intsize*knel(nlmin)
  CALL MPI_File_read_at(mpiFile, myFieldOffset,daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  
  CALL mpi_file_close(mpiFile,ierr)
  
  iEntry = 0
  DO jel=1,knel(nlmin)
  
   iel = iaux(jel)
   
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    
    iEntry = iEntry + 1
    Field1(jvt) = daux(iEntry)
    
    if (present(Field2)) then
     iEntry = iEntry + 1
     Field2(jvt) = daux(iEntry)
    end if
    
    if (present(Field3)) then
     iEntry = iEntry + 1
     Field3(jvt) = daux(iEntry)
    end if

    if (present(Field4)) then
     iEntry = iEntry + 1
     Field4(jvt) = daux(iEntry)
    end if

    if (present(Field5)) then
     iEntry = iEntry + 1
     Field5(jvt) = daux(iEntry)
    end if
    
    if (present(Field6)) then
     iEntry = iEntry + 1
     Field6(jvt) = daux(iEntry)
    end if
    
    if (present(Field7)) then
     iEntry = iEntry + 1
     Field7(jvt) = daux(iEntry)
    end if
    
    if (present(Field8)) then
     iEntry = iEntry + 1
     Field8(jvt) = daux(iEntry)
    end if
    
    if (present(Field9)) then
     iEntry = iEntry + 1
     Field9(jvt) = daux(iEntry)
    end if
    
    if (present(Field10)) then
     iEntry = iEntry + 1
     Field10(jvt) = daux(iEntry)
    end if
    
   END DO
  end do
  
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux)
  deallocate (iaux)

 end SUBROUTINE ReLoadMPIFieldQ2_NX
 
 SUBROUTINE ReLoadMPIFieldQ2_X3(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=3
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//TRIM(ADJUSTL(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nLengthV*knel(nlmin)*nF
  allocate(iaux(knel(nlmin)))
  allocate(daux(ndof))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  myFieldOffset = offset 
  CALL MPI_File_read_at(mpiFile,myFieldOffset, iaux, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)

  myFieldOffset = offset + intsize*knel(nlmin) 
  CALL MPI_File_read_at(mpiFile,myFieldOffset, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  CALL mpi_file_close(mpiFile,ierr)
  
  iEntry = 0
  DO jel=1,knel(nlmin)
   iel = iaux(jel)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    iP = nF*(jvt-1)
    Field(iP+1:iP+nF) = daux(iEntry+1:iEntry+nF)
    iEntry = iEntry + nF
   END DO
  end do
  
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux)

 END SUBROUTINE ReLoadMPIFieldQ2_X3

 SUBROUTINE CommCoordinatesWithMaster()
 implicit none
 integer LevDif
 REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
  
  IF (myid.EQ.0) THEN
    CALL CreateDumpStructures(0)
  ELSE
    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
    CALL CreateDumpStructures(LevDif)
  END IF
  
  ILEV = LinSc%prm%MGprmIn%MedLev

  nLengthV = (2**(ILEV-1)+1)**3
  nLengthE = mg_mesh%level(NLMIN)%nel

  ALLOCATE(SendVect(3,nLengthV,nLengthE))

  CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                              mg_mesh%level(ILEV)%kvert,&
                              nLengthV,&
                              nLengthE,&
                              mg_mesh%level(ILEV)%nel,&
                              mg_mesh%level(ILEV)%nvt)
  DEALLOCATE(SendVect)
  
 END SUBROUTINE CommCoordinatesWithMaster
 
 SUBROUTINE LoadMPIFieldQ2_NX(cF,nF,Field1,Field2,Field3,Field4,Field5,Field6,Field7,Field8,Field9,Field10)
 character cF*(*)
 integer :: nF
 REAL*8 Field1(*)
 REAL*8, OPTIONAL :: Field2(*)
 REAL*8, OPTIONAL :: Field3(*)
 REAL*8, OPTIONAL :: Field4(*)
 REAL*8, OPTIONAL :: Field5(*)
 REAL*8, OPTIONAL :: Field6(*)
 REAL*8, OPTIONAL :: Field7(*)
 REAL*8, OPTIONAL :: Field8(*)
 REAL*8, OPTIONAL :: Field9(*)
 REAL*8, OPTIONAL :: Field10(*)
 
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:),jaux(:)
 INTEGER ivt,jvt,jel,kel,lel,NNEL
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))
  end if
  
  IF (bParallel) THEN
   NNEL = INT(ElementOffsets(subnodes+1))
  ELSE
   NNEL = knel(nlmin)
  END IF
  ndof = nF*nLengthV!*NNEL
  allocate(daux(ndof),iaux(NNEL),jaux(NNEL))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)),MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  offset = 0
  myFieldOffset = offset
  
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_read_all(mpiFile, iaux, NNEL, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  
  do iel=1,NNEL
   jaux(iel) = iel
  end do
  
  CALL sort2D(iaux,jaux,NNEL)
  
  DO iel=1,knel(nlmin)
   jel = coarse%myELEMLINK(iel)
   kel = jaux(jel) 
   
   myFieldOffset = offset + intsize*NNEL + dblesize*(kel-1)*nLengthV*nF
   CALL MPI_File_read_at(mpiFile, myFieldOffset, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
   if (ierr.ne.0) iGlobalError = 1
  
   DO jel=1,nLengthV
    jvt = myDump%Vertices(IEL,jel)
    
    iP = nF*(jvt-1)
    iEntry = (jel-1)*nF

    iEntry = iEntry + 1
    Field1(jvt) = daux(iEntry)
    
    if (present(Field2)) then
     iEntry = iEntry + 1
     Field2(jvt) = daux(iEntry)
    end if
    
    if (present(Field3)) then
     iEntry = iEntry + 1
     Field3(jvt) = daux(iEntry)
    end if

    if (present(Field4)) then
     iEntry = iEntry + 1
     Field4(jvt) = daux(iEntry)
    end if

    if (present(Field5)) then
     iEntry = iEntry + 1
     Field5(jvt) = daux(iEntry)
    end if
    
    if (present(Field6)) then
     iEntry = iEntry + 1
     Field6(jvt) = daux(iEntry)
    end if
    
    if (present(Field7)) then
     iEntry = iEntry + 1
     Field7(jvt) = daux(iEntry)
    end if
    
    if (present(Field8)) then
     iEntry = iEntry + 1
     Field8(jvt) = daux(iEntry)
    end if
    
    if (present(Field9)) then
     iEntry = iEntry + 1
     Field9(jvt) = daux(iEntry)
    end if
    
    if (present(Field10)) then
     iEntry = iEntry + 1
     Field10(jvt) = daux(iEntry)
    end if
    
   end do
  end do
  
  CALL mpi_file_close(mpiFile,ierr)
  
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux,iaux,jaux)

 END SUBROUTINE LoadMPIFieldQ2_NX

 SUBROUTINE LoadMPIFieldQ2_X3(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:),jaux(:)
 INTEGER ivt,jvt,jel,kel,lel,NNEL
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=3
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))
  end if
  
  IF (bParallel) THEN
   NNEL = INT(ElementOffsets(subnodes+1))
  ELSE
   NNEL = knel(nlmin)
  END IF
  ndof = nF*nLengthV!*NNEL
  allocate(daux(ndof),iaux(NNEL),jaux(NNEL))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)),MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  offset = 0
  myFieldOffset = offset
  
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_read_all(mpiFile, iaux, NNEL, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
   if (ierr.ne.0) iGlobalError = 1
  
  do iel=1,NNEL
   jaux(iel) = iel
  end do
  
  CALL sort2D(iaux,jaux,NNEL)

  DO iel=1,knel(nlmin)
   jel = coarse%myELEMLINK(iel)
   kel = jaux(jel) 
   
   myFieldOffset = offset + intsize*NNEL + dblesize*(kel-1)*nLengthV*nF
   CALL MPI_File_read_at(mpiFile, myFieldOffset, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
   if (ierr.ne.0) iGlobalError = 1
   
   DO jel=1,nLengthV
    jvt = myDump%Vertices(IEL,jel)
    
    iP = nF*(jvt-1)
    iEntry = (jel-1)*nF
    
    Field(iP+1:iP+nF) = daux(iEntry+1:iEntry+nF) 
    
   end do
  end do
  
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux,iaux,jaux)

 END SUBROUTINE LoadMPIFieldQ2_X3
  
 SUBROUTINE LoadMPIFieldP1(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:),jaux(:)
 INTEGER ivt,jvt,jel,kel,lel,NNEL
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=4
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))
  end if
  
  IF (bParallel) THEN
   NNEL = INT(ElementOffsets(subnodes+1))
  ELSE
   NNEL = knel(nlmin)
  END IF
  ndof = nF*nLengthE!*NNEL
  allocate(daux(ndof),iaux(NNEL),jaux(NNEL))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)),MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  offset = 0
  myFieldOffset = offset
  
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_read_all(mpiFile, iaux, NNEL, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  
  do iel=1,NNEL
   jaux(iel) = iel
  end do
  
  CALL sort2D(iaux,jaux,NNEL)

  DO iel=1,knel(nlmin)
   jel = coarse%myELEMLINK(iel)
   kel = jaux(jel) 
   
   myFieldOffset = offset + intsize*NNEL + dblesize*(kel-1)*nLengthE*nF
   CALL MPI_File_read_at(mpiFile, myFieldOffset, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
   if (ierr.ne.0) iGlobalError = 1
      
   DO jel=1,nLengthE
    kel = myDump%Elements(IEL,jel)
    
    iP = nF*(kel-1)
    iEntry = (jel-1)*nF
    
    Field(iP+1:iP+nF) = daux(iEntry+1:iEntry+nF) 
    
   end do
  end do
  
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux,iaux,jaux)

 END SUBROUTINE LoadMPIFieldP1

 SUBROUTINE LoadMPIFieldP0_INT(cF,Field)
 character cF*(*)
 INTEGER Field(*)
 integer,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:),jaux(:)
 INTEGER ivt,jvt,jel,kel,lel,NNEL
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Loading file:"'//TRIM(ADJUSTL(cPOutFile))
  end if
  
  IF (bParallel) THEN
   NNEL = INT(ElementOffsets(subnodes+1))
  ELSE
   NNEL = knel(nlmin)
  END IF
  ndof = nLengthE!*NNEL
  allocate(daux(ndof),iaux(NNEL),jaux(NNEL))
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)),MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFile,ierr)

  offset = 0
  myFieldOffset = offset
  
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_read_all(mpiFile, iaux, NNEL, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  if (ierr.ne.0) iGlobalError = 1
  
!   if (myid.eq.1) write(*,*) daux
  do iel=1,NNEL
   jaux(iel) = iel
  end do
  
  CALL sort2D(iaux,jaux,NNEL)

  DO iel=1,knel(nlmin)
   jel = coarse%myELEMLINK(iel)
   kel = jaux(jel) 
   
   myFieldOffset = offset + intsize*NNEL + intsize*(kel-1)*nLengthE*1
   CALL MPI_File_read_at(mpiFile, myFieldOffset, daux, ndof, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
   if (ierr.ne.0) iGlobalError = 1
      
   DO jel=1,nLengthE
    kel = myDump%Elements(IEL,jel)
    
    iEntry = jel
    
    Field(kel) = daux(iEntry) 
    
   end do
  end do
  
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux,iaux,jaux)

 END SUBROUTINE LoadMPIFieldP0_INT

 SUBROUTINE SORT2D(LW,KW,N)
   INTEGER LW(N),KW(N),LWA,KWA
   INTEGER I,J,N

   DO I=2,N
   DO J=N,I,-1
   IF (LW(J).LT.LW(J-1)) THEN
     LWA     = LW(J)
     KWA     = KW(J)
     LW(J)   = LW(J-1)
     KW(J)   = KW(J-1)
     LW(J-1) = LWA
     KW(J-1) = KWA
   END IF
   END DO
   END DO

 END SUBROUTINE SORT2D

END SUBROUTINE LoadMPIDumpFiles
!
!
!
SUBROUTINE ReleaseMPIDumpFiles(iOutO,cList)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,MPI_SEEK_SET,MPI_REAL8,MPI_MODE_RDONLY,MPI_MAX
USE PP3D_MPI, ONLY:MPI_COMM_WORLD,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL,mpi_integer,mpi_status_ignore,MPI_COMM_subs,MPI_Offset_kind,mpi_seek_cur,MPI_DOUBLE_PRECISION,subnodes,comm_summn
USE var_QuadScalar,ONLY:QuadSc,LinSc,myDump,mg_mesh,Screw,MaterialDistribution,temperature,mySegmentIndicator
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:GenLinScalar,bParallel

IMPLICIT NONE
character(len=*),intent(in) :: cLIST
INTEGER iOutO
!----------------------------------------------------------------------
INTEGER  :: mpiFile,ierr
CHARACTER cPOutFile*256,cDumpFolder*256
INTEGER iOut,myCoarseElem,nLengthE,nLengthV
INTEGER iP,lCmnd
INTEGER ifilen,i
DATA ifilen/0/
real*8, allocatable :: ElementOffsets(:)

integer iFld,jFld,nFLD
character*(1),allocatable :: cFLD(:)
INTEGER iGlobalError,jGlobalError
INTEGER :: dblesize=8,intsize=4

IF (myid.ne.0) CALL CreateDumpStructures(1)

nFLD = (LEN(cLIST)+1)/2
allocate(cFLD(nFLD))
read(cLIST,*,IOSTAT=iERR) cFLD

allocate(ElementOffsets(subnodes+1))
ElementOffsets = 0d0
if (myid.ne.0) ElementOffsets(myid+1) = DBLE(KNEL(NLMIN))
CALL Comm_SummN(ElementOffsets,subnodes+1)
do i=2,subnodes+1
 ElementOffsets(i) = ElementOffsets(i) + ElementOffsets(i-1)
end do
 
IF (iOutO.LT.0) THEN
 ifilen=ifilen+1
 iOut=MOD(ifilen+insavn-1,insavn)+1
ELSE
 iOut = iOuto
END IF

if (myid.eq.0) then 
 cDumpFolder = '_dump'
 CALL CheckIfFolderIsThereCreateIfNot(cDumpFolder,-1)
 WRITE(cDumpFolder(6:),'(A,I2.2)') '/',iOut
 CALL CheckIfFolderIsThereCreateIfNot(cDumpFolder,-1)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF (myid.NE.0) THEN

 NLMAX = NLMAX+1

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-2)
 nLengthV = (2**(NLMAX-1)+1)**3

 iGlobalError = 0
 DO jFld = 1,nFLD
  if (cFLD(jFld).eq.'p'.or.cFLD(jFld).eq.'P') CALL ReleaseMPIFieldP1('pressure',LinSc%valP(NLMAX-1)%x)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','pressure')
  
  if (cFLD(jFld).eq.'m'.or.cFLD(jFld).eq.'M') CALL ReleaseMPIFieldP0_INT('material',MaterialDistribution(NLMAX-1)%x)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','material')

  if (cFLD(jFld).eq.'x'.or.cFLD(jFld).eq.'X') CALL ReleaseMPIFieldQ2_X3('coordinates',mg_mesh%level(NLMAX)%dcorvg)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','coordinates')
  
  if (cFLD(jFld).eq.'v'.or.cFLD(jFld).eq.'V') CALL ReleaseMPIFieldQ2_NX('velocity',3,QuadSc%valU,QuadSc%valV,QuadSc%valW)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','velocity')
  
  if (cFLD(jFld).eq.'d'.or.cFLD(jFld).eq.'D') CALL ReleaseMPIFieldQ2_NX('distance',1,Screw)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','distance')
  
  if (cFLD(jFld).eq.'s'.or.cFLD(jFld).eq.'S') CALL ReleaseMPIFieldQ2_NX('segment',1,mySegmentIndicator(2,:))
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','segment')
  
  if (cFLD(jFld).eq.'t'.or.cFLD(jFld).eq.'T') CALL ReleaseMPIFieldQ2_NX('temperature',1,temperature)
  CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
  if (jGlobalError.ne.0) CALL ProcessError('R','temperature')
  
  if (cFLD(jFld).eq.'q'.or.cFLD(jFld).eq.'Q') THEN
   if (allocated(GenLinScalar%Fld)) then
    DO iFld=1,GenLinScalar%nOfFields
     CALL ReleaseMPIFieldQ2_NX(adjustl(trim(GenLinScalar%prm%cField(iFld))),1,GenLinScalar%fld(iFld)%Val)
     CALL MPI_ALLREDUCE(iGlobalError,jGlobalError,1,MPI_INTEGER,MPI_MAX,MPI_COMM_SUBS,IERR)
     if (jGlobalError.ne.0) CALL ProcessError('R',adjustl(trim(GenLinScalar%prm%cField(iFld))))
     
    END DO     
   END IF
  END IF
 END DO
 NLMAX = NLMAX-1

END IF

deallocate(cFLD)
deallocate(ElementOffsets)

 CONTAINS
 
 SUBROUTINE ReleaseMPIFieldP0_INT(cF,Field)
 character cF*(*)
 INTEGER Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Writing file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nLengthE*knel(nlmin)
  allocate(iaux(ndof))
  
  iEntry = 0
  DO iel=1,knel(nlmin)
   DO jel=1,nLengthE
    kel = myDump%Elements(IEL,jel)
    iEntry = iEntry + 1
    iaux(iEntry) = Field(kel)
   END DO
  end do
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFile,ierr)

  myFieldOffset = offset + intsize*INT(ElementOffsets(myid))
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, coarse%myELEMLINK, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1
  
  myFieldOffset = offset + intsize*INT(ElementOffsets(subnodes+1)) + intsize*INT(ElementOffsets(myid))*nLengthE
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, iaux, ndof, MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (iaux)
 ! -----------------------------  Write Out the pressure field ----------------------------- ! 
 END SUBROUTINE ReleaseMPIFieldP0_INT

 SUBROUTINE ReleaseMPIFieldP1(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=4
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Writing file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nF*nLengthE*knel(nlmin)
  allocate(daux(ndof))
  
  iEntry = 0
  DO iel=1,knel(nlmin)
   DO jel=1,nLengthE
    kel = myDump%Elements(IEL,jel)
    iP = nF*(kel-1)
    daux(iEntry+1:iEntry+nF) = Field(iP+1:iP+nF)
    iEntry = iEntry + nF
   END DO
  end do
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFile,ierr)

  myFieldOffset = offset + intsize*INT(ElementOffsets(myid))
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, coarse%myELEMLINK, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1
  
  myFieldOffset = offset + intsize*INT(ElementOffsets(subnodes+1)) + dblesize*INT(ElementOffsets(myid))*nLengthE*nF
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux)

  END SUBROUTINE ReleaseMPIFieldP1
 
 SUBROUTINE ReleaseMPIFieldQ2_X3(cF,Field)
 character cF*(*)
 REAL*8 Field(*)
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 integer :: nF=3
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//TRIM(ADJUSTL(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Writing file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nF*nLengthV*knel(nlmin)
  allocate(daux(ndof))
  
  iEntry = 0
  DO iel=1,knel(nlmin)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    iP = nF*(jvt-1)
    daux(iEntry+1:iEntry+nF) = Field(iP+1:iP+nF)
    iEntry = iEntry + nF
   END DO
  end do
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFile,ierr)

  myFieldOffset = offset + intsize*INT(ElementOffsets(myid))
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, coarse%myELEMLINK, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1

  myFieldOffset = offset + intsize*INT(ElementOffsets(subnodes+1)) + dblesize*INT(ElementOffsets(myid))*nLengthV*nF
  call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
  CALL MPI_File_write(mpiFile, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
  IF (IERR.ne.0) iGlobalError = 1
  CALL mpi_file_close(mpiFile,ierr)
  IF (myid.eq.1) then
   WRITE(*,'(A)') ' ==> Done!'
  end if

  deallocate (daux)

 END SUBROUTINE ReleaseMPIFieldQ2_X3
 
 SUBROUTINE ReleaseMPIFieldQ2_NX(cF,nF,Field1,Field2,Field3,Field4,Field5,Field6,Field7,Field8,Field9,Field10)
 character cF*(*)
 integer :: nF
 REAL*8 Field1(*)
 REAL*8, OPTIONAL :: Field2(*)
 REAL*8, OPTIONAL :: Field3(*)
 REAL*8, OPTIONAL :: Field4(*)
 REAL*8, OPTIONAL :: Field5(*)
 REAL*8, OPTIONAL :: Field6(*)
 REAL*8, OPTIONAL :: Field7(*)
 REAL*8, OPTIONAL :: Field8(*)
 REAL*8, OPTIONAL :: Field9(*)
 REAL*8, OPTIONAL :: Field10(*)
 
 reaL*8,  allocatable :: daux(:)
 integer,  allocatable :: iaux(:)
 INTEGER ivt,jvt,jel,kel
 integer iEntry,ndof
 integer(kind=MPI_Offset_kind) :: offset,myFieldOffset
 
  cPOutFile = '_dump/'
  WRITE(cPOutFile(7:),'(I2.2,A,A)') iOut,'/'//ADJUSTL(TRIM(cF)),'.prf'
  offset = 0
  IF (myid.eq.1) then
   WRITE(*,'(A$)') 'Writing file:"'//TRIM(ADJUSTL(cPOutFile))//'"'
  end if
  
  ndof = nF*nLengthV*knel(nlmin)
  allocate(daux(ndof))
  
  iEntry = 0

  DO iel=1,knel(nlmin)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    
    iEntry = iEntry + 1
    daux(iEntry) = Field1(jvt)
    
    if (present(Field2)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field2(jvt)
    end if
    
    if (present(Field3)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field3(jvt)
    end if

    if (present(Field4)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field4(jvt)
    end if

    if (present(Field5)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field5(jvt)
    end if
    
    if (present(Field6)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field6(jvt)
    end if
    
    if (present(Field7)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field7(jvt)
    end if
    
    if (present(Field8)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field8(jvt)
    end if
    
    if (present(Field9)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field9(jvt)
    end if
    
    if (present(Field10)) then
     iEntry = iEntry + 1
     daux(iEntry) = Field10(jvt)
    end if
    
   END DO
  end do
  
  CALL MPI_File_open(MPI_COMM_subs, Adjustl(trim(cPOutFile)), MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFile,ierr)

  if (bParallel) then
   myFieldOffset = offset + intsize*INT(ElementOffsets(myid))
   call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
   CALL MPI_File_write(mpiFile, coarse%myELEMLINK, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
   IF (IERR.ne.0) iGlobalError = 1

   myFieldOffset = offset + intsize*INT(ElementOffsets(subnodes+1)) + dblesize*INT(ElementOffsets(myid))*nLengthV*nF
   call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
   CALL MPI_File_write(mpiFile, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
   IF (IERR.ne.0) iGlobalError = 1
   CALL mpi_file_close(mpiFile,ierr)
   IF (myid.eq.1) then
    WRITE(*,'(A)') ' ==> Done!'
   end if
  else
   if (myid.eq.1) then
    myFieldOffset = offset
    call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
    CALL MPI_File_write(mpiFile, coarse%myELEMLINK, knel(nlmin), MPI_INTEGER, MPI_STATUS_IGNORE,ierr)
    IF (IERR.ne.0) iGlobalError = 1

    myFieldOffset = offset + intsize*knel(nlmin) 
    call MPI_File_seek(mpiFile, myFieldOffset, MPI_SEEK_SET, ierr)
    CALL MPI_File_write(mpiFile, daux, ndof, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
    IF (IERR.ne.0) iGlobalError = 1
    WRITE(*,'(A)') ' ==> Done!'
   end if
  end if

  CALL mpi_file_close(mpiFile,ierr)
  
  deallocate (daux)
  
 end SUBROUTINE ReleaseMPIFieldQ2_NX
 
END SUBROUTINE ReleaseMPIDumpFiles

SUBROUTINE ProcessError(cRW,cF)
USE PP3D_MPI, ONLY:myid,subnodes
USE Sigma_User, ONLY: myErrorCode
CHARACTER*(1) cRW
CHARACTER*(*) cF

IF (cRW.eq.'R') THEN
 if (myid.eq.1) WRITE(*,'(A,I)') ' Error occured by reading the '//TRIM(cF)//' sequence!'
 CALL StopTheProgramFromReader(subnodes,myErrorCode%DUMPFILE_READER)
END IF
IF (cRW.eq.'W') THEN
 if (myid.eq.1) WRITE(*,'(A,I)') ' Error occured by writing the '//TRIM(cF)//' sequence!'
 CALL StopTheProgramFromReader(subnodes,myErrorCode%DUMPFILE_WRITER)
END IF

END SUBROUTINE ProcessError
