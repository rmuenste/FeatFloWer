module sol_out
!-------------------------------------------------------------------------------------------------
! A module for saving the solution values to 
! a file. This output(dump) is mainly done
! to a binary file.
!-------------------------------------------------------------------------------------------------
use var_QuadScalar, only: knel
contains

!------------------------------------------------------------------------------------------------- 
! Write the pressure solution to file
!-------------------------------------------------------------------------------------------------
! write_pres_sol: The structure of the pressure solution array is:
! pres(1:4*nn)
! the entries pres(4*(iel-1)+1) to pres(4*(iel-1)+4) contain
! the mean value and the dx,dy,dz derivatives in the element iel.
!
! The file format of the dump is:
! Header until the '\n' character
! After the header:
! Global coarse element idx
! #dofsInCoarseElement entries of the mean pressure values 
! #dofsInCoarseElement entries of the d/dx values 
! #dofsInCoarseElement entries of the d/dy values 
! #dofsInCoarseElement entries of the d/dz values 
! ...
!  
! 
! @param iInd 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 
subroutine write_pres_sol(iInd,iiLev,nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

integer, intent(in) :: iInd
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt

elemCoarse = KNEL(nmin)

!IF (myid.ne.0) THEN
! DO i=1,nn
!  Field2(i) = Field1(4*(i-1)+1)
!  Field3(i) = Field1(4*(i-1)+2)
!  Field4(i) = Field1(4*(i-1)+3)
!  Field5(i) = Field1(4*(i-1)+4)
! END DO
!!  WRITE(*,*) Field2(1:nn)
!END IF

! the subdivision level of an element on the 
! output level, NLMAX = 2, iiLev = 0
! 8**(2-1) = 8
! meaning on level 2 a coarse grid
! element is divided into 8 elements
dofsInCoarseElement = 8**((nmax+iiLev)-1)

if(myid.eq.1)then
  open(unit=iunit, file="test.dmp", iostat=istatus, action="write")

  do iel=1,elemCoarse
    ! a scientific format where the number before
    ! the dot is non-zero, the total width of the
    ! number of 14 with 6 digits after the dot and
    ! the ending of the number is EXXX, ie E-10
    ! write(iunit,'(ES14.6)')pres(4*(iel-1)+1)
    write(iunit,'(I0)')elemmap(iel)

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+1) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+2) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+3) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+4) ! mean
      END DO
      write(iunit,*)' '

  end do
  close(iunit)

  write(*,*)'Pres 1: ', pres(1)
  write(*,*)'Elemmap 1: ', elemmap(4)
  write(*,*)'nn 1: ', nn
  write(*,*)'elemCoarse 1: ', elemCoarse
  write(*,*)'dofsInCoarseElement 1: ', dofsInCoarseElement
  DO ivt=1,dofsInCoarseElement
    write(*,*)'edofs(1,:)', edofs(1,ivt) 
  END DO

end if

  call write_sol_pres(iiLev, nn ,elemCoarse, dofsInCoarseElement, elemmap, edofs, pres);

!CALL WriteSol(iInd,iType,iiLev,cFF,Field2,Field3,Field4,Field5)

end subroutine write_pres_sol
!------------------------------------------------------------------------------------------------- 
! Read the pressure solution from a file
!-------------------------------------------------------------------------------------------------
! read_pres_sol: The structure of the pressure solution array is:
! pres(1:4*nn)
! the entries pres(4*(iel-1)+1) to pres(4*(iel-1)+4) contain
! the mean value and the dx,dy,dz derivatives in the element iel.
! @param cInFile 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 
subroutine read_pres_sol(cInFile,iiLev,nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

character*(60) :: cInFile

integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt

elemCoarse = KNEL(nmin)

! the subdivision level of an element on the 
! output level, NLMAX = 2, iiLev = 0
! 8**(2-1) = 8
! meaning on level 2 a coarse grid
! element is divided into 8 elements
dofsInCoarseElement = 8**((nmax+iiLev)-1)

if(myid.ne.0)then

  call read_sol_pres(iiLev, nn ,elemCoarse, dofsInCoarseElement, elemmap, edofs, pres);

end if

end subroutine read_pres_sol
!
!-------------------------------------------------------------------------------------------------
! Write the pressure solution to file
!-------------------------------------------------------------------------------------------------
! write_pres_sol: The structure of the pressure solution array is:
! pres(1:4*nn)
! the entries pres(4*(iel-1)+1) to pres(4*(iel-1)+4) contain
! the mean value and the dx,dy,dz derivatives in the element iel.
! @param iInd 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 
subroutine write_vel_sol(iInd,iiLev,nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

integer, intent(in) :: iInd
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt

elemCoarse = KNEL(nmin)

!IF (myid.ne.0) THEN
! DO i=1,nn
!  Field2(i) = Field1(4*(i-1)+1)
!  Field3(i) = Field1(4*(i-1)+2)
!  Field4(i) = Field1(4*(i-1)+3)
!  Field5(i) = Field1(4*(i-1)+4)
! END DO
!!  WRITE(*,*) Field2(1:nn)
!END IF

! the subdivision level of an element on the 
! output level, NLMAX = 2, iiLev = 0
! 8**(2-1) = 8
! meaning on level 2 a coarse grid
! element is divided into 8 elements
dofsInCoarseElement = 8**((nmax+iiLev)-1)

if(myid.eq.1)then
  open(unit=iunit, file="test.dmp", iostat=istatus, action="write")

  do iel=1,elemCoarse
    ! a scientific format where the number before
    ! the dot is non-zero, the total width of the
    ! number of 14 with 6 digits after the dot and
    ! the ending of the number is EXXX, ie E-10
    ! write(iunit,'(ES14.6)')pres(4*(iel-1)+1)
    write(iunit,'(I0)')elemmap(iel)

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+1) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+2) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+3) ! mean
      END DO
      write(iunit,*)' '

      DO ivt=1,dofsInCoarseElement
        jvt = edofs(iel,ivt)
        !write(*,*)'glob', elemmap(iel) 
        write(iunit,'(A,I0, F10.4)',advance='no')' ', jvt, pres(4*(jvt-1)+4) ! mean
      END DO
      write(iunit,*)' '

  end do
  close(iunit)

  write(*,*)'Pres 1: ', pres(1)
  write(*,*)'Elemmap 1: ', elemmap(4)
  write(*,*)'nn 1: ', nn
  write(*,*)'elemCoarse 1: ', elemCoarse
  write(*,*)'dofsInCoarseElement 1: ', dofsInCoarseElement
  DO ivt=1,dofsInCoarseElement
    write(*,*)'edofs(1,:)', edofs(1,ivt) 
  END DO

end if

  call write_sol_pres(iiLev, nn ,elemCoarse, dofsInCoarseElement, elemmap, edofs, pres);

end subroutine write_vel_sol

!SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
!
!IF (myid.NE.0) THEN
!
! ILEV = NLMIN
!
! nLengthE = 8**((NLMAX+iiLev)-1)
!
! ! (2^(NLMAX+iiLev)+1)^3
! ! dofs on a cube
! nLengthV = (2**((NLMAX+iiLev))+1)**3
!
!!  WRITE(*,*)  'WRITE :::',nLengthE,nLengthV
!END IF
!
! IF (myid.ne.0) dMaxNel = DBLE(KNEL(NLMIN))
! CALL COMM_Maximum(dMaxNel)
!
! IF (myid.eq.showid) THEN
!  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
!  WRITE(MTERM,*) 'Releasing current '//TRIM(ADJUSTL(cField))//' solution into: "'//ADJUSTL(TRIM(cFile)),'"'
! END IF
!
! IF (myid.eq.0) THEN
!  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
!  OPEN(321,FILE=TRIM(ADJUSTL(cFile)))
! END IF
!
! IF (iType.eq.1) THEN
!  IF (Present(Field1)) CALL CollectVertField(Field1)
!  IF (Present(Field2)) CALL CollectVertField(Field2)
!  IF (Present(Field3)) CALL CollectVertField(Field3)
! END IF
!
! IF (iType.eq.2) THEN
!  IF (Present(Field1)) CALL CollectElemField(Field1)
!  IF (Present(Field2)) CALL CollectElemField(Field2)
!  IF (Present(Field3)) CALL CollectElemField(Field3)
!  IF (Present(Field4)) CALL CollectElemField(Field4)
! END IF
!
! IF (iType.eq.3) THEN
!   CALL CollectVertField(ViscoSc%val11)
!   CALL CollectVertField(ViscoSc%val22)
!   CALL CollectVertField(ViscoSc%val33)
!   CALL CollectVertField(ViscoSc%val12)
!   CALL CollectVertField(ViscoSc%val13)
!   CALL CollectVertField(ViscoSc%val23)
! END IF
!
! IF (myid.eq.0) THEN
!  CLOSE(321)
! END IF
!
! CONTAINS
!! -----------------------------------------------------------------
!SUBROUTINE CollectElemField(xField)
!USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
!REAL*8 xField(*)
!
! IF (myid.ne.0) THEN
!
!  CALL SENDI_myMPI(nLengthE,0)
!  CALL SENDI_myMPI(KNEL(NLMIN),0)
!  
!  ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 
!
!  DO iel = 1,KNEL(NLMIN)
!   DO ivt=1,nLengthE
!    jvt = myDump%Elements(IEL,ivt)
!    Field(ivt,iel) = xField(jvt)
!   END DO
!  END DO
! 
!  CALL SENDD_myMPI(Field,nLengthE*KNEL(NLMIN),0)
! 
! ELSE
!
!  DO pID =1,subnodes
!
!   CALL RECVI_myMPI(nLengthE,pID)
!   IF (pID.EQ.1) THEN
!    pnel = INT(dMaxNel)
!    ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 
!    ALLOCATE(auxField(nLengthE,pnel)) 
!   END IF
!   CALL RECVI_myMPI(pnel,pID)
!
!   CALL RECVD_myMPI(auxField,pnel*nLengthE,pID)
!
!   DO I=1,pnel
!   IEL = coarse%pELEMLINK(pID,I)
!    DO ivt=1,nLengthE
!     Field(ivt,iel) = auxField(ivt,I)
!    END DO
!   END DO
!  END DO
! 
!  DEALLOCATE(auxField) 
!!  WRITE(321,'(A)') '-- - ---'
!  DO iel=1,KNEL(NLMIN)
!   WRITE(321,*) Field(1:nLengthE,iel)
!   !WRITE(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
!  END DO
! END IF
!
!DEALLOCATE(Field) 
!
!END SUBROUTINE CollectElemField
!! -----------------------------------------------------------------
!
!END SUBROUTINE WriteSol

subroutine write_pres_test(fileName, nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

CHARACTER*(60) :: fileName
CHARACTER*(20) :: idName
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: i
integer :: jvt
integer :: iiLev = 0

elemCoarse = KNEL(nmin)

!IF (myid.ne.0) THEN
! DO i=1,nn
!  Field2(i) = Field1(4*(i-1)+1)
!  Field3(i) = Field1(4*(i-1)+2)
!  Field4(i) = Field1(4*(i-1)+3)
!  Field5(i) = Field1(4*(i-1)+4)
! END DO
!!  WRITE(*,*) Field2(1:nn)
!END IF

if(myid.ne.0)then
  write(idName,'(I0)')myid
  open(unit=iunit, file=trim(adjustl(fileName))//trim(adjustl(idName)), iostat=istatus, action="write")

 do i=1,elemCoarse
  write(iunit,*)pres(4*(i-1)+1)
  write(iunit,*)pres(4*(i-1)+2)
  write(iunit,*)pres(4*(i-1)+3)
  write(iunit,*)pres(4*(i-1)+4)
 end do

 close(iunit)

end if

end subroutine write_pres_test

end module sol_out

