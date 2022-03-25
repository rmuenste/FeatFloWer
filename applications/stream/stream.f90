PROGRAM pvdgen
IMPLICIT NONE
CHARACTER(LEN=10) :: cargc1,cargc2,cargc3,cargc4,cargc5,cargc6
INTEGER :: iargc,nargc
INTEGER i,j,jj,k,nL,nR,nF,mF,dAngle,nP,iPeriodicity,iCheckPoints
REAL*8 dt,f,t,dTRot
CHARACTER(LEN=200) cFileName
LOGICAL bExist
INTEGER :: NumOfFiles = 0, NumOfParticles = 0, NumOfFields = 0, ierr = 0 , maxPrtIndex = 0, NumberOfPoints = 0
INTEGER nColumns,nOtherColumns,Columns(4),nOutColumns,OutColumns(4)
INTEGER,allocatable :: OtherColumns(:)
Character*256 cHeader
character*256, allocatable :: cFields(:),cOutFields(:)
!real*8, allocatable :: dValues(:,:,:),dVals(:),dCoor(:,:,:)
real*8, allocatable :: dVals(:)
real*8 daux


type myParticle
 real*4,allocatable :: f(:,:),c(:,:)
 logical,allocatable :: b(:)
 integer :: n = 0
end type myParticle
type(myParticle),allocatable :: prt(:)

i=0
DO 

 cFileName=''
 WRITE(cFileName,'(A,I8.8,A)') '_RTD/Particles_',i,'.csv'
 INQUIRE(FILE=ADJUSTL(TRIM(cFileName)),Exist=bExist)
 if (.not.bExist) exit
 open(FILE=ADJUSTL(TRIM(cFileName)),UNIT=1)
 READ(1,*,iostat=ierr) 
 READ(1,*,iostat=ierr) daux
 close(1)
 if (ierr.ne.0) bExist=.false.
 if (.not.bExist) exit
 i = i + 1
 
END DO

NumOfFiles = i
write(*,*) NumOfFiles,' Files are available'

cFileName=''
WRITE(cFileName,'(A,I8.8,A)') '_RTD/Particles_',0,'.csv'

open(FILE=ADJUSTL(TRIM(cFileName)),UNIT=1)
read(1,iostat=ierr,fmt='(A)') cHeader

!!! Estimate the number of columns in csv
CALL getSubstring(cHeader,',',nColumns,Columns)

if (ierr.ne.0) then
 write(*,*) 'first file is empty ...'
 stop
end if

i=0
maxPrtIndex = 0
allocate(dVals(NumOfFields))

do
 read(1,*,iostat=ierr) dVals

 if (ierr.ne.0) then
  write(*,*) 'reached end of file ...'
  goto 1
 end if
 
 i=i+1
 j=int(dvals(Columns(4)))
 if (maxPrtIndex.lt.j) then
  maxPrtIndex = j
 end if
 
end do

1 close(1)

NumOfParticles = i
write(*,*) NumOfParticles,'  Particles are available', maxPrtIndex
allocate(prt(maxPrtIndex))

allocate(OtherColumns(nColumns-4))
k = 0
DO i=1,nColumns
 bExist = .false.
 do j=1,4
  if (i.eq.Columns(j)) bExist = .true.
 end do
 if (.not.bExist) then
  k = k + 1
  OtherColumns(k) = i
 end if
end do
nOtherColumns = k

write(*,*) 'Columns: ',Columns
write(*,*) 'OtherColumns: ',OtherColumns

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cFileName=''
WRITE(cFileName,'(A,I8.8,A)') '_RTD/Particles_',0,'.csv'

open(FILE=ADJUSTL(TRIM(cFileName)),UNIT=1)
read(1,iostat=ierr,fmt='(A)') cHeader

do
 i=i+1
 read(1,*,iostat=ierr) dVals

 if (ierr.ne.0) then
  write(*,*) 'reached end of file ...'
  goto 2
 end if
 
 j=int(dvals(Columns(4)))
 allocate(prt(j)%c(3,0:NumOfFiles))
 allocate(prt(j)%f(nOtherColumns,0:NumOfFiles))
 allocate(prt(j)%b(0:NumOfFiles))
 prt(j)%b = .false.
 
end do
2 continue

close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NumberOfPoints = 0

DO i=0,NumOfFiles-1

 cFileName=''
 WRITE(cFileName,'(A,I8.8,A)') '_RTD/Particles_',i,'.csv'
 open(FILE=ADJUSTL(TRIM(cFileName)),UNIT=1)
 write(*,*) ADJUSTL(TRIM(cFileName))//' is processed'
 READ(1,*) 
 j=0
 DO 
  read(1,iostat=ierr,fmt='(A)') cHeader
  if (ierr.ne.0) GOTO 3
  READ(cHeader,*) dVals
!   write(*,*) dVals,int(dvals(Columns(4)))
  NumberOfPoints = NumberOfPoints + 1
  j=int(dvals(Columns(4)))
  prt(j)%b(i) = .true.
  prt(j)%n = prt(j)%n + 1
  prt(j)%c(1:3,i) = dvals(Columns(1:3))
  prt(j)%f(1:nOtherColumns,i) = dvals(OtherColumns(1:nOtherColumns))
 END DO
3 close(1)
END DO

NumOfFiles = NumOfFiles + 1
!NumberOfPoints = NumberOfPoints + maxPrtIndex

! Set the particles to default endpoints
! i = NumOfFiles-1
! DO j=1,maxPrtIndex
!  NumberOfPoints = NumberOfPoints + 1
!  prt(j)%b(i)                 = prt(j)%b(i-1)
!  prt(j)%n                    = prt(j)%n + 1
!  prt(j)%c(1:3,i)             = prt(j)%c(1:3,i-1)
!  prt(j)%f(1:nOtherColumns,i) = prt(j)%f(1:nOtherColumns,i-1)
! END DO

OPEN(FILE='_RTD/ParticlesAtOutflow.csv',UNIT=1)
read(1,iostat=ierr,fmt='(A)') cHeader
i = NumOfFiles-1
DO k=1,maxPrtIndex
 read(1,iostat=ierr,fmt='(A)') cHeader
 if (ierr.ne.0) GOTO 4
 READ(cHeader,*) dVals
 NumberOfPoints = NumberOfPoints + 1
 j=int(dvals(Columns(4)))
 prt(j)%b(i) = .true.
 prt(j)%n = prt(j)%n + 1
 prt(j)%c(1:3,i) = dvals(Columns(1:3))
 prt(j)%f(1:nOtherColumns,i) = dvals(OtherColumns(1:nOtherColumns))
END DO
4 close(1)

iCheckPoints = 0
do j=1,maxPrtIndex
  if (allocated(prt(j)%c)) then
   do i=0,NumOfFiles-1
     if (prt(j)%b(i)) iCheckPoints = iCheckPoints + 1
   end do
  end if
end do

write(*,*)  'Check:',iCheckPoints, NumberOfPoints

CALL writeVTP()

! write(*,*) prt(55638)%f(1,0:NumOfFiles-1)
! write(*,*) prt(55638)%f(2,0:NumOfFiles-1)
! write(*,*) prt(55638)%c(3,0:NumOfFiles-1)


!
!-------------------------------------------------------
!
 contains
 
SUBROUTINE writeVTP

REAL*8 dPI,dAngle_Current,dAngle_end,dAngle_start

OPEN (UNIT=1,FILE='_RTD/stream.vtp')

write(1,'(A)') '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
write(1,'(A)') '  <PolyData>'
write(1,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NumberOfPoints,'" NumberOfVerts="0" NumberOfLines="',NumOfParticles,'" NumberOfStrips="0" NumberOfPolys="0">'
write(1,'(A)') '      <PointData Vectors="Normals">'
do k=1,nOtherColumns
 write(1,'(A,A,A)') '        <DataArray type="Float32" Name=',ADJUSTL(TRIM(cFields(OtherColumns(k)))),' format="ascii" RangeMin="-0.00018210224516" RangeMax="0.41778442264">'

 write(*,*) 'writing fields ... '
 do j=1,maxPrtIndex
  if (allocated(prt(j)%c)) then
   do i=0,NumOfFiles-1
     if (prt(j)%b(i)) write(1,*) prt(j)%f(k,i)
   end do
  end if
 end do
 write(1,'(A)') '        </DataArray>'
 write(*,*) 'done!'
end do

write(*,*) 'writing indiceIDs ... '
write(1,'(A,A,A)') '        <DataArray type="Int64" Name=','"indice"',' format="ascii" RangeMin="0" RangeMax="100000000">'
do j=1,maxPrtIndex
 if (allocated(prt(j)%c)) then
  do i=0,NumOfFiles-1
   if (prt(j)%b(i)) then
      write(1,*) j
   end if
  end do
 end if
end do
write(1,'(A)') '        </DataArray>'
write(*,*) 'done!'


write(1,'(A)') '      </PointData>'
write(1,'(A)') '      <CellData>'
write(1,'(A)') '      </CellData>'
write(1,'(A)') '      <Points>'
write(1,'(A)') '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="0" RangeMax="330.00502806">'

write(*,*) 'writing coordinates ... '
do j=1,maxPrtIndex
 if (allocated(prt(j)%c)) then
  do i=0,NumOfFiles-1
    if (prt(j)%b(i)) write(1,*) prt(j)%c(1:3,i)
  end do
 end if
end do
write(*,*) 'done!'

write(1,'(A)') '        </DataArray>'
write(1,'(A)') '      </Points>'
write(1,'(A)') '      <Verts>'
write(1,'(A)') '        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '      </Verts>'
write(1,'(A)') '      <Lines>'
write(1,'(A)') '       <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="117894">'

write(*,*) 'writing connectivity ... '
k=0
do j=1,maxPrtIndex
 if (allocated(prt(j)%c)) then
  do i=0,NumOfFiles-1
    if (prt(j)%b(i)) then
     write(1,*) k
     k=k+1
    end if
  end do
 end if
end do
write(*,*) 'done!'

write(1,'(A)') '        </DataArray>'
write(1,'(A)') '        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="554" RangeMax="117895">'

write(*,*) 'writing offsets ... '
k=0
do j=1,maxPrtIndex
 if (prt(j)%n.gt.0) then
  do i=0,NumOfFiles-1
    if (prt(j)%b(i)) then
     k=k+1
    end if
  end do
  write(1,*) k
 end if
end do
write(*,*) 'done!'

write(1,'(A)') '        </DataArray>'
write(1,'(A)') '      </Lines>'
write(1,'(A)') '      <Strips>'
write(1,'(A)') '        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="117894">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="554" RangeMax="117895">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '      </Strips>'
write(1,'(A)') '      <Polys>'
write(1,'(A)') '        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="117894">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="554" RangeMax="117895">'
write(1,'(A)') '        </DataArray>'
write(1,'(A)') '      </Polys>'
write(1,'(A)') '    </Piece>'
write(1,'(A)') '  </PolyData>'
write(1,'(A)') '</VTKFile>'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
close(1)

END SUBROUTINE writeVTP
 
SUBROUTINE getSubstring(cH,cS,nC,iC)
Character*256 cH
Character*60 cF
Character*1 cS
integer iC(4),nC

integer i,jPos1,jPos2
integer, allocatable :: iPos(:)
Character*1 cA
Integer nS

nC = 0
iC = 0

write(*,*) ",",adjustl(trim(cH)),","

nS = 0
Do i=1,256
 READ(cH(i:i),'(A1)') cA
 if (cA.eq.cS) nS = nS + 1
END DO

allocate(iPos(nS+1))
allocate(cFields(nS+1))

iPos = 0

!write(*,*) 'nS=',nS

nS = 0
Do i=1,256
 READ(cH(i:i),'(A1)') cA
 if (cA.eq.cS) THEN
  nS = nS + 1
  iPos(nS) = i
 end if
END DO

iPos(nS+1) = LEN(adjustl(trim(cH)))+1
nC = nS + 1
NumOfFields = nC

write(*,*) 'nS=',nC, 'iPos= ',iPos

jPos1 = 0
Do i=1,nS+1
 jPos2 = iPos(i)
 READ(cH(jPos1+1:jPos2-1),'(A)') cF
 cFields(i) = ADJUSTL(TRIM(cF))
 
 if (INDEX(ADJUSTL(TRIM(cF)),"coor_X").ne.0) then
  write(*,*) 'x-column is :', i,ADJUSTL(TRIM(cF))
  iC(1) = i
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"coor_Y").ne.0) then
  write(*,*) 'y-column is :', i,ADJUSTL(TRIM(cF))
  iC(2) = i
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"coor_Z").ne.0) then
  iC(3) = i
  write(*,*) 'z-column is :', i,ADJUSTL(TRIM(cF))
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"indice").ne.0) then
  iC(4) = i
  write(*,*) 'indice-column is :', i,ADJUSTL(TRIM(cF))
 end if
 jPos1=jPos2
END DO

write(*,*) 'cFields= '
do i=1,nC
 write(*,*) ADJUSTL(TRIM(cFields(i)))
end do

! pause

END SUBROUTINE getSubstring
!
!-------------------------------------------------------
!

END PROGRAM pvdgen