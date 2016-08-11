SUBROUTINE ProcessControl(MFILE,MTERM)
 USE PP3D_MPI, ONLY : myid,master,showid,ShareValueK_myMPI,&
                      Barrier_myMPI,ShareValueC_myMPI
 IMPLICIT NONE
 CHARACTER string*200,command*100,cvalue*100
 INTEGER :: myFile=989,iEnd(1),iCmnd(1),iExist(1),iPos,iFile
 INTEGER MFILE,MTERM
 LOGICAL bExist

 iExist=0
 IF (myid.EQ.showid) THEN
  INQUIRE(FILE='_data/ProcCtrl.txt',EXIST=bExist)
  IF (bExist) iExist = 1
 END IF

 CALL ShareValueK_myMPI(iExist,1,showid)

 IF (iExist(1).EQ.0) RETURN

 IF (myid.EQ.showid) OPEN(myFile,FILE='_data/ProcCtrl.txt',STATUS = 'OLD',ACTION='READWRITE')

 DO

   iCmnd = 0

   IF (myid.EQ.showid) THEN
    READ (UNIT=myFile,FMT='(A)',IOSTAT=iEnd(1)) string
    iPos = INDEX(string,"=")
    IF (iPos.NE.0) THEN
     READ(string(1:iPos-1),'(A)') command
     READ(string(iPos+1:),'(A)') cvalue
     cvalue = TRIM(ADJUSTL(cvalue))
    ELSE
     command = "Nothing to be done ... "
    END IF
   END IF

   CALL ShareValueK_myMPI(iEnd,1,showid)
   IF (iEnd(1).EQ.-1) EXIT

   IF (myid.EQ.showid) THEN
     SELECT CASE (TRIM(ADJUSTL(command)))
      CASE ("Dump_Out")
       iCmnd(1) = 1
      CASE ("Dump_In")
       iCmnd(1) = 2
      CASE ("GMV_Out")
       iCmnd(1) = 3

     END SELECT

!      IF (myid.EQ.showid)WRITE(*,*) command,"|",cvalue
   END IF

   CALL ShareValueK_myMPI(iCmnd,1,showid)

   IF (iCmnd(1).EQ.1) THEN
    CALL ShareValueC_myMPI(cvalue,100,showid)
    READ(cvalue,'(I2)') iFile
    IF (iFile.LT.0.OR.iFile.GT.99) THEN
     IF (myid.EQ.showid) WRITE(MTERM,*) "Dump output could not be done ... "
    ELSE
     CALL ReleaseSmartDumpFiles(iFile)
!     CALL Output_DUMPProfiles(iFile)
     CALL FBM_ToFile()
    END IF
   END IF
   IF (iCmnd(1).EQ.2) THEN
    CALL ShareValueC_myMPI(cvalue,100,showid)
    IF (myid.EQ.showid) write(*,*) cvalue
    CALL LoadSmartDumpFiles(cvalue,1)
!    CALL Load_DUMPProfiles(TRIM(ADJUSTL(cvalue)))
   END IF

   IF (iCmnd(1).EQ.3) THEN
    CALL ShareValueC_myMPI(cvalue,100,showid)
    READ(cvalue,'(I4)') iFile
    IF (iFile.LT.0.OR.iFile.GT.9999) THEN
     IF (myid.EQ.showid) WRITE(MTERM,*) "GMV output could not be done ... "
    ELSE
     CALL Output_Profiles(iFile)
    END IF
   END IF

   CALL Barrier_myMPI()

  END DO

  IF (myid.EQ.showid) REWIND(myFile)
  IF (myid.EQ.showid) WRITE(myFile,*) "Nothing to be done ... "
  IF (myid.EQ.showid) CLOSE(myFile)

END SUBROUTINE 


SUBROUTINE Finalize(MFILE,MTERM)
 USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
 IMPLICIT NONE
 INTEGER MFILE,MTERM
 INTEGER iERR

 ! Save the final solution vector in unformatted form
 CALL ReleaseSmartDumpFiles(-1)
! CALL Output_DUMPProfiles()

 IF (myid.eq.showid) THEN
   WRITE(MTERM,*) "PP3D_LES has successfully finished. "
   WRITE(MFILE,*) "PP3D_LES has successfully finished. "
 END IF

 CALL Barrier_myMPI()
 CALL MPI_Finalize(ierr)

END SUBROUTINE Finalize
