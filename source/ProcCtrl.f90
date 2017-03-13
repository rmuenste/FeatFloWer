SUBROUTINE ProcessControl(MFILE,MTERM)
 USE PP3D_MPI, ONLY : myid,master,showid,ShareValueK_myMPI,&
                      Barrier_myMPI,ShareValueC_myMPI
 USE Transport_UxyzP_Q2P1, ONLY : myDataFile,QuadSc,LinSc,Properties,&
                        GetVeloParameters,GetPresParameters,GetPhysiclaParameters
 IMPLICIT NONE
 CHARACTER string*200,command*100,cvalue*100,CSimPar*7
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
   END IF

   CALL ShareValueK_myMPI(iEnd,1,showid)
   IF (iEnd(1).EQ.-1) EXIT

   CALL ShareValueC_myMPI(string,200,showid)

   iPos = INDEX(string,"=")
   IF (iPos.NE.0) THEN
    READ(string(1:iPos-1),'(A)') command
    READ(string(iPos+1:),'(A)') cvalue
    cvalue = TRIM(ADJUSTL(cvalue))
   ELSE
    command = "Nothing to be done ... "
   END IF

!================================================================================
   SELECT CASE (TRIM(ADJUSTL(command)))
! -------------------------------------------------------------------------------
   CASE ("Dump_Out")

    READ(cvalue,'(I2)') iFile
    IF (iFile.LT.0.OR.iFile.GT.99) THEN
     IF (myid.EQ.showid) WRITE(MTERM,*) "Dump output could not be done ... "
    ELSE
     CALL SolToFile(iFile)
     CALL FBM_ToFile()
    END IF

! -------------------------------------------------------------------------------
   CASE ("Dump_In")

    IF (myid.EQ.showid) write(*,*) cvalue
    CALL SolFromFile(cvalue,1)
    CALL FBM_FromFile()

! -------------------------------------------------------------------------------
   CASE ("GMV_Out")

    READ(cvalue,'(I4)') iFile
    IF (iFile.LT.0.OR.iFile.GT.9999) THEN
     IF (myid.EQ.showid) WRITE(MTERM,*) "GMV output could not be done ... "
    ELSE
     CALL Output_Profiles(iFile)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Velo")
    IF (myid.EQ.showid) WRITE(MTERM,*) "Reloading Velocity parameters from file :", TRIM(ADJUSTL(cvalue))
    IF (myid.EQ.showid) WRITE(MFILE,*) "Reloading Velocity parameters from file :", TRIM(ADJUSTL(cvalue))
    INQUIRE(FILE=TRIM(ADJUSTL(cvalue)),EXIST=bExist)
    IF (bExist) THEN
     myDataFile = TRIM(ADJUSTL(cvalue))  
     CALL GetVeloParameters(QuadSc%prm,QuadSc%cName,mfile)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Pres")
    IF (myid.EQ.showid) WRITE(MTERM,*) "Reloading Pressure parameters from file :", TRIM(ADJUSTL(cvalue))
    IF (myid.EQ.showid) WRITE(MFILE,*) "Reloading Pressure parameters from file :", TRIM(ADJUSTL(cvalue))
    INQUIRE(FILE=TRIM(ADJUSTL(cvalue)),EXIST=bExist)
    IF (bExist) THEN
     myDataFile = TRIM(ADJUSTL(cvalue))
     CALL GetPresParameters(LinSc%prm,LinSc%cName,mfile)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_SimPar")
    IF (myid.EQ.showid) WRITE(MTERM,*) "Reloading Simulation parameters from file :", TRIM(ADJUSTL(cvalue))
    IF (myid.EQ.showid) WRITE(MFILE,*) "Reloading Simulation parameters from file :", TRIM(ADJUSTL(cvalue))
    INQUIRE(FILE=TRIM(ADJUSTL(cvalue)),EXIST=bExist)
    IF (bExist) THEN
     myDataFile = TRIM(ADJUSTL(cvalue))
     CSimPar = "SimPar"
     CALL  GDATNEW (CSimPar,1)
    END IF

! -------------------------------------------------------------------------------
   CASE ("Reload_Prop")
    IF (myid.EQ.showid) WRITE(MTERM,*) "Reloading Physical parameters from file :", TRIM(ADJUSTL(cvalue))
    IF (myid.EQ.showid) WRITE(MFILE,*) "Reloading Physical parameters from file :", TRIM(ADJUSTL(cvalue))
    INQUIRE(FILE=TRIM(ADJUSTL(cvalue)),EXIST=bExist)
    IF (bExist) THEN
     myDataFile = TRIM(ADJUSTL(cvalue))
     CALL GetPhysiclaParameters(Properties,Properties%cName,mfile)
    END IF

   END SELECT
!================================================================================

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
 CALL SolToFile(-1)

 IF (myid.eq.showid) THEN
   WRITE(MTERM,*) "PP3D_LES has successfully finished. "
   WRITE(MFILE,*) "PP3D_LES has successfully finished. "
 END IF

 CALL Barrier_myMPI()
 CALL MPI_Finalize(ierr)

END SUBROUTINE Finalize
