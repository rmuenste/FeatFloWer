      SUBROUTINE StatOut(time_passed,myOutFile)
      USE def_FEAT
      USE PP3D_MPI, ONLY : myid,master,showid,subnodes
      USE var_QuadScalar, ONLY : myStat,bNonNewtonian,myMatrixRenewal
      IMPLICIT NONE

      Real, intent(in) :: time_passed

      REAL*8 daux,daux1,ds
      INTEGER myFile,myOutFile,itms,istat
      LOGICAL bExist

      itms = min(itns-1,nitns-1)
      ds = DBLE(subnodes)

      IF (myid.eq.showid) THEN

      IF (myOutFile.eq.0) THEN
       myFile = 669
       OPEN (UNIT=myFile, FILE='_data/Statistics.txt',action='write',iostat=istat)
       if(istat .ne. 0)then
         write(*,*)"Could not open file for writing in StatOut(). "
       stop          
       end if
      ELSE
       myFile = myOutFile
      END IF

      WRITE(myFile,6) "Number of timesteps                ",itms
      WRITE(myFile,6) "Linear MG iterations UVW           ",myStat%iLinUVW
      WRITE(myFile,6) "Linear MG iterations   P           ",myStat%iLinP
      WRITE(myFile,6) "NonLinear iterations UVW           ",myStat%iNonLin
      WRITE(myFile,7) "Average MG-UVW iters per NonLin    ",myStat%iLinUVW/DBLE(myStat%iNonLin)
      WRITE(myFile,7) "Average MG-P   iters per timesteps ",myStat%iLinP/DBLE(itms)
      WRITE(myFile,7) "Average NonLin iters per timesteps ",myStat%iNonLin/DBLE(itms)
      WRITE(myFile,*) 
      WRITE(myFile,8) "Overall time           ",time_passed,time_passed*ds
      WRITE(myFile,*) 

      daux1=myStat%tDefUVW
      IF(bNonNewtonian.AND.myMatrixRenewal%S.EQ.0) daux1=myStat%tDefUVW-myStat%tSMat

      daux = myStat%tMGUVW+myStat%tMGP+daux1+myStat%tDefP + myStat%tCorrUVWP
      WRITE(myFile,9) "Computation time        ",daux,daux*ds,1d2*daux/time_passed
      WRITE(myFile,9) "Linear MG UVW           ",myStat%tMGUVW,myStat%tMGUVW*ds,1d2*myStat%tMGUVW/time_passed
      WRITE(myFile,9) "Linear MG   P           ",myStat%tMGP,myStat%tMGP*ds,1d2*myStat%tMGP/time_passed
      WRITE(myFile,9) "Defect computation UVW  ",daux1,daux1*ds,1d2*daux1/time_passed
      WRITE(myFile,9) "Defect computation   P  ",myStat%tDefP,myStat%tDefP*ds,1d2*myStat%tDefP/time_passed
      WRITE(myFile,9) "Correction UVWP         ",myStat%tCorrUVWP,myStat%tCorrUVWP*ds,1d2*myStat%tCorrUVWP/time_passed
      WRITE(myFile,*) 
      daux = myStat%tProlUVW+myStat%tRestUVW+myStat%tSmthUVW+myStat%tSolvUVW
      WRITE(myFile,9) "UVW Multigrid           ",daux,daux*ds,1d2*daux/myStat%tMGUVW
      WRITE(myFile,9) "Prolongation            ",myStat%tProlUVW,myStat%tProlUVW*ds,1d2*myStat%tProlUVW/myStat%tMGUVW
      WRITE(myFile,9) "Restriction             ",myStat%tRestUVW,myStat%tRestUVW*ds,1d2*myStat%tRestUVW/myStat%tMGUVW
      WRITE(myFile,9) "Smoothening             ",myStat%tSmthUVW,myStat%tSmthUVW*ds,1d2*myStat%tSmthUVW/myStat%tMGUVW
      WRITE(myFile,9) "Coarsegrid solver       ",myStat%tSolvUVW,myStat%tSolvUVW*ds,1d2*myStat%tSolvUVW/myStat%tMGUVW
      WRITE(myFile,*) 
      daux = myStat%tProlP+myStat%tRestP+myStat%tSmthP+myStat%tSolvP
      WRITE(myFile,9) "P    Multigrid          ",daux,daux*ds,1d2*daux/myStat%tMGP
      WRITE(myFile,9) "Prolongation            ",myStat%tProlP,myStat%tProlP*ds,1d2*myStat%tProlP/myStat%tMGP
      WRITE(myFile,9) "Restriction             ",myStat%tRestP,myStat%tRestP*ds,1d2*myStat%tRestP/myStat%tMGP
      WRITE(myFile,9) "Smoothening             ",myStat%tSmthP,myStat%tSmthP*ds,1d2*myStat%tSmthP/myStat%tMGP
      WRITE(myFile,9) "Coarsegrid solver       ",myStat%tSolvP,myStat%tSolvP*ds,1d2*myStat%tSolvP/myStat%tMGP
      WRITE(myFile,*) 
      daux = myStat%tKMat+myStat%tDMat+myStat%tMMat+myStat%tCMat+myStat%tSMat
      WRITE(myFile,9) "Operator assembly time  ",daux,daux*ds,1d2*daux/time_passed
      WRITE(myFile,9) "Convection matrix       ",myStat%tKMat,myStat%tKMat*ds,1d2*myStat%tKMat/time_passed
      WRITE(myFile,9) "Deformation matrix      ",myStat%tSMat,myStat%tSMat*ds,1d2*myStat%tSMat/time_passed
      WRITE(myFile,9) "Diffusion matrix        ",myStat%tDMat,myStat%tDMat*ds,1d2*myStat%tDMat/time_passed
      WRITE(myFile,9) "Mass matrix             ",myStat%tMMat,myStat%tMMat*ds,1d2*myStat%tMMat/time_passed
      WRITE(myFile,9) "Pressure-Schur matrix   ",myStat%tCMat,myStat%tCMat*ds,1d2*myStat%tCMat/time_passed
      WRITE(myFile,*) 
      daux = myStat%tGMVOut +myStat%tDumpOut
      WRITE(myFile,9) "Output time             ",daux,daux*ds,1d2*daux/time_passed
      WRITE(myFile,9) "GMV output              ",myStat%tGMVOut,myStat%tGMVOut*ds,1d2*myStat%tGMVOut/time_passed
      WRITE(myFile,9) "Dump file output        ",myStat%tDumpOut,myStat%tDumpOut*ds,1d2*myStat%tDumpOut/time_passed

      daux = myStat%tCommS + myStat%tCommV + myStat%tCommP
      WRITE(myFile,*) 
      WRITE(myFile,9) "Communication time      ",daux,daux*ds,1d2*daux/time_passed
      WRITE(myFile,9) "Maximum                 ",myStat%tCommS,myStat%tCommS*ds,1d2*myStat%tCommS/time_passed
      WRITE(myFile,9) "Velocity                ",myStat%tCommV,myStat%tCommV*ds,1d2*myStat%tCommV/time_passed
      WRITE(myFile,9) "Pressure                ",myStat%tCommP,myStat%tCommP*ds,1d2*myStat%tCommP/time_passed
      
      IF (myOutFile.eq.0) THEN
       CLOSE (myFile)
      END IF

      END IF

6     FORMAT(A35,' : ',I12  )
7     FORMAT(A35,' : ',F12.4)
8     FORMAT(A24,' : ',F12.4,'s |',F12.4,'s')
9     FORMAT(A24,' : ',F12.4,'s | ',F12.4,'s | ',F12.4,'%')

      END
