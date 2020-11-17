import os
import shutil
import sys
import getopt
import GeometryManipulator as gm

# python3 GeometryManipulator_Driver.py -f Endstueck_6_5.off -t 12.5 

def usage():
  print("Usage: heat_start.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-n, --num-processors]: defines the number of parallel jobs to be used")
  print("[-f, --in-folder]: defines the input project file (will be replaced in the q2p1-data-file)")

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:s:t:', ['help','file=', 'scale=','translation='])
    except getopt.GetoptError:
        sys.exit(2)

    dZ=0.0
    dS=1.0
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-t", "--translation"):
            dZ = float(arg)
#            print(dZ)
        elif opt in ("-s", "--scale"):
            dS = float(arg)
#            print(dZ)
        elif opt in ("-f", "--file"):
            myfile = arg
#            print(myfile)


    shortfilename = os.path.basename(myfile)
    #print(shortfilename)
    filename, file_extension = os.path.splitext(shortfilename)
    path = myfile.replace(shortfilename,"")
 
    print(myfile + " is translated by "  +str(dZ))
    myGeo = gm.GeometryObject()
    if file_extension == '.off':
      myGeo.ReadOff(myfile)
    if file_extension == '.obj':
      myGeo.ReadObj(myfile)
    if file_extension == '.ply':
      myGeo.ReadPly(myfile)
      
    myGeo.TranslateObject(0.0,0.0,dZ)
    myGeo.ScaleObject(1.0,1.0,dS)
    
    outfile = path + filename + "_Modified.off"
    print(outfile)
    myGeo.ExportOff(outfile)


if __name__ == "__main__":
    main()
