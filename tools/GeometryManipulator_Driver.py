import os
import shutil
import sys
import getopt
import GeometryManipulator as gm
import array as arr

# python3 GeometryManipulator_Driver.py -f Endstueck_6_5.off -t 12.5 

def usage():
  print("Usage: heat_start.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-n, --num-processors]: defines the number of parallel jobs to be used")
  print("[-f, --in-folder]: defines the input project file (will be replaced in the q2p1-data-file)")

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:s:t:r:', ['help','file=', 'scale=','translation=','rotation='])
    except getopt.GetoptError:
        sys.exit(2)

#    dZ=0.0
    translation=arr.array('d', [0.0, 0.0, 0.0])
    scale=arr.array('d', [1.0, 1.0, 1.0])
#    dS=1.0
    rotation=" "
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-t", "--translation"):
            translation = arg.split(',')
#            dZ = float(arg)
#            print(dZ)
        elif opt in ("-s", "--scale"):
            scale = arg.split(',')
#            dS = float(arg)
#            print(dZ)
        elif opt in ("-r", "--rotation"):
            rotation = arg
        elif opt in ("-f", "--file"):
            myfile = arg
#            print(myfile)


    shortfilename = os.path.basename(myfile)
    #print(shortfilename)
    filename, file_extension = os.path.splitext(shortfilename)
    path = myfile.replace(shortfilename,"")
 
    print(myfile + " is translated by "  + str(translation[0]) + "," + str(translation[1]) + "," + str(translation[2])) 
    print(myfile + " is scaled by "  + str(scale[0]) + "," + str(scale[1]) + "," + str(scale[2])) 

    myGeo = gm.GeometryObject()
    if file_extension == '.off':
      myGeo.ReadOff(myfile)
    if file_extension == '.obj':
      myGeo.ReadObj(myfile)
    if file_extension == '.ply':
      myGeo.ReadPly(myfile)
      
    myGeo.TranslateObject(translation[0],translation[1],translation[2])
#    myGeo.TranslateObject(0.0,0.0,dZ)
    myGeo.ScaleObject(scale[0],scale[1],scale[2])

    if rotation == 'X':
      myGeo.RotateX(90)
    if rotation == 'Y':
      myGeo.RotateY(90)
    if rotation == 'Z':
      myGeo.RotateZ(90)
    if rotation == '-X':
      myGeo.RotateX(-90)
    if rotation == '-Y':
      myGeo.RotateY(-90)
    if rotation == '-Z':
      myGeo.RotateZ(-90)
   
    outfile = path + filename + "_Modified.off"
    print(outfile)
    myGeo.ExportOff(outfile)


if __name__ == "__main__":
    main()
