import os
import shutil
import sys
import getopt
import GeometryManipulator as gm

def usage():
  print("Usage: heat_start.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-n, --num-processors]: defines the number of parallel jobs to be used")
  print("[-f, --in-folder]: defines the input project file (will be replaced in the q2p1-data-file)")

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:t:', ['help','file=', 'translation='])
    except getopt.GetoptError:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-t", "--translation"):
            dZ = float(arg)
            print(dZ)
        elif opt in ("-f", "--file"):
            myfile = arg
            print(myfile)


    myGeo = gm.GeometryObject()
    myGeo.ReadOff('Gitter/' + myfile)
    myGeo.TranslateObject(0.0,0.0,12.5)
    myGeo.ExportOff('GitterTranslated/' + myfile)


if __name__ == "__main__":
    main()
