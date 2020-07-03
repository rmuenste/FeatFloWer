import os
import subprocess
import getopt
import sys

def usage():
  print("Usage: check_off.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-p, --blender-path]: path to blender executable")
  print("[-f, --file]: the off file to check")

def main():
    blenderPath = "/home/raphael/bin/blender-2.79b-linux-glibc219-x86_64/blender" 
    fileName = "open_cylinder.off"
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:p:', ['help', 'file=', 'blender-path='])
    except getopt.GetoptError:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-p", "--blender-path"):
            blenderPath = arg
        elif opt in ("-f", "--file"):
            fileName = arg
        else:
            usage()
            sys.exit(2)

    command = "%s --background -P check_off.py " \
            "-- %s" %(blenderPath, fileName)

    subprocess.call([command], shell=True)

if __name__ == "__main__":
    main()