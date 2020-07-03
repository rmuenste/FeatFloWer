command = "/home/raphael/bin/blender-2.79b-linux-glibc219-x86_64/blender --background -P extrude.py " \
          "-- %s %s" \
          "%s %f" %(blenderBase, surfaceMeshName, outputLayerName, thickness)