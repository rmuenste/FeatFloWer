import bpy,os,re
import sys
import mathutils
import bmesh

def loadFile(fileName):
  splitName = os.path.splitext(fileName)
  ext = splitName[1]

  if ext == ".stl":
#     print("stl")
    if os.path.exists(fileName): 
      bpy.ops.import_mesh.stl(filepath=fileName)
    else:
      print("File %s was not found on the system." %fileName)
      sys.exit(2)
  elif ext == ".off":
#    print("off")
    if os.path.exists(fileName): 
      bpy.ops.import_mesh.off(filepath=fileName)
    else:
      print("File %s was not found on the system." %fileName)
      sys.exit(2)
  else:
    print("File %s has extension %s which cannot be handled." %(fileName, ext))
    sys.exit(2)

def checkManifold(fileName):
  objectName = os.path.basename(fileName)
  objectName = os.path.splitext(objectName)[0]

  count = 0
  myObjects = bpy.data.objects
  myScene = bpy.data.scenes['Scene']
  found = False
  for o in myObjects:
#    if o.type != 'MESH' and o.name != 'Cube':
#        print(o.name, o.type)
#        continue
#    if re.search(objectName, o.name, re.IGNORECASE):
    if o.type == 'MESH' and o.name != 'Cube':
        print(o.name, o.type)
        #print("Mesh name (should be the same as the file name without extension): %s" %o.name)
        found = True
        myScene = bpy.data.scenes['Scene']
        bpy.context.scene.objects.active = o
        o.select = True
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='TOGGLE')
        bpy.ops.mesh.select_non_manifold()
        bm = bmesh.from_edit_mesh(o.data)
        verts = [v.index for v in bm.verts if v.select]
        edges = [e.index for e in bm.edges if e.select]

        print("There are %s non-manifold vertices." % len(verts))
        print("There are %s non-manifold edges." % len(edges))
        if len(verts) > 0 or len(edges) > 0:
            print("This is a non-manifold mesh.")
        else:
            print("This is a manifold mesh.")

  if not found:
      print("A mesh that matches the file name could not be found in %s" %fileName)
def main():

  argv = sys.argv

  try:
      index = argv.index("--") + 1
  except ValueError:
      index = len(argv)
  argv = argv[index:]

  fileName = argv[0]
  print("Name : %s" %fileName)
  loadFile(fileName)
  checkManifold(fileName)

if __name__ == "__main__":
  main()
