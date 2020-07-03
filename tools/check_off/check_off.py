import bpy,os,re
import sys
import mathutils
import bmesh


def checkManifold(fileName):
  objectName = os.path.basename(fileName)
  objectName = os.path.splitext(objectName)[0]

  count = 0    
  myObjects = bpy.data.objects
  myScene = bpy.data.scenes['Scene']
  for o in myObjects:
#    if o.type == 'MESH':
    if re.search(objectName, o.name):
        print(o.name)                
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

def main():
  fileName = sys.argv[5]
  print("Name : %s" %fileName)
  bpy.ops.import_mesh.off(filepath=fileName)
  checkManifold(fileName)            
      
if __name__ == "__main__":
  main()  