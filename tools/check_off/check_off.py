import bpy,os,re
import mathutils
import bmesh

def checkManifold():
  count = 0    
  myObjects = bpy.data.objects
  myScene = bpy.data.scenes['Scene']
  for o in myObjects:
    if o.type == 'MESH':
        print(o.name)                
        myScene = bpy.data.scenes['Scene']
        bpy.context.scene.objects.active = o
        o.select = True
        bpy.ops.object.mode_set(mode='EDIT')
        #bpy.ops.mesh.select_all(action='SELECT')        
        bpy.ops.mesh.select_non_manifold() 
        bm = bmesh.from_edit_mesh(o.data)
        verts = [v.index for v in bm.verts if v.select]        
        edges = [e.index for e in bm.edges if e.select]                
        
        
        if len(verts) > 0 or len(edges) > 0:
            print("There are %s non-manifold vertices." % len(verts))
            print("There are %s non-manifold edges." % len(edges))                                

def main():
  checkManifold()            
      
if __name__ == "__main__":
  main()  