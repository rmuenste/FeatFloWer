import bpy
import re

def mymatch(name):
  m = re.match(r"Cube",name)
  if m != None:    
      return "Cube"
  
  m = re.match(r"Sphere",name)
  if m != None:    
      return "Sphere"
  
  m = re.match(r"Plane",name)
  if m != None:    
      return "Plane"  
  
  m = re.match(r"Cylinder",name)
  if m != None:    
      return "Cylinder"  
  

bl_info = {"name": "ODE Exporter", "category": "Object"}

class ObjectMoveX(bpy.types.Operator):
    """Write JSON Scene"""      # blender will use this as a tooltip for menu items and buttons.
    bl_idname = "object.move_x"        # unique identifier for buttons and menu items to reference.
    bl_label = "Write JSON Scene"         # display name in the interface.
    bl_options = {'REGISTER', 'UNDO'}  # enable undo for the operator.

    def execute(self, context):        # execute() is called by blender when running the operator.

        # The original script
        scene = context.scene
        first = True
        with open('/home/raphael/myspheres.json','w') as f:           
#        with open('D:/code/GitHub/FC-ode/FullC0ntact/applications/ode-test/cube.json','w') as f:            
            f.write("[\n")                      
            for obj in scene.objects:
                m = mymatch(obj.name)
                # len(obj.modifiers)
                if m != None:
                    if not first:
                        f.write(',\n')
                    else:
                        first = False
                        
                    isDyn = 1
                    if len(obj.modifiers) > 0:
                      isDyn = 0    
                      
                    print("IsDyn = " + str(isDyn))
               
                    loc, rot, scale = obj.matrix_world.decompose() 
                    loc2 = loc
                    loc2.x = 0.0   
                    loc2.y = 0.0   
                    loc2.z = 0.0       
                    
                    norm = obj.data.polygons[0].normal                    
                                                                                           
                    om = obj.matrix_world
                                        
                    wcd = om * loc2                    
#                    print(str(wc))               
                                        
#                    print(str(m))               
                    ang = obj.matrix_world.to_euler()   
                    dim = obj.dimensions                 
                    f.write('{\n')
                    f.write('"Type":'+ '"' + str(m) + '",\n')                    
                    f.write('"IsDynamic":'+ '"' + str(isDyn) + '",\n')                                        
                    f.write('"Pos":'+ '[' + str(wcd.x) + ',' + str(wcd.y) + ',' + str(wcd.z) + '],\n')
                    f.write('"Rot":'+ '[' + str(ang.x) + ',' + str(ang.y) + ',' + str(ang.z) + '],\n')                                        
                    f.write('"Vel":'+ '[0,0,0],\n')                                        
                    f.write('"AngVel":'+ '[0,0,0],\n')                                        
                    f.write('"Norm":'+ '[' + str(norm.x) + ',' + str(norm.y) + ',' + str(norm.z) + '],\n')                                                            
                    f.write('"Dim":'+ '[' + str(dim.x) + ',' + str(dim.y) + ',' + str(dim.z) + ']\n}\n')                                                                                                             
            f.write("]\n")

        return {'FINISHED'}            # this lets blender know the operator finished successfully.

def register():
    bpy.utils.register_class(ObjectMoveX)


def unregister():
    bpy.utils.unregister_class(ObjectMoveX)


# This allows you to run the script directly from blenders text editor
# to test the addon without having to install it.
if __name__ == "__main__":
    register()