
import os

import bpy
import bmesh
import mathutils as M
from mathutils.bvhtree import BVHTree
from bpy.props import (StringProperty, BoolProperty, IntProperty, FloatProperty, FloatVectorProperty, EnumProperty, PointerProperty)
from bpy.types import (Panel, Menu, Operator, PropertyGroup)

#==============================Properties ================================================================
class GeneratorProperties(PropertyGroup):
    radius: FloatProperty(name="Sphere Radius", default=6, min=0.1, max=10000)
    subDivs: IntProperty(name="Mesh Subdivisions", default=256, min=1, max=100000)
    initHeight: FloatProperty(name="Initial Noise Height", default=0.1, min=0.01, max=10000)
    noiseFreq: FloatProperty(name="Initial Noise Frequency", default=1, min=0.01, max=1000)
    
#==============================Operators ==================================================================
#Operator for initializing the terrain
class WM_OT_Initiate(Operator):
    bl_label = "Initiate Terrain"
    bl_idname = "wm.initiate_terrain"
    
    #Main function that gets called by operators
    def execute(self, context):
        scene = context.scene
        props = scene.properties
        
        #We use the ANT Landscape add on to initiate our landscape and set properties as required
        bpy.ops.mesh.landscape_add(
            refresh = True, 
            sphere_mesh = True,
            subdivision_x = props.subDivs,
            subdivision_y = props.subDivs,
            mesh_size = 2 * props.radius,
            height = props.initHeight,
            maximum = props.initHeight,
            minimum = - props.initHeight,
            noise_size = props.noiseFreq
            )
        
        return{"FINISHED"}

#===============================Panels ====================================================================
class TectonicToolsMainPanel:
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Tools"
    #bl_options = {"DEFAULT_CLOSED"}


#Main panel for tectonic tools
class OBJECT_PT_GeneratorPanel(TectonicToolsMainPanel, Panel):
    bl_idname = "OBJECT_PT_generatorPanel"
    bl_label = "Terrain Generator"

    #This is where we add components to the user interface
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

class OBJECT_PT_InitialProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Global Properties"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        layout.prop(props, "radius")
        layout.prop(props, "subDivs")
        layout.prop(props, "initHeight")
        layout.prop(props, "noiseFreq")

class OBJECT_PT_InitializationTools(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Buttons"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        #For initializing our terrain
        layout.operator("wm.initiate_terrain")



#=============================Register Classes to Blender ================================================
classes = (
    GeneratorProperties,
    WM_OT_Initiate,
    #TectonicToolsMainPanel,
    OBJECT_PT_GeneratorPanel,
    OBJECT_PT_InitialProperties,
    OBJECT_PT_InitializationTools
)

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.Scene.properties = PointerProperty(type=GeneratorProperties)
    #bpy.types.Mesh.WMTPlatePointers = PointerProperty(type=WMTPlatePointers)

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    del bpy.types.Scene.properties
    del bpy.types.Mesh.WMTRotationPointers


if __name__ == "__main__":
    register()