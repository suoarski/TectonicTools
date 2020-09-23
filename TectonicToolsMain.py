import os

import bpy
import bmesh
import mathutils as M

from bpy.props import (
    StringProperty, 
    BoolProperty, 
    IntProperty, 
    FloatProperty, 
    FloatVectorProperty, 
    EnumProperty, 
    PointerProperty
    )
from bpy.types import (
    Panel, 
    Menu, 
    Operator, 
    PropertyGroup
    )

#==============================Properties ================================================================
class GeneratorProperties(PropertyGroup):

	#Properties that are used during the planet initialization
    radius: FloatProperty(name="Sphere Radius", default=6, min=0.1, max=10000)
    subDivs: IntProperty(name="Mesh Subdivisions", default=256, min=1, max=100000)
    initHeight: FloatProperty(name="Initial Noise Height", default=0.1, min=0.01, max=10000)
    noiseFreq: FloatProperty(name="Initial Noise Frequency", default=1, min=0.01, max=1000)
    
    #Identifiers so that we can retrieve relevant data throughout code
    planetID: StringProperty(name="Planet Identifier")
    subFrontsID: StringProperty(name="Subduction Front Identifier")

    deltaTime: FloatProperty(name="Time Step Size (dt)", default=0.1, min=0.001, max=10000)
    plateSpeed: FloatProperty(name="Plate Speed", default=1, min=0.01, max=1000)

    
    #This enum property needs to be identical to the one in the ANT Landscape code, so I copy and pasted it
    noiseType: EnumProperty(
        name="Noise Type",
        default='hetero_terrain',
        description="Noise type",
        items = [
            ('multi_fractal', "Multi Fractal", "Blender: Multi Fractal algorithm", 0),
            ('ridged_multi_fractal', "Ridged MFractal", "Blender: Ridged Multi Fractal", 1),
            ('hybrid_multi_fractal', "Hybrid MFractal", "Blender: Hybrid Multi Fractal", 2),
            ('hetero_terrain', "Hetero Terrain", "Blender: Hetero Terrain", 3),
            ('fractal', "fBm Fractal", "Blender: fBm - Fractional Browninian motion", 4),
            ('turbulence_vector', "Turbulence", "Blender: Turbulence Vector", 5),
            ('variable_lacunarity', "Distorted Noise", "Blender: Distorted Noise", 6),
            ('marble_noise', "Marble", "A.N.T.: Marble Noise", 7),
            ('shattered_hterrain', "Shattered hTerrain", "A.N.T.: Shattered hTerrain", 8),
            ('strata_hterrain', "Strata hTerrain", "A.N.T: Strata hTerrain", 9),
            ('ant_turbulence', "Another Noise", "A.N.T: Turbulence variation", 10),
            ('vl_noise_turbulence', "vlNoise turbulence", "A.N.T: Real vlNoise turbulence", 11),
            ('vl_hTerrain', "vlNoise hTerrain", "A.N.T: vlNoise hTerrain", 12),
            ('distorted_heteroTerrain', "Distorted hTerrain", "A.N.T distorted hTerrain", 13),
            ('double_multiFractal', "Double MultiFractal", "A.N.T: double multiFractal", 14),
            ('rocks_noise', "Noise Rocks", "A.N.T: turbulence variation", 15),
            ('slick_rock', "Slick Rock", "A.N.T: slick rock", 16),
            ('planet_noise', "Planet Noise", "Planet Noise by: Farsthary", 17),
            ('blender_texture', "Blender Texture - Texture Nodes", "Blender texture data block", 18)])
    
#==============================Operators ==================================================================
#Operator for initializing the terrain
class WM_OT_Initiate(Operator):
    bl_label = "Initiate Terrain"
    bl_idname = "wm.initiate_terrain"
    
    #Main function that gets called by operators
    def execute(self, context):
        scene = context.scene
        props = scene.properties

        initializeProfileCurves()
        
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
            noise_size = props.noiseFreq,
            noise_type = props.noiseType
            )
        
        #Store the name of the newly created planet so that we can reference it later
        planet = bpy.context.active_object
        props.planetID = planet.name
        
        
        return{"FINISHED"}


class WM_OT_AddSubFront(Operator):
    bl_label = "Draw Subduction Front"
    bl_idname = "wm.add_sub_front"
    
    def execute(self, context):
        scene = context.scene
        props = scene.properties
        
        #We create a new Bezier Curve to represent the subduction front and prepare the curve
        #such that the user may begin to draw onto our sphere
        bpy.ops.curve.primitive_bezier_curve_add(enter_editmode=True)
        bpy.ops.curve.delete(type='VERT')
        bpy.ops.wm.tool_set_by_id(name="builtin.draw")
        bpy.data.scenes['Scene'].tool_settings.curve_paint_settings.depth_mode = "SURFACE"
        return{"FINISHED"}

class WM_OT_RunSubduction(Operator):
    bl_label = "Run Subduction"
    bl_idname = "wm.run_subduction"
    
    def execute(self, context):
        scene = context.scene
        props = scene.properties
        
        
        
        bpy.ops.object.mode_set(mode='OBJECT')
        
        #We check if the SubFronts collection already exists, and create one otherwise
        #Collections are a directory like structure for objects within the scene
        defaultCollection = bpy.data.collections["Collection"]
        if "SubFronts" not in bpy.data.collections:
            bpy.ops.collection.create(name="SubFronts")
            subFrontCollection = bpy.data.collections["SubFronts"]
            bpy.context.scene.collection.children.link(subFrontCollection)
        else:
            subFrontCollection = bpy.data.collections["SubFronts"]
        


        fronts = bpy.context.selected_objects
        for front in fronts:
            if front.type == 'CURVE':
                defaultCollection.objects.unlink(front)
                bpy.ops.object.mode_set(mode='OBJECT')
                bpy.ops.object.convert(target='MESH')
                bpy.ops.object.mode_set(mode='EDIT')
                bpy.ops.mesh.separate(type='LOOSE')
            else:
                front.select_set(False)
                print("Wrong Object Selected")

        '''
        subFronts = bpy.context.active_object
        props.subFrontsID = subFronts.name
        
        #Create Bmesh objects to allow for mesh manipulations within python code
        planet = bpy.data.objects.get(props.planetID)
        subFrontsBM = bmesh.from_edit_mesh(subFronts.data)
        planetBM = bmesh.new()
        planetBM.from_mesh(planet.data)
        '''




        '''
        #bpy.ops.outliner.id_operation(type='UNLINK')
        bpy.ops.object.mode_set(mode='EDIT')

        #Convert subfronts from type bezier curve to a mesh type object
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.convert(target='MESH')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.separate(type='LOOSE')
        


        #bpy.ops.outliner.collection_new(nested=False)
        print(fronts)
        for ob in subFrontCollection.objects:
            print(ob)

        #Store the name of the newly created mesh object representing subduction fronts
        
        

        
        print(subFrontsBM.loops)

        loops = (l for f in subFrontsBM.faces for l in f.loops)

        for l in loops:
            print(l)
        '''
        #bpy.ops.mesh.edge_face_add()
        #bpy.ops.mesh.separate(type='LOOSE')



        '''
        #Create a K-Dimensional Tree object for faster spacial searches
        subFrontsKD = M.kdtree.KDTree(len(subFrontsBM.verts))
        for i, v in enumerate(subFrontsBM.verts):
            subFrontsKD.insert(v.co, i)
        subFrontsKD.balance()
        
        

        
        for v in planetBM.verts:
            co, index, dist = subFrontsKD.find(v.co)
            
            print("VertCo = {}\nClosestSubFront = {}\nDistance = {}".format(v.co, co, dist))
        '''
        #print("Planet Name = {}, Sub Fronts Name = {}".format(props.planetID, props.subFrontsID))
        
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

class OBJECT_PT_Buttons(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Buttons"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        #For initializing our terrain
        layout.operator("wm.initiate_terrain")
        layout.operator("wm.add_sub_front")
        layout.operator("wm.run_subduction")

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
        layout.prop(props, "noiseType")

class OBJECT_PT_SubductionProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Subduction Properties"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        layout.prop(props, "deltaTime")
        layout.prop(props, "plateSpeed")

        #Node groups are used to store profile curves, 
        #and profile curves are used to allow the user to custimize functions within the code (Eg. Distance Transfer) 
        if 'ProfileCurves' not in bpy.data.node_groups:
        	bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')

        curveTree = bpy.data.node_groups['ProfileCurves'].nodes
        layout.label(text="Subduction Uplift Distance Transfer")
        if "Distance Transfer" in curveTree.keys():
        	layout.template_curve_mapping(curveTree["Distance Transfer"], "mapping")


def initializeProfileCurves():
    if 'ProfileCurves' not in bpy.data.node_groups:
        bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
    
    nods = bpy.data.node_groups['ProfileCurves'].nodes
    if "Distance Transfer" not in nods.keys():
        distCurves = nods.new('ShaderNodeRGBCurve')
        distCurves.name = "Distance Transfer"
        pnts = distCurves.mapping.curves[3].points
        pnts[0].location = [0.0, 0.8]
        pnts[1].location = [1.0, 0.0]
        pnts.new(0.2, 1.0)
        pnts.new(0.75, 0.15)


#=============================Register Classes to Blender ================================================
classes = (
    GeneratorProperties,
    WM_OT_Initiate,
    OBJECT_PT_GeneratorPanel,
    OBJECT_PT_Buttons,
    OBJECT_PT_InitialProperties,
    OBJECT_PT_SubductionProperties,
    WM_OT_AddSubFront,
    WM_OT_RunSubduction
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