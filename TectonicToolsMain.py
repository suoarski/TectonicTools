import bpy
import bmesh
import mathutils as M
from math import (acos, atan2, cos, sin, sqrt, radians)
import numpy as np

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


class GeneratorProperties(PropertyGroup):

    #Properties for the initial terrain. These will be passed onto ANT Landscape
    subDivs: IntProperty(name="Mesh Subdivisions", default=1024, min=1, max=100000)
    meshSize: FloatProperty(name="Mesh Size", default=20, min=0.1, max=100000)
    initHeight: FloatProperty(name="Initial Noise Height", default=0.2, min=0.01, max=10000)
    noiseFreq: FloatProperty(name="Initial Noise Frequency", default=1, min=0.01, max=1000)

    #Identifiers so that we can retrieve relevant data throughout code
    terrainID: StringProperty(name="Terrain Identifier")

    #Properties to be used in subduction uplift
    time: FloatProperty(name="Time (Million Years)", default=0, min=0, max=10000000)
    deltaTime: FloatProperty(name="Time Steps (dt)", default=0.05, min=0.0001, max=10000)
    iterations: IntProperty(name="Iterations Per Click", default=5, min=1, max=100000)

    baseUplift: FloatProperty(name="Base Subduction Uplift", default=10.0, min=0.0001, max=100000)
    plateSpeed: FloatProperty(name="Plate Speed", default=0.2, min=0.0001, max=1000)
    m1: FloatProperty(name="Slope Going Up From Front", default=1, min=0.01, max=1000)
    m2: FloatProperty(name="Slope Going Down", default=0.2, min=0.01, max=1000)
    plateauHeight: FloatProperty(name="Plateau Height", default=1.4, min=0.01, max=100000)

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

        bpy.ops.mesh.landscape_add(
            refresh = True, 
            subdivision_x = props.subDivs,
            subdivision_y = props.subDivs,
            mesh_size_x = props.meshSize,
            mesh_size_y = props.meshSize,
            height = props.initHeight,
            maximum = props.initHeight,
            minimum = - props.initHeight,
            noise_size = props.noiseFreq,
            noise_type = props.noiseType
            )

        terrain = bpy.context.active_object
        props.terrainID = terrain.name
        return {"FINISHED"}


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

        #If not already done so, we begin by organizing our bezier curve collections
        #Collections is blender's directory system for objects in the scene,
        #And bezier curves are used in this code to represent subduction fronts drawn by the user
        defaultCollection, subFrontCollection = getCollections()
        splitBezierCurveToMesh(defaultCollection)
        frontsBM = getSubfrontBMeshes(subFrontCollection)

        #Create a bmesh object for the terrain
        terrain = bpy.data.objects.get(props.terrainID)
        terrainBM = bmesh.new()
        terrainBM.from_mesh(terrain.data)

        #Create a numpy array of vertex coordinates
        #I will try to do all my major calculations using numpy, 
        #as it is significantly faster than using loops, and can be easily modified to use GPU acceleration

        XYZ = np.array([[i for i in v.co] for v in terrainBM.verts])

        #We use custom vertex layers to save data in blender and avoid having to recalculate these values
        #each time this operator is called. This function returns data as np arrays
        frontID, dists, speedTransfer = getCustomVertexLayers(props, terrainBM, frontsBM)

        for i in range(props.iterations):
            props.time += props.deltaTime
            distanceTransfer = plateauProfile(props.time, dists, props.m1, props.m2, props.plateauHeight)
            heightTransfer = sigmoid(XYZ[:, 2], 1.0, 0.3)
            XYZ[:, 2] += props.baseUplift * speedTransfer * heightTransfer * distanceTransfer * props.deltaTime
            XYZ[:, 2] -= (0.2 * XYZ[:, 2] / props.plateauHeight)**2

        for i, v in enumerate(terrainBM.verts):
            v.co.z = XYZ[i, 2]

        '''
        bpy.ops.object.mode_set(mode='EDIT')
        terrain.select_set(True)

        frontId = terrainBM.verts.layers.int.get('id')
        vertexSpeed = terrainBM.verts.layers.float.get('speed')

        for i, v in enumerate(terrainBM.verts):
            v.co.z += speeds[i] / (1 + dists[i])
        
        '''

        #Register changes made to blender
        terrain.select_set(True)
        bpy.ops.object.mode_set(mode='OBJECT')
        terrainBM.to_mesh(terrain.data)
        terrainBM.free()

        return {"FINISHED"}



#We use a sigmoid function for the height transfer
def sigmoid(x, mean, spread):
    return 1 / (1 + np.exp( - (x - mean) / spread))


#Here we stitch a bunch of line equations together to define the distance transfer
def plateauProfile(t, x, m1, m2, maxHeight):
    m2 = np.abs(m2)
    y = np.zeros(x.shape)
    if (m1 * t > maxHeight):
        isFirstPartOfCurve = (x * m1 < maxHeight)
        firstPart = x * m1 * isFirstPartOfCurve
        
        plateauLength = m1 * t - maxHeight
        isSecondPartOfCurve = (x * m1 >= maxHeight) * (x * m1 < (maxHeight + plateauLength))
        secondPart = np.ones(x.shape) * maxHeight * isSecondPartOfCurve
        
        xTransition = (maxHeight + plateauLength) / m1
        isThirdPartOfCurve = (x * m1 >= (maxHeight + plateauLength))
        thirdPart = ( - x * m2 + m2 * xTransition + maxHeight) * isThirdPartOfCurve
        
        y = firstPart + secondPart + thirdPart
        y *= (y >= 0)
        
    else:
        isFirstPartOfCurve = (x < t)
        firstPart = x * m1 * isFirstPartOfCurve
        
        isSecondPart = (1 - isFirstPartOfCurve)
        secondPart = ( - m2 * x + t * m1 + m2 * t) * isSecondPart
        
        y = firstPart + secondPart
        y *= (y >= 0)
    return y



#=========================Functions for Setting Up Subduction ==================================

#Collections are a directory like structure for objects within the scene
#We check if the SubFronts collection already exists, and create one otherwise
def getCollections():
    bpy.ops.object.mode_set(mode='OBJECT')
    defaultCollection = bpy.data.collections["Collection"]
    if "SubFronts" not in bpy.data.collections:
        bpy.ops.collection.create(name="SubFronts")
        subFrontCollection = bpy.data.collections["SubFronts"]
        bpy.context.scene.collection.children.link(subFrontCollection)
    else:
        subFrontCollection = bpy.data.collections["SubFronts"]
    return (defaultCollection, subFrontCollection)

#Assuming that the subduction front object is selected and is still of type CURVE,
#We convert it to a MESH and seperate the MESH by loose parts
def splitBezierCurveToMesh(defaultCollection):
    selectedObjects = bpy.context.selected_objects
    for obj in selectedObjects:
        if obj.type == 'CURVE':
            if (obj.name in defaultCollection.objects.keys()):
                defaultCollection.objects.unlink(obj)
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.convert(target='MESH')
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.separate(type='LOOSE')
    for obj in selectedObjects:
        obj.select_set(False)


#Create a list of bmesh objects for each subduction front
#Bmesh allows for mesh manipulations within python code
def getSubfrontBMeshes(subFrontCollection):
    frontsBM = []
    subFronts = subFrontCollection.objects
    for front in subFronts:
        newBM = bmesh.new()
        newBM.from_mesh(front.data)
        frontsBM.append(newBM)
    return frontsBM

#For each subfuction front, we create a kd-tree for faster spacial searches
def getKDTrees(frontsBM):
    subFrontsKDTrees = []
    for front in frontsBM:
        kd = M.kdtree.KDTree(len(front.verts))
        for i, v in enumerate(front.verts):
            kd.insert(v.co, i)
        kd.balance()
        subFrontsKDTrees.append(kd)
    return subFrontsKDTrees

def getAxii(frontsBM):
    axii = []
    for front in frontsBM:
        front.verts.ensure_lookup_table()
        axii.append(front.verts[-1].co - front.verts[0].co)
    return axii

def getCustomVertexLayers(props, terrainBM, frontsBM):
    if 'dist' not in terrainBM.verts.layers.float.keys():
        frontId = terrainBM.verts.layers.int.new('id')
        distances = terrainBM.verts.layers.float.new('dist')
        vertexSpeed = terrainBM.verts.layers.float.new('speed')
        calculateCustomVertexLayers(props, terrainBM, frontsBM, frontId, distances, vertexSpeed)
    else:
        frontId = terrainBM.verts.layers.int.get('id')
        distances = terrainBM.verts.layers.float.get('dist')
        vertexSpeed = terrainBM.verts.layers.float.get('speed')

    #Create numpy arrays for each vertex layer
    ids = np.array([v[frontId] for v in terrainBM.verts])
    dists = np.array([v[distances] for v in terrainBM.verts])
    speeds = np.array([v[vertexSpeed] for v in terrainBM.verts])

    return ids, dists, speeds

def getSubFrontTangents(frontsBM):
    tangents = []
    for front in frontsBM:
        tangs = []
        for i in range(len(front.verts) - 1):
            tangs.append(front.verts[i+1].co - front.verts[i].co)
        tangs.append(tangs[-1])
        tangents.append(tangs)
    return tangents


def calculateCustomVertexLayers(props, planetBM, frontsBM, frontId, distances, vertexSpeed):
    subFrontsKDTrees = getKDTrees(frontsBM)
    axes = getAxii(frontsBM)
    tangents = getSubFrontTangents(frontsBM)

    #We calculate and set values for our custom vertex layers
    for v in planetBM.verts:

        #List of results from the kd.find() function for each front
        kdFinds = [list(kd.find(v.co)) for kd in subFrontsKDTrees]

        #We get the distants and ID of the closest vert in subfronts
        distants = [fnd[2] for fnd in kdFinds]
        closestVertsID = [fnd[1] for fnd in kdFinds]

        #For each subduction front, we check if v is on the left hand side of it
        isLeft, dots = [], []
        for i, vertID in enumerate(closestVertsID):
            fromFrontToV = frontsBM[i].verts[vertID].co - v.co
            crossWithTangents = fromFrontToV.cross(tangents[i][vertID]).normalized()
            dot = crossWithTangents.dot(v.normal)
            isLef = (dot >= 0)
            dots.append(dot)
            isLeft.append(isLef)

        distants = [dis for i, dis in enumerate(distants) if isLeft[i]]

        if not distants:
            v[frontId] = -1
            v[distances] = 10000000
            v[vertexSpeed] = 0
        else:
            dist, idx = min((dist, idx) for (idx, dist) in enumerate(distants))
            closestFrontCo = (kdFinds[idx][0]).normalized()
            vertID = closestVertsID[idx]

            #The front ID is used to label which front a particular vertex belongs to
            v[frontId] = idx
            v[distances] = dist
            v[vertexSpeed] = np.sin(np.arccos(closestFrontCo.dot(axes[idx].normalized()))) * props.plateSpeed




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
    bl_label = "Initial World Properties"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        layout.prop(props, "meshSize")
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

        layout.prop(props, "baseUplift")
        layout.label(text="Time Elapsed: {:.2f}".format(props.time))
        layout.prop(props, "deltaTime")
        layout.prop(props, "iterations")
        layout.prop(props, "plateSpeed")

        layout.label(text="Used for Distance Transfer")
        layout.prop(props, "m1")
        layout.prop(props, "m2")
        layout.prop(props, "plateauHeight")

        layout.label(text="Used For Height Transfer")

#=============================Register Classes to Blender ================================================
classes = (
    GeneratorProperties,
    WM_OT_Initiate,
    WM_OT_AddSubFront,
    WM_OT_RunSubduction,
    OBJECT_PT_GeneratorPanel,
    OBJECT_PT_Buttons,
    OBJECT_PT_InitialProperties,
    OBJECT_PT_SubductionProperties
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