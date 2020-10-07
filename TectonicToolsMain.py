import os

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

#==============================Properties ================================================================
class GeneratorProperties(PropertyGroup):

	#Properties that are used during the planet initialization
    radius: FloatProperty(name="Sphere Radius", default=6, min=0.1, max=10000)
    subDivs: IntProperty(name="Mesh Subdivisions", default=512, min=1, max=100000)
    initHeight: FloatProperty(name="Initial Noise Height", default=0.2, min=0.01, max=10000)
    noiseFreq: FloatProperty(name="Initial Noise Frequency", default=1, min=0.01, max=1000)
    
    #Identifiers so that we can retrieve relevant data throughout code
    planetID: StringProperty(name="Planet Identifier")
    subFrontsID: StringProperty(name="Subduction Front Identifier")

    deltaTime: FloatProperty(name="Time Step Size (dt)", default=0.1, min=0.001, max=10000)
    plateSpeed: FloatProperty(name="Plate Speed", default=1.0, min=0.0001, max=1000)
    plateauHeightLim: FloatProperty(name='Plateau Height Limit', default=1, min=0.0001, max=10000)
    influenceRangeMax: FloatProperty(name="Max Range of Influence", default=10, min=0.001, max=10000)
    influenceRangeMin: FloatProperty(name="Min Range of Influence", default=1, min=0.001, max=10000)
    baseUplift: FloatProperty(name='Base Subduction Uplift', default=0.5, min=0.0001, max=10000)

    
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
        baseUp = props.baseUplift

        #If not already done so, we begin by organizing our bezier curve collections
        #Collections is blender's directory system for objects in the scene,
        #And bezier curves are used in this code to represent subduction fronts drawn by the user
        defaultCollection, subFrontCollection = getCollections()
        splitBezierCurveToMesh(defaultCollection)
        frontsBM = getSubfrontBMeshes(subFrontCollection)
        
        #Create a bmesh object for the planet
        planet = bpy.data.objects.get(props.planetID)
        planetBM = bmesh.new()
        planetBM.from_mesh(planet.data)

        #We get and/or calculate custom vertex property layers
        #This allows us to save data on every vertex within blender
        frontId, distances, vertexSpeed = getCustomVertexLayers(props, planetBM, frontsBM)

        #A list of range of influences for each subduction front
        influenceRange = getSubductionRangeOfInfluence(planetBM, frontId, props, frontsBM)

        #forTestingPurposes(planetBM)

        #Main loop for subduction uplift
        for v in planetBM.verts:
            dist = v[distances]
            index = v[frontId]
            speed = v[vertexSpeed]
            infRange = influenceRange[index]
            distTrans = getDistanceTransfer(dist, infRange)
            heightTrans = getHeightTransfer(props, v)

            #Main equation as described in paper
            uplift = baseUp * distTrans * heightTrans * speed
            v = moveAlongRadialDirection(v, uplift)
        
        #Register changes made to blender
        planetBM.to_mesh(planet.data)
        planetBM.free()
        planet.select_set(True)
        return{"FINISHED"}


def forTestingPurposes(planetBM):
    r = [np.sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2) for v in planetBM.verts]
    average = sum(r) / len(r)
    print('{}, {}, {}'.format(min(r), max(r), average))



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

#For each subduction front, we create a quaternion to represent the plate velocity of nearby vertices
def getQuaternions(frontsBM, props):
    axii, quaternions = [], []
    for front in frontsBM:
        front.verts.ensure_lookup_table()
        start = front.verts[0].co
        end = front.verts[-1].co
        axes = end - start
        axii.append(axes)
        quaternions.append(M.Quaternion(axes, radians(props.plateSpeed)))
    return axii, quaternions

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

def getCustomVertexLayers(props, planetBM, frontsBM):
    if 'dist' not in planetBM.verts.layers.float.keys():
        frontId = planetBM.verts.layers.int.new('id')
        distances = planetBM.verts.layers.float.new('dist')
        vertexSpeed = planetBM.verts.layers.float.new('speed')
        calculateCustomVertexLayers(props, planetBM, frontsBM, frontId, distances, vertexSpeed)
    else:
        frontId = planetBM.verts.layers.int.get('id')
        distances = planetBM.verts.layers.float.get('dist')
        vertexSpeed = planetBM.verts.layers.float.get('speed')
    return frontId, distances, vertexSpeed

def calculateCustomVertexLayers(props, planetBM, frontsBM, frontId, distances, vertexSpeed):
    subFrontsKDTrees = getKDTrees(frontsBM)
    axes, quaternions = getQuaternions(frontsBM, props)

    #We calculate and set values for our custom vertex layers
    for v in planetBM.verts:

        #List of results from the kd.find() function for each front
        kdFinds = [list(kd.find(v.co)) for kd in subFrontsKDTrees]

        #We find the front index corresponding to the closest subduction front
        distants = [fnd[2] for fnd in kdFinds]
        dist, idx = min((dist, idx) for (idx, dist) in enumerate(distants))

        #The front ID is used to label which front a particular vertex belongs to
        v[frontId] = idx
        v[distances] = dist

        #Calculate the vertex speed and store as a vertex custom property layer
        closestFrontCo = (kdFinds[idx][0]).normalized()
        vertToAx = (v.co - axes[idx]).normalized()
        v[vertexSpeed] = closestFrontCo.dot(vertToAx) * props.plateSpeed







#========================Functions for Subduction Uplift ====================================

def getSpeedTransfer(props):
    pSpeed = props.plateSpeed


def getHeightTransfer(props, v):
    initH = props.radius
    maxHeight = props.plateauHeightLim
    heightAboveInit = sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2) - initH + 0.01
    normalizedHeight = heightAboveInit / maxHeight
    if normalizedHeight > 1.0:
        normalizedHeight = 1.0

    heightTransCurve = bpy.data.node_groups['ProfileCurves'].nodes['Height Transfer'].mapping
    curve = heightTransCurve.curves[3]

    if (normalizedHeight <= 1.0):
        curveSamplePoint = heightTransCurve.evaluate(curve, normalizedHeight)
    else:
        #pass
        curveSamplePoint = - heightTransCurve.evaluate(curve, normalizedHeight % 1)
    hTran = curveSamplePoint * maxHeight #* 4
    #print(hTran)
    return hTran


def moveAlongRadialDirection(v, dr):
    r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
    v.co.x += sin(theta) * cos(phi) * dr
    v.co.y += sin(theta) * sin(phi) * dr
    v.co.z += cos(theta) * dr
    return v

#Coordinate transformation function from cartesian to polar
def getPolarCoords(X, Y, Z):
    R = np.sqrt(X**2 + Y**2 + Z**2)
    Theta = np.arccos(Z / R)
    Phi = np.arctan2(Y, X)    
    return (R, Theta, Phi)


#Sample a point from the distance transfer curve editor      
def getDistanceTransfer(dist, infRange):
    distTransCurve = bpy.data.node_groups['ProfileCurves'].nodes['Distance Transfer'].mapping
    curve = distTransCurve.curves[3]

    samplePoint = dist / infRange
    if (samplePoint > 1):
        samplePoint = 1
        #print('Dist: {}, InfRange: {}'.format(dist, infRange))
    else:
        pass
        #print('Dist: {}, InfRange: {}'.format(dist, infRange))

    #distTransCurve.initialize()
    return distTransCurve.evaluate(curve, samplePoint)



#Depending on how high a mountain range is, the range of subduction influence needs to be calculated accordingly
#The larger a moutain range, the larger the influence
#We use the volume of the mountain range that is above the initial world radius to calculate the range of influence
def getSubductionRangeOfInfluence(planetBM, frontId, props, frontsBM):
    
    #Get relevant properties
    initH = props.radius
    maxRange = props.influenceRangeMax
    minRange = props.influenceRangeMin
    listLength = len(frontsBM)

    volumeAboveInitialHeight = [0 for i in range(listLength)]
    totalArea = [0 for i in range(listLength)]
    for f in planetBM.faces:
        index = f.verts[0][frontId]
        area = f.calc_area()
        heightDiff = sum([(v.co.x**2 + v.co.y**2 + v.co.z**2)**0.5 - initH for v in f.verts]) / len(f.verts)
        volumeAboveInitialHeight[index] += area * heightDiff
        totalArea[index] += area

    maxHeight = props.plateauHeightLim
    maxVolume = [area * maxHeight for area in totalArea]
    proportion = [(volumeAboveInitialHeight[i] / maxVolume[i])**0.6 for i in range(listLength)]
    print(proportion)
    influenceRange = [minRange + (maxRange - minRange) * proportion[i] for i in range(listLength)]
    return influenceRange


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
        layout.prop(props, "plateauHeightLim")
        layout.prop(props, "influenceRangeMax")
        layout.prop(props, "influenceRangeMin")
        layout.prop(props, "baseUplift")

        #Node groups are used to store profile curves, 
        #and profile curves are used to allow the user to custimize functions within the code (Eg. Distance Transfer) 
        if 'ProfileCurves' not in bpy.data.node_groups:
        	bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
        curveTree = bpy.data.node_groups['ProfileCurves'].nodes
        
        #Curves for custimizing the distance transfer function
        layout.label(text="Subduction Uplift Distance Transfer")
        if "Distance Transfer" in curveTree.keys():
        	layout.template_curve_mapping(curveTree["Distance Transfer"], "mapping")

        #Curves for custimizing the height transfer function
        layout.label(text="Height Transfer")
        if "Height Transfer" in curveTree.keys():
            layout.template_curve_mapping(curveTree["Height Transfer"], "mapping")


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

    if "Height Transfer" not in nods.keys():
        distCurves = nods.new('ShaderNodeRGBCurve')
        distCurves.name = "Height Transfer"
        pnts = distCurves.mapping.curves[3].points
        pnts[0].location = [0.0, 0.0]
        pnts[1].location = [0.5, 0.25]
        pnts.new(0.75, 1.0)
        pnts.new(0.8, 0.75)
        pnts.new(0.9, 0.125)
        pnts.new(1.0, 0.0)


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