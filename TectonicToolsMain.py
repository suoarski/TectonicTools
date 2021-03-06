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
    subDivs: IntProperty(name="Mesh Subdivisions", 
        description="Defines the resolution of the mesh",
        default=512, 
        min=1,
        max=100000)  

    meshSize: FloatProperty(
        name="Mesh Size", 
        description="Defines the size a flat terrain grid will be",
        default=20, 
        min=0.1, 
        max=100000)

    radius: FloatProperty(
        name="Sphere Radius",
        description="Defines the radius of the spherical planet",
        default=6, 
        min=0.1, 
        max=10000)

    initHeight: FloatProperty(
        name="Initial Noise Height", 
        description="How high the initial noise will be",
        default=0.25, 
        min=0.01, 
        max=10000)

    noiseFreq: FloatProperty(
        name="Initial Noise Frequency", 
        description="Defines how 'close' mountains will be",
        default=1, 
        min=0.01, 
        max=1000)

    isSphere: BoolProperty( 
        name="Spherical Terrain", 
        description="Use a spherical terrain",
        default=True)

    #Identifiers so that we can retrieve relevant data throughout code
    terrainID: StringProperty(
        name="Terrain Identifier",
        )

    #Properties to be used in subduction uplift
    time: FloatProperty(
        name="Time (Million Years)", 
        default=0, 
        min=0, 
        max=10000000)

    deltaTime: FloatProperty(
        name="Time Steps",
        description="Defines how much time passes at each iteration of the simulation",
        default=0.05, 
        min=0.0001, 
        max=10000)

    iterations: IntProperty(
        name="Iterations Per Click",
        description="Defines how many iterations the simulation will run for each time the user clicks Move Tectonic Plates",
        default=20, 
        min=1, 
        max=100000)

    baseUplift: FloatProperty(
        name="Base Uplift", 
        description = "Determines how fast mountains grows.",
        default=24.0, 
        min=0.0001, 
        max=100000
        )

    m1: FloatProperty(
        name="Slope Going Up From Front",
        description="Defines how steep the slope is near the subduction front (Given as a percent)",
        default=10,
        min=0.0001,
        max=1000
        )

    m2: FloatProperty(
        name="Slope Going Down", 
        description="Defines how steep the slope is on the other side of the mountain from the subduction front (Given as a percent)",
        default=2, 
        min=0.00001, 
        max=1000
        )

    plateauHeight: FloatProperty(
        name="Plateau Height",
        description="The height after which mountains begin to form plateaus (not to scale)",
        default=0.07, 
        min=0.01, 
        max=100000
        )  

    spreadRate: FloatProperty(
        name = "Spread Rate",
        description = "How fast a mountain range spreads from the subduction boundary",
        default = 0.3,
        min = 0,
        max = 10000000
        )

    numToAverage: IntProperty(
        name="Number of Edges to Average Over",
        description="How many edges to average over when we are calculating distances that vertices are from a subduction front",
        default=4, 
        min=1, 
        max=100000)

    degrowthConst: FloatProperty(
        name="Degrowth Constant",
        description = "Defines how fast a mountain collapses under its own weight",
        default=0.5, 
        min=0, 
        max=200000)

    degrowthThreshold: FloatProperty(
        name="Degrowth Threshold",
        description="The threshold height after which mountains will collapse under their own weight",
        default=0.5
        )

    #Properties for managing tectonic plates
    plateNums: IntProperty(
        name="Number Of Plates", 
        default=0, 
        min=0, 
        max=200)

    newPlateAxes: FloatVectorProperty(
        name="Plate Rotation Axes",
        description="The axis that the new tectonic plate will rotate about",
        )

    newPlateSpeed: FloatProperty(
        name="Plate Speed",
        description="How fast the new plate will move along the sphere. Units are in minutes per iteration (There are 60 minutes in a degree)", 
        default=0.02)

    newPlateIsContinental: BoolProperty(
        name="Is Continental Plate",
        description="Mountains will only form on continental plates",
        default=True)

    randomizeProperties: BoolProperty(
        name="Randomize Plate Axes", 
        default=True)

    movePlates: BoolProperty(
        name="Move Tectonic Plates",
        description="Specifies whether to actually move the plates along the sphere. If unchecked, the mountains will still grow, but the plates will remain stationary.",
        default=True)

    #Properties for the auto plate generator
    autoCreateIsSelected: BoolProperty(
        name="Auto Create Plates",
        description="Create tectonic plates automatically using Voronoi tesselation",
        default=True)

    numOfPlatesToCreate: IntProperty(
        name="Number of Plates", 
        description="Number of Plates to create using Voronoi tesselation",
        default=16)

    boundAmp: FloatProperty(
        name="Boundary Amplitude",
        description="Amplitude of the noise used to modulate the tectonic plate boundaries with",
        default=0.6,
        min=0.001,
        max=10000)

    boundFreq: FloatProperty(
        name="Boundary Frequency",
        description="Frequency of the noise used to modulate the tectonic plate boundaries with",
        default=1,
        min=0.001,
        max=10000)

    #This enum property needs to be identical to the one in the ANT Landscape code, so I copy and pasted it
    noiseType: EnumProperty(
        name="Noise Type",
        default='hetero_terrain',
        description="The type of noise to be used for the initial landscape",
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

#Properties that each plate instance will have
class WMTPlateProperties(PropertyGroup):
    rotationAxes: FloatVectorProperty(name="Plate Rotation Axes")
    rotationSpeed: FloatProperty(name="Rotation Speed", default=0.5)
    isContinental: BoolProperty(name="Is Continental", default=True)
    centroid: FloatVectorProperty(name="Plate Centroid")

#A collection of pointers that point towards each instance of WMTPlateProperties created
class WMTPlatePointers(PropertyGroup):
    ids: bpy.props.CollectionProperty(type=WMTPlateProperties)



#==============================Operators ==================================================================
#Operator for initializing the terrain
class WM_OT_Initiate(Operator):
    bl_label = "Initiate Terrain"
    bl_idname = "wm.initiate_terrain"
    bl_description = "Creates a new A.N.T Landscape object suitable for tectonic tools"
    def execute(self, context):
        scene = context.scene
        props = scene.properties
        initializeProfileCurves()

        #Use the A.N.T Landscape operator to add an initial landscape
        bpy.ops.mesh.landscape_add(
            refresh = True, 
            sphere_mesh = props.isSphere,
            subdivision_x = props.subDivs,
            subdivision_y = props.subDivs,
            mesh_size = 2 * props.radius,
            mesh_size_x = props.meshSize,
            mesh_size_y = props.meshSize,
            height = props.initHeight,
            maximum = props.initHeight,
            minimum = - props.initHeight,
            noise_size = props.noiseFreq,
            noise_type = props.noiseType
            )

        #Save the name of the newly created object and convert it to a BMesh
        terrain = bpy.context.active_object
        props.terrainID = terrain.name

        #Create a BMesh from the terrain data
        terrainBM = bmesh.new()
        terrainBM.from_mesh(terrain.data)

        #Use the BMesh to create a new plate which initially contains all the vertices
        props.newPlateIsContinental = False
        addNewPlate(props, terrain, terrainBM)
        props.newPlateIsContinental = True
        props.newPlateAxes = M.Vector((0, 1, 0))

        #Apply changes made to the mesh
        terrainBM.to_mesh(terrain.data)
        terrainBM.free()

        #Set up the blender environment so that the user can immediately start drawing his plates
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.wm.tool_set_by_id(name="builtin.select_lasso")
        bpy.ops.mesh.select_all(action='DESELECT')
        return {"FINISHED"}


class WM_OT_AutoPlates(Operator):
    bl_label = "Auto Create Plates"
    bl_idname = "wm.auto_plates"
    bl_description = "Automatically creates tectonic plates based on Voronoi tesselation"
    def execute(self, context):
        scene = context.scene
        props = scene.properties

        #Create a bmesh object for the terrain
        bpy.ops.object.mode_set(mode='EDIT')
        terrain = bpy.data.objects.get(props.terrainID)
        terrainBM = bmesh.from_edit_mesh(terrain.data)
        plateId = terrainBM.verts.layers.int.get('id')

        props.randomizeProperties = True
        centroids = [terrain.data.WMTPlatePointers.ids[0].centroid]

        for i in range(props.numOfPlatesToCreate):
            plate = terrain.data.WMTPlatePointers.ids.add()
            plate.rotationAxes = M.Vector((
                    120 * np.random.random() - 60,
                    120 * np.random.random() - 60,
                    120 * np.random.random() - 60)).normalized()

            plate.rotationSpeed = 2 * props.newPlateSpeed * np.random.random() / 60
            plate.isContinental = np.random.choice([True, False])

            cent = np.random.choice(terrainBM.verts).co
            plate.centroid = cent
            centroids.append(cent)

        centroidsKD = M.kdtree.KDTree(len(centroids))
        for i, cent in enumerate(centroids):
            centroidsKD.insert(cent, i)
        centroidsKD.balance()

        for v in terrainBM.verts:
            noise = M.noise.noise_vector(v.co * props.boundFreq) * props.boundAmp

            co, index, dist = centroidsKD.find(v.co + noise)
            v[plateId] = index

        bmesh.update_edit_mesh(terrain.data)
        return {"FINISHED"}



class WM_OT_AddPlate(Operator):
    bl_label = "Create Plate"
    bl_idname = "wm.add_plate"
    bl_description = "Creates a new plate based on currently selected vertices"
    def execute(self, context):
        scene = context.scene
        props = scene.properties

        #Create a bmesh object for the terrain
        bpy.ops.object.mode_set(mode='EDIT')
        terrain = bpy.data.objects.get(props.terrainID)
        terrainBM = bmesh.from_edit_mesh(terrain.data)

        #Add a new plate and update changes made
        addNewPlate(props, terrain, terrainBM)
        recalculatePlateCentroids(props, terrain, terrainBM)

        bmesh.update_edit_mesh(terrain.data)
        return{"FINISHED"}


class WM_OT_MovePlates(Operator):
    bl_label = "Move Tectonic Plates"
    bl_idname = "wm.move_plates"
    bl_description = "Move all tectonic plates and run the simulation for a few steps. Note that the first time running this operator will be the   ."
    def execute(self, context):
        scene = context.scene
        props = scene.properties

        #Create BMesh object
        bpy.ops.object.mode_set(mode='EDIT')
        terrain = bpy.data.objects.get(props.terrainID)
        terrainBM = bmesh.from_edit_mesh(terrain.data)

        #Get plate vertex layers and plate properties
        plateId = terrainBM.verts.layers.int.get('id')
        plateProps = terrain.data.WMTPlatePointers.ids
        isContinent = np.array([plateProps[v[plateId]].isContinental for v in terrainBM.verts])

        #We convert the heigh transfer function from a blender profile curve to numpy arrays
        heightTransCurve = bpy.data.node_groups['ProfileCurves'].nodes['Height Transfer'].mapping
        curve = heightTransCurve.curves[3]
        xTrans = np.linspace(0, 1, num=21, endpoint=True)
        yTrans = [heightTransCurve.evaluate(curve, i) for i in xTrans]

        #Used if terrain is a sphere
        rotationQuats = [M.Quaternion(i.rotationAxes, radians(i.rotationSpeed)) for i in plateProps]

        #Used if terrain is not a sphere
        plateSpeeds = [M.Vector(i.rotationAxes) * i.rotationSpeed for i in plateProps]
        for speed in plateSpeeds:
            speed[2] = 0

        #Convert vertex coordinates to numpy arrays for faster calculations
        XYZ = np.array([[i for i in v.co] for v in terrainBM.verts])
        X = XYZ[:, 0]
        Y = XYZ[:, 1]
        Z = XYZ[:, 2]

        #To avoid having to recalculate the vertices distances and speeds to the subduction fronts,
        #We save these parameters as custom vertex layers.
        dists, speeds = getCustomVertexLayers(props, terrainBM, plateProps, XYZ)

        #Main loop for doing subduction. To speed things up, all calculations are done in numpy
        for i in range(props.iterations):
            R = np.sqrt(np.sum(XYZ**2, axis=1))
            Theta = np.arccos(Z / R)
            Phi = np.arctan2(Y, X)
            
            props.time += props.deltaTime
            heights = R - props.radius

            distanceTransfer = plateauProfile(props, dists)
            heightTransfer = np.interp(heights, xTrans, yTrans)

            #dR = 0.1 * speeds * isContinent
            
            dR = props.baseUplift * distanceTransfer * heightTransfer * isContinent * speeds * props.deltaTime
            dR -= (heights > props.degrowthThreshold) * (props.degrowthConst * (heights - 0.3 * props.degrowthThreshold))**2

            X += np.sin(Theta) * np.cos(Phi) * dR
            Y += np.sin(Theta) * np.sin(Phi) * dR
            Z += np.cos(Theta) * dR

        #Apply changes made to Bmesh
        if props.isSphere:
            for i, v in enumerate(terrainBM.verts):
                v.co.x = X[i]
                v.co.y = Y[i]
                v.co.z = Z[i]

                if props.movePlates:
                    for i in range(props.iterations):
                        v.co = rotationQuats[v[plateId]] @ v.co
        else:
            for v in terrainBM.verts:
                v.co += plateSpeeds[v[plateId]] 

        bmesh.update_edit_mesh(terrain.data)
        return{"FINISHED"}



#We use a sigmoid function for the height transfer
def sigmoid(x, mean, spread):
    return 1 / (1 + np.exp( - (x - mean) / spread))


#Here we stitch a bunch of line equations together to define the distance transfer
def plateauProfile(props, x):

    t = props.time * props.spreadRate
    m1 = props.m1 / 100
    m2 = props.m2 / 100
    maxHeight = props.plateauHeight


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



def getCustomVertexLayers(props, terrainBM, plateProps, XYZ):
    if 'dist' not in terrainBM.verts.layers.float.keys():
        distances = terrainBM.verts.layers.float.new('dist')
        vertexSpeed = terrainBM.verts.layers.float.new('speed')
        calculateCustomVertexLayers(props, terrainBM, distances, vertexSpeed, plateProps, XYZ)
    else:
        distances = terrainBM.verts.layers.float.get('dist')
        vertexSpeed = terrainBM.verts.layers.float.get('speed')
        #calculateCustomVertexLayers(props, terrainBM, distances, vertexSpeed, plateProps, XYZ)

    dists = np.array([v[distances] for v in terrainBM.verts])
    speeds = np.array([v[vertexSpeed] for v in terrainBM.verts])
    return dists, speeds


def calculateCustomVertexLayers(props, terrainBM, distances, vertexSpeed, plateProps, XYZ):
    plateId = terrainBM.verts.layers.int.get('id')
    rotationQuats = [M.Quaternion(i.rotationAxes, radians(i.rotationSpeed)) for i in plateProps]
    rotationSpeeds = [i.rotationSpeed for i in plateProps]
    plateCentroids = [i.centroid for i in plateProps]
    
    #If the two vertices at the end of an edge are not on the same plate, we have found a plate boundary
    subFrontEdges, frontEdges = [], []
    for e in terrainBM.edges:
        if (e.verts[0][plateId] != e.verts[1][plateId]) and (np.random.random() < 1):

            vert0 = e.verts[0]
            vert1 = e.verts[1]
            quat0 = rotationQuats[vert0[plateId]]   
            quat1 = rotationQuats[vert1[plateId]]

            frontSpeed = ((quat0 @ vert0.co - quat1 @ vert1.co).length - (vert0.co - vert1.co).length)

            frontHeading = (quat0 @ vert0.co - vert0.co).normalized() - (quat1 @ vert1.co - vert1.co).normalized()
            edgeTangent = (vert1.co - vert0.co).normalized()
            dotProd = frontHeading.dot(edgeTangent)

            if (frontSpeed < 0) and (dotProd > 0.7):
                e.select = True
                subFrontEdges.append(e)

    kdTree = M.kdtree.KDTree(len(subFrontEdges))
    for i, e in enumerate(subFrontEdges):
        edgeCenter = (e.verts[0].co + e.verts[1].co) / 2
        kdTree.insert(edgeCenter, i)
    kdTree.balance()

    for v in terrainBM.verts:
        edgeCos, edgeID, distancess = [], [], []
        for (edgeCo, index, dist) in kdTree.find_n(v.co, props.numToAverage):
            edgeCos.append(edgeCo)
            edgeID.append(index)
            distancess.append(dist)


        plateQuat = rotationQuats[v[plateId]]
        axis, angle = plateQuat.to_axis_angle()
        axis = axis.normalized()
        axis = np.array(axis)

        cent = plateCentroids[v[plateId]]
        vCo = np.array(v.co)
        speed = np.linalg.norm((cent - vCo) - np.dot(cent - vCo, axis) * axis)


        fromCentToPoint = vCo - cent
        crossProd = np.cross(fromCentToPoint, axis)
        dotProd = np.dot(vCo, crossProd)
        speed = sigmoid(speed * dotProd, 0, 5)

        v[distances] = sum(distancess) / len(distancess)
        v[vertexSpeed] = speed 

        





#When adding a new plate, we create a new instance of WMTPlatePointers (plate properties)
#And update the vertex layer "plateId"
def addNewPlate(props, terrain, terrainBM):
    plate = terrain.data.WMTPlatePointers.ids.add()
    plate.rotationAxes = props.newPlateAxes
    plate.rotationSpeed = props.newPlateSpeed / 60
    plate.isContinental = props.newPlateIsContinental
    plate.centroid = np.random.choice(terrainBM.verts).co

    #The vertex layer "plateId" is used to identify which plate a vertex bellongs to
    if 'id' not in terrainBM.verts.layers.int.keys():
        plateId = terrainBM.verts.layers.int.new('id')
        for v in terrainBM.verts:
            v[plateId] = 0
    else:
        plateId = terrainBM.verts.layers.int.get('id')
        for v in terrainBM.verts:
            if v.select:
                v[plateId] = props.plateNums

    #We randomize the axii for the next plate to rotate about
    if props.randomizeProperties:
        props.newPlateAxes = M.Vector((120 * np.random.random() - 60,
                                120 * np.random.random() - 60,
                                120 * np.random.random() - 60)).normalized()

    #Increase the total number of plates
    props.plateNums += 1


def recalculatePlateCentroids(props, terrain, terrainBM):
    plateProps = terrain.data.WMTPlatePointers.ids
    plateId = terrainBM.verts.layers.int.get('id')
    coordSum = np.zeros((len(plateProps), 3))
    coordCount = np.zeros(len(plateProps))

    for v in terrainBM.verts:
        xyz = np.array([v.co.x, v.co.y, v.co.z])
        coordSum[v[plateId]] += xyz
        coordCount[v[plateId]] += 1
        v.select = False

    coordAves = coordSum.T / coordCount
    coordAves = coordAves.T

    for i, plate in enumerate(plateProps):
        plate.centroid = M.Vector((coordAves[i]))



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

        layout.operator("wm.initiate_terrain")

        if props.autoCreateIsSelected:
            layout.operator("wm.auto_plates")
        else:
            layout.operator("wm.add_plate")

        layout.operator("wm.move_plates")

class OBJECT_PT_InitialProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Initial World Properties"
    bl_options = {"DEFAULT_CLOSED"}
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        #layout.prop(props, "isSphere")
        if props.isSphere:
            layout.prop(props, "radius")
        else:
            layout.prop(props, "meshSize")

        layout.prop(props, "subDivs")
        layout.prop(props, "initHeight")
        layout.prop(props, "noiseFreq")
        layout.prop(props, "noiseType")

class OBJECT_PT_PlateProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Properties For New Plates"
    bl_options = {"DEFAULT_CLOSED"}
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        layout.prop(props, "autoCreateIsSelected")

        if props.autoCreateIsSelected:
            layout.prop(props, "numOfPlatesToCreate")
            layout.prop(props, "boundAmp")
            layout.prop(props, "boundFreq")
            layout.prop(props, "newPlateSpeed")
            
        else:
            layout.prop(props, "randomizeProperties")
            layout.prop(props, "newPlateAxes")
            layout.prop(props, "newPlateSpeed")
            layout.prop(props, "newPlateIsContinental")
            layout.label(text="Total Number Of Plates: {}".format(props.plateNums))
        


class OBJECT_PT_SubductionProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Simulation Properties"  
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        layout.prop(props, "baseUplift")
        layout.label(text="Time Elapsed: {:.2f}".format(props.time))
        layout.prop(props, "deltaTime")
        layout.prop(props, "iterations")
        layout.prop(props, "movePlates")

        layout.label(text="Uplift Properties")
        layout.prop(props, "m1")
        layout.prop(props, "m2")
        layout.prop(props, "plateauHeight")
        layout.prop(props, "degrowthConst")
        layout.prop(props, "degrowthThreshold")
        layout.prop(props, "spreadRate")
        layout.prop(props, "numToAverage")
        

        #Node groups are used to store profile curves, 
        #and profile curves are used to allow the user to custimize functions within the code (Eg. Distance Transfer) 
        if 'ProfileCurves' not in bpy.data.node_groups:
            bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
        curveTree = bpy.data.node_groups['ProfileCurves'].nodes

        #Curves for custimizing the height transfer function
        layout.label(text="Height Transfer")
        if "Height Transfer" in curveTree.keys():
            layout.template_curve_mapping(curveTree["Height Transfer"], "mapping")


def initializeProfileCurves():
    if 'ProfileCurves' not in bpy.data.node_groups:
        bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
    
    nods = bpy.data.node_groups['ProfileCurves'].nodes
    if "Height Transfer" not in nods.keys():
        distCurves = nods.new('ShaderNodeRGBCurve')
        distCurves.name = "Height Transfer"
        pnts = distCurves.mapping.curves[3].points
        pnts[0].location = [0.0, 0.0]
        pnts[1].location = [0.2, 0.1]
        pnts.new(0.55, 0.75)
        pnts.new(1.0, 1.0)

#=============================Register Classes to Blender ================================================
classes = (
    GeneratorProperties,
    WMTPlateProperties,
    WMTPlatePointers,
    WM_OT_Initiate,
    WM_OT_AutoPlates,
    WM_OT_AddPlate,
    WM_OT_MovePlates,
    OBJECT_PT_GeneratorPanel,
    OBJECT_PT_InitialProperties,
    OBJECT_PT_PlateProperties,
    OBJECT_PT_SubductionProperties
)

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.Scene.properties = PointerProperty(type=GeneratorProperties)
    bpy.types.Mesh.WMTPlatePointers = PointerProperty(type=WMTPlatePointers)

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    del bpy.types.Scene.properties
    del bpy.types.Mesh.WMTPlatePointers


if __name__ == "__main__":
    register()