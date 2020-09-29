    
import bpy
import os
from bpy.props import (StringProperty, BoolProperty, IntProperty, FloatProperty, FloatVectorProperty, EnumProperty, PointerProperty)
from bpy.types import (Panel, Menu, Operator, PropertyGroup)
from math import (acos, atan2, cos, sin, sqrt, radians)

import bmesh
import mathutils as M
from mathutils.bvhtree import BVHTree
import numpy as np
from random import random, choice
import pdb

import cProfile


#================================Properties =============================================
#Note that blender deletes all python data after it has finished running an operator
#As such, global variables must be stored within blender and not as python variablesd

#Blender also converts python class definitions into some kind of C++ version of it,
#As such, class names in python need appropriate prefixes
class GeneratorProperties(PropertyGroup):
    
    #Initial planet properties
    subDivs: IntProperty(name="Sphere Subdivisions", default=5, min=1, max=50)
    radius: FloatProperty(name="Sphere Radius", default=20, min=0.1, max=10000)
    lacunarity: FloatProperty(name="Initial Lacunarity", default=2, min=0, max=100)
    noiseHeight: FloatProperty(name="Initial Noise Height", default=0.02, min=0, max=10000)
    plateNums: IntProperty(name="Number Of Plates", default=0, min=0, max=200)
    plateSpeed: FloatProperty(name="Default Plate Speed", default=0.5, min=0, max=1000)
    oceanHeight: FloatProperty(name="Ocean Height", default=0.1, min=0, max=10000)
    maxHeight: FloatProperty(name="Max Mountain Height", default=2, min=0, max=1000)
    
    #For the Voronoi plate generator
    nOfPlates: IntProperty(name="Number of Plates", default=24, min=0, max=1000)
    boundLac: FloatProperty(name="Lacunarity of Boundary", default=0.5, min=0, max=1000)
    boundNoiseAmp: FloatProperty(name="Amplitude of Boundary Noise", default=4, min=0, max=1000)
    
    #For the subduction and continental collision algorithms
    baseUplift: FloatProperty(name="Base Subduction Uplift", default=0.1, min=0, max=10000)
    maxUpliftRange: FloatProperty(name="Maximum Uplift Range", default=9, min=0, max=10000)
    minUpliftRange: FloatProperty(name="Minimum Uplift Range", default=3, min=0, max=10000)
    
    discreteCollisionCoef: FloatProperty(name="Discrete Collision Coefficient", default=0.4, min=0, max=10000)
    contColMax: FloatProperty(name="Continental Collision Radius", default=6, min=0, max=10000)
    
    #For the oceanic crust generator algorithm
    aveCrustRad: FloatProperty(name="Average Ridge Width", default=0.3, min=0, max=100000)
    crustNoiseAmp: FloatProperty(name="Crust Noise Amplitude", default=0.1, min=0, max=1)
    crustNoiseLac: FloatProperty(name="Crust Noise Lacunarity", default=2, min=0, max=100000)
    ridgeDepth: FloatProperty(name="Ridge Depth", default=0.5, min=0, max=10000)
    
    #For animations
    totalFrames: IntProperty(name="Number of Frames", default=6, min=1, max=100000)
    remeshAfter: IntProperty(name="Remesh After n Frames", default=2, min=1, max=100000)
    defaultOutputLocation = str(os.getcwd) + "\TectonicToolsAnimationOutputLocation"
    outputFolder: StringProperty(name="Output Folder", default=(os.getcwd() + "\TectonicToolsAnimationOutputLocation"), maxlen=1024)



    #Other
    contErosion: FloatProperty(name="Continental Erosion Coefficient", default=1, min=0, max=1000)

#Properties that each plate instance will have
class WMTPlateProperties(PropertyGroup):
    axes: FloatVectorProperty(name="Rotation Axes")
    angle: FloatProperty(name="Rotation Angle", default=0.1, min=0, max=1000)
    isContinental: BoolProperty(name="Is Continental")
    plateCentroid: FloatVectorProperty(name="Plate Centroid")
    plateMaterial: PointerProperty(name="PlateMaterial", type=bpy.types.Material)
    averagePlateHeight: FloatProperty(name="Average Plate Height", default=0, min=-1000, max=1000)

#A collection of pointers that point towards each instance of WMTPlateProperties created
class WMTPlatePointers(PropertyGroup):
    ids: bpy.props.CollectionProperty(type=WMTPlateProperties)







#=====================================Operators ==============================================
#Operators define behaviour for tools    

#Operator for initializing the terrain
class WM_OT_Initiate(Operator):
    bl_label = "Initiate Terrain"
    bl_idname = "wm.initiate_terrain"
    
    #Main function that gets called by operators
    def execute(self, context):
        scene = context.scene
        props = scene.properties
        
        #We create a new sphere object
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=props.subDivs, radius=props.radius)
        
        #We get the mesh data from plane and convert to bmesh for easier mesh manipulation
        sphere = bpy.context.object.data
        bm = bmesh.new()
        bm.from_mesh(sphere)
        
        #Custom vertex data layer for identifying which plate a vertex belongs to
        plateId = bm.verts.layers.int.new('id')
        isTerrane = bm.verts.layers.int.new('isTer')
        props.plateNums = 0
        createNewPlate(sphere, props, bm)
        
        #Apply initial terrain
        applyNoiseFunctions(bm, props)
        initializeProfileCurves()
        
        #We apply changes made to the sphere
        bpy.ops.object.mode_set(mode='OBJECT')
        bm.to_mesh(sphere)
        bm.free()
        bpy.ops.object.mode_set(mode='EDIT')
        
        return{"FINISHED"}
                    

#Operator for moving the selected plate
class WM_OT_MoveSelectedPlate(Operator):
    bl_label = "Move All Plates"
    bl_idname = "wm.move_selection"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        bpy.ops.object.mode_set(mode='EDIT')
        movePlates(scene, props, plateProps)
        
        #Apply changes made to blender
        bmesh.update_edit_mesh(sphere)
        return{"FINISHED"}

#A bunch of properties and lists that we will need throughout our simulation

def movePlates(scene, props, plateProps):
    sphere = bpy.data.objects.get('Icosphere')
    try:
        bpy.ops.object.mode_set(mode='EDIT')
    except:
        print("exception")
    bm = bmesh.from_edit_mesh(sphere.data)

    bm.verts.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    plateId = bm.verts.layers.int.get('id')
    isTerrane = bm.verts.layers.int.get('isTer')

    boundingEdges = getPlateBoundaries(bm)
    plateQuaternions = [M.Quaternion((r.axes[0], r.axes[1], r.axes[2]), radians(r.angle)) for r in plateProps.ids]
    centroids = [M.Vector((r.plateCentroid[0], r.plateCentroid[1], r.plateCentroid[2])) for r in plateProps.ids]
    isContinent = [r.isContinental for r in plateProps.ids]
    maxSpeed = max([np.abs(((plateQuaternions[v[plateId]] @ v.co) - v.co).length) for v in bm.verts])
    subductionFronts = getSubFronts(boundingEdges, plateQuaternions, plateId, props.plateNums, isContinent)
    subFronts = list(filter(lambda x: isSubduction(x, plateQuaternions, plateId, isContinent), boundingEdges))
    oceanContEdges = list(filter(lambda x: isOceanicContEdge(x, plateId, isContinent), boundingEdges))
    
    #Run the main function for applying tectonic phenomina
    applySubductionUplift(bm, props, subductionFronts, subFronts, plateId, plateProps, plateQuaternions, maxSpeed)
    applySlabPull(subductionFronts, plateProps, isContinent, centroids)
    applyContinentalCollision(bm, props, plateId, isTerrane, plateQuaternions, centroids, maxSpeed)
    #erodeContinents(bm, props, isContinent, plateId)
    
    #Move selected plate and centroids by specified quaternion
    bm.verts.ensure_lookup_table()
    for v in bm.verts:
        v.co = plateQuaternions[v[plateId]] @ v.co
    moveCentroids(plateQuaternions, centroids, plateProps)

#Operator for automatically generating tectonic plates based on voronoi teselations
class WM_OT_CreatePlates(Operator):
    bl_label = "Create Voronoi Plates"
    bl_idname = "wm.voronoi"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        centroids = []
        
        #Create n new plates and store their centroid location
        for i in range(props.nOfPlates):
            createNewPlate(sphere, props, bm)
            cent = plateProps.ids[i].plateCentroid
            centroids.append(M.Vector((cent[0], cent[1], cent[2])))
            props.plateNums += 1
        
        calculatePlateId(bm, centroids, props, plateId)
        initializeTerranes(bm, props, plateId, isTerrane, plateProps)
        
        return{"FINISHED"}


#We apply colors to the mesh's vertex color layer, which we then convert into a material (which the end user may edit however they like)
class WM_OT_ApplyVertexColor(Operator):
    bl_label = "Apply Vertex Colors"
    bl_idname = "wm.applyvertexcolors"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        
        #Apply vertex colors for alternative colouring scheme
        boundingEdges = getPlateBoundaries(bm)
        colorLayer = bm.loops.layers.color.new("color")
        green = (0, 1, 0, 1)
        blue = (0, 0, 1, 1)
        red = (1, 0, 0, 1)
        
        #Give each face a color based on if it is a terrain or not
        for f in bm.faces:
            for loop in f.loops:
                bool = loop.vert[isTerrane]
                loop[colorLayer] = green if bool else blue
        
        #Draw a line on plate boundaries
        for e in boundingEdges:
            for loop in e.link_loops:
                loop[colorLayer] = red
        
        #Required for assigning material to sphere
        bpy.ops.object.mode_set(mode='EDIT')
        sphre = bpy.data.objects['Icosphere']
        
        for v in bm.verts:
            v.select = True
        
        #Create a new material slot named "VertexColors" if it does not already exist
        if 'VertexColors' not in sphre.material_slots:
            vertexMaterial = bpy.ops.object.material_slot_add()
        
        #We create a new material and assign it to our sphere's material_slots
        mat = bpy.data.materials.new(name="VertexColors")
        sphre.material_slots[0].material = mat
        mat.use_nodes = True
        
        #Objects from blender for editing shader nodes trees
        tree = mat.node_tree
        nodes = tree.nodes
        links = tree.links
        
        #We link our sphere's vertex colors to the principled BSDF shader
        vertColors = nodes.new(type='ShaderNodeVertexColor')
        principledBSDF = nodes['Principled BSDF']
        links.new(vertColors.outputs['Color'], principledBSDF.inputs['Base Color'])
        
        #Apply changes made to mesh
        bmesh.update_edit_mesh(sphere)
        return{"FINISHED"}

#Creates two coloring schemes for our planet, the first creates a new material for each plate,
#The other makes use of vertex paint tools
class WM_OT_ApplyColor(Operator):
    bl_label = "Apply Plate Colors"
    bl_idname = "wm.applycolors"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        
        #Apply materials to faces of the plates
        materials = [r.plateMaterial for r in plateProps.ids]
        for f in bm.faces:
            f.material_index = f.verts[0][plateId]
        
        return{"FINISHED"}
                
#We use use remesh modifier to create a brand new mesh
#Any custom vertex data layer such as plateId needs to be recalculated and plate boundaries may have changed
class WM_OT_Remesh(Operator):
    bl_label = "Remesh"
    bl_idname = "wm.remesh"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        
        #We apply the remesh modififier to our sphere
        bpy.ops.object.mode_set(mode='OBJECT')
        sphre = bpy.context.object
        rmesh = sphre.modifiers.new(name="Remesh", type="REMESH")
        rmesh.octree_depth = props.subDivs
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Remesh")
        bpy.ops.object.mode_set(mode='EDIT')
        
        #We recalculate data that was lost due to remesh modifier
        sphere = bpy.context.object.data
        bm = bmesh.from_edit_mesh(sphere)
        plateId = bm.verts.layers.int.new('id')
        isTerrane = bm.verts.layers.int.new('isTer')
        centroids = [M.Vector((r.plateCentroid[0], r.plateCentroid[1], r.plateCentroid[2])) for r in plateProps.ids[:-1]]
        calculatePlateId(bm, centroids, props, plateId)
        initializeTerranes(bm, props, plateId, isTerrane, plateProps)
        
        #Apply Oceanic Crust
        applyOceanicCrustGeneration(bm, props, plateId, isTerrane, plateProps)
        
        #Update changes made to blender
        bmesh.update_edit_mesh(sphere)
        return{"FINISHED"}

#We use the lasso tool to select an area and press the "Split Plates" button to define a new tectonic plate
#This operator may not be compatible with remeshing
class WM_OT_SplitPlates(Operator):
    bl_label = "Split Plates"
    bl_idname = "wm.split_plates"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        
        functionToMeasurePerformanceWith(bm)
        
        #Apply changes made to mesh
        bmesh.update_edit_mesh(sphere)
        return{"FINISHED"}
















class WM_OT_Animate(Operator):
    bl_label = "Animate"
    bl_idname = "wm.animate"
    bl_context = "mesh_edit"
    def execute(self, context):

        #We initialize the sphere if it has not yet been done so
        if bpy.data.objects.get('Icosphere') is None:
            bpy.ops.wm.initiate_terrain()
            bpy.ops.wm.voronoi()
            bpy.ops.wm.applyvertexcolors()
        
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)

        #Paramateres to be used
        totFrames = props.totalFrames
        iterations = props.remeshAfter
        totalRemeshes = totFrames // iterations
        remainingFrames = totFrames % iterations
        
        #Run the simulation
        for r in range(totalRemeshes):
            for i in range(iterations):
                movePlates(scene, props, plateProps)
            bpy.ops.wm.remesh()
            bpy.context.view_layer.update()
            bpy.ops.wm.applyvertexcolors()
            renderImage(props)
        
        return {"FINISHED"}

#Render image and save it
def renderImage(props):
    bpy.ops.object.mode_set(mode='OBJECT')

    #Check if there already is a camera object in the scene, otherwise create a new camera
    if  bpy.data.objects.get('Camera') is None:
        createCamera()

    saveImage(props)
    setToEditMode()

def saveImage(props):
    i = 0
    outFold = props.outputFolder
    imageLoc = "{}\{}.png".format(outFold, i)
    while (os.path.isfile(imageLoc)):
        i += 1
        imageLoc = "{}\{}.png".format(outFold, i)
    bpy.context.scene.render.filepath = imageLoc
    bpy.ops.render.render(write_still=True)



def createCamera():
    bpy.ops.object.camera_add()
    camera = bpy.data.objects['Camera']
    camera.location = M.Vector((96, 19, 41))
    camera.rotation_euler = M.Vector((1.171, 0.013, 1.756))

def setToEditMode():
    icoSphere = bpy.data.objects.get('Icosphere')
    if icoSphere is not None:
        icoSphere.select_set(True)    
        bpy.context.view_layer.objects.active = icoSphere
        bpy.ops.object.mode_set(mode='EDIT')


        

#Whenever I have a new idea for this project, I use this operator to play and experiment with those ideas
class WM_OT_Testing(Operator):
    bl_label = "Used For Testing"
    bl_idname = "wm.testing"
    bl_context = "mesh_edit"
    def execute(self, context):
        scene, props, sphere, bm, plateProps, plateId, isTerrane = operatorSetup(context)
        
        
        
        bmesh.update_edit_mesh(sphere)
           
        return {"FINISHED"}

#These coordinate transformations are slow, I need to find a faster method
def functionToMeasurePerformanceWith(bm):
    for v in bm.verts:
        r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
        x, y, z = getCartesianCoords(r, theta, phi)
        
        rThetaPhi = cartToPolarTransformation(v.co)
        XYZ = polarToCartTransformation(rThetaPhi)

#Coordinate transformation from cartesian to polar
def cartToPolarTransformation(XYZ):
    x, y, z = XYZ.x, XYZ.y, XYZ.z
    r = sqrt(x**2 + y**2 + z**2)
    rThetaPhi = M.Vector((r, acos(z / r), atan2(y, x))) 
    return rThetaPhi

#Coordinate transformation function from polar to cartesian
def polarToCartTransformation(rThetaPhi):
    r, theta, phi = rThetaPhi.x, rThetaPhi.y, rThetaPhi.z
    sinTheta = sin(theta)
    XYZ = M.Vector((r*sinTheta*cos(phi), r*sinTheta*sin(phi), r*cos(theta)))
    return (XYZ)


#=================================Initialization and Remeshing =======================================

#Standard setup that I keep using in all my operators
def operatorSetup(context):
    scene = context.scene
    props = scene.properties

    setToEditMode()
    sphere = bpy.data.objects.get('Icosphere').data
    bm = bmesh.from_edit_mesh(sphere)
    plateId = bm.verts.layers.int.get('id')
    isTerrane = bm.verts.layers.int.get('isTer')
    plateProps = context.object.data.WMTPlatePointers
    return (scene, props, sphere, bm, plateProps, plateId, isTerrane)

#Function for creating new plate instance
def createNewPlate(sphere, props, bm):
    if hasattr(bm.verts, "ensure_lookup_table"): 
                bm.verts.ensure_lookup_table()
    
    #Create a new plateProps instance and set variables
    plate = sphere.WMTPlatePointers.ids.add()
    plate.axes = M.Vector((120 * random() - 60, 120 * random() - 60, 120 * random() - 60))
    plate.angle = props.plateSpeed
    plate.isContinental = choice([True, False, False])
    plate.plateCentroid = choice(bm.verts).co
    
    #Create a material for the plate with colors depending on if the plate is continental or not
    plate.plateMaterial = bpy.data.materials.new(name="Material")
    if plate.isContinental:
        plate.plateMaterial.diffuse_color = (0.6-0.6*random(), 1-0.3*random(), 0.1*random(), 1.0)
        plate.plateMaterial.metallic = 0.4 * random()
        plate.plateMaterial.roughness = 0.7 * random()
    else:
        plate.plateMaterial.diffuse_color = (0.1*random(), 0.3*random(), 1, 1.0)
        plate.plateMaterial.metallic = 0.7 * random()
    bpy.context.active_object.data.materials.append(plate.plateMaterial)


#Use spherical polar coordinates to apply simply noise function for initial terrain
def applyNoiseFunctions(bm, props):
    for v in bm.verts:
        r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
        
        #Some noise functions already built into blender
        noise = M.noise.noise(v.co * props.lacunarity)
        hetero = M.noise.hetero_terrain(v.co, 0.3, 2*props.lacunarity, 2, 0.2, noise_basis='VORONOI_F1')
        crackle = M.noise.hetero_terrain(v.co, 0.6, props.lacunarity, 2, 0.1, noise_basis='VORONOI_CRACKLE')
        
        #Calculate and apply dx, dy, and dz
        dr = props.noiseHeight * (noise + hetero)
        v.co.x += sin(theta) * cos(phi) * dr
        v.co.y += sin(theta) * sin(phi) * dr
        v.co.z += cos(theta) * dr

#For each vertex, we find the closest (centroid +- noise) and calculate the vertex's plate ID accordingly
#Note: This is one of the codes bottlenecks right now.
def calculatePlateId(bm, centroids, props, plateId):
    initialRadius = props.radius
    noiseAmp = props.boundNoiseAmp
    for v in bm.verts:
        
        #We need to recalculate plateId each time the sphere is being remeshed
        #We don't want the radial component of v.co to be included in the noise sample for plate boundaries
        #So we create a vector "noiseCoord" that has the same theta and phi components as v.co,
        #But has a constant radial component that is set to the initialRadius
        closestCentroid = centroids[0]
        rThetaPhi = cartToPolarTransformation(v.co)
        rThetaPhi.x = initialRadius
        noiseCoord = polarToCartTransformation(rThetaPhi)
        
        #Set the vertex's plateId to the index of the closestCentroid
        noise = M.noise.noise(noiseCoord * props.boundLac)
        for c in centroids:
            if (v.co - c).length < ((v.co - closestCentroid).length + noiseAmp * noise):
                closestCentroid = c
        v[plateId] = centroids.index(closestCentroid)
        
#Controids need to be moved along with the rest of the plate
def moveCentroids(quaternions, centroids, plateProps):
    centroids = [quaternions[i] @ cent for i, cent in enumerate(centroids)]
    for i, r in enumerate(plateProps.ids):
        r.plateCentroid[0] = centroids[i].x
        r.plateCentroid[1] = centroids[i].y
        r.plateCentroid[2] = centroids[i].z




#================================Subduction =======================================================================

#Main function for applying subduction uplift
def applySubductionUplift(bm, props, subductionFronts, subFronts, plateId, plateProps, quaternions, maxSpeed):
    
    #Properties to be used throughout subduction
    isContinental = [r.isContinental for r in plateProps.ids]
    minHeight = np.sqrt(min([(v.co.x**2 + v.co.y**2 + v.co.z**2) for v in bm.verts]))
    maxHeight = np.sqrt(max([(v.co.x**2 + v.co.y**2 + v.co.z**2) for v in bm.verts]))
    heightLimit = props.maxHeight
    aveHeights = calculateAveragePlateHeight(bm, props, plateProps, plateId)
    
    subFrontsCo = [e.verts[0].co for e in subFronts]
    kdSubFronts = createKDTree(subFrontsCo)
    maxRange = props.maxUpliftRange
    minRange = props.minUpliftRange
    baseUplift = props.baseUplift
    
    #If the list of edges "subductionFronts is not empty, iterate through all vertices in bm
    #The parameters dist, speed, and normElev correspond to d(p), v(p), and z(p) from subduction uplift equation
    if subductionFronts:
        for v in bm.verts:
            closestEdgeCo, index, dist = kdSubFronts.find(v.co)
            closestEdge = subFronts[index]
            aveHeight = aveHeights[v[plateId]]
            upRange = minRange + (maxRange - minRange) * aveHeight / heightLimit
            
            #Only bother calculating uplift if v is within the subduction region
            if (dist <= upRange) and (isContinental[v[plateId]]):
                speed = getRelativeSpeedForSubduction(v, closestEdge, plateId, quaternions)
                normElev = normalizedElevation(v, minHeight, maxHeight)
                
                #More parameters used in the subduction uplift equation
                distanceTransfer = getDistanceTransfer(dist / upRange)
                speedTransfer = speed / maxSpeed
                heightTransfer = normElev**2
                uplift = float(distanceTransfer * speedTransfer * heightTransfer)
                
                #We apply uplift to radial component of vertex v
                r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
                dr = baseUplift * uplift * (heightLimit - aveHeight)
                v.co.x += sin(theta) * cos(phi) * dr
                v.co.y += sin(theta) * sin(phi) * dr
                v.co.z += cos(theta) * dr

#Sample a point from the distance transfer curve editor      
def getDistanceTransfer(x):
    distTransCurve = bpy.data.node_groups['ProfileCurves'].nodes['Distance Transfer'].mapping
    curve = distTransCurve.curves[3]
    distTransCurve.initialize()
    return distTransCurve.evaluate(curve, x)

#Slab pull has the effect of changing the direction that a plate is moving towards
def applySlabPull(subductionFronts, plateProps, isCont, centroids):
    
    #Lists to be used in main loop
    axes = [M.Vector((r.axes[0], r.axes[1], r.axes[2])) for r in plateProps.ids]
    frontsCo = [[j.verts[0].co for j in i] for i in subductionFronts]
    
    #We apply the slab pull equation from the paper to each plate i
    #The parameters axes, centroids and fronts correspond to vectors w, c and q in the paper
    for i in range(len(frontsCo)):
        if not isCont[i]:
            sum = M.Vector((0.0, 0.0, 0.0))
            for front in frontsCo[i]:
                numerator = axes[i].cross(front)
                denominator = numerator.length
                sum += numerator / denominator
            
            #We apply the new axes to plates i
            newAxes = (axes[i] + 0.1 * sum)
            plateProps.ids[i].axes = newAxes

#subFronts is a nested list of edges which are subduction regions
#The index of each sublist corresponds to plateIds of plates in that particular edge
#While a nested list seems more complicated than necesary, this significantly helps reduce our code runtime later on
def getSubFronts(boundingEdges, plateQuaternions, plateId, plateNums, isContinent):
    subductionFronts = list(filter(lambda x: isSubduction(x, plateQuaternions, plateId, isContinent), boundingEdges))
    subFronts = [[edge for edge in subductionFronts if (edge.verts[0][plateId] == i) or (edge.verts[1][plateId] == i)] for i in range(plateNums)]
    return subFronts

#If the two vertices of a bounding edge decrease in distance, we have found a subduction front
def isSubduction(edge, quaternions, plateId, isContinent):
    quat0 = quaternions[edge.verts[0][plateId]]
    quat1 = quaternions[edge.verts[1][plateId]]
    coord0 = edge.verts[0].co
    coord1 = edge.verts[1].co
    isGettingCloser = (coord0 - coord1).length > ((quat0 @ coord0) - (quat1 @ coord1)).length
    oneEndIsOceanic = not (isContinent[edge.verts[0][plateId]] and isContinent[edge.verts[1][plateId]])
    return isGettingCloser and oneEndIsOceanic

#We loop through the list of edges "subductionFronts", and find the closest edge to specified vertex
#Note: At the moment, this function is the slowest bottleneck in our code
#As such, the function may look messy but is done so for efficiency purposes
def getClosestSubductionEdge(vertex, subductionFronts, plateId):
    id = vertex[plateId]
    if subductionFronts[id]:
        closestEdge = subductionFronts[id][0]
        closestEdgeCo = closestEdge.verts[0].co
        vertCo = vertex.co
        currentDist = vertCo - closestEdgeCo
        
        #Loop through all edges that may be on the subduction front, and find the closest edge to given vertex
        for edge in subductionFronts[id]:
            edgeCo = edge.verts[0].co
            if (vertCo - edgeCo) < currentDist:
                closestEdge = edge
                currentDist = vertCo - closestEdge.verts[0].co
        return closestEdge

def getDistanceToSubductionFront(vertex, closestEdge):
    otherVertex = closestEdge.verts[0]
    return (vertex.co - otherVertex.co).length

#Finds the relative speed that a vertex is to closest vertex on other plate of at subduction front
def getRelativeSpeedForSubduction(vertex, closestEdge, plateId, quaternions):
    
    #Each edge only has two vertices attached to them, here we find which of those vertices is on the other plate
    otherVertex = closestEdge.verts[0]
    if vertex[plateId] == otherVertex[plateId]:
        otherVertex = closestEdge.verts[1]
    
    #We calculate and return their relative speed
    quat0 = quaternions[vertex[plateId]]
    quat1 = quaternions[otherVertex[plateId]]
    velocity0 = (quat0 @ vertex.co) - vertex.co
    velocity1 = (quat1 @ otherVertex.co) - otherVertex.co
    return (velocity1 - velocity0).length

#Normalized elevation as required for subduction uplift
def normalizedElevation(v, min, max):
    vertexHeight = np.sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2)
    return (vertexHeight - min)/(max - min)

#The average plate height is used for dynamically redifining the radius of uplift influence
def calculateAveragePlateHeight(bm, props, plateProps, plateId):
    heights = [[] for i in range(len(plateProps.ids))]
    initialRadius = props.radius
    
    for v in bm.verts:
        height = np.sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2) - initialRadius
        heights[v[plateId]].append(height)
    
    aveHeights = []
    for i, p in enumerate(plateProps.ids):
        if len(heights[i]) > 0:
            p.averagePlateHeight = sum(heights[i]) / len(heights[i])
        else:
            p.averagePlateHeight = 0
        aveHeights.append(p.averagePlateHeight)
    return aveHeights

def erodeContinents(bm, props, isContinent, plateId):
    initialHeight = props.radius
    maxHeight = max([np.sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2) for v in bm.verts]) - initialHeight
    coefficient = props.contErosion
    for v in bm.verts:
        if not isContinent[v[plateId]]:
            continue
        
        r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
        height = r - initialHeight
        dr = - coefficient * height / maxHeight
        v.co.x += math.sin(theta) * math.cos(phi) * dr
        v.co.y += math.sin(theta) * math.sin(phi) * dr
        v.co.z += math.cos(theta) * dr

#=================================Continental Collision =============================================

#Main function for continental collision
def applyContinentalCollision(bm, props, plateId, isTerrane, quaternions, centroids, maxSpeed):
    
    #Bunch of properties and lists that may be needed later
    contColEdges = [e for e in bm.edges if lambda x: isContColEdge(x, plateId, quaternions)]
    contColEdgeCo = [e.verts[0].co for e in contColEdges]
    terranes = [v.co for v in bm.verts if (v[isTerrane] == 1)]
    kdColFronts = createKDTree(contColEdgeCo)
    kdTerranes = createKDTree(terranes)
    
    aveArea = getAveragePlateArea(bm, props)
    plateAreas = getPlateAreas(bm, props, plateId)
    discreteCollisionCoefficient = 0.001 * props.discreteCollisionCoef
    rMax = props.contColMax
    
    #If the list contColEdges and terranes is not empty, loop through all vertices
    if contColEdges and terranes:
        for v in bm.verts:
            
            #Calculate the radius of influence
            closestEdgeCo, index, dist = kdColFronts.find(v.co)
            closestTerraneCo, indexTerrane, terraneDistance = kdTerranes.find(v.co)
            closestEdge = contColEdges[index]
            speed = getRelativeSpeedForSubduction(v, closestEdge, plateId, quaternions)
            radiusOfInfluence = np.sqrt(speed / maxSpeed) * (plateAreas[v[plateId]] / aveArea) * rMax
            
            #If the vertex is within the radius of influence, calculate and apply elevation surge
            if terraneDistance < radiusOfInfluence:
                distanceContribution = contColDistanceTransfer(terraneDistance, radiusOfInfluence)
                area = plateAreas[v[plateId]]
                elevationSurge = discreteCollisionCoefficient * area * distanceContribution * (1 - (terraneDistance / radiusOfInfluence))
                
                #We apply elevation surve to vertex v
                r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
                dr = elevationSurge
                v.co.x += sin(theta) * cos(phi) * dr
                v.co.y += sin(theta) * sin(phi) * dr
                v.co.z += cos(theta) * dr
    

#Initiate the "isTerrane" custom layer which is used for tracking which vertices corresponds to terrane
def initializeTerranes(bm, props, plateId, isTerrane, plateProps):
    isContinent = [r.isContinental for r in plateProps.ids]
    oceanContEdgeCo = [e.verts[0].co for e in bm.edges if isOceanicContEdge(e, plateId, isContinent)]
    kdOceanCont = createKDTree(oceanContEdgeCo)
    noiseAmp = props.boundNoiseAmp
    noiseLacunarity = props.boundLac
    terraneHeight = props.radius + props.oceanHeight
    
    #Label each vertex's "isTerrane" parameter as a 1 or 0 (True or False) based on height
    for v in bm.verts:
        r = np.sqrt(v.co.x**2 + v.co.y**2 + v.co.z**2)
        if r > terraneHeight:
            v[isTerrane] = 1
        else:
            v[isTerrane] = 0

#Distance transfer function as from paper
def contColDistanceTransfer(x, r):
    return (1 - (x / (r + 1))**2 )**2
    
#Calculate the average plate area
def getAveragePlateArea(bm, props):
    totalArea = sum([f.calc_area() for f in bm.faces])
    n = props.plateNums
    return totalArea / n

#Calculate the plate surface area
def getPlateAreas(bm, props, plateId):
    areas = [0 for i in range(props.plateNums)]
    bm.verts.ensure_lookup_table()
    for i in range(len(bm.faces)):
        face = bm.faces[i]
        vert = face.verts[0]
        areas[vert[plateId]] += face.calc_area()
    return areas

#Given an edge, we check if it is a continental collision boundary
def isContColEdge(edge, plateId, quaternions):
    vOne = e.verts[0]
    vTwo = e.verts[1]
    quatOne = quaternions[vOne[plateId]]
    quatTwo = quaternions[vTwo[plateId]]
    isOnDiffPlate = (vOne[plateId] != vTwo[plateId])
    isColliding = (quatOne @ vOne - quatTwo @ vTwo).length < (vOne - vTwo).length
    return (isOnDiffPlate and isColliding)

#Used to check if a given edge is at a oceanic-continental collision front
def isOceanicContEdge(e, plateId, isContinent):
    vOne = e.verts[0]
    vTwo = e.verts[1]
    isEdge = (isContinent[vOne[plateId]] != isContinent[vTwo[plateId]])
    return isEdge



#======================================Oceanic Crust Generation ============================================

#We use the template ridge profile from the curve editor, to generate oceanic ridges
def applyOceanicCrustGeneration(bm, props, plateId, isTerrane, plateProps):
    quaternions = [M.Quaternion((r.axes[0], r.axes[1], r.axes[2]), radians(r.angle)) for r in plateProps.ids]
    boundingEdges = getPlateBoundaries(bm)
    aveRadius = props.aveCrustRad
    lacunarity = props.crustNoiseLac
    amplitude = props.crustNoiseAmp
    depth = props.ridgeDepth
    
    #Create list of ridge edges, coordinates, a ridge KDTree and another KDTree for points nearby ridge
    ridges = list(filter(lambda x: isRidge(x, plateId, quaternions, plateProps), boundingEdges))
    ridgeCo = [e.verts[0].co for e in ridges]
    ridgeKD = createKDTree(ridgeCo)
    ridgeVicinityKD = createKDTree([v.co for v in bm.verts if ridgeKD.find(v.co)[2] <= (aveRadius + amplitude)])
    
    for v in bm.verts:
        
        #We find how close v is from an oceanic ridge
        closestRidgeCo, index, distToRidge = ridgeKD.find(v.co)
        closestEdge = ridges[index]
        edgeLength = (closestEdge.verts[0].co - closestEdge.verts[1].co).length
        
        #The ridgeRadius used to define how far away points are affected by the oceanic crust generation algoritm
        noise = M.noise.noise(v.co * lacunarity)
        ridgeRadius = aveRadius + amplitude * edgeLength * (2 * noise - 1)
        
        #Break the current loop iteration if v.co is too far away from a ridge
        if distToRidge > ridgeRadius:
            continue
            
        #The interpolation factor corresponds to alpha in the tectonic planets paper
        distToPlate = ridgeRadius - distToRidge
        interpolationFactor = distToRidge / (distToRidge + distToPlate)
        
        #We find the two points that correspond to the ridge boundary
        directionToRidge = (closestRidgeCo - v.co).normalized()
        pointOne, iOne, distOne = ridgeVicinityKD.find(closestRidgeCo + ridgeRadius * directionToRidge)
        pointTwo, iTwo, distTwo = ridgeVicinityKD.find(closestRidgeCo - ridgeRadius * directionToRidge)
        ridgeRad = (pointOne - pointTwo).length / 2
        
        #Break the current loop to avoid division by zero
        if ridgeRad == 0.0:
            continue
        
        #We need a way of differentiating the left and right side of the ridge
        #So we let the left hand side be defined by the side of the ridge with a lower plateId
        isLeft = (v[plateId] < closestEdge.verts[0][plateId]) or (v[plateId] < closestEdge.verts[1][plateId])
        if isLeft:
            samplePoint = 0.5 - 0.5 * distToRidge / ridgeRad
        else:
            samplePoint = 0.5 + 0.5 * distToRidge / ridgeRad
        
        #We sample points from the user defined ridge template, 
        #And create templatePoint which represents the templates contribution in changing r 
        templateSampledPoint = sampleRidgeProfile(samplePoint)
        
        #We add the templatedSamplePoint to the radial component of v.co,
        #But preserve the old theta and phi components
        r, theta, phi = getPolarCoords(v.co.x, v.co.y, v.co.z)
        rNew = r + (templateSampledPoint + amplitude * (2 * noise - 1)) * depth
        x, y, z = getCartesianCoords(rNew, theta, phi)
        templatedPoint = M.Vector((x, y, z))
        
        #Used to calculate the linear sample point
        directionFromOneToTwo = (pointOne - pointTwo).normalized()
        distFromOneToTwo = (pointOne - pointTwo).length
        distToOne = (v.co - pointOne).length
        
        if distFromOneToTwo == 0:
            continue
        
        #The linear sample point corresponds to the linearly interpolated elevation between plates
        #Again, we disregard theta and phi changes made
        linearSamplePoint = pointOne + directionFromOneToTwo * distToOne / distFromOneToTwo
        linR, linTheta, linPhi = getPolarCoords(linearSamplePoint.x, linearSamplePoint.y, linearSamplePoint.z)
        linX, linY, linZ = getCartesianCoords(linR, theta, phi)
        linPoint = M.Vector((linX, linY, linZ))
        
        #Apply oceanic crust generation to points
        v.co = interpolationFactor * linPoint + (1 - interpolationFactor) * templatedPoint
        
#Used to find oceanic ridges
def isRidge(edge, plateId, quaternions, plateProps):
    vert1 = edge.verts[0]
    vert2 = edge.verts[1]
    quat1 = quaternions[vert1[plateId]]
    quat2 = quaternions[vert2[plateId]]
    isContinent = [r.isContinental for r in plateProps.ids]
    
    #Check if both sides of the edge are oceanic, and both sides are diverging
    isOceanOcean = (not isContinent[vert1[plateId]]) and (not isContinent[vert2[plateId]])
    isDiverging = ((quat1 @ vert1.co) - (quat2 @ vert2.co)).length > (vert1.co - vert2.co).length
    return isOceanOcean and isDiverging
    
        
#Function that maps x (from 0 to 1) to the template ridge profile as defined in curve editor
def sampleRidgeProfile(x):
    ridgeProfile = bpy.data.node_groups['ProfileCurves'].nodes[0].mapping
    curve = ridgeProfile.curves[3]
    ridgeProfile.initialize()
    return (ridgeProfile.evaluate(curve, x) - 0.5)




#==========================Random Other Functions=====================================================

#Since we are using a sphere, polar coordinates are used whenever we which to apply a change in height
#Coordinate transformation function from cartesian to polar
def getPolarCoords(X, Y, Z):
    R = np.sqrt(X**2 + Y**2 + Z**2)
    Theta = np.arccos(Z / R)
    Phi = np.arctan2(Y, X)    
    return (R, Theta, Phi)

#Coordinate transformation function from polar to cartesian
def getCartesianCoords(R, Theta, Phi):
    X = R * np.sin(Theta) * np.cos(Phi)
    Y = R * np.sin(Theta) * np.sin(Phi)
    Z = R * np.cos(Theta)
    return (X, Y, Z)

#Creates a list of all edges that correspond to a plate boundary
def getPlateBoundaries(bm):
    plateId = bm.verts.layers.int.get('id')
    boundingEdges = []
    
    #If the two vertices at the end of the edge are not on the same plate, we have found a boundary edge.
    for e in bm.edges:
        if e.verts[0][plateId] != e.verts[1][plateId]:
            boundingEdges.append(e)
    return boundingEdges

#KDTrees is a blender tool for really fast spacial searches
def createKDTree(coords):
    kd = M.kdtree.KDTree(len(coords))
    for i, c in enumerate(coords):
        kd.insert(c, i)
    kd.balance()
    return kd


#=====================Panels =================================================
#Panels are used for defining the user interface

class TectonicToolsMainPanel:
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Tools"
    #bl_options = {"DEFAULT_CLOSED"}


#Main panel for tectonic tools
class OBJECT_PT_GeneratorPanel(TectonicToolsMainPanel, bpy.types.Panel):
    bl_idname = "OBJECT_PT_generatorPanel"
    bl_label = "Terrain Generator"

    #This is where we add components to the user interface
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        
#Sub-panel for specifying initial world properties
class OBJECT_PT_InitialProperties(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Global Properties"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        layout.label(text="Initial Sphere")
        layout.prop(props, "subDivs")
        layout.prop(props, "radius")
        layout.prop(props, "lacunarity")
        layout.prop(props, "noiseHeight")
        layout.prop(props, "plateSpeed")
        layout.prop(props, "oceanHeight")
        
        layout.label(text="Voronoi Plate Generator")
        layout.prop(props, "nOfPlates")
        layout.prop(props, "boundLac")
        layout.prop(props, "boundNoiseAmp")
        
        layout.label(text="Subduction Parameters")
        layout.prop(props, "baseUplift")
        layout.prop(props, "maxHeight")
        layout.prop(props, "maxUpliftRange")
        layout.prop(props, "minUpliftRange")
        
        layout.label(text="Coninental Collision")
        layout.prop(props, "discreteCollisionCoef")
        layout.prop(props, "contColMax")
        
        layout.label(text="Oceanic Crust Generator")
        layout.prop(props, "aveCrustRad")
        layout.prop(props, "ridgeDepth")
        layout.prop(props, "crustNoiseAmp")
        layout.prop(props, "crustNoiseLac")
        
        layout.label(text="Other Parameters")
        layout.prop(props, "contErosion")
        
#Sub-panel for initilization operators
class OBJECT_PT_InitializationTools(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Buttons"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        #For initializing our terrain
        layout.operator("wm.initiate_terrain")
        layout.operator("wm.voronoi")
        layout.operator("wm.applyvertexcolors")
        layout.operator("wm.move_selection")
        layout.operator("wm.remesh")
        layout.operator("wm.split_plates")
        layout.operator("wm.applycolors")
        layout.operator("wm.testing")

#Curve editors let the user define their own transfer functions
class OBJECT_PT_CurveEditors(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Curve Editors"
    def draw(self, context):
        layout = self.layout
        curveTree = curveNodeTree()
        
        layout.label(text="Template Ridge Profile")
        if "Ridge Template" in curveTree.keys():
            layout.template_curve_mapping(curveTree["Ridge Template"], "mapping")
        
        layout.label(text="Subduction Uplift Distance Transfer")
        if "Distance Transfer" in curveTree.keys():
            layout.template_curve_mapping(curveTree["Distance Transfer"], "mapping")



#Used for automating the simulation and saving images during each iteration
class OBJECT_PT_Animate(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Animate"
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties

        layout.prop(props, "totalFrames")
        layout.prop(props, "remeshAfter")
        layout.prop(props, "outputFolder")
        layout.operator("wm.animate")
        


                    

#Sub-panels for defining plate movement
class OBJECT_PT_PlateProps(TectonicToolsMainPanel, Panel):
    bl_parent_id = "OBJECT_PT_generatorPanel"
    bl_label = "Plate Properties"
    
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        props = scene.properties
        
        if bpy.context.objects_in_mode:
            plateProps = context.object.data.WMTPlatePointers
            plateMaterials = [r.plateMaterial for r in plateProps.ids]
            for i in range(props.plateNums + 1):
                layout.label(text="Plate {}".format(i))
                layout.prop(plateProps.ids[i], "axes")
                layout.prop(plateProps.ids[i], "angle")
                layout.prop(plateProps.ids[i], "isContinental")
                layout.prop(plateMaterials[i], "diffuse_color")


#We create a node tree which can store curve data, which is used in the curve editor
def curveNodeTree():
    if 'ProfileCurves' not in bpy.data.node_groups:
        bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
    nodeTree = bpy.data.node_groups['ProfileCurves'].nodes
    return nodeTree
        
#We create a node group to store our profile curves in
#A profile curve allows the user to customize functions in the user interface
#Note that the curve editor that I'm using was originally intended for editing shaders and RGB colors,
#So some of the built in terminology refers to "shading" and "RGB" even though that's not what I'm using the curve editor for
def initializeProfileCurves():
    if 'ProfileCurves' not in bpy.data.node_groups:
        bpy.data.node_groups.new('ProfileCurves', 'ShaderNodeTree')
    
    #Create a node curve that defines the template ridge profile
    nods = bpy.data.node_groups['ProfileCurves'].nodes
    if "Ridge Template" not in nods.keys():
        ridgeCurves = nods.new('ShaderNodeRGBCurve')
        ridgeCurves.name = "Ridge Template"
        pnts = ridgeCurves.mapping.curves[3].points
        pnts[0].location = [0.0, 0.5]
        pnts[1].location = [1.0, 0.5]
        pnts.new(0.5, 0.3)
        pnts.new(0.25, 0.2)
        pnts.new(0.75, 0.2)
    
    if "Distance Transfer" not in nods.keys():
        distCurves = nods.new('ShaderNodeRGBCurve')
        distCurves.name = "Distance Transfer"
        pnts = distCurves.mapping.curves[3].points
        pnts[0].location = [0.0, 0.8]
        pnts[1].location = [1.0, 0.0]
        pnts.new(0.2, 1.0)
        pnts.new(0.75, 0.15)
        

#====================Register Classes To Blender =============================

classes = (
    GeneratorProperties,
    WMTPlateProperties,
    WMTPlatePointers, 
    WM_OT_Initiate, 
    WM_OT_SplitPlates,
    WM_OT_MoveSelectedPlate,
    WM_OT_CreatePlates,
    WM_OT_ApplyColor,
    WM_OT_ApplyVertexColor,
    WM_OT_Remesh,
    WM_OT_Animate,
    WM_OT_Testing,
    OBJECT_PT_GeneratorPanel,
    OBJECT_PT_InitialProperties,
    OBJECT_PT_InitializationTools,
    OBJECT_PT_CurveEditors,
    OBJECT_PT_Animate,
    OBJECT_PT_PlateProps,
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
    del bpy.types.Mesh.WMTRotationPointers


if __name__ == "__main__":
    register()

