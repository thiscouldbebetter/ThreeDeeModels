import bpy
import math
import mathutils

# Classes.
class Color:
    def __init__(self, name, componentsRGB):
        self.name = name
        self.componentsRGB = componentsRGB

    def Instances():
        return Color_Instances(); #hack - Figure out singletons in Python.
            
    def componentsAsTuple(self):
        components = self.componentsRGB
        return (components[0], components[1], components[2], 1)


class Color_Instances:
    def __init__(self):
        self.Black = Color("Black", [0, 0, 0])    
        self.Blue = Color("Blue", [0, 0, 1])
        self.Brown = Color("Brown", [.5, .25, 0])    
        self.Cyan = Color("Cyan", [0, 1, 1])
        self.GrayDark = Color("GrayDark", [.33, .33, .33])
        self.GrayLight = Color("GrayLight", [.67, .67, .67])
        self.Green = Color("Green", [0, 1, 0])
        self.Orange = Color("Orange", [1, .5, 0])
        self.Red = Color("Red", [1, 0, 0])
        self.Violet = Color("Violet", [1, 0, 1])
        self.White = Color("White", [1, 1, 1])    
        self.Yellow = Color("Yellow", [1, 1, 0])

        
class Coords:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def ones():
        return Coords(1, 1, 1)

    def zeroes():
        return Coords(0, 0, 0)

    def add(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def addDimensions(self, x, y, z):
        self.x += x
        self.y += y
        self.z += z
        return self
    
    def crossProduct(self, other):
        return self.overwriteWithDimensions(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        );
    
    def divideScalar(self, scalar):
        if scalar != 0:
            self.x /= scalar
            self.y /= scalar
            self.z /= scalar
            
        return self

    def magnitude(self):
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def multiply(self, other):
        self.x *= other.x
        self.y *= other.y
        self.z *= other.z
        return self
        
    def multiplyScalar(self, scalar):
        self.x *= scalar
        self.y *= scalar
        self.z *= scalar
        return self

    def normalize(self):
        self.divideScalar(self.magnitude())
        
    def overwriteWithDimensions(self, x, y, z):
        self.x = x
        self.y = y;
        self.z = z;
        return self
        
    def toStringXYZ(self):
        return str(self.x) + "x" + str(self.y) + "x" + str(self.z);

    def toSystemVertex(self):
        return (self.x, 0 - self.y, self.z)

    # Clonable.
    
    def clone(self):
        return Coords(self.x, self.y, self.z)
    
    def overwriteWith(self, other):
        self.x = other.x
        self.y = other.y
        self.z = other.z
        return self
        

class Disposition:
    def __init__(self, pos, ori):
        self.pos = pos
        self.ori = ori
       
    def fromPos(pos):
        return Disposition(
            pos,
            Orientation.default()
        )


class Material:
    def __init__(self, name, color):
        self.name = name
        self.color = color
        
        self.systemMaterial = None
    
    def renderToBpy(self, bpy):
        if self.systemMaterial == None:
            systemMaterial = bpy.data.materials.new(
                name = self.name
            )
            systemColor = self.color.componentsAsTuple()
            systemMaterial.diffuse_color = systemColor
            self.systemMaterial = systemMaterial
        return self.systemMaterial


class Mesh:
    def __init__(self, name, vertices, faces):
        self.name = name
        self.vertices = vertices
        self.faces = faces
        
        self.systemMesh = None

    # Bpy.
            
    def renderToBpy(self, bpy):

        if self.systemMesh == None:
            vertices = self.vertices
            systemVertices = []
            for vertex in vertices:
                systemVertex = vertex.toSystemVertex()
                systemVertices.append(systemVertex)
                
            faces = self.faces
            systemFaces = []
            for face in faces:
                systemFace = face.vertexIndicesAsTuple()
                systemFaces.append(systemFace)
             
            name = self.name
            systemMesh = bpy.data.meshes.new(name)

            systemMesh.from_pydata(systemVertices, [], systemFaces)
            systemMesh.update(calc_edges=True)

            self.systemMesh = systemMesh
                    
        return self.systemMesh
    

    def renderToBpyAndObject(self, bpy, systemObject):
        return self;
    
    # Clonable.
    
    def clone(self):
        return Mesh(
            self.name + "_Clone",
            list(map(lambda x : x.clone(), self.vertices)),
            list(map(lambda x :  x.clone(), self.faces))
        )
        
    # Transformable.
    
    def transform(self, transformToApply):
        for vertex in self.vertices:
            transformToApply.transformCoords(vertex)
            
        return self


class MeshBuilder:
    
    def box(self, size):
        
        x = size.x
        y = size.y
        z = size.z
        
        mesh = Mesh(
            "Box" + size.toStringXYZ(),
            #vertices
            [
                Coords(0, 0, 0),
                Coords(x, 0, 0),
                Coords(x, y, 0),
                Coords(0, y, 0),
                
                Coords(0, 0, z),
                Coords(x, 0, z),
                Coords(x, y, z),
                Coords(0, y, z),
            ],
            [
                #bottom
                MeshFace([0, 1, 2, 3]),
                
                #sides
                MeshFace([0, 1, 5, 4]),
                MeshFace([1, 2, 6, 5]),
                MeshFace([2, 3, 7, 6]),
                MeshFace([3, 0, 4, 7]),
                                        
                #top
                MeshFace([4, 5, 6, 7]),
            ]
        )
        
        return mesh
    
    def brickWithTopStuds(self, sizeInStuds):
        
        meshBody = self.box(sizeInStuds)
        
        meshesToMerge = [ meshBody ]
        
        scaleFactor = .25
        transformScale = Transform_Scale(
            Coords.ones().multiplyScalar(scaleFactor)
        )
        meshStudPrototype = self.tube(32).transform(transformScale)
        
        displacement = Coords.zeroes()
        transformDisplace = Transform_Displace(displacement)
        
        for y in range(sizeInStuds.y):
            
            for x in range(sizeInStuds.x):
                
                displacement.overwriteWithDimensions(
                    x, y, 0
                ).addDimensions(
                    .5, .5, 1 + scaleFactor
                )
                meshStudDisplaced = meshStudPrototype.clone().transform(
                    transformDisplace
                )
                
                meshesToMerge.append(meshStudDisplaced)
        
        meshName = "Brick_" + sizeInStuds.toStringXYZ()
        meshMerged = self.mergeMeshes(meshName, meshesToMerge)
        
        return meshMerged
    
    def extrudeFromVertices(self, meshName, verticesBottom):
    
        verticesBottomCount = len(verticesBottom)
        
        vertices = list(
            map(
                lambda x : x.clone(),
                verticesBottom
            )
        )
        
        displacementBottomToTop = Coords(0, 0, -1)
        for i in range(verticesBottomCount):
            vertexBottom = vertices[i]
            vertexTop = vertexBottom.clone().add(displacementBottomToTop)
            vertices.append(vertexTop)
            
        faces = []
        
        vertexIndicesBottom = []
        vertexIndicesTop = []
        for i in range(0, verticesBottomCount):
            vertexIndexBottom = i
            vertexIndexTop = i + verticesBottomCount
            vertexIndicesBottom.append(vertexIndexBottom)
            vertexIndicesTop.append(vertexIndexTop)
            
        faceBottom = MeshFace(vertexIndicesBottom)
        faceTop = MeshFace(vertexIndicesTop)
        faces.append(faceBottom)
        faces.append(faceTop)
        
        for vertexBottomIndex in range(verticesBottomCount):

            vertexTopIndex = vertexBottomIndex + verticesBottomCount

            vertexBottomNextIndex = vertexBottomIndex + 1
            if (vertexBottomNextIndex >= verticesBottomCount):
                vertexBottomNextIndex = 0

            vertexTopNextIndex = vertexTopIndex + 1
            if (vertexTopNextIndex >= verticesBottomCount * 2):
                vertexTopNextIndex = verticesBottomCount
            
            faceSide = MeshFace(
                [
                    vertexBottomIndex, vertexTopIndex,
                    vertexTopNextIndex, vertexBottomNextIndex
                ]
            )
            
            faces.append(faceSide)
                
        meshExtruded = Mesh(
            meshName,
            vertices,
            faces
        )
        
        return meshExtruded

    def gear(self, toothCount, radiusInnerOverOuter):
    
        radiansPerTurn = 2 * math.pi # "tau"
        radiansPerPoint = radiansPerTurn / toothCount
        radiansPerPointHalf = radiansPerPoint / 2
        
        vertices = []
        for i in range(0, toothCount):

            angleInRadians = i * radiansPerPoint

            vertexNearInner = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            ).multiplyScalar(
                radiusInnerOverOuter
            )
            vertices.append(vertexNearInner)

            vertexNearOuter = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            )
            vertices.append(vertexNearOuter)
            
            angleInRadians += radiansPerPointHalf

            vertexFarOuter = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            )
            vertices.append(vertexFarOuter)

            vertexFarInner = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            ).multiplyScalar(
                radiusInnerOverOuter
            )
            vertices.append(vertexFarInner)
    
        meshName = "GearWith" + str(toothCount) + "Teeth"
        mesh = self.extrudeFromVertices(meshName, vertices)
                
        return mesh

    def mergeMeshes(name, self, meshesToMerge):

        nameMerged = "Merge_"
        verticesMerged = []
        facesMerged = []

        numberOfVerticesSoFar = 0;
        
        for meshToMerge in meshesToMerge:

            nameMerged += meshToMerge.name
                        
            verticesToMerge = meshToMerge.vertices
            verticesMerged.extend(verticesToMerge)
            
            facesToMerge = meshToMerge.faces
            
            for face in facesToMerge:
                face.vertexIndicesShift(numberOfVerticesSoFar)
                facesMerged.append(face)

            numberOfVerticesSoFar += len(verticesToMerge)

        returnMesh = Mesh(
            nameMerged, verticesMerged, facesMerged
        );

        return returnMesh
    
    def tube(self, sideCount):
    
        radiansPerTurn = 2 * math.pi # "tau"
        
        vertices = []
        for i in range(0, sideCount):
            angleInRadians = i * radiansPerTurn / sideCount
            x = math.sin(angleInRadians)
            y = math.cos(angleInRadians)
            vertex = Coords(x, y, 0)
            vertices.append(vertex)
    
        meshName = "TubeWith" + str(sideCount) + "Sides"
        mesh = self.extrudeFromVertices(meshName, vertices)
                
        return mesh
    
    def star(self, pointCount, radiusInnerOverOuter):
    
        radiansPerTurn = 2 * math.pi # "tau"
        radiansPerPoint = radiansPerTurn / pointCount
        radiansPerPointHalf = radiansPerPoint / 2
        
        vertices = []
        for i in range(0, pointCount):

            angleInRadians = i * radiansPerPoint
            vertexInner = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            ).multiplyScalar(
                radiusInnerOverOuter
            )
            vertices.append(vertexInner)
            
            angleInRadians += radiansPerPointHalf
            vertexOuter = Coords(
                math.cos(angleInRadians),
                math.sin(angleInRadians),
                0
            )
            vertices.append(vertexOuter)
    
        meshName = "StarWith" + str(pointCount) + "Points"
        mesh = self.extrudeFromVertices(meshName, vertices)
                
        return mesh
                           

class MeshDisposition:
    def __init__(self, mesh, disp):
        self.mesh = mesh
        self.disp = disp

    def renderToBpy(self, bpy):

        mesh = self.mesh

        systemMesh = mesh.renderToBpy(bpy)
        systemObject = bpy.data.objects.new(
            mesh.name, systemMesh
        )

        mesh.renderToBpyAndObject(
            bpy, systemObject
        )

        disp = self.disp
        
        pos = disp.pos
        systemObject.location = pos.toSystemVertex()
        
        ori = disp.ori
        oriAsEulerRotationTuple = ori.toEulerRotationTuple()
        systemObject.rotation_euler = mathutils.Euler(
            oriAsEulerRotationTuple, "XYZ"
        )

        collection = bpy.context.scene.collection
        collection.objects.link(systemObject)
         
                 
class MeshFace:
    def __init__(self, vertexIndices):
        self.vertexIndices = vertexIndices
        
    def vertexIndicesAsTuple(self):
        return tuple(self.vertexIndices)
    
    def vertexIndicesShift(self, offset):
        
        for i in range(len(self.vertexIndices)): 
            vertexIndex = self.vertexIndices[i];
            vertexIndex += offset;
            self.vertexIndices[i] = vertexIndex;
            
        return self;
    
    # Clonable.
    
    def clone(self):
        return MeshFace(
            list(
                map(
                    lambda x : x,
                    self.vertexIndices
                )
            )
        )


class MeshWithMaterial:
    def __init__(self, mesh, material):
        self.mesh = mesh.clone()
        self.material = material
        
        self.name = mesh.name + "_" + material.name
        
    def renderToBpy(self, bpy):
        self.systemMesh = self.mesh.renderToBpy(bpy)        
        return self.systemMesh
    
    def renderToBpyAndObject(self, bpy, systemObject):
        systemMaterial = self.material.renderToBpy(bpy)
        systemObject.data.materials.append(
            systemMaterial
        )
        return systemObject

            
class Orientation:
    def __init__(self, forward, down):
        self.forward = forward
        self.down = down

        self.right = Coords.zeroes()
        self.normalize()
        self.orthogonalize()
        
    def default():
        return Orientation(
            Coords(1, 0, 0),
            Coords(0, 0, 1)
        )

    def Instances():            
         # hack - Don't know how to do singletons in Python.
        return Orientation_Instances()
    
    def normalize(self):
        self.forward.normalize()
        self.down.normalize()
        self.right.normalize()
        return self

    def orthogonalize(self):        
        self.right.overwriteWith(self.forward).crossProduct(self.down)
        return self
    
    def toEulerRotationTuple(self):
        returnValue = (0, 0, 0)

        # hack - This needs to be more general.
        if self.down.z == 1:
            forward = self.forward
            zAxisRotationInRadians = math.atan2(forward.y, forward.x)
            returnValue = (0, 0, zAxisRotationInRadians)
        
        return returnValue
    
    def toString():
        return ("Fd:" + self.forward.toStringXYZ() + ", Dn:" + self.down.toStringXYZ() )
    
    
class Orientation_Instances:
    def __init__(self):
        self.ForwardXNegativeDownZ = Orientation(
            Coords(-1, 0, 0),
            Coords(0, 0, 1)
        )
        self.ForwardYNegativeDownZ = Orientation(
            Coords(0, -1, 0),
            Coords(0, 0, 1)
        )
        self.ForwardXDownZ = Orientation(
            Coords(1, 0, 0),
            Coords(0, 0, 1)
        )
        self.ForwardYDownZ = Orientation(
            Coords(0, 1, 0),
            Coords(0, 0, 1)
        )
        

class Scene:
    def __init__(self, meshDispositions):
        self.meshDispositions = meshDispositions
        
    def renderToBpy(self, bpy):
        for meshDisp in self.meshDispositions:
            systemMesh = meshDisp.renderToBpy(bpy)
            
            
class Transform_Displace:
    def __init__(self, displacement):
        self.displacement = displacement
        
    def transformCoords(self, coordsToTransform):
        coordsToTransform.add(self.displacement)
        return coordsToTransform


class Transform_Orient:
    def __init__(self, orientation):
        self.orientation = orientation
        
        self._components = [Coords.zeroes(), Coords.zeroes(), Coords.zeroes()]
        
    def transformCoords(self, coordsToTransform):
        components = self._components;
        ori = self.orientation;

        coordsToTransform.overwriteWith
        (
            components[0].overwriteWith(ori.forward).multiplyScalar(coordsToTransform.x).add
            (
                components[1].overwriteWith(ori.right).multiplyScalar(coordsToTransform.y).add
                (
                    components[2].overwriteWith(ori.down).multiplyScalar(coordsToTransform.z)
                )
            )
        )

        return coordsToTransform


class Transform_Scale:
    def __init__(self, scaleFactors):
        self.scaleFactors = scaleFactors
        
    def transformCoords(self, coordsToTransform):
        coordsToTransform.multiply(self.scaleFactors)
        return coordsToTransform


# Main.
        
def main():
    
    print("Script begins.")
    
    context = bpy.context
    scene = context.scene

    objects = scene.objects

    cursor = scene.cursor
    cursor.location = (0.0, 0.0, 0.0)

    meshBuilder = MeshBuilder()
    
    meshSize = Coords(4, 2, 1)
    meshBrick2x4 = meshBuilder.brickWithTopStuds(meshSize)

    colors = Color.Instances()
    
    materialBlack = Material("Black", colors.Black)
    materialBlue = Material("Blue", colors.Blue)
    materialBrown = Material("Brown", colors.Brown)
    materialCyan = Material("Cyan", colors.Cyan)
    materialGrayDark = Material("GrayDark", colors.GrayDark)
    materialGrayLight = Material("GrayLight", colors.GrayLight)
    materialGreen = Material("Green", colors.Green)
    materialOrange = Material("Orange", colors.Orange)
    materialRed = Material("Red", colors.Red)
    materialViolet = Material("Violet", colors.Violet)
    materialWhite = Material("White", colors.White)
    materialYellow = Material("Yellow", colors.Yellow)

    brick2x4Black = MeshWithMaterial(meshBrick2x4, materialBlack)
    brick2x4Blue = MeshWithMaterial(meshBrick2x4, materialBlue)
    brick2x4Brown = MeshWithMaterial(meshBrick2x4, materialBrown)
    brick2x4Cyan = MeshWithMaterial(meshBrick2x4, materialCyan)
    brick2x4GrayDark = MeshWithMaterial(meshBrick2x4, materialGrayDark)
    brick2x4GrayLight = MeshWithMaterial(meshBrick2x4, materialGrayLight)
    brick2x4Green = MeshWithMaterial(meshBrick2x4, materialGreen)
    brick2x4Orange = MeshWithMaterial(meshBrick2x4, materialOrange)
    brick2x4Red = MeshWithMaterial(meshBrick2x4, materialRed)
    brick2x4Violet = MeshWithMaterial(meshBrick2x4, materialViolet)
    brick2x4White = MeshWithMaterial(meshBrick2x4, materialWhite)
    brick2x4Yellow = MeshWithMaterial(meshBrick2x4, materialYellow)

    brick2x4 = brick2x4GrayLight

    meshGear24 = meshBuilder.gear(24, 0.9).transform(
        Transform_Scale(
            Coords(2, 2, .25)
        )
    )
    gear24 = MeshWithMaterial(meshGear24, materialGrayLight)
    
    meshAxle8 = meshBuilder.star(4, 0.5).transform(
        Transform_Scale(
            Coords(.5, .5, 8)
        )
    )
    axle8 = MeshWithMaterial(meshAxle8, materialGrayDark)
            
    oris = Orientation.Instances()
    # hack - The axes are wrong here.
    # "e" = east, "n" = north, "w" = west, "s" = south.
    e = oris.ForwardXDownZ
    n = oris.ForwardYDownZ
    w = oris.ForwardXNegativeDownZ
    s = oris.ForwardYNegativeDownZ
        
    # Convenience method to cut down on typing.
    def md(mesh, x, y, z, orientation):
        return MeshDisposition(mesh, Disposition(Coords(x, y, z), orientation))
            
    meshDispositions = [

        # Tier 1.
        md(brick2x4Red, 0, 0, 0, e),
        md(brick2x4Orange, 4, 2, 0, n),
        md(brick2x4Yellow, 6, -2, 0, w),
        md(brick2x4Green, 2, -4, 0, s),

        # Tier 2.
        md(brick2x4Cyan, 2, 0, 1, e),
        md(brick2x4Blue, 4, 0, 1, n),
        md(brick2x4Violet, 4, -2, 1, w),
        md(brick2x4Black, 2, -2, 1, s),

        # Tier 3.
        md(brick2x4GrayLight, 0, 0, 2, e),
        md(brick2x4GrayDark, 4, 2, 2, n),
        md(brick2x4White, 6, -2, 2, w),
        md(brick2x4, 2, -4, 2, s),
        
        md(axle8, 3, -1, 10, e),
        md(gear24, 3, -1, 9, e),
    ]
    
    scene = Scene(meshDispositions)
    scene.renderToBpy(bpy)
    
    print("Script ends.")
    
    # To do outline shading:
    # Add a Stroke object.
    # Go into Edit mode and delete the stroke's points.
    # In the Scene Outline view, move the Stroke under the main Scene node.
    # Move all the things you want to shade under the Collection node.
    # Add a "Line Art" modifier to the Stroke.
    # Choose the "Collection" node in the Collection box.
    # Choose "Lines" in the Layer box.
    # Choose the material "Black.001" in the Material box.

main()

    