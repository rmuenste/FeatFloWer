#/usr/bin/env python
# vim: set filetype=python
"""
A module for mesh related classes and functions
"""
import operator
import sys
import numpy as np

#===============================================================================
#                          A class for a quad
#===============================================================================
class Quad:
    """
    A class for a quad element

    Attributes:
        nodeIds: A list of the node Ids (indices) of each vertex of the
                 quad
        idx: The index of the quad element
    """
    def __init__(self, nodeIds, zoneId, idx):
        self.nodeIds = nodeIds
        self.idx = idx
        self.zoneId = zoneId
        self.traits = []
        self.type = 0
        self.neighIdx = [-1] * 4
        self.hasBoundaryFace = False

#===============================================================================
#                          A class for an edge
#===============================================================================
class Edge:
    """
    A class for an edge of a quadrilateral element

    Attributes:
        nodeIds: A list of the node Ids (indices) of each vertex of the
                 edge
        layerIdx: The level that the edge belongs to
        type: The type id of the edge element
    """
    def __init__(self, nodeIds, idx, edgeType='UNINITIALIZED'):
        self.nodeIds = nodeIds
        self.idx = idx
        self.edgeType = edgeType


#===============================================================================
#                          A class for a hexa
#===============================================================================
class Hexa:
    """
    A class for a hexahedral element

    Attributes:
        nodeIds: A list of the node Ids (indices) of each vertex of the
                 hexahedron
        layerIdx: The level that the hexahedron belongs to
        type: The type id of the hexahedral element
    """
    def __init__(self, nodeIds, idx):
        self.nodeIds = nodeIds
        self.idx = idx
        self.layerIdx = 0
        self.type = 0
        self.neighIdx = [-1] * 6
        self.hasBoundaryFace = False


#===============================================================================
#                          A class for a face
#===============================================================================
class Face:
    """
    A class for face of a hexahedral element

    Attributes:
        nodeIds: A list of the node Ids (indices) of each vertex of the
                 hexahedron
        layerIdx: The level that the hexahedron belongs to
        type: The type id of the hexahedral element
    """
    def __init__(self, nodeIds, idx, faceType='UNINITIALIZED'):
        self.nodeIds = nodeIds
        self.idx = idx
        self.layerIdx = 0
        self.faceType = faceType


#===============================================================================
#                          A class for a QuadMesh
#===============================================================================
class QuadMesh:
    """
    A class for a quad mesh

    Attributes:
        quads: A list of the quad element, each element of the list is
        another list with the vertex indices of the quad

        nodes: A list of the coordinates of the nodes
    """
    def __init__(self, nodes, elems):
        self.nodes = nodes[:]
        self.elements = elems[:] 
        self.elementsAtVertex = []
        self.verticesAtBoundary = []
        self.edgesAtBoundary = []
        self.area = []
        #self.area = [0.0] * len(elems)


#===============================================================================
#                       Function checkMeshValidity1
#===============================================================================
    def checkMeshValidity1(self):
        """
          Check the connectivity of the mesh
        """

        for idx, node in enumerate(self.nodes):
            if len(self.elementsAtVertex[idx]) <= 2 and self.verticesAtBoundary[idx] == 0:
                print("Vertex index: " + str(idx) + "," + " Elements at vertex: " + 
                    str(len(self.elementsAtVertex[idx])) + " " +  
                    str(self.verticesAtBoundary[idx])
                    )
                print("An inner vertex with <= 2 incident quad elements was found. This mesh topology is invalid. Exiting") 
                sys.exit(2)


#===============================================================================
#                       Function quadArea
#===============================================================================
    def quadArea(self):
        """
          Calculate the area of each quad in the mesh
        """
        for quad in (self.elements):

            triangleA = (quad.nodeIds[0], quad.nodeIds[1], quad.nodeIds[2])
            triangleB = (quad.nodeIds[0], quad.nodeIds[3], quad.nodeIds[2])

            p0 = np.array(list(self.nodes[triangleA[0]]))  
            p1 = np.array(list(self.nodes[triangleA[1]]))  
            p2 = np.array(list(self.nodes[triangleA[2]]))  

            area1 = getTriangleArea(p0, p1, p2)

            p0 = np.array(list(self.nodes[triangleB[0]]))  
            p1 = np.array(list(self.nodes[triangleB[1]]))  
            p2 = np.array(list(self.nodes[triangleB[2]]))  

            area2 = getTriangleArea(p0, p1, p2)
            self.area.append(area1 + area2)
#            print("Idx: " + str(count) + " old idx: " + str(quad.idx))
#            print("Area1: " + str(area1))
#            print("Area2: " + str(area2))
#            count = count + 1

        minArea = np.min(self.area)
        maxArea = np.max(self.area)
        factor = maxArea/minArea

        print("=======Quad mesh element area distribution=======")
        print("Minimum Area: " + str(minArea))
        print("Maximum Area: " + str(maxArea))
        if factor < 1000.0: 
            print("Factor " + str(factor) + " < 1000.0 -> OK")
            return True
        else:
            print("Factor " + str(factor) + " >= 1000.0 -> Exiting")
            print("The input quad mesh failed the element area distribution criterion.")
            return False
            #sys.exit(2)


#===============================================================================
#                       Function generateElementsAtVertex
#===============================================================================
    def generateElementsAtVertex(self):
        """
        Compute the Elements attached to a particular vertex

        Args:
            hexMesh: The input/output hex mesh
        """
        elemAtVertIdx = []

        for node in self.nodes:
            elemAtVertIdx.append(list())

        for idx, quad in enumerate(self.elements):
            for node in quad.nodeIds:
                elemAtVertIdx[node].append(idx)


        self.elementsAtVertex = elemAtVertIdx

#===============================================================================
#                       Function generateNeighborsAtElement
#===============================================================================
    def generateNeighborsAtElement(self):
        """
        Compute the neighbors at the faces of an element
        Uses the connector data structure for a face:
        connector[4]
        connector[0-1] : the indices of the edge
        connector[2] : the idx of the quad the edge was found in
        connector[3] : the internal edge index in the quad
                       (0 for the first edge, 1 for the 2nd,...)

        Args:
            self: The input/output hex mesh
        """

        connectorList = []
        for qidx, quad in enumerate(self.elements):

            connector = []
            #=========================================================
            # first edge
            for i in range(2):
                connector.append(quad.nodeIds[i])

            # The quad idx
            connector.append(qidx)

            # The internal edge idx
            connector.append(0)

            connectorList.append(connector)

            #=========================================================
            # second edge
            connector = []
            connector.append(quad.nodeIds[1])
            connector.append(quad.nodeIds[2])

            # The quad idx
            connector.append(qidx)

            # The internal edge idx
            connector.append(1)

            connectorList.append(connector)

            #=========================================================
            # third edge
            connector = []
            connector.append(quad.nodeIds[2])
            connector.append(quad.nodeIds[3])

            # The edge idx
            connector.append(qidx)

            # The internal edge idx
            connector.append(2)

            connectorList.append(connector)

            #=========================================================
            # fourth edge
            connector = []
            connector.append(quad.nodeIds[3])
            connector.append(quad.nodeIds[0])

            # The hexa idx
            connector.append(qidx)

            # The internal edge idx
            connector.append(3)

            connectorList.append(connector)

        for connector in connectorList:
            connector[0:2] = sorted(connector[0:2])

        connectorList = sorted(connectorList, key=operator.itemgetter(1))
        connectorList = sorted(connectorList, key=operator.itemgetter(0))

        for i in range(1, len(connectorList)):
            neighA = connectorList[i-1]
            neighB = connectorList[i]
            if connectorList[i-1][0:2] == connectorList[i][0:2]:
                self.elements[neighA[2]].neighIdx[neighA[3]] = neighB[2]
                self.elements[neighB[2]].neighIdx[neighB[3]] = neighA[2]


#===============================================================================
#                       Function facesAtBoundary
#===============================================================================
    def generateEdgesAtBoundary(self):
        """
        Compute the boundary edges of the mesh

        Args:
            self: The input/output hex mesh
        """
        edgesAtBoundary = []

        edgeIndices = [[0, 1], [1, 2],
                       [2, 3], [3, 0]]

        nedges = 0
        for q in self.elements:
            for idx, item in enumerate(q.neighIdx):
                if item == -1:
                    q.hasBoundaryFace = True

                    bndryVertices = [q.nodeIds[edgeIndices[idx][0]],
                                     q.nodeIds[edgeIndices[idx][1]]]

                    bndryEdge = Edge(bndryVertices, nedges, 'boundaryEdge')
                    edgesAtBoundary.append(bndryEdge)
                    nedges = nedges + 1

        self.edgesAtBoundary = edgesAtBoundary


#===============================================================================
#                       Function verticesAtBoundary
#===============================================================================
    def generateVerticesAtBoundary(self):
        """
        Compute the boundary vertices of the mesh

        Args:
            self: The input/output hex mesh
        """
        verticesAtBoundarySet = set()
        self.verticesAtBoundary = [0] * len(self.nodes)

        for edge in self.edgesAtBoundary:
            verticesAtBoundarySet.add(edge.nodeIds[0])
            verticesAtBoundarySet.add(edge.nodeIds[1])


        for idx in range(len(self.nodes)):
            if idx in verticesAtBoundarySet:
                self.verticesAtBoundary[idx] = 1




#===============================================================================
#                    Function generateMeshStructures
#===============================================================================
    def generateMeshStructures(self):
        """
        Generate a standard set of neighborhood information structures

        """
        self.generateElementsAtVertex()
        self.generateNeighborsAtElement()
#        print(self.elements[0].neighIdx)
        self.generateEdgesAtBoundary()
        self.generateVerticesAtBoundary()





#===============================================================================
#                          A class for a HexMesh
#===============================================================================
class HexMesh:
    """
    A class for a hexa mesh

    Attributes:
        hexas: A list of the hexa elements
        nodes: A list of the coordinates of the nodes
        nodesLayer: The number of nodes on an extrusion layer
        hexasLayer: The number of hexas on an extrusion layer
    """
    def __init__(self, hexas=[], nodes=[], sliceIds=[]):
        self.hexas = hexas
        self.nodes = nodes
        self.nodesLayer = 0
        self.hexasLayer = 0
        self.elementsAtVertexIdx = []
        self.elementsAtVertex = []
        self.verticesAtBoundary = []
        self.facesAtBoundary = []
        self.nodesAtSlice = sliceIds
        self.slice = 0
        # GenElAtVert()
        # GenElAtVert()


#===============================================================================
#                       Function generateElementsAtVertex
#===============================================================================
    def generateElementsAtVertex(self):
        """
        Compute the Elements attached to a particular vertex

        Args:
            hexMesh: The input/output hex mesh
        """
        elemAtVertIdx = []

        for node in self.nodes:
            elemAtVertIdx.append(list())

        for idx, hexa in enumerate(self.hexas):
            for node in hexa.nodeIds:
                elemAtVertIdx[node].append(idx)


        self.elementsAtVertex = elemAtVertIdx


#===============================================================================
#                       Function facesAtBoundary
#===============================================================================
    def generateFacesAtBoundary(self):
        """
        Compute the boundary faces of the mesh

        Args:
            self: The input/output hex mesh
        """
        facesAtBoundary = []

        faceIndices = [[0, 1, 2, 3], [0, 1, 4, 5],
                       [1, 2, 5, 6], [3, 2, 6, 7],
                       [0, 3, 7, 4], [4, 5, 6, 7]]

        nfaces = 0
        for h in self.hexas:
            for idx, item in enumerate(h.neighIdx):
                if item == -1:
                    h.hasBoundaryFace = True

                    bndryVertices = [h.nodeIds[faceIndices[idx][0]],
                                     h.nodeIds[faceIndices[idx][1]],
                                     h.nodeIds[faceIndices[idx][2]],
                                     h.nodeIds[faceIndices[idx][3]]]

                    bndryFace = Face(bndryVertices, nfaces, 'boundaryFace')
                    bndryFace.layerIdx = h.layerIdx
                    facesAtBoundary.append(bndryFace)
                    nfaces = nfaces + 1

        self.facesAtBoundary = facesAtBoundary


#===============================================================================
#                       Function verticesAtBoundary
#===============================================================================
    def generateVerticesAtBoundary(self):
        """
        Compute the boundary vertices of the mesh

        Args:
            self: The input/output hex mesh
        """
        verticesAtBoundarySet = set()
        self.verticesAtBoundary = [0] * len(self.nodes)

        for face in self.facesAtBoundary:
            verticesAtBoundarySet.add(face.nodeIds[0])
            verticesAtBoundarySet.add(face.nodeIds[1])
            verticesAtBoundarySet.add(face.nodeIds[2])
            verticesAtBoundarySet.add(face.nodeIds[3])


        for idx in range(len(self.nodes)):
            if idx in verticesAtBoundarySet:
                self.verticesAtBoundary[idx] = 1


#===============================================================================
#                       Function generateNeighborsAtElement
#===============================================================================
    def generateNeighborsAtElement(self):
        """
        Compute the neighbors at the faces of an element
        Uses the connector data structure for a face:
        connector[6]
        connector[0-3] : the indices of the face
        connector[4] : the idx of the hexa the face was found in
        connector[5] : the internal face index in the hexa
                       (0 for the first face, 1 for the 2nd,...)

        Args:
            self: The input/output hex mesh
        """

        connectorList = []
        for hidx, hexa in enumerate(self.hexas):

            connector = []
            #=========================================================
            # first face
            for i in range(4):
                connector.append(hexa.nodeIds[i])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(0)

            connectorList.append(connector)

            #=========================================================
            # second face
            connector = []
            connector.append(hexa.nodeIds[0])
            connector.append(hexa.nodeIds[1])
            connector.append(hexa.nodeIds[4])
            connector.append(hexa.nodeIds[5])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(1)

            connectorList.append(connector)

            #=========================================================
            # third face
            connector = []
            connector.append(hexa.nodeIds[1])
            connector.append(hexa.nodeIds[2])
            connector.append(hexa.nodeIds[5])
            connector.append(hexa.nodeIds[6])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(2)

            connectorList.append(connector)

            #=========================================================
            # fourth face
            connector = []
            connector.append(hexa.nodeIds[3])
            connector.append(hexa.nodeIds[2])
            connector.append(hexa.nodeIds[6])
            connector.append(hexa.nodeIds[7])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(3)

            connectorList.append(connector)

            #=========================================================
            # fifth face
            connector = []
            connector.append(hexa.nodeIds[0])
            connector.append(hexa.nodeIds[3])
            connector.append(hexa.nodeIds[7])
            connector.append(hexa.nodeIds[4])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(4)

            connectorList.append(connector)

            #=========================================================
            # sixth face
            connector = []
            for i in range(4, 8):
                connector.append(hexa.nodeIds[i])

            # The hexa idx
            connector.append(hidx)

            # The internal face idx
            connector.append(5)

            connectorList.append(connector)

        for connector in connectorList:
            connector[0:4] = sorted(connector[0:4])

        connectorList = sorted(connectorList, key=operator.itemgetter(3))
        connectorList = sorted(connectorList, key=operator.itemgetter(2))
        connectorList = sorted(connectorList, key=operator.itemgetter(1))
        connectorList = sorted(connectorList, key=operator.itemgetter(0))

        for i in range(1, len(connectorList)):
            neighA = connectorList[i-1]
            neighB = connectorList[i]
            if connectorList[i-1][0:4] == connectorList[i][0:4]:
                self.hexas[neighA[4]].neighIdx[neighA[5]] = neighB[4]
                self.hexas[neighB[4]].neighIdx[neighB[5]] = neighA[4]


#===============================================================================
#                    Function generateMeshStructures
#===============================================================================
    def generateMeshStructures(self):
        """
        Generate a standard set of neighborhood information structures

        """
        self.generateElementsAtVertex()
        self.generateNeighborsAtElement()
        self.generateFacesAtBoundary()
        self.generateVerticesAtBoundary()


#===============================================================================
#                     Function extrudeQuadToHexMesh
#===============================================================================
def extrudeQuadMeshToHexMesh(quadMesh, dz):
    """
    Extrudes a hex mesh from the quad base layer

    Args:
        quadMesh: The input quad mesh
        dz: The extrusion length
    Returns:
        The extruded hex mesh
    """
    hexNodes = []
    hexHexas = []
    realHexas = []
    nodeSliceIds = []
    for n in quadMesh.nodes:
        coords = (n[0], n[1], n[2])
        hexNodes.append([float(coords[0]), float(coords[1]), float(coords[2])])
        nodeSliceIds.append(0)


    for q in quadMesh.elements:
        nodeIds = (q.nodeIds[0], q.nodeIds[1], q.nodeIds[2], q.nodeIds[3])
        hexHexas.append([nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3]])

    totalNodes = len(hexNodes)
    for n in quadMesh.nodes:
        coords = (n[0], n[1], n[2])
        hexNodes.append([float(coords[0]), float(coords[1]), float(coords[2])+dz])
        nodeSliceIds.append(1)

    for idx, q in enumerate(quadMesh.elements):
        nodeIdsBot = (q.nodeIds[0], q.nodeIds[1], q.nodeIds[2], q.nodeIds[3])
        hexHexas[idx].append(nodeIdsBot[0]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[1]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[2]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[3]+totalNodes)
        h = Hexa(hexHexas[idx], idx)
        h.layerIdx = 1
        h.type = q.zoneId 
        realHexas.append(h)

    return HexMesh(realHexas, hexNodes, nodeSliceIds)


#===============================================================================
#                     Function extrudeHexMeshZ
#===============================================================================
def extrudeHexMeshZ(hexMesh, offsetHex, offsetNodes, layerNodes, layerIdx, dz):
    """
    Extrudes another layer from the 'top' of the hex mesh

    Args:
        hexMesh: The input/output hex mesh
        offsetHex: The index where to add the new hexas
        offsetNodes: The index where to add the new nodes
        layerNodes: The number of nodes on an extrusion layer
        layerIdx: The index of the extrusion level
        dz: The extrusion length of the subdivision layer
    """
    hexNodes = []
    hexHexas = []
    newHexas = []

    hexMesh.slice = hexMesh.slice + 1

    for nidx in range(offsetNodes, len(hexMesh.nodes)):
        coords = (hexMesh.nodes[nidx][0], hexMesh.nodes[nidx][1], hexMesh.nodes[nidx][2])
        hexNodes.append([coords[0], coords[1], coords[2]+dz])

    for hidx in range(offsetHex, len(hexMesh.hexas)):
        nodeIds = [hexMesh.hexas[hidx].nodeIds[4], hexMesh.hexas[hidx].nodeIds[5],
                   hexMesh.hexas[hidx].nodeIds[6], hexMesh.hexas[hidx].nodeIds[7],
                   hexMesh.hexas[hidx].nodeIds[4]+layerNodes,
                   hexMesh.hexas[hidx].nodeIds[5]+layerNodes,
                   hexMesh.hexas[hidx].nodeIds[6]+layerNodes,
                   hexMesh.hexas[hidx].nodeIds[7]+layerNodes]
        h = Hexa(nodeIds, hidx)
        h.type = hexMesh.hexas[hidx].type
        h.layerIdx = layerIdx
        hexMesh.hexas.append(h)

    for n in hexNodes:
        hexMesh.nodes.append(n)
        hexMesh.nodesAtSlice.append(hexMesh.slice)


#===============================================================================
#                         Function renumberNodes
#===============================================================================
def renumberNodes(hexMesh):
    """
    Applies a renumbering algorithm to the mesh. As a side consequence
    it removes all nodes that are not connected to a hexahedron

    Args:
        hexMesh: The input/output hex mesh
    """
    nodeMap = {}
    nodeCounter = 0
    newNodes = []
    newNodesAtSlice = []

    for hidx in range(len(hexMesh.hexas)):
        newNodeIds = []
        for idx in hexMesh.hexas[hidx].nodeIds:
            if idx not in nodeMap:
                nodeMap[idx] = nodeCounter
                newNodeIds.append(nodeCounter)
                newNodes.append(hexMesh.nodes[idx])
                newNodesAtSlice.append(hexMesh.nodesAtSlice[idx])
                nodeCounter = nodeCounter + 1
            else:
                newNodeIds.append(nodeMap[idx])

        hexMesh.hexas[hidx].nodeIds = newNodeIds

    hexMesh.nodes = newNodes
    hexMesh.nodesAtSlice = newNodesAtSlice


#===============================================================================
#                         Function extractQuadMesh
#===============================================================================
def extractQuadMesh(quadMesh, zoneId):
    """
    This function extracts a quad mesh from another quad mesh.
    We are dealing with topological quad meshes that can contain
    several sub-meshes that are associated with different zoneIds(=sub-mesh ID).
    This function takes such a mesh as an input and returns a quad mesh
    that only consists of the sub-mesh with a certain ID.

    Args:
        quadMesh: The input topological quad mesh
        zoneId: The ID of the sub-mesh that should be extracted
    """
    nodeMap = {}
    nodeCounter = 0
    newNodes = []

    newQuads = []

    old_idx = 0
    for q in quadMesh.elements:
        if q.zoneId == zoneId:
            newQuad = Quad(q.nodeIds, q.zoneId, old_idx)
            newQuads.append(newQuad)
        old_idx = old_idx + 1

    newQuadMesh = QuadMesh(quadMesh.nodes, newQuads)

    for qidx in newQuadMesh.elements:
        newNodeIds = []
        for idx in qidx.nodeIds:
            if idx not in nodeMap:
                nodeMap[idx] = nodeCounter
                newNodeIds.append(nodeCounter)
                newNodes.append(newQuadMesh.nodes[idx])
                nodeCounter = nodeCounter + 1
            else:
                newNodeIds.append(nodeMap[idx])

        qidx.nodeIds = newNodeIds

    newQuadMesh.nodes = newNodes

    return newQuadMesh


#===============================================================================
#                       Function removeHexasLayer
#===============================================================================
def removeHexasLayer(hexMesh, levelIdx, typeIds):
    """
    Removes all hexas of a given type on a given level

    Args:
        hexMesh: The input/output hex mesh
        levelIdx: The level that should be processed
        typeIds: A list of the types that should be kept on the level
    """
    newHex = []

    for h in hexMesh.hexas:
        if h.layerIdx == levelIdx:
            if int(h.type) in typeIds:
                newHex.append(h)
        else:
            newHex.append(h)

    hexMesh.hexas = newHex


#===============================================================================
#                       Function parFileFromSlice
#===============================================================================
def parFileFromSlice(hexMesh, sliceId):
    """
    Builds a par file from a hex mesh and a slice

    Args:
        hexMesh: The input/output hex mesh
        sliceId: The id of the slice
    """
    layerOnePar = []
    for idx in range(len(hexMesh.nodes)):
        if hexMesh.nodesAtSlice[idx] == sliceId and hexMesh.verticesAtBoundary[idx] != 0:
            layerOnePar.append(idx)

    return layerOnePar


#===============================================================================
#                       Function getTriangleArea
#===============================================================================
def getTriangleArea(p0, p1, p2):
    
    c = p1 - p0  
    b = p2 - p0  

    cl = np.linalg.norm(c)
    bl = np.linalg.norm(b)
    
    coAlpha = np.dot(c, b)/(cl * bl)

    alpha = np.arccos(coAlpha)

    area = 0.5 * cl * bl * np.sin(alpha)

    return area


#===============================================================================
#                       Function writeSingleParFile
#===============================================================================
def writeSingleParFile(nodeIds, fileName, bndryType, dirName):
    """
    Builds a par file from a hex mesh and a slice

    Args:
        hexMesh: The input/output hex mesh
        sliceId: The id of the slice
    """

    parName = fileName
    with open(dirName + "/" + parName, "w") as parFile:
        parFile.write(str(len(nodeIds)) + " " + bndryType + "\n")
        parFile.write("' '\n")
        for nodeIdx in nodeIds:
            parFile.write(str(nodeIdx + 1) + "\n")


#===============================================================================
#                       Function writeParFiles
#===============================================================================
def writeParFiles(hexMesh, slicesOnLevel, dirName):
    """
    Writes a list of .par files from the hexa typeIds

    Args:
        hexMesh: The input/output hex mesh
    """
    typeIds = []
    parDict = {}

    print("Number of nodes: " + str(len(hexMesh.nodes)))

    for hexa in hexMesh.hexas:
        if int(hexa.type) not in typeIds:
            typeIds.append(hexa.type)
            parDict.update({int(hexa.type) : set()})

    for hexa in hexMesh.hexas:
        for node in hexa.nodeIds:
            parDict[hexa.type].add(int(node))

    #===============================================
    layerOnePar = parFileFromSlice(hexMesh, 0)

    parFileNames = []

    parName = "Inflow.par"
    parFileNames.append(parName)

    writeSingleParFile(layerOnePar, parName, "Inflow81", dirName)

    #===============================================
    lowerCylTopCapIdx = slicesOnLevel[1][0]
    lowerCylTopCap = parFileFromSlice(hexMesh, lowerCylTopCapIdx)
    parName = "plane" + str(lowerCylTopCapIdx) + ".par"
    parFileNames.append(parName)

    writeSingleParFile(lowerCylTopCap, parName, "Wall", dirName)

    #===============================================
    layerOneCyl = []
    for face in hexMesh.facesAtBoundary:
        if face.layerIdx == 1:
            sliceIds = [hexMesh.nodesAtSlice[face.nodeIds[0]],
                        hexMesh.nodesAtSlice[face.nodeIds[1]],
                        hexMesh.nodesAtSlice[face.nodeIds[2]],
                        hexMesh.nodesAtSlice[face.nodeIds[3]]]

            sliceIdLower = [0] * 4
            sliceIdUpper = [slicesOnLevel[1][0]] * 4

            if sliceIds not in (sliceIdLower, sliceIdUpper):
                for node in face.nodeIds:
                    layerOneCyl.append(node)

    parName = "cyl1.par"
    parFileNames.append(parName)
    writeSingleParFile(layerOneCyl, parName, "Wall", dirName)

    #===============================================
    cylinderLayer = []
    for face in hexMesh.facesAtBoundary:
        if face.layerIdx == 2:
            for node in face.nodeIds:
                cylinderLayer.append(node)

    parName = "allcylinders.par"
    parFileNames.append(parName)
    writeSingleParFile(cylinderLayer, parName, "Wall", dirName)

    #===============================================
    topCylBottomCapIdx = slicesOnLevel[2][0]
    topCylBottomCap = parFileFromSlice(hexMesh, topCylBottomCapIdx)

    parName = "plane" + str(topCylBottomCapIdx) + ".par"
    parFileNames.append(parName)

    writeSingleParFile(topCylBottomCap, parName, "Wall", dirName)

    #===============================================
    topCylTopCapIdx = slicesOnLevel[3][0]
    topCylTopCap = parFileFromSlice(hexMesh, topCylTopCapIdx)

    parName = "plane" + str(topCylTopCapIdx) + ".par"
    parFileNames.append(parName)

    writeSingleParFile(topCylTopCap, parName, "Wall", dirName)

    #===============================================
    topIndex = slicesOnLevel[3][len(slicesOnLevel[3])-1]
    layerOutflow = parFileFromSlice(hexMesh, topIndex)

    parName = "plane" + str(topIndex) + ".par"
    parFileNames.append(parName)

    writeSingleParFile(layerOutflow, parName, "Outflow", dirName)

    #===============================================
    parName = "cyl3.par"
    parFileNames.append(parName)

    cylinder3Layer = []
    for face in hexMesh.facesAtBoundary:
        if face.layerIdx == 3:
            sliceIds = [hexMesh.nodesAtSlice[face.nodeIds[0]],
                        hexMesh.nodesAtSlice[face.nodeIds[1]],
                        hexMesh.nodesAtSlice[face.nodeIds[2]],
                        hexMesh.nodesAtSlice[face.nodeIds[3]]]

            sliceIdLower = [slicesOnLevel[2][0]] * 4
            sliceIdUpper = [slicesOnLevel[3][0]] * 4

            if sliceIds not in (sliceIdLower, sliceIdUpper):
                for node in face.nodeIds:
                    cylinder3Layer.append(node)

    writeSingleParFile(cylinder3Layer, parName, "Wall", dirName)

    #===============================================
    parName = "profile.par"
    parFileNames.append(parName)

    profileLayer = []
    for face in hexMesh.facesAtBoundary:
        if face.layerIdx == 4:
            sliceIds = [hexMesh.nodesAtSlice[face.nodeIds[0]],
                        hexMesh.nodesAtSlice[face.nodeIds[1]],
                        hexMesh.nodesAtSlice[face.nodeIds[2]],
                        hexMesh.nodesAtSlice[face.nodeIds[3]]]

            sliceIdUpper = [topIndex] * 4

            if sliceIds != sliceIdUpper:
                for node in face.nodeIds:
                    profileLayer.append(node)

    writeSingleParFile(profileLayer, parName, "Wall", dirName)

    with open(dirName + "/file.prj", "w") as prjFile:
        prjFile.write("mesh.tri\n")
        for name in parFileNames:
            prjFile.write(name + "\n")
