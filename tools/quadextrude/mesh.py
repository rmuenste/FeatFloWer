import operator
#/usr/bin/env python
# vim: set filetype=python
"""
A module for mesh related classes and functions
"""

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
    def __init__(self, nodeIds, idx):
        self.nodeIds = nodeIds
        self.idx = idx
        self.traits = []

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
    def __init__(self, quads, nodes):
        self.quads = quads
        self.nodes = nodes

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
    def __init__(self, hexas=[], nodes=[]):
        self.hexas = hexas
        self.nodes = nodes
        self.nodesLayer = 0
        self.hexasLayer = 0
        self.elementsAtVertexIdx = []
        self.elementsAtVertex = []
        self.verticesAtBoundary = []
        # GenElAtVert()
        # GenElAtVert()


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
    for n in quadMesh.nodes:
        coords = (n[0], n[1], n[2], n[3])
        hexNodes.append([float(coords[1]), float(coords[2]), float(coords[3])])

    for q in quadMesh.quads:
        nodeIds = (int(q[5])-1, int(q[6])-1, int(q[7])-1, int(q[8])-1)
        hexHexas.append([nodeIds[0], nodeIds[1], nodeIds[2], nodeIds[3]])

    totalNodes = len(hexNodes)
    for n in quadMesh.nodes:
        coords = (n[0], n[1], n[2], n[3])
        hexNodes.append([float(coords[1]), float(coords[2]), float(coords[3])+dz])

    for idx, q in enumerate(quadMesh.quads):
        nodeIdsBot = (int(q[5])-1, int(q[6])-1, int(q[7])-1, int(q[8])-1)
        hexHexas[idx].append(nodeIdsBot[0]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[1]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[2]+totalNodes)
        hexHexas[idx].append(nodeIdsBot[3]+totalNodes)
        h = Hexa(hexHexas[idx], idx)
        h.layerIdx = 1
        h.type = int(q[4])
        realHexas.append(h)

    return HexMesh(realHexas, hexNodes)


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

    for hidx in range(len(hexMesh.hexas)):
        newNodeIds = []
        for idx in hexMesh.hexas[hidx].nodeIds:
            if not idx in nodeMap:
                nodeMap[idx] = nodeCounter
                newNodeIds.append(nodeCounter)
                newNodes.append(hexMesh.nodes[idx])
                nodeCounter = nodeCounter + 1
            else:
                newNodeIds.append(nodeMap[idx])

        hexMesh.hexas[hidx].nodeIds = newNodeIds

    hexMesh.nodes = newNodes


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
#                       Function generateElementsAtVertex
#===============================================================================
def generateElementsAtVertex(hexMesh):
    """
    Compute the Elements attached to a particular vertex

    Args:
        hexMesh: The input/output hex mesh
    """
    elemAtVertIdx = []

    for node in hexMesh.nodes:
        elemAtVertIdx.append(list())

    for idx, hex in enumerate(hexMesh.hexas):
        for node in hex.nodeIds:
            elemAtVertIdx[node].append(idx) 


    hexMesh.elementsAtVertex = elemAtVertIdx


#===============================================================================
#                       Function generateNeighborsAtElement
#===============================================================================
def generateNeighborsAtElement(hexMesh):
    """
    Compute the neighbors at the faces of an element
    Uses the connector data structure for a face:
    connector[6]
    connector[0-3] : the indices of the face
    connector[4] : the idx of the hexa the face was found in 
    connector[5] : the internal face index in the hexa 
                   (0 for the first face, 1 for the 2nd,...) 

    Args:
        hexMesh: The input/output hex mesh
    """

    connectorList = []
    for hidx, hex in enumerate(hexMesh.hexas):

        connector = []
        #=========================================================  
        # first face
        for i in range(4):  
            connector.append(hex.nodeIds[i])

        # The hexa idx
        connector.append(hidx)

        # The internal face idx
        connector.append(0)

        connectorList.append(connector)

        #=========================================================  
        # second face
        connector = []
        connector.append(hex.nodeIds[0])
        connector.append(hex.nodeIds[1])
        connector.append(hex.nodeIds[4])
        connector.append(hex.nodeIds[5])

        # The hexa idx
        connector.append(hidx)

        # The internal face idx
        connector.append(1)

        connectorList.append(connector)

        #=========================================================  
        # third face
        connector = []
        connector.append(hex.nodeIds[1])
        connector.append(hex.nodeIds[2])
        connector.append(hex.nodeIds[5])
        connector.append(hex.nodeIds[6])

        # The hexa idx
        connector.append(hidx)

        # The internal face idx
        connector.append(2)

        connectorList.append(connector)

        #=========================================================  
        # fourth face
        connector = []
        connector.append(hex.nodeIds[3])
        connector.append(hex.nodeIds[2])
        connector.append(hex.nodeIds[6])
        connector.append(hex.nodeIds[7])

        # The hexa idx
        connector.append(hidx)

        # The internal face idx
        connector.append(3)

        connectorList.append(connector)

        #=========================================================  
        # fifth face
        connector = []
        connector.append(hex.nodeIds[0])
        connector.append(hex.nodeIds[3])
        connector.append(hex.nodeIds[7])
        connector.append(hex.nodeIds[4])

        # The hexa idx
        connector.append(hidx)

        # The internal face idx
        connector.append(4)

        connectorList.append(connector)

        #=========================================================  
        # sixth face
        connector = []
        for i in range(4,8):  
            connector.append(hex.nodeIds[i])

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

    testList = connectorList[0:36] 

    for i in range(1,len(connectorList)):
        ca = connectorList[i-1]
        cb = connectorList[i]
        if connectorList[i-1][0:4] == connectorList[i][0:4]: 
            hexMesh.hexas[ca[4]].neighIdx[ca[5]] = cb[4]
            hexMesh.hexas[cb[4]].neighIdx[cb[5]] = ca[4]

 

#    print testList[13][0:4]
#    print testList[12][0:4]
#    if testList[11][0:4] == testList[13][0:4]:
#        print("ja")
#    else:
#        print("nein")

#    print("Hexas: " + str(6 * len(hexMesh.hexas)))
#    print(str(len(connectorList)))

#===============================================================================
#                       Function writeParFiles
#===============================================================================
def writeParFiles(hexMesh):
    """
    Writes a list of .par files from the hexa typeIds

    Args:
        hexMesh: The input/output hex mesh
    """
    typeIds = []
    parDict = {}

    faceIndices = [[0, 1, 2, 3], [0, 1, 4, 5], [1, 2, 5, 6], [3, 2, 6, 7], [0, 3, 7, 4], [4, 5, 6, 7]]

    hexMesh.verticesAtBoundary = [0] * len(hexMesh.nodes)

    print("Number of nodes: " + str(len(hexMesh.nodes)))

    for h in hexMesh.hexas:
        if not int(h.type) in typeIds:
            typeIds.append(h.type)
            parDict.update({int(h.type) : set()})

    for h in hexMesh.hexas:
        for node in h.nodeIds:
            parDict[h.type].add(int(node))

    for hidx, h in enumerate(hexMesh.hexas):
        for idx, item in enumerate(h.neighIdx):
            if item == -1:
                bndryVertices = [h.nodeIds[faceIndices[idx][0]], 
                                 h.nodeIds[faceIndices[idx][1]],
                                 h.nodeIds[faceIndices[idx][2]],
                                 h.nodeIds[faceIndices[idx][3]]]
                for node in bndryVertices:
                    parDict[h.type].add(int(node))
                    hexMesh.verticesAtBoundary[node] = h.type

    parFileNames = []

    for key in parDict:
        parList = parDict[key]
        parName = "parfile" + str(key) + ".par"
        parFileNames.append(parName)
        with open(parName, "w") as parFile:
            parFile.write(str(len(parList)) + " Wall")
            parFile.write("' '")
            for nodeIdx in parList:
                parFile.write(str(nodeIdx + 1) + "\n")

    with open("file.prj", "w") as prjFile:
        prjFile.write("mesh.tri\n")
        for name in parFileNames:
            prjFile.write(name + "\n")