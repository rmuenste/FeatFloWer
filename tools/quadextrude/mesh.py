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
