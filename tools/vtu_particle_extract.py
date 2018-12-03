import numpy as np
import xml.etree.ElementTree as ET
tree = ET.parse('particles.vtu')
root = tree.getroot()

print(root.tag, root.attrib)

pointsNode = ""
cellNode = ""

connectivity = ""
offset = ""

pieceNode = root.find('UnstructuredGrid').find('Piece')
for child in pieceNode:
    print(child.tag)
    if child.tag == "Points":
        pointsNode = child
    elif child.tag == "Cells": 
        cellNode = child

pointsData = pointsNode.find('DataArray')
pointsList = pointsData.text.split() 
thePoints = [float(x) for x in pointsList]

for child in cellNode:
    if child.attrib['Name'] == "connectivity":
        connectivity = child.text.split()
    elif child.attrib['Name'] == "offsets": 
        offsets = child.text.split()

connect = [int(x) for x in connectivity]
connect = zip(*[iter(connect)]*3)
connect = zip(*[iter(connect)]*224)
offsetIdx = [int(x) for x in offset]

totalCells = len(connect)
print(totalCells)

totalPoints = len(pointsList)/3
print(totalPoints)

vertices = zip(*[iter(thePoints)]*3)
vertices = zip(*[iter(vertices)]*224)

for sphere in vertices:
    center = np.zeros(3)
    for point in sphere:
        vertex = np.array([point[0], point[1], point[2]])
        center = center + vertex

    center = center * 1.0/224.0
    print(center)