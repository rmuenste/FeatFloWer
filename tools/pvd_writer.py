
import xml.etree.ElementTree as ET 
from xml.dom import minidom 
from pathlib import Path

import sys


workingDir = Path('./gmv/')
vtkList = list(workingDir.glob('u.vtk.*'))
vtkList.sort()

root = ET.Element('VTKFile') 

root.set('type', 'Collection')
root.set('version', '0.1')
root.set('byte_order', 'LittleEndian')
#root.set('compressor', 'vtkZLibDataCompressor')

collectionNode = ET.SubElement(root, 'Collection')

time = 0
dt = 0.1
for item in vtkList:
    dataSetNode = ET.SubElement(collectionNode, 'DataSet')
    dataSetNode.set('file', str(item))
    dataSetNode.set('group', '')
    dataSetNode.set('part', '0')
    dataSetNode.set('timestep', str(time))
    time = time + dt

xmlstr = minidom.parseString(ET.tostring(root, method='xml')).toprettyxml(indent=" ")
print(xmlstr)
