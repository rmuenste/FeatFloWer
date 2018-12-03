
import xml.etree.ElementTree as ET 
from xml.dom import minidom 

root = ET.Element('VTKFile') 

root.set('type', 'Collection')
root.set('version', '0.1')
root.set('byte_order', 'LittleEndian')
root.set('compressor', 'vtkZLibDataCompressor')

collectionNode = ET.SubElement(root, 'Collection')


dataSetNode = ET.SubElement(collectionNode, 'DataSet')
dataSetNode.set('file','test0001.vtp')
dataSetNode.set('group','')
dataSetNode.set('part','0')
dataSetNode.set('timestep','1')

xmlstr = minidom.parseString(ET.tostring(root, method='xml')).toprettyxml(indent=" ")
print(xmlstr)
