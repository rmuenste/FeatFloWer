 
import argparse
import math
import os
import sys
import json

import logging
from typing import List

script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_path)

import numpy as np

def getIndex(vector, steps, i):
    divide = vector[i]/steps[i]
    if divide%1 < 0.0001:
        return math.floor(divide)-0.5
    if divide%1 > 0.9999:
        return math.ceil(divide)+0.5
    return math.floor(divide)

def toIndex(divider):
    if divider%1 < 0.0001:
        return [math.floor(divider)-1, math.floor(divider)]
    elif divider%1 > 0.9999:
        return [math.ceil(divider)-1, math.ceil(divider)]
    return [math.floor(divider)]

def getVoxelIndexes(vector: List[int], startVoxel: List[int], steps: List[int]):
    dividers = (vector-startVoxel)/steps
    return list(map(toIndex, dividers))

def getVoxelIndex(vector: List[int], startVoxel: List[int], steps: List[int], dim: int):
    divider = (vector[dim]-startVoxel[dim])/steps[dim]
    if divider%1 < 0.0001:
        return [math.floor(divider)-1, math.floor(divider)]
    elif divider%1 > 0.9999:
        return [math.ceil(divider)-1, math.ceil(divider)]
    return [math.floor(divider)]



def faceArea(face):
    return 0.5*np.linalg.norm(np.cross(face[0]-face[1], face[0]-face[2]))


def common2(arr1, arr2):
    res = []
    for v1 in arr1:
        for v2 in arr2:
            if v1 == v2:
                res.append(v1)
    return res
   
def common(arrs):
    res = arrs[0]
    for arr in arrs[1:]:
        res = common2(res, arr)
    return res

def testFaceVoxel(face, startVoxel, steps, voxelData, correctionIndexes=[0,0,0], recDepth=0):
    tabs = '    '*recDepth
    #if recDepth>20:
    #    raise ValueError()
    logging.debug(tabs+'[{},{},{}],[{},{},{}],[{},{},{}]'.format(*face[0], *face[1], *face[2]))
    faceIndexes =  np.array([getVoxelIndexes(vertex, startVoxel, steps) for vertex in face])
    #print(faceIndexes)
    for dim in range(3):
        indexEdges = [[0,1,2],[0,2,1],[1,2,0]]
        for indexEdge in indexEdges:
            a1, a2, a3 = indexEdge
            i1 = faceIndexes[a1][dim]
            i2 = faceIndexes[a2][dim]
            #ix1 = getVoxelIndex(face[a1], startVoxel, steps, dim)
            #ix2 = getVoxelIndex(face[a2], startVoxel, steps, dim)
            #print(ix1, i1, ix2, i2)
            if len(common2(i1, i2))>0:
                # edge goes from one side to the other or is the same
                continue
            else: # edge goes over bound
                # compute bound
                #print(i1, i2, (face[a2][dim]-startVoxel[dim])/steps[dim])
                #print(sorted([i1,i2], key=lambda value: value[0]))
                a = i2[-1]+1 if i1[0]>i2[0] else i1[-1]+1 
                """
                a = max(sorted([i1,i2], key=lambda value: value[0])[0]) + 1
                if a != an:
                    print(i1, i2)
                    raise ValueError('{},{}'.format(an, a))
                """
                #print('split', indexEdge, dim, i1, i2, a)
                #print(face[a1][dim], steps[dim], getIndex(face[a1], steps, dim) )
                #print(face[a2][dim], steps[dim], getIndex(face[a2], steps, dim) )
                # create new vertex
                if dim == 0:
                    x = startVoxel[dim]+a*steps[dim]
                    y = (face[a1][1] + face[a2][1])/2
                    z = (face[a1][2] + face[a2][2])/2
                    
                elif dim == 1:
                    y = startVoxel[dim]+a*steps[dim]
                    x = (face[a1][0] + face[a2][0])/2
                    z = (face[a1][2] + face[a2][2])/2
                else:
                    z = startVoxel[dim]+a*steps[dim]
                    x = (face[a1][0] + face[a2][0])/2
                    y = (face[a1][1] + face[a2][1])/2
                vert = np.array([x,y,z])
                #print(vert)
                #create new faces
                nface1 = [face[a1], face[a3], vert]
                nface2 = [face[a2], face[a3], vert]
                
                logging.debug(tabs+'Split stay: {}, dim: {}, a: {}, value: {}'.format(a3, dim, a, startVoxel[dim]+a*steps[dim]))
                testFaceVoxel(nface1, startVoxel, steps, voxelData, correctionIndexes, recDepth=recDepth+1)
                testFaceVoxel(nface2, startVoxel, steps, voxelData, correctionIndexes, recDepth=recDepth+1)
                logging.debug(tabs+'Fin')
                return
    #face lies completly in a single voxel
    #indexesX = common([getVoxelIndex(face[i], startVoxel, steps, 0) for i in range(3)])
    #indexesY = common([getVoxelIndex(face[i], startVoxel, steps, 1) for i in range(3)])
    #indexesZ = common([getVoxelIndex(face[i], startVoxel, steps, 2) for i in range(3)])
    indexesX = common(faceIndexes[:,0])
    indexesY = common(faceIndexes[:,1])
    indexesZ = common(faceIndexes[:,2])
    #logging.debug(tabs, [getVoxelIndex(face[i], steps, 0) for i in range(3)])
    logging.debug(tabs+'{},{},{}'.format(indexesX, indexesY, indexesZ))
    amount = len(indexesX)*len(indexesY)*len(indexesZ)
    for iX in indexesX:
        for iY in indexesY:
            for iZ in indexesZ:
                iX -= correctionIndexes[0]
                iY -= correctionIndexes[1]
                iZ -= correctionIndexes[2]
                if iX<0 or iY<0 or iZ<0:
                    #if iZ>=0:
                    #print('boring', iX, iY, iZ)
                    continue
                if iX>=voxelData.shape[0] or iY>=voxelData.shape[1] or iZ>=voxelData.shape[2]:
                    #print('outside', iX, iY, iZ)
                    continue
                logging.debug(tabs+'{},{},{}:{}:{}'.format(iX, iY, iZ, amount, faceArea(face)))
                voxelData[iX,iY,iZ] += faceArea(face)/amount


def loadOffFile(filename, position=np.array([0.0,0.0,0.0])):
    pos    = np.array([np.inf, np.inf, np.inf])
    posMax = np.array([-np.inf, -np.inf, -np.inf])
    with open(filename, 'r')as f:
        header = f.readline()
        if header[:3] != 'OFF':
            raise ValueError('unsupported file type')
        info = f.readline()
        data = info.strip().split(' ')
        if len(data) != 3:
            raise ValueError('Off info line needs to have 3 inputs not {}'.format(len(data)))
        
        vertCount = int(data[0]) 
        faceCount = int(data[1])
        
        vertices = np.zeros((vertCount, 3))
        for i in range(vertCount):
            vertLine = f.readline()
            x, y, z = vertLine.strip().split(' ')
            vertices[i] = np.array([float(x), float(y), float(z)])-position
            
            pos[0] = min(pos[0], vertices[i][0])
            pos[1] = min(pos[1], vertices[i][1])
            pos[2] = min(pos[2], vertices[i][2])
            
            posMax[0] = max(posMax[0], vertices[i][0])
            posMax[1] = max(posMax[1], vertices[i][1])
            posMax[2] = max(posMax[2], vertices[i][2])
            
        #if not position is None:
        #    pos = np.array(position)
        
        faces = np.zeros((faceCount,3))
        faceVertexes = np.zeros((faceCount,3,3))

        for i in range(faceCount):
            faceLine = f.readline()
            polyCount, a, b, c = faceLine.strip().split(' ')
            faceIndices = [int(a),int(b),int(c)]
            faces[i] = faceIndices
            faceVertexes[i] = np.array([vertices[faceIndices[0]], vertices[faceIndices[1]], vertices[faceIndices[2]]])

    box = posMax-pos
    return pos, box, vertices, faces, faceVertexes


def createVoxelCSV(filename, x, y, z, voxelData):
    with open(filename, 'w') as outfile:
        outfile.write('x,y,z,area\n')
        for i, voxel2D in enumerate(voxelData):
            for j, voxel1D in enumerate(voxel2D):
                for k, voxel in enumerate(voxel1D):
                    #if voxel==0:
                    #    continue
                    #x = pos[0]+steps[0]/2+steps[0]*i
                    #y = pos[1]+steps[1]/2+steps[0]*j
                    #z = pos[2]+steps[2]/2+steps[0]*k
                    outfile.write('{:.4f},{:.4f},{:.4f},{:.6f}\n'.format(x[i], y[j], z[k], voxel))

def createVoxelFile(filename, x, y, z, voxelData, posOff, boxOff, pos, box, steps):
    with open(filename, 'w') as outfile:
        outfile.write('offFileBoxStart = {} {} {}\n'.format(posOff[0], posOff[1], posOff[2]))
        outfile.write('offFileBoxLength  = {} {} {}\n'.format(boxOff[0], boxOff[1], boxOff[2]))
        outfile.write('geometryStart  = {} {} {}\n'.format(pos[0], pos[1], pos[2]))
        outfile.write('geometryLength  = {} {} {}\n'.format(box[0], box[1], box[2]))
        outfile.write('voxelSize = {} {} {}\n'.format(steps[0], steps[1], steps[2]))
        outfile.write('voxelStart = {} {} {}\n'.format(pos[0], pos[1], pos[2]))
        outfile.write('voxelAmount = {} {} {}\n'.format(voxelData.shape[0], voxelData.shape[1], voxelData.shape[2]))
        outfile.write('#x, y, z, area\n')
        for i, voxel2D in enumerate(voxelData):
            for j, voxel1D in enumerate(voxel2D):
                for k, voxel in enumerate(voxel1D):
                    #if voxel==0:
                    #    continue
                    #x = pos[0]+steps[0]/2+steps[0]*i
                    #y = pos[1]+steps[1]/2+steps[0]*j
                    #z = pos[2]+steps[2]/2+steps[0]*k
                    outfile.write('{:.4f} {:.4f} {:.4f} {:.6f}\n'.format(x[i], y[j], z[k], voxel))

def createVTR(filename, x, y, z, voxelData, offset=1):
    from uvw import StructuredGrid, RectilinearGrid, DataArray
    #print(x.shape, y.shape, z.shape, voxelData.shape)
    #grid = StructuredGrid('grid.vtr', np.array((x, y, z)), voxelData.shape, compression=True)
    grid = RectilinearGrid(filename, (x, y, z))#, voxelData.shape, compression=True)
    grid.addCellData(DataArray(voxelData*offset, range(3), 'area'))
    grid.write()

def createParamFile(filename, voxelData, pos, box, steps, boundCond, refinementFraction=0.325):
    with open(filename, 'w') as outfile:
        outfile.write('geometryStart  = {} {} {}\n'.format(pos[0], pos[1], pos[2]))
        outfile.write('geometryLength  = {} {} {}\n'.format(box[0], box[1], box[2]))
        outfile.write('voxelSize = {} {} {}\n'.format(steps[0], steps[1], steps[2]))
        outfile.write('voxelStart = {} {} {}\n'.format(pos[0], pos[1], pos[2]))
        outfile.write('voxelAmount = {} {} {}\n'.format(voxelData.shape[0], voxelData.shape[1], voxelData.shape[2]))
        outfile.write('RefinementFraction = {}\n'.format(refinementFraction))
        outfile.write('XMinBC = {}\n'.format(boundCond['XMin']))
        outfile.write('XMaxBC = {}\n'.format(boundCond['XMax']))
        outfile.write('YMinBC = {}\n'.format(boundCond['YMin']))
        outfile.write('YMaxBC = {}\n'.format(boundCond['YMax']))
        outfile.write('ZMinBC = {}\n'.format(boundCond['ZMin']))
        outfile.write('ZMaxBC = {}\n'.format(boundCond['ZMax']))

def createAreaFile(filename, x, y, z, voxelData, posOff, boxOff, pos, box, steps):
    with open(filename, 'w') as outfile:
        for i, voxel2D in enumerate(voxelData):
            for j, voxel1D in enumerate(voxel2D):
                for k, voxel in enumerate(voxel1D):
                    outfile.write('{:.6f}\n'.format(voxel))

def createMeshDir(outputFolder, voxelData, pos, box, steps):
    # create folder
    logging.info('Create MeshDir Folder')
    meshDirFolder = outputFolder+'/Coarse_meshDir'
    os.makedirs(meshDirFolder, exist_ok=True)
    
    logging.info('Create PRJFile Folder')
    prjFile = meshDirFolder+'/file.prj'
    with open(prjFile, 'w') as outfile:
        outfile.write('Mesh.tri\n')
    
    logging.info('Create Tri file')
    triFile = meshDirFolder+'/Mesh.tri'
    createTRIFile(triFile, voxelData, pos, box, steps)

def createTRIFile(filename, voxelData, pos, box, steps):
    #xCount = math.ceil(box[0]/steps[0])
    #yCount = math.ceil(box[1]/steps[1])
    #zCount = math.ceil(box[2]/steps[2])
    xCount, yCount, zCount = voxelData.shape
    logging.info('writing TRI File of size: {},{},{}'.format(*voxelData.shape))
    #print(xCount, yCount, zCount)

    with open(filename, 'w') as outfile:
        outfile.write('Coarse mesh exported by VekaMesher\n')
        outfile.write('Parametrisierung PARXC, PARYC, TMAXC\n')
        outfile.write(' {} {} '.format(xCount*yCount*zCount, (xCount+1)*(yCount+1)*(zCount+1)))
        outfile.write(' 1 8 12 6     NEL,NVT,NBCT,NVE\n')

        outfile.write(' DCORVG\n')
        for i in range(xCount+1):
            for j in range(yCount+1):
                for k in range(zCount+1):
                    xx = pos[0]+steps[0]*i
                    yy = pos[1]+steps[1]*j
                    zz = pos[2]+steps[2]*k
                    outfile.write('{:.4f} {:.4f} {:.4f} \n'.format(xx, yy, zz))
                    #outfile.write('{} {} {} \n'.format(x, y, z))
        outfile.write(' KVERT\n')
        for i in range(xCount):
            for j in range(yCount):
                for k in range(zCount):
                    i1 = i*(yCount+1)*(zCount+1) + j*(zCount+1) + k + 1 
                    i2 = i1 + 1
                    i3 = i1 + (zCount+1) + 1
                    i4 = i1 + (zCount+1)
                    i5 = i1 + (yCount+1)*(zCount+1)
                    i6 = i5 + 1
                    i7 = i5 + (zCount+1) + 1
                    i8 = i5 + (zCount+1)   
                    outfile.write('{} {} {} {} {} {} {} {} \n'.format(i1, i2, i3, i4, i5, i6, i7, i8))
        outfile.write('KNPR\n')
        for i in range(xCount+1):
            for j in range(yCount+1):
                for k in range(zCount+1):
                    outfile.write('{} \n'.format(0))

def createSection(config, secName):
    config[secName] = {}
    return config[secName]

def insertE3DMaterial(config, name, material):
    print(material)
    secMaterial = createSection(config, 'E3DProcessParameters/{}'.format(name))
    secMaterial['name'] = 'LDPE'
    secMaterial['type'] = 'Polymer'
    
    materialRheo = createSection( config, 'E3DProcessParameters/{}/RheologicalData'.format(name) )
    materialRheo['calcvisco'] = str(material['viscoModel']['model'])
    materialRheo['CalcTemp']  = str(material['tempModel']['model'])
    
    if material['tempModel']['model'] == 'isotherm':
        secTempModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/Isotherm'.format(name))
    elif material['tempModel']['model'] == 'etb':
        secTempModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/etb'.format(name))
        secTempModel['ActivatingEnergy']     = str(material['tempModel']['activatingEnergy'])
        secTempModel['ReferenceTemperature'] = str(material['tempModel']['referenceTemperature'])
    elif material['tempModel']['model'] == 'tbts':
        secTempModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/TbTs'.format(name))
        secTempModel['StandardTemperature']  = str(material['tempModel']['standardTemperature'])
        secTempModel['ReferenceTemperature'] = str(material['tempModel']['referenceTemperature'])
    elif material['tempModel']['model'] == 'c1c2':
        secTempModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/C1C2'.format(name))
        secTempModel['C1'] = str(material['tempModel']['c1'])
        secTempModel['C2'] = str(material['tempModel']['c2'])
    else:
        raise ValueError('Unknown Viscosity model: {}'.format(material['tempModel']['model']))


    if material['viscoModel']['model'] == 'Carreau':
        materialViscModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/Carreau'.format(name))
        materialViscModel['zeroviscosity'] = str(material['viscoModel']['zeroViscosity'])
        materialViscModel['recipvelocity'] = str(material['viscoModel']['recipVelocity'])
        materialViscModel['exponent'] = str(material['viscoModel']['exponent'])
    elif material['viscoModel']['model'] == 'Power law':
        materialViscModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/Power law'.format(name))
        materialViscModel['Consistence'] = str(material['viscoModel']['consistence'])
        materialViscModel['Exponent']    = str(material['viscoModel']['exponent'])
    elif material['viscoModel']['model'] == 'Ellis':
        materialViscModel = createSection( config, 'E3DProcessParameters/{}/RheologicalData/Ellis'.format(name))
        materialViscModel['ZeroViscosity'] = str(material['viscoModel']['zeroViscosity'])
        materialViscModel['Gamma0']        = str(material['viscoModel']['gamma0'])
        materialViscModel['Exponent']      = str(material['viscoModel']['exponent'])
    else:
        raise ValueError('Unknown Viscosity model: {}'.format(material['viscoModel']['model']))
    
    secThermo = createSection( config, 'E3DProcessParameters/{}/ThermoData'.format(name))
    secThermo['DensityModel'] = 'DENSITY' # FIXED
    secThermo['HeatConductivity']      = '0.231'   #material['heatConductivity']
    secThermo['HeatConductivitySlope'] = '0.000017'#material['heatConductivitySlope']
    secThermo['HeatCapacity']          = '2.087'   #material['heatCapacity']
    secThermo['HeatCapacitySlope']     = '0.00431' #material['heatCapacitySlope']
    
    secThermoDensity = createSection( config, 'E3DProcessParameters/{}/ThermoData/Density'.format(name))
    secThermoDensity['density'] = '1.35'
    secThermoDensity['densityslope'] = '0.0001'
    

def createE3DFile( filename:str, infoData:dict ):
    import configparser
    config = configparser.ConfigParser()
    config.optionxform=str
    
    cadData = infoData['cad']
    inflows = cadData['inflows']
    digitalTwin = infoData['simod']['digital_twin']
    
    machine = createSection( config, 'E3DGeometryData/Machine')
    # fixed
    machine['Unit'] = 'mm'
    machine['type'] = 'die'
    machine['ScrewCylinderRendering'] = 'OFF'
    # variable
    print(math.sqrt(cadData['box'][0]**2+cadData['box'][1]**2))
    machine['barreldiameter'] = '{}'.format(math.sqrt(cadData['box'][0]**2+cadData['box'][1]**2))
    machine['noofelements']   = '1'
    machine['barrellength']   = str(cadData['box'][2])
    # @todo remove?
    machine['innerdiameter']     = '100.0'
    machine['rotationdirection'] = 'LEFT'

    element = createSection(config, 'E3DGeometryData/Machine/Element_1')
    element['Unit'] = 'mm'
    # fixed
    element['startposition'] = '0.0'
    element['type']          = 'OFF'
    element['objecttype']    = 'die'
    # variable
    element['name']         = digitalTwin['name']
    element['diameter']     = str(math.sqrt(cadData['box'][0]**2+cadData['box'][1]**2))
    element['off_filelist'] = 'surface.off'
    # remove?
    element['calculatedinflowdiammin'] = '100.0'
    element['calculatedinflowdiammax'] = '122.5'
    element['innerdiameter'] = '100.0'
    
    process = createSection(config, 'E3DProcessParameters')
    # fixed
    process['Unit'] = 'mm'
    process['processtype'] = 'THROUGHPUT'
    # variable
    process['nOfInflows']     = str(len(inflows))
    process['massthroughput'] = infoData['simod']['input_parameter']['massthroughput_global']['value']
    process['materialtemperature'] = infoData['simod']['input_parameter']['temperature']['value']
    process['barreltemperature']   = infoData['simod']['input_parameter']['temperature']['value'] # @TODO check if correct
    process['barreltemperatureadiabatic'] = 'NO'
    #process['ExtrusionSpeed_CMpS'] = 'XXX' # @TODO remove?
    process['ExtrusionGapSize_MM'] = str(infoData['cad']['extrusion']['minGap'])
    process['ExtrusionOutflowArea_MM2'] = str(infoData['cad']['extrusion']['area'])
    # remove?
    process['screwspeed']        = '430.0'
    process['screwtemperature']  = '_INVALID_'
    process['screwtemperatureadiabatic'] = 'YES'
    process['mininflowdiameter'] = '100.0'
    process['maxinflowdiameter'] = '122.5'
    
    for i, inflow in enumerate(inflows):
        inflowSimodData = infoData['simod']['input_parameter']['inflow']['value'][i]
        inflowSec = createSection( config, 'E3DProcessParameters/Inflow_{}'.format(i+1) )
        inflowSec['Unit'] = 'mm'
        # fixed
        inflowSec['material'] = '1'
        # variable
        if inflow['type'] == 'circle':
            inflowSec['type'] = 'CURVEDFLAT'
            inflowSec['InnerRadius']  = str(float(inflow['radius'])*0.8) # @todo what value?
            inflowSec['OuterRadius']  = str(inflow['radius'])
            inflowSec['massflowrate'] = '68' # @todo fix str(inflowSimodData['masstroughput'])
            inflowSec['center'] = '{},{},{}'.format(*inflow['position']) # @todo make sure it is on the bbox!!
            inflowSec['normal'] = '{},{},{}'.format(*inflow['flow_normal']) # in direction
        elif inflow['type'] == 'ring':
            inflowSec['type'] = 'CURVEDFLAT'
            inflowSec['InnerRadius']  = inflow['inner_radius']
            inflowSec['OuterRadius']  = inflow['radius']
            inflowSec['massflowrate'] = str(68.0) # @todo handle
            inflowSec['center'] = '{},{},{}'.format(*inflow['position']) # @todo make sure it is on the bbox!!
            inflowSec['normal'] = '{},{},{}'.format(*inflow['flow_normal']) # in direction
        else:
            raise ValueError('BAD')
    
    material = infoData['simod']['materials'][0]
    insertE3DMaterial(config, 'Material', material)
    
    output = createSection( config, 'E3DSimulationsettings/Output' )
    # fixed
    output['nOf1DLayers']       = '128'
    output['nOfHistogramBins']  = '16'
    output['HistogramShearMax'] = '1e4'
    output['HistogramShearMin'] = '1e-1'
    output['CutData_1D']        = '0.001' # @TODO check if correct
    
    simsetting = createSection( config, 'E3DSimulationsettings' )
    simsetting['KTPRelease'] = 'NO'
    simsetting['AutomaticTimeStepControl']  = 'yes' # CHANGE?
    simsetting['TimeStepEnlargementFactor'] = '6e0' # CHANGE?
    simsetting['PressureFBM'] = 'on'
    
    with open(filename, 'w') as configfile:
        config.write(configfile)


def main():
    parser = argparse.ArgumentParser(description='Creates a JSON file containing an XYZ-sorted three dimensional array. Each floating point value in the array corresponds with the exact area of surface present in the respective box in mm^2. \n Box dimensions in mm are to be set by the -s option. \n Example: python3 computeAreas.py data.off someresultfolder(unused so far) -s 100 100 100. This generates a mesh with 100mm x 100mm x 100mm dimensions.')
    parser.add_argument('inputFolder', metavar='inputFolder', help="The simod json file")
    
    parser.add_argument('-o,', '--output', action="store", dest='output', default=None, type=str, help='Sets the outputfolder independend of the input folder')
    parser.add_argument('--offFile',  dest='offFile',  action='store', default=None, type=str, help="The off file of the die")
    parser.add_argument('--infoData', dest='infoData', action='store', default=None, type=str, help="The json file containing several infos about the original file.")
    
    parser.add_argument('-s', '--steps', nargs='+', dest="steps", default=None, type=float, help="Set the desired step-sizes.")
    parser.add_argument('--linear', dest="linear", default=10, type=float, help="Set the desired step-sizes by defining the linear factor to the minimal gap.")
    
    parser.add_argument('--start', dest="startIndex", default=0, type=int, help="Start with this face")
    parser.add_argument('--end', dest="endIndex", default=None, type=int, help="End with this face")
    
    parser.add_argument('--ref', dest='refinementFactor', type=float, default=0.325, help="Set the desired refinment factor")
    parser.add_argument('--vtr', dest="vtr", action="store_true", help="Set if we want to write the vtr file")
    
    parser.add_argument('--debug', action="store_true", dest="debug", help="Activate debug messages")
    
    parser.add_argument('--bbXExt', action="store", dest='boundboxXExtension', type=float, default=5, help="The percentage amount to extend the bounding box x direction.")
    parser.add_argument('--bbYExt', action="store", dest='boundboxYExtension', type=float, default=5, help="The percentage amount to extend the bounding box x direction.")
    
    
    args = parser.parse_args()
    print(args)
    if args.debug:
        logging.basicConfig(format='%(asctime)s - %(levelname) 8s - %(message)s', level=logging.DEBUG, datefmt='%H:%M:%S')
    else:
        logging.basicConfig(format='%(asctime)s - %(levelname) 8s - %(message)s', level=logging.INFO, datefmt='%H:%M:%S')

    outputFolder = args.inputFolder
    if args.output:
        outputFolder = args.output
       
    offFilepath  = args.inputFolder+'/surface.off' if not args.offFile else args.offFile
    infoFilepath = args.inputFolder+'/info.json' if not args.infoData else args.infoData
    logging.info("Input off: {}".format(offFilepath))
    logging.info("Input info: {}".format(infoFilepath))
    logging.info("Output: {}".format(outputFolder))
    surfaceFilename    = outputFolder+'/surface.off'

    logging.info("Load infoFile")
    with open(infoFilepath, 'r') as infoFile:
        infoData = json.load(infoFile)

    print(infoData['cad'])
    posCAD = np.array(infoData['cad']['position'])
    boxCAD = np.array(infoData['cad']['box'])
    
    os.makedirs(outputFolder, exist_ok=True)
    

    pos = np.array(posCAD, dtype=np.float32)
    box = np.array(boxCAD, dtype=np.float32)
    logging.info("Enlarge boundbox: {}%, {}%, {}%".format(args.boundboxXExtension, args.boundboxYExtension, 0))
    logging.info("Before")
    print(pos, box)
    print(box[0]*(args.boundboxXExtension/100))
    
    pos[0] = pos[0]-box[0]*(args.boundboxXExtension/100)
    pos[1] -= box[1]*(args.boundboxXExtension/100)
    box[0] += box[0]*(args.boundboxYExtension/100)*2
    box[1] += box[1]*(args.boundboxYExtension/100)*2

    #load faces
    print(pos)
    logging.info("Load Off file")
    posOff, boxOff, _, _, faces = loadOffFile(offFilepath)
    


    logging.info("Resulting dimensions")
    print(pos, box)
    
    logging.info('Create Voxel data')
    if not args.steps:
        steps = [args.linear*infoData['cad']['extrusion']['minGap'], args.linear*infoData['cad']['extrusion']['minGap'], args.linear*infoData['cad']['extrusion']['minGap']]
    else:
        steps = args.steps
    logging.info('Box: {}, {}, {}'.format(box[0],box[1],box[2]))
    logging.info('Given steps: {}, {}, {}'.format(steps[0],steps[1],steps[2]))
    xCount = math.ceil(box[0]/steps[0])
    yCount = math.ceil(box[1]/steps[1])
    zCount = math.ceil(box[2]/steps[2])
    logging.info('{},{},{}'.format(xCount, yCount, zCount))
    steps[0] = box[0]/xCount
    steps[1] = box[1]/yCount
    steps[2] = box[2]/zCount
    logging.info('Used steps: {}, {}, {}'.format(steps[0],steps[1],steps[2]))
    logging.info('{},{},{}'.format(steps, steps[0]*xCount, steps[1]*yCount, steps[2]*zCount))
    voxelData = np.zeros((xCount, yCount, zCount))
    
    print(voxelData.shape)
    
    logging.info('Map inflows onto the simbox')
    print(infoData['cad']['inflows'])
    for i in range(len(infoData['cad']['inflows'])):
        inflow = infoData['cad']['inflows'][i]
        fnormal = np.array(inflow['face_normal'])
        if sum(abs(fnormal)) != 1:
            raise ValueError('NIMPL_YET, currently only supporting cardinal inflows')
        if fnormal[0] < -0.1:
            inflow['position'][0] = pos[0]
        elif fnormal[0] > 0.1:
            inflow['position'][0] = pos[0]+box[0]
        elif fnormal[1] < -0.1:
            inflow['position'][1] = pos[1]
        elif fnormal[1] > 0.1:
            inflow['position'][1] = pos[1]+box[1]
        elif fnormal[2] < -0.1:
            inflow['position'][2] = pos[2]
        elif fnormal[2] > 0.1:
            inflow['position'][2] = pos[2]+box[2]
        else:
            raise ValueError('Something is not right')
        infoData['cad']['inflows'][i] = inflow
    print(infoData['cad']['inflows'])
    createE3DFile(outputFolder+'/setup.e3d', infoData)
    
    
    startVoxel = pos
    print("start:", startVoxel)
    print(np.floor(pos/steps).astype(int))
    corIndexes = [-1,-1,0]
    logging.info("Process faces")
    #process faces
    for i, f in enumerate(faces[args.startIndex:]):
        #if i<args.startIndex:
        #    continue
        if i%1000==0:
            logging.info('{}/{}'.format(i+1+args.startIndex, len(faces)))
        testFaceVoxel(f, startVoxel, steps, voxelData)
        #if i==2000:
        #    break
    
    
    
    
    #print(voxelData)
    logging.info("Write data")
    data = {
        'pos'   : pos.tolist(),
        'box'   : box.tolist(),
        'steps' : steps,
        'voxels': voxelData.tolist(),
    }

    with open(outputFolder+'/area.json', 'w') as outfile:
        json.dump(data, outfile, indent = 4)
    
    x = np.linspace(pos[0]+steps[0]/2, pos[0]+steps[0]/2+steps[0]*xCount, xCount)
    y = np.linspace(pos[1]+steps[1]/2, pos[1]+steps[1]/2+steps[1]*yCount, yCount)
    z = np.linspace(pos[2]+steps[2]/2, pos[2]+steps[2]/2+steps[2]*zCount, zCount)
    
    xVert = np.linspace(pos[0], pos[0]+steps[0]*(xCount+1), xCount+1, endpoint=False)
    yVert = np.linspace(pos[1], pos[1]+steps[1]*(yCount+1), yCount+1, endpoint=False)
    zVert = np.linspace(pos[2], pos[2]+steps[2]*(zCount+1), zCount+1, endpoint=False)

    #print(x)
    #print(y)
    #print(z)
    
    createVoxelFile(outputFolder+'/voxel.txt', x,y,z, voxelData, posOff, boxOff, pos, box, steps)
    createVoxelCSV(outputFolder+'/voxel.csv', x,y,z, voxelData)
    
    

    boundCond = {
        'XMin': 'Wall',
        'XMax': 'Wall',
        'YMin': 'Wall',
        'YMax': 'Wall',
        'ZMin': 'Wall',
        'ZMax': 'Outflow',
    }
    
    inflows = infoData['cad']['inflows']
    print(pos)
    print(box)
    for i, inflow in enumerate(inflows):
        print(inflow)
        for j in range(3):
            print(inflow['position'][j] - pos[j])
            print(inflow['position'][j] - (pos[j]+box[j]))
        if abs(inflow['position'][0] - pos[0]) < box[0]*0.1: # test XMin
            if boundCond['XMin'] == 'Wall':
                boundCond['XMin'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side XMin already has an inflow')
        elif abs(inflow['position'][0] - (pos[0]+box[0])) < box[0]*0.1: # test XMax
            if boundCond['XMax'] == 'Wall':
                boundCond['XMax'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side XMax already has an inflow')
        elif abs(inflow['position'][1] - pos[1]) < box[1]*0.1: # test YMin
            if boundCond['YMin'] == 'Wall':
                boundCond['YMin'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side YMin already has an inflow')
        elif abs(inflow['position'][1] - (pos[1]+box[1])) < box[1]*0.1: # test XMax
            if boundCond['YMax'] == 'Wall':
                boundCond['YMax'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side YMax already has an inflow')
        elif abs(inflow['position'][2] - pos[2]) < box[2]*0.1: # test YMin
            if boundCond['ZMin'] == 'Wall':
                boundCond['ZMin'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side ZMin already has an inflow')
        elif abs(inflow['position'][2] - (pos[2]+box[2])) < box[2]*0.1: # test XMax
            if boundCond['ZMax'] == 'Wall':
                boundCond['ZMax'] = 'Inflow-{}'.format(i+1)
            else:
                raise ValueError('Side ZMax already has an inflow')
        else:
            raise ValueError('Could not determine simbox side for inflow: {}'.format(i))
    
    createParamFile(outputFolder+'/param.txt', voxelData, pos, box, steps, boundCond, args.refinementFactor)
    createAreaFile(outputFolder+'/area.txt', x,y,z, voxelData, posOff, boxOff, pos, box, steps)
    createMeshDir(outputFolder, voxelData, pos, box, steps)
    
    #createTRIFile(outputFolder+'/meshDir/Mesh.tri', x,y,z, voxelData, posOff, boxOff, pos, box, steps)
    #createPRJFile(outputFolder+'/meshDir/file.prj', x,y,z, voxelData, posOff, boxOff, pos, box, steps)
    
    if args.vtr:
        createVTR(outputFolder+'/voxel.vtr', xVert, yVert, zVert, voxelData)


# execute the script and catch the exception
# to allow for better traceback then what freecad provides
#import sys, traceback
#try:
main()
#except:
#   print('-'*60)
#  traceback.print_exc(file=sys.stdout)
# print('-'*60)
#exit(1)
