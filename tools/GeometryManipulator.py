from copy import deepcopy
import argparse
import os
import math

class GeometryObject(object):

  def __init__(self):
    self.NoOfVertices = None
    self.NoOfElements = None
    self.vertices = []
    self.elements = []
    self.filePresent = False

  def checkStatus(self):
    if self.filePresent == True:
      raise Exception("Error: Already loaded a File")

  def DeleteGeo(self):
    del self.vertices[:]
    del self.elements[:]
    self.NoOfVertices = None
    self.NoOfElements = None
    self.filePresent = False

#-----------------------------------------------------
# Importers
#-----------------------------------------------------

  def ReadSalomeDat(self,filename):
    self.checkStatus()
    filename = os.path.abspath(filename)
    with open(filename,'r') as f:
      firstline = f.readline()
      firstline_splitted = firstline.split()
      noOfVtx = int(firstline_splitted[0])
      self.NoOfVertices = noOfVtx
      noOfElements = int(firstline_splitted[1])
      for i in range(noOfVtx):
        line = f.readline()
        line_splitted = line.split()
        ctrlvtxnr = int(line_splitted[0])
        vtx_x = float(line_splitted[1])
        vtx_y = float(line_splitted[2])
        vtx_z = float(line_splitted[3])
        self.vertices.append([vtx_x,vtx_y,vtx_z])
      if ctrlvtxnr != noOfVtx:
        print "Warning! File might be corrupt. NoOfVtx not correct"
      for i in range(noOfElements):
        line = f.readline()
        line_splitted = line.split()
        ctrlidxnr = int(line_splitted[0])
        eltype = int(line_splitted[1])
        # Pass if it is an edge
        if eltype == 102:
          pass
        else:
          el = []
          for s in line_splitted[2::]:
            el.append(int(s)-1)
          self.elements.append(deepcopy(el))
      if ctrlidxnr != noOfElements:
        print "Warning! File might be corrupt! NoOfElements not correct!"
    self.NoOfElements = len(self.elements)
    self.filePresent = True

  def ReadOff(self,filename):
    self.checkStatus()
    filename = os.path.abspath(filename)
    with open(filename,'r') as f:
      firstline = f.readline() # Contains 'OFF'
      if firstline.strip().lower() !='off':
        raise Exception("This is not an OFF-File. Missing the magic keyword")
      secondline = f.readline()
      secondline_splitted = secondline.split()
      self.NoOfVertices = int(secondline_splitted[0])
      self.NoOfElements = int(secondline_splitted[1])
      for i in range(self.NoOfVertices):
        line = f.readline()
        line_splitted = line.split()
        vtx_x = float(line_splitted[0])
        vtx_y = float(line_splitted[1])
        vtx_z = float(line_splitted[2])
        self.vertices.append([vtx_x,vtx_y,vtx_z])
      for i in range(self.NoOfElements):
        line = f.readline()
        line_splitted = line.split()
        vtxInEl = int(line_splitted[0])
        elDef = []
        for k in range(vtxInEl):
          elDef.append(int(line_splitted[k+1]))
        self.elements.append(deepcopy(elDef))
    self.filePresent = True


  def ReadPly(self,filename):
    self.checkStatus()
    filename = os.path.abspath(filename)
    with open(filename,'r') as f:
      firstline = f.readline() # Contains 'ply'
      if firstline.strip().lower() !='ply':
        raise Exception("This is not a PLY-File. Missing the magic keyword")
      line = f.readline()
      while line.strip().lower() != 'end_header':
        line_splitted = line.split()
        if 'element' in line.lower() and 'vertex' in line.lower():
          self.NoOfVertices = int(line_splitted[2])
        if 'element' in line.lower() and 'face' in line.lower():
          self.NoOfElements = int(line_splitted[2])
        line = f.readline()
      # Now line = 'end_header'
      for i in range(self.NoOfVertices):
        line = f.readline()
        line_splitted = line.split()
        vtx_x = float(line_splitted[0])
        vtx_y = float(line_splitted[1])
        vtx_z = float(line_splitted[2])
        self.vertices.append([vtx_x,vtx_y,vtx_z])
      for i in range(self.NoOfElements):
        line = f.readline()
        line_splitted = line.split()
        vtxInEl = int(line_splitted[0])
        elDef = []
        for k in range(vtxInEl):
          elDef.append(int(line_splitted[k+1]))
        self.elements.append(deepcopy(elDef))
    self.filePresent = True
#------------------------------------------------------------------------
# Manipulation routines
#------------------------------------------------------------------------

  def TranslateObject(self,tx,ty,tz):
    tx = float(tx)
    ty = float(ty)
    tz = float(tz)
    vertices_new = []
    vertices_old = self.vertices
    for vtx in self.vertices:
      vtx_new = [vtx[0]+tx, vtx[1]+ty, vtx[2]+tz]
      vertices_new.append(deepcopy(vtx_new))
    self.vertices = vertices_new
    del vertices_old

  def ScaleObject(self,sx,sy,sz):
    sx = float(sx)
    sy = float(sy)
    sz = float(sz)
    vertices_new = []
    vertices_old = self.vertices
    for vtx in self.vertices:
      vtx_new = [vtx[0]*sx, vtx[1]*sy, vtx[2]*sz]
      vertices_new.append(deepcopy(vtx_new))
    self.vertices = vertices_new
    del vertices_old

  def RotateX(self,angle, rad=False):
    angle = float(angle)
    if rad == False:
      angle = math.radians(angle)
    
    vertices_new = []
    vertices_old = self.vertices
    s = math.sin(angle)
    c = math.cos(angle)

    for vtx in self.vertices:
      x = vtx[0]
      y = vtx[1]*c - vtx[2]*s
      z = vtx[1]*s + vtx[2]*c
      vtx_new = [x, y, z]
      vertices_new.append(deepcopy(vtx_new))
    self.vertices = vertices_new
    del vertices_old

  def RotateY(self,angle, rad=False):
    angle = float(angle)
    if rad == False:
      angle = math.radians(angle)
    
    vertices_new = []
    vertices_old = self.vertices
    s = math.sin(angle)
    c = math.cos(angle)

    for vtx in self.vertices:
      x = vtx[0]*c + s*vtx[2]
      y = vtx[1]
      z = (-1.0)*s*vtx[0] + c*vtx[2]
      vtx_new = [x, y, z]
      vertices_new.append(deepcopy(vtx_new))
    self.vertices = vertices_new
    del vertices_old

  def RotateZ(self,angle, rad=False):
    angle = float(angle)
    if rad == False:
      angle = math.radians(angle)
    
    vertices_new = []
    vertices_old = self.vertices
    s = math.sin(angle)
    c = math.cos(angle)

    for vtx in self.vertices:
      x = vtx[0]*c - vtx[1]*s
      y = vtx[0]*s + vtx[1]*c
      z = vtx[2]
      vtx_new = [x, y, z]
      vertices_new.append(deepcopy(vtx_new))
    self.vertices = vertices_new
    del vertices_old

#-------------------------------------------------------------------------
# Exporters in different format
#-------------------------------------------------------------------------

  def ExportOff(self,filename):
    filename = os.path.abspath(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
      os.path.makedirs(dirname)
    with open(filename,'w') as f:
      f.write("OFF\n")
      f.write("{} {} {}\n".format(self.NoOfVertices,self.NoOfElements,0))
      for i in range(self.NoOfVertices):
        vtx = self.vertices[i]
        #f.write('{} {} {}\n'.format(vtx[0],vtx[1],vtx[2]))
        f.write('{x:{digitsx}f} {y:{digitsy}f} {z:{digitsz}f}\n'.format(x=vtx[0],y=vtx[1],z=vtx[2], digitsx=16,digitsy=16,digitsz=16))
      for i in range(self.NoOfElements):
        el = self.elements[i]
        f.write('{:d} '.format(len(el)))
        for k in el:
          f.write('{:d} '.format(k))
        f.write('\n')
  
  def ExportPly(self,filename):
    filename = os.path.abspath(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
      os.path.makedirs(dirname)
    with open(filename,'w') as f:
      f.write("ply\n")
      f.write("format ascii 1.0\n")
      f.write("comment Exported by Maltes Geometry Manipulator\n")
      f.write("element vertex {}\n".format(self.NoOfVertices))
      f.write("property float x\nproperty float y\nproperty float z\n")
      f.write("element face {}\n".format(self.NoOfElements))
      f.write('property list uchar int vertex_indices\n')
      f.write('end_header\n')
      for i in range(self.NoOfVertices):
        vtx = self.vertices[i]
        #f.write('{} {} {}\n'.format(vtx[0],vtx[1],vtx[2]))
        f.write('{x:{digitsx}f} {y:{digitsy}f} {z:{digitsz}f}\n'.format(x=vtx[0],y=vtx[1],z=vtx[2], digitsx=16,digitsy=16,digitsz=16))
      for i in range(self.NoOfElements):
        el = self.elements[i]
        f.write('{:d} '.format(len(el)))
        for k in el:
          f.write('{:d} '.format(k))
        f.write('\n')

  def ExportStl(self,filename):
    filename = os.path.abspath(filename)
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
      os.path.makedirs(dirname)
    with open(filename,'w') as f:
      f.write('solid geometry\n')
      for i in range(self.NoOfElements):
        el = self.elements[i]
        p1 = self.vertices[el[0]]
        p2 = self.vertices[el[1]]
        p3 = self.vertices[el[2]]
        elpoints = []
        for k in el:
          elpoints.append(self.vertices[k])
        v = [p1[0] - p2[0], p1[1]-p2[1], p1[2] - p2[2]]
        w = [p2[0] - p3[0], p2[1]-p3[1], p2[2] - p3[2]]
        N = [v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]]
        NormN = N[0]**2 + N[1]**2 + N[2]**2
        f.write('  facet normal {} {} {}\n'.format(N[0]/NormN, N[1]/NormN, N[2]/NormN))
        f.write('    outer loop\n')
        for p in elpoints:
          f.write('      vertex {} {} {}\n'.format(p[0],p[1],p[2]))
        f.write('    endloop\n')
        f.write('  endfacet\n')
      f.write('endsolid geometry\n')