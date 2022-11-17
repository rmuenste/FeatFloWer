
import argparse
from math import dist

# important table for face indices
faceIDs=(
    (3,2,1,0),
    (0,1,5,4),
    (1,2,6,5),
    (2,3,7,6),
    (3,0,4,7),
    (4,5,6,7)
)

# subroutines that do the work :-)
def readTRI3D(filename):
    nodes, cells=None, None
    with open(filename,"r") as f:
        # skip the first two lines
        f.readline()
        f.readline()
        # read in the line with the mesh data
        s=f.readline().strip().split()
        nel,nvt,nbct,nve,nee,nae=map(int,s[:6])
        assert (nve==8 and nee==12 and nae==6), "This program only works with hexaeder meshes!"
        # check for nodes block
        s=f.readline().strip()
        assert s=="DCORVG", "DCORVG block was not found!"
        # read in the nodes
        nodes=tuple(tuple(map(float,f.readline().strip().split())) for i in range(nvt))
        # check for cells block
        s=f.readline().strip()
        assert s=="KVERT", "KVERT block was not found!"
        # read in the cells
        l=lambda x: int(x)-1
        cells=tuple(tuple(map(l,f.readline().strip().split())) for i in range(nel))
        # the rest will be ignored
    return (nodes,cells)

def find_boundary_node_IDs(mesh):
    eps=1e+6
    points,cells=mesh
    bnIDs=set()
    # create list of connectors over edges
    connectors=[tuple(sorted([c[faceIDs[i][j]] for j in range(4)])) for c in cells for i in range(6)]
    connectors.sort()
    # parse connector list
    a=connectors.pop()
    while connectors:
        b=connectors.pop()
        if a==b:
            if connectors:
                a=connectors.pop()
            else:
                break
        else:
            eps=min(eps,*[dist(points[a[i]],points[a[(i+1)%4]]) for i in range(4)])
            for ID in a:
               bnIDs.add(ID)
            a=b
    return (bnIDs,0.2*eps)

def write_tri3d(filename,points,cells):
    def point2str(point):
        return " ".join(f'{i}' for i in point)
    def cell2str(cell):
        return " ".join(str(i+1) for i in cell)
    with open(filename,"w") as f:
        f.write("Mesh created by fusing two meshs together\n")
        f.write("Parametrisierung PARXC, PARYC, TMAXC\n")
        f.write("%d %d 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE\n"%(len(cells),len(points)))
        f.write("DCORVG\n")
        for point in points:
            f.write(point2str(point)+"\n")
        f.write("KVERT\n")
        for cell in cells:
            f.write(cell2str(cell)+"\n")
        f.write("KNPR\n")
        f.write("0\n"*len(points))

def fuse_meshes(m1,m2,bnd1,bnd2,eps):
    points=list(m1[0])
    cells=list(m1[1])
    idx=len(points)
    mapping=dict() # map node numbers of second mesh to new numbers
    bnd1=bnd1.copy()
    n=len(m2[0])
    for (i,p) in enumerate(m2[0]):
        if i%1000==0:
           print("%6.2f%%"%(100/n*i))
        if i in bnd2:
            # try to map this point to a boundary node in the first mesh
            tmp=-1
            for j in bnd1:
                if dist(p,points[j])<eps:
                    tmp=j
                    break
            if tmp>=0: # map the node to node on old mesh
                mapping[i]=tmp
                bnd1.remove(tmp)
                continue
        mapping[i]=idx
        points.append(p)
        idx=idx+1
    print("100.00%\nDone!")
    # replace the node numbers in each cell with the new node numbers
    cells.extend(tuple(mapping[i] for i in c) for c in m2[1])
    return (tuple(points),tuple(cells))
    
# main program
parser=argparse.ArgumentParser(
    prog="Mesh Fuser",
    description="This program reads two meshes in TRI3D format and tries to fuse then at their boundaries."
)
parser.add_argument("mesh_1", help="First mesh")
parser.add_argument("mesh_2", help="Second mesh")
parser.add_argument("mesh_out", help="Name of the resulting mesh")
parser.add_argument("--eps", help="Specify distance for fusing of nodes. Will be computed by a heuristic if not given.")
args=parser.parse_args()
mesh_1=readTRI3D(args.mesh_1)
mesh_2=readTRI3D(args.mesh_2)
bnIDs_1,eps_1=find_boundary_node_IDs(mesh_1)
bnIDs_2,eps_2=find_boundary_node_IDs(mesh_2)
eps=min(eps_1,eps_2)
if args.eps!=None:
   eps=float(args.eps)
   assert eps>0.0, "Parameter 'eps' has to be larger than 0!"
mesh_out=fuse_meshes(mesh_1,mesh_2,bnIDs_1,bnIDs_2,eps)
write_tri3d(args.mesh_out,*mesh_out)
