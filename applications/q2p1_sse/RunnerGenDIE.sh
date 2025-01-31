n=128
f=11598

cp _data_BU/q2p1_paramV_DIE_0.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n $n -f $f --die-simulation

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np $n q1_scalar_multimat

cp _data_BU/q2p1_paramV_DIE_1.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n $n -f $f --die-simulation


