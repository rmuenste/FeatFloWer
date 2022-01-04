cp _data_BU/q2p1_paramV_DIE_0.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n 16 -f _ianus/DIE/RH_LOC --die-simulation

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np 16 q1_scalar 

cp _data_BU/q2p1_paramV_DIE_1.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n 16 -f _ianus/DIE/RH_LOC --die-simulation

