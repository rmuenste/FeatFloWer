#!/usr/bin/env bash

cp q2p1_sse_partfil q2p1_sse

cp _data_BU/q2p1_paramV_0.dat _data_BU/q2p1_paramV_BU.dat
python3 ./e3d_start.py -n 16 -f _ianus/TSE/ConvPF -a 0

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np 16 q1_scalar_partfil 

cp _data_BU/q2p1_paramV_1.dat _data_BU/q2p1_paramV_BU.dat
python3 ./e3d_start.py -n 16 -f _ianus/TSE/ConvPF -a 0

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np 16 q1_scalar_partfil 

cp _data_BU/q2p1_paramV_1.dat _data_BU/q2p1_paramV_BU.dat
python3 ./e3d_start.py -n 16 -f _ianus/TSE/ConvPF -a 0

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np 16 q1_scalar_partfil 

cp _data_BU/q2p1_paramV_1.dat _data_BU/q2p1_paramV_BU.dat
python3 ./e3d_start.py -n 16 -f _ianus/TSE/ConvPF -a 0

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np 16 q1_scalar_partfil 

cp _data_BU/q2p1_paramV_1.dat _data_BU/q2p1_paramV_BU.dat
python3 ./e3d_start.py -n 16 -f _ianus/TSE/ConvPF -a 0

exit 0
