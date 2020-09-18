#!/usr/bin/env bash


grep ghput _data/prot.0000.txt  | grep NaN > flow
gnuplot plot_conv_M+S.gp
rm -f flow

eog convergence_M+S_flowrate.png


exit 0
