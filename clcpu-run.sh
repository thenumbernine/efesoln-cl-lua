#!/bin/sh
# use singlethread for debugging within kernel loops
#../cl-cpu/run.lua -kernelCallMethod C-singlethread -I include -I ../../cpp/Common/include -I ../../cpp/Tensor/include ./run.lua "$@"
../cl-cpu/run.lua -kernelCallMethod C-multithread -I include -I ../../cpp/Common/include -I ../../cpp/Tensor/include ./run.lua "$@"
