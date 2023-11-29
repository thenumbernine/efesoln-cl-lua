#!/bin/sh
../cl-cpu/run.lua -kernelCallMethod C-singlethread -I include -I ../../cpp/Common/include -I ../../cpp/Tensor/include ./run.lua "$@"
