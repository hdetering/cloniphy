#!/bin/bash

# load required modules
module load gcc/6.3.0 
module load cmake/3.8.1
module load python/2.7.12

# build project (incl. debug info)
cmake -DCMAKE_BUILD_TYPE=Debug .. && make
