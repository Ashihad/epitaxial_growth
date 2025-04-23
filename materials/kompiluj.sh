#!/bin/bash
g++ -Ofast -fdefault-real-8  -fdefault-double-8  -w \
            -ffixed-form -ffixed-line-length-none -ffpe-summary=none \
           -fopenmp  -fvect-cost-model=unlimited  -ftree-vectorize  -fpic \
           -fvect-cost-model=unlimited \
           -fprefetch-loop-arrays \
           -funroll-loops \
           -march=native \
           -flto \
	   -ggdb3 \
           *.cpp \
            -lpthread -lm -ldl -o release.out        

