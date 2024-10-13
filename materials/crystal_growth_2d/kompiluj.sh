#!/bin/bash
g++ -Ofast -fdefault-real-8  -fdefault-double-8  -w \
            -ffixed-form -ffixed-line-length-none -ffpe-summary=none \
           -fopenmp  -fvect-cost-model=unlimited  -ftree-vectorize  -fpic \
           -fvect-cost-model=unlimited \
           -fprefetch-loop-arrays \
           -funroll-loops \
           -march=native \
           -flto \
           *.cpp \
            -I${HOME}/bin/bib \
            -I${MKLROOT}/include_gf \
             -Wl,--start-group \
                  ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
                  ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
                  ${MKLROOT}/lib/intel64/libmkl_core.so \
             -Wl,--end-group \
             -lgomp -lpthread -lm -ldl         

