function(append_performance_compiler_flags)
    # -Ofast                      - Allow the most aggressive performance optimizations
    # -fdefault-real-8            - Force floats to be 8 bytes long, Fortran specific...?
    # -fdefault-double-8          - Force doubles to be 8 bytes long, doubles are always 8 bytes so obsolete...?
    # -ffixed-form                - Forces the use of fixed-form source layout (as used in older versions of Fortran). Irrelevant for C++ unless legacy Fortran code is involved.
    # -ffixed-line-length-none    - Allows for lines of any length in the code when using fixed-form source layout, for compatibility with old Fortran.
    # -ffpe-summary=none          - Disables the summary of floating-point exceptions, no runtime output on overflow or illegal operations
    # -fvect-cost-model=unlimited - Enable aggressive loop vectorization (unrolling?)
    # -ftree-vectorize            - Convert loops to singular SIMD instructions when possible
    # -fpic                       - Ensures that the code can run at any address in memory without modification
    # -fprefetch-loop-arrays      - Fetch whole arrays in loops into cache whenever possible
    # -funroll-loops              - Unroll loops to increase performance
    # -march=native               - Use all hardware-specific features
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_options(-Ofast -fdefault-real-8 -fdefault-double-8 -ffixed-form -ffixed-line-length-none -ffpe-summary=none -fvect-cost-model=unlimited -ftree-vectorize -fpic -fprefetch-loop-arrays -funroll-loops -march=native -Ofast)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # TODO: Compile under Clang, some options are probably incompatible
        add_compile_options(-Ofast -fdefault-real-8 -fdefault-double-8 -ffixed-form -ffixed-line-length-none -ffpe-summary=none -fvect-cost-model=unlimited -ftree-vectorize -fpic -fprefetch-loop-arrays -funroll-loops -march=native -Ofast)
    else()
        message(ERROR "Unsupported compiler detected!")
    endif()
    get_directory_property(GLOBAL_COMPILE_OPTIONS COMPILE_OPTIONS)
    message(STATUS "Global compile options: ${GLOBAL_COMPILE_OPTIONS}")
endfunction()