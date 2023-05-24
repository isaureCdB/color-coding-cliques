from cffi import FFI
ffi = FFI()


ffi.set_source("_get_msd", open("get_msd.cpp").read(), source_extension='.cpp',
               #extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"]
               )
ffi.cdef("""
void get_msd(
  int nc1, int nc2, int natom,
  double *c1, double *c2,
  double *msd
); 
"""
)

if __name__ == "__main__":
    ffi.compile()