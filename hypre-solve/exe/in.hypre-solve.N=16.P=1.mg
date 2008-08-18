grid 0 -1 0 -4e+09 -4e+09 -4e+09  4e+09 4e+09 4e+09  0 0 0  16 16 16
grid 1 0 0 -2e+09 -2e+09 -2e+09  2e+09 2e+09 2e+09  8 8 8  16 16 16
dimension 3
domain    3 -4e9 -4e9 -4e9  4e9 4e9 4e9
boundary  dirichlet
sphere    5.993985e27 6.378137e8 0.0 0.0 0.0
point     1e43   -1.25e+08 -1.25e+08 -1.25e+08   1
point     1e43   1.25e+08 -1.25e+08 -1.25e+08   1
point     1e43   -1.25e+08 1.25e+08 -1.25e+08   1
point     1e43   1.25e+08 1.25e+08 -1.25e+08   1
point     1e43   -1.25e+08 -1.25e+08 1.25e+08   1
point     1e43   1.25e+08 -1.25e+08 1.25e+08   1
point     1e43   -1.25e+08 1.25e+08 1.25e+08   1
point     1e43   1.25e+08 1.25e+08 1.25e+08   1
discret constant
solver fac
dump_x true
dump_a true
dump_b true
