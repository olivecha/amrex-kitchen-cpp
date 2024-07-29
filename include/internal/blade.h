#include "header-reader.h"

struct boxslice{
    // Left data
    std::vector<std::vector<double>> left;
    // Left normal coord
    double nl;
    // Right data
    std::vector<std::vector<double>> right;
    // Right normal coord
    double nr;
    // Slice indexes
    std::vector<int> idxlo;
    std::vector<int> idxhi;
    // AMR level
    int lv;
};

std::vector<std::vector<int>> boxestochop(HeaderData hdr, 
                                          int norm, 
                                          float coord);

std::vector<std::vector<double>> expand(std::vector<std::vector<double>> arr,
                                        int factor);

boxslice slicebox(box bx, int lv, int mlv,
                  int cn, int cx, int cy, 
                  float coord, int varidx, float dn);


