#include <fstream>

#include "../include/internal/blade.h"


std::vector<std::vector<int>> boxestochop(HeaderData hdr, 
                                          int norm, 
                                          float coord)
{
    std::vector<std::vector<int>> indexes;
    // For each level
    for (int lv = 0; lv <= hdr.maxlevel; lv++){
        std::vector<int> lvidx;
        // For each box index
        for (int i = 0; i < hdr.nboxes[lv]; i++){
            // If plane intersect with box
            if (coord >= hdr.boxes[lv][i].geolo[norm] &&
                coord <= hdr.boxes[lv][i].geohi[norm])
            {
                // Add to list of indexes
                lvidx.push_back(i);
            }
        }
        // Add to vector with each level
        indexes.push_back(lvidx);
    }
    return indexes;
}


boxslice slicebox(box bx, int lv, int mlv,
                  int cn, int cx, int cy, 
                  float coord, int varidx, float dn)
{
    // Output struct
    boxslice bs;
    // Store level
    bs.lv = lv;
    // Factor between current level and mlv
    int factor = 1 << (mlv - lv); // 2^(mlv - lv)
    // Indexes at max level
    bs.idxlo.push_back(bx.idxlo[cx] * factor);
    bs.idxlo.push_back(bx.idxlo[cy] * factor);
    bs.idxhi.push_back(bx.idxhi[cx] * factor);
    bs.idxhi.push_back(bx.idxhi[cy] * factor);
    // Box shape and size
    std::vector<int> shape;
    int bsize = 1;
    for (int i = 0; i < 3; i++){
        shape.push_back(bx.idxhi[i] - bx.idxlo[i]);
        bsize *= bx.idxhi[i] - bx.idxlo[i];
    }
    // Grid in normal direction
    std::vector<float> ngrid;
    for (int i = 0; i < shape[cn]; i++){
        // Lower box limit + half resolution (first point)
        // + resolution * point index
        ngrid.push_back(bx.geolo[cn] + dn/2 + i*dn);
    }
    // Open the binary file
    std::ifstream bfile(bx.bpath, std::ios::binary);
    // Go to the offset
    bfile.seekg(bx.offset, std::ios::beg);
    // skip a line
    std::string temp;
    std::getline(bfile, temp);
    // Go to the field
    bfile.seekg(bsize * sizeof(double) * varidx, std::ios::cur);
    // Buffer for the bytes
    std::vector<char> buffer(bsize * sizeof(double));
    // Read the field data with the box shape
    bfile.read(buffer.data(), buffer.size());
    bfile.close();
    // Vector for the double data
    std::vector<double> vec(bsize);
    // Recast as double
    for (int i = 0; i < bsize; i++){
        double value = *(reinterpret_cast<double*>(&buffer[i * sizeof(double)]));
        vec[i] = value; 
    }
    // Reshape with Fortran order
    std::vector<std::vector<std::vector<double>>> arr(shape[0],
                                                      std::vector<std::vector<double>>(shape[1],
                                                      std::vector<double>(shape[2])));
    // Really missing numpy right now
    int index = 0;
    for (int i = 0; i < shape[2]; i++){  // Z coordinate
        for (int j = 0; j < shape[1]; j++){ // Y coordinate
            for (int k = 0; k < shape[0]; k++){ // X coordinate
                // Fortran changes X then Y then Z
                arr[k][j][i] = vec[index];
                index++;
            } 
        }
    }
    // find out which indexes to keep in the normal direction
    int idxl = -1;
    int idxr = -1;
    // Case when the plane is between the last point
    // And the box limit
    if (coord > ngrid.back()){
        idxl = shape[cn] - 1;
        bs.nl = ngrid[idxl];
    }
    // Case when the plane is between the first point
    // And the box limit
    if (coord < ngrid[0]){
        idxr = 0;
        bs.nr = ngrid[idxr];
    }
    // Case when the plane lands on a point
    for (int i = 0; i < shape[cn]; i++){
        if (std::abs(coord - ngrid[i]) < 1e-16){
            idxl = i;
            bs.nl = ngrid[i];
            idxr = i;
            bs.nr = ngrid[i];
            break;
        }
    }
    // Case when the plane is between two points
    if (idxl == -1 && idxr == -1){
        for (int i = 0; i < shape[cn] - 1; i++){
            if (coord > ngrid[i] && coord < ngrid[i + 1]){
                idxl = i;
                bs.nl = ngrid[i];
                idxr = i + 1;
                bs.nr = ngrid[i + 1];
                break;
            }
        }
    }
    // Do the mandoline
    // Normal is x dir
    if (cn ==  0){
        // Right side
        if (idxr != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[idxr][i][j];
                }
            }
            bs.right = expand(out, factor);
        }
        // Left side
        if (idxl != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[idxl][i][j];
                }
            }
            bs.left = expand(out, factor);
        }
    }
    // Normal axis is y
    else if (cn == 1){
        // Right side
        if (idxr != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[i][idxr][j];
                }
            }
            bs.right = expand(out, factor);
        }
        // Left side
        if (idxl != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[i][idxl][j];
                }
            }
            bs.left = expand(out, factor);        
        }
    }
    // Normal axis is Z
    else if (cn == 2){
        // Right side
        if (idxr != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[i][j][idxr];
                }
            }
            bs.right = expand(out, factor);
        }
        // Left side
        if (idxl != -1){
            std::vector<std::vector<double>> out(shape[cx],
                                                 std::vector<double>(shape[cy]));
            for (int i = 0; i < shape[cx]; i++){
                for (int j = 0; j< shape[cy]; j++){
                    out[i][j] = arr[i][j][idxl];
                }
            }
            bs.left = expand(out, factor);
        }
    }
    return bs;
}


std::vector<std::vector<double>> expand(std::vector<std::vector<double>> arr,
                                        int factor)
{
    if (factor > 1) {
        std::vector<std::vector<double>> newarr(arr.size()*factor,
                                                std::vector<double>(arr[0].size()*factor));
        // For every row
        for (int i = 0; i < arr.size(); i++){
            // For every col
            for (int j = 0; j < arr[0].size(); j++){
                // For every repeated value
                for (int fi = 0; fi < factor; fi++){
                    for (int fj = 0; fj < factor; fj++){
                        newarr[i*factor + fj][j*factor + fi] = arr[i][j];
                    }
                }
            }
        }
        return newarr;
    }
    else {
    return arr;
    }
}
