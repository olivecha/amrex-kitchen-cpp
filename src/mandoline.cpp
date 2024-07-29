#include <mpi.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/internal/header-reader.h"
#include "../include/internal/blade.h"

// Command line arguments
const char* optstring = "f:v:n:c:l:h";
// -h help
// -f plotfile
// -v variable
// -n normal {0,1,2}
// -c coordinate
// -l limit level

int main(int argc, char* argv[]) {
    // Declare variables that will be broadcasted
    int limitlevel;
    int norm;
    int cx;
    int cy;
    float coord;
    size_t varidx;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes and process rank
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Do this in serial
    // Parse command line arguments
    int c;
    std::string fname;
    std::string var = "density";
    norm = 0; // Normal direction
    limitlevel = -1;
    coord = 0.; // Slice coordinate
    bool input_coord = false;
    // Do the parsing
    while ((c = getopt(argc, argv, optstring)) != -1){
        switch (c){
            case 'h':
                std::cout << "Mandoline" << std::endl;
                std::cout << "_________" << std::endl;
                std::cout << "Options:" << std::endl;
                std::cout << "-h: Display this message" << std::endl;
                std::cout << "-f: Directory of the plotfile to slice"; 
                std::cout << std::endl;
                std::cout << "-v: Variables to slice" << std::endl;
                std::cout << "    (Defaults to density)" << std::endl;
                std::cout << "-n: Normal coordinate used" << std::endl;
                std::cout << "    0: x, 1: y, 2: z" << std::endl;
                std::cout << "    (Defaults to 0 when 3D)" << std::endl;
                std::cout << "-c: Slice plane coordinate" << std::endl;
                std::cout << "    (Defaults to domain center)" << std::endl;
                std::cout << "-l: Maximum level to consider" << std::endl;
                std::cout << "    (Defaults to max level)" << std::endl;
                return 0;
            case 'f':
                fname = optarg;
                break;
            case 'v':
                var = optarg;
                break;
            case 'n':
                norm = std::stoi(optarg);
                break;
            case 'c':
                coord = std::stof(optarg);
                input_coord = true;
                break;
            case 'l':
                limitlevel = std::stoi(optarg);
                break;
            case '?':
                std::cout << "Argument: " << c << std::endl;
                throw std::runtime_error("Invalid argument");
                break;
        }
    }
    // Read the header data
    HeaderData hdr(fname);
    // Check that the variable is available
    varidx = std::find(hdr.fieldnames.begin(),
                           hdr.fieldnames.end(),
                           var) - hdr.fieldnames.begin();
    if (varidx >= hdr.fieldnames.size())
    {
        if (rank == 0){
        std::string err = "Variable {" + var;
        err += "} not in plotfile";
        throw std::runtime_error(err);
        }
    }
    
    // Define the indexes of the plane coords
    cx = -1;
    if (hdr.ndims == 2){
        cx = 0;
        cy = 1;
        if (rank == 0){
        std::cout << "No need to slice its 2D" << std::endl;
        }
        return 0;
    }
    // Define the indexes of the slice plane
    // if normal = 0, then cx = 1, cy = 2, etc..
    for (int i = 0; i < 3; i++){
        if (norm != i){
            if (cx == -1){
                cx = i;
            }
            else {
                cy = i;
            }
        }
    }
    // Default coordinate to domain center
    if (!input_coord){
        float width = hdr.geohigh[norm] - hdr.geolow[norm];
        coord = hdr.geolow[norm] + width/2;
    }
    // Default limit level to max level
    if (limitlevel == -1){
        limitlevel = hdr.maxlevel;
    }
    // This would fuck stuff up
    else if (limitlevel > hdr.maxlevel){
        if (rank == 0){
            std::cout << "Max level cannot be higher ";
            std::cout << "than plotfile max level" << std::endl;
        }
        limitlevel = hdr.maxlevel;
        
    }
    // Print the slicing parameters
    if (rank == 0){
        std::cout << "Slicing " << hdr.fieldnames[varidx];
        std::cout << " at " << hdr.coordinates[norm];
        std::cout << " = " << coord << std::endl;
        std::cout << "with max level = ";
        std::cout << limitlevel << std::endl;
    }
    // Empty local matrices to perform the interpolation
    // Left side 
    std::vector<std::vector<double>> local_left(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Left normal coordinate
    std::vector<std::vector<double>> local_nleft(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Right side
    std::vector<std::vector<double>> local_right(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Right normal coordinate
    std::vector<std::vector<double>> local_nright(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Empty global matrices to perform the interpolation
    // Left side 
    std::vector<std::vector<double>> global_left(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Left normal coordinate
    std::vector<std::vector<double>> global_nleft(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Right side
    std::vector<std::vector<double>> global_right(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // Right normal coordinate
    std::vector<std::vector<double>> global_nright(hdr.grids[limitlevel][cy],
                                     std::vector<double>(hdr.grids[limitlevel][cx], -1.));

    // Find the box indexes in the slice
    std::vector<std::vector<int>> indexes = boxestochop(hdr, norm, coord);
    
    // For each level
    for (int lv = 0; lv <= limitlevel; lv++){
        // Find out how many indexes to send each process
        int nindexes = indexes[lv].size();
        int elements_per_proc = nindexes / size;
        int remainder = nindexes % size;
        if (rank == 0){
            std::cout << "dividing " << nindexes << " boxes between " << size;
            std::cout << " processes" << std::endl;
            std::cout << "Each process has " << elements_per_proc;
            std::cout << " boxes to read" << std::endl;
        }
        // Allocate space for local indexes
        std::vector<int> local_indexes(elements_per_proc + (rank < remainder ? 1 : 0));
        // Send out shares of indexes
        MPI_Scatter(&indexes[lv][0], elements_per_proc, MPI_INT,
                    &local_indexes[0], elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
        
        // Read the boxes
        for (int idx : local_indexes){
            boxslice bsl = slicebox(hdr.boxes[lv][idx], lv, limitlevel, 
                                    norm, cx, cy,
                                    coord, varidx, hdr.res[lv][norm]);
            int size_x = bsl.idxhi[0] - bsl.idxlo[0];
            int size_y = bsl.idxhi[1] - bsl.idxlo[1];
            // If there is a left side
            if (bsl.left.size() > 0){
                for (int i = 0; i < size_x; i++){
                    for (int j = 0; j < size_y; j++){
                        // Transpose here
                        local_left[bsl.idxlo[1] + j][bsl.idxlo[0] + i] = bsl.left[i][j];
                        local_nleft[bsl.idxlo[1] + j][bsl.idxlo[0] + i] = bsl.nl;
                    }
                }
            }
            // If there is a right side
            if (bsl.right.size() > 0){
                for (int i = 0; i < size_x; i++){
                    for (int j = 0; j < size_y; j++){
                        // Transpose here
                        local_right[bsl.idxlo[1] + j][bsl.idxlo[0] + i] = bsl.right[i][j];
                        local_nright[bsl.idxlo[1] + j][bsl.idxlo[0] + i] = bsl.nr;
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
    }


    //     // Output to store interpolation
    //     std::vector<std::vector<double>> output(hdr.grids[limitlevel][cy],
    //                                      std::vector<double>(hdr.grids[limitlevel][cx], -1.));
    // }
    // 
    // // For each level
    // for (int lv = 0; lv <= limitlevel; lv++){
    //     // For each box index
    //     for (int idx : indexes[lv]){
    //         bsl = slicebox(hdr.boxes[lv][idx], lv, limitlevel, 
    //                        norm, cx, cy,
    //                        coord, varidx, hdr.res[lv][norm]);
    //         // slice size
    //                 }
    //             }
    //         }
    //         std::cout << "Slice " << lv << " " << idx << " of " << indexes[lv].back() << std::endl;
    //     }
    // }

    // // For each row
    // for (int i = 0; i < output.size(); i++){
    //     // For each columns
    //     for (int j = 0; j < output[0].size(); j++){
    //         // Case when the points are the same
    //         if (std::abs(nleft[i][j] - nright[i][j]) < 1e-16){
    //             output[i][j] = (left[i][j] + right[i][j])/2;
    //         }
    //         // Linear interpolation
    //         else {
    //             output[i][j] = (left[i][j] * (nright[i][j] - coord)
    //                            + right[i][j] * (coord - nleft[i][j]))
    //                            / (nright[i][j] - nleft[i][j]);
    //         }
    //     }
    // }

    // std::cout << "Reached" << std::endl;

    // std::ofstream ofile("out.b", std::ios::binary);
    // if (!ofile.is_open()) {
    //   // Handle error opening file
    //   return 1;
    // }
    // for (const auto& inner_vector : output) {
    //   for (double value : inner_vector) {
    //     ofile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    //   }
    // }
    // ofile.close();


    // return 0;
    //}
