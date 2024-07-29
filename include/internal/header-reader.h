#ifndef HEADERDATA_H
#define HEADERDATA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


struct box {
    std::string bpath; // Path to field data
    long int offset; // offset in field data file
    std::vector<float> geolo; // Lower box bounds (ndims)
    std::vector<float> geohi; // Upper box bounds (ndims)
    std::vector<int> idxlo; // Lower indexes (ndims)
    std::vector<int> idxhi; // Upper indexes (ndims)
};


class HeaderData {
    private:
        std::string version;

    public:
        // Struct for cell meta data
        // Box data
        std::vector<int> nboxes; // number of cells at each level
        std::vector<std::vector<box>> boxes; // MultiFAB boxes data from headers
        // Defined values
        // indexes of <coordinates> match the string values
        std::vector<std::string> coordinates = {"x", "y", "z"};
        std::vector<std::string> coordsystems = {"cartesian", "cylindrical", "spherical"};
        // Values read from the plotfile
        std::string pfile;
        int nfields; // Number of fields in plotfile
        int ndims; // Number of dimensions
        float time; // Current plotfile time
        std::vector<float> lvtimes; // Per level times
        int maxlevel; // Plotfile max AMR Level
        std::vector<float> geolow; // Lower geo bounds {x, y, z}
        std::vector<float> geohigh; // Upper geo bounds {x, y, z}
        std::vector<int> refratios; // Refinement ratios (size = maxlevel + 1)
        std::vector<std::string> fieldnames; // Field names (size = nfields)
        std::vector<std::vector<int>> grids; // Grid sizes (maxlevel * ndims)
        std::vector<int> lvsteps; // current step  at each level
        std::vector<std::vector<float>> res; // Grid resolution (maxlevel * ndims)
        std::vector<std::string> lvroots; // path roots of the level data
        int coordsys; // Coordinate system index in HeaderData->coordsystems
        // Constructor
        HeaderData(const std::string& pfiledir);
        // Functions
        void ParseHeader(const std::string& pfiledir);
        void ParseBoxes();
        void plotinfo() const;
        void fieldinfo() const;
        void boxinfo() const;
};

//
//for (int lv = 0; lv <= maxlevel; lv++){
//    std::vector<box> lvboxes;
//
//    while(reading_boxes){
//        box tmpbox;
//        float tmpf;
//        for (int i = 0; i < ndims; i++){
//            hfile >> tmpf;
//            tmpbox.geolo.push_back(tmpf);
//            hfile >> tmpf;
//            tmpbox.geohi.push_back(tmpf);
//        }
//    }
//    // In Another file
//    int tint;
//    for (int i=0, i < nbox[lv]; i++){
//        ofile >> tint;
//        tmpbox.offset = tint;
//    }
//    lvboxes.push_back(tmpbox);
//}
//boxes.push_back(lvboxes)




#endif /* HEADERDATA */

