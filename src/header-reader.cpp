#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "../include/internal/header-reader.h"

// Constructor to initialize the class with data from a file
HeaderData::HeaderData(const std::string& pfiledir) {
    ParseHeader(pfiledir);
    ParseBoxes();
}

// Function to read data from a file and store it in the vector
void HeaderData::ParseHeader(const std::string& pfiledir){
    // Local temporary variables to read file
    float num; // Float entries
    std::string temp; // String entires
    int tint; // Int entries
    // Append to plofile root
    pfile = pfiledir;
    const std::string& hfilepath = pfiledir + "/Header";
    // Create the header stream
    std::ifstream hfile(hfilepath);
    // Do the Parsing
    if (hfile.is_open()) {
        // First line is version
        hfile >> version;
        // Number of available fields
        hfile >> nfields;
        // Parse the fields
        for (int i = 0; i < nfields; i++){
            hfile >> temp;
            fieldnames.push_back(temp);
        }
        // Number of dimensions
        hfile >> ndims;
        // Plotfile time
        hfile >> time;
        // Maximum AMR Level
        hfile >> maxlevel;
        // Skip a line to use getline
        std::getline(hfile, temp);
        // Parse the domain bounds data
        // Low bounds of the geometry
        std::getline(hfile, temp);
        // Create a stream for the line
        // Because we dont know ndims at compile
        std::istringstream iss(temp);
        // Store in vector
        for (int i = 0; i < ndims; i++){
            iss >> num;
            geolow.push_back(num);
        }
        // Same same with high bounds
        std::getline(hfile, temp);
        // Create a stream for the line
        std::istringstream iss2(temp);
        // Store in vector
        for (int i = 0; i < ndims; i++){
            iss2 >> num;
            geohigh.push_back(num);
        }
        // Similar process for refinement ratios
        std::getline(hfile, temp);
        // Create a stream for the line
        std::istringstream iss3(temp);
        // Store in vector
        for (int i = 0; i < ndims; i++){
            iss3 >> tint;
            refratios.push_back(tint);
        }
        // Parse the grid sizes
        // Grid sizes look like:
        // ((0,0,0) (7,7,7) (0,0,0)) ((0,0,0) (15,15,15) (0,0,0))
        // Read the line
        std::getline(hfile, temp);
        // Stream object for the line
        std::istringstream ss(temp);
        // For each grid level
        for (int lv = 0; lv <= maxlevel; lv++){
            // Vector for the grid sizes at the current lv
            std::vector<int> grid;
            // For each entry (always 3)
            for (int i = 0; i < 3; i++){
                // Get the current entry tuple
                // (Split by spaces)
                std::getline(ss, temp, ' ');
                // Remove the parentheses
                while(temp[0] == '('){
                    temp = temp.substr(1);
                }
                while(temp.back() == ')'){
                    temp.pop_back();
                }
                // Create a stream for the tuple
                std::istringstream sst(temp);
                // For each dimension
                for (int dim = 0; dim < ndims; dim++){
                    // Get an number
                    std::getline(sst, temp, ',');
                    // Convert to int
                    tint = std::stoi(temp);
                    // Case when we store the data
                    if (i == 1) {
                        grid.push_back(tint + 1);
                    }
                    // Assert other values are zero
                    else {
                        if (tint != 0) {
                            throw std::runtime_error("It seems this plotfile "
                                                     "is not cell centered data");
                        }
                    }

                } // end for each dim
            }  // end for each tuple at level
            grids.push_back(grid);
        } // end for each level

        // Read the current step at each level
        std::getline(hfile, temp);
        // Stream for the step data
        std::istringstream ss_steps(temp);
        for (int lv = 0; lv <= maxlevel; lv++){
            // Step is read into temp
            std::getline(ss_steps, temp);
            lvsteps.push_back(std::stoi(temp));
        }
        // Read the grid resolutions
        // For each AMR Level
        for (int lv = 0; lv <= maxlevel; lv++){
            // Vector for current level resolutions
            std::vector<float> lvres;
            // Line contains the data for each coordinate
            std::getline(hfile, temp);
            // Stream for the line
            std::istringstream ss_res(temp);
            // For each coordinate
            for (int i = 0; i < ndims; i++){
                // Get a string of a float
                std::getline(ss_res, temp);
                // Convert and add to vector
                lvres.push_back(stof(temp));
            }
            res.push_back(lvres);
        }
        // Get what the coordinate system is
        std::getline(hfile, temp);
        std::istringstream ss_coord(temp);
        ss_coord >> coordsys;
        // Sanity Check for the parser
        std::getline(hfile, temp);
        std::istringstream ss_ptest(temp);
        ss_ptest >> tint;
        if (tint != 0){
            throw std::runtime_error("Something went wrong, "
                                     "we are not where we are " 
                                     "supposed to be, stopping...");
        }
        // Read the boxes
        float tlo; // Lower bound
        float thi; // Higher bound
        int level; // Level read from file
        int nbox; // number of boxes at lv
        float lvtime; // time at level
        std::string lvroot;
        // For each grid level
        for (int lv = 0; lv <= maxlevel; lv++) {
            // Parse the cell data for the level
            std::getline(hfile, temp);
            std::istringstream ss_lvhead(temp);
            // Current level number of boxes in level
            ss_lvhead >> level >> nbox >> lvtime;
            // Check if we get what we expect
            if (level != lv) {
                std::cout << level << " = " << lv << std::endl;
                throw std::runtime_error("The first number was supposed "
                                         "the current AMR Level");
            }
            // Store by level data
            nboxes.push_back(nbox);
            lvtimes.push_back(lvtime);
            // Line with the step number at lv
            std::getline(hfile, temp);
            std::istringstream ss_lvstep(temp);
            ss_lvstep >> tint;
            // Check if all is good
            if (tint != lvsteps[lv]) {
                std::cout << "Current Level: " << lv << std::endl;
                std::cout << "Step read: " << tint << std::endl;
                std::cout << "Step saved: " << lvsteps[lv] << std::endl;
                throw std::runtime_error("This line was supposed to be "
                                         "the step at the current level");
            }
            // Vector of boxes for the level
            std::vector<box> lvboxes;
            // For each box
            for (int i = 0; i < nbox; i++) {
                // Box instance
                box bx;
                // For each dimension
                for (int d = 0; d < ndims; d++) {
                    // Read low and high
                    std::getline(hfile, temp);
                    std::istringstream ss_box(temp);
                    ss_box >> tlo >> thi;
                    // Add to struct
                    bx.geolo.push_back(tlo);
                    bx.geohi.push_back(thi);
                }
                // Add to the vector of boxes for the level
                lvboxes.push_back(bx);
            }
            // Add to the boxes attribute
            boxes.push_back(lvboxes);
            // Read the cell root
            getline(hfile, lvroot);
            // Add to class
            lvroots.push_back(lvroot);
        }
        hfile.close();

    } 
    else {
        throw std::runtime_error("Unable to read plotfile/Header");
    }
}

// Read the Boxes meta data from level headers
void HeaderData::ParseBoxes(){
    // Temporary variables
    std::string line;
    std::string ttup;
    std::string tstr;
    int tint;
    // For each level
    for (int lv = 0; lv <= maxlevel; lv++){
        // Construct the box header path
        const std::string& cfilepath = pfile + "/" + lvroots[lv] + "_H";
        // Stream for the box header
        std::ifstream cfile(cfilepath);
        // Skip two first lines
        getline(cfile, line);
        getline(cfile, line);
        // Check if all is good
        getline(cfile, line);
        std::istringstream ss_nfields(line);
        ss_nfields >> tint;
        if (tint != nfields) {
            std::cout << "In file: " << cfilepath << std::endl;
            std::cout << "Read: " << tint << std::endl;
            std::cout << "Expect: " << nfields << std::endl;
            throw std::runtime_error("Wrong field number in file");
        }
        // Skip line always 0
        getline(cfile, line);
        // Skip indexes start
        getline(cfile, line);
        // For each box
        for (int i = 0; i < nboxes[lv]; i++){
            // Read the indexes
            getline(cfile, line);
            // Stream the line
            std::istringstream ss_idx(line);
            // For each index tuple (always 3)
            for (int j = 0; j < 3; j++){
                // Get the tuple (space separated)
                getline(ss_idx, ttup, ' ');
                // Remove the parentheses
                while(ttup[0] == '('){
                    ttup = ttup.substr(1);
                }
                while(ttup.back() == ')'){
                    ttup.pop_back();
                }
                // Create a stream for the tuple
                std::istringstream ss_tup(ttup);
                // For each dimension
                for (int d = 0; d < ndims; d++){
                    // Get an number
                    std::getline(ss_tup, tstr, ',');
                    // Convert to int
                    tint = std::stoi(tstr);
                    // Case for low indexes
                    if (j == 0) {
                        boxes[lv][i].idxlo.push_back(tint);
                    }
                    // Case for high indexes
                    else if (j == 1){
                        boxes[lv][i].idxhi.push_back(tint+1);
                    }
                    // Could assert tint == 0 here when
                    // (Meaning data is cell centered)
                    // j == 2 but we already check in 
                    // HeaderData::ParseHeader
                } // End for each dimension
            } // End for each index tuple
        } // End for each box
       // Skip two lines 
       getline(cfile, line);
       getline(cfile, line);
       long int ltint;
       // For each box
       for (int i = 0; i < nboxes[lv]; i++){
           // Get the file info
           getline(cfile, line);
           // Stream the line
           std::istringstream ss_finfo(line);
           // Line looks like:
           // FabOnDisk: Cell_D_00001 155735
           ss_finfo >> ttup >> tstr >> ltint;
           // Store the offset
           boxes[lv][i].offset = ltint;
           // Take only the Level dir
           std::size_t pos = lvroots[lv].find("/");
           // Construct the path to cell data
           boxes[lv][i].bpath = pfile + "/" + lvroots[lv].substr(0, pos);
           boxes[lv][i].bpath += "/";
           boxes[lv][i].bpath += tstr;
       }

    } // End for each level
}  // End function

// Show some info about the plotfile
void HeaderData::plotinfo() const{
    std::cout << "Plotfile Version: " << version << std::endl;
    std::cout << "Number of fields: " << nfields << std::endl;
    std::cout << "Domain bounds:" << std::endl;
    for (int i = 0; i < ndims; i++){
        std::cout << coordinates[i] << ": " << 1000. * geolow[i] << " to ";
        std::cout << 1000. * geohigh[i] << " mm" << std::endl;
    }
    std::cout << "Refinement ratios:" << std::endl;
    for (int i = 0; i <= maxlevel; i++){
        std::cout << "Level " << i << ": " << refratios[i] << std::endl;
    }
    std::cout << "Grid sizes:" << std::endl;
    for (int lv = 0; lv <= maxlevel; lv++){
        std::cout << "Level " << lv << ": ";
        for (int d=0; d < ndims; d++){
            if (d != ndims -1){
                std::cout << grids[lv][d] << " x ";
            }
            else {
                std::cout << grids[lv][d] << std::endl;
            }
        }
    }
    std::cout << "Plotfile steps:" << std::endl;
    for (int lv = 0; lv <= maxlevel; lv++){
        std::cout << "Level " << lv << ": ";
        std::cout << lvsteps[lv] << std::endl;
    }
    std::cout << "Grid resolutions:" << std::endl;
    for (int lv = 0; lv <= maxlevel; lv++){
        std::cout << "Level " << lv << ": ";
        std::cout << 1000.*res[lv][0] << " mm"  << std::endl;
    }
    std::cout << "Coordinate system: " << coordsystems[coordsys];
    std::cout << " (" << coordsys << ")" << std::endl;
    std::cout << "Number of boxes: " << std::endl;
    for (int lv = 0; lv <= maxlevel; lv++){
        std::cout << "Level " << lv << ": ";
        std::cout << nboxes[lv] << std::endl;
    }
    boxinfo();
}

// Display stats for the boxes
void HeaderData::boxinfo() const{
    // For each AMR Level
    for (int lv = 0; lv <= maxlevel; lv++){
        std::vector<float> lvgeomax = geolow;
        std::vector<float> lvgeomin = geohigh;
        // For each box
        for (int i = 0; i < nboxes[lv]; i++){
            // For each dimension
            for (int d = 0; d < ndims; d++){
                lvgeomax[d] = std::max(lvgeomax[d],
                                       boxes[lv][i].geohi[d]);
                lvgeomin[d] = std::min(lvgeomin[d],
                                       boxes[lv][i].geolo[d]);
            }
        }
        // Display for the level
        std::cout << "Boxes bounds at level " << lv << std::endl;
        for (int d = 0; d < ndims; d++){
            std::cout << coordinates[d] << ": ";
            std::cout << lvgeomin[d] << " to " << lvgeomax[d];
            std::cout << std::endl;
            
        }
    }
}

// Function to display the available fields
void HeaderData::fieldinfo() const{
    std::cout << "Available fields:" << std::endl;
    for (const std::string& field : fieldnames) {
        std::cout << field << std::endl;
    }
}
