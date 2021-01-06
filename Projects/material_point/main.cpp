#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"
#include "mesh_query/mesh_query.h"
#include "Grid.h"


int main(int argc, char* argv[])
{
    using T = double;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;

    SimulationDriver<T, dim> driver;

    // important simulation parameters
    T youngs_modulus = 10000.0;
    T poissons_ratio = 0.4;
    T dt = 0.001;
    TV offset = TV(2.5, 2.5, 2.5); 
    std::vector<int> dimensions = {64, 64, 64};
    TV origin = TV::Zero();
    TV grid_max = TV(5.0, 5.0, 5.0);
    int grid_buffer = 5;
    T mass = 1.0;

    // deformation data
    std::vector<TM> F = std::vector<TM>();
    std::vector<TM> updateF = std::vector<TM>();

    // mesh handling
    char* filename;
    char temp_string[100];
    std::vector<T> obj_verts;
    std::vector<int> obj_mesh;

    int density = 8;

    if (argc < 2) 
    {
        std::cout << "Please input filename of an triangulated mesh (.obj) to simulate" << std::endl;
        exit(0);
    }

    // check validity of argument
    if (strlen(argv[1]) < 5) {
        std::cout << "Please input a valid filename (.obj)" << std::endl;
        exit(0);
    }

    if (argv[1][strlen(argv[1]) - 1] != 'j' || argv[1][strlen(argv[1]) - 2] != 'b' || 
        argv[1][strlen(argv[1]) - 3] != 'o' || argv[1][strlen(argv[1]) - 4] != '.') 
    {
        std::cout << "Please input a valid filename (.obj)" << std::endl;
        exit(0);
    }

    // argument is valid format
    filename = argv[1];
    FILE* obj = fopen(filename, "r");
    if (obj == NULL) {
        std::cout << "Error opening file. Does this file exist?" << std::endl;
        exit(0);
    }

    // update filename string
    filename[strlen(filename) - 4] = '\0';
    char* folder = strtok(filename, "/");
    char* file = folder;
    while (folder != NULL) {
        file = folder;
        folder = strtok(NULL, "/");
    }

    // density
    if (argv[2]) {
        try {
            density = std::stoi(argv[2]);
        }
        catch (const std::invalid_argument& ia) {}
    }

    if (argc == 4) {
        try {
            dimensions[0] = std::stoi(argv[3]);
            dimensions[1] = std::stoi(argv[3]);
            dimensions[2] = std::stoi(argv[3]);
        }
        catch (const std::invalid_argument& ia) {}
    }

    if (argc == 6) {
        try {
            dimensions[0] = std::stoi(argv[3]);
            dimensions[1] = std::stoi(argv[4]);
            dimensions[2] = std::stoi(argv[5]);
        }
        catch (const std::invalid_argument& ia) {}
    }

    // file exists: begin parsing
    int num_vertices = 0;
    int num_triangles = 0;
    while (!feof(obj)) {
        char* line = fgets(temp_string, 100, obj);
        char* token = strtok(line, " ");
        if (token != NULL && strcmp(token, "v") == 0) {
            num_vertices += 1;
            try {
                T x = std::stod(strtok(NULL, " "));
                T y = std::stod(strtok(NULL, " "));
                T z = std::stod(strtok(NULL, " "));

                obj_verts.push_back(x + offset[0]);
                obj_verts.push_back(y + offset[1]);
                obj_verts.push_back(z + offset[2]);
            }
            catch (const std::invalid_argument& ia) {
                std::cout << "Error parsing file: " << ia.what() << std::endl;
                fclose(obj);
                exit(0);
            }
        } else if (token!= NULL && strcmp(token, "f") == 0) {
            num_triangles += 1;
            try {
                int p1 = std::stoi(strtok(NULL, " ")) - 1;
                int p2 = std::stoi(strtok(NULL, " ")) - 1;
                int p3 = std::stoi(strtok(NULL, " ")) - 1;
                obj_mesh.push_back(p1);
                obj_mesh.push_back(p2);
                obj_mesh.push_back(p3);
            }
            catch (const std::invalid_argument& ia) {
                std::cout << "Error parsing file: " << ia.what() << std::endl;
                fclose(obj);
                exit(0);
            }
        }
    }

    fclose(obj);

    // construct mesh using Bridson code
    MeshObject* mesh = construct_mesh_object(num_vertices, obj_verts.data(), 
        num_triangles, obj_mesh.data());
    
    // set up grid and particle system
    Grid<T, dim> grid = Grid<T, dim>(origin, grid_max, dimensions, grid_buffer, density);
    grid.generateSamples(mesh, mass);

    // set up initial vectors needed by solver
    for (auto it = grid.parts->begin(); it != grid.parts->end(); ++it) {
        F.push_back(TM::Identity());
        updateF.push_back(TM::Zero());
    }
    std::cout << "This simulation has " << F.size() << " particles. Continue? (y/n)" << std::endl;
    char input;
    std::cin >> input;
    if (input != 'y' && input != 'Y') {
        std::cout << "Simulation aborted." << std::endl;
        exit(0);
    }

    // simulate
    driver.dt = dt;
    driver.filename = file;
    driver.mpm.grid = &grid;
    driver.mpm.parts = grid.parts;
    driver.mpm.F = F;
    driver.mpm.updateF = updateF;
    driver.mpm.youngs_modulus = youngs_modulus;
    driver.mpm.poissons_ratio = poissons_ratio;
    
    driver.mpm.computeConstants();
    driver.run(120);

    grid.parts->release();

    return 0;
}
