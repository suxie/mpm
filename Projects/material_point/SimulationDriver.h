#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MPMSystem.h"
#include <functional>


template<class T, int dim>
class SimulationDriver {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MPMSystem<T, dim> mpm;
    T dt;
    TV gravity;
    std::string filename;

    std::function<void(T, T)> helper = [&](T, T) {};

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;
    }

    void run(const int max_frame)
    {
        T accumulate_t = 0;
        mkdir("output/", 0777);
        std::string output_folder = "output/" + filename;
        mkdir(output_folder.c_str(), 0777);
        std::string output = output_folder + "/" + filename + "_0.bgeo";
        Partio::write(output.c_str(), *(mpm.grid->parts));
        for(int frame=1; frame<=max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;
            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                std::cout << "Step " << step << std::endl;
                helper(accumulate_t, dt);
                advanceOneStep();
                accumulate_t += dt;
            }
            mkdir("output/", 0777);
            std::string output_folder = "output/" + filename;
            mkdir(output_folder.c_str(), 0777);
            output = output_folder + "/" + filename + "_" + std::to_string(frame) + ".bgeo";
            Partio::write(output.c_str(), *(mpm.grid->parts));
            std::cout << std::endl;
        }
    }

    void advanceOneStep()
    {
        mpm.zeroValues();
        mpm.interpolateMassVelToGrid();
        mpm.setGridVelocities();

        applyGravity();

        mpm.computeForces();
        mpm.addForceGridVelocities(dt);
        
        mpm.interpolateVelToParticles();
        mpm.computeNewF(dt);

        // move particles
        Partio::ParticlesDataMutable* parts = mpm.parts;
        for (int idx = 0; idx < mpm.parts->numParticles(); idx++) {
            float* xp = parts->dataWrite<float>(mpm.grid->xH, idx);
            float* vp = parts->dataWrite<float>(mpm.grid->vH, idx);
            
            xp[0] += vp[0] * dt;
            xp[1] += vp[1] * dt;
            xp[2] += vp[2] * dt;
        }
        
    }

    void applyGravity() {
        for (unsigned int i = 0; i < mpm.grid->v.size(); i++) {
            mpm.grid->v[i] += gravity * dt;
        }
    }
};
