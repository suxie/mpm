#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Partio.h>
#include "mesh_query/mesh_query.h"
#include <vector>
#include <random>


template<class T, int dim>
class Grid {
public:
  using TV = Eigen::Matrix<T,dim,1>;

  TV origin;
  std::vector<int> dims;
  TV h;
  int buffer;
  int density;
  Partio::ParticlesDataMutable* parts;
  Partio::ParticleAttribute mH, xH, vH;

  std::vector<T> m;
  std::vector<TV> v;
  std::vector<TV> f;

  Grid(TV origin, TV max, std::vector<int> dims, int buffer, int density): 
    origin(origin), 
    dims(dims), 
    h(TV::Zero()), 
    buffer(buffer), 
    density(density),
    parts(Partio::create()), 
    m(std::vector<T>()), 
    v(std::vector<TV>())
  {
      TV lengths = max - origin;
      h[0] = lengths[0] / (T)dims[0];
      h[1] = lengths[1] / (T)dims[1];
      h[2] = lengths[2] / (T)dims[2];
      initialize();
  }

  void initialize() {
    for (int i = 0; i < (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1); i++) {
      m.push_back(0.0);
      v.push_back(TV::Zero());
      f.push_back(TV::Zero());
    }
  }
  
  void generateSamples(MeshObject* mesh, T mass) {
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    std::default_random_engine generator;
    std::uniform_real_distribution<T> distributionX(0.0, h[0]);
    std::uniform_real_distribution<T> distributionY(0.0, h[1]);
    std::uniform_real_distribution<T> distributionZ(0.0, h[2]);
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    generator.seed(seed);

    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    xH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    int partCount = 0;

    for (int i = 0; i < dims[0]; i++) {
      for (int j = 0; j < dims[1]; j++) {
        for (int k = 0; k < dims[2]; k++) {
          for (int n = 0; n < density; n++) {
            // initialize position
            T posX = origin[0] + i * h[0] + distributionX(generator);
            T posY = origin[1] + j * h[1] + distributionY(generator);
            T posZ = origin[2] + k * h[2] + distributionZ(generator);
            T pos[3] = {posX, posY, posZ};

            if (point_inside_mesh(pos, mesh)) {
              int idx = parts->addParticle();
              float* xp = parts->dataWrite<float>(xH, idx);
              float* vp = parts->dataWrite<float>(vH, idx);

              // initialize position and velocity              
              xp[0] = posX;
              xp[1] = posY;
              xp[2] = posZ;

              vp[0] = 0.0;
              vp[1] = 0.0;
              vp[2] = 0.0;

              partCount++;
            }
          }
        }
      }
    }

    // initialize mass distributed equally across all points
    for (int p = 0; p < partCount; p++) {
      float* mp = parts->dataWrite<float>(mH, p);
      mp[0] = mass / partCount;
    }
  }

};

#endif