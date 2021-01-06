#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include "Grid.h"
#include "Partio.h"
#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class MPMSystem {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;

    Grid<T,dim>* grid;
    Partio::ParticlesDataMutable* parts;
    
    std::vector<TM> F = std::vector<TM>();
    std::vector<TM> updateF = std::vector<TM>();
    T youngs_modulus;
    T poissons_ratio;
    T mu;
    T lambda;

    MPMSystem() {}

    void computeConstants() {
      mu = youngs_modulus / (1.0 + poissons_ratio);
      lambda = (youngs_modulus * poissons_ratio) / 
               ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
    }

    void zeroValues() {
      for (unsigned int i = 0; i < grid->m.size(); i++) {
        grid->m[i] = 0.0;
        grid->v[i].setZero();
        grid->f[i].setZero();
      }

      for (unsigned int i = 0; i < F.size(); i++) {
        updateF[i].setZero();
      }
    }

    T nHat(TV xp, TV xi, int idx) {
      T x = abs(xp[idx] - xi[idx]);
      x /= grid->h[idx];
      
      if (x < 0.5) {
        return 0.75 - (x * x);
      } else if (x < 1.5) {
        return 0.5 * (1.5 - x) * (1.5 - x);
      } else {
        return 0.0;
      }
    }

    T nHatPrime(TV xp, TV xi, int idx) {
      T x = xp[idx] - xi[idx];
      x /= grid->h[idx];
      T xabs = abs(x);
      
      if (xabs < 0.5) {
        return -2.0 * x;
      } else if (xabs < 1.5) {
        return x * (1.0 - 1.5 / xabs);
      } else {
        return 0.0;
      }
    }

    T Ni(TV xp, TV xi) {
      return nHat(xp, xi, 0) * nHat(xp, xi, 1) * nHat(xp, xi, 2);
    }

    TV NiGrad(TV xp, TV xi) {
      T x = nHatPrime(xp, xi, 0) * nHat(xp, xi, 1) * nHat(xp, xi, 2);
      T y = nHat(xp, xi, 0) * nHatPrime(xp, xi, 1) * nHat(xp, xi, 2);
      T z = nHat(xp, xi, 0) * nHat(xp, xi, 1) * nHatPrime(xp, xi, 2);
      return TV(x, y, z);
    }

    void interpolateMassVelToGrid() {
      for (int idx = 0; idx < parts->numParticles(); idx++) {
        float* xp = parts->dataWrite<float>(grid->xH, idx);
        float* mp = parts->dataWrite<float>(grid->mH, idx);
        float* vp = parts->dataWrite<float>(grid->vH, idx);
        int x = calculateBasePos(xp, 0);
        int y = calculateBasePos(xp, 1);
        int z = calculateBasePos(xp, 2);

        for (int i = 0; i < 3; i++) {
          int x_new = x + i;
          if (x_new > grid->dims[0] - grid->buffer || x_new < grid->buffer) {
            continue;
          }
          for (int j = 0; j < 3; j++) {
            int y_new = y + j;
            if (y_new > grid->dims[1] - grid->buffer || y_new < grid->buffer) {
              continue;
            }
            for (int k = 0; k < 3; k++) {
              int z_new = z + k;
              if (z_new > grid->dims[2] - grid->buffer || z_new < grid->buffer) {
                continue;
              }
              int index = getIndex(x_new, y_new, z_new);
              T weight = Ni(TV(xp[0], xp[1], xp[2]), get_pos(x_new, y_new, z_new));
              grid->m[index] += weight * mp[0];
              grid->v[index] += weight * mp[0] * TV(vp[0], vp[1], vp[2]);
            }
          }
        }
      }
    }

    void setGridVelocities() {
      for (unsigned int i = 0; i < grid->m.size() ; i++) {
        if (grid->m[i] == 0.0) {
          grid->v[i].setZero();
        } else {
          grid->v[i] /= grid->m[i]; 
        }
      }
    }

    void interpolateVelToParticles() {
      for (int idx = 0; idx < parts->numParticles(); idx++) {
        float* xp = parts->dataWrite<float>(grid->xH, idx);
        float* vp = parts->dataWrite<float>(grid->vH, idx);
        int x = calculateBasePos(xp, 0);
        int y = calculateBasePos(xp, 1);
        int z = calculateBasePos(xp, 2);

        // zero velocity
        vp[0] = 0.0;
        vp[1] = 0.0;
        vp[2] = 0.0;

        for (int i = 0; i < 3; i++) {
          int x_new = x + i;
          if (x_new > grid->dims[0] - grid->buffer || x_new < grid->buffer) {
            continue;
          }
          for (int j = 0; j < 3; j++) {
            int y_new = y + j;
            if (y_new > grid->dims[1] - grid->buffer || y_new < grid->buffer) {
              continue;
            }
            for (int k = 0; k < 3; k++) {
              int z_new = z + k;
              if (z_new > grid->dims[2] - grid->buffer || z_new < grid->buffer) {
                continue;
              }
              
              int index = getIndex(x_new, y_new, z_new);
              T weight = Ni(TV(xp[0], xp[1], xp[2]), get_pos(x_new, y_new, z_new));
              TV weightedV = grid->v[index] * weight;

              vp[0] += weightedV[0];
              vp[1] += weightedV[1];
              vp[2] += weightedV[2];
            }
          }
        }
      }
    }

    TM computeStress(int i) {
      // calculate polar SVD
      Eigen::JacobiSVD<TM>svd(F[i], Eigen::ComputeFullV | Eigen::ComputeFullU);
      TM u = TM(svd.matrixU());
      TM v = TM(svd.matrixV());
      TV sigma = TV(svd.singularValues());
      if (u.determinant() < 0) {
        for (int row = 0; row < dim; row++) {
          u(row, dim - 1) *= -1.0;
        }
        sigma[dim - 1] *= -1.0;
      }
      if (v.determinant() < 0) {
        for (int row = 0; row < dim; row++) {
          v(row, dim - 1) *= -1.0;
        }
        sigma[dim - 1] *= -1.0;
      }

      TM R = u * v.transpose();

      T J = F[i].determinant();

      TM P = 2 * mu * (F[i] - R) + lambda * (J - 1.0) * J * F[i].inverse().transpose();
      return P;
    }

    void computeForces() {
      T volume = (grid->h[0] * grid->h[1] * grid->h[2]) / grid->density;
      // std::cout << volume << std::endl;
      for (int idx = 0; idx < parts->numParticles(); idx++) {
        float* xp = parts->dataWrite<float>(grid->xH, idx);
        int x = calculateBasePos(xp, 0);
        int y = calculateBasePos(xp, 1);
        int z = calculateBasePos(xp, 2);

        TM P = computeStress(idx);
        TM PF = P * F[idx].transpose();

        for (int i = 0; i < 3; i++) {
          int x_new = x + i;
          if (x_new > grid->dims[0] - grid->buffer || x_new < grid->buffer) {
            continue;
          }
          for (int j = 0; j < 3; j++) {
            int y_new = y + j;
            if (y_new > grid->dims[1] - grid->buffer || y_new < grid->buffer) {
              continue;
            }
            for (int k = 0; k < 3; k++) {
              int z_new = z + k;
              if (z_new > grid->dims[2] - grid->buffer || z_new < grid->buffer) {
                continue;
              }
              
              int index = getIndex(x_new, y_new, z_new);
              TV weight = NiGrad(TV(xp[0], xp[1], xp[2]), get_pos(x_new, y_new, z_new));

              TV test = volume * PF * weight;

              grid->f[index] -= test;
              
            }
          }
        }
      }
    }

    void addForceGridVelocities(T dt) {
      for (unsigned int i = 0; i < grid->m.size() ; i++) {
        if (grid->m[i] > 0) {
          grid->v[i] += dt * grid->f[i] / grid->m[i];
        }
      }
    }

    void computeNewF(T dt) {
      for (int idx = 0; idx < parts->numParticles(); idx++) {
        float* xp = parts->dataWrite<float>(grid->xH, idx);
        int x = calculateBasePos(xp, 0);
        int y = calculateBasePos(xp, 1);
        int z = calculateBasePos(xp, 2);

        for (int i = 0; i < 3; i++) {
          int x_new = x + i;
          if (x_new > grid->dims[0] - grid->buffer || x_new < grid->buffer) {
            continue;
          }
          for (int j = 0; j < 3; j++) {
            int y_new = y + j;
            if (y_new > grid->dims[1] - grid->buffer || y_new < grid->buffer) {
              continue;
            }
            for (int k = 0; k < 3; k++) {
              int z_new = z + k;
              if (z_new > grid->dims[2] - grid->buffer || z_new < grid->buffer) {
                continue;
              }
              
              int index = getIndex(x_new, y_new, z_new);
              TV weight = NiGrad(TV(xp[0], xp[1], xp[2]), get_pos(x_new, y_new, z_new));

              updateF[idx] += grid->v[index] * weight.transpose();
            }
          }
        }
      }

      for (unsigned int i = 0; i < F.size(); i++) {
        F[i] = (TM::Identity() + updateF[i] * dt) * F[i];
      }
    }

private:
    int calculateBasePos(float* xp, int idx) {
      return (int) ((xp[idx] - grid->origin[idx] - 0.5 * grid->h[idx]) / (grid->h[idx]));
    }

    int getIndex(int x, int y, int z) {
      return x * (grid->dims[1] + 1) * (grid->dims[2] + 1) + y * (grid->dims[2] + 1) + z;
    }

    TV get_pos(int x, int y, int z) {
      T posX = grid->origin[0] + x * grid->h[0];
      T posY = grid->origin[1] + y * grid->h[1];
      T posZ = grid->origin[2] + z * grid->h[2];
      return TV(posX, posY, posZ);
    }
};
