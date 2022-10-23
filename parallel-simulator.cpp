#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>

#include <omp.h>
#include <numeric>

// TASK 2

// NOTE: You may modify this class definition as you see fit, as long as the
// class name, and type of simulateStep and buildAccelerationStructure remain
// the same.

bool compareParticle(Particle p1, Particle p2) {
  return p1.position.x < p2.position.x;
}

// check for 4, 8, 16, 32, 64 --> larger leaf node size works better --> shallow quad tree

const int QuadTreeLeafSize = 64;
class ParallelNBodySimulator : public INBodySimulator {
public:
  // TODO: implement a function that builds and returns a quadtree containing
  // particles. You do not have to preserve this function type.
  std::vector<int> order;
  // std::vector<Particle> original;
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {

    
    auto qtn = new QuadTreeNode();
    if (particles.size() < QuadTreeLeafSize) {
      qtn->particles = particles;
      qtn->isLeaf = true;
      return std::unique_ptr<QuadTreeNode>(qtn);
    }

    // std::vector<Particle> tobesort(particles);
    // std::sort(tobesort.begin(), tobesort.end(), compareParticle);
    // std::cout << "sorted" << std::endl;

    Vec2 centre = Vec2((bmin.x+bmax.x)/2, (bmin.y+bmax.y)/2);

    // std::cout << "step 2" <<std::endl;
  

    // push_back is not thread safe, use something else
    std::vector<Particle> p0, p1, p2, p3;

    // for (int i = 0; i < (int)particles.size(); i++) {
    //   Particle pi = particles[i];
    //   Vec2 posi = pi.position;

    //   if (posi.x < centre.x && posi.y < centre.y) {
    //     p0.push_back(pi);
    //   } else if (posi.x >= centre.x && posi.y < centre.y) {
    //     p1.push_back(pi);
    //   } else if (posi.x < centre.x && posi.y >= centre.y) {
    //     p2.push_back(pi);
    //   } else {
    //     p3.push_back(pi);
    //   }
    // }


    // OMP tasks, 1 task for each child

    // omp parallel inside recursive is not optimal, put it at the top level
    // 
    
        #pragma omp task untied if(particles.size() > 2048)
        
        {
          
          for (int i=0; i < particles.size(); i++) {
            Particle pi = particles[i];
            Vec2 posi = pi.position;
            if (posi.x < centre.x && posi.y < centre.y) {
              p0.push_back(pi);
            }
          }
          qtn->children[0] = buildQuadTree(p0, bmin, centre);
        }
        
        #pragma omp task untied if(particles.size() > 2048)
        
        {
          for (int i=particles.size()-1; i > -1; i--) {
            Particle pi = particles[i];
            Vec2 posi = pi.position;
            
            if (posi.x >= centre.x && posi.y < centre.y) {
              p1.push_back(pi);
            }
          }
          qtn->children[1] = buildQuadTree(p1, Vec2((bmin.x+bmax.x)/2, bmin.y), Vec2(bmax.x, (bmin.y+bmax.y)/2)); 
        }
        
        #pragma omp task untied if(particles.size() > 2048)
        
        {
          for (int i=0; i < particles.size(); i++) {
            Particle pi = particles[i];
            Vec2 posi = pi.position;
            
            if (posi.x < centre.x && posi.y >= centre.y) {
              p2.push_back(pi);
            }
            
          }
          qtn->children[2] = buildQuadTree(p2, Vec2(bmin.x, (bmin.y+bmax.y)/2), Vec2((bmin.x+bmax.x)/2, bmax.y)); 
        }
        
        #pragma omp task untied if(particles.size() > 2048)
        
        {
          
          for (int i=0; i < particles.size(); i++) {
            Particle pi = particles[i];
            Vec2 posi = pi.position;
            
            if (posi.x >= centre.x && posi.y >= centre.y) {
              p3.push_back(pi);
            }
          }
          qtn->children[3] = buildQuadTree(p3, centre, bmax); 
        }

    // OMP 4 threads for 
    // #pragma omp parallel 
    // {
    //   // std::cout << "numthreads: " << omp_get_num_threads() << "size: " << particles.size() << std::endl;
      
    //   volatile bool flag=false;
    //   std::vector<Particle> priv0;

    //   #pragma omp for nowait
    //   for (int i = 0; i < (int)tobesort.size(); i++) {
    //     if (flag) continue;
    //     Particle pi = tobesort[i];
    //     Vec2 posi = pi.position;
    //     if (posi.x >= centre.x) {
    //       flag=true;
    //     }
        
    //     if (posi.x < centre.x && posi.y < centre.y) {
    //       priv0.push_back(pi);
    //     } 
    //   }
    //   #pragma omp critical
    //   p0.insert(p0.end(), priv0.begin(), priv0.end());
    // }
    // // for (auto &p: p0) {
    // //   std::cout << p.position.x << " " << p.position.y << std::endl;
    // // }
    // // std::cout << "-------" << std::endl;
    // qtn->children[0] = buildQuadTree(p0, bmin, centre);
    // #pragma omp parallel 
    // {
    //   volatile bool flag=false;
    //   std::vector<Particle> priv1;

    //   #pragma omp for nowait
    //   for (int i = tobesort.size()-1; i > -1; i--) {
    //     if (flag) continue;
    //     Particle pi = tobesort[i];
    //     Vec2 posi = pi.position;
    //     if (posi.x < centre.x) {
    //       flag=true;
    //     }
        
    //     if (posi.x >= centre.x && posi.y < centre.y) {
    //       priv1.push_back(pi);
    //     } 
    //   }
    //   #pragma omp critical
    //   p1.insert(p1.end(), priv1.begin(), priv1.end());
    // }
    // qtn->children[1] = buildQuadTree(p1, Vec2((bmin.x+bmax.x)/2, bmin.y), Vec2(bmax.x, (bmin.y+bmax.y)/2)); 

    // #pragma omp parallel 
    // {
    //   volatile bool flag=false;
    //   std::vector<Particle> priv2;

    //   #pragma omp for nowait
    //   for (int i = 0; i < (int)tobesort.size(); i++) {
    //     if (flag) continue;
    //     Particle pi = tobesort[i];
    //     Vec2 posi = pi.position;
    //     if (posi.x >= centre.x) {
    //       flag=true;
    //     }
        
    //     if (posi.x < centre.x && posi.y >= centre.y) {
    //       priv2.push_back(pi);
    //     }
    //   }
    //   #pragma omp critical
    //   p2.insert(p2.end(), priv2.begin(), priv2.end());
    // }
    // qtn->children[2] = buildQuadTree(p2, Vec2(bmin.x, (bmin.y+bmax.y)/2), Vec2((bmin.x+bmax.x)/2, bmax.y)); 
    
    // #pragma omp parallel 
    // {
    //   volatile bool flag=false;
    //   std::vector<Particle> priv3;

    //   #pragma omp for nowait
    //   for (int i = tobesort.size()-1; i > -1; i--) {
    //     if (flag) continue;;
    //     Particle pi = tobesort[i];
    //     Vec2 posi = pi.position;
    //     if (posi.x < centre.x) {
    //       flag=true;
    //     }
        
    //     if (posi.x >= centre.x && posi.y >= centre.y) {
    //       priv3.push_back(pi);
    //     }
    //   }
    //   #pragma omp critical
    //   p3.insert(p3.end(), priv3.begin(), priv3.end());
    // }
    // qtn->children[3] = buildQuadTree(p3, centre, bmax); 


    // printf("step 1 %d %d %d %d\n", p0.size(), p1.size(), p2.size(), p3.size());
    
    // #pragma omp parallel
    // {
    //   qtn->children[0] = buildQuadTree(p0, bmin, centre);
    //   qtn->children[1] = buildQuadTree(p1, Vec2((bmin.x+bmax.x)/2, bmin.y), Vec2(bmax.x, (bmin.y+bmax.y)/2)); 
    //   qtn->children[2] = buildQuadTree(p2, Vec2(bmin.x, (bmin.y+bmax.y)/2), Vec2((bmin.x+bmax.x)/2, bmax.y)); 
    //   qtn->children[3] = buildQuadTree(p3, centre, bmax); 
    // }

    // qtn->children[0] = buildQuadTree(p0, bmin, centre);
    // qtn->children[1] = buildQuadTree(p1, Vec2((bmin.x+bmax.x)/2, bmin.y), Vec2(bmax.x, (bmin.y+bmax.y)/2)); 
    // qtn->children[2] = buildQuadTree(p2, Vec2(bmin.x, (bmin.y+bmax.y)/2), Vec2((bmin.x+bmax.x)/2, bmax.y)); 
    // qtn->children[3] = buildQuadTree(p3, centre, bmax); 
    
    
    // #pragma omp parallel for
    // for (int i = 0; i < (int)particles.size(); i++) {
    //   Particle pi = particles[i];
    //   Vec2 posi = pi.position;

    //   if (posi.x < centre.x && posi.y < centre.y) {
    //     p0.push_back(pi);
    //   } else if (posi.x >= centre.x && posi.y < centre.y) {
    //     p1.push_back(pi);
    //   } else if (posi.x < centre.x && posi.y >= centre.y) {
    //     p2.push_back(pi);
    //   } else {
    //     p3.push_back(pi);
    //   }
    // }

    // std::cout << "step 3" <<std::endl;

    std::unique_ptr<QuadTreeNode> topleft, topright, btmleft, btmright;
    // #pragma omp parallel
    // {
    //   qtn->children[0] = buildQuadTree(p0, bmin, centre);
    //   qtn->children[1] = buildQuadTree(p1, Vec2((bmin.x+bmax.x)/2, bmin.y), Vec2(bmax.x, (bmin.y+bmax.y)/2)); 
    //   qtn->children[2] = buildQuadTree(p2, Vec2(bmin.x, (bmin.y+bmax.y)/2), Vec2((bmin.x+bmax.x)/2, bmax.y)); 
    //   qtn->children[3] = buildQuadTree(p3, centre, bmax); 
    // }
    qtn->isLeaf=false;

    // std::cout << "step 4" << std::endl;
    
    return std::unique_ptr<QuadTreeNode>(qtn);
  }

  // Do not modify this function type.
  virtual std::unique_ptr<AccelerationStructure>
  buildAccelerationStructure(std::vector<Particle> &particles) {
    // build quad-tree
    auto quadTree = std::make_unique<QuadTree>();

    // find bounds
    Vec2 bmin(1e30f, 1e30f);
    Vec2 bmax(-1e30f, -1e30f);

    for (auto &p : particles) {
      bmin.x = fminf(bmin.x, p.position.x);
      bmin.y = fminf(bmin.y, p.position.y);
      bmax.x = fmaxf(bmax.x, p.position.x);
      bmax.y = fmaxf(bmax.y, p.position.y);
    }

    quadTree->bmin = bmin;
    quadTree->bmax = bmax;

    // build nodes
    // std::sort(particles.begin(), particles.end(), compareParticle);

    // original = std::vector<Particle>(particles);
    // copy(particles.begin(), particles.end(), back_inserter(original));
    // order.reserve(particles.size());
    // for (int i=0; i<particles.size(); i++) {
    //   order[i] = i;
    // }
    // std::sort(order.begin(), order.end(),
    //    [particles](size_t i1, size_t i2) { return particles[i1].position.x < particles[i2].position.x; });
    // for (int i=0; i<particles.size(); i++) {
    //   if (order[i] != i) {
    //     printf("%d:%d ", i, order[i]);
    //   }
    // }
    // printf("\n");

    // put omp parallel single here, untied if (...) particles size threadshold 
    #pragma omp parallel
    {
      #pragma omp single
      {
        quadTree->root = buildQuadTree(particles, bmin, bmax);
      }
    }
    
    

    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }

  // Do not modify this function type.
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement parallel version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure
    // std::cout << "inside simulate step" << std::endl;
    auto quadTree = buildAccelerationStructure(particles);

    #pragma omp parallel for schedule(dynamic, 64) // static is not optimal
    for (int i =0; i < particles.size(); i++) {
      auto p = particles[i];
      std::vector<Particle> nearby;
      quadTree->getParticles(nearby, p.position, params.cullRadius);
      Vec2 force = Vec2(0.0f, 0.0f);
      #pragma omp parallel for schedule(auto) 
      for (auto &near: nearby) {
        force += computeForce(p, near, params.cullRadius);
      }
      newParticles[i] = updateParticle(p, force, params.deltaTime);
    }
  }
};

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelNBodySimulator() {
  return std::make_unique<ParallelNBodySimulator>();
}
