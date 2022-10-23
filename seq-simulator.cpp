#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>

// TASK 1

// NOTE: You may modify any of the contents of this file, but preserve all
// function types and names. You may add new functions if you believe they will
// be helpful.

const int QuadTreeLeafSize = 8;
class SequentialNBodySimulator : public INBodySimulator {
public:
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {
    // TODO: implement a function that builds and returns a quadtree containing
    // particles.

    //create a root node pointer
    auto root = new QuadTreeNode();
    // check if cur particles are able to fit one tree node
    if (particles.size() < QuadTreeLeafSize){
      root -> particles = particles;
      root -> isLeaf = true;
      return std::unique_ptr<QuadTreeNode>(root);
    }

    // cannot fit --> not leaf node  -> split
    Vec2 middle = Vec2((bmin.x+bmax.x)/2,(bmin.y+bmax.y)/2);

    // split into 4 regions
    //create particles lists for each partition --> check if need further split
  
    std::vector<Particle> child_0;
    std::vector<Particle> child_1;
    std::vector<Particle> child_2;
    std::vector<Particle> child_3;

    // assign each particle to specific region/partition

    // four child nodes are stored in following order:
    //  x0, y0 --------------- x1, y0
    //    |           |           |
    //    |children[0]|children[1]|
    //    | ----------+---------  |
    //    |children[2]|children[3]|
    //    |           |           |
    //  x0, y1 ----------------- x1, y1
    // where x0 < x1 and y0 < y1.

    for (int i = 0; i < (int)particles.size(); i++){
      
      Particle p = particles[i];
      Vec2 p_pos = p.position;
      // Debug boundary cases cannot pass --> particle can locate on boundary --> will miss count 
      /*
      if (p_pos.x > bmin.x && p_pos.x < middle.x && p_pos.y > bmin.y && p_pos.y < middle.y)
      {
        child_0.push_back(p);
      }
      else if (p_pos.x > middle.x && p_pos.x < bmax.x && p_pos.y > bmin.y && p_pos.y < middle.y)
      {
        child_1.push_back(p);
      }
      else if (p_pos.x > bmin.x && p_pos.x < middle.x && p_pos.y > middle.y && p_pos.y < bmax.y)
      {
        child_2.push_back(p);
      }
      else if (p_pos.x < bmax.x && p_pos.x > middle.x && p_pos.y > middle.y && p_pos.y < bmax.y)
      {
        child_3.push_back(p);
      }
    }
    */
    if (p_pos.x < middle.x && p_pos.y < middle.y)
      {
        child_0.push_back(p);
      }
      else if (p_pos.x >= middle.x && p_pos.y < middle.y)
      {
        child_1.push_back(p);
      }
      else if (p_pos.x < middle.x && p_pos.y >= middle.y)
      {
        child_2.push_back(p);
      }
      else if (p_pos.x >= middle.x && p_pos.y >= middle.y )
      {
        child_3.push_back(p);
      }
    }










    Vec2 up_middle = Vec2((bmin.x+bmax.x)/2, bmin.y);
    Vec2 left_middle = Vec2(bmin.x, (bmin.y+bmax.y)/2);
    Vec2 right_middle = Vec2(bmax.x, (bmin.y+bmax.y)/2);
    Vec2 bottom_middle = Vec2((bmin.x+bmax.x)/2, bmax.y);


    // recursive steps for dividing into 4 regions
    root ->children[0] = buildQuadTree(child_0, bmin, middle);
    root ->children[1] = buildQuadTree(child_1, up_middle, right_middle);
    root ->children[2] = buildQuadTree(child_2, left_middle, bottom_middle);
    root ->children[3] = buildQuadTree(child_3, middle, bmax);

    // since this cannot fit a leaf node can it should be splitted
    root -> isLeaf = false;

    return std::unique_ptr<QuadTreeNode>(root);
  }



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
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement sequential version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure

    // first build quadTree  
    std::unique_ptr<AccelerationStructure> quadTree = buildAccelerationStructure(particles);


    // iterate all particles and check its pos in the quad tree built
    for (int i = 0; i < (int)particles.size(); i++) {
      auto pi = particles[i];
      // create a particle list to recording particles within radius
      std::vector<Particle> neighbors;
      quadTree->getParticles(neighbors, pi.position, params.cullRadius);


      Vec2 force = Vec2(0.0f, 0.0f);
      // accumulate attractive forces in neighbors to apply to particle i
      for (auto &neighbor: neighbors) {
        force += computeForce(pi, neighbor, params.cullRadius);
      }
      // update particle state using the computed force
      newParticles[i] = updateParticle(pi, force, params.deltaTime);
    }

  }
};

std::unique_ptr<INBodySimulator> createSequentialNBodySimulator() {
  return std::make_unique<SequentialNBodySimulator>();
}


