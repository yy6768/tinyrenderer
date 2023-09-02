#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>

#include "geometry.h"
/**
 * 模型类
 */
class Model {
 private:
  std::vector<Vec3f> verts_;              // 定点集合
  std::vector<std::vector<int> > faces_;  // 面集合

 public:
  Model(const char *filename);
  ~Model();
  int nverts();
  int nfaces();
  Vec3f vert(int i);
  std::vector<int> face(int idx);
};

#endif  //__MODEL_H__