#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>

#include "geometry.h"
#include "tgaimage.h"
/**
 * 模型类
 */
class Model {
 private:
  std::vector<Vec3f> verts_;              // 定点集合
  std::vector<std::vector<Vec3i> > faces_;  // 面集合
  std::vector<Vec2f> tcoords_;            // 纹理坐标集合
  TGAImage texture_;                      // 纹理
  void load_texture(std::string, const char*, TGAImage&);
 public:
  Model(const char *filename);
  ~Model();
  int nverts();
  int nfaces();
  int ntcoords();
  Vec3f vert(int i);
  std::vector<Vec3i> face(int idx);
  Vec2i tcoords(int i);
  TGAColor get_texture(Vec2i uv);
};

#endif  //__MODEL_H__