#include <cmath>
#include <limits>
#include <vector>

#include "geometry.h"
#include "model.h"
#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
Model *model = NULL;
const int width = 800;
const int height = 800;

inline int vec2idx(int x, int y) {
  return x + y * width;
}

inline Vec2i idxvec2(const int& id){
  return {id % width, id / width};
}

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
  bool steep = false;
  if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
    std::swap(p0.x, p0.y);
    std::swap(p1.x, p1.y);
    steep = true;
  }
  if (p0.x > p1.x) {
    std::swap(p0, p1);
  }

  for (int x = p0.x; x <= p1.x; x++) {
    float t = (x - p0.x) / (float)(p1.x - p0.x);
    int y = p0.y * (1. - t) + p1.y * t;
    if (steep) {
      image.set(y, x, color);
    } else {
      image.set(x, y, color);
    }
  }
}

/**
 * @brief compute barycentric coords for triangles
 */
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
  Vec3f s[2];
  for (int i = 2; i--;) {
    s[i][0] = C[i] - A[i];
    s[i][1] = B[i] - A[i];
    s[i][2] = A[i] - P[i];
  }
  Vec3f u = cross(s[0], s[1]);
  if (std::abs(u[2]) > 1e-2)  // dont forget that u[2] is integer. If it is zero
                              // then triangle ABC is degenerate
    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
  return Vec3f(-1, 1, 1);  // in this case generate negative coordinates, it
                           // will be thrown away by the rasterizator
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage & image,
                TGAColor color) {
    if (t0.y == t1.y && t0.y == t2.y)
      return;  // I dont care about degenerate triangles
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);

    // Rasteriaze
    float alpha = 0.0;
    float beta = 0.0;
    float da = (float) 1 / (t2.y - t0.y + 1);
    float db = (float) 1 / (t1.y - t0.y + 1);
    for(int y = t0.y; y <= t1.y; ++y) {
      Vec2i A = t0 + (t2 - t0) * alpha;
      Vec2i B = t0 + (t1 - t0) * beta;
      alpha += da, beta += db;
      if(A.x > B.x) std::swap(A.x, B.x);
      for(int j = A.x; j < B.x; ++ j){
        image.set(j, y, color);
      }
    }
    beta = 0.0;
    db = (float) 1 / (t2.y - t1.y + 1);
    for(int y = t1.y; y <= t2.y; ++y){
      Vec2i A = t0 + (t2 - t0) * alpha;
      Vec2i B = t1 + (t2 - t1) * beta;
      alpha += da; beta += db;
      if(A.x > B.x) std::swap(A.x, B.x);
      for(int j = A.x; j < B.x; ++ j){
        image.set(j, y, color);
      }
    }
}

void triangle(Vec3f *pts, Vec2i *uvs, float* zbuffer, TGAImage &image, float intensity) {
  Vec2f bboxmin(std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max());
  Vec2f bboxmax(-std::numeric_limits<float>::max(),
                -std::numeric_limits<float>::max());
  Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
      bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
    }
  }
  Vec3f P;
  Vec2i Puv;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      Vec3f bc_screen = barycentric(pts[0],pts[1], pts[2], P);
      if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
      P.z = 0;
      Puv.x = 0, Puv.y = 0;
      for (int i = 0; i < 3; ++i) {
        P.z += pts[i][2] * bc_screen[i];
        Puv = Puv + Vec2i(uvs[i][0] * bc_screen[i], uvs[i][1] * bc_screen[i]);
      }
      // std::cout<<uvs[1]<<std::endl;
      if (zbuffer[vec2idx(P.x, P.y)] < P.z){
        zbuffer[vec2idx(P.x, P.y)] = P.z;
        TGAColor texture = model->get_texture(Puv);
        texture.set_intensity(intensity);
        // std::cout<<Puv<<std::endl;
        image.set(P.x, P.y, texture);
      } 
    }
  }
}

/**
 * 等价于model矩阵转移
*/
inline Vec3f world2screen(Vec3f v) {
  return Vec3f(int((v.x + 1.) * width / 2. + .5),
               int((v.y + 1.) * height / 2. + .5), v.z);
}

int main(int argc, char **argv) {
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("obj/african_head/african_head.obj");
  }

  float *zbuffer = new float[width * height];
  for (int i = width * height; i--;
       zbuffer[i] = -std::numeric_limits<float>::max())
    ;
  TGAImage image(width, height, TGAImage::RGB);
  Vec3f light_dir(0, 0, -1);
  for (int i = 0; i < model->nfaces(); i++) {
    std::vector<Vec3i> face = model->face(i);
    Vec3f world_coords[3];
    Vec3f pts[3];
    Vec2i uvs[3];
    for (int i = 0; i < 3; i++){
      Vec3f v = model->vert(face[i].x);
      pts[i] = world2screen(v);
      world_coords[i] = v;
      uvs[i] = model->tcoords(face[i].y);
      // std::cout<<uvs[i] << std::endl;
    }
    // 叉乘生成法向量
    Vec3f n = cross((world_coords[2] - world_coords[0]),
              (world_coords[1] - world_coords[0]));
    n.normalize();
    float intensity = n * light_dir;  // 夹角越小，强度越高
    if (intensity > 0.) {              // 判断是否遮挡
      triangle(pts, uvs, zbuffer, image, intensity);
    }
  }

  image.flip_vertically();  // i want to have the origin at the left bottom
                            // corner of the image
  image.write_tga_file("output.tga");
  delete model;
  return 0;
}