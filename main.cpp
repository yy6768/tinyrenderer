#include <cmath>
#include <vector>

#include "geometry.h"
#include "model.h"
#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = NULL;
const int width = 800;
const int height = 800;

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
Vec3f barycentric(Vec2i *pts, Vec2i P) {
  Vec3f u = Vec3f(pts[2].x - pts[0].x, pts[1].x - pts[0].x, pts[0].x - P.x) 
          ^ Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y);
  if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);
  return {1.f - (u.x +u.y) / u.z, u.y / u.z , u.x / u.z};
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

void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
  Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
  Vec2i bboxmax(0, 0);
  Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
  for (int i = 0; i < 3; i++) {
    bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
    bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

    bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
    bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
  }
  Vec2i P;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      Vec3f bc_screen = barycentric(pts, P);
      if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
      image.set(P.x, P.y, color);
    }
  }
}

int main(int argc, char **argv) {
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("obj/african_head.obj");
  }

  TGAImage image(width, height, TGAImage::RGB);
  Vec3f light_dir(0, 0, -1);
  for (int i = 0; i < model->nfaces(); i++) {
    std::vector<int> face = model->face(i);
    Vec2i screen_coords[3];
    Vec3f world_coords[3];
    for (int j = 0; j < 3; j++) {
      Vec3f v = model->vert(face[j]);
      screen_coords[j] =
          Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
      world_coords[j] = v;
    }
    // 叉乘生成法向量
    Vec3f n = (world_coords[2] - world_coords[0]) ^
              (world_coords[1] - world_coords[0]);
    n.normalize();
    float intensity = n * light_dir; // 夹角越小，强度越高
    if (intensity > 0) {
      triangle(
          screen_coords[0], screen_coords[1], screen_coords[2], image,
          TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    }
  }

  image.flip_vertically();  // i want to have the origin at the left bottom
                            // corner of the image
  image.write_tga_file("output.tga");
  delete model;
  return 0;
}