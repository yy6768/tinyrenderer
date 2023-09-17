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
const int depth = 255;

Vec3f eye(50., 10., 50.);
Vec3f center(0, 0, 0);

Vec3f m2v(Matrix m) {
  return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
}

Matrix v2m(Vec3f v) {
  Matrix m(4, 1);
  m[0][0] = v.x;
  m[1][0] = v.y;
  m[2][0] = v.z;
  m[3][0] = 1.f;
  return m;
}

Matrix viewport(float x, float y, float w, float h) {
  Matrix m = Matrix::identity(4);
  m[0][3] = x + w / 2.f;
  m[1][3] = y + h / 2.f;
  m[2][3] = depth / 2.f;

  m[0][0] = w / 2.f;
  m[1][1] = h / 2.f;
  m[2][2] = depth / 2.f;
  return m;
}

/**
 * @brief 
 * @param eye the direction which camera points to
 * @param center the center of object
 * @param up defined by camera
 * @return the lookat Matrix
*/
Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
  Vec3f z = (eye - center).normalize(); // 指向摄像机朝向-z,所以z轴方向是摄像机反方向
  Vec3f x = cross(up, z).normalize(); // 指向摄影机右边
  Vec3f y = cross(z, x).normalize(); // 摄影上方
  Matrix Minv = Matrix::identity(4); // 基变换矩阵的逆矩阵
  Matrix Tr = Matrix::identity(4); // 平移矩阵
  for (int i = 0; i < 3;  ++i){
    Minv[0][i] = x[i];
    Minv[1][i] = y[i];
    Minv[2][i] = z[i];
    Tr[i][3] = -eye[i];
  }
  return Minv * Tr; 
}

inline int vec2idx(int x, int y) {
  return x + y * width;
}

inline Vec2i idxvec2(const int& id){
  return {id % width, id / width};
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
  if (std::abs(u[2]) > 0.)  // dont forget that u[2] is integer. If it is zero
                              // then triangle ABC is degenerate
    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
  return Vec3f(-1, 1, 1);  // in this case generate negative coordinates, it
                           // will be thrown away by the rasterizator
}



void triangle(Vec3f *pts, Vec2i *uvs, int* zbuffer, TGAImage &image, float intensity) {
  Vec2f bboxmin(std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max());
  Vec2f bboxmax(-std::numeric_limits<float>::max(),
                -std::numeric_limits<float>::max());
  Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
  for (int i = 0; i < 3; i++) {
    // std::cout << pts[i] << std::endl;
      bboxmin.x = std::max(0.f, std::min(bboxmin.x, pts[i].x));
      bboxmin.y = std::max(0.f, std::min(bboxmin.y, pts[i].y));
      bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
      bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
  }
  // std::cout << bboxmax << std::endl;
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
        Puv = Puv + Vec2i(uvs[i].x * bc_screen[i], uvs[i].y * bc_screen[i]);
      }
      if (zbuffer[vec2idx(P.x, P.y)] < P.z){
        zbuffer[vec2idx(P.x, P.y)] = P.z;
        TGAColor texture = model->get_texture(Puv);
        texture.set_intensity(intensity);
        image.set(P.x, P.y, texture);
      } 
    }
  }
}

int main(int argc, char **argv) {
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("obj/african_head/african_head.obj");
  }

  int *zbuffer = new int[width * height];
  for (int i = width * height; i--;
       zbuffer[i] = -std::numeric_limits<int>::max())
       ;
  

  TGAImage image(width, height, TGAImage::RGB);
  Vec3f light_dir(0, 0, -1);
  Matrix Projection = Matrix::identity(4);
  Matrix ViewPort   = viewport(width/8.0, height/8.0, width*3./4., height*3./4.);
  Matrix ModelView  = lookat(eye, center, Vec3f(0, 1, 0));
  Projection[3][2] = -1.f/(eye - center).norm();
  for (int i = 0; i < model->nfaces(); i++) {
    std::vector<Vec3i> face = model->face(i);
    Vec3f world_coords[3];
    Vec3f pts[3];
    Vec2i uvs[3];
    
    for (int i = 0; i < 3; i++){
      Vec3f v = model->vert(face[i].x);
      // pts[i] = world2screen(v);
      pts[i] = m2v(ViewPort * Projection * ModelView * v2m(v));  
      pts[i].x = int(pts[i].x + 0.5);
      pts[i].y = int(pts[i].y + 0.5);
      world_coords[i] = v;
      uvs[i] = model->tcoords(face[i].y);
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