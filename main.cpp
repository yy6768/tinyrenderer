#include <iostream>
#include <vector>
#include <limits>

#include "geometry.h"
#include "model.h"
#include "our_gl.h"
#include "tgaimage.h"

float *shadowbuffer = NULL;
Model *model = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir(1, 1, 1);
Vec3f eye(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

struct GouraudShader : public IShader {
  Vec3f varying_intensity;  // written by vertex shader, read by fragment shader

  virtual Vec4f vertex(int iface, int nthvert) {
    varying_intensity[nthvert] =
        std::max(0.f, model->normal(iface, nthvert) *
                          light_dir);  // get diffuse lighting intensity
    Vec4f gl_Vertex = embed<4>(
        model->vert(iface, nthvert));  // read the vertex from .obj file
    return Viewport * Projection * ModelView *
           gl_Vertex;  // transform it to screen coordinates
  }

  virtual bool fragment(Vec3f bar, TGAColor &color) {
    float intensity =
        varying_intensity * bar;  // interpolate intensity for the current pixel
    color = TGAColor(255, 255, 255) * intensity;  // well duh
    return false;  // no, we do not discard this pixel
  }
};

struct TextureShader : public IShader {
  Vec3f varying_intensity;
  mat<2, 3, float> varying_uv;

  virtual Vec4f vertex(int iface, int nthvert) {
    varying_uv.set_col(nthvert, model->uv(iface, nthvert));
    varying_intensity[nthvert] = std::max (0.f, model->normal(iface, nthvert) * light_dir);
    Vec4f gl_vertex = embed<4>(model->vert(iface, nthvert));
    return Viewport * Projection * ModelView * gl_vertex;
  }

  virtual bool fragment(Vec3f bar, TGAColor &color) {
    float intensity = varying_intensity * bar;
    Vec2f uv = varying_uv * bar;
    color = model->diffuse(uv) * intensity;
    return false;
  }
};

struct DarbouxShader : public IShader {
  mat<2, 3, float> varying_uv;
  mat<4, 4, float> uniform_matrix;
  mat<4, 4, float> uniform_invert_transpose;

  virtual Vec4f vertex(int iface, int nthvert) {
    varying_uv.set_col(nthvert, model->uv(iface, nthvert)); 
    Vec4f gl_vertex = embed<4>(model->vert(iface, nthvert));
    return Viewport * Projection * ModelView * gl_vertex;
  }

  virtual bool fragment(Vec3f bar, TGAColor& color) {
    Vec2f uv = varying_uv * bar;
    Vec3f n = proj<3>(uniform_invert_transpose * embed<4>(model->normal(uv))).normalize();
    Vec3f l = proj<3>(uniform_matrix * embed<4>(light_dir)).normalize();
    float intensity = std::max(0.f, n * l);
    color = model->diffuse(uv) * intensity;
    return false; 
  } 
};

struct Shader : public IShader {
  mat<4, 4, float> uniform_M;    //  Projection*ModelView
  mat<4, 4, float> uniform_MIT;  // (Projection*ModelView).invert_transpose()
  mat<4, 4, float> uniform_Mshadow;  // transform framebuffer screen coordinates
                                     // to shadowbuffer screen coordinates
  mat<2, 3, float> varying_uv;       // triangle uv coordinates, written by the
                                // vertex shader, read by the fragment shader
  mat<3, 3, float> varying_tri;  // triangle coordinates before Viewport
                                 // transform, written by VS, read by FS

  Shader(Matrix M, Matrix MIT, Matrix MS)
      : uniform_M(M),
        uniform_MIT(MIT),
        uniform_Mshadow(MS),
        varying_uv(),
        varying_tri() {}

  virtual Vec4f vertex(int iface, int nthvert) {
    varying_uv.set_col(nthvert, model->uv(iface, nthvert));
    Vec4f gl_Vertex = Viewport * Projection * ModelView *
                      embed<4>(model->vert(iface, nthvert));
    varying_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
    return gl_Vertex;
  }

  virtual bool fragment(Vec3f bar, TGAColor &color) {
    Vec4f sb_p = uniform_Mshadow *
                 embed<4>(varying_tri *
                          bar);  // corresponding point in the shadow buffer
    sb_p = sb_p / sb_p[3];
    int idx =
        int(sb_p[0]) + int(sb_p[1]) * width;  // index in the shadowbuffer array
    float shadow =
        .3 + .7 * (shadowbuffer[idx] <
                   sb_p[2] + 43.34);     // magic coeff to avoid z-fighting
    Vec2f uv = varying_uv * bar;         // interpolate uv for the current pixel
    Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv)))
                  .normalize();  // normal
    Vec3f l =
        proj<3>(uniform_M * embed<4>(light_dir)).normalize();  // light vector
    Vec3f r = (n * (n * l * 2.f) - l).normalize();  // reflected light
    float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
    float diff = std::max(0.f, n * l);
    TGAColor c = model->diffuse(uv);
    for (int i = 0; i < 3; i++)
      color[i] =
          std::min<float>(20 + c[i] * shadow * (1.2 * diff + .6 * spec), 255);
    return false;
  }
};

struct PhongShader : public IShader {
    mat<2,3,float> varying_uv;  
    mat<4, 4, float> uniform_matrix;
    mat<4, 4, float> uniform_invert_transpose;

    virtual Vec4f vertex(int iface, int nthvert) {
      varying_uv.set_col(nthvert, model->uv(iface, nthvert));
      Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); 
      return Viewport*Projection*ModelView*gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        Vec2f uv = varying_uv*bar;
        Vec3f n = proj<3>(uniform_invert_transpose *embed<4>(model->normal(uv))).normalize();
        Vec3f l = proj<3>(uniform_matrix  *embed<4>(light_dir)).normalize();
        Vec3f r = (n*(n*l*2.f) - l).normalize();  
        float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
        float diff = std::max(0.f, n*l);
        TGAColor c = model->diffuse(uv);
        color = c;
        for (int i=0; i<3; i++) color[i] = std::min<float>(5 + c[i]*(diff + .6*spec), 255);
        return false;
  }
};

struct DepthShader : public IShader{
  mat<3, 3, float> varying_tri;

  const float depth = 255.0; // 透视深度
  DepthShader() : varying_tri() {}

  virtual Vec4f vertex(int iface, int nthvert) {
    Vec4f gl_vertex = Viewport * ModelView * embed<4>(model->vert(iface, nthvert));    
    varying_tri.set_col(nthvert, proj<3>(gl_vertex / gl_vertex[3]));
    return gl_vertex;
  }

  virtual bool fragment(Vec3f bar, TGAColor& color) {
    Vec3f p = varying_tri * bar;
    color = TGAColor(255, 255, 255) * (p.z/depth);
    return false;
  }
};

int main(int argc, char **argv) {
  if (2 > argc) {
    std::cerr << "Usage: " << argv[0] << "obj/model.obj" << std::endl;
    return 1;
  }

  float *zbuffer = new float[width * height];
  shadowbuffer = new float[width * height];
  for (int i = width * height; --i;) {
    zbuffer[i] = shadowbuffer[i] = -std::numeric_limits<float>::max();
  }

  model = new Model(argv[1]);
  light_dir.normalize();

  {  // rendering the shadow buffer
    TGAImage depth(width, height, TGAImage::RGB);
    lookat(light_dir, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(0);

    DepthShader depthshader;
    Vec4f screen_coords[3];
    for (int i = 0; i < model->nfaces(); i++) {
      for (int j = 0; j < 3; j++) {
        screen_coords[j] = depthshader.vertex(i, j);
      }
      triangle(screen_coords, depthshader, depth, shadowbuffer);
    }
    depth.flip_vertically();  // to place the origin in the bottom left corner
                              // of the image
    depth.write_tga_file("depth.tga");
  }

  Matrix M = Viewport * Projection * ModelView;

  {  // rendering the frame buffer
    TGAImage frame(width, height, TGAImage::RGB);
    lookat(eye, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(-1.f / (eye - center).norm());

    Shader shader(ModelView, (Projection * ModelView).invert_transpose(),
                  M * (Viewport * Projection * ModelView).invert());
    Vec4f screen_coords[3];
    for (int i = 0; i < model->nfaces(); i++) {
      for (int j = 0; j < 3; j++) {
        screen_coords[j] = shader.vertex(i, j);
      }
      triangle(screen_coords, shader, frame, zbuffer);
    }
    frame.flip_vertically();  // to place the origin in the bottom left corner
                              // of the image
    frame.write_tga_file("framebuffer.tga");
  }

  delete model;
  delete[] zbuffer;
  delete[] shadowbuffer;

  return 0;
}