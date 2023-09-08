#include "model.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

Model::Model(const char *filename) : verts_(), faces_() {
  std::ifstream in;
  in.open(filename, std::ifstream::in);
  if (in.fail()) return;
  std::string line;
  while (!in.eof()) {
    std::getline(in, line);
    std::istringstream iss(line.c_str());
    char trash;
    if (!line.compare(0, 2, "v ")) {
      iss >> trash;
      Vec3f v;
      for (int i = 0; i < 3; i++) iss >> v[i];
      verts_.push_back(v);
    } else if (!line.compare(0, 2, "f ")) {
      std::vector<Vec3i> f;
      int vidx, tidx, nidx;
      iss >> trash;
      while (iss >> vidx >> trash >> tidx >> trash >> nidx) {
        vidx--, tidx--, nidx--;  // in wavefront obj all indices start at 1, not zero
        f.emplace_back(vidx, tidx, nidx);
        
      }
      faces_.push_back(f);
    } else if (!line.compare(0, 3, "vt ")) {
      float u, v;
      iss >> trash >> trash >> u >> v;
      tcoords_.emplace_back(u, v);
    }
  }
  std::cerr << "# v# " << verts_.size() << " f# " << faces_.size() << "vt# " << tcoords_.size() << std::endl;
  load_texture(filename, "_diffuse.tga", texture_);
}

Model::~Model() {}

int Model::nverts() { return (int)verts_.size(); }

int Model::nfaces() { return (int)faces_.size(); }

int Model::ntcoords() { return (int)tcoords_.size(); }

std::vector<Vec3i> Model::face(int idx) { return faces_[idx]; }

Vec3f Model::vert(int i) { return verts_[i]; }

TGAColor Model::get_texture(Vec2i uv) { return texture_.get(uv.x, uv.y); }

Vec2i Model::tcoords(int i) { 

  return Vec2i(tcoords_[i].x * texture_.get_width(), tcoords_[i].y * texture_.get_height()); 
}

void Model::load_texture(std::string filename, const char* suffix, TGAImage& image) {
  std::string texfile(filename);
  size_t dot = texfile.find_last_of(".");
  if (dot != std::string::npos) {
    texfile = texfile.substr(0, dot) + std::string(suffix);
    std::cout << "texture loading: " << texfile<< std::endl;
    if (image.read_tga_file(texfile.c_str())) {
      std::cout<< "Ok!" << std::endl;
    } else {
      std::cout<< "Loading failed: not exists";
    }
    image.flip_vertically();
    std::cout<<image.get_width() << " "<< image.get_height()<< std::endl;
  }
}