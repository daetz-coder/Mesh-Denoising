// include/denoise_obj.hpp

#ifndef DENOISE_OBJ_HPP
#define DENOISE_OBJ_HPP

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>

// 使用 OpenMesh 内置的三角网格类
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

// 定义一个3D向量
struct Vec3 {
    float x, y, z;

    Vec3();
    Vec3(float x, float y, float z);

    Vec3 operator+(const Vec3& other) const;
    Vec3 operator-(const Vec3& other) const;
    Vec3 operator*(float scalar) const;
    
    float dot(const Vec3& other) const;
    Vec3 cross(const Vec3& other) const;
    
    float length() const;
    Vec3 normalize() const;

    void print() const;
};

// 计算高斯权重
float gaussian_weight(float r, float sigma);

// 计算面的质心
Vec3 compute_face_centroid(const MyMesh& mesh, const MyMesh::FaceHandle& face);

// 计算面的真实法线
Vec3 compute_face_normal(const MyMesh& mesh, const MyMesh::FaceHandle& face);

// 载入OBJ文件
bool load_obj(const std::string& file_path, MyMesh& mesh);

// 保存OBJ文件
bool save_obj(const std::string& file_path, const MyMesh& mesh);

// 顶点投影到切平面
// Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal);
// 顶点投影到切平面
Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal, const Vec3& base_point);

// 对网格进行去噪
std::vector<Vec3> smooth_mesh(const MyMesh& mesh, float sigma_f = 1.0f, float sigma_g = 1.0f);

// 设置日志记录
void setup_logging(std::ofstream& log_file);

// 记录顶点变化
void log_vertex_changes(std::ofstream& log_file, const Vec3& original_vertex, const Vec3& smoothed_vertex, size_t index);

// 主去噪函数
void denoise_obj(const std::string& input_obj, const std::string& output_obj);

#endif // DENOISE_OBJ_HPP
