// src/denoise_obj.cpp

#include "denoise_obj.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>

// Vec3 实现
Vec3::Vec3() : x(0), y(0), z(0) {}

Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

Vec3 Vec3::operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z);
}

Vec3 Vec3::operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z);
}

Vec3 Vec3::operator*(float scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);
}

float Vec3::dot(const Vec3& other) const {
    return x * other.x + y * other.y + z * other.z;
}

Vec3 Vec3::cross(const Vec3& other) const {
    return Vec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

float Vec3::length() const {
    return std::sqrt(x * x + y * y + z * z);
}

Vec3 Vec3::normalize() const {
    float len = length();
    if(len == 0) return Vec3(0,0,0);
    return Vec3(x / len, y / len, z / len);
}

void Vec3::print() const {
    std::cout << "(" << x << ", " << y << ", " << z << ")";
}

// 计算高斯权重
float gaussian_weight(float r, float sigma) {
    return std::exp(- (r * r) / (2 * sigma * sigma));
}

// 计算面的质心
Vec3 compute_face_centroid(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    Vec3 centroid;
    int count = 0;
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it);
        centroid = centroid + Vec3(point[0], point[1], point[2]);
        count++;
    }
    if(count > 0)
        return centroid * (1.0f / count);
    else
        return Vec3(0, 0, 0);
}

// 计算面的真实法线
Vec3 compute_face_normal(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    // 获取当前面的所有顶点
    std::vector<Vec3> face_vertices;
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it);
        face_vertices.emplace_back(point[0], point[1], point[2]);
    }

    if (face_vertices.size() < 3) {
        return Vec3(0, 0, 0); // 面片顶点不足，法线为0向量
    }

    // 取前三个顶点计算法线（假设是三角形）
    Vec3 A = face_vertices[0];
    Vec3 B = face_vertices[1];
    Vec3 C = face_vertices[2];

    Vec3 AB = B - A;
    Vec3 AC = C - A;
    Vec3 normal = AB.cross(AC).normalize();

    return normal;
}

// 计算面的面积
float compute_face_area(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    // 获取当前面的所有顶点
    std::vector<Vec3> face_vertices;
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it);
        face_vertices.emplace_back(point[0], point[1], point[2]);
    }

    if (face_vertices.size() < 3) {
        return 0.0f; // 面片顶点不足，面积为0
    }

    // 取前三个顶点计算面积（假设是三角形）
    Vec3 A = face_vertices[0];
    Vec3 B = face_vertices[1];
    Vec3 C = face_vertices[2];

    Vec3 AB = B - A;
    Vec3 AC = C - A;
    Vec3 cross_product = AB.cross(AC);
    float area = 0.5f * cross_product.length();

    return area;
}

// 载入OBJ文件
bool load_obj(const std::string& file_path, MyMesh& mesh) {
    if (!OpenMesh::IO::read_mesh(mesh, file_path)) {
        std::cerr << "Error loading OBJ file: " << file_path << std::endl;
        return false;
    }
    return true;
}

// 保存OBJ文件
bool save_obj(const std::string& file_path, const MyMesh& mesh) {
    if (!OpenMesh::IO::write_mesh(mesh, file_path)) {
        std::cerr << "Error saving OBJ file: " << file_path << std::endl;
        return false;
    }
    return true;
}

// 顶点投影到切平面
// Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal) {
//     return vertex - normal.dot(vertex) * normal;
// }
// Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal) {
//     return vertex - normal * normal.dot(vertex);
// }
// 顶点投影到切平面
Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal, const Vec3& base_point) {
    // 计算向量(vertex - base_point)
    Vec3 vertex_to_base = vertex - base_point;
    // 计算投影
    return vertex - normal * normal.dot(vertex_to_base);
}


// 对网格进行去噪
std::vector<Vec3> smooth_mesh(const MyMesh& mesh, float sigma_f, float sigma_g) {
    std::vector<Vec3> smoothed_vertices;
    smoothed_vertices.reserve(mesh.n_vertices());

    // 计算每个面的质心和真实法线，并计算面积
    std::vector<std::pair<Vec3, Vec3>> face_centroids_normals;
    std::vector<float> face_areas; // 存储每个面的面积
    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        Vec3 centroid = compute_face_centroid(mesh, *f_it);
        Vec3 normal = compute_face_normal(mesh, *f_it); // 使用真实法线
        face_centroids_normals.emplace_back(centroid, normal);

        float area = compute_face_area(mesh, *f_it); // 计算面面积
        face_areas.push_back(area);
    }

    // 遍历每个顶点进行去噪
    size_t total_vertices = mesh.n_vertices();
    size_t current_vertex = 0;
    std::cout << "Smoothing vertices: " << std::endl;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++current_vertex) {
        const auto& point = mesh.point(*v_it);
        Vec3 vertex(point[0], point[1], point[2]);

        std::vector<float> f_weights;
        f_weights.reserve(face_centroids_normals.size());
        std::vector<float> g_weights;
        g_weights.reserve(face_centroids_normals.size());
        std::vector<Vec3> projections;
        projections.reserve(face_centroids_normals.size());

        std::vector<float> triangle_areas;
        triangle_areas.reserve(face_centroids_normals.size());

        // 计算每个面对当前顶点的权重和投影
        for (size_t i = 0; i < face_centroids_normals.size(); ++i) {
            const Vec3& centroid = face_centroids_normals[i].first;
            const Vec3& normal = face_centroids_normals[i].second;

            float distance = (vertex - centroid).length();
            float f_weight = gaussian_weight(distance, sigma_f);
            f_weights.push_back(f_weight);

            // 投影到真实法线的切平面
            Vec3 projection = project_to_tangent_plane(vertex, normal, centroid);
            projections.push_back(projection);

            float projection_dist = (projection - vertex).length();
            float g_weight = gaussian_weight(projection_dist, sigma_g);
            g_weights.push_back(g_weight);

            // 获取当前面的面积
            float area = face_areas[i];
            triangle_areas.push_back(area);
        }

        // 计算加权和，考虑面积
        float total_weight = 0.0f;
        Vec3 weighted_sum(0, 0, 0);
        for (size_t i = 0; i < projections.size(); ++i) {
            float weight = triangle_areas[i] * f_weights[i] * g_weights[i];
            weighted_sum = weighted_sum + (projections[i] * weight);
            total_weight += weight;
        }

        if(total_weight != 0.0f){
            smoothed_vertices.push_back(weighted_sum * (1.0f / total_weight));
        }
        else{
            smoothed_vertices.push_back(vertex); // 保持原样
        }

        // 显示进度
        if(current_vertex % 1000 == 0){
            std::cout << "\rProcessed " << current_vertex << " / " << total_vertices << " vertices." << std::flush;
        }
    }
    std::cout << "\rProcessed " << total_vertices << " / " << total_vertices << " vertices." << std::endl;

    return smoothed_vertices;
}

// 设置日志记录
void setup_logging(std::ofstream& log_file) {
    log_file.open("denoise_log.txt");
    if (!log_file.is_open()) {
        std::cerr << "Error opening log file." << std::endl;
        exit(1);
    }
}

// 记录顶点变化
void log_vertex_changes(std::ofstream& log_file, const Vec3& original_vertex, const Vec3& smoothed_vertex, size_t index) {
    log_file << "Vertex " << index << " - Original: (" << original_vertex.x << ", " << original_vertex.y << ", " << original_vertex.z 
             << ") | Smoothed: (" << smoothed_vertex.x << ", " << smoothed_vertex.y << ", " << smoothed_vertex.z << ")\n";
}

// 主去噪函数
void denoise_obj(const std::string& input_obj, const std::string& output_obj) {
    MyMesh mesh;

    // 载入OBJ文件
    if (!load_obj(input_obj, mesh)) {
        return;
    }

    // 去噪
    std::vector<Vec3> smoothed_vertices = smooth_mesh(mesh, 1.0f, 1.0f);

    // 记录日志
    std::ofstream log_file;
    setup_logging(log_file);
    size_t index = 0;
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++index) {
        const auto& original = mesh.point(*v_it);
        Vec3 original_vertex(original[0], original[1], original[2]);
        Vec3 smoothed_vertex = smoothed_vertices[index];
        log_vertex_changes(log_file, original_vertex, smoothed_vertex, index);
    }
    log_file.close();

    // 更新网格顶点
    index = 0;
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++index) {
        mesh.set_point(*v_it, OpenMesh::Vec3f(smoothed_vertices[index].x, smoothed_vertices[index].y, smoothed_vertices[index].z));
    }

    // 保存去噪后的OBJ文件
    if (!save_obj(output_obj, mesh)) {
        return;
    }

    std::cout << "Denoised OBJ saved to " << output_obj << std::endl;
}
