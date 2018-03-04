#include <VVRScene/canvas.h>
#include <VVRScene/mesh.h>
#include <VVRScene/settings.h>
#include <VVRScene/utils.h>
#include <MathGeoLib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <set>
#include "symmetriceigensolver3x3.h"

#define FLAG_SHOW_AXES       1
#define FLAG_SHOW_WIRE       2
#define FLAG_SHOW_SOLID      4
#define FLAG_SHOW_NORMALS    8
#define FLAG_SHOW_PLANE     16
#define FLAG_SHOW_AABB      32
#define INTERSECTION_POINTS 64

void Task_1_FindCenterMass(std::vector<vec> &vertices, vec &cm);
void Task_2_FindAABB(std::vector<vec> &vertices, vvr::Box3D &aabb);
void Task_3_AlignOriginTo(std::vector<vec> &vertices, const vec &cm);
void Task_4_Draw_PCA(vec &center, vec &dir);
void Task_5_Intersect(std::vector<vvr::Triangle>& triangles, Plane &plane, std::vector<int> &intersection_indices, std::vector<vec> &m_i_points);
void Task_5_Split(vvr::Mesh &mesh, Plane &plane, std::vector<vec> &m_i_points, std::vector<vvr::Triangle>& new_triangles, vvr::Mesh &m_greenmesh);
void pca(std::vector<vec>& vertices, vec &center, vec &dir);

void erwthma_2(std::vector<vvr::Triangle>& triangles, std::vector<int> &m_intersections, std::vector<vec> &m_i_points, Plane &plane);
void just_draw(std::vector<vvr::Triangle>& triangles, std::vector<int> &m_intersections);
class Mesh3DScene : public vvr::Scene
{
public:
    Mesh3DScene();
    const char* getName() const { return "3D Scene"; }
    void keyEvent(unsigned char key, bool up, int modif) override;
    void arrowEvent(vvr::ArrowDir dir, int modif) override;

private:
    void draw() override;
    void reset() override;
    void resize() override;
    void Tasks();
    vec Task_3_Pick_Origin();

private:
    int m_style_flag;
    float m_plane_d;
    vvr::Canvas2D m_canvas;
    vvr::Colour m_obj_col;
    vvr::Mesh m_model_original, m_model,m_greenmesh;
    vvr::Box3D m_aabb;
    math::vec m_center_mass;
    math::vec m_pca_cen;
    math::vec m_pca_dir;
    math::Plane m_plane;
    std::vector<int> m_intersections;
	std::vector<vec> m_i_points;
	std::vector<vvr::Triangle> new_triangles;
};
