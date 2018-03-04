#ifndef SCENE_TRIANGULATION_H
#define SCENE_TRIANGULATION_H

#include <VVRScene/scene.h>
#include <VVRScene/canvas.h>
#include <VVRScene/utils.h>
#include <GeoLib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <ctime>

struct Tri;

class TriangulationScene : public vvr::Scene
{
public:
    TriangulationScene();

    const char* getName() const override {
        return "UNIVERSITY OF PATRAS - VVR GROUP - COMPUTATIONAL GEOMETRY LAB";
    }

protected:
    void draw() override;
    void reset() override;
    void resize() override;
    void keyEvent(unsigned char key, bool up, int modif) override;
    void mousePressed(int x, int y, int modif) override;
    void mouseMoved(int x, int y, int modif) override;
    void mouseWheel(int dir, int modif) override;

private:
    void processPoint(C2DPoint* const p);

private:
    vvr::Canvas2D m_canvas;
    C2DPointSet m_pts;
    std::vector<Tri> m_triangles;
    int m_style_flag;
    float m_lw_canvas;
    float m_lw_tris;
    float m_sz_pt;
	vvr::Canvas2D m_canvas_1;
};

/**
* Struct representing a triangle with pointers to its 3 vertices
*/
struct Tri
{
    C2DPoint *v1;
    C2DPoint *v2;
    C2DPoint *v3;
    float area;

    Tri(C2DPoint *v1, C2DPoint *v2, C2DPoint *v3) : v1(v1), v2(v2), v3(v3) {
        area = to_C2D().GetArea();
    }

    C2DTriangle to_C2D() const { return C2DTriangle(*v1, *v2, *v3); }

    vvr::Triangle2D to_vvr(vvr::Colour col = vvr::Colour::black, bool filled = false) const {
        vvr::Triangle2D t(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y, col);
        t.setSolidRender(filled);
        return t;
    }

    bool operator < (const Tri& other) const {
        if (area != other.area) return (area < other.area);
        else if (v1 != other.v1) return v1 < other.v1;
        else if (v2 != other.v2) return v2 < other.v2;
        else if (v3 != other.v3) return v3 < other.v3;
    }
};

C2DCircle GetCircumCircle(const C2DTriangle &t);
C2DPoint isOnSeg(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3);
bool IsDelaunay(const C2DTriangle &t, const C2DPointSet &pset);
bool TriIsDelaunay(std::vector<Tri> &m_triangles, int i, const C2DPointSet &m_pts);
bool PtIsInCircle(const C2DTriangle &t, const C2DPoint *pt);
bool isIntersect(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3, C2DPoint i);

bool FindAdjofTri(std::vector<Tri> &tris,int num_of_tri, C2DPoint *p1, C2DPoint *p2, unsigned &tri_adj_index, C2DPoint* &opposite_vertex);

bool FindAdjacentTriangle(std::vector<Tri> &tris, C2DPoint *p1, C2DPoint *p2, unsigned &tri_adj_index, C2DPoint* &opposite_vertex);

void FindViolations(std::vector<Tri> &tris, const C2DPointSet &ptset, std::vector<unsigned> &violations);

void ShowViolations(std::vector<Tri> &tris, const std::vector<unsigned> &violations, vvr::Canvas2D &canvas, const vvr::Colour &col);

#endif // SCENE_TRIANGULATION_H
