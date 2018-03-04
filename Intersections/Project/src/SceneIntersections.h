#ifndef SCENE_INTERSECTIONS_H
#define SCENE_INTERSECTIONS_H

#include <VVRScene/scene.h>
#include <VVRScene/canvas.h>
#include <VVRScene/utils.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <GeoLib.h>

class SceneIntersections : public vvr::Scene
{
public:
    SceneIntersections();

    const char* getName() const override {
        return "UNIVERSITY OF PATRAS - VVR GROUP - COMPUTATIONAL GEOMETRY LAB";
    }

protected: // Overriden methods
    void draw() override;
    void reset() override;
    void mousePressed(int x, int y, int modif) override;
    void mouseMoved(int x, int y, int modif) override;

private: // Methods
    void Task1(const C2DPoint &p);
    void Task2(const C2DPoint &p);
    void Task3(const C2DPoint &p);

private: // Data
    C2DLine         m_bound_horizontal;
    C2DLine         m_bound_vertical;
    vvr::Canvas2D   m_canvas_0;
    vvr::Canvas2D   m_canvas_1;
    vvr::Canvas2D   m_canvas_2;
    vvr::Canvas2D   m_canvas_3;
    C2DTriangle     m_triangle_1;
    C2DTriangle     m_triangle_2;
    C2DCircle       m_circle_1;
    C2DCircle       m_circle_2;
    C2DLine         m_line_1;
    C2DLine         m_line_2;
	C2DLine m_trig_1;
	C2DLine m_trig_2;
	C2DPolygon m_mesa;
	C2DPolygon m_tri;
	
};

bool isIntersect(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3, C2DPoint i);
C2DPoint isOnSeg(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3);
int isOnLeftC2D(C2DLine line, C2DPoint t);

#endif // SCENE_INTERSECTIONS_H
