#ifndef SCENE_CONVEX_HULL_H
#define SCENE_CONVEX_HULL_H

#include <VVRScene/scene.h>
#include <VVRScene/canvas.h>
#include <VVRScene/utils.h>
#include <GeoLib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>

class ConvexHullScene : public vvr::Scene
{
public:
	ConvexHullScene();
	const char* getName() const override {
		return "UNIVERSITY OF PATRAS - "
			"VVR GROUP - COMPUTATIONAL GEOMETRY LAB";
	}

protected:
	void draw() override;
	void reset() override;
	void createRandomPts();
	void mousePressed(int x, int y, int modif) override;
	void mouseMoved(int x, int y, int modif) override;
	void mouseReleased(int x, int y, int modif) override;

private:
	vvr::Canvas2D m_canvas;
	vvr::Canvas2D m_canvas_1;
	vvr::Canvas2D m_canvas_2;
	C2DPointSet m_ptset;
	C2DPolygon m_convex_hull_polygon;
	C2DPolygon m_magic_polygon;
	vvr::Point2D* m_mouse_pos;
	vvr::Point2D p1;
	vvr::Point2D p2;
	vvr::Line2D line2;
	C2DPointSet m_setaki;
	C2DLine grammh;

};

//int isOnRight(vvr::Point2D p1, vvr::Point2D p2, vvr::Point2D m_mouse_pos);
int isOnRight(vvr::Line2D line2, vvr::Point2D m_mouse_pos);
int isOnLeftC2D(C2DLine line, C2DPoint t);

C2DPolygon ConvexHull_Fast(C2DPointSet &ptset);
C2DPolygon ConvexHull_Slow(C2DPointSet &ptset);

#endif // SCENE_CONVEX_HULL_H