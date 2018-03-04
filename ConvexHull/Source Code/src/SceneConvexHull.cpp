#include "SceneConvexHull.h"

#define NUM_PTS 50

ConvexHullScene::ConvexHullScene()
{
	m_bg_col = vvr::Colour(0x44, 0x44, 0x44);
	m_hide_sliders = true;
	m_hide_log = true;
	p1 = vvr::Point2D(30, 30);
	p2 = vvr::Point2D(100, 130);
	line2 = vvr::Line2D(p1.x, p1.y, p2.x, p2.y, vvr::Colour::blue);
	C2DBaseSet m_setaki;

	reset();

	//...
}

void ConvexHullScene::reset()
{
	Scene::reset(); // Call base class reset
	m_canvas.clear();
	createRandomPts();
	// Compute convex hull.
	std::cout << "Computing Convex Hull..." << std::endl;
	float duration = vvr::getSeconds();

	m_convex_hull_polygon = ConvexHull_Slow(m_ptset);

	m_setaki.DeleteAll();
	// ta vazw se BaseSet gia na ginontai delete kathe fora pou kanw reset
	// kanw copy ola ta points tou polygwnou se auto 
	m_convex_hull_polygon.GetPointsCopy(m_setaki);

	m_canvas.clear();
	m_canvas_1.clear();
	m_canvas_2.clear();

	duration = vvr::getSeconds() - duration;
	std::cout << "Duration: " << duration << " sec" << std::endl;
}

void ConvexHullScene::createRandomPts()
{
	const int BW = 100, BH = 100;
	C2DRect bound(-BW, BH, BW, -BH);
	C2DPolygon rand_poly;
	m_ptset.DeleteAll();
	rand_poly.CreateRandom(bound, NUM_PTS, NUM_PTS);
	rand_poly.GetPointsCopy(m_ptset);
}

void ConvexHullScene::mousePressed(int x, int y, int modif)
{
	m_canvas.clear();
	m_canvas_1.clear();
	m_canvas_2.clear();

	m_mouse_pos = new vvr::Point2D;
	m_mouse_pos->x = x; //diavazei se pointer
	m_mouse_pos->y = y;
	m_canvas_2.add(m_mouse_pos);
	const C2DPoint p(x, y);// diavazei se C2D
	// dhmiourgw topiko LineSet gia na valw ola ta euth tmhmata pou vlepoun to p
	C2DLineSet m_magic; 
	C2DPoint n1;
	C2DPoint n2;
	vvr::Point2D mouse = *m_mouse_pos; //kanei to pointer se vvr
	
	C2DPointSet seten; //twra molis
	//Vazw to BaseSet se PointSet gia na to epeksergastw pio eukola
	seten.AddCopy(m_setaki);
	int n = seten.size();
	// ta bools gia na elegxw tis 2 periptwseis sto 2.2 erwthma
	bool pame = false;
	bool pame1 = false;
	bool pame2 = false;
	echo(n);
	int sum = 0;
	int mesa_ola;// gia na vrw to size tou LineSet sthn deuterh periptwsh
	for (int i = 0; i < seten.size() - 1; i++) {
		n1 = seten[i];
		n2 = seten[i + 1];
		float shmeio1 = (n2.x - n1.x)*(mouse.y - n1.y) - (n2.y - n1.y)*(mouse.x - n1.x);
		if (pame2 == true) {
			if (pame == true) {
			mesa_ola = m_magic.size();
			}
		}
		pame2 = true;
		if (shmeio1 > 0) {
			sum++;
			grammh = C2DLine(n1, n2);
			m_magic.AddCopy(grammh);
			if (i == 0) {
				pame = true;
			}
			pame2 = false;

		}
	}
	//elegxw to teleutaio me to prwto
	n1 = seten[n - 1]; 
	n2 = seten[0];
	float shmeio1 = (n2.x - n1.x)*(mouse.y - n1.y) - (n2.y - n1.y)*(mouse.x - n1.x);
	if (shmeio1 > 0) {
		sum++;
		grammh = C2DLine(n1, n2);
		m_magic.AddCopy(grammh);
		pame1 = true;

	}
	// an einai mesa kai to prwto kai to teleutaio pare ta shmeia pou spane kai ftiakse tis grammes
	//1h periptwsh
	if (pame == true && pame1 == true){
			n1 = m_magic[mesa_ola].GetPointFrom();
			n2 = m_magic[mesa_ola-1].GetPointTo();
			grammh = C2DLine(n1, p);
			m_canvas.add(grammh, vvr::Colour::red);
			grammh = C2DLine(n2, p);
			m_canvas.add(grammh, vvr::Colour::red);

	}
	else if ( sum >0){
		//2h periptwsh
		int magic_size = m_magic.size();
		n1 = m_magic[0].GetPointFrom();
		n2 = m_magic[magic_size - 1].GetPointTo();
		grammh = C2DLine(n1, p);
		m_canvas.add(grammh, vvr::Colour::red);
		grammh = C2DLine(p, n2);
		m_canvas.add(grammh, vvr::Colour::red);
	}
	echo(sum);
	if (sum == 0) {
		m_mouse_pos->colour = vvr::Colour::red;

	}
	else {
		m_mouse_pos->colour = vvr::Colour::cyan;


	}
	echo(m_magic.size());

	for (int i=0; i < m_magic.size(); i++) {
		C2DPoint n1 = m_magic[i].GetPointFrom();
		C2DPoint n2 = m_magic[i].GetPointTo();
		grammh = C2DLine(n1, n2);
		m_canvas.add(grammh, vvr::Colour::red);
	}


}

void ConvexHullScene::mouseMoved(int x, int y, int modif)
{
	m_canvas.clear();
	m_canvas_1.clear();
	m_canvas_2.clear();
	m_mouse_pos = new vvr::Point2D;
	m_mouse_pos->x = x; //diavazei se pointer
	m_mouse_pos->y = y;
	int result = isOnRight(line2, *m_mouse_pos);

	if (result > 0) {
	m_mouse_pos->colour = vvr::Colour::green;
	}
	else {
	m_mouse_pos->colour = vvr::Colour::red;
	}
	m_canvas_1.add(m_mouse_pos);
	mousePressed(x, y, modif);

}

void ConvexHullScene::mouseReleased(int x, int y, int modif)
{

}

void ConvexHullScene::draw()
{

	enterPixelMode();

	vvr::draw(m_convex_hull_polygon, vvr::Colour::green);
	vvr::draw(m_ptset, vvr::Colour::black);
	
	m_canvas.draw();
	//m_canvas_1.draw();
	m_canvas_2.draw();

	
	//p1.draw();
	//p2.draw();
	//line2.draw();

	//vvr::LineSeg2D line(p1.x, p1.y, p2.x, p2.y, vvr::Colour::red);
	//line.draw();
	returnFromPixelMode();
}

C2DPolygon ConvexHull_Slow(C2DPointSet &ptset)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Breite to kyrto polygwno twn simeinwn `ptset`.
	//!
	//! HINTS:
	//!
	//!  - An `p` einai C2DPoint* tote gia na parete to `i` simeio tou point cloud
	//!    grafete to eksis:
	//!      p = ptset.GetAt(i);
	//!
	//!  - Dimiourgia euthigrammou tmimatos apo 2 simeia:
	//!      C2DLine line(p1,  p2);
	//!
	//!  - To idio, alla me pointers sta simeia:
	//!      C2DLine line(*p1, *p2);
	//!
	//!  - Stin metavliti `lineset` prepei na apothikeysete ta eythigramma tmimata
	//!    pou vriskete oti anikoun sto Convex Hull.
	//!
	//!  - An `pq` einai ena C2DLine, etsi to prosthetete sto lineset:
	//!      lineset.AddCopy(pq);
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	const int num_pts = ptset.size(); // To plithos twn simeiwn.
	C2DLineSet lineset;
	C2DPoint *p, *q, *t;
	int valid;
	for (int pi = 0; pi < num_pts; pi++)
	{
		p = ptset.GetAt(pi);
		for (int qi = 0; qi < num_pts; qi++)
		{
			q = ptset.GetAt(qi);
			if (p != q) {
				C2DLine line(*p, *q);
				valid = 1;
				for (int ti = 0; ti < num_pts; ti++)
				{
					t = ptset.GetAt(ti);
					if (t != q) {
						if (t != p) {

							int result = isOnLeftC2D(line, *t);
							if (result > 0) {
								valid = -1;
							}
						}
					}
				}
				if (valid > 0) {
					lineset.AddCopy(line);
				}
			}
		}
	}

	// LineSet --> Polygon --> Convex Hull
	C2DPoint *pts = new C2DPoint[lineset.size() * 2];
	for (int i = 0; i < lineset.size(); i++) {
		pts[i] = lineset.GetAt(i)->GetPointFrom();

	}

	C2DPolygon convex_hull_polygon;
	convex_hull_polygon.Create(pts, lineset.size(), true);
	convex_hull_polygon.RemoveNullLines();
	delete[] pts;
	return convex_hull_polygon;
}

C2DPolygon ConvexHull_Fast(C2DPointSet &ptset)
{
	// Copy the points to a temporary array in order to
	// be compatible with GeoLib's functions.
	C2DPoint *pts = new C2DPoint[ptset.size()];
	for (int i = 0; i < ptset.size(); i++) {
		pts[i] = *ptset.GetAt(i);
	}

	// Create convex hull
	C2DPolygon cloud_polygon, convex_hull_polygon;
	cloud_polygon.Create(pts, ptset.size());
	convex_hull_polygon.CreateConvexHull(cloud_polygon);
	delete[] pts;
	return convex_hull_polygon;
}

//int isOnRight(vvr::Point2D p1, vvr::Point2D p2, vvr::Point2D m_mouse_pos) {
//float shmeio = m_mouse_pos.y - p1.y - ((p1.y - p2.y) / (p1.x - p2.x))*(m_mouse_pos.x - p1.x);
//
//if (shmeio > 0) {
//  return 1;
//}
//  else {
//  return -1;
//  }
//
//}
int isOnLeftC2D(C2DLine line, C2DPoint t) {
	C2DPoint end = line.GetPointTo();
	C2DPoint start = line.GetPointFrom();
	//float shmeio1 = t.y - start.y - ((start.y - end.y) / (start.x - end.x))*(t.x - start.x);
	float shmeio1 = (end.x - start.x)*(t.y - start.y) - (end.y - start.y)*(t.x - start.x);

	if (shmeio1 > 0) {
		return 1; //aristera
	}
	else {
		return -1; // deksia
	}

}

int isOnRight(vvr::Line2D line2, vvr::Point2D m_mouse_pos) {
	float shmeio = m_mouse_pos.y - line2.y1 - ((line2.y1 - line2.y2) / (line2.x1 - line2.x2))*(m_mouse_pos.x - line2.x1);

	if (shmeio > 0) {
		return 1;
	}
	else {
		return -1;
	}

}

int main(int argc, char* argv[])
{
	return vvr::mainLoop(argc, argv, new ConvexHullScene);
}