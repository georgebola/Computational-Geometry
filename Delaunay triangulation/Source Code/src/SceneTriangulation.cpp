#include "SceneTriangulation.h"
#include <VVRScene/utils.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <ctime>

#define FLAG_SHOW_CANVAS 1

using namespace std;
using namespace vvr;

TriangulationScene::TriangulationScene()
{
	m_bg_col = Colour(0x44, 0x44, 0x44);
	m_create_menus = false;
	m_lw_canvas = 2;
	m_lw_tris = 3;
	m_sz_pt = 10;
}

void TriangulationScene::reset()
{
	Scene::reset();
	m_style_flag = 0;
	m_style_flag |= FLAG_SHOW_CANVAS;
	m_canvas.clear();
	m_canvas_1.clear();
	m_triangles.clear();
	m_pts.DeleteAll();

	const int W = getViewportWidth() * 0.9;
	const int H = getViewportHeight()  * 0.9;

	C2DPoint* A = new C2DPoint(-W / 2, -H / 2);
	C2DPoint* B = new C2DPoint(W / 2, -H / 2);
	C2DPoint* C = new C2DPoint(0, H / 2);

	m_pts.Add(A);
	m_pts.Add(B);
	m_pts.Add(C);
	m_triangles.push_back(Tri(A, B, C));
}

void TriangulationScene::resize()
{
	static bool first_pass = true;
	if (first_pass) reset();
	first_pass = false;
}

void TriangulationScene::mousePressed(int x, int y, int modif)
{
	Scene::mousePressed(x, y, modif);
	processPoint(new C2DPoint(x, y));
}

void TriangulationScene::mouseMoved(int x, int y, int modif)
{
	Scene::mouseMoved(x, y, modif);
	if (m_pts.GetLast()->Distance(C2DPoint(x, y)) > 80)
		processPoint(new C2DPoint(x, y));
}

void TriangulationScene::keyEvent(unsigned char key, bool up, int modif)
{
	Scene::keyEvent(key, up, modif);
	const bool ctrl_down = ctrlDown(modif);
	const bool alt_down = altDown(modif);
	const bool shift_down = shiftDown(modif);
	key = tolower(key);

	switch (key)
	{
	case 'c': m_style_flag ^= FLAG_SHOW_CANVAS; break;
	}
}

void TriangulationScene::mouseWheel(int dir, int modif)
{
	if (shiftDown(modif))
	{
		m_lw_canvas *= pow(1.2f, dir);
		m_lw_tris *= pow(1.2f, dir);
		m_sz_pt *= pow(1.2f, dir);
	}
	else
	{
		Scene::mouseWheel(dir, modif);
	}
}

void TriangulationScene::draw()
{
	enterPixelMode();
	
	//! Draw violations and anything else added to canvas
	if (m_style_flag & FLAG_SHOW_CANVAS) {
		Shape::DEF_LINE_WIDTH = m_lw_canvas;
		m_canvas.draw();
	}

	//! Draw triangles
	Shape::DEF_LINE_WIDTH = m_lw_tris;

	vector<Tri>::const_iterator tri_iter;
	for (tri_iter = m_triangles.begin(); tri_iter != m_triangles.end(); ++tri_iter) {
		tri_iter->to_vvr(Colour::black).draw();
	}

	//! Draw points
	Shape::DEF_POINT_SIZE = m_sz_pt;
	vvr::draw(m_pts, Colour::red);
	m_canvas_1.draw();
	returnFromPixelMode();
}

int main(int argc, char* argv[])
{
	try
	{
		return vvr::mainLoop(argc, argv, new TriangulationScene);
	}
	catch (std::string exc)
	{
		cerr << exc << endl;
		return 1;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//! Task Methods / Functions
/////////////////////////////////////////////////////////////////////////////////////////

void TriangulationScene::processPoint(C2DPoint* const p)
{
	//! [Check 1]
	//! Check whether a point already exists in the same coords.

	for (size_t i = 0; i < m_pts.size(); i++)
	{
		if (m_pts.GetAt(i)->x == p->x &&
			m_pts.GetAt(i)->y == p->y)
		{
			delete p;
			return;
		}
	}

	//! Find enclosing triangle.
	unsigned i_enclosing;
	unsigned count_enclosing = 0;
	for (int i = 0; i < m_triangles.size(); i++) {
		if (m_triangles[i].to_C2D().Contains(*p)) {
			count_enclosing++;
			i_enclosing = i;
		}
	}
	//! [Check 2]
	//! If no enclosing triangle was found.
	//! Or if more than one were found.

	if (count_enclosing != 1) {
		delete p;
		return;
	}

	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!  - Create the 3 subdivision triangles.
	//!  - Check the 3 new triangles for Delaunay violations.
	//!  - Vreite to geitoniko (adjacent) trigwno tou tri_check.
	//!  - Kante to flip.
	//!
	//! HINTS:
	//!  - Remove the enclosing triangle
	//!
	//!  - Store the 3 new triangles which will replace old one.
	//!    Store them *AFTER* the flipping step, in order to avoid finding one of these
	//!    as adjacent triangle.
	//!
	//!  - To delete a triangle from `m_triangles` :
	//!      m_triangles.erase(m_triangles.begin() + index_trigwnou_pros_diagrafh);
	//!
	//!  - To add a filled `Tri` object to the canvas:
	//!      m_canvas.add(new Triangle2D(tri.to_vvr(color, true)));
	//!
	//!  - Kathe trigwnno exei eite 3 geitonika, eite 2 an prokeitai gia trigwno
	//!    pou vrisketai sta oria (=> anikei merikws sto Convex Hull)
	//!
	//!  - bool adj_exists = FindAdjacentTriangle(m_triangles, edge_start, edge_end, &tri_adjacent, &v_opposite);
	//!
	//! VARIABLES:
	//!
	//!  - C2DPoint* p;                // New point
	//!  - Tri &tri_enclosing          // Enclosing triangle of the new point
	//!  - Tri tris_new[3];            // The 3 triangles that will replace the enclosing.
	//!  - unsigned tri_adjacent_index // The adjacent triangle index, i.e. m_triangles[tri_adjacent_index] is what we want
	//!  - C2DPoint *v_opposite        // The opposite vertex of the adjacent triangle
	//!
	//!//////////////////////////////////////////////////////////////////////////////////

	vector<Tri> new_triangles;
	Tri &tri_enclosing = m_triangles[i_enclosing];
	new_triangles.push_back(Tri(p, tri_enclosing.v1, tri_enclosing.v2));
	new_triangles.push_back(Tri(p, tri_enclosing.v2, tri_enclosing.v3));
	new_triangles.push_back(Tri(p, tri_enclosing.v3, tri_enclosing.v1));

	//! [Check 3]
	//! Check if any of the 3 triangles are colinear. (Use GeoLib's function)

	if (new_triangles[0].to_C2D().Collinear() ||
		new_triangles[1].to_C2D().Collinear() ||
		new_triangles[2].to_C2D().Collinear())
	{
		delete p;
		return;
	}

	//! HERE: We have a valid point, and we can proceed

	m_canvas.clear();
	m_canvas_1.clear();
	m_pts.Add(p);
	m_triangles.erase(m_triangles.begin() + i_enclosing);

	for (int i = 0; i < 3; i++)
	{
		bool did_flip = false;

		if (!IsDelaunay(new_triangles[i].to_C2D(), m_pts))
		{
			C2DPoint *p1 = new_triangles[i].v1;
			C2DPoint *p2 = new_triangles[i].v2;
			C2DPoint *p3 = new_triangles[i].v3;

			C2DPoint *v_opposite = NULL;
			unsigned tri_adjacent_index = NULL;
			bool adj_exists = FindAdjacentTriangle(m_triangles, p2, p3, tri_adjacent_index, v_opposite);
			
			if (adj_exists) // Adjacent triangle exists
			{
				// prepei to euthigrammo tmhma v_opposite-p1 na temnei to p2,p3
				C2DPoint v_opp, po1, po2, po3;
				v_opp.x = v_opposite->x;
				v_opp.y = v_opposite->y;
				po1 = new_triangles[i].to_C2D().GetPoint1();
				po2 = new_triangles[i].to_C2D().GetPoint2();
				po3 = new_triangles[i].to_C2D().GetPoint3();
				C2DPoint Seg = isOnSeg(v_opp, po1, po2, po3);
				bool vlepei = isIntersect(v_opp, po1, po2, po3, Seg);
				if (vlepei) {
					m_triangles.push_back(Tri(p3, v_opposite, p1));
					m_triangles.push_back(Tri(p2, v_opposite, p1));
					//m_canvas_1.add(Tri(p3, v_opposite, p1).to_C2D(), vvr::Colour::red);
					//m_canvas_1.add(Tri(p2, v_opposite, p1).to_C2D(), vvr::Colour::red);
					m_triangles.erase(m_triangles.begin() + tri_adjacent_index);
				}
				else // No adjacent found with a straight view of i triangle
				{
					m_triangles.push_back(new_triangles[i]);
				}
			}
			else // No adjacent found, i.e. side on boundary
			{
				m_triangles.push_back(new_triangles[i]);
			}
		}
		else // triangle is delaunay, so we add it
		{
			m_triangles.push_back(new_triangles[i]);
		}

	}

	bool vrhke;
	// trexoume mia apo katw pros ta panw kai mia apo panw pros ta katw
	//giati mporei na petuxei kapoio pou den ginetai na allaksei kai na mhn mporei na allaksei meta kai ta ypoloipa
	// otan petuxei vevaia kapoio pou den allazei ta elegxei ola k fores
	//giauto meiwnontas to k trexei pio grhgora eidika otan exoun mazeutei polla trigwna
	for (int k = 0; k < 5; k++) {
		vrhke = false;
		int cnt=0;
		for (int i = 0; i< m_triangles.size() ; i++) {
			vrhke = TriIsDelaunay(m_triangles, i, m_pts);
			if (!vrhke) {
				cnt++;
			}
		}
		if (cnt == m_triangles.size()) {
			break;			
		}
		for (int i = m_triangles.size()-1; i > -1; i--) {
			vrhke = TriIsDelaunay(m_triangles, i, m_pts);
			if (!vrhke) {
				
				cnt++;
			}
		}
		if (cnt == m_triangles.size()) {
			break;
		}
			
	}

	vector<unsigned> violations;
	FindViolations(m_triangles, m_pts, violations);
	ShowViolations(m_triangles, violations, m_canvas, Colour::magenta);
	C2DPoint pos;
	pos.x = p->x;
	pos.y = p->y;
	//m_canvas_1.add(pos, vvr::Colour::green);

}

C2DCircle GetCircumCircle(const C2DTriangle &t)
{
	C2DCircle circle;
	C2DPoint p1;
	C2DPoint p2;
	C2DPoint p3;


	p1 = t.GetPoint1();
	p2 = t.GetPoint2();
	p3 = t.GetPoint3();

	float l1 = (p2.y - p1.y) / (p2.x - p1.x);

	float meso_1_x = (p1.x + p2.x) / 2;
	float meso_1_y = (p1.y + p2.y) / 2;
	C2DPoint meso_point(meso_1_x, meso_1_y);

	float l1_k = -(1 / l1); //lamda mesokathetou 1
	float deutero_y = l1_k*(p1.x + 1 - meso_point.x) + meso_point.y; //xnew = x+1
	C2DPoint deutero_shmeio(p1.x + 1, deutero_y);

	float l2 = (p2.y - p3.y) / (p2.x - p3.x);

	float meso_2_x = (p3.x + p2.x) / 2;
	float meso_2_y = (p3.y + p2.y) / 2;
	C2DPoint meso_point2(meso_2_x, meso_2_y);

	float l2_k = -(1 / l2); //lamda mesokathetou 2
	float deutero_y2 = l2_k*(p3.x + 1 - meso_point2.x) + meso_point2.y; //xnew = x+1
	C2DPoint deutero_shmeio2(p3.x + 1, deutero_y2);

	C2DPoint tomh = isOnSeg(meso_point, deutero_shmeio, meso_point2, deutero_shmeio2);
	//aktina= apostash tomhs apo p1
	float r = sqrt(pow((p1.x - tomh.x), 2) + pow((p1.y - tomh.y), 2));

	circle.Set(tomh, r);
	return circle;
}

bool PtIsInCircle(const C2DTriangle &t, const C2DPoint *pt)
{
	//an to shmeio periexete ston kyklo epistrefei true
	C2DCircle c = GetCircumCircle(t);
	if (c.Contains(*pt)) {
		return true;
	}
	return false;
}
bool IsDelaunay(const C2DTriangle &t, const C2DPointSet &pset)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Check if this triangle is a Delaunay triangle.
	//!
	//! HINTS:
	//!
	//!  - Shrink the circle a bit in order to exclude points of its circumference.
	//!
	//!//////////////////////////////////////////////////////////////////////////////////

	for (int i = 0; i < pset.size(); i++) {
		C2DCircle c = GetCircumCircle(t);
		c.SetRadius(c.GetRadius()*0.9);
		if (c.Contains(*pset.GetAt(i))) {
			return false;
		}
	}

	return true;
}

bool TriIsDelaunay(vector<Tri> &m_triangles, int i, const C2DPointSet &m_pts) {
	Tri &tri = m_triangles[i];
	C2DTriangle t(*tri.v1, *tri.v2, *tri.v3);

	if (!IsDelaunay(t, m_pts)) {
		C2DPoint *v_opposite1 = NULL;
		C2DPoint *v_opposite2 = NULL;
		C2DPoint *v_opposite3 = NULL;

		unsigned tri_adjacent_index1 = NULL;
		unsigned tri_adjacent_index2 = NULL;
		unsigned tri_adjacent_index3 = NULL;

		//vriskw ta stoixeia tou trigwnou ths kathe pleuras pou akoumpaei
		//to trigwno pou den einai delanay
		bool adj_exists1 = FindAdjofTri(m_triangles, i, m_triangles[i].v1, m_triangles[i].v2, tri_adjacent_index1, v_opposite1);
		bool adj_exists2 = FindAdjofTri(m_triangles, i, m_triangles[i].v1, m_triangles[i].v3, tri_adjacent_index2, v_opposite2);
		bool adj_exists3 = FindAdjofTri(m_triangles, i, m_triangles[i].v2, m_triangles[i].v3, tri_adjacent_index3, v_opposite3);
		bool Pt_In_Circle1 = false;
		bool Pt_In_Circle2 = false;
		bool Pt_In_Circle3 = false;
		//C2DPoint opp;
		if (adj_exists1) {
			Pt_In_Circle1 = PtIsInCircle(m_triangles[i].to_C2D(), v_opposite1);
		}
		if (adj_exists2) {
			Pt_In_Circle2 = PtIsInCircle(m_triangles[i].to_C2D(), v_opposite2);
		}
		if (adj_exists3) {
			Pt_In_Circle3 = PtIsInCircle(m_triangles[i].to_C2D(), v_opposite3);
		}
		if (adj_exists1 == true & Pt_In_Circle1 == true) // Adjacent triangle exists and his opposite point is inside the circle
		{
			m_triangles.push_back(Tri(m_triangles[i].v1, v_opposite1, m_triangles[i].v3));
			m_triangles.push_back(Tri(m_triangles[i].v2, v_opposite1, m_triangles[i].v3));
			//m_canvas_1.add(t.GetPoint1(), vvr::Colour::green);
			//m_canvas_1.add(t.GetPoint2(), vvr::Colour::green);
			if (i > tri_adjacent_index1) {
				m_triangles.erase(m_triangles.begin() + i);
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index1);
			}
			else {
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index1);
				m_triangles.erase(m_triangles.begin() + i);
			}


		}
		else if (adj_exists2 == true & Pt_In_Circle2 == true) // Adjacent triangle exists and his opposite point is inside the circle
		{
			m_triangles.push_back(Tri(m_triangles[i].v1, v_opposite2, m_triangles[i].v2));
			m_triangles.push_back(Tri(m_triangles[i].v3, v_opposite2, m_triangles[i].v2));
			if (i > tri_adjacent_index2) {
				m_triangles.erase(m_triangles.begin() + i);
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index2);

			}
			else {
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index2);
				m_triangles.erase(m_triangles.begin() + i);
			}

		}
		else if (adj_exists3 == true & Pt_In_Circle3 == true) // Adjacent triangle exists and his opposite point is inside the circle
		{
			m_triangles.push_back(Tri(m_triangles[i].v2, v_opposite3, m_triangles[i].v1));
			m_triangles.push_back(Tri(m_triangles[i].v3, v_opposite3, m_triangles[i].v1));
			//m_canvas_1.add(t.GetPoint2(), vvr::Colour::green);
			//m_canvas_1.add(t.GetPoint3(), vvr::Colour::green);
			if (i > tri_adjacent_index3) {
				m_triangles.erase(m_triangles.begin() + i);
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index3);
			}
			else {
				m_triangles.erase(m_triangles.begin() + tri_adjacent_index3);
				m_triangles.erase(m_triangles.begin() + i);

			}

		}
		return true;
	}
	else {
		return false;
	}

}
void FindViolations(vector<Tri> &tris, const C2DPointSet &ptset, vector<unsigned> &violations)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Check in `tris` for Delaunay violations.
	//!
	//! HINTS:
	//!
	//!  - If triangle i causes violation add it to violations vector like this:
	//!      violations.push_back(i);
	//!
	//!//////////////////////////////////////////////////////////////////////////////////

	violations.clear();

	for (int i = 0; i < tris.size(); i++)
	{
		Tri &tri = tris[i];
		C2DTriangle t(*tri.v1, *tri.v2, *tri.v3);
		if (!IsDelaunay(t, ptset)) {
			violations.push_back(i);
		}
	}
}

void ShowViolations(vector<Tri> &tris, const vector<unsigned> &violations, Canvas2D &canvas, const Colour &col)
{
	for (int i = 0; i < violations.size(); i++) {
		Tri &tri = tris[violations[i]];
		C2DTriangle t(*tri.v1, *tri.v2, *tri.v3);
		canvas.add(GetCircumCircle(t), col, false);
	}
}

bool FindAdjacentTriangle(vector<Tri> &tris, C2DPoint *p1, C2DPoint *p2, unsigned &tri_adj_index, C2DPoint* &opp_ver)
{
	for (int i = 0; i < tris.size(); i++)
	{
		Tri &tri = tris[i];
		C2DPoint v1 = *tri.v1;
		C2DPoint v2 = *tri.v2;
		C2DPoint v3 = *tri.v3;

		if ((*p1 == v1 && *p2 == v2) || (*p1 == v2 && *p2 == v1) ){
			opp_ver = tris[i].v3;
			tri_adj_index = i;
			return true;
		}
		if ((*p1 == v1 && *p2 == v3) || (*p1 == v3 && *p2 == v1)) {
			tri_adj_index = i;

			opp_ver = tris[i].v2;
			return true;
		}
		if ((*p1 == v2 && *p2 == v3) || (*p1 == v3 && *p2 == v2)) {
			tri_adj_index = i;

			opp_ver = tris[i].v1;
			return true;
		}

	}

	return false;
}
bool FindAdjofTri(vector<Tri> &tris,int num_of_tri, C2DPoint *p1, C2DPoint *p2, unsigned &tri_adj_index, C2DPoint* &opp_ver)
{
	for (int i = 0; i < tris.size(); i++)
	{
		if (i != num_of_tri) {
			Tri &tri = tris[i];
			C2DPoint v1 = *tri.v1;
			C2DPoint v2 = *tri.v2;
			C2DPoint v3 = *tri.v3;

			if ((*p1 == v1 && *p2 == v2) || (*p1 == v2 && *p2 == v1)) {
				opp_ver = tris[i].v3;
				tri_adj_index = i;
				return true;
			}
			if ((*p1 == v1 && *p2 == v3) || (*p1 == v3 && *p2 == v1)) {
				tri_adj_index = i;

				opp_ver = tris[i].v2;
				return true;
			}
			if ((*p1 == v2 && *p2 == v3) || (*p1 == v3 && *p2 == v2)) {
				tri_adj_index = i;

				opp_ver = tris[i].v1;
				return true;
			}
		}
	}

	return false;
}
C2DPoint isOnSeg(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3) {
	C2DPoint i;
	// eutheia etoimh
	float a_1 = (p2.y - p3.y) / (p2.x - p3.x);
	float b_1 = -(a_1*p3.x) + p3.y;
	// eutheia me to xeri
	float a_2 = (p.y - p1.y) / (p.x - p1.x);
	float b_2 = -(a_2*p1.x) + p1.y;
	i.x = (b_2 - b_1) / (a_1 - a_2);
	i.y = a_2*i.x + b_2;

	return i;
}
bool isIntersect(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3, C2DPoint i) {
	bool seg_intersect = false;

	if (i.x > std::min(p3.x, p2.x) && i.x < std::max(p3.x, p2.x)) {
		if (i.x > std::min(p1.x, p.x) && i.x < std::max(p1.x, p.x)) {
			seg_intersect = true;
		}
	}

	return seg_intersect;
}

//////
