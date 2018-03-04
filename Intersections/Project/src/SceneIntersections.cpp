#include "SceneIntersections.h"

using namespace std;
using namespace vvr;

/* Construct - Load  - Setup */

SceneIntersections::SceneIntersections()
{
	m_bg_col = Colour(0x44, 0x44, 0x44);
	m_hide_log = false;
	reset();

}

void SceneIntersections::reset()
{
	Scene::reset();

	// Clear everything
	m_canvas_0.clear();
	m_canvas_1.clear();
	m_canvas_2.clear();
	m_canvas_3.clear();

	// Divide window to Tasks
	m_bound_vertical.Set(C2DPoint(0, -3000), C2DPoint(0, 3000));
	m_bound_horizontal.Set(C2DPoint(4000, 0), C2DPoint(-4000, 0));
	m_canvas_0.add(m_bound_horizontal, Colour::black);
	m_canvas_0.add(m_bound_vertical, Colour::black);

	// Setup Task 1:
	{
		C2DPoint a1(-300, 100);
		C2DPoint a2(-100, 200);
		C2DPoint b1(-350, 230);
		C2DPoint b2(-50, 50);
		m_line_1 = C2DLine(a1, a2);
		m_line_2 = C2DLine(b1, b2);
		m_canvas_0.add(a1, Colour::orange);
		m_canvas_0.add(a2, Colour::orange);
		m_canvas_0.add(m_line_1, Colour::orange);
		m_canvas_1.add(b1, Colour::cyan);
		m_canvas_1.add(b2, Colour::cyan);
		m_canvas_1.add(m_line_2, Colour::cyan);
	}

	// Task 2:
	{
		C2DPoint t1a(-300, -50);
		C2DPoint t1b(-40, -45);
		C2DPoint t1c(-70, -170);
		m_triangle_1 = C2DTriangle(t1a, t1b, t1c);
		C2DPoint t2a(-197, -266);
		C2DPoint t2b(-368, -136);
		C2DPoint t2c(-108, -76);
		m_triangle_2 = C2DTriangle(t2a, t2b, t2c);
		m_canvas_0.add(m_triangle_1, Colour::orange);
		m_canvas_2.add(m_triangle_2, Colour::cyan);
	}

	// Setup Task 3:
	{
		C2DPoint c1(166, 112);
		C2DPoint c2(290, 150);
		m_circle_1 = C2DCircle(c1, 80);
		m_circle_2 = C2DCircle(c2, 60);
		m_canvas_0.add(c1, Colour::orange);
		m_canvas_0.add(m_circle_1, Colour::orange);
		m_canvas_3.add(c2, Colour::cyan);
		m_canvas_3.add(m_circle_2, Colour::cyan);
	}

	Task1(m_line_2.GetPointTo());
	Task2(m_triangle_2.GetPoint3());
	Task3(m_circle_2.GetCentre());
	
}

/* UI Handling */

void SceneIntersections::mousePressed(int x, int y, int modif)
{
	Scene::mousePressed(x, y, modif);
	echo(x);
	echo(y);
	const C2DPoint p(x, y);
	if (m_bound_horizontal.IsOnRight(p) && !m_bound_vertical.IsOnRight(p)) Task1(p);
	if (m_bound_horizontal.IsOnRight(p) && m_bound_vertical.IsOnRight(p)) Task3(p);
	if (!m_bound_horizontal.IsOnRight(p) && !m_bound_vertical.IsOnRight(p)) Task2(p);
}

void SceneIntersections::mouseMoved(int x, int y, int modif)
{
	Scene::mouseMoved(x, y, modif);
	mousePressed(x, y, modif);
}

/* Tasks */

void SceneIntersections::Task1(const C2DPoint &p)
{
	C2DPoint p1 = m_line_2.GetPointFrom();   // To arxiko simeio paremenei idio.
	m_line_2 = C2DLine(p1, p);   // To teliko simeio tis grammis akolouthei to mouse.


	m_canvas_1.clear();
	m_canvas_1.add(p, Colour::cyan);
	m_canvas_1.add(p1, Colour::cyan);
	m_canvas_1.add(m_line_2, Colour::cyan);


	C2DPoint p2 = m_line_1.GetPointFrom();
	C2DPoint p3 = m_line_1.GetPointTo();


	C2DPoint i;

	i = isOnSeg(p1, p, p2,p3);
	bool seg_intersect = isIntersect(p1, p, p2, p3, i);

	m_canvas_1.add(i, seg_intersect ? Colour::green : Colour::red);
}

void SceneIntersections::Task2(const C2DPoint &p)
{
	const C2DPoint &p1 = m_triangle_2.GetPoint1();
	const C2DPoint &p2 = m_triangle_2.GetPoint2();
	//stathero
	const C2DPoint &p4 = m_triangle_1.GetPoint1();
	const C2DPoint &p5 = m_triangle_1.GetPoint2();
	const C2DPoint &p6 = m_triangle_1.GetPoint3();

	m_triangle_2.Set(p1, p2, p);
	m_canvas_2.clear();
	m_canvas_2.add(m_triangle_2, Colour::cyan);

	C2DPoint i1, i2, i3, i4, i5, i6;

	i1 = isOnSeg(p1, p, p4, p5);
	bool seg_intersect1 = isIntersect(p1, p, p4, p5, i1);
	if (seg_intersect1 == true) {
		m_canvas_2.add(i1,vvr::Colour::red);
	}

	i2 = isOnSeg(p1, p, p5, p6);
	bool seg_intersect2 = isIntersect(p1, p, p5, p6, i2);
	if (seg_intersect2 == true) {
		m_canvas_2.add(i2, vvr::Colour::red);
	}

	i3 = isOnSeg(p1, p, p4, p6);
	bool seg_intersect3 = isIntersect(p1, p, p4, p6, i3);
	if (seg_intersect3 == true) {
		m_canvas_2.add(i3, vvr::Colour::red);
	}
	// h panw grammh me oles tou trigwnou
	i4 = isOnSeg(p2, p, p4, p5);
	bool seg_intersect4 = isIntersect(p2, p, p4, p5, i4);
	if (seg_intersect4 == true) {
		m_canvas_2.add(i4, vvr::Colour::red);
	}

	i5 = isOnSeg(p2, p, p5, p6);
	bool seg_intersect5 = isIntersect(p2, p, p5, p6, i5);
	if (seg_intersect5 == true) {
		m_canvas_2.add(i5, vvr::Colour::red);
	}
	i6 = isOnSeg(p2, p, p4, p6);
	bool seg_intersect6 = isIntersect(p2, p, p4, p6, i6);
	if (seg_intersect6 == true) {
		m_canvas_2.add(i6, vvr::Colour::red);
	}

	C2DPointSet pts_tri;
	pts_tri.AddCopy(p1);
	pts_tri.AddCopy(p2);
	pts_tri.AddCopy(p); 
	pts_tri.AddCopy(p4);
	pts_tri.AddCopy(p5);
	pts_tri.AddCopy(p6);

	int size = pts_tri.size();
	echo(size);

	C2DPoint *pts = new C2DPoint[size];
	for (int i = 0; i < size; i++) {
		pts[i] = *pts_tri.GetAt(i);
	}

	// Create convex hull
	C2DPolygon cloud_polygon;
	cloud_polygon.Create(pts, size);
	m_tri.CreateConvexHull(cloud_polygon);
	delete[] pts;

	//kuklika oles
	C2DLine line1(p4, p5);
	C2DLine line2(p5, p6);
	C2DLine line3(p6, p4);
	C2DLine line4(p2, p4);
	C2DLine line5(p6, p1);
	C2DLine line6(p1, p2);;
	C2DLineSet lines_tri;
	bool aristera1 = false;
	bool aristera2 = false;
	bool aristera3 = false;

	bool aristera4 = false;
	bool aristera5 = false;
	bool aristera6 = false;
	if (isOnLeftC2D(line1, p) > 0) {
		aristera1 = true;
	}
	if (isOnLeftC2D(line2, p) > 0) {
		aristera2 = true;
	}
	if (isOnLeftC2D(line3, p) > 0) {
		aristera3 = true;
	}

	if (isOnLeftC2D(line4, p) > 0) {
		aristera4 = true;
	}
	if (isOnLeftC2D(line5, p) > 0) {
		aristera5 = true;
	}
	if (isOnLeftC2D(line6, p) > 0) {
		aristera6 = true;
	}
	C2DTriangle triangle_mesa1;
	C2DTriangle triangle_mesa2;
	C2DTriangle triangle_mesa3;
	if (aristera1 == false & aristera2 == false & aristera3==false) {
		triangle_mesa1.Set(i6, p4, p2);
		triangle_mesa2.Set(i3, p6, p1);

	}
	//an einai mesa sto parallhlogramo p2,p4,p6,p1
	


	if (aristera4 == false & aristera5==false & aristera6 == false & aristera3 == true) {
		triangle_mesa1.Set(p, p4, p2);
		triangle_mesa2.Set(p, p4, p6);
		triangle_mesa3.Set(p, p6, p1);

	}
	//aristero trigwno
	if (aristera6 == false & aristera4 == true & aristera1 == false) {
		triangle_mesa1.Set(p, p4, p1);
		triangle_mesa2.Set(p1, p4, p6);

	}
	//deksi trigwno////////////////////...//// & seg_intersect1 == false
	if (aristera6 == false & aristera2 == false & aristera5 == true & seg_intersect1 == false ) {
		triangle_mesa1.Set(p, p6, p2);
		triangle_mesa2.Set(p2, p4, p6);
	}
	//pio aristero apo to aristero
	if (aristera3 == true & aristera1 == true & aristera6==false) {
		triangle_mesa1.Set(p, p5, p4);
		triangle_mesa2.Set(p, p4, p1);
		triangle_mesa3.Set(p6, p4, p1);


	}
	// pio deksi apo to deksi
	if (aristera3 == true & aristera2 == true & aristera6==false & seg_intersect1 == false) {
		triangle_mesa1.Set(p, p5, p6);
		triangle_mesa2.Set(p, p2, p6);
		triangle_mesa3.Set(p6, p4, p2);


	}
	//apo pisw(mono aristera ths line6(p1,p2))
	C2DTriangle triangle_mesa4;
	if (aristera6 == true ) {
		triangle_mesa1.Set(p4, p2, p1);
		triangle_mesa2.Set(p1, p4, p6);
		//piso aristera
		if (aristera1 == false & aristera4 == true) {
			triangle_mesa3.Set(p, p4, p2);
		}
		//piso deksia
		if (aristera2 == false & aristera5 == true) {
			triangle_mesa3.Set(p, p1, p6);
		}
		//piso pio aristera
		if (aristera1 == true & aristera4 == true) {
			triangle_mesa3.Set(p, p4, p5);
			triangle_mesa4.Set(p, p4, p2);
		}
		//piso pio deksia
		if (aristera2 == true & aristera5 == true) {
			triangle_mesa3.Set(p, p5, p6);
			triangle_mesa4.Set(p, p1, p6);
		}
	}
	// dio tomes sthn line1
	C2DTriangle triangle_mesa5;
	if (seg_intersect1 == true & seg_intersect4 == true ) {
		triangle_mesa1.Set(p, i4, p4);
		triangle_mesa2.Set(p, i1, p5);
		triangle_mesa3.Set(p2, i6, p4);
		triangle_mesa4.Set(p1, i3, p6);
		if (aristera2 == true) {
			triangle_mesa5.Set(p, p5, p6);
			
		}

	}
	//mia tomh sthn line 1
	if (seg_intersect1 == true & seg_intersect4 == false) {
		triangle_mesa1.Set(p, i1, p5);
		triangle_mesa2.Set(p1, i3, p6);
	}
	//panw apo thn line3 kai line1 alla den temnete me to trigwno katholou
	if (seg_intersect1 == false & seg_intersect4 == false & aristera3 == false & aristera1 == true) {
		triangle_mesa1.Set(p, p4, p5);
		triangle_mesa2.Set(p, p4, p1);
		triangle_mesa3.Set(p1, p4, p6);
	}
	//2 tomes sthn line 2
	if (seg_intersect2 == true & seg_intersect5 == true) {
		triangle_mesa1.Set(p, i5, p5);
		triangle_mesa2.Set(p, i2, p6);
		triangle_mesa3.Set(p1, i3, p6);
		triangle_mesa4.Set(p2, i6, p4);
		if (aristera1 == true) {
			triangle_mesa5.Set(p, p4, p5);
		}
	}
	//mia tomh sthn line 2
	if (seg_intersect5 == true & seg_intersect2 == false) {
		triangle_mesa1.Set(p, i5, p5);
		triangle_mesa2.Set(p2, i6, p4);
	}
	//panw apo thn line2 kai line3 alla den temnete me to trigwno katholou
	if (seg_intersect5 == false & seg_intersect2 == false & seg_intersect1 == false & aristera3 == false & aristera2 == true) {
		triangle_mesa1.Set(p, p6, p5);
		triangle_mesa2.Set(p, p6, p2);
		triangle_mesa3.Set(p2, p4, p6);
	}
	//mia tomh se line1 kai mia tomh se line2
	if (seg_intersect4 == true & seg_intersect2 == true ) {
		triangle_mesa1.Set(p, i4, p4);
		triangle_mesa2.Set(p, i2, p6);
		triangle_mesa3.Set(p2, p4, i6);
		triangle_mesa4.Set(p1, p6, i3);
	}
	m_canvas_2.add(triangle_mesa1, Colour::green, true);
	m_canvas_2.add(triangle_mesa2, Colour::green, true);
	m_canvas_2.add(triangle_mesa3, Colour::green, true);
	m_canvas_2.add(triangle_mesa4, Colour::green, true);
	m_canvas_2.add(triangle_mesa5, Colour::green, true);

}

void SceneIntersections::Task3(const C2DPoint &p)
{
	m_circle_2.SetCentre(p);
	m_canvas_3.clear();
	m_canvas_3.add(p, Colour::cyan);
	m_canvas_3.add(m_circle_2, Colour::cyan);

	const double x1 = m_circle_1.GetCentre().x;
	const double y1 = m_circle_1.GetCentre().y;
	const double r1 = m_circle_1.GetRadius();
	const double x2 = m_circle_2.GetCentre().x;
	const double y2 = m_circle_2.GetCentre().y;
	const double r2 = m_circle_2.GetRadius();

	float a1 = 2 * (x2 - x1);
	float a2 = 2 * (y2 - y1);
	float a3 = pow(x1, 2) - pow(x2, 2) + pow(y1, 2) - pow(y2, 2) + pow(r2, 2) - pow(r1, 2);

	float b1 = a2 / a1;
	float b2 = a3 / a1;

	float c1 = pow(x1, 2) + pow(y1, 2) - pow(r1, 2);
	float c2 = pow(b2, 2) + 2*b2*x1 +c1;

	float c3 = 2 * b1*b2 + 2 * b1*x1 - 2 * y1;
	float c4 = 1+ pow(b1, 2);

	float D = pow(c3, 2) - 4 * c4*c2;

	float y_1_found = (-c3 + sqrt(D)) / (2 * c4);
	float y_2_found = (-c3 - sqrt(D)) / (2 * c4);
	float x_1_found = -b1*y_1_found - b2;
	float x_2_found = -b1*y_2_found - b2;

	C2DPoint i1(x_1_found, y_1_found);
	C2DPoint i2(x_2_found, y_2_found);
	m_canvas_3.add(i1, Colour::red);
	m_canvas_3.add(i2, Colour::red);

	//pts_kinoumenou = new C2DPoint[60];
	C2DPointSet pts_kinoumenou;

	float deigmato = 2*r2;
	float deigma = deigmato/200;

	//float x_new = pts_kinoumenou[0].x;
	float y_new,y_new_2;
	float z1, z2, Dz;
	int i = 0;

	float x_new = x2-r2;
	
	while(x_new < x2+r2){
		x_new = deigma + x_new;
		z1 = -2 * y2;
		z2 = pow(y2, 2) + pow(x_new, 2) - 2 * x_new*x2 + pow(x2, 2) - pow(r2, 2);
		Dz = pow(z1, 2) - 4 * z2;
		y_new = (-z1 + sqrt(Dz)) / (2);
		y_new_2 = (-z1 - sqrt(Dz)) / (2);
		C2DPoint i_sample(x_new, y_new);
		C2DPoint i_sample1(x_new, y_new_2);
		float sample_apost = sqrt(pow((i_sample.x - x1), 2) + pow((i_sample.y - y1), 2));
		float sample_apost1 = sqrt(pow((i_sample1.x - x1), 2) + pow((i_sample1.y - y1), 2));
		// pare mono ta shmeia pou exoun vriskontai mesa ston stathero kuklo
		if (sample_apost < r1) {
			pts_kinoumenou.AddCopy(i_sample);
			//pts_kinoumenou[i] = i_sample;
			//m_canvas_3.add(pts_kinoumenou[i], Colour::red);
		}
		if (sample_apost1 < r1) {
			//pts_kinoumenou[i + 1] = i_sample1;
			pts_kinoumenou.AddCopy(i_sample1);
			//m_canvas_3.add(pts_kinoumenou[i + 1], Colour::red);
		}
		i = i + 2;
	}
	//deigmatolipsia staherou kuklou
	float deigmato1 = 2 * r1;
	float deigma1 = deigmato1 / 200;

	//float x_new = pts_kinoumenou[0].x;

	x_new = x1 - r1;

	while (x_new < x1 + r1) {
		x_new = deigma1 + x_new;
		z1 = -2 * y1;
		z2 = pow(y1, 2) + pow(x_new, 2) - 2 * x_new*x1 + pow(x1, 2) - pow(r1, 2);
		Dz = pow(z1, 2) - 4 * z2;
		y_new = (-z1 + sqrt(Dz)) / (2);
		y_new_2 = (-z1 - sqrt(Dz)) / (2);
		C2DPoint i_sample(x_new, y_new);
		C2DPoint i_sample1(x_new, y_new_2);
		float sample_apost = sqrt(pow((i_sample.x - x2), 2) + pow((i_sample.y - y2), 2));
		float sample_apost1 = sqrt(pow((i_sample1.x - x2), 2) + pow((i_sample1.y - y2), 2));
		// pare mono ta shmeia pou exoun vriskontai mesa ston stathero kuklo
		if (sample_apost < r2) {
			pts_kinoumenou.AddCopy(i_sample);
			//pts_kinoumenou[i] = i_sample;
			//m_canvas_3.add(pts_kinoumenou[i], Colour::red);
		}

		if (sample_apost1 < r2) {
			pts_kinoumenou.AddCopy(i_sample1);
			//pts_kinoumenou[i + 1] = i_sample1;
			//m_canvas_3.add(pts_kinoumenou[i + 1], Colour::red);
		}
		i = i + 2;
	}
	int size = pts_kinoumenou.size();
	echo(size);

	C2DPoint *pts = new C2DPoint[size];
	for (int i = 0; i < size; i++) {
		pts[i] = *pts_kinoumenou.GetAt(i);
	}

	// Create convex hull
	C2DPolygon cloud_polygon;
	cloud_polygon.Create(pts, size);
	m_mesa.CreateConvexHull(cloud_polygon);
	delete[] pts;

}


void SceneIntersections::draw()
{
	enterPixelMode();
	vvr::draw(m_mesa, vvr::Colour::green,true);
	//vvr::draw(m_tri, vvr::Colour::green);
	m_canvas_0.draw();
	m_canvas_1.draw();
	m_canvas_2.draw();
	m_canvas_3.draw();
}

/* Application Entry Point */

int main(int argc, char* argv[])
{
	return vvr::mainLoop(argc, argv, new SceneIntersections);
}


C2DPoint isOnSeg(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3) {
	//vvr::Point2D *im;
	//im = new vvr::Point2D;
	//vvr::Point2D i = *im; //kanei to pointer se vvr
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
bool isIntersect(C2DPoint p1, C2DPoint p, C2DPoint p2, C2DPoint p3,C2DPoint i){
	bool seg_intersect = false;

	if (i.x > std::min(p3.x, p2.x) && i.x < std::max(p3.x, p2.x)) {
		if (i.x > std::min(p1.x, p.x) && i.x < std::max(p1.x, p.x)) {
			seg_intersect = true;
		}
	}
	
	return seg_intersect;
}
int isOnLeftC2D(C2DLine line, C2DPoint t) {
	C2DPoint end = line.GetPointTo();
	C2DPoint start = line.GetPointFrom();
	float shmeio1 = (end.x - start.x)*(t.y - start.y) - (end.y - start.y)*(t.x - start.x);

	if (shmeio1 > 0) {
		return 1; //aristera
	}
	else {
		return -1; // deksia
	}

}