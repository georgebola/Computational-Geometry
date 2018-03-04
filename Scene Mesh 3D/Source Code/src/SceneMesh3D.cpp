#include "SceneMesh3D.h"

#define SPLIT_INSTEAD_OF_INTERSECT 1

using namespace std;
using namespace vvr;

Mesh3DScene::Mesh3DScene()
{
	//! Load settings.
	vvr::Shape::DEF_LINE_WIDTH = 4;
	vvr::Shape::DEF_POINT_SIZE = 10;
	m_perspective_proj = true;
	m_bg_col = Colour("768E77");
	m_obj_col = Colour("454545");
	const string objDir = getBasePath() + "resources/obj/";
	const string objFile = objDir + "armadillo_low_low.obj";
	m_model_original = vvr::Mesh(objFile);
	reset();
}
void Mesh3DScene::reset()
{
	Scene::reset();

	//! Define plane
	m_plane_d = 0;
	m_plane = Plane(vec(0, 1, 1).Normalized(), m_plane_d);

	//! Define what will be vissible by default
	m_style_flag = 0;
	m_style_flag |= FLAG_SHOW_SOLID;
	m_style_flag |= FLAG_SHOW_WIRE;
	m_style_flag |= FLAG_SHOW_AXES;
	m_style_flag |= FLAG_SHOW_AABB;
	//m_style_flag |= FLAG_SHOW_PLANE;
}

void Mesh3DScene::resize()
{
	//! By Making `first_pass` static and initializing it to true,
	//! we make sure that the if block will be executed only once.

	static bool first_pass = true;

	if (first_pass)
	{
		m_model_original.setBigSize(getSceneWidth() / 2);
		m_model_original.update();
		m_model = m_model_original;
		Tasks();
		first_pass = false;
	}
}

void Mesh3DScene::Tasks()
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! Task 1
	//!//////////////////////////////////////////////////////////////////////////////////

	m_center_mass = vec(-10, 0, 0);
	Task_1_FindCenterMass(m_model.getVertices(), m_center_mass);

	//!//////////////////////////////////////////////////////////////////////////////////
	//! Task 2
	//!//////////////////////////////////////////////////////////////////////////////////

	Task_2_FindAABB(m_model.getVertices(), m_aabb);

	//!//////////////////////////////////////////////////////////////////////////////////
	//! Task 3
	//!//////////////////////////////////////////////////////////////////////////////////

	Task_3_AlignOriginTo(m_model.getVertices(), Task_3_Pick_Origin());
	Task_1_FindCenterMass(m_model.getVertices(), m_center_mass);
	Task_2_FindAABB(m_model.getVertices(), m_aabb);

	//!//////////////////////////////////////////////////////////////////////////////////
	//! Task 4
	//!//////////////////////////////////////////////////////////////////////////////////

	pca(m_model.getVertices(), m_pca_cen, m_pca_dir);

	//!//////////////////////////////////////////////////////////////////////////////////
	//! Task 5
	//!//////////////////////////////////////////////////////////////////////////////////

	m_model_original = m_model;

	if (SPLIT_INSTEAD_OF_INTERSECT == false) {
		m_intersections.clear();
		m_i_points.clear();
		new_triangles.clear();

		Task_5_Intersect(m_model.getTriangles(), m_plane, m_intersections, m_i_points);
	}
	else {
		m_i_points.clear();
		m_model = Mesh(m_model_original);
		Task_5_Split(m_model, m_plane, m_i_points, new_triangles,m_greenmesh);
	}
}

void Mesh3DScene::arrowEvent(ArrowDir dir, int modif)
{
	math::vec n = m_plane.normal;
	if (dir == UP) m_plane_d += 1;
	if (dir == DOWN) m_plane_d -= 1;
	else if (dir == LEFT) n = math::float3x3::RotateY(DegToRad(1)).Transform(n);
	else if (dir == RIGHT) n = math::float3x3::RotateY(DegToRad(-1)).Transform(n);
	m_plane = Plane(n.Normalized(), m_plane_d);

	if (SPLIT_INSTEAD_OF_INTERSECT == false) {
		m_intersections.clear();
		m_i_points.clear();
		Task_5_Intersect(m_model.getTriangles(), m_plane, m_intersections, m_i_points);

	}
	else {
		m_i_points.clear();
		new_triangles.clear();
		m_model = Mesh(m_model_original);
		Task_5_Split(m_model, m_plane, m_i_points,new_triangles, m_greenmesh);
	}

}

void Mesh3DScene::keyEvent(unsigned char key, bool up, int modif)
{
	Scene::keyEvent(key, up, modif);
	key = tolower(key);

	switch (key)
	{
	case 's': m_style_flag ^= FLAG_SHOW_SOLID; break;
	case 'w': m_style_flag ^= FLAG_SHOW_WIRE; break;
	case 'n': m_style_flag ^= FLAG_SHOW_NORMALS; break;
	case 'a': m_style_flag ^= FLAG_SHOW_AXES; break;
	case 'p': m_style_flag ^= FLAG_SHOW_PLANE; break;
	case 'b': m_style_flag ^= FLAG_SHOW_AABB; break;
	case 'i': m_style_flag ^= INTERSECTION_POINTS; break;
	}
}

void Mesh3DScene::draw()
{
	//! Draw plane
	if (m_style_flag & FLAG_SHOW_PLANE) {
		vvr::Colour colPlane(0x41, 0x14, 0xB3);
		float u = 20, v = 20;
		math::vec p0(m_plane.Point(-u, -v, math::vec(0, 0, 0)));
		math::vec p1(m_plane.Point(-u, v, math::vec(0, 0, 0)));
		math::vec p2(m_plane.Point(u, -v, math::vec(0, 0, 0)));
		math::vec p3(m_plane.Point(u, v, math::vec(0, 0, 0)));
		math2vvr(math::Triangle(p0, p1, p2), colPlane).draw();
		math2vvr(math::Triangle(p2, p1, p3), colPlane).draw();
	}

	if (m_style_flag & FLAG_SHOW_SOLID) m_model.draw(m_obj_col, SOLID);
	if (m_style_flag & FLAG_SHOW_WIRE) m_model.draw(Colour::black, WIRE);
	if (m_style_flag & FLAG_SHOW_NORMALS) m_model.draw(Colour::black, NORMALS);
	if (m_style_flag & FLAG_SHOW_AXES) m_model.draw(Colour::black, AXES);

	//! Draw pca line
	Task_4_Draw_PCA(m_pca_cen, m_pca_dir);

	//! Draw center mass
	Point3D(m_center_mass.x, m_center_mass.y, m_center_mass.z, Colour::red).draw();

	//! Draw AABB
	if (m_style_flag & FLAG_SHOW_AABB) {
		m_aabb.setColour(Colour::black);
		m_aabb.setTransparency(1);
		m_aabb.draw();
	}

	//! Draw intersecting triangles of model
	//just_draw(m_model.getTriangles(), m_intersections);

	///! ERWTHMA 2
	erwthma_2(m_model.getTriangles(),m_intersections, m_i_points, m_plane );

	///!ERWTHMA 3
	// Draw my new green split mesh, solid and wireframe
	m_greenmesh.draw(Colour::green, SOLID);
	m_greenmesh.draw(Colour::black, WIRE);

	///PRESS i to draw inter_points
	//Draw Intersection_points 
	if (m_style_flag & INTERSECTION_POINTS) {
		for (int i = 0; i < m_i_points.size(); i++) {
			Point3D(m_i_points[i].x, m_i_points[i].y, m_i_points[i].z, Colour::red).draw();
		}
	}

}

void pca(vector<vec>& vertices, vec &center, vec &dir)
{
	const int count = vertices.size();

	float w0 = 0;
	float x0 = 0, y0 = 0, z0 = 0;
	float x2 = 0, y2 = 0, z2 = 0, xy = 0, yz = 0, xz = 0;
	float dx2, dy2, dz2, dxy, dxz, dyz;
	float det[9];

	for (int i = 0; i < count; i++)
	{
		float x = vertices[i].x;
		float y = vertices[i].y;
		float z = vertices[i].z;

		x2 += x * x;
		xy += x * y;
		xz += x * z;
		y2 += y * y;
		yz += y * z;
		z2 += z * z;
		x0 += x;
		y0 += y;
		z0 += z;
	}
	w0 = (float)count;

	x2 /= w0;
	xy /= w0;
	xz /= w0;
	y2 /= w0;
	yz /= w0;
	z2 /= w0;

	x0 /= w0;
	y0 /= w0;
	z0 /= w0;

	dx2 = x2 - x0 * x0;
	dxy = xy - x0 * y0;
	dxz = xz - x0 * z0;
	dy2 = y2 - y0 * y0;
	dyz = yz - y0 * z0;
	dz2 = z2 - z0 * z0;

	det[0] = dz2 + dy2;
	det[1] = -dxy;
	det[2] = -dxz;
	det[3] = det[1];
	det[4] = dx2 + dz2;
	det[5] = -dyz;
	det[6] = det[2];
	det[7] = det[5];
	det[8] = dy2 + dx2;

	/* Searching for a eigenvector of det corresponding to the minimal eigenvalue */
	gte::SymmetricEigensolver3x3<float> solver;
	std::array<float, 3> eval;
	std::array<std::array<float, 3>, 3> evec;
	solver(det[0], det[1], det[2], det[4], det[5], det[8], true, 1, eval, evec);

	center.x = x0;
	center.y = y0;
	center.z = z0;

	dir.x = evec[0][0];
	dir.y = evec[0][1];
	dir.z = evec[0][2];
}

int main(int argc, char* argv[])
{
	try {
		return vvr::mainLoop(argc, argv, new Mesh3DScene);
	}
	catch (std::string exc) {
		cerr << exc << endl;
		return 1;
	}
	catch (...)
	{
		cerr << "Unknown exception" << endl;
		return 1;
	}
}

//! LAB Tasks

void Task_1_FindCenterMass(vector<vec> &vertices, vec &cm)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Breite to kentro mazas twn simeiwn `vertices`.
	//!  - Apothikeyste to apotelesma stin metavliti `cm`.
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	int sumx = 0;
	int sumy = 0;
	int sumz = 0;
	for (int i = 0; i < vertices.size(); i++) {
		sumx = sumx + vertices[i].x;
		sumy = sumy + vertices[i].y;
		sumz = sumz + vertices[i].z;
	}
	cm.x = sumx / float(vertices.size());
	cm.y = sumy / float(vertices.size());
	cm.z = sumz / float(vertices.size());

}

void Task_2_FindAABB(vector<vec> &vertices, Box3D &aabb)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Breite to Axis Aligned Bounding Box tou montelou
	//!
	//! HINTS:
	//!
	//!  - To `aabb` orizetai apo 2 gwniaka simeia. (V_min, V_max)
	//!  - V_min: { aabb.x1, aabb.y1, aabb.z1 }
	//!  - V_max: { aabb.x2, aabb.y2, aabb.z2 }
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	float minx = vertices[0].x;
	float miny = vertices[0].y;
	float minz = vertices[0].z;
	float maxx = vertices[0].x;
	float maxy = vertices[0].y;
	float maxz = vertices[0].z;
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].x < minx) {
			minx = vertices[i].x;
		}
		if (vertices[i].y < miny) {
			miny = vertices[i].y;
		}
		if (vertices[i].z < minz) {
			minz = vertices[i].z;
		}
		if (vertices[i].x > maxx) {
			maxx = vertices[i].x;
		}
		if (vertices[i].y > maxy) {
			maxy = vertices[i].y;
		}
		if (vertices[i].z > maxz) {
			maxz = vertices[i].z;
		}

	}
	aabb.x1 = minx;
	aabb.x2 = maxx;
	aabb.y1 = miny;
	aabb.y2 = maxy;
	aabb.z1 = minz;
	aabb.z2 = maxz;

}

vec Mesh3DScene::Task_3_Pick_Origin()
{
	return m_center_mass;
}

void Task_3_AlignOriginTo(vector<vec> &vertices, const vec &cm)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Metatopiste to montelo esti wste to simeio `cm` na erthei sto (0,0,0).
	//!
	//!//////////////////////////////////////////////////////////////////////////////////


	for (int i = 0; i < vertices.size(); i++) {
		vertices[i].x = vertices[i].x - cm.x;
		vertices[i].y = vertices[i].y - cm.y;
		vertices[i].z = vertices[i].z - cm.z;
	}


}

void Task_4_Draw_PCA(vec &center, vec &dir)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Apeikoniste to kentro mazas kai ton Principal Axis tou PCA.
	//!    Gia tin apeikonisi, xreiazeste ena simeio kai mia eytheia.
	//!
	//! HINTS:
	//!  - Auti i synartisi kaleitai mesa apo tin `Mesh3DScene::draw()`.
	//!    Ara, mporeite na kalesete amesa tis metodous draw ton diaforwn antikeimenwn
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	float xnew = center.x + dir.x * 30;
	float ynew = center.y + dir.y * 30;
	float znew = center.z + dir.z * 30;
	LineSeg3D line(center.x, center.y, center.z, xnew, ynew, znew, vvr::Colour::magenta);
	line.draw();
}

void Task_5_Intersect(vector<vvr::Triangle> &triangles, Plane &plane, vector<int> &intersection_indices, vector<vec> &intersection_points)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Brete ta trigwna pou temnontai me to epipedo `plane`.
	//!  - Kante ta push_back sto vector intersection_indices.
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < triangles.size(); i++) {
		float distance1 = triangles[i].v1().Dot(plane.normal) - plane.d;
		float distance2 = triangles[i].v2().Dot(plane.normal) - plane.d;
		float distance3 = triangles[i].v3().Dot(plane.normal) - plane.d;
		if (distance1*distance2 < 0 || distance1*distance3 < 0) {
			intersection_indices.push_back(i);
			//triangles.erase(triangles.begin() + i--);

			if (distance1*distance2 < 0) {
				//lamda is the direction
				//t is the distance 
				//triangles[i].v1() is the starting point
				vec lamda = triangles[i].v1() - triangles[i].v2();
				float distance = plane.normal.Dot(lamda);
				float t = -(plane.normal.Dot(triangles[i].v1()) + plane.d) / distance;
				
				float ix = triangles[i].v1().x + t * lamda.x;
				float iy = triangles[i].v1().y + t * lamda.y;
				float iz = triangles[i].v1().z + t * lamda.z;

				intersection_points.push_back(vec(ix, iy, iz));
			}
			if (distance1*distance3 < 0) {
				vec lamda = triangles[i].v1() - triangles[i].v3();

				float distance = plane.normal.Dot(lamda);
				float t = -(plane.normal.Dot(triangles[i].v1()) + plane.d) / distance;

				float ix = triangles[i].v1().x + t * lamda.x;
				float iy = triangles[i].v1().y + t * lamda.y;
				float iz = triangles[i].v1().z + t * lamda.z;
				intersection_points.push_back(vec(ix, iy, iz));
			}
			if (distance2*distance3 < 0) {
				vec lamda = triangles[i].v2() - triangles[i].v3();

				float distance = plane.normal.Dot(lamda);
				float t = -(plane.normal.Dot(triangles[i].v2()) + plane.d) / distance;

				float ix = triangles[i].v2().x + t * lamda.x;
				float iy = triangles[i].v2().y + t * lamda.y;
				float iz = triangles[i].v2().z + t * lamda.z;
				intersection_points.push_back(vec(ix, iy, iz));
			}
		}
		
	}

}

void Task_5_Split(Mesh &mesh, Plane &plane,  vector<vec> &intersection_points, vector<vvr::Triangle> &triangles,Mesh &m_greenmesh)
{
	//!//////////////////////////////////////////////////////////////////////////////////
	//! TASK:
	//!
	//!  - Kopste to antikeimeno sta 2. (Odigies tin ora tou ergasthriou)
	//!
	//!//////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < mesh.getTriangles().size(); i++) {
		float distance1 = mesh.getTriangles()[i].v1().Dot(plane.normal) - plane.d;
		float distance2 = mesh.getTriangles()[i].v2().Dot(plane.normal) - plane.d;
		float distance3 = mesh.getTriangles()[i].v3().Dot(plane.normal) - plane.d;
		if (distance1*distance2 < 0) {
			vec lamda = mesh.getTriangles()[i].v1() - mesh.getTriangles()[i].v2();
			float distance = plane.normal.Dot(lamda);
			float t = -(plane.normal.Dot(mesh.getTriangles()[i].v1())+ plane.d) / distance;
			//plane.IntersectLinePlane(plane.normal, plane.d, mesh.getTriangles()[i].v1(), lamda, t);
			
			//lamda is the direction
			//t is the distance 
			//triangles[i].v1() is the starting point
			float ix = mesh.getTriangles()[i].v1().x + t * lamda.x;
			float iy = mesh.getTriangles()[i].v1().y + t * lamda.y;
			float iz = mesh.getTriangles()[i].v1().z + t * lamda.z;
			//ta pernaw kai panw kai katw
			intersection_points.push_back(vec(ix + 3 * plane.normal.x, iy + 3 * plane.normal.y, iz + 3 * plane.normal.z));

			intersection_points.push_back(vec(ix - 3 * plane.normal.x, iy - 3 * plane.normal.y, iz - 3 * plane.normal.z));


		}
		if (distance1*distance3 < 0) {
			vec lamda = mesh.getTriangles()[i].v1() - mesh.getTriangles()[i].v3();

			float distance = plane.normal.Dot(lamda);
			float t = -(plane.normal.Dot(mesh.getTriangles()[i].v1()) + plane.d) / distance;
			//plane.IntersectLinePlane(plane.normal, plane.d, mesh.getTriangles()[i].v1(), lamda, t);

			float ix = mesh.getTriangles()[i].v1().x + t * lamda.x;
			float iy = mesh.getTriangles()[i].v1().y + t * lamda.y;
			float iz = mesh.getTriangles()[i].v1().z + t * lamda.z;

			intersection_points.push_back(vec(ix + 3 * plane.normal.x, iy + 3 * plane.normal.y, iz + 3 * plane.normal.z));
			intersection_points.push_back(vec(ix - 3 * plane.normal.x, iy - 3 * plane.normal.y, iz - 3 * plane.normal.z));
			
		}
		if (distance2*distance3 < 0) {
			vec lamda = mesh.getTriangles()[i].v2() - mesh.getTriangles()[i].v3();
			float distance = plane.normal.Dot(lamda);
			float t = -(plane.normal.Dot(mesh.getTriangles()[i].v2()) + plane.d) / distance;
			//plane.IntersectLinePlane(plane.normal, plane.d, mesh.getTriangles()[i].v2(), lamda, t);

			float ix = mesh.getTriangles()[i].v2().x + t * lamda.x;
			float iy = mesh.getTriangles()[i].v2().y + t * lamda.y;
			float iz = mesh.getTriangles()[i].v2().z + t * lamda.z;

			intersection_points.push_back(vec(ix + 3 * plane.normal.x, iy + 3 * plane.normal.y, iz + 3 * plane.normal.z));
			intersection_points.push_back(vec(ix - 3 * plane.normal.x, iy - 3 * plane.normal.y, iz - 3 * plane.normal.z));
			
		}
		if (distance1*distance2< 0 || distance1*distance3 < 0) {
			//vvr::Triangle &t = mesh.getTriangles[i];
			
			triangles.push_back(mesh.getTriangles()[i]);
			mesh.getTriangles().erase(mesh.getTriangles().begin() + i--);
			
		}

	}
	for (int i = 0; i < mesh.getVertices().size(); i++) {
		float over = mesh.getVertices()[i].Dot(plane.normal) - plane.d;
		if (over > 0) {
			mesh.getVertices()[i] +=3 * plane.normal;
		}
		if (over < 0) {
			mesh.getVertices()[i] -= 3 * plane.normal;
		}
	}
	int j = 0;
	for (int i = 0; i < triangles.size(); i++) {

		//exw m_i_points[j] kai m_i_poits[j+1] gia to kathe trigwno
		float distance1 = triangles[i].v1().Dot(plane.normal) - plane.d;
		float distance2 = triangles[i].v2().Dot(plane.normal) - plane.d;
		float distance3 = triangles[i].v3().Dot(plane.normal) - plane.d;
		vvr::Triangle &t = triangles[i];
		vvr::Point3D p1, p2, p3;
		float inner;
		if (distance1*distance2 > 0) { // an dhladh einai sthn idia meria

			p1.x = t.v1().x;
			p1.y = t.v1().y;
			p1.z = t.v1().z;
			p2.x = t.v2().x;
			p2.y = t.v2().y;
			p2.z = t.v2().z;
			p3.x = t.v3().x;
			p3.y = t.v3().y;
			p3.z = t.v3().z;
			inner = distance1;
		}


		else if (distance1*distance3 > 0) { // an dhladh einai sthn idia meria
			p1.x = t.v1().x;
			p1.y = t.v1().y;
			p1.z = t.v1().z;
			p2.x = t.v3().x;
			p2.y = t.v3().y;
			p2.z = t.v3().z;
			p3.x = t.v2().x;
			p3.y = t.v2().y;
			p3.z = t.v2().z;

			inner = distance1;
		}
		else if (distance2*distance3 > 0) { // an dhladh einai sthn idia meria
			p1.x = t.v2().x;
			p1.y = t.v2().y;
			p1.z = t.v2().z;
			p2.x = t.v3().x;
			p2.y = t.v3().y;
			p2.z = t.v3().z;
			p3.x = t.v1().x;
			p3.y = t.v1().y;
			p3.z = t.v1().z;

			inner = distance2;
		}
		if (inner>0) { //an einai apo panw dld

			m_greenmesh.getVertices().push_back(vec(p1.x, p1.y, p1.z)); // 1st vert
			m_greenmesh.getVertices().push_back(vec(p2.x, p2.y, p2.z)); // 2nd vert
			m_greenmesh.getVertices().push_back(intersection_points[j]); // 3rd vert

			int index1 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index1, index1 + 1, index1 + 2));
			/// 2o TRIGWNO
			m_greenmesh.getVertices().push_back(vec(p2.x, p2.y, p2.z)); // 4th vert
			m_greenmesh.getVertices().push_back(intersection_points[j]); // 5th vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 2]); //6th vert

			int index2 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index2, index2 + 1, index2 + 2));
			///3o trigwno
			m_greenmesh.getVertices().push_back(vec(p3.x, p3.y, p3.z)); // 7th vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 1]); // 8th vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 3]); // 9th vert

			int index3 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index3, index3 + 1, index3 + 2));

		}
		else {//an einai apo katw to 1o kai 2o shmeio dld
			m_greenmesh.getVertices().push_back(vec(p1.x, p1.y, p1.z)); // 1st vert
			m_greenmesh.getVertices().push_back(vec(p2.x, p2.y, p2.z)); // 2nd vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 1]); // 3rd vert

			int index1 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index1, index1 + 1, index1 + 2));
			/// 2o TRIGWNO
			m_greenmesh.getVertices().push_back(vec(p2.x, p2.y, p2.z)); // 4st vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 1]); // 5nd vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 3]); //6rd vert

			int index2 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index2, index2 + 1, index2 + 2));
			///3o trigwno
			m_greenmesh.getVertices().push_back(vec(p3.x, p3.y, p3.z)); // 7st vert
			m_greenmesh.getVertices().push_back(intersection_points[j]); // 8nd vert
			m_greenmesh.getVertices().push_back(intersection_points[j + 2]); // 9rd vert

			int index3 = m_greenmesh.getVertices().size() - 3;
			m_greenmesh.getTriangles().push_back(vvr::Triangle(&m_greenmesh.getVertices(), index3, index3 + 1, index3 + 2));
		}
		j = j + 4;
	}
	
}

void just_draw(vector<vvr::Triangle>& triangles, vector<int> &m_intersections) {
	for (int i = 0; i < m_intersections.size(); i++) {
		vvr::Triangle &t = triangles[m_intersections[i]];
		Triangle3D t3d(
			t.v1().x, t.v1().y, t.v1().z,
			t.v2().x, t.v2().y, t.v2().z,
			t.v3().x, t.v3().y, t.v3().z,
			Colour::green);
		t3d.draw();

	}
}

void erwthma_2(vector<vvr::Triangle>& triangles, vector<int> &m_intersections, vector<vec> &m_i_points, Plane &m_plane) {
	int j = 0;
	for (int i = 0; i < m_intersections.size(); i++) {

		float distance1 = triangles[m_intersections[i]].v1().Dot(m_plane.normal) - m_plane.d;
		float distance2 = triangles[m_intersections[i]].v2().Dot(m_plane.normal) - m_plane.d;
		float distance3 = triangles[m_intersections[i]].v3().Dot(m_plane.normal) - m_plane.d;
		vvr::Triangle &t = triangles[m_intersections[i]];

		float a2 = sqrt(pow(t.v1().x - t.v2().x, 2) + pow(t.v1().y - t.v2().y, 2) + pow(t.v1().z - t.v2().z, 2));
		float b2 = sqrt(pow(t.v1().x - t.v3().x, 2) + pow(t.v1().y - t.v3().y, 2) + pow(t.v1().z - t.v3().z, 2));
		float c2 = sqrt(pow(t.v2().x - t.v3().x, 2) + pow(t.v2().y - t.v3().y, 2) + pow(t.v2().z - t.v3().z, 2));
		float s2 = (a2 + b2 + c2) / 2;
		float A2 = sqrt(s2*(s2 - a2)*(s2 - b2)*(s2 - c2));
		float inner;
		float A1;
		vvr::Point3D opp;
		if (distance1*distance2 > 0) { // an dhladh einai sthn idia meria to 1 kai 2
			opp.x = t.v3().x;
			opp.y = t.v3().y;
			opp.z = t.v3().z;
			inner = distance3;
		}
		else if (distance1*distance3 > 0) { // an dhladh einai sthn idia meria to 1 kai 3
			opp.x = t.v2().x;
			opp.y = t.v2().y;
			opp.z = t.v2().z;
			inner = distance2;
		}
		else if (distance2*distance3 > 0) { // an dhladh einai sthn idia meria to 2 kai 3
			opp.x = t.v1().x;
			opp.y = t.v1().y;
			opp.z = t.v1().z;
			inner = distance1;
		}

		float a1 = sqrt(pow(opp.x - m_i_points[j].x, 2) + pow(opp.y - m_i_points[j].y, 2) + pow(opp.z - m_i_points[j].z, 2));
		float b1 = sqrt(pow(opp.x - m_i_points[j + 1].x, 2) + pow(opp.y - m_i_points[j + 1].y, 2) + pow(opp.z - m_i_points[j + 1].z, 2));
		float c1 = sqrt(pow(m_i_points[j + 1].x - m_i_points[j].x, 2) + pow(m_i_points[j + 1].y - m_i_points[j].y, 2) + pow(m_i_points[j + 1].z - m_i_points[j].z, 2));
		float s1 = (a1 + b1 + c1) / 2;
		A1 = sqrt(s1*(s1 - a1)*(s1 - b1)*(s1 - c1));

		vvr::Colour color;
		if ((inner>0 && A1>(A2 / 2)) || (inner < 0 && A1 < (A2 / 2))) {
			color = vvr::Colour::yellow;
		}
		else {
			color = vvr::Colour::green;
		}
		Triangle3D t4d(
			t.v1().x, t.v1().y, t.v1().z,
			t.v2().x, t.v2().y, t.v2().z,
			t.v3().x, t.v3().y, t.v3().z,
			color);
		t4d.draw();
		j = j + 2;
	}
	
}

