#include "SceneKDTrees.h"

#define AUTOPLAY false
#define DIMENSIONS 3
#define NUM_PTS_DEFAULT 100
#define SPHERE_RAD 2.5
#define POINT_SIZE 6
#define GND_WIDTH 40
#define GND_TOP 10
#define GND_BOTTOM -10
#define GND_DEPTH 10
#define SEC_PER_FLOOR 10

//! KDTreeScene::

KDTreeScene::KDTreeScene()
{
	vvr::Shape::DEF_LINE_WIDTH = 1;
	vvr::Shape::DEF_POINT_SIZE = 10;
	m_bg_col = vvr::Colour::grey;
	m_perspective_proj = true;
	m_hide_log = false;
	m_hide_sliders = false;
	m_fullscreen = false;
	if (AUTOPLAY) m_anim.update(true);
	m_KDTree = NULL;
	reset();
}

void KDTreeScene::reset()
{
	Scene::reset();
	m_kn = 1;
	m_current_tree_level = 0;

	//! Define what will be vissible by default
	m_flag = 0;
	m_flag |= FLAG(SHOW_NN);
	m_flag |= FLAG(SHOW_KDTREE);
	m_flag |= FLAG(SHOW_SPHERE);
	m_flag |= FLAG(SHOW_PTS_KDTREE);
	m_flag |= FLAG(SHOW_PTS_ALL);
	m_flag |= FLAG(SHOW_PTS_IN_SPHERE);
	m_flag |= FLAG(BRUTEFORCE);

	//! Define scene objects
	m_sphere = vvr::Sphere3D(-GND_WIDTH / 2, 0, 0, SPHERE_RAD, vvr::Colour::white);

	//! Create random points
	const float mw = getSceneWidth() * 0.3;
	const float mh = getSceneHeight() * 0.3;
	const float mz = std::min(mw, mh);

	if (FLAG_ON(m_flag, POINTS_ON_SURFACE)) {
		createSurfacePts(m_pts.empty() ? NUM_PTS_DEFAULT : m_pts.size());
	}
	else {
		createRandomPts(m_pts.empty() ? NUM_PTS_DEFAULT : m_pts.size());
	}

	delete m_KDTree;
	m_KDTree = new KDTree(m_pts);
	m_tree_invalidation_sec = -1;

	//! Reset animation
	m_anim.setTime(0);
}

void KDTreeScene::resize()
{
	static bool first_pass = true;
	if (first_pass) printKeyboardShortcuts();
	first_pass = false;
}

void KDTreeScene::createRandomPts(int num_pts)
{
	m_pts.clear();
	for (int i = 0; i < num_pts; ++i) {
		float x = GND_WIDTH * 0.8 * (m_lcg.Float() - 0.5f);
		float y = GND_BOTTOM + (GND_TOP - GND_BOTTOM) * 0.8 * (m_lcg.Float() + 0.1);
		float z = GND_DEPTH * 0.8 * (m_lcg.Float() - 0.5f);
		m_pts.push_back(vec(x, y, z));
	}
	m_pts.shrink_to_fit();
}

void KDTreeScene::createSurfacePts(int num_pts)
{
	num_pts /= 10;
	num_pts += 10; // To avoid zero points which would cause crash.

	m_pts.clear();

	const int nx = math::Sqrt(num_pts) * GND_WIDTH / GND_DEPTH / 2;
	const int nz = num_pts / nx;

	const float x0 = GND_WIDTH * 0.8 * -0.5;
	const float x1 = GND_WIDTH * 0.8 * +0.5;
	const float dx = (x1 - x0) / nx;

	const float z0 = GND_DEPTH * 0.8 * -0.5f;
	const float z1 = GND_DEPTH * 0.8 * +0.5f;
	const float dz = (z1 - z0) / nz;

	const float y0 = (GND_TOP + GND_BOTTOM) * 0.5 - 1;
	const float y1 = (GND_TOP + GND_BOTTOM) * 0.5 + 1;

	for (size_t xi = 0; xi < nx; xi++)
	{
		for (size_t zi = 0; zi < nz; zi++)
		{
			float x = x0 + dx * xi;
			float z = z0 + dz * zi;
			float y = y0 + sin(x / 2) * sin(z / 2) * (y1 - y0);
			m_pts.push_back(vec(x, y, z));
		}
	}

	m_pts.shrink_to_fit();
}

bool KDTreeScene::idle()
{
	if (m_tree_invalidation_sec > 0 &&
		vvr::getSeconds() - m_tree_invalidation_sec > 0.8)
	{
		delete m_KDTree;
		m_KDTree = new KDTree(m_pts);
		m_tree_invalidation_sec = -1;
	}
	m_anim.update();
	return true;
}

void KDTreeScene::draw()
{
	const float POINT_SIZE_SAVE = vvr::Shape::DEF_POINT_SIZE;
	float duration1,duration2;
	//! Draw axes
	if (m_flag & FLAG(SHOW_AXES)) drawAxes();

	//! Draw ground
	Ground ground(GND_WIDTH, GND_DEPTH, GND_BOTTOM, GND_TOP, vvr::Colour(35, 45, 55));
	ground.draw();

	//! Animate sphere
	float t = m_anim.t;
	vvr::Sphere3D sphere_moved(m_sphere);
	sphere_moved.x += t * ((float)GND_WIDTH / SEC_PER_FLOOR);
	vec sc = vec(sphere_moved.x, sphere_moved.y, sphere_moved.z);
	Sphere sphere(sc, sphere_moved.rad);
	if (sphere_moved.x > GND_WIDTH / 2) m_anim.setTime(0); // Bring back to start

														   //! Draw points
	if (FLAG_ON(m_flag, SHOW_PTS_ALL)) {
		for (size_t i = 0; i < m_pts.size(); i++) {
			vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
			math2vvr(m_pts[i], vvr::Colour::white).draw();
			vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;
		}
	}

	//! Draw points in sphere
	if (FLAG_ON(m_flag, SHOW_PTS_IN_SPHERE)) {
		if (FLAG_ON(m_flag, SHOW_SPHERE)) {
			sphere_moved.draw();
		}
		VecVector pts_in;
		duration1 = vvr::getSeconds();
		if (FLAG_ON(m_flag, BRUTEFORCE)) {
			for (size_t i = 0; i < m_KDTree->pts.size(); i++)
				if (sphere.Contains(m_KDTree->pts.at(i))) pts_in.push_back(m_KDTree->pts.at(i));
		}
		else {
			Task_03_InSphere(sphere, m_KDTree->root(), pts_in);
		}
		duration2 = vvr::getSeconds();
		vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
		for (size_t i = 0; i < pts_in.size(); i++) {
			math2vvr(pts_in[i], vvr::Colour::magenta).draw();
		}
		vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;

		
	}

	//! Find and Draw Nearest Neighbour
	if (FLAG_ON(m_flag, SHOW_NN)) {
		float dist;
		dist = std::numeric_limits<float>::max();
		const KDNode *nearest = NULL;
		vec nn;
		if (FLAG_ON(m_flag, BRUTEFORCE))
		{
			dist = std::numeric_limits<float>::max();
			for (auto &d : m_KDTree->pts)
			{
				if (d.Distance(sc) < dist)
				{
					nn = d;
					dist = d.DistanceSq(sc);
				}
			}
		}
		else
		{
			Task_02_Nearest(sc, m_KDTree->root(), nearest, dist);
			if (nearest) nn = nearest->split_point;
		}
		vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
		math2vvr(sc, vvr::Colour::blue).draw();
		math2vvr(nn, vvr::Colour::green).draw();
		vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;
	}

	//! Find and Draw K Nearest Neighbour
	const float duration_knn_before = vvr::getSeconds();
	if (FLAG_ON(m_flag, SHOW_KNN)) {
		float dist;
		const KDNode **nearests = new const KDNode*[m_kn];
		memset(nearests, NULL, m_kn * sizeof(KDNode*));
		Task_04_NearestK(m_kn, sc, m_KDTree->root(), *nearests, dist);

		for (int i = 0; i < m_kn; i++) {
			if (!nearests[i]) continue;
			vec nn = nearests[i]->split_point;
			vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
			math2vvr(sc, vvr::Colour::blue).draw();
			math2vvr(nn, vvr::Colour::green).draw();
			vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;
		}
		delete[] nearests;
	}
	else if (FLAG_ON(m_flag, BRUTEFORCE)) {
		std::vector<float> dists;
		std::vector<vec> p_d;
		//const KDNode **nearests = new const KDNode*[m_kn];
		//memset(nearests, NULL, m_kn * sizeof(KDNode*));
		for (auto &d : m_KDTree->pts)
		{
			float distance = d.DistanceSq(sc);
			dists.push_back(distance);
			p_d.push_back(d);
		}
		std::vector<std::size_t> index_vec;
		assert(p_d.size() == dists.size());
		for (std::size_t i = 0; i != p_d.size(); ++i) { index_vec.push_back(i); }

		std::sort(
			index_vec.begin(), index_vec.end(),
			[&](std::size_t a, std::size_t b) { return dists[a] < dists[b]; });

	
		for (int i = 0; i < m_kn; i++) {
			vec nn = p_d[index_vec[i]];
			vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
			math2vvr(sc, vvr::Colour::blue).draw();
			math2vvr(nn, vvr::Colour::green).draw();
			vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;
		}
	}
	const float duration_knn_after = vvr::getSeconds();
	//! Draw KDTree
	if (FLAG_ON(m_flag, SHOW_KDTREE)) {
		for (int level = m_current_tree_level; level <= m_current_tree_level; level++) {
			std::vector<KDNode*> levelNodes = m_KDTree->getNodesOfLevel(level);
			for (int i = 0; i < levelNodes.size(); i++) {
				if (m_flag & FLAG(SHOW_PTS_KDTREE)) {
					VecVector pts;
					Task_01_FindPtsOfNode(levelNodes[i], pts);
					vvr::Shape::DEF_POINT_SIZE = vvr::Shape::DEF_POINT_SIZE = POINT_SIZE;
					for (int pi = 0; pi < pts.size(); pi++) {
						math2vvr(pts[pi], Pallete[i % 6]).draw();
					}
					vvr::Shape::DEF_POINT_SIZE = POINT_SIZE_SAVE;
				}
				vec c1 = levelNodes[i]->aabb.minPoint;
				vec c2 = levelNodes[i]->aabb.maxPoint;
				vvr::Box3D box(c1.x, c1.y, c1.z, c2.x, c2.y, c2.z);
				box.setTransparency(0.9);
				box.setColour(vvr::Colour::cyan);
				box.draw();
				
			}
		}
	}

	//! Compute & Display FPS
	static float last_update = 0;
	static float last_show = 0;
	const float sec = vvr::getSeconds();
	const float dt = sec - last_update;
	const float dt_show = sec - last_show;
	int FPS = 1.0 / dt;
	last_update = sec;
	if (FLAG_ON(m_flag, SHOW_FPS) && dt_show >= 2) {
		echo(FPS);
		last_show = sec;
		//std::cout << "Duration is  " << duration2-duration1 << " sec" << std::endl;
		echo(m_kn);
		std::cout << "Duration knn is  " << duration_knn_after - duration_knn_before << " sec" << std::endl;
	}
	
}

void KDTreeScene::arrowEvent(vvr::ArrowDir dir, int modif)
{
	if (dir == vvr::UP)
	{
		++m_current_tree_level;
		if (m_current_tree_level > m_KDTree->depth())
			m_current_tree_level = m_KDTree->depth();
	}
	else if (dir == vvr::DOWN)
	{
		--m_current_tree_level;
		if (m_current_tree_level < 0)
			m_current_tree_level = 0;
	}
}

void KDTreeScene::keyEvent(unsigned char key, bool up, int modif)
{
	Scene::keyEvent(key, up, modif);
	key = tolower(key);

	switch (key)
	{
		FLAG_TOGGLE(m_flag, 'b', BRUTEFORCE);
		FLAG_TOGGLE(m_flag, 'n', SHOW_NN);
		FLAG_TOGGLE(m_flag, 'k', SHOW_KNN);
		FLAG_TOGGLE(m_flag, 'f', SHOW_FPS);
		FLAG_TOGGLE(m_flag, 'a', SHOW_AXES);
		FLAG_TOGGLE(m_flag, 's', SHOW_SPHERE);
		FLAG_TOGGLE(m_flag, 't', SHOW_KDTREE);
		FLAG_TOGGLE(m_flag, 'p', SHOW_PTS_ALL);
		FLAG_TOGGLE(m_flag, 'd', SHOW_PTS_KDTREE);
		FLAG_TOGGLE(m_flag, 'c', SHOW_PTS_IN_SPHERE);
		FLAG_TOGGLE(m_flag, 'u', POINTS_ON_SURFACE);
	}

	if (key == ' ')
	{
		if (m_anim.paused()) m_anim.update(true); else m_anim.pause();
	}
	else if (key == '?')
	{
		printKeyboardShortcuts();
	}
	else if (key == 'u')
	{
		if (FLAG_ON(m_flag, POINTS_ON_SURFACE)) {
			createSurfacePts(m_pts.size());
		}
		else {
			createRandomPts(m_pts.size());
		}
		m_tree_invalidation_sec = vvr::getSeconds();
	}
}

void KDTreeScene::printKeyboardShortcuts()
{
	std::cout << "Keyboard shortcuts:"
		<< std::endl << "'?' => This shortcut list:"
		<< std::endl << "'b' => BRUTEFORCE"
		<< std::endl << "'n' => SHOW_NN"
		<< std::endl << "'k' => SHOW_KNN"
		<< std::endl << "'f' => SHOW_FPS"
		<< std::endl << "'a' => SHOW_AXES"
		<< std::endl << "'s' => SHOW_SPHERE"
		<< std::endl << "'t' => SHOW_KDTREE"
		<< std::endl << "'p' => SHOW_PTS_ALL"
		<< std::endl << "'d' => SHOW_PTS_KDTREE"
		<< std::endl << "'c' => SHOW_PTS_IN_SPHERE"
		<< std::endl << "'u' => POINTS_ON_SURFACE"
		<< std::endl << "'l' => BRUTEFORCE_2"
		<< std::endl << std::endl;
}

void KDTreeScene::mouseWheel(int dir, int modif)
{
	if (false) // Placeholder
	{
	}
	else
	{
		Scene::mouseWheel(dir, modif);
	}
}

void KDTreeScene::sliderChanged(int slider_id, float v)
{
	switch (slider_id)
	{
	case 0:
		m_anim.setSpeed(v);
		break;
	case 1:
	{
		int num_pts = NUM_PTS_DEFAULT * (SQUARE(100 * v) + 1);
		if (num_pts <= 3) num_pts = 3;
		if (FLAG_ON(m_flag, POINTS_ON_SURFACE)) {
			createSurfacePts(num_pts);
		}
		else {
			createRandomPts(num_pts);
		}
		echo(m_pts.size());
	}
	m_tree_invalidation_sec = vvr::getSeconds();
	break;
	case 2:
		m_sphere.rad = v * 30 * SPHERE_RAD;
		break;
	case 3:
		m_kn = v * 100;
		break;
	}
}

int main(int argc, char* argv[])
{
	try {
		return vvr::mainLoop(argc, argv, new KDTreeScene);
	}
	catch (std::string exc) {
		std::cerr << exc << std::endl;
		return 1;
	}
	catch (...)
	{
		std::cerr << "Unknown exception" << std::endl;
		return 1;
	}
}

//! Ground::

Ground::Ground(const float W, const float D, const float B, const float T, const vvr::Colour &colour)
	: m_col(colour)
{
	const vec vA(-W / 2, B, -D / 2);
	const vec vB(+W / 2, B, -D / 2);
	const vec vC(+W / 2, B, +D / 2);
	const vec vD(-W / 2, B, +D / 2);
	const vec vE(-W / 2, T, -D / 2);
	const vec vF(+W / 2, T, -D / 2);

	m_floor_tris.push_back(math::Triangle(vB, vA, vD));
	m_floor_tris.push_back(math::Triangle(vB, vD, vC));
	m_floor_tris.push_back(math::Triangle(vF, vE, vA));
	m_floor_tris.push_back(math::Triangle(vF, vA, vB));
}

void Ground::draw() const
{
	for (int i = 0; i < m_floor_tris.size(); i++)
	{
		vvr::Triangle3D floor_tri = vvr::math2vvr(m_floor_tris.at(i), m_col);
		floor_tri.setSolidRender(true);
		floor_tri.draw();
	}
}

//! KDTree::

KDTree::KDTree(VecVector &pts)
	: pts(pts)
{
	const float t = vvr::getSeconds();
	m_root = new KDNode();
	m_depth = makeNode(m_root, pts, 0);
	const float KDTree_construction_time = vvr::getSeconds() - t;
	echo(KDTree_construction_time);
	echo(m_depth);
}

KDTree::~KDTree()
{
	const float t = vvr::getSeconds();
	delete m_root;
	const float KDTree_destruction_time = vvr::getSeconds() - t;
	echo(KDTree_destruction_time);
}

int KDTree::makeNode(KDNode *node, VecVector &pts, const int level)
{
	//! Sort along the appropriate axis, find median point and split.
	const int axis = level % DIMENSIONS;
	std::sort(pts.begin(), pts.end(), VecComparator(axis));
	const int i_median = pts.size() / 2;

	//! Set node members
	node->level = level;
	node->axis = axis;
	node->split_point = pts[i_median];
	node->aabb.SetFrom(&pts[0], pts.size());

	//! Continue recursively or stop.
	if (pts.size() <= 1)
	{
		return level;
	}
	else
	{
		int level_left = 0;
		int level_right = 0;

		VecVector pts_left(pts.begin(), pts.begin() + i_median);
		VecVector pts_right(pts.begin() + i_median + 1, pts.end());

		if (!pts_left.empty())
		{
			node->child_left = new KDNode();
			level_left = makeNode(node->child_left, pts_left, level + 1);

		}
		if (!pts_right.empty())
		{
			node->child_right = new KDNode();
			level_right = makeNode(node->child_right, pts_right, level + 1);
		}

		int max_level = std::max(level_left, level_right);
		return max_level;
	}
}

void KDTree::getNodesOfLevel(KDNode *node, std::vector<KDNode*> &nodes, int level)
{
	if (!level)
	{
		nodes.push_back(node);
	}
	else
	{
		if (node->child_left) getNodesOfLevel(node->child_left, nodes, level - 1);
		if (node->child_right) getNodesOfLevel(node->child_right, nodes, level - 1);
	}
}

std::vector<KDNode*> KDTree::getNodesOfLevel(const int level)
{
	std::vector<KDNode*> nodes;
	if (!m_root) return nodes;
	getNodesOfLevel(m_root, nodes, level);
	return nodes;
}

//! Tasks

void Task_01_FindPtsOfNode(const KDNode* root, VecVector &pts)
{
	// Push all points under "root" to "pts" by traversing the KD-tree
	//...
	//pts.push_back(root->split_point);

	if (root->child_left) {
		Task_01_FindPtsOfNode(root->child_left, pts);
	}
	if (root->child_right) {
		Task_01_FindPtsOfNode(root->child_right, pts);
	}



}

void Task_02_Nearest(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist)
{
	// Find the Nearest Neighbour to "test_pt" and store it to "nn"
	// Traverse the tree from "root" down to search the best distance to "test_pt"
	// and then traverse the tree bottom up to search for alternatives,
	// store the best distance to "best_dist"
	// ...
	
	int axis = root->axis;

	//float d = sqrt(pow(test_pt.x - root->split_point.x, 2) + pow(test_pt.y - root->split_point.y, 2) + pow(test_pt.z - root->split_point.z, 2));
	float d = test_pt.DistanceSq(root->split_point);
	if (d < best_dist) {
		best_dist = d;
		nn = root;
	}

	float distance = test_pt.ptr()[axis] - root->split_point.ptr()[axis];

	if (distance < 0) {
		if (root->child_left) {
			Task_02_Nearest(test_pt, root->child_left, nn, best_dist);
		}

	}
	else {
		if (root->child_right) {
			Task_02_Nearest(test_pt, root->child_right, nn, best_dist);
		}
	}

	if (distance > 0) {
		if (root->child_left) {
			Task_02_Nearest(test_pt, root->child_left, nn, best_dist);
		}

	}
	else {
		if (root->child_right) {
			Task_02_Nearest(test_pt, root->child_right, nn, best_dist);
		}
	}
	
}

void Task_03_InSphere(const Sphere &sphere, const KDNode *root, VecVector &pts)
{
	// Find all the points of the KD-tree that are inside "sphere" and store them to "pts"
	//...
	float r = sphere.Diameter() / 2;

	//sphere.Contains(root->split_point)
	if (r > sphere.Centroid().Distance(root->split_point)) {
		pts.push_back(root->split_point);
	}

	if (root->child_left) {
		Task_03_InSphere(sphere, root->child_left, pts);
	}
	if (root->child_right) {
		Task_03_InSphere(sphere, root->child_right, pts);
	}
	
}

void Task_04_NearestK(const int k, const vec& test_pt, const KDNode* root, const KDNode* &knn, float& best_dist)
{
	// Find the K-nearest neighbours of "test_pt" and store them to "knn"
	float old_distance;
	old_distance = 0;
	
	for (int i = 0; i < k; i++)
	{
		const KDNode *nearest = root;
		float dist;
		dist = std::numeric_limits<float>::max();

		Nearest(test_pt, root, nearest, dist, old_distance);

		*(&knn + i) = nearest;		// store it in knn array

		old_distance = dist;
	}
}

void Nearest(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist, float& old_distance)
{
	int axis = root->axis;

	//float d = sqrt(pow(test_pt.x - root->split_point.x, 2) + pow(test_pt.y - root->split_point.y, 2) + pow(test_pt.z - root->split_point.z, 2));
	float d = test_pt.DistanceSq(root->split_point);
	if (d < best_dist) {
		if (d > old_distance) {
			best_dist = d;
			nn = root;
		}
	}

	float distance = test_pt.ptr()[axis] - root->split_point.ptr()[axis];

	if (distance < 0) {
		if (root->child_left) {
			Nearest(test_pt, root->child_left, nn, best_dist, old_distance);
		}

	}
	else {
		if (root->child_right) {
			Nearest(test_pt, root->child_right, nn, best_dist, old_distance);
		}
	}

	if (distance > 0) {
		if (root->child_left) {
			Nearest(test_pt, root->child_left, nn, best_dist, old_distance);
		}

	}
	else {
		if (root->child_right) {
			Nearest(test_pt, root->child_right, nn, best_dist, old_distance);
		}
	}

}
