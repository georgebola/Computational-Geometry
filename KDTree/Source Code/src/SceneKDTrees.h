#include <VVRScene/scene.h>

//! MACROS used for toggling and testing bitwise flags.

#define FLAG(x) (1<<(x))
#define FLAG_ON(v,f) (v & FLAG(f))
#define FLAG_TOGGLE(v,c,f) case c: v ^= FLAG(f); std::cout << #f << " = " << (FLAG_ON(v,f) ? "ON" : "OFF") << std::endl; break

/**
 * To save typing.
 */
typedef std::vector<vec> VecVector;

/**
 * Array with 6 predefined colours.
 */
static const vvr::Colour Pallete[6] = {
    vvr::Colour::red, vvr::Colour::green, vvr::Colour::blue, vvr::Colour::magenta,
    vvr::Colour::orange, vvr::Colour::yellow,
};

/**
 * A node of a KD-Tree
 */
struct KDNode
{
    vec split_point;
    int axis;
    int level;
    AABB aabb;
    KDNode *child_left;
    KDNode *child_right;
    KDNode() : child_left(NULL), child_right(NULL) {}
    ~KDNode() { delete child_left; delete child_right; }
};

/**
 * KD-Tree wrapper. Holds a ptr to tree root.
 */
class KDTree
{
public:
    KDTree(VecVector &pts);
    ~KDTree();
    std::vector<KDNode*> getNodesOfLevel(int level);
    int depth() const { return m_depth; }
    const KDNode* root() const { return m_root; }
    const VecVector &pts;

private:
    static int makeNode(KDNode *node, VecVector &pts, const int level);
    static void getNodesOfLevel(KDNode *node, std::vector<KDNode*> &nodes, int level);

private:
    KDNode *m_root;
    int m_depth;
};

/**
 * Scene
 */
class KDTreeScene : public vvr::Scene
{
    enum {
        BRUTEFORCE, POINTS_ON_SURFACE, SHOW_AXES, SHOW_FPS, SHOW_NN, SHOW_KNN, SHOW_SPHERE,
        SHOW_KDTREE, SHOW_PTS_ALL, SHOW_PTS_KDTREE,
        SHOW_PTS_IN_SPHERE
    };

public:
    KDTreeScene();
    const char* getName() const { return "KD Tree Scene"; }
    void keyEvent(unsigned char key, bool up, int modif) override;
    void arrowEvent(vvr::ArrowDir dir, int modif) override;
    void mouseWheel(int dir, int modif) override;
    void sliderChanged(int slider_id, float val);


private:
    void draw() override;
    void reset() override;
    void resize() override;
    bool idle() override;
    void createRandomPts(int num_pts);
    void createSurfacePts(int num_pts);
    void printKeyboardShortcuts();

private:
    KDTree *m_KDTree;
    VecVector m_pts;
    vvr::Sphere3D m_sphere;
    vvr::Animation m_anim;
    int m_flag;
    math::LCG m_lcg;
    int m_current_tree_level;
    float m_tree_invalidation_sec;
    int m_kn;
};

/**
 * Class Representing scene floor and background wall.
 */
class Ground : public vvr::IRenderable
{
public:
    Ground(const float W, const float D, const float B, const float T, const vvr::Colour &colour);
    void draw() const override;

private:
    std::vector<math::Triangle> m_floor_tris;
    vvr::Colour m_col;
};

/**
 * Function object to compare 2 3D-vecs in the specified axis.
 */
struct VecComparator {
    unsigned axis;
    VecComparator(unsigned axis) : axis(axis % 3) {}
    virtual inline bool operator() (const vec& v1, const vec& v2) {
        return (v1.ptr()[axis] < v2.ptr()[axis]);
    }
};

//! Task function prototypes

/**
 * Find all the points under `root` node of the tree.
 */
void Task_01_FindPtsOfNode(const KDNode* root, VecVector &pts);

/**
 * Find the nearest neighbour of `test_pt` inside `root`.
 */
void Task_02_Nearest(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist);

/**
 * Find the points of `kdtree` that are contained inside `sphere`.
 */
void Task_03_InSphere(const Sphere &sphere, const KDNode *root, VecVector &pts);

/**
 * Find the `k` nearest neighbours of `test_pt` inside `root`.
 */
void Task_04_NearestK(const int k, const vec& test_pt, const KDNode* root, const KDNode* &knn, float& best_dist);


void Nearest(const vec& test_pt, const KDNode* root, const KDNode* &nn, float& best_dist,float& old_distance);