//
// Created by ylzha on 2025/2/25.
//
#ifndef slic3r_FillBridge_hpp_
#define slic3r_FillBridge_hpp_

#include "../libslic3r.h"
#include <boost/geometry.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "FillBase.hpp"


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


namespace Slic3r{

class Surface;

// 定义点类型
typedef bg::model::d2::point_xy<double> Point_t;

//定义线段
typedef bg::model::segment<Point_t> Segment;

// 定义多边形类型
typedef bg::model::polygon<Point_t> Polygon_t;

// 多边形集合
typedef bg::model::multi_polygon<Polygon_t> MultiPolygon;

// 定义环类型
typedef bg::model::ring<Point_t> Ring;
typedef bg::model::ring<Point_t> Ring_t;

using BBox = bg::model::box<Point_t>;

// 定义折线类型（用于表示直线）
typedef bg::model::linestring<Point_t> Linestring;
typedef bg::model::box<Point_t> Box_t;

struct EdgeProperty {
    bool is_temp = false;  // 是否为临时添加的边
    bool visited = false;  // 边是否已被访问
};

// 2. 图结构定义（Boost.Graph）
typedef boost::adjacency_list<
    boost::vecS,               // 边存储：动态数组
    boost::vecS,               // 顶点存储：动态数组
    boost::undirectedS,        // 无向图（桥接点连接无方向）
    Point_t,                   // 顶点属性：存储 Point_t（坐标）
    EdgeProperty               // 边属性：包含访问标记
> ShapeGraph;

typedef boost::graph_traits<ShapeGraph>::edge_iterator EdgeIter;  // 边迭代器
typedef boost::graph_traits<ShapeGraph>::vertex_iterator VertexIter; // 顶点迭代器
typedef boost::graph_traits<ShapeGraph>::vertex_descriptor Vertex; // 顶点句柄
typedef boost::graph_traits<ShapeGraph>::edge_descriptor   Edge;   // 边句柄
typedef boost::graph_traits<ShapeGraph>::out_edge_iterator OutEdgeIter;

// 定义路径类型（顶点序列）
typedef std::vector<Vertex> Path;

struct PointCompare {
    bool operator()(const Point_t& a, const Point_t& b) const {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    }
};

struct MidPoints {
    size_t index; //边索引
    Point_t mid_point; //中点
    double distance; //某点离边的距离
    MidPoints(size_t index, Point_t mid_point, double distance) :index(index)
        , mid_point(mid_point), distance(distance) {
    }
};

struct IdIndex {
    size_t id{};
    size_t index{};
    // 定义小于运算符
    bool operator<(const IdIndex& other) const {
        return id < other.id || (id == other.id && index < other.index);
    }
    IdIndex() = default;
    IdIndex(const size_t id, const size_t index) : id(id), index(index) {}

    bool operator==(const IdIndex& ii) const {
        return id == ii.id && index == ii.index;
    }
};

//环路径树结构中节点类型
struct RingNode {
    IdIndex id;
    Ring ring; //当前环
    int orientation; //方向 -1 向内 1 向外
    std::vector<RingNode> children; //子节点集
    RingNode* parent; //父节点
    bool isHide{ false }; //在最终路径中是否呈现 false 呈现  true 不呈现
    // 构造函数，方便初始化
    RingNode(const IdIndex id, const Ring& ring, const int orientation = 0)
        : id(id), ring(ring), orientation(orientation), parent(nullptr) {
    }

    // 重载==运算符方便比较
    bool operator==(const RingNode& other) const {
        return id.id == other.id.id && id.index == other.id.index;
    }

};

//桥接映射
struct BridgeMap {
    Point_t from;
    Point_t to;
    Point_t from2;
    Point_t to2;
    IdIndex from_ii;
    IdIndex to_ii;
    // 默认构造函数
    BridgeMap()
        : from(0, 0), to(0, 0), from2(0, 0), to2(0, 0),  // 给 Point_t 传默认参数
        from_ii({}), to_ii({}) {
    }                       

    BridgeMap(Point_t from, Point_t to, 
        Point_t from2, Point_t to2, IdIndex from_ii
        , IdIndex to_ii
    ) :
        from(from), to(to), 
        from2(from2), to2(to2), 
        from_ii(from_ii), to_ii(to_ii){
    }
};

//合并映射
struct MergeMap {
    IdIndex ii1; //环1
    IdIndex ii2; //环2
    std::vector<IdIndex> nodes; //产生的环集
    // 构造函数
    MergeMap(const IdIndex& i1, const IdIndex& i2, const std::vector<IdIndex>& nodeList)
        : ii1(i1), ii2(i2), nodes(nodeList) {
    }
};

struct MergeMap2 {
    IdIndex outer_ii; //外多边形
    std::vector<IdIndex> inner_iis; //多个内多边形
    std::vector<IdIndex> merged_iis; //合并产生的多边形
    MergeMap2(const IdIndex& outer_ii, const std::vector<IdIndex>& inner_iis,
        const std::vector<IdIndex> merged_iis)
        :outer_ii(outer_ii), inner_iis(inner_iis), merged_iis(merged_iis) {
    }
};

class FillBridge : public Fill {

public:
    ~FillBridge() override = default;
    bool is_self_crossing() override { return false; }
    Fill* clone() const override { return new FillBridge(*this); }
    //Polylines fill_surface(const Surface* surface, const FillParams& params) override;
    void _fill_surface_single(
        const FillParams& params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>& direction,
        ExPolygon                        expolygon,
        Polylines& polylines_out);

private:


    std::map<size_t, std::vector<RingNode>> ringNodes;
    std::vector<BridgeMap> bridges;

    std::map<IdIndex, std::vector<IdIndex>> offsetMap;
    std::vector<MergeMap> containMap;
    std::vector<MergeMap2> mergeMap2;
    size_t maxRid = 1;
    bool isFront = false;
    double _offset = 40000; //偏移距离
    double _line_spacing = 0;
    double area_threshold = sqr(scaled<double>(0.001)); //面积阈值

    Polygon_t b_polygon;
    std::vector<Polygon_t> o_polygons;

private:

    // 判断两个环是否相交（不包括包含关系）
    bool polyIntersect(const Polygon_t poly1, const Polygon_t poly2);

    //获取环偏移后的环
    std::vector<Ring> offsetRing(const Ring& ring,
        double distance,
        double area_threshold
    );

    std::vector<Ring> getFallbackRings(
        const Ring& original_ring,
        double original_area,
        const MultiPolygon& intersectionMP
    );

    //生成环集
    void generateRings();
    //形成节点
    RingNode formatNode(size_t id, size_t index, const Ring& ring, int orientation);
    void addNode(RingNode& node);
    void removeNode(IdIndex id);
    //查找节点
    RingNode& findNode(IdIndex ii);

    //形成树
    void formatTree();
    // 深度优先搜索遍历树
    void dfs(RingNode& node, std::vector<IdIndex>& visited);
    Polygon_t safe_union(const Polygon_t& a, const Polygon_t& b);

    // 在环上按顺时针方向查找距离给定点d的另一个点
    Point_t find_point_at_distance_clockwise(Ring& ring, const Point_t& start_point,
        size_t _index, double d, size_t& e_index);

    bool equal(Point_t p1, Point_t p2);


    void handleBridge(IdIndex o_ii, IdIndex i_ii);



    std::vector<RingNode> compute_complex_polygon_merge(
        const std::vector<RingNode>& outers
        , const std::vector<RingNode>& inners
    );
    bool isPolygonContained(const Polygon_t& polyA, const Polygon_t& polyB);

    Point_t closest_point_on_segment(const Point_t& p, const Point_t& seg_start, const Point_t& seg_end);
    Point_t find_closest_point_on_ring_edges(Ring& ring, const Point_t& p0, size_t& e_index1);
    int findPointIndex(const Ring& ring, const Point_t& p0);
    bool does_segment_cross_ring(const Segment& seg, const Ring& ring);


    auto get_bbox(const Polygon_t& poly);
    void polygonsMerge(
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inputs,
        std::vector<std::pair<Polygon_t, IdIndex>>& outputs, size_t orientation);


    //环所有边按照边长遍历边的中点
    std::vector<MidPoints> find_mid_points_on_ring(const Ring& ring, double& perimeter) {
        // 存储边长度和索引的结构
        std::vector<MidPoints> points;
        perimeter = 0;
        // 遍历环的边（注意闭合环最后一点与起点相连）
        for (size_t i = 0; i < ring.size(); ++i) {
            const Point_t& p1 = ring[i];
            const Point_t& p2 = ring[(i + 1) % ring.size()]; // 处理闭合边
            // 计算线段中点
            Point_t midpoint(
                (bg::get<0>(p1) + bg::get<0>(p2)) / 2,
                (bg::get<1>(p1) + bg::get<1>(p2)) / 2
            );
            double len = bg::distance(p1, p2);
            MidPoints m(i, midpoint, len);
            points.push_back(m);
            perimeter += len;
        }
        //根据距离降序排序
        std::sort(points.begin(), points.end(),
            [](const MidPoints& a, const MidPoints& b) {
                return a.distance > b.distance; // 降序：a的distance大于b时，a排在前面
            });
        return points;

    }

    //距离点point最远的长边中点
    void find_further_point(std::vector<MidPoints> midpoints,
        Point_t point, size_t _index, Point_t& point2, size_t& point2_index) {
        double max_distance = 0.0;
        for (auto& mp : midpoints) {
            double d = bg::distance(mp.mid_point, point);
            if (max_distance < d && mp.index != _index) {
                point2 = mp.mid_point;
                point2_index = mp.index;
                max_distance = d;
            }
        }
    }


    void insertPointIntoRing(Ring& r, const Point_t& p0) {
        if (r.empty()) {
            r.push_back(p0);
            return;
        }

        if (r.size() == 1) {
            r.push_back(p0);
            return;
        }

        std::size_t n = r.size();
        double min_distance = std::numeric_limits<double>::max();
        std::size_t insert_index = 0;

        for (std::size_t i = 0; i < n; ++i) {
            const Point_t& p1 = r[i];
            const Point_t& p2 = r[(i + 1) % n];

            // 创建线段
            Segment seg(p1, p2);

            // 计算点到线段的距离
            double dist = bg::distance(p0, seg);

            if (dist < min_distance) {
                min_distance = dist;
                insert_index = i + 1;
            }
        }

        // 在找到的位置插入点
        if (insert_index == n) {
            r.push_back(p0);
        }
        else {
            r.insert(r.begin() + insert_index, p0);
        }
    }


    // 修复后的凸包计算函数
    Polygon_t computeConvexHull(const std::vector<Point_t>& points) {
        bg::model::linestring<Point_t> line;
        for (const auto& point : points) {
            bg::append(line, point);
        }

        Polygon_t hull;
        bg::convex_hull(line, hull);
        return hull;
    }



    // 计算向量叉积：(b - a) × (c - b)
    double cross_product(const Point_t& a, const Point_t& b, const Point_t& c) {
        double dx1 = b.x() - a.x();
        double dy1 = b.y() - a.y();
        double dx2 = c.x() - b.x();
        double dy2 = c.y() - b.y();
        return dx1 * dy2 - dy1 * dx2; // 叉积结果：正→逆时针，负→顺时针，0→共线
    }

    // 手动判断多边形是否为凸多边形
    bool is_convex(const Polygon_t& poly) {
        const auto& ring = poly.outer(); // 取多边形外环
        if (ring.size() < 4) return true; // 3个点以下视为凸（退化多边形）

        // 计算所有连续边的叉积符号（忽略共线点的0值）
        int sign = 0;
        for (size_t i = 0; i < ring.size() - 1; ++i) {
            const Point_t& a = ring[i];
            const Point_t& b = ring[i + 1];
            const Point_t& c = ring[(i + 2) % (ring.size() - 1)]; // 注意：外环最后一点与第一点重复，需取模

            double cross = cross_product(a, b, c);
            if (std::abs(cross) < 1e-9) continue; // 共线点跳过

            int current_sign = cross > 0 ? 1 : -1;
            if (sign == 0) {
                sign = current_sign; // 初始化符号
            }
            else if (current_sign != sign) {
                return false; // 符号不一致，为凹多边形
            }
        }
        return true; // 所有转向符号一致，为凸多边形
    }

    // 计算三点叉积（用于判断转向）
    double crossProduct(const Point_t& a, const Point_t& b, const Point_t& c) {
        return (b.x() - a.x()) * (c.y() - b.y()) -
            (b.y() - a.y()) * (c.x() - b.x());
    }

    // 方法3：手动检查凹性（基于叉积符号变化）
    bool isConcaveManual(const Polygon_t& polygon) {
        const auto& outer_ring = polygon.outer();
        int n = outer_ring.size();

        // 三角形总是凸的
        if (n <= 4) return false; // 3个点+1个重复的闭合点

        int sign = 0;

        for (int i = 0; i < n - 1; ++i) {
            const Point_t& a = outer_ring[i];
            const Point_t& b = outer_ring[(i + 1) % n];
            const Point_t& c = outer_ring[(i + 2) % n];

            double cross = crossProduct(a, b, c);

            // 忽略接近零的叉积（共线情况）
            if (std::abs(cross) > 1e-10) {
                if (sign == 0) {
                    // 确定初始符号
                    sign = (cross > 0) ? 1 : -1;
                }
                else {
                    // 检查符号是否变化
                    if ((cross > 0 && sign < 0) || (cross < 0 && sign > 0)) {
                        return true; // 发现凹顶点
                    }
                }
            }
        }

        return false;
    }
   
    // 使用凸包方法提取分离的多边形
    std::vector<Polygon_t> extractPolygonsWithoutThreshold(
        const std::vector<Point_t>& all_points, Polygon_t& polygon) {
        std::vector<Polygon_t> polygons;

        // 如果点集为空，直接返回
        if (all_points.empty()) {
            return polygons;
        }
        Polygon_t _poly;
        if (!isConcaveManual(polygon)) {
            // 计算整个点集的凸包
            _poly = computeConvexHull(all_points);
        }
        else {
            bg::append(_poly.outer(), all_points);
            polygons.push_back(_poly);
            return polygons;
        }

        // 如果凸包内没有孔，说明只有一个多边形
        if (_poly.inners().empty()) {
            polygons.push_back(_poly);
            return polygons;
        }

        // 如果有孔，则每个孔对应一个分离的多边形
        // 主多边形是第一个分离的多边形
        polygons.push_back(_poly);

        // 处理每个孔（分离的多边形）
        for (const auto& inner_ring : _poly.inners()) {
            std::vector<Point_t> inner_points(inner_ring.begin(), inner_ring.end());
            Polygon_t inner_poly = computeConvexHull(inner_points);
            polygons.push_back(inner_poly);
        }
        return polygons;
    }



    // 在向量中查找点的索引
    size_t findPointIndex(const std::vector<Point_t>& points, const Point_t& target) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (equal(points[i], target)) {
                return i;
            }
        }
        return points.size(); // 返回无效索引
    }



    // 查找欧拉回路的函数（基于Hierholzer算法）
    std::vector<Vertex> findEulerCircuit(ShapeGraph& graph) {
        std::vector<Vertex> circuit;
        std::stack<Vertex> vertexStack;

        // 选择起始顶点（任意顶点）
        Vertex start = *boost::vertices(graph).first;
        vertexStack.push(start);

        while (!vertexStack.empty()) {
            Vertex v = vertexStack.top();
            bool hasUnvisitedEdge = false;

            // 查找未访问的边
            std::pair<OutEdgeIter, OutEdgeIter> edges = boost::out_edges(v, graph);
            for (OutEdgeIter eit = edges.first; eit != edges.second; ++eit) {
                Edge e = *eit;
                if (!graph[e].visited) {
                    graph[e].visited = true;  // 标记为已访问
                    hasUnvisitedEdge = true;

                    // 获取相邻顶点
                    Vertex u = boost::target(e, graph);
                    if (u == v) u = boost::source(e, graph);

                    vertexStack.push(u);
                    break;
                }
            }

            if (!hasUnvisitedEdge) {
                vertexStack.pop();
                circuit.push_back(v);
            }
        }

        // 反转得到正确的路径顺序
        std::reverse(circuit.begin(), circuit.end());
        return circuit;
    }

    // 环形 vector：精准处理 v1 和 v2 的前后关系，删除较短路径的中间顶点
    void erase_between_vertices_ring(ShapeGraph& graph, std::vector<Vertex>& vertices, Vertex v1, Vertex v2) {
        std::vector<size_t> _vertex_index_list;
        size_t h, k;
        size_t v_count = vertices.size();
        for (size_t i = 0; i < v_count; i++) {
            if (equal(graph[v1], graph[vertices[i]])) {
                k = i;
            }
            else if (equal(graph[v2], graph[vertices[i]])) {
                h = i;
            }
        }
        size_t start = k < h ? k : h;
        size_t end = k > h ? k : h;
        if (end - start < v_count - end + start) {
            for (size_t j = start + 1; j < end; j++) {
                _vertex_index_list.push_back(j);
            }
        }
        else {
            for (size_t j = 0; j < start; j++) {
                _vertex_index_list.push_back(j);
            }
            for (size_t j = end + 1; j < v_count; j++) {
                _vertex_index_list.push_back(j);
            }
        }

        // 将删除索引转换为集合便于快速查找
        std::unordered_set<size_t> remove_set(_vertex_index_list.begin(),
            _vertex_index_list.end());

        // 创建新向量，只保留不需要删除的顶点
        std::vector<Vertex> new_vertices;
        new_vertices.reserve(vertices.size() - remove_set.size());

        for (size_t i = 0; i < vertices.size(); ++i) {
            if (remove_set.find(i) == remove_set.end()) {
                new_vertices.push_back(vertices[i]);
            }
        }

        // 交换内容
        vertices.swap(new_vertices);
    }

    // 删除两个顶点之间的边
    bool remove_edge_between_vertices(ShapeGraph& graph, Vertex v1, Vertex v2) {
        // 调用 BGL 的 remove_edge 函数，返回是否成功删除边
        // 对于无向图，v1 和 v2 的顺序不影响结果
        std::pair<boost::graph_traits<ShapeGraph>::edge_descriptor, bool> result =
            boost::edge(v1, v2, graph);

        if (result.second) {
            // 边存在，执行删除
            boost::remove_edge(result.first, graph);
            return true; // 删除成功
        }
        else {
            // 边不存在，无需删除
            return false; // 删除失败（边不存在）
        }
    }

    Vertex findVertexByPoint(const ShapeGraph& graph, const Point_t& from,
        const std::vector<Vertex>& from_ring_vertex) {
        auto it = std::find_if(from_ring_vertex.begin(), from_ring_vertex.end(),
            [&graph, &from](Vertex v) {
                return graph[v].x() == from.x() && graph[v].y() == from.y();
            });
        return *it;
    }

    
};


}; // namespace Slic3r
#endif // slic3r_FillBridge_hpp_
