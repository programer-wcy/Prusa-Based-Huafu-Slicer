#include "../ClipperUtils.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"
#include "../ShortestPath.hpp"
#include <windows.h>
#include <sstream>
#include <string>
#include <filesystem>

#include "FillBridge.hpp"

namespace Slic3r {

    Ring polygonToRing(const Polygon polygon) {
        Points points = polygon.points;
        Point first = points[0];
        points.push_back(first); // 闭合
        Ring inner_ring;
        inner_ring.reserve(points.size());
        for (const Point& p : points) {
            inner_ring.emplace_back(p.x(), p.y());
        }
        return inner_ring;
    }

    void FillBridge::_fill_surface_single(
        const FillParams& params,
        unsigned int                     thickness_layers,
        const std::pair<float, Point>& direction,
        ExPolygon                        expolygon,
        Polylines& polylines_out)
    {
        double area = expolygon.contour.area();
        if (area < 9000000000) {
            return;// polylines_out;
        }

        BoundingBox bounding_box = expolygon.contour.bounding_box();

        coord_t min_spacing = scale_(this->spacing);
        coord_t distance = coord_t(min_spacing / params.density);
        if (params.density > 0.9999f && !params.dont_adjust) {
            distance = this->_adjust_solid_spacing(bounding_box.size()(0), distance);
            this->spacing = unscale<double>(distance);
            min_spacing = scale_(this->spacing);
        }

        // commented by wangcy
        // _offset = min_spacing
        _offset = 1000000*(this->print_config->carbon_fiber_parallel_line_spacing.value); //added by wangcy


        o_polygons.clear();
        b_polygon.inners().clear();

        //构造外边界（闭合）
        Points outer_points = expolygon.contour.points;
        if (!outer_points.empty()) {
            Point first = outer_points[0];
            outer_points.push_back(first); // 闭合

            std::vector<Point_t> outer_boundary;
            outer_boundary.reserve(outer_points.size());
            for (const Point& p : outer_points) {
                outer_boundary.emplace_back(p.x(), p.y());
            }

            Polygon_t poly;
            bg::append(poly.outer(), outer_boundary);

            // 提取多边形
            o_polygons = extractPolygonsWithoutThreshold(outer_boundary, poly);
        }

        //  处理内孔（每个内孔都作为独立的 inner ring）
        for (const Polygon& hole : expolygon.holes) {
            Ring inner_ring = polygonToRing(hole);
            if (inner_ring.size() > 0) {
                b_polygon.inners().push_back(inner_ring);
            }
        }


        //初始化
        maxRid = 1;
        ringNodes.clear();
        offsetMap.clear();
        containMap.clear();
        mergeMap2.clear();

        generateRings();
        formatTree();
        std::vector<IdIndex> iis;
        for (auto& rn : ringNodes) {
            for (auto& node : rn.second) {
                if (node.parent == nullptr && !node.isHide) {
                    iis.push_back(node.id);
                }
            }
        }
        for (auto& id : iis) {
            bridges.clear();
            path.clear();
            all_rings.clear();
            std::vector<IdIndex> visited;
            dfs(findNode(id), visited);
            std::reverse(bridges.begin(), bridges.end());
            RingNode& root = findNode(id);
            if (bridges.size() > 0) {
                Point_t start = bridges[0].from;   // 第一个桥接点
                Point_t end = bridges[0].from2;  // 第二个桥接点
                traverseRing(root, id, start, end, true);
            }
            else {
                path.insert(path.begin(), root.ring.begin(), root.ring.end());
            }
            // 转为 Polylines
            for (size_t i = 0; i < path.size(); ++i) {
                size_t j = (i + 1) % path.size();
                polylines_out.push_back({ { path[i].x(),path[i].y()},
                    { path[j].x() ,path[j].y() } });
            }
        }
    }


    //判断两个点是否相等
    bool FillBridge::equal(Point_t p1, Point_t p2) {
        const double eps = 1e-6;
        return std::abs(p1.x() - p2.x()) < eps && std::abs(p1.y() - p2.y()) < eps;
    }

    //将节点加入处理节点集合中
    RingNode FillBridge::formatNode(const size_t id, const size_t index, const Ring& ring, const int orientation) {
        Ring ringNode;
        ringNode.assign(ring.begin(), ring.end());
        RingNode node({ id, index }, ringNode, orientation);
        return node;
    }
    void FillBridge::addNode(RingNode& node) {

        // 检查键是否存在
        const auto it = ringNodes.find(node.id.id);
        if (it != ringNodes.end()) {
            // 键已存在，将新节点添加到对应向量的末尾
            it->second.emplace_back(node);
        }
        else {
            // 键不存在，创建一个新的向量并插入到 map 中
            ringNodes[node.id.id] = { node };
        }
    }

    void FillBridge::removeNode(IdIndex id) {
        auto it = ringNodes.find(id.id);
        if (it != ringNodes.end()) {
            std::vector<RingNode>& nodes = it->second;
            nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                [id](RingNode x) { return x.id.id == id.id && x.id.index == id.index; }),
                nodes.end());
        }
    }

    //查找节点
    RingNode& FillBridge::findNode(const IdIndex ii) {
        const auto it = ringNodes.find(ii.id);
        std::vector<RingNode>& nodes = it->second;
        const auto it1 = std::find_if(nodes.begin(), nodes.end(), [ii](const RingNode& node) {
            return node.id.index == ii.index;
            });
        return *it1;
    }
    // 判断多边形A是否包含多边形B
    bool FillBridge::isPolygonContained(const Polygon_t& polyA, const Polygon_t& polyB) {

        // 检查B的所有顶点是否在A内部且不在A的内环上
        const auto& a_outer = bg::exterior_ring(polyA);
        const auto& b_outer = bg::exterior_ring(polyB);

        for (const auto& p : b_outer) {
            if (!bg::within(p, polyA)) return false;
            for (const auto& a_inner : polyA.inners()) {
                if (bg::within(p, a_inner)) return false;
            }
        }

        // 检查A与B的交集是否等于B的面积（确保完全包含）
        MultiPolygon intersection;
        bg::intersection(polyA, polyB, intersection);
        double b_area = bg::area(polyB);
        double inter_area = 0.0;
        for (const auto& p : intersection) inter_area += bg::area(p);
        return std::abs(inter_area - b_area) < 1e-6;

    }

    bool FillBridge::polyIntersect(const Polygon_t poly1, const Polygon_t poly2)
    {
        // 首先检查是否相交（包括边界相交）
        if (!bg::intersects(poly1, poly2)) {
            return false;
        }

        if (bg::area(poly1) == bg::area(poly2)) {
            return false;
        }
        // 如果存在包含关系，则不算作相交
        if (isPolygonContained(poly1, poly2) || isPolygonContained(poly2, poly1)) {
            return false;
        }

        return true;
    }


    //获取环偏移后的环
    std::vector<Ring> FillBridge::offsetRing(
        const Ring& ring,
        const double distance,
        double area_threshold
    ) {

        Polygon_t inputPolygon;
        inputPolygon.outer() = ring;
        bg::correct(inputPolygon);

        double original_area = std::abs(bg::area(ring));
        if (original_area < area_threshold) {
            return {};
        }

        MultiPolygon offsettedMP; // 偏移结果

        auto distanceStrategy = bg::strategy::buffer::distance_symmetric<double>(distance);
        auto sideStrategy = bg::strategy::buffer::side_straight();
        auto joinStrategy = bg::strategy::buffer::join_miter(10); // 斜接策略（防尖角）
        auto endStrategy = bg::strategy::buffer::end_flat();
        auto pointStrategy = bg::strategy::buffer::point_circle(8); // 点拟合精度

        boost::geometry::buffer(
            inputPolygon, offsettedMP,
            distanceStrategy, sideStrategy,
            joinStrategy, endStrategy, pointStrategy
        );
        // 处理偏移后为空的情况
        if (offsettedMP.empty()) {
            return {};
        }

        // 计算原多边形与偏移后几何的交集（得到相交区域）
        MultiPolygon intersectionMP;
        bg::intersection(inputPolygon, offsettedMP, intersectionMP);

        // 从相交区域提取所有有效环
        std::vector<Ring> resultRings;
        //if (distance < 0) {
        //    for (const auto& poly : offsettedMP) {
        //        const Ring& outer_ring = poly.outer();
        //        double offset_area = std::abs(bg::area(outer_ring));
        //        if (
        //            offset_area > 0.04 * original_area &&  // 最小面积限制offset_area > area_threshold &&
        //            offset_area < 5 * original_area) {   // 最大面积限制
        //            resultRings.push_back(outer_ring);
        //        }
        //    }
        //}

        bool isI = false;
        for (const auto& poly : offsettedMP) {
            if (polyIntersect(inputPolygon, poly)) {
                isI = true;
                break;
            }
        }
        if (!isI) {
            for (const auto& poly : offsettedMP) {
                const Ring& outer_ring = poly.outer();
                double offset_area = std::abs(bg::area(outer_ring));
                if (
                    offset_area > 0.03 * original_area &&  // 最小面积限制offset_area > area_threshold &&
                    offset_area < 5 * original_area) {   // 最大面积限制
                    resultRings.push_back(outer_ring);
                }
            }
        }
        else {
            for (const auto& poly : intersectionMP) {
                const Ring& outer_ring = poly.outer();
                double offset_area = std::abs(bg::area(outer_ring));
                if (offset_area > area_threshold) {
                    resultRings.push_back(outer_ring);
                }
            }
        }
        return resultRings;
    }


    auto FillBridge::get_bbox(const Polygon_t& poly) {
        bg::model::box<Point_t> bbox;
        bg::envelope(poly, bbox);
        return bbox;
    }

    // 改进的安全合并函数
    Polygon_t FillBridge::safe_union(const Polygon_t& a, const Polygon_t& b) {
        const double tolerance = 1e-6;

        try {
            // 确保多边形有效
            Polygon_t valid_a = a;
            Polygon_t valid_b = b;

            if (!bg::is_valid(valid_a)) {
                bg::correct(valid_a);
            }

            if (!bg::is_valid(valid_b)) {
                bg::correct(valid_b);
            }

            // 首先检查两个多边形是否相交
            if (!polyIntersect(valid_a, valid_b)) {
                // 如果不相交，尝试合并它们
                MultiPolygon union_result;
                bg::union_(valid_a, valid_b, union_result);

                // 如果合并成功，返回结果
                if (union_result.size() == 1) {
                    return union_result[0];
                }

                // 如果合并产生多个多边形，尝试使用凸包
                if (union_result.size() > 1) {
                    Polygon_t hull;
                    bg::convex_hull(union_result, hull);

                    if (bg::is_valid(hull)) {
                        return hull;
                    }
                }

                // 如果所有方法都失败，返回第一个多边形
                return valid_a;
            }

            // 如果多边形相交，尝试合并它们
            MultiPolygon union_result;
            bg::union_(valid_a, valid_b, union_result);

            // 如果合并成功，返回结果
            if (union_result.size() == 1) {
                return union_result[0];
            }

            // 如果合并产生多个多边形，尝试使用交集和差异来合并
            if (union_result.size() > 1) {
                // 尝试找到最大的多边形
                double max_area = 0;
                Polygon_t largest_poly;

                for (const auto& poly : union_result) {
                    double area = bg::area(poly);
                    if (area > max_area) {
                        max_area = area;
                        largest_poly = poly;
                    }
                }

                return largest_poly;
            }

            // 如果所有方法都失败，返回第一个多边形
            return valid_a;
        }
        catch (const std::exception& e) {
            std::cerr << "Precise union failed: " << e.what() << std::endl;
            return a;
        }
    }

    // 改进的多边形合并函数
    void FillBridge::polygonsMerge(
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inputs,
        std::vector<std::pair<Polygon_t, IdIndex>>& outputs, size_t orientation) {

        // 使用R-tree加速多边形合并
        typedef std::pair<bg::model::box<Point_t>, size_t> BoxValue;
        bgi::rtree<BoxValue, bgi::quadratic<16>> rtree;

        for (size_t i = 0; i < inputs.size(); i++) {
            rtree.insert(std::make_pair(get_bbox(inputs[i].first), i));
        }

        std::vector<bool> merged(inputs.size(), false);
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> combined_inner_polys;

        // 尝试合并所有明显相交的多边形
        for (size_t i = 0; i < inputs.size(); i++) {
            if (merged[i]) continue;

            auto current = inputs[i];
            merged[i] = true;

            // 查找所有可能相交的多边形
            std::vector<BoxValue> candidates;
            auto bbox = get_bbox(current.first);
            rtree.query(bgi::intersects(bbox), std::back_inserter(candidates));

            // 尝试合并所有相交的多边形
            for (const auto& candidate : candidates) {
                size_t j = candidate.second;
                if (merged[j] || i == j) continue;

                if (polyIntersect(current.first, inputs[j].first)) {
                    Polygon_t union_result = safe_union(current.first, inputs[j].first);

                    // 检查合并后的多边形是否有效
                    if (bg::is_valid(union_result)) {
                        current.first = union_result;
                        current.second.insert(current.second.end(),
                            inputs[j].second.begin(),
                            inputs[j].second.end());
                        merged[j] = true;
                    }
                }
            }

            combined_inner_polys.push_back(current);
        }

        //检查合并后的多边形之间是否还有相交
        bool has_intersection = true;
        while (has_intersection) {
            has_intersection = false;
            std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> new_combined;
            std::vector<bool> processed(combined_inner_polys.size(), false);

            for (size_t i = 0; i < combined_inner_polys.size(); i++) {
                if (processed[i]) continue;

                auto current = combined_inner_polys[i];
                processed[i] = true;

                // 检查当前多边形是否与其他多边形相交
                for (size_t j = i + 1; j < combined_inner_polys.size(); j++) {
                    if (processed[j]) continue;

                    if (polyIntersect(current.first, combined_inner_polys[j].first)) {
                        Polygon_t union_result = safe_union(current.first, combined_inner_polys[j].first);

                        if (bg::is_valid(union_result)) {
                            current.first = union_result;
                            current.second.insert(current.second.end(),
                                combined_inner_polys[j].second.begin(),
                                combined_inner_polys[j].second.end());
                            processed[j] = true;
                            has_intersection = true;
                        }
                    }
                }

                new_combined.push_back(current);
            }

            combined_inner_polys = std::move(new_combined);
        }

        // 最终检查，确保没有多边形相交
        for (size_t i = 0; i < combined_inner_polys.size(); i++) {
            for (size_t j = i + 1; j < combined_inner_polys.size(); j++) {
                if (polyIntersect(combined_inner_polys[i].first, combined_inner_polys[j].first)) {
                    // 尝试最后一次合并
                    Polygon_t union_result = safe_union(combined_inner_polys[i].first, combined_inner_polys[j].first);

                    if (bg::is_valid(union_result)) {
                        combined_inner_polys[i].first = union_result;
                        combined_inner_polys[i].second.insert(combined_inner_polys[i].second.end(),
                            combined_inner_polys[j].second.begin(),
                            combined_inner_polys[j].second.end());
                        combined_inner_polys.erase(combined_inner_polys.begin() + j);
                        j--; // 调整索引
                    }
                }
            }
        }

        // 创建合并后的内多边形并记录映射
        for (const auto& combined : combined_inner_polys) {
            if (combined.second.size() > 1) {
                Ring outer_ring = bg::exterior_ring(combined.first);
                if (bg::area(outer_ring) > area_threshold) {
                    RingNode node = formatNode(++maxRid, 1, outer_ring, orientation);
                    addNode(node);

                    std::vector<IdIndex> result_ids{ node.id };
                    MergeMap2 mm2{ IdIndex{0, 0}, combined.second, result_ids };
                    mergeMap2.push_back(mm2);

                    outputs.emplace_back(combined.first, node.id);
                }
            }
            else {
                outputs.emplace_back(combined.first, combined.second[0]);
            }
        }
    }


    std::vector<RingNode> FillBridge::compute_complex_polygon_merge(
        const std::vector<RingNode>& outers
        , const std::vector<RingNode>& inners
    ) {

        // 准备内多边形
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inner_polys;
        for (const auto& inner : inners) {
            Polygon_t poly;
            poly.outer() = inner.ring;
            bg::correct(poly);
            inner_polys.emplace_back(poly, std::vector<IdIndex>{inner.id});
        }

        std::vector<std::pair<Polygon_t, IdIndex>> final_inners;
        polygonsMerge(inner_polys, final_inners, 1);

        //外多边形与内多边形合并
        std::vector<IdIndex> result_iis;
        for (const auto& outer : outers) {
            Polygon_t outer_poly;
            outer_poly.outer() = outer.ring;
            bg::correct(outer_poly);
            std::vector<Polygon_t> intersected_polys;
            std::vector<IdIndex> intersected_iis;

            // 与内多边形依次检查交集
            for (const auto& combined : final_inners) {
                if (polyIntersect(outer_poly, combined.first)) {
                    intersected_polys.emplace_back(combined.first);
                    intersected_iis.emplace_back(combined.second);
                }
                else {
                    result_iis.emplace_back(combined.second);
                }
            }
            //去掉相交的内多边形
            final_inners.erase(
                std::remove_if(
                    final_inners.begin(),
                    final_inners.end(),
                    [&intersected_iis](const auto& pair) {
                        return std::find(
                            intersected_iis.begin(),
                            intersected_iis.end(),
                            pair.second
                        ) != intersected_iis.end();
                    }
                ),
                final_inners.end()
            );

            if (intersected_polys.size() > 0) {  // 有相交的内多边形
                MultiPolygon result;
                MultiPolygon intersected_mp;
                for (const auto& poly : intersected_polys) {
                    intersected_mp.push_back(poly);
                }
                bg::difference(outer_poly, intersected_mp, result);
                for (const auto& poly : result) {
                    Ring _ring = bg::exterior_ring(poly);
                    if (bg::area(_ring) > area_threshold) {
                        RingNode node = formatNode(++maxRid, 1, _ring, -1);
                        addNode(node);
                        IdIndex _ii(node.id.id, 1);
                        std::vector<IdIndex> id_indices;
                        id_indices.emplace_back(_ii);
                        MergeMap2 mm2(outer.id, intersected_iis, id_indices);
                        mergeMap2.emplace_back(mm2);
                        result_iis.emplace_back(_ii);
                    }
                }
            }
            else {  // 与外多边形无相交
                result_iis.emplace_back(outer.id);
            }
        }

        //外多边形合并
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> outer_polys;
        for (IdIndex ii : result_iis) {
            RingNode rn = findNode(ii);
            Polygon_t poly;
            poly.outer() = rn.ring;
            bg::correct(poly);
            outer_polys.emplace_back(poly, std::vector<IdIndex>{ii});
        }
        std::vector<std::pair<Polygon_t, IdIndex>> final_polys;
        polygonsMerge(outer_polys, final_polys, -1);

        std::vector<RingNode> result;
        for (auto& pair : final_polys) {
            result.emplace_back(findNode(pair.second));
        }
        return result;
    }


    //生成环集
    void FillBridge::generateRings() {
        std::vector<RingNode> nodeHandles;

        for (auto& poly : o_polygons) {
            Ring ring = bg::exterior_ring(poly);
            bg::correct(ring);
            RingNode node = formatNode(maxRid++, 1, ring, -1);
            nodeHandles.emplace_back(node);
            addNode(node);
        }
        for (auto& inner : b_polygon.inners()) {
            bg::correct(inner);
            RingNode node = formatNode(maxRid++, 1, inner, 1);
            nodeHandles.emplace_back(node);
            addNode(node);
        }
        double t = 1;
        //遍历当前用于判断互交的节点集
        while (!nodeHandles.empty()) {
            std::vector<RingNode> outers;
            std::vector<RingNode> inners;
            for (const auto& rn : nodeHandles) {
                if (rn.orientation == -1) {
                    outers.emplace_back(rn);
                }
                else {
                    inners.emplace_back(rn);
                }
            }
            std::vector<RingNode> nodes = compute_complex_polygon_merge(
                outers, inners);
            std::vector<RingNode> offsetNodes;
            //偏移
            for (auto node : nodes) {
                std::vector<Ring> ring0 = offsetRing(
                    node.ring,
                    static_cast<double>(node.orientation) * _offset * t, area_threshold
                );

                if (ring0.size() == 1) {
                    RingNode node0 = formatNode(node.id.id, ++node.id.index, ring0[0], node.orientation);
                    if (bg::area(node0.ring) > area_threshold) {
                        addNode(node0);
                        offsetNodes.emplace_back(node0);
                    }
                }
                else if (ring0.size() > 1) {
                    std::vector<IdIndex> id_indices;
                    for (auto r : ring0) {
                        RingNode node0 = formatNode(++maxRid, 1, r, node.orientation);
                        if (bg::area(node0.ring) > area_threshold) {
                            addNode(node0);
                            offsetNodes.emplace_back(node0);
                            IdIndex ii(maxRid, 1);
                            id_indices.emplace_back(ii);
                        }
                    }
                    //分裂映射
                    offsetMap.insert(std::pair<IdIndex, std::vector<IdIndex>>({ node.id.id,node.id.index }, id_indices));
                }
            }
            t = 1;
            //清空
            nodeHandles.clear();
            if (!offsetNodes.empty()) {
                if (offsetNodes.size() == 2) {
                    Polygon_t poly1, poly2;
                    bg::append(poly1.outer(), offsetNodes[0].ring);
                    bg::append(poly2.outer(), offsetNodes[1].ring);
                    std::vector<IdIndex> id_indices;
                    if (offsetNodes[0].orientation == 1 && isPolygonContained(poly1, poly2)) {
                        MergeMap mm(offsetNodes[0].id, offsetNodes[1].id, id_indices);
                        containMap.emplace_back(mm);
                        break;
                    }
                    if (offsetNodes[1].orientation == 1 && isPolygonContained(poly2, poly1)) {
                        MergeMap mm(offsetNodes[1].id, offsetNodes[0].id, id_indices);
                        containMap.emplace_back(mm);
                        break;
                    }
                }
                //重填
                nodeHandles.insert(nodeHandles.end(), offsetNodes.begin(), offsetNodes.end());
            }
        }
    }

    //形成树
    void FillBridge::formatTree() {
        //遍历不同环类型的节点集合
        for (auto& ringNode : ringNodes) {

            // 使用反向迭代器进行反序遍历
            std::vector<RingNode>& vec = ringNode.second;
            const auto it0 = ringNode.second.begin();
            if (it0->orientation == 1) { // 向外  倒序
                for (size_t i = ringNode.second.size() - 1; i > 0; i--) {
                    ringNode.second[i].children.emplace_back(ringNode.second[i - 1]);
                    ringNode.second[i - 1].parent = &ringNode.second[i];
                }
            }
            else { //向内 正序
                for (int i = 0; i < vec.size() - 1; i++) {
                    ringNode.second[i].children.emplace_back(ringNode.second[i + 1]);
                    ringNode.second[i + 1].parent = &ringNode.second[i];
                }
            }
        }
        //遍历分裂
        for (auto& map : offsetMap) { //偏移产生的分裂
            RingNode& it1 = findNode(map.first);
            for (const auto& ii : map.second) {
                RingNode& it2 = findNode(ii);
                if (it1.orientation == -1) {
                    it1.children.emplace_back(it2);
                    it2.parent = &it1;
                }
                else {
                    it1.parent = &it2;
                    it2.children.emplace_back(it1);
                }
            }
        }
        for (auto& map : containMap) {
            RingNode& it1 = findNode(map.ii1); //内
            RingNode& it2 = findNode(map.ii2);  //外
            it2.children.insert(it2.children.end(), it1.children.begin(), it1.children.end());
            it1.isHide = true;
            it1.parent = nullptr;
        }
        for (auto& map : mergeMap2) {
            if (map.outer_ii.id == 0 && map.outer_ii.index == 0) {  //内多边形合并
                for (const auto& merged : map.merged_iis) {
                    RingNode& merge_rn = findNode(merged);
                    for (const auto& inner_ii : map.inner_iis) {
                        RingNode& inner_rn = findNode(inner_ii);
                        merge_rn.children.insert(merge_rn.children.end(),
                            inner_rn.children.begin(), inner_rn.children.end());
                        inner_rn.isHide = true;
                    }
                }
            }
            else {  //外多边形  合并
                RingNode& outer = findNode(map.outer_ii);
                for (const auto& merged : map.merged_iis) {
                    RingNode& merge_rn = findNode(merged);

                    outer.parent->children.emplace_back(merge_rn);
                    merge_rn.parent = outer.parent;
                    outer.isHide = true;
                    IdIndex _ii = outer.id;
                    //将outer从父节点的子节点集中移除
                    outer.parent->children.erase(std::remove_if(outer.parent->children.begin(), outer.parent->children.end(),
                        [_ii](RingNode x) { return x.id.id == _ii.id && x.id.index == _ii.index; }),
                        outer.parent->children.end());
                    //将内多边形的子节点作为merge_rn的子节点
                    for (const auto& inner_ii : map.inner_iis) {
                        RingNode& inner_rn = findNode(inner_ii);
                        merge_rn.children.insert(merge_rn.children.end(),
                            inner_rn.children.begin(), inner_rn.children.end());
                        inner_rn.isHide = true;
                    }
                }

            }
        }

    }
    //树遍历
    void FillBridge::dfs(
        RingNode& node,
        std::vector<IdIndex>& visited
    ) {
        for (auto& one : visited) {
            if (one == node.id) {
                return;
            }
        }
        visited.emplace_back(node.id);
        all_rings.push_back(node.ring);
        for (auto& child : node.children) {
            dfs(findNode(child.id), visited);
            handleBridge(node.id, child.id);
        }
    }

    // 计算点到线段的最近点
    Point_t FillBridge::closest_point_on_segment(const Point_t& p, const Point_t& seg_start, const Point_t& seg_end) {
        const double x = bg::get<0>(p);
        const double y = bg::get<1>(p);
        const double x1 = bg::get<0>(seg_start);
        const double y1 = bg::get<1>(seg_start);
        const double x2 = bg::get<0>(seg_end);
        const double y2 = bg::get<1>(seg_end);

        // 线段向量
        const double dx = x2 - x1;
        const double dy = y2 - y1;

        // 如果线段长度为零，返回起点
        if (dx == 0 && dy == 0) {
            return seg_start;
        }

        // 计算投影参数t
        double t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy);

        // 限制 t 在 [0,1] 范围内，确保投影点在线段上
        t = std::max(0.0, std::min(1.0, t));

        // 计算投影点坐标
        return Point_t(bg::get<0>(seg_start) + t * dx,
            bg::get<1>(seg_start) + t * dy);
    }

    // 在环的边上查找离给定点最近的点
    Point_t FillBridge::find_closest_point_on_ring_edges(Ring& ring, const Point_t& p0, size_t& e_index1) {
        Point_t closest = ring[0]; // 初始化为Ring的第一个点
        double minDist = bg::distance(p0, closest);

        // 遍历Ring的每条边
        for (size_t i = 0; i < ring.size(); ++i) {
            // 获取当前边的两个端点
            Point_t a = ring[i];
            Point_t b = ring[(i + 1) % ring.size()]; // 处理闭合边

            // 计算点到边的最近点
            Point_t projection = closest_point_on_segment(p0, a, b);

            // 计算距离
            double dist = bg::distance(p0, projection);

            // 更新最近点
            if (dist < minDist) {
                minDist = dist;
                closest = projection;
                e_index1 = i;
            }
        }
        return closest;
    }

    Point_t FillBridge::find_point_at_distance_clockwise(Ring& ring, const Point_t& start_point,
        size_t _index, double d, size_t& e_index)
    {
        // 确保环是闭合的（首尾点相同）
        Ring normalized_r1 = ring;
        double remaining_distance = d;
        size_t current_edge = _index;
        Point_t current_point = start_point;
        while (remaining_distance > 0) {
            // 获取当前边的终点
            Point_t next_point = normalized_r1[(current_edge + 1) % normalized_r1.size()];

            // 计算当前边剩余长度
            double edge_length = bg::distance(current_point, next_point);

            if (edge_length >= remaining_distance) {
                // 目标点在当前边上
                double ratio = remaining_distance / edge_length;
                double x = current_point.x() + ratio * (next_point.x() - current_point.x());
                double y = current_point.y() + ratio * (next_point.y() - current_point.y());
                e_index = current_edge;
                return Point_t(x, y);
            }
            else {
                // 目标点在下一条边上，更新剩余距离
                remaining_distance -= edge_length;
                current_point = next_point;
                current_edge = (current_edge + 1) % normalized_r1.size();
            }
        }
        e_index = _index;
        return start_point;
    }


    //点在环上的位置索引
    int FillBridge::findPointIndex(const Ring& ring, const Point_t& p0) {
        for (size_t i = 0; i < ring.size(); ++i) {
            if (equal(ring[i], p0)) {
                return static_cast<int>(i);
            }
        }
        return -1; // 未找到
    }


    void FillBridge::handleBridge(IdIndex o_ii, IdIndex i_ii) {
        // 检查是否已存在相同桥接，存在则直接返回
        for (const auto& b : bridges) {
            if ((b.to_ii == i_ii && b.from_ii == o_ii) || (b.to_ii == o_ii && b.from_ii == i_ii)) {
                return;
            }
        }

        RingNode& inner = findNode(i_ii);
        RingNode& outer = findNode(o_ii);


        // 计算内环中点及周长，同时计算内环到外环的平均距离
        double perimeter = 0;
        std::vector<MidPoints> inner_mid_points = find_mid_points_on_ring(inner.ring, perimeter);
        if (inner_mid_points.empty() || perimeter <= 0) { // 异常处理：内环无有效点
            return;
        }

        // 计算内环中点到外环的平均距离（环间平均距离）
        double total_ring_dist = 0.0;
        int valid_points = 0;
        for (const auto& mid : inner_mid_points) {
            size_t dummy_idx;
            Point_t closest_outer = find_closest_point_on_ring_edges(outer.ring, mid.mid_point, dummy_idx);
            double dist = bg::distance(mid.mid_point, closest_outer);
            if (dist > 1e-6) {
                total_ring_dist += dist;
                valid_points++;
            }
        }
        if (valid_points == 0) {
            std::cout << "Warning: No valid distance between inner and outer rings." << std::endl;
            return;
        }
        double _avg_ring_dist = total_ring_dist / valid_points; // 环间平均距离（核心基准）
        if (_avg_ring_dist > 1.5 * _offset) {
            _avg_ring_dist = 1.5 * _offset;
        }

        // 初始化参数
        Point_t p0, p1, p2, p3;
        size_t p0_idx = 0, p1_idx = 0, p2_idx = 0, p3_idx = 0;
        double total_len = 0;
        int start_idx = 0;
        const int max_attempts = 300;
        int attempt_count = 0;
        bool find_bridge = false;

        // 初始化p0（首次或基于最后一个桥接）
        if (bridges.empty()) {
            p0 = inner_mid_points[0].mid_point;
            p0_idx = inner_mid_points[0].index;
        }
        else {
            const BridgeMap& last_bridge = bridges.back();
            size_t last_idx = findPointIndex(inner.ring, last_bridge.from);
            find_further_point(inner_mid_points, last_bridge.from, last_idx, p0, p0_idx);
        }
        bool b = false;
        // 辅助函数：调整p0位置和基准距离
        auto adjust_params = [&]() {
            if (total_len > perimeter && start_idx < inner_mid_points.size()) {
                _avg_ring_dist *= 0.9;
                total_len = 0;
                start_idx = (start_idx + 1) % inner_mid_points.size();
                p0 = inner_mid_points[start_idx].mid_point;
                p0_idx = inner_mid_points[start_idx].index;
            }
            else {
                _avg_ring_dist = total_ring_dist / valid_points;
                if (_avg_ring_dist > 1.5 * _offset) {
                    _avg_ring_dist = 1.5 * _offset;
                }
                p0 = p1;
                p0_idx = p1_idx;
            }
            };

        // 核心桥接逻辑（基于环间平均距离）
        while (!find_bridge && attempt_count < max_attempts) {
            attempt_count++;

            // 找外环上与p0、p1最近的点p2、p3
            p2 = find_closest_point_on_ring_edges(outer.ring, p0, p2_idx);
            // 在内环上按环间平均距离找p1（顺时针方向）
            p1 = find_point_at_distance_clockwise(inner.ring, p0, p0_idx, _avg_ring_dist, p1_idx);
            //p3 = find_point_at_distance_clockwise(outer.ring, p2, p2_idx, _avg_ring_dist, p3_idx);
            p3 = find_closest_point_on_ring_edges(outer.ring, p1, p3_idx);
            //p1 = find_closest_point_on_ring_edges(inner.ring, p3, p1_idx);

            total_len += _avg_ring_dist;

            // 检查点对距离是否过近（基于环间平均距离的95%）
            const double min_dist = 0.95 * _avg_ring_dist;
            if (bg::distance(p0, p1) < min_dist || bg::distance(p2, p3) < min_dist) {
                adjust_params();
                continue;
            }

            // 检查桥接线段是否相交（seg1: p2-p1 与 seg2: p3-p0）
            Segment seg1(p2, p1), seg2(p3, p0);
            std::vector<Point_t> intersections;
            bg::intersection(seg1, seg2, intersections);
            if (intersections.empty() || intersections.size() > 1) { // 无交点或多交点（无效）
                adjust_params();
                continue;
            }

            // 检查桥接长度是否过长（不超过环间平均距离的1.5倍）
            const double max_len = 1.5 * _avg_ring_dist;
            if (bg::distance(p1, p2) > max_len || bg::distance(p0, p3) > max_len) {
                adjust_params();
                continue;
            }

            // 检查与已有桥接是否过近（距离小于环间平均距离）
            bool too_close = false;
            for (const auto& existing : bridges) {
                if (bg::distance(p0, existing.from) < _avg_ring_dist ||
                    bg::distance(p1, existing.from2) < _avg_ring_dist ||
                    bg::distance(p0, existing.from2) < _avg_ring_dist ||
                    bg::distance(p1, existing.from) < _avg_ring_dist) {
                    too_close = true;
                    break;
                }
            }
            if (too_close) {
                adjust_params();
                continue;
            }

            // 插入点并添加桥接
            insertPointIntoRing(inner.ring, p0);
            insertPointIntoRing(inner.ring, p1);
            insertPointIntoRing(outer.ring, p2);
            insertPointIntoRing(outer.ring, p3);

            bridges.emplace_back(p2, p1, p3, p0, o_ii, i_ii, false);
            find_bridge = true;
        }
    }


    int FillBridge::findIndex(std::vector<Point_t> points, const Point_t& p0) {
        for (int i = 0; i < points.size(); ++i) {
            if (equal(points[i], p0)) {
                return static_cast<int>(i);
            }
        }
        return -1; // 未找到
    }

    void FillBridge::traverseRing(
        RingNode& node,
        IdIndex parent,
        Point_t& start,
        Point_t& end,
        bool isOutermostLayer
    ) {
        int size = node.ring.size();
        int s_index = findPointIndex(node.ring, start);
        int e_index = findPointIndex(node.ring, end);
        if (s_index == -1 || e_index == -1) {
            return;
        }
        std::vector<Point_t> c_points;
        std::vector<Point_t> cc_points;
        std::vector<Point_t> points;
        for (int i = s_index; !equal(node.ring[i], end); i = (i + 1) % size) {
            c_points.emplace_back(node.ring[i]);
        }
        c_points.emplace_back(end);
        for (int i = s_index; !equal(node.ring[i], end); i = (i - 1 + size) % size) {
            cc_points.emplace_back(node.ring[i]);
        }
        cc_points.emplace_back(end);

        if (cc_points.size() > c_points.size()) {
            points.assign(cc_points.begin(), cc_points.end());
        }
        else
        {
            points.assign(c_points.begin(), c_points.end());
        }
        int index = 0;
        do {
            if (findIndex(path, points[index]) == -1) {
                path.emplace_back(points[index]);
                BridgeMap* bridge = nullptr;
                //判断该点是否为桥接点
                for (auto& bm : bridges) {
                    if (equal(bm.from, points[index]) || equal(bm.from2, points[index])) {
                        bridge = &bm;
                        break;
                    }
                }
                if (bridge != nullptr) {
                    RingNode& rn = findNode(bridge->to_ii); //to 环
                    if (equal(bridge->from, points[index])) {
                        traverseRing(rn, node.id, bridge->to, bridge->to2, false);
                        index = findIndex(points, bridge->from2);
                    }
                    else if (equal(bridge->from2, points[index])) {
                        traverseRing(rn, node.id, bridge->to2, bridge->to, false);
                        index = findIndex(points, bridge->from);
                    }
                    if (index == -1) {
                        std::cout << bridge->from.x() << "," << bridge->from.y() << std::endl;
                        break;
                    }
                    path.emplace_back(points[index]);
                }
            }
            index += isOutermostLayer ? -1 : 1;
        } while (index > 0 && index < points.size());
    }


} // namespace Slic3r


