#include "LocalNet.h"

void LocalNet::initGridBoxes() {
    // rangeSearch & slice routeGuides
    vector<vector<db::GridBoxOnLayer>> guides(database.getLayerNum());
    for (auto& guide : routeGuides) {
        guides[guide.layerIdx].push_back(database.rangeSearch(guide));
    }
    routeGuides.clear();
    gridRouteGuides.clear();
    for (unsigned layerIdx = 0; layerIdx != database.getLayerNum(); ++layerIdx) {
        db::GridBoxOnLayer::sliceGridPolygons(guides[layerIdx]);
        for (const auto& guide : guides[layerIdx]) {
            gridRouteGuides.push_back(guide);
            routeGuides.push_back(database.getLoc(guide));
        }
    }

    // init gridPinAccessBoxes
    database.getGridPinAccessBoxes(dbNet, gridPinAccessBoxes);
    replaceMacroPabs(gridPinAccessBoxes);
    for (unsigned pinIdx = 0; pinIdx != numOfPins(); ++pinIdx) {
        pinAccessBoxes[pinIdx].clear(); 
        for (auto& box : gridPinAccessBoxes[pinIdx]) {
            pinAccessBoxes[pinIdx].push_back(database.getLoc(box));
        }
    }

    getRouteGuideMapping();
}

void LocalNet::getRouteGuideMapping() {
    dbRouteGuideIdxes.resize(gridRouteGuides.size());
    for (unsigned i = 0; i < dbRouteGuideIdxes.size(); i++) {
        const db::BoxOnLayer& routeGuide = routeGuides[i];
        boostBox query_box(boostPoint(routeGuide.x.low, routeGuide.y.low), boostPoint(routeGuide.x.high, routeGuide.y.high));
        std::vector<std::pair<boostBox, int>> results;
        dbNet.routeGuideRTrees[routeGuide.layerIdx].query(bgi::intersects(query_box), std::back_inserter(results));

        dbRouteGuideIdxes.reserve(results.size());

        for (const auto& v : results) {
            dbRouteGuideIdxes[i].push_back(v.second);
        }
    }
}

void LocalNet::initNumOfVertices() {
    estimatedNumOfVertices = 0;
    for (unsigned b1 = 0; b1 < gridRouteGuides.size(); b1++) {
        for (unsigned b2 : guideConn[b1]) {
            if (b2 <= b1) {
                continue;
            }

            const db::GridBoxOnLayer& box1 = gridRouteGuides[b1];
            const db::GridBoxOnLayer& box2 = gridRouteGuides[b2];

            if (!database.isValid(box1) || !database.isValid(box2)) {
                continue;
            }

            const int upperIdx = box1.layerIdx > box2.layerIdx ? b1 : b2;
            const int lowerIdx = box1.layerIdx < box2.layerIdx ? b1 : b2;

            const db::ViaBox viaBox = database.getViaBoxBetween(gridRouteGuides[lowerIdx], gridRouteGuides[upperIdx]);
            if (!database.isValid(viaBox)) {
                continue;
            }
            estimatedNumOfVertices += ((viaBox.lower.trackRange.high + 1 - viaBox.lower.trackRange.low) *
                                       (viaBox.upper.trackRange.high + 1 - viaBox.upper.trackRange.low));
        }
    }
}

bool LocalNet::checkPin() const {
    for (const auto& gridAccessBoxes : gridPinAccessBoxes) {
        int isPinValid = false;
        for (const auto& gridBox : gridAccessBoxes) {
            if (database.isValid(gridBox)) {
                isPinValid = true;
                break;
            }
        }
        if (!isPinValid) {
            return false;
        }
    }
    return true;
}

void LocalNet::print() const {
    db::NetBase::print();

    log() << "dbRouteGuideIdxes" << std::endl;
    for (unsigned i = 0; i < dbRouteGuideIdxes.size(); i++) {
        log() << "guide " << i << ": " << dbRouteGuideIdxes[i] << std::endl;
    }

    RouteGuideGraph::print();
}

int LocalNet::getViaPenalty(int guideIdx1, int trackIdx1, int cpIdx1, int guideIdx2, int trackIdx2, int cpIdx2) const {
    bool isContain = false;
    for (auto idx1 : dbRouteGuideIdxes[guideIdx1]) {
        const db::GridBoxOnLayer& box1 = dbNet.gridRouteGuides[idx1];
        if (box1.trackRange.Contain(trackIdx1) && box1.crossPointRange.Contain(cpIdx1)) {
            isContain = true;
            break;
        }
    }
    if (!isContain) return 1;

    for (auto idx2 : dbRouteGuideIdxes[guideIdx2]) {
        const db::GridBoxOnLayer& box2 = dbNet.gridRouteGuides[idx2];
        if (box2.trackRange.Contain(trackIdx2) && box2.crossPointRange.Contain(cpIdx2)) return 0;
    }
    return 1;
}

double LocalNet::getWireSegmentPenalty(int guideIdx, int trackIdx, int cpIdx1, int cpIdx2) const {
    vector<utils::IntervalT<DBU>> intersections;

    db::GridBoxOnLayer segmentGridBox(gridRouteGuides[guideIdx].layerIdx, {trackIdx, trackIdx}, {cpIdx1, cpIdx2});
    db::BoxOnLayer segmentBox = database.getLoc(segmentGridBox);

    utils::IntervalT<DBU> segment = segmentBox.x.range() == 0 ? segmentBox.y : segmentBox.x;

    for (auto idx : dbRouteGuideIdxes[guideIdx]) {
        const db::GridBoxOnLayer& box = dbNet.gridRouteGuides[idx];

        if (box.trackRange.Contain(trackIdx)) {
            if (box.crossPointRange.Contain(cpIdx1) && box.crossPointRange.Contain(cpIdx2)) {
                return 0;
            } else {
                const utils::IntervalT<DBU>& boxSegment =
                    segmentBox.y.range() == 0 ? dbNet.routeGuides[idx].x : dbNet.routeGuides[idx].y;
                utils::IntervalT<DBU> intersection = boxSegment.IntersectWith(segment);
                if (intersection.IsValid()) {
                    intersections.push_back(intersection);
                }
            }
        }
    }

    if (intersections.empty()) {
        return 1;
    } else {
        sort(intersections.begin(),
             intersections.end(),
             [](const utils::IntervalT<DBU>& a, const utils::IntervalT<DBU>& b) {
                 return a.low < b.low || (a.low == b.low && a.high < b.high);
             });

        DBU length = 0, low = 0, high = 0;
        int sz = intersections.size();
        for (int i = 0; i < sz; ++i) {
            low = intersections[i].low;
            high = intersections[i].high;
            while (i + 1 < sz && high >= intersections[i + 1].low) {
                ++i;
                high = intersections[i].high > high ? intersections[i].high : high;
            }
            length += (high - low);
        }

        return 1 - length * 1.0 / segment.range();
    }
}

int LocalNet::getCrossPointPenalty(int guideIdx, int trackIdx, int cpIdx) const {
    for (auto idx : dbRouteGuideIdxes[guideIdx]) {
        const db::GridBoxOnLayer& box = dbNet.gridRouteGuides[idx];
        if (box.trackRange.Contain(trackIdx) && box.crossPointRange.Contain(cpIdx)) return 0;
    }
    return 1;
}

// PARITAL RIPUP
void LocalNet::traverse(std::shared_ptr<db::GridSteiner> node) {
    bool isPin = node->pinIdx >= 0;
    node->distance = -1;
    for (auto c : node->children) {
        if (c->isVio)
            traverse(c);
        else
            isPin = true;
    }
    if (isPin && node != pnetPins[0])
        pnetPins.push_back(node);
}

void LocalNet::initPNetPins() {
    gridTopo.clear();
    pnetPins.push_back(dbNet.pnets[pnetIdx]);
    traverse(pnetPins[0]);

    // do not reroute small pnets
    int dist = 0;
    for (size_t i=1; i<pnetPins.size(); i++) {
        dist += (abs(pnetPins[0]->trackIdx - pnetPins[i]->trackIdx) +
                abs(pnetPins[0]->crossPointIdx - pnetPins[i]->crossPointIdx));
    }
    if (dist < 5) {
        // pnetPins[0]->isVio = false;
        pnetPins.clear();
        return;
    }

    pinAccessBoxes.clear();
    pinAccessBoxes.resize(pnetPins.size());
    for (size_t i=0; i<pnetPins.size(); i++) {
        auto p = pnetPins[i].get();
        bool notfix = true;
        if (p->pinIdx < 0 || (p->parent && !p->parent->isVio)) notfix = false;
        else {
            for (auto c : p->children) {
                if (!c->isVio)
                    notfix = false;
            }
        }
        if (notfix)
            pinAccessBoxes[i] = dbNet.pinAccessBoxes[p->pinIdx];
        pinAccessBoxes[i].emplace_back(p->layerIdx, database.getLoc(*(p)));
    }
}

inline void printndoes(std::shared_ptr<db::GridSteiner> node) {
    // const auto &np = database.getLoc(*node);
    // std::cout << "  " << node->layerIdx << ", " << np << ";  pinIdx : " << node->pinIdx
    //             << "; vio : " << int(node->isVio) << std::endl;
    std::cout << *node << "; " << node->distance << std::endl;
    for (auto c : node->children) {
        if (c->isVio)
            printndoes(c);
    }
}

void LocalNet::creatLocalRouteGuides() {
    db::GridBoxOnLayer gbox;
    vector<std::shared_ptr<db::GridSteiner>> current, other;
    other.push_back(pnetPins[0]);
    while (!other.empty()) {
        auto node = other.back();
        other.pop_back();
        gbox.layerIdx = node->layerIdx;
        gbox.trackRange.Set(node->trackIdx);
        gbox.crossPointRange.Set(node->crossPointIdx);
        current.push_back(node);
        while (!current.empty()) {
            node = current.back();
            current.pop_back();
            for (auto c : node->children) {
                if (!c->isVio) continue;
                if (c->layerIdx == node->layerIdx) {
                    gbox.trackRange.FastUpdate(c->trackIdx);
                    gbox.crossPointRange.FastUpdate(c->crossPointIdx);
                    current.push_back(c);
                }
                else other.push_back(c);
            }
        }
        gbox.trackRange.low -= 7;
        if (gbox.trackRange.low < 0) gbox.trackRange.low = 0;
        gbox.trackRange.high = min(database.getTrackEnd(gbox.layerIdx),
                                    gbox.trackRange.high + 7);
        routeGuides.push_back(database.getLoc(gbox));
    }
    // if (routeGuides.size() == 1 && routeGuides[0].layerIdx < (database.getLayerNum() - 1)) {
    //     routeGuides.emplace_back(routeGuides[0].layerIdx + 1, routeGuides[0]);
    //     db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_1, 1);
    // }
}

inline bool inGuide(db::GridSteiner *node, const vector<db::BoxOnLayer> &guides) {
    const auto &nloc = database.getLoc(*node);
    for (size_t i=0; i<guides.size(); i++) {
        if (guides[i].layerIdx == node->layerIdx && guides[i].Contain(nloc)) {
            node->distance = int(i);
            return true;
        }
    }
    return false;
}

inline void selectGuides(db::GridSteiner *node,
                            int preOut,
                            vector<int> &selected,
                            vector<std::pair<int, int>> &breaks,
                            const vector<db::BoxOnLayer> &guides) {
    if (inGuide(node, guides) || node->distance >= 0) {
        selected[node->distance]++;
        // if (node->parent && node->parent->distance < 0)
        if (preOut > 0)
            breaks.emplace_back(preOut, node->distance);
    }
    // else if (node->parent && node->parent->distance >= 0)
    //     breaks.emplace_back(node->parent->distance, -1);
    else if (preOut == -1) preOut = node->parent->distance;
    
    for (auto c : node->children) {
        if (c->isVio) {
            selectGuides(c.get(), preOut, selected, breaks, guides);
        }
    }
}

void LocalNet::creatAdaptiveRouteGuides() {
    // const auto &mguides = dbNet.mergedGuides;
    // size_t msize = mguides.size();
    size_t msize = routeGuides.size();
    for (size_t i=0; i<msize; i++) {
        if (dbNet.routeGuideVios[i] >= 4)
            database.expandBox(routeGuides[i], 8);
        else
            database.expandBox(routeGuides[i], 5);
    }
    vector<int> selected(msize, 0), svios;
    vector<std::pair<int, int>> breaks;
    vector<db::BoxOnLayer> pickedGuides;

    // select route guides
    for (size_t i=0; i<pnetPins.size(); i++) {
        if (pnetPins[i]->pinIdx >=0 && !inGuide(pnetPins[i].get(), routeGuides)) {
            printf(" WARNING : net %d pin %d out of guide\n", idx, pnetPins[i]->pinIdx);
            for (const auto &pab : pinAccessBoxes[i]) {
                for (size_t j=0; j<msize; j++) {
                    auto &g = routeGuides[j];
                    if (g.layerIdx == pab.layerIdx && g.HasIntersectWith(pab)) {
                        g.Set(g.layerIdx, g.UnionWith(pab));
                        pnetPins[i]->distance = int(j);
                        break;
                    }
                }
            }
            if (pnetPins[i]->distance < 0) {
                printf(" ERROR : net %d pin %d totally out of guide\n", idx, pnetPins[i]->pinIdx);
                std::cout << pinAccessBoxes[i] << std::endl;
                std::cout << dbNet.pinAccessBoxes[pnetPins[i]->pinIdx] << std::endl;
                std::cout << routeGuides << std::endl;
            }
        }
    }
    selectGuides(pnetPins[0].get(), -1, selected, breaks, routeGuides);
    for (size_t i=0; i<pnetPins.size(); i++) {
        if (pnetPins[i]->distance < 0)
            printf(" ERROR : net %d pnet %d pin %d (%d) out of guide\n", idx, pnetIdx, pnetPins[i]->pinIdx, int(i));
    }
    
    // ensure connectivity
    if (!breaks.empty()) {
        vector<int> curr, next;
        std::unordered_set<int> visitedGuides;
        std::unordered_map<int, int> parent;
        for (const auto &pair : breaks) {
            visitedGuides.insert(pair.first);
            curr.push_back(pair.first);
            while (!curr.empty()) {
                for (int id : curr) {
                    for (int i=0; i<msize; i++) {
                        if (visitedGuides.count(i)) continue;
                        if (abs(routeGuides[id].layerIdx-routeGuides[i].layerIdx) < 2 &&
                            routeGuides[id].HasIntersectWith(routeGuides[i])) {
                            next.push_back(i);
                            visitedGuides.insert(i);
                            parent.insert({i, id});
                        }
                    }
                }
                if (visitedGuides.count(pair.second))
                    next.clear();
                curr = move(next);
            }
            if (visitedGuides.count(pair.second)) {
                int currId = pair.second;
                while (parent.count(currId)) {
                    currId = parent[currId];
                    selected[currId]++;
                }
            }
            else {
                printf(" ERROR : net %d pnet %d open : %d, %d\n", idx, pnetIdx, pair.first, pair.second);
                printndoes(pnetPins[0]);
                for (size_t i=0; i<msize; i++) {
                    std::cout << i << " : " << routeGuides[i] << std::endl;
                }
            }
            visitedGuides.clear();
            parent.clear();
        }
    }

    for (size_t i=0; i<msize; i++) {
        if (selected[i] > 0) {
            // pickedGuides.push_back(mguides[i]);
            pickedGuides.push_back(routeGuides[i]);
            // svios.push_back(dbNet.mergedVios[i]);
            svios.push_back(dbNet.routeGuideVios[i]);
        }
    }

    // slice route guide
    for (auto p : pnetPins) {
        db::GridSteiner *nei = nullptr;
        if (p->parent && !(p->parent->isVio)) nei = p->parent.get();
        else {
            for (auto c : p->children) {
                if (!c->isVio) {
                    nei = c.get();
                    break;
                }
            }
        }
        if (nei == nullptr || nei->layerIdx != p->layerIdx) continue;
        const auto &pp = database.getLoc(*p);
        const auto &np = database.getLoc(*nei);
        auto &box = routeGuides[p->distance];
        if ((box.x.range() > 20000 || box.y.range() > 20000) && box.Contain(np)) {
            int dir = 1 - database.getLayerDir(box.layerIdx);
            if (np[dir] < pp[dir]) box[dir].low = pp[dir];
            else if (np[dir] > pp[dir]) box[dir].high = pp[dir];
        }
    }
    routeGuides = move(pickedGuides);
    // addDiffLayerRouteGuides();
    int lnum = database.getLayerNum() - 1;
    size_t rsize = svios.size();
    for (size_t i=0; i<rsize; i++) {
        int gl = routeGuides[i].layerIdx;
        if (svios[i] > db::setting.diffLayerGuideVioThres) {
            auto g = routeGuides[i];
            if (gl > 2) routeGuides.emplace_back(gl - 1, g);
            if (gl < lnum) routeGuides.emplace_back(gl + 1, g);
            db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_1, 1);
        }
        if (svios[i] > db::setting.diffLayerGuideVioThres * 2) {
            auto g = routeGuides[i];
            if (gl > 3) routeGuides.emplace_back(gl - 2, g);
            if (gl < lnum-1) routeGuides.emplace_back(gl + 2, g);
            db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_2, 1);
        }
    }
}

void LocalNet::getGridBoxes() {
    // pin access boxes
    vector<vector<db::GridBoxOnLayer>> gridBoxes;
    database.getGridPinAccessBoxes(dbNet, gridBoxes);
    replaceMacroPabs(gridBoxes);
    gridPinAccessBoxes.resize(pnetPins.size());
    for (size_t i=0; i<pnetPins.size(); i++) {
        // if (pnetPins[i]->children.empty() || pnetPins[i]->parent == nullptr) {
        //     if (pnetPins[i]->pinIdx < 0)
        //         printf("net %d pnetIdx %d pnet %ld is not pin but have no children\n", idx, pnetIdx, i);
        if (pinAccessBoxes[i][0].area() == 0) {
            gridPinAccessBoxes[i].push_back(database.rangeSearch(pinAccessBoxes[i][0]));
        }
        else {
            gridPinAccessBoxes[i] = move(gridBoxes[pnetPins[i]->pinIdx]);
            pinAccessBoxes[i].clear();
            for (const auto &gbox : gridPinAccessBoxes[i])
                pinAccessBoxes[i].push_back(database.getLoc(gbox));
        }
    }
    
    // route guides
    gridBoxes.clear();
    gridBoxes.resize(database.getLayerNum());
    for (const auto& guide : routeGuides) {
        const auto &gbox = database.rangeSearch(guide);
        if (database.isValid(gbox))
            gridBoxes[guide.layerIdx].push_back(gbox);
    }
    routeGuides.clear();
    gridRouteGuides.clear();
    for (unsigned layerIdx = 0; layerIdx != database.getLayerNum(); ++layerIdx) {
        db::GridBoxOnLayer::sliceGridPolygons(gridBoxes[layerIdx]);
        for (const auto& guide : gridBoxes[layerIdx]) {
            gridRouteGuides.push_back(guide);
            routeGuides.push_back(database.getLoc(guide));
        }
    }
    
    getRouteGuideMapping();
}

void LocalNet::addDiffLayerRouteGuides() {
    size_t gsize = routeGuides.size();
    unsigned int lnum = database.getLayerNum();
    for (size_t i=0; i<gsize; i++) {
        const auto &box = routeGuides[i];
        boostBox qbox(boostPoint(box.x.low, box.y.low),
                        boostPoint(box.x.high, box.y.high));
        vector<std::pair<boostBox, int>> results;
        dbNet.routeGuideRTrees[box.layerIdx].query(bgi::intersects(qbox), std::back_inserter(results));
        int count = 0;
        for (const auto &pair : results) {
            count += pair.second;
        }
        int bl = box.layerIdx;
        if (count > db::setting.diffLayerGuideVioThres) {
            if (bl > 2) routeGuides.emplace_back(bl - 1, box);
            if (bl < (lnum - 1)) routeGuides.emplace_back(bl + 1, box);
            db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_1, 1);
        }
    }
}

void LocalNet::replaceMacroPabs(vector<vector<db::GridBoxOnLayer>> &gboxes) {
    for (size_t i=0; i<dbNet.pinAccessBoxes.size(); i++) {
        if (!dbNet.rsynPins[i].isMacroBlockPin()) continue;
        gboxes[i].clear();
        for (auto pab : dbNet.pinAccessBoxes[i]) {
            DBU halfWidth = database.getLayer(pab.layerIdx).width / 2,
                pitch = database.getLayer(pab.layerIdx).pitch;
            for (const auto &obs : database.obsBoxes[pab.layerIdx]) {
                if (obs.HasIntersectWith(pab)) {
                    for (int d=0; d<2; d++) {
                        if (pab[d].low == obs[d].low)
                            pab[d].low -= pitch * 2;
                        else pab[d].low += halfWidth;

                        if (pab[d].high == obs[d].high)
                            pab[d].high += pitch * 2;
                        else pab[d].high -= halfWidth;
                    }
                    break;
                }
            }
            const auto &g = database.rangeSearch(pab);
            if (!db::rrrIterSetting.fullyRoute && !database.isValid(g)) {
                printf(" net %d pnet %d pin %d : ", idx, pnetIdx, int(i));
                std::cout << g << "; " << pab << std::endl;
            }
            gboxes[i].push_back(g);
        }
    }
}