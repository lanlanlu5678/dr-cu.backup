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
        dist += abs(pnetPins[0]->trackIdx - pnetPins[i]->trackIdx) +
                abs(pnetPins[0]->crossPointIdx - pnetPins[i]->crossPointIdx);
    }
    if (pnetPins.size() == 1 || dist < 10) {
        pnetPins[0]->isVio = false;
        return;
    }

    pinAccessBoxes.clear();
    pinAccessBoxes.resize(pnetPins.size());
    for (size_t i=0; i<pnetPins.size(); i++) {
        if (pnetPins[i]->pinIdx < 0 || pnetPins[i]->children.size() > 0) {
            pinAccessBoxes[i].emplace_back(pnetPins[i]->layerIdx, database.getLoc(*(pnetPins[i])));
        }
        else {
            pinAccessBoxes[i] = dbNet.pinAccessBoxes[pnetPins[i]->pinIdx];
        }
    }
}

inline void printndoes(std::shared_ptr<db::GridSteiner> node) {
    const auto &np = database.getLoc(*node);
    std::cout << "  " << node->layerIdx << ", " << np << ";  pinIdx : " << node->pinIdx << std::endl;
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
    size_t rsize = routeGuides.size();
    int lnum = database.getLayerNum() - 1;
    for (size_t i=0; i<rsize; i++) {
        const auto &g = routeGuides[i];
        int l = g.layerIdx;
        // if (l > 2) routeGuides.emplace_back(l-1, g);
        if (l < lnum) routeGuides.emplace_back(l+1, g);
        else routeGuides.emplace_back(l-1, g);
        db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_1, 1);
    }
}

void selectGuides(db::GridSteiner *node,
                    vector<int> &selected,
                    std::set<size_t> &breaks,
                    const vector<db::BoxOnLayer> &guides) {
    const auto &np = database.getLoc(*node);
    bool currentOut = true;
    for (size_t i=0; i<guides.size(); i++) {
        if (guides[i].layerIdx == node->layerIdx && guides[i].Contain(np)) {
            selected[i] += 1;
            node->distance = int(i);
            currentOut = false;
            break;
        }
    }

    if (currentOut) {
        node->distance = -1;
        if (node->parent) {
            const auto &pp = database.getLoc(*(node->parent));
            for (size_t i=0; i<guides.size(); i++) {
                if (guides[i].layerIdx == node->parent->layerIdx &&
                    guides[i].Contain(pp)) {
                    breaks.insert(i);
                    break;
                }
            }
        }
        for (auto c : node->children) {
            const auto &cp = database.getLoc(*c);
            for (size_t i=0; i<guides.size(); i++) {
                if (guides[i].layerIdx == c->layerIdx &&
                    guides[i].Contain(cp)) {
                    breaks.insert(i);
                    break;
                }
            }
        }
    }

    for (auto c : node->children) {
        if (c->isVio) selectGuides(c.get(), selected, breaks, guides);
    }
}

void LocalNet::creatAdaptiveRouteGuides() {
    const auto &mguides = dbNet.mergedGuides;
    size_t msize = mguides.size();
    vector<int> selected(mguides.size(), 0), svios;
    std::set<size_t> breaks;
    vector<db::BoxOnLayer> pickedGuides;
    selectGuides(pnetPins[0].get(), selected, breaks, mguides);

    if (!breaks.empty()) {
        size_t root = *(breaks.begin());
        breaks.erase(root);
        vector<size_t> curr, next, ends;
        curr.push_back(root);
        std::unordered_set<size_t> visitedGuides;
        std::unordered_map<size_t, size_t> parent;

        if (breaks.empty()) {
            for (auto p : pnetPins) {
                if (p->distance > -1 || p->pinIdx < 0) continue;
                DBU mindist = 10000, dist = 0;
                size_t minid = 0;
                const auto &ploc = database.getLoc(*p);
                for (size_t i=0; i<mguides.size(); i++) {
                    if (mguides[i].layerIdx != p->layerIdx) continue;
                    dist = utils::Dist(mguides[i], ploc);
                    if (dist < mindist) {
                        mindist = dist;
                        minid = i;
                    }
                }
                breaks.insert(minid);
                pickedGuides.push_back(mguides[minid]);
                pickedGuides.back().Update(ploc);
            }
        }

        while (!breaks.empty() && !curr.empty()) {
            for (size_t i : curr) {
                const auto &gi = mguides[i];
                for (size_t j=0; j<msize; j++) {
                    if (visitedGuides.count(j)) continue;
                    const auto &gj = mguides[j];
                    if (abs(gi.layerIdx-gj.layerIdx) == 1 &&
                        gi.HasIntersectWith(gj)) {
                        parent.insert({j, i});
                        visitedGuides.insert(j);
                        if (breaks.count(j)) {breaks.erase(j); ends.push_back(j);}
                        next.push_back(j);
                    }
                }
            }
            curr = move(next);
        }
        for (auto i : ends) {
            while (i != root) {
                selected[parent[i]]++;
                i = parent[i];
            }
        }
    }
    for (size_t i=0; i<mguides.size(); i++) {
        if (selected[i] > 0) {
            pickedGuides.push_back(mguides[i]);
            svios.push_back(dbNet.mergedVios[i]);
        }
    }
    // ensure connectivity
    for (size_t i=0; i<pnetPins.size(); i++) {
        auto p = pnetPins[i];
        if (p->distance > -1 || p->pinIdx >= 0) continue;
        int lid = p->layerIdx;
        db::BoxOnLayer patch;
        const auto &pp = database.getLoc(*p);
            printf(" Warning : net %d pnet %d pin %d not real pin but out of guide\n", idx, pnetIdx, int(i));
            for (const auto &g : pickedGuides) {
                if (g.Contain(pp) && abs(g.layerIdx-p->layerIdx) < 2) {
                    patch = g;
                    break;
                }
            }
            const auto &layer = database.getLayer(lid);
            DBU margin = layer.pitch * 20;
            int dir = 1 - layer.direction;
            patch.layerIdx = lid;
            patch[dir].low = max<DBU>(patch[dir].low, pp[dir]-margin);
            patch[dir].high = min<DBU>(patch[dir].high, pp[dir]+margin);
        // }
        pickedGuides.push_back(patch);
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
        for (auto &box : pickedGuides) {
            if (box.layerIdx == p->layerIdx && box.Contain(pp)) {
                if (box.Contain(np)) {
                    int dir = 1 - database.getLayerDir(box.layerIdx);
                    // DBU dist = pp[dir] - box[dir].low;
                    // if (dist > box[dir].range() * 0.3 || dist < box[dir].range() * 0.7)
                    //     break;
                    if (np[dir] < pp[dir]) box[dir].low = pp[dir];
                    else if (np[dir] > pp[dir]) box[dir].high = pp[dir];
                }
                break;
            }
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
    gridPinAccessBoxes.resize(pnetPins.size());
    for (size_t i=0; i<pnetPins.size(); i++) {
        // if (pnetPins[i]->children.empty() || pnetPins[i]->parent == nullptr) {
        //     if (pnetPins[i]->pinIdx < 0)
        //         printf("net %d pnetIdx %d pnet %ld is not pin but have no children\n", idx, pnetIdx, i);
        if (pnetPins[i]->pinIdx < 0 || pnetPins[i]->children.size() > 0)
            gridPinAccessBoxes[i].push_back(database.rangeSearch(pinAccessBoxes[i][0]));
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