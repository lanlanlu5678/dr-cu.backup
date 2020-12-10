#include <fstream>
#include "LocalNet.h"

void LocalNet::initGridBoxes() {
    // rangeSearch & slice routeGuides
    vector<vector<db::GridBoxOnLayer>> guides(database.getLayerNum());
    for (auto& guide : routeGuides) {
        const auto &gbox = database.rangeSearch(guide);
        if (database.isValid(gbox))
            guides[guide.layerIdx].push_back(gbox);
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
    // PARTIAL RIPUP
    initPNetPA();
    for (unsigned pinIdx = 0; pinIdx != numOfPins(); ++pinIdx) {
        pinAccessBoxes[pinIdx].clear(); 
        for (auto& box : gridPinAccessBoxes[pinIdx]) {
            pinAccessBoxes[pinIdx].push_back(database.getLoc(box));
        }
    }

    if (db::rrrIterSetting.constrainInGuide)
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

// PARTIAL RIPUP
void LocalNet::traverse(std::shared_ptr<db::GridSteiner> node) {
    bool isPin = false;
    for (auto c : node->children) {
        if (c->isVio)
            traverse(c);
        else
            isPin = true;
    }
    if ((node->children.empty() || isPin) && node != pnetPins[0])
        pnetPins.push_back(node);
}

void initIntersection(const db::GridSteiner &p,
                        vector<db::BoxOnLayer> &picked,
                        const vector<db::BoxOnLayer> &guides,
                        size_t oriSize) {

    const auto &np = database.getLoc(p);
    bool out = true;
    for (const auto &box : picked) {
        if (box.layerIdx == p.layerIdx && box.Contain(np)) {
            out = false;
            break;
        }
    }
    if (out) {
        std::set<int> otherGuides;
        for (size_t i=0; i<guides.size(); i++) {
            if (guides[i].layerIdx == p.layerIdx && guides[i].Contain(np)) {
                otherGuides.insert(i);
            }
        }
        for (int i : otherGuides) {
            for (size_t j=0; j<oriSize; j++) {
                const auto &its = guides[i].IntersectWith(picked[j]);
                if (its.IsValid())
                    picked.emplace_back(guides[i].layerIdx, its);
            }
        }
    }

    for (auto c : p.children) {
        if (c->isVio) initIntersection(*c, picked, guides, oriSize);
    }
}

void LocalNet::quickRG() {
    db::GridBoxOnLayer box;
    // vector<vector<int>> connList;
    vector<std::shared_ptr<db::GridSteiner>> cbox, obox;

    obox.push_back(pnetPins[0]);
    while (!obox.empty()) {
        auto node = obox.back();
        obox.pop_back();
        box.layerIdx = node->layerIdx;
        box.trackRange.Set(node->trackIdx);
        box.crossPointRange.Set(node->crossPointIdx);
        cbox.push_back(node);
        while (!cbox.empty()) {
            node = cbox.back();
            cbox.pop_back();
            for (auto c : node->children) {
                if (!c->isVio) continue;
                if (c->layerIdx == node->layerIdx) {
                    box.trackRange.FastUpdate(c->trackIdx);
                    box.crossPointRange.FastUpdate(c->crossPointIdx);
                    cbox.push_back(c);
                }
                else obox.push_back(c);
            }
        }
        box.trackRange.low -= db::rrrIterSetting.extraTracks;
        if (box.trackRange.low < 0) box.trackRange.low = 0;
        box.trackRange.high += db::rrrIterSetting.extraTracks;
        if (box.trackRange.high > database.getTrackLimit(box.layerIdx))
            box.trackRange.high = database.getTrackLimit(box.layerIdx);
        routeGuides.push_back(database.getLoc(box));
    }

    // if (db::rrrIterSetting.addDiffLayerGuides) {
    //     int oriSize = routeGuides.size();
    //     for (int i=0; i<oriSize; i++) {
    //         int l = routeGuides[i].layerIdx;
    //         if (l > 0)
    //             routeGuides.emplace_back(l-1, routeGuides[i]);
    //         if (l < database.getLayerNum()-1)
    //             routeGuides.emplace_back(l+1, routeGuides[i]);
    //     }
    // }
}

void LocalNet::guideRipup() {
    // check this pseudo net inside which group of ripup nets
    vector<db::BoxOnLayer> pickedGuides;
    std::set<int> gids;
    getGroups(*(pnetPins[0]), gids, pickedGuides);
    for (int guideIdx : gids) {
        // debug
        if (guideIdx < 0) {
            database.debugLock.lock();
            log() << "NET " << idx << "_" << pnetIdx << " has vio node out of groups" << std::endl;
            database.debugLock.unlock();
            continue;
        }
        for (int id : dbNet.ripupGuides[guideIdx]) {
            const auto &box = routeGuides[id];
            pickedGuides.push_back(box);
            // diff layer
            int l = box.layerIdx;
            if (dbNet.routeGuideVios[id] >= db::setting.diffLayerGuideVioThres) {
                if (l > 2) pickedGuides.emplace_back(l - 1, box);  // do not add to layers 0, 1
                if ((l + 1) < database.getLayerNum()) pickedGuides.emplace_back(l + 1, box);
                db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_1, 1);
            }
            if (dbNet.routeGuideVios[id] >= db::setting.diffLayerGuideVioThres * 2) {
                if (l > 3) pickedGuides.emplace_back(l - 2, box);  // do not add to layers 0, 1
                if ((l + 2) < database.getLayerNum()) pickedGuides.emplace_back(l + 2, box);
                db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::ADD_DIFF_LAYER_GUIDE_2, 1);
            }
        }
    }
    // if (gids.size() > 1) fixLongEdge(*(pnetPins[0]), pickedGuides); // if long edge
    // initIntersection(pnetPins, pickedGuides, routeGuides);
    initIntersection(*(pnetPins[0]), pickedGuides, routeGuides, pickedGuides.size());
    // debugBox = pickedGuides; // debug
    routeGuides = std::move(pickedGuides);
}

// void LocalNet::fixLongEdge(const db::GridSteiner &node, vector<db::BoxOnLayer> &picked) {
//     const auto &np = database.getLoc(node);
//     int nid = getGroup(np);
//     for (auto c : node.children) {
//         if (!c->isVio || (c->layerIdx != node.layerIdx))
//             continue;
//         const auto &cp = database.getLoc(*c);
//         int cid = getGroup(cp);
//         if (nid != cid) {
//             // db::BoxOnLayer box(node.layerIdx, np);
//             // box.FastUpdate(cp);
//             for (const auto &g : routeGuides)
//                 // if (g.HasStrictIntersectWith(box))
//                 if (g.Contain(np) && g.Contain(cp))
//                     picked.push_back(g);
//         }
//     }
//     for (auto c : node.children)
//         if (c->isVio) fixLongEdge(*c, picked);
// }

void LocalNet::getGroups(const db::GridSteiner &node,
                            std::set<int> &gids,
                            vector<db::BoxOnLayer> &picked) {
    const auto &np = database.getLoc(node);
    int nid = getGroup(np);
    gids.insert(nid);
    for (auto c : node.children) {
        if (!c->isVio || (c->layerIdx != node.layerIdx))
            continue;
        const auto &cp = database.getLoc(*c);
        int cid = getGroup(cp);
        if (nid != cid || cid == -1) {
            db::BoxOnLayer box(node.layerIdx, np);
            box.Update(cp);
            for (const auto &g : routeGuides) {
                if (g.layerIdx != node.layerIdx) continue;
                const auto &its = g.IntersectWith(box);
                if (its.IsValid() && (its.x.IsStrictValid() || its.y.IsStrictValid()))
                    picked.push_back(g);
            }
        }
    }
    
    for (auto c : node.children)
        if (c->isVio) getGroups(*c, gids, picked);
}

int LocalNet::getGroup(const utils::PointT<DBU> &p) {
    for (size_t i=0; i<dbNet.ripupGuides.size(); i++) {
        for (int id : dbNet.ripupGuides[i]) {
            if (routeGuides[id].Contain(p)) {
                return int(i);
            }
        }
    }
    return -1;
}

void LocalNet::initPNet(int numPitch) {

    pnetPins.push_back(dbNet.pnets[pnetIdx]);
    traverse(pnetPins[0]);

    pinAccessBoxes.clear();
    pinAccessBoxes.resize(pnetPins.size());
    if (pnetPins.size() == 1) {
        pnetPins[0]->isVio = false;
        return;   // one vio node won't ripup, safe
    }
    for (int i=0; i<pnetPins.size(); i++) {
        pinAccessBoxes[i].emplace_back(pnetPins[i]->layerIdx, database.getLoc(*(pnetPins[i])));
    }

    // routeGuides.clear();
    if (db::rrrIterSetting.quickFixMode) {
        routeGuides.clear();
        quickRG();
    }
    else {
        guideRipup();
    }
    // // debug
    // if (idx == 179558) {
    //     database.debugLock.lock();
    //     debugBox = routeGuides;
    //     printDebug();
    //     database.debugLock.unlock();
    // }

    gridTopo.clear();
}

void LocalNet::initPNetPA() {
    if (pnetIdx<0) return;

    vector<vector<db::GridBoxOnLayer>> gridPPins;

    gridPPins.resize(pnetPins.size());

    gridPPins[0].push_back({pnetPins[0]->layerIdx,
                        utils::IntervalT<int>(pnetPins[0]->trackIdx),
                        utils::IntervalT<int>(pnetPins[0]->crossPointIdx)});

    for (int i=1; i<pnetPins.size(); i++) {
        int pidx = pnetPins[i]->pinIdx;
        if (pnetPins[i]->children.empty() && pidx>0) {
            gridPPins[i] = move(gridPinAccessBoxes[pidx]);
        }
        else {
            gridPPins[i].push_back({pnetPins[i]->layerIdx,
                                utils::IntervalT<int>(pnetPins[i]->trackIdx),
                                utils::IntervalT<int>(pnetPins[i]->crossPointIdx)});
        }
    }
    gridPinAccessBoxes = move(gridPPins);
}

void LocalNet::printDebug() const {
    std::ofstream ofs(std::to_string(idx)+"_"+std::to_string(pnetIdx)+".txt");

    // ori route guides
    ofs << "ori route guides : " << dbNet.routeGuides.size() << std::endl;
    for (auto g : dbNet.routeGuides) {
        database.expandBox(g, db::rrrIterSetting.defaultGuideExpand);
        ofs << g << std::endl;
    }

    // picked route guides
    ofs << "picked route guides : " << debugBox.size() << std::endl;
    for (const auto &g : debugBox) {
        ofs << g << std::endl;
    }

    // pnet pin nodes
    ofs << "pnet pin nodes : " << pnetPins.size() << std::endl;
    for (size_t i=0; i<pnetPins.size(); i++) {
        const auto &np = database.getLoc(*(pnetPins[i]));
        ofs << pnetPins[i]->layerIdx << " " << np.x << " " << np.y << std::endl;
    }

    // vio nodes
    ofs << "vio nodes : " << dbNet.vioNodes.size() << std::endl;
    for (size_t i=0; i<dbNet.vioNodes.size(); i++) {
        const auto &np = database.getLoc(*(dbNet.vioNodes[i]));
        ofs << dbNet.vioNodes[i]->layerIdx << " " << np.x << " " << np.y << std::endl;
    }

    // ori grid points
    ofs << "ori grid points " << std::endl;
    dbNet.postOrderVisitGridTopo([&ofs](std::shared_ptr<db::GridSteiner> node){
    // postOrderVisitGridTopo([&ofs](std::shared_ptr<db::GridSteiner> node){
        const auto &p = database.getLoc(*node);
        ofs << node->layerIdx << " " << p.x << " " << p.y << std::endl;
    });

    // // pnet edges
    // ofs << std::endl;
    // ofs << "pnet edges" << std::endl;
    // dbNet.postOrderVisitGridTopo([&ofs](std::shared_ptr<db::GridSteiner> node){
    //     if (node->isVio) {
    //         const auto &np = database.getLoc(*node);
    //         const auto &pp = database.getLoc(*(node->parent));
    //         ofs << node->layerIdx << " " << np.x << " " << np.y << " ; "
    //             << node->parent->layerIdx << " " << pp.x << " " << pp.y << std::endl;
    //     }
    // });

    // vio groups
    ofs << std::endl;
    for (const auto &set : dbNet.ripupGuides) {
        for (int id : set) {
            ofs << id << ",";
        }
        ofs << std::endl;
    }
}