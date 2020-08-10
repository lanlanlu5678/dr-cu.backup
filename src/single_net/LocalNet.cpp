#include "LocalNet.h"
#include "PartialRipup.h"

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
    // partial ripup
    if (pseudoNetIdx >= 0) initPPinBoxes();
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

/*
    in case of repeatly push same guide, use vector<int> to mark
*/
void LocalNet::initGuidesAndPins() {

    vector<int> guideVio(routeGuides.size(), 0);
    int d;
    // auto getOverlap = [](db::BoxOnLayer &path, db::BoxOnLayer &guide, db::BoxOnLayer &ol)->bool {
    //     DBU lx = max(path.lx(), guide.lx());
    //     DBU hx = min(path.hx(), guide.hx());
    //     DBU ly = max(path.ly(), guide.ly());
    //     DBU hy = min(path.hy(), guide.hy());
    //     if (lx<=hx && ly<=hy) {
    //         auto dir = database.getLayerDir(guide.layerIdx);
    //         DBU margin = database.getPitch(guide.layerIdx) * db::rrrIterSetting.defaultGuideExpand;
    //         if (dir == 0) {
    //             ol.y.low = ly - margin;
    //             ol.y.high = hy + margin;
    //             if (ol.y.low < guide.y.low) ol.y.low = guide.y.low;
    //             if (ol.y.high > guide.y.high) ol.y.high = guide.y.high;
    //             ol.x.low = guide.x.low;
    //             ol.x.high = guide.x.high;
    //         }
    //         else {
    //             ol.x.low = lx - margin;
    //             ol.x.high = hx + margin;
    //             if (ol.x.low < guide.x.low) ol.x.low = guide.x.low;
    //             if (ol.x.high > guide.x.high) ol.x.high = guide.x.high;
    //             ol.y.low = guide.y.low;
    //             ol.y.high = guide.y.high;
    //         }
    //         return true;
    //     }
    //     else return false;
    // };

    pins.push_back(dbNet.pinsOfPNets[pseudoNetIdx]);

    // find end pins
    if (pins[0]->isVio) PartialRipup::checkGuidesAndPins(pins[0], pins, routeGuides, guideVio);
    else
        for (auto c : pins[0]->children) {
            if (c->isVio) PartialRipup::checkGuidesAndPins(c, pins, routeGuides, guideVio);
        }
    
    // find route guides of start pin
    auto p = database.getLoc(*pins[0]);
    for (int i=0; i<routeGuides.size(); i++) {
        d = PartialRipup::inGuide(pins[0]->layerIdx, p.x, p.y, routeGuides[i]);
        if (d > -3) guideVio[i] |= (1<<(d+2));
    }

    // find end pins' route guides
    for (int i=1; i<pins.size(); i++){
        p = database.getLoc(*(pins[i]));
        for (int j=0; j<routeGuides.size(); j++) {
            d = PartialRipup::inGuide(pins[i]->layerIdx, p.x, p.y, routeGuides[j]);
            if (d > -3) guideVio[j] |= (1<<(d+2));
        }
    }
    // update route guides
    for (int i=0; i<guideVio.size(); i++) {
        int j = routeGuides[i].layerIdx;
        if ((guideVio[i] & 24) && (j+2) < database.getLayerNum()) relatedGuides.emplace_back(j+2, routeGuides[i]);
        if ((guideVio[i] & 28) && (j+1) < database.getLayerNum())relatedGuides.emplace_back(j+1, routeGuides[i]);
        if (guideVio[i] & 14) relatedGuides.emplace_back(j, routeGuides[i]);
        if ((guideVio[i] & 7) && (j-1) > 0) relatedGuides.emplace_back(j-1, routeGuides[i]);
        if ((guideVio[i] & 3) && (j-2) > 0) relatedGuides.emplace_back(j-2, routeGuides[i]);
        if (guideVio[i] > 0) pickedGuides.push_back(routeGuides[i]);
    }
    // for (auto &path : paths) {
    //     for (auto &guide : routeGuides) {
    //         l = path.layerIdx;
    //         diff = l - guide.layerIdx;
    //         if (abs(diff)>=db::rrrIterSetting.diffLayerBound || diff<-1) continue;
    //         if (getOverlap(path, guide, overlap)) {
    //             if (diff<2 && l+1<L) relatedGuides.emplace_back(l+1, overlap);
    //             if (diff==2 || l+1==L) relatedGuides.emplace_back(l-1, overlap);
    //             relatedGuides.emplace_back(l, overlap);
    //             pickedGuides.emplace_back(guide.layerIdx, overlap);
    //         }
    //     }
    // }
    // routeGuides = std::move(relatedGuides);
    std::swap(routeGuides, relatedGuides);
    gridTopo.clear();

    PartialRipup::checkPaths(pins[0], paths);
}

/*
    start pin's position must be fix to keep edges to it parent & no-vio children valid.
    posistion of end pins' with children must be fixed
*/
void LocalNet::initPPinBoxes() {

    vector<vector<db::GridBoxOnLayer>> ppins;
    int pidx;

    ppins.resize(pins.size());
    pinAccessBoxes.clear();
    pinAccessBoxes.resize(pins.size());

    ppins[0].push_back({pins[0]->layerIdx,
                        utils::IntervalT<int>(pins[0]->trackIdx),
                        utils::IntervalT<int>(pins[0]->crossPointIdx)});

    for (int i=1; i<pins.size(); i++) {
        if (!(pins[i]->children.empty())) {
            ppins[i].push_back({pins[i]->layerIdx,
                                utils::IntervalT<int>(pins[i]->trackIdx),
                                utils::IntervalT<int>(pins[i]->crossPointIdx)});
        }
        else {
            pidx = pins[i]->pinIdx;
            if (pidx < 0 || pidx > gridPinAccessBoxes.size())
                log() << "ERROR : TERMINAL " << i << " OF PNET " << idx << "_" << pseudoNetIdx << " IS NOT VALID PIN\n"
                        << "       PINIDX : " << pidx << ";     CHILDREN : " << pins[i]->children.size() << "\n";
            if (gridPinAccessBoxes[pidx].empty())
                log() << "ERROR : A PIN WITHOUT ACCESS BOXES\n";
            ppins[i] = move(gridPinAccessBoxes[pidx]);
        }
    }
    gridPinAccessBoxes = move(ppins);
}

void LocalNet::writePNet(const std::string &prefix) {
    std::ofstream ofs;
    db::BoxOnLayer box;
    std::string filename = "failed_nets/"+prefix+std::to_string(idx)+"_"+std::to_string(pseudoNetIdx)+".log";
    ofs.open(filename);
    ofs << "Net " << getName() << " (idx = " << idx << ") with " << pins.size() << " pins " << std::endl;
    for (int i=0; i<pins.size(); i++) {
        ofs << "pin " << i << " p; isVio " << pins[i]->isVio << "; child " << !(pins[i]->children.empty()) << std::endl;
        auto p = database.getLoc(*(pins[i]));
        box.Set(pins[i]->layerIdx, p);
        ofs << box << std::endl;
    }

    ofs << relatedGuides.size() << " route guides" << std::endl;
    for (auto &guide : relatedGuides) {
        ofs << guide << std::endl;
    }

    ofs << paths.size() << " paths" << std::endl;
    for (auto &path : paths) {
        ofs << path << std::endl;
    }

    ofs << pickedGuides.size() << " related guides" << std::endl;
    for (auto &guide : pickedGuides) {
        ofs << guide << std::endl;
    }

    ofs << std::endl;
    PartialRipup::plotPNet(ofs, pins[0]);
    ofs << std::endl;
    ofs.close();
}