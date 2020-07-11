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

void LocalNet::makePseudo() {

    vector<vector<db::GridBoxOnLayer>> pGridBoxes;
    vector<int> relatedGuides(gridRouteGuides.size(), 0);
    vector<db::GridBoxOnLayer> trimedGuides;

    // INIT PINS & RELATED GUIDES
    pins.push_back(dbNet.pinsOfPNets[pseudoNetIdx]);
    if (pins[0]->isVio) {
        // PartialRipup::checkStartPin(pins[0]);
        PartialRipup::purge(pins[0], pins, relatedGuides, gridRouteGuides);
    }
    else {
        for (int i=0; i<gridRouteGuides.size(); i++) {
            if (gridRouteGuides[i].includePoint(*(pins[0]))) {
                relatedGuides[i] = 1;
                break;
            }
        }
        for (auto c : pins[0]->children)
            if (c->isVio) PartialRipup::purge(c, pins, relatedGuides, gridRouteGuides);
    }

    for (int i=0; i<relatedGuides.size(); i++) {
        if (relatedGuides[i] > 0) trimedGuides.push_back(gridRouteGuides[i]);
    }
    gridRouteGuides = std::move(trimedGuides);

    pGridBoxes.resize(pins.size());
    if (pins[0]->parent) pGridBoxes[0].push_back({pins[0]->layerIdx,
                            utils::IntervalT<int>(pins[0]->trackIdx),
                            utils::IntervalT<int>(pins[0]->crossPointIdx)});
    else pGridBoxes[0] = std::move(gridPinAccessBoxes[pins[0]->pinIdx]);
    for (int i=1; i<pins.size(); i++) {
        if (pins[i]->children.empty()) {
            // PartialRipup::checkEndPin(i,pseudoNetIdx,pins,dbNet);
            pGridBoxes[i] = std::move(gridPinAccessBoxes[pins[i]->pinIdx]);
        }
        else pGridBoxes[i].push_back({pins[i]->layerIdx,
                utils::IntervalT<int>(pins[i]->trackIdx),
                utils::IntervalT<int>(pins[i]->crossPointIdx)});
    }
    gridPinAccessBoxes = std::move(pGridBoxes);
    pinAccessBoxes.clear();
    for (auto &Boxes : gridPinAccessBoxes) {
        pinAccessBoxes.emplace_back();
        auto &pBoxes = pinAccessBoxes.back();
        for (auto box : Boxes) pBoxes.push_back(database.getLoc(box));
    }

    gridTopo.clear();
}

/*
    in case of repeatly push same guide, use vector<int> to mark
*/
void LocalNet::initPseudo() {

    vector<int> guideIdx(routeGuides.size());
    vector<db::BoxOnLayer> relatedGuides;

    // traverse pseudo net
    // TODO : Diff-layer guides & trim
    pins.push_back(dbNet.pinsOfPNets[pseudoNetIdx]);
        
    if (pins[0]->isVio) PartialRipup::traversePNet(pins[0], pins, routeGuides, guideIdx);
    else
        for (auto c : pins[0]->children) {
            if (c->isVio) PartialRipup::traversePNet(c, pins, routeGuides, guideIdx);
        }
    
    // UPDATE GUIDES
    for (auto p : pins) {
        auto pbox = database.getLoc(*p);
        for (int i=0; i<routeGuides.size(); i++) {
            auto &box = routeGuides[i];
            if (box.layerIdx == p->layerIdx &&
                box.x.low < pbox.x && box.x.high > pbox.x &&
                box.y.low < pbox.y && box.y.high > pbox.y) {
                guideIdx[i] = 1;
                break;
            }
        }
    }
    for (int i=0; i<guideIdx.size(); i++)
        if (guideIdx[i] > 0) relatedGuides.push_back(routeGuides[i]);
    routeGuides = std::move(relatedGuides);

    gridTopo.clear();
}

/*
    pins[0]'s position must be fix to keep edge valid.
*/
void LocalNet::initPPinBoxes() {

    vector<vector<db::GridBoxOnLayer>> ppins;

    ppins.resize(pins.size());
    pinAccessBoxes.clear();
    pinAccessBoxes.resize(pins.size());

    ppins[0].push_back({pins[0]->layerIdx,
                        utils::IntervalT<int>(pins[0]->trackIdx),
                        utils::IntervalT<int>(pins[0]->crossPointIdx)});

    for (int i=1; i<pins.size(); i++) {
        if (pins[i]->pinIdx < 0) {
            ppins[i].push_back({pins[i]->layerIdx,
                                utils::IntervalT<int>(pins[i]->trackIdx),
                                utils::IntervalT<int>(pins[i]->crossPointIdx)});
        }
        else ppins[i] = move(gridPinAccessBoxes[pins[i]->pinIdx]);
    }
    gridPinAccessBoxes = move(ppins);
}