#include "PreRoute.h"
#include "PartialRipup.h"

db::RouteStatus PreRoute::run(int numPitchForGuideExpand) {
    // expand guides uniformally
    auto& guides = localNet.routeGuides;
    for (int i = 0; i < guides.size(); ++i) {
        int expand = localNet.dbNet.routeGuideVios[i] ? numPitchForGuideExpand : db::setting.defaultGuideExpand;
        database.expandBox(guides[i], numPitchForGuideExpand);
    }

    // expand guides by cross layer connection
    expandGuidesToMargin();

    db::RouteStatus status = db::RouteStatus::SUCC_NORMAL;
    if (localNet.numOfPins() < 2) {
        status = db::RouteStatus::SUCC_ONE_PIN;
    } else {
        // expand guides to cover pin
        status = expandGuidesToCoverPins();
        if (db::isSucc(status)) {
            // init localNet and check
            // partial ripup
            if (localNet.pseudoNetIdx>=0 && !localNet.gridTopo.empty()) localNet.initGuidesAndPins();
            localNet.initGridBoxes();
            localNet.initConn(localNet.gridPinAccessBoxes, localNet.gridRouteGuides);
            localNet.initNumOfVertices();

            if (!localNet.checkPin()) {
                status = db::RouteStatus::FAIL_PIN_OUT_OF_GRID;
            } else if (!localNet.checkPinGuideConn()) {
                status = db::RouteStatus::FAIL_DETACHED_PIN;
            } else if (!checkGuideConnTrack()) {
                status = db::RouteStatus::FAIL_DETACHED_GUIDE;
            }
        }
    }
    printWarnMsg(status, localNet.dbNet);
    return status;
}

db::RouteStatus PreRoute::runIterative() {
    db::RouteStatus status = run(db::rrrIterSetting.defaultGuideExpand);

    int iter = 0;
    int numPitchForGuideExpand = db::rrrIterSetting.defaultGuideExpand;
    utils::timer singleNetTimer;
    while (status == +db::RouteStatus::FAIL_DETACHED_GUIDE && iter < db::setting.guideExpandIterLimit) {
        iter++;
        numPitchForGuideExpand += iter;

        status = run(numPitchForGuideExpand);
    }

    const std::string name = localNet.getName();
    if (iter >= 1) {
        log() << "Warning: Net " << name << " expands " << iter << " iterations"
              << ", which takes " << singleNetTimer.elapsed() << " s in total." << std::endl;
    }

    switch (status) {
        case +db::RouteStatus::FAIL_DETACHED_GUIDE:
            database.debugPrintLock.lock();
            log() << "Error: Exceed the guideExpandIterLimit, but Net " << localNet.idx
            << "-" << localNet.pseudoNetIdx << " still FAIL_DETACHED_GUIDE\n";
            if (--db::rrrIterSetting.writeNetNum > 0)
                localNet.writePNet("guide_");
            database.debugPrintLock.unlock();
            break;
        case +db::RouteStatus::FAIL_DETACHED_PIN:
            database.debugPrintLock.lock();
            log() << "Error: Net " << localNet.idx
            << "-" << localNet.pseudoNetIdx << " FAIL_DETACHED_PIN\n";
            if (--db::rrrIterSetting.writeNetNum > 0)
                localNet.writePNet("pin_");
            database.debugPrintLock.unlock();
            break;
        default:
            break;
    }

    db::routeStat.increment(db::RouteStage::PRE, status);
    return status;
}

void PreRoute::expandGuidesToMargin() {
    vector<db::BoxOnLayer>& guides = localNet.routeGuides;
    vector<vector<int>> crossLayerConn(guides.size());
    for (unsigned g1 = 0; g1 < guides.size(); g1++) {
        db::BoxOnLayer& box1 = guides[g1];
        for (unsigned g2 = g1 + 1; g2 < guides.size(); g2++) {
            db::BoxOnLayer& box2 = guides[g2];
            if (abs(box1.layerIdx - box2.layerIdx) == 1 && box1.HasIntersectWith(box2)) {
                crossLayerConn[g1].push_back(g2);
            }
        }
    }

    for (unsigned g1 = 0; g1 < guides.size(); g1++) {
        for (auto g2 : crossLayerConn[g1]) {
            Dimension dir1 = database.getLayerDir(guides[g1].layerIdx);
            Dimension dir2 = database.getLayerDir(guides[g2].layerIdx);

            guides[g1][dir2].Update(guides[g2][dir2].low);
            guides[g1][dir2].Update(guides[g2][dir2].high);
            guides[g2][dir1].Update(guides[g1][dir1].low);
            guides[g2][dir1].Update(guides[g1][dir1].high);
        }
    }
}

db::RouteStatus PreRoute::expandGuidesToCoverPins() {
    db::RouteStatus status = db::RouteStatus::SUCC_NORMAL;
    for (int i = 0; i < localNet.numOfPins(); ++i) {
        int bestAB = -1;
        int bestGuide = -1;
        DBU bestDist = std::numeric_limits<DBU>::max();
        for (int j = 0; j < localNet.pinAccessBoxes[i].size(); ++j) {
            const auto& accessBox = localNet.pinAccessBoxes[i][j];
            for (int k = 0; k < localNet.routeGuides.size(); ++k) {
                const auto& guide = localNet.routeGuides[k];
                if (abs(guide.layerIdx - accessBox.layerIdx) > 1 || !guide.IsValid()) continue;
                DBU dist = Dist(accessBox, guide);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestAB = j;
                    bestGuide = k;
                    if (bestDist == 0) {
                        break;
                    }
                }
            }
        }
        if (bestDist > 0) {
            if (bestAB == -1) {
                return db::RouteStatus::FAIL_DETACHED_PIN;
            } else {
                localNet.routeGuides[bestGuide] = {
                    localNet.routeGuides[bestGuide].layerIdx,
                    localNet.routeGuides[bestGuide].UnionWith(localNet.pinAccessBoxes[i][bestAB])};
                    db::routeStat.increment(db::RouteStage::PRE, db::MiscRouteEvent::FIX_DETACHED_PIN, 1);
            }
        }
    }
    return status;
}

bool PreRoute::checkGuideConnTrack() const {
    // init to all false
    vector<char> pinVisited(localNet.pinGuideConn.size(), false);
    vector<vector<char>> guideVisited(localNet.guideConn.size());
    for (int i = 0; i < localNet.guideConn.size(); i++) {
        guideVisited[i].resize(localNet.gridRouteGuides[i].trackRange.range() + 1, false);
    }

    // define recursive std::functions VisitPin & VisitGuide
    std::function<void(int)> VisitPin;
    std::function<void(int, utils::IntervalT<int>)> VisitGuide;
    VisitPin = [&](int pinIdx) {
        pinVisited[pinIdx] = true;
        for (const auto& ga : localNet.pinGuideConn[pinIdx]) {
            int guideIdx = ga.first;
            utils::IntervalT<int> trackInterval = localNet.gridRouteGuides[ga.first].trackRange.IntersectWith(
                localNet.gridPinAccessBoxes[pinIdx][ga.second].trackRange);
            bool isVisited = true;
            int trackLowIdx = localNet.gridRouteGuides[ga.first].trackRange.low;
            for (int t = trackInterval.low - trackLowIdx; t <= trackInterval.high - trackLowIdx; t++) {
                if (guideVisited[ga.first][t] == false) {
                    isVisited = false;
                    guideVisited[ga.first][t] = true;
                }
            }
            if (!isVisited) {
                VisitGuide(guideIdx, trackInterval);
            }
        }
    };
    VisitGuide = [&](int guideIdx, utils::IntervalT<int> trackRange) {
        for (int adjGuideIdx : localNet.guideConn[guideIdx]) {
            // to upper layer
            if (localNet.gridRouteGuides[guideIdx].layerIdx < localNet.gridRouteGuides[adjGuideIdx].layerIdx) {
                utils::IntervalT<int> upperCpInterval =
                    database.getLayer(localNet.gridRouteGuides[guideIdx].layerIdx).getUpperCrossPointRange(trackRange);
                if (localNet.gridRouteGuides[adjGuideIdx].crossPointRange.IntersectWith(upperCpInterval).IsValid()) {
                    db::ViaBox viaBox = database.getViaBoxBetween(localNet.gridRouteGuides[guideIdx],
                                                                  localNet.gridRouteGuides[adjGuideIdx]);
                    bool isVisited = true;
                    int trackLowIdx = localNet.gridRouteGuides[adjGuideIdx].trackRange.low;
                    for (int t = viaBox.upper.trackRange.low - trackLowIdx;
                         t <= viaBox.upper.trackRange.high - trackLowIdx;
                         t++) {
                        if (guideVisited[adjGuideIdx][t] == false) {
                            isVisited = false;
                            guideVisited[adjGuideIdx][t] = true;
                        }
                    }
                    if (!isVisited) {
                        VisitGuide(adjGuideIdx, viaBox.upper.trackRange);
                    }
                }
            }
            if (localNet.gridRouteGuides[guideIdx].layerIdx > localNet.gridRouteGuides[adjGuideIdx].layerIdx) {
                utils::IntervalT<int> lowerCpInterval =
                    database.getLayer(localNet.gridRouteGuides[guideIdx].layerIdx).getLowerCrossPointRange(trackRange);
                if (localNet.gridRouteGuides[adjGuideIdx].crossPointRange.IntersectWith(lowerCpInterval).IsValid()) {
                    db::ViaBox viaBox = database.getViaBoxBetween(localNet.gridRouteGuides[adjGuideIdx],
                                                                  localNet.gridRouteGuides[guideIdx]);
                    bool isVisited = true;
                    int trackLowIdx = localNet.gridRouteGuides[adjGuideIdx].trackRange.low;
                    for (int t = viaBox.lower.trackRange.low - trackLowIdx;
                         t <= viaBox.lower.trackRange.high - trackLowIdx;
                         t++) {
                        if (guideVisited[adjGuideIdx][t] == false) {
                            isVisited = false;
                            guideVisited[adjGuideIdx][t] = true;
                        }
                    }
                    if (!isVisited) {
                        VisitGuide(adjGuideIdx, viaBox.lower.trackRange);
                    }
                }
            }
        }
        for (const auto& pa : localNet.guidePinConn[guideIdx]) {
            utils::IntervalT<int> trackInterval =
                trackRange.IntersectWith(localNet.gridPinAccessBoxes[pa.first][pa.second].trackRange);
            if (trackInterval.IsValid() && !pinVisited[pa.first]) {
                VisitPin(pa.first);
            }
        }
    };

    // DFS from pin 0
    VisitPin(0);

    return all_of(pinVisited.begin(), pinVisited.end(), [](bool visited) { return visited; });
}
