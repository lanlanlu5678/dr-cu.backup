#include "PostRoute.h"
#include "PinTapConnector.h"

inline bool isEnclosed(const utils::BoxT<DBU> &viaBox, const db::BoxOnLayer &pbox) {
    return !((viaBox.x.low < pbox.x.low) ||
            (viaBox.y.low < pbox.y.low) ||
            (viaBox.x.high > pbox.x.high) ||
            (viaBox.y.high > pbox.y.high));
}

bool PostRoute::handlePinVia(const db::GridSteiner &tap) {
    const db::ViaType *type = nullptr;
    if (tap.viaType) type = tap.viaType;
    else if (tap.parent == nullptr)
        type = tap.children[0]->viaType;
    if (tap.layerIdx > 0 || type == nullptr)
        return false;

    const auto &pabs = dbNet.pinAccessBoxes[tap.pinIdx];
    const auto &tapLoc = database.getLoc(tap);
    const auto &viaBox = type->getShiftedBotMetal(tapLoc);
    utils::BoxT<DBU> connBox;
    vector<size_t> closeIdx;
    DBU space = database.getLayer(0).getParaRunSpace(0), margin = 0,
        width = database.getLayer(0).width;
    int dir = 0;

    for (size_t i=0; i<pabs.size(); i++) {
        const auto &pbox = pabs[i];
        if (pbox.layerIdx > 0) continue;
        if (utils::Dist(pbox, viaBox) < space) {
            closeIdx.push_back(i);
            if (viaBox.HasIntersectWith(pbox)) {
                if (isEnclosed(viaBox, pbox)) return true;
                connBox.x.low = std::min<DBU>(connBox.x.low, pbox.x.low);
                connBox.x.high = std::max<DBU>(connBox.x.high, pbox.x.high);
                connBox.y.low = std::min<DBU>(connBox.y.low, pbox.y.low);
                connBox.y.high = std::max<DBU>(connBox.y.high, pbox.y.high);
            }
        }
    }
    if (closeIdx.empty()) return false;

    if (safeShift(tap, closeIdx) ||
        ((database.countFixedMetals(viaBox, dbNet.idx) > 0) && shiftPinVia(tap, closeIdx))) {
        db::routeStat.increment(db::RouteStage::POST, db::MiscRouteEvent::MOVE_PIN_VIA, 1);
        return true;
    }

    if (!connBox.IsValid()) return false;
    const auto &its = viaBox.IntersectWith(connBox);
    if (its.x.range() < width && its.y.range() < width) {
        if (viaBox.x.range() > viaBox.y.range()) {
            margin = its.x.range();
        }
        else {
            dir = 1;
            margin = its.y.range();
        }
        margin = width - margin;
        if (margin > 0) {
            connBox = viaBox;
            if (its[dir].low == viaBox[dir].low)
                connBox[dir].low -= margin;
            else
                connBox[dir].high += margin;
            database.writeDEFFillRect(dbNet, connBox, 0);
            pinViaPatches[tap.pinIdx] = connBox;
        }
    }
    return true;
}

bool PostRoute::shiftPinVia(const db::GridSteiner &tap, const vector<size_t> &closeIdx) {
    const auto &pabs = dbNet.pinAccessBoxes[tap.pinIdx];
    const auto &tapLoc = database.getLoc(tap);
    std::shared_ptr<db::GridSteiner> nei1, nei2;
    int dir = 1 - database.getLayerDir(1), ntdir = 0;
    std::multimap<int, size_t> canPior;
    vector<utils::PointT<DBU>> canLocs;
    vector<const db::ViaType *> canTypes;
    const db::ViaType *type = nullptr;

    // get candidates
    for (size_t id : closeIdx) {
        const auto &pbox = pabs[id];
        type = database.getBestViaTypeForShift(pbox);
        if (type == nullptr) continue;
        const auto &bot = type->getShiftedBotMetal(tapLoc);
        auto cand = tapLoc;
        for (int tdir=0; tdir<2; tdir++) {
            if (bot[tdir].low < pbox[tdir].low)
                cand[tdir] += (pbox[tdir].low - bot[tdir].low);
            else if (bot[tdir].high > pbox[tdir].high)
                cand[tdir] += (pbox[tdir].high - bot[tdir].high);
        }
        int pior = abs(cand[1-dir] - tapLoc[1-dir]);
        if (pior) {
            canLocs.push_back(cand);
            canTypes.push_back(type);
            canPior.insert({pior, canLocs.size()-1});
        }
    }
    // if (dbNet.getName() == "net12727") {
    //     printf("\n net12727 tap loc : %ld, %ld\n", tapLoc.x, tapLoc.y);
    //     printf(" close pin boxes :\n");
    //     for (size_t i : closeIdx)
    //         std::cout << pabs[i] << std::endl;
    //     printf(" can locs : \n");
    //     for (const auto &loc : canLocs)
    //         std::cout << loc << std::endl;
    // }
    if (canLocs.empty()) return false;

    // check topo
    if (tap.parent) {
        nei1 = tap.parent;
        nei2 = nei1->parent;
    }
    else {
        nei1 = tap.children[0];
        nei2 = nei1->children[0];
    }
    ntdir = nei1->trackIdx == nei2->trackIdx ? 1-dir : dir;
    const auto &n1 = database.getLoc(*nei1);
    const auto &n2 = database.getLoc(*nei2);
    const auto &n12 = PinTapConnector::getLinkMetal({n1, n2}, 1);

    // get confict boxes
    const auto &tapu = database.getUpper(tap);
    db::GridBoxOnLayer qGrid(tapu.layerIdx, {tapu.trackIdx-2, tapu.trackIdx+2},
                        {tapu.crossPointIdx-2, tapu.crossPointIdx+2});
    const auto &routedBoxes = database.getRoutedBox(qGrid, dbNet.idx);
    const auto &fixedBoxes = database.getFixedBox(database.getLoc(qGrid), dbNet.idx);

    // check vio
    vector<db::BoxOnLayer> boxesToCheck;
    vector<utils::SegmentT<DBU>> link;
    for (const auto &pair : canPior) {
        size_t i = pair.second;
        boxesToCheck.clear();
        boxesToCheck.emplace_back(1, canTypes[i]->getShiftedTopMetal(canLocs[i]));

        // get link to new via
        link.clear();
        if (pair.first < 3) 
            link.emplace_back(tapLoc, canLocs[i]);
        else {
            auto turn = tapLoc;
            turn[dir] = canLocs[i][dir];
            link.emplace_back(tapLoc, turn);
            link.emplace_back(turn, canLocs[i]);
        }
        for (const auto &wire : link) {
            boxesToCheck.emplace_back(1, PinTapConnector::getLinkMetal(wire, 1));
        }
        
        // patch
        if (canLocs[i][ntdir] != n1[ntdir]) {
            auto its = boxesToCheck[0].IntersectWith(n12);
            its[ntdir] = its[ntdir].UnionWith(n12[ntdir]);
            boxesToCheck.emplace_back(1, its);
        }

        int vio = 0;
        for (const auto &box : boxesToCheck) {
            const auto &regions = database.getAccurateMetalRectForbidRegions(box);
            vio += database.countOvlp(box, regions, routedBoxes);
            vio += database.countOvlp(box, regions, fixedBoxes);
            if (vio > 0) {
                // if (dbNet.getName() == "net12727") {
                //     printf("\n net12727 tap loc : %ld, %ld ; with neighbours = %d\n", tapLoc.x, tapLoc.y, neighbours);
                //     printf(" can loc : %ld, %ld and via upper box :", canLocs[i].x, canLocs[i].y);
                //     std::cout << canTypes[i]->top;
                //     printf("\n box checking : ");
                //     std::cout << box << std::endl;
                //     printf(" checking box forbid regions :\n");
                //     for (const auto &rbox : regions)
                //         std::cout << rbox << std::endl;
                //     printf(" routed boxes :\n");
                //     for (const auto &rbox : routedBoxes)
                //         std::cout << rbox << std::endl;
                //     printf(" routed box forbid regions :\n");
                //     for (const auto &rbox : routedBoxes) {
                //         const auto &rfbr = database.getAccurateMetalRectForbidRegions({1, rbox});
                //         for (const auto &rfb : rfbr)
                //             std::cout << rfb << std::endl;
                //     }
                //     // printf(" fixed boxes :\n");
                //     // for (const auto &fbox : fixedBoxes)
                //     //     std::cout << fbox << std::endl;
                //     printf("\n");
                // }
                break;
            }
        }
        // succ
        if (vio == 0) {
            // update topo
            if (tap.parent) {
                for (auto it=nei1->children.begin(); it!=nei1->children.end(); it++) {
                    if ((*it).get() == &tap) {
                        *it = nei1->children.back();
                        nei1->children.pop_back();
                        break;
                    }
                }
            }
            else {
                nei1->parent = nullptr;
                for (auto it=dbNet.gridTopo.begin(); it!=dbNet.gridTopo.end(); it++) {
                    if ((*it).get() == &tap) {
                        *it = nei1;
                        break;
                    }
                }
            }
            nei1->pinIdx = tap.pinIdx;
            nei1->fakePin = true;

            // store changes
            linkToPins[nei1] = move(link);
            linkViaToPins[nei1] = std::make_pair(0, canLocs[i]);
            linkViaTypes[nei1] = canTypes[i];
            return true;
        }
    }
    return false;
}

bool PostRoute::safeShift(const db::GridSteiner &tap, const vector<size_t> &closeIdx) {
    int dir = 1 - database.getLayerDir(1), mindist = 1000, dist = 0;
    const auto &tapLoc = database.getLoc(tap);
    const auto &cutLayer = database.getCutLayer(0);
    const auto &pabs = dbNet.pinAccessBoxes[tap.pinIdx];
    const db::ViaType *canType = nullptr;
    utils::PointT<DBU> canLoc(tapLoc.x, tapLoc.y);

    for (size_t id : closeIdx) {
        const auto &pbox = pabs[id];
        for (const auto &type : cutLayer.allViaTypes) {
            if (type.topDir != dir) break;
            const auto &bot = type.getShiftedBotMetal(tapLoc);
            if ((bot[1-dir].low < pbox[1-dir].low) || (bot[1-dir].high > pbox[1-dir].high) ||
                (bot[dir].range() > pbox[dir].range()))
                continue;
            if (bot[dir].low < pbox[dir].low)
                dist = pbox[dir].low - bot[dir].low;
            else if (bot[dir].high > pbox[dir].high)
                dist = pbox[dir].high - bot[dir].high;
            if (abs(dist) < abs(mindist)) {
                mindist = dist;
                canType = &type;
            }
        }
    }
    if (canType == nullptr)
        return false;
    
    canLoc[dir] += mindist;
    auto qgr = database.getSurroundingCrossPoint(1, canLoc[dir]);
    const auto &gu = database.getUpper(tap);
    qgr.Update(gu.crossPointIdx);
    if (database.hasOtherRoutedWireOnTrack(dbNet.idx, gu.layerIdx, gu.trackIdx, qgr) ||
        database.hasOtherRoutedViaOnTrack(dbNet.idx, gu.layerIdx, gu.trackIdx, qgr))
        return false;

    // update topo
    std::shared_ptr<db::GridSteiner> nei1;
    if (tap.parent) {
        nei1 = tap.parent;
        for (auto it=nei1->children.begin(); it!=nei1->children.end(); it++) {
            if ((*it).get() == &tap) {
                *it = nei1->children.back();
                nei1->children.pop_back();
                break;
            }
        }
    }
    else {
        nei1 = tap.children[0];
        nei1->parent = nullptr;
        for (auto it=dbNet.gridTopo.begin(); it!=dbNet.gridTopo.end(); it++) {
            if ((*it).get() == &tap) {
                *it = nei1;
                break;
            }
        }
    }
    nei1->pinIdx = tap.pinIdx;
    nei1->fakePin = true;

    // store changes
    linkToPins[nei1].emplace_back(tapLoc, canLoc);
    linkViaToPins[nei1] = std::make_pair(0, canLoc);
    linkViaTypes[nei1] = canType;
    return true;
}