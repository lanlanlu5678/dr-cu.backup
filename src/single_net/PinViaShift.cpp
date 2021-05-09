#include "PostRoute.h"
#include "PinTapConnector.h"

inline bool isEnclosed(const utils::BoxT<DBU> &viaBox, const db::BoxOnLayer &pbox) {
    return !((viaBox.x.low < pbox.x.low) ||
            (viaBox.y.low < pbox.y.low) ||
            (viaBox.x.high > pbox.x.high) ||
            (viaBox.y.high > pbox.y.high));
}

void patchPinVia(const utils::BoxT<DBU> &pin, utils::BoxT<DBU> viabot,
                    vector<utils::BoxT<DBU>> &pvp, db::Net &net) {
    DBU margin = 0, width = database.getLayer(0).width;
    int dir = 0;
    const auto &its = viabot.IntersectWith(pin);
    if (its.x.range() < width && its.y.range() < width) {
        if (viabot.x.range() > viabot.y.range()) {
            margin = its.x.range();
        }
        else {
            dir = 1;
            margin = its.y.range();
        }
        margin = width - margin;
        if (margin > 0) {
            if (its[dir].low == viabot[dir].low)
                viabot[dir].low -= margin;
            else
                viabot[dir].high += margin;
            database.writeDEFFillRect(net, viabot, 0);
            pvp.push_back(viabot);
        }
    }
}

bool PostRoute::handlePinVia(const db::GridSteiner &tap) {
    const db::ViaType *type = nullptr;
    if (tap.parent) {
        type = tap.viaType;
    }
    else if (tap.children.size() > 0) {
        type = tap.children[0]->viaType;
    }
    if (tap.layerIdx > 0 || type == nullptr) return false;

    const auto &pboxes = dbNet.pinAccessBoxes[tap.pinIdx];
    const auto &vialoc = database.getLoc(tap);
    const auto &viabot = type->getShiftedBotMetal(vialoc);
    utils::BoxT<DBU> connBox;
    vector<size_t> closeIdx;
    DBU patchSpace = database.getLayer(0).getParaRunSpace(0);
    // // debug
    // if (dbNet.getName() == "net9187" && vialoc.x > 1335100)
    //     printf("    enter handlePinVia\n");
    for (size_t i=0; i<pboxes.size(); i++) {
        const auto &pbox = pboxes[i];
        if (pbox.layerIdx > 0 ||
            utils::Dist(pbox, viabot) >= patchSpace)
            continue;
        closeIdx.push_back(i);
        if (pbox.HasIntersectWith(viabot)) {
            if (isEnclosed(viabot, pbox)) return true;
            connBox.x.low = std::min<DBU>(connBox.x.low, pbox.x.low);
            connBox.x.high = std::max<DBU>(connBox.x.high, pbox.x.high);
            connBox.y.low = std::min<DBU>(connBox.y.low, pbox.y.low);
            connBox.y.high = std::max<DBU>(connBox.y.high, pbox.y.high);
        }
    }

    if (shiftPinVia(tap, type, closeIdx)) {
        db::routeStat.increment(db::RouteStage::POST, db::MiscRouteEvent::MOVE_PIN_VIA, 1);
        return true;
    }

    if (!connBox.IsValid()) return false;
    patchPinVia(connBox, viabot, pinViaPatches[tap.pinIdx], dbNet);
    int forbidSide = 2;
    if (connBox.x.low < viabot.x.low && connBox.x.high > viabot.x.high) {
        // coveredSide = connBox.y.Contain(viabot.y.low) ? 1 : 0;
        forbidSide = 1;
    }
    else if (connBox.y.low <viabot.y.low && connBox.y.high > viabot.y.high) {
        // coveredSide = connBox.x.Contain(viabot.x.low) ? 2 : 3;
        forbidSide = 0;
    }
    if (forbidSide < 2) {
        for (size_t i : closeIdx) {
            const auto &patchIts = pboxes[i].IntersectWith(viabot);
            if (connBox.HasIntersectWith(pboxes[i]) ||
                patchIts[forbidSide].IsValid() ||
                !patchIts[1-forbidSide].IsValid())
                continue;
            patchPinVia(pboxes[i], viabot, pinViaPatches[tap.pinIdx], dbNet);
        }
    }
    return true;
}

bool PostRoute::safeShift(const db::GridSteiner &tap) {
    const db::ViaType *newtype = nullptr;
    int dir = 1 - database.getLayerDir(1);
    const auto &vialoc = database.getLoc(tap);
    const auto &pboxes = dbNet.pinAccessBoxes[tap.pinIdx];
    const auto &cutlayer = database.getCutLayer(0);
    DBU disp = 0, mindisp = 500;
    // // debug
    // if (dbNet.getName() == "net9187" && vialoc.x > 1335100)
    //     printf("    enter safeShift\n");
    for (const auto &pbox : pboxes) {
        for (const auto &type : cutlayer.allViaTypes) {
            if (type.top[dir].range() <= type.top[1-dir].range() ||
                type.bot[dir].range() > pbox[dir].range())
                continue;
            const auto &viabot = type.getShiftedBotMetal(vialoc);
            if (pbox[1-dir].low > viabot[1-dir].low ||
                pbox[1-dir].high < viabot[1-dir].high)
                continue;
            if (viabot[dir].low < pbox[dir].low)
                disp = pbox[dir].low - viabot[dir].low;
            else if (viabot[dir].high > pbox[dir].high)
                disp = pbox[dir].high - viabot[dir].high;
            if (abs(disp) < abs(mindisp)) {
                mindisp = disp;
                newtype = &type;
            }
        }
    }
    // // debug
    // if (dbNet.getName() == "net9187" && vialoc.x > 1335100)
    //     printf("    new type == nullptr : %d;   mindisp : %ld\n", int(newtype==nullptr), mindisp);
    if (newtype == nullptr || mindisp == 0) return false;

    auto newloc = vialoc;
    newloc[dir] += mindisp;
    const auto &ug = database.getUpper(tap);
    const auto &newtop = newtype->getShiftedTopMetal(newloc);
    // // debug
    // if (dbNet.getName() == "net9187" && vialoc.x > 1335100) {
    //     printf("    newtop : ");
    //     std::cout << newtop << std::endl;
    //     database.debugHasVioRoutedMetalOnTrack(dbNet.idx, ug.trackIdx, newtop[dir].low, newtop[dir].high);
    // }
    if (database.hasVioRoutedMetalOnTrack(dbNet.idx, 1, ug.trackIdx, newtop[dir].low, newtop[dir].high))
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
    linkToPins[nei1].emplace_back(vialoc, newloc);
    linkViaToPins[nei1] = std::make_pair(0, newloc);
    linkViaTypes[nei1] = newtype;
    return true;
}

bool PostRoute::shiftPinVia(const db::GridSteiner &tap, const db::ViaType *type, const vector<size_t> &cidx) {
    const auto &vialoc = database.getLoc(tap);
    // const auto &viabot = type->getShiftedBotMetal(vialoc);
    const auto &pboxes = dbNet.pinAccessBoxes[tap.pinIdx];
    const auto &layer = database.getLayer(0);
    auto viabot = type->getShiftedBotMetal(vialoc);
    vector<utils::BoxT<DBU>> neiMetals, connBoxes;
    DBU wSpace = layer.paraRunSpaceForLargerWidth,
        // sSpace = layer.getParaRunSpace(0),
        eSpace = layer.maxEolSpace,
        eolwithin = layer.maxEolWithin,
        wWidth = layer.parallelWidth[1];
    int dir = 1 - layer.direction;
    DBU disp = 0, maxdisp = 0;

    // nei metals
    db::BoxOnLayer qBox(0, viabot);
    for (int dir=0; dir<2; dir++) {
        qBox[dir].low -= wSpace;
        qBox[dir].high += wSpace;
    }
    auto qGrid = database.rangeSearch(qBox);
    qGrid.trackRange.low -= 1;
    qGrid.trackRange.high += 1;
    qGrid.crossPointRange.low -= 1;
    qGrid.crossPointRange.high += 1;
    database.getFixedBox(qBox, neiMetals, dbNet.idx);
    database.getRoutedBox(qGrid, neiMetals, dbNet.idx);

    if (neiMetals.empty()) return false;
    else if (safeShift(tap)) return true;

    // check wider width metal
    connBoxes.push_back(viabot);
    for (size_t i : cidx) {
        connBoxes.push_back(pboxes[i]);
        auto its = pboxes[i].IntersectWith(viabot);
        if (its.x.range() > 0 && its.y.range() < 0) {
            std::swap(its.y.low, its.y.high);
            connBoxes.push_back(its);
        }
        else if (its.x.range() < 0 && its.y.range() > 0) {
            std::swap(its.x.low, its.x.high);
            connBoxes.push_back(its);
        }
    }
    vector<DBU> xlocs;
    vector<utils::IntervalT<DBU>> yRanges;
    for (const auto &box : connBoxes) {
        xlocs.push_back(box.x.low);
        xlocs.push_back(box.x.high);
    }
    sort(xlocs.begin(), xlocs.end());
    xlocs.erase(std::unique(xlocs.begin(), xlocs.end()), xlocs.end());
    size_t numx = xlocs.size() - 1;
    yRanges.resize(numx);
    for (const auto &box : connBoxes) {
        for (size_t i=0; i<yRanges.size(); i++) {
            if (xlocs[i] < box.x.low) continue;
            else if (xlocs[i] < box.x.high) {
                yRanges[i].low = min<DBU>(yRanges[i].low, box.y.low);
                yRanges[i].high = max<DBU>(yRanges[i].high, box.y.high);
            }
            else break;
        }
    }
    connBoxes.clear();
    for (size_t i=0; i<numx; i++) {
        if (yRanges[i].range() < wWidth) continue;
        size_t j = i+1;
        while (j < numx) {
            const auto &temp = yRanges[i].IntersectWith(yRanges[j]);
            if (temp.range() < wWidth)
                break;
            yRanges[i] = temp;
            j++;
        }

                // if (dbNet.getName() == "net35533") {
                //     printf(" sSpace : %ld", layer.getParaRunSpace(0));
                //     printf(" vialoc : ");
                //     std::cout << vialoc << std::endl;
                //     for (DBU x : xlocs) printf(" %ld,", x);
                //     printf("\n");
                //     for (const auto &y : yRanges) printf(" %ld,%ld; ", y.low, y.high);
                //     printf("\n");
                //     printf(" %ld, %ld, %ld, %ld\n", xlocs[i], xlocs[j], yRanges[i].low, yRanges[i].high);
                // }

        if (xlocs[j] - xlocs[i] >= wWidth &&
            ((yRanges[i].low <= viabot.y.low && yRanges[i].high >= viabot.y.high) ||
            (xlocs[i] <= viabot.x.low && xlocs[j] >= viabot.x.high))) {
            connBoxes.emplace_back(xlocs[i]-wSpace, viabot.y.low-wSpace, xlocs[j]+wSpace, viabot.y.high+wSpace);
            break;
        }
    }
    dir = viabot.x.range() > viabot.y.range() ? 0 : 1; // length dir
    connBoxes.push_back(viabot);
    auto &backbox = connBoxes.back();
    backbox[dir].low -= eSpace;
    backbox[dir].high += eSpace;
    backbox[1-dir].low -= eolwithin;
    backbox[1-dir].high += eolwithin;

    dir = 1 - database.getLayerDir(1); // metal 2 pref dir
    for (const auto &nbox : neiMetals) {
        for (const auto &qbox : connBoxes) {
            const auto &its = nbox.IntersectWith(qbox);
            if (!its.IsValid() || its[dir].range() <= maxdisp) continue;
            maxdisp = its[dir].range();
            if (nbox[dir].high > viabot[dir].high) disp = maxdisp * -1;
            else disp = maxdisp;
        }
    }

    if (disp == 0 || abs(disp) > 80) {
        return false;
    }

    auto newloc = vialoc;
    newloc[dir] += disp;
    viabot[dir].ShiftBy(disp);
    const auto &ug = database.getUpper(tap);
    const auto &newtop = type->getShiftedTopMetal(newloc);
    if (database.getPinLinkVio({0, viabot}, dbNet.idx, false) > 0 ||
        database.hasVioRoutedMetalOnTrack(dbNet.idx, 1, ug.trackIdx, newtop[dir].low, newtop[dir].high))
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
    linkToPins[nei1].emplace_back(vialoc, newloc);
    linkViaToPins[nei1] = std::make_pair(0, newloc);
    linkViaTypes[nei1] = type;
    // // check connected
    // disp = 1000;
    // for (size_t i : cidx) {
    //     maxdisp = utils::Dist(pboxes[i], viabot);
    //     if (maxdisp < disp) {
    //         disp = maxdisp;
    //         dir = int(i);
    //     }
    // }
    // if (disp > 0) {
    //     viabot.Set(pboxes[dir].x, pboxes[dir].y);
    // }
    return true;
}

bool PostRoute::handleMacroPins(const db::GridSteiner &tap) {
    int lid = tap.layerIdx;
    if (lid == 0 || tap.fakePin) return false;

    std::shared_ptr<db::GridSteiner> neiPtr;
    if (tap.parent == nullptr) neiPtr = tap.children[0];
    else neiPtr = tap.parent;

    utils::BoxT<DBU> legalBox, obsBox;
    for (const auto &pab : dbNet.pinAccessBoxes[tap.pinIdx]) {
        if (pab.layerIdx == lid) {
            legalBox = pab;
            break;
        }
    }
    for (const auto &obs : database.obsBoxes[lid]) {
        if (obs.HasIntersectWith(legalBox)) {
            obsBox = obs;
            break;
        }
    }
    if (!obsBox.IsValid()) return false;

        // printf(" net %d pin %d abs num : %d\n", dbNet.idx, tap.pinIdx, int(dbNet.pinAccessBoxes[tap.pinIdx].size()));

    const auto &taploc = database.getLoc(tap);

    bool fix = false;
    if (neiPtr->layerIdx == lid) {
        DBU w = database.getLayer(lid).width / 2,
            thre = database.getLayer(lid).paraRunSpaceForLargerWidth + w;
        if (legalBox.x.low == obsBox.x.low) legalBox.x.low = -1;
        else if (legalBox.x.high == obsBox.x.high) legalBox.x.high = std::numeric_limits<DBU>::max();
        if (legalBox.y.low == obsBox.y.low) legalBox.y.low = -1;
        else if (legalBox.y.high == obsBox.y.high) legalBox.y.high = std::numeric_limits<DBU>::max();
        legalBox.x.low += w;
        legalBox.x.high -= w;
        legalBox.y.low += w;
        legalBox.y.high -= w;
        if (utils::Dist(obsBox, taploc) < thre && !legalBox.Contain(taploc)) {
            const auto &u = database.getLoc(*neiPtr);
            const auto &t = legalBox.GetNearestPointTo(u);
            const auto &v = obsBox.GetNearestPointTo(t);
            // database.writeDEFWireSegment(dbNet, u, t, lid);
            // database.writeDEFWireSegment(dbNet, t, v, lid);
            // std::cout << dbNet.getName() << ", " << u << ", " << t << ", " << v << "; " << legalBox << "; " << obsBox << std::endl;
            linkToPins[neiPtr].emplace_back(u, t);
            linkToPins[neiPtr].emplace_back(t, v);
            fix = true;
        }
    }
    else {
        neiPtr->pinIdx = tap.pinIdx;
        utils::PointT<DBU> vx(legalBox.cx(), taploc.y), vy(legalBox.cx(), legalBox.cy());
        linkToPins[neiPtr].emplace_back(taploc, vx);
        linkToPins[neiPtr].emplace_back(vx, vy);
        // database.writeDEFWireSegment(dbNet, taploc, vx, p->layerIdx);
        // database.writeDEFWireSegment(dbNet, vx, vy, p->layerIdx);
        linkViaToPins[neiPtr] = std::make_pair(min(lid, neiPtr->layerIdx), vy);
        linkViaTypes[neiPtr] = tap.viaType;
        fix = true;
    }

    if (fix) {
        if (tap.parent == nullptr) {
            neiPtr->parent = nullptr;
            for (size_t i=0; i<dbNet.gridTopo.size(); i++) {
                if (dbNet.gridTopo[i].get() == &tap)
                    dbNet.gridTopo[i] = neiPtr;
            }
        }
        else {
            for (size_t i=0; i<neiPtr->children.size(); i++) {
                if (neiPtr->children[i].get() == &tap) {
                    neiPtr->children[i] = neiPtr->children.back();
                    neiPtr->children.pop_back();
                    break;
                }
            }
        }
        neiPtr->pinIdx = tap.pinIdx;
        neiPtr->fakePin = true;
    }
    return fix;
}