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
            if (pboxes[i].HasIntersectWith(connBox) ||
                pboxes[i][forbidSide].HasIntersectWith(viabot[forbidSide]))
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
    vector<utils::BoxT<DBU>> neiMetals;
    DBU wSpace = layer.paraRunSpaceForLargerWidth,
        sSpace = layer.getParaRunSpace(0),
        // eSpace = layer.maxEolSpace,
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

    // if (neiMetals.empty()) return false;

    if (neiMetals.size() > 0 && safeShift(tap)) return true;

    bool wider = false;
    vector<utils::IntervalT<DBU>> touched(4);   // yhi(up), ylo(down), xlo(left), xhi(right)
    for (size_t i : cidx) {
        const auto &pbox = pboxes[i];
        const auto &its = pbox.IntersectWith(viabot);
        // if (its.x.range() > 0) {
        //     if (its.x.range() >= wWidth &&
        //         (max<DBU>(pbox.y.high, viabot.y.high)-min<DBU>(pbox.y.low, viabot.y.low)) >= wWidth)
        //         wider = true;
        //     if (pbox.x.low < viabot.x.low) touched[2].UnionWith(its.x);
        //     else touched[3].UnionWith(its.x);
        // }
        // if (its.y.range() > 0) {
        //     if (its.y.range() >= wWidth &&
        //         (max<DBU>(pbox.x.high, viabot.x.high)-min<DBU>(pbox.x.low, viabot.x.low)) >= wWidth)
        //         wider = true;
        //     if (pbox.y.low < viabot.y.low) touched[1].UnionWith(its.y);
        //     else touched[0].UnionWith(its.y);
        // }
        if (its.x.range() > 0) {
            if (its.x.range() >= wWidth &&
                (max<DBU>(pbox.y.high, viabot.y.high)-min<DBU>(pbox.y.low, viabot.y.low)) >= wWidth)
                wider = true;
            if (pbox.y.low <= viabot.y.low) touched[1] = touched[1].UnionWith(its.x);
            if (pbox.y.high >= viabot.y.high) touched[0] = touched[0].UnionWith(its.x);
        }
        if (its.y.range() > 0) {
            if (its.y.range() >= wWidth &&
                (max<DBU>(pbox.x.high, viabot.x.high)-min<DBU>(pbox.x.low, viabot.x.low)) >= wWidth)
                wider = true;
            if (pbox.x.low <= viabot.x.low) touched[2] = touched[2].UnionWith(its.y);
            if (pbox.x.high >= viabot.x.high) touched[3] = touched[3].UnionWith(its.y);
        }
    }
    int coveredSides = 0;
    for (size_t i=0; i<2; i++) {
        if (touched[i] == viabot.x) coveredSides++;
    }
    for (size_t i=2; i<4; i++) {
        if (touched[i] == viabot.y) coveredSides++;
    }

    qBox.Set(0, viabot);
    dir = viabot.x.range() > viabot.y.range() ? 0 : 1; // length dir
    if (wider || coveredSides > 1) {
        qBox[1-dir].high += wSpace;
        qBox[1-dir].low -= wSpace;
        qBox[dir].high += wSpace;
        qBox[dir].low -= wSpace;
        // qBox[dir].high += max<DBU>(wSpace, eSpace);
        // qBox[dir].low -= max<DBU>(wSpace, eSpace);
    }
    else {
        qBox[1-dir].high += sSpace;
        qBox[1-dir].low -= sSpace;
        qBox[dir].high += sSpace;
        qBox[dir].low -= sSpace;
    }

    dir = 1 - database.getLayerDir(1); // metal 2 pref dir
    for (const auto &nbox : neiMetals) {
        const auto &its = nbox.IntersectWith(qBox);
        if (!its.IsValid() || its[dir].range() <= maxdisp) continue;
        maxdisp = its[dir].range();
        if (nbox[dir].high > viabot[dir].high) disp = maxdisp * -1;
        else disp = maxdisp;
    }

    if (disp == 0 || abs(disp) > 60) {
        // if (coveredSides > 0) {
        //     for (size_t i : cidx) {
        //         const auto &pbox = pboxes[i];
        //         const auto &its = pbox.IntersectWith(viabot);
        //         if (its.IsValid()) continue;
        //         else if (its[dir].IsValid() && its[dir].range() <= 10) {
        //             if (pbox[dir].high > viabot[dir].high) disp = -60;
        //             else disp = 60;
        //         }
        //         else if (its[1-dir].IsValid() && its[1-dir].range() <= 10) {
        //             disp = viabot[dir].range() > viabot[1-dir].range() ? eSpace : sSpace;
        //             disp += its[dir].range();
        //             if (pbox[dir].low > viabot[dir].high) disp = disp * -1;
        //         }
        //     }
        // }
        // else
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
    return true;
}