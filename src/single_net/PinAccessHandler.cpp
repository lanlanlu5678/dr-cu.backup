#include "PinAccessHandler.h"

inline bool isEnclosed(const utils::BoxT<DBU> &viaBox, const db::BoxOnLayer &pbox) {
    return !((viaBox.x.low < pbox.x.low) ||
            (viaBox.y.low < pbox.y.low) ||
            (viaBox.x.high > pbox.x.high) ||
            (viaBox.y.high > pbox.y.high));
}

void patch(const utils::BoxT<DBU> &pin, utils::BoxT<DBU> viabot,
                    vector<utils::BoxT<DBU>> &pvp, db::Net &net) {
    DBU margin = 0, width = database.getLayer(0).width;
    int dir = 0;
    const auto &its = viabot.IntersectWith(pin);
    // if (its.x.range() < width && its.y.range() < width) {
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
    // }
}

void PinAccessHandler::run() {
    if (tap->layerIdx > 0) {
        handleMacroPin();
        return;
    }
    const db::ViaType *type = nullptr;
    if (tap->parent) type = tap->viaType;
    else if (!(tap->children.empty())) type = tap->children[0]->viaType;

    if (type != nullptr)
        handlePinVia(type);
}

void PinAccessHandler::handlePinVia(const db::ViaType *type) {
    const auto &vialoc = database.getLoc(*tap);
    const auto &viabot = type->getShiftedBotMetal(vialoc);
    vector<utils::BoxT<DBU>> neiMetals;
    // get nei metals
    DBU widerSpace = database.getLayer(0).paraRunSpaceForLargerWidth;
    db::BoxOnLayer qbox(0, viabot.x.low-widerSpace, viabot.y.low-widerSpace, viabot.x.high+widerSpace, viabot.y.high+widerSpace);
    db::GridBoxOnLayer qgrid(0, {tap->trackIdx-1, tap->trackIdx+1}, {tap->crossPointIdx-1, tap->crossPointIdx+1});
    database.getFixedBox(qbox, neiMetals, net.idx);
    database.getRoutedBox(qgrid, neiMetals, net.idx);
    if (!neiMetals.empty() && (safeShift() || handleVioPinVia(type, neiMetals)))
        return;
    patchPinVia(type);
}

bool PinAccessHandler::handleVioPinVia(const db::ViaType *type, const vector<utils::BoxT<DBU>> &neiMetals) {
    const auto &pabs = net.pinAccessBoxes[tap->pinIdx];
    const auto &vialoc = database.getLoc(*tap);
    const auto &viabot = type->getShiftedBotMetal(vialoc);
    
    const auto &layer = database.getLayer(0);
    DBU minSpace = layer.getParaRunSpace(0),
        powSpace = 0, fixdist0 = 0, fixdist1 = 0;

    vector<DBU> disp(2, 0);
    utils::BoxT<DBU> eolRegion;
    int lendir = viabot.x.range() > viabot.y.range() ? 0 : 1,
        shiftdir = layer.direction;

    if (neiMetals.empty()) return false;

    bool eol = true;
    if (isWider(viabot, minSpace, layer.parallelWidth[1])) {
        powSpace = std::pow(layer.paraRunSpaceForLargerWidth, 2);
        eol = false;
    }
    else {
        powSpace = std::pow(minSpace, 2);
        eolRegion = viabot;
        eolRegion[lendir].low -= layer.maxEolSpace;
        eolRegion[lendir].high += layer.maxEolSpace;
        eolRegion[1-lendir].low -= layer.maxEolWithin;
        eolRegion[1-lendir].high += layer.maxEolWithin;
    }
    for (const auto &box : neiMetals) {
        fixdist0 = std::pow(utils::Dist(viabot[1-shiftdir], box[1-shiftdir]), 2);
        fixdist1 = std::pow(utils::Dist(viabot[shiftdir], box[shiftdir]), 2);
        if (fixdist0 + fixdist1 < powSpace) {
            fixdist1 = std::ceil(std::sqrt(powSpace - fixdist0));
            disp[0] = max<DBU>(disp[0], box[shiftdir].high + fixdist1 - viabot[shiftdir].low);
            disp[1] = min<DBU>(disp[1], box[shiftdir].low - fixdist1 - viabot[shiftdir].high);
        }
        if (eol) {
            const auto &its = box.IntersectWith(eolRegion);
            if (its.IsValid()) {
                disp[0] = max<DBU>(disp[0], box[shiftdir].high - eolRegion[shiftdir].low);
                disp[1] = min<DBU>(disp[1], box[shiftdir].low - eolRegion[shiftdir].high);
            }
        }
    }
        // if (net.getName() == "net28877") {
        //     printf(" eol %d, disp %d, %d; powSpace : %d; num nei %d\n", int(eol), int(disp[0]), int(disp[1]), int(powSpace), int(neiMetals.size()));
        //     std::cout << viabot << std::endl;
        //     std::cout << neiMetals[0] << std::endl;
        // }
    if (abs(disp[0]) > abs(disp[1])) std::swap(disp[0], disp[1]);
    const auto &upper = database.getUpper(*tap);
    for (DBU displacement : disp) {
        if (displacement == 0 || abs(displacement) > 160) continue;
        bool notConn = true;
        auto newloc = vialoc;
        newloc[shiftdir] += displacement;
        const auto &newbot = type->getShiftedBotMetal(newloc);
        for (const auto &pab : pabs) {
            if (pab.HasIntersectWith(newbot)) {
                notConn = false;
                break;
            }
        }
        if (notConn) continue;
        const auto &newtop = type->getShiftedTopMetal(newloc);
        if (!database.hasVioRoutedMetalOnTrack(net.idx, 1, upper.trackIdx, newtop[shiftdir].low, newtop[shiftdir].high)) {
            editGridTopo(vialoc, newloc, type);
            yes = true;
            break;;
        }
    }
    return yes;
}

bool PinAccessHandler::isWider(const utils::BoxT<DBU> &viabot, DBU space, DBU wWidth) {
    vector<utils::BoxT<DBU>> connBoxes {viabot};
    for (const auto &pab : net.pinAccessBoxes[tap->pinIdx]) {
        if (utils::Dist(pab, viabot) >= space) continue;
        connBoxes.push_back(pab);
        auto its = pab.IntersectWith(viabot);
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
        if (xlocs[j] - xlocs[i] >= wWidth &&
            ((yRanges[i].low <= viabot.y.low && yRanges[i].high >= viabot.y.high) ||
            (xlocs[i] <= viabot.x.low && xlocs[j] >= viabot.x.high)))
            return true;
    }
    return false;
}

void PinAccessHandler::patchPinVia(const db::ViaType *type) {
    const auto &pabs = net.pinAccessBoxes[tap->pinIdx];
    const auto &vialoc = database.getLoc(*tap);
    const auto &viabot = type->getShiftedBotMetal(vialoc);
    DBU patchDist = database.getLayer(0).getParaRunSpace(0);
    utils::BoxT<DBU> touchBox;
    vector<size_t> closeIdx;

    for (size_t i=0; i<pabs.size(); i++) {
        const auto &pab = pabs[i];
        if (pab.layerIdx > 0 || utils::Dist(pab, viabot) >= patchDist)
            continue;
        closeIdx.push_back(i);
        if (pab.HasIntersectWith(viabot)) {
            if (isEnclosed(viabot, pab)) {
                yes = true;
                return;
            }
            touchBox = touchBox.UnionWith(pab);
        }
    }
    if (!touchBox.IsValid()) return;
    patch(touchBox, viabot, pr->pinViaPatches[tap->pinIdx], net);
    int forbidSide = 2;
    if (touchBox.x.low < viabot.x.low && touchBox.x.high > viabot.x.high)
        forbidSide = 1;
    else if (touchBox.y.low <viabot.y.low && touchBox.y.high > viabot.y.high)
        forbidSide = 0;
    if (forbidSide < 2) {
        for (size_t i : closeIdx) {
            const auto &patchIts = pabs[i].IntersectWith(viabot);
            if (touchBox.HasIntersectWith(pabs[i]) ||
                patchIts[forbidSide].IsValid() ||
                !patchIts[1-forbidSide].IsValid())
                continue;
            patch(pabs[i], viabot, pr->pinViaPatches[tap->pinIdx], net);
        }
    }
    yes = true;
}

bool PinAccessHandler::safeShift() {
    const db::ViaType *newType = nullptr;
    int dir = 1 - database.getLayerDir(1);
    const auto &vialoc = database.getLoc(*tap);
    const auto &cutlayer = database.getCutLayer(0);
    const auto &pabs = net.pinAccessBoxes[tap->pinIdx];
    DBU disp = 0, mindisp = 500;

    for (const auto &pab : pabs) {
        for (const auto &type : cutlayer.allViaTypes) {
            if (type.top[dir].range() <= type.top[1-dir].range() ||
                type.bot[dir].range() > pab[dir].range())
                continue;
            const auto &viabot = type.getShiftedBotMetal(vialoc);
            if (pab[1-dir].low > viabot[1-dir].low ||
                pab[1-dir].high < viabot[1-dir].high)
                continue;
            if (viabot[dir].low < pab[dir].low)
                disp = pab[dir].low - viabot[dir].low;
            else if (viabot[dir].high > pab[dir].high)
                disp = pab[dir].high - viabot[dir].high;
            if (abs(disp) < abs(mindisp)) {
                mindisp = disp;
                newType = &type;
            }
        }
    }

    if (newType == nullptr || mindisp == 0) return false;

    auto newloc = vialoc;
    newloc[dir] += mindisp;
    const auto &upper = database.getUpper(*tap);
    const auto &newtop = newType->getShiftedTopMetal(newloc);

    if (database.hasVioRoutedMetalOnTrack(net.idx, 1, upper.trackIdx, newtop[dir].low, newtop[dir].high))
        return false;
    
    editGridTopo(vialoc, newloc, newType);
    yes = true;
    return true;
}

void PinAccessHandler::editGridTopo(const utils::PointT<DBU> &vialoc,
                                    const utils::PointT<DBU> &newloc,
                                    const db::ViaType *type) {
    std::shared_ptr<db::GridSteiner> nei1;
    if (tap->parent) {
        nei1 = tap->parent;
        for (auto it=nei1->children.begin(); it!=nei1->children.end(); it++) {
            if ((*it).get() == tap) {
                *it = nei1->children.back();
                nei1->children.pop_back();
                break;
            }
        }
    }
    else {
        nei1 = tap->children[0];
        nei1->parent = nullptr;
        for (auto it=net.gridTopo.begin(); it!=net.gridTopo.end(); it++) {
            if ((*it).get() == tap) {
                *it = nei1;
                break;
            }
        }
    }
    nei1->pinIdx = tap->pinIdx;
    nei1->fakePin = true;

    pr->linkToPins[nei1].emplace_back(vialoc, newloc);
    pr->linkViaToPins[nei1] = std::make_pair(0, newloc);
    pr->linkViaTypes[nei1] = type;
}

void PinAccessHandler::handleMacroPin() {
    if (!net.rsynPins[tap->pinIdx].isMacroBlockPin() || tap->fakePin) return;
    int lid = tap->layerIdx;
    std::shared_ptr<db::GridSteiner> neiPtr;
    if (tap->parent == nullptr) neiPtr = tap->children[0];
    else neiPtr = tap->parent;

    utils::BoxT<DBU> legalBox, obsBox;
    for (const auto &pab : net.pinAccessBoxes[tap->pinIdx]) {
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

    const auto &taploc = database.getLoc(*tap);

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
            pr->linkToPins[neiPtr].emplace_back(u, t);
            pr->linkToPins[neiPtr].emplace_back(t, v);
            yes = true;
        }
    }
    else {
        neiPtr->pinIdx = tap->pinIdx;
        utils::PointT<DBU> vx(legalBox.cx(), taploc.y), vy(legalBox.cx(), legalBox.cy());
        pr->linkToPins[neiPtr].emplace_back(taploc, vx);
        pr->linkToPins[neiPtr].emplace_back(vx, vy);
        pr->linkViaToPins[neiPtr] = std::make_pair(min(lid, neiPtr->layerIdx), vy);
        pr->linkViaTypes[neiPtr] = tap->viaType;
        yes = true;
    }

    if (yes) {
        if (tap->parent == nullptr) {
            neiPtr->parent = nullptr;
            for (size_t i=0; i<net.gridTopo.size(); i++) {
                if (net.gridTopo[i].get() == tap)
                    net.gridTopo[i] = neiPtr;
            }
        }
        else {
            for (size_t i=0; i<neiPtr->children.size(); i++) {
                if (neiPtr->children[i].get() == tap) {
                    neiPtr->children[i] = neiPtr->children.back();
                    neiPtr->children.pop_back();
                    break;
                }
            }
        }
        neiPtr->pinIdx = tap->pinIdx;
        neiPtr->fakePin = true;
    }
    return;
}

void PinAccessHandler::handlePinWire() {
    // const db::ViaType *type;
    // std::shared_ptr<db::GridSteiner> u = tap->parent;
    // if (u) {
    //     while (u->parent != nullptr && u->layerIdx == 0 && u->children.size() == 1)
    //         u = u->parent;
    // }
    // else {
    //     u = tap->children[0];
    //     while (u->layerIdx == 0 && u->children.size() == 1)
    //         u = u->children[0];
    // }
    // if (u->viaType) type = u->viaType;
    // else type = u->children[0]->viaType;
    // if (u->layerIdx != 1 || database.getViaVioCost(database.getLower(*u), net.idx, true, type) == 0)
    //     return;

    // DBU margin = database.getLayer(0).pitch * 3;
    // int dir = database.getLayerDir(0);
    // const auto &vialoc = database.getLoc(*u);
    // const auto &cutLayer = database.getCutLayer(0);
    // vector<std::pair<const db::ViaType *, utils::PointT<DBU>>> candidates;
    // for (auto &box : net.pinAccessBoxes[tap->pinIdx]) {
    //     if (utils::Dist(box, vialoc) > margin) continue;
    //     for (const auto &type : cutLayer.allViaTypes) {
    //         if (box.x.range() >= type.bot.x.range() &&
    //             box.y.range() >= type.bot.y.range()) {
    //             auto p = box.GetNearestPointTo(vialoc);
    //             if (p.x != vialoc.x && p.y != vialoc.y) continue;
    //             if (p.x == box.x.low) p.x += type.bot.x.high;
    //             else p.x += type.bot.x.low;
    //             if (p.y == box.y.low) p.y += type.bot.y.high;
    //             else p.y += type.bot.y.low;
    //             candidates.emplace_back(&type, p);
    //         }
    //     }
    // }
    // for (const auto &pair : candidates) {
    //     const auto &newtop = pair.first->getShiftedTopMetal(pair.second);
    //     if (database.hasVioGridlessRoutedMetal(net.idx, {1, newtop}) ||
    //         database.hasVioGridlessRoutedMetal(net.idx, {1, database.getWireBox(1, vialoc, pair.second)}))
    //         continue;
    //     if (tap->parent) {
    //         for (auto it=u->children.begin(); it!=u->children.end(); it++) {
    //             if ((*it)->viaType != nullptr) {
    //                 *it = u->children.back();
    //                 u->children.pop_back();
    //                 break;
    //             }
    //         }
    //     }
    //     else {
    //         u->parent = nullptr;
    //         for (auto it=net.gridTopo.begin(); it!=net.gridTopo.end(); it++) {
    //             if ((*it).get() == tap) {
    //                 *it = u;
    //                 break;
    //             }
    //         }
    //     }
    //     u->pinIdx = tap->pinIdx;
    //     u->fakePin = true;
    //     pr->linkToPins[u].emplace_back(vialoc, pair.second);
    //     pr->linkViaToPins[u] = std::make_pair(0, pair.second);
    //     pr->linkViaTypes[u] = pair.first;
    //     yes = true;
    //     return;
    // }
}