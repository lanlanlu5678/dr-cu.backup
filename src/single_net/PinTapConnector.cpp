#include "PinTapConnector.h"

db::RouteStatus PinTapConnector::run(const db::ViaType *tapViaType) {
    // 1 Get bestBox
    auto tapXY = database.getLoc(tap);
    db::BoxOnLayer bestBox;
    db::RouteStatus status = getBestPinAccessBox(tapXY, tap.layerIdx, dbNet.pinAccessBoxes[pinIdx], bestBox);
    if (status != +db::RouteStatus::SUCC_CONN_EXT_PIN) {
        return status;
    }

    // 2 Get tapsOnPin
    const auto &layer = database.getLayer(tap.layerIdx);
    if (tap.layerIdx == 0 && tapViaType != nullptr) {
        if (noNeedExtraLink(tapViaType, tapXY, dbNet.pinAccessBoxes[pinIdx]))
            return db::RouteStatus::SUCC_NORMAL;
    }
    std::vector<utils::PointT<DBU>> tapsOnPin;
    if (bestBox.layerIdx != tap.layerIdx) {
        tapsOnPin.emplace_back(bestBox.cx(), bestBox.cy());
        // also need a via
        bestLinkVia.first = min(bestBox.layerIdx, tap.layerIdx);
        bestLinkVia.second = tapsOnPin[0];
    } else {
        // // a tap within a pin access box
        // if (min(bestBox.width(), bestBox.height()) > layer.widthForSuffOvlp) {
        //     auto bestBoxTmp = bestBox;
        //     shrinkBox(bestBoxTmp, layer.shrinkForSuffOvlp);
        //     tapsOnPin.push_back(bestBoxTmp.GetNearestPointTo(tapXY));
        // }
        // route "into" a pin access box
        shrinkBox(bestBox, std::max<DBU>(layer.shrinkForSuffOvlp, layer.width * 0.5));
        tapsOnPin.push_back(bestBox.GetNearestPointTo(tapXY));
    }

    // 3 Get candidateLinks
    vector<vector<utils::SegmentT<DBU>>> candidateLinks;
    for (auto tapOnPin : tapsOnPin) {
        utils::PointT<DBU> turn1(tapXY.x, tapOnPin.y), turn2(tapOnPin.x, tapXY.y);
        if (layer.direction == Y) {
            // prefer routing from tap on track first (to avoid uncontrolled violations with vias)
            std::swap(turn1, turn2);
        }
        candidateLinks.push_back(getLinkFromPts({tapXY, turn1, tapOnPin}));
        candidateLinks.push_back(getLinkFromPts({tapXY, turn2, tapOnPin}));
    }

    // 4 Get bestLink
    bestVio = std::numeric_limits<int>::max();
    for (auto &candidateLink : candidateLinks) {
        int vio = getLinkPinSpaceVio(candidateLink, tap.layerIdx);
        if (vio < bestVio) {
            bestLink = move(candidateLink);
            bestVio = vio;
            if (bestVio == 0) break;
        }
    }

    return db::RouteStatus::SUCC_CONN_EXT_PIN;
}

void PinTapConnector::shrinkInterval(utils::IntervalT<DBU> &interval, DBU margin) {
    if (interval.range() < margin * 2) {
        DBU center = interval.center();
        interval.Set(center);
    } else {
        interval.low += margin;
        interval.high -= margin;
    }
}

void PinTapConnector::shrinkBox(db::BoxOnLayer &box, DBU margin) {
    shrinkInterval(box.x, margin);
    shrinkInterval(box.y, margin);
}

vector<utils::SegmentT<DBU>> PinTapConnector::getLinkFromPts(const vector<utils::PointT<DBU>> &linkPts) {
    vector<utils::SegmentT<DBU>> link;
    for (int i = 0; (i + 1) < linkPts.size(); ++i) {
        if (linkPts[i] != linkPts[i + 1]) {
            link.emplace_back(linkPts[i], linkPts[i + 1]);
        }
    }
    return link;
}

int PinTapConnector::getLinkPinSpaceVio(const vector<utils::SegmentT<DBU>> &link, int layerIdx) {
    int numVio = 0;
    for (const auto &linkSeg : link) {
        auto linkMetal = getLinkMetal(linkSeg, layerIdx);
        // numVio += database.getFixedMetalVio({layerIdx, linkMetal}, dbNet.idx);
        numVio += database.getPinLinkVio({layerIdx, linkMetal}, dbNet.idx);
        // if (dbNet.getName() == "net27817" && pinIdx == 1) {
        //     printf("\n      net %d pin %d getSpaceVio\n", dbNet.idx, pinIdx);
        //     int n = database.getPinLinkVio({layerIdx, linkMetal}, dbNet.idx, true);
        //     printf("        vio : %d\n\n", n);
        // }
    }
    return numVio;
}

db::RouteStatus PinTapConnector::getBestPinAccessBox(const utils::PointT<DBU> &tapXY,
                                                     int layerIdx,
                                                     const std::vector<db::BoxOnLayer> &pinAccessBoxes,
                                                     db::BoxOnLayer &bestBox) {
    DBU minDist = std::numeric_limits<DBU>::max();

    // 1 Get bestBox
    // 1.1 try same-layer boxes
    for (auto &box : pinAccessBoxes) {
        if (box.layerIdx != layerIdx) {
            continue;
        }
        DBU dist = utils::Dist(box, tapXY);
        if (dist < minDist) {
            minDist = dist;
            bestBox = box;
        }
        if (dist == 0) {
            DBU halfWidth = database.getLayer(layerIdx).halfWidth;
            if ((tapXY.x >= box.x.low + halfWidth &&
                tapXY.x <= box.x.high - halfWidth) ||
                (tapXY.y >= box.y.low + halfWidth &&
                tapXY.y <= box.y.high - halfWidth)) {
                return db::RouteStatus::SUCC_NORMAL;
            }
        }
    }
    // 1.2 try diff-layer boxes
    if (minDist == std::numeric_limits<DBU>::max()) {
        DBU maxArea = std::numeric_limits<DBU>::min();
        for (auto &box : pinAccessBoxes) {
            if (abs(box.layerIdx - layerIdx) == 1 && maxArea < box.area()) {
                maxArea = box.area();
                bestBox = box;
            }
        }
        if (maxArea == std::numeric_limits<DBU>::min()) {
            return db::RouteStatus::FAIL_CONN_EXT_PIN;
        }
    }

    return db::RouteStatus::SUCC_CONN_EXT_PIN;
}

utils::BoxT<DBU> PinTapConnector::getLinkMetal(const utils::SegmentT<DBU> &link, int layerIdx) {
    utils::BoxT<DBU> box(link.x, link.y);  // copy
    int dir = (box[0].range() == 0) ? 0 : 1;
    if (box[1 - dir].low > box[1 - dir].high) {
        std::swap(box[1 - dir].low, box[1 - dir].high);
    }
    DBU halfWidth = database.getLayer(layerIdx).width / 2;
    for (int d = 0; d < 2; ++d) {
        box[d].low -= halfWidth;
        box[d].high += halfWidth;
    }
    return box;
}

bool PinTapConnector::noNeedExtraLink(const db::ViaType *type,
                                        const utils::PointT<DBU> &tapLoc,
                                        const vector<db::BoxOnLayer> &pinABs) {
    const auto &viaBox = type->getShiftedBotMetal(tapLoc);
    size_t numVio = 0, dir = 0, bsize = pinABs.size();
    utils::BoxT<DBU> connBox;
    vector<size_t> closeBoxIdx;
    DBU temp = 0, width = database.getLayer(0).width, space = database.getLayer(0).getParaRunSpace(0);
    for (size_t i=0; i<bsize; i++) {
        const auto &pbox = pinABs[i];
        if (pbox.layerIdx > 0) continue;
        if (utils::Dist(viaBox, pbox) < space) {
            const auto &its = viaBox.IntersectWith(pbox);
            if (its.IsValid()) {
                connBox.x.low = std::min<DBU>(connBox.x.low, pbox.x.low);
                connBox.x.high = std::max<DBU>(connBox.x.high, pbox.x.high);
                connBox.y.low = std::min<DBU>(connBox.y.low, pbox.y.low);
                connBox.y.high = std::max<DBU>(connBox.y.high, pbox.y.high);
            }
            else if (its.x.IsValid() || its.y.IsValid())
                closeBoxIdx.push_back(i);   // except connected box
        }
    }

    // if (!closeBoxIdx.empty()) {
    //     for (size_t i=0; i<bsize; i++) {
    //         if (viaBox.HasIntersectWith(pinABs[i]))
    //             closeBoxIdx.push_back(i);
    //     }
    //     legalPinVia(type, tapLoc, closeBoxIdx);
    //     if (pinViaMoved) {
    //         bestLinkVia.first = 0;
    //         db::routeStat.increment(db::RouteStage::POST, db::MiscRouteEvent::MOVE_PIN_VIA, 1);
    //         // log() << " *** " << dbNet.getName();
    //         // printf("    net %d move pin %d via ***\n", dbNet.idx, pinIdx);
    //         return true;
    //     }
    // }

    if (!connBox.IsValid()) return false;

    const auto &its = viaBox.IntersectWith(connBox);
    if (its.x.range() < width && its.y.range() < width) {
        if (viaBox.x.range() > viaBox.y.range()) {
            temp = its.x.range();
        }
        else {
            dir = 1;
            temp = its.y.range();
        }
        temp = width - temp;
        if (temp > 0) {
            connBox = viaBox;
            if (its[dir].low == viaBox[dir].low)
                connBox[dir].low -= temp;
            else
                connBox[dir].high += temp;
            database.writeDEFFillRect(database.nets[dbNet.idx], connBox, 0);
        }
    }
    return true;
}

void PinTapConnector::legalPinVia(const db::ViaType *type,
                                    const utils::PointT<DBU> &tapLoc,
                                    const vector<size_t> &closeBoxIdx) {
    const auto &pabs = dbNet.pinAccessBoxes[pinIdx];
    DBU minX = 100, minY = 200, distX = 0, distY = 0;
    const auto &oriBox = type->getShiftedBotMetal(tapLoc);
    vector<std::pair<int, int>> targets(2, {0, -1});
    for (auto i : closeBoxIdx) {
        const auto &pbox = pabs[i];
        distX = abs(tapLoc.x - pbox.cx());
        distY = abs(tapLoc.y - pbox.cy());
        if (distX < minX &&
            (oriBox.x.low < pbox.x.low || oriBox.x.high > pbox.x.high) &&
            oriBox.y.HasIntersectWith(pbox.y)) {
            minX = distX;
            targets[0].second = i;
        }
        if (distY < minY &&
            (oriBox.y.low < pbox.y.low || oriBox.y.high > pbox.y.high) &&
            oriBox.x.HasIntersectWith(pbox.x)) {
            minY = distY;
            targets[1].second = i;
        }
    }

    const auto &routedBox = database.getRoutedBox(database.getUpper(tap), dbNet.idx);
    // if (dbNet.getName() == "net9679" || dbNet.getName() == "net9699") {
    //     log() << std::endl;
    //     log() << " ------ " << dbNet.getName();
    //     printf(" (%d) pin %d legalPinVia ------\n", dbNet.idx, pinIdx);
    //     log() << " pin loc : " << tapLoc << std::endl;
    //     printf(" X dir target box : %d    ", targets[0].second);
    //     if (targets[0].second > -1) std::cout << pabs[targets[0].second];
    //     printf("\n");
    //     printf(" Y dir target box : %d    ", targets[1].second);
    //     if (targets[1].second > -1) std::cout << pabs[targets[1].second];
    //     printf("\n");
    //     printf("  routed box :\n");
    //     for (const auto &rbox : routedBox)
    //         log() << "      " << rbox << std::endl;
    //     printf("\n");
    // }

    if (database.getLayerDir(1) == X) {
        std::swap(targets[0], targets[1]);
        targets[0].first = 1;
    }
    else
        targets[1].first = 1;

    for (const auto &i : targets) {
        if (i.second < 0) continue;
        const auto &pbox = pabs[i.second];
        DBU disPlacement = pbox[i.first].center() - tapLoc[i.first],
            diff = abs(pbox[i.first].range() - oriBox[i.first].range()) * 0.5;
        if (disPlacement > 0) disPlacement -= diff;
        else disPlacement += diff;
        utils::PointT<DBU> viaCan = tapLoc;
        viaCan[i.first] += disPlacement;
        const auto &viaBox = type->getShiftedTopMetal(viaCan);
        db::BoxOnLayer newViaBox(1, viaBox.UnionWith(type->getShiftedTopMetal(tapLoc)));
        const auto &regions = database.getAccurateMetalRectForbidRegions(newViaBox);
        if (database.countOvlp(newViaBox, regions, routedBox) == 0) {
            bestLinkVia.second = viaCan;
            pinViaMoved = true;
            break;
        }
    }
}