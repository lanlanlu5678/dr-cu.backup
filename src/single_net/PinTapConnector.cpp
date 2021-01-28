#include "PinTapConnector.h"

db::RouteStatus PinTapConnector::run() {
    // 0 Try shift via
    const auto &layer = database.getLayer(tap.layerIdx);
    const auto &tapXY = database.getLoc(tap);
    // if (tap.layerIdx == 0 && tapViaType != nullptr) {
    //     if (noNeedExtraLink(tapViaType, tapXY, dbNet.pinAccessBoxes[pinIdx]))
    //         return db::RouteStatus::SUCC_NORMAL;
    // }

    // 1 Get bestBox
    db::BoxOnLayer bestBox;
    db::RouteStatus status = getBestPinAccessBox(tapXY, tap.layerIdx, dbNet.pinAccessBoxes[pinIdx], bestBox);
    if (status != +db::RouteStatus::SUCC_CONN_EXT_PIN) {
        return status;
    }

    // 2 Get tapsOnPin
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

// bool PinTapConnector::noNeedExtraLink(const db::ViaType *type,
//                                         const utils::PointT<DBU> &tapLoc,
//                                         const vector<db::BoxOnLayer> &pinABs) {
//     const auto &viaBox = type->getShiftedBotMetal(tapLoc);
//     int numVio = 0, dir = 0;
//     utils::BoxT<DBU> connBox;
//     vector<size_t> connIdx;
//     DBU temp = 0, width = database.getLayer(0).width, space = database.getLayer(0).getParaRunSpace(0);

//     for (size_t i=0; i<pinABs.size(); i++) {
//         const auto &pbox = pinABs[i];
//         if (pbox.layerIdx > 0) continue;
//         if (utils::Dist(viaBox, pbox) < space) {
//             connIdx.push_back(i);
//             if (viaBox.HasIntersectWith(pbox)) {
//                 connBox.x.low = std::min<DBU>(connBox.x.low, pbox.x.low);
//                 connBox.x.high = std::max<DBU>(connBox.x.high, pbox.x.high);
//                 connBox.y.low = std::min<DBU>(connBox.y.low, pbox.y.low);
//                 connBox.y.high = std::max<DBU>(connBox.y.high, pbox.y.high);
//             }
//         }
//     }

//     // try shift pin via
//     if (connIdx.empty()) return false;
//     tryShift(tapLoc, connIdx);
//     if (pinViaShifted) {
//         db::routeStat.increment(db::RouteStage::POST, db::MiscRouteEvent::MOVE_PIN_VIA, 1);
//         return true;
//     }

//     if (!connBox.IsValid()) return false;
//     const auto &its = viaBox.IntersectWith(connBox);
//     if (its.x.range() < width && its.y.range() < width) {
//         if (viaBox.x.range() > viaBox.y.range()) {
//             temp = its.x.range();
//         }
//         else {
//             dir = 1;
//             temp = its.y.range();
//         }
//         temp = width - temp;
//         if (temp > 0) {
//             connBox = viaBox;
//             if (its[dir].low == viaBox[dir].low)
//                 connBox[dir].low -= temp;
//             else
//                 connBox[dir].high += temp;
//             database.writeDEFFillRect(database.nets[dbNet.idx], connBox, 0);
//         }
//     }
//     return true;
// }

// void PinTapConnector::tryShift(const utils::PointT<DBU> &tapLoc, const vector<size_t> &connIdx) {
//     std::shared_ptr<db::GridSteiner> nei1;
//     const auto &pabs = dbNet.pinAccessBoxes[pinIdx];
//     int dir = 1 - database.getLayerDir(1);
//     std::multimap<int, size_t> canPior;
//     vector<utils::PointT<DBU>> canLocs;
//     vector<const db::ViaType *> canTypes;
//     DBU hwidth = database.getLayer(1).width * 0.5;
//     // check branches
//     if (tapPtr->parent) {
//         nei1 = tapPtr->parent;
//     }
//     else {
//         nei1 = tapPtr->children[0];
//     }
    
//     // get candidates
//     for (size_t id : connIdx) {
//         const auto &pbox = pabs[id];
//         bestType = database.getBestViaTypeForShift(pbox);
//         if (bestType == nullptr) continue;
//         const auto &bot = bestType->getShiftedBotMetal(tapLoc);
//         auto cand = tapLoc;
//         for (int tdir=0; tdir<2; tdir++) {
//             if (bot[tdir].low < pbox[tdir].low)
//                 cand[tdir] += (pbox[tdir].low - bot[tdir].low);
//             else if (bot[tdir].high > pbox[tdir].high)
//                 cand[tdir] += (pbox[tdir].high - bot[tdir].high);
//         }
//         int pior = 0;
//         if (cand[dir] != tapLoc[dir]) pior += 1;
//         if (cand[1-dir] != tapLoc[dir]) pior += 2;
//         if (pior) {
//             canLocs.push_back(cand);
//             canTypes.push_back(bestType);
//             canPior.insert({pior, canLocs.size()-1});
//         }
//     }
//     if (canLocs.empty()) return;
    
//     // get surrounding routed box
//     const auto &routedBoxes = database.getRoutedBox(*nei1, dbNet.idx);

//     // try candidates
//     utils::BoxT<DBU> oriMetal(tapLoc.x-hwidth, tapLoc.y-hwidth, tapLoc.x+hwidth, tapLoc.y+hwidth);
//     for (const auto &pair : canPior) {
//         size_t i = pair.second;
//         const auto &newViaUpperBox = canTypes[i]->getShiftedTopMetal(canLocs[i]);
//         db::BoxOnLayer newUpperBox(1, newViaUpperBox.UnionWith(oriMetal));
//         if (database.getFixedMetalVio(newUpperBox, dbNet.idx) > 0) continue;
//         const auto &regions = database.getAccurateMetalRectForbidRegions(newUpperBox);
//         if (database.countOvlp(newUpperBox, regions, routedBoxes) == 0) {
//             // update grid topo
//             if (tapPtr->parent) {
//                 for (auto it=nei1->children.begin(); it!=nei1->children.end(); it++) {
//                     if (*it == tapPtr) {
//                         *it = nei1->children.back();
//                         nei1->children.pop_back();
//                         break;
//                     }
//                 }
//             }
//             else
//                 nei1->parent = nullptr;
//             nei1->pinIdx = pinIdx;
//             nei1->fakePin = true;
//             tapPtr = nei1;
//             // get extra link
//             auto turn = tapLoc;
//             turn[dir] = canLocs[i][dir];
//             bestLink = getLinkFromPts({tapLoc, turn, canLocs[i]});
//             bestLinkVia.first = 0;
//             bestLinkVia.second = canLocs[i];
//             bestType = canTypes[i];
//             pinViaShifted = true;
//             break;
//         }
//     }
// }