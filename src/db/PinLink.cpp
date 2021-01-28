#include "RouteGrid.h"

namespace db {

utils::BoxT<DBU> RouteGrid::getWireBox(int layerIdx, utils::PointT<DBU> pu, utils::PointT<DBU> pv) const {
    DBU width = layers[layerIdx].width * 0.5;
    if (pu.x > pv.x || pu.y > pv.y)
        std::swap(pu, pv);
    return utils::BoxT<DBU>(pu.x-width, pu.y-width, pv.x+width, pv.y+width);
}

int RouteGrid::countOvlp(const BoxOnLayer &box,
                            const vector<utils::BoxT<DBU>> &regions,
                            const vector<utils::BoxT<DBU>> &neiMetals) const {
    int numOvlp = 0;
    for (const auto& neiMetal : neiMetals) {
        // getOvlpFixedMetals
        for (auto forbidRegion : regions) {
            const auto &ovlp = forbidRegion.IntersectWith(neiMetal);
            if (ovlp.IsValid()) {
                numOvlp += (ovlp.area() > 0);
            }
        }

        // getOvlpFixedMetalForbidRegions
        const auto &forbidRegions = getAccurateMetalRectForbidRegions({box.layerIdx, neiMetal});
        for (auto forbidRegion : forbidRegions) {
            const auto &ovlp = forbidRegion.IntersectWith(box);
            if (ovlp.IsValid()) {
                numOvlp += (ovlp.area() > 0);
            }
        }

        // getOvlpC2CMetals
        if (!layers[box.layerIdx].isEolDominated(neiMetal)) {
            DBU space = layers[box.layerIdx].getParaRunSpace(neiMetal);
            numOvlp += (utils::L2Dist(box, neiMetal) < space);
        }
    }
    return numOvlp;
}

const ViaType * RouteGrid::getBestViaTypeForShift(const BoxOnLayer &accessBox) const {
    const auto &cutLayer = cutLayers[0];
    // const ViaType *bestType = nullptr;
    for (const auto &type : cutLayer.allViaTypes) {
        if (accessBox.x.range() >= type.bot.x.range() &&
            accessBox.y.range() >= type.bot.y.range()) {
            return &type;
        }
    }
    return nullptr;
}

int RouteGrid::getPinLinkVio(const BoxOnLayer& box, int netIdx, bool debug) const {
    int lid = box.layerIdx;
    const auto &regions = getAccurateMetalRectForbidRegions(box);
    utils::BoxT<DBU> queryBox = box;
    queryBox.x.low -= layers[lid].fixedMetalQueryMargin;
    queryBox.x.high += layers[lid].fixedMetalQueryMargin;
    queryBox.y.low -= layers[lid].fixedMetalQueryMargin;
    queryBox.y.high += layers[lid].fixedMetalQueryMargin;

    for (const auto& region : regions) {
        queryBox = queryBox.UnionWith(region);
    }

    // fixed metal
    boostBox rtreeQueryBox(boostPoint(queryBox.x.low, queryBox.y.low), boostPoint(queryBox.x.high, queryBox.y.high));
    vector<std::pair<boostBox, int>> queryResults;
    fixedMetals[lid].query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    vector<utils::BoxT<DBU>> neiMetals;
    for (const auto& queryResult : queryResults) {
        if (queryResult.second != netIdx) {
            const auto& b = queryResult.first;
            neiMetals.emplace_back(bg::get<bg::min_corner, 0>(b),
                                    bg::get<bg::min_corner, 1>(b),
                                    bg::get<bg::max_corner, 0>(b),
                                    bg::get<bg::max_corner, 1>(b));
        }
    }

    // routed metal (via only)
    GridBoxOnLayer queryGridBox = rangeSearch(BoxOnLayer(lid, queryBox));
    GridPoint via;
    int lc=0, uc=0;
    if (lid > 0) {
        int low = lid - 1;
        via.layerIdx = low;
        const auto &loBox = getLower(queryGridBox);
        lc = max(0, loBox.crossPointRange.low - 1);
        uc = min(layers[low].numCrossPoints() - 1, loBox.crossPointRange.high + 1);
        for (int t=loBox.trackRange.low-1; t<=loBox.trackRange.high; t++) {
            via.trackIdx = t;
            auto itLo = routedViaMap[low][t].lower_bound(lc);
            auto itHi = routedViaMap[low][t].upper_bound(uc);
            for (auto it=itLo; it!=itHi; it++) {
                if (it->second == netIdx) continue;
                via.crossPointIdx = it->first;
                const auto &viaLoc = getLoc(via);
                const auto viaType = getViaType(via);
                neiMetals.push_back(viaType->getShiftedTopMetal(viaLoc));
            }
        }
    }
    lc = max(0, queryGridBox.crossPointRange.low - 1);
    uc = min(layers[lid].numCrossPoints() - 1, queryGridBox.crossPointRange.high + 1);
    via.layerIdx = lid;
    for (int t=queryGridBox.trackRange.low-1; t<=queryGridBox.trackRange.high; t++) {
        via.trackIdx = t;
        auto itLo = routedViaMap[lid][t].lower_bound(lc);
        auto itHi = routedViaMap[lid][t].upper_bound(uc);
        for (auto it=itLo; it!=itHi; it++) {
            if (it->second == netIdx) continue;
            via.crossPointIdx = it->first;
            const auto &viaLoc = getLoc(via);
            const auto viaType = getViaType(via);
            neiMetals.push_back(viaType->getShiftedBotMetal(viaLoc));
        }
    }

    if (debug) {
        log() << " ------ getPinLinkVio ------ " << std::endl;
        log() << "   forbid regions for : " << box << std::endl;
        for (const auto &region : regions)
            log() << "      " << region << std::endl;
        log() << "   neiMetals : " << std::endl;
        for (const auto &nm : neiMetals)
            log() << "      " << nm << std::endl;
    }

    return countOvlp(box, regions, neiMetals);
}

vector<utils::BoxT<DBU>> RouteGrid::getRoutedBox(GridBoxOnLayer &queryGrid, int netIdx, bool debug) const {
    int lid = queryGrid.layerIdx;
    const auto &layer = layers[lid];
    vector<utils::BoxT<DBU>> neiMetals;
    auto dir = layer.direction;
    utils::PointT<DBU> pu, pv;
    int lc = 0, uc = 0;
    GridPoint via;

    queryGrid.trackRange.low = max(0, queryGrid.trackRange.low);
    queryGrid.crossPointRange.low = max(0, queryGrid.crossPointRange.low);
    queryGrid.trackRange.high = min(layers[lid].numTracks()-1, queryGrid.trackRange.high);
    queryGrid.crossPointRange.high = min(layers[lid].numCrossPoints()-1, queryGrid.crossPointRange.high);

    int low = lid - 1;  // lid shoud be 1
    via.layerIdx = low;
    const auto &lowerGrid = getLower(queryGrid);
    lc = max(0, lowerGrid.crossPointRange.low),
    uc = min(layers[low].numCrossPoints() - 1, lowerGrid.crossPointRange.high);
    for (int t=lowerGrid.trackRange.low; t<=lowerGrid.trackRange.high; t++) {
        via.trackIdx = t;
        auto itLo = routedViaMap[low][t].lower_bound(lc);
        auto itHi = routedViaMap[low][t].upper_bound(uc);
        for (auto it=itLo; it!=itHi; it++) {
            if (it->second == netIdx) continue;
            via.crossPointIdx = it->first;
            const auto &viaLoc = getLoc(via);
            const auto viaType = getViaType(via);
            neiMetals.push_back(viaType->getShiftedTopMetal(viaLoc));
        }
    }
    via.layerIdx = lid;
    lc = max(0, queryGrid.crossPointRange.low ),
    uc = min(layer.numCrossPoints() - 1, queryGrid.crossPointRange.high);
    for (int t=queryGrid.trackRange.low; t<=queryGrid.trackRange.high; t++) {
        via.trackIdx = t;
        auto itLo = routedViaMap[lid][t].lower_bound(lc);
        auto itHi = routedViaMap[lid][t].upper_bound(uc);
        for (auto it=itLo; it!=itHi; it++) {
            if (it->second == netIdx) continue;
            via.crossPointIdx = it->first;
            const auto &viaLoc = getLoc(via);
            const auto viaType = getViaType(via);
            neiMetals.push_back(viaType->getShiftedBotMetal(viaLoc));
        }
        auto queryInterval = boost::icl::interval<int>::closed(lc, uc);
        auto intervals = routedWireMap[lid][t].equal_range(queryInterval);
        for (auto it = intervals.first; it != intervals.second; ++it) {
            int usage = it->second.size();
            if (it->second.count(netIdx))
                usage--;
            if (usage > 0) {
                pu[dir] = layer.tracks[t].location;
                pv[dir] = layer.tracks[t].location;
                pu[1-dir] = layer.crossPoints[first(it->first)].location;
                pv[1-dir] = layer.crossPoints[last(it->first)].location;
                neiMetals.push_back(getWireBox(lid, pu, pv));
            }
        }
    }
    if (debug) {
        log() << std::endl;
        log() << " ------- get new via & routed metal vio ------ " << std::endl;
        log() << "    --- target box : " << queryGrid << std::endl;
        log() << "    --- routed metals :" << std::endl;
        for (const auto &box : neiMetals) {
            log() << "    " << box <<  std::endl;
        }
        log() << std::endl;
    }

    return neiMetals;
}

vector<utils::BoxT<DBU>> RouteGrid::getFixedBox(const BoxOnLayer &queryBox, int netIdx) const {
    boostBox rtreeQueryBox(boostPoint(queryBox.x.low, queryBox.y.low),
                            boostPoint(queryBox.x.high, queryBox.y.high));
    vector<std::pair<boostBox, int>> queryResults;
    fixedMetals[queryBox.layerIdx].query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));

    vector<utils::BoxT<DBU>> results;
    for (const auto &pair : queryResults) {
        if (pair.second != netIdx) {
            const auto &b = pair.first;
            results.emplace_back(bg::get<bg::min_corner, 0>(b),
                                    bg::get<bg::min_corner, 1>(b),
                                    bg::get<bg::max_corner, 0>(b),
                                    bg::get<bg::max_corner, 1>(b));
        }
    }
    return results;
}

int RouteGrid::countFixedMetals(const utils::BoxT<DBU> &viaBox, int netIdx) const {
    DBU space = layers[0].paraRunSpaceForLargerWidth;
    boostBox rtreeQueryBox(boostPoint(viaBox.x.low-space, viaBox.y.low-space),
                            boostPoint(viaBox.x.high+space, viaBox.y.high+space));
    vector<std::pair<boostBox, int>> queryResults;
    fixedMetals[0].query(bgi::intersects(rtreeQueryBox), std::back_inserter(queryResults));
    
    int num = 0;
    for (const auto &pair : queryResults) {
        if (pair.second != netIdx) num++;
    }
    return num;
}

bool RouteGrid::hasOtherRoutedWireOnTrack(int netIdx, int layerIdx, int trackIdx, const utils::IntervalT<int> &range) const {
    auto queryItvl = boost::icl::interval<int>::closed(range.low-1, range.high+1);
    auto itvls = routedWireMap[layerIdx][trackIdx].equal_range(queryItvl);
    for (auto it=itvls.first; it!=itvls.second; it++) {
        for (int id : it->second) {
            if (id != netIdx)
                return true;
        }
    }
    return false;
}

bool RouteGrid::hasOtherRoutedViaOnTrack(int netIdx, int layerIdx, int trackIdx, const utils::IntervalT<int> &range) const {
    auto lb = routedViaMap[layerIdx][trackIdx].lower_bound(range.low-1);
    auto ub = routedViaMap[layerIdx][trackIdx].upper_bound(range.high+1);
    for (auto it=lb; it!=ub; it++) {
        if (it->second != netIdx)
            return true;
    }
    lb = routedViaMapUpper[layerIdx][trackIdx].lower_bound(range.low-1);
    ub = routedViaMapUpper[layerIdx][trackIdx].upper_bound(range.high+1);
    for (auto it=lb; it!=ub; it++) {
        if (it->second != netIdx)
            return true;
    }
    return false;
}


}   // namespace db