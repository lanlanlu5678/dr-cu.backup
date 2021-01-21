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

void RouteGrid::getRoutedViaBox(const std::multimap<int, int> &viamap,
                                int lc,
                                int uc,
                                int netIdx,
                                GridPoint &via,
                                vector<utils::BoxT<DBU>> &neiMetals) const {
    auto itLo = viamap.lower_bound(lc);
    auto itHi = viamap.upper_bound(uc);
    for (auto it=itLo; it!=itHi; it++) {
        if (it->second == netIdx) continue;
        via.crossPointIdx = it->first;
        const auto &viaLoc = getLoc(via);
        const auto viaType = getViaType(via);
        neiMetals.push_back(viaType->getShiftedTopMetal(viaLoc));
    }
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

vector<utils::BoxT<DBU>> RouteGrid::getRoutedBox(const GridPoint &tap, int netIdx, bool debug) const {
    GridBoxOnLayer queryGrid(tap.layerIdx, {tap.trackIdx-1, tap.trackIdx+1},
                            {tap.crossPointIdx-1, tap.crossPointIdx+1});
    int lid = tap.layerIdx;
    const auto &layer = layers[lid];
    vector<utils::BoxT<DBU>> neiMetals;
    auto dir = layer.direction;
    utils::PointT<DBU> pu, pv;
    int lc = 0, uc = 0;
    GridPoint via;
    int low = lid - 1;  // lid shoud be 1
    via.layerIdx = low;
    const auto &lowerGrid = getLower(queryGrid);
    lc = max(0, lowerGrid.crossPointRange.low - 1),
    uc = min(layers[low].numCrossPoints() - 1, lowerGrid.crossPointRange.high + 1);
    for (int t=lowerGrid.trackRange.low-1; t<=lowerGrid.trackRange.high; t++) {
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
    lc = max(0, queryGrid.crossPointRange.low - 1),
    uc = min(layer.numCrossPoints() - 1, queryGrid.crossPointRange.high + 1);
    for (int t=queryGrid.trackRange.low-1; t<=queryGrid.trackRange.high; t++) {
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
        log() << "    --- target grid : " << tap.layerIdx << "," << tap.trackIdx
                        << "," << tap.crossPointIdx << std::endl;
        log() << "    --- routed metals :" << std::endl;
        for (const auto &box : neiMetals) {
            log() << "    " << box <<  std::endl;
        }
        log() << std::endl;
    }

    return neiMetals;
}


}   // namespace db