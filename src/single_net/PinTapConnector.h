#pragma once

#include "db/Database.h"

class PinTapConnector {
public:
    vector<utils::SegmentT<DBU>> bestLink;
    std::pair<int, utils::PointT<DBU>> bestLinkVia = {-1, {}};

    int bestVio = 0;
    bool pinViaMoved = false;

    PinTapConnector(const db::GridPoint& pinTap, const db::Net& databaseNet, int pinIndex)
        : tap(pinTap), dbNet(databaseNet), pinIdx(pinIndex) {}

    db::RouteStatus run(const db::ViaType *tapViaType);
    bool noNeedExtraLink(const db::ViaType *type, const utils::PointT<DBU> &tapLoc, const vector<db::BoxOnLayer> &pinABs);
    void legalPinVia(const db::ViaType *type, const utils::PointT<DBU> &tapLoc, const vector<size_t> &closeBoxIdx);

    static db::RouteStatus getBestPinAccessBox(const utils::PointT<DBU>& tapXY,
                                               int layerIdx,
                                               const std::vector<db::BoxOnLayer>& pinAccessBoxes,
                                               db::BoxOnLayer& bestBox);

    static utils::BoxT<DBU> getLinkMetal(const utils::SegmentT<DBU>& link, int layerIdx);

private:
    const db::GridPoint& tap;
    const db::Net& dbNet;
    int pinIdx;

    void shrinkInterval(utils::IntervalT<DBU>& interval, DBU margin);
    void shrinkBox(db::BoxOnLayer& box, DBU margin);
    vector<utils::SegmentT<DBU>> getLinkFromPts(const vector<utils::PointT<DBU>>& linkPts);
    int getLinkPinSpaceVio(const vector<utils::SegmentT<DBU>>& link, int layerIdx);
};