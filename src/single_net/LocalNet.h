#pragma once

#include "RouteGuideGraph.h"

class LocalNet : public RouteGuideGraph, public db::NetBase {
public:
    LocalNet(const db::Net& databaseNet) : dbNet(databaseNet), db::NetBase(databaseNet) {}

    const db::Net& dbNet;

    // note: routeGuides & pinAccessBoxes (inherited from db::NetBase) may be different from db::Net
    // (pinIdx, accessBoxIdx) -> (Grid)BoxOnLayer
    vector<vector<db::GridBoxOnLayer>> gridPinAccessBoxes;
    // guideIdx -> dbGuideIdxes
    vector<vector<int>> dbRouteGuideIdxes;

    int estimatedNumOfVertices = 0;

    void initGridBoxes();
    void initNumOfVertices();
    bool checkPin() const;

    int getViaPenalty(int guideIdx1, int trackIdx1, int cpIdx1, int guideIdx2, int trackIdx2, int cpIdx2) const;
    double getWireSegmentPenalty(int guideIdx, int trackIdx, int cpIdx1, int cpIdx2) const;
    int getCrossPointPenalty(int guideIdx, int trackIdx, int cpIdx) const;

    void print() const;  // db::NetBase::print() + RouteGuideGraph::print()
    void printDebug() const;
    void printPnetDetail() const;

    // PARTIAL RIPUP
    int pnetIdx = -1, repeat = 0;
    // DBU maxLength;
    vector<std::shared_ptr<db::GridSteiner>> pnetPins;
    vector<db::BoxOnLayer> debugBox;
    void initPNetPA();

    void initPNet(int numPitch);
    void traverse(std::shared_ptr<db::GridSteiner> node);
    void quickRG();
    void guideRipup();
    void fixLongEdge(const db::GridSteiner &node, vector<db::BoxOnLayer> &picked);
    void getGroups(const db::GridSteiner &node, std::set<int> &gids, vector<db::BoxOnLayer> &picked);
    int getGroup(const utils::PointT<DBU> &p);

    void initMaxLength();

private:
    void getRouteGuideMapping();
};