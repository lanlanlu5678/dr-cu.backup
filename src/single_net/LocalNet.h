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

    // PARTIAL RIPUP
    int pnetIdx = -1;
    vector<std::shared_ptr<db::GridSteiner>> pnetPins;
    void traverse(std::shared_ptr<db::GridSteiner> node);
    void initPNetPins();
    void creatLocalRouteGuides();
    void creatAdaptiveRouteGuides();
    void addDiffLayerRouteGuides();
    void getGridBoxes();
    void replaceMacroPabs(vector<vector<db::GridBoxOnLayer>> &gboxes);
    // // debug
    // DBU ripupArea = 0;

private:
    void getRouteGuideMapping();
};