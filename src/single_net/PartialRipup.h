#pragma once

#include "db/Database.h"

class PartialRipup {
    public:

    // static bool inDiffGuides(db::GridSteiner &node1, db::GridSteiner &node2, int currGuideIdx, vector<BoxOnLayer> &guides);
    // static void checkRouteGuides(vector<rgNodes> &rgGraph, db::Net &dbNet);
    static int contains(std::shared_ptr<db::GridSteiner> &node, vector<db::BoxOnLayer> &guides);

    static void pureMark(std::shared_ptr<db::GridSteiner> &node, vector<int> &pnets);
    static void mark(std::shared_ptr<db::GridSteiner> &node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num);
    static void purge(std::shared_ptr<db::GridSteiner> &node, vector<std::shared_ptr<db::GridSteiner>> &pins, int netIdx);

    // static void clearOuterGuide(db::Net &dbNet);
    
    static void extractPseudoNets(db::Net &dbnet, int &num, int &ignore);
};

struct rgNodes {
    int parent=0, breakPinIdx=0;
    vector<int> children;
    std::shared_ptr<db::GridSteiner> enterPoint = nullptr;
    vector<std::shared_ptr<db::GridSteiner>> leavePoints;

    rgNodes() {}
    rgNodes(int p, int b) : parent(p), breakPinIdx(b) {}
};
