#pragma once

#include <fstream>
#include "db/Database.h"

class PartialRipup {
    public:
    static void mergeRouteGuides(db::Net &net);

    static void extractPNets(vector<int> &vioNets);
    static void markLocalRipup(db::Net &net);
    static void markAdaptiveRipup(db::Net &net);

    static void merge(std::shared_ptr<db::GridSteiner> node, int thre);
    static void removeDBEdges(std::shared_ptr<db::GridSteiner> node, int idx);

    static void removeSmallLayerSwitch(db::Net &net);
    static void removeCorners(std::shared_ptr<db::GridSteiner> node, int idx);
    static void shiftPinGrid(db::Net &net);
    static void handlePinSplitVias(db::Net &net);
};