#pragma once

#include "SingleNetRouter.h"

class SeqRouter {
public:
    int numV = -1;
    vector<int> routerIds;
    vector<utils::BoxT<DBU>> bboxes;

    SeqRouter(int id, int estimatedV) {
        numV = estimatedV;
        routerIds.push_back(id);
        bboxes.resize(database.getLayerNum());
    }

};