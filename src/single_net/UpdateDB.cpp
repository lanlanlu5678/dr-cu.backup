#include "UpdateDB.h"

void UpdateDB::commitRouteResult(LocalNet &localNet, db::Net &dbNet) {

    if (dbNet.pinsOfPseudoNets.empty()) {
        // update db::Net
        dbNet.gridTopo = move(localNet.gridTopo);
        // update RouteGrid
        dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
            if (node->parent) {
                database.useEdge({*node, *(node->parent)}, dbNet.idx);
            }
            if (node->extWireSeg) {
                database.useEdge(*(node->extWireSeg), dbNet.idx);
            }
        });
    }
    else {
        // partial ripup : update dbNet
        // original GridSteiner should be replaced
        auto &pins = dbNet.pinsOfPseudoNets[localNet.pseudoNetIdx];
        int pidx;

        // merge end pins
        localNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node){
            pidx = node->pinIdx;
            if (pidx > 0) {
                for (auto &c : pins[pidx]->children) {
                    c->parent = node;
                    node->children.push_back(c);
                }
                node->pinIdx = pins[pidx]->pinIdx;
            }

            // update RouteGrid
            if (node->parent) {
                database.useEdge({*node, *(node->parent)}, dbNet.idx);
            }
            if (node->extWireSeg) {
                database.useEdge(*(node->extWireSeg), dbNet.idx);
            }
        });

        // merge root
        auto root = localNet.gridTopo[0];
        pidx = pins[0]->pinIdx;
        if (pidx > -1) {
            for (auto &n : dbNet.gridTopo) {
                if (n == pins[0]) n = root;
            }
        }
        if (pins[0]->parent) {
            root->parent = pins[0]->parent;
            for (auto &c : pins[0]->parent->children)
                if (c == pins[0]) c = root;
        }
        for (auto &c : pins[0]->children) {
            if (!c->isVio) root->children.push_back(c);
        }
        root->pinIdx = pidx;
        if (root->extWireSeg) database.useEdge(*(root->extWireSeg), dbNet.idx);
    }
}

void UpdateDB::clearRouteResult(db::Net &dbNet) {
    // update RouteGrid
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (node->parent) {
            database.removeEdge({*node, *(node->parent)}, dbNet.idx);
        }
        if (node->extWireSeg) {
            database.removeEdge(*(node->extWireSeg), dbNet.idx);
        }
    });
    // update db::Net
    dbNet.clearResult();
}

void UpdateDB::commitMinAreaRouteResult(db::Net& dbNet) {
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (node->extWireSeg) {
            database.useEdge(*(node->extWireSeg), dbNet.idx);
        }
    });
};

void UpdateDB::clearMinAreaRouteResult(db::Net& dbNet) {
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (node->extWireSeg) {
            database.removeEdge(*(node->extWireSeg), dbNet.idx);
        }
    });
};

void UpdateDB::commitViaTypes(db::Net& dbNet) {
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (!(node->parent)) return;
        db::GridEdge edge(*node, *(node->parent));
        if (!edge.isVia())  return;
        database.markViaType(edge.lowerGridPoint(), node->viaType);
    });
};

bool UpdateDB::checkViolation(db::Net &dbNet) {
    bool hasVio = false;

    auto checkEdge = [&](const db::GridEdge& edge, bool &isVio) {
        if (database.getEdgeVioCost(edge, dbNet.idx, false)) {
            hasVio = true;
            auto uLoc = database.getLoc(edge.u);
            auto vLoc = database.getLoc(edge.v);
            std::vector<std::pair<boostBox, int>> relatedGuides;
            auto checkBox = [&](const db::BoxOnLayer& boxOnLayer) {
                auto box = boxOnLayer;
                database.expandBox(box, db::rrrIterSetting.defaultGuideExpand);
                boostBox query_box(boostPoint(box.x.low, box.y.low), boostPoint(box.x.high, box.y.high));
                dbNet.routeGuideRTrees[box.layerIdx].query(bgi::intersects(query_box), std::back_inserter(relatedGuides));
            };
            if (edge.u.layerIdx != edge.v.layerIdx) {
                checkBox({edge.u.layerIdx, uLoc, uLoc});
                checkBox({edge.v.layerIdx, vLoc, vLoc});
            }
            else {
                checkBox({edge.u.layerIdx, uLoc, vLoc});
            }
            for (const auto& guide : relatedGuides) {
                ++dbNet.routeGuideVios[guide.second];
            }

            // partial ripup
            isVio = true;
        }
    };
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (node->parent) checkEdge({*node, *(node->parent)}, node->isVio);
        if (node->extWireSeg) checkEdge(*(node->extWireSeg), node->isVio);
        
        // partial ripup
        node->pNetIdx = -1;
        node->distance = 4;
    });
    return hasVio;
}