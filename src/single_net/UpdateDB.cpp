#include "UpdateDB.h"

void UpdateDB::commitRouteResult(LocalNet &localNet, db::Net &dbNet) {
    // update RouteGrid, PARTIAL RIPUP
    localNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        if (node->parent) {
            database.useEdge({*node, *(node->parent)}, dbNet.idx);
        }
        if (node->extWireSeg) {
            database.useEdge(*(node->extWireSeg), dbNet.idx);
        }
    });
    if (localNet.pnetIdx < 0) {
        // update db::Net
        dbNet.gridTopo = move(localNet.gridTopo);
    }
    else {
        // // debug
        // if (localNet.idx == 179558 && localNet.pnetIdx == 1)
        //     localNet.printDebug();

        localNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
            int pidx = node->pinIdx;
            if (pidx > 0) {
                for (auto c : localNet.pnetPins[pidx]->children) {
                    if (!(c->isVio)) {
                        c->parent = node;
                        node->children.push_back(c);
                        // c = nullptr;
                    }
                }
                node->pinIdx = localNet.pnetPins[pidx]->pinIdx;
            }
        });

        for (auto root : localNet.gridTopo) {
            if (root->pinIdx == 0) {
                auto p = localNet.pnetPins[0]->parent;
                if (p) {
                    root->parent = p;
                    for (int i=0; i<p->children.size(); i++)
                        if (p->children[i] == localNet.pnetPins[0])
                            p->children[i] = root;
                }
                else {
                    for (int i=0; i<dbNet.gridTopo.size(); i++)
                        if (dbNet.gridTopo[i] == localNet.pnetPins[0])
                            dbNet.gridTopo[i] = root;
                }
                for (auto c : localNet.pnetPins[0]->children) {
                    if (!(c->isVio)) {
                        c->parent = root;
                        root->children.push_back(c);
                        c = nullptr;
                    }
                }
                root->pinIdx = localNet.pnetPins[0]->pinIdx;
            }
            else dbNet.gridTopo.push_back(root);
        }

        localNet.gridTopo.clear();
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

// PARTIAL RIPUP
bool UpdateDB::checkViolation(db::Net &dbNet) {
    bool hasVio = false, currVio = false;

    auto checkEdge = [&](const db::GridEdge& edge) {
        if (database.getEdgeVioCost(edge, dbNet.idx, false)) {
            hasVio = true;
            currVio = true;
            if (db::rrrIterSetting.addDiffLayerGuides) {
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
            }
        }
    };

    dbNet.pnets.clear();
    dbNet.vioNodes.clear();
    dbNet.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        currVio = false;
        if (node->parent) checkEdge({*node, *(node->parent)});
        if (node->extWireSeg) checkEdge(*(node->extWireSeg));
        if (currVio && !(db::rrrIterSetting.oriMode)) {
            node->isVio = true;
            dbNet.vioNodes.push_back(node);
        }
        else node->isVio = false;
    });
    return hasVio;
}