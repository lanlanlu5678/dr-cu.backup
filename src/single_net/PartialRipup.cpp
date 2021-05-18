#include "PartialRipup.h"

inline void printndoes(std::shared_ptr<db::GridSteiner> node) {
    const auto &np = database.getLoc(*node);
    std::cout << "  " << node->layerIdx << ", " << np << ";  pinIdx : " << node->pinIdx
                << "; vio : " << int(node->isVio) << std::endl;
    for (auto c : node->children) {
        // if (c->isVio)
            printndoes(c);
    }
}

void PartialRipup::mergeRouteGuides(db::Net &net) {
    vector<int> clusters;
    vector<vector<utils::BoxT<DBU>>> guides(database.getLayerNum());
    for (const auto &g : net.routeGuides) {
        guides[g.layerIdx].push_back(g);
    }
    for (size_t l=0; l<database.getLayerNum(); l++) {
        size_t lsize = guides[l].size();
        if (lsize < 2) continue;
        int count = 0;
        clusters.clear();
        clusters.resize(lsize, -1);
        // merge connected clusters
        for (size_t i=0; i<lsize; i++) {
            const auto &iGuide = guides[l][i];
            for (size_t j=0; j<i; j++) {
                const auto &jGuide = guides[l][j];
                if (iGuide.HasIntersectWith(jGuide)) {
                    if (clusters[j] < 0) {
                        if (clusters[i] < 0) {
                            clusters[j] = count;
                            clusters[i] = count++;
                        }
                        else
                            clusters[j] = clusters[i];
                    }
                    else {
                        if (clusters[i] >= 0) {
                            for (size_t k=0; k<j; k++) {
                                if (clusters[k] == clusters[i])
                                    clusters[k] = clusters[j];
                            }
                        }
                        clusters[i] = clusters[j];
                    }
                }
            }
        }
        for (int c=0; c<count; c++) {
            utils::BoxT<DBU> ubox;
            for (size_t i=0; i<lsize; i++) {
                if (clusters[i] == c) {
                    ubox = ubox.UnionWith(guides[l][i]);
                    guides[l][i].Set(1, 1, 0, 0);
                }
            }
            if (ubox.IsValid()) {
                guides[l].push_back(ubox);
            }
        }
    }
    for (int l=0; l<database.getLayerNum(); l++) {
        if (guides[l].empty()) continue;
        for (const auto &guide : guides[l]) {
            if (guide.IsValid()) {
                net.mergedGuides.emplace_back(l, guide);
                int vio = 0;
                for (size_t i=0; i<net.routeGuides.size(); i++) {
                    if (net.routeGuides[i].layerIdx != l) continue;
                    if (guide.HasIntersectWith(net.routeGuides[i]))
                        vio += net.routeGuideVios[i];
                }
                net.mergedVios.push_back(vio);
            }
        }
    }
}

void PartialRipup::extractPNets(vector<int> &vioNets) {
    // if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE)
    //     log() << "Extracting pseudo nets from violated nets " << std::endl;

    if (db::rrrIterSetting.localRipup) {
        runJobsMT(vioNets.size(), [&](int id) {
            auto &net = database.nets[vioNets[id]];
            markLocalRipup(net);
        });
    }
    else if (db::rrrIterSetting.adaptiveRipup) {
        runJobsMT(database.nets.size(), [&](int id) {
            mergeRouteGuides(database.nets[id]);
        });
        runJobsMT(vioNets.size(), [&](int id) {
            auto &net = database.nets[vioNets[id]];
            // mergeRouteGuides(net);
            for (auto &g : net.mergedGuides)
                database.expandBox(g, db::rrrIterSetting.defaultGuideExpand);
            markAdaptiveRipup(net);
        });
    }

    if (!db::rrrIterSetting.greedyRoute) {
        for (int id : vioNets) {
            for (auto proot : database.nets[id].pnets)
                removeDBEdges(proot, id);
        }
    }
}

void PartialRipup::markLocalRipup(db::Net &net) {
    for (auto vio : net.vioNodes) {
        int layerIdx = vio->layerIdx;
        db::GridBoxOnLayer gbox(layerIdx, {vio->trackIdx-7, vio->trackIdx+7},
                                {vio->crossPointIdx-20, vio->crossPointIdx+20});
        if (gbox.trackRange.low < 0) gbox.trackRange.low = 0;
        if (gbox.crossPointRange.low < 0) gbox.crossPointRange.low = 0;
        gbox.trackRange.high = min(database.getTrackEnd(layerIdx), gbox.trackRange.high);
        gbox.crossPointRange.high = min(database.getCrossPointEnd(layerIdx), gbox.crossPointRange.high);
        auto dbuBox = database.getLoc(gbox);

        std::function<void(std::shared_ptr<db::GridSteiner>)> markBranches;
        markBranches = [&](std::shared_ptr<db::GridSteiner> node) {
            for (auto c : node->children) {
                if (c->isVio) continue;
                const auto &cp = database.getLoc(*c);
                if (dbuBox.Contain(cp)) {
                    c->isVio = true;
                    markBranches(c);
                }
                else if (dbuBox.StrictlyContain(database.getLoc(*node))) {
                    // split long edge
                    std::shared_ptr<db::GridSteiner> p = std::make_shared<db::GridSteiner>(db::GridPoint(node->layerIdx,0,0));
                    dbuBox.layerIdx = node->layerIdx;
                    const auto &qBoxOnLayer = database.rangeSearch(dbuBox);
                    if (node->trackIdx == c->trackIdx) {
                        p->trackIdx = node->trackIdx;
                        p->crossPointIdx = qBoxOnLayer.crossPointRange.GetNearestPointTo(c->crossPointIdx);
                    }
                    else {
                        p->trackIdx = qBoxOnLayer.trackRange.GetNearestPointTo(c->trackIdx);
                        p->crossPointIdx = node->crossPointIdx;
                    }
                    for (size_t i=0; i<node->children.size(); i++) {
                        if (node->children[i] == c) node->children[i] = p;
                    }
                    p->parent = node;
                    p->children.push_back(c);
                    c->parent = p;
                    p->isVio = true;
                }
            }
        };

        // downwards
        markBranches(vio);

        // upwards
        std::shared_ptr<db::GridSteiner> up, curr = vio;
        while (curr->parent) {
            up = curr->parent;
            if (up->isVio) break;
            const auto &upp = database.getLoc(*up);
            if (dbuBox.Contain(upp)) {
                up->isVio = true;
                markBranches(up);
                curr = up;
            }
            else if (dbuBox.StrictlyContain(database.getLoc(*curr))) {
                // split long edge
                std::shared_ptr<db::GridSteiner> p = std::make_shared<db::GridSteiner>(db::GridPoint(curr->layerIdx,0,0));
                dbuBox.layerIdx = curr->layerIdx;
                const auto &qBoxOnLayer = database.rangeSearch(dbuBox);
                if (curr->trackIdx == up->trackIdx) {
                    p->trackIdx = curr->trackIdx;
                    p->crossPointIdx = qBoxOnLayer.crossPointRange.GetNearestPointTo(up->crossPointIdx);
                }
                else {
                    p->trackIdx = qBoxOnLayer.trackRange.GetNearestPointTo(up->trackIdx);
                    p->crossPointIdx = curr->crossPointIdx;
                }
                for (size_t i=0; i<up->children.size(); i++) {
                    if (up->children[i] == curr) up->children[i] = p;
                }
                p->parent = up;
                curr->parent = p;
                p->children.push_back(curr);
                p->isVio = true;
            }
            else
                break;
        }
    }

    for (auto root : net.gridTopo)
        mergeByGrids(root, 8);
    net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
        if ((node->isVio) && (!(node->parent) || !(node->parent->isVio)))
            net.pnets.push_back(node);
    });
}

inline void handleOutOfGuide(std::shared_ptr<db::GridSteiner> node,
                                vector<int> &selected,
                                const vector<db::BoxOnLayer> &guides) {
    // upwards
    auto up = node;
    bool out = true;
    size_t gsize = guides.size();
    while (out && up->parent) {
        up = up->parent;
        up->isVio = true;
        const auto &upp = database.getLoc(*up);
        for (size_t i=0; i<gsize; i++) {
            if (guides[i].layerIdx == up->layerIdx &&
                guides[i].Contain(upp)) {
                selected[i]++;
                out = false;
                break;
            }
        }
    }

    // downwards
    vector<db::GridSteiner *> children, grandChildren;
    children.push_back(node.get());
    while (!children.empty()) {
        for (auto p : children) {
            p->isVio = true;
            for (auto c : p->children) {
                out = true;
                const auto &cp = database.getLoc(*c);
                for (size_t i=0; i<gsize; i++) {
                    if (guides[i].layerIdx == c->layerIdx &&
                        guides[i].Contain(cp)) {
                        selected[i]++;
                        out = false;
                        break;
                    }
                }
                if (out)
                    grandChildren.push_back(c.get());
            }
        }
        children = move(grandChildren);
    }
}

inline void dfs(db::GridSteiner *node,
                db::GridSteiner *exit,
                bool preIn,
                const db::BoxOnLayer &g) {
    const auto &np = database.getLoc(*node);
    if (g.Contain(np)) {
        node->isVio = true;
        if (!preIn && exit != nullptr) {
            node = node->parent.get();
            while (node != exit) {
                node->isVio = true;
                node = node->parent.get();
            }
        }
        preIn = true;
    }
    else {
        if (preIn) exit = node;
        preIn = false;
    }

    for (auto c : node->children) {
        dfs(c.get(), exit, preIn, g);
    }
}

void PartialRipup::markAdaptiveRipup(db::Net &net) {
    size_t guideSize = net.mergedGuides.size();
    vector<int> selected(guideSize, 0);
    for (auto vio : net.vioNodes) {
        const auto &vp = database.getLoc(*vio);
        bool out = true;
        for (size_t i=0; i<guideSize; i++) {
            const auto &gbox = net.mergedGuides[i];
            if (gbox.layerIdx == vio->layerIdx &&
                gbox.Contain(vp)) {
                selected[i]++;
                out = false;
                break;
            }
        }
        if (out) {
            handleOutOfGuide(vio, selected, net.mergedGuides);
        }
    }
    // net.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
    //     const auto &np = database.getLoc(*node);
    //     for (size_t i=0; i<guideSize; i++) {
    //         if (selected[i] > 0 && net.mergedGuides[i].Contain(np)) {
    //             node->isVio = true;
    //             break;
    //         }
    //     }
    // });
    for (size_t i=0; i<guideSize; i++) {
        if (selected[i] == 0) continue;
        for (auto root : net.gridTopo) {
            dfs(root.get(), nullptr, false, net.mergedGuides[i]);
        }
    }

    for (auto root : net.gridTopo)
        mergeByGrids(root, 40);
    net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
        if ((node->isVio) && (!(node->parent) || !(node->parent->isVio)))
            net.pnets.push_back(node);
    });
}

void mergePNetsInSubtrees(std::shared_ptr<db::GridSteiner> node, int thre) {
    for (auto c : node->children)
        if (!(c->isVio) && c->distance < thre) {
            c->isVio = true;
            mergePNetsInSubtrees(c, thre);
        }
}

void PartialRipup::mergeByGrids(std::shared_ptr<db::GridSteiner> node, int thre) {
    for (auto c : node->children) mergeByGrids(c, thre);

    int d = 64, d1 = 64, d2 = 64;
    db::GridSteiner *danger = nullptr;
    
    if (node->isVio) {
        mergePNetsInSubtrees(node, thre);
        node->distance = -1;
    }
    else {
        for (auto c : node->children) {
            d = c->distance;

            if (d1 > d) {
                d2 = d1;
                d1 = d;
                danger = c.get();
            }
            else if (d2 > d) {
                d2 = d;
            }
        }

        if (d1 + d2 < thre) {
            mergePNetsInSubtrees(node, thre);
            node->distance = 0;
        }
        else {
            if (danger == nullptr)
                node->distance = d1;
            else
                node->distance = abs(node->crossPointIdx - danger->crossPointIdx) +
                                    abs(node->trackIdx - danger->trackIdx) +
                                    abs(node->layerIdx - danger->layerIdx) + d1;
        }
    }
}

void PartialRipup::mergeByNodes(std::shared_ptr<db::GridSteiner> node, int thre) {
    for (auto c : node->children) mergeByNodes(c, thre);

    int d = 64, d1 = 64, d2 = 64;
    
    if (node->isVio) {
        mergePNetsInSubtrees(node, thre);
        node->distance = -1;
    }
    else {
        for (auto c : node->children) {
            d = c->distance;

            if (d1 > d) {
                d2 = d1;
                d1 = d;
            }
            else if (d2 > d) {
                d2 = d;
            }
        }

        if (d1 + d2 < thre) {
            mergePNetsInSubtrees(node, thre);
            node->distance = 0;
        }
        else {
            node->distance = d1 + 1;
        }
    }
}

void PartialRipup::removeDBEdges(std::shared_ptr<db::GridSteiner> node, int idx) {
    for (auto c : node->children) {
        if (c->isVio) {
            database.removeEdge({*c, *node}, idx);
            if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), idx);
            removeDBEdges(c, idx);
        }
    }
}
