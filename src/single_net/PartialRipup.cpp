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
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE)
        log() << "Extracting pseudo nets from violated nets " << std::endl;

    if (db::rrrIterSetting.localRipup) {
        runJobsMT(vioNets.size(), [&](int id) {
            auto &net = database.nets[vioNets[id]];
            markLocalRipup(net);
        });
    }
    else if (db::rrrIterSetting.adaptiveRipup) {
        runJobsMT(vioNets.size(), [&](int id) {
            auto &net = database.nets[vioNets[id]];
            mergeRouteGuides(net);
            for (auto &g : net.mergedGuides)
                database.expandBox(g, db::rrrIterSetting.defaultGuideExpand);
            markAdaptiveRipup(net);
        });
    }

    for (int id : vioNets) {
        for (auto proot : database.nets[id].pnets)
            removeDBEdges(proot, id);
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
        merge(root, 8);
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
    net.postOrderVisitGridTopo([&](std::shared_ptr<db::GridSteiner> node) {
        const auto &np = database.getLoc(*node);
        for (size_t i=0; i<guideSize; i++) {
            if (selected[i] > 0 && net.mergedGuides[i].Contain(np)) {
                node->isVio = true;
                break;
            }
        }
    });

    for (auto root : net.gridTopo)
        merge(root, 40);
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

void PartialRipup::merge(std::shared_ptr<db::GridSteiner> node, int thre) {
    for (auto c : node->children) merge(c, thre);

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
            if ((danger == nullptr) || (node->layerIdx != danger->layerIdx))
                node->distance = d1;
            else
                node->distance = abs(node->crossPointIdx - danger->crossPointIdx) +
                                    abs(node->trackIdx - danger->trackIdx) + d1;
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


void printgrids(std::shared_ptr<db::GridSteiner> node) {
    // const auto &np = database.getLoc(*node);
    // std::cout << "  " << node->layerIdx << ", " << np << std::endl;
    printf(" %d, %d, %d\n", node->layerIdx, node->trackIdx, node->crossPointIdx);
    for (auto c : node->children) {
            printgrids(c);
    }
}

void extractClosedVias(db::GridSteiner *node,
                        vector<std::pair<db::GridSteiner *, db::GridSteiner *>> &pairs) {
    if (!node->parent) node->distance = 6;
    else {
        auto parent = node->parent.get();
        if (node->viaType) {
            node->distance = 0;
            if (parent->distance < 4) {
                while (parent->distance > 0) parent = parent->parent.get();
                pairs.push_back({parent, node});
            }
        }
        else {
            node->distance = abs(node->crossPointIdx - parent->crossPointIdx) + 
                                abs(node->trackIdx - parent->trackIdx) +
                                parent->distance;
        }
    }

    for (auto c : node->children) extractClosedVias(c.get(), pairs);
}

inline bool connSameLayerFailed(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx, bool debug) {
    std::shared_ptr<db::GridSteiner> turn;
    utils::IntervalT<int> trange, cprange;
    int lid = u->layerIdx, vio = 0;
    if (u->trackIdx > v->trackIdx) trange.Set(v->trackIdx, u->trackIdx);
    else trange.Set(u->trackIdx, v->trackIdx);
    if (u->crossPointIdx > v->crossPointIdx) cprange.Set(v->crossPointIdx, u->crossPointIdx);
    else cprange.Set(u->crossPointIdx, v->crossPointIdx);

    if (trange.range() == 0 || cprange.range() == 0) {
        u->children.push_back(v);
        v->parent = u;
    }
    else {
        db::GridPoint t(lid, u->trackIdx, u->crossPointIdx);
        vio += database.getWireSegmentVioCost({lid, t.trackIdx, cprange}, netIdx, false);
        if (debug) printf(" u -> %d,%d; vio : %d\n", t.trackIdx, t.crossPointIdx, vio);
        vio += database.getWrongWayWireSegmentVioCost({lid, trange, v->crossPointIdx}, netIdx, false);
        if (debug) printf(" %d,%d -> v; vio : %d\n", t.trackIdx, t.crossPointIdx, vio);
        if (vio == 0)
            t.crossPointIdx = v->crossPointIdx;
        else {
            vio = 0;
            vio += database.getWrongWayWireSegmentVioCost({lid, trange, t.crossPointIdx}, netIdx, false);
            if (debug) printf(" u -> %d,%d; vio : %d\n", t.trackIdx, t.crossPointIdx, vio);
            vio += database.getWireSegmentVioCost({lid, v->trackIdx, cprange}, netIdx, false);
            if (debug) printf(" %d,%d -> v; vio : %d\n", t.trackIdx, t.crossPointIdx, vio);
            t.trackIdx = v->trackIdx;
        }
        if (vio > 0) return true;
        turn = std::make_shared<db::GridSteiner>(t);
        u->children.push_back(turn);
        turn->parent = u;
        v->parent = turn;
        turn->children.push_back(v);
        database.useEdge({*u, *turn}, netIdx);
        database.useEdge({*turn, *v}, netIdx);
        // std::cout << *u << "," << *turn << "," << *v << std::endl;
    }
    v->viaType = nullptr;
    return false;
}

void PartialRipup::removeSmallLayerSwitch(db::Net &net) {
    vector<std::pair<db::GridSteiner *, db::GridSteiner *>> pairs;
    for (auto root : net.gridTopo) {
        extractClosedVias(root.get(), pairs);
    }
    // if (net.getName() == "net48841") {
    //     printf(" net 48841 %d\n\n\n\n", net.idx);
    //     for (const auto &p : pairs) {
    //         const auto &first = *p.first;
    //         const auto &second = *p.second;
    //         printf(" f : %d, %d, %d;  s : %d, %d, %d\n", first.layerIdx, first.trackIdx, first.crossPointIdx,
    //                                                         second.layerIdx, second.trackIdx, second.crossPointIdx);
    //         const auto &up = database.getLoc(first);
    //         const auto &vp = database.getLoc(second);
    //         printf(" u : %d, %ld, %ld;  v : %d, %ld, %ld\n", first.layerIdx, up.x/2, up.y/2, second.layerIdx, vp.x/2, vp.y/2);
    //     }
    // }
    for (auto &pair : pairs) {
        std::shared_ptr<db::GridSteiner> u = pair.first->parent, v;
        for (auto c : pair.second->parent->children) if (c.get() == pair.second) v = c;

        // // // debug
        // if (net.getName() == "net48841" && pair.first->trackIdx == 7931) {
        //     printf(" net 48841 %d\n\n\n\n", net.idx);
        // //     // const auto &up = database.getLoc(*u);
        // //     // const auto &vp = database.getLoc(*v);
        // //     // printf(" net %d vio : %d\n", net.idx, vio);
        // //     printf(" net %d; unull : %d; vnull : %d\n", net.idx, int(u==nullptr), int(v==nullptr));
        //     const auto &first = *pair.first;
        //     const auto &second = *pair.second;
        //     printf(" f : %d, %d, %d;  s : %d, %d, %d  ", first.layerIdx, first.trackIdx, first.crossPointIdx,
        //                                                     second.layerIdx, second.trackIdx, second.crossPointIdx);
        //     printf(" u : %d, %d, %d;  v : %d, %d, %d\n\n\n\n", u->layerIdx, u->trackIdx, u->crossPointIdx,
        //                                                     v->layerIdx, v->trackIdx, v->crossPointIdx);
            // connSameLayerFailed(u, v, net.idx, true);
        // //     if (v == nullptr) {
        // //         auto tt = pair.second->parent;
        // //         for (auto c : tt->children) printf(" %d,%d,%d; ", c->layerIdx, c->trackIdx, c->crossPointIdx);
        // //         printf("\n\n\n\n");
        // //     }
        // //     if (v->trackIdx == 2250 && v->crossPointIdx == 304) {
        // //         printgrids(u);
        // //     }
        // //     // printf(" u : %d, %ld, %ld;  v : %d, %ld, %ld\n", u->layerIdx, up.x/2, up.y/2, v->layerIdx, vp.x/2, vp.y/2);
        // //     // printf("    turn : %d, %d, %d\n", t.layerIdx, t.trackIdx, t.crossPointIdx);
        // //     auto temp = v;
        // //     while (temp != u) {printf(" %d,%d,%d; ", temp->layerIdx, temp->trackIdx, temp->crossPointIdx); temp = temp->parent;}
        // //     printf("\n\n");
        // }

        auto up = v->parent;
        if (u->layerIdx != v->layerIdx ||
            abs(u->trackIdx-v->trackIdx) > 1 ||
            abs(pair.first->trackIdx-pair.second->parent->trackIdx) > 1 ||
            connSameLayerFailed(u, v, net.idx, false))
            continue;

        // the latter via should always be removed
        database.removeEdge({*v, *up}, net.idx);
        if (up->extWireSeg) database.removeEdge(*(up->extWireSeg), net.idx);
        // if (v->extWireSeg) database.removeEdge(*(v->extWireSeg), net.idx);
        // v->extWireSeg = nullptr;
        auto &cs = up->children;
        for (auto it=cs.begin(); it!=cs.end(); it++) {
            if ((*it)==v) {
                *it = cs.back();
                cs.pop_back();
                break;
            }
        }

            // printf(" %d,", net.idx);

            // if (net.getName() == "net48841" && pair.first->trackIdx == 7931) {
            //     printf(" net 48841 %d\n\n\n\n", net.idx);
            //     auto c = u->children[0];
            //     printf(" uc : %d, %d ,%d\n", c->layerIdx, c->trackIdx, c->crossPointIdx);
            // }
    }

    net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
        if (node->pinIdx < 0 && node->children.empty()) {
            database.removeEdge({*node, *node->parent}, net.idx);
            if (node->extWireSeg) database.removeEdge(*(node->extWireSeg), net.idx);
            auto &cs = node->parent->children;
            for (auto it=cs.begin(); it!=cs.end(); it++) {
                if ((*it)==node) {
                    *it = cs.back();
                    cs.pop_back();
                    break;
                }
            }
            node->parent = nullptr;
        }
    });
}

void PartialRipup::handlePinSplitVias(db::Net &net) {
    vector<vector<std::shared_ptr<db::GridSteiner>>> pinTaps(net.numOfPins());
    std::shared_ptr<db::GridSteiner> u, v, temp;
    // int numOfOrphal = 0;

    net.postOrderVisitGridTopo([&pinTaps](std::shared_ptr<db::GridSteiner> node) {
        if (node->pinIdx >=0 && node->layerIdx == 0) {
            pinTaps[node->pinIdx].push_back(node);
        }
    });
    for (size_t i=0; i<pinTaps.size(); i++) {
        if (pinTaps[i].size() != 2) continue;
        u = pinTaps[i][0];
        v = pinTaps[i][1];

            // printf("net %d pin %d (%d,%d,%d) (%d,%d,%d) \n",
            //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
            //                             v->layerIdx, v->trackIdx, v->crossPointIdx);

            if ((u->parent!=nullptr && u->children.size() > 0) ||
                (v->parent!=nullptr && v->children.size() > 0)) {
                // printf("Warning : net %d pin %d (%d,%d,%d) (%d,%d,%d) not terminal\n\n",
                //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
                //                             v->layerIdx, v->trackIdx, v->crossPointIdx);
                continue;
            }
            if ((u->parent==nullptr && u->children.empty()) ||
                (v->parent==nullptr && v->children.empty())) {
                // printf("Warning : net %d pin %d (%d,%d,%d) (%d,%d,%d) orphal\n",
                //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
                //                             v->layerIdx, v->trackIdx, v->crossPointIdx);
                // numOfOrphal++;
                continue;
            }

        if (abs(u->trackIdx-v->trackIdx) > 1 ||
            abs(u->crossPointIdx-v->crossPointIdx) > 1)
            continue;

        // both taps are arriving at this pin, conn new gird at v.upper to u, removed via at v
        if (u->children.empty() && v->children.empty()) {
            u = u->parent;
            temp = std::make_shared<db::GridSteiner>(database.getUpper(*v));
            if (u->layerIdx != 1 ||
                v->parent->layerIdx != 1 ||
                connSameLayerFailed(u, temp, net.idx, false))
                continue;
            database.removeEdge({*v, *(v->parent)}, net.idx);
            auto &cs = v->parent->children;
            for (auto it=cs.begin(); it!=cs.end(); it++) {
                if ((*it)==v) {
                    *it = cs.back();
                    cs.pop_back();
                    break;
                }
            }
        }
        // keep v the children of a root, conn v to u, removed via at v
        else {
            if (v->children.empty())
                u.swap(v);
            temp = v;
            v = v->children[0];
            u = u->parent == nullptr ? u->children[0] : u->parent;
            if (u->layerIdx != 1 ||
                v->layerIdx != 1 ||
                connSameLayerFailed(u, v, net.idx, false))
                continue;
            database.removeEdge({*v, *temp}, net.idx);
            auto &cs = net.gridTopo;
            for (auto it=cs.begin(); it!=cs.end(); it++) {
                if ((*it)==temp) {
                    *it = cs.back();
                    cs.pop_back();
                    break;
                }
            }
        }

            // printf(" %d:%d;", net.idx, int(i));
    }

    // printf(" number of orphal : %d\n", numOfOrphal);
}

void PartialRipup::handleMacroPins(db::Net &net) {
    return;
}

inline bool notCorner(const db::GridPoint &u, const db::GridPoint &v) {
    return u.trackIdx == v.trackIdx || u.crossPointIdx == v.crossPointIdx;
}

void PartialRipup::removeCorners(std::shared_ptr<db::GridSteiner> node, int idx) {
    int step = 3;
    bool remove = true;
    auto v = node;
    // MYTODO : last node can have multiple children
    if (node->parent) {
        while (step--) {
            if (v->children.size() != 1 ||
                v->layerIdx != node->layerIdx ||
                notCorner(*(v->parent), *(v->children[0]))) {
                remove = false;
                break;
            }
            v = v->children[0];
        }
        int dist = abs(node->trackIdx-v->trackIdx) + abs(node->crossPointIdx-v->crossPointIdx);
        if (remove && dist < 3) {
            // printf(" u : %d,%d,%d; v : %d,%d,%d\n", node->layerIdx,node->trackIdx,node->crossPointIdx,
            //                                         v->layerIdx, v->trackIdx,v->crossPointIdx);
            auto t = v;
            while (t != node) {
                // printf("  %d,%d,%d\n", t->layerIdx,t->trackIdx,t->crossPointIdx);
                database.removeEdge({*t, *(t->parent)}, idx);
                if (t->extWireSeg) database.removeEdge(*(t->extWireSeg), idx);
                t = t->parent;
            }
            if (dist == 0) {
                node->children = move(v->children);
                for (auto c : node->children) c->parent = node;
            }
            else {
                node->children[0] = v;
                v->parent = node;
            }
        }
    }

    for (auto c : node->children) removeCorners(c, idx);
}

void PartialRipup::fixMAR(db::Net &net) {

    vector<db::GridSteiner *> otherPaths, bfsCurr, bfsNext;
    vector<std::pair<db::GridSteiner *, db::GridSteiner *>> unSafePaths;

    for (auto root : net.gridTopo)
        otherPaths.push_back(root.get());
    
    while (!otherPaths.empty()) {
        db::GridSteiner *begin = otherPaths.back(), *end = nullptr;
        otherPaths.pop_back();
        bool unsafe = begin->extWireSeg == nullptr && begin->parent != nullptr;
        int layerIdx = begin->layerIdx;

        // bfs for a path on 1 layer
        bfsCurr.push_back(begin);
        while (!bfsCurr.empty()) {
            for (auto n : bfsCurr) {
                if (n->children.empty()) {
                    unsafe = false; // assume all leaves are pin conns
                    continue;
                }
                for (auto c : n->children) {
                    if (c->layerIdx == layerIdx)
                        bfsNext.push_back(c.get());
                    else {
                        otherPaths.push_back(c.get());
                        if (end == nullptr) {
                            end = c.get();
                            unsafe = unsafe && n->extWireSeg == nullptr;
                        }
                        else
                            unsafe = false;
                    }
                }
            }
            bfsCurr = move(bfsNext);
        }
        if (unsafe)
            unSafePaths.emplace_back(begin, end);
    }

    // for (const auto &pair : unSafePaths)
    //     if (pair.first->parent == nullptr || pair.second->parent == nullptr)
    //         printf(" ERROR : unSafe paths has nullptr\n");
    // printf("net %d extracted\n", net.idx);
    // printf(" num of unsafe paths : %d\n", int(unSafePaths.size()));

    for (const auto &pair : unSafePaths) {
        const auto &u = *pair.first;
        const auto &v = *pair.second;   // node on adj layer
        const auto &p = *(v.parent);
        const auto &uloc = database.getLoc(u);
        const auto &vloc = database.getLoc(v);
        const auto &layer = database.getLayer(u.layerIdx);
        int lid = u.layerIdx, dir = 1 - layer.direction, tempCP = -1,
            numGrids = abs(u.trackIdx-p.trackIdx) + abs(u.crossPointIdx-p.crossPointIdx);   // wirelength estimation by hpwl
        DBU mina = layer.minArea, width = layer.width, patchArea = 0, patchLen = 0;

        // get patch area, available patch intervals
        utils::BoxT<DBU> patchBox;
        vector<utils::BoxT<DBU>> viaBoxes;
        vector<utils::IntervalT<DBU>> availItvls;
        if (u.parent->layerIdx > lid)
            viaBoxes.emplace_back(u.viaType->getShiftedBotMetal(uloc));
        else
            viaBoxes.emplace_back(u.viaType->getShiftedTopMetal(uloc));
        if (v.layerIdx > lid)
            viaBoxes.emplace_back(v.viaType->getShiftedBotMetal(vloc));
        else
            viaBoxes.emplace_back(v.viaType->getShiftedTopMetal(vloc));
        availItvls.emplace_back(database.getEmptyRange(lid, u.trackIdx, u.crossPointIdx, net.idx));
        availItvls.emplace_back(database.getEmptyRange(lid, p.trackIdx, p.crossPointIdx, net.idx));
        if (numGrids == 0) {
            const auto &its = viaBoxes[0].IntersectWith(viaBoxes[1]);
            patchArea = mina - viaBoxes[0].area() - viaBoxes[1].area() + its.area();
        }
        if (numGrids > 0) {
            patchArea = mina - viaBoxes[0].area() - viaBoxes[1].area() - numGrids * layer.pitch * width;
            tempCP = u.children[0]->crossPointIdx;
            if (tempCP == u.crossPointIdx) patchArea += viaBoxes[0][1-dir].range() * width / 2;
            else {
                patchArea += viaBoxes[0][dir].range() * width / 2;
                if (tempCP < u.crossPointIdx) availItvls[0].low = uloc[dir];
                else availItvls[0].high = uloc[dir];
            }
            tempCP = p.parent->crossPointIdx;
            if (tempCP == p.crossPointIdx) patchArea += viaBoxes[1][1-dir].range() * width / 2;
            else {
                patchArea += viaBoxes[1][dir].range() * width / 2;
                if (tempCP < p.crossPointIdx) availItvls[1].low = vloc[dir];
                else availItvls[1].high = vloc[dir];
            }
        }

            // printf(" u : %d,%d,%d;  v : %d,%d,%d\n", u.layerIdx, u.trackIdx, u.crossPointIdx,
            //                                             p.layerIdx, p.trackIdx, p.crossPointIdx);
            // printf(" u : %ld,%ld;  v : %ld,%ld\n", uloc.x, uloc.y, vloc.x, vloc.y);
            // printf("%ld, %ld, %d\n", mina, patchArea, numGrids);

        // patch
        if (patchArea > 0) {
            for (size_t i=0; i<2; i++) {
                if (viaBoxes[i][dir].range() < viaBoxes[i][1-dir].range()) continue;
                patchLen = patchArea / viaBoxes[i][1-dir].range();

                    // printf(" %ld, %ld\n", patchLen, availItvls[i].range());
                    // std::cout << viaBoxes[i] << " ; " << availItvls[i] << std::endl;

                if (patchLen + viaBoxes[i][dir].range() > availItvls[i].range()) continue;
                if (viaBoxes[i][dir].low - availItvls[i].low >= patchLen)
                    viaBoxes[i][dir].low -= patchLen;
                else {
                    viaBoxes[i][dir].low = availItvls[i].low;
                    viaBoxes[i][dir].high += (patchLen - viaBoxes[i][dir].low + availItvls[i].low);
                }
                database.writeDEFFillRect(net, viaBoxes[i], lid);
                break;
            }
        }
    }
}
