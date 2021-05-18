#include "PostProcess.h"

void printgrids(std::shared_ptr<db::GridSteiner> node) {
    // const auto &np = database.getLoc(*node);
    // std::cout << "  " << node->layerIdx << ", " << np << std::endl;
    printf(" %d, %d, %d\n", node->layerIdx, node->trackIdx, node->crossPointIdx);
    for (auto c : node->children) {
            printgrids(c);
    }
}

void extractVias(int treeId,
                    std::shared_ptr<db::GridSteiner> node,
                    vector<std::shared_ptr<db::GridSteiner>> &allVias,
                    vector<std::pair<int, int>> &closeVias,
                    vector<std::map<int, std::map<int, int>>> &vias) {
    node->distance = treeId;
    if (node->viaType != nullptr) {
        auto lower = node->layerIdx > node->parent->layerIdx ? node->parent : node;
        int viaId = allVias.size();
        allVias.push_back(node);
        auto &Vias = vias[lower->layerIdx];
        std::map<int, std::map<int, int>>::iterator it;
        for (int t=lower->trackIdx-1; t<lower->trackIdx+2; t++) {
            it = Vias.find(t);
            if (it == Vias.end()) continue;
            for (int c=lower->crossPointIdx-1; c<lower->crossPointIdx+2; c++) {
                if (it->second.count(c))
                    closeVias.push_back(std::make_pair(it->second[c], viaId));
            }
        }
        if (!Vias.count(lower->trackIdx))
            Vias.insert({lower->trackIdx, std::map<int, int>()});
        Vias[lower->trackIdx].insert({lower->crossPointIdx, viaId});
    }

    for (auto c : node->children)
        extractVias(treeId, c, allVias, closeVias, vias);
}

inline bool connSameLayer(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx) {
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
        vio += database.getWrongWayWireSegmentVioCost({lid, trange, v->crossPointIdx}, netIdx, false);
        if (vio == 0)
            t.crossPointIdx = v->crossPointIdx;
        else {
            vio = 0;
            vio += database.getWrongWayWireSegmentVioCost({lid, trange, t.crossPointIdx}, netIdx, false);
            vio += database.getWireSegmentVioCost({lid, v->trackIdx, cprange}, netIdx, false);
            t.trackIdx = v->trackIdx;
        }
        if (vio > 0) return false;
        turn = std::make_shared<db::GridSteiner>(t);
        u->children.push_back(turn);
        turn->parent = u;
        v->parent = turn;
        turn->children.push_back(v);
        database.useEdge({*u, *turn}, netIdx);
        database.useEdge({*turn, *v}, netIdx);
    }
    v->viaType = nullptr;
    return true;
}

void PostProcess::removeSameNetVioVias(db::Net &net) {
    vector<std::shared_ptr<db::GridSteiner>> allVias;
    vector<std::pair<int, int>> closeVias;
    vector<std::map<int, std::map<int, int>>> vias(database.getLayerNum());
    for (int i=0; i<net.gridTopo.size(); i++) {
        extractVias(i, net.gridTopo[i], allVias, closeVias, vias);
    }

    for (const auto &pair : closeVias) {
        auto u = allVias[pair.first], v = allVias[pair.second];
        if (u->viaType == nullptr || v->viaType == nullptr) continue;
        if (u->distance == v->distance)
            removeSameTreeVias(u, v, net.idx);
        else
            removeDiffTreeVias(u, v, net.idx);
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
    vector<std::shared_ptr<db::GridSteiner>> temp;
    for (size_t i=0; i<net.gridTopo.size(); i++) {
        if (net.gridTopo[i]->children.size() > 0)
            temp.push_back(net.gridTopo[i]);
    }
    net.gridTopo = std::move(temp);
}

void PostProcess::removeSameTreeVias(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx) {
    // assume : u->lid == v->lid --> in different sub-trees
    //          u->parent->lid = v->lid --> u is in sub-tree of v
    db::GridSteiner *toRemove = v->parent.get();
    bool fixed = false;
    if (u->layerIdx == v->layerIdx)
        fixed = connSameLayer(u, v, netIdx);
    else if (u->parent->layerIdx == v->layerIdx)
        fixed = connSameLayer(u->parent, v, netIdx);
    if (fixed) {
        database.removeEdge({*v, *toRemove}, netIdx);
        if (toRemove->extWireSeg) database.removeEdge(*(toRemove->extWireSeg), netIdx);
        for (size_t i=0; i<toRemove->children.size(); i++) {
            if (toRemove->children[i] == v) {
                toRemove->children[i] = toRemove->children.back();
                toRemove->children.pop_back();
                break;
            }
        }
    }
    else {
        std::cout << database.nets[netIdx].getName();
        printf(" remove same tree vias failed\n");
    }
}

void PostProcess::removeDiffTreeVias(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx) {
    // assume : close to metal 1 pins
    // children -> metal 1 --> ahead for pin
    // children -> metal 2 --> away from pin
    // remove the away via if any
    if (u->layerIdx > 1 || v->layerIdx > 1) return;
    if (u->layerIdx == 0 && v->layerIdx == 0) {
        std::shared_ptr<db::GridSteiner> node = std::make_shared<db::GridSteiner>(database.getUpper(*v));
        if (connSameLayer(u->parent, node, netIdx)) {
            auto p = v->parent.get();
            p->fakePin = v->pinIdx > -1;
            p->pinIdx = v->pinIdx;
            node->fakePin = p->fakePin;
            node->pinIdx = v->pinIdx;
            database.removeEdge({*v, *p}, netIdx);
            if (p->extWireSeg) database.removeEdge(*(p->extWireSeg), netIdx);
            v->parent = nullptr;
            for (size_t i=0; i<p->children.size(); i++) {
                if (p->children[i] == v) {
                    p->children[i] = p->children.back();
                    p->children.pop_back();
                    break;
                }
            }
        }
        else {
            std::cout << database.nets[netIdx].getName();
            printf(" remove diff tree vias (%d, %d) failed\n", u->layerIdx, v->layerIdx);
        }
    }
    else {
        if (u->layerIdx == 1) u.swap(v);  // make sure v is a children on metal 2
        if (u->layerIdx == 0) u = u->parent;
        auto toRemove = v->parent.get();
        if (connSameLayer(u, v, netIdx)) {
            database.removeEdge({*v, *toRemove}, netIdx);
            if (toRemove->extWireSeg) database.removeEdge(*(toRemove->extWireSeg), netIdx);
            for (size_t i=0; i<toRemove->children.size(); i++) {
                if (toRemove->children[i] == v) {
                    toRemove->children[i] = toRemove->children.back();
                    toRemove->children.pop_back();
                    break;
                }
            }
        }
        else {
            std::cout << database.nets[netIdx].getName();
            printf(" remove diff tree vias (%d, %d) failed\n", u->layerIdx, v->layerIdx);
        }
    }
}

inline bool notCorner(const db::GridPoint &u, const db::GridPoint &v) {
    return u.trackIdx == v.trackIdx || u.crossPointIdx == v.crossPointIdx;
}

void PostProcess::removeCorners(std::shared_ptr<db::GridSteiner> node, int idx) {
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

void PostProcess::fixMAR(db::Net &net) {

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

// void getSameTreeVias(db::GridSteiner *node,
//                         vector<std::pair<db::GridSteiner *, db::GridSteiner *>> &pairs) {
//     if (!node->parent) node->distance = 6;
//     else {
//         auto parent = node->parent.get();
//         if (node->viaType) {
//             node->distance = 0;
//             if (parent->distance < 4) {
//                 while (parent->distance > 0) parent = parent->parent.get();
//                 pairs.push_back({parent, node});
//             }
//         }
//         else {
//             node->distance = abs(node->crossPointIdx - parent->crossPointIdx) + 
//                                 abs(node->trackIdx - parent->trackIdx) +
//                                 parent->distance;
//         }
//     }

//     for (auto c : node->children) getSameTreeVias(c.get(), pairs);
// }

// void PostProcess::removeSameNetVioVias(db::Net &net) {
//     vector<std::pair<db::GridSteiner *, db::GridSteiner *>> pairs;
//     for (auto root : net.gridTopo) {
//         getSameTreeVias(root.get(), pairs);
//     }
//     for (auto &pair : pairs) {
//         std::shared_ptr<db::GridSteiner> u = pair.first->parent, v;
//         for (auto c : pair.second->parent->children) if (c.get() == pair.second) v = c;

//         auto up = v->parent;
//         if (u->layerIdx != v->layerIdx ||
//             abs(u->trackIdx-v->trackIdx) > 1 ||
//             abs(pair.first->trackIdx-pair.second->parent->trackIdx) > 1 ||
//             connSameLayerFailed(u, v, net.idx))
//             continue;

//         // the latter via should always be removed
//         database.removeEdge({*v, *up}, net.idx);
//         if (up->extWireSeg) database.removeEdge(*(up->extWireSeg), net.idx);
//         auto &cs = up->children;
//         for (auto it=cs.begin(); it!=cs.end(); it++) {
//             if ((*it)==v) {
//                 *it = cs.back();
//                 cs.pop_back();
//                 break;
//             }
//         }
//     }

//     net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
//         if (node->pinIdx < 0 && node->children.empty()) {
//             database.removeEdge({*node, *node->parent}, net.idx);
//             if (node->extWireSeg) database.removeEdge(*(node->extWireSeg), net.idx);
//             auto &cs = node->parent->children;
//             for (auto it=cs.begin(); it!=cs.end(); it++) {
//                 if ((*it)==node) {
//                     *it = cs.back();
//                     cs.pop_back();
//                     break;
//                 }
//             }
//             node->parent = nullptr;
//         }
//     });
// }

// void PostProcess::removeDiffTreeVioVias(db::Net &net) {
//     vector<vector<std::shared_ptr<db::GridSteiner>>> pinTaps(net.numOfPins());
//     std::shared_ptr<db::GridSteiner> u, v, temp;
//     // int numOfOrphal = 0;

//     net.postOrderVisitGridTopo([&pinTaps](std::shared_ptr<db::GridSteiner> node) {
//         if (node->pinIdx >=0 && node->layerIdx == 0) {
//             pinTaps[node->pinIdx].push_back(node);
//         }
//     });
//     for (size_t i=0; i<pinTaps.size(); i++) {
//         if (pinTaps[i].size() != 2) continue;
//         u = pinTaps[i][0];
//         v = pinTaps[i][1];

//             // printf("net %d pin %d (%d,%d,%d) (%d,%d,%d) \n",
//             //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
//             //                             v->layerIdx, v->trackIdx, v->crossPointIdx);

//             if ((u->parent!=nullptr && u->children.size() > 0) ||
//                 (v->parent!=nullptr && v->children.size() > 0)) {
//                 // printf("Warning : net %d pin %d (%d,%d,%d) (%d,%d,%d) not terminal\n\n",
//                 //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
//                 //                             v->layerIdx, v->trackIdx, v->crossPointIdx);
//                 continue;
//             }
//             if ((u->parent==nullptr && u->children.empty()) ||
//                 (v->parent==nullptr && v->children.empty())) {
//                 // printf("Warning : net %d pin %d (%d,%d,%d) (%d,%d,%d) orphal\n",
//                 //         net.idx, int(i), u->layerIdx, u->trackIdx, u->crossPointIdx,
//                 //                             v->layerIdx, v->trackIdx, v->crossPointIdx);
//                 // numOfOrphal++;
//                 continue;
//             }

//         if (abs(u->trackIdx-v->trackIdx) > 1 ||
//             abs(u->crossPointIdx-v->crossPointIdx) > 1)
//             continue;

//         // both taps are arriving at this pin, conn new gird at v.upper to u, removed via at v
//         if (u->children.empty() && v->children.empty()) {
//             u = u->parent;
//             temp = std::make_shared<db::GridSteiner>(database.getUpper(*v));
//             if (u->layerIdx != 1 ||
//                 v->parent->layerIdx != 1 ||
//                 connSameLayerFailed(u, temp, net.idx))
//                 continue;
//             database.removeEdge({*v, *(v->parent)}, net.idx);
//             auto &cs = v->parent->children;
//             for (auto it=cs.begin(); it!=cs.end(); it++) {
//                 if ((*it)==v) {
//                     *it = cs.back();
//                     cs.pop_back();
//                     break;
//                 }
//             }
//         }
//         // keep v the children of a root, conn v to u, removed via at v
//         else {
//             if (v->children.empty())
//                 u.swap(v);
//             temp = v;
//             v = v->children[0];
//             u = u->parent == nullptr ? u->children[0] : u->parent;
//             if (u->layerIdx != 1 ||
//                 v->layerIdx != 1 ||
//                 connSameLayerFailed(u, v, net.idx))
//                 continue;
//             database.removeEdge({*v, *temp}, net.idx);
//             auto &cs = net.gridTopo;
//             for (auto it=cs.begin(); it!=cs.end(); it++) {
//                 if ((*it)==temp) {
//                     *it = cs.back();
//                     cs.pop_back();
//                     break;
//                 }
//             }
//         }

//             // printf(" %d:%d;", net.idx, int(i));
//     }

//     // printf(" number of orphal : %d\n", numOfOrphal);
// }