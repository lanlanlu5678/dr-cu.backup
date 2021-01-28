#include "PartialRipup.h"

void PartialRipup::preProcessRouteGuides() {
    // std::mutex plock;
    // int expandTimes = 0;
    auto expand = [](int id) {
        auto &net = database.nets[id];
        auto &guides = net.routeGuides;
        int oriSize = guides.size(), localTimes = 0;
        vector<vector<utils::BoxT<DBU>>> gol(database.getLayerNum());
        for (int i=0; i<oriSize; i++) {
            int l = guides[i].layerIdx;
            gol[l].push_back(guides[i]);
            auto tdir = database.getLayerDir(l);
            for (int j=i+1; j<oriSize; j++) {
                if (guides[j].layerIdx!=l) continue;
                auto box = guides[i].IntersectWith(guides[j]);
                if ((box[tdir].range() == 0) && (box[1-tdir].range() > 0)) {
                    int longer = i, shorter = j;
                    if (guides[i][1-tdir].range() < guides[j][1-tdir].range()) std::swap(longer, shorter);
                    box[tdir] = guides[longer][tdir];
                    box[1-tdir] = guides[shorter][1-tdir];
                    gol[l].push_back(box);
                    localTimes++;
                }
            }
        }
        guides.clear();
        for (int l=0; l<database.getLayerNum(); l++) {
            utils::SlicePolygons(gol[l], 1-database.getLayerDir(l));
            for (auto &box : gol[l])
                guides.emplace_back(l, box);
        }
        net.routeGuideVios.resize(guides.size());
        // plock.lock();
        // expandTimes += localTimes;
        // plock.unlock();
    };

    runJobsMT(database.nets.size(), expand);
    // log() << std::endl;
    // log() << "Expand " << expandTimes << " times" << std::endl;
}

// void removeDbEdges(std::shared_ptr<db::GridSteiner> node, int netIdx) {
//     for (auto c : node->children) {
//         if (c->isVio) {
//             database.removeEdge({*c, *node}, netIdx);
//             if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), netIdx);
//             removeDbEdges(c, netIdx);
//         }
//     }
// }

// void mergePNetsInSubtrees(std::shared_ptr<db::GridSteiner> node) {
//     for (auto c : node->children)
//         if (!(c->isVio) && c->distance<5) {
//             c->isVio = true;
//             mergePNetsInSubtrees(c);
//         }
// }

// void merge(std::shared_ptr<db::GridSteiner> node) {
//     for (auto c : node->children) merge(c);

//     int d = 6, d1 = 6, d2 = 6;
//     db::GridSteiner *danger = nullptr;
    
//     if (node->isVio) {
//         mergePNetsInSubtrees(node);
//         node->distance = -1;
//     }
//     else {
//         for (auto c : node->children) {
//             d = c->distance;

//             if (d1 > d) {
//                 d2 = d1;
//                 d1 = d;
//                 danger = c.get();
//             }
//             else if (d2 > d) {
//                 d2 = d;
//             }
//         }

//         if (d1 + d2 < 4) {
//             mergePNetsInSubtrees(node);
//             node->distance = 0;
//         }
//         else {
//             if ((danger == nullptr) || (node->layerIdx != danger->layerIdx))
//                 node->distance = d1;
//             else
//                 node->distance = abs(node->crossPointIdx - danger->crossPointIdx) + d1;
//         }
//     }
// }

// void mergeDistEdge(std::shared_ptr<db::GridSteiner> node) {
//     for (auto c : node->children) mergeDistEdge(c);

//     int d = 6, d1 = 6, d2 = 6;
    
//     if (node->isVio) {
//         mergePNetsInSubtrees(node);
//         node->distance = -1;
//     }
//     else {
//         for (auto c : node->children) {
//             d = c->distance;

//             if (d1 > d) {
//                 d2 = d1;
//                 d1 = d;
//             }
//             else if (d2 > d) {
//                 d2 = d;
//             }
//         }

//         if (d1 + d2 < 4) {
//             mergePNetsInSubtrees(node);
//             node->distance = 0;
//         }
//         else {
//             node->distance = d1 + 1;
//         }
//     }
// }

// void quickFix(db::Net &dbNet) {

//     dbNet.gridRouteGuides.clear();

//     for (auto vio : dbNet.vioNodes) {
//         int l = vio->layerIdx, lt = vio->trackIdx-db::rrrIterSetting.extraTracks,
//             ut = vio->trackIdx+db::rrrIterSetting.extraTracks, lc = vio->crossPointIdx-db::rrrIterSetting.extraCPs,
//             uc = vio->crossPointIdx+db::rrrIterSetting.extraCPs;
//         std::map<int, db::GridBoxOnLayer> boxes;
//         db::GridBoxOnLayer gbox(l,utils::IntervalT<int>(lt,ut),utils::IntervalT<int>(lc,uc));
//         if (!database.isValid(gbox)) {
//             if (gbox.trackRange.low < 0) gbox.trackRange.low = 0;
//             if (gbox.trackRange.high > database.getTrackLimit(l)) gbox.trackRange.high = database.getTrackLimit(l);
//             if (gbox.crossPointRange.low < 0) gbox.crossPointRange.low = 0;
//             if (gbox.crossPointRange.high > database.getCPLimit(l)) gbox.crossPointRange.high = database.getCPLimit(l);
//         }
//         boxes.insert({l, gbox});
//         auto box = database.getLoc(boxes[l]);

//         std::function<void(std::shared_ptr<db::GridSteiner>)> markBranches;
//         markBranches = [&](std::shared_ptr<db::GridSteiner> node) {
//             for (auto c : node->children) {
//                 if (c->isVio) continue;
//                 if (!boxes.count(c->layerIdx)) {
//                     box.layerIdx = c->layerIdx;
//                     boxes.insert({c->layerIdx, database.rangeSearch(box)});
//                 }
//                 if (boxes[c->layerIdx].crossPointRange.Contain(c->crossPointIdx)) {
//                     c->isVio = true;
//                     markBranches(c);
//                 }
//                 else if (c->layerIdx == node->layerIdx &&
//                         abs(c->crossPointIdx-node->crossPointIdx) > 1 &&
//                         boxes[c->layerIdx].crossPointRange.StrictlyContain(node->crossPointIdx)) { // long edge
//                     std::shared_ptr<db::GridSteiner> p;
//                     if (c->crossPointIdx < node->crossPointIdx)
//                         p = std::make_shared<db::GridSteiner>(db::GridPoint(c->layerIdx, c->trackIdx,
//                                                                             boxes[c->layerIdx].crossPointRange.low));
//                     else
//                         p = std::make_shared<db::GridSteiner>(db::GridPoint(c->layerIdx, c->trackIdx,
//                                                                             boxes[c->layerIdx].crossPointRange.high));
//                     p->parent = node;
//                     for (int i=0; i<node->children.size(); i++)
//                         if (node->children[i] == c) node->children[i] = p;
//                     c->parent = p;
//                     p->children.push_back(c);
//                     p->isVio = true;
//                 }
//             }
//         };

//         // upwards
//         std::shared_ptr<db::GridSteiner> up, curr = vio;
//         while (curr->parent) {
//             up = curr->parent;
//             if (!boxes.count(up->layerIdx)) {
//                 box.layerIdx = up->layerIdx;
//                 boxes.insert({up->layerIdx, database.rangeSearch(box)});
//             }
//             if (boxes[up->layerIdx].crossPointRange.Contain(up->crossPointIdx)) {
//                 up->isVio = true;
//                 if (up->children.size() > 1)
//                     markBranches(up);
//             }
//             else {
//                 if (up->layerIdx == curr->layerIdx &&
//                     abs(up->crossPointIdx-curr->crossPointIdx) > 1 &&
//                     boxes[up->layerIdx].crossPointRange.StrictlyContain(curr->crossPointIdx)) { // long edge
//                     std::shared_ptr<db::GridSteiner> p;
//                     if (up->crossPointIdx < curr->crossPointIdx)
//                         p = std::make_shared<db::GridSteiner>(db::GridPoint(up->layerIdx, up->trackIdx,
//                                                                             boxes[up->layerIdx].crossPointRange.low));
//                     else
//                         p = std::make_shared<db::GridSteiner>(db::GridPoint(up->layerIdx, up->trackIdx,
//                                                                             boxes[up->layerIdx].crossPointRange.high));
//                     p->parent = up;
//                     for (int i=0; i<up->children.size(); i++)
//                         if (up->children[i] == curr) up->children[i] = p;
//                     curr->parent = p;
//                     p->children.push_back(curr);
//                     p->isVio = true;
//                 }
//                 break;
//             }
//             curr = up;
//         }

//         // downwards
//         markBranches(vio);
//     }
// }

// void guideRipup(db::Net &dbNet) {
//     vector<int> guideMark(dbNet.routeGuides.size()), ripupGuides;
//     vector<db::BoxOnLayer> gboxes = dbNet.routeGuides, ripupBoxes;
//     auto markGuide = [&](const utils::PointT<DBU> &vioP)->bool {
//         bool out = true;
//         for (size_t i=0; i<guideMark.size(); i++) {
//             if (gboxes[i].Contain(vioP)) {
//                 guideMark[i]++;
//                 out = false;
//             }
//         }
//         return out;
//     };
//     dbNet.ripupGuides.clear();
//     // expand
//     for (auto &g : gboxes) {
//         database.expandBox(g, db::rrrIterSetting.defaultGuideExpand);
//     }
//     // get guides cover vio nodes
//     for (auto vio : dbNet.vioNodes) {
//         markGuide(database.getLoc(*(vio->parent)));
//         bool outOfGuide = markGuide(database.getLoc(*vio));
//         if (outOfGuide) {
//             auto up = vio;
//             while (up->parent) {
//                 up = up->parent;
//                 if (up->isVio) break;
//                 up->isVio = true;
//                 if (!(markGuide(database.getLoc(*up))))
//                     break;
//             }
//             vector<db::GridSteiner *> downs, ddowns;
//             for (auto c : vio->children) {
//                 if (!c->isVio)
//                     downs.push_back(c.get());
//             }
//             while (!downs.empty()) {
//                 for (auto dp : downs) {
//                     if (dp->isVio) continue;
//                     dp->isVio = true;
//                     if (!(markGuide(database.getLoc(*dp))))
//                         continue;
//                     for (auto c : dp->children) {
//                         if (!c->isVio)
//                             ddowns.push_back(c.get());
//                     }
//                 }
//                 downs = std::move(ddowns);
//             }
//         }
//     }
//     for (int i=0; i<guideMark.size(); i++) {
//         if (guideMark[i] > 0) {
//             ripupGuides.push_back(i);
//             ripupBoxes.push_back(gboxes[i]);
//         }
//     }
//     // ripup nodes inside marked guides
//     dbNet.postOrderVisitGridTopo([&ripupBoxes](std::shared_ptr<db::GridSteiner> node) {
//         if (node->isVio) return;
//         const auto &nodeP = database.getLoc(*node);
//         for (auto &gbox : ripupBoxes) {
//             if (gbox.Contain(nodeP)) {
//                 node->isVio = true;
//                 break;
//             }
//         }
//     });
//     // check connectivity of marked guides
//     std::map<int, int> bosses;
//     vector<std::set<int>> unMergedConns;
//     size_t size = ripupGuides.size(), connNum = 0;
//     for (size_t i=0; i<size-1; i++) {
//         int first = ripupGuides[i];
//         for (size_t j=0; j<size; j++) {
//             int second = ripupGuides[j];
//             if (gboxes[first].HasIntersectWith(gboxes[second])) {
//                 if (bosses.count(first) && bosses.count(second)) {
//                     if (bosses[first] == bosses[second]) continue;
//                     int old = bosses[second], nbs = bosses[first];
//                     for (int id : unMergedConns[old]) {
//                         bosses[id] = nbs;
//                         unMergedConns[nbs].insert(id);
//                     }
//                     unMergedConns[old].clear();
//                 }
//                 else if (bosses.count(first)) {
//                     bosses.emplace(second, bosses[first]);
//                     unMergedConns[bosses[first]].insert(second);
//                 }
//                 else if (bosses.count(second)) {
//                     bosses.emplace(first, bosses[second]);
//                     unMergedConns[bosses[second]].insert(first);
//                 }
//                 else {
//                     bosses.emplace(first, connNum);
//                     bosses.emplace(second, connNum);
//                     unMergedConns.push_back({first, second});
//                     connNum++;
//                 }
//             }
//         }
//     }
//     for (size_t i=0; i<unMergedConns.size(); i++) {
//         if (!unMergedConns[i].empty()) {
//             dbNet.ripupGuides.emplace_back();
//             std::swap(dbNet.ripupGuides.back(), unMergedConns[i]);
//         }
//     }
//     for (int i : ripupGuides) {
//         if (!bosses.count(i))
//             dbNet.ripupGuides.push_back({i});
//     }
// }

// void PartialRipup::extractPseudoNets(vector<int> &nets) {
//     vector<int> rnets;
//     auto qc = [&](int id) {
//         auto &net = database.nets[nets[id]];
//         if (net.vioNodes.size() > 5) return;
//         quickFix(net);
//         for (auto r : net.gridTopo)
//             merge(r);
//         net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
//             if ((node->isVio) && (!(node->parent) || !(node->parent->isVio))) {
//                 net.pnets.push_back(node);
//             }
//         });
//     };
//     auto gr = [&](int id) {
//         auto &net = database.nets[nets[id]];
//         guideRipup(net);
//         net.postOrderVisitGridTopo([&net](std::shared_ptr<db::GridSteiner> node) {
//             if ((node->isVio) && (!(node->parent) || !(node->parent->isVio))) {
//                 net.pnets.push_back(node);
//             }
//         });
//     };

//     if (db::rrrIterSetting.quickFixMode)
//         runJobsMT(nets.size(), qc);
//     else
//         runJobsMT(nets.size(), gr);
//     for (int id : nets) {
//         if (database.nets[id].pnets.empty()) continue;
//         for (auto r : database.nets[id].pnets)
//             removeDbEdges(r, id);
//         rnets.push_back(id);
//     }
//     if (!rnets.empty())
//         nets = std::move(rnets);
// }