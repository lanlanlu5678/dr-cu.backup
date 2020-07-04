#include "PartialRipup.h"

// void PartialRipup::checkRouteGuides(vector<rgNodes> &rgGraph, db::Net &dbNet) {
//     int currGuideIdx, nextGuideIdx, currLayerIdx;
//     vector<std::shared_ptr<db::GridSteiner>> queue1;
//     vector<std::shared_ptr<db::GridSteiner>> queue2;

//     rgGraph.clear();
//     rgGraph.resize(dbNet.routeGuides.size());
    
//     for (auto node : dbNet.gridTopo) {
//         currGuideIdx = contains(node, dbNet.routeGuides);
//         currLayerIdx = node->layerIdx;
//         if (node->pinIdx == 0) {
//             rgGraph[currGuideIdx].parent = -1;
//             rgGraph[currGuideIdx].enterPoint = node;
//         }
//         else rgGraph[currGuideIdx].breakPinIdx = node->pinIdx;

//         while (true) {
//             for (auto &c : node->children) {
//                 if (c->layerIdx != currLayerIdx) {
//                     nextGuideIdx = contains(c, dbNet.routeGuides);
//                     // ignore nodes from/to expanded & diff-layer guides
//                     if (nextGuideIdx > -1 && nextGuideIdx != currGuideIdx) {
//                         rgGraph[currGuideIdx].children.push_back(nextGuideIdx);
//                         rgGraph[currGuideIdx].leavePoints.push_back(c);
//                         queue1.push_back(c);
//                         continue;
//                     }
//                 }
//                 queue2.push_back(c);
//             }
//             if (!queue2.empty()) {
//                 node = queue2.back();
//                 queue2.pop_back();
//             }
//             else if (!queue1.empty()) {
//                 node = queue1.back();
//                 currGuideIdx = contains(node, dbNet.routeGuides);
//                 currLayerIdx = node->layerIdx;
//                 queue1.pop_back();
//             }
//             else break;
//         }
//     }
// }

int PartialRipup::contains(std::shared_ptr<db::GridSteiner> &node, vector<db::BoxOnLayer> &guides) {
    auto loc = database.getLoc(*(node));
    for (int i=0; i<guides.size(); i++) {
        auto &box = guides[i];
        if (node->layerIdx == box.layerIdx) {
            if ((loc.x >= box.x.low && loc.x <= box.x.high) && (loc.y >= box.y.low && loc.y <= box.y.high))
                return i;
        }
    }
    return -1;
}

void PartialRipup::pureMark(std::shared_ptr<db::GridSteiner> &node, vector<int> &pnets) {
    for (auto &c : node->children) {
        if (c->isVio) pnets.push_back(c->pNetIdx);
        else if (c->distance < 3) {
            c->isVio = true;
            pureMark(c, pnets);
        }
    }
}

void PartialRipup::mark(std::shared_ptr<db::GridSteiner> &node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num) {

    for (auto &c : node->children) mark(c, roots, num);

    int d = 4, d1 = 4, d2 = 4, danger = -1; // d1,d2 is the 2 smallest distance of children
    vector<int> adjPnets; // adjacent pseudo nets that need to be merge

    if (node->isVio) {

        // current node is a terminal of vio edges, mark all children with d<3
        // as propagating, every marked child is another terminal
        pureMark(node, adjPnets);

        // if any, delete the roots of adjacent pseudo nets, merge to current node
        for (int i : adjPnets) roots.erase(i);
        roots[num] = node;
        node->pNetIdx = num;
        num++;

        // update distance
        node->distance = -1;
    }
    else {
        for (int i=0; i<node->children.size(); i++) {
            d = node->children[i]->distance;
            if (d < 3) {
                if (d < d1) {
                    d2 = d1;
                    d1 = d;
                    danger = i;
                }
                else if (d < d2) d2 = d;
            }
        }
        
        if (d1 + d2 < 2) {

            // current node is a terminal of vio edges, mark all children with d<3
            // as propagating, ever child has been marked is another terminal
            pureMark(node, adjPnets);

            // merge pseudo nets in subtree to current node
            for (int i : adjPnets) roots.erase(i);
            roots[num] = node;
            node->pNetIdx = num;
            num++;

            node->isVio = true;
            node->distance = -1;
        }
        else if (d1 < 3) {

            // has dangerous children, but no need to mark
            // record the distance to the nearest vio edge
            node->distance = d1 + 1;
        }
        // else a pseudo net is completely found
    }
}

void PartialRipup::purge(std::shared_ptr<db::GridSteiner> &node, vector<std::shared_ptr<db::GridSteiner>> &pins, int netIdx) {

    // every marked node except pins will be purged
    // un-marked children & real pins are end pins
    for (auto &c : node->children) {
        database.removeEdge({*c, *node}, netIdx);
        if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), netIdx);

        if (c->isVio) purge(c, pins, netIdx);
        else {
            c->parent = nullptr;
            pins.push_back(c);
        }
    }

    if (node->children.empty()) { // although vio, but an end pin
        // node->isVio = false;
        // node->distance = 4;
        // node->pNetIdx = -1;
        // node->extWireSeg = nullptr;
        pins.push_back(node);
    }

    node->children.clear();
    node->parent = nullptr;
}

void PartialRipup::extractPseudoNets(db::Net &dbNet, int &num, int &ignore) {
    
    std::unordered_map<int, std::shared_ptr<db::GridSteiner>> rootsOfPnets;
    int n=0, guideIdx=-1;

    dbNet.pinsOfPseudoNets.clear();

    // mark originally not violated edges to be violated, when they are near violated edges
    // after marking, each pseudo net is a set of consecutive violated edges, and is at least 4 edges away from another
    // real-time update the root of pseudo nets
    for (auto &node : dbNet.gridTopo) {
        mark(node, rootsOfPnets, n);
    }
    
    // boundaries nodes are not-vio nodes & real pins
    // position of nodes remains, but GridSteiner will be replaced
    num += rootsOfPnets.size();
    for (auto &p : rootsOfPnets) {
        auto &node = p.second;

        // // if start pin out of guide, ignore
        // guideIdx = contains(node, dbNet.routeGuides);
        // if (guideIdx > 0) {
            dbNet.pinsOfPseudoNets.emplace_back();
            auto &pins = dbNet.pinsOfPseudoNets.back();

            // TODO: extend reroute space

            if (node->parent) {
                auto &root = node->parent;
                database.removeEdge({*node, *root}, dbNet.idx);
                if (root->extWireSeg) database.removeEdge(*(root->extWireSeg),dbNet.idx);
                pins.push_back(root);
                // real pin access boxes may be different, so remove a few more nodes,
                // when commit reroute result only need to replace the start pin
                // ***********HOWEVER : the mis-match happens between the start pin and its root
                // for (auto &c : root->children) {
                //     if (!c->isVio) {
                //         database.removeEdge({*c, *root}, dbNet.idx);
                //         if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), dbNet.idx);
                //         pins.push_back(c);
                //         c->parent = nullptr;
                //     }
                // }
                // root->children.clear();
            }
            else pins.push_back(node);

            purge(node, pins, dbNet.idx);
        // }
        // else ignore++;
    }
}




// void PartialRipup::extractPseudoNet(db::Net &dbNet) {
//     std::unordered_map<int, std::shared_ptr<db::GridSteiner>> rootsOfPnets;
//     vector<std::shared_ptr<db::GridSteiner>> tempChildren;
//     int num = 0;

//     dbNet.pinsOfPseudoNets.clear();

//     for (auto &root : dbNet.gridTopo) {
//         // mark originally not violated edges to be violated, when they are near violated edges
//         // after marking, each pseudo net is a set of consecutive violated edges, and is at least 4 edges away from another
//         // real-time update the root of pseudo nets
//         mark(root, rootsOfPnets, num);
//     }

//     numOfPseudoNets += rootsOfPnets.size();
//     for (auto &p : rootsOfPnets) {
//         auto root = p.second;
//         dbNet.pinsOfPseudoNets.emplace_back();
//         auto &pnet = dbNet.pinsOfPseudoNets.back();

//         // root of pseudo net should reset all status of partial ripup
//         if (root->parent) { // one step beyond, get a safe node
//             root->isVio = true;
//             root = root->parent;
//         }
//         else { // reset the origin root
//             root->isVio = false;
//             root->distance = 4;
//             root->pNetIdx = -1;
//             if (root->extWireSeg) { // origin root may have a vio extseg
//                 database.removeEdge(*(root->extWireSeg), dbNet.idx);
//                 root->extWireSeg = nullptr;
//             }
//         }
//         pnet.push_back(root);

//         tempChildren = root->children;
//         root->children.clear();

//         for (auto &c : tempChildren) {
//             if (c->isVio) {
//                 // purge() will check all marked nodes and remove any edge/extSeg attached
//                 // clear all the shared_ptr a marked node stores (parents & children)
//                 // end pins' extSeg will be removed
//                 purge(c, pnet, dbNet.idx);
//                 database.removeEdge({*c, *root}, dbNet.idx);
//                 if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), dbNet.idx);
//             }
//             else root->children.push_back(c); // vio children will be clear
//         }

//         if (pnet.size() < 2) {
//             log() << "pseudo net " << p.first << " of dbNet " << dbNet.idx << " has only one pin\n";
//             for (auto &c : root->children) {
//                 database.removeEdge({*c, *root}, dbNet.idx);
//                 if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), dbNet.idx);
//                 c->extWireSeg = nullptr;
//                 c->parent = nullptr;
//                 pnet.push_back(c);
//             }
//         }
//     }

// }

// void PartialRipup::partialRipUp(const vector<int>& netsToRoute) {

//     numOfPseudoNets = 0;

//     for (int netIdx : netsToRoute) {
        
//         extractPseudoNet(database.nets[netIdx]);

//         // if (db::setting.multiNetVerbose >= +db::VerboseLevelT::HIGH) {
//         //     log() << "After extract net :" << netIdx << " , " << numOfPseudoNets << " pNets found.\n";
//         // }

//         allNetStatus[netIdx] = db::RouteStatus::FAIL_UNPROCESSED;
//     }

//     if (db::setting.multiNetVerbose >= +db::VerboseLevelT::HIGH)
//         log() << "All pNets are extracted.\n";
// }

// void PartialRipup::extractPseudoNet(db::Net &dbNet, int &num) {

//     std::vector<std::shared_ptr<db::GridSteiner>> oQueue;
//     std::vector<std::shared_ptr<db::GridSteiner>> iQueue;
//     auto &pins = dbNet.pinsOfPseudoNets;
//     auto &guides = dbNet.guidesOfPseudoNets;
//     bool inVio = false;
//     int guideIdx, n=0;

//     /*  Mark boundaries of pNet
//         assume paths in different route guides are all connected by vias
//         pins of pNets are via nodes or real pins
//         if a root in gridTopo is in a vio guide of previous pNet
//         still breaks, get all pin access boxes of real pins!
//     */

//     oQueue = dbNet.gridTopo;
//     auto node = oQueue.back();
//     oQueue.pop_back();

//     while (true)
//     {
//         // check route guide
//         if (!node->parent || node->viaType) {
//             guideIdx = contains(node, dbNet.routeGuides);
//             if (guideIdx > 0 && dbNet.vioGuides[guideIdx] > 0) {
//                 if (inVio) guides[n].push_back(guideIdx); // record the consecutive vio guides
//                 else {
//                     // find a new pNet
//                     inVio = true;
//                     pins.emplace_back();
//                     guides.emplace_back();
//                     n++;
//                     pins[n].push_back(node);
//                     guides[n].push_back(guideIdx);
//                 }
//             }
//             if (guideIdx > 0 && dbNet.vioGuides[guideIdx] < 1 && inVio) {
//                 // find an pseudo pin
//                 inVio = false;
//                 pins[n].push_back(node);
//             }
//         }

//         if (inVio) {
//             if (node->parent) database.removeEdge({*node, *(node->parent)}, dbNet.idx);
//             if (node->extWireSeg) database.removeEdge(*(node->extWireSeg), dbNet.idx);

//             if (node->children.empty()) {
//                 // a real pin, wouldn't be the start pin
//                 pins[n].push_back(node);
//                 if (iQueue.empty()) {
//                     inVio = false;
//                     node = oQueue.back();
//                     oQueue.pop_back();
//                 }
//             } 
//             for (auto &c : node->children) iQueue.push_back(c);
            
//             node->parent = nullptr;
//             node->children.clear();

//             node = iQueue.back();
//             iQueue.pop_back();
//         }
//         else {
//             if (iQueue.empty()) {
//                 if (node->children.empty()) {
//                     if (oQueue.empty()) break;
//                     node = oQueue.back();
//                     oQueue.pop_back();
//                 }
//                 for (int i=1; i<node->children.size(); i++)
//                     oQueue.push_back(node->children[i]);
//                 node = node->children[0];
//             }
//             else {
//                 inVio = true;
//                 for (auto &c : node->children)
//                     oQueue.push_back(c);
//                 node = iQueue.back();
//                 iQueue.pop_back();
//             }
//         }
//     }
    
//     num += ++n;
// }