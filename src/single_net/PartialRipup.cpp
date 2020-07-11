#include "PartialRipup.h"

void PartialRipup::pureMark(std::shared_ptr<db::GridSteiner> node, vector<int> &pnets) {
    for (auto c : node->children) {
        if (c->isVio) pnets.push_back(c->pNetIdx);
        else if (c->distance < 3) {
            c->isVio = true;
            pureMark(c, pnets);
        }
    }
}

void PartialRipup::mark(std::shared_ptr<db::GridSteiner> node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num) {

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

void PartialRipup::purge(std::shared_ptr<db::GridSteiner> node, vector<std::shared_ptr<db::GridSteiner>> &pins, vector<int> &guides, vector<db::GridBoxOnLayer> &gridGuides) {

    // find related grid route guides
    if (node->parent && node->layerIdx != node->parent->layerIdx) {
        for (int i=0; i<gridGuides.size(); i++) {
            if (gridGuides[i].includePoint(*node)) {
                guides[i] = 1;
                break;
            }
        }
    }

    for (auto c : node->children) {
        if (c->isVio) purge(c, pins, guides, gridGuides);
        else {
            c->parent = nullptr;
            pins.push_back(c);
        }
    }
    if (node->children.empty()) {
        // if (node->pinIdx < 0) {
        //     std::cout << "*********** ERROR : IN PURGE INVALID END PIN\n";
        // }
        pins.push_back(node);
    }
    node->children.clear();
    node->parent = nullptr;
}

void PartialRipup::traversePNet(std::shared_ptr<db::GridSteiner> node, 
                                vector<std::shared_ptr<db::GridSteiner>> &pins,
                                vector<db::BoxOnLayer> &oriGuides,
                                vector<int> &guideIdx) {

    if (node->parent && node->layerIdx != node->parent->layerIdx) {
        db::BoxOnLayer box;
        auto p = database.getLoc(*node);
        for (int i=0; i<oriGuides.size(); i++) {
            box = oriGuides[i];
            if (box.layerIdx == node->layerIdx &&
                box.x.low < p.x && box.x.high > p.x &&
                box.y.low < p.y && box.y.high > p.y) {
                guideIdx[i] = 1;
                break;
            }
        }
    }

    for (auto c : node->children) {
        if (c->isVio) traversePNet(c, pins, oriGuides, guideIdx);
        else {
            c->parent = nullptr;
            pins.push_back(c);
        }
    }
    if (node->children.empty()) pins.push_back(node);
    node->children.clear();
    node->parent = nullptr;
}

void PartialRipup::removeDbEdges(std::shared_ptr<db::GridSteiner> node, int netIdx) {
    for (auto c : node->children) {
        database.removeEdge({*c,*node}, netIdx);
        if (c->extWireSeg) database.removeEdge(*(c->extWireSeg), netIdx);
        if (c->isVio) removeDbEdges(c, netIdx);
    }
}

void PartialRipup::extractPseudoNets(db::Net &dbNet, int &num) {

    // V0.3 mark nodes too near to any vio-nodes & remove all marked nodes.

    int n=0;
    std::unordered_map<int, std::shared_ptr<db::GridSteiner>> roots;

    dbNet.pinsOfPNets.clear();
    for (auto &node : dbNet.gridTopo) {
        mark(node, roots, n);
    }
    num += roots.size();
    for (auto &r : roots) {
        auto rp = r.second;
        if (rp->parent) {
            database.removeEdge({*rp, *(rp->parent)}, dbNet.idx);
            if (rp->extWireSeg) database.removeEdge(*(rp->extWireSeg), dbNet.idx);
            dbNet.pinsOfPNets.push_back(rp->parent);
        }
        else {
            dbNet.pinsOfPNets.push_back(rp);
        }
        removeDbEdges(rp, dbNet.idx);
    }
}

void PartialRipup::printPNet(std::shared_ptr<db::GridSteiner> node) {
    vector<std::shared_ptr<db::GridSteiner>> q;
    vector<std::shared_ptr<db::GridSteiner>> tq;
    vector<std::shared_ptr<db::GridSteiner>> all;
    bool loop = false;

    auto findin = [&all](std::shared_ptr<db::GridSteiner> node) ->bool {
        for (auto n : all) {
            if (n == node) return true;
        }
        return false;
    };

    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
    std::cout << (int)node->isVio << "\n";

    q.push_back(node);
    all.push_back(node);

    while (!q.empty()) {
        for (auto n : q) {
            for (auto c : n->children) {
                std::cout << (int)c->isVio << ",";
                if (findin(c)) loop = true;
                tq.push_back(c);
                all.push_back(c);
            }
            std::cout << ";    ";
        }
        if (loop) {
            std::cout << " A LOOP \n";
            return;
        }
        q = std::move(tq);
        std::cout << "\n";
    }
    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
}

// void PartialRipup::checkStartPin(std::shared_ptr<db::GridSteiner> sp) {
//     if (sp->pinIdx < 0) {
//         database.debugPrintLock.lock();
//         std::cout << "*************** ERROR : INVALID BREAK POINT IN ROOT\n";
//         database.debugPrintLock.unlock();
//     }
// }

// void PartialRipup::checkEndPin(int ppinIdx, int pNetIdx, vector<std::shared_ptr<db::GridSteiner>> &pins, const db::Net &dbNet) {
//     if (pins[ppinIdx]->pinIdx == -1) {
//         database.debugPrintLock.lock();
//         log() << "ERROR : net " << dbNet.idx << " pseudo net " << pNetIdx 
//             << " pin " << ppinIdx << " is not real pin but without children" << std::endl;
//         log() << "ERROR : pseudo pins number : " << pins.size() << "\n";
//         log() << "ERROR : pseudo pin is vio  : " << pins[ppinIdx]->isVio << "\n";
//         log() << "ERROR : pseudo nets number : " << dbNet.pinsOfPNets.size() << "\n";
//         if (pins[ppinIdx] == database.breakpoint) log() << "is that breakpoint\n";
//         database.debugPrintLock.unlock();
//     }
// }