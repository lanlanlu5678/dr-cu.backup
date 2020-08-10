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

void PartialRipup::checkGuidesAndPins(std::shared_ptr<db::GridSteiner> node,
                            vector<std::shared_ptr<db::GridSteiner>> &pins,
                            vector<db::BoxOnLayer> &guides,
                            vector<int> &vios) {

    for (int i=0; i<guides.size(); i++) {
        auto p = database.getLoc(*node);
        int d = inGuide(node->layerIdx, p.x, p.y, guides[i]);
        if (d > -3) vios[i] |= (1<<(d+2));
    }

    for (auto c : node->children) {
        if (c->isVio) checkGuidesAndPins(c, pins, guides, vios);
        else pins.push_back(c);
    }
    if (node->children.empty()) pins.push_back(node);
}

void PartialRipup::purge(std::shared_ptr<db::GridSteiner> node) {
    for (auto c : node->children) {
        if (c->isVio) purge(c);
        else c->parent = nullptr;
    }
    node->children.clear();
    node->parent = nullptr;
}

int PartialRipup::inGuide(int l, DBU x, DBU y, db::BoxOnLayer &guide) {
    int d = l - guide.layerIdx;
    if (abs(d)<db::rrrIterSetting.diffLayerBound &&
        guide.x.low < x && guide.x.high > x &&
        guide.y.low < y && guide.y.high > y) return d;
    
    else return -3;
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

void PartialRipup::mergeUnique(std::shared_ptr<db::GridSteiner> old, std::shared_ptr<db::GridSteiner> nld) {
    bool notdul;
    for (auto c : old->children) {
        if (!c->isVio) {
            notdul = true;
            for (auto nc : nld->children) {
                if (c->layerIdx == nc->layerIdx &&
                    c->trackIdx == nc->trackIdx &&
                    c->crossPointIdx == nc->crossPointIdx) {
                    mergeUnique(c, nc);
                    c->parent = nullptr;
                    notdul = false;
                }
            }
            if (notdul) {
                c->parent = nld;
                nld->children.push_back(c);
            }
        }
    }
    old->children.clear();
}

void PartialRipup::checkPaths(std::shared_ptr<db::GridSteiner> pin0, vector<db::BoxOnLayer> &paths) {
    vector<std::shared_ptr<db::GridSteiner>> newPaths;
    vector<std::shared_ptr<db::GridSteiner>> samePath;
    utils::PointT<DBU> cp, p;
    db::BoxOnLayer box;
    std::shared_ptr<db::GridSteiner> node, next;

    if (pin0->isVio) newPaths.push_back(pin0);
    else
        for (auto c : pin0->children)
            if (c->isVio) {
                newPaths.push_back(c);
                // if (c->viaType) {
                //     p = database.getLoc(*(pin0));
                //     paths.emplace_back(pin0->layerIdx, p);
                // }
                p = database.getLoc(*(pin0));
                box.Set(pin0->layerIdx, p);
                cp = database.getLoc(*c);
                box.FastUpdate(cp);
                paths.push_back(box);
            }

    while (!newPaths.empty()) {
        node = newPaths.back();
        newPaths.pop_back();

        p = database.getLoc(*node);
        box.Set(node->layerIdx, p);
        
        while (node != nullptr) {
            p = database.getLoc(*node);
            box.FastUpdate(p);
            for (auto c : node->children) {
                if (c->isVio) {
                    if (c->viaType)
                        newPaths.push_back(c);
                    else
                        samePath.push_back(c);
                }
                else {
                    cp = database.getLoc(*c);
                    if (c->viaType)
                        paths.emplace_back(c->layerIdx, cp);
                    else
                        box.FastUpdate(cp);
                }
            }
            if (samePath.empty())
                node = nullptr;
            else {
                node = samePath.back();
                samePath.pop_back();
            }
        }
        paths.push_back(box);
    }
}

void PartialRipup::checkPaths(vector<std::shared_ptr<db::GridSteiner>> &pins, vector<db::BoxOnLayer> &paths) {
    vector<std::shared_ptr<db::GridSteiner>> newPaths;
    vector<std::shared_ptr<db::GridSteiner>> samePath;
    utils::PointT<DBU> cp, p;
    db::BoxOnLayer box;
    std::shared_ptr<db::GridSteiner> node, next;

    if (pins[0]->isVio) newPaths.push_back(pins[0]);
    else
        for (auto c : pins[0]->children)
            if (c->isVio) {
                newPaths.push_back(c);
                // if (c->viaType) {
                //     p = database.getLoc(*(pins[0]));
                //     paths.emplace_back(pins[0]->layerIdx, p);
                // }
                p = database.getLoc(*(pins[0]));
                box.Set(pins[0]->layerIdx, p);
                cp = database.getLoc(*c);
                box.FastUpdate(cp);
                paths.push_back(box);
            }

    while (!newPaths.empty()) {
        node = newPaths.back();
        newPaths.pop_back();

        p = database.getLoc(*node);
        box.Set(node->layerIdx, p);
        
        while (node != nullptr) {
            p = database.getLoc(*node);
            box.FastUpdate(p);
            if (node->children.empty())
                pins.push_back(node);
            for (auto c : node->children) {
                if (c->isVio) {
                    if (c->viaType)
                        newPaths.push_back(c);
                    else
                        samePath.push_back(c);
                }
                else {
                    pins.push_back(c);
                    cp = database.getLoc(*c);
                    if (c->viaType)
                        paths.emplace_back(c->layerIdx, cp);
                    else
                        box.FastUpdate(cp);
                }
            }
            if (samePath.empty())
                node = nullptr;
            else {
                node = samePath.back();
                samePath.pop_back();
            }
        }
        paths.push_back(box);
    }
}

void PartialRipup::plotPNet(std::ofstream &ofs, std::shared_ptr<db::GridSteiner> node) {
    utils::PointT<DBU> p;
    db::BoxOnLayer box;
    for (auto c : node->children) {
        if (c->isVio) plotPNet(ofs, c);
        else {
            p = database.getLoc(*node);
            box.Set(node->layerIdx, p);
            ofs << box << std::endl;
            p = database.getLoc(*c);
            box.Set(c->layerIdx, p);
            ofs << box << std::endl;
        }
    }
    if (node->parent) {
        p = database.getLoc(*(node->parent));
        box.Set(node->parent->layerIdx, p);
    }
    else {
        p = database.getLoc(*node);
        box.Set(node->layerIdx, p);
    }
    ofs << box << std::endl;
    p = database.getLoc(*node);
    box.Set(node->layerIdx, p);
    ofs << box << std::endl;
}