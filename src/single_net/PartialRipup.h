#pragma once

#include <fstream>
#include "db/Database.h"

class PartialRipup {
    public:

    static void pureMark(std::shared_ptr<db::GridSteiner> node, vector<int> &pnets);
    static void mark(std::shared_ptr<db::GridSteiner> node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num);

    static void removeDbEdges(std::shared_ptr<db::GridSteiner> node, int netIdx);
    static void checkGuidesAndPins(std::shared_ptr<db::GridSteiner> node, vector<std::shared_ptr<db::GridSteiner>> &pins, vector<db::BoxOnLayer> &guides, vector<int> &vios);
    static void purge(std::shared_ptr<db::GridSteiner> node);
    static int inGuide(int l, DBU x, DBU y, db::BoxOnLayer &box);

    static void extractPseudoNets(db::Net &dbNet, int &num);

    static void mergeUnique(std::shared_ptr<db::GridSteiner> old, std::shared_ptr<db::GridSteiner> nld);

    static void checkPaths(std::shared_ptr<db::GridSteiner> root, vector<db::BoxOnLayer> &paths);
    static void checkPaths(vector<std::shared_ptr<db::GridSteiner>> &pins, vector<db::BoxOnLayer> &paths);

    static void plotPNet(std::ofstream &ofs, std::shared_ptr<db::GridSteiner> node);
};