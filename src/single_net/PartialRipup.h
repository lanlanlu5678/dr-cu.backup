#pragma once

#include "db/Database.h"

class PartialRipup {
    public:

    static void pureMark(std::shared_ptr<db::GridSteiner> node, vector<int> &pnets);
    static void mark(std::shared_ptr<db::GridSteiner> node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num);
    // static void purge(std::shared_ptr<db::GridSteiner> &node, vector<std::shared_ptr<db::GridSteiner>> &pins, int netIdx);

    static void removeDbEdges(std::shared_ptr<db::GridSteiner> node, int netIdx);
    // static void purge(std::shared_ptr<db::GridSteiner> node, vector<std::shared_ptr<db::GridSteiner>> &pins);
    static void purge(std::shared_ptr<db::GridSteiner> node, vector<std::shared_ptr<db::GridSteiner>> &pins, vector<int> &guides, vector<db::GridBoxOnLayer> &gridGuides);
    static void traversePNet(std::shared_ptr<db::GridSteiner> node, vector<std::shared_ptr<db::GridSteiner>> &pins, vector<db::BoxOnLayer> &oriGuides, vector<int> &guideIdx);

    static void extractPseudoNets(db::Net &dbNet, int &num);

    static void printPNet(std::shared_ptr<db::GridSteiner> node);
    static void checkStartPin(std::shared_ptr<db::GridSteiner> sp);
    static void checkEndPin(int ppinIdx, int pNetIdx, vector<std::shared_ptr<db::GridSteiner>> &pins, const db::Net &dbNet);
};