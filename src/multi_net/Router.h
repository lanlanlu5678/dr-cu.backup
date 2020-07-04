#pragma once

#include "db/Database.h"
#include "single_net/SingleNetRouter.h"

class Router {
public:
    void run();

private:
    int iter;
    vector<db::RouteStatus> allNetStatus;

    vector<int> getNetsToRoute();
    void ripup(const vector<int>& netsToRoute);
    void updateCost(const vector<int>& netsToRoute);
    void route(const vector<int>& netsToRoute);
    void finish();
    void unfinish();

    void printStat(bool major = false);

    // partial ripup
    int numOfPseudoNets = 0;
    void partialRipUp(const vector<int>& netsToRoute);
    // V0.1
    // void pureMark(std::shared_ptr<db::GridSteiner> &node, vector<int> &pnets);
    // void mark(std::shared_ptr<db::GridSteiner> &node, std::unordered_map<int,std::shared_ptr<db::GridSteiner>> &roots, int &num);
    // void purge(std::shared_ptr<db::GridSteiner> &node, vector<std::shared_ptr<db::GridSteiner>> &pnet, int netIdx);
    // void extractPseudoNet(db::Net &dbNet);
};
