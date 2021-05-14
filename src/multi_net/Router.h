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

    // PARTIAL RIPUP
    int numPNets;
    void partialRipup(vector<int> &netsToRoute);
    void greedyRoute();
    // void assignSeqRouter(const vector<SingleNetRouter> &routers, vector<SeqRouter> &seqRouters)
};
