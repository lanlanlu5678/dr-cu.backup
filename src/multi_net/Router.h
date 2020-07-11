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
    int numOfPNets = 0;
    // vector<db::mutex_wrapper> netsLock;
    void partialRipUp(const vector<int>& netsToRoute);
};
