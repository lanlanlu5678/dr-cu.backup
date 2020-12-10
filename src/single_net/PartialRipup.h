#pragma once

#include <fstream>
#include "db/Database.h"

class PartialRipup {
    public:
    static void preProcessRouteGuides();

    static void extractPseudoNets(vector<int> &nets);

    static void reOrderAndGuideMark(const vector<int> &nets);

    static void markOfGuides(const vector<int> &nets);
    static void expandRouteGuides(const vector<int> &nets, bool oriMode);

    static void reduceNets(vector<int> &netsToRoute);
};
