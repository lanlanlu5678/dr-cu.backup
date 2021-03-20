#include "Setting.h"
#include "Database.h"

namespace db {

void Setting::makeItSilent() {
    singleNetVerbose = VerboseLevelT::LOW;
    multiNetVerbose = VerboseLevelT::LOW;
    dbVerbose = VerboseLevelT::LOW;
}

void Setting::adapt() {
    if (database.nets.size() < 10000) {
        ++rrrIterLimit;
    }
    else if (database.nets.size() > 800000) {
        --rrrIterLimit;
    }
}

Setting setting;

void RrrIterSetting::update(int iter) {
    if (iter == 0) {
        defaultGuideExpand = setting.defaultGuideExpand;
        wrongWayPointDensity = setting.wrongWayPointDensity;
        addDiffLayerGuides = false;
    } else {
        defaultGuideExpand += iter * 2;
        wrongWayPointDensity = std::min(1.0, wrongWayPointDensity + 0.1);
        // if (database.nets.size() < 200000) {
        //     // high-effort mode (exclude million-net test case)
        //     addDiffLayerGuides = true;
        // }

        // PARTIAL RIPUP
        if (iter == 1) {
            fullyRoute = false;
            constrainInGuides = false;
            localRipup = true;
        }
        else if (iter == 2) {
            localRipup = false;
            adaptiveRipup = true;
        }
        // else if (iter == 3) {
        //     adaptiveRipup = false;
        //     guideRipup = true;
        // }
        // else if (iter == 4) {
        //     guideRipup = false;
        //     fullyRoute = true;
        // }
        else {
            fullyRoute = true;
        }
    }
    converMinAreaToOtherVio = ((iter + 1) < setting.rrrIterLimit);
}

void RrrIterSetting::print() const {
    printlog("defaultGuideExpand =", defaultGuideExpand);
    printlog("wrongWayPointDensity =", wrongWayPointDensity);
    printlog("addDiffLayerGuides =", addDiffLayerGuides);
}

RrrIterSetting rrrIterSetting;

}  // namespace db
