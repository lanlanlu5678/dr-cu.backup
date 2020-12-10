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
        // defaultGuideExpand += iter;
        wrongWayPointDensity = std::min(1.0, wrongWayPointDensity + 0.1);

        // PARTIAL RIPUP
        constrainInGuide = false;
        // if (iter < 3) {
        //     quickFixMode = true;
        //     extraCPs += 5;
        //     extraTracks += 2;
        //     oriMode = false;
        // }
        // else if (iter == 3) {
        //     quickFixMode = false;
        //     guideRipup = true;
        //     addDiffLayerGuides = true;
        //     // setting.numThreads = 0;
        // }
        // else {
        //     // guideRipup = false;
        //     // oriMode = true;
        //     defaultGuideExpand -= iter * 2;
        // }
        addDiffLayerGuides = true;
        oriMode = false;
        guideRipup = true;
        defaultGuideExpand = 14;
    }
    converMinAreaToOtherVio = ((iter + 1) < setting.rrrIterLimit);
}

void RrrIterSetting::print() const {
    printlog("defaultGuideExpand =", defaultGuideExpand);
    printlog("wrongWayPointDensity =", wrongWayPointDensity);
    printlog("addDiffLayerGuides =", addDiffLayerGuides);
    if (quickFixMode) {
        printlog("extraTracks = ", extraTracks);
        printlog("extraCrossPoints = ", extraCPs);
    }
}

RrrIterSetting rrrIterSetting;

}  // namespace db
