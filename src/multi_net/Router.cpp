#include "Router.h"
#include "Scheduler.h"
#include "single_net/PartialRipup.h"
#include "single_net/PostProcess.h"

const MTStat& MTStat::operator+=(const MTStat& rhs) {
    auto dur = rhs.durations;
    std::sort(dur.begin(), dur.end());
    if (durations.size() < dur.size()) {
        durations.resize(dur.size(), 0.0);
    }
    for (int i = 0; i < dur.size(); ++i) {
        durations[i] += dur[i];
    }
    return *this;
}

ostream& operator<<(ostream& os, const MTStat mtStat) {
    double minDur = std::numeric_limits<double>::max(), maxDur = 0.0, avgDur = 0.0;
    for (double dur : mtStat.durations) {
        minDur = min(minDur, dur);
        maxDur = max(maxDur, dur);
        avgDur += dur;
    }
    avgDur /= mtStat.durations.size();
    os << "#threads=" << mtStat.durations.size() << " (dur: min=" << minDur << ", max=" << maxDur << ", avg=" << avgDur
       << ")";
    return os;
}

inline void printnodes(std::shared_ptr<db::GridSteiner> node) {
    // const auto &np = database.getLoc(*node);
    std::cout << *node << "; ";
    if (node->parent) std::cout << *(node->parent);
    std::cout << std::endl;
    for (auto c : node->children) {
        // if (c->isVio)
            printnodes(c);
    }
}

void Router::run() {
    allNetStatus.resize(database.nets.size(), db::RouteStatus::FAIL_UNPROCESSED);
    for (iter = 0; iter < db::setting.rrrIterLimit; iter++) {
        log() << std::endl;
        log() << "################################################################" << std::endl;
        log() << "Start RRR iteration " << iter << std::endl;
        log() << std::endl;
        db::routeStat.clear();
        // PARTIAL RIPUP
        db::rrrIterSetting.update(iter);
        vector<int> netsToRoute = getNetsToRoute();
        if (netsToRoute.empty()) {
            if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
                log() << "No net is identified for this iteration of RRR." << std::endl;
                log() << std::endl;
            }
            break;
        }
        if (iter > 0) {
            // updateCost should before ripup, otherwise, violated nets have gone
            updateCost(netsToRoute);
            if (db::rrrIterSetting.fullyRoute)
                ripup(netsToRoute);
            else
                partialRipup(netsToRoute);
        }
        database.statHistCost();
        if (db::setting.rrrIterLimit > 1) {
            double step = (1.0 - db::setting.rrrInitVioCostDiscount) / (db::setting.rrrIterLimit - 1);
            database.setUnitVioCost(db::setting.rrrInitVioCostDiscount + step * iter);
        }
        if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
            db::rrrIterSetting.print();
        }
        route(netsToRoute);
        log() << std::endl;
        log() << "Finish RRR iteration " << iter << std::endl;
        log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
              << std::endl;
        if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
            printStat(db::setting.rrrWriteEachIter);
        }
        // if (db::setting.rrrWriteEachIter) {
        //     std::string fn = "iter" + std::to_string(iter) + "_" + db::setting.outputFile;
        //     printlog("Write result of RRR iter", iter, "to", fn, "...");
        //     finish();
        //     database.writeDEF(fn);
        //     unfinish();
        // }
    }
    // greedyRoute();
    finish();
    log() << std::endl;
    log() << "################################################################" << std::endl;
    log() << "Finish all RRR iterations and PostRoute" << std::endl;
    log() << "MEM: cur=" << utils::mem_use::get_current() << "MB, peak=" << utils::mem_use::get_peak() << "MB"
          << std::endl;
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printStat(true);
    }
}

vector<int> Router::getNetsToRoute() {
    vector<int> netsToRoute;
    if (iter == 0) {
        for (int i = 0; i < database.nets.size(); i++) {
            netsToRoute.push_back(i);
        }
    } else {
        for (auto& net : database.nets) {
            if (UpdateDB::checkViolation(net)) {
                netsToRoute.push_back(net.idx);
            }
        }
    }

    return netsToRoute;
}

void Router::ripup(const vector<int>& netsToRoute) {
    for (auto netIdx : netsToRoute) {
        // UpdateDB::clearRouteResult(database.nets[netIdx]);
        allNetStatus[netIdx] = db::RouteStatus::FAIL_UNPROCESSED;
    }
}

// PARTIAL RIPUP
void Router::partialRipup(vector<int> &nets) {
    PartialRipup::extractPNets(nets);
    numPNets = 0;
    for (int i : nets) {
        numPNets += database.nets[i].pnets.size();
        // allNetStatus[i] = db::RouteStatus::FAIL_UNPROCESSED;
    }
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE)
        log() << "From " << nets.size() << " nets, " << numPNets << " PNets found" << std::endl;
}

void Router::updateCost(const vector<int>& netsToRoute) {
    database.addHistCost();
    database.fadeHistCost(netsToRoute);
}

void Router::route(const vector<int>& netsToRoute) {
    vector<double> stages(6, 0.0);
    utils::timer t;
    // init SingleNetRouters
    t.start();
    vector<SingleNetRouter> routers;
    if (db::rrrIterSetting.fullyRoute) numPNets = netsToRoute.size();
    routers.reserve(numPNets);
    for (int i : netsToRoute) {
        if (database.nets[i].pnets.empty())
            routers.emplace_back(database.nets[i]);
        else {
            for (int pnet=0; pnet<int(database.nets[i].pnets.size()); pnet++) {
                routers.emplace_back(database.nets[i]);
                routers.back().localNet.pnetIdx = pnet;
            }
        }
    }
    stages[0] = t.elapsed();

    // pre route
    t.start();
    auto preMT = runJobsMT(numPNets, [&](int netIdx) { routers[netIdx].preRoute(); });
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("preMT", preMT);
        printStat();
    }
    stages[1] = t.elapsed();

    // schedule
    t.start();
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "Start multi-thread scheduling. There are " << numPNets << " nets to route." << std::endl;
    }
    Scheduler scheduler(routers);
    const vector<vector<int>>& batches = scheduler.schedule();
    stages[2] = t.elapsed();
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        log() << "Finish multi-thread scheduling" << ((db::setting.numThreads == 0) ? " using simple mode" : "")
              << ". There will be " << batches.size() << " batches." << std::endl;
        log() << std::endl;
    }

    // maze route and commit DB by batch
    int iBatch = 0;
    MTStat allMazeMT, allCommitMT, allGetViaTypesMT, allCommitViaTypesMT;
    for (const vector<int>& batch : batches) {
        // 1 maze route
        t.start();
        auto mazeMT = runJobsMT(batch.size(), [&](int jobIdx) {
            auto& router = routers[batch[jobIdx]];
            router.mazeRoute();
            allNetStatus[router.dbNet.idx] = router.status;
        });
        stages[3] += t.elapsed();
        allMazeMT += mazeMT;
        // 2 commit nets to DB
        t.start();
        auto commitMT = runJobsMT(batch.size(), [&](int jobIdx) {
            auto& router = routers[batch[jobIdx]];
            if (!db::isSucc(router.status)) return;
            router.commitNetToDB();
            router.localNet.pnetPins.clear();
        });
        allCommitMT += commitMT;
        stages[4] += t.elapsed();
        // 3 get via types
        t.start();
        allGetViaTypesMT += runJobsMT(batch.size(), [&](int jobIdx) {
            auto& router = routers[batch[jobIdx]];
            if (!db::isSucc(router.status)) return;
            PostRoute postRoute(router.dbNet);
            postRoute.getViaTypes();
        });
        allCommitViaTypesMT += runJobsMT(batch.size(), [&](int jobIdx) {
            auto& router = routers[batch[jobIdx]];
            if (!db::isSucc(router.status)) return;
            UpdateDB::commitViaTypes(router.dbNet);
        });
        stages[5] += t.elapsed();
        // 4 stat
        if (db::setting.multiNetVerbose >= +db::VerboseLevelT::HIGH && db::setting.numThreads != 0) {
            int maxNumVertices = 0;
            for (int i : batch) maxNumVertices = std::max(maxNumVertices, routers[i].localNet.estimatedNumOfVertices);
            log() << "Batch " << iBatch << " done: size=" << batch.size() << ", mazeMT " << mazeMT << ", commitMT "
                  << commitMT << ", peakM=" << utils::mem_use::get_peak() << ", maxV=" << maxNumVertices << std::endl;
        }
        iBatch++;
    }
    if (iter > 0) {
        runJobsMT(netsToRoute.size(), [&](int id) {
            auto &net = database.nets[netsToRoute[id]];
            PostMazeRoute(net).run();
            PostProcess::removeSameNetVioVias(net);
            // for (auto root : net.gridTopo)
            //     PostProcess::removeCorners(root, net.idx);
        });
    }
    else {
        runJobsMT(netsToRoute.size(), [&](int id) {
            auto &net = database.nets[netsToRoute[id]];
            PartialRipup::mergeRouteGuides(net);
            net.routeGuideVios.resize(net.routeGuides.size(), 0);
            net.gridRouteGuides.clear();
            for (const auto &g : net.routeGuides)
                net.gridRouteGuides.push_back(database.rangeSearch(g));
        });
    }
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("allMazeMT", allMazeMT);
        printlog("allCommitMT", allCommitMT);
        printlog("allGetViaTypesMT", allGetViaTypesMT);
        printlog("allCommitViaTypesMT", allCommitViaTypesMT);
        printf(" init     | preRoute    | schedule   | mazeRoute   | commit    | viaType\n");
        double total = 0.0;
        for (double time : stages) {
            total += time;
            printf(" %f   ", time);
        }
        printf("\n");
        for (double time : stages) {
            printf(" %f   ", time / total * 100);
        }
        printf(" %% \n");
    }
}

void Router::finish() {
    PostScheduler postScheduler(database.nets);
    const vector<vector<int>>& batches = postScheduler.schedule();
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("There will be", batches.size(), "batches for getting via types.");
    }
    // 1. redo min area handling
    MTStat allPostMaze2MT;
    for (const vector<int>& batch : batches) {
        runJobsMT(batch.size(), [&](int jobIdx) {
            int netIdx = batch[jobIdx];
            if (!db::isSucc(allNetStatus[netIdx])) return;
            UpdateDB::clearMinAreaRouteResult(database.nets[netIdx]);
        });
        allPostMaze2MT += runJobsMT(batch.size(), [&](int jobIdx) {
            int netIdx = batch[jobIdx];
            if (!db::isSucc(allNetStatus[netIdx])) return;
            PostMazeRoute(database.nets[netIdx]).run2();
        });
        runJobsMT(batch.size(), [&](int jobIdx) {
            int netIdx = batch[jobIdx];
            if (!db::isSucc(allNetStatus[netIdx])) return;
            UpdateDB::commitMinAreaRouteResult(database.nets[netIdx]);
        });
    }
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("allPostMaze2MT", allPostMaze2MT);
    }
    // 1.5 post processing
    auto ppMt = runJobsMT(database.nets.size(), [&](int netIdx) {
        PostProcess::removeSameNetVioVias(database.nets[netIdx]);
        PostProcess::fixMAR(database.nets[netIdx]);
    });

    // 2. get via types again
    for (int iter = 0; iter < db::setting.multiNetSelectViaTypesIter; iter++) {
        MTStat allGetViaTypesMT, allCommitViaTypesMT;
        for (const vector<int>& batch : batches) {
            allGetViaTypesMT += runJobsMT(batch.size(), [&](int jobIdx) {
                int netIdx = batch[jobIdx];
                if (!db::isSucc(allNetStatus[netIdx])) return;
                PostRoute postRoute(database.nets[netIdx]);
                if (iter == 0) postRoute.considerViaViaVio = false;
                postRoute.getViaTypes();
            });
            allCommitViaTypesMT += runJobsMT(batch.size(), [&](int jobIdx) {
                int netIdx = batch[jobIdx];
                if (!db::isSucc(allNetStatus[netIdx])) return;
                UpdateDB::commitViaTypes(database.nets[netIdx]);
            });
        }
        if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
            printlog("allGetViaTypesMT", allGetViaTypesMT);
            printlog("allCommitViaTypesMT", allCommitViaTypesMT);
        }
    }
    // 3. post route
    auto postMT = runJobsMT(database.nets.size(), [&](int netIdx) {
        if (!db::isSucc(allNetStatus[netIdx])) return;
        PostRoute postRoute(database.nets[netIdx]);
        postRoute.run();
    });
    if (db::setting.multiNetVerbose >= +db::VerboseLevelT::MIDDLE) {
        printlog("postMT", postMT);
    }
    // final open fix
    if (db::setting.fixOpenBySST) {
        int count = 0;
        for (auto& net : database.nets) {
            if (net.defWireSegments.empty() && net.numOfPins() > 1) {
                // log() << " connect net " << net.idx << " by stt for " << allNetStatus[net.idx] << std::endl;
                connectBySTT(net);
                count++;
            }
        }
        if (count > 0) log() << "#nets connected by STT: " << count << std::endl;
    }
}

void Router::unfinish() {
    runJobsMT(database.nets.size(), [&](int netIdx) { database.nets[netIdx].clearPostRouteResult(); });
}

void Router::printStat(bool major) {
    log() << std::endl;
    log() << "----------------------------------------------------------------" << std::endl;
    db::routeStat.print();
    if (major) {
        database.printAllUsageAndVio();
    }
    log() << "----------------------------------------------------------------" << std::endl;
    log() << std::endl;
}
