#include "db/Database.h"

class PostProcess {
public:
    static void removeSameNetVioVias(db::Net &net);
    static void removeSameTreeVias(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx);
    static void removeDiffTreeVias(std::shared_ptr<db::GridSteiner> u, std::shared_ptr<db::GridSteiner> v, int netIdx);
    static void removeCorners(std::shared_ptr<db::GridSteiner> node, int idx);
    static void fixMAR(db::Net &net);
};