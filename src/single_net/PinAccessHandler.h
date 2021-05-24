#include "PostRoute.h"

class PinAccessHandler {
public:
    PinAccessHandler(db::Net &dbNet, db::GridSteiner *tapPtr, PostRoute *prData) :
                    net(dbNet), tap(tapPtr), pr(prData) {}

    void run();

    bool yes = false;

    db::Net &net;
    db::GridSteiner *tap;
    PostRoute *pr;
private:
    void editGridTopo(const utils::PointT<DBU> &vialoc,
                        const utils::PointT<DBU> &newloc,
                        const db::ViaType *type);

    void handleMacroPin();

    void handlePinWire();
    void offGridReroute();

    void handlePinVia(const db::ViaType *type);
    bool safeShift();
    // bool aggrShift();
    bool handleVioPinVia(const db::ViaType *type, const vector<utils::BoxT<DBU>> &neiMetals);
    bool isWider(const utils::BoxT<DBU> &viabot, DBU space, DBU wWidth);
    void patchPinVia(const db::ViaType *type);
};