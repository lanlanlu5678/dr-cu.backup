#pragma once

#include "db/Net.h"

class PseudoNet : public db::NetBase {
public:
    PseudoNet(const db::Net& _dbNet) : db::Net(_dbNet), db::NetBase(_dbNet) {}

    const db::Net& dbNet;
}