/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef EHH_HPP
#define EHH_HPP

#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include "hapmap.hpp"

struct HapStats
{
    HapStats()
        : probs{}
        , probsNot{} {}
    double probs;
    double probsNot;
};

struct EHH
{
    EHH()
        : upstream(std::vector<HapStats>())
        , downstream(std::vector<HapStats>())
        , index{}
        , num{}
        , numNot{}
        , iHH_0{}
        , iHH_1{}
        {}
    std::vector<HapStats> upstream;
    std::vector<HapStats> downstream;

    unsigned long long index;

    int num;
    int numNot;

    double iHH_0;
    double iHH_1;

    ~EHH() {}

    void printEHH() const {
        for (auto it = upstream.crbegin(); it != upstream.crend(); ++it) {
            std::cout << it->probs << " " << it->probsNot << std::endl;
        }
        std::cout << "----" << std::endl;
        for (auto it = downstream.cbegin(); it != downstream.cend(); ++it) {
            std::cout << it->probs << " " << it->probsNot << std::endl;
        }
    }

    void printEHH(HapMap* hm)
    {
        if (upstream.size() == 0 && downstream.size() == 0)
        {
            std::cout << "EHH not calculated (<MAF)" << std::endl;
            return;
        }
        std::size_t i = upstream.size();
        while(i != 0)
        {
            --i;
            std::cout << hm->lineToId(index-i-1) << " "
                      << hm->geneticPosition(index-i-1) << " "
                      << hm->physicalPosition(index-i-1) << " "
                      << upstream.at(i).probs << " "
                      << upstream.at(i).probsNot << std::endl;
        }
        std::cout << hm->lineToId(index) << " "
                  << hm->geneticPosition(index) << " "
                  << hm->physicalPosition(index) << " "
                  << 1.0 << " " << 1.0 << std::endl;
        for (i = 0; i < downstream.size(); ++i)
        {
            std::cout << hm->lineToId(index+i+1) << " "
                      << hm->geneticPosition(index+i+1) << " "
                      << hm->physicalPosition(index+i+1) << " "
                      << downstream.at(i).probs << " "
                      << downstream.at(i).probsNot << std::endl;
        }
    }
};

struct XPEHH
{
    XPEHH()
        : index(0ULL)
        , xpehh(0.0)
        , numA(0)
        , numB(0)
        , numNotA(0)
        , numNotB(0)
        , iHH_A1(0.0)
        , iHH_B1(0.0)
        , iHH_P1(0.0)
    {}
    std::size_t index;

    double xpehh;

    int numA;
    int numB;
    int numNotA;
    int numNotB;

    double iHH_A1;
    double iHH_B1;
    double iHH_P1;
};

struct IhsScore
{
    IhsScore() : iHS(0.0), iHH_0(0.0), iHH_1(0.0), freq(0.0) {}
    IhsScore(double s, double a, double d, double f) : iHS(s), iHH_0(a), iHH_1(d), freq(f) {}
    double iHS;
    double iHH_0;
    double iHH_1;
    double freq;
};

#endif // EHH_HPP
