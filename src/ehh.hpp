/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
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
#include <tuple>
#include "hapmap.hpp"

struct HapStats
{
    HapStats()
        : probs{}
        , probsNot{} {}
    double probs;
    double probsNot;
};

struct HapPairStats
{
    HapPairStats()
        : probs{}
        , probsNot{}
        {}
    HapPairStats(const double (&p)[4], const double (&pn)[4])
        : probs{p[0],p[1],p[2],p[3]}
        , probsNot{pn[0],pn[1],pn[2],pn[3]}
        {}
    double probs[4];
    double probsNot[4];
};

struct EHHPair
{
    std::vector<HapPairStats> upstream;
    std::vector<HapPairStats> downstream;
    std::size_t focus[2];
    int num[4];
    int numNot[4];

    int count[4];
    double af[4];
    double ihh[4];
    double ihhNot[4];
    double sl[4];
    double slNot[4];

    void printEHH() const {
        for (auto it = upstream.crbegin(); it != upstream.crend(); ++it) {
            std::cout << it->probs[0b00] << " " << it->probs[0b01] << " " << it->probs[0b10] << " " << it->probs[0b11]
                << it->probsNot[0b00] << " " << it->probsNot[0b01] << " " << it->probsNot[0b10] << " " << it->probsNot[0b11] << std::endl;
        }
        for (auto it = downstream.cbegin(); it != downstream.cend(); ++it) {
            std::cout << it->probs[0b00] << " " << it->probs[0b01] << " " << it->probs[0b10] << " " << it->probs[0b11]
                << it->probsNot[0b00] << " " << it->probsNot[0b01] << " " << it->probsNot[0b10] << " " << it->probsNot[0b11] << std::endl;
        }
    }

    void printEHH(const HapMap *hm) const {
        if (upstream.size() == 0 && downstream.size() == 0)
        {
            std::cout << "EHH not calculated (<MAF)" << std::endl;
            return;
        }
        std::size_t i = upstream.size();
        while(i != 0)
        {
            --i;
            std::cout << hm->indexToId(focus[0]-i) << " "
                      << hm->geneticPosition(focus[0]-i) << " "
                      << hm->physicalPosition(focus[0]-i) << " "
                      << upstream.at(i).probs[0b00] << " "
                      << upstream.at(i).probs[0b01] << " "
                      << upstream.at(i).probs[0b10] << " "
                      << upstream.at(i).probs[0b11] << " "
                      << upstream.at(i).probsNot[0b00] << " "
                      << upstream.at(i).probsNot[0b01] << " "
                      << upstream.at(i).probsNot[0b10] << " "
                      << upstream.at(i).probsNot[0b11] << std::endl;
        }
        /*std::cout << "------" << std::endl;
        std::cout << hm->indexToId(focus[0]) << " "
                  << hm->geneticPosition(focus[0]) << " "
                  << hm->physicalPosition(focus[0]) << " "
                  << 1.0 << " " << 1.0 << " " << 1.0 << " " << 1.0 << " "
                  << 1.0 << " " << 1.0 << " " << 1.0 << " " << 1.0 << std::endl;
        std::cout << hm->indexToId(focus[1]) << " "
                  << hm->geneticPosition(focus[1]) << " "
                  << hm->physicalPosition(focus[1]) << " "
                  << 1.0 << " " << 1.0 << " " << 1.0 << " " << 1.0 << " "
                  << 1.0 << " " << 1.0 << " " << 1.0 << " " << 1.0 << std::endl;
        std::cout << "------" << std::endl;*/
        for (i = 0; i < downstream.size(); ++i)
        {
            std::cout << hm->indexToId(focus[1]+i) << " "
                      << hm->geneticPosition(focus[1]+i) << " "
                      << hm->physicalPosition(focus[1]+i) << " "
                      << downstream.at(i).probs[0b00] << " "
                      << downstream.at(i).probs[0b01] << " "
                      << downstream.at(i).probs[0b10] << " "
                      << downstream.at(i).probs[0b11] << " "
                      << downstream.at(i).probsNot[0b00] << " "
                      << downstream.at(i).probsNot[0b01] << " "
                      << downstream.at(i).probsNot[0b10] << " "
                      << downstream.at(i).probsNot[0b11] << std::endl;
        }
    }
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
        , sl_0{}
        , sl_1{}
        {}
    std::vector<HapStats> upstream;
    std::vector<HapStats> downstream;

    std::size_t index;

    int num;
    int numNot;

    double iHH_0;
    double iHH_1;

    double sl_0;
    double sl_1;

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

    void printEHH(const HapMap* hm)
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
            std::cout << hm->indexToId(index-i-1) << " "
                      << hm->geneticPosition(index-i-1) << " "
                      << hm->physicalPosition(index-i-1) << " "
                      << upstream.at(i).probs << " "
                      << upstream.at(i).probsNot << std::endl;
        }
        std::cout << hm->indexToId(index) << " "
                  << hm->geneticPosition(index) << " "
                  << hm->physicalPosition(index) << " "
                  << 1.0 << " " << 1.0 << std::endl;
        for (i = 0; i < downstream.size(); ++i)
        {
            std::cout << hm->indexToId(index+i+1) << " "
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
        , sl_A1(0.0)
        , sl_B1(0.0)
        , sl_P1(0.0)
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

    double sl_A1;
    double sl_B1;
    double sl_P1;
};

struct IhsScore
{
    IhsScore()
        : iHS{}
        , iHH_0{}
        , iHH_1{}
        , nSL{}
        , sl_0{}
        , sl_1{}
        {}
    IhsScore(double s, double a, double d, double n, double na, double nd)
        : iHS{s}
        , iHH_0{a}
        , iHH_1{d}
        , nSL{n}
        , sl_0{na}
        , sl_1{nd}
        {}
    double iHS;
    double iHH_0;
    double iHH_1;
    double nSL;
    double sl_0;
    double sl_1;
};

struct StdScorePair
{
    StdScorePair()
        : freq{}
        , score_0{}
        , score_1{}
        , score{}
        {}
    StdScorePair(const double (&r)[4],
                 const double (&s)[4],
                 const double (&s0)[4],
                 const double (&s1)[4]
                )
        : freq{r[0],r[1],r[2],r[3]}
        , score_0{s0[0],s0[1],s0[2],s0[3]}
        , score_1{s1[0],s1[1],s1[2],s0[3]}
        , score{s[0],s[1],s[2],s[3]}
        {}
    //std::size_t focus[2];
    double freq[4];
    double score_0[4];
    double score_1[4];
    double score[4];
};

#endif // EHH_HPP
