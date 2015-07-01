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
        , num{}
        , numNot{}
        , iHH_a{}
        , iHH_d{}
        , index{}
        {}
    std::vector<HapStats> upstream;
    std::vector<HapStats> downstream;
    
    unsigned long long index;
    
    int num;
    int numNot;
    
    double iHH_a;
    double iHH_d;
    
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

#endif // EHH_HPP
