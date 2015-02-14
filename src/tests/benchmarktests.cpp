/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2014  Colin MacLean <colin@colin-maclean.com>
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
#include "benchmarktests.hpp"
#include "../hapmap.hpp"
#include "../snpset.hpp"
#include "../calcselect.hpp"

#include <QtTest/QtTest>

QTEST_MAIN(BenchmarkTests);

void BenchmarkTests::ihs()
{
    QString hap = QFINDTESTDATA("../../data/chr22.500.hap");
    QString map = QFINDTESTDATA("../../data/test_chr22.map");
    QBENCHMARK {
        calcIhsNoMpi(hap.toLatin1().constData(), map.toLatin1().constData(), "out.txt", true, 0.05, 0.05);
    }
}

#include "benchmarktests.moc"
