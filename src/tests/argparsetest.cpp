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
#include "argparsetest.hpp"

#include <QtTest>
#include <QString>
#include <QCoreApplication>
#include <QDebug>

#include "../argparse.hpp"

QTEST_MAIN(ArgParseTest);

void ArgParseTest::initTestCase()
{
    // Called before the first testfunction is executed
}

void ArgParseTest::cleanupTestCase()
{
    // Called after the last testfunction was executed
}

void ArgParseTest::init()
{
    // Called before each testfunction is executed
}

void ArgParseTest::cleanup()
{
    // Called after every testfunction
}

void ArgParseTest::testFlag()
{
    int argc = 2;
    const char* argv[2] = { QCoreApplication::arguments().at(0).toLatin1(), "--test"};
    Argument<bool> flagarg('t', "test", "Test flag", true, false);
    ArgParse parse({&flagarg});
    parse.parseArguments(argc, (char**) argv);
    if (!flagarg.wasFound())
        QFAIL("Flag was not found");
    if (!flagarg.value())
        QFAIL("Flag not set");
}

void ArgParseTest::testFlagShort()
{
    int argc = 2;
    const char* argv[2] = { QCoreApplication::arguments().at(0).toLatin1(), "-t"};
    Argument<bool> flagarg('t', "test", "Test flag", true, false);
    ArgParse parse({&flagarg});
    parse.parseArguments(argc, (char**) argv);
    if (!flagarg.wasFound())
        QFAIL("Flag was not found");
    if (!flagarg.value())
        QFAIL("Flag not set");
}

void ArgParseTest::testFlagNotFound()
{
    int argc = 2;
    const char* argv[2] = { QCoreApplication::arguments().at(0).toLatin1(), "-h"};
    Argument<bool> flagarg('t', "test", "Test flag", true, false);
    ArgParse parse({&flagarg});
    parse.parseArguments(argc, (char**) argv);
    if (flagarg.wasFound())
        QFAIL("Flag was found when it shouldn't be!");
}

void ArgParseTest::testInt()
{
    int argc = 3;
    const char* argv[3] = { QCoreApplication::arguments().at(0).toLatin1(), "--test", "123"};
    Argument<int> arg1('t', "test", "Test flag", false, false, 500);
    ArgParse parse({&arg1});
    parse.parseArguments(argc, (char**) argv);
    if (!arg1.wasFound())
        QFAIL("Argument was not found");
    if (arg1.value() != 123)
        QFAIL("Incorrect value");
}

void ArgParseTest::testIntShort()
{
    int argc = 3;
    const char* argv[3] = { QCoreApplication::arguments().at(0).toLatin1(), "-t", "123"};
    Argument<int> arg1('t', "test", "Test flag", false, false, 500);
    ArgParse parse({&arg1});
    parse.parseArguments(argc, (char**) argv);
    if (!arg1.wasFound())
        QFAIL("Argument was not found");
    if (arg1.value() != 123)
    {
        std::cout << arg1.value() << std::endl;
        QFAIL("Incorrect value");
    }
}

void ArgParseTest::testIntParseError()
{
    int argc = 3;
    const char* argv[3] = { QCoreApplication::arguments().at(0).toLatin1(), "-t", "4564561545456545465"};
    Argument<int> arg1('t', "test", "Test flag", false, false, 500);
    ArgParse parse({&arg1});
    bool fail = parse.parseArguments(argc, (char**) argv);
    if (!arg1.wasFound())
        QFAIL("Argument was not found");
    if (!fail)
    {
        QFAIL("This should fail");
    }
}

void ArgParseTest::testStdString()
{
    int argc = 3;
    char* test1 = new char[7];
    char* test2 = new char[6];
    strcpy(test1, "--test");
    strcpy(test2, "foo");
    char* argv[3];
    argv[0] = QCoreApplication::arguments().at(0).toLatin1().data();
    argv[1] = test1;
    argv[2] = test2;
    Argument<std::string> arg1('t', "test", "Test flag", false, false, "bar");
    ArgParse parse({&arg1});
    parse.parseArguments(argc, (char**) argv);
    if (!arg1.wasFound())
        QFAIL("Argument was not found");
    if (arg1.value() != std::string("foo"))
        QFAIL("Incorrect value");
}

void ArgParseTest::testRemoveBackslashes()
{
    char test[] = "\\-\\-test";
    removeBackslashes(test);
    QVERIFY(strcmp(test, "--test") == 0);
}

void ArgParseTest::testRemoveBackslashes2()
{
    char test[] = "test";
    removeBackslashes(test);
    QVERIFY(strcmp(test, "test") == 0);
}

#include "argparsetest.moc"
