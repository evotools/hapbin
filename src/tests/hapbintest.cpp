#include "hapbintest.hpp"
#include "../hapbin.hpp"
#include "../hapmap.hpp"
#include "../snpset.hpp"
#include <QtTest/QTest>
#include <list>
#include <QtCore/QDebug>
#include <QTextStream>
#include <QFile>
#include <QCryptographicHash>
#include "../calcselect.hpp"

QDebug& operator<<(QDebug d, const std::string& s)
{
    d.nospace() << QString(s.c_str());
    return d.space();
}

void HapbinTest::initTestCase()
{
    QFile f("test1.hap");
    if (f.open(QFile::ReadWrite | QFile::Truncate))
    {
        QTextStream out(&f);
        out << "0 0 0 0 0 0 0 0\n"
            << "0 0 0 0 1 1 1 1\n"
            << "0 0 1 1 0 0 0 0\n"
            << "0 0 0 0 0 0 1 1\n"
            << "1 0 0 0 0 0 0 0\n"
            << "0 0 1 0 0 0 0 0\n"
            << "0 0 0 0 1 0 0 0\n"
            << "0 0 0 0 0 0 1 0\n";
    }
}

void HapbinTest::cleanupTestCase()
{
    QFile f1("test1.hap"), f2("test1.bin");
    f1.remove();
    f2.remove();
}

void HapbinTest::save()
{
    HapMap map;
    if (!map.loadHapAscii("test1.hap"))
    {
        QFAIL("Test file test1.hap not found");
    }
    
    map.save("test1.bin");
}

void HapbinTest::testIHS()
{
    QString hap = QFINDTESTDATA("test.500.bin");
    QString map = QFINDTESTDATA("test.map");
    QString corr = QFINDTESTDATA("ihs.txt");
    calcIhsNoMpi(hap.toLatin1().constData(), map.toLatin1().constData(), "out.txt", false, 0.05, 0.05);
    QCryptographicHash hash1( QCryptographicHash::Sha1 );
    QCryptographicHash hash2( QCryptographicHash::Sha1 );
    QFile file1("out.txt");
    QFile file2(corr.toLatin1().constData());
    if (file1.open(QIODevice::ReadOnly))
    {
        hash1.addData(file1.readAll());
        file1.close();
        if (file2.open(QIODevice::ReadOnly))
        {
            hash2.addData(file2.readAll());
            file2.close();
            QCOMPARE(hash1.result(), hash2.result());
            return;
        }
        QFAIL("Could not open ihs.txt");
        return;
    }
    QFAIL("Could not open out.txt");
}


QTEST_APPLESS_MAIN(HapbinTest)
#include "hapbintest.moc"
