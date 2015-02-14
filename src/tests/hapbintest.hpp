#ifndef HAPBINTEST_HPP
#define HAPBINTEST_HPP

#include <QtCore/QObject>
#include <deque>

class HapbinTest : public QObject
{
Q_OBJECT
private slots:
    void save();
    void testIHS();
    
    void initTestCase();
    void cleanupTestCase();
    
};

class HapMapNode;
void printHapMapTree(std::list<HapMapNode*> tiers);

QDebug& operator<<(QDebug d, const std::string &s);

Q_DECLARE_METATYPE(std::string)

#endif // HAPBINTEST_HPP
