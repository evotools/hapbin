#include "mpirpctest.hpp"

#include "../parameterstream.hpp"
#include <QDebug>
#include <type_traits>

template<typename T>
T testParamStream(T t)
{
    std::vector<char> buffer;
    ParameterStream s(&buffer);
    s << t;
    s.seek(0);
    T tmp;
    s >> tmp;
    return tmp;
}

void MpirpcTest::stream_double_test_data()
{
    QTest::addColumn<double>("val");
    
    QTest::newRow("2.0") << 2.0;
    QTest::newRow("3.14") << 3.14;
}

void MpirpcTest::stream_double_test()
{
    QFETCH(double, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_string_test_data()
{
    QTest::addColumn<std::string>("val");
    
    QTest::newRow("test") << std::string("test");
}

void MpirpcTest::stream_string_test()
{
    QFETCH(std::string, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_charp_test_data()
{
    QTest::addColumn<std::string>("val");
    QTest::newRow("test") << std::string("test");
}

void MpirpcTest::stream_charp_test()
{
    QFETCH(std::string, val);
    std::vector<char> buffer;
    ParameterStream s(&buffer);;
    s << val.c_str();
    s.seek(0);
    char* tmp;
    s >> tmp;
    QCOMPARE(strcmp(tmp, val.c_str()), 0);
}

void MpirpcTest::stream_uint64_t_test_data() {
    QTest::addColumn<uint64_t>("val");
    
    QTest::newRow("1234567890") << (uint64_t) 1234567890UL;
}
void MpirpcTest::stream_uint64_t_test() {
    QFETCH(uint64_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_uint32_t_test_data() {
    QTest::addColumn<uint32_t>("val");
    
    QTest::newRow("1234567890") << (uint32_t) 1234567890U;
}
void MpirpcTest::stream_uint32_t_test() {
    QFETCH(uint32_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_uint16_t_test_data() {
    QTest::addColumn<uint16_t>("val");
    
    QTest::newRow("12345") << (uint16_t) 12345;
}
void MpirpcTest::stream_uint16_t_test() {
    QFETCH(uint16_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_uint8_t_test_data() {
    QTest::addColumn<uint8_t>("val");
    
    QTest::newRow("200") << (uint8_t) 200;
}
void MpirpcTest::stream_uint8_t_test() {
    QFETCH(uint8_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_int64_t_test_data() {
    QTest::addColumn<int64_t>("val");
    
    QTest::newRow("1234567890") << (int64_t) 123456789L;
}
void MpirpcTest::stream_int64_t_test() {
    QFETCH(int64_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_int32_t_test_data() {
    QTest::addColumn<int32_t>("val");
    
    QTest::newRow("1234567890") << (int32_t) 1234567890;
}
void MpirpcTest::stream_int32_t_test() {
    QFETCH(int32_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_int16_t_test_data() {
    QTest::addColumn<int16_t>("val");
    
    QTest::newRow("12345") << (int16_t) 12345;
}
void MpirpcTest::stream_int16_t_test() {
    QFETCH(int16_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_int8_t_test_data() {
    QTest::addColumn<int8_t>("val");
    
    QTest::newRow("100") << (int8_t) 100;
}
void MpirpcTest::stream_int8_t_test() {
    QFETCH(int8_t, val);
    QCOMPARE(testParamStream(val), val);
}

void MpirpcTest::stream_combo() {
    std::vector<char> buffer;
    ParameterStream s(&buffer);
    float a = 1.1f, a2;
    double b = 2.2, b2;
    uint32_t c = 1234, c2;
    int64_t d = 2345542341L, d2;
    const char * charp = "blah blah";
    char* charp2;
    std::string str("testtest!");
    std::string str2;
    s << a << str << b << c << charp << d;
    s >> a2 >> str2 >> b2 >> c2 >> charp2 >> d2;
    bool ok = true;
    if (a != a2)
        ok = false;
    if (str != str2)
        ok = false;
    if (b != b2)
        ok = false;
    if (c != c2)
        ok = false;
    if (d != d2)
        ok = false;
    if (strcmp(charp, charp2) != 0)
        ok = false;
    QVERIFY(ok == true);
}

QTEST_APPLESS_MAIN(MpirpcTest)
