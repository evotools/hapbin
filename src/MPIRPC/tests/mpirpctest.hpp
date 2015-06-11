#ifndef MPIRPCTEST_H
#define MPIRPCTEST_H

#include <QtCore/QObject>
#include <QtTest/QTest>
#include <string>

class MpirpcTest : public QObject
{
Q_OBJECT
private slots:
    void stream_uint8_t_test();
    void stream_uint8_t_test_data();
    
    void stream_uint16_t_test();
    void stream_uint16_t_test_data();
    
    void stream_uint32_t_test();
    void stream_uint32_t_test_data();
    
    void stream_uint64_t_test();
    void stream_uint64_t_test_data();
    
    void stream_int8_t_test();
    void stream_int8_t_test_data();
    
    void stream_int16_t_test();
    void stream_int16_t_test_data();
    
    void stream_int32_t_test();
    void stream_int32_t_test_data();
    
    void stream_int64_t_test();
    void stream_int64_t_test_data();
    
    void stream_double_test();
    void stream_double_test_data();
    
    void stream_string_test();
    void stream_string_test_data();

    void stream_charp_test();
    void stream_charp_test_data();

    void stream_combo();
};

Q_DECLARE_METATYPE(std::string)

#endif
