#include "XYZParser.h"
#include "unity.h"

void setUp(void) {}

void tearDown(void) {}

void test_parse_comment(void)
{
    char comment_string[] = "Lattice=\"128.000000 -128.000000 -128.000000 128.000000 128.000000 "
                            "128.000000 -128.000000 128.000000 -128.000000\" "
                            "Origin=\"0.000000 0.000000 128.000000\" "
                            "pbc=\"T T T\" "
                            "Properties=id:I:1:species:S:1:pos:R:3 time=3.507347e-07 "
                            "temperature=293.000000 potential=0.900000 energy=140594.660000";

    // typedef struct KV TKV;
    struct KV expected_kvpairs[] = {
        {"lattice", "128.000000 -128.000000 -128.000000 128.000000 128.000000 128.000000 "
                    "-128.000000 128.000000 -128.000000"},
        {"origin", "0.000000 0.000000 128.000000"},
        {"pbc", "t t t"},
        {"properties", "id:i:1:species:s:1:pos:r:3"},
        {"time", "3.507347e-07"},
        {"temperature", "293.000000"},
        {"potential", "0.900000"},
        {"energy", "140594.660000"},
    };

    size_t kvpairs_cnt;
    struct KV *kvpairs; // pointer to array of KVs
    int pc = parse_comment(comment_string, &kvpairs, &kvpairs_cnt);
    TEST_ASSERT_EQUAL_INT(0, pc);
    TEST_ASSERT_EQUAL_size_t_MESSAGE(7, kvpairs_cnt, "Unexpected number of key-value pairs");

    for (int i = 0; i < (int)kvpairs_cnt; i++) {
        TEST_ASSERT_EQUAL_STRING(expected_kvpairs[i].key, kvpairs[i].key);
        TEST_ASSERT_EQUAL_STRING(expected_kvpairs[i].value, kvpairs[i].value);
    }
}

void test_parse_properties_value(void)
{
    char propval[] = "id:I:1:species:S:1:pos:R:3";

    PropertyDesc expected_properties[] = {
        {"id", 'I', 1},
        {"species", 'S', 1},
        {"pos", 'R', 3},
    };

    PropertyDesc *properties = NULL;
    int properties_cnt;
    int pp = parse_properties_value(propval, &properties, &properties_cnt);
    TEST_ASSERT_EQUAL_INT(0, pp);
    TEST_ASSERT_EQUAL_size_t_MESSAGE(3, properties_cnt, "Unexpected number of key-value pairs");

    for (int i = 0; i < properties_cnt; i++) {
        TEST_ASSERT_EQUAL_STRING(expected_properties[i].name, properties[i].name);
        TEST_ASSERT_EQUAL_CHAR(expected_properties[i].type, properties[i].type);
        TEST_ASSERT_EQUAL_INT(expected_properties[i].ncols, properties[i].ncols);
    }
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_parse_comment);
    RUN_TEST(test_parse_properties_value);
    UNITY_END();

    // return 0 else makefile throws error
    return 0;
}
