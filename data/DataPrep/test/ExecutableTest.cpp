#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include "boost/process.hpp"

namespace pr = boost::process;
using namespace std;

class ExecutableTest : public::testing::Test {
protected:
    virtual ~ExecutableTest() {
    }

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }
};

TEST_F(ExecutableTest, CountAlleles) {
    string executable = "../bin/CountAlleles";
    vector<string> args = { "--in", "all_samples.test_input.vcf.gz", "--out", "all_samples.counts.test_output.vcf.gz" };
    int result = pr::system(executable, args);
    ASSERT_EQ(result, 0u);

    ifstream ok_file("all_samples.counts.test_ok.vcf.gz");
    ifstream test_output_file("all_samples.counts.test_output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));

};

TEST_F(ExecutableTest, Histograms) {
    string executable = "../bin/ComputeHistograms";
    vector<string> args = { "--in", "all_samples.test_input.vcf.gz", "--out", "all_samples.histograms.test_output.vcf.gz" };
    int result = pr::system(executable, args);
    ASSERT_EQ(result, 0u);

    ifstream ok_file("all_samples.histograms.test_ok.vcf.gz");
    ifstream test_output_file("all_samples.histograms.test_output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));
};
