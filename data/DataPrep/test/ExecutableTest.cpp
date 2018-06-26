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

TEST_F(ExecutableTest, CountAllelesAll) {
    string executable = "../bin/ComputeAlleleCounts";
    vector<string> args = { "--in", "input.vcf.gz", "--out", "allelecounts.output.vcf.gz" };
    ASSERT_EQ(pr::system(executable, args), 0u);
    ifstream ok_file("allelecounts.output_ok.vcf.gz");
    ifstream test_output_file("allelecounts.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));

}

TEST_F(ExecutableTest, CountAllelesSubset) {
    string executable = "../bin/ComputeAlleleCounts";
    vector<string> args1 = { "--in", "input.vcf.gz", "--out", "allelecounts.samples_option.output.vcf.gz", "--samples", "samples.txt" };
    vector<string> args2 = { "--in", "input.subset.vcf.gz", "--out", "allelecounts.subset.output.vcf.gz" };
    ASSERT_EQ(pr::system(executable, args1), 0u);
    ASSERT_EQ(pr::system(executable, args2), 0u);
    ifstream file1("allelecounts.samples_option.output.vcf.gz");
    ifstream file2("allelecounts.subset.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(file1.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(file2.rdbuf())));

}

TEST_F(ExecutableTest, CountAllelesRegion) {
    string executable = "../bin/ComputeAlleleCounts";
    vector<string> args = { "--in", "input.vcf.gz", "--out", "allelecounts.region.output.vcf.gz", "--region", "22:16050739-16051721" };
    ASSERT_EQ(pr::system(executable, args), 0u);
    ifstream ok_file("allelecounts.region.output_ok.vcf.gz");
    ifstream test_output_file("allelecounts.region.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));

}

TEST_F(ExecutableTest, HistogramsAll) {
    string executable = "../bin/ComputeHistograms";
    vector<string> args = { "--in", "input.vcf.gz", "--out", "histograms.output.vcf.gz" };
    ASSERT_EQ(pr::system(executable, args), 0u);
    ifstream ok_file("histograms.output_ok.vcf.gz");
    ifstream test_output_file("histograms.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));
}

TEST_F(ExecutableTest, HistogramsSubset) {
    string executable = "../bin/ComputeHistograms";
    vector<string> args1 = { "--in", "input.vcf.gz", "--out", "histograms.samples_option.output.vcf.gz", "--samples", "samples.txt" };
    vector<string> args2 = { "--in", "input.subset.vcf.gz", "--out", "histograms.subset.output.vcf.gz" };
    ASSERT_EQ(pr::system(executable, args1), 0u);
    ASSERT_EQ(pr::system(executable, args2), 0u);
    ifstream file1("histograms.samples_option.output.vcf.gz");
    ifstream file2("histograms.subset.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(file1.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(file2.rdbuf())));
}

TEST_F(ExecutableTest, HistogramsRegion) {
    string executable = "../bin/ComputeHistograms";
    vector<string> args = { "--in", "input.vcf.gz", "--out", "histograms.region.output.vcf.gz", "--region", "22:16050739-16051721" };
    ASSERT_EQ(pr::system(executable, args), 0u);
    ifstream ok_file("histograms.region.output_ok.vcf.gz");
    ifstream test_output_file("histograms.region.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));
}

TEST_F(ExecutableTest, AlleleCountsAndHistograms) {
    string executable = "../bin/ComputeAlleleCountsAndHistograms";
    vector<string> args = { "--in", "input.vcf.gz", "--out", "allelecounts_and_histograms.output.vcf.gz" };
    ASSERT_EQ(pr::system(executable, args), 0u);
    ifstream ok_file("allelecounts_and_histograms.output_ok.vcf.gz");
    ifstream test_output_file("allelecounts_and_histograms.output.vcf.gz");
    ASSERT_TRUE(equal(istreambuf_iterator<char>(ok_file.rdbuf()), istreambuf_iterator<char>(), istreambuf_iterator<char>(test_output_file.rdbuf())));
}

