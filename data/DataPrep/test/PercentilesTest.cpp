#include <gtest/gtest.h>
#include <limits>
#include "../src/Percentiles.h"

using namespace std;

class PercentilesTest : public::testing::Test {
protected:
	virtual ~PercentilesTest() {
	}

	virtual void SetUp() {
	}

	virtual void TearDown() {
	}
};

TEST_F(PercentilesTest, CreatePercentiles) {
	Percentiles p;
	ASSERT_EQ(0u, p.n());
   ASSERT_EQ(0u, p.n_pass());
   ASSERT_EQ(0u, p.n_nonpass());
   double min = -1.0;
   double max = 1.0;
   p.min_max(min, max);
   ASSERT_EQ(-1.0, min);
   ASSERT_EQ(1.0, max);
   ASSERT_EQ(0, p.percentiles.size());
}

TEST_F(PercentilesTest, PopulatePercentiles1) {
   Percentiles p;
   double min = -1.0;
   double max = 1.0;
   p.add(0.9, false);
   p.add(0.1, false);
   p.add(0.4, false);
   p.add(1.2, true);
   p.add(0.5, true);
   ASSERT_EQ(5, p.n());
   ASSERT_EQ(2, p.n_pass());
   ASSERT_EQ(3, p.n_nonpass());
   p.min_max(min, max);
   ASSERT_EQ(0.1, min);
   ASSERT_EQ(1.2, max);
   ASSERT_EQ(0, p.percentiles.size());
}

TEST_F(PercentilesTest, ComputePercentiles1) {
   Percentiles p;
   vector<double> values = { 50 };
   vector<double> probabilities = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   vector<double> expected_percentiles = { 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0 };
   for (double v: values) {
      p.add(v, false);
   }
   p.compute(probabilities);
   ASSERT_EQ(probabilities.size(), p.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], p.percentiles[i].value, 0.0000001));
   }
}

TEST_F(PercentilesTest, ComputePercentiles2) {
   Percentiles p;
   vector<double> values = { 0.0, 1.0 };
   vector<double> probabilities = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   vector<double> expected_percentiles = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   for (double v: values) {
      p.add(v, false);
   }
   p.compute(probabilities);
   ASSERT_EQ(probabilities.size(), p.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], p.percentiles[i].value, 0.0000001));
   }
}

TEST_F(PercentilesTest, ComputePercentiles3) {
   Percentiles p;
   vector<double> values = { 25, 66, 58, -64, 16, -48, 40, 62, -8, 9, -12, 55, 22, -14, -57, 53, 82, -56, 63, 29 };
   vector<double> probabilities = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   vector<double> expected_percentiles = { -64.0, -56.1, -49.2, -20.8, -12.5, -9.2, 23.5,  58.8, 62.15, 63.3, 66.8, 82.0 };
   for (double v: values) {
      p.add(v, false);
   }
   p.compute(probabilities);
   ASSERT_EQ(probabilities.size(), p.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], p.percentiles[i].value, 0.0000001));
   }
}

TEST_F(PercentilesTest, ComputePassProportion) {
   Percentiles p;
   vector<double> values = { 25, 66, 58, -64, 16, -48, 40, 62, -8, 9, -12, 55, 22, -14, -57, 53, 82, -56, 63, 29 };
   vector<double> pass = { true, true, true, false, true, false, false, true, true, true, false, false, true, true, true, true, false, false, true, false };
   vector<double> probabilities = { 0.25, 0.5, 0.75, 1.0};
   vector<double> expected_percentiles = { -12.50, 23.50, 55.75, 82.0 };
   vector<unsigned int> expected_pass = { 2, 6, 8, 12 };
   vector<unsigned int> expected_total = { 5, 10, 15, 20 }; 
   ASSERT_EQ(values.size(), pass.size());
   for (unsigned int i = 0u; i < values.size(); ++i) {
      p.add(values[i], pass[i]);
   }
   p.compute(probabilities);
   ASSERT_EQ(probabilities.size(), p.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], p.percentiles[i].value, 0.0000001));
      ASSERT_EQ(expected_pass[i], p.percentiles[i].n_pass);
      ASSERT_EQ(expected_total[i], p.percentiles[i].n);
   }
}

TEST_F(PercentilesTest, ComputeProbability1) {
   Percentiles p;
   vector<double> values = { 25, 66, 58, -64, 16, -48, 40, 62, -8, 9, -12, 55, 22, -14, -57, 53, 82, -56, 63, 29 };
   double min = 1.0, max = 1.0;
   for (double v: values) {
      p.add(v, false);
   }
   
   p.probability(-65, min, max);
   ASSERT_EQ(0.0, min);
   ASSERT_EQ(0.0, max);

   p.probability(-64, min, max);
   ASSERT_EQ(0.05, min);
   ASSERT_EQ(0.05, max);

   p.probability(-58, min, max);
   ASSERT_EQ(0.05, min);
   ASSERT_EQ(0.05, max);

   p.probability(-57, min, max);
   ASSERT_EQ(0.10, min);
   ASSERT_EQ(0.10, max);

   p.probability(67, min, max);
   ASSERT_EQ(0.95, min);
   ASSERT_EQ(0.95, max);

   p.probability(82, min, max);
   ASSERT_EQ(1.0, min);
   ASSERT_EQ(1.0, max);

   p.probability(83, min, max);
   ASSERT_EQ(1.0, min);
   ASSERT_EQ(1.0, max);
}

TEST_F(PercentilesTest, ComputeProbability2) {
   Percentiles p;
   vector<double> values = { 25, 66, 58, -64, 16, -57, 40, 62, -8, 9, -12, 55, 22, -14, -57, 53, 82, -57, 63, 29 };
   vector<double> probabilities = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   vector<double> expected_percentiles = { -64.0, -57.0, -57.0, -22.6, -12.5, -9.2, 23.5, 58.8, 62.15, 63.3, 66.8, 82.0 };
   double min = 1.0, max = 1.0;
   for (double v: values) {
      p.add(v, false);
   }
   p.compute(probabilities);
   ASSERT_EQ(probabilities.size(), p.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], p.percentiles[i].value, 0.0000001));
   }
   p.probability(-57, min, max);
   ASSERT_EQ(0.1, min);
   ASSERT_EQ(0.2, max);
   p.probability(-14, min, max);
   ASSERT_EQ(0.25, min);
   ASSERT_EQ(0.25, max);
}

TEST_F(PercentilesTest, ComputeProbability3) {
   Percentiles p;
   vector<double> values = { 25, 66, 58, -57, 16, -57, 40, 62, -8, 9, -12, 55, 22, -14, -57, 53, 82, -57, 63, 29 };
   double min = 1.0, max = 1.0;
   for (double v: values) {
      p.add(v, false);
   }
   p.probability(-57, min, max);
   ASSERT_EQ(0.0, min);
   ASSERT_EQ(0.2, max);
}

TEST_F(PercentilesTest, ComputeProbability4) {
   Percentiles p;
   vector<double> values = { -57, -57, -57, -57 };
   double min = 1.0, max = 1.0;
   for (double v: values) {
      p.add(v, false);
   }
   p.probability(-57, min, max);
   ASSERT_EQ(0.0, min);
   ASSERT_EQ(1.0, max);
}

TEST_F(PercentilesTest, MergePercentiles) {
   vector<Percentiles> percentiles;
   percentiles.emplace_back();
   percentiles.emplace_back();
   percentiles.emplace_back();
   vector<double> values1 = { 25, 66, 58, -64 };
   vector<double> values2 = { 16, -57, 40, 62, -8, 9, -12 };
   vector<double> values3 = { 55, 22, -14, -57, 53, 82, -57, 63, 29 };
   vector<double> probabilities = { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0 };
   vector<double> expected_percentiles = { -64.0, -57.0, -57.0, -22.6, -12.5, -9.2, 23.5, 58.8, 62.15, 63.3, 66.8, 82.0 };
   for (double v: values1) {
      percentiles[0].add(v, false);
   }
   for (double v: values2) {
      percentiles[1].add(v, false);
   }
   for (double v: values3) {
      percentiles[2].add(v, false);
   }
   Percentiles percentile("", "", percentiles);
   ASSERT_EQ(20, percentile.n());
   percentile.compute(probabilities);
   ASSERT_EQ(probabilities.size(), percentile.percentiles.size());
   for (unsigned int i = 0u; i < expected_percentiles.size(); ++i) {
      ASSERT_EQ(0, fcmp(expected_percentiles[i], percentile.percentiles[i].value, 0.0000001));
   }
}
