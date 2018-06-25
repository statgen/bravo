#include <gtest/gtest.h>
#include <limits>
#include "../src/Histogram.h"

using namespace std;

class HistogramTest : public::testing::Test {
protected:
	virtual ~HistogramTest() {
	}

	virtual void SetUp() {
	}

	virtual void TearDown() {
	}
};

TEST_F(HistogramTest, CreateHistogram) {
	vector<double> expected_borders = { 0, 11, 21, 31, 41, 51, 61, 71, 81, 91,  numeric_limits<double>::max() };
	vector<double> borders = { 61, 0, 11, 31, 21, 41, numeric_limits<double>::max(), 51, 81, 71, 91 };

	Histogram histogram(borders);

	vector<pair<double, unsigned int>> observed_bins = histogram.get_bins();

	ASSERT_EQ(expected_borders.size(), observed_bins.size());
	for (unsigned int i = 0u; i < expected_borders.size(); ++i) {
		ASSERT_EQ(expected_borders[i], observed_bins[i].first);
		ASSERT_EQ(observed_bins[i].second, 0u);
	}
}

TEST_F(HistogramTest, PopulateHistogram1) {
	vector<pair<double, unsigned int>> expected_bins = {
			pair<double, unsigned int>(0, 2u),
			pair<double, unsigned int>(11, 1u),
			pair<double, unsigned int>(21, 4u),
			pair<double, unsigned int>(31, 0u),
			pair<double, unsigned int>(41, 2u),
			pair<double, unsigned int>(51, 0u),
			pair<double, unsigned int>(61, 0u),
			pair<double, unsigned int>(71, 0u),
			pair<double, unsigned int>(81, 0u),
			pair<double, unsigned int>(91, 3u),
//			pair<double, unsigned int>(101, 0u)
			pair<double, unsigned int>(100, 0u)
	};

//	vector<double> bins = { 61, 0, 11, 31, 21, 41, 101, 51, 81, 71, 91 };
	vector<double> bins = { 61, 0, 11, 31, 21, 41, 100, 51, 81, 71, 91 };
	Histogram histogram(bins);

	// check extreme values
	histogram.add(-1);
	histogram.add(101);

	// check values that fall into boundaries
	histogram.add(0);
	histogram.add(10);
	histogram.add(11);
	histogram.add(21);
	histogram.add(22);
	histogram.add(100);

	// check values tha fall in the middle
	histogram.add(45);
	histogram.add(46);
	histogram.add(25);
	histogram.add(26);
	histogram.add(92);
	histogram.add(93);

	vector<pair<double, unsigned int>> observed_bins = histogram.get_bins();
    ASSERT_EQ(expected_bins.size(), observed_bins.size());
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(expected_bins[i].second, observed_bins[i].second);
	}

	ASSERT_EQ(histogram.get_n(), 14u);
	ASSERT_EQ(histogram.get_total(), 591.0);
	ASSERT_DOUBLE_EQ(histogram.get_average(), 42.21428571428571);

	histogram.clear();
	observed_bins = histogram.get_bins();
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(0u, observed_bins[i].second);
	}

    ASSERT_EQ(histogram.get_n(), 0u);
    ASSERT_EQ(histogram.get_total(), 0.0);
}

TEST_F(HistogramTest, PopulateHistogram2) {
	vector<pair<double, unsigned int>> expected_bins = {
			pair<double, unsigned int>(-numeric_limits<double>::max(), 3u),
			pair<double, unsigned int>(11, 1u),
			pair<double, unsigned int>(21, 4u),
			pair<double, unsigned int>(31, 0u),
			pair<double, unsigned int>(41, 2u),
			pair<double, unsigned int>(51, 0u),
			pair<double, unsigned int>(61, 0u),
			pair<double, unsigned int>(71, 0u),
			pair<double, unsigned int>(81, 0u),
			pair<double, unsigned int>(91, 3u),
			pair<double, unsigned int>(numeric_limits<double>::max(), 0u)
	};

	vector<double> bins = { 61, -numeric_limits<double>::max(), 11, 31, 21, 41, 91, 51, 81, 71, numeric_limits<double>::max() };
	Histogram histogram(bins);

	// check extreme values
	histogram.add(-1);
	histogram.add(101);

	// check values that fall into boundaries
	histogram.add(0);
	histogram.add(10);
	histogram.add(11);
	histogram.add(21);
	histogram.add(22);

	// check values that fall in the middle
	histogram.add(45);
	histogram.add(46);
	histogram.add(25);
	histogram.add(26);
	histogram.add(92);
	histogram.add(93);

    ASSERT_EQ(histogram.get_n(), 13u);
    ASSERT_EQ(histogram.get_total(), 491.0);
    ASSERT_DOUBLE_EQ(histogram.get_average(), 37.76923076923077);

	vector<pair<double, unsigned int>> observed_bins = histogram.get_bins();
    ASSERT_EQ(expected_bins.size(), observed_bins.size());
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(expected_bins[i].second, observed_bins[i].second);
	}

	histogram.clear();
	observed_bins = histogram.get_bins();
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(0u, observed_bins[i].second);
	}

    ASSERT_EQ(histogram.get_n(), 0u);
    ASSERT_EQ(histogram.get_total(), 0.0);
}

TEST_F(HistogramTest, PopulateHistogram3) {
	vector<pair<double, unsigned int>> expected_bins = {
			pair<double, unsigned int>(10, 3u),
	};

	vector<double> bins = { 10, 10 };
	Histogram histogram(bins);

	// check extreme values
	histogram.add(-1);
	histogram.add(101);

	// check values that fall into boundaries
	histogram.add(10);
	histogram.add(10);
	histogram.add(10);

    ASSERT_EQ(histogram.get_n(), 5u);
    ASSERT_EQ(histogram.get_total(), 130.0);
    ASSERT_DOUBLE_EQ(histogram.get_average(), 26.0);

	vector<pair<double, unsigned int>> observed_bins = histogram.get_bins();
    ASSERT_EQ(expected_bins.size(), observed_bins.size());
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(expected_bins[i].second, observed_bins[i].second);
	}

	histogram.clear();
	observed_bins = histogram.get_bins();
	for (unsigned int i = 0u; i < expected_bins.size(); ++i) {
		ASSERT_EQ(expected_bins[i].first, observed_bins[i].first);
		ASSERT_EQ(0u, observed_bins[i].second);
	}

    ASSERT_EQ(histogram.get_n(), 0u);
    ASSERT_EQ(histogram.get_total(), 0.0);
}


