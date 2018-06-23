#ifndef PERCENTILES_H_
#define PERCENTILES_H_

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <iterator>
#include "aux.h"

using namespace std;
using namespace aux;

class Percentiles {
private:
   struct raw_data_entry {
      double value;
      unsigned char filter;
      raw_data_entry(double value, unsigned char filter): value(value), filter(filter) {}
   };

   vector<raw_data_entry> raw_data;

   void sort_raw_data();
   bool raw_data_sorted;

   string name;
   string description;

public:
   struct percentile {
      double probability;
      double value;
      unsigned int n;
      unsigned int n_pass;
      percentile(double probability, double value, unsigned int n, unsigned int n_pass):
         probability(probability), value(value), n(n), n_pass(n_pass) {}
   };

   vector<percentile> percentiles;

   Percentiles();
   Percentiles(const string& name, const string& description);
   Percentiles(const string& name, const string& description, vector<Percentiles>& many_percentiles);
   virtual ~Percentiles() noexcept;

   void add(double value, bool pass_filter);
   unsigned int n() noexcept;
   unsigned int n_pass() noexcept;
   unsigned int n_nonpass() noexcept;
   void min_max(double& min, double& max) noexcept;
   void probability(double value, double& min, double& max) noexcept;
   void compute(vector<double> sorted_probabilities) noexcept;
   void write(BGZF* f) throw (runtime_error);
};

#endif
