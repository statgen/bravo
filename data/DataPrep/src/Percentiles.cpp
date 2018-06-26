#include "Percentiles.h"

Percentiles::Percentiles(): raw_data_sorted(false) {
}

Percentiles::Percentiles(const string& name, const string& description): raw_data_sorted(false), name(name), description(description) {
}

Percentiles::Percentiles(const string& name, const string& description, vector<Percentiles>& many_percentiles): raw_data_sorted(false), name(name), description(description) {
   for (Percentiles& p : many_percentiles) {
      p.raw_data.shrink_to_fit();
   }
   for (Percentiles& p : many_percentiles) {
      raw_data.reserve(raw_data.size() + p.raw_data.size());
      raw_data.insert(raw_data.end(), p.raw_data.begin(), p.raw_data.end());
      p.raw_data.clear();
      p.raw_data.shrink_to_fit();
   }
}

Percentiles::~Percentiles() noexcept {
}

void Percentiles::sort_raw_data() {
   if (raw_data_sorted != true) {
      random_device rd;
      mt19937 g(rd());
      shuffle(raw_data.begin(), raw_data.end(), g); // we shuffle to have random order of equal elements after sort
      sort(raw_data.begin(), raw_data.end(), [](const raw_data_entry& e1, const raw_data_entry& e2) {
         return fcmp(e1.value, e2.value, 0.00000001) < 0;
      });
      raw_data_sorted = true;
   }
}

void Percentiles::add(double value, bool pass_filter) {
   raw_data.emplace_back(value, pass_filter ? 1u : 0u);
}

unsigned int Percentiles::n() noexcept {
   return raw_data.size();
}

unsigned int Percentiles::n_pass() noexcept {
   unsigned int n_pass = 0u;
   for (auto&& raw_data_entry: raw_data) {
      if (raw_data_entry.filter == 1u) {
         ++n_pass;
      }
   }
   return n_pass;
}

unsigned int Percentiles::n_nonpass() noexcept {
   return  raw_data.size() - n_pass();
}

void Percentiles::min_max(double& min, double& max) noexcept {
   if (raw_data.size() > 0) {
      min = max = raw_data[0].value;
      for (auto&& raw_data_entry: raw_data) {
         if (fcmp(raw_data_entry.value, max, 0.00000001) > 0) {
            max = raw_data_entry.value;
         }
         if (fcmp(min, raw_data_entry.value, 0.00000001) > 0) {
            min = raw_data_entry.value;
         }
      }
   }
}

void Percentiles::compute(vector<double> sorted_probabilities) noexcept {
   percentiles.clear();

   if ((raw_data.size() == 0) || (sorted_probabilities.size() == 0)) {
      return;
   }
   auto out_of_scope = find_if(sorted_probabilities.begin(), sorted_probabilities.end(), [](double p) {
      return ((fcmp(p, 0.0, 0.00000001) < 0) || (fcmp(p, 1.0, 0.00000001) > 0));
   });
   if (out_of_scope != sorted_probabilities.end()) {
      return;
   }
   if (raw_data.size() == 1) {
      for (double& p: sorted_probabilities) {
         percentiles.emplace_back(p, raw_data[0].value, 1u, raw_data[0].filter == 0u ? 0u : 1u );
      }
   } else {
      sort_raw_data();

      // use method #7 from R quantile() function
      double m = 0.0, j = 0.0, gamma = 0.0, percentile = 0.0;
      unsigned int i = 0u;
      unsigned int n = 0u, n_pass = 0u;
      for (double& p: sorted_probabilities) {
         m = 1.0 - p;
         j = floor(p * raw_data.size() + m);
         gamma = p * raw_data.size() + m - j;
         percentile = (1.0 - gamma) * raw_data[j - 1].value + gamma * raw_data[j].value;
         for (; (i < raw_data.size()) && (i < j) && (fcmp(raw_data[i].value, percentile, 0.00000001) <= 0); ++i, ++n) {
             if (raw_data[i].filter == 1u) {
                ++n_pass;
             }
         }
         percentiles.emplace_back(p, percentile, n, n_pass);
      } 
   }
}

void Percentiles::probability(double value, double& min, double& max) noexcept {
   if (raw_data.size() == 0) {
      min = max = 0.0;
      return;
   }
   sort_raw_data();
   raw_data_entry e(value, 0u);
   auto upper = upper_bound(raw_data.begin(), raw_data.end(), e, [](const raw_data_entry& e1, const raw_data_entry& e2) {
      return fcmp(e1.value, e2.value, 0.00000001) < 0;
   });
   auto lower = lower_bound(raw_data.begin(), raw_data.end(), e, [](const raw_data_entry& e1, const raw_data_entry& e2) {
      return fcmp(e1.value, e2.value, 0.00000001) < 0;
   });
   int i = upper - raw_data.begin();
   int j = lower - raw_data.begin();
   max = i / (double)raw_data.size();
   if (i - j == 1) {
      min = max;
   } else {
      if ((j > 0) && (fcmp(value, raw_data[j].value, 0.00000001) == 0)) {
         ++j;
      }
      min = j / (double)raw_data.size();
   }
}

void Percentiles::write(BGZF* f) throw (runtime_error) {
    double min = 0.0, max = 0.0;
    min_max(min, max);
    aux::write(f, "{");
    aux::write(f, "\"name\":\"%s\",", name.c_str());
    aux::write(f, "\"description\":\"%s\",", description.c_str());
    aux::write(f, "\"type\":\"percentiles\",");
    aux::write(f, "\"n\":%d,", n());
    aux::write(f, "\"n_pass\":%d,", n_pass());
    aux::write(f, "\"min\":%g,", min);
    aux::write(f, "\"max\":%g,", max);
    aux::write(f, "\"percentiles\":[");
    if (percentiles.size() > 0) {
        aux::write(f, "{\"probability\":%g,\"value\":%g,\"n\":%d,\"n_pass\":%d}", percentiles[0].probability, percentiles[0].value, percentiles[0].n, percentiles[0].n_pass);
        for (unsigned int i = 1u; i < percentiles.size(); ++i) {
            aux::write(f, ",{\"probability\":%g,\"value\":%g,\"n\":%d,\"n_pass\":%d}", percentiles[i].probability, percentiles[i].value, percentiles[i].n, percentiles[i].n_pass);
        }
    }
    aux::write(f, "]");
    aux::write(f, "}");
}
