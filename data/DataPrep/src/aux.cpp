#include "aux.h"

string aux::read_samples(const char* samples_file) throw (runtime_error) {
   string line;
   stringstream buffer;
   ifstream if_samples;
   if_samples.exceptions(ifstream::failbit | ifstream::badbit);
   try {
      if_samples.open(samples_file);
      if (getline(if_samples, line)) {
         buffer << line;
      } 
      while (getline(if_samples, line)) {
         buffer << "," << line;
      }
      if_samples.close();
   } catch (exception& e) {
      if (!if_samples.eof()) {
         throw runtime_error("Error while reading samples file! Check if file exists and has read permissions.");
      }
      if_samples.close();
   }
   buffer.flush();
   return buffer.str();
}
