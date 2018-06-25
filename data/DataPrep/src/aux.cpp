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

void aux::write(BGZF* f, const char* format, ...) throw (runtime_error) {
   va_list arguments;
   long int n = 0;
   const unsigned int BUFFER_SIZE = 32768u;
   unsigned int max_string_length = BUFFER_SIZE - 1u;
   char buffer[BUFFER_SIZE];

   va_start(arguments, format);
   if ((n = vsnprintf(buffer, max_string_length, format, arguments)) < 0) {
      throw runtime_error("Error while writing to memory buffer!");
   } else if (n > max_string_length) {
      throw runtime_error("Too small memory buffer size for writing!");
   }
   va_end(arguments);

   if (bgzf_write(f, buffer, n) < n) {
      throw runtime_error("Error while writing to output file!");
   }
}

