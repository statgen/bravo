#ifndef SRC_GZIPWRITER_H_
#define SRC_GZIPWRITER_H_

#include <iostream>
#include <cstdarg>
#include <cstring>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <new>
#include <zlib.h>

using namespace std;

class GzipWriter {
private:
	char* buffer;
	unsigned int max_string_length;

	gzFile gzfile;

public:
	static const unsigned int DEFAULT_BUFFER_SIZE;

	GzipWriter(unsigned int buffer_size = DEFAULT_BUFFER_SIZE) throw (runtime_error);
	virtual ~GzipWriter();

	void open(const char* file_name) throw (runtime_error);
	void close() throw (runtime_error);

	void write(const char* format, ...) throw (runtime_error);

};

#endif
