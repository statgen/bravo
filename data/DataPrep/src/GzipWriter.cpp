#include "GzipWriter.h"

const unsigned int GzipWriter::DEFAULT_BUFFER_SIZE = 33554432u;

GzipWriter::GzipWriter(unsigned int buffer_size) throw (runtime_error) :
		buffer(nullptr), max_string_length(0u), gzfile(nullptr) {
	if (strcmp(zlibVersion(), ZLIB_VERSION) != 0) {
		throw runtime_error("Incompatible ZLIB version!");
	}

	try {
		buffer = new char[buffer_size];
		max_string_length = buffer_size - 1u;
	} catch (bad_alloc& e) {
		throw runtime_error("Error in memory allocation for gzip writer!");
	}
}

GzipWriter::~GzipWriter() {
	delete[] buffer;
	buffer = nullptr;

	if (gzfile != nullptr) {
		gzclose(gzfile);
		gzfile = nullptr;
	}
}

void GzipWriter::open(const char* file_name) throw (runtime_error) {
	if (gzfile == nullptr) {
		gzfile = gzopen(file_name, "wb");
		if (gzfile == nullptr) {
			throw runtime_error("Error while opening file for writing in gzip writer!");
		}
	}
}

void GzipWriter::close() throw (runtime_error) {
	if (gzfile != nullptr) {
		int gzerrno = 0;
		gzerrno = gzclose(gzfile);
		if (gzerrno != Z_OK) {
			throw runtime_error("Error while closing file in gzip writer!");
		}
		gzfile = nullptr;
	}
}

void GzipWriter::write(const char* format, ...) throw (runtime_error) {
	va_list arguments;
	long int n = 0;

	va_start(arguments, format);
	if ((n = vsnprintf(buffer, max_string_length, format, arguments)) < 0) {
		throw runtime_error("Error while writing to memory buffer in gzip writer!");
	} else if (n > max_string_length) {
		throw runtime_error("Too small buffer size in gzip writer!");
	}
	va_end(arguments);

	if (gzputs(gzfile, buffer) < 0) {
		throw runtime_error("Error while writing to file in gzip writer!");
	}
}
