#include "TypeSwitcher.h"

TypeSwitcher::TypeSwitcher() : fmt(nullptr), fmt_p(nullptr), read(nullptr) {

}

TypeSwitcher::~TypeSwitcher() noexcept {

}

template<typename Tbcf, Tbcf (*convert)(const uint8_t*), Tbcf missing>
void TypeSwitcher::get_read(vector<int32_t>& v) {
	v.clear();
	for (int j = 0; j < fmt->n; ++j) {
		Tbcf p = convert(fmt_p + j * sizeof(Tbcf));
		if (p == missing) {
			v.push_back(bcf_int32_missing);
		} else {
			v.push_back((int32_t)p);
		}
	}
	fmt_p += fmt->size;
}

void TypeSwitcher::init(bcf_fmt_t* fmt) throw (runtime_error) {
	this->fmt = fmt;
	this->fmt_p = fmt->p;
	switch (fmt->type) {
		case BCF_BT_INT8:
			read = &TypeSwitcher::get_read<int8_t, le_to_i8, bcf_int8_missing>;
			break;
		case BCF_BT_INT16:
			read = &TypeSwitcher::get_read<int16_t, le_to_i16, bcf_int16_missing>;
			break;
		case BCF_BT_INT32:
			read = &TypeSwitcher::get_read<int32_t, le_to_i32, bcf_int32_missing>;
			break;
		default:
			throw runtime_error("Error while resolving data types in FORMAT!");
	}
}

