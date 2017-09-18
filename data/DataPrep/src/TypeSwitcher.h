#ifndef TYPESWITCHER_H_
#define TYPESWITCHER_H_

#include <iostream>
#include <functional>
#include <vector>
#include "vcf.h"
#include "hts.h"

using namespace std;

class TypeSwitcher {
private:
	template<typename Tbcf, Tbcf (*convert)(const uint8_t*), Tbcf missing>
	void get_read(vector<int32_t>& v);

public:
	bcf_fmt_t* fmt;
	uint8_t* fmt_p;

	void (TypeSwitcher::*read)(vector<int32_t>&);

	TypeSwitcher();
	virtual ~TypeSwitcher() noexcept;

	void init(bcf_fmt_t* fmt) throw (runtime_error);
};

#endif
