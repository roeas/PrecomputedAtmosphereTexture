#pragma once

#include <cassert>

#define IN(x) const x&
#define OUT(x) x&
#define TEMPLATE(x) template<class x>
#define TEMPLATE_ARGUMENT(x) <x>

float clamp(const float num, const float min, const float max) {
	assert(min <= max);
	return (num < min) ? (min) : ((num > max) ? (max) : (num));
}

Number ClampCosine(Number mu) {
	return clamp(mu, Number(-1.0), Number(1.0));
}
Length ClampDistance(Length d) {
	return std::max(d, Length(0.0 * m));
}
Length ClampRadius(IN(AtmosphereParameters) atmosphere, Length r) {
	return clamp(r, atmosphere.bottom_radius, atmosphere.top_radius);
}
Length SafeSqrt(Area a) {
	return sqrt(std::max(a, Length(0.0 * m2)));
}
