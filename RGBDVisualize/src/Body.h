#pragma once

#pragma unmanaged

struct Body
{
public:
	double x, y, z;
	double vx, vy, vz;
	double mass;
	bool isDisabled;
};