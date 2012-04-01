#pragma once
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#pragma unmanaged

struct CloudPoint
{
public:
	double x, y, z;
};


struct Cloud
{
public:
	vector<CloudPoint> positions;

	void Save(char* filename)
	{
		ofstream file(filename, ios::binary); 
		if(!file) { 
			return; 
		} 
		unsigned int size = positions.size();
		file.write((char*)&size, sizeof(size));
		file.write((char*)&positions[0], sizeof(CloudPoint) * size);

		file.close();	
	}
};