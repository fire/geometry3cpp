#pragma once

#include "g3types.h"
#include <list>

class ISampledCurve3d {
private:
	std::list<Wm5::Vector3d> Vertices;

public:
	virtual int VertexCount() = 0;
	virtual int SegmentCount() = 0;
	virtual bool Closed() = 0;

	virtual Wm5::Vector3d GetVertex(int i) = 0;
	virtual Wm5::Segment3d GetSegment(int i) = 0;
};




