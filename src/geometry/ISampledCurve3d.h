#pragma once

#include "DCurve3.h"
#include "g3types.h"
#include "core/math/vector3.h"
#include "core/templates/list.h"

namespace g3 {
struct Segment3d;
class ISampledCurve3d {
private:
	List<Vector3> Vertices;

public:
	virtual int VertexCount() = 0;
	virtual int SegmentCount() = 0;
	virtual bool Closed() = 0;

	virtual Vector3 GetVertex(int i) = 0;
	virtual g3::Segment3d GetSegment(int i) = 0;
};

} // namespace g3