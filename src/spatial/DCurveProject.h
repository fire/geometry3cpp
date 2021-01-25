#pragma once

#include "src/geometry/DCurve3.h"
#include "src/spatial/SpatialInterfaces.h"

#include "g3types.h"

#include <limits>
#include <string>

namespace g3 {
	class DCurveProjectionTarget : public IProjectionTarget {
	public:
		DCurve3Ptr Curve;

		DCurveProjectionTarget(DCurve3Ptr curve) {
			Curve = curve;
		}

		virtual Vector3d Project(const Vector3d &vPoint, int identifier = -1) override {
			Vector3d vNearest;
			double fNearestSqr = std::numeric_limits<double>::max();

			int N = Curve->VertexCount();
			int NStop = (Curve->Closed()) ? N : N - 1;
			for (int i = 0; i < NStop; ++i) {
				Segment3d seg = Segment3d(Curve->GetVertex(i), Curve->GetVertex((i + 1) % N));
				Vector3d pt = seg.NearestPoint(vPoint);
				double dsqr = (vPoint - pt).squaredNorm();
				if (dsqr < fNearestSqr) {
					fNearestSqr = dsqr;
					vNearest = pt;
				}
			}

			return (fNearestSqr < std::numeric_limits<double>::max()) ? vNearest : vPoint;
		}
	};
}
