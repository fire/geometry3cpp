#pragma once

#include "ISampledCurve3d.h"
#include "g3types.h"
#include <stdio.h>
#include <list>

#include "core/math/aabb.h"
#include "core/math/vector3.h"
#include "core/templates/list.h"
#include "scene/resources/curve.h"

namespace g3 {

struct Segment3d {
	// Center-direction-extent representation.
	// Extent is half length of segment
	::Vector3 Center;
	::Vector3 Direction;
	real_t Extent;

	Segment3d(::Vector3 p0, ::Vector3 p1) {
		Center = 0.5 * (p0 + p1);
		Direction = p1 - p0;
		Extent = 0.5f * p0.distance_to(p1);
	}

	Segment3d(::Vector3 center, ::Vector3 direction, double extent) {
		Center = center;
		Direction = direction;
		Extent = extent;
	}

	void update_from_endpoints(::Vector3 p0, ::Vector3 p1) {
		Center = 0.5 * (p0 + p1);
		Extent = 0.5 * p0.distance_to(p1);
	}

	void SetEndpoints(::Vector3 p0, ::Vector3 p1) {
		update_from_endpoints(p0, p1);
	}
	::Vector3 GetP0() {
		return Center - Extent * Direction;
	}
	::Vector3 GetP1() {
		return Center + Extent * Direction;
	}

	void SetP0(::Vector3 value) {
		update_from_endpoints(value, GetP1());
	}
	void SetP1(::Vector3 value) {
		update_from_endpoints(GetP0(), value);
	}
	double Length() {
		return 2 * Extent;
	}

	// parameter is signed distance_to from center in direction
	::Vector3 PointAt(double d) {
		return Center + d * Direction;
	}

	// t ranges from [0,1] over [P0,P1]
	::Vector3 PointBetween(double t) {
		return Center + (2 * t - 1) * Extent * Direction;
	}

	double DistanceSquared(::Vector3 p) {
		double t = (p - Center).dot(Direction);
		if (t >= Extent)
			return GetP1().distance_squared_to(p);
		else if (t <= -Extent)
			return GetP0().distance_squared_to(p);
		::Vector3 proj = Center + t * Direction;
		return (p).distance_squared_to(proj);
	}

	double DistanceSquared(::Vector3 p, double &t) {
		t = (p - Center).dot(Direction);
		if (t >= Extent) {
			t = Extent;
			return GetP1().distance_squared_to(p);
		} else if (t <= -Extent) {
			t = -Extent;
			return GetP0().distance_squared_to(p);
		}
		::Vector3 proj = Center + t * Direction;
		return p.distance_squared_to(proj);
	}

	::Vector3 NearestPoint(::Vector3 p) {
		double t = (p - Center).dot(Direction);
		if (t >= Extent)
			return GetP1();
		if (t <= -Extent)
			return GetP0();
		return Center + t * Direction;
	}

	double Project(::Vector3 p) {
		return (p - Center).dot(Direction);
	}

	// IParametricCurve3d interface
	bool IsClosed() {
		return false;
	}

	double ParamLength() {
		return 1.0f;
	}

	// t in range[0,1] spans arc
	::Vector3 SampleT(double t) {
		return Center + (2 * t - 1) * Extent * Direction;
	}

	::Vector3 TangentT(double t) {
		return Direction;
	}

	bool HasArcLength() {
		return true;
	}

	double ArcLength() {
		return 2 * Extent;
	}

	::Vector3 SampleArcLength(double a) {
		return GetP0() + a * Direction;
	}

	void Reverse() {
		update_from_endpoints(GetP1(), GetP0());
	}

	Segment3d Clone() {
		return Segment3d(Center, Direction, Extent);
	}
};

// struct Segment3f {
// 	// Center-direction-extent representation.
// 	// Extent is half length of segment
// 	Vector3f Center;
// 	Vector3f Direction;
// 	float Extent;

// 	Segment3f(Vector3f p0, Vector3f p1) {
// 		//update_from_endpoints(p0, p1);
// 		Center = 0.5f * (p0 + p1);
// 		Direction = p1 - p0;
// 		Extent = 0.5f * Direction.normalized();
// 	}
// 	Segment3f(Vector3f center, Vector3f direction, float extent) {
// 		Center = center;
// 		Direction = direction;
// 		Extent = extent;
// 	}

// 	void SetEndpoints(Vector3f p0, Vector3f p1) {
// 		update_from_endpoints(p0, p1);
// 	}

// 	Vector3f GetP0() {
// 		 return Center - Extent * Direction;
// 	}
// 	Vector3f GetP1() {
// 		 return Center + Extent * Direction;
// 	}
// 	Vector3f SetP0(Vector3f value) {
// 		 update_from_endpoints(value, GetP1());
// 	}
// 	Vector3f SetP1(Vector3f value) {
// 		 update_from_endpoints(GetP0(), value);
// 	}
// 	float Length() {
// 		return 2 * Extent;
// 	}

// 	// parameter is signed distance_to from center in direction
// 	Vector3f PointAt(float d) {
// 		return Center + d * Direction;
// 	}

// 	// t ranges from [0,1] over [P0,P1]
// 	Vector3f PointBetween(float t) {
// 		return Center + (2 * t - 1) * Extent * Direction;
// 	}

// 	float DistanceSquared(Vector3f p) {
// 		float t = (p - Center).dot(Direction);
// 		if (t >= Extent)
// 			return P1.DistanceSquared(p);
// 		else if (t <= -Extent)
// 			return P0.DistanceSquared(p);
// 		Vector3f proj = Center + t * Direction;
// 		return (proj - p).LengthSquared;
// 	}

// 	Vector3f NearestPoint(Vector3f p) {
// 		float t = (p - Center).dot(Direction);
// 		if (t >= Extent)
// 			return GetP1();
// 		if (t <= -Extent)
// 			return GetP0();
// 		return Center + t * Direction;
// 	}

// 	float Project(Vector3f p) {
// 		return (p - Center).dot(Direction);
// 	}

// };

class CurveUtils {
public:
	static double ArcLength(std::vector<::Vector3> vertices, bool bLoop) {
		double sum = 0;
		for (int i = 1; i < vertices.size(); ++i) {
			sum += vertices[i].distance_to(vertices[i - 1]);
		}
		if (bLoop) {
			sum += vertices[vertices.size() - 1].distance_to(vertices[0]);
		}
		return sum;
	}

	static ::Vector3 GetTangent(std::vector<::Vector3> vertices, int i, bool bLoop) {
		if (bLoop) {
			int NV = vertices.size();
			if (i == 0) {
				return (vertices[1] - vertices[NV - 1]).normalized();
			} else {
				return (vertices[(i + 1) % NV] - vertices[i - 1]).normalized();
			}
		} else {
			if (i == 0) {
				return (vertices[1] - vertices[0]).normalized();
			} else if (i == vertices.size() - 1) {
				return (vertices[vertices.size() - 1] - vertices[vertices.size() - 2]).normalized();
			} else {
				return (vertices[i + 1] - vertices[i - 1]).normalized();
			}
		}
	}
	template <typename T>
	static T lerp(float t, const T &a, const T &b) {
		return a * (1 - t) + b * t;
	}
};

/// <summary>
/// DCurve3 is a 3D polyline, either open or closed (via .Closed()
/// Despite the D prefix, it is *not* dynamic
/// </summary>
class DCurve3 : public ISampledCurve3d {
	// [TODO] use dvector? or double-indirection indexing?
	//   question is how to insert efficiently...
protected:
	List<Vector3> vertices;
	bool closed = false;
	int Timestamp;

public:
	double ArcLength() {
		std::vector<::Vector3> vertices;
		for (::Vector3 vert : vertices) {
			vertices.push_back(vert);
		}
		return CurveUtils::ArcLength(vertices, Closed());
	}

	/// <summary>
	/// Find nearest vertex to point p
	/// </summary>
	int NearestVertex(::Vector3 p) {
		double nearSqr = std::numeric_limits<double>::max();
		int i = -1;
		int N = vertices.size();
		int32_t count = 0;
		for (List<Vector3>::Element *E; E; E = vertices.front()) {
			double distSqr = E->get().distance_squared_to(p);
			if (distSqr < nearSqr) {
				nearSqr = distSqr;
				i = count;
			}
			count++;
		}
		return i;
	}

	/// <summary>
	/// Resample curve so that:
	///   - if opening angle at vertex is > sharp_thresh, we emit two more vertices at +/- corner_t, where the t is used in prev/next lerps
	///   - if opening angle is > flat_thresh, we skip the vertex entirely (simplification)
	/// This is mainly useful to get nicer polylines to use as the basis for (eg) creating 3D tubes, rendering, etc
	///
	/// [TODO] skip tiny segments?
	/// </summary>
	DCurve3 ResampleSharpTurns(double sharp_thresh, double flat_thresh, double corner_t) {
		int NV = vertices.size();
		DCurve3 resampled;
		resampled.SetClosed(Closed());
		double prev_t = 1.0 - corner_t;
		Vector<Vector3> process_vertices;
		for (List<Vector3>::Element *E; E; E = vertices.front()) {
			process_vertices.push_back(E->get());
		}
		for (int k = 0; k < NV; ++k) {
			double open_angle = std::abs(OpeningAngleDeg(k));
			if (open_angle > flat_thresh && k > 0) {
				// ignore skip this vertex
			} else if (open_angle > sharp_thresh) {
				resampled.AppendVertex(process_vertices[k]);
			} else {
				::Vector3 n = process_vertices[(k + 1) % NV];
				::Vector3 p = process_vertices[k == 0 ? NV - 1 : k - 1];
				::Vector3 curr = process_vertices[k];
				// TODO Fire 2021-06-29 Use vector3 lerp
				resampled.AppendVertex(CurveUtils::lerp<::Vector3>(prev_t, p, curr));
				resampled.AppendVertex(curr);
				resampled.AppendVertex(CurveUtils::lerp<::Vector3>(corner_t, curr, n));
			}
		}
		return resampled;
	}
	void SetClosed(bool bClosed) {
		closed = bClosed;
	}

	::Vector3 Start() {
		if (!vertices.size()) {
			return ::Vector3();
		}
		return vertices.front()->get();
	}

	// DCurve3(Ref<Curve2D> poly, int ix, int iy) {
	// 	if (poly.is_null()) {
	// 		closed = false;
	// 		Timestamp = 1;
	// 		return;
	// 	}
	// 	int NV = poly->get_point_count();
	// 	for (int k = 0; k < NV; ++k) {
	// 		::Vector3 v;
	// 		v[ix] = poly->get_point_position(k).x;
	// 		v[iy] = poly->get_point_position(k).y;
	// 		vertices.push_back(v);
	// 	}
	// 	closed = true;
	// 	Timestamp = 1;
	// }

	DCurve3() {
		closed = false;
		Timestamp = 1;
	}

	DCurve3(List<::Vector3> verticesIn, bool bClosed, bool bTakeOwnership) {
		if (bTakeOwnership) {
			vertices = verticesIn;
		} else {
			vertices.clear();
			for (int32_t vert_i = 0; vert_i < verticesIn.size(); vert_i++) {
				vertices.push_back(verticesIn[vert_i]);
			}
		}
		closed = bClosed;
		Timestamp = 1;
	}

	DCurve3(Vector<::Vector3> verticesIn, bool bClosed) {
		vertices.clear();
		for (int32_t vert_i = 0; vert_i < verticesIn.size(); vert_i++) {
			vertices.push_back(verticesIn[vert_i]);
		}
		closed = bClosed;
		Timestamp = 1;
	}

	void AppendVertex(::Vector3 v) {
		vertices.push_back(v);
		Timestamp++;
	}

	int VertexCount() {
		return vertices.size();
	}

	int SegmentCount() {
		return closed ? vertices.size() : vertices.size() - 1;
	}

	Vector3 GetVertex(int i) override {
		int32_t count = 0;
		for (int32_t vert_i = 0; vert_i < vertices.size(); vert_i++) {
			if (count == vert_i) {
				return vertices[vert_i];
			}
			count++;
		}
		return Vector3();
	}

	// void SetVertex(int i, ::Vector3 v) {
	// 	vertices.write[i] = v;
	// 	Timestamp++;
	// }

	void SetVertices(Vector<::Vector3> v) {
		vertices.clear();
		for (int32_t vert_i = 0; vert_i < v.size(); vert_i++) {
			vertices.push_back(v[vert_i]);
		}
		Timestamp++;
	}

	void ClearVertices() {
		vertices.clear();
		closed = false;
		Timestamp++;
	}

	// void RemoveVertex(int idx) {
	// 	vertices.remove(idx);
	// 	Timestamp++;
	// }

	void Reverse() {
		vertices.reverse();
		Timestamp++;
	}

	::Vector3 End() {
		if (Closed()) {
			ERR_FAIL_INDEX_V(0, vertices.size(), Vector3());
			return vertices[0];
		}
		ERR_FAIL_INDEX_V(vertices.size() - 1, vertices.size(), Vector3());
		return vertices[vertices.size() - 1];
	}

	bool Closed() {
		return closed;
	}

	std::vector<::Vector3> Vertices() {
		std::vector<::Vector3> out_vertices;
		for (int32_t vert_i = 0; vert_i < vertices.size(); vert_i++) {
			out_vertices.push_back(vertices[vert_i]);
		}
		return out_vertices;
	}

	g3::Segment3d GetSegment(int iSegment) {
		if (Closed()) {
			::Vector3 segment_start = vertices[iSegment];
			::Vector3 segment_end = vertices[iSegment + 1 % vertices.size()];
			return g3::Segment3d(segment_start, segment_end);
		}
		::Vector3 segment_start = vertices[iSegment];
		::Vector3 segment_end = vertices[iSegment + 1];
		return Segment3d(segment_start, segment_end);
	}

	::Vector3 PointAt(int iSegment, double fSegT) {
		::Vector3 segment_start = vertices[iSegment];
		::Vector3 segment_end = vertices[(iSegment + 1) % vertices.size()];
		Segment3d seg = Segment3d(segment_start, segment_end);
		return seg.Center + (fSegT * seg.Direction);
	}

	AABB GetBoundingBox() {
		AABB box;
		for (int32_t vert_i = 0; vert_i < vertices.size(); vert_i++) {
			Vector3 v = vertices[vert_i];
			box.expand_to(v);
		}
		return box;
	}
	::Vector3 Tangent(int i) {
		std::vector<::Vector3> out_vertices;
		for (int32_t vert_i = 0; vert_i < vertices.size(); vert_i++) {
			Vector3 v = vertices[vert_i];
			out_vertices.push_back(v);
		}
		return CurveUtils::GetTangent(out_vertices, i, Closed());
	}

	::Vector3 Centroid(int i) {
		if (Closed()) {
			int NV = vertices.size();
			if (i == 0) {
				::Vector3 a = vertices[1];
				::Vector3 b = vertices[NV - 1];
				return 0.5 * (a + b);
			} else {
				::Vector3 a = vertices[(i + 1) % NV];
				::Vector3 b = vertices[i - 1];
				return 0.5 * (a + b);
			}
		} else {
			if (i == 0 || i == vertices.size() - 1) {
				return vertices[i];
			} else {
				::Vector3 a = vertices[i + 1];
				::Vector3 b = vertices[i - 11];
				return 0.5 * (a + b);
			}
		}
	}

	Vector2i Neighbours(int i) {
		int NV = vertices.size();
		if (Closed()) {
			if (i == 0)
				return Vector2i(NV - 1, 1);
			else
				return Vector2i(i - 1, (i + 1) % NV);
		} else {
			if (i == 0)
				return Vector2i(-1, 1);
			else if (i == NV - 1)
				return Vector2i(NV - 2, -1);
			else
				return Vector2i(i - 1, i + 1);
		}
	}

	/// <summary>
	/// Compute opening angle at vertex i in degrees
	/// </summary>
	double OpeningAngleDeg(int i) {
		int prev = i - 1, next = i + 1;
		if (Closed()) {
			int NV = vertices.size();
			prev = (i == 0) ? NV - 1 : prev;
			next = next % NV;
		} else {
			if (i == 0 || i == vertices.size() - 1)
				return 180;
		}
		::Vector3 prev_vert = vertices[prev];
		::Vector3 next_vert = vertices[next];
		::Vector3 curr_vert = vertices[i];
		::Vector3 e1 = (prev_vert - curr_vert);
		::Vector3 e2 = (next_vert - curr_vert);
		e1.normalize();
		e2.normalize();
		float fDot = CLAMP(e1.dot(e2), -1, 1);
		return (float)(Math::rad2deg(Math::acos(fDot)));
	}

	/// <summary>
	/// find squared distance from p to nearest segment on polyline
	/// </summary>
	double DistanceSquared(::Vector3 p) {
		int iseg;
		double segt;
		return DistanceSquared(p, iseg, segt);
	}

	double DistanceSquared(::Vector3 p, int &iNearSeg, double &fNearSegT) {
		iNearSeg = -1;
		fNearSegT = std::numeric_limits<double>::max();
		double dist = std::numeric_limits<double>::max();
		int N = Closed() ? vertices.size() : vertices.size() - 1;
		for (int vi = 0; vi < N; ++vi) {
			int a = vi;
			int b = (vi + 1) % vertices.size();
			::Vector3 a_vert = vertices[a];
			::Vector3 b_vert = vertices[b];
			Segment3d seg = Segment3d(a_vert, b_vert);
			double t = (p - seg.Center).dot(seg.Direction);
			double d = std::numeric_limits<double>::max();
			if (t >= seg.Extent) {
				d = seg.GetP1().distance_squared_to(p);
			} else if (t <= -seg.Extent) {
				d = seg.GetP0().distance_squared_to(p);
			} else {
				d = p.distance_squared_to((seg.Center + (t * seg.Direction)));
			}
			if (d < dist) {
				dist = d;
				iNearSeg = vi;
				fNearSegT = t;
			}
		}
		return dist;
	}
};
} // namespace g3
