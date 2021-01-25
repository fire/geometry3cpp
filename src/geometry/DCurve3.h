#pragma once

#include "ISampledCurve3d.h"
#include "Wm5Segment3.h"
#include "g3types.h"
#include <stdio.h>
#include <list>

namespace g3 {

struct Segment3d {
	// Center-direction-extent representation.
	// Extent is half length of segment
	Vector3d Center;
	Vector3d Direction;
	double Extent;

	Segment3d(Vector3d p0, Vector3d p1) {
		//update_from_endpoints(p0, p1);
		Center = 0.5 * (p0 + p1);
		Direction = p1 - p0;
		Extent = 0.5 * Direction.normalize();
	}
	Segment3d(Vector3d center, Vector3d direction, double extent) {
		Center = center;
		Direction = direction;
		Extent = extent;
	}

	void update_from_endpoints(Vector3f p0, Vector3f p1) {
		Center = 0.5f * (p0 + p1);
		Direction = p1 - p0;
		Extent = 0.5f * Direction.normalize();
	}

	void SetEndpoints(Vector3d p0, Vector3d p1) {
		update_from_endpoints(p0, p1);
	}
	Vector3d GetP0() {
		return Center - Extent * Direction;
	}
	Vector3d GetP1() {
		return Center + Extent * Direction;
	}

	Vector3d SetP0(Vector3d value) {
		update_from_endpoints(value, GetP1());
	}
	Vector3d SetP1(Vector3d value) {
		update_from_endpoints(GetP0(), value);
	}
	double Length() {
		return 2 * Extent;
	}

	// parameter is signed distance from center in direction
	Vector3d PointAt(double d) {
		return Center + d * Direction;
	}

	// t ranges from [0,1] over [P0,P1]
	Vector3d PointBetween(double t) {
		return Center + (2 * t - 1) * Extent * Direction;
	}

	double DistanceSquared(Vector3d p) {
		double t = (p - Center).Dot(Direction);
		if (t >= Extent)
			return GetP1().DistanceSquared(p);
		else if (t <= -Extent)
			return GetP0().DistanceSquared(p);
		Vector3d proj = Center + t * Direction;
		return (proj - p).LengthSquared;
	}

	double DistanceSquared(Vector3d p, double &t) {
		t = (p - Center).dot(Direction);
		if (t >= Extent) {
			t = Extent;
			return GetP1().DistanceSquared(p);
		} else if (t <= -Extent) {
			t = -Extent;
			return GetP0().DistanceSquared(p);
		}
		Vector3d proj = Center + t * Direction;
		return (proj - p).LengthSquared;
	}

	Vector3d NearestPoint(Vector3d p) {
		double t = (p - Center).dot(Direction);
		if (t >= Extent)
			return GetP1();
		if (t <= -Extent)
			return GetP0();
		return Center + t * Direction;
	}

	double Project(Vector3d p) {
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
	Vector3d SampleT(double t) {
		return Center + (2 * t - 1) * Extent * Direction;
	}

	Vector3d TangentT(double t) {
		return Direction;
	}

	bool HasArcLength {
		get { return true; }
	}

	double ArcLength {
		get { return 2 * Extent; }
	}

	Vector3d SampleArcLength(double a) {
		return P0 + a * Direction;
	}

	void Reverse() {
		update_from_endpoints(P1, P0);
	}

	IParametricCurve3d Clone() {
		return new Segment3d(this.Center, this.Direction, this.Extent);
	}

};

struct Segment3f {
	// Center-direction-extent representation.
	// Extent is half length of segment
	Vector3f Center;
	Vector3f Direction;
	float Extent;

	Segment3f(Vector3f p0, Vector3f p1) {
		//update_from_endpoints(p0, p1);
		Center = 0.5f * (p0 + p1);
		Direction = p1 - p0;
		Extent = 0.5f * Direction.normalized();
	}
	Segment3f(Vector3f center, Vector3f direction, float extent) {
		Center = center;
		Direction = direction;
		Extent = extent;
	}

	void SetEndpoints(Vector3f p0, Vector3f p1) {
		update_from_endpoints(p0, p1);
	}

	Vector3f GetP0() {
		 return Center - Extent * Direction; 
	}
	Vector3f GetP1() {
		 return Center + Extent * Direction; 
	}
	Vector3f SetP0(Vector3f value) {
		 update_from_endpoints(value, GetP1()); 
	}
	Vector3f SetP1(Vector3f value) {
		 update_from_endpoints(GetP0(), value); 
	}
	float Length() {
		return 2 * Extent;
	}

	// parameter is signed distance from center in direction
	Vector3f PointAt(float d) {
		return Center + d * Direction;
	}

	// t ranges from [0,1] over [P0,P1]
	Vector3f PointBetween(float t) {
		return Center + (2 * t - 1) * Extent * Direction;
	}

	float DistanceSquared(Vector3f p) {
		float t = (p - Center).dot(Direction);
		if (t >= Extent)
			return P1.DistanceSquared(p);
		else if (t <= -Extent)
			return P0.DistanceSquared(p);
		Vector3f proj = Center + t * Direction;
		return (proj - p).LengthSquared;
	}

	Vector3f NearestPoint(Vector3f p) {
		float t = (p - Center).Dot(Direction);
		if (t >= Extent)
			return P1;
		if (t <= -Extent)
			return P0;
		return Center + t * Direction;
	}

	float Project(Vector3f p) {
		return (p - Center).Dot(Direction);
	}

}

class CurveUtils {
public:
	static double ArcLength(std::vector<Vector3d> vertices, bool bLoop) {
		double sum = 0;
		for (int i = 1; i < vertices.size(); ++i) {
			sum += (vertices[i - 1] - vertices[i]).norm();
		}
		if (bLoop) {
			sum += (vertices[0] - vertices[vertices.size() - 1]).norm();
		}
		return sum;
	}

	static Vector3d GetTangent(std::vector<Wm5::Vector3d> vertices, int i, bool bLoop) {
		if (bLoop) {
			int NV = vertices.size();
			if (i == 0) {
				return (vertices[1] - vertices[NV - 1]).Normalized();
			} else {
				return (vertices[(i + 1) % NV] - vertices[i - 1]).Normalized();
			}
		} else {
			if (i == 0) {
				return (vertices[1] - vertices[0]).Normalized();
			} else if (i == vertices.size() - 1) {
				return (vertices[vertices.size() - 1] - vertices[vertices.size() - 2]).Normalized();
			} else {
				return (vertices[i + 1] - vertices[i - 1]).Normalized();
			}
		}
	}
};

class DCurve3;
typedef std::shared_ptr<DCurve3> DCurve3Ptr;

/// <summary>
/// DCurve3 is a 3D polyline, either open or closed (via .Closed()
/// Despite the D prefix, it is *not* dynamic
/// </summary>
class DCurve3 : public ISampledCurve3d {
	// [TODO] use dvector? or double-indirection indexing?
	//   question is how to insert efficiently...
protected:
	std::list<Vector3d> vertices;
	bool closed = false;
	int Timestamp;

public:
	double ArcLength() {
		std::vector<Vector3d> vertices;
		for (Vector3d vert : vertices) {
			vertices.push_back(vert);
		}
		return CurveUtils::ArcLength(vertices, Closed());
	}

	/// <summary>
	/// Find nearest vertex to point p
	/// </summary>
	int NearestVertex(Vector3d p) {
		double nearSqr = std::numeric_limits<double>::max();
		int i = -1;
		int N = vertices.size();
		for (int vi = 0; vi < N; ++vi) {
			double distSqr = (p - vertices[vi]).squaredNorm();
			if (distSqr < nearSqr) {
				nearSqr = distSqr;
				i = vi;
			}
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
		for (int k = 0; k < NV; ++k) {
			double open_angle = Wm5::Math::FAbs(OpeningAngleDeg(k));
			if (open_angle > flat_thresh && k > 0) {
				// ignore skip this vertex
			} else if (open_angle > sharp_thresh) {
				resampled.AppendVertex(vertices[k]);
			} else {
				Vector3d n = vertices[(k + 1) % NV];
				Vector3d p = vertices[k == 0 ? NV - 1 : k - 1];
				resampled.AppendVertex(Vector3d::Lerp(p, vertices[k], prev_t));
				resampled.AppendVertex(vertices[k]);
				resampled.AppendVertex(Vector3d::Lerp(vertices[k], n, corner_t));
			}
		}
		return resampled;
	}
	void SetClosed(bool bClosed) {
		closed = bClosed;
	}

	Vector3d Start() {
		if (!vertices.size()) {
			return Vector3d();
		}
		return *vertices.begin();
	}

	DCurve3(Wm5::Polygon2d poly, int ix, int iy) {
		int NV = poly.VertexCount();
		vertices.resize(NV);
		for (int k = 0; k < NV; ++k) {
			Vector3d v = Vector3d::Zero;
			v[ix] = poly.GetVertex(k).X();
			v[iy] = poly.GetVertex(k).Y();
			vertices.push_back(v);
		}
		closed = true;
		Timestamp = 1;
	}

	DCurve3() {
		closed = false;
		Timestamp = 1;
	}

	DCurve3(std::list<Vector3d> verticesIn, bool bClosed, bool bTakeOwnership) {
		if (bTakeOwnership) {
			vertices = verticesIn;
		} else {
			vertices.clear();
			for (Vector3d vec : verticesIn) {
				vertices.push_back(vec);
			}
		}
		closed = bClosed;
		Timestamp = 1;
	}

	DCurve3(std::vector<Vector3d> verticesIn, bool bClosed) {
		vertices.clear();
		for (Vector3d vert : verticesIn) {
			vertices.push_back(vert);
		}
		closed = bClosed;
		Timestamp = 1;
	}

	void AppendVertex(Vector3d v) {
		vertices.push_back(v);
		Timestamp++;
	}

	int VertexCount() {
		return vertices.size();
	}

	int SegmentCount() {
		return closed ? vertices.size() : vertices.size() - 1;
	}

	Vector3d GetVertex(int i) {
		int32_t count = 0;
		for (Vector3d vert : vertices) {
			if (count == i) {
				return vert;
			}
			count++;
		}
		return Vector3d();
	}

	void SetVertex(int i, Vector3d v) {
		std::list<Vector3d>::iterator it = vertices.begin();
		std::advance(it, i);
		*it = v;
		Timestamp++;
	}

	void SetVertices(std::vector<Vector3d> v) {
		vertices.clear();
		for (Vector3d vert : v) {
			vertices.push_back(vert);
		}
		Timestamp++;
	}

	void ClearVertices() {
		vertices.clear();
		closed = false;
		Timestamp++;
	}

	void RemoveVertex(int idx) {
		std::list<Vector3d>::iterator it = vertices.begin();
		std::advance(it, idx);
		vertices.erase(it);
		Timestamp++;
	}

	void Reverse() {
		std::reverse(vertices.begin(), vertices.end());
		Timestamp++;
	}

	Vector3d End() {
		if (Closed()) {
			std::list<Vector3d>::iterator it = vertices.begin();
			std::advance(it, 0);
			return *it;
		}
		std::list<Vector3d>::iterator it = vertices.end();
		return *it;
	}

	bool Closed() {
		return closed;
	}

	std::vector<Vector3d> Vertices() {
		std::vector<Vector3d> out_vertices;
		for (Vector3d vert : vertices) {
			out_vertices.push_back(vert);
		}
		return out_vertices;
	}

	Wm5::Segment3d GetSegment(int iSegment) {
		if (Closed()) {
			std::list<Vector3d>::iterator it = vertices.begin();
			std::advance(it, iSegment);
			Vector3d segment_start = *it;
			it = vertices.begin();
			std::advance(it, iSegment + 1 % vertices.size());
			Vector3d segment_end = *it;
			return Wm5::Segment3d(segment_start, segment_end);
		}
		std::list<Vector3d>::iterator it = vertices.begin();
		std::advance(it, iSegment);
		Vector3d segment_start = *it;
		it = vertices.begin();
		std::advance(it, iSegment + 1);
		Vector3d segment_end = *it;
		return Wm5::Segment3d(segment_start, segment_end);
	}

	Vector3d PointAt(int iSegment, double fSegT) {
		std::list<Vector3d>::iterator it = vertices.begin();
		std::advance(it, iSegment);
		Vector3d segment_start = *it;
		it = vertices.begin();
		std::advance(it, (iSegment + 1) % vertices.size());
		Vector3d segment_end = *it;
		Wm5::Segment3d seg = Wm5::Segment3d(segment_start, segment_end);
		return seg.Center + (fSegT * seg.Direction);
	}

	AxisAlignedBox3d GetBoundingBox() {
		AxisAlignedBox3d box;
		for (Vector3d v : vertices) {
			box.Contain(v);
		}
		return box;
	}
	Vector3d Tangent(int i) {
		std::vector<Vector3d> out_vertices;
		for (Vector3d vert : vertices) {
			out_vertices.push_back(vert);
		}
		return CurveUtils::GetTangent(out_vertices, i, Closed());
	}

	Vector3d Centroid(int i) {
		if (Closed()) {
			int NV = vertices.size();
			if (i == 0) {
				std::list<Vector3d>::iterator it = vertices.begin();
				std::advance(it, 1);
				Vector3d a = *it;
				it = vertices.begin();
				std::advance(it, NV - 1);
				Vector3d b = *it;
				return 0.5 * (a + b);
			} else {
				std::list<Vector3d>::iterator it = vertices.begin();
				std::advance(it, (i + 1) % NV);
				Vector3d a = *it;
				it = vertices.begin();
				std::advance(it, i - 1);
				Vector3d b = *it;
				return 0.5 * (a + b);
			}
		} else {
			if (i == 0 || i == vertices.size() - 1) {
				std::list<Vector3d>::iterator it = vertices.begin();
				std::advance(it, i);
				return *it;
			} else {
				std::list<Vector3d>::iterator it = vertices.begin();
				std::advance(it, i + 1);
				Vector3d a = *it;
				it = vertices.begin();
				std::advance(it, i - 1);
				Vector3d b = *it;
				return 0.5 * (a + b);
			}
		}
	}

	Index2i Neighbours(int i) {
		int NV = vertices.size();
		if (Closed()) {
			if (i == 0)
				return Index2i(NV - 1, 1);
			else
				return Index2i(i - 1, (i + 1) % NV);
		} else {
			if (i == 0)
				return Index2i(-1, 1);
			else if (i == NV - 1)
				return Index2i(NV - 2, -1);
			else
				return Index2i(i - 1, i + 1);
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
		std::list<Vector3d>::iterator it = vertices.begin();
		std::advance(it, prev);
		Vector3d prev_vert = *it;
		it = vertices.begin();
		std::advance(it, next);
		Vector3d next_vert = *it;
		it = vertices.begin();
		std::advance(it, i);
		Vector3d curr_vert = *it;
		Vector3d e1 = (prev_vert - curr_vert);
		Vector3d e2 = (next_vert - curr_vert);
		e1.normalize();
		e2.normalize();
		float fDot = Mathd::Clamp(e1.dot(e2), -1, 1);
		return (float)(Mathd::ACos(fDot) * Mathd::RAD_TO_DEG);
	}

	/// <summary>
	/// find squared distance from p to nearest segment on polyline
	/// </summary>
	double DistanceSquared(Vector3d p) {
		int iseg;
		double segt;
		return DistanceSquared(p, iseg, segt);
	}

	double DistanceSquared(Vector3d p, int &iNearSeg, double &fNearSegT) {
		iNearSeg = -1;
		fNearSegT = std::numeric_limits<double>::max();
		double dist = std::numeric_limits<double>::max();
		int N = Closed() ? vertices.size() : vertices.size() - 1;
		for (int vi = 0; vi < N; ++vi) {
			int a = vi;
			int b = (vi + 1) % vertices.size();
			std::list<Vector3d>::iterator it = vertices.begin();
			std::advance(it, a);
			Vector3d a_vert = *it;
			it = vertices.begin();
			std::advance(it, b);
			Vector3d b_vert = *it;
			Wm5::Segment3d seg = Wm5::Segment3d(a_vert, b_vert);
			double t = (p - seg.Center).dot(seg.Direction);
			double d = std::numeric_limits<double>::max();
			if (t >= seg.Extent) {
				d = (p - seg.P1).squaredNorm();
			} else if (t <= -seg.Extent) {
				d = (p - seg.P0).squaredNorm();
			} else {
				d = ((seg.Center + (t * seg.Direction)) - p).SquaredLength();
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
