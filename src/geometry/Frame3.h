#pragma once

#include <g3types.h>

#include <VectorUtil.h>

#include "core/math/transform_3d.h"

/*
* Frame3D represents an orthonormal 3D coordinate frame, ie 3 orthogonal axes and an origin point.
* Although it is possible to pack all that information into a Transform, or a Basis+Vector3, it
* is not very convenient to work with. Making the axes available explicitly allows for much
* clearer code.
*
* Frame3D also provides functions to:
*    - manipulate the frame by (eg) setting one axis to point in a particular direction
*    - transform vectors and frames to/from frame-relative coordinates
*    - (TODO: ray-intersection, project to plane, ...)
*
* Internally, Frame3D stores the inverse of the rotation matrix that would take the XYZ axes
* to the Frame axes. This is so that the axes are the rows of the matrix, and hence in
* row-major storage the 3 floats of each axis are contiguous.
*
*/
namespace g3 {

class Frame3D : public Transform3D {
public:
	//! create a reference frame at a specific origin
	Frame3D(const Vector3 &origin = Vector3()) {
		set_origin(origin);
	}

	//! create an orthogonal frame at an origin with a given z axis vector
	Frame3D(const Vector3 &vOrigin,
			const Vector3 &vNormalizedAxis, int nAxis,
			bool bUseFastPerpVectors) {
		origin = vOrigin;
		if (bUseFastPerpVectors && nAxis == 2) {
			Vector3 kTangent0, kTangent1;
			ComputePerpVectors(vNormalizedAxis, kTangent0, kTangent1, true);
			basis.get_axis(0) = kTangent0;
			basis.get_axis(1) = kTangent1;
			basis.get_axis(2) = vNormalizedAxis;
		} else {
			AlignAxis(nAxis, vNormalizedAxis);
		}
	}

	//! copy a reference frame
	Frame3D(const Frame3D &copy) :
			origin(copy.origin), elements(copy.get_basis()) {}

	~Frame3D();

	//! get an axis of the frame
	Vector3 Axis(unsigned int nAxis) const {
		return get_basis().get_axis(nAxis);
	}

	//! access the axes of the frame
	Vector3 X() const {
		return get_basis().get_axis(0);
	}
	Vector3 Y() const {
		return get_basis().get_axis(1);
	}
	Vector3 Z() const {
		return get_basis().get_axis(2);
	}

	void AlignAxis(int nAxis, const Vector3 &toAxis, bool bNormalized) {
		Basis matAlign;
		ComputeAlignAxisMatrix(
				Axis(nAxis), (bNormalized) ? toAxis : toAxis.normalized(), matAlign);
		Rotate(matAlign, false);
	}

	//! matrix that rotates canonical/unit XYZ axes to frame axes
	Basis GetRotation() const {
		return get_basis().transpose();
	}

	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	void SetToWorldCoords(Vector3 &v) const {
		v = get_origin() + get_basis().transpose().xform(v);
	}
	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	Vector3 ToWorldCoords(const Vector3 &v) const {
		return get_origin() + get_basis().transposed().xform(v);
	}

	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	void SetToAxisCoords(Vector3 &v) const {
		v = get_basis().xform(v - origin);
	}
	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	Vector3 ToAxisCoords(const Vector3 &v) const {
		return get_basis().xform(v - origin);
	}

	//! treat f as being relative to this frame, transform to absolute/world frame
	Frame3D ToWorldCoords(const Frame3D &f, bool bPreserveOrigin = false) const {
		Vector3 x = get_basis().transposed().xform(f.X());
		Vector3 y = get_basis().transposed().xform(f.Y());
		Vector3 z = get_basis().transposed().xform(f.Z());
		Frame3D l;
		l.origin = ((bPreserveOrigin) ? f.origin : ToWorldCoords(f.origin));
		l.SetFrame(x, y, z);
		return l;
	}

	//! treat f as being an absolute/world frame, transform to be relative to this frame
	Frame3D ToAxisCoords(const Frame3D &f, bool bPreserveOrigin = false) const {
		Vector3 x = get_basis().xform(f.X());
		Vector3 y = get_basis().xform(f.Y());
		Vector3 z = get_basis().xform(f.Z());
		Frame3D w;
		w.origin = (bPreserveOrigin) ? f.origin : ToAxisCoords(f.origin);
		w.SetFrame(x, y, z);
		return w;
	}

	//! translate the reference frame origin
	void Translate(const Vector3 &vTranslate, bool bRelative) {
		if (bRelative) {
			origin += vTranslate;
		} else {
			origin = vTranslate;
		}
	}
	//! rotate the reference frame
	void Rotate(const Basis &mRotation, bool bReNormalize) {
		rotate(Quat(mRotation));
		if (bReNormalize) {
			normalize()
		}
	}

	//! make axes perpendicular. can "preserve" one axis by setting nPreserveAxis (0=X,1=Y,2=Z)

	void ReNormalize(int nPreserveAxis) {
		Vector3 rowX(get_basis().get_axis(0));
		Vector3 rowY(get_basis().get_axis(1));
		Vector3 rowZ(get_basis().get_axis(2));

		switch (nPreserveAxis) {
			default:
			case -1: {
				Vector3 vCrossY = rowZ.cross(rowX);
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
				Vector3 vCrossX = rowZ.cross(rowY);
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3 vCrossZ = rowX.cross(rowY);
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 0: {
				rowX.normalize();
				Vector3 vCrossY(rowX.cross(rowZ));
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
				Vector3 vCrossZ(rowX.cross(rowY));
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 1: {
				rowY.normalize();
				Vector3 vCrossX(rowY.cross(rowZ));
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3 vCrossZ(rowX.cross(rowY));
				vCrossZ.normalize();
				rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
			} break;

			case 2: {
				rowZ.normalize();
				Vector3 vCrossX(rowY.cross(rowZ));
				vCrossX.normalize();
				rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
				Vector3 vCrossY(rowX.cross(rowZ));
				vCrossY.normalize();
				rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
			} break;
		}

		get_basis().get_axis(0) = rowX;
		get_basis().get_axis(1) = rowY;
		get_basis().get_axis(2) = rowZ;
	}

	//! intersect ray with plane defined by one of the axis vectors
	bool IntersectRay(const Vector3 &vRayOrigin,
			const Vector3 &vRayDirection,
			Vector3 &vRayHit, int nPlaneNormalAxis) {
		Vector3 N = Axis(nPlaneNormalAxis);
		real_t d = -(origin.dot(N));
		real_t fDenom = vRayDirection.dot(N);
		if (Math::is_zero_approx(fDenom)) {
			return false;
		}

		real_t t = -(vRayOrigin.dot(N) + d) / fDenom;
		vRayHit = vRayOrigin + t * vRayDirection;
		return true;
	}
	Vector3 IntersectRay(const Vector3 &vRayOrigin,
			const Vector3 &vRayDirection,
			int nPlaneNormalAxis) {
		Vector3 vHit;
		bool bOK = IntersectRay(vRayOrigin, vRayDirection, vHit, nPlaneNormalAxis);
		// [TODO] return invalid vector here instead of ZERO
		return (bOK) ? vHit : Vector3();
	}

	//! create frame at an origin with a given x/y/z vectors
	Frame3D(const Vector3 &vOrigin, const Vector3 &vXAxis,
			const Vector3 &vYAxis, const Vector3 &vZAxis) :
			origin(vOrigin) {
		SetFrame(vXAxis, vYAxis, vZAxis);
	}

	//! create frame with a given rotation
	Frame3D(const Basis &mRotation) {
		set_basis(mRotation.transposed());
	}

	//! create frame at an origin with a given rotation
	Frame3D(const Vector3 &vOrigin,
			const Basis &mRotation) :
	{
		set_origin(vOrigin);
		set_basis(mRotation);
	}

	//! allow external setting of matrix...
	void SetFrame(const Basis &matRotate) {
		basis = matRotate.transposed();
	}

	//! set the reference frame axes
	void SetFrame(const Vector3 &vAxisX,
			const Vector3 &vAxisY,
			const Vector3 &vAxisZ) {
		basis.set_axis(0, vAxisX.normalized());
		basis.set_axis(1, vAxisY.normalized());
		basis.set_axis(2, vAxisZ.normalized());
	}

	void ComputePerpVectors(const Vector3 &vIn,
			Vector3 &vOut1, Vector3 &vOut2, bool bInIsNormalized) {
		Vector3 vPerp(vIn);
		if (!bInIsNormalized)
			vPerp.normalize();

		if (fabs(vPerp.x) >= fabs(vPerp.y) && fabs(vPerp.x) >= fabs(vPerp.z)) {
			vOut1.x = -vPerp.y;
			vOut1.y = vPerp.x;
			vOut1.z = 0.0;
		} else {
			vOut1.x = 0.0;
			vOut1.y = vPerp.z;
			vOut1.z = -vPerp.y;
		}

		vOut1.normalize();
		vOut2 = vPerp.cross(vOut1);
	}

	void ComputeAlignAxisMatrix(const Vector3 &vInitial,
			const Vector3 &vAlignWith, Basis &matrix) {
		// compute cosine of angle between vectors
		real_t axisDot = vAlignWith.dot(vInitial);

		// compute rotation axis
		Vector3 axisCross(vInitial.cross(vAlignWith));

		// apply rotation if necessary
		if (!Math::is_zero_approx(origin.distance_squared_to(axisCross))) {
			// compute normalized axis and angle, then create rotation around axis
			axisCross.normalize();
			real_t fAngle = Math::acos(axisDot / vAlignWith.norm());
			matrix = Basis(axisCross, fAngle);

		} else if (axisDot < (real_t)0) {
			// find some perpendicular vectors
			Vector3 vPerp1, vPerp2;
			ComputePerpVectors(vInitial, vPerp1, vPerp2, false);

			matrix = Basis(vPerp1, Math::deg2rad(180.0f));
		} else {
			matrix = Basis();
		}
	}

	//! rotate selected axis of this frame into toAxis
	void AlignAxis(int nAxis, const Vector3 &toAxis,
			bool bNormalized) {
		Basis matAlign;
		ComputeAlignAxisMatrix(
				Axis(nAxis), (bNormalized) ? toAxis : toAxis.normalized(), matAlign);
		Rotate(matAlign, false);
	}

	//! compute matrix that rotates this frame into toFrame
	void ComputeAlignmentMatrix(const Frame3D &toFrame,
			Basis &matRotate) {
		// align vCurFrame.Z() with vDestFrame.Z()
		Basis matAlignZ;
		ComputeAlignAxisMatrix(this->Z(), toFrame.Z(), matAlignZ);
		Frame3D vCopy(*this);
		vCopy.Rotate(matAlignZ, false);

		// compute rotation angle around vDestFrame.Z()
		Vector3 vX1(toFrame.X());
		Vector3 vX2(vCopy.X());
		Basis matAlignX;
		ComputeAlignAxisMatrix(vX2, vX1, matAlignX);

		matRotate = matAlignX * matAlignZ;

		vCopy = Frame3D(*this);
		vCopy.Rotate(matRotate, false);
		real_t fDotX = vCopy.X().dot(toFrame.X());
		real_t fDotY = vCopy.Y().dot(toFrame.Y());
		real_t fDotZ = vCopy.Z().dot(toFrame.Z());

		// check if Y is flipped - if it is, flip Y...
		// if ( vCopy.Y().dot( toFrame.Y() ) < 0 ) {
		//	Basis matFlip( 1,0,0, 0,-1,0, 0,0,1 );
		//	matRotate = matRotate * matFlip;
		//}
		//// check if Z is flipped - if it is, flip Z...
		if (vCopy.Z().dot(toFrame.Z()) < 0) {
			Basis matFlip;
			matFlip.get_axis(0) = Vector3(1, 0, 0);
			matFlip.get_axis(1) = Vector3(0, 1, 0);
			matFlip.get_axis(2) = Vector3(0, 0, -1);
			matRotate = matRotate * matFlip;
		}

		// [RMS] does this do anything??
		vCopy = Frame3D(*this);
		vCopy.Rotate(matRotate, false);
		fDotX = vCopy.X().dot(toFrame.X());
		fDotY = vCopy.Y().dot(toFrame.Y());
		fDotZ = vCopy.Z().dot(toFrame.Z());
	}
};
} // namespace g3