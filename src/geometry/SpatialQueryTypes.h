#pragma once

#include <g3types.h>

namespace g3 {

class IPixelHitRadius {
public:
	virtual ~IPixelHitRadius() {}
	virtual float GetWorldHitRadius(const Vector3 &vHit) const = 0;
};

class HitTestRay {
public:
	Vector3 vOrigin;
	Vector3 vDirection;
	IPixelHitRadius *pHitThresh;

	HitTestRay() {
		vOrigin = vDirection = Vector3::Zero();
		pHitThresh = nullptr;
	}
	HitTestRay(const Vector3 &o, const Vector3 &d) {
		vOrigin = o;
		vDirection = d;
		pHitThresh = nullptr;
	}
	HitTestRay(const Vector3 &o, const Vector3 &d, IPixelHitRadius *pPixelRadius) {
		vOrigin = o;
		vDirection = d;
		pHitThresh = pPixelRadius;
	}
};

} // namespace g3
