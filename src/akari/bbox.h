#ifndef _BBOX_H_
#define _BBOX_H_

#include <algorithm>

#include "vec.h"
#include "constant.h"
#include "ray.h"

class BBox {
public:
	Vec pmin_, pmax_;

	BBox() {
		pmin_ = Vec(kINF, kINF, kINF);
		pmax_ = -1.0 * pmin_;
	}

	BBox(const Vec &p) : pmin_(p), pmax_(p) {}

	BBox(const Vec &p1, const Vec &p2) {
		pmin_ = Vec(std::min(p1.x_, p2.x_), std::min(p1.y_, p2.y_), std::min(p1.z_, p2.z_));
		pmax_ = Vec(std::max(p1.x_, p2.x_), std::max(p1.y_, p2.y_), std::max(p1.z_, p2.z_));
	}

	bool inside(const Vec &pt) const {
		return (pmin_.x_ <= pt.x_ && pt.x_ <= pmax_.x_ &&
				pmin_.y_ <= pt.y_ && pt.y_ <= pmax_.y_ &&
				pmin_.z_ <= pt.z_ && pt.z_ <= pmax_.z_);
	}

	Vec &operator[](int i) {
		if (i == 0)
			return pmin_;
		return pmax_;
	}	
	const Vec &operator[](int i) const {
		if (i == 0)
			return pmin_;
		return pmax_;
	}

	void expand(float delta) {
		const Vec v(delta, delta, delta);
		pmin_ = pmin_ - v;
		pmax_ = pmax_ + v;
	}

	float surface_area() {
		const Vec d = pmax_ - pmin_;
		return 2.0f * (d.x_ * d.y_ + d.x_ * d.z_ + d.y_ * d.z_);
	}
	float volume() {
		const Vec d = pmax_ - pmin_;
		return d.x_ * d.y_ * d.z_;
	}

	enum LongestAxis {
		AxisX,
		AxisY,
		AxisZ,
	};

	LongestAxis maximum_extent() const {
		const Vec diag = pmax_ - pmin_;
		if (diag.x_ > diag.y_ && diag.x_ > diag.z_)
			return AxisX;
		else if (diag.y_ > diag.z_)
			return AxisY;
		else
			return AxisZ;
	}
	
//tigra: need to speedup
	bool check_intersect(const Ray &ray, float *hitt0, float *hitt1) {
		float t0 = 0.0, t1 = kINF;
		for (int i = 0; i < 3; ++i) {
			// Update interval for _i_th bounding box slab
			float invRayDir = 1.f / ray.dir_[i];
			float tNear = (pmin_[i] - ray.org_[i]) * invRayDir;
			float tFar  = (pmax_[i] - ray.org_[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$s
			if (tNear > tFar) std::swap(tNear, tFar);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar  < t1 ? tFar  : t1;
			if (t0 > t1) return false;
		}
		if (hitt0) *hitt0 = t0;
		if (hitt1) *hitt1 = t1;
		return true;
	}
};

/*
//pp можно vec* чтоб не копировать значение min, max
//vec *pp[2];
//pp[0]=&min;
//pp[1]=&max;

//pp[0]->x

class aabb
{
    public:
        aabb() {}
        aabb(const vec3& a, const vec3& b) { pp[0] = a; pp[1] = b; }
		
        vec3 min() const { return pp[0]; }
        vec3 max() const { return pp[1]; }
        vec3 pp[2];
};


inline bool aabb::hit14(const ray& r, float tmin, float tmax) const {	
int posneg[6];	
		
		if(r.direction()[0]<0.f)
		{
			if(r.origin()[0]<min()[0])
				return false;
			posneg[0]= 1;
			posneg[1]= 0;
		}
		else
		{
			if(r.origin()[0]>max()[0])
				return false;
			posneg[0]= 0;
			posneg[1]= 1;
		}			
		
	const float rd0 = 1.f / r.direction()[0];
		
    float t0 = (pp[posneg[0]].x() - r.o.x() ) * rd0 ;
    float t1 = (pp[posneg[1]].x() - r.o.x() ) * rd0 ;	

    float interval_min = tmin;
    float interval_max = tmax;
	
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
	
	
		
		if(r.direction()[1]<0.f)
		{
			if(r.origin()[1]<min()[1])
				return false;
			posneg[2]= 1;
			posneg[3]= 0;
		}
		else
		{
			if(r.origin()[1]>max()[1])
				return false;
			posneg[2]= 0;
			posneg[3]= 1;
		}
	
	
	const float rd1 = 1.f / r.direction()[1];
		
    t0 = (pp[posneg[2]].y() - r.o.y() ) * rd1 ;
    t1 = (pp[posneg[3]].y() - r.o.y() ) * rd1 ;	
	
    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    if (interval_min > interval_max) return false;
	
		
		if(r.direction()[2]<0.f)
		{
			if(r.origin()[2]<min()[2])
				return false;
			posneg[4]= 1;
			posneg[5]= 0;
		}
		else
		{
			if(r.origin()[2]>max()[2])
				return false;
			posneg[4]= 0;
			posneg[5]= 1;
		}
	
	const float rd2 = 1.f / r.direction()[2];
		
    t0 = (pp[posneg[4]].z() - r.o.z() ) * rd2 ;
    t1 = (pp[posneg[5]].z() - r.o.z() ) * rd2 ;

    if (t0 > interval_min) interval_min = t0;
    if (t1 < interval_max) interval_max = t1;
    return (interval_min <= interval_max);
}
*/


inline BBox Union(const BBox &b, const Vec &p) {
	BBox ret = b;
	ret.pmin_.x_ = std::min(b.pmin_.x_, p.x_);
	ret.pmin_.y_ = std::min(b.pmin_.y_, p.y_);
	ret.pmin_.z_ = std::min(b.pmin_.z_, p.z_);
		
	ret.pmax_.x_ = std::max(b.pmax_.x_, p.x_);
	ret.pmax_.y_ = std::max(b.pmax_.y_, p.y_);
	ret.pmax_.z_ = std::max(b.pmax_.z_, p.z_);

	return ret;
}

inline BBox Union(const BBox &b1, const BBox &b2) {
	BBox ret = b1;
	ret.pmin_.x_ = std::min(b1.pmin_.x_, b2.pmin_.x_);
	ret.pmin_.y_ = std::min(b1.pmin_.y_, b2.pmin_.y_);
	ret.pmin_.z_ = std::min(b1.pmin_.z_, b2.pmin_.z_);
	
	ret.pmax_.x_ = std::max(b1.pmax_.x_, b2.pmax_.x_);
	ret.pmax_.y_ = std::max(b1.pmax_.y_, b2.pmax_.y_);
	ret.pmax_.z_ = std::max(b1.pmax_.z_, b2.pmax_.z_);

	return ret;
}

#endif