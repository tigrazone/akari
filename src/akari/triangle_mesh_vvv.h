#ifndef _TRIANGLE_MESH_VVV_
#define _TRIANGLE_MESH_VVV_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "qbvh.h"
#include "vec.h"

#include <time.h>

struct VVVVertex {
	Vec position_;
	Vec color_;
	Vec normal_;
};

struct VVVFace {
	unsigned int v0, v1, v2;
};

class VVVTriangle {
private:
public:
	VVVVertex *v_[3];
	Vec *p_[3]; // 位置

	VVVTriangle(VVVVertex *v0, VVVVertex *v1, VVVVertex *v2) {
		v_[0] = v0;
		v_[1] = v1;
		v_[2] = v2;
		
		p_[0] = &v0->position_;
		p_[1] = &v1->position_;
		p_[2] = &v2->position_;
	}

	BBox ObjectBound() const {
		return Union(BBox(*p_[0], *p_[1]), *p_[2]);
	}
};


class VVVTriangleMesh {	
private:
	std::vector<VVVVertex> vertex;
	std::vector<VVVFace> face;

	BBox bbox;

	QBVHSSE<VVVTriangle> bvh;

public:
	std::vector<VVVTriangle*> triangles;

	VVVTriangleMesh(std::string filename, const float scale, const Vec &origin) {
		
		std::ifstream ifs(filename, std::ios::in | std::ios::binary);
		if (!ifs) {
			std::cout << "Triangle Mesh Error: " << filename << std::endl;
			return;
		}
		
		/*
		vvv format:
		2	0x76767676	signature
		2				vsize
		vsize*9 float	position, color, nornal
		2				fsize
		fsize*3 int		3 indexes
		*/
		
		/*
		bool exportOBJ = true;
		FILE *objf;
		*/
		
		unsigned int header = 0;

		ifs.read((char *)&header, sizeof(unsigned int));
		if (header != 0x76767676)
			return;

		unsigned int vsize = 0;
		ifs.read((char *)&vsize, sizeof(unsigned int));

		/*
		if (exportOBJ)
		{
			size_t position = filename.find(".");
			std::string extractName = (std::string::npos == position) ? filename : filename.substr(0, position);
			//std::cout << "extractName=" << extractName.c_str() << std::endl;

			std::string fn_obj = extractName + ".obj";

			objf = fopen(fn_obj.c_str(), "wt");

			fprintf(objf, "# exported from akari vvv loader\n");
		}
		*/
		

		for (unsigned int i = 0; i < vsize; ++i) {
			float f[9];
			ifs.read((char *)f, sizeof(float) * 9);

			VVVVertex v;

			v.position_.x_ = f[0];
			v.position_.y_ = f[1];
			v.position_.z_ = f[2];
//			put(v.position_);
			v.position_ = origin + scale * v.position_;

			v.color_.x_ = f[3];
			v.color_.y_ = f[4];
			v.color_.z_ = f[5];

			v.normal_.x_ = f[6];
			v.normal_.y_ = f[7];
			v.normal_.z_ = f[8];

			/*
			if (exportOBJ)
			{
				fprintf(objf, "v %f %f %f\n", v.position_.x_, v.position_.y_, v.position_.z_);
				fprintf(objf, "vn %f %f %f\n", v.normal_.x_, v.normal_.y_, v.normal_.z_);
			}
			*/

			vertex.push_back(v);
			
			//std::cout << f[0] << " " << f[1] << " " << f[2] << std::endl;
		}
		
		unsigned int fsize = 0;
		ifs.read((char *)&fsize, sizeof(unsigned int));

		for (unsigned int i = 0; i < fsize; ++i) {
			unsigned int v[3];
			ifs.read((char *)&v, sizeof(unsigned int) * 3);
			VVVFace f;
			f.v0 = v[0];
			f.v1 = v[1];
			f.v2 = v[2];

			/*
			if (exportOBJ)
			{
				fprintf(objf, "f v%d//vn%d v%d//vn%d v%d//vn%d\n", f.v0+1, f.v0+1, f.v1+1, f.v1+1, f.v2+1, f.v2+1);
			}
			*/

			face.push_back(f);
		}

		/*
		
		if (exportOBJ)
		{
			fclose(objf);
		}
		*/
		

		std::cout << "=== QBVH ===" << std::endl;
		std::cout << "vertex: " << vertex.size() << std::endl;
		std::cout << "face: " << face.size() << std::endl;
		
		for (int i = 0; i < face.size(); i ++) {
			triangles.push_back(new VVVTriangle(&vertex[face[i].v0], &vertex[face[i].v1], &vertex[face[i].v2]));
		}

		for (int i = 0; i < vertex.size(); i ++) {
			bbox = Union(bbox, vertex[i].position_);
		}

		put(bbox.pmin_); std::cout << std::endl;
		put(bbox.pmax_); std::cout << std::endl;
		
		
		std::cout << "create BVH" << std::endl;

		clock_t t0 = clock();

		bvh.CreateBVHFromTriangle2s(triangles);

		clock_t t1 = clock();
		
		std::cout << "done BVH ";

		printf("%.2fs\n\n", (float)(t1 - t0) / CLOCKS_PER_SEC);
	}

	virtual ~VVVTriangleMesh() {
		for (int i = 0; i < triangles.size(); i ++)
			delete triangles[i];
	}

	BBox ObjectBound() const {
		return bbox;
	}

	bool intersect(const Ray &ray, const Hitpoint &old_hitpoint, Hitpoint* hitpoint, const bool using_interpolated_normal) {
		bool check = bbox.check_intersect(ray, NULL, NULL);
		if (check) {
			const bool ret = bvh.intersect(ray, old_hitpoint, hitpoint, using_interpolated_normal);
			return ret;
		}

		return false;
	}
};

#endif