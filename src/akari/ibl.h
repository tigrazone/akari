#ifndef _IBL_H_
#define _IBL_H_

#include "hdr.h"
#include "constant.h"
#include "ray.h"
#include "random.h"

#include <algorithm>
#include <array>

class IBL {
private:
	HDRImage hdr_;

	std::vector<float> luminance_;
	
	//not used
	std::vector<float> normal_;

	//not used
	//std::vector<Color> importance_map_;

	std::vector<Vec> direction_table_;

	static const int importance_map_size_ = 64;

public:
	int width_, height_;
	float width_1;

	inline float luminance(const Color &color) {
		return 0.298912f * color.x_ + 0.586611f * color.y_ + 0.114478f * color.z_;
	}
	

	IBL(const std::string &filename, const bool importance = false) {
		hdr_.load_unsafe(filename);
		width_ = hdr_.width();
		height_ = hdr_.height();
		
		width_1 = width_ / (kPI2);
		
//		importance_map_size_ = importance_map_size;
		
		if (importance) {
			normal_.resize(width_ * height_);

			// インポータンスマップ作成
			// cos項考慮した環境マップにワープする処理
			std::vector<Color> warped(width_ * height_);
			for (int iy = 0; iy < height_; ++iy) {
				for (int ix = 0; ix < width_; ++ix) {
					const float u = (ix + 0.5f) / width_;
					const float v = 1.0f - (iy + 0.5f) / height_;
					
					const float phi = u * kPI2;
					const float y = v +v - 1.0f;
					
					const float yy = sqrt(1.0f - y*y);

					const Vec vec(yy * cos(phi), y, yy * sin(phi));

					warped[iy * width_ + ix] = sample(vec);
				}
			}

			// 縮小
			const int scale = width_ / importance_map_size_;
			luminance_.resize(importance_map_size_ * importance_map_size_);
			direction_table_.resize(importance_map_size_ * importance_map_size_);
			
//			HDRImage tmphdr(importance_map_size_, importance_map_size_);

			for (int iy = 0; iy < importance_map_size_; ++iy) {
				for (int ix = 0; ix < importance_map_size_; ++ix) {
					const float u = (ix + 0.5f) / importance_map_size_;
					const float v = (iy + 0.5f) / importance_map_size_;
					
					const float phi = u * 2.0f * kPI;
					const float y   = v * 2.0f - 1.0f;
					const Vec vec(sqrt(1.0f - y*y) * cos(phi), y, sqrt(1.0f - y*y) * sin(phi));

					direction_table_[iy * importance_map_size_ + ix] = vec;

					Color accum;
					for (int iiy = iy * scale; iiy < (iy + 1) * scale; ++iiy) {
						for (int iix = ix * scale; iix < (ix + 1) * scale; ++iix) {
							accum = accum + warped[iiy * width_ + iix];
						}
					}
//					*tmphdr.image_ptr(ix, iy) = accum;
					luminance_[(importance_map_size_ - iy - 1) * importance_map_size_ + ix] = luminance(accum);
				}
			}
		}
	}
	
	// DebevecのオリジナルHDR
	inline Color sample_sphere(const Vec &dir) {
		const float r = (kPI_1) * acos(dir.z_) / sqrt(dir.x_ * dir.x_ + dir.y_ * dir.y_);

		float u = (dir.x_ * r + 1.0f) * 0.5f;
		float v = 1.0f - (dir.y_ * r + 1.0f) * 0.5f;
		
		if (u < 0.0f)
			u += 1.0f;
		if (v < 0.0f)
			v += 1.0f;

		const int x = (int)(u * width_) % width_;
		const int y = height_ - 1 - (int)(v * height_) % height_;

		return *hdr_.image_ptr(x, y);
	}

	// phi = [0, 2pi], theta = [0, pi]
	// theta = 0 が真上
	// v = 1が真上
	// u = [0, 1), v = [0, 1)が直接phi,thetaにマッピングされてるサンプリング方法
	inline Color sample(const Vec &dir) {
		const float theta = acos(dir.y_);
		float phi = acos(dir.x_ / sqrt(dir.x_ * dir.x_ + dir.z_ * dir.z_));
		
		if (dir.z_ < 0.0f)
			phi = kPI + kPI - phi;

		//const float u = phi / (2.0f * kPI);
		const float v = 1.0f - theta * kPI_1;

		//const int x = (int)(u * width_) % width_;
		const int x = (int)(phi * width_1) % width_;
		const int y = (int)(v * height_) % height_;

		return *hdr_.image_ptr(x, y);
	}
	struct Sample {
		Vec dir_;
		float weight_;
		float weight_1;
	};
	

	void importance_sampling_unsafe_fast(const int usamples, const int vsamples, const Vec &normal, std::vector<Sample> &vecs, Random *rnd) {
		struct Bin {
			int index_;
			float probability_;
		};

		const long ims2 = importance_map_size_ * importance_map_size_;

		const float ims_1 = 1.0 / importance_map_size_;

		float is1;
		float ttl1;

		Bin bin[ims2];
		// std::vector<Bin> bin(importance_map_size_ * importance_map_size_);

		int num_bin = 0;
		// bin.reserve(importance_map_size_ * importance_map_size_);
		float total = 0.0f;
		for (int i = 0; i < ims2; ++i) {
			const float lumi = luminance_[i] * std::max(dot(normal, direction_table_[i]), 0.0f);
			if (lumi > 0.0f) {
				bin[num_bin].index_ = i;
				bin[num_bin].probability_ = lumi;
				num_bin ++;
				total += lumi;
			}
		}
		

		ttl1 = 1.0f / total;

		for (int i = 0; i < num_bin; ++i) {
			bin[i].probability_ = bin[i].probability_ * ttl1;
			const float samples = usamples * vsamples * bin[i].probability_;

			int int_samples = (int)samples;
			const float fraction = samples - int_samples;
			if (rnd->next01() < fraction)
				int_samples ++;
			
			const int ix = bin[i].index_ % importance_map_size_;
			const int iy = bin[i].index_ / importance_map_size_;

			is1 = 1.0f / int_samples;

			// QMC
			// hammersley
			float p, qu, qv;
			int k, kk, pos;
			for (k=0, pos=0; k < int_samples ; ++k) {
				qu = 0;
				for (p=0.5f, kk = k; kk; p *= 0.5f, kk>>=1)
					if (kk & 1) // kk mod 2 == 1
						qu += p;

				qv = (k + 0.5f) ;
				
				const float u = (ix + qu) * ims_1;
				const float v = (iy + qv) * ims_1;
				
				const float phi = u * kPI2;
				const float t = v + v - 1.0f;

				const float sqqq = sqrt(1.0f - t * t);

					//const float theta = sqqq * kPI;
			
				const Vec vec(sqqq * cos(phi), t, sqqq * sin(phi));

				Sample s;
				s.dir_ = vec;
				
				/*
				s.weight_ = bin[i].probability_ * (importance_map_size_ * importance_map_size_ / (4.0f * kPI));
				s.weight_1 = 1.0 / (s.weight_ * kPI);
				*/
				
				s.weight_ = bin[i].probability_ * (ims2);
				s.weight_1 = 4.0f / s.weight_;
				
				vecs.push_back(s);

				pos ++;
			}
		}
	}


};

#endif