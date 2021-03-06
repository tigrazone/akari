﻿#include <iostream>
#include <omp.h>

//#include "bmpexporter.h"

#include "scene.h"
#include "setting.h"
#include "vec.h"
#include "hdr.h"
#include "random.h"
#include "ibl.h"
#include "triangle_mesh_vvv.h"
#include "render.h"

#include <time.h>

#pragma comment(lib, "winmm.lib")


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "bmpexporter.h"

#pragma pack(push,2)

typedef struct tagBITMAPFILEHEADER {
  unsigned short bfType;
  unsigned long  bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned long  bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER{
	unsigned long  biSize;
	long           biWidth;
	long           biHeight;
	unsigned short biPlanes;
	unsigned short biBitCount;
	unsigned long  biCompression;
	unsigned long  biSizeImage;
	long           biXPixPerMeter;
	long           biYPixPerMeter;
	unsigned long  biClrUsed;
	unsigned long  biClrImporant;
} BITMAPINFOHEADER;

#pragma pack(pop)

int exportToBmp(
	const char* fileName, 
	unsigned char* pixel, 
	unsigned int width,
	unsigned int height )
{
	//
	int success = 0;
	unsigned int x = 0;
	unsigned int y = 0;
	FILE* file = 0;
	BITMAPFILEHEADER* header = 0;
	BITMAPINFOHEADER* infoHeader = 0;
	unsigned int scanLineLengthInBytes = 0;
	int totalFileLength = 0;
	unsigned char* bmpMemoryStart = 0;
	unsigned char* bmpMemoryCursor = 0;
	//
	// スキャンライン長計算
	scanLineLengthInBytes = width*3;
	// ファイル全体の長さを計算
	totalFileLength = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+scanLineLengthInBytes*height;
	// メモリを確保
	bmpMemoryStart = (unsigned char*)malloc(totalFileLength);
	
	if( !bmpMemoryStart )
	{ goto EXIT; }

	bmpMemoryCursor = bmpMemoryStart;
	// BITMAPFILEHEADERを作成
	header = (BITMAPFILEHEADER*)bmpMemoryCursor; //
	header->bfType = 'B' | ('M' << 8); // ファイルタイプ
	header->bfSize = totalFileLength;// ファイルサイズ (byte)
	header->bfReserved1 = 0; // 予約領域
	header->bfReserved2 = 0; // 予約領域
	header->bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER); // ファイル先頭から画像データまでのオフセット (byte)
	bmpMemoryCursor += sizeof(BITMAPFILEHEADER);
	// BITMAPINFOHEADERを作成
	infoHeader = (BITMAPINFOHEADER*)bmpMemoryCursor;
	infoHeader->biSize = sizeof(BITMAPINFOHEADER); // 情報ヘッダのサイズ
	infoHeader->biWidth = width; // 画像の幅 (ピクセル)
	infoHeader->biHeight = height; // 画像の高さ (ピクセル)
	infoHeader->biPlanes = 1; // プレーン数
	infoHeader->biBitCount = 24; // 1 画素あたりのデータサイズ (bit)
	infoHeader->biCompression = 0; // 圧縮形式(無圧縮)
	infoHeader->biSizeImage = width*height*3; // 画像データ部のサイズ (byte)
	infoHeader->biXPixPerMeter = 3780; // 横方向解像度 (1mあたりの画素数)
	infoHeader->biYPixPerMeter = 3780; // 縦方向解像度 (1mあたりの画素数)
	infoHeader->biClrUsed = 0; // 格納されているパレット数 (使用色数)	
	infoHeader->biClrImporant = 0; // 重要なパレットのインデックス
	bmpMemoryCursor += sizeof(BITMAPINFOHEADER);
	// 全てのデータを書き込む
	for(y=0;y<height;++y)
	{
		for(x=0;x<width;++x)
		{
			const unsigned int srcBase = (x+(height-y-1)*width)*3;
			const unsigned int dstBase = x*3;
			bmpMemoryCursor[dstBase+0] = pixel[srcBase+0];
			bmpMemoryCursor[dstBase+1] = pixel[srcBase+1];
			bmpMemoryCursor[dstBase+2] = pixel[srcBase+2];
		}
		bmpMemoryCursor += scanLineLengthInBytes;
	}
	// ファイルに書き込む
	
	#if defined(__GNUC__)
	file=fopen(fileName,"wb");
#elif defined(__MSVC__)
	fopen_s(&file,fileName,"wb");
#else
#error Unsupported compiler!
#endif

	
	if(!file)
	{ goto EXIT; }
	if(fwrite(bmpMemoryStart,totalFileLength,1,file)!=1)
	{ goto EXIT; }
	if(fclose(file)==EOF)
	{ goto EXIT; }
	// 最後まで到達したので成功
	success = 1;
EXIT:
	if(bmpMemoryStart)
	{
		free(bmpMemoryStart);
	}
	if(file)
	{
		fclose(file);
	}
	return success;
}


int width0, height0;

inline float luminance(const Color &color) {
	return 0.298912f * color.x_ + 0.586611f * color.y_ + 0.114478f * color.z_;
}
class LessColor {
public:
	bool operator()(const Color* riLeft, const Color* riRight) const {
		return luminance(*riLeft) < luminance(*riRight);
	}
};

unsigned char to_bmp_value(const float value, const float display_gamma) {
	const unsigned int v = (unsigned int)(pow(value, 1.0f / display_gamma) * 255.0f + 0.5f);
	if (v >= 256)
		return (unsigned char)255;
	return (unsigned char)v;
}

void to_ldr(HDRImage *hdr_image, const float fstop) {
	const float scale = pow(2.0f, fstop);
	for (int iy = 0; iy < hdr_image->height(); ++iy) {
		for (int ix = 0; ix < hdr_image->width(); ++ix) {
			Color col = *hdr_image->image_ptr(ix, iy) * scale;
			if (col.x_ >= 1.0)
				col.x_ = 1.0;
			if (col.y_ >= 1.0)
				col.y_ = 1.0;
			if (col.z_ >= 1.0)
				col.z_ = 1.0;
			*hdr_image->image_ptr(ix, iy) = col;
		}
	}
}

//tigra
void flo2hdr(const char *fn,const std::vector<float> mColor)
{	
	std::ofstream hdr(fn, std::ios::binary);

	const static int mResY = height0;
	const static int mResX = width0;
		
        hdr << "#?RADIANCE" << '\n';
        hdr << "# SmallVCM" << '\n';
        hdr << "FORMAT=32-bit_rle_rgbe" << '\n' << '\n';
        hdr << "-Y " << mResY << " +X " << mResX << '\n';

		for (int y = mResY-1; y >=0; y--)
        {
            for(int x=0; x<mResX; x++)
            {
                typedef unsigned char byte;
                byte rgbe[4] = {0,0,0,0};

                const Vec rgbF = Vec(mColor[x + y*mResX],mColor[x + y*mResX],mColor[x + y*mResX]);
                float v = std::max(rgbF.x_, std::max(rgbF.y_, rgbF.z_));

                //if(v >= 1e-32f)
                {
                    int e;
                    //v = float(frexp(v, &e) * 256.f *rcp(v));
                    v = float(frexp(v, &e) * 256.f / v);
                    rgbe[0] = byte(rgbF.x_ * v);
                    rgbe[1] = byte(rgbF.y_ * v);
                    rgbe[2] = byte(rgbF.z_ * v);
                    rgbe[3] = byte(e + 128);
                }

                hdr.write((char*)&rgbe[0], 4);
            }
        }
	
}

void int2hdr(const char *fn, const std::vector<int> mColor)
{
	std::ofstream hdr(fn, std::ios::binary);

	const static int mResY = height0;
	const static int mResX = width0;

	hdr << "#?RADIANCE" << '\n';
	hdr << "# SmallVCM" << '\n';
	hdr << "FORMAT=32-bit_rle_rgbe" << '\n' << '\n';
	hdr << "-Y " << mResY << " +X " << mResX << '\n';

	for (int y = mResY - 1; y >= 0; y--)
	{
		for (int x = 0; x<mResX; x++)
		{
			typedef unsigned char byte;
			byte rgbe[4] = { 0, 0, 0, 0 };

			const Vec rgbF = Vec(mColor[x + y*mResX], mColor[x + y*mResX], mColor[x + y*mResX]);
			float v = std::max(rgbF.x_, std::max(rgbF.y_, rgbF.z_));

			//if (v >= 1e-32f)
			{
				int e;
				//v = float(frexp(v, &e) * 256.f *rcp(v));
				v = float(frexp(v, &e) * 256.f / v);
				rgbe[0] = byte(rgbF.x_ * v);
				rgbe[1] = byte(rgbF.y_ * v);
				rgbe[2] = byte(rgbF.z_ * v);
				rgbe[3] = byte(e + 128);
			}

			hdr.write((char*)&rgbe[0], 4);
		}
	}

}


void tone_curve(HDRImage *hdr_image, const Color input, const Color output) {
	const Color z = Color(log(output.x_) / log(input.x_), log(output.y_) / log(input.y_), log(output.z_) / log(input.z_));

	for (int iy = 0; iy < hdr_image->height(); ++iy) {
		for (int ix = 0; ix < hdr_image->width(); ++ix) {
			Color col = *hdr_image->image_ptr(ix, iy);
			*hdr_image->image_ptr(ix, iy) = Color(pow(col.x_, z.x_), pow(col.y_, z.y_), pow(col.z_, z.z_));
		}
	}
}

void vignetting(HDRImage *hdr_image, const float radius) {
	for (int iy = 0; iy < hdr_image->height(); ++iy) {
		for (int ix = 0; ix < hdr_image->width(); ++ix) {

			int cx = ix - hdr_image->width() / 2;
			int cy = iy - hdr_image->height() / 2;
			float length = sqrt((float)(cx * cx + cy * cy));

			float vig = 1.0f;

			if (length - radius >= 0.0)
				vig = 1.0f - pow((length - radius) / (0.5f * radius), 1.0f);
			if (vig < 0.0f)
				vig = 0.0f;

			*hdr_image->image_ptr(ix, iy) = *hdr_image->image_ptr(ix, iy) * vig;
		}
	}
}

template <typename T>
void blur(std::vector<T> &arr, const std::vector<float> &depth_arr, const int width, const int height, const float radius, const float depth_threashold) {
	std::vector<T> tmp = arr;

	// カーネル分解した方が良いと思われる
	#pragma omp parallel for schedule(dynamic, 1)
		for (int y = 0; y < height; y ++) {
			for (int x = 0; x < width; x ++) {
				T accum = T(0.0);
				float num = 0;

				const float depth = depth_arr[y * width + x];

				for (int ix = -radius; ix <= radius; ++ix) {
					for (int iy = -radius; iy <= radius; ++iy) {
						float weight = sqrt((float)(ix * ix + iy * iy));
						if (weight > radius)
							continue;

						weight = 1.0 / (weight + 1.0);

						const int nx = x + ix;
						const int ny = y + iy;
						if (nx < 0 || width <= nx)
							continue;
						if (ny < 0 || height <= ny)
							continue;

						if (fabs(depth - depth_arr[ny * width + nx]) >= depth_threashold)
							continue;

						num += weight;
						accum = accum + weight * tmp[ny * width + nx];
					}
				}

				arr[y * width + x] = accum / num;
			}
		}
}

int main() {


	clock_t beginTime000 = clock();

	clock_t endTime1;
	float now1;

	Setting setting("setting.txt");
	Render render(&setting);

	// スレッド数設定
	omp_set_num_threads(setting.int_value("num_threads", 1));

	// 絞りファイル読み込み
	HDRImage apertureImage;
	apertureImage.load_unsafe(setting.string_value("aperture"));
	std::vector<float> lens_cdf(apertureImage.width() * apertureImage.height());
	float total = 0.0f;
	for (int iy = 0; iy < apertureImage.height(); ++iy) {
		for (int ix = 0; ix < apertureImage.width(); ++ix) {
			total += apertureImage.image_ptr(ix, iy)->x_;
			lens_cdf[iy * apertureImage.width() + ix] = total;
		}
	}
	for (int i = 0; i < lens_cdf.size(); ++i) {
		lens_cdf[i] /= total;
	}

	// 各種パラメータ設定
	const int width  = setting.int_value("width"); 
	const int height = setting.int_value("height"); 

	const int direct_light_samples   = setting.int_value("direct_light_samples");
	const int indirect_light_samples = setting.int_value("indirect_light_samples");

	const int iteration_num = setting.int_value("iteration_num");
	const int sample_per_pixel_per_iteration = setting.int_value("sample_per_pixel_per_iteration");
	int total_samples = 0;


	const Vec camera_position = Vec(setting.float_value("camera_position_x"), setting.float_value("camera_position_y"), setting.float_value("camera_position_z"));
	const Vec camera_dir      = normalize(Vec(setting.float_value("camera_dir_x"), setting.float_value("camera_dir_y"), setting.float_value("camera_dir_z")));
	const Vec camera_up       = normalize(Vec(setting.float_value("camera_up_x"), setting.float_value("camera_up_y"), setting.float_value("camera_up_z")));


	// ワールド座標系でのスクリーンの大きさ
	const float screen_width = setting.float_value("screen_width");
	const float screen_height= setting.float_value("screen_height");
	// スクリーンまでの距離
	const float screen_dist  = setting.float_value("screen_dist");
	// スクリーンを張るベクトル
	const Vec screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
	const Vec screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;

	const Vec screen_center = camera_position + camera_dir * screen_dist;

	const float lens_param_a = setting.float_value("lens_param_a");
	const float lens_param_b = setting.float_value("lens_param_b");
	const float lens_radius   = setting.float_value("lens_radius");
	const float focal_distance = setting.float_value("focal_distance");
	const float blur_depth_threashold = setting.float_value("blur_depth_threashold");
	const float blur_total_samples_coeff = setting.float_value("blur_total_samples_coeff");
	const int blur_max_radius = setting.int_value("blur_max_radius");
	const float fstop = setting.float_value("fstop");
	const float vignetting_value = setting.float_value("vignetting_value");

	const int tone_curve_input_r = setting.int_value("tone_curve_input_r");
	const int tone_curve_input_g = setting.int_value("tone_curve_input_g");
	const int tone_curve_input_b = setting.int_value("tone_curve_input_b");
	
	const int tone_curve_output_r = setting.int_value("tone_curve_output_r");
	const int tone_curve_output_g = setting.int_value("tone_curve_output_g");
	const int tone_curve_output_b = setting.int_value("tone_curve_output_b");

	const float output_interval = setting.float_value("output_interval");
	const float display_gamma = setting.float_value("display_gamma");

	HDRImage hdr(width, height);
	HDRImage depth_hdr(width, height);
	HDRImage result_hdr(width, height);
	HDRImage tmp_hdr(width, height);
	
	width0	= width;
	height0	= height;

	int count = 0;
	
	clock_t beginTime = clock(); 


	float ttt = (float)(beginTime - beginTime000) / (float)CLOCKS_PER_SEC;
	
	printf("pre-render time %.2fs\n", ttt);

	clock_t beginTime0 = clock();

	//timeBeginPeriod(1);
	{
		
		std::vector <int> sample_map(width * height);
		std::vector <float> sample_mapf(width * height);
		
		for (int i = 0; i < sample_map.size(); ++i) {
			sample_map[i] = sample_per_pixel_per_iteration;
		}
		
		std::vector<float> total_samples(width * height);
		std::vector<float> luminance_map(width * height);
		std::vector<float> luminance2_map(width * height);
		std::vector<float> variance_map(width * height);
		std::vector<float> importance_map(width * height);
		std::vector<float> importance_map1(width * height);
		std::vector<float> depth_arr(width * height);

		std::vector<unsigned char> bmp_arr(width * height * 3);

		std::cout << "Image size " << width << "x" << height << std::endl;
		std::cout << "direct_light_samples=" << direct_light_samples << std::endl;
		std::cout << "indirect_light_samples=" << indirect_light_samples << std::endl;
		std::cout << "sample_per_pixel_per_iteration=" << sample_per_pixel_per_iteration << std::endl;

		Random grnd(0);
		
		for (int iteration = 0; iteration < iteration_num; ++iteration) {
			std::cout << std::endl << "iteration: " << iteration << std::endl;
			
			for (int y = 0; y < height; y ++) {
				if (y % 64 == 0)
					std::cerr << y << " ";

			#pragma omp parallel for schedule(dynamic, 1)
				for (int x = 0; x < width; x ++) {
					Random rnd((y * width + x) * iteration_num + iteration);
					Color accumulated_radiance = Color();
					float accumulated_depth = 0.0f;
					for (int ss = 0; ss < sample_map[y * width + x]; ss ++) {
						const float r1 = rnd.next01();
						const float r2 = rnd.next01();
						// スクリーン上の位置
						const float wx = (r1 + x) / width - 0.5f;
						const float wy = (r2 + y) / height - 0.5f;

						const float length = sqrt(wx * wx + wy * wy);

						const float wx2 = wx * pow(lens_param_a, length + lens_param_b);
						const float wy2 = wy * pow(lens_param_a, length + lens_param_b);

						const Vec screen_position = 
							screen_center + 
							screen_x * wx2 +
							screen_y * wy2;

						// レイを飛ばす方向
						const Vec dir = normalize(screen_position - camera_position);
						Ray ray(camera_position, dir);

						// レイがレンズにあたって屈折するシミュレーション
						float lensU, lensV;
						const float r = rnd.next01();
						
						std::vector<float>::iterator it = std::lower_bound(lens_cdf.begin(), lens_cdf.end(), r);
						const int idx = (int)(it - lens_cdf.begin());
						const int lens_ix = idx % apertureImage.width();
						const int lens_iy = idx / apertureImage.width();
						
						lensU = ((float)lens_ix / apertureImage.width()) * 2.0f - 1.0f;
						lensV = ((float)lens_iy / apertureImage.height()) * 2.0f - 1.0f;
						lensU *= lens_radius;
						lensV *= lens_radius;
						const float ft = fabs(focal_distance / ray.dir_.z_);
						const Vec Pfocus = ray.org_ + ray.dir_ * ft;
						ray.org_ = ray.org_ + Vec(lensU, lensV, 0.f);
						ray.dir_ = normalize(Pfocus - ray.org_);

						// 直接光の影響計算
						Intersection intersection;
						const Color direct_light = render.radiance_direct_light_env(ray, &rnd, direct_light_samples, direct_light_samples, Intersection(), &intersection, 0);

						if (valid(direct_light)) {
							accumulated_radiance = accumulated_radiance + direct_light;
							if (intersection.hitpoint_.distance_ <= -kINF)
								intersection.hitpoint_.distance_ = kINF;
							accumulated_depth += intersection.hitpoint_.distance_;
						}
						
						// 間接光の影響計算
						
						const float ils1 = 1.0 / (float)indirect_light_samples;
						
						for (int v = 0; v < indirect_light_samples; ++ v) {
							const Color indirect_light = render.radiance_indirect_light_env(ray, &rnd, Intersection(), NULL, 0) * ils1;
							if (valid(indirect_light))
								accumulated_radiance = accumulated_radiance + indirect_light;
						}

						const float lumi =  luminance(accumulated_radiance);
						luminance_map[y * width + x] += lumi;
						luminance2_map[y * width + x] += lumi * lumi;
						total_samples[y * width + x] ++;
					}
					(*hdr.image_ptr(x, y)) = (*hdr.image_ptr(x, y)) + accumulated_radiance;
					(*depth_hdr.image_ptr(x, y)) = (*depth_hdr.image_ptr(x, y)) + Color(accumulated_depth, accumulated_depth, accumulated_depth);
				}

				// 経過時間チェック
				clock_t endTime = clock();
				const float now = (float)(endTime - beginTime) / (float)CLOCKS_PER_SEC;

				if (now >= output_interval) 
				{
					std::cout << std::endl;
					beginTime = clock();

					// 出力
					for (int y = 0; y < height; y ++) {
						for (int x = 0; x < width; x ++) {
							// 画像計算
							(*tmp_hdr.image_ptr(x, y)) = (*hdr.image_ptr(x, y)) / total_samples[y * width + x];
						}
					}

					
					// 最終結果をぼかす
					//std::cout << "blur rendered image" << std::endl;
					for (int y = 0; y < height; y ++) {
				// #pragma omp parallel for schedule(dynamic)
				
				// #pragma omp parallel for shared(depth_arr,y,result_hdr,tmp_hdr)
				
				
				#pragma omp parallel for schedule(dynamic, 1)
						for (int x = 0; x < width; x ++) {
							const float depth = depth_arr[y * width + x];
							
							//tigra: off blur
							//if (fabs(depth - focal_distance) >= blur_depth_threashold) 
							if (0)
							{
								float radius = (fabs(depth - focal_distance) / (blur_depth_threashold * 2.0f) - total_samples[y * width + x] / blur_total_samples_coeff);
								if (radius > blur_max_radius)
									radius = blur_max_radius;
								if (radius < 0)
									radius = 0;
								Color accum;
								float num = 0;
								for (int ix = -radius; ix <= radius; ++ix) {
									for (int iy = -radius; iy <= radius; ++iy) {
										float weight = sqrt((float)(ix * ix + iy * iy));
										if (weight > radius)
											continue;

										weight = 1.0f / (weight + 1.0f);

										const int nx = x + ix;
										const int ny = y + iy;
										if (nx < 0 || width <= nx)
											continue;
										if (ny < 0 || height <= ny)
											continue;

										if (fabs(depth_arr[ny * width + nx] - focal_distance) < blur_depth_threashold)
											continue;

										num += weight;
										accum = accum + weight * (*tmp_hdr.image_ptr(nx, ny));
									}
								}
								(*result_hdr.image_ptr(x, y)) = accum / num;
							} else {
								(*result_hdr.image_ptr(x, y)) = (*tmp_hdr.image_ptr(x, y));
							}
						}
					}
//#pragma omp critical
					// LDR化
					to_ldr(&result_hdr, fstop);

					// トーンカーブ調整
					tone_curve(&result_hdr, 
						Color((float)tone_curve_input_r, (float)tone_curve_input_g, (float)tone_curve_input_b) / 255.0f, 
						Color((float)tone_curve_output_r, (float)tone_curve_output_g, (float)tone_curve_output_b) / 255.0f);

					// 周辺光量
					vignetting(&result_hdr, height * vignetting_value);
			
					char str[256];
								
					sprintf(str, "result_%03d__%03d.hdr", iteration, count);
					result_hdr.save(str);
					
					sprintf(str, "tmp_hdr_%03d__%03d.hdr", iteration, count);
					tmp_hdr.save(str);
					
					sprintf(str, "hdr_%03d__%03d.hdr", iteration, count);
					hdr.save(str);
					

					// bmpに出力
					for (int iy = 0; iy < height; ++iy) {
						for (int ix = 0; ix < width; ++ix) {
							const unsigned char r = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->x_, display_gamma);
							const unsigned char g = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->y_, display_gamma);
							const unsigned char b = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->z_, display_gamma);
							bmp_arr[(iy * width + ix) * 3 + 0] = b;
							bmp_arr[(iy * width + ix) * 3 + 1] = g;
							bmp_arr[(iy * width + ix) * 3 + 2] = r;

//							printf("%f %d\n", result_hdr.image_ptr(ix, height - iy - 1)->x_, r);
						}
					}
					sprintf(str, "%03d__%03d.bmp", iteration, count);
					exportToBmp(str, &bmp_arr[0], width, height);

					std::cout << "output: " << str << std::endl;


					endTime1 = clock();
					now1 = (float)(endTime1 - beginTime0) / (float)CLOCKS_PER_SEC;

					printf("%.2fsec.\n", now1);

					++count;
				}
			}
			std::cout << std::endl;

			// 深度マップ計算
			for (int y = 0; y < height; y ++) {
				for (int x = 0; x < width; x ++) {
					depth_arr[y * width + x] =  (*depth_hdr.image_ptr(x, y)).x_ / total_samples[y * width + x];
				}
			}
			
			
			for (int y = 0; y < height; y ++) {
				for (int x = 0; x < width; x ++) {
					depth_arr[y * width + x] =  depth_arr[y * width + x] * 100.0;
				}
			}
					
			char str[256];

			sprintf(str, "depth_%03d__%03d.hdr", iteration, count);
					//depth_hdr.save(str);
					flo2hdr(str, depth_arr);

			// 次のサンプルマップを計算
			std::cout << "calculate variance_map" << std::endl;
			for (int y = 0; y < height; y ++) {
				for (int x = 0; x < width; x ++) {
					const float tts1 = 1.0f / total_samples[y * width + x];

					const float lumi = luminance_map[y * width + x] * tts1;
					float lumi0 = lumi;
					float lumi1 = lumi;
					if (x + 1 < width)
						lumi0 = luminance_map[y * width + (x + 1)] / total_samples[y * width + (x + 1)];
					if (y + 1 < height)
						lumi1 = luminance_map[(y + 1) * width + x] / total_samples[(y + 1) * width + x];

					const float dfdx = (lumi - lumi0);
					const float dfdy = (lumi - lumi1);

					variance_map[y * width + x] = sqrt(dfdx * dfdx + dfdy * dfdy) + tts1;
				}
			}

			// importance_mapをぼかす
			std::cout << "blur importance_map" << std::endl;
		#pragma omp parallel for schedule(dynamic, 1)
			for (int y = 0; y < height; y ++) {
				for (int x = 0; x < width; x ++) {
					float accum = 0.0;
					int num = 0;
					int radius = 5;
					for (int ix = -radius; ix <= radius; ++ix) {
						for (int iy = -radius; iy <= radius; ++iy) {
							const float weight = sqrt((float)(ix * ix + iy * iy));
							if (weight > radius)
								continue;

							const int nx = x + ix;
							const int ny = y + iy;
							if (nx < 0 || width <= nx)
								continue;
							if (ny < 0 || height <= ny)
								continue;
							num ++;
							accum = accum + variance_map[ny * width + nx];
						}
					}
					importance_map[y * width + x] = accum / num;
				}
			}
				
			float total = 0.0f;
			for (int i = 0; i < importance_map.size(); ++i)
				total += importance_map[i];

			const float ttl1 = 1.0f / total;

			for (int i = 0; i < importance_map.size(); ++i) {
				importance_map1[i] = importance_map[i];
				importance_map[i] *= ttl1;
			}
			
			std::cout << "calculate sample_map" << std::endl;
			// 各ピクセルのサンプル数を適当に決定
			const int geta = sample_per_pixel_per_iteration - 1;
			//const int sample_num = sample_per_pixel_per_iteration * width * height - (geta * width * height);
			const int sample_num = width * height;
			

			long ttl_s = 0;


			int samples_min, samples_max;
			
			for (int i = 0; i < sample_map.size(); ++i) {
				sample_map[i] = geta+1;

				const int now = (int)(sample_num * importance_map[i]);
				const float nowf = ((float)sample_num * importance_map[i] - now);

				sample_map[i] += now;
				if (grnd.next01() < nowf) {
					sample_map[i] ++;
				}
				
				ttl_s += sample_map[i];
				
				if(i)
				{
					if(sample_map[i]>samples_max) samples_max=sample_map[i];
					if(sample_map[i]<samples_min) samples_min=sample_map[i];
				}
				else
				{
					samples_max=sample_map[i];
					samples_min=sample_map[i];
				}
			}

			for (int i = 0; i < sample_map.size(); ++i)
			{
				sample_mapf[i] = (float)100.0f*sample_map[i]/samples_max;
			}
			
			printf("samples_min=%d samples_max=%d\n--------------------\n",samples_min,samples_max);

			//tigra: здесь можно выводить variance_map, importance_map, luminance_map - float,  sample_map - int

			//char str[256];

			sprintf(str, "variance_map_%03d__%03d.hdr", iteration, count);
			flo2hdr(str, variance_map);

			sprintf(str, "importance_map_%03d__%03d.hdr", iteration, count);
			flo2hdr(str, importance_map1);

			sprintf(str, "luminance_map_%03d__%03d.hdr", iteration, count);
			flo2hdr(str, luminance_map);

			sprintf(str, "sample_map_%03d__%03d.hdr", iteration, count);
			flo2hdr(str, sample_mapf);
			
			
					// LDR化
					to_ldr(&result_hdr, fstop);

					// トーンカーブ調整
					tone_curve(&result_hdr, 
						Color((float)tone_curve_input_r, (float)tone_curve_input_g, (float)tone_curve_input_b) / 255.0f, 
						Color((float)tone_curve_output_r, (float)tone_curve_output_g, (float)tone_curve_output_b) / 255.0f);

					// 周辺光量
					vignetting(&result_hdr, height * vignetting_value);
			
								
					sprintf(str, "result_%03d__%03d.hdr", iteration, count);
					result_hdr.save(str);
					
					sprintf(str, "tmp_hdr_%03d__%03d.hdr", iteration, count);
					tmp_hdr.save(str);
					
					sprintf(str, "hdr_%03d__%03d.hdr", iteration, count);
					hdr.save(str);
					

					// bmpに出力
					for (int iy = 0; iy < height; ++iy) {
						for (int ix = 0; ix < width; ++ix) {
							const unsigned char r = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->x_, display_gamma);
							const unsigned char g = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->y_, display_gamma);
							const unsigned char b = to_bmp_value(result_hdr.image_ptr(ix, height - iy - 1)->z_, display_gamma);
							bmp_arr[(iy * width + ix) * 3 + 0] = b;
							bmp_arr[(iy * width + ix) * 3 + 1] = g;
							bmp_arr[(iy * width + ix) * 3 + 2] = r;

//							printf("%f %d\n", result_hdr.image_ptr(ix, height - iy - 1)->x_, r);
						}
					}
					sprintf(str, "%03d__%03d.bmp", iteration, count);
					exportToBmp(str, &bmp_arr[0], width, height);
			
			
					endTime1 = clock();
					now1 = (float)(endTime1 - beginTime0) / (float)CLOCKS_PER_SEC;

					printf("Iteration time %.2fsec.\n", now1);
			
		}
	}

	//timeEndPeriod(1);
	
	
					endTime1 = clock();
					now1 = (float)(endTime1 - beginTime000) / (float)CLOCKS_PER_SEC;

					printf("all time %.2fsec.\n", now1);

	return 0;
}