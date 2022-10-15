#pragma once

#include "math/vec.h"

#define Length                    double
#define Wavelength                double
#define Angle                     double
#define SolidAngle                double
#define Power                     double
#define LuminousPower             double

#define Number                    double
#define InverseLength             double
#define Area                      double
#define Volume                    double
#define NumberDensity             double
#define Irradiance                double
#define Radiance                  double
#define SpectralPower             double
#define SpectralIrradiance        double
#define SpectralRadiance          double
#define SpectralRadianceDensity   double
#define ScatteringCoefficient     double
#define InverseSolidAngle         double
#define LuminousIntensity         double
#define Luminance                 double
#define Illuminance               double

// 波长到其他类型的通用函数
#define AbstractSpectrum          Vec3d
// 波长到数字的函数				  
#define DimensionlessSpectrum     Vec3d
// 波长到光谱 Power 的函数		  
#define PowerSpectrum             Vec3d
// 波长到光谱 Irradiance 的函数	  
#define IrradianceSpectrum        Vec3d
// 波长到光谱 Radiance 的函数	  
#define RadianceSpectrum          Vec3d
// 波长到光谱 Radiance 密度的函数 
#define RadianceDensitySpectrum   Vec3d
// 波长到散射系数				   
#define ScatteringSpectrum        Vec3d

#define Position                  Vec3d
#define Direction                 Vec3d
#define Luminance3                Vec3d
#define Illuminance3              Vec3d

#define sampler2D int
#define sampler3D int
#define TransmittanceTexture      sampler2D
#define AbstractScatteringTexture sampler3D
#define ReducedScatteringTexture  sampler3D
#define ScatteringTexture         sampler3D
#define ScatteringDensityTexture  sampler3D
#define IrradianceTexture         sampler2D

// 物理单位
constexpr Length m = 1.0;
constexpr Wavelength nm = 1.0;
constexpr Angle rad = 1.0;
constexpr SolidAngle sr = 1.0;
constexpr Power watt = 1.0;
constexpr LuminousPower lm = 1.0;

// 衍生出的物理单位
constexpr double PI = 3.14159265358979323846;
constexpr Length km = 1000.0 * m;
constexpr Area m2 = m * m;
constexpr Volume m3 = m * m * m;
constexpr Angle pi = PI * rad;
constexpr Angle deg = pi / 180.0;
constexpr Irradiance watt_per_square_meter = watt / m2;
constexpr Radiance watt_per_square_meter_per_sr = watt / (m2 * sr);
constexpr SpectralIrradiance watt_per_square_meter_per_nm = watt / (m2 * nm);
constexpr SpectralRadiance watt_per_square_meter_per_sr_per_nm = watt / (m2 * sr * nm);
constexpr SpectralRadianceDensity watt_per_cubic_meter_per_sr_per_nm = watt / (m3 * sr * nm);
constexpr LuminousIntensity cd = lm / sr;
constexpr LuminousIntensity kcd = 1000.0 * cd;
constexpr Luminance cd_per_square_meter = cd / m2;
constexpr Luminance kcd_per_square_meter = kcd / m2;

// 宽度为 width 的大气层，
// 密度由 'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term' 定义，
// clamp 至 [0, 1]，h 为海拔高度
// 只有臭氧会用到 linear_term 和 constant_term
struct DensityProfileLayer {
	Length width;
	// 非臭氧大气粒子为 1，臭氧为 0
	Number exp_term;
	// 负标准海拔高度分之一
	InverseLength exp_scale;
	InverseLength linear_term;
	Number constant_term;
};

// 一个从下到上由若干层叠加而成的大气密度分布，
// 最后一层的宽度被忽略，也就是说总是延伸到大气顶部的边界，
// 密度分布值在 [0, 1] 之间变化
// 这里双层的设计仅仅服务于臭氧，其他的大气粒子无论什么海拔都只使用 layers[1] 的数据
struct DensityProfile {
	DensityProfileLayer layers[2];
};

struct AtmosphereParameters {
	// 大气顶部太阳的 Irradiance
	IrradianceSpectrum solar_irradiance;
	// 太阳的角半径
	// 注意：实现中使用的近似值只有当该角度小于 0.1 时才有效
	Angle sun_angular_radius;
	// 地球半径
	Length bottom_radius;
	// 地球半径 + 大气高度
	Length top_radius;

	// 空气分子的密度分布，即由海拔高度得到 [0, 1] 之间值的函数
	DensityProfile rayleigh_density;
	// 空气分子在其密度最大的海拔高度（通常是大气底部）的散射系数，是波长的函数
	// 海拔高度 h 处的散射系数为 'rayleigh_scattering' 乘以 'rayleigh_density'
	ScatteringSpectrum rayleigh_scattering;

	// 气溶胶的密度分布，即由海拔高度得到 [0, 1] 之间值的函数
	DensityProfile mie_density;
	// 气溶胶在其密度最大的海拔高度（通常是大气底部）的散射系数，是波长的函数
	// 海拔高度 h 处的散射系数为 'mie_scattering' 乘以 'mie_density'
	ScatteringSpectrum mie_scattering;
	// 气溶胶在其密度最大的海拔高度（通常是大气底部）的消光系数，是波长的函数
	// 海拔高度 h 处的消光系数为 'mie_extinction' 乘以 'mie_density'
	ScatteringSpectrum mie_extinction;
	// Cornette-Shanks 气溶胶相函数的非对称性参数
	Number mie_phase_function_g;

	// 会吸收光的气体分子（比如臭氧）的密度分布，即由海拔高度得到 [0, 1] 之间值的函数
	DensityProfile absorption_density;
	// 吸收光的分子（比如臭氧）在其密度最大的海拔高度时的消光系数，时波长的函数
	// 海拔高度 h 处的消光系数为 'absorption_extinction' 乘以 'absorption_density'
	ScatteringSpectrum absorption_extinction;
	
	// 地表的平均 albedo
	DimensionlessSpectrum ground_albedo;
	// 必须预计算大气散射的最大太阳天顶角
	// 为了获得最大精度，使用产生可忽略的天光 radiance 的最小太阳天顶角。
	// 例如，对于地球的情况，102 度是一个很好的选择 - 产生 mu_s_min = -0.2
	Number mu_s_min;
};
