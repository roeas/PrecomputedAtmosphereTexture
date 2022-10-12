#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>
#include "definitions.h"

// 在当前的实现中，该类仅仅是为了获得用于初始化 AtmosphereParameters 的正确参数
class Model {
public:
    Model(
        // 提供太阳辐射度、雷利散射、米氏散射、米氏消光和地面消光样本的波长值，单位是纳米，并按递增顺序排序
        const std::vector<double> &wavelengths,
        // The solar irradiance at the top of the atmosphere, in W/m^2/nm. This
        // 大气顶部的太阳 irradiance，该 vetor 必须与 wavelengths 相同大小
        const std::vector<double> &solar_irradiance,
        // 太阳的角半径，单位是弧度
        double sun_angular_radius,
        // 地球半径
        double bottom_radius,
        // 地球半径 + 大气高度
        double top_radius,
        // 空气分子密度分布
        const std::vector<DensityProfileLayer> &rayleigh_density,
        // 空气分子散射系数
        const std::vector<double> &rayleigh_scattering,
        // 气溶胶密度分布
        const std::vector<DensityProfileLayer> &mie_density,
        // 气溶胶散射系数
        const std::vector<double> &mie_scattering,
        // 气溶胶消光系数
        const std::vector<double> &mie_extinction,
        // 气溶胶像函数非对称性系数
        double mie_phase_function_g,
        // 臭氧吸密度分布
        const std::vector<DensityProfileLayer> &absorption_density,
        // 臭氧吸光系数
        const std::vector<double> &absorption_extinction,
        // 地面的平均 albedo，该 vetor 必须与 wavelengths 相同大小
        const std::vector<double> &ground_albedo,
        // 最大的地平线天顶角
        double max_sun_zenith_angle,
        // 模型/坐标系的单位，即一个整体的缩放
        double length_unit_in_meters,
        // 小于等于3使用 {kLambdaR, kLambdaG, kLambdaB} 进行渲染，否则拆分光谱
        unsigned int num_precomputed_wavelengths,
        // 打包 single Mie.r 和 Rayleigh and multiple scattering
        bool combine_scattering_textures,
        // 半精度
        bool half_precision);

    static constexpr double kLambdaR = 680.0;
    static constexpr double kLambdaG = 550.0;
    static constexpr double kLambdaB = 440.0;

    // 将生成的 GLSL 代码打印出来，手动初始化 AtmosphereParameters
    void PrintAtmParameter();

private:
    typedef std::array<double, 3> vec3;
    typedef std::array<float, 9> mat3;

    unsigned int num_precomputed_wavelengths_;
    bool half_precision_;
    std::function<std::string(const vec3 &)> glsl_header_factory_;
    int transmittance_texture_;
    int scattering_texture_;
    int optional_single_mie_scattering_texture_;
    int irradiance_texture_;
    int atmosphere_shader_;
    int full_screen_quad_vao_;
    int full_screen_quad_vbo_;
};
