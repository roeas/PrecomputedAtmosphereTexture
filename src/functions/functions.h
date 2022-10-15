#pragma once

#include "atmosphereParameters/constants.h"
#include "atmosphereParameters/definitions.h"
#include "util.h"

// 距离大气顶部的距离
Length DistanceToTopAtmosphereBoundary(IN(AtmosphereParameters) atmosphere, Length r, Number mu) {
    assert(r <= atmosphere.top_radius);
    assert(mu >= -1.0 && mu <= 1.0);
    const Area discriminant = r * r * (mu * mu - 1.0) + atmosphere.top_radius * atmosphere.top_radius;
    return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

// 距离地表的距离
Length DistanceToBottomAtmosphereBoundary(IN(AtmosphereParameters) atmosphere, Length r, Number mu) {
    assert(r >= atmosphere.bottom_radius);
    assert(mu >= -1.0 && mu <= 1.0);
    const Area discriminant = r * r * (mu * mu - 1.0) + atmosphere.bottom_radius * atmosphere.bottom_radius;
    return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

// 射线是否与地球相交
bool RayIntersectsGround(IN(AtmosphereParameters) atmosphere, Length r, Number mu) {
    assert(r >= atmosphere.bottom_radius);
    assert(mu >= -1.0 && mu <= 1.0);
    return (mu < 0.0) && (r * r * (mu * mu - 1.0) + atmosphere.bottom_radius * atmosphere.bottom_radius >= 0.0f * m2);
}

// 返回指定大气层级的大气密度
// 在计算非臭氧大气粒子时该公式退化为 Number = exp(layer.exp_scale * altitude);
// 在计算臭氧时该公式退化为 Numer = layer.linear_term * altitude + layer.constant_term;
Number GetLayerDensity(IN(DensityProfileLayer) layer, Length altitude) {
    const Number density = layer.exp_term * exp(layer.exp_scale * altitude) + layer.linear_term * altitude + layer.constant_term;
    return clamp(density, Number(0.0), Number(1.0));
}

// 多层大气，根据海拔选择不同的大气层级
Number GetProfileDensity(IN(DensityProfile) profile, Length altitude) {
    return altitude < profile.layers[0].width ?
        GetLayerDensity(profile.layers[0], altitude) :
        GetLayerDensity(profile.layers[1], altitude);
}

// 从射线起点到其与大气顶部交点之间的光学距离
Length ComputeOpticalLengthToTopAtmosphereBoundary(IN(AtmosphereParameters) atmosphere, IN(DensityProfile) profile, Length r, Number mu) {
    assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
    assert(mu >= -1.0 && mu <= 1.0);
    // 采样次数
    const int SAMPLE_COUNT = 500;
    // 积分的长度
    const Length dx = DistanceToTopAtmosphereBoundary(atmosphere, r, mu) / Number(SAMPLE_COUNT);
    Length result = 0.0 * m;
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
        const Length d_i = Number(i) * dx;
        // 当前采样点与行星中心的距离
        const Length r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
        // 采样点的大气密度（除以大气底部的大气密度）
        const Number y_i = GetProfileDensity(profile, r_i - atmosphere.bottom_radius);
        // 采样权重
        const Number weight_i = i == 0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
        result += y_i * weight_i * dx;
    }
    return result;
}

// 从射线起点到其与大气顶部交点之间的透射率，考虑瑞利散射、米氏散射、吸收
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundary(IN(AtmosphereParameters) atmosphere, Length r, Number mu) {
    assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
    assert(mu >= -1.0 && mu <= 1.0);
    // rayleigh_scattering == rayleigh_extinction，瑞利散射不吸收光
    const DimensionlessSpectrum rayleighTerm = atmosphere.rayleigh_scattering * ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.rayleigh_density, r, mu);
    const DimensionlessSpectrum mieTerm = atmosphere.mie_extinction * ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.mie_density, r, mu);
    const DimensionlessSpectrum ozoneTerm = atmosphere.absorption_extinction * ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.absorption_density, r, mu);
    return exp(-(rayleighTerm + mieTerm + ozoneTerm));
}

// 参数化重映射
Number GetUnitRangeFromTextureCoord(Number u, int texture_size) {
    return (u - 0.5 / Number(texture_size)) / (1.0 - 1.0 / Number(texture_size));
}

 // UV -> RMu
void GetRMuFromTransmittanceTextureUv(IN(AtmosphereParameters) atmosphere, IN(Vec2d) uv, OUT(Length) r, OUT(Number) mu) {
    assert(uv.x >= 0.0 && uv.x <= 1.0);
    assert(uv.y >= 0.0 && uv.y <= 1.0);
    Number x_mu = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    Number x_r = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);
    // 射向地平线的射线，从地表到大气顶部的距离
    Length H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    // 射向地平线的射线，从起点到地表的距离
    Length rho = H * x_r;
    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    // 射线起点到顶部大气边界的距离，及其在所有 mu 上的最小值（r, 1）和最大值（r, mu_horizon）
    Length d_min = atmosphere.top_radius - r;
    Length d_max = rho + H;
    Length d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 * m ? Number(1.0) : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

// TRANSMITTANCE
DimensionlessSpectrum ComputeTransmittanceToTopAtmosphereBoundaryTexture(IN(AtmosphereParameters) atmosphere, IN(Vec2d) frag_coord) {
    const Vec2d TRANSMITTANCE_TEXTURE_SIZE = Vec2d(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    Length r;
    Number mu;
    GetRMuFromTransmittanceTextureUv(atmosphere, frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, mu);
    return ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, mu);
}
