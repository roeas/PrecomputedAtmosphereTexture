#include "model.h"

// 使用目标波长（wavelength）在所有波长（wavelengths）中匹配对应的数据（wavelength_function）
// TODO：太麻烦了，我们只用RGB，回头得把这套光谱渲染剔除掉
double Interpolate(
    const std::vector<double> &wavelengths,
    const std::vector<double> &wavelength_function,
    double wavelength) {
    assert(wavelength_function.size() == wavelengths.size());
    if (wavelength < wavelengths[0]) {
        return wavelength_function[0];
    }
    for (unsigned int i = 0; i < wavelengths.size() - 1; ++i) {
        if (wavelength < wavelengths[i + 1]) {
            // 全光谱表（wavelengths）默认以 10nm 为单位步进，这里以具体的波长在表中做插值
            double u = (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
            return wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
        }
    }
    return wavelength_function[wavelength_function.size() - 1];
}

Model::Model(
    const std::vector<double> &wavelengths,
    const std::vector<double> &solar_irradiance,
    const double sun_angular_radius,
    double bottom_radius,
    double top_radius,
    const std::vector<DensityProfileLayer> &rayleigh_density,
    const std::vector<double> &rayleigh_scattering,
    const std::vector<DensityProfileLayer> &mie_density,
    const std::vector<double> &mie_scattering,
    const std::vector<double> &mie_extinction,
    double mie_phase_function_g,
    const std::vector<DensityProfileLayer> &absorption_density,
    const std::vector<double> &absorption_extinction,
    const std::vector<double> &ground_albedo,
    double max_sun_zenith_angle,
    double length_unit_in_meters,
    unsigned int num_precomputed_wavelengths,
    bool combine_scattering_textures,
    bool half_precision) :

    num_precomputed_wavelengths_(num_precomputed_wavelengths){

    auto to_string = [&wavelengths](const std::vector<double> &v, const vec3 &lambdas, double scale) {
            double r = Interpolate(wavelengths, v, lambdas[0]) * scale;
            double g = Interpolate(wavelengths, v, lambdas[1]) * scale;
            double b = Interpolate(wavelengths, v, lambdas[2]) * scale;
            return "vec3(" + std::to_string(r) + "," + std::to_string(g) + "," + std::to_string(b) + ")";
    };
    auto density_layer =
        [length_unit_in_meters](const DensityProfileLayer &layer) {
        return "DensityProfileLayer(" +
            std::to_string(layer.width / length_unit_in_meters) + "," +
            std::to_string(layer.exp_term) + "," +
            std::to_string(layer.exp_scale * length_unit_in_meters) + "," +
            std::to_string(layer.linear_term * length_unit_in_meters) + "," +
            std::to_string(layer.constant_term) + ")";
    };
    auto density_profile =
        [density_layer](std::vector<DensityProfileLayer> layers) {
        constexpr int kLayerCount = 2;
        while (layers.size() < kLayerCount) {
            layers.insert(layers.begin(), DensityProfileLayer());
        }
        std::string result = "DensityProfile(DensityProfileLayer[" +
            std::to_string(kLayerCount) + "](";
        for (int i = 0; i < kLayerCount; ++i) {
            result += density_layer(layers[i]);
            result += i < kLayerCount - 1 ? "," : "))";
        }
        return result;
    };

    // 原代码中使用该 Lambda 生成 GLSL 代码
    glsl_header_factory_ = [=](const vec3 &lambdas) {
        return
            "const AtmosphereParameters ATMOSPHERE = AtmosphereParameters(\n" +
            to_string(solar_irradiance, lambdas, 1.0) + ",\n" +
            std::to_string(sun_angular_radius) + ",\n" +
            std::to_string(bottom_radius / length_unit_in_meters) + ",\n" +
            std::to_string(top_radius / length_unit_in_meters) + ",\n" +
            density_profile(rayleigh_density) + ",\n" +
            to_string(rayleigh_scattering, lambdas, length_unit_in_meters) + ",\n" +
            density_profile(mie_density) + ",\n" +
            to_string(mie_scattering, lambdas, length_unit_in_meters) + ",\n" +
            to_string(mie_extinction, lambdas, length_unit_in_meters) + ",\n" +
            std::to_string(mie_phase_function_g) + ",\n" +
            density_profile(absorption_density) + ",\n" +
            to_string(absorption_extinction, lambdas, length_unit_in_meters) + ",\n" +
            to_string(ground_albedo, lambdas, 1.0) + ",\n" +
            std::to_string(cos(max_sun_zenith_angle)) + ");\n";
    };
}

void Model::PrintAtmParameter() {
    std::string str = glsl_header_factory_({ kLambdaR , kLambdaG , kLambdaB });
    std::cout << str << std::endl;
}