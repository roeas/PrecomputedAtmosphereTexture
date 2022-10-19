#include <vector>
#include <cmath>
#include <fstream>

#include "functions/functions.h"
#include "atmosphereParameters/model.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

constexpr double kEpsilon = 1e-3;
constexpr SpectralIrradiance kSolarIrradiance = 123.0 * watt_per_square_meter_per_nm;
constexpr Length kBottomRadius = 1000.0 * km;
constexpr Length kTopRadius = 1500.0 * km;
constexpr Length kScaleHeight = 60.0 * km;
constexpr Length kRayleighScaleHeight = 60.0 * km;
constexpr Length kMieScaleHeight = 30.0 * km;
constexpr ScatteringCoefficient kRayleighScattering = 0.001 / km;
constexpr ScatteringCoefficient kMieScattering = 0.0015 / km;
constexpr ScatteringCoefficient kMieExtinction = 0.002 / km;
constexpr Number kGroundAlbedo = 0.1;

constexpr double kPi = 3.1415926;
constexpr double kSunAngularRadius = 0.00935 / 2.0;
constexpr double kSunSolidAngle = kPi * kSunAngularRadius * kSunAngularRadius;
constexpr double kLengthUnitInMeters = 1000.0;

void InitModel() {
    // Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
    // (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
    // summed and averaged in each bin (e.g. the value for 360nm is the average
    // of the ASTM G-173 values for all wavelengths between 360 and 370nm).
    // Values in W.m^-2.
    constexpr int kLambdaMin = 360;
    constexpr int kLambdaMax = 830;
    constexpr double kSolarIrradiance[48] = {
      1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
      1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
      1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
      1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
      1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
      1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
    };
    // Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
    // referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
    // each bin (e.g. the value for 360nm is the average of the original values
    // for all wavelengths between 360 and 370nm). Values in m^2.
    constexpr double kOzoneCrossSection[48] = {
      1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
      8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
      1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
      4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
      2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
      6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
      2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
    };
    // From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
    constexpr double kDobsonUnit = 2.687e20;
    // Maximum number density of ozone molecules, in m^-3 (computed so at to get
    // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
    // the ozone density profile defined below, which is equal to 15km).
    constexpr double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
    // Wavelength independent solar irradiance "spectrum" (not physically
    // realistic, but was used in the original implementation).
    constexpr double kConstantSolarIrradiance = 1.5;
    constexpr double kBottomRadius = 6360000.0;
    constexpr double kTopRadius = 6420000.0;
    constexpr double kRayleigh = 1.24062e-6;
    constexpr double kRayleighScaleHeight = 8000.0;
    constexpr double kMieScaleHeight = 1200.0;
    constexpr double kMieAngstromAlpha = 0.0;
    constexpr double kMieAngstromBeta = 5.328e-3;
    constexpr double kMieSingleScatteringAlbedo = 0.9;
    constexpr double kMiePhaseFunctionG = 0.8;
    constexpr double kGroundAlbedo = 0.1;
    const double max_sun_zenith_angle = 102.0 / 180.0 * kPi;

    DensityProfileLayer rayleigh_layer{ 0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0 };
    DensityProfileLayer mie_layer{ 0.0, 1.0, -1.0 / kMieScaleHeight, 0.0, 0.0 };
    // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
    // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
    // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
    // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
    std::vector<DensityProfileLayer> ozone_density;
    ozone_density.push_back(DensityProfileLayer{ 25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0 });
    ozone_density.push_back(DensityProfileLayer{ 0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0 });

    std::vector<double> wavelengths;
    std::vector<double> solar_irradiance;
    std::vector<double> rayleigh_scattering;
    std::vector<double> mie_scattering;
    std::vector<double> mie_extinction;
    std::vector<double> absorption_extinction;
    std::vector<double> ground_albedo;
    for (int l = kLambdaMin; l <= kLambdaMax; l += 10) {
        double lambda = static_cast<double>(l) * 1e-3;  // micro-meters
        double mie = kMieAngstromBeta / kMieScaleHeight * pow(lambda, -kMieAngstromAlpha);
        wavelengths.push_back(l);
        solar_irradiance.push_back(kConstantSolarIrradiance);

        rayleigh_scattering.push_back(kRayleigh * pow(lambda, -4));
        mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
        mie_extinction.push_back(mie);
        absorption_extinction.push_back(kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10]);
        ground_albedo.push_back(kGroundAlbedo);
    }

    Model model = Model(wavelengths, solar_irradiance, kSunAngularRadius,
        kBottomRadius, kTopRadius, { rayleigh_layer }, rayleigh_scattering,
        { mie_layer }, mie_scattering, mie_extinction, kMiePhaseFunctionG,
        ozone_density, absorption_extinction, ground_albedo, max_sun_zenith_angle,
        kLengthUnitInMeters, 3,
        false, false);
    model.PrintAtmParameter();
}

float data[TRANSMITTANCE_TEXTURE_WIDTH * TRANSMITTANCE_TEXTURE_HEIGHT * 3];

int main(int argc, char **argv) {
    // 初始化 Model 并打印 AtmosphereParameters 的初始化参数
    InitModel();

    const AtmosphereParameters ATMOSPHERE = AtmosphereParameters{
        Vec3d(1.500000, 1.500000, 1.500000),
        0.004675,
        6360.000000,
        6420.000000,
        DensityProfile{
            DensityProfileLayer{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 },
            DensityProfileLayer{ 0.000000, 1.000000, -0.125000, 0.000000, 0.000000 }
        },
        Vec3d(0.005802, 0.013558, 0.033100),
        DensityProfile{
            DensityProfileLayer{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 },
            DensityProfileLayer{ 0.000000, 1.000000, -0.833333, 0.000000, 0.000000 }
        },
        Vec3d(0.003996, 0.003996, 0.003996),
        Vec3d(0.004440, 0.004440, 0.004440),
        0.800000,
        DensityProfile{
            DensityProfileLayer{ 25.000000, 0.000000, 0.000000, 0.066667, -0.666667 },
            DensityProfileLayer{ 0.000000, 0.000000, 0.000000, -0.066667, 2.666667 }
        },
        Vec3d(0.000650, 0.001881, 0.000085),
        Vec3d(0.100000, 0.100000, 0.100000),
        -0.207912 };

    // 只计算 Transmittance 并将结果保存为二维的贴图
    for (int i = 0; i < TRANSMITTANCE_TEXTURE_HEIGHT; i++) {
        for (int j = 0; j < TRANSMITTANCE_TEXTURE_WIDTH; j++) {
            const Vec2d UV = { static_cast<double>(j), static_cast<double>(i) };
            const Vec3d trans = ComputeTransmittanceToTopAtmosphereBoundaryTexture(ATMOSPHERE, UV);
            //std::cout << trans.x << " " << trans.y << " " << trans.z << std::endl;
            const int pixelIndex = i * TRANSMITTANCE_TEXTURE_WIDTH + j;
            data[pixelIndex * 3 + 0] = static_cast<float>(trans.x);
            data[pixelIndex * 3 + 1] = static_cast<float>(trans.y);
            data[pixelIndex * 3 + 2] = static_cast<float>(trans.z);
        }
    }
    //stbi_flip_vertically_on_write(true);
    std::string outPutPath(argv[1]);
    outPutPath += "/LUT.hdr";
    stbi_write_hdr(outPutPath.c_str(), TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT, 3, data);

	return 0;
}

/*
const AtmosphereParameters ATMOSPHERE = AtmosphereParameters(
vec3(1.500000,1.500000,1.500000),
0.004675,
6360.000000,
6420.000000,
DensityProfile(DensityProfileLayer[2](DensityProfileLayer(0.000000,0.000000,0.000000,0.000000,0.000000),DensityProfileLayer(0.000000,1.000000,-0.125000,0.000000,0.000000))),
vec3(0.005802,0.013558,0.033100),
DensityProfile(DensityProfileLayer[2](DensityProfileLayer(0.000000,0.000000,0.000000,0.000000,0.000000),DensityProfileLayer(0.000000,1.000000,-0.833333,0.000000,0.000000))),
vec3(0.003996,0.003996,0.003996),
vec3(0.004440,0.004440,0.004440),
0.800000,
DensityProfile(DensityProfileLayer[2](DensityProfileLayer(25.000000,0.000000,0.000000,0.066667,-0.666667),DensityProfileLayer(0.000000,0.000000,0.000000,-0.066667,2.666667))),
vec3(0.000650,0.001881,0.000085),
vec3(0.100000,0.100000,0.100000),
-0.207912);
*/
