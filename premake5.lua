workspace("lut_gen")
	configurations { "Debug", "Release" }
	
project("Transmittance")
	kind("ConsoleApp")
	language("C++")
	cppdialect("C++17")
	
	targetdir("bin")
	files {
		"src/*.*"
	}
	