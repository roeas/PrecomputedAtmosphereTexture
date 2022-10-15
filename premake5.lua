workspace("lut_gen")
	configurations { "Debug", "Release" }
	
	architecture "x64"
	
	-- Debug is a strict debug mode. No optimization will be performed.
	filter "configurations:Debug"
		defines { "DEBUG" }
		symbols("On")
		optimize("Off")
	-- Full optimization.
	filter "configurations:Release"
		defines { "NDEBUG" }
		symbols("Off")
		optimize("Full")
	filter {}
	
project("Transmittance")
	kind("ConsoleApp")
	language("C++")
	cppdialect("C++17")
	
	location("build")
	
	targetdir("bin")
	files {
		"src/**.**"
	}
	
	vpaths {
		["AtmosphereParameters/*"] = {
			"atmosphereParameters/**.*"
		},
		["Functions/*"] = { 
			"functions/**.*",
		},
		["Math/*"] = {
			"math/**.*"
		},
		["stb/*"] = { 
			"stb/**.*",
		},
	}
	
	includedirs {
		"src"
	}
	
	cwd = os.getcwd()
	postbuildcommands {
		"{MKDIR} "..path.join(cwd, "precomputed")
	}
	postbuildmessage ("creat precomputed path")
