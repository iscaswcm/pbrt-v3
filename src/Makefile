###########################################################################
# user-configurable section
###########################################################################

ECHO=echo
# common locations for the OpenEXR libraries; may need to be updated
# for unusual installation locations
#

OPENEXR_INC    =-Iext/openexr/IlmBase \
				-Iext/openexr/IlmBase/Iex \
			    -Iext/openexr/IlmBase/Half \
			    -Iext/openexr/IlmBase/IexMath \
			    -Iext/openexr/IlmBase/Imath \
			    -Iext/openexr/IlmBase/IlmThread \
	    	    -Iext/openexr/IlmBase/config \
                -Iext/openexr/OpenEXR \
			    -Iext/openexr/OpenEXR/IlmImf \
		        -Iext/openexr/OpenEXR/config \
	    	    -Iext/openexr/OpenEXR/IlmImfUtil \
				-Iext/openexr/include
GLOG_INC       =-Iext/glog/src \
				-Iext/glog/src/glog \
				-Iext/glog/src/base \
				-Iext/glog/src/windows
ZLIB_INC       =-Iext/zlib
LIB_GENERAL_INC=-Iext

THIRD_PARTY_INC= $(ZLIB_INC) $(LIB_JPEG_INC) $(OPENEXR_INC) \
				 $(OPENEXR_INC) $(LIB_TARGA_INC) $(GLOG_INC) $(LIB_GENERAL_INC)
#libsolarmodel.a
LIB_SOLARMODEL_HEADERS= ext/ArHosekSkyModelData_CIEXYZ.h ext/ArHosekSkyModelData_RGB.h ext/ArHosekSkyModelData_Spectral.h ext/ArHosekSkyModel.h
LIB_SOLARMODEL_SRCS=ext/ArHosekSkyModel.c
LIB_SOLARMODEL_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_SOLARMODEL_SRCS:.c=.o))) 

#liblodepng
LIB_LODEPNG_HEADERS=ext/lodepng.h
LIB_LODEPNG_SRCS=ext/lodepng.cpp
LIB_LODEPNG_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_LODEPNG_SRCS:.cpp=.o))) 

#librply
LIB_RPLY_HEADERS=ext/rply.h
LIB_RPLY_SRCS=ext/rply.cpp
LIB_RPLY_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_RPLY_SRCS:.cpp=.o))) 
#zlib
LIB_ZLIB_HEADERS=$(wildcard ext/zlib/*.h )
LIB_ZLIB_SRCS=$(wildcard ext/zlib/*.c)
LIB_ZLIB_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ZLIB_SRCS:.c=.o))) 


#libjpeg
LIB_JPEG_HEADERS=$(wildcard ext/jpeg/*.h )
LIB_JPEG_SRCS=$(wildcard ext/jpeg/*.c)
LIB_JPEG_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_JPEG_SRCS:.c=.o))) 
#ilmbase
#ilmbase Half
LIB_ILMBASE_HALF_HEADERS=$(wildcard ext/openexr/IlmBase/Half/*.h )
LIB_ILMBASE_HALF_SRCS=$(wildcard ext/openexr/IlmBase/Half/*.cpp)
LIB_ILMBASE_HALF_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ILMBASE_HALF_SRCS:.cpp=.o))) 

#ilmbase libIex
LIB_ILMBASE_IEX_HEADERS=$(wildcard ext/openexr/IlmBase/Iex/*.h )
LIB_ILMBASE_IEX_SRCS=$(wildcard ext/openexr/IlmBase/Iex/*.cpp)
LIB_ILMBASE_IEX_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ILMBASE_IEX_SRCS:.cpp=.o))) 

#ilmbase IexMath
LIB_ILMBASE_IEXMATH_HEADERS=$(wildcard ext/openexr/IlmBase/IexMath/*.h )
LIB_ILMBASE_IEXMATH_SRCS=$(wildcard ext/openexr/IlmBase/IexMath/*.cpp)
LIB_ILMBASE_IEXMATH_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ILMBASE_IEXMATH_SRCS:.cpp=.o))) 

#ilmbase libIlmThread
LIB_ILMBASE_ILMTHREAD_HEADERS=$(wildcard ext/openexr/IlmBase/IlmThread/*.h )
LIB_ILMBASE_ILMTHREAD_SRCS=$(wildcard ext/openexr/IlmBase/IlmThread/*.cpp)
LIB_ILMBASE_ILMTHREAD_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ILMBASE_ILMTHREAD_SRCS:.cpp=.o))) 

#ilmbase Imath
LIB_ILMBASE_IMATH_HEADERS=$(wildcard ext/openexr/IlmBase/Imath/*.h )
LIB_ILMBASE_IMATH_SRCS=$(wildcard ext/openexr/IlmBase/Imath/*.cpp)
LIB_ILMBASE_IMATH_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_ILMBASE_IMATH_SRCS:.cpp=.o))) 
#ilmbase

#openexr
#openexr IlmImf
LIB_OPENEXR_ILMIMF_SRCS=$(wildcard ext/openexr/OpenEXR/IlmImf/*.cpp)
LIB_OPENEXR_ILMIMF_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_OPENEXR_ILMIMF_SRCS:.cpp=.o))) 
LIB_OPENEXR_ILMIMF_HEADERS=$(wildcard ext/openexr/OpenEXR/IlmImf/*.h )

#openexr IlmImfUTIL
LIB_OPENEXR_ILMIMFUTIL_HEADERS=$(wildcard ext/openexr/OpenEXR/IlmImfUtil/*.h )
LIB_OPENEXR_ILMIMFUTIL_SRCS=$(wildcard ext/openexr/OpenEXR/IlmImfUtil/*.cpp)
LIB_OPENEXR_ILMIMFUTIL_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_OPENEXR_ILMIMFUTIL_SRCS:.cpp=.o))) 

LIB_OPENEXR_OBJS = $(LIB_ILMBASE_HALF_OBJS) \
				  $(LIB_ILMBASE_IEX_OBJS) \
				  $(LIB_ILMBASE_IEXMATH_OBJS) \
				  $(LIB_ILMBASE_ILMTHREAD_OBJS) \
				  $(LIB_ILMBASE_IMATH_OBJS) \
				  $(LIB_OPENEXR_ILMIMF_OBJS) \
                  $(LIB_OPENEXR_ILMIMFUTIL_OBJS) 

#libtarga.a
LIB_TARGA_HEADERS=ext/targa.h
LIB_TARGA_SRCS=ext/targa.cpp
LIB_TARGA_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_TARGA_SRCS:.cpp=.o))) 

#glog
LIB_GLOG_SRCS= ext/glog/src/demangle.cc \
			   ext/glog/src/logging.cc \
			   ext/glog/src/raw_logging.cc \
			   ext/glog/src/signalhandler.cc \
			   ext/glog/src/symbolize.cc \
			   ext/glog/src/utilities.cc \
			   ext/glog/src/vlog_is_on.cc
#LIB_GLOG_SRCS=$(wildcard ext/glog/src/*.cc)
LIB_GLOG_OBJS=$(addprefix ../objs/, $(subst /,_,$(LIB_GLOG_SRCS:.cc=.o))) 
LIB_GLOG_HEADERS=$(wildcard ext/glog/*.h )

# 32 bit
#MARCH=-m32 -msse2 -mfpmath=sse

OS_TYPE =	$(shell uname -m | grep 64)
ifeq "$(strip $(OS_TYPE))" ""
	MARCH=-m32
else
	MARCH=-m64
endif
# 64 bit

# change this to -g3 for debug Buildings
OPT=-O3
# comment out this line to enable assertions at runtime

#########################################################################
# nothing below this line should need to be changed (usually)
#########################################################################

ARCH = $(shell uname)
#make -D parameter
#macro OPEXR_PLATFORM_WINDOWS only for MICROSOFT WINDOWS PLATFORM
ifeq ($(ARCH),Darwin)
	DEFS = -DMACOS_APPLE -DGOOGLE_GLOG_DLL_DECL="" 
else
	DEFS = -DGOOGLE_GLOG_DLL_DECL="" -DOPEXR_PLATFORM_WINDOWS0 -DNDEBUG=""
endif
LEX=flex
YACC=bison -d -v -t
SDLFLAGS = -I/usr/local/sdl2-2.0.5/include  -L/usr/local/sdl2-2.0.5/lib
CC = gcc  
CXX= g++  -std=c++11  -fPIC

INCLUDE=-I. -Icore $(THIRD_PARTY_INC) 
WARN= -w 

CXXFLAGS=$(OPT) $(MARCH) $(INCLUDE) $(WARN) $(DEFS) $(SDLFLAGS)
CCFLAGS=$(CXXFLAGS)
LINK_LIBS= -lpthread -lSDL2  -ll
LIB_CXXSRCS  = $(wildcard core/*.cpp) core/pbrtlex.cpp core/pbrtparse.cpp
LIB_CXXSRCS += $(wildcard accelerators/*.cpp cameras/*.cpp film/*.cpp filters/*.cpp )
LIB_CXXSRCS += $(wildcard integrators/*.cpp lights/*.cpp materials/*.cpp media/*.cpp renderers/*.cpp )
LIB_CXXSRCS += $(wildcard samplers/*.cpp shapes/*.cpp textures/*.cpp volumes/*.cpp)

LIBOBJS  = $(addprefix ../objs/, $(subst /,_,$(LIB_CSRCS:.c=.o)))
LIBOBJS += $(addprefix ../objs/, $(subst /,_,$(LIB_CXXSRCS:.cpp=.o)))

HEADERS  = $(wildcard */*.h)

PBRT_TEST_CXXSRCS = $(wildcard tests/*.cpp)
PBRT_TEST_CCSRCS  = $(wildcard tests/gtest/*.cc)
PBRT_TEST_OBJS    = $(addprefix ../objs/, $(subst /,_,$(PBRT_TEST_CXXSRCS:.cpp=.o)))
PBRT_TEST_OBJS   += $(addprefix ../objs/, $(subst /,_,$(PBRT_TEST_CCSRCS:.cc=.o)))
#############################################pbrt_test############################################
THIRD_PARTY_LIBS :=../libs/libopenexr.a \
		           ../libs/libjpeg.a \
		           ../libs/libzlib.a  \
		           ../libs/liblodepng.a \
		           ../libs/libtarga.a \
		           ../libs/librply.a \
		           ../libs/libsolarmodel.a \
		           ../libs/libglog.a 
###############################################pbrt###############################################
default: dirs ../bin/pbrt  ../bin/bsdftest ../bin/cyhair2pbrt ../bin/obj2pbrt #../bin/imgtool ../bin/pbrt_test
../bin/%: dirs
dirs:
	@echo "\033[0;32;33mMakeing directories needed....done\033[0m"
	mkdir -p ../bin ../objs  ../libs

.PHONY: dirs

###############################################pbrt###############################################
../bin/pbrt: ../objs/main_pbrt.o ../libs/libpbrt.a  $(THIRD_PARTY_LIBS)
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK_LIBS) 

../bin/pbrt_test: $(PBRT_TEST_OBJS) ../libs/libpbrt.a $(THIRD_PARTY_LIBS) 
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK_LIBS)

../bin/bsdftest: ../objs/tools_bsdftest.o ../libs/libpbrt.a $(THIRD_PARTY_LIBS)
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK_LIBS)

../bin/cyhair2pbrt: ../objs/tools_cyhair2pbrt.o ../libs/libpbrt.a $(THIRD_PARTY_LIBS)
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"

../bin/obj2pbrt: ../objs/tools_obj2pbrt.o ../libs/libpbrt.a $(THIRD_PARTY_LIBS) 
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK_LIBS)

../bin/imgtool: ../objs/tools_imgtool.o ../libs/libpbrt.a  $(THIRD_PARTY_LIBS)
	@echo "\033[44;37mLinking object\033[0m \033[32;39;1m$@\033[0m with \033[0;32;33m$^\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK_LIBS)

../objs/tools_%.o: tools/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../libs/libpbrt.a: $(LIBOBJS) 
	@echo "\033[0;32;33mBuildinging the core rendering library (libpbrt.a)\033[0m"
	@ar rcs $@ $(LIBOBJS)

../objs/media_%.o: media/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/accelerators_%.o: accelerators/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/cameras_%.o: cameras/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/core_%.o: core/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/core_%.o: core/%.c
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CC) $(CCFLAGS) -o $@ -c $<

../objs/film_%.o: film/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/filters_%.o: filters/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/integrators_%.o: integrators/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/lights_%.o: lights/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/main_%.o: main/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/materials_%.o: materials/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/renderers_%.o: renderers/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/samplers_%.o: samplers/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/shapes_%.o: shapes/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/textures_%.o: textures/%.cpp
	@echo "\033[32;49;1mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

core/pbrtlex.cpp: core/pbrtlex.ll core/pbrtparse.cpp
	@echo "\033[32;49;1mLEX'ing pbrtlex.ll\033[0m"
	@$(LEX) -o$@ core/pbrtlex.ll

core/pbrtparse.cpp: core/pbrtparse.y
	@echo "\033[32;49;1mYACC'ing pbrtparse.y\033[0m"
	@$(YACC) -o $@ core/pbrtparse.y
	@if [ -e core/pbrtparse.cpp.h ]; then /../bin/mv core/pbrtparse.cpp.h core/pbrtparse.h; fi
	@if [ -e core/pbrtparse.hpp ]; then /../bin/mv core/pbrtparse.hpp core/pbrtparse.h; fi
#for pbrt_test
../objs/tests_%.o: tests/%.cpp
	@echo "\033[32;49;34mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

../objs/tests_gtest_%.o: tests/gtest/%.cc
	@echo "\033[32;49;34mBuilding object $@\033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#######################################################################################
#librply
$(LIB_RPLY_OBJS): $(LIB_RPLY_HEADERS)
../libs/librply.a: $(LIB_RPLY_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core rply library librply.a\033[0m"
	@ar rcs $@ $(LIB_RPLY_OBJS)
../objs/ext_rply.o: ext/rply.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#libzlib
$(LIB_ZLIB_OBJS): $(LIB_ZLIB_HEADERS)
../libs/libzlib.a: $(LIB_ZLIB_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core libzlib library libzlib.a\033[0m"
	@ar rcs $@ $(LIB_ZLIB_OBJS)
../objs/ext_zlib_%.o: ext/zlib/%.c
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CC) $(CXXFLAGS) -o $@ -c $<

#libjpeg
$(LIB_JPEG_OBJS): $(LIB_JPEG_HEADERS)
../libs/libjpeg.a: $(LIB_JPEG_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core jpeg library libjpeg.a\033[0m"
	@ar rcs $@ $(LIB_JPEG_OBJS)
../objs/ext_jpeg_%.o: ext/jpeg/%.c
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase
#ilmbase Half
../objs/ext_openexr_IlmBase_Half_%.o: ext/openexr/IlmBase/Half/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase Iex
../objs/ext_openexr_IlmBase_Iex_%.o: ext/openexr/IlmBase/Iex/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase IexMath
../objs/ext_openexr_IlmBase_IexMath_%.o: ext/openexr/IlmBase/IexMath/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase libIlmThread
../objs/ext_openexr_IlmBase_IlmThread_%.o: ext/openexr/IlmBase/IlmThread/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase Imath
../objs/ext_openexr_IlmBase_Imath_%.o: ext/openexr/IlmBase/Imath/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#ilmbase

#libopenexr
../objs/ext_openexr_OpenEXR_IlmImfUtil_%.o: ext/openexr/OpenEXR/IlmImfUtil/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#openexr libopenexr.a
$(LIB_OPENEXR_ILMIMF_OBJS): $(LIB_OPENEXR_ILMIMF_HEADERS)
../libs/libopenexr.a: $(LIB_OPENEXR_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core OpenEXR library libopenexr.a\033[0m"
	@ar rcs $@ $(LIB_OPENEXR_OBJS)
../objs/ext_openexr_OpenEXR_IlmImf_%.o: ext/openexr/OpenEXR/IlmImf/%.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#libtarga
$(LIB_TARGA_OBJS): $(LIB_TARGA_HEADERS)
../libs/libtarga.a: $(LIB_TARGA_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core targa library libtarga.a\033[0m"
	@ar rcs $@ $(LIB_TARGA_OBJS)
../objs/ext_targa.o: ext/targa.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#liblodepng
$(LIB_LODEPNG_OBJS): $(LIB_LODEPNG_HEADERS)
../libs/liblodepng.a: $(LIB_LODEPNG_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core lodepng library liblodepng.a\033[0m"
	@ar rcs $@ $(LIB_LODEPNG_OBJS)
../objs/ext_lodepng.o: ext/lodepng.cpp
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

#libglog
$(LIB_GLOG_OBJS): $(LIB_GLOG_HEADERS)
../libs/libglog.a: $(LIB_GLOG_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core glog library libglog.a\033[0m"
	@ar rcs $@ $(LIB_GLOG_OBJS)
../objs/ext_glog_src_%.o: ext/glog/src/%.cc
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<
#libarhoseskymodel
$(LIB_SOLARMODEL_OBJS): $(LIB_SOLARMODEL_HEADERS)
../libs/libsolarmodel.a: $(LIB_SOLARMODEL_OBJS)
	@$(ECHO) "\033[32;49;1mBuilding the core ArHosekSkyModel library libsolarmodel.a\033[0m"
	@ar rcs $@ $(LIB_SOLARMODEL_OBJS)
	
../objs/ext_ArHosekSkyModel.o: ext/ArHosekSkyModel.c
	@$(ECHO) "\033[0;32;33mBuilding object $@ \033[0m"
	@$(CXX) $(CXXFLAGS)  -o $@ -c $<
#######################################################################################

clean:
	@$(ECHO) "\033[32;49;1mRemoving directories needed......done\033[0m"
	rm -f ../objs/* ../bin/* ../libs/libpbrt.a core/pbrtlex.[ch]* core/pbrtparse.[ch]*
