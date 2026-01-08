# MultistageOptimizer Makefile (wraps CMake)
#
# Usage examples:
#   make build
#   make test
#   make clean
#
# Override variables:
#   make BUILD_DIR=build CONFIG=Release
#   make configure GENERATOR="Ninja"
#   make configure GENERATOR="Visual Studio 17 2022" ARCH=x64

CMAKE ?= cmake
CTEST ?= ctest

BUILD_DIR ?= build
CONFIG ?= Release

# Optional: set CMake generator and platform/arch (mainly for Windows VS generators).
GENERATOR ?=
ARCH ?=

# For single-config generators (Ninja/Unix Makefiles). Multi-config generators ignore this.
CMAKE_BUILD_TYPE ?= $(CONFIG)

CMAKE_GEN_ARGS :=
ifneq ($(strip $(GENERATOR)),)
  CMAKE_GEN_ARGS += -G "$(GENERATOR)"
endif
ifneq ($(strip $(ARCH)),)
  CMAKE_GEN_ARGS += -A $(ARCH)
endif

.PHONY: all configure build test clean distclean help

all: build

configure:
	$(CMAKE) -S . -B $(BUILD_DIR) $(CMAKE_GEN_ARGS) -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE)

build: configure
	$(CMAKE) --build $(BUILD_DIR) --config $(CONFIG)

test: build
	$(CTEST) --test-dir $(BUILD_DIR) -C $(CONFIG) --output-on-failure

clean:
	-$(CMAKE) --build $(BUILD_DIR) --target clean --config $(CONFIG)

distclean:
	$(CMAKE) -E rm -rf $(BUILD_DIR)

help:
	@echo "Targets: build, test, clean, distclean"
	@echo "Vars: BUILD_DIR (default: build), CONFIG (default: Release)"
	@echo "Optional vars: GENERATOR, ARCH"

