
# Read README.md for information about compilation and running

# Configure both root and direct module builds without requiring shell-profile
# variables. The optional install-shell-env target remains for older workflows.
include flick-env.mk

ACCURT_KERNEL_DIR ?= $(FLICK_PATH)/external/accurt-kernel
ACCURT_KERNEL_AVAILABLE := $(if $(wildcard $(ACCURT_KERNEL_DIR)/kernel.f95),1,)

export ACCURT_KERNEL_DIR ACCURT_KERNEL_AVAILABLE

MODULE_DIRS := \
	environment \
	astronomy \
	numeric/linalg \
	numeric \
	numeric/legendre \
	numeric/wigner \
	numeric/spherical_harmonics \
	mie \
	geometry \
	polarization \
	component \
	material \
	coating \
	material/gas \
	material/gas/smooth_input \
	material/snow_impurity \
	material/aerosols \
	material/water \
	material/water/refractive_index \
	material/ice \
	material/marine_cdom \
	material/marine_particles \
	accurt_api \
	$(if $(ACCURT_KERNEL_AVAILABLE),accurt_kernel_api) \
	transporter \
	radiator \
	model \
	main \
	main/commands \
	Example/single_layer_slab

TEST_DIRS := $(filter-out main Example/% material/gas/smooth_input,$(MODULE_DIRS))
CLEAN_DIRS := $(MODULE_DIRS) 

BINDIR ?= /usr/local/bin

.PHONY: all with-python build test clean test-python install install-eigen update-eigen install-shell-env install-accurt-kernel update-accurt-kernel binary-mac

all:	install-eigen build test

with-python:	all python

build:
	@set -e; for dir in $(MODULE_DIRS); do \
		$(MAKE) -C "$$dir" obj link; \
	done

test:
	@set -e; for dir in $(TEST_DIRS); do \
		$(MAKE) -C "$$dir" test; \
	done

clean:
	@set -e; for dir in $(CLEAN_DIRS); do \
		$(MAKE) -C "$$dir" clean; \
	done
	@rm -rf external/eigen
	@rm -f *~

test-python:	
	@echo ''
	@echo 'Testing all python scripts. May take a while ...'
	cd Example/python_plots; python3 test_all.py
	cd Example/accurt_calls/logo; python3 test_all.py
	cd Example/accurt_calls/atmosphere_ocean; python3 test_all.py
	cd Example/accurt_calls/snow_albedo; python3 test_all.py
	cd Example/accurt_calls/arctic_surface_radiation; python3 test_all.py
	cd Example/accurt_calls/surface_reference; python3 test_all.py
	cd Example/accurt_calls/conical_radiation; python3 test_all.py

install: build
	@mkdir -p "$(BINDIR)"
	@ln -sfn "$(CURDIR)/flick" "$(BINDIR)/flick"
	@echo "Installed Flick launcher in $(BINDIR)/flick"


EIGEN_DIR = external/eigen
EIGEN_REPO = https://gitlab.com/libeigen/eigen.git
ACCURT_KERNEL_REPO = https://git.app.uib.no/optics/model/accurt-kernel.git

install-accurt-kernel:
	@if [ ! -d "$(ACCURT_KERNEL_DIR)" ]; then \
		echo "AccuRT kernel not found. Cloning from $(ACCURT_KERNEL_REPO)..."; \
		mkdir -p "$(dir $(ACCURT_KERNEL_DIR))"; \
		git clone $(ACCURT_KERNEL_REPO) "$(ACCURT_KERNEL_DIR)"; \
	else \
		echo "AccuRT kernel already exists in $(ACCURT_KERNEL_DIR)."; \
	fi

update-accurt-kernel:
	@if [ -d "$(ACCURT_KERNEL_DIR)/.git" ]; then \
		echo "Updating AccuRT kernel in $(ACCURT_KERNEL_DIR)..."; \
		git -C "$(ACCURT_KERNEL_DIR)" pull; \
	else \
		echo "AccuRT kernel repository not found in $(ACCURT_KERNEL_DIR). Run 'make install-accurt-kernel' first."; \
		exit 1; \
	fi

# Build a relocatable macOS archive containing Flick, its run-time input data,
# Eigen, the AccuRT kernel sources, and all non-system dynamic libraries used by
# the executable.  The resulting archive is written below dist/.
binary-mac: install-eigen install-accurt-kernel build
	@./scripts/package_binary_mac.sh

install-eigen:
	@if [ ! -d "$(EIGEN_DIR)" ]; then \
		echo "Eigen not found. Cloning from $(EIGEN_REPO)..."; \
		git clone --depth 1 $(EIGEN_REPO) $(EIGEN_DIR); \
	else \
		echo "Eigen already exists in $(EIGEN_DIR)."; \
	fi

update-eigen:
	@if [ -d "$(EIGEN_DIR)/.git" ]; then \
		echo "Updating Eigen in $(EIGEN_DIR)..."; \
		git -C "$(EIGEN_DIR)" pull; \
	else \
		echo "Eigen repository not found in $(EIGEN_DIR). Run 'make install-eigen' first."; \
		exit 1; \
	fi

install-shell-env:
	@./update_shell.sh
