
# Read README.md for information about compilation and running

# Keep the first build usable before the shell profile has been reloaded.
# update_shell.sh persists these values for future shells, while these exports
# make them available to recursive makes and tests in the current invocation.
FLICK_ENV_MISSING := $(if $(strip $(FLICK_PATH)),,1)

ifeq ($(FLICK_ENV_MISSING),1)
FLICK_PATH := $(CURDIR)
endif

ifeq ($(strip $(FLICK_COMPILER)),)
ifeq ($(shell uname -s),Darwin)
FLICK_COMPILER := clang++ -std=c++20
else
FLICK_COMPILER := g++ -std=c++20
endif
endif

ifeq ($(strip $(EIGEN_PATH)),)
EIGEN_PATH := $(FLICK_PATH)/external/eigen
endif

CPLUS_INCLUDE_PATH := $(FLICK_PATH):$(EIGEN_PATH):$(CPLUS_INCLUDE_PATH)
PATH := $(FLICK_PATH)/main:$(PATH)

export FLICK_PATH FLICK_COMPILER EIGEN_PATH CPLUS_INCLUDE_PATH PATH

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
	transporter \
	radiator \
	model \
	main \
	main/commands \
	Example/single_layer_slab

TEST_DIRS := $(filter-out main Example/% material/gas/smooth_input,$(MODULE_DIRS))
CLEAN_DIRS := $(MODULE_DIRS) 

.PHONY: all with-python build test clean python check-eigen check-env

all:	check-env check-eigen build test

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

python:	
	@echo ''
	@echo 'Testing all python scripts. May take a while ...'
	cd Example/python_plots; python3 test_all.py
	cd Example/accurt_calls/logo; python3 test_all.py
	cd Example/accurt_calls/atmosphere_ocean; python3 test_all.py
	cd Example/accurt_calls/snow_albedo; python3 test_all.py


EIGEN_DIR = external/eigen
EIGEN_REPO = https://gitlab.com/libeigen/eigen.git

check-eigen:
	@if [ ! -d "$(EIGEN_DIR)" ]; then \
		echo "Eigen not found. Cloning from $(EIGEN_REPO)..."; \
		git clone --depth 1 $(EIGEN_REPO) $(EIGEN_DIR); \
	else \
		echo "Eigen already exists in $(EIGEN_DIR)."; \
	fi

check-env:
	@if [ "$(FLICK_ENV_MISSING)" = "1" ]; then \
		./update_shell.sh; \
	fi
