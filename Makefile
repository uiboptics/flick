
# Read README.md for information on compilation and running

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

TEST_DIRS := $(filter-out main,$(MODULE_DIRS))
CLEAN_DIRS := $(MODULE_DIRS) material/gas/smooth_input

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
	@echo 'Testing all python scripts. May take an hour ...'
	cd Example/python_plots; python3 test_all.py
	cd Example/accurt_calls/logo; python3 test_all.py
	cd Example/accurt_calls/atmosphere_ocean; python3 test_all.py


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
	@if ! env | grep -q '^FLICK_PATH=' || [ -z "$${FLICK_PATH:-}" ]; then \
		./update_shell.sh; \
		false; \
	fi

