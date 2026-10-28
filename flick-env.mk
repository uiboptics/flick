# Shared build environment for both root and direct module invocations.
ifeq ($(strip $(FLICK_PATH)),)
FLICK_PATH := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
endif

ifeq ($(strip $(FLICK_COMPILER)),)
ifeq ($(shell uname -s),Darwin)
FLICK_COMPILER := clang++ -O2 -std=c++20
else
FLICK_COMPILER := g++ -O2 -std=c++20
endif
endif

EIGEN_PATH ?= $(FLICK_PATH)/external/eigen

CPLUS_INCLUDE_PATH := $(FLICK_PATH):$(EIGEN_PATH):$(CPLUS_INCLUDE_PATH)
PATH := $(FLICK_PATH):$(PATH)

export FLICK_PATH FLICK_COMPILER EIGEN_PATH CPLUS_INCLUDE_PATH PATH
