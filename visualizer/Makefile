# Attempt to load a config.make file.
# If none is found, project defaults in config.project.make will be used.
ifneq ($(wildcard config.make),)
	include config.make
endif

OF_ROOT=/home/missinguser/Downloads/of_v0.12.0_linux64gcc6_release
# make sure the the OF_ROOT location is defined
OF_ROOT?=$(realpath ../../..)


# call the project makefile!
include $(OF_ROOT)/libs/openFrameworksCompiled/project/makefileCommon/compile.project.mk
