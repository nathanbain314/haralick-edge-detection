# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home1/03726/nbain/har/edges

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home1/03726/nbain/har/edges

# Include any dependencies generated for this target.
include src/CMakeFiles/iio.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/iio.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/iio.dir/flags.make

src/CMakeFiles/iio.dir/iio.c.o: src/CMakeFiles/iio.dir/flags.make
src/CMakeFiles/iio.dir/iio.c.o: src/iio.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home1/03726/nbain/har/edges/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object src/CMakeFiles/iio.dir/iio.c.o"
	cd /home1/03726/nbain/har/edges/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/iio.dir/iio.c.o   -c /home1/03726/nbain/har/edges/src/iio.c

src/CMakeFiles/iio.dir/iio.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/iio.dir/iio.c.i"
	cd /home1/03726/nbain/har/edges/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home1/03726/nbain/har/edges/src/iio.c > CMakeFiles/iio.dir/iio.c.i

src/CMakeFiles/iio.dir/iio.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/iio.dir/iio.c.s"
	cd /home1/03726/nbain/har/edges/src && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home1/03726/nbain/har/edges/src/iio.c -o CMakeFiles/iio.dir/iio.c.s

src/CMakeFiles/iio.dir/iio.c.o.requires:
.PHONY : src/CMakeFiles/iio.dir/iio.c.o.requires

src/CMakeFiles/iio.dir/iio.c.o.provides: src/CMakeFiles/iio.dir/iio.c.o.requires
	$(MAKE) -f src/CMakeFiles/iio.dir/build.make src/CMakeFiles/iio.dir/iio.c.o.provides.build
.PHONY : src/CMakeFiles/iio.dir/iio.c.o.provides

src/CMakeFiles/iio.dir/iio.c.o.provides.build: src/CMakeFiles/iio.dir/iio.c.o

# Object files for target iio
iio_OBJECTS = \
"CMakeFiles/iio.dir/iio.c.o"

# External object files for target iio
iio_EXTERNAL_OBJECTS =

src/libiio.a: src/CMakeFiles/iio.dir/iio.c.o
src/libiio.a: src/CMakeFiles/iio.dir/build.make
src/libiio.a: src/CMakeFiles/iio.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libiio.a"
	cd /home1/03726/nbain/har/edges/src && $(CMAKE_COMMAND) -P CMakeFiles/iio.dir/cmake_clean_target.cmake
	cd /home1/03726/nbain/har/edges/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/iio.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/iio.dir/build: src/libiio.a
.PHONY : src/CMakeFiles/iio.dir/build

src/CMakeFiles/iio.dir/requires: src/CMakeFiles/iio.dir/iio.c.o.requires
.PHONY : src/CMakeFiles/iio.dir/requires

src/CMakeFiles/iio.dir/clean:
	cd /home1/03726/nbain/har/edges/src && $(CMAKE_COMMAND) -P CMakeFiles/iio.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/iio.dir/clean

src/CMakeFiles/iio.dir/depend:
	cd /home1/03726/nbain/har/edges && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home1/03726/nbain/har/edges /home1/03726/nbain/har/edges/src /home1/03726/nbain/har/edges /home1/03726/nbain/har/edges/src /home1/03726/nbain/har/edges/src/CMakeFiles/iio.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/iio.dir/depend
