# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /home/liangzi/download/clion-2018.3.4/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/liangzi/download/clion-2018.3.4/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/liangzi/CLionProjects/DanChunXingFa

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/DanChunXingFa.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DanChunXingFa.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DanChunXingFa.dir/flags.make

CMakeFiles/DanChunXingFa.dir/demo.cpp.o: CMakeFiles/DanChunXingFa.dir/flags.make
CMakeFiles/DanChunXingFa.dir/demo.cpp.o: ../demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DanChunXingFa.dir/demo.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DanChunXingFa.dir/demo.cpp.o -c /home/liangzi/CLionProjects/DanChunXingFa/demo.cpp

CMakeFiles/DanChunXingFa.dir/demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DanChunXingFa.dir/demo.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/liangzi/CLionProjects/DanChunXingFa/demo.cpp > CMakeFiles/DanChunXingFa.dir/demo.cpp.i

CMakeFiles/DanChunXingFa.dir/demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DanChunXingFa.dir/demo.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/liangzi/CLionProjects/DanChunXingFa/demo.cpp -o CMakeFiles/DanChunXingFa.dir/demo.cpp.s

# Object files for target DanChunXingFa
DanChunXingFa_OBJECTS = \
"CMakeFiles/DanChunXingFa.dir/demo.cpp.o"

# External object files for target DanChunXingFa
DanChunXingFa_EXTERNAL_OBJECTS =

DanChunXingFa: CMakeFiles/DanChunXingFa.dir/demo.cpp.o
DanChunXingFa: CMakeFiles/DanChunXingFa.dir/build.make
DanChunXingFa: CMakeFiles/DanChunXingFa.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable DanChunXingFa"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DanChunXingFa.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DanChunXingFa.dir/build: DanChunXingFa

.PHONY : CMakeFiles/DanChunXingFa.dir/build

CMakeFiles/DanChunXingFa.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DanChunXingFa.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DanChunXingFa.dir/clean

CMakeFiles/DanChunXingFa.dir/depend:
	cd /home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/liangzi/CLionProjects/DanChunXingFa /home/liangzi/CLionProjects/DanChunXingFa /home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug /home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug /home/liangzi/CLionProjects/DanChunXingFa/cmake-build-debug/CMakeFiles/DanChunXingFa.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DanChunXingFa.dir/depend

