# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/makovey/CLionProjects/mat_vector

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/makovey/CLionProjects/mat_vector/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/mat_vector.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mat_vector.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mat_vector.dir/flags.make

CMakeFiles/mat_vector.dir/main.cpp.o: CMakeFiles/mat_vector.dir/flags.make
CMakeFiles/mat_vector.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/makovey/CLionProjects/mat_vector/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/mat_vector.dir/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mat_vector.dir/main.cpp.o -c /Users/makovey/CLionProjects/mat_vector/main.cpp

CMakeFiles/mat_vector.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mat_vector.dir/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/makovey/CLionProjects/mat_vector/main.cpp > CMakeFiles/mat_vector.dir/main.cpp.i

CMakeFiles/mat_vector.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mat_vector.dir/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/makovey/CLionProjects/mat_vector/main.cpp -o CMakeFiles/mat_vector.dir/main.cpp.s

CMakeFiles/mat_vector.dir/Matrix.cpp.o: CMakeFiles/mat_vector.dir/flags.make
CMakeFiles/mat_vector.dir/Matrix.cpp.o: ../Matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/makovey/CLionProjects/mat_vector/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/mat_vector.dir/Matrix.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mat_vector.dir/Matrix.cpp.o -c /Users/makovey/CLionProjects/mat_vector/Matrix.cpp

CMakeFiles/mat_vector.dir/Matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mat_vector.dir/Matrix.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/makovey/CLionProjects/mat_vector/Matrix.cpp > CMakeFiles/mat_vector.dir/Matrix.cpp.i

CMakeFiles/mat_vector.dir/Matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mat_vector.dir/Matrix.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/makovey/CLionProjects/mat_vector/Matrix.cpp -o CMakeFiles/mat_vector.dir/Matrix.cpp.s

CMakeFiles/mat_vector.dir/Vector.cpp.o: CMakeFiles/mat_vector.dir/flags.make
CMakeFiles/mat_vector.dir/Vector.cpp.o: ../Vector.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/makovey/CLionProjects/mat_vector/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/mat_vector.dir/Vector.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mat_vector.dir/Vector.cpp.o -c /Users/makovey/CLionProjects/mat_vector/Vector.cpp

CMakeFiles/mat_vector.dir/Vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mat_vector.dir/Vector.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/makovey/CLionProjects/mat_vector/Vector.cpp > CMakeFiles/mat_vector.dir/Vector.cpp.i

CMakeFiles/mat_vector.dir/Vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mat_vector.dir/Vector.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/makovey/CLionProjects/mat_vector/Vector.cpp -o CMakeFiles/mat_vector.dir/Vector.cpp.s

# Object files for target mat_vector
mat_vector_OBJECTS = \
"CMakeFiles/mat_vector.dir/main.cpp.o" \
"CMakeFiles/mat_vector.dir/Matrix.cpp.o" \
"CMakeFiles/mat_vector.dir/Vector.cpp.o"

# External object files for target mat_vector
mat_vector_EXTERNAL_OBJECTS =

mat_vector: CMakeFiles/mat_vector.dir/main.cpp.o
mat_vector: CMakeFiles/mat_vector.dir/Matrix.cpp.o
mat_vector: CMakeFiles/mat_vector.dir/Vector.cpp.o
mat_vector: CMakeFiles/mat_vector.dir/build.make
mat_vector: CMakeFiles/mat_vector.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/makovey/CLionProjects/mat_vector/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable mat_vector"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mat_vector.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mat_vector.dir/build: mat_vector

.PHONY : CMakeFiles/mat_vector.dir/build

CMakeFiles/mat_vector.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mat_vector.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mat_vector.dir/clean

CMakeFiles/mat_vector.dir/depend:
	cd /Users/makovey/CLionProjects/mat_vector/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/makovey/CLionProjects/mat_vector /Users/makovey/CLionProjects/mat_vector /Users/makovey/CLionProjects/mat_vector/cmake-build-debug /Users/makovey/CLionProjects/mat_vector/cmake-build-debug /Users/makovey/CLionProjects/mat_vector/cmake-build-debug/CMakeFiles/mat_vector.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mat_vector.dir/depend

