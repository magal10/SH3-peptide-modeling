# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/APP/jetbrains/clion/2021.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /usr/local/APP/jetbrains/clion/2021.2/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/peptide.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/peptide.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/peptide.dir/flags.make

CMakeFiles/peptide.dir/main.cpp.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/peptide.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/main.cpp.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/main.cpp

CMakeFiles/peptide.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/main.cpp > CMakeFiles/peptide.dir/main.cpp.i

CMakeFiles/peptide.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/main.cpp -o CMakeFiles/peptide.dir/main.cpp.s

CMakeFiles/peptide.dir/Atom.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Atom.cc.o: ../Atom.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/peptide.dir/Atom.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Atom.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Atom.cc

CMakeFiles/peptide.dir/Atom.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Atom.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Atom.cc > CMakeFiles/peptide.dir/Atom.cc.i

CMakeFiles/peptide.dir/Atom.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Atom.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Atom.cc -o CMakeFiles/peptide.dir/Atom.cc.s

CMakeFiles/peptide.dir/Match.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Match.cc.o: ../Match.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/peptide.dir/Match.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Match.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Match.cc

CMakeFiles/peptide.dir/Match.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Match.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Match.cc > CMakeFiles/peptide.dir/Match.cc.i

CMakeFiles/peptide.dir/Match.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Match.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Match.cc -o CMakeFiles/peptide.dir/Match.cc.s

CMakeFiles/peptide.dir/Matrix3.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Matrix3.cc.o: ../Matrix3.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/peptide.dir/Matrix3.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Matrix3.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Matrix3.cc

CMakeFiles/peptide.dir/Matrix3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Matrix3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Matrix3.cc > CMakeFiles/peptide.dir/Matrix3.cc.i

CMakeFiles/peptide.dir/Matrix3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Matrix3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Matrix3.cc -o CMakeFiles/peptide.dir/Matrix3.cc.s

CMakeFiles/peptide.dir/numerics.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/numerics.cc.o: ../numerics.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/peptide.dir/numerics.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/numerics.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/numerics.cc

CMakeFiles/peptide.dir/numerics.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/numerics.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/numerics.cc > CMakeFiles/peptide.dir/numerics.cc.i

CMakeFiles/peptide.dir/numerics.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/numerics.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/numerics.cc -o CMakeFiles/peptide.dir/numerics.cc.s

CMakeFiles/peptide.dir/PDB.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/PDB.cc.o: ../PDB.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/peptide.dir/PDB.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/PDB.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/PDB.cc

CMakeFiles/peptide.dir/PDB.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/PDB.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/PDB.cc > CMakeFiles/peptide.dir/PDB.cc.i

CMakeFiles/peptide.dir/PDB.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/PDB.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/PDB.cc -o CMakeFiles/peptide.dir/PDB.cc.s

CMakeFiles/peptide.dir/RigidTrans3.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/RigidTrans3.cc.o: ../RigidTrans3.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/peptide.dir/RigidTrans3.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/RigidTrans3.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/RigidTrans3.cc

CMakeFiles/peptide.dir/RigidTrans3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/RigidTrans3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/RigidTrans3.cc > CMakeFiles/peptide.dir/RigidTrans3.cc.i

CMakeFiles/peptide.dir/RigidTrans3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/RigidTrans3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/RigidTrans3.cc -o CMakeFiles/peptide.dir/RigidTrans3.cc.s

CMakeFiles/peptide.dir/Rotation3.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Rotation3.cc.o: ../Rotation3.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/peptide.dir/Rotation3.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Rotation3.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Rotation3.cc

CMakeFiles/peptide.dir/Rotation3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Rotation3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Rotation3.cc > CMakeFiles/peptide.dir/Rotation3.cc.i

CMakeFiles/peptide.dir/Rotation3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Rotation3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Rotation3.cc -o CMakeFiles/peptide.dir/Rotation3.cc.s

CMakeFiles/peptide.dir/Triangle.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Triangle.cc.o: ../Triangle.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/peptide.dir/Triangle.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Triangle.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Triangle.cc

CMakeFiles/peptide.dir/Triangle.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Triangle.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Triangle.cc > CMakeFiles/peptide.dir/Triangle.cc.i

CMakeFiles/peptide.dir/Triangle.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Triangle.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Triangle.cc -o CMakeFiles/peptide.dir/Triangle.cc.s

CMakeFiles/peptide.dir/Vector3.cc.o: CMakeFiles/peptide.dir/flags.make
CMakeFiles/peptide.dir/Vector3.cc.o: ../Vector3.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/peptide.dir/Vector3.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/peptide.dir/Vector3.cc.o -c /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Vector3.cc

CMakeFiles/peptide.dir/Vector3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/peptide.dir/Vector3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Vector3.cc > CMakeFiles/peptide.dir/Vector3.cc.i

CMakeFiles/peptide.dir/Vector3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/peptide.dir/Vector3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/Vector3.cc -o CMakeFiles/peptide.dir/Vector3.cc.s

# Object files for target peptide
peptide_OBJECTS = \
"CMakeFiles/peptide.dir/main.cpp.o" \
"CMakeFiles/peptide.dir/Atom.cc.o" \
"CMakeFiles/peptide.dir/Match.cc.o" \
"CMakeFiles/peptide.dir/Matrix3.cc.o" \
"CMakeFiles/peptide.dir/numerics.cc.o" \
"CMakeFiles/peptide.dir/PDB.cc.o" \
"CMakeFiles/peptide.dir/RigidTrans3.cc.o" \
"CMakeFiles/peptide.dir/Rotation3.cc.o" \
"CMakeFiles/peptide.dir/Triangle.cc.o" \
"CMakeFiles/peptide.dir/Vector3.cc.o"

# External object files for target peptide
peptide_EXTERNAL_OBJECTS =

peptide: CMakeFiles/peptide.dir/main.cpp.o
peptide: CMakeFiles/peptide.dir/Atom.cc.o
peptide: CMakeFiles/peptide.dir/Match.cc.o
peptide: CMakeFiles/peptide.dir/Matrix3.cc.o
peptide: CMakeFiles/peptide.dir/numerics.cc.o
peptide: CMakeFiles/peptide.dir/PDB.cc.o
peptide: CMakeFiles/peptide.dir/RigidTrans3.cc.o
peptide: CMakeFiles/peptide.dir/Rotation3.cc.o
peptide: CMakeFiles/peptide.dir/Triangle.cc.o
peptide: CMakeFiles/peptide.dir/Vector3.cc.o
peptide: CMakeFiles/peptide.dir/build.make
peptide: CMakeFiles/peptide.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable peptide"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/peptide.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/peptide.dir/build: peptide
.PHONY : CMakeFiles/peptide.dir/build

CMakeFiles/peptide.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/peptide.dir/cmake_clean.cmake
.PHONY : CMakeFiles/peptide.dir/clean

CMakeFiles/peptide.dir/depend:
	cd /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug /cs/usr/magal.gafni/Desktop/3D_HACKATHON/extractClosestProtein/cmake-build-debug/CMakeFiles/peptide.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/peptide.dir/depend

