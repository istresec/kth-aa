set(BINARY ${CMAKE_PROJECT_NAME})
#file(GLOB_RECURSE BLOSSOM5 "blossom5-v2.05.src/.cpp" "blossom5-v2.05.src/.h")
#add_executable(${BINARY}_kattis main_kattis.cpp TSP.h TSP.cpp utility.h utility.cpp ${BLOSSOM5})
#add_executable(${BINARY}_test main.cpp TSP.h TSP.cpp utility.h utility.cpp ${BLOSSOM5})
add_executable(${BINARY}_kattis main_kattis.cpp TSP.h utility.h utility.cpp blossom.h)
add_executable(${BINARY}_test main.cpp TSP.h utility.h utility.cpp blossom.h)


# For adding all recursively into one executable
#set(BINARY ${CMAKE_PROJECT_NAME})
#file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)
#set(SOURCES ${SOURCES})
#add_executable(${BINARY}_run ${SOURCES})
#add_library(${BINARY}_lib STATIC ${SOURCES})