cmake_minimum_required(VERSION 3.21)
project(Opulator)

set(CMAKE_CXX_STANDARD 14)

add_executable(Opulator main.cpp)

#target_link_libraries(Opulator -lconio)
target_link_libraries(Opulator -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_ttf -lSDL2_mixer)