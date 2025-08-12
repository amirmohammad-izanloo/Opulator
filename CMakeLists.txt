cmake_minimum_required(VERSION 3.21)
project(Opulator)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Opulator main.cpp save_load.cpp)

#target_link_libraries(Opulator -lconio)
target_link_libraries(Opulator -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_ttf -lSDL2_mixer)
target_link_libraries(Opulator comdlg32)

include_directories(${CMAKE_SOURCE_DIR}/cereal/include)
