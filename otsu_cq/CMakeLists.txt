project(otsu)
find_package(OpenCV REQUIRED)
add_executable(otsu main.cpp)
target_link_libraries(otsu ${OpenCV_LIBS})
