cmake_minimum_required(VERSION 2.8.3)
project(flatten)

find_package(catkin REQUIRED
  COMPONENTS
  grid_map_core
  grid_map_ros
)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(yaml-cpp REQUIRED)
catkin_package(
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARIES_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(${PROJECT_NAME}_node src/flatten.cpp src/helper.cpp)
# add_executable(curb_data src/curb_from_outline.cpp)

target_link_libraries(${PROJECT_NAME}_node
  yaml-cpp
  ${catkin_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${PCL_LIBRARIES}
)

# target_link_libraries(curb_data
#   ${catkin_LIBRARIES}
#   ${OpenCV_LIBRARIES}
#   ${PCL_LIBRARIES}
# )
