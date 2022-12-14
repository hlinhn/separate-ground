#ifndef FLATTEN_HELPER_H_
#define FLATTEN_HELPER_H_

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

typedef pcl::PointXYZI Point;
typedef pcl::PointCloud<pcl::PointXYZI> Cloud;

struct ProcessMapCloudParam
{
  std::string map_name;
  bool separate_ground;
  double piece_size;
  int min_processing_size;
  bool save_separate_clouds;
  std::string save_ground_name;
  std::string save_nonground_name;
  double ground_resolution;
  double ground_max_diff;
  double ground_max_height;
  double cluster_resolution;
  double cluster_distance;
  double ground_near_min;
  double ground_near_max;
  double ground_erase_fudge;
  double min_distance_considered;
  double ground_highest_diff;

  bool create_mean_intensity;
  bool create_height_image;
  double image_resolution;
  double height_im_max;
};

struct HeightInd
{
  float height;
  int ind;
  HeightInd(float h, int i)
    : height(h)
    , ind(i)
  {
  }
};

ProcessMapCloudParam readParamFile(std::string file_name);
Cloud::Ptr readCloud(std::string name);
std::pair<Point, Point> cloudMinmaxPoint(Cloud::Ptr points);
Cloud::Ptr cutCloud(Cloud::Ptr points, Point cmin, Point cmax);
std::vector<int> slidingWindow(std::vector<HeightInd> ordered, double cluster_res, double cluster_dist);
Cloud::Ptr getGround(Cloud::Ptr points, Point cmin, Point cmax, ProcessMapCloudParam param);
Cloud::Ptr getNonground(Cloud::Ptr points, Cloud::Ptr ground, Point cmin, Point cmax, ProcessMapCloudParam param);

#endif
