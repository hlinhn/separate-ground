#include <flatten/helper.h>

#include <grid_map_core/GridMap.hpp>
#include <grid_map_core/GridMapMath.hpp>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <yaml-cpp/yaml.h>

ProcessMapCloudParam
readParamFile(std::string file_name)
{
  ProcessMapCloudParam param;
  const auto yaml = YAML::LoadFile(file_name);
  param.map_name = yaml["map_name"].as<std::string>();
  param.separate_ground = yaml["separate"]["enable"].as<bool>();
  param.piece_size = yaml["separate"]["size"].as<double>();
  param.min_processing_size = yaml["separate"]["min_size"].as<int>();
  param.save_separate_clouds = yaml["separate"]["save"].as<bool>();
  param.save_ground_name = yaml["separate"]["ground_name"].as<std::string>();
  param.save_nonground_name = yaml["separate"]["nonground_name"].as<std::string>();

  param.ground_resolution = yaml["separate"]["resolution"].as<bool>();
  param.ground_max_diff = yaml["separate"]["max_diff"].as<double>();
  param.ground_max_height = yaml["separate"]["max_height"].as<double>();
  param.cluster_resolution = yaml["separate"]["cluster_resolution"].as<double>();
  param.cluster_distance = yaml["separate"]["cluster_distance"].as<double>();
  param.ground_near_min = yaml["separate"]["near_min"].as<double>();
  param.ground_near_max = yaml["separate"]["near_max"].as<double>();
  param.ground_erase_fudge = yaml["separate"]["erase_fudge"].as<double>();
  param.min_distance_considered = yaml["separate"]["min_distance_considered"].as<double>();
  param.ground_highest_diff = yaml["separate"]["highest_diff"].as<double>();

  param.create_mean_intensity = yaml["image"]["intensity"].as<bool>();
  param.create_height_image = yaml["image"]["height"].as<bool>();
  param.image_resolution = yaml["image"]["resolution"].as<double>();
  param.height_im_max = yaml["image"]["max_height"].as<double>();
  param.mean_intensity_image_name = yaml["image"]["intensity_name"].as<std::string>();
  param.height_image_name = yaml["image"]["height_name"].as<std::string>();
  return param;
}

Cloud::Ptr
readCloud(std::string name)
{
  Cloud::Ptr points(new Cloud);

  bool read_success = false;
  std::string type = name.substr(name.size() - 3);
  if (type.compare(std::string("ply")) == 0)
  {
    read_success = pcl::io::loadPLYFile(name, *points) == 0;
  }
  else if (type.compare(std::string("pcd")) == 0)
  {
    read_success = pcl::io::loadPCDFile(name, *points) == 0;
  }

  if (!read_success)
  {
    std::cerr << "PCL failed to load file: " << name << std::endl;
    throw "Failed to load file";
  }
  return points;
}

std::pair<Point, Point>
cloudMinmaxPoint(Cloud::Ptr points)
{
  if (points->size() == 0)
  {
    std::cerr << "Empty cloud\n";
    throw std::invalid_argument("No point in cloud");
  }

  Point min, max;
  std::pair<Point, Point> res;
  min.x = points->at(0).x;
  min.y = points->at(0).y;
  max = min;

  for (const auto p : *points)
  {
    if (p.x < min.x)
    {
      min.x = p.x;
    }
    if (p.y < min.y)
    {
      min.y = p.y;
    }
    if (p.x > max.x)
    {
      max.x = p.x;
    }
    if (p.y > max.y)
    {
      max.y = p.y;
    }
  }

  res.first = min;
  res.second = max;
  return res;
}

Cloud::Ptr
cutCloud(Cloud::Ptr points, Point cmin, Point cmax)
{
  Cloud::Ptr box(new Cloud);
  pcl::ConditionalRemoval<Point> condrem;
  pcl::ConditionAnd<Point>::Ptr cond(new pcl::ConditionAnd<Point>);
  cond->addComparison(
      pcl::FieldComparison<Point>::ConstPtr(new pcl::FieldComparison<Point>("x", pcl::ComparisonOps::LE, cmax.x)));
  cond->addComparison(
      pcl::FieldComparison<Point>::ConstPtr(new pcl::FieldComparison<Point>("x", pcl::ComparisonOps::GE, cmin.x)));
  cond->addComparison(
      pcl::FieldComparison<Point>::ConstPtr(new pcl::FieldComparison<Point>("y", pcl::ComparisonOps::LE, cmax.y)));
  cond->addComparison(
      pcl::FieldComparison<Point>::ConstPtr(new pcl::FieldComparison<Point>("y", pcl::ComparisonOps::GE, cmin.y)));

  condrem.setCondition(cond);
  condrem.setInputCloud(points);
  condrem.setKeepOrganized(true);
  condrem.filter(*box);

  return box;
}

bool
operator<(HeightInd& lhs, HeightInd& rhs)
{
  return (lhs.height < rhs.height);
}

std::vector<int>
slidingWindow(std::vector<HeightInd> ordered, double cluster_res, double cluster_dist)
{
  std::vector<HeightInd> cluster_count;
  cluster_count.reserve(ordered.size());

  float last_height = ordered.front().height;
  cluster_count.push_back(HeightInd(last_height, 0));

  for (auto& h : ordered)
  {
    float dist = h.height - last_height;
    if (dist < cluster_res)
    {
      continue;
    }
    cluster_count.push_back(HeightInd(h.height, 0));
    last_height = h.height;
  }

  size_t last_start = 0;
  size_t largest_count = 0;
  size_t largest_ind = 0;

  for (auto& h : ordered)
  {
    bool first_point = true;
    for (size_t i = last_start; i < cluster_count.size(); i++)
    {
      auto& cpair = cluster_count.at(i);
      float dist = cpair.height - h.height;

      if (std::fabs(dist) <= cluster_dist)
      {
        if (first_point)
        {
          last_start = i;
          first_point = false;
        }

        cpair.ind++;

        if (cpair.ind > largest_count)
        {
          largest_count = cpair.ind;
          largest_ind = i;
        }
      }

      else if (dist > cluster_dist)
      {
        break;
      }
    }
  }

  float best_height = cluster_count.at(largest_ind).height;
  std::vector<int> chosen;
  for (auto& h : ordered)
  {
    if (std::fabs(h.height - best_height) <= cluster_dist)
      chosen.push_back(h.ind);
  }
  return chosen;
}

Cloud::Ptr
getGround(Cloud::Ptr points, Point cmin, Point cmax, ProcessMapCloudParam param)
{
  grid_map::GridMap map;
  map.setGeometry(grid_map::Length(cmax.x - cmin.x, cmax.y - cmin.y), param.ground_resolution);
  map.add("lowest");

  std::vector<std::vector<HeightInd>> cell_points;
  cell_points.resize(map.getSize().x() * map.getSize().y());
  std::vector<int> chosen;

  Point adjust;
  adjust.x = (cmax.x + cmin.x) / 2;
  adjust.y = (cmax.y + cmin.y) / 2;
  for (const auto& p : *points)
  {
    grid_map::Position pos(p.x - adjust.x, p.y - adjust.y);
    if (!map.isInside(pos))
    {
      continue;
    }
    float val = map.atPosition("lowest", pos);
    if (std::isnan(val))
    {
      map.atPosition("lowest", pos) = p.z;
      continue;
    }

    if (p.z < val)
    {
      map.atPosition("lowest", pos) = p.z;
    }
  }

  for (int i = 0; i < points->size(); i++)
  {
    auto p = points->at(i);
    grid_map::Position pos(p.x - adjust.x, p.y - adjust.y);
    if (!map.isInside(pos))
    {
      continue;
    }
    float val = map.atPosition("lowest", pos);

    if (p.z >= val + param.ground_max_diff || p.z >= param.ground_max_height)
    {
      continue;
    }

    grid_map::Index index;
    map.getIndex(pos, index);
    size_t linear_index = grid_map::getLinearIndexFromIndex(index, map.getSize());
    cell_points.at(linear_index).push_back(HeightInd(p.z, i));
  }

  for (size_t i = 0; i < cell_points.size(); i++)
  {
    auto& ordered = cell_points.at(i);
    if (ordered.empty())
    {
      continue;
    }
    if (ordered.size() == 1)
    {
      chosen.push_back(ordered.at(0).ind);
      continue;
    }
    std::sort(ordered.begin(), ordered.end());
    std::vector<int> populated = slidingWindow(ordered, param.cluster_resolution, param.cluster_distance);
    chosen.insert(chosen.end(), populated.begin(), populated.end());
  }

  Cloud::Ptr ground(new Cloud);

  pcl::ExtractIndices<Point> extract;
  extract.setInputCloud(points);
  pcl::IndicesPtr chosen_ptr(new std::vector<int>());
  *chosen_ptr = chosen;
  extract.setIndices(chosen_ptr);
  extract.filter(*ground);

  return ground;
}

Cloud::Ptr
getNonground(Cloud::Ptr points, Cloud::Ptr ground, Point cmin, Point cmax, ProcessMapCloudParam param)
{
  grid_map::GridMap map;
  map.setGeometry(grid_map::Length(cmax.x - cmin.x, cmax.y - cmin.y), param.ground_resolution);
  map.add("highest");
  map.add("erase");

  Point adjust;
  adjust.x = (cmax.x + cmin.x) / 2;
  adjust.y = (cmax.y + cmin.y) / 2;
  for (const auto& p : *ground)
  {
    grid_map::Position pos(p.x - adjust.x, p.y - adjust.y);
    if (!map.isInside(pos))
    {
      continue;
    }
    float val = map.atPosition("highest", pos);
    if (std::isnan(val))
    {
      map.atPosition("highest", pos) = p.z;
      continue;
    }

    if (p.z > val)
    {
      map.atPosition("highest", pos) = p.z;
    }
  }

  for (const auto& p : *points)
  {
    grid_map::Position pos(p.x - adjust.x, p.y - adjust.y);
    if (!map.isInside(pos))
    {
      continue;
    }
    float val = map.atPosition("highest", pos);
    if (std::isnan(val))
    {
      continue;
    }
    if (p.z > val + param.ground_near_min && p.z < val + param.ground_near_max)
    {
      float curval = map.atPosition("erase", pos);
      if (std::isnan(curval))
      {
        map.atPosition("erase", pos) = 1;
      }
      else
      {
        map.atPosition("erase", pos) = curval + param.ground_erase_fudge;
      }
    }
  }

  std::vector<int> not_chosen;
  for (int i = 0; i < points->size(); i++)
  {
    auto p = points->at(i);
    grid_map::Position pos(p.x - adjust.x, p.y - adjust.y);
    if (!map.isInside(pos))
    {
      continue;
    }
    float val = map.atPosition("highest", pos);
    if (std::isnan(val))
    {
      continue;
    }
    float inter = map.atPosition("erase", pos);
    if (!std::isnan(inter) && inter > param.min_distance_considered)
    {
      continue;
    }
    if (p.z <= val + param.ground_highest_diff)
    {
      not_chosen.push_back(i);
    }
  }

  Cloud::Ptr nonground(new Cloud);

  pcl::ExtractIndices<Point> extract;
  extract.setInputCloud(points);
  pcl::IndicesPtr chosen_ptr(new std::vector<int>());
  *chosen_ptr = not_chosen;
  extract.setIndices(chosen_ptr);
  extract.setNegative(true);
  extract.filter(*nonground);

  return nonground;
}
