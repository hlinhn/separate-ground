#include <flatten/flatten.h>
#include <flatten/helper.h>
#include <flatten/image_generation.h>
#include <pcl/io/pcd_io.h>

void
separateGround(Cloud::Ptr points, ProcessMapCloudParam param)
{
  auto minmax = cloudMinmaxPoint(points);
  Point min, max;
  min = minmax.first;
  max = minmax.second;

  auto size_x = max.x - min.x;
  auto size_y = max.y - min.y;

  int num_x = std::max(1, static_cast<int>(std::floor(size_x / param.piece_size)));
  int num_y = std::max(1, static_cast<int>(std::floor(size_y / param.piece_size)));

  Cloud::Ptr all_ground(new Cloud);
  Cloud::Ptr all_nonground(new Cloud);

  for (int i = 0; i < num_x; i++)
  {
    for (int j = 0; j < num_y; j++)
    {
      Point cmax, cmin;
      cmin.x = min.x + i * param.piece_size;
      cmin.y = min.y + j * param.piece_size;
      cmax.x = cmin.x + param.piece_size;
      cmax.y = cmin.y + param.piece_size;
      if (i == num_x - 1)
      {
        cmax.x = max.x;
      }
      if (j == num_y - 1)
      {
        cmax.y = max.y;
      }

      auto cut_cloud = cutCloud(points, cmin, cmax);
      if (cut_cloud->size() < param.min_processing_size)
      {
        continue;
      }
      auto current_ground = getGround(cut_cloud, cmin, cmax, param);
      auto current_nonground = getNonground(cut_cloud, current_ground, cmin, cmax, param);
      *all_ground += *current_ground;
      *all_nonground += *current_nonground;
    }
  }

  if (param.save_separate_clouds)
  {
    pcl::io::savePCDFileBinary(param.save_ground_name, *all_ground);
    pcl::io::savePCDFileBinary(param.save_nonground_name, *all_nonground);
  }
}

int
main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cout << "USAGE: flatten config_file\n";
    return 1;
  }

  auto param = readParamFile(std::string(argv[1]));
  Cloud::Ptr points = readCloud(param.map_name);

  if (param.separate_ground)
  {
    separateGround(points, param);
  }
  if (param.create_mean_intensity)
  {
    auto im = createMeanIntensityImage(points, param.image_resolution);
    cv::imwrite(param.mean_intensity_image_name, im);
  }
  if (param.create_height_image)
  {
    auto im = createHeightImage(points, param.image_resolution, param.height_im_max);
    cv::imwrite(param.height_image_name, im);
  }
  return 0;
}
