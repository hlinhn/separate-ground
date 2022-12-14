#include <flatten/flatten.h>
#include <flatten/helper.h>
#include <flatten/image_generation.h>
#include <pcl/io/pcd_io.h>

// void
// pickCurb(CloudT::Ptr points, double resolution, double max_val, double min_val)
// {
//   std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
//   cv::Point2d min, max;
//   min = minmax.first;
//   max = minmax.second;

//   unsigned int sizex = std::ceil((max.x - min.x) / resolution);
//   unsigned int sizey = std::ceil((max.y - min.y) / resolution);

//   cv::Mat minIm(sizex, sizey, CV_64FC1, cv::Scalar(0));
//   cv::Mat maxIm(sizex, sizey, CV_64FC1, cv::Scalar(0));
//   cv::Mat countIm(sizex, sizey, CV_8UC1, cv::Scalar(0));

//   for (const auto p : *points)
//   {
//     int x = std::ceil((p.x - min.x) / resolution);
//     int y = std::ceil((p.y - min.y) / resolution);
//     if (x < 0 || x >= sizex || y < 0 || y >= sizey)
//       continue;
//     if (countIm.at<char>(x, y) == 0)
//     {
//       minIm.at<double>(x, y) = p.z;
//       maxIm.at<double>(x, y) = p.z;
//       countIm.at<unsigned char>(x, y) = 1;
//       continue;
//     }
//     if (p.z < minIm.at<double>(x, y))
//       minIm.at<double>(x, y) = p.z;
//     if (p.z > maxIm.at<double>(x, y))
//       maxIm.at<double>(x, y) = p.z;
//   }

//   cv::Mat bindiffIm(sizex, sizey, CV_8UC1, cv::Scalar(255));

//   for (unsigned int i = 0; i < minIm.rows; i++)
//   {
//     for (unsigned int j = 0; j < minIm.cols; j++)
//     {
//       double d = maxIm.at<double>(i, j) - minIm.at<double>(i, j);
//       if (d < max_val && d > min_val)
//         bindiffIm.at<unsigned char>(i, j) = 255 - (d - min_val) / (max_val - min_val) * 255;
//     }
//   }

//   std::ostringstream out;
//   out << "bindiff_" << std::setprecision(6) << min.x << "_" << min.y << "_.png";
//   cv::imwrite(out.str().c_str(), bindiffIm);

//   /*
//   CloudT::Ptr curb(new CloudT());
//   for (const auto p: *points) {
//     int x = std::ceil((p.x - min.x) / resolution);
//     int y = std::ceil((p.y - min.y) / resolution);
//     if (x < 0 || x >= sizex || y < 0 || y >= sizey) continue;
//     double cell_diff = diffIm.at<double>(x, y);
//   */
// }

// void
// chooseDenseGround(CloudT::Ptr points, CloudT::Ptr ground, double resolution = 0.1)
// {
//   std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
//   cv::Point2d min, max;
//   min = minmax.first;
//   max = minmax.second;

//   unsigned int sizex = std::ceil((max.x - min.x) / resolution);
//   unsigned int sizey = std::ceil((max.y - min.y) / resolution);

//   cv::Mat groundIm(sizex, sizey, CV_64FC1, cv::Scalar(0));

//   for (const auto p : *ground)
//   {
//     int x = std::ceil((p.x - min.x) / resolution);
//     int y = std::ceil((p.y - min.y) / resolution);
//     if (x < 0 || x >= sizex || y < 0 || y >= sizey)
//       continue;
//     groundIm.at<int>(x, y) = p.z;
//   }

//   CloudT::Ptr dense_ground(new CloudT());
//   for (const auto p : *points)
//   {
//     int x = std::ceil((p.x - min.x) / resolution);
//     int y = std::ceil((p.y - min.y) / resolution);
//     if (x < 0 || x >= sizex || y < 0 || y >= sizey)
//       continue;
//     double z = groundIm.at<double>(x, y);
//     if (std::fabs(p.z - z) <= 0.1)
//     {
//       dense_ground->points.push_back(p);
//     }
//   }
//   pcl::io::savePCDFileBinary("dense.pcd", *dense_ground);
// }

// void
// cutCloud(CloudT::Ptr points, cv::Mat mask, double resolution, cv::Point2d min)
// {
//   // std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
//   // cv::Point2d min;
//   // min = minmax.first;

//   CloudT::Ptr white(new CloudT());
//   CloudT::Ptr black(new CloudT());

//   for (const auto p : *points)
//   {
//     if (p.z < -3.0 || p.z > 2.0)
//       continue;
//     int x = std::ceil((p.x - min.x) / resolution);
//     int y = std::ceil((p.y - min.y) / resolution);
//     if (x < 0 || x >= mask.rows || y < 0 || y >= mask.cols)
//     {
//       black->push_back(p);
//       continue;
//     }
//     if (mask.at<unsigned char>(x, y) >= 253)
//     {
//       white->push_back(p);
//     }
//     else
//     {
//       black->push_back(p);
//     }
//   }
//   pcl::io::savePCDFileBinary("white.pcd", *white);
//   pcl::io::savePCDFileBinary("black.pcd", *black);
// }

// void
// createCurbImage(CloudT::Ptr curb_image, double size, double resolution = 0.1)
// {
//   int edge_size = std::ceil(size / resolution);
//   cv::Mat countIm(edge_size, edge_size, CV_32SC1, cv::Scalar(0));
//   double min_pt = size / 2.0;
//   for (const auto p : *curb_image)
//   {
//     int x = std::ceil((p.x + min_pt) / resolution);
//     int y = std::ceil((p.y + min_pt) / resolution);
//     if (x < 0 || x >= edge_size || y < 0 || y >= edge_size)
//       continue;
//     countIm.at<int>(edge_size - x, edge_size - y) += 1;
//   }

//   for (unsigned int i = 0; i < countIm.rows; i++)
//   {
//     for (unsigned int j = 0; j < countIm.cols; j++)
//     {
//       if (countIm.at<int>(i, j) > 255)
//       {
//         countIm.at<int>(i, j) = 255;
//       }
//     }
//   }

//   std::ostringstream out_height;
//   out_height << "image_" << std::setprecision(6) << edge_size << "_" << resolution << ".png";
//   cv::imwrite(out_height.str().c_str(), countIm);
// }

// void
// voxelize(CloudT::Ptr points, double size = 0.2)
// {
//   CloudT::Ptr filtered(new CloudT());

//   /*
//   pcl::VoxelGrid<PointT> sor;
//   sor.setInputCloud(points);
//   sor.setLeafSize(size, size, size);
//   sor.filter(*filtered);
//   pcl::io::savePCDFileBinary("vox.pcd", *filtered);
//   */

//   std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
//   cv::Point2d min, max;
//   min = minmax.first;

//   std::map<std::pair<int, int>, std::vector<double>> zmap;
//   for (const auto p : *points)
//   {
//     int x = std::ceil((p.x - min.x) / size);
//     int y = std::ceil((p.y - min.y) / size);

//     std::pair<int, int> key(x, y);
//     if (zmap.find(key) == zmap.end())
//     {
//       std::vector<double> zs;
//       zs.push_back(p.z);
//       zmap[key] = zs;
//     }
//     else
//     {
//       zmap[key].push_back(p.z);
//     }
//   }

//   for (const auto& kv : zmap)
//   {
//     PointT p;
//     p.x = (kv.first.first + 0.5) * size + min.x;
//     p.y = (kv.first.second + 0.5) * size + min.y;

//     std::vector<double> copy = kv.second;
//     sort(copy.begin(), copy.end());
//     // p.z = copy[copy.size() / 2];
//     p.z = std::accumulate(copy.begin(), copy.end(), 0.0) / copy.size();

//     filtered->points.push_back(p);
//   }
//   pcl::io::savePCDFileBinary("vox.pcd", *filtered);
// }

// void
// flattenCloud(CloudT::Ptr points, double radius = 0.5)
// {
//   pcl::KdTreeFLANN<PointT> kdtree;
//   CloudT::Ptr flatPoints(new CloudT());
//   CloudT::Ptr averPoints(new CloudT());

//   for (const auto p : *points)
//   {
//     auto pflat = p;
//     pflat.z = 0;
//     flatPoints->points.push_back(pflat);
//   }
//   kdtree.setInputCloud(flatPoints);
//   int countPts = 0;
//   int empty = 0;
//   double stdevsum = 0;
//   double diffsum = 0;
//   for (const auto p : *points)
//   {
//     auto pSearch = p;
//     pSearch.z = 0;

//     /*
//     std::vector<int> pointIdxClose;
//     std::vector<float> pointDistanceClose;

//     int numpts = kdtree.radiusSearch(pSearch, radius * 2 / 5, pointIdxClose, pointDistanceClose);
//     if (numpts < 2) {
//       continue;
//     }
//     std::vector<double> zclose;
//     for (auto i: pointIdxClose)
//       zclose.push_back(points->at(i).z);

//     double zcloseMean = std::accumulate(zclose.begin(), zclose.end(), 0.0) / zclose.size();
//     if (fabs(p.z - zcloseMean) < 0.0025) {
//       averPoints->points.push_back(p);
//       empty++;
//       continue;
//     }
//     */

//     std::vector<int> pointIdx;
//     std::vector<float> pointDistance;

//     if (kdtree.radiusSearch(pSearch, radius, pointIdx, pointDistance) > 5)
//     {
//       countPts += 1;
//       std::vector<double> zvec;
//       for (auto ind : pointIdx)
//         zvec.push_back(points->at(ind).z);
//       auto pfix = p;
//       std::sort(zvec.begin(), zvec.end());
//       double farMean = std::accumulate(zvec.begin(), zvec.end(), 0.0) / zvec.size();
//       std::vector<double> diff(zvec.size());
//       std::transform(zvec.begin(), zvec.end(), diff.begin(), [farMean](double x) { return fabs(x - farMean); });
//       double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
//       double stdev = std::sqrt(sq_sum / zvec.size());
//       stdevsum += stdev;
//       diffsum += fabs(p.z - farMean);
//       /*
//       if (fabs(p.z - farMean) >= 0.1 * stdev) {
//         empty++;
//         continue;
//       }
//       */
//       // pfix.z = farMean;

//       pfix.z = zvec[zvec.size() / 2];

//       averPoints->points.push_back(pfix);
//     }
//     else
//     {
//       empty++;
//       averPoints->points.push_back(p);
//     }
//   }
//   std::cout << stdevsum / countPts << std::endl;
//   std::cout << diffsum / countPts << std::endl;
//   std::cout << empty << std::endl;
//   std::cout << countPts << std::endl;
//   pcl::io::savePCDFileBinary("aver.pcd", *averPoints);
// }

// cv::Point2f
// initPos(std::string name)
// {
//   int x_start = name.find("_");
//   int y_start = name.find("_", x_start + 1);
//   int y_end = name.find("_", y_start + 1);
//   std::string x_pos = name.substr(x_start + 1, y_start);
//   std::string y_pos = name.substr(y_start + 1, y_end);
//   cv::Point2f init(std::stof(x_pos), std::stof(y_pos));
//   return init;
// }

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
    createMeanIntensityImage(points, param.image_resolution);
  }
  if (param.create_height_image)
  {
    createHeightImage(points, param.image_resolution, param.height_im_max);
  }
  // if (argc == 3)
  // {
  //   // double resolution = std::stod(argv[2]);
  //   // pickCurb(points, 0.1);
  //   // flattenCloud(points);
  //   // voxelize(points);
  //   // chooseDenseGround(points, ground);
  //   createCurbImage(points, std::stod(argv[2]));
  //   return 0;
  //   cv::Mat mask = cv::imread(std::string(argv[2]), cv::IMREAD_GRAYSCALE);
  //   cv::Point2f min_points = initPos(std::string(argv[2]));
  //   // double minv, maxv;
  //   // cv::minMaxLoc(mask, &minv, &maxv);
  //   // std::cout << minv << " " << maxv << std::endl;
  //   cutCloud(points, mask, 0.05, min_points);
  // }
  // else
  // {
  //   /*
  //   cv::Mat mask = cv::imread(std::string(argv[2]), cv::IMREAD_GRAYSCALE);
  //   double minv, maxv;
  //   cv::minMaxLoc(mask, &minv, &maxv);
  //   std::cout << minv << " " << maxv << std::endl;
  //   cutCloud(points, mask, 0.2);
  //   */
  //   double resolution = std::stod(argv[2]);
  //   double maxval = std::stod(argv[3]);
  //   double minval = std::stod(argv[4]);
  //   pickCurb(points, resolution, maxval, minval);
  // }
  return 0;
}
