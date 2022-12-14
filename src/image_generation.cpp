#include <flatten/helper.h>
#include <flatten/image_generation.h>

cv::Mat
createMeanIntensityImage(Cloud::Ptr points, double resolution)
{
  auto minmax = cloudMinmaxPoint(points);
  auto min = minmax.first;
  auto max = minmax.second;

  unsigned int size_x = std::ceil((max.x - min.x) / resolution);
  unsigned int size_y = std::ceil((max.y - min.y) / resolution);

  cv::Mat sum_image(size_x, size_y, CV_64FC1, cv::Scalar(0));
  cv::Mat count_image(size_x, size_y, CV_16UC1, cv::Scalar(0));

  for (const auto p : *points)
  {
    int x = std::ceil((p.x - min.x) / resolution);
    int y = std::ceil((p.y - min.y) / resolution);
    if (x < 0 || x >= size_x || y < 0 || y >= size_y)
    {
      continue;
    }
    if (count_image.at<unsigned int>(x, y) == 0)
    {
      sum_image.at<double>(x, y) = p.intensity;
      count_image.at<unsigned int>(x, y) = 1;
      continue;
    }
    sum_image.at<double>(x, y) += p.intensity;
    count_image.at<unsigned int>(x, y) += 1;
  }

  for (unsigned int i = 0; i < count_image.rows; i++)
  {
    for (unsigned int j = 0; j < count_image.cols; j++)
    {
      if (count_image.at<unsigned int>(i, j) == 0)
      {
        continue;
      }
      sum_image.at<double>(i, j) = sum_image.at<double>(i, j) / count_image.at<unsigned int>(i, j);
    }
  }

  double hmin, hmax;
  cv::minMaxLoc(sum_image, &hmin, &hmax);

  sum_image -= hmin;
  sum_image /= (hmax + hmin);
  sum_image *= 255;
  return sum_image;
}

cv::Mat
createHeightImage(Cloud::Ptr points, double resolution, double max_height)
{
  auto minmax = cloudMinmaxPoint(points);
  auto min = minmax.first;
  auto max = minmax.second;

  unsigned int size_x = std::ceil((max.x - min.x) / resolution);
  unsigned int size_y = std::ceil((max.y - min.y) / resolution);

  cv::Mat min_image(size_x, size_y, CV_64FC1, cv::Scalar(0));
  cv::Mat max_image(size_x, size_y, CV_64FC1, cv::Scalar(0));
  cv::Mat count_image(size_x, size_y, CV_64FC1, cv::Scalar(0));

  for (const auto p : *points)
  {
    int x = std::ceil((p.x - min.x) / resolution);
    int y = std::ceil((p.y - min.y) / resolution);
    if (x < 0 || x >= size_x || y < 0 || y >= size_y)
    {
      continue;
    }

    if (count_image.at<double>(x, y) < 1)
    {
      min_image.at<double>(x, y) = p.z;
      max_image.at<double>(x, y) = p.z;
      count_image.at<double>(x, y) = 1;
      continue;
    }
    if (p.z < min_image.at<double>(x, y))
    {
      min_image.at<double>(x, y) = p.z;
    }
    if (p.z > max_image.at<double>(x, y))
    {
      max_image.at<double>(x, y) = p.z;
    }
  }

  cv::Mat diff_image(size_x, size_y, CV_8UC1, cv::Scalar(255));

  // double max_val = std::log(max_height + 8);
  double max_val = max_height;
  for (unsigned int i = 0; i < min_image.rows; i++)
  {
    for (unsigned int j = 0; j < min_image.cols; j++)
    {
      if (count_image.at<double>(i, j) < 1)
      {
        continue;
      }
      double d = max_image.at<double>(i, j) - min_image.at<double>(i, j);
      if (d < 0.05 || d > 0.3)
      {
        continue;
      }
      double v = d; // std::log((d + 8) / (max_height + 8)) + 1;

      if (d > max_val)
      {
        v = max_val;
      }
      v = (v - 0.05) / (max_val - 0.05) * 255;
      diff_image.at<unsigned char>(i, j) = 255 - std::floor(v);
    }
  }

  return diff_image;
}
