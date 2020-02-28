#include <flatten/flatten.h>

std::pair<cv::Point2d, cv::Point2d> cloudMinmax(CloudT::Ptr points) {
  if (points->size() == 0) {
    std::cerr << "Empty cloud\n";
    throw std::invalid_argument("No point in cloud");
  }

  cv::Point2d min, max;
  std::pair<cv::Point2d, cv::Point2d> res;
  min = cv::Point2d(points->at(0).x, points->at(0).y);
  max = min;
  
  for (const auto p: *points) {
    if (p.x < min.x) min.x = p.x;
    if (p.y < min.y) min.y = p.y;
    if (p.x > max.x) max.x = p.x;
    if (p.y > max.y) max.y = p.y;
  }

  res.first = min;
  res.second = max;
  return res;
}

cv::Mat createMeanHeight(CloudT::Ptr points, double resolution) {
  std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
  cv::Point2d min, max;
  min = minmax.first;
  max = minmax.second;
  
  unsigned int sizex = std::ceil((max.x - min.x) / resolution);
  unsigned int sizey = std::ceil((max.y - min.y) / resolution);
  
  cv::Mat sumIm(sizex, sizey, CV_64FC1, cv::Scalar(0));
  cv::Mat countIm(sizex, sizey, CV_8UC1, cv::Scalar(0));

  for (const auto p: *points) {
    int x = std::ceil((p.x - min.x) / resolution);
    int y = std::ceil((p.y - min.y) / resolution);
    if (x < 0 || x >= sizex || y < 0 || y >= sizey) continue;
    if (countIm.at<char>(x, y) == 0) {
      sumIm.at<double>(x, y) = p.intensity;
      countIm.at<unsigned char>(x, y) = 1;
      continue;
    }
    sumIm.at<double>(x, y) += p.intensity;
    countIm.at<unsigned char>(x, y) += 1;
  }

  cv::Mat diffIm(sizex, sizey, CV_8UC1, cv::Scalar(0));
  for (unsigned int i = 0; i < countIm.rows; i++) {
    for (unsigned int j = 0; j < countIm.cols; j++) {
      if (countIm.at<unsigned char>(i, j) == 0) continue;
      sumIm.at<double>(i, j) = sumIm.at<double>(i, j) / countIm.at<unsigned char>(i, j);
    }
  }

  double hmin, hmax;
  cv::minMaxLoc(sumIm, &hmin, &hmax);

  sumIm -= hmin;
  sumIm /= hmax;
  sumIm *= 255;
  
  cv::imwrite("aver.png", sumIm);
  return sumIm;

}

cv::Mat createImage(CloudT::Ptr points, double resolution, double max_height = 0.5) {
  std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
  cv::Point2d min, max;
  min = minmax.first;
  max = minmax.second;
  
  unsigned int sizex = std::ceil((max.x - min.x) / resolution);
  unsigned int sizey = std::ceil((max.y - min.y) / resolution);
  
  cv::Mat minIm(sizex, sizey, CV_64FC1, cv::Scalar(0));
  cv::Mat maxIm(sizex, sizey, CV_64FC1, cv::Scalar(0));
  cv::Mat countIm(sizex, sizey, CV_8UC1, cv::Scalar(0));

  for (const auto p: *points) {
    int x = std::ceil((p.x - min.x) / resolution);
    int y = std::ceil((p.y - min.y) / resolution);
    if (x < 0 || x >= sizex || y < 0 || y >= sizey) continue;
    if (countIm.at<char>(x, y) == 0) {
      minIm.at<double>(x, y) = p.z;
      maxIm.at<double>(x, y) = p.z;
      countIm.at<unsigned char>(x, y) = 1;
      continue;
    }
    if (p.z < minIm.at<double>(x, y))
      minIm.at<double>(x, y) = p.z;
    if (p.z > maxIm.at<double>(x, y))
      maxIm.at<double>(x, y) = p.z;
  }

  cv::Mat diffIm(sizex, sizey, CV_8UC1, cv::Scalar(0));

  //double max_val = std::log(max_height + 8);
  double max_val = max_height;
  for (unsigned int i = 0; i < minIm.rows; i++) {
    for (unsigned int j = 0; j < minIm.cols; j++) {
      if (countIm.at<unsigned char>(i, j) == 0) continue;
      double d = maxIm.at<double>(i, j) - minIm.at<double>(i, j);
      double v = d;//std::log((d + 8) / (max_height + 8)) + 1;
      
      if (d > max_val) v = max_val;
      v = v / max_val * 255;
      diffIm.at<unsigned char>(i, j) = 255 - std::floor(v);
    }
  }

  cv::imwrite("diff.png", diffIm);
  return diffIm;
}

void cutCloud(CloudT::Ptr points, cv::Mat mask, double resolution) {
  std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
  cv::Point2d min;
  min = minmax.first;
  CloudT::Ptr white(new CloudT());
  CloudT::Ptr black(new CloudT());

  for (const auto p: *points) {
    int x = std::ceil((p.x - min.x) / resolution);
    int y = std::ceil((p.y - min.y) / resolution);
    if (x < 0 || x >= mask.rows || y < 0 || y >= mask.cols) continue;
    if (mask.at<unsigned char>(x, y) == 0) {
      black->points.push_back(p);
    } else {
      white->points.push_back(p);
    }
  }
  pcl::io::savePCDFileBinary("white.pcd", *white);
  pcl::io::savePCDFileBinary("black.pcd", *black);
}

void voxelize(CloudT::Ptr points, double size = 0.2) {
  CloudT::Ptr filtered(new CloudT());

  /*
  pcl::VoxelGrid<PointT> sor;
  sor.setInputCloud(points);
  sor.setLeafSize(size, size, size);
  sor.filter(*filtered);
  pcl::io::savePCDFileBinary("vox.pcd", *filtered);  
  */

  std::pair<cv::Point2d, cv::Point2d> minmax = cloudMinmax(points);
  cv::Point2d min, max;
  min = minmax.first;

  std::map<std::pair<int, int>, std::vector<double> > zmap;
  for (const auto p: *points) {
    int x = std::ceil((p.x - min.x) / size);
    int y = std::ceil((p.y - min.y) / size);

    std::pair<int, int> key(x, y);
    if (zmap.find(key) == zmap.end()) {
      std::vector<double> zs;
      zs.push_back(p.z);
      zmap[key] = zs;
    } else {
      zmap[key].push_back(p.z);
    }
  }

  for (const auto &kv: zmap) {
    PointT p;
    p.x = (kv.first.first + 0.5) * size + min.x;
    p.y = (kv.first.second + 0.5) * size + min.y;

    std::vector<double> copy = kv.second;
    sort(copy.begin(), copy.end());
    //p.z = copy[copy.size() / 2];
    p.z = std::accumulate(copy.begin(), copy.end(), 0.0) / copy.size();

    filtered->points.push_back(p);
  }
  pcl::io::savePCDFileBinary("vox.pcd", *filtered);  

  
}

void flattenCloud(CloudT::Ptr points, double radius = 0.5) {
  pcl::KdTreeFLANN<PointT> kdtree;
  CloudT::Ptr flatPoints(new CloudT());
  CloudT::Ptr averPoints(new CloudT());

  for (const auto p: *points) {
    auto pflat = p;
    pflat.z = 0;
    flatPoints->points.push_back(pflat);
  }
  kdtree.setInputCloud(flatPoints);
  int countPts = 0;
  int empty = 0;
  double stdevsum = 0;
  double diffsum = 0;
  for (const auto p: *points) {
    auto pSearch = p;
    pSearch.z = 0;

    /*
    std::vector<int> pointIdxClose;
    std::vector<float> pointDistanceClose;

    int numpts = kdtree.radiusSearch(pSearch, radius * 2 / 5, pointIdxClose, pointDistanceClose);
    if (numpts < 2) {
      continue;
    }
    std::vector<double> zclose;
    for (auto i: pointIdxClose)
      zclose.push_back(points->at(i).z);

    double zcloseMean = std::accumulate(zclose.begin(), zclose.end(), 0.0) / zclose.size();
    if (fabs(p.z - zcloseMean) < 0.0025) {
      averPoints->points.push_back(p);
      empty++;
      continue;
    }
    */
    
    std::vector<int> pointIdx;
    std::vector<float> pointDistance;
    
    if (kdtree.radiusSearch(pSearch, radius, pointIdx, pointDistance) > 5) {
      countPts += 1;
      std::vector<double> zvec;
      for (auto ind: pointIdx)
	zvec.push_back(points->at(ind).z);
      auto pfix = p;
      std::sort(zvec.begin(), zvec.end());
      double farMean = std::accumulate(zvec.begin(), zvec.end(), 0.0) / zvec.size();
      std::vector<double> diff(zvec.size());
      std::transform(zvec.begin(), zvec.end(), diff.begin(), [farMean](double x) { return fabs(x - farMean); });
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      double stdev = std::sqrt(sq_sum / zvec.size());
      stdevsum += stdev;
      diffsum += fabs(p.z - farMean);
      /*
      if (fabs(p.z - farMean) >= 0.1 * stdev) {
	empty++;
	continue;
      }
      */
      //pfix.z = farMean;
      
      pfix.z = zvec[zvec.size() / 2];
      
      averPoints->points.push_back(pfix);
    } else {
      empty++;
      averPoints->points.push_back(p);
    }
  }
  std::cout << stdevsum / countPts << std::endl;
  std::cout << diffsum / countPts << std::endl;
  std::cout << empty << std::endl;
  std::cout << countPts << std::endl;
  pcl::io::savePCDFileBinary("aver.pcd", *averPoints);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "USAGE: flatten filename\n";
    return 1;
  }

  CloudT::Ptr points(new CloudT());
  const bool read_success = pcl::io::loadPLYFile(std::string(argv[1]), *points) == 0;
  if(!read_success)
  {
    std::cerr << "PCL failed to load .ply file: " << argv[1] << std::endl;
    return 1;
  }

  if (argc == 2) {
    createImage(points, 0.2);
    //createMeanHeight(points, 0.05);
    //flattenCloud(points);
    //voxelize(points);
  } else {
    cv::Mat mask = cv::imread(std::string(argv[2]), cv::IMREAD_GRAYSCALE);
    double minv, maxv;
    cv::minMaxLoc(mask, &minv, &maxv);
    std::cout << minv << " " << maxv << std::endl;
    cutCloud(points, mask, 0.2);
  }
  return 0;
}
