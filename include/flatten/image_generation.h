#ifndef FLATTEN_IMAGE_GENERATION_H_
#define FLATTEN_IMAGE_GENERATION_H_

#include <opencv2/opencv.hpp>

cv::Mat createMeanIntensityImage(Cloud::Ptr points, double resolution);
cv::Mat createHeightImage(Cloud::Ptr points, double resolution, double max_height);

#endif
