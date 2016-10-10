#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H

#include <vector>

#include <CGAL/Data_classification/Color.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Segmentation attribute based on colorimetric information.

    If the input point cloud comes with colorimetric information, it
    can be useful for classification purposes. This attribute computes
    a distance between the color of a point and a user-defined color
    region in the HSV color-space, defined by a Gaussian
    distribution with user-defined mean and standard deviation values.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointPMap Property map to access the colors of the input points.
  */
template <typename Kernel, typename RandomAccessIterator, typename ColorPMap>
class Attribute_color : public Attribute
{
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  
public:
  /*!
    \brief Constructs an attribute based on the given color.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the colors of the input points
    \param weight Weight of the attribute
    \param mean_h Mean hue of the selected color
    \param mean_s Mean saturation of the selected color
    \param mean_v Mean value of the selected color
    \param sd_h Standard deviation of the hue of the selected color
    \param sd_s Standard deviation of the saturation of the selected color
    \param sd_v Standard deviation of the value of the selected color

    \note The default values describe a gray color region
    corresponding to the color of a concrete road.
  */
  Attribute_color (RandomAccessIterator begin,
                   RandomAccessIterator end,
                   ColorPMap color_pmap,
                   double weight = 1.,
                   double mean_h = 156., double mean_s = 5., double mean_v = 76.,
                   double sd_h = 70., double sd_s = 12., double sd_v = 8.4)
  {
    this->weight = weight;
    for(std::size_t i = 0; i < (std::size_t)(end - begin);i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (get(color_pmap, begin[i]));
        color_attribute.push_back (std::exp (-(c[0] - mean_h) * (c[0] - mean_h) / (2. * sd_h * sd_h))
                                   * std::exp (-(c[1] - mean_s) * (c[1] - mean_s) / (2. * sd_s * sd_s))
                                   * std::exp (-(c[2] - mean_v) * (c[2] - mean_v) / (2. * sd_v * sd_v)));
      }
    this->compute_mean_max (color_attribute, this->mean, this->max);
  }

  /// \cond SKIP_IN_MANUAL
  virtual double value (std::size_t pt_index)
  {
    return color_attribute[pt_index];
  }

  virtual std::string id() { return "color"; }
  /// \endcond
};


  /// \cond SKIP_IN_MANUAL
template <typename Kernel, typename RandomAccessIterator, typename ColorPMap>
class Attribute_hsv : public Attribute
{
  typedef typename Data_classification::RGB_Color RGB_Color;
  typedef typename Data_classification::HSV_Color HSV_Color;
  
  std::vector<double> color_attribute;
  std::string m_id;
  
public:
  
  Attribute_hsv (RandomAccessIterator begin,
                 RandomAccessIterator end,
                 ColorPMap color_pmap,
                 std::size_t channel,
                 double mean, double sd,
                 double weight = 1.)
  {
    this->weight = weight;
    for(std::size_t i = 0; i < (std::size_t)(end - begin);i++)
      {
        HSV_Color c = Data_classification::rgb_to_hsv (get(color_pmap, begin[i]));
        color_attribute.push_back (std::exp (-(c[channel] - mean) * (c[channel] - mean) / (2. * sd * sd)));
      }
    this->compute_mean_max (color_attribute, this->mean, this->max);

    std::ostringstream oss;
    if (channel == 0) oss << "hue";
    else if (channel == 1) oss << "saturation";
    else if (channel == 2) oss << "value";
    oss << "_" << mean;
    m_id = oss.str();
  }

  virtual double value (std::size_t pt_index)
  {
    return color_attribute[pt_index];
  }

  virtual std::string id() { return m_id; }

};
/// \endcond

} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTE_COLOR_H
