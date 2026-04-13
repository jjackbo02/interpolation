#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

namespace Interpolation
{
namespace Chebyshev
{

double get_chebyshev_lobatto_point(size_t j, size_t p)
{
   if (j > p) {
      throw std::domain_error("[get_chebyshev_lobatto_point]: j>p is not valid");
   }
   double jd = static_cast<double>(j);
   double pd = static_cast<double>(p);

   return cos(jd * M_PI / pd);
}

StandardGrid::StandardGrid(size_t p)
    : _p(p), _tj(_p + 1, 0), _betaj(_p + 1, 1), _Dij(_p + 1, vector_d(_p + 1, 0))
{
   _betaj[0]  = 0.5;
   _betaj[_p] = 0.5;
   for (size_t j = 0; j <= _p; j++) {
      _tj[j]     = get_chebyshev_lobatto_point(j, _p);
      _betaj[j] *= (j % 2 == 0 ? 1 : -1);
   }

   _Dij[0][0]   = (2 * _p * _p + 1) / 6.0;
   _Dij[_p][_p] = -((2.0 * _p * _p + 1) / 6.0);
   for (size_t j = 1; j < _p; j++) {
      _Dij[j][j] = (-_tj[j]) / (2.0 * (1.0 - _tj[j] * _tj[j]));
   }
   for (size_t j = 0; j <= _p; j++) {
      for (size_t k = 0; k < j; k++) {
         _Dij[j][k] = -(_betaj[j] / _betaj[k]) / (_tj[j] - _tj[k]);
      }
      for (size_t k = j + 1; k <= _p; k++) {
         _Dij[j][k] = -(_betaj[j] / _betaj[k]) / (_tj[j] - _tj[k]);
      }
   }
}

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
   if (t < -1 || t > 1 || (end - start) != _p) {
      throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
                              + " \\notin [-1, +1] OR view "
                                "into fj of wrong size: ["
                              + std::to_string(start) + ", " + std::to_string(end) + "]");
   }

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return fj[l + start];
      den += _betaj[l] / (t - _tj[l]);
   }

   double sum = 0;
   for (size_t i = start, j = 0; i <= end; i++, j++) {
      sum += poli_weight(t, j, den) * fj[i];
   }
   return sum;
}

double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{
   if (t < -1 || t > 1 || (end - start) != _p) {
      throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
                              + " \\notin [-1, +1] OR view "
                                "into fj of wrong size: ["
                              + std::to_string(start) + ", " + std::to_string(end) + "]");
   }

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) {
         double sum = 0;
         for (size_t i = start, j = 0; i <= end; i++, j++) {
            sum += fj[i] * _Dij[j][l];
         }
         return sum;
      }
      den += _betaj[l] / (t - _tj[l]);
   }

   double sum = 0;
   for (size_t i = start, j = 0; i <= end; i++, j++) {
      sum += poli_weight_der(t, j, den) * fj[i];
   }
   return sum;
}

double StandardGrid::poli_weight(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   if (fabs(t - _tj[j]) < 1.0e-15) return 1.0;

   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) return 0.0;

      den += _betaj[l] / (t - _tj[l]);
   }

   return _betaj[j] / den / (t - _tj[j]);
}

double StandardGrid::poli_weight(double t, size_t j, double den) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   if (fabs(t - _tj[j]) < 1.0e-15) return 1.0;

   return _betaj[j] / den / (t - _tj[j]);
}

double StandardGrid::poli_weight_der(double t, size_t j) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   double res = 0;
   for (size_t i = 0; i <= _p; i++) {
      if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
      res += _Dij[j][i] * poli_weight(t, i);
   }
   return res;
}

double StandardGrid::poli_weight_der(double t, size_t j, double den) const
{
   if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }
   double res = 0;
   for (size_t i = 0; i <= _p; i++) {
      if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
      res += _Dij[j][i] * poli_weight(t, i, den);
   }
   return res;
}

void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{
   if (end - start != _p) {
      throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply "
                                  "derivative matrix to partial vector.");
   }
   vector_d temp((end - start + 1), 0.);
   for (size_t i = 0; i <= _p; i++) {
      for (size_t j = 0, k = start; k <= end; k++, j++) {
         temp[i] += fj[k] * _Dij[j][i];
      }
   }

   for (size_t i = start; i <= end; i++) {
      fj[i] = temp[i - start];
   }
}

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
   vector_d result(_p + 1, 0.);
   for (size_t i = 0; i <= _p; i++) {
      result[i] = fnc(_tj[i]);
   }
   return result;
}

} // namespace Chebyshev
} // namespace Interpolation