#include <Interpolation/interpolation.hh>

using namespace Interpolation;

template <class PIM>
requires cpt::isPIM<PIM>
static SingleDiscretizationInfo make_discretization_info(std::vector<double> inter,
                                                         std::vector<size_t> g_size)
{
   return SingleDiscretizationInfo(inter, g_size, PIM::tis, PIM::tis_d, PIM::tps, PIM::tps_d);
};

double foo(double x)
{
   return log(x);
}

int main()
{
   Grid1D grid(make_discretization_info<details::log_0_maps>({0., 0.2, 1.0}, {16, 16}));
   std::vector<double> fj = Discretize<std::vector<double>, double>(grid, foo, [](size_t n) {
      return std::vector<double>(n, 0.);
   });

   double y = 0.56661151;

   double tmp = grid.interpolate<double, std::vector<double>>(y, fj, []() -> double {
      return 0.;
   });

   double exact = foo(y);

   std::printf("%.16e\n", tmp - exact);
   std::printf("%.16e\n", tmp);
   std::printf("%.16e\n", exact);
}