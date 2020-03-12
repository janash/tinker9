#include "energy.h"
#include "io_text.h"
#include "tinker_rt.h"
#include <tinker/detail/inform.hh>


TINKER_NAMESPACE_BEGIN
void x_testgrad(int argc, char** argv)
{
   TINKER_RT(initial)();
   TINKER_RT(getxyz)();
   TINKER_RT(mechanic)();


   auto out = stdout;
   int digits = inform::digits;


   int flags = calc::xyz;
   flags += (calc::energy + calc::grad);


   rc_flag = flags;
   initialize();
   energy(rc_flag);


   std::vector<double> gdx(n), gdy(n), gdz(n);
   copy_gradient(calc::grad, gdx.data(), gdy.data(), gdz.data());


   const int len_e = 20 + digits;
   const char* fmt_e = "\n Total Potential Energy :{0:{1}.{2}f} Kcal/mole\n\n";
   print(out, fmt_e, esum, len_e, digits);

   std::string fmt_t;
   std::string fmt;
   if (digits == 8) {
      fmt_t =
         "\n  Type{0:4s}Atom{0:10s}dE/dX{0:11s}dE/dY{0:11s}dE/dZ{0:11s}Norm\n";
      fmt = "\n Anlyt{:>8d} {:16.8f}{:16.8f}{:16.8f}{:16.8f}";
   } else if (digits == 6) {
      fmt_t =
         "\n  Type{0:6s}Atom{0:11s}dE/dX{0:9s}dE/dY{0:9s}dE/dZ{0:11s}Norm\n";
      fmt = "\n Anlyt{:>10d}   {:14.6f}{:14.6f}{:14.6f}  {:14.6f}";
   } else {
      fmt_t =
         "\n  Type{0:6s}Atom{0:14s}dE/dX{0:7s}dE/dY{0:7s}dE/dZ{0:10s}Norm\n";
      fmt = "\n Anlyt{:>10d}       {:12.4f}{:12.4f}{:12.4f}  {:12.4f}";
   }
   print(out, fmt_t, "");
   auto do_print = [](int i, int n, int top_m) {
      if (n <= 2 * top_m)
         return true;
      else if (i < top_m)
         return true;
      else if (i >= n - top_m)
         return true;
      else
         return false;
   };
   int print_top_n = 15;
   for (int i = 0; i < n; ++i) {
      if (!do_print(i, n, print_top_n))
         continue;


      real x1 = gdx[i];
      real y1 = gdy[i];
      real z1 = gdz[i];
      real norm = x1 * x1 + y1 * y1 + z1 * z1;
      norm = REAL_SQRT(norm);
      print(out, fmt, i + 1, x1, y1, z1, norm);
   }


   double totnorm = 0;
   for (int i = 0; i < n; ++i) {
      totnorm += gdx[i] * gdx[i] + gdy[i] * gdy[i] + gdz[i] * gdz[i];
   }
   totnorm = std::sqrt(totnorm);
   double rms = totnorm / std::sqrt(n);
   print(out, "\n\n Total Gradient Norm and RMS Gradient per Atom :\n");
   const char* fmt3 = "\n Anlyt      {0:<30s}{1:{2}.{3}f}\n";
   const int len3 = 13 + digits;
   print(out, fmt3, "Total Gradient Norm Value", totnorm, len3, digits);
   print(out, fmt3, "RMS Gradient over All Atoms", rms, len3, digits);


   finish();
   TINKER_RT(final)();
}
TINKER_NAMESPACE_END
