#include "util_mdstate.h"

TINKER_NAMESPACE_BEGIN
void sum_energy_acc_impl_(real* ebuf, int top) {
  #pragma acc serial deviceptr(ebuf,esum)
  {
    for (int i = 1; i < top; ++i)
      *esum += ebuf[i];
  }
}

void sum_virial_acc_impl_(real* vbuf, int top, int virlen) {
  #pragma acc serial deviceptr(vbuf,vir)
  {
    for (int i = 1; i < top; ++i)
      for (int j = 0; j < 9; ++j)
        vir[j] += vbuf[i * virlen + j];
  }
}
TINKER_NAMESPACE_END
