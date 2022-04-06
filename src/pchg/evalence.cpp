#include "tool/error.h"
#include "tool/platform.h"

namespace tinker {
extern void evalence_cu(int vers);
void evalence(int vers)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      evalence_cu(vers);
   else
#endif
      TINKER_THROW("Combined valence energy term should not have been called.\n");
}
}
