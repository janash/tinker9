#pragma once
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>
#include "mdi.h"

namespace tinker {

class MDIEngine
{
protected:
   static std::FILE* outputFile;
public:
   static void initialize(std::FILE* o);
   static void run_mdi();

   // temporary function for MDI testing purposes
   template <class F, class... Ts>
   static void mdiprint(const F& f, Ts&&... args){
       print(outputFile, f, args...);
   };
};

}
