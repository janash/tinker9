#include "md/mdiengine.h"
#include "md/integrator.h"
#include "ff/energy.h"
#include "md/lflpiston.h"
#include "md/misc.h"
#include "md/pq.h"
#include "md/rattle.h"
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>
#include <string.h>
#include "mdi.h"

namespace tinker {

std::FILE* MDIEngine::outputFile;

void MDIEngine::initialize(std::FILE* o)
{
   outputFile = o;

   // the call to MDI_Init and all other initialization should go here
}

void MDIEngine::run_mdi()
{
  int ret;
  int is_initialized = 0;
  ret = MDI_Initialized(&is_initialized);
  if ( ret ) {
    TINKER_THROW(format("MDI  --  Error in MDI_Initialized\n"));
  }
  if (not is_initialized) { return; }

  MDI_Comm mdi_comm;
  ret = MDI_Accept_communicator(&mdi_comm);
  if ( ret ) {
    TINKER_THROW(format("MDI  --  Error in MDI_Accept_communicator\n"));
  }

  /* Exit flag for the main MDI loop */
  bool exit_flag = false;
 
  /* MDI command from the driver */
  char* command = new char[MDI_COMMAND_LENGTH];
 
  /* Main MDI loop */
  while (not exit_flag) {
    /* Receive a command from the driver */
    ret = MDI_Recv_command(command, mdi_comm);
    if ( ret ) {
       TINKER_THROW(format("MDI  --  Error in MDI_Recv_command\n"));
    }

    /* Respond to the received command */
    if ( strcmp(command, "EXIT") == 0 ) {
      exit_flag = true;
    }
    else {
      TINKER_THROW(format("MDI  --  Received unsupported command: %s\n",
         command));
    }
  }

  // Free any memory allocations
  delete [] command;

}

}
