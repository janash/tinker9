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

MDI_Comm MDIEngine::mdi_comm;
int MDIEngine::target_node_id;
int MDIEngine::next_node_id;
int MDIEngine::default_node_id;
int MDIEngine::forces_node_id;
int MDIEngine::coords_node_id;

//void MDIEngine::initialize(std::FILE* o)
void MDIEngine::initialize(int* argc_ptr, char*** argv_ptr)
{
   next_node_id = 1;
   default_node_id = 2;
   forces_node_id = 3;
   coords_node_id = 4;
   target_node_id = next_node_id;

   mdiprint("Intializing MDI\n");

   // Initialize the MDI Library
   MDI_Init(argc_ptr, argv_ptr);
}

void MDIEngine::run_mdi(int node_id)
{
  int ret;
  int is_initialized = 0;
  ret = MDI_Initialized(&is_initialized);
  if ( ret ) {
    TINKER_THROW(format("MDI  --  Error in MDI_Initialized\n"));
  }
  if (not is_initialized) { return; }

  mdiprint("Accepting MDI communicator\n");
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
    mdiprint("Received command: %s\n",command);
    if ( ret ) {
       TINKER_THROW(format("MDI  --  Error in MDI_Recv_command\n"));
    }

    /* Respond to the received command */
    if ( strcmp(command, "<@") == 0 ) {
      if ( node_id == default_node_id ) {
          MDI_Send("@DEFAULT", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
      }
      else if ( node_id == forces_node_id ) {
          MDI_Send("@FORCES", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
      }
      else if ( node_id == coords_node_id ) {
          MDI_Send("@COORDS", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
      }
    }
    else if ( strcmp(command, "EXIT") == 0 ) {
      exit_flag = true;
    }
    else {
      TINKER_THROW(format("MDI  --  Received unsupported command: %s\n",
         command));
      exit_flag = true;
    }
  }

  // Free any memory allocations
  delete [] command;

}

}
