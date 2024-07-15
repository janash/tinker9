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
int MDIEngine::initmd_node_id;
int MDIEngine::forces_node_id;
int MDIEngine::coords_node_id;
bool MDIEngine::terminate_node;
bool MDIEngine::exit_mdi;

//void MDIEngine::initialize(std::FILE* o)
void MDIEngine::initialize(int* argc_ptr, char*** argv_ptr)
{
   int ret;

   next_node_id = 1;
   default_node_id = 2;
   initmd_node_id = 3;
   forces_node_id = 4;
   coords_node_id = 5;
   target_node_id = next_node_id;
   terminate_node = false;
   exit_mdi = false;

   mdiprint("Intializing MDI\n");

   // Initialize the MDI Library
   ret = MDI_Init(argc_ptr, argv_ptr);
   if ( ret ) {
      TINKER_THROW(format("MDI  --  Error in MDI_Init\n"));
   }

   // Establish a connection with the driver
   mdiprint("Accepting MDI communicator\n");
   ret = MDI_Accept_communicator(&mdi_comm);
   if ( ret ) {
      TINKER_THROW(format("MDI  --  Error in MDI_Accept_communicator\n"));
   }
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

   // Check if this is a node MDI should listen from
   if ( exit_mdi ) { return; }
   if ( node_id != target_node_id and target_node_id != next_node_id ) { return; }
   terminate_node = false;

   /* Exit flag for the main MDI loop */
   bool exit_flag = false;

   /* MDI command from the driver */
   char* command = new char[MDI_COMMAND_LENGTH];

   /* Get conversion factors */
   double kcal_to_hartree;
   ret = MDI_Conversion_factor("kilocalorie_per_mol", "hartree", &kcal_to_hartree);
   if ( ret ) {
      TINKER_THROW(format("MDI  --  Error in MDI_Conversion_factor\n"));
      exit_mdi = true;
   }

   /* Main MDI loop */
   while (not terminate_node and not exit_mdi) {
      /* Receive a command from the driver */
      ret = MDI_Recv_command(command, mdi_comm);
      mdiprint("Received command: %s\n",command);
      if ( ret ) {
         TINKER_THROW(format("MDI  --  Error in MDI_Recv_command\n"));
         exit_mdi = true;
         delete [] command;
         return;
      }

      /* Respond to the received command */
      if ( strcmp(command, "<@") == 0 ) {
         if ( node_id == default_node_id ) {
            ret = MDI_Send("@DEFAULT", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
            if ( ret ) {
               TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
               exit_mdi = true;
            }
         }
         else if ( node_id == initmd_node_id ) {
            ret = MDI_Send("@INIT_MD", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
            if ( ret ) {
               TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
               exit_mdi = true;
            }
         }
         else if ( node_id == forces_node_id ) {
            ret = MDI_Send("@FORCES", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
            if ( ret ) {
               TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
               exit_mdi = true;
            }
         }
         else if ( node_id == coords_node_id ) {
            ret = MDI_Send("@COORDS", MDI_COMMAND_LENGTH, MDI_CHAR, mdi_comm);
            if ( ret ) {
               TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
               exit_mdi = true;
            }
         }
      }
      else if ( strcmp(command, "@") == 0 ) {
         terminate_node = true;
         target_node_id = next_node_id;
      }
      else if ( strcmp(command, "@INIT_MD") == 0 ) {
         terminate_node = true;
         target_node_id = initmd_node_id;
      }
      else if ( strcmp(command, "@COORDS") == 0 ) {
         terminate_node = true;
         target_node_id = coords_node_id;
      }
      else if ( strcmp(command, "@FORCES") == 0 ) {
         terminate_node = true;
         target_node_id = forces_node_id;
      }
      else if ( strcmp(command, "<NATOMS") == 0 ) {
         ret = MDI_Send(&n, 1, MDI_INT, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<COORDS") == 0 ) {
         double* coords = new double[3 * n];
         double angstrom_to_bohr;
         ret = MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Conversion_factor\n"));
            exit_mdi = true;
         }
         for (int iatom = 0; iatom < n; iatom++) {
            coords[(3 * iatom) + 0] = x[iatom] * angstrom_to_bohr;
            coords[(3 * iatom) + 1] = y[iatom] * angstrom_to_bohr;
            coords[(3 * iatom) + 2] = z[iatom] * angstrom_to_bohr;
         }
         ret = MDI_Send(coords, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
         delete [] coords;
      }
      else if ( strcmp(command, ">COORDS") == 0 ) {
         double* coords = new double[3 * n];
         double bohr_to_angstrom;
         ret = MDI_Conversion_factor("bohr", "angstrom", &bohr_to_angstrom);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Conversion_factor\n"));
            exit_mdi = true;
         }
         ret = MDI_Recv(coords, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Recv\n"));
            exit_mdi = true;
         }
         for (int iatom = 0; iatom < n; iatom++) {
            x[iatom] = coords[(3 * iatom) + 0] * bohr_to_angstrom;
            y[iatom] = coords[(3 * iatom) + 1] * bohr_to_angstrom;
            z[iatom] = coords[(3 * iatom) + 2] * bohr_to_angstrom;
         }
         delete [] coords;
      }
      else if ( strcmp(command, "<KE") == 0 ) {
         double mdi_ke = eksum * kcal_to_hartree;
         ret = MDI_Send(&mdi_ke, 1, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<PE") == 0 ) {
         double mdi_e;
         energy_prec tinker_e;
         copyEnergy(calc::v0, &tinker_e);
         mdi_e = tinker_e * kcal_to_hartree;
         ret = MDI_Send(&mdi_e, 1, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<FORCES") == 0 ) {
         double angstrom_to_bohr;
         double* forces = new double[3 * n];
         double* tinker_gx = new double[n];
         double* tinker_gy = new double[n];
         double* tinker_gz = new double[n];
         ret = MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Conversion_factor\n"));
            exit_mdi = true;
         }
         double conv_factor = - kcal_to_hartree / angstrom_to_bohr;
         energy_prec tinker_e;
         copyGradient(calc::v4, tinker_gx, tinker_gy, tinker_gz);
         for (int iatom = 0; iatom < n; iatom++) {
            forces[(3 * iatom) + 0] = tinker_gx[iatom] * conv_factor;
            forces[(3 * iatom) + 1] = tinker_gy[iatom] * conv_factor;
            forces[(3 * iatom) + 2] = tinker_gz[iatom] * conv_factor;
         }
         ret = MDI_Send(forces, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<MASSES") == 0 ) {
         ret = MDI_Send(mass, n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "EXIT") == 0 ) {
         exit_mdi = true;
      }
      else {
         TINKER_THROW(format("MDI  --  Received unsupported command: %s\n",
            command));
         exit_mdi = true;
      }
   }

   // Free any memory allocations
   delete [] command;

}

}
