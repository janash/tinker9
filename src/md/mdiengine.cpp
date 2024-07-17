#include "md/mdiengine.h"
#include "md/integrator.h"
#include "ff/energy.h"
#include "md/lflpiston.h"
#include "md/misc.h"
#include "md/pq.h"
#include "md/rattle.h"
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/atoms.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>
#include <string.h>
#include "mdi.h"

namespace tinker {

MDI_Comm MDIEngine::mdi_comm;
int MDIEngine::nrespa_mdi;
int MDIEngine::target_node_id;
int MDIEngine::next_node_id;
int MDIEngine::default_node_id;
int MDIEngine::initmd_node_id;
int MDIEngine::forces_node_id;
int MDIEngine::coords_node_id;
bool MDIEngine::terminate_node;
bool MDIEngine::exit_mdi;

/// Converts a floating-point value \c f to fixed-point value on host.
template <class F>
inline fixed MDIEngine::toFixedPoint(F f)
{
    static_assert(std::is_same<F, float>::value or std::is_same<F, double>::value, "");
    return static_cast<fixed>(static_cast<long long>(f * 0x100000000ull));
}

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

void MDIEngine::set_nrespa(int nrespa_in)
{
   nrespa_mdi = nrespa_in;
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
   double angstrom_to_bohr;
   ret = MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
   if ( ret ) {
      TINKER_THROW(format("MDI  --  Error in MDI_Conversion_factor\n"));
      exit_mdi = true;
   }
   double time_conv;
   ret = MDI_Conversion_factor("picosecond", "atomic_unit_of_time", &time_conv);
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
         darray::copyout(g::q0, n, atoms::x, xpos);
         darray::copyout(g::q0, n, atoms::y, ypos);
         darray::copyout(g::q0, n, atoms::z, zpos);
         waitFor(g::q0);
         for (int iatom = 0; iatom < n; iatom++) {
            coords[(3 * iatom) + 0] = atoms::x[iatom] * angstrom_to_bohr;
            coords[(3 * iatom) + 1] = atoms::y[iatom] * angstrom_to_bohr;
            coords[(3 * iatom) + 2] = atoms::z[iatom] * angstrom_to_bohr;
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
         double conv_factor = 1.0 / angstrom_to_bohr;
         ret = MDI_Recv(coords, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Recv\n"));
            exit_mdi = true;
         }/*
         for (int iatom = 0; iatom < n; iatom++) {
            x[iatom] = coords[(3 * iatom) + 0] * bohr_to_angstrom;
            y[iatom] = coords[(3 * iatom) + 1] * bohr_to_angstrom;
            z[iatom] = coords[(3 * iatom) + 2] * bohr_to_angstrom;
         }*/

         std::vector<pos_prec> x_mdi(n), y_mdi(n), z_mdi(n);
         for (int iatom = 0; iatom < n; ++iatom) {
            x_mdi[iatom] = coords[(3 * iatom) + 0] * conv_factor;
            y_mdi[iatom] = coords[(3 * iatom) + 1] * conv_factor;
            z_mdi[iatom] = coords[(3 * iatom) + 2] * conv_factor;
         }
         darray::copyin(g::q0, n, xpos, x_mdi.data());
         darray::copyin(g::q0, n, ypos, y_mdi.data());
         darray::copyin(g::q0, n, zpos, z_mdi.data());
         waitFor(g::q0);

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
      else if ( strcmp(command, "<ENERGY") == 0 ) {
         double mdi_e;
         energy_prec tinker_e;
         copyEnergy(calc::v0, &tinker_e);
         tinker_e += eksum;
         mdi_e = tinker_e * kcal_to_hartree;
         ret = MDI_Send(&mdi_e, 1, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<FORCES") == 0 ) {
         double* forces = new double[3 * n];
         double conv_factor = - kcal_to_hartree / angstrom_to_bohr;
         if ( nrespa_mdi == 1 ) {
            double* tinker_gx = new double[n];
            double* tinker_gy = new double[n];
            double* tinker_gz = new double[n];
            copyGradient(calc::v4, tinker_gx, tinker_gy, tinker_gz);
            for (int iatom = 0; iatom < n; iatom++) {
               forces[(3 * iatom) + 0] = tinker_gx[iatom] * conv_factor;
               forces[(3 * iatom) + 1] = tinker_gy[iatom] * conv_factor;
               forces[(3 * iatom) + 2] = tinker_gz[iatom] * conv_factor;
            }
            delete [] tinker_gx;
            delete [] tinker_gy;
            delete [] tinker_gz;
         }
         else {
            for (int iatom = 0; iatom < n; iatom++) {
               forces[(3 * iatom) + 0] = (gx1[iatom] + gx2[iatom]) * conv_factor;
               forces[(3 * iatom) + 1] = (gy1[iatom] + gy2[iatom]) * conv_factor;
               forces[(3 * iatom) + 2] = (gz1[iatom] + gz2[iatom]) * conv_factor;
            }
         }
         ret = MDI_Send(forces, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
         delete [] forces;
      }
      else if ( strcmp(command, ">FORCES") == 0 ) {
         double* forces = new double[3 * n];
         double conv_factor = - angstrom_to_bohr / kcal_to_hartree;
         energy_prec tinker_e;
         ret = MDI_Recv(forces, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
         std::vector<grad_prec> gx_mdi(n), gy_mdi(n), gz_mdi(n);
#if TINKER_DETERMINISTIC_FORCE
         for (int iatom = 0; iatom < n; iatom++) {
            gx_mdi[iatom] = toFixedPoint<double>( forces[(3 * iatom) + 0] * conv_factor );
            gy_mdi[iatom] = toFixedPoint<double>( forces[(3 * iatom) + 1] * conv_factor );
            gz_mdi[iatom] = toFixedPoint<double>( forces[(3 * iatom) + 2] * conv_factor );
         }
#else
         for (int iatom = 0; iatom < n; iatom++) {
            gx_mdi[iatom] = forces[(3 * iatom) + 0] * conv_factor;
            gy_mdi[iatom] = forces[(3 * iatom) + 1] * conv_factor;
            gz_mdi[iatom] = forces[(3 * iatom) + 2] * conv_factor;
         }
#endif
         darray::copyin(g::q0, n, gx, gx_mdi.data());
         darray::copyin(g::q0, n, gy, gy_mdi.data());
         darray::copyin(g::q0, n, gz, gz_mdi.data());
         waitFor(g::q0);
         // Only adjust the RESPA forces if using RESPA
         if ( nrespa_mdi != 1 ) {
            for (int iatom = 0; iatom < n; iatom++) {
               // fast RESPA forces
               gx1[iatom] = forces[(3 * iatom) + 0] * conv_factor;
               gy1[iatom] = forces[(3 * iatom) + 1] * conv_factor;
               gz1[iatom] = forces[(3 * iatom) + 2] * conv_factor;

               // slow RESPA forces
               gx2[iatom] = 0.0;
               gy2[iatom] = 0.0;
               gz2[iatom] = 0.0;
            }
         }
         delete [] forces;
      }
      else if ( strcmp(command, "<MASSES") == 0 ) {
         ret = MDI_Send(mass, n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
      }
      else if ( strcmp(command, "<VELOCITIES") == 0 ) {
         std::vector<vel_prec> vvx(n), vvy(n), vvz(n);
         darray::copyout(g::q0, n, vvx.data(), vx);
         darray::copyout(g::q0, n, vvy.data(), vy);
         darray::copyout(g::q0, n, vvz.data(), vz);
         waitFor(g::q0);

         double* velocities = new double[3 * n];
         double conv_factor = angstrom_to_bohr / time_conv;
         for (int iatom = 0; iatom < n; iatom++) {
            velocities[(3 * iatom) + 0] = vvx[iatom] * conv_factor;
            velocities[(3 * iatom) + 1] = vvy[iatom] * conv_factor;
            velocities[(3 * iatom) + 2] = vvz[iatom] * conv_factor;
         }
         ret = MDI_Send(velocities, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Send\n"));
            exit_mdi = true;
         }
         delete [] velocities;
      }
      else if ( strcmp(command, ">VELOCITIES") == 0 ) {
         double* velocities = new double[3 * n];
         double conv_factor = time_conv / angstrom_to_bohr;
         ret = MDI_Recv(velocities, 3 * n, MDI_DOUBLE, mdi_comm);
         if ( ret ) {
            TINKER_THROW(format("MDI  --  Error in MDI_Recv\n"));
            exit_mdi = true;
         }

         std::vector<vel_prec> vvx(n), vvy(n), vvz(n);
         for (int iatom = 0; iatom < n; ++iatom) {
            vvx[iatom] = velocities[(3 * iatom) + 0] * conv_factor;
            vvy[iatom] = velocities[(3 * iatom) + 1] * conv_factor;
            vvz[iatom] = velocities[(3 * iatom) + 2] * conv_factor;
         }
         darray::copyin(g::q0, n, vx, vvx.data());
         darray::copyin(g::q0, n, vy, vvy.data());
         darray::copyin(g::q0, n, vz, vvz.data());
         waitFor(g::q0);

         delete [] velocities;
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
