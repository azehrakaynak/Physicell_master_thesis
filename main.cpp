/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	bool XML_status = false;
	char copy_command[1024];
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] );
		sprintf(copy_command, "cp %s %s", argv[1], PhysiCell_settings.folder.c_str());
	}
	else
	{
		XML_status = load_PhysiCell_config_file("./config/PhysiCell_settings.xml");
		sprintf(copy_command, "cp ./config/PhysiCell_settings.xml %s", PhysiCell_settings.folder.c_str());
	}
	if( !XML_status )
	{
		std::cout << "Error: Failed to load PhysiCell_settings.xml" << std::endl;
		exit(-1);
	}
	std::cerr << "[CONFIG] progression_threshold = " << parameters.ints("progression_threshold") << std::endl;
	// copy config file to output directory
	system(copy_command);

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// Microenvironment setup
	std::cout << "STEP 1: Loading microenvironment" << std::endl;
	setup_microenvironment();

	// Cell container setup
	double mechanics_voxel_size = 30;
	Cell_Container* cell_container = create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);

	// Cell and tissue setup
	std::cout << "STEP 2: Creating cell types" << std::endl;
	create_cell_types();
	std::cout << "STEP 3: Setting up tissue" << std::endl;
	setup_tissue();

	// Save initial state
	char filename[1024];
	sprintf(filename, "%s/initial", PhysiCell_settings.folder.c_str());
	save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);

	sprintf(filename, "%s/initial.svg", PhysiCell_settings.folder.c_str());
	SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cancer_immune_coloring_function, paint_by_density_percentage);
	sprintf(filename, "%s/legend.svg", PhysiCell_settings.folder.c_str());
	create_plot_legend(filename, cancer_immune_coloring_function);

	display_citations();

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	std::ofstream report_file;
	if (PhysiCell_settings.enable_legacy_saves)
	{
		sprintf(filename, "%s/simulation_report.txt", PhysiCell_settings.folder.c_str());
		report_file.open(filename);
		report_file << "simulated time\tnum cells\tnum division\tnum death\twall time" << std::endl;
	}

	// --- Main simulation loop ---
	static bool chemo_on = false;
	static bool immuno_on = false;

	while (PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1 * diffusion_dt)
	{
		// Therapy control logic
		int living_cancer_cells = 0;
		for (Cell* pCell : *all_cells)
		{
			if ((pCell->type == 0 || pCell->type == 1) && !pCell->phenotype.death.dead)
				living_cancer_cells++;
		}
			// Detect strategy mode from XML
		bool arm_b_mode = (parameters.strings("strategy") == "ArmB");

		// Chemotherapy activation logic
		if (!chemo_on &&
			((arm_b_mode && immuno_on && living_cancer_cells >= parameters.ints("progression_threshold")) ||  // Arm B: progression-based
			(!arm_b_mode && std::abs(PhysiCell_globals.current_time - parameters.ints("chemo_start")) < 0.01 * diffusion_dt))) // Arm A / Early: time-based
		{
			std::cout << "Chemo ON at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = true;
			chemo_on = true;
		}

		// Chemo OFF logic stays as-is
		if (chemo_on && std::abs(PhysiCell_globals.current_time - parameters.ints("chemo_end")) < 0.01 * diffusion_dt)
		{
			std::cout << "Chemo OFF at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = false;
			chemo_on = false;
		}
		if (!chemo_on && std::abs(PhysiCell_globals.current_time - parameters.ints("chemo_start")) < 0.01 * diffusion_dt)
		{
			std::cout << "Chemo ON at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = true;
			chemo_on = true;
		}
		if (chemo_on && std::abs(PhysiCell_globals.current_time - parameters.ints("chemo_end")) < 0.01 * diffusion_dt)
		{
			std::cout << "Chemo OFF at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = false;
			chemo_on = false;
		}

		if (!immuno_on &&
			(PhysiCell_globals.current_time >= parameters.ints("immunotherapy_start") ||
			living_cancer_cells >= parameters.ints("progression_threshold") ||
			PhysiCell_globals.current_time >= parameters.ints("early_switch_time")))
		{
			std::cout << "Immunotherapy ON at " << PhysiCell_globals.current_time << std::endl;
			for (Cell* pCell : *all_cells)
			{
				if (pCell->type == 2)
				{
					pCell->custom_data["oncoprotein_threshold"] = 0.8;
				}
			}
			immuno_on = true;
		}

		// Save data if it's time
		if (PhysiCell_globals.current_time > PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt)
		{
			display_simulation_status(std::cout);
			if (PhysiCell_settings.enable_legacy_saves)
				log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);

			if (PhysiCell_settings.enable_full_saves)
			{
				sprintf(filename, "%s/output%08u", PhysiCell_settings.folder.c_str(), PhysiCell_globals.full_output_index);
				save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
			}

			PhysiCell_globals.full_output_index++;
			PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
		}

		// Save SVG if needed
		if (PhysiCell_globals.current_time > PhysiCell_globals.next_SVG_save_time - 0.5 * diffusion_dt)
		{
			if (PhysiCell_settings.enable_SVG_saves)
			{
				sprintf(filename, "%s/snapshot%08u.svg", PhysiCell_settings.folder.c_str(), PhysiCell_globals.SVG_output_index);
				SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cancer_immune_coloring_function, paint_by_density_percentage);

				PhysiCell_globals.SVG_output_index++;
				PhysiCell_globals.next_SVG_save_time += PhysiCell_settings.SVG_save_interval;
			}
		}

		// Run simulation step
		microenvironment.simulate_diffusion_decay(diffusion_dt);
		((Cell_Container*)microenvironment.agent_container)->update_all_cells(PhysiCell_globals.current_time);

		PhysiCell_globals.current_time += diffusion_dt;
	}

	if (PhysiCell_settings.enable_legacy_saves)
	{
		log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
		report_file.close();
	}

	// Final save
	sprintf(filename, "%s/final", PhysiCell_settings.folder.c_str());
	save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);

	sprintf(filename, "%s/final.svg", PhysiCell_settings.folder.c_str());
	SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cancer_immune_coloring_function, paint_by_density_percentage);

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value(std::cout, BioFVM::runtime_stopwatch_value());

	return 0;
}
