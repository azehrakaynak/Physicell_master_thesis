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
#include <thread>
#include <chrono>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;
bool chemo_just_started = false;
bool immuno_just_started = false;
bool immuno_permanently_off = false;


#include <thread>
#include <chrono>
int counter = 0;
int main( int argc, char* argv[] )
{
	bool XML_status = false;
	char copy_command[1024];
	static int initial_tumor_cell_count = -1; // ✅ Declared here — only once

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
	system(copy_command);

	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	setup_microenvironment();

	double mechanics_voxel_size = 30;
	Cell_Container* cell_container = create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);

	std::cout << "STEP 2: Creating cell types" << std::endl;
	create_cell_types();
	std::cout << "STEP 3: Setting up tissue" << std::endl;
	setup_tissue();
	
	initial_tumor_cell_count = all_cells->size(); 



	
	initial_tumor_cell_count = all_cells->size();

	char filename[1024];
	sprintf(filename, "%s/initial", PhysiCell_settings.folder.c_str());
	save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);

	sprintf(filename, "%s/initial", PhysiCell_settings.folder.c_str());
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

	// Safety fallback
	if (diffusion_dt <= 0.0)
	{
		std::cerr << "[WARNING] diffusion_dt was 0. Setting to default value 0.01\n";
		diffusion_dt = 0.01;
	}

	std::cout << "[DEBUG] Starting simulation loop. max_time = " 
	          << PhysiCell_settings.max_time 
	          << ", diffusion_dt = " << diffusion_dt << std::endl;

	static bool chemo_on = false;
	static bool immuno_on = false;

	const double min_progression_check_time = 100.0; // Prevent false trigger at t=0

	while (PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1 * diffusion_dt)
	{
		int living_cancer_cells = 0;
		for (Cell* pCell : *all_cells)
		{
			if ((pCell->type == 0 || pCell->type == 1) && !pCell->phenotype.death.dead)
				living_cancer_cells++;
		}

		// Store initial tumor size
		if (initial_tumor_cell_count < 0)
			initial_tumor_cell_count = living_cancer_cells;

		// Define progression dynamically
		bool tumor_has_progressed = false;
		if (PhysiCell_globals.current_time > min_progression_check_time &&
			living_cancer_cells >= 1.1 * initial_tumor_cell_count)
		{
			tumor_has_progressed = true;
		}

		bool arm_b_mode = (parameters.strings("strategy") == "ArmB");

		// --- Chemo logic ---
		bool early_switch_mode = (parameters.strings("strategy") == "early_switch");

		if (!chemo_on &&
			(
				(arm_b_mode && immuno_on && tumor_has_progressed) ||
				(parameters.strings("strategy") == "ArmA" &&
				std::abs(PhysiCell_globals.current_time - parameters.ints("chemo_start")) < 0.01 * diffusion_dt) ||
				(early_switch_mode &&
				PhysiCell_globals.current_time >= parameters.ints("early_switch_time"))
			))
		{
			std::cout << "Chemo ON at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = true;
			chemo_on = true;
			chemo_just_started = true;

			// Always turn off immunotherapy if it is on
			if (immuno_on)
			{
				immuno_on = false;

				if (arm_b_mode)  // ✅ only turn off permanently in ArmB
					immuno_permanently_off = true;

				std::cout << "Immunotherapy OFF at " << PhysiCell_globals.current_time << std::endl;

				for (Cell* pCell : *all_cells)
				{
					if (pCell->type == 2)
						pCell->custom_data["oncoprotein_threshold"] = 0.5;
				}
			}

		}
		// Switch chemo off due to progression
		if (chemo_on && tumor_has_progressed)
		{
			std::cout << "Chemo OFF due to progression at " << PhysiCell_globals.current_time << std::endl;
			parameters.bools("treatment") = false;
			chemo_on = false;
		}


		// --- Immunotherapy logic ---
		if (!immuno_on && !chemo_on && tumor_has_progressed)

		{
			std::cout << "Immunotherapy ON at " << PhysiCell_globals.current_time
					<< " due to progression (tumor grew from " << initial_tumor_cell_count
					<< " to " << living_cancer_cells << ")" << std::endl;

			for (Cell* pCell : *all_cells)
			{
				if (pCell->type == 2)
					pCell->custom_data["oncoprotein_threshold"] = 0.2;
			}

			immuno_on = true;
			immuno_just_started = true;
		}

		// --- Rest of loop: Save, plot, advance ---
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

		if (PhysiCell_globals.current_time > PhysiCell_globals.next_SVG_save_time - 0.5 * diffusion_dt &&
			PhysiCell_settings.enable_SVG_saves)
		{
			sprintf(filename, "%s/snapshot%08u.svg", PhysiCell_settings.folder.c_str(), PhysiCell_globals.SVG_output_index);
			SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cancer_immune_coloring_function, paint_by_density_percentage);

			std::string marker_text;
			if (chemo_on)
				marker_text += "<text x='300' y='40' font-size='28' fill='red'>CHEMO ON</text>\n";
			if (immuno_on)
				marker_text += "<text x='300' y='80' font-size='28' fill='blue'>IMMUNO ON</text>\n";

			if (!marker_text.empty())
			{
				std::string svg_content;
				bool file_loaded = false;
				for (int attempts = 0; attempts < 10; ++attempts)
				{
					std::ifstream in_svg(filename);
					if (in_svg.is_open())
					{
						std::stringstream buffer;
						buffer << in_svg.rdbuf();
						svg_content = buffer.str();
						in_svg.close();

						if (!svg_content.empty())
						{
							file_loaded = true;
							break;
						}
					}
					std::this_thread::sleep_for(std::chrono::milliseconds(50));
				}

				if (file_loaded)
				{
					size_t pos = svg_content.rfind("</svg>");
					if (pos != std::string::npos)
					{
						svg_content.insert(pos, marker_text);
						std::ofstream out_svg(filename);
						out_svg << svg_content;
						out_svg.close();
					}
				}
				else
				{
					std::cerr << "[WARNING] Could not read SVG for annotation.\n";
				}
			}

			PhysiCell_globals.SVG_output_index++;
			PhysiCell_globals.next_SVG_save_time += PhysiCell_settings.SVG_save_interval;
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

	sprintf(filename, "%s/final", PhysiCell_settings.folder.c_str());
	save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);

	sprintf(filename, "%s/final.svg", PhysiCell_settings.folder.c_str());
	SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cancer_immune_coloring_function, paint_by_density_percentage);

	std::cout << "\nTotal simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value(std::cout, BioFVM::runtime_stopwatch_value());

	return 0;
}

