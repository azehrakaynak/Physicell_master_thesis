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
#include "./custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

bool immuno_on = false; // global flag for immunotherapy

void introduce_immune_cells()
{
    Cell_Definition* pCD = find_cell_definition("T_cell");
    for (int i = 0; i < 100; i++)
    {
        std::vector<double> pos = {
            UniformRandom() * 500 - 250,
            UniformRandom() * 500 - 250,
            0.0 };
        Cell* pC = create_cell(*pCD);
        pC->assign_position(pos);
    }
}

int main(int argc, char* argv[])
{
    bool XML_status = false;
    char copy_command[1024];
    if (argc > 1)
    {
        XML_status = load_PhysiCell_config_file(argv[1]);
        sprintf(copy_command, "cp %s %s", argv[1], PhysiCell_settings.folder.c_str());
    }
    else
    {
        XML_status = load_PhysiCell_config_file("./config/PhysiCell_settings.xml");
        sprintf(copy_command, "cp ./config/PhysiCell_settings.xml %s", PhysiCell_settings.folder.c_str());
    }
    if (!XML_status) { exit(-1); }

    system(copy_command);
    omp_set_num_threads(PhysiCell_settings.omp_num_threads);

    setup_microenvironment();
    Cell_Container* cell_container = create_cell_container_for_microenvironment(microenvironment, 30);
    create_cell_types();
    setup_tissue();

    std::ofstream tumor_volume_file("output/tumor_volume.txt");
    tumor_volume_file << "time\tsensitive_volume\tresistant_volume\ttotal_volume\tchemo\timmuno\n";

    bool chemo_on = false;
    bool immune_cells_introduced = false;

    int resistance_threshold_on = parameters.ints("resistance_threshold_on");
    int resistance_threshold_off = parameters.ints("resistance_threshold_off");
    int immune_threshold_on = parameters.ints("immune_activation_threshold_on");
    int immune_threshold_off = parameters.ints("immune_activation_threshold_off");
    int immune_activation_time= parameters.ints("immune_activation_time");

    int doxo_index = microenvironment.find_density_index("doxorubicin");

    while (PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1 * diffusion_dt)
    {
        microenvironment.simulate_diffusion_decay(diffusion_dt);
        ((Cell_Container*)microenvironment.agent_container)->update_all_cells(PhysiCell_globals.current_time);

        double sensitive_volume = 0.0;
        double resistant_volume = 0.0;
        double total_volume = 0.0;
        
        for (Cell* pC : *all_cells)
        {
            if (!pC->phenotype.death.dead)
            {
                if (pC->type == 0) // sensitive cells
                    sensitive_volume += 1.0;
                else if (pC->type == 1) // resistant cells
                    resistant_volume += 1.0;
            }
        }
        total_volume = sensitive_volume + resistant_volume;

        // Adaptive chemo logic based on total tumor volume
        if (!chemo_on && total_volume > resistance_threshold_on)
        {
            chemo_on = true;
            std::cout << "[CHEMO ON] time = " << PhysiCell_globals.current_time 
                     << " | total volume = " << total_volume 
                     << " | sensitive = " << sensitive_volume 
                     << " | resistant = " << resistant_volume << "\n";
        }
        else if (chemo_on && total_volume < resistance_threshold_off)
        {
            chemo_on = false;
            std::cout << "[CHEMO OFF] time = " << PhysiCell_globals.current_time 
                     << " | total volume = " << total_volume 
                     << " | sensitive = " << sensitive_volume 
                     << " | resistant = " << resistant_volume << "\n";
        }

        // Adaptive immunotherapy logic based on resistant cell population
        if (!immuno_on && (resistant_volume > immune_threshold_on || 
            (PhysiCell_globals.current_time >= immune_activation_time && resistant_volume > 0)))
        {
            immuno_on = true;
            std::cout << "[IMMUNO ON] time = " << PhysiCell_globals.current_time 
                     << " | resistant volume = " << resistant_volume << "\n";
        }
        else if (immuno_on && resistant_volume < immune_threshold_off)
        {
            immuno_on = false;
            std::cout << "[IMMUNO OFF] time = " << PhysiCell_globals.current_time 
                     << " | resistant volume = " << resistant_volume << "\n";
        }

        // Apply chemotherapy
        for (int i = 0; i < microenvironment.number_of_voxels(); i++)
            microenvironment(i)[doxo_index] = chemo_on ? 1.0 : 0.0;

        // Introduce immune cells when immunotherapy is activated
        if (immuno_on && !immune_cells_introduced)
        {
            introduce_immune_cells();
            immune_cells_introduced = true;
        }

        tumor_volume_file << PhysiCell_globals.current_time << "\t"
                         << sensitive_volume << "\t"
                         << resistant_volume << "\t"
                         << total_volume << "\t"
                         << (chemo_on ? 1 : 0) << "\t"
                         << (immuno_on ? 1 : 0) << "\n";

        PhysiCell_globals.current_time += diffusion_dt;
    }

    tumor_volume_file.close();
    return 0;
}

