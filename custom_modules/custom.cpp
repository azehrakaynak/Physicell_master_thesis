#include "physicell_mac_patch.h"
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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"
extern bool immuno_on;

void custom_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    if (pCell->type == 2) // T cell logic
    {
        if (!immuno_on) return;

        for (Cell* neighbor : pCell->cells_in_my_container())
        {
            if (neighbor->type != 2 && !neighbor->phenotype.death.dead)
            {
                if (dist(pCell->position, neighbor->position) < 20.0)
                {
                    if (UniformRandom() < 0.01 * dt)
                        neighbor->start_death(neighbor->phenotype.death.find_death_model_index("apoptosis"));
                }
            }
        }
        return;
    }

    int doxo_index = microenvironment.find_density_index("doxorubicin");
    double drug = pCell->nearest_density_vector()[doxo_index];
    double oncoprotein = get_single_signal(pCell, "custom:oncoprotein");
    
    // Handle sensitive tumor cells (type 0)
    if (pCell->type == 0 && drug > 0.1)
    {
        double death_prob = 0.001 * drug * dt;
        if (UniformRandom() < death_prob)
            pCell->start_death(pCell->phenotype.death.find_death_model_index("apoptosis"));
    }
    // Handle resistant tumor cells (type 1)
    else if (pCell->type == 1 && drug > 0.5) // Resistant cells need higher drug concentration
    {
        double death_prob = 0.0001 * drug * dt; // Lower death probability for resistant cells
        if (UniformRandom() < death_prob)
            pCell->start_death(pCell->phenotype.death.find_death_model_index("apoptosis"));
    }
}

void tumor_cell_phenotype_with_oncoprotein(Cell* pCell, Phenotype& phenotype, double dt)
{
    update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
    if (pCell->phenotype.death.dead)
    {
        pCell->functions.update_phenotype = NULL;
        return;
    }

    double cycle_rate = get_single_behavior(pCell, "cycle entry");
    cycle_rate *= get_single_signal(pCell, "custom:oncoprotein");
    set_single_behavior(pCell, "cycle entry", cycle_rate);
}

void create_cell_types(void)
{
    initialize_default_cell_definition();
    cell_defaults.functions.custom_cell_rule = custom_function;
    cell_defaults.functions.update_phenotype = NULL;

    initialize_cell_definitions_from_pugixml();
    
    // Set up sensitive tumor cells
    Cell_Definition* pCD_sensitive = find_cell_definition("tumor_cell_sensitive");
    pCD_sensitive->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein;
    
    // Set up resistant tumor cells
    Cell_Definition* pCD_resistant = find_cell_definition("tumor_cell_resistant");
    pCD_resistant->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein;

    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    display_cell_definitions(std::cout);
}

void setup_microenvironment(void)
{
    initialize_microenvironment();
}

void setup_tissue(void)
{
    load_cells_from_pugixml(); // optional, from XML file

    // Example tumor seeding (if not using load_cells)
    /*
    Cell_Definition* pCD_sensitive = find_cell_definition("tumor_cell_sensitive");
    Cell_Definition* pCD_resistant = find_cell_definition("tumor_cell_resistant");
    
    // Create sensitive tumor cells
    for (int i = 0; i < 80; i++)
    {
        std::vector<double> pos = {
            UniformRandom() * 200 - 100,
            UniformRandom() * 200 - 100,
            0.0 };
        Cell* pC = create_cell(*pCD_sensitive);
        pC->assign_position(pos);
    }
    
    // Create resistant tumor cells
    for (int i = 0; i < 20; i++)
    {
        std::vector<double> pos = {
            UniformRandom() * 200 - 100,
            UniformRandom() * 200 - 100,
            0.0 };
        Cell* pC = create_cell(*pCD_resistant);
        pC->assign_position(pos);
    }
    */
}

std::vector<std::string> heterogeneity_coloring_function(Cell* pCell)
{
    double p = get_single_signal(pCell, "custom:oncoprotein");
    static double p_min = parameters.doubles("oncoprotein_min");
    static double p_max = parameters.doubles("oncoprotein_max");

    std::vector<std::string> output(4, "black");
    if (pCell->type == 1) return output;

    if (!pCell->phenotype.death.dead)
    {
        int val = (int)round((p - p_min) / (p_max - p_min) * 255.0);
        char color[128];
        sprintf(color, "rgb(%d,%d,%d)", val, val, 255 - val);
        output[0] = color;
        output[1] = color;
        output[2] = color;
        return output;
    }

    output[0] = "rgb(250,138,38)";
    output[2] = "rgb(139,69,19)";
    return output;
}
