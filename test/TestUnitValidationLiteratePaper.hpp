/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef TESTUNITVALIDATION_HPP_
#define TESTUNITVALIDATION_HPP_
/*
 * = Validate a simple three cell simulation in parallel and serial (Section 3.1) =
 *
 *
 * On this wiki page we describe in detail some code that is used to compare the results of a simple three
 *  cell simulation in parallel and serial.  The results should be the same (to machine output precision)
 *  as described in Section 3.1.
 *
 * == Use ==
 *
 * Both tests in this file are designed to be run twice:
 *  {{{
 *  # in serial
 *  scons build=GccOptNative projects/Harvey2015/test/TestUnitValidationLiteratePaper.hpp
 *  # In parallel
 *  scons build=GccOptNative_2 projects/Harvey2015/test/TestUnitValidationLiteratePaper.hpp
 *  }}}
 *
 *  After this the positional output may be checked to machine output precision:
 *  {{{
 *  cd /tmp/$USER/testoutput
 *  # first test comparison
 *  diff ValidateTwoCells_1_Procs/results_from_time_0/results.viznodes ValidateTwoCells_2_Procs/results_from_time_0/results.viznodes
 *  # second test
 *  diff ValidateTwoCellsOneProc_1_Procs/results_from_time_0/results.viznodes ValidateTwoCellsOneProc_2_Procs/results_from_time_0/results.viznodes
 *  }}}
 *  The `results.viznode` file show one line per timestep with the x,y,z coordinates of each cell listed in order.
 *
 *  Note that the output is ''not'' given to machine precision (only C++ `stdio` precision) unless you amend the trunk with `setprecision(..)` in
 *  the  `NodeLocationWriter` (see wiki:PaperTutorials/Harvey2015/ValidateSimulation).
 *
 *  VTK files will contain full machine precision position information, together with process ownership.
 *  {{{
 *  paraview --data=/tmp/$USER/testoutput/ValidateTwoCells_2_Procs/results_from_time_0/results.pvd
 *  # View the cells by adding the Glyph filter and rotating the z-axis
 *  }}}
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 */

// For any extra output
#include <iostream>

// The testing framework we use
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// The main implementation used in this publication
#include "NodeBasedCellPopulation.hpp"

// Other Chaste cell_based code
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

// Required to run the simulation in parallel
#include "PetscSetupAndFinalize.hpp"

/*
 * == The test suite ==
 */
class TestUnitValidation : public AbstractCellBasedTestSuite
{
public:

    /*
     *  == First unit test ==
     *
     *  This test places the three cells at
     *  1. z=0.0
     *  1. z=0.4
     *  1. z=0.6
     *
     *  The natural strip size of the simulation is 0.5 so, when run on 2 or more process,
     *  the start configuration is
     *  1. z=0.0 on Process 0
     *  1. z=0.4 on Process 0
     *  1. z=0.6 on Process 1
     *  All other processes are assigned no cells.
     *
     *  The forces between cells are such that the cells are in repulsion
     */
    void TestTwoCellOnTwoProcesses() throw (Exception)
    {
        /*
         *  Here we make define the node positions.  When run in parallel this is done on every process.
         */
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.4));
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.6));

        /*
         * The nodes are used to construct the `NodesOnlyMesh` object with a "maximum interaction distance" of 0.5.
         * This parameter informs strip size for the parallelisation.  The bounding geometric region
         */
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 0.5);
        /*
         * The following is to confirm that when run parallel in the distributed mesh
         *  * process 0 owns the first two nodes
         *  * process 1 owns the 3rd node
         *  * other processes do not take part
         */
        std::cout << "Node/cell ownership on process "<<PetscTools::GetMyRank();
        for (unsigned i=0; i<3; i++)
        {
            std::cout<<"\t"<<mesh.rGetInitiallyOwnedNodes()[i];
        }
        std::cout<<"\n";

        /* Cells are constructed for each of the 3 nodes */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        /* A population is constructed.  The absolute movement threshold is there to alter us
         * if a cell moves further than expected (since such a cell might move over multiple strips).
         */
        NodeBasedCellPopulation<3> population(mesh, cells);
        population.SetAbsoluteMovementThreshold(0.5);

        /* Name the output folder based on how many processes we are running.
         * This folder is to be found relative to `CHASTE_TEST_OUTPUT` which by default is
         * `/tmp/$USER/testoutput`
         */
        std::ostringstream num_procs;
        num_procs << PetscTools::GetNumProcs();
        std::string output_directory = "ValidateTwoCells_" + num_procs.str() + "_Procs";

        /* A simulation is set up to run for one hour (of simulated time) in 120 discrete steps. */
        OffLatticeSimulation<3> simulator(population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(1.0);
        simulator.SetDt(1.0/120.0);

        /* Create a force law and pass it to the `OffLatticeSimulation`.  The cells have been positioned
         * such that the force will repel them from each other. */
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(0.5);
        simulator.AddForce(p_force);

        /* Run the simulation */
        simulator.Solve();
    }

    /*
     *  == Second unit test ==
     *
     *  This test places the three cells at
     *  1. z=0.0
     *  1. z=0.4
     *  1. z=0.45
     *
     *  The natural strip size of the simulation is 0.5 so, when run on any number of processes,
     *  the start configuration has
     *   * All cells on Process 0
     *   * No cells on Process 1, 2, ...
     *
     *  The forces between cells are such that the cells are in repulsion.  Cells will migrate during the
     *  simulation.
     */
    void TestTwoCellOnOneProcesses() throw (Exception)
    {
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.4));
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.45));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 0.5);
        std::cout << "Node/cell ownership on process "<<PetscTools::GetMyRank();
        for (unsigned i=0; i<3; i++)
        {
            std::cout<<"\t"<<mesh.rGetInitiallyOwnedNodes()[i];
        }
        std::cout<<"\n";

        // Cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Population
        NodeBasedCellPopulation<3> population(mesh, cells);
        population.SetAbsoluteMovementThreshold(0.5);

        // Name the output file based on how many procs we are running.
        std::ostringstream num_procs;
        num_procs << PetscTools::GetNumProcs();
        std::string output_directory = "ValidateTwoCellsOneProc_" + num_procs.str() + "_Procs";

        // Simulation
        OffLatticeSimulation<3> simulator(population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(1.0);
        simulator.SetDt(1.0/120.0);

        // Create a force law and pass it to the OffLatticeSimulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(0.5);
        simulator.AddForce(p_force);

        // Run
        simulator.Solve();
    }
};

#endif /*TESTUNITVALIDATION_HPP_*/
