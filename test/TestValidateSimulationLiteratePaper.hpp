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

#ifndef TESTVALIDATESIMULATION_HPP_
#define TESTVALIDATESIMULATION_HPP_
/*
 * = Validate a 256 cell simulation in parallel and serial (Figure 3) =
 *
 * This class was used to generate the results in Figure 3.
 * A script `CompareParallelResults.py` is provided to aid comparison
 * of the output results of the simulation.
 *
 * '''Note:  before compiling this code you need to alter the precision of the output.'''
 * You can do so by adding single line `setprecision(..)` to `NodeLocationWriter` in the main code base
 *
{{{
 Index: cell_based/src/population/writers/population_writers/NodeLocationWriter.cpp
===================================================================
--- cell_based/src/population/writers/population_writers/NodeLocationWriter.cpp (revision 21854)
+++ cell_based/src/population/writers/population_writers/NodeLocationWriter.cpp (working copy)
@@ -50,6 +50,7 @@
 template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
 void NodeLocationWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
 {
+    *this->mpOutStream << setprecision(20);
     for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
          node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
          ++node_iter)
}}}
 * == Use ==
 *
 * This test suite is designed to be run twice.  Each run will take roughly a minute (depending on your machine configuration).
 *
 *  {{{
 *  # in serial
 *  scons build=GccOptNative projects/Harvey2014/test/TestValidateSimulationLiteratePaper.hpp
 *  # In parallel
 *  scons build=GccOptNative_2 projects/Harvey2014/test/TestValidateSimulationLiteratePaper.hpp
 *  }}}
 *
 *  After this the positional output may be checked to machine output precision:
 *  {{{
 *  # May need to see this to $CHASTE_TEST_OUTPUT
 *  export OUTPUT=/tmp/$USER/testoutput
 *  ./projects/Harvey2014/test/CompareParallelResults.py 2 $OUTPUT/ValidateSimulation3Rand1/results_from_time_0/results.viznodes $OUTPUT/ValidateSimulation3Rand2/results_from_time_0/results.viznodes
 *  }}}
 * == Code overview ==
 *
 * The first thing to do is to include the necessary header files.
 *
 * === Include header files ===
 */
// The testing framework
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Cell-based Chaste include files
#include "NodeBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "ParallelCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

// Needed for parallel code
#include "PetscSetupAndFinalize.hpp"

/**
 *
 * == The test suite ==
 */
class TestValidateSimulation : public AbstractCellBasedTestSuite
{
public:

    void Test2dSimulation() throw (Exception)
    {
        NodesOnlyMesh<2> mesh;
        mesh.SetCalculateNodeNeighbours(false);
        mesh.SetMaximumInteractionDistance(1.6);

        std::vector<CellPtr> cells;

        /*
         * Here we call a parallel helper method which reads the cell locations
         * from a file on disk.
         *
         * The proliferative type ensures that the cells are not growing.
         */
        ParallelCellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> generator;

        generator.GenerateParallelCells("projects/Harvey2014/test/data/2DCellsCircle.dat",
                                        cells,
                                        mesh,
                                        CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);

        /*
         * Output from this simulation is to be found relative to `CHASTE_TEST_OUTPUT` which by default is
         * `/tmp/$USER/testoutput`
         * The folder is suffixed by the number of processes involved in this calculation.
         */
        std::ostringstream procs;
        procs << PetscTools::GetNumProcs();
        std::string output_directory = "ValidateSimulation3Rand" + procs.str();
        simulator.SetOutputDirectory(output_directory);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);

        simulator.AddForce(p_force);

        simulator.SetDt(1.0/240.0);
        simulator.SetSamplingTimestepMultiple(240);
        simulator.SetEndTime(100.0);
        simulator.Solve();

        /*
         * Report on the time taken to run the simulation
         */
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();
    }
};

#endif /*TESTVALIDATESIMULATION_HPP_*/
