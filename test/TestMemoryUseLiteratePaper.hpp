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

#ifndef TESTMEMORYUSE_HPP_
#define TESTMEMORYUSE_HPP_
/*
 * = Measure memory use per process (Figure 4) =
 * This class was used to produce the results in Figure 4.
 * It constructs a large cell population and measures the amount of memory
 * in use. By running on larger numbers of processes, a proportionally
 * smaller amount of memory is used by each process.
 *
 * The geometry for the construction of the population is contained in
 * {{{
 * projects/Harvey2014/test/data/1024000_2d_cells.dat
 * }}}
 *
 * == Use ==
 *
 * This test suite should be run in parallel.  It should be run on several numbers of processes (from 1 to 32 processes in figure
 * in the paper.)
 *
 * A useful for loop (in `bash`) would be
 * {{{
 * for i in {1..32}; do echo $i "processes ===";scons build=GccOptNative_$i projects/Harvey2014/test/TestMemoryUseLiteratePaper.hpp | grep memory; done
 * }}}
 *
 * This runs the test on increasingly larger number of processes and only outputs the lines which state
 * the individual memory use.  Some post-processing will be needed in order to select the maximum memory use to
 * plot on the vertical axis.
 */

// The testing framework
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes for printing out memory
#include <unistd.h>
#include <sys/resource.h>

// Cell-based includes
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ParallelCellsGenerator.hpp"
#include "CellPropertyRegistry.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

// Needed to run in parallel
#include "PetscSetupAndFinalize.hpp"


/*
 * This function prints the memory usage at the time it is called.
 * `rPrefix` contains the rank (identifier) of the process
 */
void PrintMemoryUsage(const std::string& rPrefix)
{
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );

    double max_res = (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB

    std::cout << rPrefix << ": memory use = " << max_res <<  " MB.\n";
}

/*
 * == The test suite ==
 * This class was used to produce the results in Figure 4.
 * It constructs a large cell population and measures the amount of memory
 * in use. By running on larger numbers of processes, a proportionally
 * smaller amount of memory is used by each process.
 */
class TestMemoryUse : public AbstractCellBasedTestSuite
{
public:

    void TestProfile2dSimulation() throw (Exception)
    {
        /*
         * Record the rank of each process so that it can be output
         */
        std::ostringstream rank;
        rank << PetscTools::GetMyRank();

        /*
         *  Make a `NodesOnlyMesh` in which to store the nodes which are to be read from file
         */
        NodesOnlyMesh<2> mesh;
        mesh.SetCalculateNodeNeighbours(false);
        /*
         *  There is a 1024 y-range over which this mesh is split.  This means that each of 32 processes will take a y-range of 32.
         *  The maximum interaction distance will define the size of the indexing boxes and thus define
         *  the size of the strips which are the unit of parallel distribution
         */
        mesh.SetMaximumInteractionDistance(4);

        /*
         * Here we call a parallel helper method which reads the cell locations
         * from a file on disk.  This is a uniformly spread rectangle which is 1000 cells in the
         * x direction  and 1024 cells wide in the y direction.
         *
         * The proliferative type ensures that the cells are not growing, although this information is never
         * used.  No simulations are run.
         */
        std::vector<CellPtr> cells;
        ParallelCellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> generator;
        generator.GenerateParallelCells("projects/Harvey2014/test/data/1024000_2d_cells.dat",
                                        cells,
                                        mesh,
                                        CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        /*
         * Print the approximate memory use
         */
        PrintMemoryUsage(rank.str());
    }
};

#endif /*TESTMEMORYUSE_HPP_*/
