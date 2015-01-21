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

#ifndef TESTLOADBALANCESENSITIVITY_HPP_
#define TESTLOADBALANCESENSITIVITY_HPP_
/*
 * = Test the load-balancing algorithm, and its sensitivity to the rebalancing frequency  (Figure 6).=
 * 
 * This test suites demonstrates the use of the load-balancing algorithm
 * for parallel cell-based simulations.
 *
 * = Use =
 * 
 * The test in this class is designed to be run using multiple processes to
 * recreate figure 6.
 * 
 * For example:
 * {{{
 * scons build=GccOptNative_2 test_suite=projects/Harvey2014/tests/TestLoadBalanceSensitivityLiteratePaper.hpp
 * }}}
 * or, to produce the 8-way simulation of Figure 6c and Figure 7,
 * {{{
 * scons build=GccOptNative_8 projects/Harvey2014/tests/TestLoadBalanceSensitivityLiteratePaper.hpp
 * }}}
 */

/*
 * === Header files ===
 */
// Used for the test-suite set up
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Cell-based headers
#include "OffLatticeSimulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellBasedEventHandler.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"

// Used to set up the simulation in a parallel environment
#include "PetscSetupAndFinalize.hpp"

/**
 * == The test suite ==
 * This class was used to produce the results of Figure 6.
 *
 * It runs a cell population simulation with different load-balancing
 * frequencies and is able to record the relative speed of the simulation
 * as the frequency is varied.
 */
class TestLoadBalanceSensitivity : public AbstractCellBasedTestSuite
{
public:

    /*
     * == The first unit test ==
     *
     * This test starts a population from a seed of 25 cells, and allows
     * the population to grow for 100 hours. This test should be run on 
     * progressively more processes to experience the speed-up
     * show in Figure 6(a).
     */
    void TestParallelLoadBalance () throw (Exception)
    {
        /*
         * '''Warning note: this test may take of the order of an hour to run''' (without the accompanying test in the suite) and
         *  will not output to the screen during that time.
         *
         * This test should be run in parallel,
         * so we enforce this condition with an exception.
         */
        if (PetscTools::GetNumProcs() == 1)
        {
          EXCEPTION ("This code should be run in parallel.");
        }
        else
        {
            if (PetscTools::AmMaster())
            {
                std::cout << "Running load balancing test on "
                          << PetscTools::GetNumProcs()
                          << " processes" << std::endl;
            }
        }
        /*
         * The initial geometry is 25 cells in a 'honeycomb' configuration,
         * such that the cells are close to a stable equilibrium.
         */
        std::vector< Node<2>* > nodes;
        for (unsigned i = 0; i < 5; ++i)
        {
            for (unsigned j = 0; j < 5; ++j)
            {
                unsigned index = j + 5 * i;
                c_vector<double, 2> location;
                location[0] = static_cast<double> (i) / 2.0 + 0.25 * (j % 2);
                location[1] = static_cast<double> (j) / 2.0;
                nodes.push_back(new Node<2>(index, location, false));
            }
        }

        /*
         * Using the cell locations configured above, we create a mesh (a collection of nodes).
         * that defines the geometry of the population.
         */
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        mesh.SetCalculateNodeNeighbours(false);

        /*
         * Having created the geometry, we then assign one cell to each node in the
         * mesh. Each cell within the Chaste library has two states which must be 
         * specified in the initial condition: the mutation state and the proliferative type.
         *
         * Using `WildTypeCellMutationState` simply denotes that each cell has no defined
         * mutation, and behaves 'normally'.
         * 
         * Using `TransitCellProliferativeType` gives the cell a proliferative type such that
         * it can carry on dividing up to a maximum number of transit generations. For the purpose
         * of this test, it allows us to ensure that the cells keep dividing, and so the
         * population will keep growing.
         */
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProliferativeType> p_type(new TransitCellProliferativeType);

        for (unsigned i = 0; i < mesh.GetNumNodes(); ++i)
        {
            /*
             * The cells are given a cell cycle model that causes them to divide after
             * a fixed duration has elapsed.
             */
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetStemCellG1Duration(4.0);
            p_model->SetMaxTransitGenerations(20u);

            CellPtr p_cell(new Cell(p_state, p_model));

            /*
             * Cells are configured to be born at some random time in the past,
             * so that cell division is not synchronised between all cells.
             */
            double birth_time = -1.0 * RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(p_type);
            cells.push_back(p_cell);
        }

        /* Set up a cell population from the geometry defined in the mesh and the set of cells created above. */
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.0);
        cell_population.SetOutputResultsForChasteVisualizer(false);

        /*
         * Here we set a flag to make the population rebalance the distribution of 
         * cells between processes during the simulation. To see the scaling
         * characteristics of this test without load-balancing, these two lines
         * should be removed or commented out.
         */
        cell_population.SetLoadBalanceMesh(true);
        cell_population.SetLoadBalanceFrequency(1000);

        /* 
         * Make a simulation object from the cell population, and set the end time to 99 hours,
         * this allows us to reset the timers and measure the speed-up in the final hour
         * of simulation.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetEndTime(99.0);
        simulator.SetDt(1.0/1200.0);
        simulator.SetSamplingTimestepMultiple(1200);

        /*
         * Set output directory
         * This folder is to be found relative to `CHASTE_TEST_OUTPUT` which by default is
         * `/tmp/$USER/testoutput`
         */
        std::ostringstream exp_num_str;
        exp_num_str << PetscTools::GetNumProcs();
        std::string output_directory = "LoadBalance_" + exp_num_str.str();
        simulator.SetOutputDirectory(output_directory);

        /* 
         * Create a force law and pass it to the `OffLatticeSimulation`. 
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        /*
         * We first solve the system up to 99 hours.
         */
        simulator.Solve();

        /*
         * And then reset the timers to measure the final hour.
         */
        CellBasedEventHandler::Report();
        CellBasedEventHandler::Reset();

        simulator.SetEndTime(100.0);
        simulator.Solve();

        CellBasedEventHandler::Report();
        SimulationTime::Destroy();

        /*
         * Print the number of local cells to screen.
         */
        if (PetscTools::AmMaster())
        {
            std::cout << "Number of cells on process "
                  << PetscTools::GetMyRank () << ": "
                  << cell_population.GetNumNodes()
                  << std::endl;
        }
    }

    /*
     * == The second unit test ==
     *
     * This test starts a population from a seed of 25 cells, and allows
     * the population to grow for 100 hours. The entire simulation is 
     * run for different values of the load balancing frequency, which 
     * is the rate at which the distribution of cells is rebalanced between
     * each of the processes.
     *
     * Some post-processing is required to produce the data for plot in Figure 6c.
     */
    void TestParallelLoadBalanceSensitivity () throw(Exception)
    {
        /*
         * '''Warning note: each iteration of the main "rebalancing_frequency" loop may take of the order of 30 minutes to run''' (without the accompanying test in the suite)
         *  as can be seen from the vertical axis of Figure 6c.  This test will output intermittently and that output may well be buffered.
         *
          * This test should be run in parallel,
         * so we enforce this condition with an exception.
         * Figure 6c was generated using 8 processes.
         */
        if (PetscTools::IsSequential())
        {
          EXCEPTION ("This code should be run in parallel.");
        }

        /*
         * We run the test for rebalancing frequencies from every 10 time-steps,
         * up to every 100,000 time-steps.
         *
         * We multiply by 5 or 2 on alternating iterations - see bottom of loop.
         */
        bool even_loop = false;
        for (unsigned rebalancing_frequency = 10; rebalancing_frequency <= 100000; even_loop = !even_loop /*rebalancing_frequency*=10*/)
        {
            /*
             * Print the rebalancing frequency for reference in final results.
             */
            if (PetscTools::AmMaster())
            {
                std::cout << "Rebalancing frequency: " << rebalancing_frequency << "\n";
            }

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            CellBasedEventHandler::Reset();

            /*
             * The initial geometry 25 cells in a 'honeycomb' configuration,
             * such that the cells are close to a stable equilibrium.
             */
            std::vector< Node<2>* > nodes;
            for (unsigned i = 0; i < 5; ++i)
            {
                for (unsigned j = 0; j < 5; ++j)
                {
                    unsigned index = j + 5 * i;
                    c_vector<double, 2> location;
                    location[0] = static_cast<double> (i) / 2.0 + 0.25 * (j % 2);
                    location[1] = static_cast<double> (j) / 2.0;
                    nodes.push_back(new Node<2>(index, location, false));
                }
            }

            /*
             * Using the cell locations configured above, we create a mesh (a collection of nodes).
             * that defines the geometry of the population.
             */
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);
            mesh.SetCalculateNodeNeighbours(false);

            /*
             * Having created the geometry, we then assign one cell to each node in the
             * mesh.
             * See the previous test for a brief overview of the states assigned to cells.
             */
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
            boost::shared_ptr<AbstractCellProliferativeType> p_type(new TransitCellProliferativeType);

            for (unsigned i = 0; i < mesh.GetNumNodes(); ++i)
            {
                /*
                 * The cells are given a cell cycle model that causes them to divide after
                 * a fixed duration has elapsed.
                 */
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetDimension(2);
                p_model->SetStemCellG1Duration(4.0);
                p_model->SetMaxTransitGenerations(20u);

                CellPtr p_cell(new Cell(p_state, p_model));

                /* 
                 * Cells are configured to be born at some random time in the past,
                 * so that cell division is not synchronised between all cells.
                 */
                double birth_time = -1.0 * RandomNumberGenerator::Instance()->ranf()*18.0;
                p_cell->SetBirthTime(birth_time);
                p_cell->SetCellProliferativeType(p_type);
                cells.push_back(p_cell);
            }

            /* Set up a cell population from the geometry defined in the mesh, and the set of cells. */
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.SetAbsoluteMovementThreshold(1.0);
            cell_population.SetOutputResultsForChasteVisualizer(false);
            cell_population.SetLoadBalanceMesh(true);
            cell_population.SetLoadBalanceFrequency(rebalancing_frequency);

            /* 
             * Make a simulation object from the cell population, and set the end time to 99 hours,
             * this allows us to reset the timers and measure the speed-up in the final hour
             * of simulation.
             */
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetEndTime(99.0);
            simulator.SetDt(1.0/1200.0);
            simulator.SetSamplingTimestepMultiple(1200);

            /*
             * Set output directory
             * This folder is to be found relative to `CHASTE_TEST_OUTPUT` which by default is
             * `/tmp/$USER/testoutput`
             */
            std::ostringstream exp_num_str;
            exp_num_str << PetscTools::GetNumProcs() << "_" << rebalancing_frequency;
            std::string output_directory = "LoadBalSensitivity_" + exp_num_str.str(); 

            simulator.SetOutputDirectory(output_directory);

            /* 
             * Create a force law and pass it to the `OffLatticeSimulation`. 
             */
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            p_force->SetCutOffLength(1.5);
            simulator.AddForce(p_force);

            /*
             * We first solve the system up to 99 hours.
             */
            simulator.Solve();

            /*
             * Report, and then reset the timers to measure the final hour.
             */
            if (PetscTools::AmMaster())
            {
                std::cout<<"Report for bulk of simulation (99 simulation hours)\n";
            }
            CellBasedEventHandler::Report();
            CellBasedEventHandler::Reset();

            simulator.SetEndTime(100.0);
            simulator.Solve();
            if (PetscTools::AmMaster())
            {
                std::cout<<"Report for final simulation hour\n";
            }
            CellBasedEventHandler::Report();

            /*
             * We multiply by 5 or 2 on alternating iterations
             */
            if (even_loop)
            {
                rebalancing_frequency *= 2;
            }
            else
            {
                rebalancing_frequency *= 5;
            }
        }
    }
};

#endif /*TESTLOADBALANCESENSITIVITY_HPP_*/
