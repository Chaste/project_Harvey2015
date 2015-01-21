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

#ifndef SEMCELLSGENERATOR_HPP_
#define SEMCELLSGENERATOR_HPP_

#include "Cell.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileFinder.hpp"
#include "Node.hpp"
#include "RandomNumberGenerator.hpp"
#include "SemMesh.hpp"
#include "UblasVectorInclude.hpp"
#include "WildTypeCellMutationState.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cmath>
#include <vector>

/**
 * @class SemCellsGenerator
 *
 * This class generates Sem type cells from an archive of element
 * locations for a template cell.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class SemCellsGenerator
{
protected:
    /** Test is a friend to allow easier testing. */
    friend class TestSemCellsGenerator;

    /** The number of cells to be generated. */
    unsigned mNumCells;

    /** The locations of the centre of mass of each cell */
    std::vector<c_vector<double, DIM> > mCellLocations;

    /** The relative size of each cell. */
    std::vector<double> mCellScalings;

    /**
     * Generate a basic Chaste cell to attach with a node.
     */
    CellPtr GenerateCell();

    /**
     * Given an archive location, get the bounding box that surrounds all the node locations.
     * @param archivePath the path to the archive.
     * @return the boundaing box.
     */
    c_vector<double, 2*DIM> GetArchiveBoundingBox(std::string archivePath);

public:

    /**
     * Default constructor.
     */
    SemCellsGenerator();

    /**
     * Set the total number of ScEM cells to be created.
     *
     * @param num_cells the number of cells.
     */
    void SetNumCells(unsigned num_cells);

    /**
     * Set the centre of mass for each Sem cell.
     *
     * @param semIndex the Sem index of the cell
     * @param location its centre of mass (assuming all elements have equal 'mass')
     */
    void SetCellLocation(unsigned semIndex, c_vector<double, DIM> location);

    /**
     * Set the scale factor for the cell area / volume. Defaults to 1.0.
     *
     * @param semIndex the index of the cell to scale.
     * @param The scale factor. Must be < 1.0.
     */
    void SetCellScaleFactor(unsigned semIndex, double factor);

    /**
     * Generate nodes and cells for a given decomposition of space, dictated by mesh.
     *
     * @param archiveLocation the name of the file containing the template element locations for a cell
     * @param cells a (typically empty) vector into which cells can be inserted.
     * @param mesh a pointer to the mesh that will store the nodes.
     */
    void GenerateSemCells(std::string archiveLocation, std::vector<CellPtr>& cells, SemMesh<DIM>& mesh);
};

#endif /* SEMCELLSGENERATOR_HPP_ */
