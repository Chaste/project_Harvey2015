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

#include "SemCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

#include <fstream>

template<class CELL_CYCLE_MODEL, unsigned DIM>
SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::SemCellsGenerator()
    : mNumCells(0u)
{
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::SetNumCells(unsigned num_cells)
{
    mNumCells = num_cells;
    mCellLocations.resize(num_cells, zero_vector<double>(DIM));
    mCellScalings.resize(num_cells, 1.0);
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::SetCellLocation(unsigned semIndex, c_vector<double, DIM> location)
{
    assert(semIndex < mNumCells);
    mCellLocations[semIndex] = location;
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::SetCellScaleFactor(unsigned semIndex, double factor)
{
    assert(semIndex < mNumCells);
    assert(0.0 < factor);
    assert(!(factor > 1.0));
    mCellScalings[semIndex] = factor;
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::GenerateSemCells(std::string archiveLocation, std::vector<CellPtr>& cells, SemMesh<DIM>& mesh)
{
    // Get a bounding box for the archived nodes
    c_vector<double, 2*DIM> base_bounding_box = GetArchiveBoundingBox(archiveLocation);

    // Enlarge the bounding box to the size of the enscribing circle, to allow for rotations
    c_vector<double, DIM> min_corner;
    c_vector<double, DIM> max_corner;
    for (unsigned i=0; i<DIM; i++)
    {
        min_corner[i] = base_bounding_box[2*i];
        max_corner[i] = base_bounding_box[2*i+1];
    }
    double max_mod = norm_2(max_corner);
    double min_mod = norm_2(min_corner);

    double box_size = std::max(max_mod, min_mod);

    for (unsigned i=0; i<DIM; i++)
    {
        base_bounding_box[2*i] = -box_size;
        base_bounding_box[2*i+1] = box_size;
    }

    c_vector<double, 2*DIM> bounding_box;
    for (unsigned i=0; i<DIM; i++)
    {
        bounding_box[2*i] = DBL_MAX;
        bounding_box[2*i+1] = -DBL_MAX;
    }

    for (unsigned i=0; i<mNumCells; i++)
    {
        c_vector<double, 2*DIM> cell_bounding_box;

        for (unsigned j=0; j<DIM; j++)
        {
            cell_bounding_box[2*j] = base_bounding_box[2*j] + mCellLocations[i][j];
            cell_bounding_box[2*j+1] = base_bounding_box[2*j+1] + mCellLocations[i][j];
        }
        for (unsigned j=0; j<DIM; j++)
        {
            bounding_box[2*j] = std::min(cell_bounding_box[2*j], bounding_box[2*j]);
            bounding_box[2*j+1] = std::max(cell_bounding_box[2*j+1], bounding_box[2*j+1]);
        }
    }

    // Construct box collection of Sem mesh
    mesh.SetInitialBoxCollection(bounding_box, mesh.GetMaximumInteractionDistance());

    unsigned node_index = 0;
    for (unsigned i=0; i<mNumCells; i++)
    {
        std::ifstream infile(archiveLocation.c_str(), std::ios::in);

        std::string line;
        std::getline(infile, line);

        std::istringstream iss(line);

        unsigned num_lines;
        unsigned file_dimension;

        iss >> num_lines >> file_dimension;

        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            c_vector<double, DIM> location;

            for (unsigned k=0; k<DIM; k++)
            {
                iss >> location[k];
            }

            // Apply scaling to volume. Assumes a square initial cell in region -0.5, 0.5 for now.
            double scale_factor = mCellScalings[i];
            bool accepted_node = true;

            double cell_width = base_bounding_box[2*DIM-1] - base_bounding_box[2*DIM-2];
            double cell_mean = 0.5 * (base_bounding_box[2*DIM-1] + base_bounding_box[2*DIM-2]);

            if (fabs(location[1] - cell_mean) > scale_factor*0.5*cell_width)
            {
                accepted_node = false;
            }

            // Apply translation
            location += mCellLocations[i];

            Node<DIM>* p_node = new Node<DIM>(node_index, location);
            if (mesh.IsNodeLocallyOwned(p_node) && accepted_node)
            {
                p_node->SetRegion(i);

                mesh.AddNode(p_node);
                CellPtr p_cell = GenerateCell();
                cells.push_back(p_cell);
            }
            else
            {
                delete p_node;
            }
            node_index++;
        }

        infile.close();
    }

    NodeMap map(node_index+1);
    mesh.ReMesh(map);
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
c_vector<double, 2*DIM> SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::GetArchiveBoundingBox(std::string archivePath)
{
    c_vector<double, 2*DIM> bounding_box;
    for (unsigned i=0; i<DIM; i++)
    {
        bounding_box[2*i] = DBL_MAX;
        bounding_box[2*i + 1] = -DBL_MAX;
    }

    std::ifstream infile(archivePath.c_str());

    std::string line;
    std::getline(infile, line);

    std::istringstream iss(line);

    unsigned num_lines;
    unsigned file_dimension;

    iss >> num_lines >> file_dimension;

    if (file_dimension != DIM)
    {
        EXCEPTION("Space dimension of SemCellsGenerator and archive file do not match");
    }

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        for (unsigned i=0; i<DIM; i++)
        {
            double point;
            iss >> point;
            bounding_box[2*i] = (point < bounding_box[2*i]) ? point : bounding_box[2*i];
            bounding_box[2*i+1] = (point > bounding_box[2*i+1]) ? point : bounding_box[2*i+1];
        }
    }

    infile.close();

    return bounding_box;
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
CellPtr SemCellsGenerator<CELL_CYCLE_MODEL, DIM>::GenerateCell()
{
    CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
    p_cell_cycle_model->SetDimension(DIM);

    boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
    p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

    return p_cell;
}

// Explicit instantiation
template class SemCellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1>;
template class SemCellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2>;
template class SemCellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3>;
