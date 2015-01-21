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


#include "SemMesh.hpp"

template<unsigned SPACE_DIM>
SemMesh<SPACE_DIM>::SemMesh()
    : NodesOnlyMesh<SPACE_DIM>()
{
}

template<unsigned SPACE_DIM>
void SemMesh<SPACE_DIM>::SetInitialBoxCollection(const c_vector<double, 2*SPACE_DIM> domainSize, double maxInteractionDistance)
{
    this->SetUpBoxCollection(maxInteractionDistance, domainSize);
}

template<unsigned SPACE_DIM>
bool SemMesh<SPACE_DIM>::IsNodeLocallyOwned(Node<SPACE_DIM>* pNode)
{
    assert(this->GetBoxCollection());

    return this->GetBoxCollection()->IsOwned(pNode);
}

template<unsigned SPACE_DIM>
void SemMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes, double maxInteractionDistance)
{
    /*
     * Call the base class method, and then copy the region indices from the nodes.
     */
    NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(rNodes, maxInteractionDistance);
    std::vector<bool> initial_nodes = this->rGetInitiallyOwnedNodes();

    unsigned j = 0;
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        if (initial_nodes[i])
        {
            unsigned sem_id = rNodes[i]->GetRegion();
            this->mNodes[j]->SetRegion(sem_id);
            bool is_particle = rNodes[i]->IsParticle();
            this->mNodes[j]->SetIsParticle(is_particle);
            j++;
        }
    }
}

// Explicit instantiation
template class SemMesh<1>;
template class SemMesh<2>;
template class SemMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemMesh)
