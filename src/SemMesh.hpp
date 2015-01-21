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

#ifndef SEMMESH_HPP_
#define SEMMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>

#include "NodesOnlyMesh.hpp"

/**
 * Mesh class for storing lists of nodes (no elements) that additionally
 * stores an integer index with each node.
 */
template<unsigned SPACE_DIM>
class SemMesh: public NodesOnlyMesh<SPACE_DIM>
{
public:

    /** 
     * Default contructor 
     *
     * See documentation for NodesOnlyMesh::NodesOnlyMesh
     */
    SemMesh();

    /** 
     * Overridden ConstructNodesWithoutMesh to keep the 
     * integer index associated with nodes.
     * 
     * See documentation for NodesOnlyMesh::ConstructNodesWithoutMesh 
     */
    void ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes,
                                   double maxInteractionDistance);

    /**
     * Set the initial distribution of boxes for the mesh.
     * This is a helper method for constructing large
     * parallel meshes.
     *
     * @param [in] domainSize              the axis aligned bounding box
     *                                     for the domain size, specified
     *                                     by its extreme corners.
     * @param [in] maxInteractionDistance  the maximum distance between
     *                                     two interacting cells.
     */
    void SetInitialBoxCollection(c_vector<double, 2*SPACE_DIM> domainSize, 
                                 double maxInteractionDistance);

    /**
     * A helper method for constructing nodes
     * and identifying whether a node is in the local spatial
     * region of this process. 
     */
    bool IsNodeLocallyOwned(Node<SPACE_DIM>* pNode);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemMesh)

#endif /*SEMMESH_HPP_*/
