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

#ifndef SEMFORCE_HPP_
#define SEMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SemParameterScaler.hpp"

#include <boost/serialization/base_object.hpp>

/**
 * @class SemForce
 * 
 * A force class used for a fine-grained cell centre simulation
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class SemForce : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /** Needed for testing. */
    friend class TestSemForce;

    /** Spring constant for interactions between elements. */
    double mSpringStiffness;

    /** Scaling constant for force between two nuclear elements. */
    double mNuclearSpringFactor;

    /** A map of the relative spring stiffness for different ScEM cells */
    std::map<unsigned, double> mRelativeSpringStiffnessMapping;

    /** A cache of the rest length between nodes */
    double mRestLength;

    /** A cache of squared rest length */
    double mRestLengthSquared;

    /** A cache of Newman scaling factor */
    double mNewmanScalingFactor;

    /** 
     * Get the force magnitude between two nodes.
     *
     * @param [in] distanceBetweenNodes       the distance between the cells.
     * @param [in] springStiffness            stiffness of the spring between the cells.
     * @return the magnitude of the force.
     */
    double GetMagnitudeOfForce(double distanceBetweenNodes, double springStiffness);

protected:

    /**
     * Get scaled stiffness between elements so that bulk cell
     * properties are independent of number of elements.
     * @param semIndex the index of the cell
     * @return the scaled value of u_0.
     */
    double GetScaledSpringStiffness(unsigned semIndex);

public:

    /**
     * Constructor.
     */
    SemForce();

    /**
     * Destructor.
     */
    virtual ~SemForce();

    /**
     * WE DO NOT USE THIS VIRTUAL METHOD IN THE CLASS FOR EFFICIENY.
     * We avoid the case of getting a node, finding its index, passing
     * index to this method, then getting the node again...
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM>
    CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                               unsigned nodeBGlobalIndex,
                               AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Set the spring stiffness k_0 in Sandersius and Newman Phys. Biol. 5:015002, (2008).
     * @param stiffness the value of the stiffness.
     */
    void SetSpringStiffness(double stiffness);

    /**
     * Set the value of the relative spring stiffness for different cell types.
     *
     * @param semIndex the index of the cell type to set
     * @param relativeStiffness the multiplication factor for the spring stiffness
     */
    void SetRelativeSpringStiffness(unsigned semIndex, double relativeStiffness);

    /**
     * @param semIndex the index of the ScEM cell
     * @return the spring stiffness.
     */
    double GetRelativeSpringStiffness(unsigned semIndex);

    /**
     * Set the nuclear spring factor for cell nuclei.
     * @param factor the value to set it to.
     */
    void SetNuclearSpringFactor(double factor);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SemForce)

#endif /*SEMFORCE_HPP_*/
