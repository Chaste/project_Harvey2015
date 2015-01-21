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

#include "SemForce.hpp"
#include "CellPropertyCollection.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SemForce<ELEMENT_DIM,SPACE_DIM>::SemForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM,SPACE_DIM>(),
    mSpringStiffness(1e-8),    // 1e-9 Newtons per micron in Chaste units
    mNuclearSpringFactor(1.0)  // Defaults to homogeneous cytoplasm
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SemForce<ELEMENT_DIM,SPACE_DIM>::~SemForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SemForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Not used in this class for efficiency
    return zero_vector<double>(SPACE_DIM);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SemForce<ELEMENT_DIM,SPACE_DIM>::GetMagnitudeOfForce(double distanceBetweenNodes, double springStiffness)
{
    double exponent = exp(mNewmanScalingFactor * (1.0 - distanceBetweenNodes*distanceBetweenNodes / mRestLengthSquared));

    double coefficient = 4.0 * mNewmanScalingFactor * springStiffness * distanceBetweenNodes / mRestLengthSquared;

    return -1.0 * coefficient * exponent*(exponent - 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SemForce<ELEMENT_DIM, SPACE_DIM>::GetScaledSpringStiffness(unsigned semIndex)
{
    return SemParameterScaler<SPACE_DIM>::Instance()->ScaleSpringConstant(semIndex, mSpringStiffness);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SemForce<ELEMENT_DIM, SPACE_DIM>::SetSpringStiffness(double stiffness)
{
    assert(!(stiffness<0.0));
    mSpringStiffness = stiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SemForce<ELEMENT_DIM, SPACE_DIM>::SetRelativeSpringStiffness(unsigned semIndex, double relativeStiffness)
{
    if (relativeStiffness < 0.0)
    {
        EXCEPTION("Relative spring stiffness constant must be >= 0.0");
    }
    mRelativeSpringStiffnessMapping[semIndex] = relativeStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SemForce<ELEMENT_DIM, SPACE_DIM>::GetRelativeSpringStiffness(unsigned semIndex)
{
    std::map<unsigned, double>::iterator iter = mRelativeSpringStiffnessMapping.find(semIndex);

    if (iter == mRelativeSpringStiffnessMapping.end())
    {
        EXCEPTION("No relative spring stiffness set for sem cell");
    }

    return iter->second;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SemForce<ELEMENT_DIM, SPACE_DIM>::SetNuclearSpringFactor(double factor)
{
    assert(factor > 0.0);
    mNuclearSpringFactor = factor;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SemForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

    std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

    // Cache some constants to reduce multiple calls to SemParameterScaler
    mRestLength = SemParameterScaler<SPACE_DIM>::Instance()->GetRestLength();
    mRestLengthSquared = mRestLength * mRestLength;
    mNewmanScalingFactor = SemParameterScaler<SPACE_DIM>::Instance()->GetNewmanScalingFactor();

    for (typename std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >::iterator iter = r_node_pairs.begin();
        iter != r_node_pairs.end();
        iter++)
    {
        Node<SPACE_DIM>* node_a = iter->first;
        Node<SPACE_DIM>* node_b = iter->second;

        // Calculate the distance between nodes
        double distance_between_nodes = norm_2(node_b->rGetLocation() - node_a->rGetLocation());

        if (distance_between_nodes >= this->GetCutOffLength())
        {
            // Do nothing
        }
        else
        {
            unsigned sem_id_a = node_a->GetRegion();
            unsigned sem_id_b = node_b->GetRegion();

            double spring_stiffness;

            if (sem_id_a == sem_id_b)
            {
                double relative_stiffness;

                // If both elements are nuclear the force should be scaled using mNuclearSpringFactor
                if ((node_a->IsParticle()) && (node_b->IsParticle()))
                {
                    relative_stiffness = mNuclearSpringFactor;
                }
                else
                {
                    relative_stiffness = GetRelativeSpringStiffness(sem_id_a);
                }
                // Scaled u_0 based on number of ScEM elements per cell and relative scaling. u_0 in Sandersius and Newman Phys. Biol. 5:015002, (2008)
                spring_stiffness = relative_stiffness * GetScaledSpringStiffness(sem_id_a);
            }
            else
            {
                double relative_stiffness = 0.5; // Interactions between different cell types are weaker

                // Scaled u_0 based on number of ScEM elements per cell. u_0 in Sandersius and Newman Phys. Biol. 5:015002, (2008)
                double scaled_stiffness_a = GetScaledSpringStiffness(sem_id_a);
                double scaled_stiffness_b = GetScaledSpringStiffness(sem_id_b);

                spring_stiffness = relative_stiffness * (0.5 * (scaled_stiffness_a + scaled_stiffness_b) );
            }
            double magnitude_of_force = GetMagnitudeOfForce(distance_between_nodes, spring_stiffness);

            // For rheology plates, double the force with a non-rheology plate to ensure that they adhere to the cell
            if ((sem_id_a > UINT_MAX -3) != (sem_id_b > UINT_MAX -3))
            {
                magnitude_of_force *= 2.0;
            }

            // We divide by distance_between_nodes to normalise (node_b->rGetLocation - node_a->rGetLocation) which saves creating an extra temp c_vector
            c_vector<double, SPACE_DIM> force = (magnitude_of_force / distance_between_nodes) * (node_b->rGetLocation() - node_a->rGetLocation());
            node_a->AddAppliedForceContribution(force);
            force *= -1.0;
            node_b->AddAppliedForceContribution(force);
        }
    }
}

// Explicit instantiation
template class SemForce<1,1>;
template class SemForce<1,2>;
template class SemForce<2,2>;
template class SemForce<1,3>;
template class SemForce<2,3>;
template class SemForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SemForce)
