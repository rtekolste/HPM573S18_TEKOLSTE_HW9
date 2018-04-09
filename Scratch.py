import ParameterClasses as P
import SupportMarkov as SupportMarkov
import HW9MarkovClasses as Cls
import ParameterClassesTreatment as Ptreatment
import InputDataNoTreatment as Data

# print the outcomes of this simulated cohort
#print("Problem 3")
#cohort = Cls.Cohort(
#    id=0,
#    therapy=P.Therapies.NOTREATMENT)

P._Parameters.get_transition_prob(P.Therapies.NOTREATMENT)
