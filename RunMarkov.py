import ParameterClasses as P
import SupportMarkov as SupportMarkov
import HW9MarkovClasses as Cls
import ParameterClassesTreatment as Ptreatment

# print the outcomes of this simulated cohort
print("Problem 3")
cohort = Cls.Cohort(
    id=0,
    therapy=P.Therapies.NOTREATMENT)

#ProbMatrix = P.calculate_prob_matrix()
#print(ProbMatrix)
#P.add_background_mortality(ProbMatrix)

# simulate the cohort
simOutputs = cohort.simulate()
SupportMarkov.print_outcomes(simOutputs, 'No Treatment')


print("Problem 4")

print("TRANS_MATRIX = [    [0.75,    0.15,    0,       0.1],   #Well",
   " [0,      0,      1,       0],    #Stroke",
   " [0.0000,      0.1625,  0.701,   0.1365],",
   " [0,      0,      0,       1],#.65*.25 .2*1.05 ",
    "]")


print("Problem 5")
cohort2 = Cls.Cohort(
    id=1,
    therapy=P.Therapies.ANTICOAG)

# simulate the cohort
simOutputs2 = cohort2.simulate()
SupportMarkov.print_outcomes(simOutputs2, 'Anticoagulation Treatment')


print("Problem 6")
SupportMarkov.print_comparative_outcomes(simOutputs_none=simOutputs, simOutputs_anticoag=simOutputs2)
SupportMarkov.draw_survival_curves_and_histograms(simOutputs_none=simOutputs, simOutputs_anticoag=simOutputs2)




print("Problem 7")
simOutputs.print_mean_stroke_outcomes(therapyName="No Therapy")
simOutputs2.print_mean_stroke_outcomes(therapyName="Anticoagulation Therapy")
SupportMarkov.draw_stroke_histograms(simOutputs_none=simOutputs, simOutputs_anticoag=simOutputs2)
