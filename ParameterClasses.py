from enum import Enum
import numpy as np
import scipy.stats as stat
import math as math
import InputDataNoTreatment as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.ProbDistParEst as Est


class HealthStats(Enum):
    """ health states of patients with Atrial Fibrillation """
    WELL = 0
    STROKE = 1
    POST_STROKE = 2
    STROKE_DEATH = 3
    BACKGROUND_DEATH = 4

class Therapies(Enum):
    """ no vs. anticoagualtion therapy """
    NOTREATMENT = 0
    ANTICOAG = 1


class _Parameters:

    def __init__(self, therapy):

        # selected therapy
        self._therapy = therapy

        # simulation time step
        self._delta_t = Data.DELTA_T

        # calculate the adjusted discount rate
        self._adjDiscountRate = Data.DISCOUNT*Data.DELTA_T

        # initial health state
        self._initialHealthState = HealthStats.WELL

        # annual treatment cost
        if self._therapy == Therapies.NOTREATMENT:
            self._annualTreatmentCost = Data.NO_TREATMENT_COST
        else:
            self._annualTreatmentCost = Data.NO_TREATMENT_COST + Data.ANTI_COAG_COST

        # transition probability matrix of the selected therapy
        self._prob_matrix = []
        # treatment relative risk
        self._treatmentRR = 0

        # annual state costs and utilities
        self._annualStateCosts = []
        self._annualStateUtilities = []

    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t

    def get_adj_discount_rate(self):
        return self._adjDiscountRate

    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]

    def get_annual_state_cost(self, state):
        if state == HealthStats.STROKE_DEATH or state == HealthStats.BACKGROUND_DEATH:
            return 0
        else:
            return self._annualStateCosts[state.value]

    def get_annual_state_utility(self, state):
        if state == HealthStats.STROKE_DEATH or state == HealthStats.BACKGROUND_DEATH:
            return 0
        else:
            return self._annualStateUtilities[state.value]

    def get_annual_treatment_cost(self):
        return self._annualTreatmentCost


class ParametersFixed(_Parameters):
    def __init__(self, therapy):

        # initialize the base class
        _Parameters.__init__(self, therapy)

        # calculate transition probabilities between stroke states
        self._prob_matrix = calculate_prob_matrix()
        # add background mortality if needed
        if Data.ADD_BACKGROUND_MORT:
            add_background_mortality(self._prob_matrix)

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == Therapies.ANTICOAG:
            # treatment relative risk
            self._treatmentRR = Data.TREATMENT_RR
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = calculate_prob_matrix_anticoag(
                matrix_none=self._prob_matrix, anticoag_rr=Data.TREATMENT_RR)

        # annual state costs and utilities
        self._annualStateCosts = Data.ANNUAL_STATE_COST
        self._annualStateUtilities = Data.ANNUAL_STATE_UTILITY


class ParametersProbabilistic(_Parameters):
    def __init__(self, seed, therapy):

        # initializing the base class
        _Parameters.__init__(self, therapy)

        self._rng = Random.RNG(seed)    # random number generator to sample from parameter distributions
        self._strokeProbMatrixRVG = []  # list of dirichlet distributions for transition probabilities
        self._lnRelativeRiskRVG = None  # random variate generator for the treatment relative risk
        self._annualStateCostRVG = []       # list of random variate generators for the annual cost of states
        self._annualStateUtilityRVG = []    # list of random variate generators for the annual utility of states

        # Stroke transition probabilities
        j = 0
        for prob in Data.TRANS_MATRIX:
            self._strokerobMatrixRVG.append(Random.Dirichlet(prob[j:]))
            j += 1

        # treatment relative risk
        # find the mean and st_dev of the normal distribution assumed for ln(RR)
        sample_mean_lnRR = math.log(Data.TREATMENT_RR)
        sample_std_lnRR = (Data.TREATMENT_RR_CI[1]-Data.TREATMENT_RR_CI[0])/(2*stat.norm.ppf(1-0.05/2))
        self._lnRelativeRiskRVG = Random.Normal(mean=sample_mean_lnRR, st_dev=sample_std_lnRR)

        # annual state cost
        for cost in Data.ANNUAL_STATE_COST:
            # find shape and scale of the assumed gamma distribution
            shape, scale = Est.get_gamma_parameters(mean=cost, st_dev=cost/ 4)
            # append the distribution
            self._annualStateCostRVG.append(Random.Gamma(shape, scale))

        # annual state utility
        for utility in Data.ANNUAL_STATE_UTILITY:
            # find alpha and beta of the assumed beta distribution
            a, b = Est.get_beta_parameters(mean=utility, st_dev=utility/5)
            # append the distribution
            self._annualStateUtilityRVG.append(Random.Beta(a, b))

        # resample parameters
        self.__resample()

    def __resample(self):

        # calculate transition probabilities
        # create an empty matrix populated with zeroes
        self._prob_matrix = []
        for s in HealthStats:
            self._prob_matrix.append([0] * len(HealthStats))

        # for all health states
        for s in HealthStats:
            # if the current state is death
            if s in [HealthStats.STROKE_DEATH, HealthStats.BACKGROUND_DEATH]:
                # the probability of staying in this state is 1
                self._prob_matrix[s.value][s.value] = 1
            else:
                # sample from the dirichlet distribution to find the transition probabilities between stroke states
                dist = self._strokeProbMatrixRVG[s.value]
                sample = dist.sample(self._rng)
                for j in range(len(sample)):
                    self._prob_matrix[s.value][s.value+j] = sample[j]

        # add background mortality if needed
        if Data.ADD_BACKGROUND_MORT:
            add_background_mortality(self._prob_matrix)

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == Therapies.ANTICOAG:
            # treatment relative risk
            self._treatmentRR = math.exp(self._lnRelativeRiskRVG.sample(self._rng))
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = calculate_prob_matrix_anticoag(
                matrix_none=self._prob_matrix, anticoag_rr=self._treatmentRR)

        # sample from gamma distributions that are assumed for annual state costs
        self._annualStateCosts = []
        for dist in self._annualStateCostRVG:
            self._annualStateCosts.append(dist.sample(self._rng))

        # sample from beta distributions that are assumed for annual state utilities
        self._annualStateUtilities = []
        for dist in self._annualStateUtilityRVG:
            self._annualStateUtilities.append(dist.sample(self._rng))


def calculate_prob_matrix():
    """ :returns transition probability matrix for stroke states under no therapy"""

    # create an empty matrix populated with zeroes
    prob_matrix = []
    for s in HealthStats:
        prob_matrix.append([0] * len(HealthStats))

    # for all health states
    for s in HealthStats:
        # if the current state is death
        if s in [HealthStats.STROKE_DEATH, HealthStats.BACKGROUND_DEATH]:
            # the probability of staying in this state is 1
            prob_matrix[s.value][s.value] = 1
        else:
            # calculate total counts of individuals
            sum_counts = sum(Data.TRANS_MATRIX[s.value])
            # calculate the transition probabilities out of this state
            for j in range(s.value, HealthStats.BACKGROUND_DEATH.value):
                prob_matrix[s.value][j] = Data.TRANS_MATRIX[s.value][j] / sum_counts

    return prob_matrix


def add_background_mortality(prob_matrix):

    # find the transition rate matrix
    rate_matrix = MarkovCls.discrete_to_continuous(prob_matrix, 1)
    # add mortality rates
    for s in HealthStats:
        if s not in [HealthStats.STROKE_DEATH, HealthStats.BACKGROUND_DEATH]:
            rate_matrix[s.value][HealthStats.BACKGROUND_DEATH.value] \
                = -np.log(1 - Data.ANNUAL_PROB_BACKGROUND_MORT)

    # convert back to transition probability matrix
    prob_matrix[:], p = MarkovCls.continuous_to_discrete(rate_matrix, Data.DELTA_T)
    # print('Upper bound on the probability of two transitions within delta_t:', p)


def calculate_prob_matrix_anticoag(matrix_none, anticoag_rr):
    """
    :param matrix_none: (list of lists) transition probability matrix under no therapy
    :param anticoag_rr: relative risk of the combination treatment
    :returns (list of lists) transition probability matrix under combination therapy """

    # create an empty list of lists
    matrix_anticoag = []
    for l in matrix_none:
        matrix_anticoag.append([0] * len(l))

    # populate the anticoag matrix
    # first non-diagonal elements
    for s in HealthStats:
        for next_s in range(s.value + 1, len(HealthStats)):
            matrix_anticoag[s.value][next_s] = anticoag_rr * matrix_none[s.value][next_s]

    # diagonal elements are calculated to make sure the sum of each row is 1
    for s in HealthStats:
        if s not in [HealthStats.STROKE_DEATH, HealthStats.BACKGROUND_DEATH]:
            matrix_anticoag[s.value][s.value] = 1 - sum(matrix_none[s.value][s.value + 1:])

    return matrix_anticoag
