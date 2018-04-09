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
        self._prob_matrix = Data.TRANS_MATRIX
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

    def get_transition_prob(self):
        return self._prob_matrix

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
