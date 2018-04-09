
# simulation settings
POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 50    # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
DISCOUNT = 0.0     # annual discount rate

ADD_BACKGROUND_MORT = True  # if background mortality should be added
DELTA_T = 1       # years. I did not change this

PSA_ON = False      # if probabilistic sensitivity analysis is on

# transition matrix
TRANS_MATRIX = [
    [.75,  .15,    0,    .10],   # Well
    [0,     0,     1.00,     0],   # Stroke
    [0,    .25,    .55,  .20],  # Post Stroke
    [0,    0,      0,    1.00],
    ]


# annual cost of each health state - No data yet
ANNUAL_STATE_COST = [
    0,   # Well
    0,   # Stroke
    0    # Post-Stroke
    ]

# annual health utility of each health state - No data yet
ANNUAL_STATE_UTILITY = [
    0,   # Well
    0,   # Stroke
    0    # Post-Stroke
    ]

# annual drug costs - No data yet
NO_TREATMENT_COST = 0
ANTI_COAG_COST = 0

# treatment relative risk - #4
TREATMENT_RR = 1
#TREATMENT_RR_CI = 0, 0  # lower 95% CI, upper 95% CI

# annual probability of background mortality (number per year per 1,000 population) - keeping the same
ANNUAL_PROB_BACKGROUND_MORT = 0

