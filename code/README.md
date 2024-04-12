# PANDA RCT Go/NoGo Task Data Analysis

This project encompasses MATLAB scripts for modeling the task and Python notebooks for analyzing the randomized control trial (RCT) data and parameter estimates.

## matlab Folder: Modeling the Go/NoGo Task
- **emfit/**: A set of Matlab scripts which performs fixed-effects empirical Bayes batch EM inference by first fitting priors for a set of models from each model class, plotting the inferred posterior parameters, performing approximate model comparison, generating surrogate data and performing some visual comparisons of the true and surrogate data.
- **emfit/mAffectiveGoNogo/**: Contains different models fitted to the Go/NoGo task data.
- **powersimulation/**: Code for power simulation, a crucial step for preregistration.
- **parameter_recovery/**: Code for running parameter recovery.

## python Folder: RCT Data and Parameter Estimates Analysis

The Python folder comprises Jupyter notebooks, each tailored for a specific task. Here's an overview of the notebooks:

### Notebooks Overview

- **'prep_trial_data.ipynb'**: Load RCT data and select variables for further analysis.
- **'basic_task_characteristics.ipynb'**: Explore dataset characteristics and task performance.
- **'parameter_investigation.ipynb'**: Examine correlation, test-retest reliability of parameter estimates, recoverability of parameters, and their relation to basic task characteristics.
- **'main_analyses.ipynb'**: Conduct the main analysis investigating non-informative data, sertraline effects, and associations with anxiety and depressive symptoms.
- **'plots_for_paper.ipynb'**: Code for generating plots shown in the paper.
- **'addressing_reviewer_comments_additional_variables.ipynb'**: Investigate variables from additional questionnaires, including medication and other interventions based on reviewers' suggestions.
- **'addressing_reviewer_comments_additional_analyses.ipynb'**: Run additional analyses investigating subsamples and controlling for suggested factors.

## Prerequisites

Ensure you have Conda installed. You can create the necessary environment using the provided `environment.yml` file.

## Data Access

Please note that the data is not publicly available.
