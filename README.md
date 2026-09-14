# drug-toxicity-qsar-prediction

Internship Project: Predicting Drug Toxicity Endpoints using QSAR Modelling. Multi-stage ML pipeline for toxicity prediction with mathematical modeling and parallel computing optimization.

## High Level Workflow

![High Level Workflow](./assets/high-level-workflow.png)

## Overview

This project implements a comprehensive three-stage pipeline for predicting drug toxicity endpoints using Quantitative Structure-Activity Relationship (QSAR) modeling:

### Stage 1: Model Training & Classification
- Dataset acquisition from Therapeutic Data Commons
- Data preprocessing (standardization, handling missing values)
- Classifier selection (SVC, Random Forest)
- Cross-validation and hyperparameter tuning
- Optional hybrid classifier approach

### Stage 2: QSAR Mathematical Modeling
- QSAR mathematical model development
- Lead optimization
- Pharmacokinetics exposure data analysis
- Biological activity and lipophilicity prediction
- Bioavailability assessment using Dorst's Law
- Pfizer 3/75 Rule validation

### Stage 3: Parallel Computing Architecture
- Implementation using OMP, SMP, and MSI parallel computing frameworks
- Optimization for high-performance computation

## Installation

```bash
pip install -r requirements.txt
```

## Requirements

- `setuptools<=81.0.0`
- `pytdc-nextml`

## Contributing

This is an internship project. Please reach out for contribution guidelines.
