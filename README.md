# Modeling of binary and time-to-event data using R and Stan

This content is intended for a course which spans four 3-4 hour sessions.

## Learning Objectives

1. Confounded exposure-response
2. Binary data (exploratory analysis and modeling using R and Stan)
3. Time-to-event data (exploratory analysis, semi-parametric and parametric modeling using R and Stan)
4. Advanced topics
    * TTE models with a continuously time-varying hazard
    * Repeated time-to-event data
    * Markov models for analysis of longitudinal discrete state data
    * Joint modeling of continuous and time-to-event data (TD and TTE)
    * Modeling of treatment related AEs (neutropenia, diarrhea, etc)
    


## Agenda

1. Day 1
   a. Study designs and confounded exposure-response
   b. General theory / background 
   c. Binary data
   d. Models for binary data
   e. Bayesian models for binary data
   
2. Day 2
   a. Time-to-event (TTE) data
      i. What makes it different?
         * Describing a distribution w/o censoring (density, CDF)
         * Describing a distribution w/o censoring (hazard, CDF, survival function)
         * Non-parametric estimation of S(t), H(t) and h(t)
         * Study design implications: number of events vs number of subjects
      ii. Visualizing TTE data vs predictors: K-M plot (session 1)
         * How to interpret it in general?
         * Considerations for exposure metrics in TTE analyses
         * Utility and pitfalls of exposure quartiles
         * Hands-on example: visualizing TTE endpoint vs treatment or exposure


   b. TTE semi-parametric modeling
   c. Introduction to parametric Survival analysis
   
3. Day 3
   a. TTE parametric modeling (Bayesian)
   b. Additional topics in TTE modeling
      1. TTE models with a continuously time-varying hazard
      2. Repeated time-to-event data
      3. Markov models for analysis of longitudinal discrete state data
