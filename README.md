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

1. Day 1 ( [Link to recording](https://metrumrg.zoom.us/rec/share/axJg0GVfaA9PIac1_L-NDeR3MGk2zWKSb76DPEvydq-odFnX4cPJpKjCKld2vgOS.dbKlnJTeOoVcryRr); Password: `Nw738$$X` )
    * Study designs and confounded exposure-response
    * General theory / background 
    * Binary data
    * Models for binary data
    * Bayesian models for binary data
   
2. Day 2 ( [Link to recording](https://metrumrg.zoom.us/rec/share/nP4CNLd_hH-wTC8ULII3oF9nJBeIOi6G2xNHTmkwKubYgzX_g4Zr_HJXKeaeqhsC.KBM5j5EmX_xa0znd); Password: `*Fi#j2zW` )
    * Time-to-event (TTE) data
        * What makes it different?
             * Describing a distribution w/o censoring (density, CDF)
             * Describing a distribution w/o censoring (hazard, CDF, survival function)
             * Non-parametric estimation of S(t), H(t) and h(t)
             * Study design implications: number of events vs number of subjects
       * Visualizing TTE data vs predictors: K-M plot (session 1)
           * How to interpret it in general?
           * Considerations for exposure metrics in TTE analyses
           * Utility and pitfalls of exposure quartiles
           * Hands-on example: visualizing TTE endpoint vs treatment or exposure
    * TTE semi-parametric modeling
    * Introduction to parametric Survival analysis
   
3. Day 3 ( [Link to recording](https://metrumrg.zoom.us/rec/share/6gd3F7JLlAo5LuUGTvs-6jvXorPxny_l874ms7YFNLYcCqz5m1K5hbo0O_IpfF5C.XlvqcChk3SXhMxwG); Password: `8&p!Fd#t`)
    * TTE parametric modeling (Bayesian)
    * Additional topics in TTE modeling
  
4. Day 4
        * TTE models with a continuously time-varying hazard
        * Repeated time-to-event data
