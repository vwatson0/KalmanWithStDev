# KalmanWithStDev
This is a Kalamn Filter Based tolld developped to track and filter times series related to ECR ion sources. I recursively calculate the posterior mean estimator of the state but also estimate the standard deviation


The goal of this tool, is to track the statistics of a time series to be able to have at each instant the best estimate of it's mean standard deviation and 1st order dynamic.

The posterior mean (PM) estimate of the state gives the estimate of the beam current value, but also the slope that it is following.
There is a second layer added to the filter to recursively estimate the standard deviation of the time series using Bayesian statistics.

Initialization - the object is created at the begining of the monitoring with:
- import KalmanFilterStdEst as KFlib
- KF = KFlib.KFobject(y0, cMeas)
  where y0 is the first measure and cMeas is a scalarused to build the covariance matrix of the measurements.
  the more noisy the data is the highest we need cMeas.
  If you realize that the noise is not enough filtered -> increase cMeas 
  If you realize the filter takes for ever to get back on track after a transition -> decrease cMeas

Update - At any new measurement update the Filtwer with:

- KF.EstimateState(y[k], t[k] - t[k-1])
  where y[k] is the new measure taken at t = t[k] and t[k-1] is the time when the Kalman Filter got updated last

Outputs - To access the last state of the Filter:

Current PM estimate of the time series value = KF.X[0]
Current estimate of the slope followed by the time series = KF.X[1]
Current estimated value of the standard deviation = KF.Sig[0]

These estimate necessitate some time for the filter to converge. When a change is detected or ordered between two stable values, wait for the slope KF.X[1] to increase and and then decrease so that the slope estimate is almost null. Then you can grab the estimate X[0] and Sig[0] as instant estimate of the time series statistics.

If the slope KF.X[1] is high on purpose, the estimate of the standard deviation will be affected and become non-reliable. The average value estimate should be reliable as long as the slope has been constant (or almost) for the convergence to happen.
