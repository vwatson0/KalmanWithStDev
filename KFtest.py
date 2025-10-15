import numpy as np
import matplotlib.pylab as plt
import KalmanFilterStdEst as KFlib

# tuning test function parameters
Amp = .0001 # amplitude of the arctangent -> depth of the transition
sigNoise = .00001 # standard deviation of the added gaussian noise
stretch = .5 # violence of the transition

# tuning the Kalman filter
cMeas = 1E+1 # covariance of measure -> the highest the more filter the lowest the more reactive to trajectory changes



t = np.linspace(0, 60, 250) # time vector
x = np.arctan(stretch * (t-30)) * Amp # real value
y = x + np.random.randn(len(x)) * sigNoise # measure

KF = KFlib.KFobject(y[0], cMeas) # Initialization of the Kalman Filter object

Xest = np.zeros(len(x)) # for storage Estimated Value
Xest[0] = y[0]
SlopeEst =  np.zeros(len(x)) # for storage estimated slope
StdEst =  np.zeros(len(x)) # for storage estimated noise standard deviation

for k in range(len(y)-1): # going through the time series as if we were measuring live

    KF.EstimateState(y[k+1], t[k+1] - t[k]) # updating the KF with last measure

    # storing values before they get replaced at the next estimate
    Xest[k+1] = KF.X[0]
    SlopeEst[k+1] = KF.X[1]
    StdEst[k+1] = KF.Sig[0]



fig = plt.figure()
ax = fig.add_subplot(311)
ax.plot(t, x)
ax.plot(t, y)
ax.plot(t, Xest)
plt.tick_params('x', labelbottom = False)
plt.ylabel('X')
ax.legend(['true value', 'measure', 'Filtered'])
ax1 = fig.add_subplot(312, sharex = ax)
ax1.plot(t, StdEst)
ax1.plot(t, sigNoise * np.ones(len(t)))
ax1.legend(['Estimated StDev', 'Model stDev'])
plt.tick_params('x', labelbottom = False)
plt.ylabel('Sig')
ax2 = fig.add_subplot(313, sharex = ax)
ax2.plot(t, SlopeEst)
ax2.plot(t[1::],(x[1::]-x[0:(len(x)-1)])/(t[1::]-t[0:(len(t)-1)]))
ax2.legend(['Estimated Slope', 'Model Slope'])
plt.subplots_adjust(hspace =0)
plt.ylabel('Slope')
plt.show()
