import numpy as np
import matplotlib.pylab as plt
import GenKFlib as KFlib

# tuning test function parameters
Amp = .0001 # amplitude of the arctangent -> depth of the transition
sigNoise = .00001 # standard deviation of the added gaussian noise
stretch = .5 # violence of the transition

# tuning the Kalman filter
cMeas = 1E+1 # covariance of measure -> the highest the more filter the lowest the more reactive to trajectory changes



t = np.linspace(0, 60, 250) # time vector
x1 = np.arctan(stretch * (t-30)) * Amp # real value
x2 = - np.arctan(stretch * (t-30)) * Amp # real value
y = np.asarray([x1 + np.random.randn(len(x1)) * sigNoise, x2 + np.random.randn(len(x2)) * sigNoise]).T # measure

KF = KFlib.KFobject(y[0], cMeas) # Initialization of the Kalman Filter object
print(y[0])
Xest = np.zeros(np.shape(y)) # for storage Estimated Value
Xest[0] = y[0]
SlopeEst =  np.zeros(np.shape(y)) # for storage estimated slope
StdEst =  np.zeros(np.shape(y)) # for storage estimated noise standard deviation

for k in range(len(y)-1): # going through the time series as if we were measuring live

    KF.EstimateState(y[k+1], t[k+1] - t[k]) # updating the KF with last measure

    # storing values before they get replaced at the next estimate
    Xest[k+1] = KF.X[0:len(y[k+1])]
    SlopeEst[k+1] = KF.X[len(y[k+1])::]
    StdEst[k+1] = KF.Sig[0:len(y[k+1])]



fig = plt.figure()
ax = fig.add_subplot(411)
ax.plot(t, x1)
ax.plot(t, y[:, 0])
ax.plot(t, Xest[:,0])
plt.tick_params('x', labelbottom = False)
plt.ylabel('X1')
ax.legend(['true value', 'measure', 'Filtered'])
ax4 = fig.add_subplot(412)
ax4.plot(t, x2)
ax4.plot(t, y[:, 1])
ax4.plot(t, Xest[:,1])
plt.tick_params('x', labelbottom = False)
plt.ylabel('X2')
ax.legend(['true value', 'measure', 'Filtered'])
ax1 = fig.add_subplot(413, sharex = ax)
ax1.plot(t, StdEst)
ax1.plot(t, sigNoise * np.ones(len(t)))
plt.tick_params('x', labelbottom = False)
plt.ylabel('Sigmas')
ax2 = fig.add_subplot(414, sharex = ax)
ax2.plot(t, SlopeEst)
plt.subplots_adjust(hspace =0)
plt.ylabel('Slopes')
plt.show()
