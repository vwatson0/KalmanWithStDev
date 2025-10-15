
import numpy as np

"This is a Kalman filter tool including the estimate of the instant standard deviation of noise at any time." \
"This can be used to track the statistics of any time series under the assumption of smooth non linearities " \
"and gaussian noise." \
"The object KFobject is initialized with the First measure and a large value for the measure covariance cMeas" \
"The object is updated using the method Estimate state with arguments the current measure and the elapsed time " \
"since the last one." \
"At any time the last etimate of the state KFobject.X and the PM (posterior mean) estimate of the noise standard deviation" \
" KFobject.Sig are accessible with: " \
"- Estimates posterior mean of the state = KFobject.X[0]" \
"- Estimated current dynamic of the state dx/dt = KFobject.X[1]" \
"- Estimated posterior mean of the noise standard deviation = KFobject.Sig[0]" \
"- Estimated current dynamic of the noise standard deviation = KFobject.Sig[1]"

class KFobject:
    def __init__(self, FirstMeasure,  cMeas, cLin=10E-3) :
        ''' Initialize:
        Initialize with a large value for cMeas (10E+3 for example), if you need to filter more noise increase cMeas,
        if you see the filter lagging in tracking your time series decrease cMeas.
        cLin is relative to cMeas no need to
        The first measure will initialize X[0]
        '''
        self.X = np.zeros(2)  # X[0] contains the current estimate of the Beam mean X[1] is the estimate of the slope
        self.X[0] = FirstMeasure
        self.PX = np.zeros([2, 2])  # Covariance Matrix of X (Calculated by the filter)
        self.Sig = np.zeros(2)  # Sig[0] contains the current estimate of Beam std
        self.PS = np.zeros([2, 2])  # Covariance Matrix of Sig (Calculated by the filter)
        self.Q = np.ones(
            [2, 2]) * cLin  # relative confidence in the linear dynamic if increased less noise slower convergence
        self.R = np.array(
            [[cMeas]])  # relative confidence in the measurement if increased faster convergence more noise
        self.F = np.array([[1., 1.], [0, 1.]])  # self dynamic of the system
        # [[link measure to state estimate, link fist order to state estimate],[link state to 1st order, propagate 1st order]]
        self.H = np.array([1.,
                           0])  # link from the state space to the measures space (here the transformation from the measure to X[0] is 1)
        # and we do not measure directly dx/dt but deduce it with F

    def EstimateState(self, measure, deltaT) :
        # extracting current values from the filter object
        "This updates the Kalman filter object with the next measure and the elapsed time since the previous one"
        PoldX = self.PX
        PoldSig = self.PS
        oldx = self.X
        oldSig = self.Sig

        F = self.F
        Q = self.Q
        R = self.R
        H = self.H
        F[0, 1] = deltaT  # updating F

        # predictions
        xPred = np.dot(F,
                       oldx)  # predicting the state xPred[k,0] = xEst[k-1,0] + (dx/dt)_estimate | xPred[k,1] = (dx/dt)_est
        pPred = np.dot(np.dot(F, PoldX),
                       F.T) + Q  # Covariance matrix of the prediction (the bigger, the less confident we are)

        SigPred = np.dot(F, oldSig)  # Same thing but with standard deviation of the beam current
        SigpPred = np.dot(np.dot(F, PoldSig), F.T) + Q

        Inow = measure

        # updates

        y = Inow - np.dot(H,
                          xPred)  # Calculating the innovation (diff between measure and prediction in the measure space)
        S = np.dot(np.dot(H, pPred),
                   H.T) + R  # Calculating the Covariance of the measure (the bigger the less confident in the measure)

        K = np.dot(np.array([np.dot(PoldX, H.T)]).T, np.linalg.inv(S))  # Setting the Kalman optimal gain
        newX = xPred + np.dot(K, np.atleast_1d(y)).T  # Estimating the state at this instant
        PnewX = np.dot((np.eye(len(pPred)) - np.dot(K, np.array([H]))), pPred)  # Covariance matrix of the state

        # same steps followed for the standard deviation
        y = np.sqrt((Inow - newX[0]) ** 2) - np.dot(H, SigPred)  # Innovation of the standard deviation
        # this is an additional drawer to the Kalman filter and it is rather uncommon to estimate another variable that is
        # statistically dependent there may be better solutions, but this one works
        S = np.dot(np.dot(H, SigpPred), H.T) + R
        K = np.dot(np.array([np.dot(PoldSig, H.T)]).T, np.linalg.inv(S))
        newSig = (SigPred + np.dot(K, np.atleast_1d(y)).T)
        if newSig[0] < 0 :  # VW mod Aug9 floor sigma
            newSig[0] = 0
        PnewSig = np.dot((np.eye(len(pPred)) - np.dot(K, np.array([H]))), SigpPred)

        # Updating the Kalman filter object
        self.PX = PnewX
        self.X = newX
        self.Sig = newSig
        self.PS = PnewSig
