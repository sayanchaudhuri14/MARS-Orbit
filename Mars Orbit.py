import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from math import sin, cos, tan, sqrt, pi
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.optimize import basinhopping
from scipy.optimize import brute
import sympy


if __name__ == "__main__":

    # Import oppositions data from the CSV file provided

    mars = pd.read_csv(
        r'C:\Users\sayan\OneDrive\Desktop\Data Analysis Projects\Data\01_data_mars_opposition_updated.csv'
    )

    # Extract times from the data in terms of number of days.
    # "times" is a numpy array of length 12. The first time is the reference
    # time and is taken to be "zero". That is times[0] = 0.0
    times = [0]
    for i in range(1, 12):
        times.append(
            (dt.datetime(mars['Year'][i], mars['Month'][i], mars['Day'][i],
                         mars['Hour'][i], mars['Minute'][i]) -
             dt.datetime(mars['Year'][i - 1], mars['Month'][i - 1],
                         mars['Day'][i - 1], mars['Hour'][i - 1],
                         mars['Minute'][i - 1])).total_seconds() / (24 * 3600))
    mars['times'] = times
    mars['times'] = mars['times'].cumsum()
    times = np.array(mars['times'])
    assert len(times) == 12, "times array is not of length 12"

    # Extract angles from the data in degrees. "oppositions" is
    # a numpy array of length 12.
    mars['oppositions'] = (mars["ZodiacIndex"]) * 30 + mars[
        'Degree'] + mars['Minute.1'] / 60 + mars['Second'] / 3600  # degrees
    oppositions = np.array(mars['oppositions'])
    assert len(oppositions) == 12, "oppositions array is not of length 12"

    # Call the top level function for optimization
    ###########################################################################################################################
    ###########################################################################################################################
    # Sridharacharya Formula


    def quad_root(a, b, c):
        root1 = (-b + sqrt(b**2 - 4 * a * c)) / (2 * a)
        root2 = (-b - sqrt(b**2 - 4 * a * c)) / (2 * a)

        return [root1, root2]

    ###########################################################################################################################
    ###########################################################################################################################
    def Plotting(c, e1, e2, z, r, s, times, oppositions):

        Ex = e1 * cos(e2)  # Ex and Ey are Equant locations with Sun at Origin
        Ey = e1 * sin(e2)
        Cx = cos(c)
        Cy = sin(c)
        Z = (z + s * times) % (
            2 * pi
        )  # Z is the angle made by line joining Equant and Mars with Aries
        calcXs = [Ex]
        calcYs = [Ey]
        obsXs = [0]
        obsYs = [0]

        plt.figure()
        circle = plt.Circle((Cx, Cy), r, color='red', fill=False, ec="black")
        for i in range(len(Z)):
            K_line = Ey - tan(Z[i]) * Ex
            a = (1 + tan(Z[i])**2)
            b = -2 * Cx + 2 * tan(Z[i]) * (Ey - tan(Z[i]) * Ex - Cy)
            c = Cx**2 + (Ey - tan(Z[i]) * Ex - Cy)**2 - r**2

            if b**2 < 4 * a * c:
                break
            Root = quad_root(a, b, c)

            if Z[i] < 3 * pi / 2 and Z[i] > pi / 2:
                Rx = Root[1]
            else:
                Rx = Root[0]

            Ry = tan(Z[i]) * Rx + K_line

            Ry = float(
                Ry
            )  # if I were to use the sympy approach, this step wud be crucial to determine arctan2(Ry,Rx)
            Rx = float(Rx)
            calcXs.append(Rx)
            calcYs.append(Ry)

        for i in range(len(oppositions)):
            a = 1 + tan(np.radians(oppositions[i]))**2
            b = -2 * Cx - 2 * Cy * tan(np.radians(oppositions[i]))
            c = Cx**2 + Cy**2 - r**2

            Root = quad_root(a, b, c)

            if np.radians(oppositions[i]) < 3 * pi / 2 and np.radians(
                    oppositions[i]) > pi / 2:
                Rx = Root[1]
            else:
                Rx = Root[0]

            Ry = tan(np.radians(oppositions[i])) * Rx

            obsXs.append(Rx)
            obsYs.append(Ry)

        ax = plt.gca()
        for i in range(1, len(calcXs)):
            plt.plot([calcXs[i], calcXs[0]], [calcYs[i], calcYs[0]],
                     'o--',
                     markersize=5,
                     mfc='black',
                     mec='black')
            plt.plot([obsXs[i], obsXs[0]], [obsYs[i], obsYs[0]],
                     'o-',
                     markersize=5,
                     mfc='black',
                     mec='black')

        ax.add_patch(circle)
        plt.show()
        return None

    ###########################################################################################################################
    ############################################***ANSWER 1***#################################################################
    ###########################################################################################################################


    def MarsEquantModel(
            Variables, r, s, times,
            oppositions):  # here the angles are in radians except oppositions

        c, e1, e2, z = Variables[0], Variables[1], Variables[2], Variables[
            3]  # Unboxing the Variables

        Ex = e1 * cos(e2)  # Ex and Ey are Equant locations with Sun at Origin
        Ey = e1 * sin(e2)
        Cx = cos(c)
        Cy = sin(c)  # Cx and Cy are center of the Mars's orbit
        error = []
        Z = (z + s * times) % (
            2 * pi
        )  # Z is the angle made by line joining Equant and Mars with Aries

        ###########################################################################################################################

        #  x = sympy.symbols('x')             #Tried using sympy library. The function was evaluating but the minimize function was
        # taking too long to optimize. Also, the number of iterations required were higher.
        ###########################################################################################################################
        # Finding the intersection of the circle described the center C and Radius 'r' and the line joining Equant and Mars

        for i in range(len(Z)):
            K_line = Ey - tan(
                Z[i]
            ) * Ex  # Constant of the line joining equant and the orbit if mars
            ###########################################################################################################################
            #         y = tan(Z[i]) * x + K_line                            #sympy approach. Can evaluate them by commenting in the following
            # four lines and commenting out the next 3 lines
            #         expr = (x-Cx)**2 + (y-Cy)**2 - r**2
            #         smpl = sympy.simplify(expr)
            #         (a,b,c) = (smpl.coeff(x,2),smpl.coeff(x,1),smpl.coeff(x,0))
            ###########################################################################################################################
            a = (1 + tan(Z[i])**2)
            b = -2 * Cx + 2 * tan(Z[i]) * (Ey - tan(Z[i]) * Ex - Cy)
            c = Cx**2 + (Ey - tan(Z[i]) * Ex - Cy)**2 - r**2

            if b**2 < 4 * a * c:
                break
            Root = quad_root(a, b, c)

            if 3 * pi / 2 > Z[i] > pi / 2:
                Rx = Root[1]
            else:
                Rx = Root[0]

            Ry = tan(Z[i]) * Rx + K_line

            Ry = float(
                Ry
            )  # if I were to use the sympy approach, this step wud be crucial to determine arctan2(Ry,Rx)
            Rx = float(Rx)
            Calc_long = (
                (np.arctan2(float(Ry), float(Rx)) + 4 * pi) % (2 * pi)
            ) * 180 / pi  # Calculated Heliocentric Longitude at oppositions

            error.append(
                Calc_long -
                oppositions[i])  # Oppositions Discrepancy Measured in Degrees

        error = np.array(error)
        error = (error + 720) % 360

        for x in range(len(error)):

            if error[x] > 180:
                error[x] = error[x] - 360
            elif error[x] < -180:
                error[x] = error[x] + 360

        if len(error) == 0:
            maxError = 1000
        else:
            maxError = max(np.abs(np.array(error)))

        return error, maxError

    ###########################################################################################################################
    ###########################################################################################################################
    # Sub_routines

    def bestc(c, e1, e2, z, r, s, times, oppositions):
        maxError = MAXERROR(c, e1, e2, z, r, s, times, oppositions)
        return maxError

    def beste1(e1, c, e2, z, r, s, times, oppositions):
        maxError = MAXERROR(c, e1, e2, z, r, s, times, oppositions)
        return maxError

    def beste2(e2, c, e1, z, r, s, times, oppositions):
        maxError = MAXERROR(c, e1, e2, z, r, s, times, oppositions)
        return maxError

    def bestz(z, c, e1, e2, r, s, times, oppositions):
        maxError = MAXERROR(c, e1, e2, z, r, s, times, oppositions)
        return maxError

    def MAXERROR(c, e1, e2, z, r, s, times, oppositions):
        error, maxError = MarsEquantModel((c, e1, e2, z), r, s, times,
                                          oppositions)
        return maxError

    def MAXERROR_net(Vars, r, s, times, oppositions):
        (c, e1, e2, z) = Vars
        error, maxError = MarsEquantModel((c, e1, e2, z), r, s, times,
                                          oppositions)

        return maxError

    ###########################################################################################################################
    ############################################***ANSWER 2***#################################################################
    ###########################################################################################################################
    def bestOrbitInnerParams(r, s, times, oppositions):
        c = np.radians(100)
        z = np.radians(60)
        e2 = np.radians(70)
        e1 = 1

        # TRIALS      #This was my other option of the order in which I can optimize the variables.
        # The orders were obtained after trying all fact(4)=24 combinations.

        #     res=minimize_scalar(beste1,e1,args=(c,e2,z,r,s,times,oppositions),method='Bounded',options={'disp': True,'maxiter': 20000},bounds=(0.9,1.7))
        #     e1=res.x

        #     res=minimize_scalar(beste2,e2,args=(c,e1,z,r,s,times,oppositions),method='Bounded',options={'disp': True,'maxiter': 20000},bounds=(-pi,pi))
        #     e2=res.x
        #     res=minimize_scalar(bestz,z,args=(c,e1,e2,r,s,times,oppositions),method='Bounded',options={'disp': True,'maxiter': 20000},bounds=(-pi,pi))
        #     z=res.x
        #     res=minimize_scalar(bestc,c,args=(e1,e2,z,r,s,times,oppositions),method='Bounded',options={'disp': True,'maxiter': 20000},bounds=(-pi,pi))
        #     c=res.x

        # BEST SET    #This order was obtained after trying all fact(4)=24 combinations.

        res = minimize_scalar(beste1,
                              e1,
                              args=(c, e2, z, r, s, times, oppositions),
                              method='Bounded',
                              options={
                                  'disp': False,
                                  'maxiter': 20000
                              },
                              bounds=(0.9, 1.7))
        e1 = res.x

        res = minimize_scalar(bestz,
                              z,
                              args=(c, e1, e2, r, s, times, oppositions),
                              method='Bounded',
                              options={
                                  'disp': False,
                                  'maxiter': 20000
                              },
                              bounds=(-pi, pi))
        z = res.x

        res = minimize_scalar(bestc,
                              c,
                              args=(e1, e2, z, r, s, times, oppositions),
                              method='Bounded',
                              options={
                                  'disp': False,
                                  'maxiter': 20000
                              },
                              bounds=(-pi, pi))
        c = res.x

        res = minimize_scalar(beste2,
                              e2,
                              args=(c, e1, z, r, s, times, oppositions),
                              method='Bounded',
                              options={
                                  'disp': False,
                                  'maxiter': 20000
                              },
                              bounds=(-pi, pi))
        e2 = res.x

        res = minimize(MAXERROR_net, (c, e1, e2, z),
                       args=(r, s, times, oppositions),
                       method='Nelder-Mead',
                       options={
                           'disp': False,
                           'maxiter': 20000
                       },
                       bounds=((-pi, pi), (0.9, 1.7), (-pi, pi), (-pi, pi)))
        c, e1, e2, z = res.x[0], res.x[1], res.x[2], res.x[3]
        variables = [c, e1, e2, z]

        error, maxError = MarsEquantModel(variables, r, s, times, oppositions)

        return c, e1, e2, z, error, maxError

    ###########################################################################################################################
    ##########################################***ANSWER 3***###################################################################
    ###########################################################################################################################
    def bestS(r, times, oppositions):

        S = np.arange((1 - 0.002) * 2 * pi / 687, (1 + 0.002) * 2 * pi / 687,
                      pi / 10000000)
        maxError = 10000
        error = []
        for test_s in range(len(S)):
            c, e1, e2, z, err, maxErr = bestOrbitInnerParams(
                r, S[test_s], times, oppositions)

            if maxError > maxErr:
                maxError = maxErr
                error = err
                best_s = S[test_s]

        return best_s, error, maxError

    ###########################################################################################################################
    ############################################***ANSWER 4***#################################################################
    ###########################################################################################################################
    def bestR(s, times, oppositions):

        R = np.arange(7.5, 9, 0.005)
        error = []
        maxError = 10000
        for test_r in range(len(R)):

            c, e1, e2, z, err, maxErr = bestOrbitInnerParams(
                R[test_r], s, times, oppositions)

            if maxError > maxErr:
                maxError = maxErr
                error = err
                best_r = R[test_r]

        return best_r, error, maxError

    ###########################################################################################################################
    #############################################***ANSWER 5***################################################################
    ###########################################################################################################################
    def bestMarsOrbitParams(times, oppositions):

        s_init = 2 * pi / 687

        BEST_S = s_init
        BEST_R = 10
        best_s = 2 * pi / 687
        maxError = 10000
        error = [
        ]  # BruteForcing my way into the perfect combination of R ans S

        for r_test in np.arange(8, 9, 0.005):
            for s_test in np.arange((1 - 0.002) * 2 * pi / 687,
                                    (1 + 0.002) * 2 * pi / 687, pi / 10000000):

                c, e1, e2, z, err, maxErr = bestOrbitInnerParams(
                    r_test, s_test, times, oppositions)
                # print(r_test, s_test, maxError, maxErr)
                if maxErr < maxError:
                    maxError = maxErr
                    error = err
                    goodVar = (c, e1, e2, z)
                    BEST_R = r_test
                    BEST_S = s_test

        (c, e1, e2, z) = (goodVar[0], goodVar[1], goodVar[2], goodVar[3])

        maxError = MAXERROR(c, e1, e2, z, BEST_R, BEST_S, times, oppositions)
        Plotting(c, e1, e2, z, BEST_R, BEST_S, times, oppositions)
        (c, e1, e2, z) = (goodVar[0]*180/pi, goodVar[1], goodVar[2]*180/pi, goodVar[3]*180/pi)

        return BEST_R, BEST_S, c, e1, e2, z, error, maxError

    # The angles are all in degrees
    r, s, c, e1, e2, z, errors, maxError = bestMarsOrbitParams(
        times, oppositions)

    assert max(list(map(
        abs, errors))) == maxError, "maxError is not computed properly!"
    print(
        "Fit parameters: r = {:.4f}, s = {:.4f}, c = {:.4f}, e1 = {:.4f}, e2 = {:.4f}, z = {:.4f}"
        .format(r, s, c, e1, e2, z))
    print("The maximum angular error = {:2.4f}".format(maxError))
