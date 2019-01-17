import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
# Define some universal params
rho0 = 2.3769e-3
g = 32.26

# Start with some functions
# Density and Temperature
def temp_ft(h):
    # Calculates atmospheric temp under standard atmosphere model
    # imperial units
    if h <= 36089:
        return 288.16 - 1.9812e-3 * h
    else:
        return 288.16 - 1.9812e-3 * 36089

def dens_imp(h):
    # Calculate air density in imperial units
    if h <= 36089:
        return 0.0023769 * (temp_ft(h) / 288.16) ** (-1 - 32.15 / (-1.9812e-3 * 3.0892e3))
    else:
        return dens_imp(36089) * np.exp(-32.15 * (h - 36089) / (3.0892e3 * temp_ft(36089)))

def long_stab(plane, rho, u0, CL1, CD1, Cm1, CDu, CLu, Cmu, CDalpha, CLalpha, CMalpha,CLalphadot, CMalphadot, CLq, CMq):
    '''
    Computes longitudinal stability matrix based on the longitudinal linear model

    :param plane: aircraft object
    :param rho: density at performance altitude
    :param u0: airspeed
    :param CL1: steady-state lift coefficient
    :param CD1: steady state drag coefficient
    :param Cm1: steady-state moment coefficient

    NON-DIMENSIONAL STABILITY DERIVATIVES:
    :param CDu:
    :param CLu:
    :param Cmu:

    :param CDalpha:
    :param CLalpha:
    :param CMalpha:

    :param CLalphadot:
    :param CMalphadot:

    :param CLq:
    :param CMq:
    :return: The longitudinal stability matrix
    '''

    # Additional / Derived non-dimensional stability derivatives
    CXu = -(CDu + 2 * CD1)
    CXalpha = -(CDalpha - CL1)

    CXalphadot = 0
    CXq = 0
    Czu = -(CLu + 2 * CL1)
    Czalpha = -(CLalpha + CD1)
    Czalphadot = -CLalphadot
    Czq = -CLq

    Cw0 = plane.W / (1 / 2 * rho * u0 ** 2 * plane.S)

    # Dimensionalized derivatives
    Xu = rho * u0 * plane.S * Cw0 * np.sin(plane.thta0) + 1 / 2 * rho * u0 * plane.S * CXu
    Xw = 1 / 2 * rho * u0 * plane.S * CXalpha
    Xq = 1 / 4 * rho * u0 * plane.c * plane.S * CXq
    Xwdot = 1 / 4 * rho * plane.c * plane.S * CXalphadot

    Zu = -rho * u0 * plane.S * Cw0 * np.cos(plane.thta0) + 1 / 2 * rho * u0 * plane.S * Czu
    Zw = 1 / 2 * rho * u0 * plane.S * Czalpha
    Zq = 1 / 4 * rho * u0 * plane.c * plane.S * Czq
    Zwdot = 1 / 4 * rho * plane.c * plane.S * Czalphadot

    Mu = 1 / 2 * rho * u0 * plane.c * plane.S * Cmu
    Mw = 1 / 2 * rho * u0 * plane.c * plane.S * CMalpha
    Mq = 1 / 4 * rho * u0 * plane.c ** 2 * plane.S * CMq
    Mwdot = 1 / 4 * rho * plane.c ** 2 * plane.S * CMalphadot

    A1 = [Xu / plane.m, Xw / plane.m, 0, -g * np.cos(plane.thta0)]
    A2 = [Zu / (plane.m - Zwdot), Zw / (plane.m - Zwdot), (Zq + plane.m * u0) / (plane.m - Zwdot), -plane.m * g * np.sin(plane.thta0) / (plane.m - Zwdot)]
    A3 = 1 / plane.Jyy * np.array(
        [Mu + Mwdot * Zw / (plane.m - Zwdot), Mw + Mwdot * Zw / (plane.m - Zwdot), Mq + Mwdot * (Zq + plane.m * u0) / (plane.m - Zwdot),
         -Mwdot * plane.m * g * np.sin(plane.thta0) / (plane.m - Zwdot)])
    A4 = [0, 0, 1, 0]
    A_long = np.array([A1, A2, A3, A4])
    return A_long


def lat_stab(plane, u0, rho, CYbeta, Clbeta, Cnbeta, CYp, Clp, Cnp, CYr, Clr, Cnr):

    scale = 1/2*rho*u0*plane.S*plane.b
    Ixprim = (plane.Jxx*plane.Jzz - plane.Jxz**2)/plane.Jzz
    Izprim = (plane.Jxx*plane.Jzz - plane.Jxz**2)/plane.Jxx
    Izxprim = plane.Jxz/(plane.Jxx*plane.Jzz - plane.Jxz**2)

    scale = 1 / 2 * rho * u0 * plane.S * plane.b # dimensional scaling parameter

    Yv = scale/plane.b*CYbeta
    Yp = 1/2*scale*CYp
    Yr = 1/2*scale*CYr

    Lv = scale*Clbeta
    Lp = 1/2*scale*plane.b*Clp
    Lr = 1/2*scale*plane.b*Clr

    Nv = scale*Cnbeta
    Np = 1/2*scale*plane.b*Cnp
    Nr = 1/2*scale*plane.b*Cnr

    A1 = [Yv/plane.m, Yp/plane.m, (Yr/plane.m - u0), g*np.cos(plane.thta0)]
    A2 = [(Lv/Ixprim + Izxprim*Nv), (Lp/Ixprim + Izxprim*Np),
        (Lr/Ixprim + Izxprim*Nr), 0]
    A3 = [(Izxprim*Lv + Nv/Izprim), (Izxprim*Lp + Np/Izprim),
        (Izxprim*Lr + Nr/Izprim), 0]
    A4 = [0, 1, np.tan(plane.thta0), 0]

    A_lat = np.array([A1, A2, A3, A4])

    return A_lat


class Aircraft:
    def __init__(self, h_cg, h_ac, S, St, lt, c, b, a, at, ae, de_dalpha, Cmac, Clow, Cd0, K, it, e0, W, T_sl, Jxx, Jyy, Jzz, Jxz):
        '''

        Structural parameters
        :param h_cg: chord fraction of CG
        :param h_ac: chord fraction of AC
        :param S: wing area
        :param St: tail area
        :param lt: distance from wing AC to tail AC
        :param c: mean chord length

        Aerodynamic parameters
        :param a: wing lift curve slope
        :param at: tail lift curve slope
        :param ae: elevator effectiveness
        :param de_dalpha: downwash derivative
        :param Cmac: natural wing pitching moment
        :param Clow: wing lift at 0 aoa
        :param Cd0: Constant Drag term
        :param K: Induced drag coefficient
        :param it: tail incidence
        :param e0: initial downwash angle


        :param W: Weight
        :param T_sl: Thrust at sea-level
        '''
        # Set variables
        self.h_cg = h_cg
        self.h_ac = h_ac
        self.S = S
        self.St = St
        self.lt = lt
        self.c = c
        self.b = b
        self.a = a
        self.at = at
        self.ae = ae
        self.de_dalpha = de_dalpha
        self.Cmac = Cmac
        self.Clow = Clow
        self.Cd0 = Cd0
        self.K = K
        self.it = it
        self.e0 = e0
        self.W = W
        self.T_sl = T_sl
        self.Jxx = Jxx
        self.Jyy = Jyy
        self.Jzz = Jzz
        self.Jxz = Jxz

        # Derived variables
        self.L_D_max = math.sqrt(1/(4*self.Cd0*K))
        self.W_S = W/S
        self.Vh = lt / c * St / S
        a_bar = a + at * (1 - de_dalpha) * St / S
        self.h_np= h_ac + at / a_bar * (1 - de_dalpha) * self.Vh
        self.derivs = self.stability_derivs()
        self.m = W/g  # mass
        self.thta0 = 0
    def describe_performance(self,h, V):
        '''
        Describes the performance of the aircraft at a given altitude and velocity
        :param h: altitude
        :param V: Velocity
        :return: Outputs various performance descriptions
        '''
        print('Describing performance at',h,'ft')
        print('Max L/D :',self.L_D_max)
        rho =dens_imp(h)
        T = self.T_sl *  rho/rho0

        z = 1 + math.sqrt(1 + 3/(self.L_D_max**2*T/self.W))
        RC_max = math.sqrt(self.W/self.S*z/(3*rho*self.Cd0)) * (T/self.W)**(3/2) * \
                 (1 - z/6 - 3/(2*(T/self.W)**2 * self.L_D_max**2))
        print('Max Rate of Climb:', RC_max, 'ft/s')

    def neutral_point(self):
        a_bar = self.a + self.at*(1-self.de_dalpha)*self.St/self.S
        return self.h_ac + self.at / a_bar * (1 - self.de_dalpha)*self.Vh

    def req_incidence(self, alpha_target):
        Cm_dalpha = (self.a + self.at * self.St / self.S * (1 - self.de_dalpha)) * (self.h_cg - self.h_ac) \
                    - self.at * self.Vh * (1 - self.de_dalpha);

        return (self.Cmac + self.Clow * (self.h_cg - self.h_ac) + Cm_dalpha * alpha_target) /\
             (self.at * self.Vh * (1 - (self.h_cg - self.h_ac) * self.c / self.lt)) + self.e0

    def stability_derivs(self):
        Cl_0 = self.Clow + self.St / self.S * self.at * (self.it - self.e0)
        Cl_alpha = self.a + self.at * (1 - self.de_dalpha) * self.St / self.S
        Cl_delta = self.St / self.S * self.ae
        Cm_0 = self.Cmac + self.Clow * (self.h_cg - self.h_ac) - self.at * \
                                    (self.it - self.e0) * self.Vh * (1 - (self.h_cg - self.h_ac) * self.c / self.lt)

        Cm_alpha = (self.a + self.at * self.St / self.S * (1 - self.de_dalpha)) * \
                   (self.h_cg - self.h_ac) - self.at * self.Vh * (1 - self.de_dalpha)
        Cm_delta = Cl_delta * (self.h_cg - self.h_ac) - self.ae * self.Vh
        return [Cl_0, Cl_alpha, Cl_delta, Cm_0, Cm_alpha, Cm_delta]

    def static_stability(self, output=False, plot_vals=False):

        derivs = self.stability_derivs()
        Cm_0 = derivs[3]
        Cm_alpha = derivs[4]
        if output:
            print('Cm_0:', Cm_0, '\nCm_alpha:', Cm_alpha)
        if Cm_alpha >= 0:
            print('Aircraft is not statically stable\n')
            return None
        else:
            alpha_trim = -Cm_0 / Cm_alpha
            if output:
                print('Trim alpha = ', alpha_trim * 180 / np.pi, 'degrees\n');

        if plot_vals:
            alpha_span = np.linspace(-10, 10)
            alpha_span = alpha_span * np.pi / 180
            Cm_plot = Cm_0 + Cm_alpha * alpha_span
            plt.plot(alpha_span * 180 / np.pi, Cm_plot)
            plt.grid()
            plt.plot([-10, 10], [0, 0], 'k--')
            plt.plot([0, 0], [np.min(Cm_plot), np.max(Cm_plot)], 'k--')
            plt.xlabel(r'$\alpha$', fontsize=20)
            plt.ylabel('$C_m$', fontsize=20)
            plt.title('Static Stability Analysis', fontsize=20)

        return alpha_trim

    def elevator_analysis(self, deltas, alt=0, plot_vals=False):
        rho = dens_imp(alt)
        trim_alphas = []
        trim_Vs = []
        alpha_span = np.arange(-10, 11, .1) * np.pi / 180
        if plot_vals:
            fig = plt.figure(figsize=(6, 8))

        stability_vals = self.stability_derivs()
        Cl_0 = stability_vals[0]
        Cl_alpha = stability_vals[1]
        Cl_delta = stability_vals[2]
        Cm_0 = stability_vals[3]
        Cm_alpha = stability_vals[4]
        Cm_delta = stability_vals[5]

        for delta in deltas:

            delta = delta * np.pi / 180
            alpha_trim = (-Cm_delta*delta - Cm_0)/Cm_alpha
            trim_alphas.append(alpha_trim)
            Cl = Cl_0 + Cl_alpha*alpha_trim + Cl_delta * delta
            if Cl<=0:
                print(f'Warning, negative lift calculated for delta = {round(delta*180/np.pi,2)} degrees '
                      f'and trim alpha = {round(alpha_trim*180/np.pi,2)} degrees')
                V_trim = 0
            elif alpha_trim*180/np.pi > 20:
                print(f'Warning, trim alpha for delta = {round(delta*180/np.pi,2)} degrees is above stall angle at'
                      f'alpha = {round(alpha_trim*180/np.pi,2)} degrees')
                V_trim = 0
            else:
                V_trim = np.sqrt(self.W/(1/2*rho*Cl*self.S))

            trim_Vs.append(V_trim)
            if plot_vals:
                plt.subplot(2, 1, 1)
                plt.plot(alpha_span * 180 / np.pi, Cm_0 + Cm_alpha * alpha_span + Cm_delta * delta)
                plt.subplot(2, 1, 2)
                plt.plot(alpha_span * 180 / np.pi, Cl_0 + Cl_alpha * alpha_span + Cl_delta * delta)
        if plot_vals:
            plt.subplot(2, 1, 1)
            plt.title('Elevator Analysis', fontsize=20)
            plt.ylabel('$C_m$', fontsize=16)
            plt.grid()

            plt.subplot(2, 1, 2)
            ax = fig.gca()
            plt.ylabel('$C_L$', fontsize=16)
            plt.xlabel(r'$\alpha$', fontsize=20)
            leg_entries = deltas.astype(str)
            for i in range(0, len(leg_entries)):
                leg_entries[i] = leg_entries[i] + r'$^\circ$'

            plt.legend(leg_entries)

            plt.subplot(2, 1, 1)
            plt.plot([0, 0], [-.4, .4], 'k--')
            plt.plot([-10, 10], [0, 0], 'k--')

            plt.subplot(2, 1, 2)
            plt.plot([0, 0], [-.4, 2], 'k--')
            plt.plot([-10, 10], [0, 0], 'k--')
            plt.grid()

        return trim_alphas, trim_Vs

    def launch_sim(self, V0, thta0, h0, plot_vals=True):
        '''
        Models launch as 2D rigid body and integrates equations of motion
        :param V0: initial launch velocity
        :param thta0: initial launch angle relative to ground
        :param h0: launch height
        :return: time and state simulation
        '''

        def launch_dot(state, params):
            '''

            :param t: time
            :param state: 2 position values, 2 velocity values, 2 angular values
            :param params: additional parameters for state transition
            :return:
            '''
            # State Variables
            x = state[0]  # x position, ground-fixed
            z = state[1]  # z position, ground-fixed
            thta = state[2]  # flight path angle
            u = state[3]  # x velocity, body fixed
            w = state[4]  # z velocity, body fixed
            q = state[5]  # angular velocity of flight path angle, body fixed

            # Stability derivatives
            stability_vals = self.stability_derivs()
            Cl_0 = stability_vals[0]
            Cl_alpha = stability_vals[1]
            Cl_delta = stability_vals[2]
            Cm_0 = stability_vals[3]
            Cm_alpha = stability_vals[4]
            Cm_delta = stability_vals[5]

            # other params
            m = params[0]  # mass
            rho = params[1]  # air density
            Cd0 = params[2]  # constant drag term
            K = params[3]  # Induced drag coeff
            Jy = params[4]  # moment of inertia
            delta = params[5]

            # Thrust interpolation
            V = math.sqrt(u ** 2 + w ** 2)
            alpha = math.atan(w / u)  # calculate angle of attack to make sure it does not exceed stall

            # Force Calcs
            Cl_total = Cl_0 + Cl_alpha * alpha + Cl_delta * delta
            Cm_total = Cm_0 + Cm_alpha * alpha + Cm_delta * delta
            Cd_total = Cd0 + K * Cl_total ** 2 + .08

            L = 1 / 2 * rho * u ** 2 * self.S * Cl_total
            M = 1 / 2 * rho * u ** 2 * self.S * self.c * Cm_total
            D = 1 / 2 * rho * u ** 2 * self.S * Cd_total

            # summing aerodynamic force
            X = self.T_sl - D
            Z = -L

            # State derivative
            xdot = u * math.cos(thta) + w * math.sin(thta)
            zdot = -u * math.sin(thta) + w * math.cos(thta)
            thtadot = q
            udot = X / m - g * math.sin(thta) - q * w
            wdot = Z / m + g * math.cos(thta) + q * u
            qdot = M / Jy

            return np.array([xdot, zdot, thtadot, udot, wdot, qdot])

        delta = 0
        # Passing extra parameters to the ODE
        params = [self.W/32.26, dens_imp(0), self.Cd0, self.K, self.Jyy, delta]


        # initial state
        state0 = [0, -h0, thta0, V0, 0, 0]
        # initialize data frame to hold values each loop
        state_sim = pd.DataFrame(columns=['x(t)', 'z(t)', 'theta(t)', 'u(t)', 'w(t)', 'q(t)'])
        state = state0
        state_sim = state_sim.append({'x(t)': state[0], 'z(t)': -state[1], 'theta(t)': state[2] * 180 / math.pi,
                                      'u(t)': state[3], 'w(t)': state[4], 'q(t)': state[5]}, ignore_index=True)
        dt = 0.01
        ts = np.arange(0, 8, dt)
        # RK4 integration scheme
        for i in range(1, len(ts)):
            t = ts[i]
            k1 = launch_dot(state, params)
            k2 = launch_dot(state + .5 * dt * k1, params)
            k3 = launch_dot(state + .5 * dt * k2, params)
            k4 = launch_dot(state + dt * k3, params)

            state = state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            state_sim = state_sim.append({'x(t)': state[0], 'z(t)': -state[1], 'theta(t)': state[2] * 180 / math.pi,
                                          'u(t)': state[3], 'w(t)': state[4], 'q(t)': state[5]}, ignore_index=True)

        v_cruise = math.sqrt(self.W / (1 / 2 * rho0 * self.Clow * self.S))

        # calculate airspeed and alpha as function of time
        V_sim = np.sqrt(state_sim['u(t)'].values ** 2 + state_sim['w(t)'].values ** 2)
        alpha_sim = np.arctan(state_sim['w(t)'] / state_sim['u(t)']) * 180 / math.pi

        # plots
        if(plot_vals):
            plt.figure(figsize=(10, 8))
            plt.subplot(3, 1, 1)
            plt.title('Launch Simulation', fontsize=22)
            plt.plot(ts, state_sim['z(t)'], lw=2)
            plt.plot([0, ts[-1]], [0, 0], 'k--')
            plt.grid()
            plt.ylabel('Height (ft)', fontsize=20)

            plt.subplot(3, 1, 2)
            plt.plot(ts, V_sim)
            plt.plot([0, ts[-1]], [v_cruise, v_cruise], 'k--', label='Cruise Speed')
            plt.legend(fontsize=16)
            plt.grid()
            plt.ylabel('Airspeed (ft/s)', fontsize=20)

            plt.subplot(3, 1, 3)
            plt.plot(ts, alpha_sim)
            plt.ylabel(r'$\alpha$ ($\circ$)', fontsize=20)
            plt.xlabel('time (s)', fontsize=20)
            plt.grid()

        return ts, state_sim