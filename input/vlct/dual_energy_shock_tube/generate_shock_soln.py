# this produces a tables of solutions to test against the Sod Shock
# it's a little over-engineered but can be reused to prepare solutions to other
# shock tube test problems

# I have confirmed that the wavespeeds are consistent with the ones used in
# Athena for the Sod Shock Tube test Problem

import numpy as np

def sound_speed(p,rho, gamma):
    return np.sqrt(gamma*p/rho)

class RarefactionFanSoln:
    """ 
    Solution to the Rarefaction fan. 

    K is a variable that should be substituted for L or R. It is specified
    by passing True or False to left
    """
    def __init__(self, gamma, rhoK, pK, uK, p_star, u_star, rho_starK,
                 left = True):
        self.gamma     = float(gamma)
        self.rhoK      = float(rhoK)
        self.pK        = float(pK)
        self.uK        = float(uK)
        self.p_star    = float(p_star)
        self.u_star    = float(u_star)
        self.rho_starK = float(rho_starK)
        self.left = left

        self.aK = sound_speed(self.pK,self.rhoK, self.gamma)
        # eqns 4.54/4.61 from Toro
        self.a_starK = self.aK * (self.p_star/self.pK)**(0.5 - 0.5/self.gamma)

    def head_speed(self):
        # eqns 4.55/4.62 from Toro
        if self.left:
            return self.uK - self.aK
        else:
            return self.uK + self.aK

    def tail_speed(self):
        # eqns 4.55/4.62 from Toro
        if self.left:
            return self.u_star - self.a_starK
        else:
            return self.u_star + self.a_starK

    def _func(self, xdivt):
        signed_aK = self.aK
        if not self.left:
            signed_aK *= -1.

        inv_gp1 = 1./(self.gamma+1)
        gm1 = self.gamma - 1.

        return ( 2.0 * inv_gp1 +
                 (gm1 * inv_gp1 * (self.uK - xdivt) / signed_aK) )**(2. / gm1)

    def rho(self, x, t):
        # adapted from 4.56/4.63 in Toro
        return self.rhoK * self._func(x/t)

    def p(self, x, t):
        # adapted from 4.56/4.63 in Toro
        return self.pK * (self._func(x/t) ** self.gamma)

    def u(self, x, t):
        # adapted from 4.56/4.63 in Toro
        signed_aK = self.aK
        if not self.left:
            signed_aK *= -1.
        return 2./(self.gamma+1) * (signed_aK + 0.5*(self.gamma - 1.)*self.uK
                                    + x/t)

    def dict_of_state_funcs(self):
        return dict(p   = lambda x,t : self.p(x,t),
                    rho = lambda x,t : self.rho(x,t),
                    u   = lambda x,t : self.u(x,t))

def shock_speed(gamma, rhoK, pK, uK, p_star, u_star, rho_starK, left = True):
    """
    Returns the speed of the shock transmission.

    K is a variable that should be substituted for L or R. It is specified
    by passing True or False to left

    Implements equations 4.52 and 4.59 from Toro
    """

    signed_aK = sound_speed(pK, rhoK, gamma)
    if left:
        signed_aK*=-1.

    half_invgamma= 0.5/gamma
    p_ratio = p_star/pK

    return uK + signed_aK * np.sqrt( half_invgamma * (gamma + 1.) * p_ratio +
                                     half_invgamma * (gamma - 1.))

    


class RiemannSoln:
    def __init__(self, gamma, x0,
                 rhoL, pL, uL, rhoR, pR, uR,
                 p_star, u_star, rho_starL, rho_starR):
        self.gamma     = gamma
        self.x0        = x0 # the initial discontinuity

        self.rhoL      = rhoL
        self.pL        = pL
        self.uL        = uL

        self.rhoR      = rhoR
        self.pR        = pR
        self.uR        = uR

        self.p_star    = p_star
        self.u_star    = u_star
        self.rho_starL = rho_starL
        self.rho_starR = rho_starR

        self._wave_speeds = []
        self._solns = []

        self._setup_soln()

    def _setup_soln(self):

        # from left to right:

        self._solns.append(dict(p   = lambda x,t : self.pL,
                                rho = lambda x,t : self.rhoL,
                                u   = lambda x,t : self.uL))
        # first wavespeed
        if self.p_star <= self.pL:
            # there is a Rarefaction Wave to the left
            lfan = RarefactionFanSoln(self.gamma, self.rhoL, self.pL, self.uL,
                                      self.p_star, self.u_star, self.rho_starL,
                                      left = True)
            # append speeds of heads and tails of rarefaction wave
            self._wave_speeds.append(lfan.head_speed())
            self._wave_speeds.append(lfan.tail_speed())
            # append the solution between the head and tail
            self._solns.append(lfan.dict_of_state_funcs())
        else:
            # there is a shock
            S = shock_speed(self.gamma, self.rhoL, self.pL, self.uL,
                            self.p_star, self.u_star, self.rho_starL,
                            left = True)
            self._wave_speeds.append(S)

        # now for the state between the left wave (rarefaction or shock) and
        # the contact wave
        self._solns.append(dict(p   = lambda x,t : self.p_star,
                                rho = lambda x,t : self.rho_starL,
                                u   = lambda x,t : self.u_star))
        # append speed of contact wave
        self._wave_speeds.append(self.u_star)
        # now for the state between contact wave and right wave
        self._solns.append(dict(p   = lambda x,t : self.p_star,
                                rho = lambda x,t : self.rho_starR,
                                u   = lambda x,t : self.u_star))

        # now for the right wave
        if self.p_star <= self.pR:
            # there is a Rarefaction Wave to the left
            lfan = RarefactionFanSoln(self.gamma, self.rhoR, self.pR, self.uR,
                                      self.p_star, self.u_star, self.rho_starR,
                                      left = False)
            # append speeds of heads and tails of rarefaction wave
            # (the order of speeds is inverted from before)
            self._wave_speeds.append(lfan.tail_speed())
            self._wave_speeds.append(lfan.head_speed())
            # append the solution between the head and tail
            self._solns.append(lfan.dict_of_state_funcs())
        else:
            # there is a Shock Wave
            S = shock_speed(self.gamma, self.rhoR, self.pR, self.uR,
                            self.p_star, self.u_star, self.rho_starR,
                            left = False)
            self._wave_speeds.append(S)

        # now for the rightmost state
        self._solns.append(dict(p   = lambda x,t : self.pR,
                                rho = lambda x,t : self.rhoR,
                                u   = lambda x,t : self.uR))

    def get_speeds(self):
        return self._wave_speeds

    def evaluate(self, x_vals, t):

        x = x_vals - x0

        # identify which wavespeed to use

        condlist = [x <= self._wave_speeds[0]*t]
        for i,elem in enumerate(self._wave_speeds[:-1]):
            loc = elem * t
            loc_next = self._wave_speeds[i+1]*t
            condlist.append(np.logical_and(x > loc, x<=loc_next))
        condlist.append(x > self._wave_speeds[-1]*t) 

        out = {}
        for quantity in ['rho','u','p']:
            funclist = [fdict[quantity] for fdict in self._solns]
            out[quantity] = np.piecewise(x,condlist,funclist,t)
        return out

def plot(x, val_dict):
    import matplotlib.pyplot as plt

    fig,ax_arr = plt.subplots(3,1,sharex=True,figsize=(3,8))

    for i,name in enumerate(['rho','u','p']):
        ax = ax_arr[i]
        ax.plot(x,val_dict[name])
        ax.set_xlabel('x')
        ax.set_ylabel(name)
    plt.show()

def write_table(f, time, left_edge, right_edge, x_vals, field_vals):
    f.write('#time = {0!r}\n'.format(       float(time)      ) )
    f.write('#left_edge = {0!r}\n'.format(  float(left_edge) ) )
    f.write('#right_edge = {0!r}\n'.format( float(right_edge)) )

    prepared_dict = {'x'          : x_vals,
                     'density'    : field_vals['rho'],
                     'velocity_x' : field_vals['u'],
                     'pressure'   : field_vals['p']}
    colnames = ['x','density','velocity_x','pressure']
    f.write('{:s}\n'.format(','.join(colnames)))

    cols = [prepared_dict[name].tolist() for name in colnames]
    for row in zip(*cols):
        f.write('{:s}\n'.format(','.join([repr(elem) for elem in row])))
    
if __name__ == '__main__':
    gamma = 1.4

    # 4 characteristic speeds from left to right:
    # head of rarefaction wave, tail of rarefaction wave, contact discontinuity
    # and the shock

    # Happens to the Left of the Head rarefaction wave
    rhoL = 1.0
    pL   = 1.0
    uL   = 0.0

    # Describes to the state to the right of the shock
    rhoR = 0.125
    pR   = 0.1
    uR   = 0.0

    # according to Table 4.3 from Sod:
    # p_star and u_star describe the state everywhere between the tail of the
    # rarefeaction wave and the shock 
    p_star = 0.30313
    u_star = 0.92745
    # rho_starL describes the density to the left of the contact wave (but to
    # the right of the tail of the rarefaction wave)
    rho_starL = 0.42632
    # rho_starR describes the density between contact wave and shock
    rho_starR = 0.26557

    x0 = 0.5
    tfinal = 0.25

    rsoln = RiemannSoln(gamma, x0, rhoL, pL, uL, rhoR, pR, uR,
                        p_star, u_star, rho_starL, rho_starR)
    #print(len(rsoln._wave_speeds))
    #print(len(rsoln._solns))

    N_vals = 128
    x_min = 0.
    x_max = 1.0

    x_vals = np.linspace(x_min,x_max,num=(2*N_vals+1))[1::2]
    result = rsoln.evaluate(x_vals, tfinal)

    with open('sod_shock_tube_t0.25_res128.csv','w') as f:
        write_table(f, tfinal, x_min, x_max, x_vals, result)

    

    

    
