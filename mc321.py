import random
import math

#CONSTANTS
PI = 3.1415926
LIGHTSPEED = 2.997925E10 # in vacuo speed of light [cm/s] */
ALIVE = 1   		# if photon not yet terminated */
DEAD = 0    		# if photon is to be terminated */
THRESHOLD =0.01		# used in roulette */
CHANCE=0.1  		# used in roulette */
COS90D = 0.000001
     # If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
ONE_MINUS_COSZERO=0.000000000001 
     # If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
     # If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */


def InitRandomGen():
    return float(random.uniform(0,1))

     # Initializes the seed for the random number generator. */     
def RandomNum():
    return float(random.uniform(0,1))
     #Calls for a random number from the randum number generator. */

def SIGN(x):
    if x>=0:
        return 1
    else:
        return 0
def main():
    # dummy variables */
    rnd=0        # assigned random value 0-1 */
      #*** INPUT
    #   Input the optical properties
    #   Input the bin and array sizes
    #   Input the number of photons
    #*****/

    mua = 1.0   # cm^-1 */
    mus = 0.0  # cm^-1 */
    g = 0.90
    nt = 1.33
    Nphotons = 500000 # set number of photons in simulation */
    radial_size = 3.05 # cm, total range over which bins extend */
    NR = 100	 # set number of bins.  */
       # IF NR IS ALTERED, THEN USER MUST ALSO ALTER THE ARRAY DECLARATION TO A SIZE = NR + 1. */
    dr = radial_size/NR  # cm */
    albedo = mus/(mus + mua)
    Csph=[]
    Ccyl=[]
    Cpla=[]

    #*** INITIALIZATIONS
    #*****/
    i_photon = 0

    for ir in range(0,NR+1,1):
       Csph.append(0)
       Ccyl.append(0)
       Cpla.append(0)

    #print Ccyl
    #*** RUN
    #   Launch N photons, initializing each one before progation.
    #*****/
    while (i_photon < Nphotons):
        #*** LAUNCH
        #   Initialize photon position and trajectory.
        #   Implements an isotropic point source.
        #*****/
        i_photon += 1	# increment photon count */
        W = 1.0                    # set photon weight to one */
        photon_status = ALIVE      # Launch an ALIVE photon */

        x = 0                      # Set photon position to origin. */
        y = 0
        z = 0

        # Randomly set photon trajectory to yield an isotropic source. */
        costheta = 2.0*RandomNum() - 1.0
        sintheta = math.sqrt(1.0 - costheta*costheta)	# sintheta is always positive */
        psi = 2.0*PI*RandomNum()
        ux = sintheta*math.cos(psi)
        uy = sintheta*math.sin(psi)
        uz = costheta


        # HOP_DROP_SPIN_CHECK
        #   Propagate one photon until it dies as determined by ROULETTE.
        #*******/

        while (photon_status == ALIVE):


            while(rnd <= 0.0):
                rnd = RandomNum()  # yields 0 < rnd <= 1 */
            #*** HOP
            #   Take step to new position
            #   s = stepsize
            #   ux, uy, uz are cosines of current photon trajectory
            #*****/

            s = -math.log(rnd,math.e)/(mua + mus)          # Step size.  Note: log() is base e */
            x += s * ux                        # Update positions. */
            y += s * uy
            z += s * uz


            #*** DROP
            #   Drop photon weight (W) into local bin.
            #*****/
            absorb = W*(1 - albedo)      # photon weight absorbed at this step */
            W -= absorb                  # decrement WEIGHT by amount absorbed */

            # spherical */
            r = math.sqrt(x*x + y*y + z*z)    # current spherical radial position */
            ir = int(r/dr)           # ir = index to spatial bin */
            if (ir >= NR): ir = NR        # last bin is for overflow */

            Csph[ir] += absorb           # DROP absorbed weight into bin */

            # cylindrical */
            r = math.sqrt(x*x + y*y)          # current cylindrical radial position */
            ir = r/dr       # ir = index to spatial bin */
            if (ir >= NR): ir = NR        # last bin is for overflow */
            Ccyl[int(ir)] += absorb           # DROP absorbed weight into bin */

            # planar */
            r = math.fabs(z)                  # current planar radial position */
            ir = r/dr           # ir = index to spatial bin */
            if (ir >= NR): ir = NR        # last bin is for overflow */
            Cpla[int(ir)] += absorb           # DROP absorbed weight into bin */


            #*** SPIN
            #   Scatter photon into new trajectory defined by theta and psi.
            #   Theta is specified by cos(theta), which is determined
            #   based on the Henyey-Greenstein scattering function.
            #   Convert theta and psi into cosines ux, uy, uz.
            #*****/
            # Sample for costheta */
            rnd = RandomNum()
            if (g == 0.0):
                costheta = 2.0*rnd - 1.0
            else:
                temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd)
                costheta = (1.0 + g*g - temp*temp)/(2.0*g)

            sintheta = math.sqrt(1.0 - costheta*costheta) # sqrt() is faster than sin(). */

                # Sample psi. */
            psi = 2.0*PI*RandomNum()
            cospsi = math.cos(psi)
            if (psi < PI):
                sinpsi = math.sqrt(1.0 - cospsi*cospsi)     # sqrt() is faster than sin(). */
            else:
                sinpsi = -math.sqrt(1.0 - cospsi*cospsi)

                # New trajectory. */
            if (1 - math.fabs(uz) <= ONE_MINUS_COSZERO):      # close to perpendicular. */
                uxx = sintheta * cospsi
                uyy = sintheta * sinpsi
                uzz = costheta * SIGN(uz)   # SIGN() is faster than division. */

            else:				# usually use this option */
                temp = math.sqrt(1.0 - uz * uz)
                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
                uzz = -sintheta * cospsi * temp + uz * costheta


            # Update trajectory */
            ux = uxx
            uy = uyy
            uz = uzz


            #*** CHECK ROULETTE
            #If photon weight below THRESHOLD, then terminate photon using Roulette technique.
            #Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
            #and 1-CHANCE probability of terminating.
            #*****/
            if (W < THRESHOLD):
                if (RandomNum() <= CHANCE):
                    W /= CHANCE
                else: photon_status = DEAD
            # end STEP_CHECK_HOP_SPIN */
            # If photon dead, then launch new photon. */
            # end RUN */



        #*** SAVE
        #Convert data to relative fluence rate [cm^-2] and save to file called "mcmin321.out".
    #*****/
    target = open("mc321.txt", "w")

    # print header */
    target.writelines("%s %s \n" % ("number of photons", Nphotons))
    #print Nphotons

    target.writelines("%s %s \n" %  ("bin size = ",str(dr) + "[cm]" ))
    target.writelines("%s \n" % ("last row is overflow. Ignore.",))

    # print column titles */
    target.writelines("%s \n" % ("r [cm] \t Fsph [1/cm2] \t Fcyl [1/cm2] \t Fpla [1/cm2]\n",))

    # print data:  radial position, fluence rates for 3D, 2D, 1D geometries */
    for ir in range(0,NR,1) :
        # r = sqrt(1.0/3 - (ir+1) + (ir+1)*(ir+1))*dr; */
        r = (ir + 0.5)*dr
        shellvolume = 4.0*PI*r*r*dr # per spherical shell */
        Fsph = Csph[ir]/Nphotons/shellvolume/mua
        shellvolume = 2.0*PI*r*dr   # per cm length of cylinder */
        Fcyl = Ccyl[ir]/Nphotons/shellvolume/mua
        shellvolume = dr            # per cm2 area of plane */
        Fpla =Cpla[ir]/Nphotons/shellvolume/mua
        target.writelines("%s %s %s %s \n" % ( r, Fsph, Fcyl, Fpla))


    target.close()

    import matplotlib.pyplot as plt

    lines = open('mc321.txt').read().split("\n")
    x = []
    y = [[], [], []]

    for i in range(5, len(lines) - 1, 1):
        line = lines[i].split()
        #print line
        x.append(line[0])
        y[0].append(line[1])
        y[1].append(line[2])
        y[2].append(line[3])
    labels = ["Fesf", "Fcil", "Fpla"]

    for i in range(0, 3, 1):
        plt.scatter(x, y[i], s=7)
        plt.plot(x, y[i], label=labels[i])

    plt.title("MC321")
    plt.xlabel("r(cm)")
    plt.legend()
    plt.ylabel("Ritmo de fluencia(W/m2)")
    plt.yscale("log")
    plt.show()


if __name__ == "__main__":
    main()


