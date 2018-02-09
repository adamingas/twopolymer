from __future__ import division

import sys
import os
import math
import linecache
import config
import argparse
import numpy as np

"""
Version 3.0.0Alpha
Changed to two polymer chains. Simulate run for a number of constants and find maximum separation squared.
"""

def walk(k,max_file,chi_element):
    """
    Performs the random walk.
    There are 4 possibilities of a walk, which are combinations of the binary states of hydrodynamic
    interactions and thermal noise.
    k is the weissenberg number
    """
    # No hydrodynamic interaction and no noise

    # Hydrodynamic interactions but not noise
    #if str.upper(noise) in ("NO","N","NO NOISE"):

        #for j in xrange(runs):
        #     out = open("Run{}_Wi{}_chi{}".format(j, k,chi_element), "w")
        #     z = zinitial
        #     x = xinitial
        #     y = yinitial
        #     rseparation = x * x + y * y + z*z
        #     out.write("{} {} {} {}\n".format(0, x, y, z))
        #     maxr = rseparation
        #     tmax = 0
        #     for i in xrange(steps):
        #         rseparation = x * x + y * y +z*z
        #         ynew = y + time_step*(1/(1+3/(4*math.sqrt(rseparation)*ar_ratio)))*(-y/(1-rseparation) +
        #             (1/(ar_ratio*math.sqrt(rseparation)*(4/3) + 2))*(-y*y*x*k/rseparation +
        #             (y*y*y + y*x*x + y*z*z)/((1-rseparation)*rseparation)))
        #
        #         xnew = x + time_step * (1 / (1 + 3 / (4 * math.sqrt(rseparation) * ar_ratio))) * (k*y
        #         -x / (1 - rseparation) +(1 / (ar_ratio * math.sqrt(rseparation)*(4/3) + 2)) * (-y * x * x * k / rseparation
        #         + (x * x * x + x * y * y + x*z*z) / ((1 - rseparation) * rseparation)))
        #
        #         znew = z + time_step * (1 / (1 + 3 / (4 * math.sqrt(rseparation) * ar_ratio))) * (
        #                 -z / (1 - rseparation) +(1 / (ar_ratio * math.sqrt(rseparation) * (4 / 3) + 2)) * (-z*y*x*k/rseparation +
        #                 (z*y*y + z*x*x + z*z*z) / ((1 - rseparation) * rseparation)))
        #         x, y, z = xnew, ynew, znew
        #         # Finding out largest separation squared
        #         if maxr <= x * x + y * y + z*z and args.max:
        #             maxr = x * x + y * y +z*z
        #             tmax = (i+1) * time_step
        #         out.write("{} {} {} {}\n".format(time_step * (i + 1), x, y, z))
        #     if args.max:
        #         max_file.write("{} {} {}\n".format(k, maxr / (init_separation ** 2), tmax))
        #     update_progress(j / (runs))
        #     out.close()

    if str.upper(noise) in ("YES","Y","NOISE"):
        for j in xrange(runs):
            out = open("Run{}_Wi{}_chi{}".format(j, k,chi_element), "w")
            particles = [np.array(initial_positions[0:3]),np.array(initial_positions[3:6]),
                         np.array(initial_positions[6:9]),np.array(initial_positions[9:12])]
            particles_interm = [None]*4
            polyvec1 = particles[1] - particles[0]
            polyvec2 = particles[3] - particles[2]
            out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (1), polyvec1[0], polyvec1[1], polyvec1[2],
            np.linalg.norm(polyvec1),polyvec2[0], polyvec2[1], polyvec2[2],np.linalg.norm(polyvec2)))
            for i in xrange(steps):


                # Constructing the oseen tensor (stokeslet matrix) and a list of the separation vectors between particles
                # The list of vectors goes like: pointing from particle ()- to particle ()
                # 1-2, 1-3, 2-3, 1-4, 2-4, 3-4
                oseen_tensor, sep_lst = oseen(particles)
                # Mxy = - (3*ar_ratio/4)*x*y/(rseparation + epsilon_squared)**(1.5)
                # Mxz = - (3* ar_ratio/4)*x*z/(rseparation + epsilon_squared)**(1.5)
                # Myy = 1- (3* ar_ratio/4)*(rseparation + 2*epsilon_squared + y*y)/(rseparation + epsilon_squared)**(1.5)
                # Myz = - (3* ar_ratio/4)*y*z/(rseparation + epsilon_squared)**(1.5)
                # Mzz = 1- (3* ar_ratio/4)*((rseparation + 2*epsilon_squared) + z*z)/(rseparation + epsilon_squared)**(1.5)
                noise_vector = noise_producer(oseen_tensor, chi_element)
                # Updating new position of particles
                particles_interm[0] = particles[0] + time_step*function(Part1,Part2,Part3,Part4,sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part1]
                particles_interm[1] = particles[1] + time_step*function(Part2,Part1,Part3,Part4,-sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part2]
                particles_interm[2] = particles[2] + time_step*function(Part3,Part4,Part1,Part2,sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part3]
                particles_interm[3] = particles[3] + time_step*function(Part4,Part3,Part1,Part2,-sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part4]

                oseen_tensor, sep_lst = oseen(particles_interm)

                particles[0] = particles_interm[0] + time_step*function(Part1,Part2,Part3,Part4,sep_lst[0],sep_lst[5],particles,oseen_tensor,k) + noise_vector[Part1]
                particles[1] = particles_interm[1] + time_step*function(Part2,Part1,Part3,Part4,-sep_lst[0],sep_lst[5],particles,oseen_tensor,k) + noise_vector[Part2]
                particles[2] = particles_interm[2] + time_step*function(Part3,Part4,Part1,Part2,sep_lst[5],sep_lst[0],particles,oseen_tensor,k) + noise_vector[Part3]
                particles[3] = particles_interm[3] + time_step*function(Part4,Part3,Part1,Part2,-sep_lst[5],sep_lst[0],particles,oseen_tensor,k) + noise_vector[Part4]

                particles[0] = particles[0]+ np.array([time_step*k*particles[0][1],0,0]) + (time_step *
                    (sep_lst[0]) / (2 * (1 - np.linalg.norm(sep_lst[0]) ** 2))) + time_step * (
                    -oseen_tensor[3:6,0:3].dot(sep_lst[0]) + oseen_tensor[6:9, 0:3].dot(sep_lst[5])
                    -oseen_tensor[9:12,0:3].dot(sep_lst[5])) + noise_vector[0:3]

                particles[1] = particles[1]+ np.array([time_step*k*particles[1][1],0,0]) + (time_step *
                                 -(sep_lst[0]) / (2 * (1 - np.linalg.norm(sep_lst[0]) ** 2))) + time_step * (
                                oseen_tensor[3:6,0:3].dot(sep_lst[0]) + oseen_tensor[6:9, 3:6].dot(sep_lst[5])
                                -oseen_tensor[9:12,3:6].dot(sep_lst[5])) + noise_vector[3:6]

                particles[2] = particles[2]+ np.array([time_step*k*particles[2][1],0,0]) + (time_step *
                                                                                            (sep_lst[5]) / (2 * (1 - np.linalg.norm(sep_lst[5]) ** 2))) + time_step * (
                    -oseen_tensor[6:9,3:6].dot(sep_lst[0]) + oseen_tensor[6:9, 0:3].dot(sep_lst[0])
                    -oseen_tensor[9:12,6:9].dot(sep_lst[5])) + noise_vector[6:9]

                particles[3] = particles[3]+ np.array([time_step*k*particles[3][1],0,0]) + (time_step *
                                                                                            -(sep_lst[5]) / (2 * (1 - np.linalg.norm(sep_lst[5]) ** 2))) + time_step * (
                                oseen_tensor[9:12,0:3].dot(sep_lst[0]) - oseen_tensor[9:12, 3:6].dot(sep_lst[0])
                                +oseen_tensor[9:12,6:9].dot(sep_lst[5])) + noise_vector[9:12]
                polyvec1 = particles[1] - particles[0]
                polyvec2 = particles[3] - particles[2]
                out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (i + 1),polyvec1[0],polyvec1[1],polyvec1[2],np.linalg.norm(polyvec1),
                                                    polyvec2[0], polyvec2[1], polyvec2[2], np.linalg.norm(polyvec2)))
            update_progress(j / (runs))
            out.close()

def function(refparticle,refparticlepair,particle3,particle4,sepvec,sepvec2, particles, oseen_tensor,k):

    """

    :param refparticle: The particle on which the force acts on
    :param refparticlepair: The particle on the same chain of the reference particle
    :param particle3: The 3rd particle
    :param particle4: The 4th particle

    :param sepvec: Separation vector from reference particle to its pair on the same chain
    :param sepvec2: Separation vector between the other two particles from smallest particle number
                    to largest (1 to 2, or, 3 to 4)
    :param particles:
    :return: The velocity of the reference particle at the current time step
    """
    vector = np.array([ k * particles[int(refparticle[0]/3)][1], 0, 0]) + (
            (sepvec) / (2 * (1 - np.linalg.norm(sepvec) ** 2))) +  (
        -oseen_tensor[refparticlepair, refparticle].dot(sepvec) + oseen_tensor[particle3, refparticle].dot(sepvec2)
        - oseen_tensor[particle4, refparticle].dot(sepvec2))
    return vector
def oseen(particles):
    # Construct list of separations between particles
    sep_lst = []
    oseen_tensor = np.zeros((12,12))
    for i in range(4):
        for j in range(i):
            # separation vector
            svec = particles[i] - particles[j]
            # separation magnitude
            smagn = np.linalg.norm(svec)**2
            sep_lst.append(svec)
            # The units might not be correct
            oseen_tensor[i*3:(i+1)*3,j*3:(j+1)*3] = (3*ar_ratio/8)*(np.identity(3)*(smagn + 2*epsilon_squared)+ (np.outer(svec,svec)))/(smagn +epsilon_squared)**(1.5)

    # Filling in the diagonal elements of the oseen tensor and returning the symmetric matrix

    np.fill_diagonal(oseen_tensor,1/2)
    return np.maximum(oseen_tensor,oseen_tensor.T),sep_lst
def noise_producer(matrix,chi_element):
    """
    This method calculates the noise exerted on the 2 spheres and returns a 3 dimensional vector.
    :param matrix: Stokeslet matrix (oseen tensor)
    :return: Random noise vector
    """
    # Creating the sigma matrix
    sigma_m = np.zeros((12,12))
    sigma_m[0,0] = math.sqrt(matrix[0,0])
    for i in range(12):
        for j in range(i):
            sigma_m[i,j] = matrix[i,j]
            k = 0
            while k < j:
                sigma_m[i,j] -= sigma_m[i,k]*sigma_m[j,k]
                k += 1

            sigma_m[i,j] = sigma_m[i,j]/sigma_m[j,j]
        print (sigma_m)
        print ("/n")
        sigma_m[i,i] = math.sqrt(matrix[i,i] - np.sum(sigma_m[i,0:i]**2))


    return sigma_m.dot(np.random.normal(0,math.sqrt(chi_element* time_step),12))
def analyse():
    """
    This function goes to the folder of the previously run simulation, and averages over
    all the runs (which, because of thermal noise, are different). It writes the mean
    squared separation and mean separation of each constant simulated in a file. It also
    creates a file with the maximum separation of every constant and at which time it occurred
    """

    os.chdir("Shear_constants:{}-{}_numc:{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],
    constant[-1],len(constant),hydro, steps,time_step, ar_ratio,noise))
    max_file = open("AMax_File_con"
                    ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
        constant[0], constant[-1], len(constant), hydro, steps, time_step, ar_ratio, noise), "w")
    max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")

    for v in constant:
        #Rfile = open(os.getcwd() + "/Results.out", "a")

        out = open("MSS_{}_{}_{}_{}_{}_{}.out".format(v,hydro, steps,time_step,ar_ratio,noise), "w")
        nout = open("MS_{}_{}_{}_{}_{}_{}.out".format(v,hydro, steps,time_step,ar_ratio,noise), "w")

        #The following way of reading files is because of a limitation
        #in the number of files a computer can have open at the same time
        thousands_of_runs = int(math.ceil(runs / 1000))
        ms_list = []
        mss_list = []
        # Reads every thousand runs of a simulation
        for k in range(thousands_of_runs):
            # Opens the first 1000 runs in a dictionary, then opens the next 1000 and so on.
            filedata = {i: open("Run{}_Wi{}_chi{}".format(i,v), "r") for i in xrange(k * 1000, min(runs, (k + 1) * 1000))}
            # Mean separation and Mean square separation lists that contain temporary files
            # with the respective values for every thousand runs. They are deleted afterwards
            ms_list.append(open("ms_{}th_thousand.tmp".format(k), "w"))
            mss_list.append(open("mss_{}th_thousand.tmp".format(k), "w"))

            # Adding squared separation and separation together
            # to average noise
            for lines in xrange(steps + 1):
                s1 = 0
                ssq1 = 0
                s2 = 0
                ssq2 = 0

                for file in filedata.values():
                    token = str.split(file.readline())
                    # This convenion will most likely change in the 3rd version of the program
                    t = float(token[0])
                    x1 = float(token[1])
                    y1 = float(token[2])
                    z1 = float(token[3])
                    x2 = float(token[5])
                    y2 = float(token[6])
                    z2 = float(token[7])
                    rsepparation1 = float(token[4])
                    rsepparation2 = float(token[8])

                    s1 += rsepparation1**2
                    ssq1 += rsepparation1
                    s2 += rsepparation2 ** 2
                    ssq2 += rsepparation2
                mss_list[k].write("{} {} {}\n".format(t, s1 / runs,s2 / runs))
                ms_list[k].write("{} {} {}\n".format(t, (ssq1 / runs),(ssq2 / runs)))
                update_progress(lines / (steps))
            for fruns in filedata.values():
                fruns.close()
            ms_list[k].close()
            mss_list[k].close()
            ms_list[k] = open("ms_{}th_thousand.tmp".format(k), "r")
            mss_list[k] = open("mss_{}th_thousand.tmp".format(k), "r")

        #THIS HAS NOT BEEN UPDATED TO WORK WITH TWO POLYMERS YET
        # This loop goes through the temporary file in ms_list and mss_list and finds the
        # largest separation. It also finds the mean separation and separation squared if
        # the number of runs was more than 1000. If its under 1000 runs then this loop will
        # slow down the computation by a bit.
        # ~~~~~~~~~ NOTE: If computation time is an issue then modify this ~~~~~~~~~~~~~~~~~~~~~~~~~~
        print "~~~~~~~Merging and finding Max value~~~~~~~~~"
        maxr = 0
        tmax = 0
        for j in xrange(steps + 1):
            mean_mss = 0
            mean_ms = 0

            for k in range(thousands_of_runs):
                mstoken = str.split(ms_list[k].readline())
                msstoken = str.split(mss_list[k].readline())
                t = float(mstoken[0])
                mssn = float(msstoken[1])
                msn = float(mstoken[1])
                mean_mss += mssn
                mean_ms += msn

            out.write("{} {}\n".format(t, mean_mss))
            nout.write("{} {}\n".format(t, mean_ms))
            if maxr <= mean_mss:
                maxr = mean_mss
                tmax = t
        # Max separation squared over initial separation squared is stored in a max file for
        # every constant
        # The loop deletes the unnecessary temporary files
        init_separation = np.linalg.norm(initial_positions[3:6] - initial_positions[0:3])
        max_file.write("{} {} {}\n".format(v, maxr / (init_separation ** 2), tmax))
        for k in range(thousands_of_runs):
            os.remove(mss_list[k].name)
            os.remove(ms_list[k].name)
        out.close()
        nout.close()
        meansqsep = float(str.split(linecache.getline(out.name, steps + 1))[1])
        meansep = float(str.split(linecache.getline(nout.name, steps + 1))[1])
        #print("Mean squared separation over {} runs: {} ".format(runs, meansqsep))
        #print("Root Mean squared separation over {} runs: {} ".format(runs, math.sqrt(meansqsep)))
        #print ("Mean separation over {} runs : {}".format(runs, meansep))
        # Appending the results at the end of the Results.out file
        # Rfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~R E S U L T S~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
        # Rfile.write(
        #     "Time-Step:{} Steps:{} runs:{} constant:{} Initial separation:{} hydro: {} time&date: {} \n".format(time_step,
        #                                                                                                         steps, runs,
        #                                                                                                         constant[j],
        #                                                                                                         init_separation,
        #                                                                                                         hydro,
        #                                                                                                         time.strftime(
        #                                                                                                             "%c")))
        # Rfile.write("Mean Squared separation {}\n".format(meansqsep))
        # Rfile.write("Root Mean squared separation {}\n".format(math.sqrt(meansqsep)))
        # Rfile.write("Mean separation {}\n".format(meansep))
        # Rfile.close()
        # Mean squared displacement. Each row has a colour of the rainbow.
        # if args.walk:
        #     os.chdir("Max_Separation_constants:{}-{}_numc:{}".format(constant[0],constant[-1],len(constant)))
        #     max_file = open("Max_Separation_constants:{}-{}_numc:{}".format(constant[0],constant[-1],len(constant)),"w")
        #     for n,j in enumerate(constant):

def average_sepparation():
    os.chdir("Shear_Wi:{}-{}_chi:{}-{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0], constant[-1],
    chi[0], chi[-1], hydro, steps,time_step, ar_ratio, noise))

    for j,chi_element in enumerate(chi):
        avsep_file = open(("Av.Sep_hydro:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro, chi_element,
                                                                constant[0],constant[-1],len(constant)),"w")
        for i,wi in enumerate(constant):
            wisep_array = np.loadtxt("Run{}_Wi{}_chi{}".format(0, wi,chi_element))
            avsep_file.write("{} {}\n".format(wi,np.mean(wisep_array[500:,-1])))
            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))

def angle():
    os.chdir("Shear_Wi:{}-{}_chi:{}-{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0], constant[-1],
    chi[0], chi[-1], hydro, steps,time_step, ar_ratio, noise))

    for j,chi_element in enumerate(chi):
        avsep_file = open(("Angle_hydro:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro, chi_element,
                                                                constant[0],constant[-1],len(constant)),"w")
        for i,wi in enumerate(constant):
            file_array = np.loadtxt("Run{}_Wi{}_chi{}".format(0, wi,chi_element))
            angle_array = np.arccos(np.absolute(file_array[:,1])/file_array[:,-1])
            avsep_file.write("{} {}\n".format(wi,np.mean(angle_array[500:])))
            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))

def simulate():
    """
    For version 1 the creation of folder is inside the for loop.
    For version 2 all simulations are created in a single folder but separate files.
    :return:
    """
    #This code simulates the walk of the polymer and stores the max separation in a file

    try:
        os.mkdir("Shear_Wi:{}-{}_chi:{}-{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],constant[-1],
        chi[0],chi[-1],hydro,steps,time_step,ar_ratio,noise))
    except OSError as e:
        if e.errno != 17:
            raise
    os.chdir("Shear_Wi:{}-{}_chi:{}-{}_hydro:{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0],constant[-1],
        chi[0],chi[-1],hydro,steps,time_step,ar_ratio,noise))
    if args.max:
        max_file = open("Max_File_con"
                        ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
            constant[0],constant[-1],len(constant),hydro,steps,time_step,ar_ratio,noise),"w")
        max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")
    else:
        max_file = None
    for o,k in enumerate(chi):
        for n,j in enumerate(constant):

            walk(j,max_file,k)
            print("wi number {} done, {} left".format(j,len(constant)-n -1))
        print("chi number {} done, {} left".format(k, len(chi) -o-1))


def update_progress(progress):
    barLength = 10  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__ == "__main__":
    """
    Program gets variables from config file. 
    To run the program, arguments have to be passed while initialising it in the form of flags.
    The -walk argument runs the random walk simulation using the parameters of the config file
    The -analyse argument finds the folder that matches the parameters of the config file
    and subsequently runs through the files to produce a mean separation squared and max separation
    files for every constant.
    The -max argument is used only when there is no noise in the simulation so as to avoid running
    the computationally expensive -analyse method 
    """
    #Constants used in the simulation
    steps = int(config.steps)
    runs = int(config.runs)
    constant = config.constant #Weissenberg numbers in a numpy array
    time_step = float(config.time_step)
    initial_positions = config.initial_positions
    # init_separation = math.sqrt(xinitial**2 + yinitial**2 + zinitial**2)
    hydro = str(config.hydro)
    ar_ratio = float(config.ar_ratio)
    noise = str(config.noise)
    chi = config.chi
    epsilon_squared = 4*ar_ratio*ar_ratio
    # Lists that signify particle position in the Oseen tensor to make code more readable
    Part1 = range(0,3)
    Part2 = range(3,6)
    Part3 = range(6,9)
    Part4 = range(9,12)
    parser = argparse.ArgumentParser(description="Program version 2"
                                                 "The program simulates the motion of a polymer in shear flow.\n "
                                                 "The model is of a finite extensibility non-linear elastic spring"
                                                 "(FENE).\n"
                                                 "Parameters of the simulation can be found in the config.py file.\n"
                                                 "\nVersion 1: Has no thermal fluctuations so analyser doesn't do anything\n"
                                     "Version 2: There are hydrodynamic interactions between the two ends of the polymer chain.\n"
                                     "Version 3: There is thermal noise acting on the polymer")
    # parser.add_argument("echo", help="echo the string you use here")
    parser.add_argument("-a", "--analyse", help="Run analyser. For version 1 this does nothing (v1) since there "
                                                "are no thermal fluctuations.",
                        action="store_true")
    parser.add_argument("-w", "--walk", help="Simulate walks with parameters of config.py file", action="store_true")

    parser.add_argument("-m", "--max", help="Store max separation squared in a separate file", action="store_true")

    parser.add_argument("-av", "--average", help="Averages separation for each weisenberg number", action="store_true")

    parser.add_argument("-an", "--angle", help="Averages angle between separation vector and x-axis for each weisenberg number", action="store_true")

    args = parser.parse_args()
    # print args.echo
    if args.walk:
        simulate()
    if args.analyse:
        analyse()
    if args.average:
        average_sepparation()
    if args.angle:
        angle()