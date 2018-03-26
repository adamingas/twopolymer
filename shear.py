from __future__ import division

import sys
import os
import math
import linecache
import config
import argparse
import numpy as np
import timeit
import matplotlib.pyplot as plt

"""
Version 4.0.1
Bug fix in randomiser.
Version 4.0
The position of the second polymer is randomised around a sphere of radius l of the first one. 
This happens for every run. When the polymer passes outside this sphere, its position is randomised again.
Version 3.0
Changed to two polymer chains. Simulate run for a number of constants and find maximum separation squared.
"""

def walk(k,max_file,chi_element,l):
    """
    Performs the random walk.
    There are 4 possibilities of a walk, which are combinations of the binary states of hydrodynamic
    interactions and thermal noise.
    k is the weissenberg number
    """


    if str.upper(noise) in ("YES","Y","NOISE"):
        for j in xrange(runs):
            out = open("Run{}_Wi{}_chi{}_l{}".format(j, k,chi_element,l), "w")

            p3,p4 = config.randomiser(l)
            particles = [config.p1,config.p2,p3,p4]
            particles_interm = [None]*4
            #polyvec1 = particles[1] - particles[0]
            #polyvec2 = particles[3] - particles[2]
            # out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (1), polyvec1[0], polyvec1[1], polyvec1[2],
            # np.linalg.norm(polyvec1),polyvec2[0], polyvec2[1], polyvec2[2],np.linalg.norm(polyvec2)))
            out.write("{} {} {} {} {}\n".format(0, str(particles[0]).strip("array([,])"),
                        str(particles[1]).strip("array([,])"), str(particles[2]).strip("array([,])"),
                        str(particles[3]).strip("array([,])")))

            for i in xrange(steps):


                # Constructing the oseen tensor (stokeslet matrix) and a list of the separation vectors between particles
                # The list of vectors goes like: pointing from particle ()- to particle ()
                # 1-2, 1-3, 2-3, 1-4, 2-4, 3-4
                oseen_tensor, sep_lst = oseen(particles)
                # Chi_element is divided by two because the time step is half in the intermediate step
                noise_vector = noise_producer(oseen_tensor, chi_element/2)
                # Updating new position of particles
                particles_interm[0] = particles[0] + time_step*function(0,1,2,3,sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part1]
                particles_interm[1] = particles[1] + time_step*function(1,0,2,3,-sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part2]
                particles_interm[2] = particles[2] + time_step*function(2,3,0,1,sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part3]
                particles_interm[3] = particles[3] + time_step*function(3,2,0,1,-sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part4]

                oseen_tensor, sep_lst = oseen(particles_interm)
                noise_vector = noise_producer(oseen_tensor, chi_element)

                particles[0] = particles[0] + time_step*function(0,1,2,3,sep_lst[0],sep_lst[5],particles_interm,oseen_tensor,k) + noise_vector[Part1]
                particles[1] = particles[1] + time_step*function(1,0,2,3,-sep_lst[0],sep_lst[5],particles_interm,oseen_tensor,k) + noise_vector[Part2]
                particles[2] = particles[2] + time_step*function(2,3,0,1,sep_lst[5],sep_lst[0],particles_interm,oseen_tensor,k) + noise_vector[Part3]
                particles[3] = particles[3] + time_step*function(3,2,0,1,-sep_lst[5],sep_lst[0],particles_interm,oseen_tensor,k) + noise_vector[Part4]

                com1 = (particles[0] + particles[1]) /2
                com2 = (particles[3] + particles[2]) / 2
                comsep = np.linalg.norm(com2 - com1)
                if comsep >= l:

                    particles[2],particles[3] = randomiser(l,com1,particles[2],particles[3])

                # polyvec1 = particles[1] - particles[0]
                # polyvec2 = particles[3] - particles[2]
                # out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (i + 1),polyvec1[0],polyvec1[1],polyvec1[2],np.linalg.norm(polyvec1),
                #                                     polyvec2[0], polyvec2[1], polyvec2[2], np.linalg.norm(polyvec2)))
                out.write("{} {} {} {} {}\n".format(time_step * (i), str(particles[0]).strip("array([,])"),
                str(particles[1]).strip("array([,])"),str(particles[2]).strip("array([,])"),str(particles[3]).strip("array([,])")))
            update_progress(j / (runs))
            out.close()

def function(ref, refpair, p3, p4, sepvec, sepvec2, particles, oseen_tensor, k):

    """

    :param ref: The particle on which the force acts on
    :param refpair: The particle on the same chain of the reference particle
    :param p3: The 3rd particle
    :param p4: The 4th particle

    :param sepvec: Separation vector from reference particle to its pair on the same chain
    :param sepvec2: Separation vector between the other two particles from smallest particle number
                    to largest (1 to 2, or, 3 to 4)
    :param particles:
    :return: The velocity of the reference particle at the current time step
    """
    magvec = np.linalg.norm(sepvec)**2
    magvec2 = np.linalg.norm(sepvec2)**2

    hydro_forces = ((sepvec) / (2 * (1 - np.linalg.norm(sepvec) ** 2))) + (
        -oseen_tensor[refpair*3:(refpair+1)*3, ref*3:(ref+1)*3].dot(sepvec)/(1-magvec)) + (oseen_tensor[3*p3:3*(p3+1), ref*3:(ref+1)*3].dot(sepvec2)
        - oseen_tensor[3*p4:3*(p4+1), ref*3:(ref+1)*3].dot(sepvec2))/(1-magvec2)

    vector = np.array([k * particles[int(ref)][1], 0, 0]) + hydro_forces



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
        sigma_m[i,i] = math.sqrt(matrix[i,i] - np.sum(sigma_m[i,0:i]**2))


    return sigma_m.dot(np.random.normal(0,math.sqrt(chi_element* time_step),12))
def randomiser(l,com1,p3,p4):
    radius = np.random.uniform() * l
    costh = np.random.uniform(-1, 1)
    fi = np.random.uniform() * 2 * math.pi
    sinth = np.sqrt(1 - costh ** 2)

    extvec2 = p4 -p3
    particle3 = com1 + np.array([radius * math.cos(fi) * sinth, radius * math.sin(fi) * sinth, radius * costh]) -extvec2/2
    particle4 = particle3+ extvec2
    return particle3,particle4
def analyse():
    """
    This function goes to the folder of the previously run simulation, and averages over
    all the runs (which, because of thermal noise, are different). It writes the mean
    squared separation and mean separation of each constant simulated in a file. It also
    creates a file with the maximum separation of every constant and at which time it occurred
    """


    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    # max_file = open("AMax_File_con"
    #                 ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
    #     constant[0], constant[-1], len(constant), hydro, steps, time_step, ar_ratio, noise), "w")
    # max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")

    for l in ldensity:
        for x in chi:
            for v in constant:
                #Rfile = open(os.getcwd() + "/Results.out", "a")
                # results is the file where all the results will be printed out
                results = open("RES_{}_{}_{}_{}_{}_{}_{}.out".format(l,v, x, steps, time_step, ar_ratio, noise), "w")

                results.write("# Time-step, Mean Separation1, Mean Squared Separation1, Mean Separation2, Mean Squared Separation2, Mean COM separation, Angle 1, Angle 2 \n")

                #The following way of reading files is because of a limitation
                #in the number of files a computer can have open at the same time
                thousands_of_runs = int(math.ceil(runs / 1000))
                averages_list = []

                # Reads every thousand runs of a simulation
                for k in range(thousands_of_runs):
                    # Opens the first 1000 runs in a dictionary, then opens the next 1000 and so on.
                    filedata = {i: open("Run{}_Wi{}_chi{}_l{}".format(i,v,x,l), "r") for i in xrange(k * 1000, min(runs, (k + 1) * 1000))}
                    # Mean separation and Mean square separation lists that contain temporary files
                    # with the respective values for every thousand runs. They are deleted afterwards
                    # Averages_list contains temporary files which hold the average values of various parameters for every
                    # thousand runs of the simulation
                    averages_list.append(open("averages_{}th_thousand.tmp".format(k), "w"))


                    # Adding squared separation and separation together
                    # to average noise
                    for lines in xrange(steps + 1):
                        # sepsq1 and 2 are the separations squared of the polymers. sep1 and 2 are the separations.
                        # comsep is the separation between the centres of mass of the two polymers.
                        sepsq1 = 0
                        sep1 = 0
                        sepsq2 = 0
                        sep2 = 0
                        comsep = 0
                        totangle1 = 0
                        totangle2 = 0

                        for file in filedata.values():
                            token = str.split(file.readline())
                            t = float(token[0])
                            particles = []
                            for i in range(4):
                                particles.append(np.array([float(token[3*i +j +1]) for j in range(3)]))
                            # Polymer separation vector and centre of mass vecotr
                            polyvec1 = particles[1] - particles[0]
                            polymagn1 = np.linalg.norm(polyvec1)
                            polyvec2 = particles[3] - particles[2]
                            polymagn2 = np.linalg.norm(polyvec2)
                            comvec = (particles[3] + particles[2] - particles[1] - particles[0])/2
                            commag = np.linalg.norm(comvec)
                            # Angle with shear axis
                            angle1 = np.degrees(np.arccos(np.clip(abs(np.dot(polyvec1,np.array([1,0,0])))/polymagn1,0,1.0)))
                            angle2 = np.degrees(np.arccos(np.clip(abs(np.dot(polyvec2, np.array([1, 0, 0]))) / polymagn2,0, 1.0)))


                            sepsq1 += polymagn1 ** 2
                            sep1 += polymagn1
                            sepsq2 += polymagn2 ** 2
                            sep2 += polymagn2
                            comsep += commag
                            totangle1 += angle1
                            totangle2 += angle2
                        averages_list[k].write("{} {} {} {} {} {} {} {}\n".format(t,sep1/runs, sepsq1 / runs,sep2/runs, sepsq2 / runs,comsep/runs,totangle1/runs,totangle2/runs))

                        update_progress(lines / (steps))
                    for fruns in filedata.values():
                        fruns.close()
                    averages_list[k].close()

                    averages_list[k] = open("averages_{}th_thousand.tmp".format(k), "r")


                # This loop goes through the temporary file in averages_list and finds the
                # mean separation and separation squared of the two polymers, and the average separation between them
                # if the number of runs was more than 1000. If its under 1000 runs then this loop will
                # slow down the computation by a bit.
                # ~~~~~~~~~ NOTE: If computation time is an issue then modify this ~~~~~~~~~~~~~~~~~~~~~~~~~~
                print "~~~~~~~ Merging ~~~~~~~~~"
                #maxr = 0
                #tmax = 0
                for j in xrange(steps + 1):

                    mean_sep1 =0
                    mean_sep2 = 0
                    mean_sepsq1 = 0
                    mean_sepsq2 = 0
                    mean_comsep = 0
                    mean_angle1 = 0
                    mean_angle2 = 0

                    for k in range(thousands_of_runs):
                        # Reading the averaged parameters and averages the again for every thousand runs
                        token = str.split(averages_list[k].readline())
                        # mstoken = str.split(ms_list[k].readline())
                        # msstoken = str.split(mss_list[k].readline())
                        t = float(token[0])
                        mean_sep1 += float(token[1])
                        mean_sep2 += float(token[3])
                        mean_sepsq1 += float(token[2])
                        mean_sepsq2 += float(token[4])
                        mean_comsep += float(token[5])
                        mean_angle1 += float(token[6])
                        mean_angle2 += float(token[7])

                    results.write("{} {} {} {} {} {} {} {}\n".format(t,mean_sep1,mean_sepsq1,mean_sep2,mean_sepsq2,mean_comsep, mean_angle1,mean_angle2))

                    # if maxr <= mean_mss:
                    #     maxr = mean_mss
                    #     tmax = t
                # Max separation squared over initial separation squared is stored in a max file for
                # every constant
                # The loop deletes the unnecessary temporary files
                # init_separation = np.linalg.norm(initial_positions[3:6] - initial_positions[0:3])
                # max_file.write("{} {} {}\n".format(v, maxr / (init_separation ** 2), tmax))
                for k in range(thousands_of_runs):
                    os.remove(averages_list[k].name)

                results.close()
    os.chdir("..")
    print os.getcwd()
def maximum():
    """
    This method calculates the maximum separation between the polymers and the maximum separation
    between particles of the same polymer
    :return:
    """
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    for j,chi_element in enumerate(chi):
        max_file = open("MAX_chi:{}_Wi:{}-{}_numWi:{}".format( chi_element,constant[0], constant[-1], len(constant)),"w")
        max_file.write("Weisenberg, Max mean squared separation 1, Max mean squared separation 2, Max COM separation \n")
        for i,wi in enumerate(constant):
            results_array = np.loadtxt("RES_{}_{}_{}_{}_{}_{}.out".format(wi, chi_element, steps, time_step, ar_ratio, noise),skiprows=1)
            initial_sepsq1 = results_array[0,2]
            initial_sepsq2 = results_array[0, 4]
            max_sepsq1 = np.amax(results_array[:,2])
            max_sepsq2 = np.amax(results_array[:, 4])
            max_comsep = np.amax(results_array[:, 5])

            max_file.write("{} {} {} {}\n".format(wi,max_sepsq1/initial_sepsq1,max_sepsq2/initial_sepsq2,max_comsep))
            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
    max_file.close()
    os.chdir("..")
def average_sepparation():
    print ("Average Sepparation calculation")
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))


    for m,l in enumerate(ldensity):
        for j,chi_element in enumerate(chi):
            avsep_file = open(("Av.Sep_hydro:{}_l:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro,l, chi_element,
                                                                    constant[0],constant[-1],len(constant)),"w")
            for i,wi in enumerate(constant):
                wisep_array = np.loadtxt("RES_{}_{}_{}_{}_{}_{}_{}.out".format(l,wi, chi_element, steps, time_step, ar_ratio, noise))


                average1 = np.mean(wisep_array[500:,1])
                #average2 = np.mean(wisep_array[500:, 3])
                mean = average1
                avsep_file.write("{} {}\n".format(wi,mean))
                print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
            print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
        print("l number {} done, {} left".format(l, len(ldensity) - m - 1))

    os.chdir("..")

def angle():
    print ("Angle Calculation")
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    for m,l in enumerate(ldensity):
        for j,chi_element in enumerate(chi):
            avsep_file = open(("Angle_hydro:{}_l:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro,l, chi_element,
                                                                    constant[0],constant[-1],len(constant)),"w")
            for i,wi in enumerate(constant):
                file_array = np.loadtxt("RES_{}_{}_{}_{}_{}_{}_{}.out".format(l,wi, chi_element, steps, time_step, ar_ratio, noise))

                angle_array1 = np.mean(file_array[500:,6])
                #angle_array2 = np.mean(file_array[500:,7])
                avsep_file.write("{} {}\n".format(wi,angle_array1))
                print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
            print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
        print("l number {} done, {} left".format(l, len(ldensity) - m - 1))

    os.chdir("..")
def plot():
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    for j,chi_element in enumerate(chi):

        plots = []
        color = iter(plt.cm.rainbow(np.linspace(0, 1,len(ldensity))))

        for m,l in enumerate(ldensity):


            avsep = np.loadtxt("Av.Sep_hydro:{}_l:{}_chi:{}_Wi:{}-{}_numWi:{}".format(hydro,l, chi_element,
                                                                    constant[0],constant[-1],len(constant)))

            c = next(color)
            plt.plot(avsep[:,0],avsep[:,1], c=c, label = "L: {}".format(l))
        plt.legend()
        plt.ylabel("Average Extension")
        plt.xlabel("Weisenberg")
        plt.title("Average Extension vs Weisenberg number for l {}-{} hydro {}".format(ldensity[0],ldensity[-1],hydro))
        plt.savefig("Av.Ext:{}_steps:{}_runs:{}_l:{}-{}.png".format(time_step,steps,runs,ldensity[0],ldensity[-1]))
def simulate():
    """
    For version 1 the creation of folder is inside the for loop.
    For version 2 all simulations are created in a single folder but separate files.
    :return:
    """
    #This code simulates the walk of the polymer and stores the max separation in a file

    try:
        os.mkdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))
    except OSError as e:
        if e.errno != 17:
            raise
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    if args.max:
        max_file = open("Max_File_con"
                        ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
            constant[0],constant[-1],len(constant),hydro,steps,time_step,ar_ratio,noise),"w")
        max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")
    else:
        max_file = None
    for m,l in enumerate(ldensity):
        for o,k in enumerate(chi):
            for n,j in enumerate(constant):

                walk(j,max_file,k,l)
                print("wi number {} done, {} left".format(j,len(constant)-n -1))
            print("chi number {} done, {} left".format(k, len(chi) -o-1))
        print("l number {} done, {} left".format(l, len(ldensity) -m-1))
    os.chdir("..")

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
def distribution():
    os.chdir("Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0],ldensity[-1],constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    for l in ldensity:
        for x in chi:
            for w in d_constant:

                comb_file = open("Comb.ext_Wi{}_chi{}_l{}_ts{}_step{}".format(w,x,l,time_step,steps),"w")
                comb_angle = open("Comb.ang_Wi{}_chi{}_l{}_ts{}_step{}".format(w,x,l,time_step,steps),"w")

                for i in range(runs):
                    file = np.loadtxt("Run{}_Wi{}_chi{}_l{}".format(i, w, x, l))
                    polymag1 = np.sqrt(np.square(file[:,4] -file[:,1]) +np.square(file[:,5] -file[:,2]) +np.square(file[:,6] -file[:,3]))
                    polyvec1x = file[:,4]-file[:,1]
                    angle1 = np.degrees(
                        np.arccos(np.clip(polyvec1x / polymag1, -1.0, 1.0)))

                    comb_file.write("\n".join(map(str,polymag1)))
                    comb_file.write("\n")
                    comb_angle.write("\n".join(map(str, angle1)))
                    comb_angle.write("\n")
                comb_file.close()
                comb_angle.close()

    os.chdir("..")

def histogram():
    os.chdir(
        "Shear_l:{}-{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(ldensity[0], ldensity[-1], constant[0],
                                                                              constant[-1],
                                                                              chi[0], chi[-1], steps, time_step,
                                                                              ar_ratio, noise))
    for l in ldensity:
        for x in chi:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(1, 1, 1)
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(d_constant))))

            for w in d_constant:
                comb_file = np.loadtxt("Comb.ext_Wi{}_chi{}_l{}_ts{}_step{}".format(w, x, l, time_step, steps))
                values, bins = np.histogram(comb_file, density=True, bins=config.bins)
                centre = (bins[:-1] + bins[1:]) / 2


                c = next(color)
                ax.bar(centre, values, color=c,width=(centre[1] - centre[0]),alpha = 0.5 ,label = "Wi {}".format(w))
                ax.legend()
                comb_angle = np.loadtxt("Comb.ang_Wi{}_chi{}_l{}_ts{}_step{}".format(w,x,l,time_step,steps))
                angles , abins = np.histogram(comb_angle,density=True,bins = 180)
                acentre = (abins[:-1] + abins[1:]) / 2


                ax2.bar(acentre, angles,color =c, width=(acentre[1] - acentre[0]),alpha = 0.5, label = "Wi {}".format(w))
                ax2.legend()
            ax.set_title("Normalised Distribution of Extension for l{} w{}".format(l, w))
            ax.set_xlabel("Extension split in {} bins of width {}".format(len(bins), centre[1] - centre[0]))
            ax.set_ylabel("Normalised Frequency")
            fig.savefig("Dist.ext_l{}_Wi{}_chi{}_ts{}_bins{}_steps{}.png".format(l, d_constant , x, time_step, len(bins), steps))

            plt.close(fig)

            ax2.set_title("Normalised Distribution of Angle with x-axis for l{} w{}".format(l, w))
            ax2.set_xlabel("Angle split in {} bins of width {}".format(len(abins), acentre[1] - acentre[0]))
            ax2.set_ylabel("Normalised Frequency")
            fig2.savefig(
                "Dist.ang_l{}_Wi{}_chi{}_ts{}_bins{}_steps{}.png".format(l, d_constant, x, time_step, len(bins), steps))

            plt.close(fig2)
    os.chdir("..")
def fileshistogram(data):
    extdata = data[0]
    angdata = data[1]
    ldata = data[2]


    fig = plt.figure()
    lineax = fig.add_subplot(1, 1, 1)
    #test = plt.figure()
    #lineax = test.add_subplot(1,1,1)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1, 1, 1)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(ldata))))
    for i,l in enumerate(ldata):

        comb_file = np.loadtxt(extdata[i])
        values, bins = np.histogram(comb_file, density=True, bins=config.bins)
        centre = (bins[:-1] + bins[1:]) / 2


        c = next(color)
        #ax.bar(centre, values, color=c,width=(centre[1] - centre[0]),alpha = 0.5 ,label = l)
        lineax.plot(centre, values, color=c,label = l, marker = ".", ms = 8,markevery = (config.ext_markevery))
        lineax.legend()
        #ax.legend()
        comb_angle = np.loadtxt(angdata[i])
        angles , abins = np.histogram(comb_angle,density=True,bins = 180)
        acentre = (abins[:-1] + abins[1:]) / 2


        ax2.plot(acentre, angles,color =c, label = l, marker = ".", ms = 8, markevery = config.ang_markevery)
        ax2.legend()
    lineax.set_title("Normalised Distribution of Extension for l{}".format(ldata))
    lineax.set_xlabel("Extension split in {} bins of width {}".format(len(bins), centre[1] - centre[0]))
    lineax.set_ylabel("Normalised Frequency")
    fig.savefig("Dist.ext_l{}_bins{}.png".format(ldata, len(bins)))

    plt.close(fig)

    ax2.set_title("Normalised Distribution of Angle with x-axis for l{}".format(ldata))
    ax2.set_xlabel("Angle split in {} bins of width {}".format(len(abins), acentre[1] - acentre[0]))
    ax2.set_ylabel("Normalised Frequency")
    fig2.savefig(
        "Dist.ang_l{}.png".format(ldata))

    plt.close(fig2)



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
    # init_separation = math.sqrt(xinitial**2 + yinitial**2 + zinitial**2)
    hydro = str(config.hydro)
    ar_ratio = float(config.ar_ratio)
    noise = str(config.noise)
    chi = config.chi
    ldensity = config.ldensity
    epsilon_squared = 4*ar_ratio*ar_ratio
    d_constant = config.distr_constant
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

    parser.add_argument("-p", "--plot", help="Plot configurations", action="store_true")

    parser.add_argument("-d", "--distribution", help="Finds the distribution of angles and extension", action="store_true")

    parser.add_argument("-hi", "--histogram", help="Plots histogram of angle and extension", action="store_true")

    parser.add_argument("-f",'--files', nargs='+',help = "Takes file names as input. This file names are used to draw histograms")

    parser.add_argument("-app",'--append',nargs = "+", action='append')
    args = parser.parse_args()

    print args.append


    if args.walk:
        simulate()
    if args.analyse:
        print timeit.timeit(analyse, number=1)
    if args.average:
        average_sepparation()
    if args.angle:
        angle()
    if args.max:
        maximum()
    if args.plot:
        plot()
    if args.distribution:
        distribution()
    if args.histogram:
        histogram()
    if args.append:
        fileshistogram(args.append)