#!/usr/bin/env python
__author__ = "I-Ta Hsieh"
__copyright__ = "Copyright 2017, I-Ta Hsieh(ASIAA)"
__credits__ = ["I-Ta Hsieh"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "I-Ta Hsieh"
__email__ = "ita.hsieh@gamil.com"
__status__ = "Developing"

#######################################################
# set up the parameters of the genetic algorithm
GoalOfError = 1e-2
MAX_GEN = 100
TemperatureFactor = 1.0  
cooling_factor = 1.0
# First ancestor
# the parameters of the searching domain
params = {'density_cutoff':0.0002,
          'density_scale': 1e15,
          'temp1':71.0,
          'temp2':60.7}
params_fixed = {'density_powerlaw':-2.2}

#######################################################
from numpy import zeros
from math import log,exp,tan,pi
import random
import time
from subprocess import call

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

nparams = len(params)
npop = size
max_nselect = npop/2
npop_proc = npop /size

#######################################################
import target as tg
import run_chromosome

# run a specific population
def run_pop(par,chromosome_id):
    global params_fixed
    result = run_chromosome.run( par, params_fixed, chromosome_id, tg)
    # return the value of weighted difference 
    return compare_error(result)


ndata=0
for i in range(len(tg.luminance)):
    ndata += len(tg.luminance[i])

# compare the result to the observational data
def compare_error(simulation):
    error=0.
    global ndata
    for i in range(len(tg.luminance)):
        for j in range(len(tg.luminance[i])):
            ratio = simulation[i][j] / tg.luminance[i][j][2]
            if ratio > 0.0:
                rel_diff = abs( log( ratio ) )
            else: 
                rel_diff = 10.0
            error += rel_diff
    return exp( error / ndata ) - 1.0


# mutate from the reproduced chromosome
def mutate(par,FurnaceTemperature):
    global TemperatureFactor
    new_params = par.copy()
    for key in par.keys():
        RandomDistribution = 1e3
        while abs(RandomDistribution) > 1e2:
            RandomFloat = random.uniform(-1.0, 1.0) 
            RandomDistribution = FurnaceTemperature * tan( 0.5 * pi * RandomFloat )
        factor = exp( RandomDistribution )
        new_params[key] = TemperatureFactor * factor * par[key]
    return new_params


def GatherErrorAndPopulation_Selection(critical_error):
    global error_proc, pops_proc
    global pops_current
    global npop_proc
    global rank
    
    if rank == 0:
        error = zeros(npop)
        error[0:npop_proc] = error_proc[0:npop_proc]
        pops_current = [params.copy() for i in range(npop)]
        for i in range(npop_proc):
            pops_current[i] = pops_proc[i].copy()

    comm.Barrier()

    for proc_id in range(1,size):
        if rank == proc_id:

            comm.send( error_proc, dest=0, tag=100)
            comm.send( pops_proc, dest=0, tag=101)

        elif rank == 0:
            error_proc = comm.recv( source = proc_id, tag=100)
            error[ proc_id*npop_proc: (proc_id+1)*npop_proc] = error_proc
            pops_proc = comm.recv( source = proc_id, tag=101)
            for i in range(npop_proc):
                pops_current[proc_id*npop_proc+i] = pops_proc[i].copy()

    global sorted_error,sorted_index
    if rank == 0:
        # sort the error
        sorted_index = sorted(range(npop), key=lambda k: error[k])
        sorted_error = error[ [sorted_index[i] for i in range(npop)] ]
        # find the children who are better than the father
        i = 0
        while sorted_error[i] < critical_error:    
            i += 1
            if i == npop: break
        nselect = i
    else: 
        nselect = None

    return nselect

  
def message_output(gen_id, NumberOfSurvivor, error):
    hour,minute,second = GetElapsedTime()
    print 'GEN %6d. %3d suvivors. Minimum error:%10.5e. Elapsed:%5d:%02d:%02d' %( gen_id, NumberOfSurvivor, error, hour, minute, second)


def GetElapsedTime():
    global t0
    elapsed_time = time.time() - t0
    hour = int(elapsed_time/3600)
    minute = int(elapsed_time/60%60)
    second = int(elapsed_time%60)
    return hour,minute,second

#######################################################





GEN_ID = 0

if rank == 0:
    # create the working directory if it doesn't exist
    call(['rm','-rf','working_dir'])
    call(['mkdir','working_dir'])

    # use timer measure execution time
    t0 = time.time()
 
    # test the first ancestor
    error_ancestor = run_pop(params,rank)

    hour,minute,second = GetElapsedTime()
    print 'GEN %6d.   The error of the ancestor:%10.5e. Elapsed:%5d:%02d:%02d' %(GEN_ID,error_ancestor,hour,minute,second)
else:
    error_ancestor = None 

comm.Barrier()
error_ancestor = comm.bcast(error_ancestor,root=0)

while 1:
    # generate the first generation
    pops_proc = []
    error_proc = zeros(npop_proc)
    for i in range(npop_proc):
        # self-mutation and test the first generation
        pops_proc.append( mutate( params, error_ancestor) )
        error_proc[i] = run_pop( pops_proc[i], i * size + rank )
    # Gather all error and population
    nselect = GatherErrorAndPopulation_Selection(error_ancestor)
    
    comm.Barrier()
    nselect = comm.bcast(nselect,root=0)

    if nselect == 0: 
        TemperatureFactor *= cooling_factor
        if rank == 0:
            hour,minute,second = GetElapsedTime()
            print 'GEN %6d.   No survivor apears, multiply again...  Elapsed:%5d:%02d:%02d' %(GEN_ID,hour,minute,second)
    else:
        break

GEN_ID = 1
if rank == 0:
    # the suvivors are the selected children plus the ancestor 
    if nselect < max_nselect:
        nselect += 1
        
        pops_select = [params.copy() for i in range(nselect)]
        for i in range(nselect-1):
            pops_select[i] = pops_current[sorted_index[i]].copy()
        pops_select[nselect-1] = params.copy()
        
        error_select = zeros(nselect)
        error_select[0:nselect-1] = sorted_error[0:nselect-1]
        error_select[nselect-1] = error_ancestor
    else:
        nselect = max_nselect
        
        pops_select = [params.copy() for i in range(nselect)]
        for i in range(nselect):
            pops_select[i] = pops_current[sorted_index[i]].copy()

        error_select = zeros(nselect)
        error_select[0:nselect] = sorted_error[0:nselect]

    minimum_error = error_select[0]
    message_output(GEN_ID,nselect,minimum_error)



###############################################################
# the first generation of survivors appears now
# Let's go genetic algorithm
for GEN_ID in range(2,MAX_GEN+1):
    # store the current population to the previous
    if rank == 0:
        nselect_old = nselect
        error_old = error_select
        pops_old = pops_select
        minimum_error_old = minimum_error
    else:
        nselect_old = None
        error_old = None
        pops_old = None
        minimum_error_old = None

    nselect_old = comm.bcast(nselect_old,root=0)
    error_old = comm.bcast(error_old,root=0)
    pops_old = comm.bcast(pops_old,root=0)
    minimum_error_old = comm.bcast(minimum_error_old,root=0)

    params_crossover = params.copy()
    for i in range(npop_proc):
        # crossover
        random_combination = random.sample(xrange(nselect_old), 2)
        random_keys_combination = random.sample(params_crossover.keys(), nparams/2)
        for key in params_crossover.keys():
            if key in random_keys_combination:
                params_crossover[key] = pops_old[ random_combination[0] ][ key ]
            else:
                params_crossover[key] = pops_old[ random_combination[1] ][ key ]
        error_crossover = min( error_old[random_combination[0]],
                               error_old[random_combination[1]])

        # mutation
        pops_proc[i] = mutate( params_crossover, error_crossover)
        
        # test
        child_id = i * size + rank
        error_proc[i] = run_pop( pops_proc[i], child_id )

    # Gather all error and population, then select
    nselect = GatherErrorAndPopulation_Selection(error_old[nselect_old-1])

    # the suvivors are the selected children plus the previous generation 
    if rank == 0:
        # combine the new and old selected list
        nselect_new = nselect
        nselect += nselect_old

        pops_select = []
        for i in range( nselect_new):
            pops_select.append( pops_current[sorted_index[i]].copy() )
        for i in range( nselect_old):
            pops_select.append( pops_old[i].copy() )

        error_select = zeros(nselect)
        error_select[0:nselect_new] = sorted_error[0:nselect_new]
        error_select[nselect_new:nselect] = error_old[0:nselect_old]
        
        # sort angain
        sorted_index = sorted(range(nselect), key=lambda k: error_select[k])
        sorted_error = error_select[ [sorted_index[i] for i in range(nselect)] ]
        sorted_pops = []
        for i in range(nselect):
            sorted_pops.append( pops_select[sorted_index[i]].copy() )
        
        # if list is longer than the default length of selection, reset
        if nselect > max_nselect:
            nselect = max_nselect
            error_select = sorted_error[0:max_nselect]
            pops_select = sorted_pops[0:max_nselect]
        else:
            error_select = sorted_error
            pops_select = sorted_pops
        minimum_error = error_select[0]
        
        message_output(GEN_ID,nselect,minimum_error)
    else:
        minimum_error = None

    minimum_error = comm.bcast(minimum_error,root=0)
    
    if minimum_error == minimum_error_old:
        TemperatureFactor *= cooling_factor
    
    if minimum_error < GoalOfError : break

# Final output message
if rank == 0:
    print 'The fixed parameters:'
    for key in params_fixed:
        print key,params_fixed[key]
        
    print 'The chosen parameters:'
    for key in pops_select[0]:
        print key,pops_select[0][key]
    
    print 'The error:',minimum_error



