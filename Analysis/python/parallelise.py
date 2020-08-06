import sys, os, time, glob
import subprocess
import multiprocessing as mp

from Utils import *
from fileManipulation import *


def MultiProcDecorator(f,*args, **kw):
    def g(*args, **kw):
        result = f(*args, **kw)
        return result
    g.__module__ = "__main__"
    return g

@timeit
def MultiProcess(func,start,stop,step,*args):
    # func = MultiProcDecorator(func)
    nproc = 20
    print "Running on ", ((stop-start)/step), "processes in step of ", nproc
    pool = mp.Pool(processes=nproc)
    result = [pool.apply_async(func, args=(x,x+step)+args) for x in range(start,stop,step)]
    result = [p.get() for p in result if p.get() is not None]
    pool.close()
    pool.join()
    return result


@timeit
def parallelise(list_processes, MaxProcess=10, list_logfiles=[], cwd=None):
  ntotal = len(list_processes)
  processes = []
  logfiles = []
  condition = len(list_logfiles)>0 and len(list_logfiles)==len(list_processes)
  for index, process in enumerate(list_processes):
    wait = True
    while wait:
      nrunning = 0
      ncompleted = 0
      idx=0
      for proc in processes:
        if proc.poll() == None :
          nrunning += 1
        else:
          ncompleted += 1
          if not logfiles[idx].closed:
            logfiles[idx].close()
            print 'Job "%s" has finished.' % logfiles[idx].name
        idx += 1
      if nrunning >= MaxProcess:
        percentage = float(ncompleted)/float(ntotal)*100
        sys.stdout.write( 'Already completed '+str(ncompleted)+' out of '+str(ntotal)+' jobs --> '+str(percentage)+'%. Currently running: '+str(nrunning)+' \r')
        sys.stdout.flush()
        time.sleep(5)
      else:
        print 'only %i jobs are running, going to spawn new ones.' % nrunning
        wait = False
    if condition:
      f = open(list_logfiles[index],'w')
    else:
      f = open("log_"+str(index)+".txt",'w')
    logfiles.append(f)
    if cwd:
      processes.append(subprocess.Popen(process[1:], stdout=f, cwd=process[0]))
      time.sleep(1)
    else:
      processes.append(subprocess.Popen(process, stdout=f))

  for proc in processes:
    proc.wait()
  for file in logfiles:
    file.close()
  if not condition:
    # os.remove("log.txt")
    a = map(os.remove, glob.glob("log_*.txt"))
