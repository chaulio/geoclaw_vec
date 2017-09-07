#!/usr/bin/python

#---  CONFIG  ----
list_versions="1 2 3"
list_threading="singlethread multithread"
list_machines="host mic"
#-----------------

list_versions=list_versions.split()
list_threading=list_threading.split()
list_machines=list_machines.split()


import re # used by grep function
def grep(fileObj, pattern):
  r=[]
  for line in fileObj:
    if re.search(pattern,line):
      r.append(line)
  return r
#end def

class Simulation:
    version=""
    machine=""
    threading=""
    log_file=""
    
    total_time=0.0
    rpn2_time=0.0
    rpt2_time=0.0
    rpn2_calls=0
    rpt2_calls=0
    
    def __init__(self, version, threading, machine):
        self.version=version
        self.threading=threading
        self.machine=machine
        self.log_file="logs/log-" + version + "-" + threading + "." + machine
        
        self.read()
    
    def read(self):
        try:
            in_file = open(self.log_file)
        except:
            print "Warning: File not found:", self.log_file, "Ignoring."
            return
        
        # get wall time
        grep_result = grep(in_file, "Total time:")
        self.total_time = float(grep_result[-1].split()[2])
        
        # get rpn2 stats
        in_file.seek(0)
        grep_result  = grep(in_file, "rpn2 \(calls/time\):")
        self.rpn2_time = float(grep_result[-1].split()[3])
        self.rpn2_calls = int(grep_result[-1].split()[2])
        
        # get rpt2 stats
        in_file.seek(0)
        grep_result  = grep(in_file, "rpt2 \(calls/time\):")
        self.rpt2_time = float(grep_result[-1].split()[3])
        self.rpt2_calls = int(grep_result[-1].split()[2])
    
#end class Simulation


#-----------------------#
# Execution starts here #
#-----------------------#

# list of simulations
sim_list = []

# load all and print table
print "%-30s \ttotal_time \trpn2Time \trpt2_time \trpn2_calls \trpt2_calls" % ("Simulation")
for threading in list_threading:
    for machine in list_machines:
        for version in list_versions:
            sim = Simulation(version, threading, machine)
            sim_list.append(sim)
            print "%-30s \t%f \t%f \t%f \t%10d \t%10d " % (sim.log_file, sim.total_time, sim.rpn2_time, sim.rpt2_time, sim.rpn2_calls, sim.rpt2_calls)
            


