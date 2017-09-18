#!/usr/bin/python

import re # used by grep function
def grep(fileObj, pattern):
  r=[]
  for line in fileObj:
    if re.search(pattern,line):
      r.append(line)
  return r
#end def

class Simulation:
    log_file=""
    machine=""
    num_threads=-1
    version=""
    
    total_time=-1.0
    rpn2_time=-1.0
    rpt2_time=-1.0
    rpn2_calls=-1
    rpt2_calls=-1
    
    def __init__(self, filename):
        self.log_file=filename
        try:
            filename=filename.split('/')[1]
            self.machine=filename.split('_')[0]
            self.num_threads=int(filename.split('_')[1])
            self.version=filename.split('_')[2]
        except:
            print "Warning: something wrong with filename ", self.log_file, ". Ignoring."
            
        self.read()
    
    def read(self):
        try:
            in_file = open(self.log_file)
        except:
            print "Warning: File not found:", self.log_file, ". Ignoring."
            return
        
        try:
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
        except:
            print "Warning: something missing in ", self.log_file, ". Ignoring."
    
#end class Simulation


#-----------------------#
# Execution starts here #
#-----------------------#

# find *.log files in ./logs
import glob
list_files = glob.glob("logs/*.log")
list_files.sort()

    # list of simulations
list_simulations = []

# load all and print table
print "%-30s \ttotal_time \trpn2Time \trpt2_time \trpn2_calls \trpt2_calls" % ("Simulation")
for file in list_files:
    sim = Simulation(file)
    list_simulations.append(sim)
    
for sim in list_simulations:
    print "%-30s \t%f \t%f \t%f \t%10d \t%10d " % (sim.log_file, sim.total_time, sim.rpn2_time, sim.rpt2_time, sim.rpn2_calls, sim.rpt2_calls)

