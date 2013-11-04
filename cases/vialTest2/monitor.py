import os

tspan = 180

class Task(object):
    def __init__(self, name, time, pid):
        self.name = name
        self.time = time
        self.pid = [pid]
        self.nprocs = len(self.pid)
        
    def __eq__(self,other):
        return self.name == other.name and abs(self.time - other.time) < tspan
        
    def __neq__(self,other):
        return not self == other
        
    def addPID(self, newPID, newTime):
        self.time = (self.time*len(self.pid) + newTime) / (len(self.pid)+1)
        self.pid.append(newPID)
        self.nprocs = len(self.pid)

    def __str__(self):
        return "%s: Running for %4.1f hours (%d processors)" % (self.name, self.time/3600., self.nprocs)


f = os.popen('ps -u tvoskuil')
lines = f.readlines()



tasks = []

for i,l in enumerate(lines):
    if i > 0:
        pid = int(l[0:5].strip())
        timeStr = l[15:23]
        time = int(timeStr[0:2])*3600 + int(timeStr[3:5])*60 + int(timeStr[6:8])
        name = l[24:].strip()
        
        newTask = Task(name, time, pid)
        
        if newTask not in tasks:
            tasks.append(newTask)
        else:
            ti = tasks.index(newTask)
            tasks[ti].addPID(pid, time)
        

for t in tasks:
    print t
