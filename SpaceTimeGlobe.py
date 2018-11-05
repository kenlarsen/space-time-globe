import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.path as path
import matplotlib.patches as patches
import numpy as np
import sys
from numpy.core.umath import sqrt
import time as time


class LorentzGraph:
    """
        Draw a spacetime diagram
        The diagram has a coordinate system with a time cone, proper time curves and a set of space-time events
        The diagram can be boosted to any fraction of lightspeed (-1..1)
    """

    def lorentz(self, x, t, v):
        'Lorentz transform'
        c = 1.0 # speed of light is 1 for convenience
        gamma = 1 / sqrt(1 - (v*v/c/c))
        xp = gamma * (x - v*np.asarray(t))
        tp = gamma * (t - v*np.asarray(x)/(c*c))
        return xp,tp

    def __init__(self, ax, Xmax, Tmax, options):
        self.ax      = ax
        self.Xmax    = Xmax
        self.Tmax    = Tmax
        self.options = options
        self.infoG = ax.text(-self.Xmax+0.5, self.Tmax-0.5, "", fontsize=20)
        self.meshVerts, self.meshCodes = self.__generatemesh()
        self.zmeshVerts, tmp = self.__generatemesh()
        pth = path.Path(self.meshVerts, self.meshCodes)
        patch = patches.PathPatch(pth, facecolor='none', edgecolor='blue', lw=1)
        self.ax.add_patch(patch)
        self.meshG = patch

    def __drawcoord(self):
        'draw the spacetime coordinate system'
    
        # the axis, grid, labels etc.
        self.ax.set(xlim=(-self.Xmax,self.Xmax), ylim=(-1, self.Tmax))
        self.ax.set_xticks(np.arange(-self.Xmax, self.Xmax, 1))
        self.ax.set_yticks(np.arange(-1, self.Tmax, 1))
        self.ax.set_title(title) # "Lorentz Transformation")
        self.ax.set_xlabel("distance")
        self.ax.set_ylabel("time")
        self.ax.grid(True, ls = '-', color = '#a0a0a0', which = 'both')   # turn on grid lines
    
        # light cone
        self.ax.plot([-1,self.Xmax],[-1,self.Xmax], 'g,-', c='black', lw=2)
        self.ax.plot([-self.Xmax,1],[self.Xmax,-1], 'g,-', c='black', lw=2)
        
        # regular proper time curves
        self.__drawpropertimes(np.linspace(1, self.Tmax, self.Tmax))
        
        self.ax.legend()
        
    def __propertime(self, x, t):
        return sqrt(t*t - x*x)

    def __drawpropertime(self, s, color, lw, ls):
        if s > 0:
            x = np.linspace(0, Xmax, 40)
            t = sqrt(s*s + x*x)
            self.ax.plot(x, t, 'b,:', c=color, linewidth=lw, markersize=0, linestyle=ls)
            self.ax.plot(-x, t, 'b,:', c=color, linewidth=lw, markersize=0, linestyle=ls)

    def __drawpropertimes(self, s):
        if not options["pt"]: return
        for r in s:
            self.__drawpropertime(r, '#bbbbbb', 1, '-')

    def __generatemesh(self):
        x = range(-self.Xmax, self.Xmax+1)
        t = range(-1, self.Tmax+1)
        verts = [(i,j) for i in x for j in [-1, self.Tmax]]
        verts.extend([(j,i) for i in t for j in [-self.Xmax, self.Xmax]])
        verts.append((0,0)) # closepoly
        codes = [path.Path.MOVETO, path.Path.LINETO] * (len(verts)//2);
        codes.append(path.Path.CLOSEPOLY)

        # verts must be a mutable array in  order to perform animation 
        averts = np.zeros((len(verts), 2)) # 
        averts[:] = verts
        return averts, codes
        
    def __updatemesh(self, v):
        if self.options["grid"]:
            for i in range(0, len(self.meshVerts)):
                self.meshVerts[i] = self.lorentz(self.zmeshVerts[i][0], self.zmeshVerts[i][1], v)
        else:
            for i in range(0, len(self.meshVerts)):
                self.meshVerts[i] = (0,0)
        return self.meshG

    def propertimeevent(self, evt):
        'return the proper time length of a given event'
        dx = self.events[evt]["x"][1] - self.events[evt]["x"][0]
        dt = self.events[evt]["t"][1] - self.events[evt]["t"][0]
        s  = self.__propertime(dx, dt)
        return s
    
    def addevents(self, events):
        'add a set of space-time events (name=name,x=x,t=t,c=color)'
        self.events  = events
        
        # add some derived event info
        for e in self.events: 
            e["dx"] = e["x"][1]-e["x"][0]
            e["dt"] = e["t"][1]-e["t"][0]
            e["s"]  = self.__propertime(e["dx"], e["dt"])
        
        # draw each event as a line with markers + the simultaneity lines at the ends
        self.eventsG = []
        i = 0
        for e in self.events:
            self.eventsG.append( self.ax.plot([], [], e["c"], lw=3, markersize=10, 
                                              label="%d %s: pt=%0.2f" % (i+1, e["name"], e["s"]))[0] ) # event
            self.eventsG.append( self.ax.plot([], [], e["c"], lw=3, linestyle='dotted')[0] ) # sim line at event start
            self.eventsG.append( self.ax.plot([], [], e["c"], lw=3, linestyle='dotted')[0] ) # sim line at event end
            i += 1
        
        # draw proper time of events
        for e in self.events: 
            self.__drawpropertime(e["s"], 'black', 3, ':')
        
    def __updateevents(self, v):
        # draw Lorentz transformed events
        i = 0
        for e in self.events:
            x,t = self.lorentz(e["x"], e["t"], v)
            self.eventsG[i].set_data(x, t) # the event
            dx,dt = self.lorentz(e["dx"], e["dt"], v)
            if options["sim"]: # draw sim lines
                for j in [0,1]:
                    slope = dx / dt # sim line slope
                    self.eventsG[i+j+1].set_data([-Xmax,Xmax], [t[j] + (-Xmax-x[j])*slope, t[j] + (Xmax-x[j])*slope])
            else:
                for j in [0,1]:
                    self.eventsG[i+j+1].set_data([],[]) # does not work ???
            i += 3
        return self.eventsG

    def boostevent(self, evt):
        'return the boost required to set a given event as the reference'
        dx = self.events[evt]["x"][1] - self.events[evt]["x"][0]
        dt = self.events[evt]["t"][1] - self.events[evt]["t"][0]
        boost = dx/dt
        print("Boost: ",dx,"/",dt)
        return boost
    
    def __updateinfo(self, v):
        self.infoG.set_text("Boost=%.2f" % (v))
        return self.infoG
    
    def init(self):
        'initialise everything'
        self.__drawcoord()
        
    def update(self, v):
        'update and return all dynamic graphics'
        return (self.__updateinfo(v), *self.__updateevents(v), self.__updatemesh(v))


class Booster:
    """
        Manages easing the boost factor for animation
    """
    
    def __init__(self):
        self.boost = self.boostStart = self.boostLen = self.boostN = 0.0
        
    def setBoost(self, newboost):
        'Set a target boost to be eased into'
        self.boostStart = self.boost
        self.boostLen   = newboost - self.boostStart
        self.boostN     = 0
        
    def setBoostNow(self, newboost):
        'Set the boost immediately'
        self.boostStart = self.boost = newboost
        self.boostLen   = 0
        self.boostN     = 0
        
    def getBoost(self):
        return self.boost
    
    def __easeInOut(self, t):
        if t < 0.5:
            return 2 * t * t
        return (-2 * t * t) + (4 * t) - 1

    def next(self):
        if self.boostN < 10:
            self.boostN += 1
            self.boost = self.boostStart + self.__easeInOut(self.boostN / 10)*self.boostLen
        return self.boost


class UserInput:
    """
        Handles user input from the mouse and keyboard
    """

    def __init__(self, boost, ax, options):
        self.boost   = boost
        self.ax      = ax
        self.options = options
        self.press = None
        ax.figure.canvas.mpl_connect('button_press_event'  , self.on_press)
        ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        ax.figure.canvas.mpl_connect('motion_notify_event' , self.on_motion)
        ax.figure.canvas.mpl_connect('key_press_event'     , self.on_keypress)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != ax: return
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
        self.press = event.xdata, event.ydata, boost.getBoost()
    
    def on_release(self, event):
        'on button release we clear the press data'
        if event.inaxes != ax: return
        print('%s release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
        self.press = None
    
    def on_motion(self, event):
        'on motion we will drag the boost'
        if self.press is None: return
        if event.inaxes != self.ax: return
        xpress, ypress, bpress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        v = bpress - dx/20.0 # boost delta if proportional to dx/20
        if v > 1.0: v = 1.0
        if v < -1.0: v = -1.0
        boost.setBoostNow(v)
#         print('xpress=%f, event.xdata=%f, dx=%f boost=%f' %
#               (xpress, event.xdata, dx, boost.getBoost())
    
    def on_keypress(self, event):
        'on key press we react accordingly'
#         print('keypress: key=%s' % (event.key))
        if event.key == "'": 
            options["grid"] = not options["grid"]
            print("grid:",options["grid"])
            return
        if event.key == ';': 
            options["sim"] = not options["sim"]
            print("sim:",options["sim"])
            return
        if event.key == '.': 
            options["pt"] = not options["pt"]
            print("pt:",options["pt"])
            return
        if event.key == 'escape': sys.exit()
        if   event.key == '1': evt = 0
        elif event.key == '2': evt = 1
        elif event.key == '3': evt = 2
        else: return
        boost.setBoost(lorentz.boostevent(evt))
       
#
# defaulTs
#    
params = {'axes.labelsize' : 15,
          'axes.titlesize' : 25,
          'legend.fontsize': 20,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': False}
plt.rcParams.update(params);
options = {'grid':True, 'sim':True, 'pt':True}

# Define some space-time events
#events = [ {"x":[0,0],"t":[0,6],"c":'bo-'}, {"x":[0,2],"t":[0,3],"c":'ro-'}, {"x":[2,0],"t":[3,6],"c":'go-'} ]
title = 'Twin Paradox'
events = [ {"name":"earth"    ,"x":[0,0],"t":[0,10],"c":'bo-'}, 
           {"name":"twin out" ,"x":[0,4],"t":[0,5],"c":'ro-'}, 
           {"name":"twin back","x":[4,0],"t":[5,10],"c":'go-'} ]

# set the boundaries to reflect the extend of the events
Xmax = max([j for i in [events[k]["x"] for k in range(0,3)] for j in i])*2
Tmax = max([j for i in [events[k]["t"] for k in range(0,3)] for j in i])*2
if Tmax > Xmax*3/4: Xmax = Tmax*3//4
print(Xmax, Tmax)

# allocate the plot area
fig, ax = plt.subplots(figsize=(Xmax*2+1,Tmax+2), dpi=30)

# Create Lorentz graphics, the booster handler and the user input
lorentz = LorentzGraph(ax, Xmax, Tmax, options)
boost   = Booster()
userio  = UserInput(boost, ax, options)

# Add the space-time events
lorentz.addevents(events)

# Animate
def init():
    print("init")
    lorentz.init()
    return lorentz.update(boost.getBoost())
 
lastboost = -1  
laststate = [] 
def animate(frame):
    global lastboost, laststate
    if boost.next() != lastboost:
        lastboost = boost.getBoost();
        laststate = lorentz.update(boost.getBoost())
        return laststate
    time.sleep(0.25) # save CPU from running berserk when we have nothing to animate
    return laststate

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=1000, interval=50, blit=True)

plt.show()
