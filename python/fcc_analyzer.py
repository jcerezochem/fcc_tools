#!/usr/bin/python
"""
DESCRIPTION
-----------
Analyze FCclasses transitions stored in fort.21, by plotting the sticks
by classes and indicating the information of the transition. Picking on
a stick, it raises the transition description on the plot and more info
on the console. The TI convoluted spectrum is available also read from  
fort.18 and plotted with the intensity scaled to the stick range.

Instructions:
-Get info from a transition: right-mouse-click on a stick
-Set the next/previous indexed transiton: +/- keys
-Place a label: left-mouse-click on a stick
-Move a label: hold the label with the left-mouse-button
-Remove one label: on the label, right-mouse-click
-Deactivate a class: mouse-click on the legend
-Clean info about transitons: left-mouse-click on the title
-Clean all labels: right-mouse-click on the title
-Export to xmgrace: central-mouse-click on the title

HISTORY
-------
v0              : First stable version
v0.1  (14/12/2014): Improved layout
v0.1-1(14/12/2014): Using 'fort.18' instead of 'spec_Int_TI.dat'
                  tidy up: Function onpick_legend put out on main
v0.1-2(08/12/2014): Merge analyzer and labeler (exploiting that mouse
                  have three buttons...)
v0.1-3(22/04/2015): Solved some bugs on fcc_analyzer: now properly reads 
                    incomplete fort.21 and fort.18 with negative intensities
v0.2  (13/04/2015): Added export_xmgrace function to export the plot to xmgrace
                    including labels

NOTES AND BUGS
---------------
* Tested on Python 2.7.5+ (from ubuntu repo). Will not work on python3 
  Required packages:
   - numpy (tested with version: 1.9.1)
   - matplotlib (tested with version: 1.4.2)

* If the 0-0 transion has an intensity above 10, it is not printed 
  correctly by FCclasses (format issue). In such a case, an error
  is raised. It can be solved by manually setting the right value
  on fort.21 (getting it, e.g., from fort.8)
"""
import numpy as np

def helptext():
    print """
    Description:
    ------------
    This python scripts parses fort.21 retrieving the informaition about the 
    transitions and plot them using matplotlib. An interactive plot is generated
    which allow the used to assign the stick transition by clicking on them.
    The script can also export the plot to xmgrace format.
    
    Instructions:
    ------------
    * Interacting with the plot:
      -Get info from a transition: right-mouse-click on a stick
      -Set the next/previous indexed transiton: +/- keys
      -Place a label: left-mouse-click on a stick
      -Move a label: hold the label with the left-mouse-button
      -Remove one label: on the label, right-mouse-click
      -Deactivate a class: mouse-click on the legend
      -Clean info about transitons: left-mouse-click on the title
      -Clean all labels: right-mouse-click on the title
      -Export to xmgrace: central-mouse-click on the title
    
    * Command line flags
      By default the script plots the transitions for all classes in absolute intensity
      but this behaviour can be tuned with command line flags
      
       Flag      Action
      ---------------------
       -fc       Plot FCfactors
       -fc-abs   Plot absolute value of FCfactors
       -fc-sqr   Plot square values of FCfactors
       -maxC     Maximum Class to show
       -v        Print version info and quit
       -h        This help
       
    """

class spectral_transition:
    """
    A class for spectral transitions
    """
    def __init__(self):
        self.motherstate = 0
        self.fcclass = 0
        self.position = 0.0
        self.fcfactor = 0.0
        self.intensity = 0.0
        self.einit = 0.0
        self.efin = 0.0
        self.DE = 0.0
        self.DE00cm = 0.0
        self.init = [0]
        self.final = [0]
        self.qinit = [0]
        self.qfinal = [0]
        self.index = 0
    def def_transitions(self):
        #DEFINE TRANSITIONS
        #initial modes
        modesI=''
        for i in range(0,self.init.count(0)):
            self.init.remove(0)
        if len(self.init) > 0:
            for i in range(0,len(self.init)):
                modesI = modesI+str(self.qinit[i])+'('+str(self.init[i])+'),'
            #Remove trailing comma
            modesI = modesI[0:-1]
        else:
            modesI = '0'
        #final modes
        modesF=''
        for i in range(0,self.final.count(0)):
            self.final.remove(0)
        if len(self.final) > 0:
            for i in range(0,len(self.final)):
                modesF = modesF+str(self.qfinal[i])+'('+str(self.final[i])+'),'
            #Remove trailing comma
            modesF = modesF[0:-1]
        else:
            modesF = '0'
        #Define attribute transition
        self.transition = modesI+" --> "+modesF
    def summary(self):
        print "\nINDEX:        %s" % self.index
        print "Mother State: %s" % self.motherstate
        print "FC class:     %s" % self.fcclass
        print "Einitial:     %s" % self.einit
        print "Efinal:       %s" % self.efin
        print "DE:           %s" % self.DE
        print "Intensity:    %s" % self.intensity
        print "FCfactor:     %s" % self.fcfactor
        self.def_transitions()
        print "Transition:   %s\n" % self.transition


# Class PointBrowser will take care of 
        
class PointBrowser:
    """
    DESCRIPTION
    ------------
    Functions to highlight a transition and identify it with the init->final modes
    
    NOTES
    -----
    Based on matplotlib example: data_browser.py
    http://matplotlib.org/examples/event_handling/data_browser.html
    
    *Comments on the original class:
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    
    *Updates in this addaptation
    This class handles vlines (not plot) to represent the stick spectrum
    and we use +/- instead of n/p
    """
    def __init__(self):
        self.lastind = 0

        #tranform: set text coordinates to viewport(grace analog) instead of world, which is the default
        self.text = ax.text(0.5, 0.95, '',
                            transform=ax.transAxes, va='top')
        
        self.selected  = ax.vlines([ener[0]], [zero[0]], [intens[0]], linewidths=3,
                                  color='yellow', visible=False)

    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('+', '-'): return
        if event.key=='+': inc = 1
        else:  inc = -1

        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(ener)-1)
        self.update()

    def onpick(self, event):
        if event.mouseevent.button != 3: return True
        picked_stick = False
        for spectrum in stickspc:
            if event.artist==spectrum: picked_stick = True
        if not picked_stick: return True
        
        """
        Since array indexes for each artist (sticksX) do not correspond to
        the whole array index (becouse we are using parts of the whole array)
        we need to readjust the index for each one
        We use the ind[], a dictionary that associates the index corresponding 
        to the begining of each class.
        """
        indplus = ind[event.artist]

        N = len(event.ind)
        if not N: return True

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        #print "TESTING:"
        #print x, y, event.ind
        #print ener[event.ind], intens[event.ind], zero[event.ind]

        # Compute the distance between the clicked point (x,y) and 
        # the 
        distances = np.hypot(x-ener[event.ind], y-intens[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin] + indplus

        self.lastind = dataind
        self.update()

    def update(self):
        if self.lastind is None: return

        dataind = self.lastind

        self.selected.set_visible(False)
        
        self.selected  = ax.vlines([ener[dataind]], [zero[dataind]], [intens[dataind]], linewidths=3,
                                  color='yellow', visible=True, alpha=0.7)

        #Summary on console (also defines the transition)
        tr[dataind].summary()

        self.text.set_text('Index: %d\nTransition: %s'%((dataind+1),tr[dataind].transition))
        
        
        fig.canvas.draw()
        
    def reset(self,event):
        # Clean the view when clicking on the title
        if event.artist!=title: return True
        self.text.set_text('')
        self.selected.set_visible(False)
        fig.canvas.draw()
        
class LabelSet:
    """
    DESCRIPTION
    ------------
    Functions to create and manage labels over the sticks
    
    NOTES
    -----
    Based on matplotlib example: looking_glass.py
    http://matplotlib.org/examples/event_handling/looking_glass.html
    """
    def __init__(self):
        self.pressevent = None
        self.x0 = ener[0]
        self.y0 = intens[0]
        self.lab = ax.annotate('',xy=(0.,0.))
        return
    
    def onpick(self, event):
        if event.mouseevent.button != 1: return True
        picked_stick = False
        for spectrum in stickspc:
            if event.artist==spectrum: picked_stick = True
        if not picked_stick: return True
        
        """
        Since array indexes for each artist (sticksX) do not correspond to
        the whole array index (because we are using parts of the whole array)
        we need to readust the index for each one
        """
        indplus = ind[event.artist]

        N = len(event.ind)
        if not N: return True

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        #print "TESTING:"
        #print x, y, event.ind
        #print ener[event.ind], intens[event.ind], zero[event.ind]

        distances = np.hypot(x-ener[event.ind], y-intens[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin] + indplus

        self.lastind = dataind
        #self.update()
        self.put_label()
   
    def put_label(self):
        if self.lastind is None: return

        dataind = self.lastind

        #Summary
        tr[dataind].summary()
        #Now set data and label positions
        xd = ener[dataind]
        yd = intens[dataind]
        xl = ener[dataind] - 0.00
        #Since the intensity may change several orders of magnitude, 
        # we better don't shift it
        yl = intens[dataind] + 0.0
        
        #Define label as final modes description
        if dataind == 0:
            label='0-0'
            agrlabel='0-0'
        else:
            label='$'
            agrlabel=''
            for i in range(0,len(tr[dataind].final)):
                # Labels changed to Mode^quata
                label = label+str(tr[dataind].final[i])+'^{'+str(tr[dataind].qfinal[i])+'},'
                agrlabel = agrlabel+str(tr[dataind].final[i])+'\S'+str(tr[dataind].qfinal[i])+'\N,'
            #Remove trailing comma
            label = label[0:-1]+'$'
            agrlabel = agrlabel[0:-1]

        #In labelref we get the annotation class corresponding to the 
        #current label. labelref is Class(annotation)
        labelref = ax.annotate(label, xy=(xd, yd), xytext=(xl, yl),picker=1,
                            arrowprops=dict(arrowstyle="-",
                                            color='grey'))

        #Check whether the label was already assigned or not
        set_lab = True
        for labref in labs:
            lab = labs[labref]
            if lab == agrlabel:
                print "This label was already defined"
                set_lab = False
        if set_lab:
            # The dictionary labs relates each labelref(annotation) to 
            # the agrlabel. Note that the mpl label can be retrieved
            # from labelref.get_name()
            labs[labelref] = agrlabel
            fig.canvas.draw()
        else:
            labelref.remove()
    
    def getlabel(self,event):
        picked_lab = False
        for labref in labs:
            if event.artist==labref: picked_lab = True
        if not picked_lab: return True
    
        self.lab = event.artist
        self.x0, self.y0 = self.lab.get_position()
    
    def onpress(self,event):
        if event.button != 1: return
        if event.inaxes!=ax:
         return
        if not self.lab.contains(event)[0]:
            return
        self.pressevent = event
    
    def onrelease(self,event):
        if event.button != 1: return
        if event.inaxes!=ax:
         return
        self.pressevent = None

    def onmove(self,event):
        if event.button != 1: return
        if self.pressevent is  None or event.inaxes!=self.pressevent.inaxes: return
    
        dx = event.xdata - self.pressevent.xdata
        dy = event.ydata - self.pressevent.ydata
        self.lab.set_position((self.x0+dx,self.y0+dy))
        fig.canvas.draw() 
        
    def delete(self,event):
        if event.button != 3: return
        if not self.lab.contains(event)[0]:
            return
        #Remove the label
        self.lab.remove()
        #And substract the corresponding entry from the dict
        labs.pop(self.lab)
        #We deactivate the lab. Otherwise, if we reclikc on the same point, it raises an error
        self.lab.set_visible(False)
        fig.canvas.draw()
        
    def reset(self,event):
        # Clean the view when clicking on the title
        if event.artist!=title: return True
        if event.mouseevent.button!=3: return True
        #We need a local copy of labs to iterate while popping
        labs_local = dict.copy(labs)
        for lab in labs_local:
            lab.remove()
            labs.pop(lab)
        fig.canvas.draw()
        
def onpick_legend(event):
    """
    DESCRIPTION
    ------------
    Function to activate/deactivate a plot, clicking on the legend 
    
    NOTES
    -----
    Based on matplotlib example: legend_picking.py
    http://matplotlib.org/examples/event_handling/legend_picking.html
    """
    # Avoid pick even when clicking on the stick spectrum
    picked_leg = False
    for legline in leg.get_lines():
        if event.artist==legline: picked_leg = True
    if not picked_leg: return True
    # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
    legline = event.artist
    origline = lined[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)
    # Change the alpha on the line in the legend so we can see what lines
    # have been toggled
    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()
    
    
def export_xmgrace(event):
    """
    DESCRIPTION
    ------------
    Function to convert the current data into a xmgrace plot, including
    labels currently on the screen
    
    Public variables taken from __main__
    * labs,  dict     : labels in xmgr format (arg: labels as annotation Class)
    * ind,   dict     : array indices (arg: stick spectra as 
                        matplotlib.collections.LineCollection Class)
    * ax,    mpl.axes : graphic info
    * ener,  np.array : vector with the energy of all transitions
    * intens,np.array : vector with the intensity of all transitions
    * spc,   np.array : convoluted spectrum. x is spc[:,0], y is spc[:,1]
    """

    if event.artist!=title: return True
    if event.mouseevent.button!=2: return True

    print "\nExporting plot to xmgrace (fcc_analyzer.agr)..."

    f = open('fcc_analyzer.agr','w')
    
    print >> f, "# XMGRACE CREATED BY FCC_ANALYZER"
    print >> f, "# Only data and labels. Format will"
    print >> f, "# be added by your default xmgrace"
    print >> f, "# defaults (including colors, fonts...)"
    print >> f, "# Except the followins color scheme:"
    print >> f, '@map color 0 to (255, 255, 255), "white"'
    print >> f, '@map color 1 to (0, 0, 0), "black0"'
    print >> f, '@map color 2 to (0, 0, 0), "black"'
    print >> f, '@map color 3 to (0, 0, 255), "blue"'
    print >> f, '@map color 4 to (255, 0, 0), "red"'
    print >> f, '@map color 5 to (0, 139, 0), "green4"'
    print >> f, '@map color 6 to (0, 255, 255), "cyan"'
    print >> f, '@map color 7 to (255, 0, 255), "magenta"'
    print >> f, '@map color 8 to (255, 255, 0), "yellow"'
    print >> f, '@map color 9 to (188, 143, 143), "brown"'
    print >> f, '@map color 10 to (220, 220, 220), "grey"'
    print >> f, '@map color 11 to (0, 255, 0), "green"'
    print >> f, '@map color 12 to (148, 0, 211), "violet"'
    print >> f, '@map color 13 to (255, 165, 0), "orange"'
    print >> f, '@map color 14 to (114, 33, 188), "indigo"'
    print >> f, '@map color 15 to (103, 7, 72), "maroon"'
    print >> f, '@map color 16 to (64, 224, 208), "turquoise"'
    print >> f, '@map color 17 to (50, 50, 50), "gris2"'
    print >> f, '@map color 18 to (100, 100, 100), "gris3"'
    print >> f, '@map color 19 to (150, 150, 150), "gris4"'
    print >> f, '@map color 20 to (200, 200, 200), "gris5"'
    print >> f, '@map color 21 to (255, 150, 150), "red2"'
    print >> f, '@map color 22 to (150, 255, 150), "green2"'
    print >> f, '@map color 23 to (150, 150, 255), "blue2"'  
    # Without the @version, it makes auto-zoom (instead of taking world coords) 
    print >> f, "@version 50123"
    print >> f, "@page size 792, 612"
    print >> f, "@default symbol size 0.010000"
    print >> f, "@default char size 0.800000"
    for lab in labs:
        print >> f, "@with line"
        print >> f, "@    line on"
        print >> f, "@    line g0"
        print >> f, "@    line loctype world"
        print >> f, "@    line color 20"
        print >> f, "@    line ",lab.xy[0],",",lab.xy[1],",",lab.xyann[0],",",lab.xyann[1]
        print >> f, "@line def"
        print >> f, "@with string"
        print >> f, "@    string on"
        print >> f, "@    string g0"
        print >> f, "@    string loctype world"
        print >> f, "@    string ", lab.xyann[0],",",lab.xyann[1]
        print >> f, "@    string def \"",labs[lab],"\""
    print >> f, "@with g0"
    # Set a large view
    print >> f, "@    view 0.150000, 0.150000, 1.2, 0.92"
    #Get plotting range from mplt
    x=ax.get_xbound()
    y=ax.get_ybound()
    print >> f, "@    world ",x[0],",",y[0],",",x[1],",",y[1]
    #Get xlabel from mplt
    print >> f, "@    xaxis  label \""+ax.get_xlabel()+"\""
    #Get tick spacing from mplt
    x=ax.get_xticks()
    y=ax.get_yticks()
    print >> f, "@    xaxis  tick major", x[1]-x[0]
    print >> f, "@    yaxis  tick major", y[1]-y[0]
    #Legend
    print >> f, "@    legend loctype view"
    print >> f, "@    legend 0.95, 0.9"
    #Now include data
    counter=-1
    if (spc[:,:].size != 0):
        counter+=1
        print >> f, "# Spect"
        print >> f, "@    s"+str(counter),"line type 1"
        print >> f, "@    s"+str(counter),"line linestyle 3"
        print >> f, "@    s"+str(counter),"legend  \"Spec\""
        for i in range(0,spc[:,0].size):
            print >> f, spc[i,0], spc[i,1]
    if True:
        counter+=1
        print >> f, "& C0"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"0-0\""
        print >> f, ener[0], intens[0]
    if nclass1 != 0:
        counter+=1
        print >> f, "& C1"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C1\""
        for i in range(ind[sticks1],ind[sticks1]+nclass1):
            print >> f, ener[i], intens[i]
    if nclass2 != 0:
        counter+=1
        print >> f, "& C2"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C2\""
        for i in range(ind[sticks2],ind[sticks2]+nclass2):
            print >> f, ener[i], intens[i]
    if nclass3 != 0:
        counter+=1
        print >> f, "& C3"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C3\""
        for i in range(ind[sticks3],ind[sticks3]+nclass3):
            print >> f, ener[i], intens[i]
    if nclass4 != 0:
        counter+=1
        print >> f, "& C4"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C4\""
        for i in range(ind[sticks4],ind[sticks4]+nclass4):
            print >> f, ener[i], intens[i]
    if nclass5 != 0:
        counter+=1
        print >> f, "& C5"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C5\""
        for i in range(ind[sticks5],ind[sticks5]+nclass5):
            print >> f, ener[i], intens[i]
    if nclass6 != 0:
        counter+=1
        print >> f, "& C6"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C5\""
        for i in range(ind[sticks6],ind[sticks6]+nclass6):
            print >> f, ener[i], intens[i]
    if nclass7 != 0:
        counter+=1
        print >> f, "& C7"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"C7\""
        for i in range(ind[sticks7],ind[sticks7]+nclass7):
            print >> f, ener[i], intens[i]
    if nhot != 0:
        counter+=1
        print >> f, "& Hot"
        print >> f, "@type bar"
        print >> f, "@    s"+str(counter),"line type 0"
        print >> f, "@    s"+str(counter),"legend  \"Hot\""
        for i in range(ind[sticksH],ind[sticksH]+nhot):
            print >> f, ener[i], intens[i]
            
    f.close()
    print "Done\n"


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import sys
    try:
        import version_tag
    except:
        class version_tag:
            COMMIT="Untracked"
            DATE="No date"
    
    
    print """
    =====================================================
                     FCC_ANALYZER:                  
          a python tool to analyze FCCLASSES (TI) output
          
          VERSION INFO:
            Git commit: %s
            Date      : %s
    =====================================================
    """%(version_tag.COMMIT,version_tag.DATE)
    
    #====================================================
    # Command line options
    #====================================================
    # Defaults
    plot_intensity=True
    plot_fcabs=False
    plot_fcsqr=False
    skip_next=False
    MaxClass=7
    for i,option in enumerate(sys.argv):
        if skip_next:
            skip_next=False
        elif option == "-fc":
            plot_intensity=False
        elif option == "-fc-abs":
            plot_intensity=False
            plot_fcabs=True
        elif option == "-fc-sqr":
            plot_intensity=False
            plot_fcsqr=True
        elif option == "-maxC":
            MaxClass = int(sys.argv[i+1])
            skip_next = True
        elif option == "-v":
            sys.exit()
        elif option == "-h":
            helptext()
            sys.exit()
        elif option == sys.argv[0]:
            pass
        else:
            exit("ERROR: Unknown command line option: %s"%(option))

    #====================================================
    #Analyze fort.21 (this is the core of the script)
    #====================================================
    tr=[]
    print "Loading transitions (fort.21)..."
    try:
        f = open('fort.21','r')
    except:
        try:
            f = open('Assignments.dat','r')
        except:
            exit("ERROR: Cannot open file 'fort.21' or 'Assignments.dat'")
        
    #First read 0-0 transition
    itrans = 0
    for line in f:
        if "INDEX" in line:
            line = f.next()
            data = line.split()
            tr.append(spectral_transition())
            tr[itrans].motherstate = 1
            tr[itrans].fcclass     = 0
            tr[itrans].index       = float(data[0])
            tr[itrans].einit       = float(data[1])
            tr[itrans].efin        = float(data[2])
            tr[itrans].DE          = float(data[3])
            tr[itrans].DE00cm      = float(data[4])
            tr[itrans].fcfactor    = float(data[5])
            """
            Sometimes, intens for 0-0 is ******
            Could be read from fort.8 although 
            the cleanest is to fix it fcclasses
            """ 
            try: tr[itrans].intensity   = float(data[6])
            except: exit("ERROR: Check 0-0 transition") 
            #For this case, modes are assigned manually
            tr[itrans].init = [0] 
            tr[itrans].final = [0] 
            #Now break to continue with other Classes
            break

    nhot = 0
    nclass1 = 0
    nclass2 = 0
    nclass3 = 0
    nclass4 = 0
    nclass5 = 0
    nclass6 = 0
    nclass7 = 0
    for line in f:
        if "MOTHER STATE" in line:
            motherstate = int(line.split("N.")[1])
        elif "C1 " in line:
            fcclass = 1
            if motherstate == 1:
                nclass = 0
        elif "C2 " in line:
            fcclass = 2
            if motherstate == 1:
                nclass1 = nclass
                nclass = 0
        elif "C3 " in line:
            fcclass = 3
            if motherstate == 1:
                nclass2 = nclass
                nclass = 0
        elif "C4 " in line:
            fcclass = 4
            if motherstate == 1:
                nclass3 = nclass
                nclass = 0
        elif "C5 " in line:
            fcclass = 5
            if motherstate == 1:
                nclass4 = nclass
                nclass = 0
        elif "C6 " in line:
            fcclass = 6
            if motherstate == 1:
                nclass5 = nclass
                nclass = 0
        elif "C7 " in line:
            fcclass = 7
            if motherstate == 1:
                nclass6 = nclass
                nclass = 0
        elif "M-0 TRANSITION" in line and motherstate == 2:
            if fcclass == 1:
                nclass1 = nclass
            elif fcclass == 2:
                nclass2 = nclass
            elif fcclass == 3:
                nclass3 = nclass
            elif fcclass == 4:
                nclass4 = nclass
            elif fcclass == 5:
                nclass5 = nclass
            elif fcclass == 6:
                nclass6 = nclass
            elif fcclass == 7:
                nclass7 = nclass
            nclass = 0
            fcclass = 0
        elif ("E+" in line) | ("E-" in line):
            nclass += 1
            if fcclass<=MaxClass:
                data = line.split()
                itrans += 1
                tr.append(spectral_transition())
                tr[itrans].motherstate = motherstate
                tr[itrans].fcclass     = fcclass
                tr[itrans].index       = float(data[0])
                tr[itrans].efin        = float(data[1])
                tr[itrans].einit       = float(data[2])
                tr[itrans].DE          = float(data[3])
                tr[itrans].DE00cm      = float(data[4])
                tr[itrans].fcfactor    = float(data[5])
                tr[itrans].intensity   = float(data[6])
        #Indentify modes in the transition
        #Initial
        elif 'state 1 = GROUND' in line and fcclass<=MaxClass:
            tr[itrans].init = [0]
        elif 'Osc1=' in line:
            A = line.split('Osc1=')
            del A[0]
            for i in range(0,len(A)):
                A[i] = int(A[i])
            tr[itrans].init = A 
        #Number of quanta involved
        elif 'Nqu1=' in line and fcclass<=MaxClass:
            A = line.split('Nqu1=')
            del A[0]
            for i in range(0,len(A)):
                A[i] = int(A[i])
            tr[itrans].qinit = A 
        #Final
        elif 'state 1 = GROUND' in line and fcclass<=MaxClass:
            tr[itrans].final = [0]
        elif 'Osc2=' in line and fcclass<=MaxClass:
            A = line.split('Osc2=')
            del A[0]
            for i in range(0,len(A)):
                A[i] = int(A[i])
            tr[itrans].final = A 
        #Number of quanta involved
        elif 'Nqu2=' in line and fcclass<=MaxClass:
            A = line.split('Nqu2=')
            del A[0]
            for i in range(0,len(A)):
                A[i] = int(A[i])
            tr[itrans].qfinal = A 
    #If there were mother states>1, no nclass7 was already updated, otherwise:
    if tr[itrans].motherstate == 1: nclass7 = nclass
    else: nhot = nclass

    f.close()
    loadC=True
    nclass_list = [nclass1,nclass2,nclass3,nclass4,nclass5,nclass6,nclass7]
    print 'Transitions read:'
    print ' Class     N. trans.         Load?  '
    for i,nclass in enumerate(nclass_list):
        if MaxClass<i+1: loadC=False
        print ' C{0}        {1:5d}             {2}   '.format(i+1,nclass,loadC)
        if MaxClass<i+1: nclass_list[i]=0
    print     ' Hot       {0:5d}                   '.format(nhot)
    print ''
    print 'Loaded transitions: ',(itrans)
    nclass1,nclass2,nclass3,nclass4,nclass5,nclass6,nclass7 = nclass_list
    #========== Done with fort.21 ====================================

    #==================================================
    #Now construct the spectra and plot them
    #==================================================
    #Plot specifications (labels and so)
    fig, ax = plt.subplots()
    title = ax.set_title('TI stick spectrum from $\mathcal{FC}classes$',fontsize=18,picker=1)
    ax.set_xlabel('Energy (eV)',fontsize=16)
    ax.set_ylabel('Inensity',fontsize=16)
    ax.tick_params(direction='out',top=False, right=False)
    
    #Inialize variables
    intens   = []
    ener     = []
    stickspc = []
    ind = dict()
    for i in range(0,len(tr)):
        if plot_intensity:
            intens.append(tr[i].intensity)
        elif plot_fcsqr:
            intens.append((tr[i].fcfactor)**2)
        elif plot_fcabs:
            intens.append(abs(tr[i].fcfactor))
        else: #plot FC factors
            intens.append(tr[i].fcfactor)
        ener.append(tr[i].DE)

    # Build the np.arrays needed for plotting
    intens = np.array(intens)
    ener = np.array(ener)
    zero = np.zeros(len(intens))

    #Read fort.18 (from "Total spectrum at the chosen temperature")
    x = []
    y = []
    print "Loading convoluted spectrum (fort.18)..."
    try:
        f = open('fort.18','r')
        read_spc = False
        for line in f:
            #We ensure that it do not stop at intermediate steps
            if "Total spectrum at the chosen temperature up to now" in line:
                continue
            elif "Total spectrum at the chosen temperature" in line:
                read_spc = True
                continue
            #Now read the file
            if read_spc:
                data = line.split()
                if "integral" in line: break
                x.append(float(data[0]))
                y.append(float(data[1]))
        f.close()
        spc = np.array((x,y))
        #Tranpose to use the same index order as in previous version
        spc = spc.transpose()
        try: 
            spcmax = spc[:,1].max()
            spcmin = spc[:,1].min()
            if abs(spcmin) > spcmax:
                spc[:,1] = spc[:,1]/spcmin*intens.min()
            else:
                spc[:,1] = spc[:,1]/spcmax*intens.max()
            #Plot
            spctot, = ax.plot(spc[:,0],spc[:,1],'--',color='k',label="spec")
            stickspc.append(spctot)
        except:
            print "Spectra on fort.18 could not be loaded"
    except:
        print "Spectra on fort.18 could not be loaded"
        spc = np.array((x,y))
        read_spc = False
    
    #Now start plotting sticks (should be done with a list over the classes)
    # NOTE: ind[] stores the array indexes (from 0), not the fort.21 INDEX (from 1)
    #0-0 (always present)
    sticks0 = ax.vlines(ener[0],zero[0],intens[0],linewidths=1,color='k',label="0-0",picker=5)
    stickspc.append(sticks0)
    ind[sticks0] = 0
    indini = 0
    indfin = 1
    if nclass1 != 0:
        #Class 1
        indini = 1
        indfin = 1 + nclass1
        ind1 = indini
        sticks1 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='b',label="C1",picker=5)
        stickspc.append(sticks1)
        ind[sticks1] = indini
    if nclass2 != 0:
        #Class 2
        indini = 1 + nclass1
        indfin = 1 + nclass1 + nclass2
        ind2 = indini
        sticks2 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='r',label="C2",picker=5)
        stickspc.append(sticks2)
        ind[sticks2] = indini
    if nclass3 != 0:
        #Class 3
        indini = 1 + nclass1 + nclass2
        indfin = 1 + nclass1 + nclass2 + nclass3
        ind3 = indini
        sticks3 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='g',label="C3",picker=5)
        stickspc.append(sticks3)
        ind[sticks3] = indini
    if nclass4 != 0:
        #Class 4
        indini = 1 + nclass1 + nclass2 + nclass3
        indfin = 1 + nclass1 + nclass2 + nclass3 + nclass4
        sticks4 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='c',label="C4",picker=5)
        stickspc.append(sticks4)
        ind[sticks4] = indini
    if nclass5 != 0:
        #Class 5
        indini = 1 + nclass1 + nclass2 + nclass3 + nclass4
        indfin = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5
        sticks5 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='m',label="C5",picker=5)
        stickspc.append(sticks5)
        ind[sticks5] = indini
    if nclass6 != 0:
        #Class 6
        indini = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5
        indfin = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5 + nclass6
        sticks6 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='brown',label="C6",picker=5)
        stickspc.append(sticks6)
        ind[sticks6] = indini
    if nclass7 != 0:
        #Class 7
        indini = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5 + nclass6
        indfin = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5 + nclass6 + nclass7
        sticks7 = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='pink',label="C7",picker=5)
        stickspc.append(sticks7)
        ind[sticks7] = indini
    if nhot != 0:
        #Hot bands (we put all them together in the spectra)
        indini = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5 + nclass6 + nclass7
        indfin = 1 + nclass1 + nclass2 + nclass3 + nclass4 + nclass5 + nclass6 + nclass7 + nhot
        sticksH = ax.vlines(ener[indini:indfin],zero[indini:indfin],intens[indini:indfin],linewidths=1,color='orange',label="Hot",picker=5)
        stickspc.append(sticksH)
        ind[sticksH] = indini

    #Legends management
    leg = ax.legend(loc='upper right', fancybox=True, shadow=True)
    #leg.set_bbox_to_anchor([1.1,1.])
    leg.get_frame().set_alpha(0.4)
    #leg.get_frame().set_x(100.)
    #leg.get_frame().set_y(100.)
    leg.set_picker(5)
    # we will set up a dict mapping legend line to orig line, and enable
    # picking on the legend line (from legend_picking.py)
    lined = dict()
    for legline, origline in zip(leg.get_lines(), stickspc):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline
        
    #Also create a dictionary for the labels. Used by both browser and labset
    #Being in __main__, it is public to all Classes and Functions
    labs = dict()
    
    #Instance PointBrowser functions
    browser = PointBrowser()
    labset  = LabelSet()
    
    #Make the connection of events with callbacks
    #Highlight sticks
    fig.canvas.mpl_connect('pick_event', browser.onpick)
    fig.canvas.mpl_connect('pick_event', browser.reset)
    fig.canvas.mpl_connect('key_press_event', browser.onpress)
    #On/Off spectra
    fig.canvas.mpl_connect('pick_event', onpick_legend)
    #Manage labels
    fig.canvas.mpl_connect('pick_event', labset.onpick)
    fig.canvas.mpl_connect('pick_event', labset.getlabel)
    fig.canvas.mpl_connect('button_press_event', labset.onpress)
    fig.canvas.mpl_connect('button_release_event', labset.onrelease)
    fig.canvas.mpl_connect('motion_notify_event', labset.onmove)
    fig.canvas.mpl_connect('button_press_event', labset.delete)
    #Delete all labels
    fig.canvas.mpl_connect('pick_event', labset.reset)
    #Export to xmgrace function
    fig.canvas.mpl_connect('pick_event', export_xmgrace)
    
    plt.show()
    
    #print(intens)

