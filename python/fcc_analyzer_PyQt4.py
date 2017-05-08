#!/usr/bin/python
# -*- coding: utf8 -*-

"""
This demo demonstrates how to embed a matplotlib (mpl) plot 
into a PyQt4 GUI application, including:
* Using the navigation toolbar
* Adding data to the plot
* Dynamically modifying the plot's properties
* Processing mpl events

* Saving the plot to a file from a menu
The main goal is to serve as a basis for developing rich PyQt GUI
applications featuring mpl plots (using the mpl OO API).
Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 19.01.2009
"""
import sys, os, random
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtCore, QtGui

import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
#from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib
import re

try:
    import version_tag
except:
    class version_tag:
        COMMIT="Untracked"
        DATE="No date"

stick_type = 'fc'

def helptext():
    print """
    Description:
    ------------
    This python script parses fort.21 retrieving the informaition about the 
    transitions and plot them using matplotlib. An interactive plot is generated
    which allow to assign the stick transition by clicking on them.
    The script can also export the plot to xmgrace format. The convoluted spectrum
    is generated from fort.22 (otherwise, it is read from fort.18)
    
    Instructions:
    ------------
    General call:
    
        fcc_analyzer_PyQt4 [options]
    
    The data to be analyzed (fort.21) should be on the folder from where the script is being called
    If no fort.21 file is present, a window will pop-up to select the appropriate path. It the type
    of spectrum was not specified by command line, a new window will pop-up to select the type of
    spectrum.
    
    * Interacting with the plot, once it is loaded:
      -Get info from a transition: right-mouse-click on a stick
      -Seach a transition: use the serach text box
      -Set the next/previous indexed transiton: +/- keys
      -Place a label: left-mouse-click on a stick
      -Move a label: hold the label with the left-mouse-button
      -Remove one label: on the label, right-mouse-click
      -Deactivate a class: mouse-click on the legend
      -Clean info about transitons: push "Clean(panel)" button
      -Clean all labels: push "Clean(labels)" button
      -Export to xmgrace: use File->Export to xmgrace
    """

class SpcConstants:
    exp = {"abs":1,"ecd":1,"emi":3,"cpl":3}
    factor = {"abs":703.30,"ecd":20.5288,"emi":1063.055,"cpl":4252.216}
    # The factors already include the conversion between eV <-> au
    # to handle some issues in the emission Lineshape (need to used
    # (27.2116) factor instead of (27.2116)^3)
    #factor = {"abs":25.8459,"ecd":0.75441,"emi":39.0662,"cpl":156.265}


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
                modesI = modesI+str(self.init[i])+'('+str(self.qinit[i])+'),'
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
                modesF = modesF+str(self.final[i])+'('+str(self.qfinal[i])+'),'
            #Remove trailing comma
            modesF = modesF[0:-1]
        else:
            modesF = '0'
        #Define attribute transition
        return modesI+" --> "+modesF
    
    def info(self):
        transition = self.def_transitions()
        msg = """ Transition:   \t%s 
  =========================
  MotherState:\t%s
  FC class:     \t%s 
  Einit(eV):    \t%s 
  Efin (eV):    \t%s 
  DE   (eV):    \t%s 
  Intensity:    \t%s 
  FCfactor:     \t%s 
  INDEX:        \t%s
        """%(transition,
             self.motherstate,
             self.fcclass,
             self.einit,
             self.efin,
             self.DE,
             self.intensity,
             self.fcfactor,
             self.index)
        return msg

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        """
        The _init__ function takes care of:
        1) Building the app window 
        2) Load data (fort.21, fort.22)
        3) Initialize some variables
        """
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('FCclasses analyzer')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        # Initialize additional items
        # Broadening info
        self.broadening="Gau"
        self.update_hwhm_from_slider(UpdateConvolute=False)
        # Data type
        self.data_type ="Intensity"
        # Line markers
        self.selected  = None
        # Label dictionary
        self.labs = dict()
        # Other values that need initialization
        self.active_label = None
        self.active_tr = None
        self.spectrum_ref = None
        self.spectrum_sim = None
        self.AnnotationType = None
        
        # Get command line arguments
        cml_args = get_args()
        
        # First locate the fort.21 file
        if os.path.isfile('fort.21'):
            path=""
        else:
            file_choices = r"FCclasses aux (fort.21) (fort.21);; All files (*)"
            path = unicode(QFileDialog.getOpenFileName(self, 
                            'Set the location of FClasses output', '', 
                            file_choices))
            if not path or not "fort.21" in path:
                return
            path = path.replace("fort.21","")
        # Parse command line args
        MaxClass = cml_args.get("-maxC")
        MaxClass = int(MaxClass)
        self.spc_type = cml_args.get("-type")
        self.spc_type = str(self.spc_type).lower()
        calc_type_list = ("abs","emi","ecd","cpl")
        if self.spc_type not in calc_type_list:
            # Get spc_type if not given on input or was wrong
            get_option,self.spc_type = self.start_assistant()
            if not get_option:
                sys.exit()
        # Data load
        # Stick transitions (fort.21)
        self.fcclass_list = read_fort21(path+'fort.21',MaxClass)
        # Bins to convolute spectrum (only in new versions of FCclasses)
        if os.path.isfile(path+'fort.22'):
            self.with_fort22 = True
            x,y = read_spc_xy(path+'fort.22')
            # If there are hot bands, the fort.22 file
            # repeats the bins for each MotherState.
            # Here we just want them all together (so as
            # to properly handle the "Input bins" checkbox option
            nMotherStates = x.count(x[0])
            nbins = len(x)/nMotherStates
            self.xbin = np.array(x[:nbins])
            self.ybin = np.array(y[:nbins])
            # Sum all MotherStates into the same bins
            for i in range(nMotherStates-1):
                n = (i+1)*nbins
                self.ybin += np.array(y[n:n+nbins])
            # ybin is in intensity and atomic inits. Changed to "experimental" units
            factor = SpcConstants.factor[self.spc_type]
            self.ybin = self.ybin*factor
        elif os.path.isfile(path+'fort.18'):
            self.with_fort22 = False
            x,y = read_spc_xy(path+'fort.18',fromsection="Total spectrum at the chosen temperature")
            # Deactivate convolution controls
            self.inputBins_cb.setEnabled(False)
            self.broadbox.setText("-")
            self.broadbox.setEnabled(False)
            self.slider.setEnabled(False)
            self.select_broad.setEnabled(False)
        else:
            sys.exit("No convoluted spectrum could be loaded")
        
        # This is the load driver
        self.load_sticks(cml_args.get("-stick"))
        self.set_axis_labels()
        if os.path.isfile(path+'fort.22'):
            self.load_convoluted()
        else:
            self.load_fixconvolution(x,y)
        self.load_legend()        
        
        if cml_args.get("--test"):
            sys.exit()
        
        
    #==========================================================
    # FUNCTIONS CONNECTED WITH THE MENU OPTIONS
    #==========================================================
    def save_plot(self):
        file_choices = r"Portable Network Graphics (*.png) (*.png);; All files (*)"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=100)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
            
    def open_plot(self):
        file_choices = r"Data file (*.dat) (*.dat);; All files (*)"
        
        path = unicode(QFileDialog.getOpenFileName(self, 
                        'Open spectrum', '', 
                        file_choices))
        if path:
            self.statusBar().showMessage('Opened %s' % path, 2000)    
            x,y = read_spc_xy(path)
            x = np.array(x)
            y = np.array(y)
            # Pop-up a selector to indicate units of Xaxis
            ok, units = self.spc_import_assitant_Xaxis()
            if not ok:
                self.statusBar().showMessage('Load data aborted', 2000)
                return
            # Transform X axis if needed
            if units == "cm^-1":
                x = x * 1.23981e-4
            elif units == "nm":
                x = 1239.81/x
                
            # Pop-up a selector to indicate units of Yaxis
            ok, data_type = self.spc_import_assitant_Yaxis()
            if not ok:
                self.statusBar().showMessage('Load data aborted', 2000)
                return
            if data_type != self.data_type:
                # If data_type is not the same as the graph, transform the data
                n = SpcConstants.exp[self.spc_type]
                factor = SpcConstants.factor[self.spc_type]
                if self.data_type == "Lineshape":
                    # Division x/27.2116 could be included in the factor
                    y /= (x/27.2116)**n * factor
                elif self.data_type == "Intensity":
                    # Division x/27.2116 could included in the factor
                    y *= (x/27.2116)**n * factor
            
            # Update Table
            self.refspc_table.setItem(1,1, QTableWidgetItem(path.split('/')[-1]))
            cell = self.refspc_table.item(1,1)
            cell.setTextColor(Qt.black)
            cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
                
            self.load_experiment_spc(x,y)
            
            # Initialize shift an scale
            self.ref_shift = 0.0
            self.ref_scale = 1.0
            self.refspc_table.setItem(2,1, QTableWidgetItem(str(self.ref_shift)))
            self.refspc_table.setItem(3,1, QTableWidgetItem(str(self.ref_scale)))
            
            # Enable manipulations
            self.shiftref_action.setEnabled(True)
            self.scaleref_action.setEnabled(True)


    #==========================================================
    # POP-UP WINDOWS
    #==========================================================
    def start_assistant(self):
        calc_type_list = ("abs","emi","ecd","cpl")
                 
        spc_type, ok = QInputDialog.getItem(self, "Select type of calculation", 
                            "Type of calculation", calc_type_list, 0, False)
        
        return ok, str(spc_type)
    
    
    def spc_import_assitant_Xaxis(self):
        unit_list = ("eV", "cm^-1", "nm")
                 
        units, ok = QInputDialog.getItem(self, "Import Assistant", 
                    "X-Units", unit_list, 0, False)
                         
        if not ok or not units:
            units = ""
        
        return ok, units
        
        
    def spc_import_assitant_Yaxis(self):
        data_type_list = ("Intensity","Lineshape")
                 
        data_type, ok = QInputDialog.getItem(self, "Import Assistant", 
                            "Data type", data_type_list, 0, False)
                         
        if not ok or not data_type:
            data_type = ""
                
        return ok, data_type
        
            
    def xmgr_export(self):
        # Export mode Dialog
        export_mode_list = ("Overlaid graphs","Same graph")
                 
        export_mode, ok = QInputDialog.getItem(self, "Select export mode", 
                            "Organization of the plots", export_mode_list, 0, False)
        if not ok:
            return
        # File Dialog
        file_choices = "xmgrace graph (*.agr) (*.agr);; All files (*)"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Export to file', '', 
                        file_choices))
        if path:
            if self.spectrum_ref:
                spc = self.spectrum_sim+self.spectrum_ref
            else:
                spc = self.spectrum_sim
            if export_mode == "Overlaid graphs":
                ax2 = self.axes2
            else:
                ax2 = None
            export_xmgrace(path,self.axes,self.fcclass_list,self.labs,ax2=ax2,specs=spc)
    
    def on_about(self):
        msg = """
       A python application to analyze FCclasses TI spectra
       J.Cerezo, May 2016
       
       Version info 
        Git commit: %s
        Date: %s
        
       This program is free software; you can redistribute it and/or modify  
       it under the terms of the GNU General Public License as published by 
       the Free Software Foundation; either version 2 of the License, or
       (at your option) any later version.
       
       Write comments to: j.cerezo@pi.iccom.cnr.it
        """%(version_tag.COMMIT,version_tag.DATE)
        QMessageBox.about(self, "About the app", msg.strip())
        
    def on_man(self):
        msg = """ 
       SHORT GUIDE:
        *Mouse/keyboard interactions with the plot:
         -Right-mouse-click on a stick: get info from a transition
         -Keyboard "+"/"-" keys: set the next/previous indexed transiton: 
         -Left-mouse-click on a stick: place a label
         -Left-mouse-button (press and hold) on a label: move a label 
         -Right-mouse-click on the label: remove the label
        
        *File operations ("File" menu)
         -Save plot: save png figure (also available from matplotlib toolbar)
         -Export to xmgrace: export current plot (including labels) to xmgrace
         -Import plot: import reference plot
         
        *Manipulation of the reference plot:
         -In "Manipulation" menu: shift and scale with respect to simulated
         -In "Reference spectrum" table:
           [X] button: clear spectrum
           [T] buttons: reset to current values
           Scale/Shift cells: manually change the values
           
        *Search transition/progression box
         Select a given transition of progression
         The different syntax accepted are
         
          Mode1(Quanta1),Mode2(Quanta2)... 
           e.g. 1(1),2(1) : select the transiton from ground to 1(1),2(1)
           
          Mode1(Quanta1),Mode2(P)... 
           e.g. 1(1),2(P) : select the transiton from ground to 
                            1(1),2(1); 1(1),2(2); 1(1),2(3)...
           Note: 
            (P) can only be specified on one mode
           
          Mode1'(Quanta1'),Mode2'(Quanta2')... --> Mode1(Quanta1),Mode2(Quanta2)... 
           e.g. 1(1) --> 1(1) : select hot transition from 1'(1) to 1(1)
           Notes:
            (P) can be specified (on the final modes only)
            The ground state is specified as 0. E.g. 0-->0 is the 0-0 transition,
            1(1)-->0 is the M-0 transition and 0-->8(1) is equivalent to 8(1)
        """
        QMessageBox.about(self, "Instructions", msg.strip())
        
    #==========================================================
    # ANALYSIS FUNCTIONS
    #==========================================================
    def compute_moments(self):
        """
        Compute First and Second moments for the convoluted spectrum
        """
        x = self.spectrum_sim[0].get_xdata()
        y = self.spectrum_sim[0].get_ydata()
        
        # Zero
        m0 = np.trapz(y, x)
        y /= m0
        # First
        m1 = np.trapz(y*x, x)
        # Second
        m2 = np.trapz(y*x**2, x)
        # Sigma
        sgm = np.sqrt(m2-m1**1)
        
        result = """  MOMENTA ANALYSIS
  ====================
  
  * Covoluted Spectrum
  ---------------------------------
  1st Moment (eV) = %.3f
  
  2nd Moment (eV^2) = %.3f
  
  Sigma (eV) = %.3f
  ---------------------------------
        """ % (m1, m2, sgm)
        
        if self.spectrum_ref:
            x = self.spectrum_ref[0].get_xdata()
            y = self.spectrum_ref[0].get_ydata()
            
            # Zero
            m0 = np.trapz(y, x)
            y /= m0
            # First
            m1 = np.trapz(y*x, x)
            # Second
            m2 = np.trapz(y*x**2, x)
            # Sigma
            sgm = np.sqrt(m2-m1**1)
            # Add to message
            result = result+"""

  * Reference Spectrum
  ---------------------------------
  1st Moment (eV) = %.3f
  
  2nd Moment (eV^2) = %.3f
  
  Sigma (eV) = %.3f
  ---------------------------------
        """ % (m1, m2, sgm)
        
        self.analysis_box.setText(result)
        
        
    def shift_to_simulated(self):
        """
        Shift the reference expectrum to match the first moment or Emax of the simulated one
        The type of match is set through a pop-up question
        """
        match_options = ("First Moment","Peak Maximum")
        match, ok = QInputDialog.getItem(self, "Shift assistant", 
                    "Reference for shifting", match_options, 0, False)
        if not ok:
            return
        if match == "First Moment":
            x = self.spectrum_sim[0].get_xdata()
            y = self.spectrum_sim[0].get_ydata()
            # Load simulated data
            # Zero
            m0 = np.trapz(y, x)
            # First
            x0_sim = np.trapz(y*x/m0, x)
            
            # Load reference data
            if self.spectrum_ref:
                x = self.spectrum_ref[0].get_xdata()
                y = self.spectrum_ref[0].get_ydata()
            else:
                return
            # Zero
            m0 = np.trapz(y, x)
            # First
            x0_ref = np.trapz(y*x/m0, x)
        elif match == "Peak Maximum":
            x = self.spectrum_sim[0].get_xdata()
            y = self.spectrum_sim[0].get_ydata()
            x0_sim = x[y.argmax()]
            
            # Load reference data
            if self.spectrum_ref:
                x = self.spectrum_ref[0].get_xdata()
                y = self.spectrum_ref[0].get_ydata()
            else:
                return
            x0_ref = x[y.argmax()]
        
        # Update global data (this triggers the shift)
        ref_shift_new = self.ref_shift + x0_sim-x0_ref
        self.refspc_table.setItem(2,1, QTableWidgetItem(str(ref_shift_new)))
        
        
    def scale_to_simulated(self):
        """
        Scale the reference expectrum to match the maximum of the simulated one
        """
        y = self.spectrum_sim[0].get_ydata()
        y0_sim = max(abs(y.max()),abs(y.min()))
        
        # Load reference data
        if self.spectrum_ref:
            y = self.spectrum_ref[0].get_ydata()
        else:
            return
        y0_ref = max(abs(y.max()),abs(y.min()))
        
        # Update global data (this triggers the shift)
        ref_scale_new = self.ref_scale * y0_sim/y0_ref
        self.refspc_table.setItem(3,1, QTableWidgetItem(str(ref_scale_new)))
        
        
    #==========================================================
    # LOAD DATA AND ELEMENTS
    #==========================================================
    def load_sticks(self,stick_type):
        """ 
        Load stick spectra for all classes (C0,C1...,CHot)
        - The spectra objects are stored in a list: self.stickspc
        """
        # clear the axes and redraw the plot anew
        # 
        self.axes.set_title('TI stick spectrum from $\mathcal{FC}classes$',fontsize=18)
        self.axes.set_xlabel('Energy (eV)',fontsize=16)
        self.axes.set_ylabel('Stick Intensity',fontsize=16)
        self.axes.tick_params(direction='out',top=False, right=False)
        
        
        #Plotting sticks and store objects
        # Set labels and colors
        label_list = ['0-0']+[ 'C'+str(i) for i in range(1,8) ]+['Hot']
        color_list = ['k', 'b', 'r', 'g', 'c', 'm', 'brown', 'pink', 'orange' ]


        #Inialize variables
        self.stickspc = []
        xmin =  999.
        xmax = -999
        for iclass in range(9):
            x = np.array([ self.fcclass_list[iclass][i].DE        for i in range(len(self.fcclass_list[iclass])) ])
            if stick_type == "fc":
                # Get intensity as FC^2
                for i in range(len(self.fcclass_list[iclass])):
                    self.fcclass_list[iclass][i].intensity = self.fcclass_list[iclass][i].fcfactor**2
            y = np.array([ self.fcclass_list[iclass][i].intensity for i in range(len(self.fcclass_list[iclass])) ])
            z = np.zeros(len(x))
            if len(x) == 0:
                self.stickspc.append(None)
            else:
                self.stickspc.append(self.axes.vlines(x,z,y,linewidths=1,color=color_list[iclass],
                                                      label=label_list[iclass],picker=5))
                xmin = min([xmin,min(x)])
                xmax = max([xmax,max(x)])
                # Getting the type from the objects when loaded for future use
                self.LineCollectionType = type(self.stickspc[-1])
                
        self.axes2.set_xlim([xmin-0.15,xmax+0.15])
        
        self.canvas.draw()
        
        
    def load_legend(self):
        """
        Plot legend, which is pickable so as to turn plots on/off
        """
        
        #Legends management
        # First get lines from both axes
        lns  = list(filter(lambda x:x,self.stickspc))+self.spectrum_sim
        labs = [l.get_label() for l in lns]
        # Set the location according to the type of calculation
        if self.spc_type in ["abs","ecd"]:
            position="upper right"
        else:
            position="upper left"
        self.legend = self.axes.legend(lns,labs,loc=position, fancybox=True, shadow=True)
        self.LegendType = type(self.legend)
        self.legend.get_frame().set_alpha(0.4)
        self.legend.set_picker(5)
        # we will set up a dict mapping legend line to orig line, and enable
        # picking on the legend line (from legend_picking.py)
        self.legend_lines = dict()
        # Note: mechanism to get rid out of None elements in the list
        #  filter(lambda x:x,lista)) evaluates every member in the list, and only take if
        #  the result of the evaluation is True. None gives a False
        for legline, origline in zip(self.legend.get_lines(), lns):
            legline.set_picker(5)  # 5 pts tolerance
            self.legend_lines[legline] = origline
            # Get type from the object for future use
            self.Line2DType = type(legline)
            
            
    def set_axis_labels(self):
        if self.data_type == "Lineshape":
            self.axes2.set_ylabel(r'Lineshape (a.u.)',fontsize=16)
        elif self.spc_type == 'abs':
            self.axes2.set_ylabel(r'$\varepsilon$ (dm$^3$mol$^{-1}$cm$^{-1}$)',fontsize=16)
        elif self.spc_type == 'ecd':
            self.axes2.set_ylabel(r'$\Delta\varepsilon$ (dm$^3$mol$^{-1}$cm$^{-1}$)',fontsize=16)
        elif self.spc_type == 'emi':
            self.axes2.set_ylabel(r'I (molecule$^{-1}$ns$^{-1}$]',fontsize=16)
        elif self.spc_type == 'cpl':
            self.axes2.set_ylabel(r'$\Delta$I (molec$^{-1}$ns$^{-1}$]',fontsize=16)
        self.axes2.set_xlabel('Energy (eV)',fontsize=16)
        #self.axes.tick_params(direction='out',top=False, right=False)
            

    def load_fixconvolution(self,x,y):
        self.spectrum_sim = self.axes2.plot(x,y,'--',color='k',label="Conv")
        
        self.canvas.draw()
        
            
    def load_convoluted(self):
        str = unicode(self.broadbox.text())
        hwhm = float(str)
        fixaxes = self.fixaxes_cb.isChecked()
        
        #Convolution
        xc,yc = convolute([self.xbin,self.ybin],hwhm=hwhm,broad=self.broadening,input_bins=self.inputBins_cb.isChecked())
        # Plot convoluted
        self.spectrum_sim = self.axes2.plot(xc,yc,'--',color='k',label="Conv")
        if not fixaxes:
            self.rescale_yaxis()
        
        self.canvas.draw()
        
        
    def update_convolute(self):
        str = unicode(self.broadbox.text())
        hwhm = float(str)
        fixaxes = self.fixaxes_cb.isChecked()
        
        #Convolution (in energy(eV))
        xc,yc = convolute([self.xbin,self.ybin],hwhm=hwhm,broad=self.broadening,input_bins=self.inputBins_cb.isChecked())
        # Re-Plot convoluted
        #self.spectrum_sim[0].remove()
        self.spectrum_sim[0].set_xdata(xc)
        self.spectrum_sim[0].set_ydata(yc)
        if not fixaxes:
            self.rescale_yaxis()

        self.canvas.draw()
        
        
    def load_experiment_spc(self,x,y):
        fixaxes = self.fixaxes_cb.isChecked()
        
        x,y = [np.array(x), np.array(y)]
        # Plot experiment (only one experiment is allowed)
        if self.spectrum_ref:
            self.spectrum_ref[0].remove()
        self.spectrum_ref = self.axes2.plot(x,y,'-',color='gray',label="Ref")
        if not fixaxes:
            self.rescale_yaxis()
        
        self.canvas.draw()
        
    def rescale_yaxis(self):
        """"
        Set the range so as to keep the same zero as in the case of the sticks
        getting the maximum between exp and sim
        """
        ysim = self.spectrum_sim[0].get_ydata()
        if self.spectrum_ref:
            yexp = self.spectrum_ref[0].get_ydata()
            ymax2 = max(ysim.max(),yexp.max())
            ymin2 = min(ysim.min(),yexp.min())
        else:
            ymax2 = ysim.max()
            ymin2 = ysim.min()
            
        ymin,ymax = self.axes.get_ylim()
        if abs(ymin2) > abs(ymax2):
            ymin2 *= 1.05
            ymax2 = ymin2/ymin * ymax
        else:
            ymax2 *= 1.05
            ymin2 = ymax2/ymax * ymin
            
        self.axes2.set_ylim([ymin2,ymax2])
        
        return
    
    def shift_spectrum(self,x,y,shift):
        # If Intensity, we need to pass to LS before shifting
        if self.data_type == "Intensity":
            n = SpcConstants.exp[self.spc_type]
            factor = SpcConstants.factor[self.spc_type]
            # Division x/27.2116 could be included in the factor
            y /= (x/27.2116)**n * factor
        x = x + shift
        # If Intensity, set back from Lineshape
        if self.data_type == "Intensity":
            # Division x/27.2116 could be included in the factor
            y *= (x/27.2116)**n * factor
        
        return x,y
        
        
    #==========================================================
    # FUNCTIONS TO INTERACT WITH THE PLOT
    #==========================================================
    def on_pick(self, event):
        """"
        This is a common manager for picked objects
        
        The event received here is of the type
        matplotlib.backend_bases.PickEvent, they
        differeciate the artist
        
        "Pickable" opjects are labels, sticks and legend lines
        but they have different attributes and, thus,
        need to be handled differently
        NOTE: The matplotlib.collections.LineCollection requires 
        to load matplotlib apart from matplotlib.pyplot (which is 
        al ready used). Not a problem as it does not affect the weight
        of the final app
        * with matplotlib: 223MB
        * only matplotlib.pyplot: 223MB
        Anyway, it might be useful to get the types from the objects
        so as to avoid the need of knowing them
        So, we get the type from the object when it is called
        """
        if type(event.artist) == self.LineCollectionType:
            self.select_stick(event)
        elif type(event.artist) == self.AnnotationType:
            self.active_label = event.artist
        elif type(event.artist) == self.Line2DType:
            self.del_stick_marker()
            self.interact_with_legend(event)
        #elif type(event.artist) == self.LegendType:
            #self.legend.draggable()
            
            
    def select_stick(self,event):
        """
        Select a given stick transiton by clicking with the mouse, and do:
         1) Put a label if the left-button is clicked
         2) Highlight is right button is clicked
        """
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        iclass = self.stickspc.index(event.artist)

        #Now set data and label positions
        stick_x = [ self.fcclass_list[iclass][i].DE        for i in event.ind ]
        stick_y = [ self.fcclass_list[iclass][i].intensity for i in event.ind ]
        distances = np.hypot(x-stick_x, y-stick_y)
        indmin = distances.argmin()
        dataind = event.ind[indmin]
        
        # Highlight or add label
        tr = self.fcclass_list[iclass][dataind]
        if event.mouseevent.button == 3:
            self.active_tr = tr
            # Clear seach_box as it would be outdated
            self.search_box.setText('')
            self.set_stick_marker()
        elif event.mouseevent.button == 1:
            # The transition info need to be gathered
            self.fcclass_list[iclass][dataind].info()
            self.add_label(tr)


    def on_press_key(self, event):
        """"
        Manage keyword interactions with the graph
        """
        if self.active_tr is None: 
            return
        if event.key not in ('+', '-'): 
            return
        if event.key=='+': 
            inc = 1
        else:  
            inc = -1
            
        # Clear seach_box as it would be outdated
        self.search_box.setText('')

        # Get the active class
        if self.active_tr.motherstate == 1:
            iclass = self.active_tr.fcclass
        else:
            iclass = 8 #Hot bands class
        index = self.fcclass_list[iclass].index(self.active_tr)
        index += inc
        
        # Manage transition between classes
        if index > len(self.fcclass_list[iclass])-1:
            iclass += 1
            index = 0
        elif index < 0:
            iclass -= 1
            try:
                index = len(self.fcclass_list[iclass])-1
            except:
                pass
        
        # Get next transition. If out-of-bounds, simply 
        # clean the panel
        try:
            self.active_tr = self.fcclass_list[iclass][index]
            self.set_stick_marker()
        except IndexError:
            self.del_stick_marker()
        
        
    def on_press_mouse(self,event):
        self.past_event = None
        if self.active_label is None:
            return
        self.past_event = event
        
        
    def move_label(self,event):
        if self.active_label is None: return
        if event.button != 1: return
        if self.past_event is  None or event.inaxes!=self.past_event.inaxes: return
    
        # The label is picked even if the connecting line
        # is pressed. Since, this is not a nice behabiour 
        # we reject the press events done on the line
        x0, y0 = self.active_label.get_position()
        xe, ye = self.past_event.xdata,self.past_event.ydata
        x_range=self.axes.get_xlim()
        # Set the presision based on the plot scale
        eps_x = (x_range[1]-x_range[0])/50.
        y_range=self.axes.get_ylim()
        eps_y = (y_range[1]-y_range[0])/50.
        if abs(x0-xe) > eps_x and abs(x0-xe) > eps_y:
            return
        # Get the plot distance
        
        dx = event.xdata - self.past_event.xdata
        dy = event.ydata - self.past_event.ydata
        self.active_label.set_position((x0+dx,y0+dy))
        self.canvas.draw()
        # Update the pressevent
        self.past_event = event
        
        
    def release_label(self,event):
        if self.active_label is None: return
        if event.button != 1: return
        if event.inaxes!=self.axes:
         return
        self.past_event = None
        self.active_label = None
        
        
    def delete_label(self,event):
        if self.active_label is None: return
        if event.button != 3: return
    
        # As in the move_label, only remove if we click
        # on the text (not on the line)
        x0, y0 = self.active_label.get_position()
        xe, ye = self.past_event.xdata,self.past_event.ydata
        x_range=self.axes.get_xlim()
        # Set the presision based on the plot scale
        eps_x = (x_range[1]-x_range[0])/70.
        y_range=self.axes.get_ylim()
        eps_y = (y_range[1]-y_range[0])/70.
        if abs(x0-xe) > eps_x and abs(x0-xe) > eps_y:
            return
    
        #Remove the label
        self.active_label.remove()
        #And substract the corresponding entry from the dict
        self.labs.pop(self.active_label)
        #We deactivate the lab. Otherwise, if we reclikc on the same point, it raises an error
        self.active_label.set_visible(False)
        self.canvas.draw()


    def interact_with_legend(self,event):
        """
        DESCRIPTION
        ------------
        Function to activate/deactivate a plot, clicking on the legend 
        
        NOTES
        -----
        Based on matplotlib example: legend_picking.py
        http://matplotlib.org/examples/event_handling/legend_picking.html
        """
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        class_spc = self.legend_lines[legline]
        vis = not class_spc.get_visible()
        class_spc.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
            
        self.canvas.draw()
            
    # FUNCTIONS WITHOUT EVENT
    def add_label(self,tr):
        
        stick_x = tr.DE
        stick_y = tr.intensity
        
        xd = stick_x
        yd = stick_y
        xl = stick_x - 0.00
        #Since the intensity may change several orders of magnitude, 
        # we better don't shift it
        yl = stick_y + 0.0
        
        #Define label as final modes description
        if tr.fcclass == 0 and tr.motherstate == 1:
            label='0-0'
            agrlabel='0-0'
        else:
            label='$'
            agrlabel=''
            if tr.motherstate > 1:
                # For Hot bands, write also initial state
                for i in range(0,len(tr.init)):
                    label = label+str(tr.init[i])+'^{'+str(tr.qinit[i])+'},'
                    agrlabel = agrlabel+str(tr.init[i])+'\S'+str(tr.qinit[i])+'\N,'
                if len(tr.init) == 0:
                    label=label+"0,"
                    agrlabel=agrlabel+"0,"
                label = label[0:-1]+'-'
                agrlabel = agrlabel[0:-1]+'-'
            for i in range(0,len(tr.final)):
                label = label+str(tr.final[i])+'^{'+str(tr.qfinal[i])+'},'
                agrlabel = agrlabel+str(tr.final[i])+'\S'+str(tr.qfinal[i])+'\N,'
            if len(tr.final) == 0:
                label=label+"0,"
                agrlabel=agrlabel+"0,"
            #Remove trailing comma
            label = label[0:-1]+'$'
            agrlabel = agrlabel[0:-1]

        #In labelref we get the annotation class corresponding to the 
        #current label. labelref is Class(annotation)
        labelref = self.axes.annotate(label, xy=(xd, yd), xytext=(xl, yl),picker=1,
                            arrowprops=dict(arrowstyle="-",
                                            color='grey'))
        self.AnnotationType = type(labelref)

        #Check whether the label was already assigned or not
        set_lab = True
        for labref in self.labs:
            lab = self.labs[labref]
            if lab == agrlabel:
                print "This label was already defined"
                set_lab = False
        if set_lab:
            # The dictionary labs relates each labelref(annotation) to 
            # the agrlabel. Note that the mpl label can be retrieved
            # from labelref.get_name()
            self.labs[labelref] = agrlabel
            self.canvas.draw()
        else:
            labelref.remove()
    
    
    # RESPONSES TO SIGNALS
    def update_fixlegend(self):
        fixlegend = not self.fixlegend_cb.isChecked()
        self.legend.draggable(fixlegend)
        
    
    def update_hwhm_from_slider(self,UpdateConvolute=True):
        hwhmmin = 0.01
        hwhmmax = 0.1
        slidermin = 1   # this is not changed
        slidermax = 100 # this is not changed
        hwhm = float((hwhmmax-hwhmmin)/(slidermax-slidermin) * (self.slider.value()-slidermin) + hwhmmin)
        hwhm = round(hwhm,3)
        self.broadbox.setText(str(hwhm))
        if (UpdateConvolute):
            self.update_convolute()
        
        
    def update_hwhm_from_textbox(self):
        hwhmmin = 0.01
        hwhmmax = 0.1
        slidermin = 1   # this is not changed
        slidermax = 100 # this is not changed
        str = unicode(self.broadbox.text())
        hwhm = float(str)
        sliderval = int((slidermax-slidermin)/(hwhmmax-hwhmmin) * (hwhm-hwhmmin) + slidermin)
        sliderval = min(sliderval,slidermax)
        sliderval = max(sliderval,slidermin)
        self.slider.setValue(sliderval)
        self.update_convolute()
        
        
    def set_stick_marker(self):
        if self.active_tr is None: return
        
        tr = self.active_tr 
        stick_x = tr.DE
        stick_y = tr.intensity
        
        # Add transition info to analysis_box
        self.analysis_box.setText(self.active_tr.info())
        
        if self.selected:
            self.selected.remove()
        self.selected  = self.axes.vlines([stick_x], [0.0], [stick_y], linewidths=3,
                                  color='yellow', visible=True, alpha=0.7)
        self.canvas.draw()
        
    
    def del_stick_marker(self):
        
        self.analysis_box.setText("")
        self.search_box.setText('')
        self.active_tr = None 
        if self.selected:
            self.selected.remove()
            self.selected = None
            self.canvas.draw()
        
        
    def reset_labels(self):
        #We need a local copy of labs to iterate while popping
        labs_local = [ lab for lab in self.labs ]
        for lab in labs_local:
            lab.remove()
            self.labs.pop(lab)
        self.canvas.draw()
        
        
    def update_broad_function(self):
        self.broadening = self.select_broad.currentText()
        self.update_convolute()
        
        
    def update_data_type(self):
        current_data_type = self.data_type
        self.data_type = self.select_data_type.currentText()
        
        if current_data_type != self.data_type:
            factor = SpcConstants.factor[self.spc_type]
            n = SpcConstants.exp[self.spc_type]
            if self.with_fort22:
                if self.data_type == "Lineshape":
                    # Division x/27.2116 could be included in the factor
                    self.ybin /= (self.xbin/27.2116)**n * factor
                elif self.data_type == "Intensity":
                    # Division x/27.2116 could be included in the factor
                    self.ybin *= (self.xbin/27.2116)**n * factor
            else:
                if self.spectrum_sim:
                    x = self.spectrum_sim[0].get_xdata()
                    y = self.spectrum_sim[0].get_ydata()
                    if self.data_type == "Lineshape":
                        # Division x/27.2116 could be included in the factor
                        y /= (x/27.2116)**n * factor
                    elif self.data_type == "Intensity":
                        # Division x/27.2116 could be included in the factor
                        y *= (x/27.2116)**n * factor
                    self.spectrum_sim[0].set_ydata(y)
            
            if self.spectrum_ref:
                x = self.spectrum_ref[0].get_xdata()
                y = self.spectrum_ref[0].get_ydata()
                if self.data_type == "Lineshape":
                    # Division x/27.2116 could be included in the factor
                    y /= (x/27.2116)**n * factor
                elif self.data_type == "Intensity":
                    # Division x/27.2116 could be included in the factor
                    y *= (x/27.2116)**n * factor
                self.spectrum_ref[0].set_ydata(y)
                
            self.set_axis_labels()
            if self.with_fort22:
                self.update_convolute()
            else:
                fixaxes = self.fixaxes_cb.isChecked()
                if not fixaxes:
                    self.rescale_yaxis()
                self.canvas.draw()
                
                
    def table_buttons_action(self,i,j):
        # "X" button: clear spectrum
        # "T" reset shift/scale
        if not self.spectrum_ref:
            return
        if (i,j) == (1,2):
            clear_msg = "Clear reference spectrum?"
            reply = QtGui.QMessageBox.question(self, 'Clear Spectrum', 
                         clear_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return            
            self.spectrum_ref[0].remove()
            self.spectrum_ref = None
            self.canvas.draw()
            celllabel = ["No","-","-"]
            for i,j in [(1,1),(2,1),(3,1)]:
                self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
                cell = self.refspc_table.item(i,j)
                cell.setTextColor(Qt.black)
                cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
            # Disable manipulations
            self.shiftref_action.setEnabled(False)
            self.scaleref_action.setEnabled(False)

        elif (i,j) == (2,2):
            # Tare the shift
            msg = "Reset shift to current value?"
            reply = QtGui.QMessageBox.question(self, 'Reset Shift', 
                         msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return      
            self.ref_shift = 0.0
            self.refspc_table.setItem(2,1, QTableWidgetItem(str(self.ref_shift)))
            
        elif (i,j) == (3,2):
            # Tare the scale
            msg = "Reset scale to current value?"
            reply = QtGui.QMessageBox.question(self, 'Reset Scale', 
                         msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return      
            self.ref_scale = 1.0
            self.refspc_table.setItem(3,1, QTableWidgetItem(str(self.ref_scale)))
            
            
    def change_refspc(self,i,j):
        fixaxes = self.fixaxes_cb.isChecked()
        if not self.spectrum_ref:
            return
        cell = self.refspc_table.item(i,j)
        if (i,j) == (2,1):
            # Shift
            #========
            new_shift = float(cell.text())
            # Get the relative shift from the current global shift
            shift = new_shift - self.ref_shift
            self.ref_shift = new_shift
            x = self.spectrum_ref[0].get_xdata()
            y = self.spectrum_ref[0].get_ydata()
            x,y = self.shift_spectrum(x,y,shift)
            self.spectrum_ref[0].set_xdata(x)
            self.spectrum_ref[0].set_ydata(y)
        elif (i,j) == (3,1):
            # Scale
            #========
            new_scale = float(cell.text())
            # Get the relative scale from the current global scaling
            scale = new_scale/self.ref_scale
            self.ref_scale = new_scale
            y = self.spectrum_ref[0].get_ydata() * scale
            self.spectrum_ref[0].set_ydata(y)
            
        if not fixaxes:
            self.rescale_yaxis()
        self.canvas.draw()
        
        
    def search_transitions(self):
        """
        Read the content in the seach box and interpret
        it to a search action to pick the selected transitions
        """
        command = str(self.search_box.text())
        
        # First, check if the used is only asking for help
        if command == 'help':
            msg = """ 
Examples
---------
1(1),2(1)   - From ground to 1(1),2(1)
1(1),2(P)   - Progression over mode 2
1(1)-->1(1) - Hot transition
1(1)-->0    - M-0 transition
0-->0       - 0-0 transition
        """
            self.analysis_box.setText(msg)
            return
        
        # Start processing the command
        self.statusBar().showMessage('Searching transition '+command, 2000)
        # Here, use regexp to set the validity of the command
        pattern = r'(( )*(([0-9pP]+\([0-9]+\)( )*\,( )*)*( )*[0-9pP]+\([0-9]+\)|0)( )*\-\->( )*){0,1}(( )*([0-9]+\([0-9pP]+\)( )*\,( )*)*( )*[0-9]+\([0-9pP]+\)|0)( )*'
        match = re.match(pattern,command)
        if not match or match.group() != command:
        #     [m'1(q'1),m'2(q'2)... -->] m1(q1),m2(q2)...
        # or: 0 --> 0 
            if self.selected:
                self.selected.remove()
                self.selected = None
                self.canvas.draw()
            self.analysis_box.setText("")
            self.statusBar().showMessage('Invalid syntax', 2000)
            return
        transition_list=command.split('-->')
        if len(transition_list) == 2:
            state_ini,state_fin = transition_list
        elif len(transition_list) == 1:
            state_ini  = "0"
            state_fin, = transition_list
            
            
        # Initial state
        imodes  = []
        iquanta = []
        if not "(" in state_ini and int(state_ini) == 0:
            transitions = []
        else:
            transitions = state_ini.split(',')
        for tr in transitions:
            m,q = tr.split('(')
            # If re.match is done well (as it is?) these try/except
            # should no be needed. To be removed
            try:
                imodes.append(int(m))
            except:
                self.statusBar().showMessage('Invalid syntax', 2000)
                return
            try:
                iquanta.append(int(q.replace(')','')))
            except:
                self.statusBar().showMessage('Invalid syntax', 2000)
                return
        
        
        # Final state
        fmodes  = []
        fquanta = []
        pindex = []
        if not "(" in state_fin and int(state_fin) == 0:
            transitions = []
        else:
            transitions = state_fin.split(',')
        for tr in transitions:
            m,q = tr.split('(')
            try:
                fmodes.append(int(m))
            except:
                self.statusBar().showMessage('Invalid syntax', 2000)
                return
            if q.strip().upper() == 'P)':
                progression=True
                pindex.append(len(fquanta))
                fquanta.append(-1)
            else:
                try:
                    fquanta.append(int(q.replace(')','')))
                except:
                    self.statusBar().showMessage('Invalid syntax', 2000)
                    return

        if (len(pindex) > 1):
            self.statusBar().showMessage('Invalid syntax', 2000)
            return
        elif (len(pindex) == 1):
            pindex = pindex[0]
        elif (len(pindex) == 0):
            progression = False

        # Set FCclass
        fclass = len(fmodes)
        iclass = len(imodes)
        if iclass != 0:
            fcclass = self.fcclass_list[8]
        else:
            fcclass = self.fcclass_list[fclass]
        
        tr_select = []
        i  = 0
        take_next_tr = True
        # To handle the cases where the even are active but not the odd ones we use the imissing counter
        imissing = 0
        # Sorting two list (apply change of one into the other)
        # http://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
        while take_next_tr:
            take_next_tr = False
            if progression:
                i += 1
                fquanta[pindex] = i
            for tr in fcclass:
                if iclass != 0:
                    mi_srch,qi_srch = zip(*sorted(zip(imodes,iquanta)))
                    mi_tr,qi_tr     = zip(*sorted(zip(tr.init[:iclass],tr.qinit[:iclass])))
                else:
                    mi_srch,qi_srch = 0,0
                    mi_tr,qi_tr     = 0,0
                if fclass != 0:
                    # HotClass include X->0, which has len(tr.final)=0
                    # and crash with the zip. This is becase in tr.info()
                    # the lists tr.init and tr.final are cut to delete the 
                    # zeroes
                    if len(tr.final) == 0:
                        continue
                    mf_srch,qf_srch = zip(*sorted(zip(fmodes,fquanta)))
                    mf_tr,qf_tr     = zip(*sorted(zip(tr.final[:fclass],tr.qfinal[:fclass])))
                else:
                    mf_srch,qf_srch = 0,0
                    mf_tr,qf_tr     = 0,0
                if mf_srch==mf_tr and qf_srch==qf_tr and mi_srch==mi_tr and qi_srch==qi_tr:
                    tr_select.append(tr)
                    imissing -= 2
                    break
            imissing += 1
            if progression and imissing<2:
                take_next_tr = True
        
        if progression and tr_select:
            # Get the origin of the transition
            fclass -= 1
            if iclass == 0:
                fcclass = self.fcclass_list[fclass]
            if fclass == 0:
                tr = fcclass[0]
                tr_select.insert(0,tr)
            else:
                fmodes.pop(pindex)
                fquanta.pop(pindex)
                for tr in fcclass:
                    if iclass != 0:
                        mi_srch,qi_srch = zip(*sorted(zip(imodes,iquanta)))
                        mi_tr,qi_tr     = zip(*sorted(zip(tr.init[:iclass],tr.qinit[:iclass])))
                    else:
                        mi_srch,qi_srch = 0,0
                        mi_tr,qi_tr     = 0,0
                    if fclass != 0:
                        mf_srch,qf_srch = zip(*sorted(zip(fmodes,fquanta)))
                        mf_tr,qf_tr     = zip(*sorted(zip(tr.final[:fclass],tr.qfinal[:fclass])))
                    else:
                        mf_srch,qf_srch = 0,0
                        mf_tr,qf_tr     = 0,0
                    if mf_srch==mf_tr and qf_srch==qf_tr and mi_srch==mi_tr and qi_srch==qi_tr:
                        tr_select.insert(0,tr)
                        break
        
        # If only one tr selected, set it to active_tr
        # This allows to browse with +/- from this tr
        if (len(tr_select)) == 1:
            self.active_tr
            
        # Remove stick if there was already 
        if self.selected:
            self.selected.remove()
            self.selected = None
            self.canvas.draw()
        self.analysis_box.setText("")
                
        if tr_select:
            self.statusBar().showMessage('Found', 2000)
            msg=""
            stick_x = []
            stick_y = []
            for tr in tr_select:
                stick_x.append(tr.DE)
                stick_y.append(tr.intensity)
                msg = msg+"\n"+tr.info()
            zero = np.zeros(len(stick_x))
            # Add transition info to analysis_box
            self.analysis_box.setText(msg.strip())
            self.selected  = self.axes.vlines(stick_x, zero, stick_y, linewidths=3,
                                      color='yellow', visible=True, alpha=0.7)
            self.canvas.draw()
        else:
            self.statusBar().showMessage('Not found', 2000)

        
        
    #==========================================================
    # MAIN FRAME, CONNECTIONS AND MENU ACTIONS
    #==========================================================    
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        #self.dpi = 100
        #self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        # (lo de arriba, estaba en el original)
        # 
        # Usamos la figura/ejes creados com subplots
        self.fig, self.axes2 = plt.subplots()
        # Second axis for the convoluted graph
        self.axes = self.axes2.twinx()
        # The canvas is the graphical object managed by PyQt
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        #self.canvas.setMaximumWidth(1000)
        # The following allows the interaction with the keyboard
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        #self.axes = self.fig.add_subplot(111)
        
        # Bind events to interact with the graph
        # Picking maganement
        self.canvas.mpl_connect('pick_event', self.on_pick)
        #Manage labels
        self.canvas.mpl_connect('button_press_event', self.on_press_mouse)
        self.canvas.mpl_connect('button_release_event', self.release_label)
        self.canvas.mpl_connect('motion_notify_event', self.move_label)
        self.canvas.mpl_connect('button_press_event', self.delete_label)
        # Manage highlighting
        self.canvas.mpl_connect('key_press_event', self.on_press_key)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI controls
        # 
        self.broadbox = QLineEdit()
        self.broadbox.setMinimumWidth(50)
        self.broadbox.setMaximumWidth(50)
        self.connect(self.broadbox, SIGNAL('editingFinished ()'), self.update_hwhm_from_textbox)
        
        clean_button1 = QPushButton("&Clean(Panel)")
        self.connect(clean_button1, SIGNAL('clicked()'), self.del_stick_marker)
        clean_button2 = QPushButton("&Clean(Labels)")
        self.connect(clean_button2, SIGNAL('clicked()'), self.reset_labels)
        
        self.select_broad = QComboBox()
        self.select_broad.addItems(["Gau","Lor"])
        self.select_broad.currentIndexChanged.connect(self.update_broad_function)
        
        self.select_data_type = QComboBox()
        self.select_data_type.addItems(["Intensity","Lineshape"])
        self.select_data_type.currentIndexChanged.connect(self.update_data_type)
        
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMaximumWidth(200)
        self.slider.setRange(1, 100)
        self.slider.setValue(30)
        self.slider.setTracking(True)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.connect(self.slider, SIGNAL('valueChanged(int)'), self.update_hwhm_from_slider)
        
        self.fixaxes_cb = QCheckBox("Fix y-axis")
        self.fixaxes_cb.setChecked(False)
        self.fixaxes_cb.setMaximumWidth(100)
        
        self.fixlegend_cb = QCheckBox("Fix legend")
        self.fixlegend_cb.setChecked(True)
        self.fixlegend_cb.setMaximumWidth(100)
        self.connect(self.fixlegend_cb, SIGNAL('stateChanged(int)'), self.update_fixlegend)
        
        self.inputBins_cb = QCheckBox("Input bins")
        self.inputBins_cb.setChecked(True)
        self.inputBins_cb.setMaximumWidth(100)
        # Update when clicking
        self.connect(self.inputBins_cb, SIGNAL('stateChanged(int)'), self.update_convolute)
        
        # Labels
        hwhm_label  = QLabel('   HWHM')
        eV_label    = QLabel('(eV)')
        broad_label = QLabel('Broadening')
        datatype_label = QLabel('Data Type')
        search_label = QLabel('Select transition/progression')
        # Splitters
        vline = QFrame()
        vline.setFrameStyle(QFrame.VLine)
        vline.setLineWidth(1)
        vline2 = QFrame()
        vline2.setFrameStyle(QFrame.VLine)
        vline2.setLineWidth(1)
        
        # Analysis box
        self.analysis_box = QTextEdit(self.main_frame)
        self.analysis_box.setFontFamily('Arial')
        self.analysis_box.setReadOnly(True)
        self.analysis_box.setMinimumWidth(200)
        self.analysis_box.setMaximumWidth(250)

        # Search box
        self.search_box = QLineEdit(self.main_frame)
        self.search_box.setMinimumWidth(200)
        self.search_box.setMaximumWidth(250)
        self.search_box.returnPressed.connect(self.search_transitions)
        
        # Table for the reference spectrum
        # Ids ordered by column
        cellids = [[(1,0),(2,0),(3,0)] ,[(1,1),(2,1),(3,1)],[(1,2),(2,2),(3,2)]]
        self.refspc_table = QTableWidget(self.main_frame)
        self.refspc_table.setRowCount(4)
        self.refspc_table.setColumnCount(3)
        self.refspc_table.setMinimumWidth(238)
        self.refspc_table.setMaximumWidth(238)
        self.refspc_table.setMinimumHeight(126)
        self.refspc_table.setMaximumHeight(126)
        # Set data
        ## Title row
        self.refspc_table.setSpan(0,0,1,3)
        self.refspc_table.setItem(0,0, QTableWidgetItem("Reference Spectrum"))
        font = QFont()
        font.setBold(False)
        #font.setWeight(75)
        font.setPointSize(12)
        title = self.refspc_table.item(0,0)
        title.setBackgroundColor(Qt.lightGray)
        title.setTextColor(Qt.black)
        title.setFont(font)
        title.setTextAlignment(Qt.AlignCenter)
        title.setFlags(title.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
        ## Fist column
        celllabel = ["Loaded?","Shift(eV)","Y-scale"]
        font = QFont()
        font.setBold(True)
        font.setWeight(75)
        for i,j in cellids[0]:
            self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
            cell = self.refspc_table.item(i,j)
            cell.setBackgroundColor(Qt.gray)
            cell.setTextColor(Qt.white)
            cell.setFont(font)
            # Set non editable. See: http://stackoverflow.com/questions/2574115/how-to-make-a-column-in-qtablewidget-read-only
            #cell.setFlags(cell.flags() ^ Qt.ItemIsEditable)
            cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
        ## Second column
        celllabel = ["No","-","-"]
        for i,j in cellids[1]:
            self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
            cell = self.refspc_table.item(i,j)
            cell.setTextColor(Qt.black)
            # Set non editable. See: http://stackoverflow.com/questions/2574115/how-to-make-a-column-in-qtablewidget-read-only
            cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
        ## Last column
        celllabel = ["X","T","T"]
        for i,j in cellids[2]:
            self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
            cell = self.refspc_table.item(i,j)
            # Set non editable. See: http://stackoverflow.com/questions/2574115/how-to-make-a-column-in-qtablewidget-read-only
            cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEditable ^ QtCore.Qt.ItemIsSelectable)
        ## Tune the X button
        xbutton = self.refspc_table.item(1,2)
        xbutton.setBackgroundColor(Qt.red)
        ## Tune the Tare buttons
        tbutton = self.refspc_table.item(2,2)
        tbutton.setBackgroundColor(Qt.blue)
        tbutton.setTextColor(Qt.white)
        tbutton = self.refspc_table.item(3,2)
        tbutton.setBackgroundColor(Qt.blue)
        tbutton.setTextColor(Qt.white)
        # Connecting cellPressed(int,int), passing its arguments to the called function
        # My function definition uses the irow and icol pressed:
        # self.table_buttons_action(self,irow,icol)
        # irow and icol are hold (someway) by the event and the fed to cellPressed:
        # cellPressed(irow,icol)
        # The connect call below passes the arguments of cellPressed to my function:
        self.refspc_table.cellPressed.connect(self.table_buttons_action)
        self.refspc_table.cellChanged.connect(self.change_refspc)
        # If we try with code below, I see no way to make it pass the args
        #self.connect(self.refspc_table, SIGNAL('cellPressed(int,int)'), lambda: self.table_buttons_action())
        # Table format
        self.refspc_table.horizontalHeader().hide()
        self.refspc_table.verticalHeader().hide()
        self.refspc_table.setColumnWidth(0,80)
        self.refspc_table.setColumnWidth(1,137)
        self.refspc_table.setColumnWidth(2,15)
        self.refspc_table.setShowGrid(False)
        
        #
        # Layout with box sizers
        # 
        # Broad selector
        vbox_select = QVBoxLayout()
        vbox_select.addWidget(broad_label)
        vbox_select.addWidget(self.select_broad)
        
        # HWHM slider with textbox
        hbox_slider = QHBoxLayout()
        hbox_slider.addWidget(self.slider)
        hbox_slider.addWidget(self.broadbox)
        hbox_slider.addWidget(eV_label)
        # HWHM label with inputBins checkbox
        hbox_hwhmlab = QHBoxLayout()
        hbox_hwhmlab.addWidget(hwhm_label)
        hbox_hwhmlab.setAlignment(hwhm_label, Qt.AlignLeft)
        hbox_hwhmlab.addWidget(self.inputBins_cb)
        hbox_hwhmlab.setAlignment(self.inputBins_cb, Qt.AlignLeft)
        # Complete HWHM widget merging all box here
        vbox_slider = QVBoxLayout()
        vbox_slider.addLayout(hbox_hwhmlab)
        vbox_slider.setAlignment(hbox_hwhmlab, Qt.AlignLeft)
        vbox_slider.addLayout(hbox_slider)
        vbox_slider.setAlignment(hbox_slider, Qt.AlignLeft)
        # Clean button
        vbox_cleaner = QVBoxLayout()
        vbox_cleaner.addWidget(clean_button1)
        vbox_cleaner.addWidget(clean_button2)
        # DataType sector
        vbox_datatype = QVBoxLayout()
        vbox_datatype.addWidget(datatype_label)
        vbox_datatype.setAlignment(datatype_label, Qt.AlignLeft)
        vbox_datatype.addWidget(self.select_data_type)
        vbox_datatype.setAlignment(self.select_data_type, Qt.AlignLeft)
        ## MAIN LOWER BOX
        hbox = QHBoxLayout()
        hbox.addLayout(vbox_cleaner)
        hbox.addWidget(vline)
        hbox.addLayout(vbox_select)
        hbox.setAlignment(vbox_select, Qt.AlignTop)
        hbox.addLayout(vbox_slider)
        hbox.setAlignment(vbox_slider, Qt.AlignTop)
        hbox.setAlignment(vbox_slider, Qt.AlignLeft)
        hbox.addWidget(vline2)
        hbox.setAlignment(vline2, Qt.AlignLeft)
        hbox.addLayout(vbox_datatype)
        hbox.setAlignment(vbox_datatype, Qt.AlignTop)
        
        # Hbox below plot
        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.fixaxes_cb)
        hbox2.setAlignment(self.fixaxes_cb, Qt.AlignLeft)
        hbox2.addWidget(self.fixlegend_cb)
        hbox2.setAlignment(self.fixlegend_cb, Qt.AlignLeft)
        hbox2.addWidget(self.mpl_toolbar)
        hbox2.setAlignment(self.mpl_toolbar, Qt.AlignLeft)
        
        # Search widget
        search_widget = QVBoxLayout()
        search_widget.addWidget(search_label)
        search_widget.addWidget(self.search_box)
            
        #hbox_plot = QHBoxLayout()
        #hbox_plot.addWidget(self.canvas)
        #hbox_plot.addWidget(self.analysis_box)
        #hbox_plot.addStretch(1)
        
        #VMainBox = QVBoxLayout()
        #VMainBox.addLayout(hbox_plot)
        #VMainBox.addWidget(self.mpl_toolbar)
        #VMainBox.setAlignment(self.mpl_toolbar, Qt.AlignLeft)
        #VMainBox.addLayout(hbox)

        grid = QGridLayout()
        grid.setSpacing(10)
        
        #                   (row-i,col-i,row-expand,col-expand)     
        grid.addWidget(self.canvas,      0,0 ,2,1)
        grid.addWidget(self.analysis_box,0,1, 1,1)
        grid.addLayout(search_widget,    1,1, 1,1)
        grid.addWidget(self.refspc_table,2,1, 2,1)
        grid.addLayout(hbox2,            2,0 ,1,1)
        grid.addLayout(hbox,             3,0 ,1,1)
        grid.setAlignment(hbox,   Qt.AlignLeft)
        grid.setAlignment(hbox2,  Qt.AlignLeft)
        
        self.main_frame.setLayout(grid)
        
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        # /File
        self.file_menu = self.menuBar().addMenu("&File")
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        xmgr_export_action = self.create_action("&Export to xmgrace", slot=self.xmgr_export, 
            shortcut="Ctrl+E", tip="Export to xmgrace format")
        spc_import_action = self.create_action("&Import plot", slot=self.open_plot, 
            shortcut="Ctrl+N", tip="Import spectrum")
        # Now place the actions in the menu
        self.add_actions(self.file_menu, 
            (load_file_action, xmgr_export_action, spc_import_action, None, quit_action))
        
        # /Analyze
        self.anlyze_menu = self.menuBar().addMenu("&Analyze")
        momenta_action = self.create_action("&Momenta", 
            slot=self.compute_moments, 
            tip='Compute moments')
        self.add_actions(self.anlyze_menu, (momenta_action,))
        
        # /Manipulate
        self.manip_menu = self.menuBar().addMenu("&Manipulate")
        self.shiftref_action = self.create_action("&Shift to simulated", 
            slot=self.shift_to_simulated, 
            tip='Shift reference spectrum to match the simulated one')
        self.scaleref_action = self.create_action("&Scale to simulated", 
            slot=self.scale_to_simulated, 
            tip='Scale reference spectrum to match the simulated one')
        self.add_actions(self.manip_menu, (self.shiftref_action,self.scaleref_action))
        # Initially, the reference spectrum is not available
        self.shiftref_action.setEnabled(False)
        self.scaleref_action.setEnabled(False)
        
        # /About
        self.help_menu = self.menuBar().addMenu("&Help")
        man_action = self.create_action("&Instructions", 
            shortcut='F1', slot=self.on_man, 
            tip='Short manual')
        about_action = self.create_action("&About", 
            shortcut='F2', slot=self.on_about, 
            tip='About the app')
        self.add_actions(self.help_menu, (man_action,about_action))
        
    # The following two wrapper functions came from:
    # https://github.com/eliben/code-for-blog/blob/master/2009/qt_mpl_bars.py
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)
                

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
    
    
#==========================================================
# GLOBAL FUNCTIONS
#==========================================================
# FILE READERS
def read_fort21(fort21file,MaxClass):
    """
    Function to extract transtion infor from fort.21 file
    The content is taken from the standard version and may
    need some polish
    
    Arguments:
     MaxClass (int): maximum class to be loaded (0 to 7)
                     Hot bands are always loaded
                     
    Returns:
     List of list of transitions: [ [C0], [C1], .., [Chot] ]
    """
    # Open and read file
    tr=[]
    print "Loading transitions (fort.21)..."
    try:
        f = open(fort21file,'r')
    except:
        exit("ERROR: Cannot open file 'fort.21'")
        
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
        elif "*************************************************" in line:
            nclass0 = len(tr)
            itrans  = nclass0-1
            break

    nhot = 0
    nclass1 = 0
    nclass2 = 0
    nclass3 = 0
    nclass4 = 0
    nclass5 = 0
    nclass6 = 0
    nclass7 = 0
    fcclass = 0
    nclass  = 0
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
        elif 'state 2 = GROUND' in line and fcclass<=MaxClass:
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
    
    # Load filter transitions if required
    loadC=True
    total_transitions = 0
    nclass_list = [nclass0,nclass1,nclass2,nclass3,nclass4,nclass5,nclass6,nclass7]
    print 'Transitions read:'
    print ' Class     N. trans.         Load?  '
    for i,nclass in enumerate(nclass_list):
        if MaxClass<i: 
            loadC=False
        total_transitions += nclass
        print ' C{0}        {1:5d}             {2}   '.format(i,nclass,loadC)
        if MaxClass<i: 
            nclass_list[i]=0
    print     ' Hot       {0:5d}             {1}   '.format(nhot,True)
    nclass_list.append(nhot)
    total_transitions += nhot
    print 'Total transitions : ',(total_transitions)
    print 'Loaded transitions: ',(itrans+1)
    print ''
    #========== Done with fort.21 ====================================

    # This is a conversion from old stile reader to class_list
    # maybe it'd be better to change the code above
    class_list = []
    for nclass in nclass_list:
        class_list.append([tr.pop(0) for i in range(nclass)])

    return class_list


def read_spc_xy(filename,fromsection=None):
    """
    Function to read fort.22, which contains the bins to reconstruct
    the convoluted spectrum_sim.
    This is a simple [x,y] file
    An alternative would be to use the np.loadtxt. But for the 
    moment we leave it like that
    """
    # Open and read file
    print "Loading spectral data from '"+filename+"'..."
    try:
        f = open(filename,'r')
    except:
        sys.exit("ERROR: Cannot open file '"+filename+"'")
    i = 0
    x = []
    y = []
    
    if fromsection:
        for line in f:
            if fromsection in line:
                break
    
    for line in f:
        data = line.split()
        try:
            x.append(float(data[0]))
            y.append(float(data[1]))
        except:
            continue

    f.close()

    return x,y

# CONVOLUTION
def convolute(spc_stick,npoints=1000,hwhm=0.1,broad="Gau",input_bins=False):
    """
    Make a Gaussian convolution of the stick spectrum
    The spectrum must be in energy(eV) vs Intens (LS?)
    
    Arguments:
    spc_stick  list of list  stick spectrum as [x,y]
               list of array
    npoints    int           number of points (for the final graph)
                             Can be a bit more if input_bins is False
    hwhm       float         half width at half maximum
    
    Retunrs a list of arrays [xconv,yconv]
    """
    x = spc_stick[0]
    y = spc_stick[1]
   
    # ------------------------------------------------------------------------
    # Convert discrete sticks into a continuous function with an histogram
    # ------------------------------------------------------------------------
    # (generally valid, but the exact x values might not be recovered)
    # Make the histogram for an additional 20% (if the baseline is not recovered, enlarge this)
    extra_factor = 0.2
    recovered_baseline=False
    sigma = hwhm / np.sqrt(2.*np.log(2.))
    while not recovered_baseline:
        if input_bins:
            # Backup npoints
            npts = npoints
            npoints = len(x)
            xhisto = x
            yhisto = y
            width = (x[1] - x[0])
        else:
            extra_x = (x[-1] - x[0])*extra_factor
            yhisto, bins =np.histogram(x,range=[x[0]-extra_x,x[-1]+extra_x],bins=npoints,weights=y)
            # Use bin centers as x points
            width = (bins[1] - bins[0])
            xhisto = bins[0:-1] + width/2
        
        # ----------------------------------------
        # Build Gaussian (centered around zero)
        # ----------------------------------------
        dxgau = width
        # The same range as xhisto should be used
        # this is bad. We can get the same using 
        # a narrower range and playing with sigma.. (TODO)
        if npoints%2 == 1:
            # Zero is included in range
            xgau_min = -dxgau*(npoints/2)
            xgau_max = +dxgau*(npoints/2)
        else:
            # Zero is not included
            xgau_min = -dxgau/2. - dxgau*((npoints/2)-1)
            xgau_max = +dxgau/2. + dxgau*((npoints/2)-1)
        xgau = np.linspace(xgau_min,xgau_max,npoints)
        if broad=="Gau":
            ygau = np.exp(-xgau**2/2./sigma**2)/sigma/np.sqrt(2.*np.pi)
        elif broad=="Lor":
            ygau = hwhm/(xgau**2+hwhm**2)/np.pi
        else:
            sys.exit("ERROR: Unknown broadening function: "+broad)
        
        # ------------
        # Convolute
        # ------------
        # with mode="same", we get the original xhisto range.
        # Since the first moment of the Gaussian is zero, 
        # xconv is exactly xhisto (no shifts)
        yconv = np.convolve(yhisto,ygau,mode="same")
        xconv = xhisto

        # Check baseline recovery (only with automatic bins
        if yconv[0] < yconv.max()/100.0 and yconv[-1] < yconv.max()/100.0:
            recovered_baseline=True
        if input_bins:
            recovered_baseline=True
            # If the input_bins are larger than npts, then reduce the grid to npts
            if (len(xconv) > npts):
                skip = len(xconv)/npts + 1
                x = xconv[0::skip]
                y = yconv[0::skip]
                xconv = x
                yconv = y

        extra_factor = extra_factor + 0.05

    return [xconv,yconv]
    
# GENERATE XMGR
def export_xmgrace(filename,ax,sticks,labs,ax2=None,specs=None):
    """
    DESCRIPTION
    ------------
    Function to convert the current data into a xmgrace plot, including
    labels currently on the screen
    
    Variables
    * filename string : path to the file to save the export
    * ax,         mpl.axes        : graphic info
    * class_list  list of vlines  : stick spectra (classes)
    * specs       Line2D          : convoluted spectrum
    * labs        dict            : labels in xmgr format (arg: labels as annotation Class)
    """

    f = open(filename,'w')
    
    print >> f, "# XMGRACE CREATED BY FCC_ANALYZER"
    print >> f, "# Only data and labels. Format will"
    print >> f, "# be added by your default xmgrace"
    print >> f, "# defaults (including colors, fonts...)"
    print >> f, "# Except the followins color scheme:"
    print >> f, '@map color 0  to (255, 255, 255), "white"'
    print >> f, '@map color 1  to (0, 0, 0), "black"'
    print >> f, '@map color 2  to (0, 0, 255), "blue"'
    print >> f, '@map color 3  to (255, 0, 0), "red"'
    print >> f, '@map color 4  to (0, 139, 0), "green4"'
    print >> f, '@map color 5  to (0, 255, 255), "cyan"'
    print >> f, '@map color 6  to (255, 0, 255), "magenta"'
    print >> f, '@map color 7  to (188, 143, 143), "brown"'
    print >> f, '@map color 8  to (100, 0, 100), "pink"'
    print >> f, '@map color 9  to (255, 165, 0), "orange"'
    print >> f, '@map color 10 to (255, 255, 0), "yellow"'
    print >> f, '@map color 11 to (220, 220, 220), "grey"'
    print >> f, '@map color 12 to (0, 255, 0), "green"'
    print >> f, '@map color 13 to (148, 0, 211), "violet"'
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
    print >> f, "@    view 0.180000, 0.150000, 1.15, 0.92"
    #Get plotting range from mplt
    x=ax.get_xbound()
    y=ax.get_ybound()
    print >> f, "@    world ",x[0],",",y[0],",",x[1],",",y[1]
    #Get xlabel from mplt
    print >> f, "@    xaxis  label \""+ax.get_xlabel()+"\""
    print >> f, "@    yaxis  label \""+ax.get_ylabel()+"\""
    if ax2:
        # Position
        print >> f, "@    yaxis  label place opposite"
        print >> f, "@    yaxis  ticklabel place opposite"
        print >> f, "@    yaxis  tick place opposite"
    # Char sizes
    print >> f, "@    xaxis  ticklabel char size 1.250000"
    print >> f, "@    yaxis  ticklabel char size 1.250000"
    print >> f, "@    xaxis  label char size 1.500000"
    print >> f, "@    yaxis  label char size 1.500000"
    #Get tick spacing from mplt
    x=ax.get_xticks()
    y=ax.get_yticks()
    print >> f, "@    xaxis  tick major", x[1]-x[0]
    print >> f, "@    yaxis  tick major", y[1]-y[0]
    #Legend
    print >> f, "@    legend char size 1.250000"
    print >> f, "@    legend loctype view"
    print >> f, "@    legend 0.95, 0.9"
    #Now include data
    label_list = ['0-0']+[ 'C'+str(i) for i in range(1,8) ]+['Hot']
    color_list = [ 1, 2, 3, 4, 5, 6, 7, 8, 9]
    k=-1
    ymax = -999.
    ymin =  999.
    for iclass in range(9):
        if (len(sticks[iclass]) == 0):
            continue
        k += 1
        x = np.array([ sticks[iclass][i].DE        for i in range(len(sticks[iclass])) ])
        y = np.array([ sticks[iclass][i].intensity for i in range(len(sticks[iclass])) ])
        print >> f, "& %s"%(label_list[iclass])
        print >> f, "@type bar"
        print >> f, "@    s"+str(k)," line type 0"
        print >> f, "@    s"+str(k)," legend  \"%s\""%(label_list[iclass])
        print >> f, "@    s"+str(k)," symbol color %s"%(color_list[iclass])
        for i in range(len(x)):
            print >> f, x[i], y[i]
        ymax = max(ymax,y.max())
        ymin = max(ymin,y.min())
            
    # Convoluted spectrum. scaled TODO: use the alternative axis
    # This can be done either with another Graph or using alt-x axis
    xs_max = max([ max(max(s.get_xdata()),abs(min(s.get_xdata()))) for s in specs ])
    ys_max = max([ max(max(s.get_ydata()),abs(min(s.get_ydata()))) for s in specs ])
    if not ax2:
        k += 1
        if abs(ymax) > abs(ymin):
            scale_factor = ymax/ys_max
        else:
            scale_factor = abs(ymax)/abs(ys_max)
        print "Scale factor applied to convoluted spectrum %s"%(scale_factor)
    else:
        k = 0
        print >> f, "@g1 on"
        print >> f, "@g1 hidden false"
        print >> f, "@with g1"
        # Set a large view
        print >> f, "@    view 0.180000, 0.150000, 1.15, 0.92"
        #Get plotting range from mplt
        x=ax2.get_xbound()
        y=ax2.get_ybound()
        print >> f, "@    world ",x[0],",",y[0],",",x[1],",",y[1]
        #Get labels from mplt
        print >> f, "@    yaxis  label \""+latex2xmgrace(ax2.get_ylabel())+"\""
        # Char sizes
        print >> f, "@    xaxis  ticklabel char size 1.250000"
        print >> f, "@    yaxis  ticklabel char size 1.250000"
        print >> f, "@    xaxis  label char size 1.500000"
        print >> f, "@    yaxis  label char size 1.500000"
        #Set axis ticks on the left
        print >> f, "@    xaxis  off"
        print >> f, "@    yaxis  tick out"
        print >> f, "@    yaxis  label place normal"
        print >> f, "@    yaxis  ticklabel place normal"
        print >> f, "@    yaxis  tick place normal"
        #Get tick spacing from mplt
        y=ax2.get_yticks()
        print >> f, "@    yaxis  tick major", y[1]-y[0]
        #Legend
        print >> f, "@    legend char size 1.250000"
        print >> f, "@    legend loctype view"
        print >> f, "@    legend 0.8, 0.9"
        # No scale
        scale_factor = 1.0
        
    # Note that y holds the adress of the data! So, if we change it, we change the original data!
    #y *= scale_factor
    color_list = [ 1, 18]
    style_list = [ 3, 1 ]
    for i,s in enumerate(specs):
        leg_label=specs[i].get_label()
        x = s.get_xdata()
        y = s.get_ydata()
        print >> f, "& "+leg_label
        print >> f, "@type xy"
        print >> f, "@    s"+str(k)," line type 1"
        print >> f, "@    s"+str(k)," line linestyle %s"%(style_list[i])
        print >> f, "@    s"+str(k)," line color %s"%(color_list[i])
        print >> f, "@    s"+str(k)," legend  \""+leg_label+"\""
        for j in range(len(x)):
            print >> f, x[j], y[j]*scale_factor
        k += 1
            
    f.close()

def latex2xmgrace(string):
    """
    Transform latex mathenv entries to xmgrace format string
    """
    if not '$' in string:
        return string
    
    # Dictionary of latex commands and xmgrace translation
    ltx2xmgr=dict()
    ltx2xmgr[r'\alpha']     =r'\xa\f{}'
    ltx2xmgr[r'\beta']      =r'\xb\f{}'
    ltx2xmgr[r'\gamma']     =r'\xg\f{}'
    ltx2xmgr[r'\epsilon']   =r'\xe\f{}'
    ltx2xmgr[r'\varepsilon']=r'\xe\f{}'
    ltx2xmgr[r'\Delta']     =r'\xD\f{}'
    
    str_parts = string.split('$')
    str_math = str_parts[1::2]
    str_text = str_parts[0::2]
    
    for i,item in enumerate(str_math):
        if '\\' in item:
            # There is a command to replace
            pattern=r'\\[A-Za-z_]+'
            for ltx_cmd in re.findall(pattern,item):
                xmgr_cmd=ltx2xmgr[ltx_cmd]
                item = item.replace(ltx_cmd,xmgr_cmd)
        if '^' in item:
            # There is a supperscript
            pattern=r'\^[\{\}\-0-9]'
            for exp in re.findall(pattern,item):
                item = item.replace('^','\\S')
                item = item.replace('{','')
                item = item.replace('}','')
                item += '\\N'
            
        str_math[i] = item
        
    # Join list
    for i,v in enumerate(str_math):
        str_text.insert(2*i+1,v)
        
    string=""
    for char in str_text:
        string += char
        
    return string

    

# INPUT PARSER
def get_args():
    
    # Options and their defaults 
    final_arguments = dict()
    final_arguments["-maxC"]="7"
    final_arguments["-type"]="(interactive)"
    final_arguments["--test"]=False
    final_arguments["-h"]=False
    final_arguments["-stick"]="int"
    # Description of the options
    arg_description = dict()
    arg_description["-maxC"] ="Maximum class to load"
    arg_description["-type"] ="Type of calculation [abs|emi|ecd|cpl]"
    arg_description["--test"]="Load spectra and quit"
    arg_description["-stick"]="Load intensity or square FC [int|fc]"
    arg_description["-h"]    ="Show this help"
    # Type for arguments
    arg_type = dict()
    arg_type["-maxC"] ="int"
    arg_type["-type"] ="char"
    arg_type["--test"]="-"
    arg_type["-h"]    ="-"
    arg_type["-stick"]="char"
    
    # Get list of input args
    input_args_list = []
    iarg = -1
    for s in sys.argv[1:]:
        # get -flag [val] arguments 
        if s[0]=="-":
            iarg=iarg+1
            input_args_list.append([s])
        else:
            input_args_list[iarg].append(s)
            
    # Transform into dict. Associtaing lonely flats to boolean   
    input_args_dict=dict()
    for input_arg in input_args_list:
        if len(input_arg) == 1:
            # Boolean option. Can be -Bool or -noBool
            input_arg.append(True)
            if input_arg[0][1:3] == "no":
                input_arg[0] = "-" + input_arg[0][3:]
                input_arg[1] = not input_arg[1]
        elif len(input_arg) != 2:
            raise BaseException("Sintax error. Too many arguments")

        input_args_dict[input_arg[0]] = input_arg[1]
    
    for key,value in input_args_dict.iteritems():
        # Check it is allowed
        isValid = final_arguments.get(key,None)
        if isValid is None:
            raise BaseException("Sintax error. Unknown label: " + key)
        # If valid, update final argument
        final_arguments[key]=value
        
    if final_arguments.get("-h"):
        
        print """
 ----------------------------------------
           FCclasses analyzer
   A GUI to analyze FCclasses output 
 ----------------------------------------
        """
        helptext()
        print "    Options:"
        print "    --------"
        print '      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format("Flag","Type","Description","Value")
        print '      {0:-<10}  {1:-^4}  {2:-<41}  {3:-<7}'.format("","","","")
        for key,value in final_arguments.iteritems():
            descr = arg_description[key]
            atype = arg_type[key]
            #atype=str(type(value)).replace("<type '","").replace("'>","")
            print '      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format(key, atype, descr, str(value))
        print ""
        
        sys.exit()
        
    return final_arguments
        
    
    
def main():
    app = QtGui.QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()


