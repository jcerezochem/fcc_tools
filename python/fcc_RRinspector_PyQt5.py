#!/usr/bin/env python
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
# PyQt4
# from PyQt4.QtCore import *
# from PyQt4.QtGui import *
# from PyQt4 import QtCore, QtGui
# PyQt5
from PyQt5 import QtGui, QtCore, uic, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication, QAction, QWidget, QLineEdit, QPushButton, QComboBox, \
                            QSlider, QCheckBox, QLabel, QFrame, QTextEdit, QTableWidget, QTableWidgetItem,  \
                            QVBoxLayout, QHBoxLayout, QGridLayout, QInputDialog, QFileDialog, QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.QtCore import *

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
#from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib
import re

try:
    import fcc_version_tag as version_tag
except:
    class version_tag:
        COMMIT="Untracked"
        DATE="No date"

stick_type = 'fc'

def helptext():
    print("""
    Description:
    ------------
    This python script parses fort.21/Assignments.dat retrieving the informaition 
    about the transitions and plot them using matplotlib. An interactive plot is 
    generated which allow to assign the stick transition by clicking on them.
    The script can also export the plot to xmgrace format. The convoluted spectrum
    is generated from fort.22/Bin_Spectrum.dat (otherwise, it is read from fort.18)
    
    Instructions:
    ------------
    General call:
    
        fcc_analyzer_PyQt5 [options]
    
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
    """)

class SpcConstants:
    exp = {"opa":1,"ecd":1,"emi":3,"cpl":3,'mcd':1,'tpa':2,'tpcd':2}
    factor = {"opa":703.30,"ecd":20.5288,"emi":1063.055,"cpl":4252.216,'mcd':5.98442e-3,'tpa':8.35150e-4,'tpcd':8.35150e-4}
    # JC: 'mcd':-5.98442e-3 (should be negative), but then it seems not consistent with fcc3
    # The factors already include the conversion between eV <-> au
    # to handle some issues in the emission Lineshape (need to used
    # (27.2116) factor instead of (27.2116)^3)
    #factor = {"opa":25.8459,"ecd":0.75441,"emi":39.0662,"cpl":156.265}


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
        msg = """ Transition: %s 
  =========================
  MotherState: %s
  FC class:    %s 
  Einit(eV):   %s 
  Efin (eV):   %s 
  DE   (eV):   %s 
  DE-00(cm-1): %s 
  Intensity:   %s 
  FCfactor:    %s 
  INDEX:       %s
        """%(transition,
             self.motherstate,
             self.fcclass,
             self.einit,
             self.efin,
             self.DE,
             self.DE00cm,
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
        self.setWindowTitle('FCclasses RR inspector')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        # Initialize additional items
        # Broadening info
        self.broadening="Lor"
        self.update_hwhm_from_slider(UpdateConvolute=False)
        # UpdateW (required to get the right value from textbox)
        self.updateW = True
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
        if os.path.isfile('RR_Spectrum_2D.dat'):
            path=""
            RR2D_file='RR_Spectrum_2D.dat'
        else:
            file_choices = r"RR_Spectrum_2D.dat (RR_Spectrum_2D.dat);; All files (*)"
            path = QFileDialog.getOpenFileName(self, 
                            'Set the location of FClasses output', '', 
                            file_choices)
            # Management of QFileDialog output is different in PyQt5
            #  * Do not use unicode() to wrap the call
            #  * It is now an array. Take first value
            path = path[0]
            if not path or not 'RR_Spectrum_2D.dat' in path:
                return
            path = path.replace('RR_Spectrum_2D.dat',"")
            RR2D_file='RR_Spectrum_2D.dat'
        ## Parse command line args
        #MaxClass = cml_args.get("-maxC")
        #MaxClass = int(MaxClass)
        #self.spc_type = cml_args.get("-type")
        #self.spc_type = str(self.spc_type).lower()
        #calc_type_list = ("opa","emi","ecd","cpl","mcd","tpa","tpcd")
        #if self.spc_type not in calc_type_list:
            ## Get spc_type if not given on input or was wrong
            #get_option,self.spc_type = self.start_assistant()
            #if not get_option:
                #sys.exit()
                
        # Get Ev from RR_Spectrum_VertE.dat
        if os.path.isfile(path+'RR_Spectrum_VertE.dat'):
            with open(path+'RR_Spectrum_VertE.dat') as f:
                line = f.readline()
                line = f.readline()
                line = line.split('maximum')[-1]
                Ev = float(line.split()[0]) / 1.23981e-4
        else:
            Ev = -1

        # Read all data from file
        data = np.loadtxt(path+RR2D_file)
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]
        # Get size of x(vib bands) and y(indicend freqs)
        nx = len(np.where(y==y[0])[0])
        ny = len(np.where(x==x[0])[0])
        # Expose all data to the class methods
        self.X = x.reshape(nx,ny)
        self.Y = y.reshape(nx,ny)
        self.Z = z.reshape(nx,ny)
        
        # Get vector with wI
        self.wI = self.Y[0,:]
        if Ev>0:
            w_ind = (np.abs(self.wI - Ev)).argmin()
        else:
            w_ind = int(len(self.wI)/2)
        w = self.wI[w_ind]
        self.winc_box.setText(str(round(w,2)))
        self.winc_slider.setRange(1, len(self.wI))
        #self.update_incident_freq_from_slider()
        
        # First selections
        i0=1
        self.xbin = self.X[i0:,100]
        self.ybin = self.Z[i0:,100]
        
        # This is the load driver
        #self.load_sticks()
        self.load_convoluted()
        self.set_axis_labels()
        #self.load_legend()
        
        # We call update to get fix the slider
        self.update_incident_freq_from_textbox()
        
        if cml_args.get("--test"):
            sys.exit()
        
        
    #==========================================================
    # FUNCTIONS CONNECTED WITH THE MENU OPTIONS
    #==========================================================
    def save_plot(self):
        file_choices = r"Portable Network Graphics (*.png) (*.png);; All files (*)"
        
        path = QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices)
        # Management of QFileDialog output is different in PyQt5
        #  * Do not use unicode() to wrap the call
        #  * It is now an array. Take first value
        path = path[0]
        if path:
            self.canvas.print_figure(path, dpi=100)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
            
    def open_plot(self):
        file_choices = r"Data file (*.dat) (*.dat);; All files (*)"
        
        path = QFileDialog.getOpenFileName(self, 
                        'Open spectrum', '', 
                        file_choices)
        # Management of QFileDialog output is different in PyQt5
        #  * Do not use unicode() to wrap the call
        #  * It is now an array. Take first value
        path = path[0]
        if path:
            self.statusBar().showMessage('Opened %s' % path, 2000)    
            x,y = read_spc_xy(path)
            x = np.array(x)
            y = np.array(y)
            
            # Update Table
            self.refspc_table.setItem(1,1, QTableWidgetItem(path.split('/')[-1]))
            cell = self.refspc_table.item(1,1)
            # PyQt4
            #cell.setTextColor(Qt.black)
            # PyQt5
            cell.setForeground(Qt.black)
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
        calc_type_list = ("opa","emi","ecd","cpl","mcd","tpa","tpcd")
                 
        spc_type, ok = QInputDialog.getItem(self, "Select type of calculation", 
                            "Type of calculation", calc_type_list, 0, False)
        
        return ok, str(spc_type)
    
    
    def spc_import_assitant_Xaxis(self):
        unit_list = ("cm^-1", "eV", "nm")
                 
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
        
        path = QFileDialog.getSaveFileName(self, 
                        'Export to file', '', 
                        file_choices)
        # Management of QFileDialog output is different in PyQt5
        #  * Do not use unicode() to wrap the call
        #  * It is now an array. Take first value
        path = path[0]
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
            
    def dat_export(self):
        # File Dialog
        file_choices = "dat file (*.dat) (*.dat);; All files (*)"
        
        path = QFileDialog.getSaveFileName(self, 
                        'Export to file', '', 
                        file_choices)
        # Management of QFileDialog output is different in PyQt5
        #  * Do not use unicode() to wrap the call
        #  * It is now an array. Take first value
        path = path[0]
        if path:
            # Export simulated spectrum only:
            spc = self.spectrum_sim
            f = open(path,'w')
            x = spc[0].get_xdata()
            y = spc[0].get_ydata()
            for j in range(len(x)):
                print(x[j], y[j], file=f)
            f.close()
    
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
    def load_sticks(self,wi_index=0):
        """ 
        Load stick spectra for all classes (C0,C1...,CHot)
        - The spectra objects are stored in a list: self.stickspc
        """
        # clear the axes and redraw the plot anew
        # 
        self.axes.set_title('RR stick spectrum from $\mathcal{FC}classes$',fontsize=18)
        self.axes.set_xlabel('Energy (cm-1)',fontsize=16)
        self.axes.set_ylabel('Stick Intensity',fontsize=16)
        self.axes.tick_params(direction='out',top=False, right=False)
        
        
        #Plotting sticks and store objects
        # Set labels and colors
        label_list = ['Rayleigh','Stokes']
        color_list = ['k', 'b']
        
        # Get data

        #Inialize variables
        self.stickspc = []
        xmin =  999.
        xmax = -999
        # Rayleigh
        x = self.X[0,wi_index]
        y = self.Y[0,wi_index]
        z = self.Z[0,wi_index]
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
        lns  = self.spectrum_sim
        labs = [l.get_label() for l in lns]
        # Set the location according to the type of calculation
        position="upper right"
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
        self.axes.set_title('vRR spectrum from $\mathcal{FC}classes$',fontsize=18)
        self.axes2.set_ylabel(r'$d\sigma/d\Omega$ (cm$^2$/sr)',fontsize=16)
        self.axes2.set_xlabel('Frequency (cm$^{-1}$)',fontsize=16)
        #self.axes.tick_params(direction='out',top=False, right=False)
        
            
    def load_convoluted(self):
        txt = str(self.broadbox.text())
        hwhm = float(txt)
        fixaxes = self.fixaxes_cb.isChecked()
        #Convolution
        xc,yc = convolute([self.xbin,self.ybin],hwhm=hwhm,broad=self.broadening)
        # Plot convoluted
        self.spectrum_sim = self.axes2.plot(xc,yc,color='k',label="RR spectrum")
        if not fixaxes:
            self.rescale_yaxis()
        
        self.canvas.draw()
        
        
    def update_convolute(self):
        txt = str(self.broadbox.text())
        hwhm = float(txt)
        fixaxes = self.fixaxes_cb.isChecked()
        
        #Convolution (in energy(eV))
        xc,yc = convolute([self.xbin,self.ybin],hwhm=hwhm,broad=self.broadening)
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
        n = 1 # w_s * w_I³ ?
        y /= x**n 
        x = x + shift
        # If Intensity, set back from Lineshape
        y *= x**n
        
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
        #if type(event.artist) == self.LineCollectionType:
            #self.select_stick(event)
        if type(event.artist) == self.AnnotationType:
            self.active_label = event.artist
        elif type(event.artist) == self.Line2DType:
            #self.del_stick_marker()
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
                    agrlabel = agrlabel+str(tr.init[i])+'\\S'+str(tr.qinit[i])+'\\N,'
                if len(tr.init) == 0:
                    label=label+"0,"
                    agrlabel=agrlabel+"0,"
                label = label[0:-1]+'-'
                agrlabel = agrlabel[0:-1]+'-'
            for i in range(0,len(tr.final)):
                label = label+str(tr.final[i])+'^{'+str(tr.qfinal[i])+'},'
                agrlabel = agrlabel+str(tr.final[i])+'\\S'+str(tr.qfinal[i])+'\\N,'
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
                print("This label was already defined")
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
        #fixlegend = not self.fixlegend_cb.isChecked()
        #self.legend.set(draggable=fixlegend)
        pass
        
    
    def update_hwhm_from_slider(self,UpdateConvolute=True):
        hwhmmin = 1
        hwhmmax = 26
        slidermin = 1   # this is not changed
        slidermax = 101 # this is not changed
        hwhm = float((hwhmmax-hwhmmin)/(slidermax-slidermin) * (self.slider.value()-slidermin) + hwhmmin)
        hwhm = round(hwhm,3)
        self.broadbox.setText(str(hwhm))
        if (UpdateConvolute):
            self.update_convolute()
            
        
    def update_incident_freq_from_slider(self):
        wmin = self.wI.min()
        wmax = self.wI.max()
        slidermin = 1   # this is not changed
        slidermax = len(self.wI) # this is not changed
        w = float((wmax-wmin)/(slidermax-slidermin) * (self.winc_slider.value()-slidermin) + wmin)
        # Get the closest value
        w_ind = (np.abs(self.wI - w)).argmin()
        w = self.wI[w_ind]
        # Update slider
        sliderval = int((slidermax-slidermin)/(wmax-wmin) * (w-wmin) + slidermin)
        sliderval = min(sliderval,slidermax)
        sliderval = max(sliderval,slidermin)
        self.winc_slider.setValue(sliderval)
        # Update textbox (this will trigger the update)
        if self.updateW:
            self.winc_box.setText(str(round(w,2)))
            # Get new bins
            if self.showRayleigh_cb.isChecked():
                i0 = 0
            else:
                i0 = 1
            self.xbin = self.X[i0:,w_ind] * float(self.scale_factor_box.text())
            self.ybin = self.Z[i0:,w_ind]
            self.update_convolute()
            
        
    def update_hwhm_from_textbox(self):
        hwhmmin = 1
        hwhmmax = 25
        slidermin = 1   # this is not changed
        slidermax = 100 # this is not changed
        msg = str(self.broadbox.text())
        hwhm = float(msg)
        sliderval = int((slidermax-slidermin)/(hwhmmax-hwhmmin) * (hwhm-hwhmmin) + slidermin)
        sliderval = min(sliderval,slidermax)
        sliderval = max(sliderval,slidermin)
        self.slider.setValue(sliderval)
        self.update_convolute()
        
        
    def update_incident_freq_from_textbox(self):
        wmin = self.wI.min()
        wmax = self.wI.max()
        slidermin = 1   # this is not changed
        slidermax = len(self.wI) # this is not changed
        msg = str(self.winc_box.text())
        w = float(msg)
        # Get the closest value
        w_ind = (np.abs(self.wI - w)).argmin()
        w = self.wI[w_ind]
        self.winc_box.setText(str(round(w,2)))
        # Update slider
        sliderval = int((slidermax-slidermin)/(wmax-wmin) * (w-wmin) + slidermin)
        sliderval = min(sliderval,slidermax)
        sliderval = max(sliderval,slidermin)
        self.updateW = False
        self.winc_slider.setValue(sliderval)
        self.updateW = True
        # Get new bins
        if self.showRayleigh_cb.isChecked():
            i0 = 0
        else:
            i0 = 1
        self.xbin = self.X[i0:,w_ind] * float(self.scale_factor_box.text())
        self.ybin = self.Z[i0:,w_ind]
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
        
        
    #def update_incident_freq(self):
        #current_incident_freq = self.incident_freq
        #self.data_type = self.select_data_type.currentText()
        
        #if current_data_type != self.data_type:
            #factor = SpcConstants.factor[self.spc_type]
            #n = SpcConstants.exp[self.spc_type]
            #if self.with_fort22:
                #if self.data_type == "Lineshape":
                    ## Division x/27.2116 could be included in the factor
                    #self.ybin /= (self.xbin/27.2116)**n * factor
                #elif self.data_type == "Intensity":
                    ## Division x/27.2116 could be included in the factor
                    #self.ybin *= (self.xbin/27.2116)**n * factor
            #else:
                #if self.spectrum_sim:
                    #x = self.spectrum_sim[0].get_xdata()
                    #y = self.spectrum_sim[0].get_ydata()
                    #if self.data_type == "Lineshape":
                        ## Division x/27.2116 could be included in the factor
                        #y /= (x/27.2116)**n * factor
                    #elif self.data_type == "Intensity":
                        ## Division x/27.2116 could be included in the factor
                        #y *= (x/27.2116)**n * factor
                    #self.spectrum_sim[0].set_ydata(y)
            
            #if self.spectrum_ref:
                #x = self.spectrum_ref[0].get_xdata()
                #y = self.spectrum_ref[0].get_ydata()
                #if self.data_type == "Lineshape":
                    ## Division x/27.2116 could be included in the factor
                    #y /= (x/27.2116)**n * factor
                #elif self.data_type == "Intensity":
                    ## Division x/27.2116 could be included in the factor
                    #y *= (x/27.2116)**n * factor
                #self.spectrum_ref[0].set_ydata(y)
                
            #self.set_axis_labels()
            #if self.with_fort22:
                #self.update_convolute()
            #else:
                #fixaxes = self.fixaxes_cb.isChecked()
                #if not fixaxes:
                    #self.rescale_yaxis()
                #self.canvas.draw()
                
                
    def table_buttons_action(self,i,j):
        # "X" button: clear spectrum
        # "T" reset shift/scale
        if not self.spectrum_ref:
            return
        if (i,j) == (1,2):
            clear_msg = "Clear reference spectrum?"
            reply = QMessageBox.question(self, 'Clear Spectrum', 
                         clear_msg, QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                return            
            self.spectrum_ref[0].remove()
            self.spectrum_ref = None
            self.canvas.draw()
            celllabel = ["No","-","-"]
            for i,j in [(1,1),(2,1),(3,1)]:
                self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
                cell = self.refspc_table.item(i,j)
                # PyQt4 (function deprecated in PyQt5)
                #cell.setTextColor(Qt.black)
                # PyQt5
                cell.setForeground(Qt.black)
                cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
            # Disable manipulations
            self.shiftref_action.setEnabled(False)
            self.scaleref_action.setEnabled(False)

        elif (i,j) == (2,2):
            # Tare the shift
            msg = "Reset shift to current value?"
            reply = QMessageBox.question(self, 'Reset Shift', 
                         msg, QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                return      
            self.ref_shift = 0.0
            self.refspc_table.setItem(2,1, QTableWidgetItem(str(self.ref_shift)))
            
        elif (i,j) == (3,2):
            # Tare the scale
            msg = "Reset scale to current value?"
            reply = QMessageBox.question(self, 'Reset Scale', 
                         msg, QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
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
        return
        #command = str(self.search_box.text())
        
        ## First, check if the used is only asking for help
        #if command == 'help':
            #msg = """ 
#Examples
#---------
#1(1),2(1)   - From ground to 1(1),2(1)
#1(1),2(P)   - Progression over mode 2
#1(1)-->1(1) - Hot transition
#1(1)-->0    - M-0 transition
#0-->0       - 0-0 transition
        #"""
            #self.analysis_box.setText(msg)
            #return
        
        ## Start processing the command
        #self.statusBar().showMessage('Searching transition '+command, 2000)
        ## Here, use regexp to set the validity of the command
        #pattern = r'(( )*(([0-9pP]+\([0-9]+\)( )*\,( )*)*( )*[0-9pP]+\([0-9]+\)|0)( )*\-\->( )*){0,1}(( )*([0-9]+\([0-9pP]+\)( )*\,( )*)*( )*[0-9]+\([0-9pP]+\)|0)( )*'
        #match = re.match(pattern,command)
        #if not match or match.group() != command:
        ##     [m'1(q'1),m'2(q'2)... -->] m1(q1),m2(q2)...
        ## or: 0 --> 0 
            #if self.selected:
                #self.selected.remove()
                #self.selected = None
                #self.canvas.draw()
            #self.analysis_box.setText("")
            #self.statusBar().showMessage('Invalid syntax', 2000)
            #return
        #transition_list=command.split('-->')
        #if len(transition_list) == 2:
            #state_ini,state_fin = transition_list
        #elif len(transition_list) == 1:
            #state_ini  = "0"
            #state_fin, = transition_list
            
            
        ## Initial state
        #imodes  = []
        #iquanta = []
        #if not "(" in state_ini and int(state_ini) == 0:
            #transitions = []
        #else:
            #transitions = state_ini.split(',')
        #for tr in transitions:
            #m,q = tr.split('(')
            ## If re.match is done well (as it is?) these try/except
            ## should no be needed. To be removed
            #try:
                #imodes.append(int(m))
            #except:
                #self.statusBar().showMessage('Invalid syntax', 2000)
                #return
            #try:
                #iquanta.append(int(q.replace(')','')))
            #except:
                #self.statusBar().showMessage('Invalid syntax', 2000)
                #return
        
        
        ## Final state
        #fmodes  = []
        #fquanta = []
        #pindex = []
        #if not "(" in state_fin and int(state_fin) == 0:
            #transitions = []
        #else:
            #transitions = state_fin.split(',')
        #for tr in transitions:
            #m,q = tr.split('(')
            #try:
                #fmodes.append(int(m))
            #except:
                #self.statusBar().showMessage('Invalid syntax', 2000)
                #return
            #if q.strip().upper() == 'P)':
                #progression=True
                #pindex.append(len(fquanta))
                #fquanta.append(-1)
            #else:
                #try:
                    #fquanta.append(int(q.replace(')','')))
                #except:
                    #self.statusBar().showMessage('Invalid syntax', 2000)
                    #return

        #if (len(pindex) > 1):
            #self.statusBar().showMessage('Invalid syntax', 2000)
            #return
        #elif (len(pindex) == 1):
            #pindex = pindex[0]
        #elif (len(pindex) == 0):
            #progression = False

        ## Set FCclass
        #fclass = len(fmodes)
        #iclass = len(imodes)
        #if iclass != 0:
            #fcclass = self.fcclass_list[8]
        #else:
            #fcclass = self.fcclass_list[fclass]
        
        #tr_select = []
        #i  = 0
        #take_next_tr = True
        ## To handle the cases where the even are active but not the odd ones we use the imissing counter
        #imissing = 0
        ## Sorting two list (apply change of one into the other)
        ## http://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
        #while take_next_tr:
            #take_next_tr = False
            #if progression:
                #i += 1
                #fquanta[pindex] = i
            #for tr in fcclass:
                #if iclass != 0:
                    #mi_srch,qi_srch = list(zip(*sorted(zip(imodes,iquanta))))
                    #mi_tr,qi_tr     = list(zip(*sorted(zip(tr.init[:iclass],tr.qinit[:iclass]))))
                #else:
                    #mi_srch,qi_srch = 0,0
                    #mi_tr,qi_tr     = 0,0
                #if fclass != 0:
                    ## HotClass include X->0, which has len(tr.final)=0
                    ## and crash with the zip. This is becase in tr.info()
                    ## the lists tr.init and tr.final are cut to delete the 
                    ## zeroes
                    #if len(tr.final) == 0:
                        #continue
                    #mf_srch,qf_srch = list(zip(*sorted(zip(fmodes,fquanta))))
                    #mf_tr,qf_tr     = list(zip(*sorted(zip(tr.final[:fclass],tr.qfinal[:fclass]))))
                #else:
                    #mf_srch,qf_srch = 0,0
                    #mf_tr,qf_tr     = 0,0
                #if mf_srch==mf_tr and qf_srch==qf_tr and mi_srch==mi_tr and qi_srch==qi_tr:
                    #tr_select.append(tr)
                    #imissing -= 2
                    #break
            #imissing += 1
            #if progression and imissing<2:
                #take_next_tr = True
        
        #if progression and tr_select:
            ## Get the origin of the transition
            #fclass -= 1
            #if iclass == 0:
                #fcclass = self.fcclass_list[fclass]
            #if fclass == 0:
                #tr = fcclass[0]
                #tr_select.insert(0,tr)
            #else:
                #fmodes.pop(pindex)
                #fquanta.pop(pindex)
                #for tr in fcclass:
                    #if iclass != 0:
                        #mi_srch,qi_srch = list(zip(*sorted(zip(imodes,iquanta))))
                        #mi_tr,qi_tr     = list(zip(*sorted(zip(tr.init[:iclass],tr.qinit[:iclass]))))
                    #else:
                        #mi_srch,qi_srch = 0,0
                        #mi_tr,qi_tr     = 0,0
                    #if fclass != 0:
                        #mf_srch,qf_srch = list(zip(*sorted(zip(fmodes,fquanta))))
                        #mf_tr,qf_tr     = list(zip(*sorted(zip(tr.final[:fclass],tr.qfinal[:fclass]))))
                    #else:
                        #mf_srch,qf_srch = 0,0
                        #mf_tr,qf_tr     = 0,0
                    #if mf_srch==mf_tr and qf_srch==qf_tr and mi_srch==mi_tr and qi_srch==qi_tr:
                        #tr_select.insert(0,tr)
                        #break
        
        ## If only one tr selected, set it to active_tr
        ## This allows to browse with +/- from this tr
        #if (len(tr_select)) == 1:
            #self.active_tr
            
        ## Remove stick if there was already 
        #if self.selected:
            #self.selected.remove()
            #self.selected = None
            #self.canvas.draw()
        #self.analysis_box.setText("")
                
        #if tr_select:
            #self.statusBar().showMessage('Found', 2000)
            #msg=""
            #stick_x = []
            #stick_y = []
            #for tr in tr_select:
                #stick_x.append(tr.DE)
                #stick_y.append(tr.intensity)
                #msg = msg+"\n"+tr.info()
            #zero = np.zeros(len(stick_x))
            ## Add transition info to analysis_box
            #self.analysis_box.setText(msg.strip())
            #self.selected  = self.axes.vlines(stick_x, zero, stick_y, linewidths=3,
                                      #color='yellow', visible=True, alpha=0.7)
            #self.canvas.draw()
        #else:
            #self.statusBar().showMessage('Not found', 2000)

        
        
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
        self.broadbox.setMaximumWidth(150)
        # PyQt4 (old-style signals)
        #self.connect(self.broadbox, SIGNAL('editingFinished ()'), self.update_hwhm_from_textbox)
        # PyQt5 (new-style signals)
        self.broadbox.editingFinished.connect(self.update_hwhm_from_textbox)
        
        self.winc_box = QLineEdit()
        self.winc_box.setMinimumWidth(100)
        self.winc_box.setMaximumWidth(150)
        # PyQt4 (old-style signals)
        #self.connect(self.winc_box, SIGNAL('editingFinished ()'), self.update_hwhm_from_textbox)
        # PyQt5 (new-style signals)
        self.winc_box.editingFinished.connect(self.update_incident_freq_from_textbox)
        
        self.scale_factor_box = QLineEdit()
        self.scale_factor_box.setMinimumWidth(100)
        self.scale_factor_box.setMaximumWidth(100)
        self.scale_factor_box.setText('1.000')
        self.scale_factor_box.editingFinished.connect(self.update_incident_freq_from_textbox)
        
        clean_button1 = QPushButton("&Clean(Panel)")
        # PyQt4 (old-style signals)
        #self.connect(clean_button1, SIGNAL('clicked()'), self.del_stick_marker)
        # PyQt5 (new-style signals)
        clean_button1.clicked.connect(self.del_stick_marker)
        clean_button2 = QPushButton("&Clean(Labels)")
        # PyQt4 (old-style signals)
        #self.connect(clean_button2, SIGNAL('clicked()'), self.reset_labels)
        # PyQt5 (new-style signals)
        clean_button2.clicked.connect(self.reset_labels)
        
        self.select_broad = QComboBox()
        self.select_broad.addItems(["Lor","Gau"])
        self.select_broad.currentIndexChanged.connect(self.update_broad_function)
        
        self.select_data_type = QComboBox()
        self.select_data_type.addItems(["Disabled"])
        #self.select_data_type.currentIndexChanged.connect(self.update_data_type)
        
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMaximumWidth(200)
        self.slider.setRange(1, 101)
        self.slider.setValue(30)
        self.slider.setTracking(True)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        # PyQt4 (old-style signals)
        #self.connect(self.slider, SIGNAL('valueChanged(int)'), self.update_hwhm_from_slider)
        # PyQt5 (new-style signals)
        self.slider.valueChanged.connect(self.update_hwhm_from_slider)
        
        self.winc_slider = QSlider(Qt.Horizontal)
        self.winc_slider.setMaximumWidth(200)
        self.winc_slider.setValue(30)
        self.winc_slider.setTracking(True)
        self.winc_slider.setTickPosition(QSlider.TicksBothSides)
        # PyQt4 (old-style signals)
        #self.connect(self.winc_slider, SIGNAL('valueChanged(int)'), self.update_hwhm_from_slider)
        # PyQt5 (new-style signals)
        self.winc_slider.valueChanged.connect(self.update_incident_freq_from_slider)
        
        self.fixaxes_cb = QCheckBox("Fix y-axis")
        self.fixaxes_cb.setChecked(False)
        self.fixaxes_cb.setMaximumWidth(100)

        self.fixlegend_cb = QCheckBox("Fix legend")
        self.fixlegend_cb.setChecked(True)
        self.fixlegend_cb.setMaximumWidth(100)
        # PyQt4 (old-style signals)
        #self.connect(self.fixlegend_cb, SIGNAL('stateChanged(int)'), self.update_fixlegend)
        # PyQt5 (new-style signals)
        self.fixlegend_cb.stateChanged.connect(self.update_fixlegend)
        
        self.showRayleigh_cb = QCheckBox("Show Rayleigh")
        self.showRayleigh_cb.setChecked(False)
        self.showRayleigh_cb.setMaximumWidth(100)
        # Update when clicking
        # PyQt4 (old-style signals)
        #self.connect(self.show_Rayleigh, SIGNAL('stateChanged(int)'), self.update_convolute)
        # PyQt5 (new-style signals)
        self.showRayleigh_cb.stateChanged.connect(self.update_incident_freq_from_slider)
        
        # Labels
        hwhm_label  = QLabel('   HWHM')
        cm_label    = QLabel('(cm-1)')
        cm_label2   = QLabel('(cm-1)')
        broad_label = QLabel('Broadening')
        datatype_label = QLabel('Data Type')
        scalefactor_label = QLabel('Scale')
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
        self.analysis_box.setFontFamily('Courier') #use a monospaced typeface to ensure good aligment of the text
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
        # PyQt4
        #title.setBackgroundColor(Qt.lightGray)
        # PyQt5
        title.setBackground(Qt.lightGray)
        # PyQt4 (function deprecated in PyQt5)
        #title.setTextColor(Qt.black)
        # PyQt5
        title.setForeground(Qt.black)
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
            # PyQt4
            #cell.setBackgroundColor(Qt.gray)
            # PyQt5
            cell.setBackground(Qt.gray)
            # PyQt4
            #cell.setTextColor(Qt.white)
            # PyQt5
            cell.setForeground(Qt.white)
            cell.setFont(font)
            # Set non editable. See: http://stackoverflow.com/questions/2574115/how-to-make-a-column-in-qtablewidget-read-only
            #cell.setFlags(cell.flags() ^ Qt.ItemIsEditable)
            cell.setFlags(cell.flags() ^ QtCore.Qt.ItemIsEnabled ^ QtCore.Qt.ItemIsEditable)
        ## Second column
        celllabel = ["No","-","-"]
        for i,j in cellids[1]:
            self.refspc_table.setItem(i,j, QTableWidgetItem(celllabel[i-1]))
            cell = self.refspc_table.item(i,j)
            # PyQt4
            #cell.setTextColor(Qt.black)
            # PyQt5
            cell.setForeground(Qt.black)
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
        # PyQt4
        #xbutton.setBackgroundColor(Qt.red)
        # PyQt5
        xbutton.setBackground(Qt.red)
        ## Tune the Tare buttons
        tbutton = self.refspc_table.item(2,2)
        # PyQt4
        #tbutton.setBackgroundColor(Qt.blue)
        # PyQt5
        tbutton.setBackground(Qt.blue)
        # PyQt4
        #tbutton.setTextColor(Qt.white)
        # PyQt5
        tbutton.setForeground(Qt.white)
        tbutton = self.refspc_table.item(3,2)
        # PyQt4
        #tbutton.setBackgroundColor(Qt.blue)
        # PyQt5
        tbutton.setBackground(Qt.blue)
        # PyQt4
        #tbutton.setTextColor(Qt.white)
        # PyQt5
        tbutton.setForeground(Qt.white)
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
        hbox_slider.addWidget(cm_label)
        # WInc slider with textbox
        hbox_wslider = QHBoxLayout()
        hbox_wslider.addWidget(self.winc_slider)
        hbox_wslider.addWidget(self.winc_box)
        hbox_wslider.addWidget(cm_label2)
        # HWHM label with inputBins checkbox
        hbox_hwhmlab = QHBoxLayout()
        hbox_hwhmlab.addWidget(hwhm_label)
        hbox_hwhmlab.setAlignment(hwhm_label, Qt.AlignLeft)
        hbox_hwhmlab.addWidget(self.showRayleigh_cb)
        hbox_hwhmlab.setAlignment(self.showRayleigh_cb, Qt.AlignLeft)
        # Complete HWHM widget merging all box here
        vbox_slider = QVBoxLayout()
        vbox_slider.addLayout(hbox_hwhmlab)
        vbox_slider.setAlignment(hbox_hwhmlab, Qt.AlignLeft)
        vbox_slider.addLayout(hbox_slider)
        vbox_slider.setAlignment(hbox_slider, Qt.AlignLeft)
        vbox_slider.addLayout(hbox_wslider)
        vbox_slider.setAlignment(hbox_wslider, Qt.AlignLeft)
        # Clean button
        vbox_cleaner = QVBoxLayout()
        vbox_cleaner.addWidget(clean_button1)
        vbox_cleaner.addWidget(clean_button2)
        # DataType sector
        vbox_scalefactor = QVBoxLayout()
        vbox_scalefactor.addWidget(scalefactor_label)
        vbox_scalefactor.setAlignment(scalefactor_label, Qt.AlignLeft)
        vbox_scalefactor.addWidget(self.scale_factor_box)
        vbox_scalefactor.setAlignment(self.scale_factor_box, Qt.AlignLeft)
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
        hbox.addLayout(vbox_scalefactor)
        hbox.setAlignment(vbox_scalefactor, Qt.AlignTop)
        
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
        dat_export_action = self.create_action("&Export to dat", slot=self.dat_export, 
            shortcut="Ctrl+D", tip="Export to text file (dat)")
        spc_import_action = self.create_action("&Import plot", slot=self.open_plot, 
            shortcut="Ctrl+N", tip="Import spectrum")
        # Now place the actions in the menu
        self.add_actions(self.file_menu, 
            (load_file_action, xmgr_export_action, dat_export_action, spc_import_action, None, quit_action))
        
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
                

    # PyQt4
    #def create_action(  self, text, slot=None, shortcut=None, 
    #                    icon=None, tip=None, checkable=False, 
    #                    signal="triggered()"):
    # PyQt5
    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            # PyQt4 (old-style signals)
            #self.connect(action, SIGNAL(signal), slot)
            # PyQt5 (new-style signals)
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action
    
    
#==========================================================
# GLOBAL FUNCTIONS
#==========================================================
# FILE READERS
def read_spc_xy(filename,fromsection=None):
    """
    Function to read fort.22, which contains the bins to reconstruct
    the convoluted spectrum_sim.
    This is a simple [x,y] file
    An alternative would be to use the np.loadtxt. But for the 
    moment we leave it like that
    """
    # Open and read file
    print("Loading spectral data from '"+filename+"'...")
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
def convolute(spc_stick,res=1,hwhm=0.1,broad="Gau",input_bins=False):
    """
    Make a Gaussian convolution of the stick spectrum
    The spectrum must be in energy(eV) vs Intens (LS?)
    
    Arguments:
    spc_stick  list of list  stick spectrum as [x,y]
               list of array
    res        float         resolution (for the final graph)
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
    n_extensions = 0
    sigma = hwhm / np.sqrt(2.*np.log(2.))
    while not recovered_baseline:
        n_extensions += 1
        if input_bins:
            # Backup npoints
            npts = npoints
            npoints = len(x)
            xhisto = x
            yhisto = y
            width = (x[1] - x[0])
        else:
            extra_x = (x[-1] - x[0])*extra_factor
            npoints = int((x[-1]+extra_x-x[0]-extra_x)/res)
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
        ymax = max(abs(yconv.max()),abs(yconv.min()))
        if abs(yconv[0]) < ymax/100.0 and abs(yconv[-1]) < ymax/100.0:
            recovered_baseline=True
        # If exteded 100 times exit anyway
        if n_extensions > 100:
            print('WARNING in convolution: Baseline cannot be recovered!')
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
    
    print("# XMGRACE CREATED BY FCC_ANALYZER", file=f)
    print("# Only data and labels. Format will", file=f)
    print("# be added by your default xmgrace", file=f)
    print("# defaults (including colors, fonts...)", file=f)
    print("# Except the followins color scheme:", file=f)
    print('@map color 0  to (255, 255, 255), "white"', file=f)
    print('@map color 1  to (0, 0, 0), "black"', file=f)
    print('@map color 2  to (0, 0, 255), "blue"', file=f)
    print('@map color 3  to (255, 0, 0), "red"', file=f)
    print('@map color 4  to (0, 139, 0), "green4"', file=f)
    print('@map color 5  to (0, 255, 255), "cyan"', file=f)
    print('@map color 6  to (255, 0, 255), "magenta"', file=f)
    print('@map color 7  to (188, 143, 143), "brown"', file=f)
    print('@map color 8  to (100, 0, 100), "pink"', file=f)
    print('@map color 9  to (255, 165, 0), "orange"', file=f)
    print('@map color 10 to (255, 255, 0), "yellow"', file=f)
    print('@map color 11 to (220, 220, 220), "grey"', file=f)
    print('@map color 12 to (0, 255, 0), "green"', file=f)
    print('@map color 13 to (148, 0, 211), "violet"', file=f)
    print('@map color 14 to (114, 33, 188), "indigo"', file=f)
    print('@map color 15 to (103, 7, 72), "maroon"', file=f)
    print('@map color 16 to (64, 224, 208), "turquoise"', file=f)
    print('@map color 17 to (50, 50, 50), "gris2"', file=f)
    print('@map color 18 to (100, 100, 100), "gris3"', file=f)
    print('@map color 19 to (150, 150, 150), "gris4"', file=f)
    print('@map color 20 to (200, 200, 200), "gris5"', file=f)
    print('@map color 21 to (255, 150, 150), "red2"', file=f)
    print('@map color 22 to (150, 255, 150), "green2"', file=f)
    print('@map color 23 to (150, 150, 255), "blue2"', file=f)  
    # Without the @version, it makes auto-zoom (instead of taking world coords) 
    print("@version 50123", file=f)
    print("@page size 792, 612", file=f)
    print("@default symbol size 0.010000", file=f)
    print("@default char size 0.800000", file=f)
    for lab in labs:
        print("@with line", file=f)
        print("@    line on", file=f)
        print("@    line g0", file=f)
        print("@    line loctype world", file=f)
        print("@    line color 20", file=f)
        print("@    line ",lab.xy[0],",",lab.xy[1],",",lab.xyann[0],",",lab.xyann[1], file=f)
        print("@line def", file=f)
        print("@with string", file=f)
        print("@    string on", file=f)
        print("@    string g0", file=f)
        print("@    string loctype world", file=f)
        print("@    string ", lab.xyann[0],",",lab.xyann[1], file=f)
        print("@    string def \"",labs[lab],"\"", file=f)
    print("@with g0", file=f)
    # Set a large view
    print("@    view 0.180000, 0.150000, 1.15, 0.92", file=f)
    #Get plotting range from mplt
    x=ax.get_xbound()
    y=ax.get_ybound()
    print("@    world ",x[0],",",y[0],",",x[1],",",y[1], file=f)
    #Get xlabel from mplt
    print("@    xaxis  label \""+ax.get_xlabel()+"\"", file=f)
    print("@    yaxis  label \""+ax.get_ylabel()+"\"", file=f)
    if ax2:
        # Position
        print("@    yaxis  label place opposite", file=f)
        print("@    yaxis  ticklabel place opposite", file=f)
        print("@    yaxis  tick place opposite", file=f)
    # Char sizes
    print("@    xaxis  ticklabel char size 1.250000", file=f)
    print("@    yaxis  ticklabel char size 1.250000", file=f)
    print("@    xaxis  label char size 1.500000", file=f)
    print("@    yaxis  label char size 1.500000", file=f)
    #Get tick spacing from mplt
    x=ax.get_xticks()
    y=ax.get_yticks()
    print("@    xaxis  tick major", x[1]-x[0], file=f)
    print("@    yaxis  tick major", y[1]-y[0], file=f)
    #Legend
    print("@    legend char size 1.250000", file=f)
    print("@    legend loctype view", file=f)
    print("@    legend 0.95, 0.9", file=f)
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
        print("& %s"%(label_list[iclass]), file=f)
        print("@type bar", file=f)
        print("@    s"+str(k)," line type 0", file=f)
        print("@    s"+str(k)," legend  \"%s\""%(label_list[iclass]), file=f)
        print("@    s"+str(k)," symbol color %s"%(color_list[iclass]), file=f)
        for i in range(len(x)):
            print(x[i], y[i], file=f)
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
        print("Scale factor applied to convoluted spectrum %s"%(scale_factor))
    else:
        k = 0
        print("@g1 on", file=f)
        print("@g1 hidden false", file=f)
        print("@with g1", file=f)
        # Set a large view
        print("@    view 0.180000, 0.150000, 1.15, 0.92", file=f)
        #Get plotting range from mplt
        x=ax2.get_xbound()
        y=ax2.get_ybound()
        print("@    world ",x[0],",",y[0],",",x[1],",",y[1], file=f)
        #Get labels from mplt
        print("@    yaxis  label \""+latex2xmgrace(ax2.get_ylabel())+"\"", file=f)
        # Char sizes
        print("@    xaxis  ticklabel char size 1.250000", file=f)
        print("@    yaxis  ticklabel char size 1.250000", file=f)
        print("@    xaxis  label char size 1.500000", file=f)
        print("@    yaxis  label char size 1.500000", file=f)
        #Set axis ticks on the left
        print("@    xaxis  off", file=f)
        print("@    yaxis  tick out", file=f)
        print("@    yaxis  label place normal", file=f)
        print("@    yaxis  ticklabel place normal", file=f)
        print("@    yaxis  tick place normal", file=f)
        #Get tick spacing from mplt
        y=ax2.get_yticks()
        print("@    yaxis  tick major", y[1]-y[0], file=f)
        #Legend
        print("@    legend char size 1.250000", file=f)
        print("@    legend loctype view", file=f)
        print("@    legend 0.8, 0.9", file=f)
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
        print("& "+leg_label, file=f)
        print("@type xy", file=f)
        print("@    s"+str(k)," line type 1", file=f)
        print("@    s"+str(k)," line linestyle %s"%(style_list[i]), file=f)
        print("@    s"+str(k)," line color %s"%(color_list[i]), file=f)
        print("@    s"+str(k)," legend  \""+leg_label+"\"", file=f)
        for j in range(len(x)):
            print(x[j], y[j]*scale_factor, file=f)
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
    
    for key,value in input_args_dict.items():
        # Check it is allowed
        isValid = final_arguments.get(key,None)
        if isValid is None:
            raise BaseException("Sintax error. Unknown label: " + key)
        # If valid, update final argument
        final_arguments[key]=value
        
    if final_arguments.get("-h"):
        
        print("""
 ----------------------------------------
           FCclasses analyzer
   A GUI to analyze FCclasses output 
 ----------------------------------------
 Version info 
        Git commit: %s
        Date: %s
        """%(version_tag.COMMIT,version_tag.DATE))
        helptext()
        print("    Options:")
        print("    --------")
        print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format("Flag","Type","Description","Value"))
        print('      {0:-<10}  {1:-^4}  {2:-<41}  {3:-<7}'.format("","","",""))
        for key,value in final_arguments.items():
            descr = arg_description[key]
            atype = arg_type[key]
            #atype=str(type(value)).replace("<type '","").replace("'>","")
            print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format(key, atype, descr, str(value)))
        print("")
        
        sys.exit()
        
    return final_arguments
        
    
    
def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()


