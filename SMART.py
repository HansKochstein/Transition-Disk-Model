""""""""""""""" PROBLEM SETUP """""""""""""""
# This code creates the protoplanetary disk structures and outputs as a set of files, 
# which will be input for later radiative transfer calculations by RADMC-3D.


import tkinter as tk
from tkinter import filedialog
import numpy as np
import os
import tools
import threading
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# ----------------------------------
#   Physical constants in cgs unit
# ----------------------------------

# Define useful constants.  
# Note that cgs units are commonly used in astronomy.
# Length is measured in [cm], mass is measured in [g], and time is measured in [s].


au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
ss  = 5.6703e-5      # Stefan-Boltzmann const  [erg/cm^2/K^4/s]
kb  = 1.3807e-16     # Bolzmann's constant     [erg/K]
mp  = 1.6726e-24     # Mass of proton          [g]
GG  = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
pi  = np.pi          # Pi


''''''''''''' INPUT PARAMETER '''''''''''''''

# OFFSET = False

nr       = 120                       # Nr of radial grid points
nphi     = 100                       # Nr of azimuthual grid points 
ntheta   = 100                       # Nr of polar grid points

rin      = 1*au                     # Inner radius
rout     = 120*au                   # Outer radius

#INNER DISK
sigmad0   = 0.1                     # Sigma dust (in the middle of the Ring)  
ringcau1  = 7                       # R_{ring, center}: the radial location of the ring center in [au]
ringwau1  = 1                       # R_{ring, width}: the radial width of the ring in [au]
elevD     = 20                      # inner DISK elevation in DEGREES
rotaD1    = 210                     # inner DISK rotation in DEGREES (azimuthal rotation)

#TESTING
INNERASSYMETRIC = False
azipos1   = pi                      # inner DISK rotation of cresendent structur in RAD
azimwid1  = 0.45                    # inner DISK azimuthual width of cresendent structure in RAD

#OUTER DISK
sigmadc2  = 10                      # Sigma_{dust, center}: dust surface density at the center of the ring [g/cm^2]
ringcau2  = 88                      # R_{ring, center}: the radial location of the ring center in [au]
ringwau2  = 9                       # R_{ring, width}: the radial width of the ring in [au]
rotaD2    = 95                      # outer DISK rotation in (azimuthal rotation) & outer DISK rotation of cresendent structur in DEG

OUTERASSYMETRIC = True
azimwid2  = 0.75                   # outer DISK azimuthual width of cresendent structure in RAD


ELLIPSE = True
a_semi  = ringcau2 * au            # [cm] semi-major axis
e       = 0.2                      # eccentricity

nphot    = 1000000                 # Photones for Monte Carlo Simulation  



''' Star parameters '''

mstar    = ms * 2.2
rstar    = rs * 1.8
tstar    = ts + 3700  # K 
pstar    = np.array([0.,0.,0.])


class DiskModelGUI:
    def __init__(self, master):
        self.master = master
        master.title("Transitdisk Disk Model")
   

        #Create output frame and text widget
        # Create main frames for left and right sides
        left_frame = tk.Frame(master)
        left_frame.grid(row=0, column=0, sticky="nsew")
        right_frame = tk.Frame(master)
        right_frame.grid(row=0, column=1, sticky="nsew")

        # Configure the grid to expand
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure(0, weight=1)
        master.grid_columnconfigure(1, weight=1)
        
        # Create input frames within left_frame
        self.create_grid_frame(left_frame)
        self.create_inner_disk_frame(left_frame)
        self.create_outer_disk_frame(left_frame)
        self.create_plot_para_frame(left_frame)
        self.running = False  # Flag to track if simulation is running

        # # Create output elements within right_frame
        # self.output_text = tk.Text(right_frame, wrap=tk.WORD)  
        # self.output_text.grid(row=0, column=0, sticky="nsew")

        self.loading_label = tk.Label(left_frame, text="Please Input Parameters", font=("Courier", 12))
        self.loading_label.grid(row=5, column=0, sticky="nsew")

        # Store output_frame as an attribute
        self.output_frame = right_frame
        self.input_frame  = left_frame

        # Create placeholder frame for the plot
        self.plot_frame = tk.LabelFrame(self.output_frame, text='PLOT')
        self.plot_frame.grid(row=0, column=0, rowspan=4, columnspan=6,padx=5, pady=5, sticky="nsew")
        self.plot_canvases = []

        # Create a placeholder frame for the plot
        fig = plt.Figure(figsize=(8, 6))
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.plot_canvases.append(canvas)


        # Create abort/close button
        abort_button = tk.Button(self.input_frame, text="Close", command=self.abort_and_close)
        abort_button.grid(row=6, column=0,columnspan=2, pady=10, sticky="e")  # Place next to loading label

        # Create button to run simulation
        button = tk.Button(self.input_frame, text="Run Simulation", command=self.update_mode)
        button.grid(row=6, column=0, columnspan=2, pady=10, sticky='w')

        # Create a Checkbutton to toggle between scattering and thermal re-emission of starlight by dust grains
        self.scattering_var = tk.IntVar(value=1)  # This will hold the state of the checkbox
        scattering_checkbutton = tk.Checkbutton(self.input_frame, text="Scattering(on)/thermal(off)", variable=self.scattering_var)  # Create the checkbox
        scattering_checkbutton.grid(row=0, column=0, columnspan=2, pady=10, sticky="w")  # Place the checkbox in the input frame


        # Create a Checkbutton to toggle mask on/off
        self.mask_var = tk.IntVar(value=1)  # This will hold the state of the checkbox
        mask_checkbutton = tk.Checkbutton(self.output_frame, text="Apply Mask", variable=self.mask_var)
        mask_checkbutton.grid(row=7, column=0, columnspan=1, pady=10, sticky="n")

        # Create a Checkbutton to toggle between single plot and grid of plots
        self.single_plot_var = tk.IntVar(value=1)
        single_plot_checkbutton = tk.Checkbutton(self.output_frame, text="Single Plot/ 3x3 Grid", variable=self.single_plot_var, command=self.update_plots)
        single_plot_checkbutton.grid(row=7, column=2, columnspan=1, pady=10, sticky="n")

         # Add a Save button
        save_button = tk.Button(self.output_frame, text="Save Image", command=self.save_plot)
        save_button.grid(row=7, column=4, columnspan=1, pady=10, sticky="n")  # Place at the bottom of the plot frame


    def update_plots(self):
        # Clear the plot frame
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        self.plot_canvases = []

        if self.single_plot_var.get() == True:
            # Create a single plot
            fig = plt.Figure(figsize=(8, 6))
            canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            self.plot_canvases.append(canvas)
        else:
            # Create a grid of plots
            for i in range(9):
                plot_frame = tk.Frame(self.plot_frame)
                plot_frame.grid(row=i // 3, column=i % 3, padx=2, pady=2, sticky="nsew")

                fig = plt.Figure(figsize=(3.2, 2.6))
                canvas = FigureCanvasTkAgg(fig, master=plot_frame)
                canvas.draw()
                canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
                self.plot_canvases.append(canvas)

    '''Create Frame for input parameters'''

    # Create a frame for grid parameters
    def create_grid_frame(self,parent_frame):
        grid_frame = tk.LabelFrame(parent_frame, text="Grid Parameters")
        grid_frame.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        # Input labels and entries
        labels = ["Radial grid number", "Azimuthal grid number", "Polar grid number", "Grid inner radius [au]", "Grid outer radius [au]", "Photon number"]
        variable_names = ["nr", "nphi", "ntheta", "rin", "rout", "nphot"]
        defaults = [40, 40, 40, 1, 120, 100000]

        for i, label in enumerate(labels):
            tk.Label(grid_frame, text=label).grid(row=i, column=0, sticky="w")
            entry = tk.Entry(grid_frame)
            entry.insert(0, str(defaults[i]))
            entry.grid(row=i, column=1, padx=5)
            setattr(self, variable_names[i], entry)

    # Create a frame for inner disk parameters
    def create_inner_disk_frame(self,parent_frame):
        inner_disk_frame = tk.LabelFrame(parent_frame, text="Inner Disk Parameters")
        inner_disk_frame.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")

        # Input labels and entries
        labels = ["Peak Surface Denisty [g/cm^2]", "Ring distance [au]", "Radial width [au]","Inclination [deg]","Rotation [deg]"]
        variable_names = ["sigmad0", "ringcau1", "ringwau1","elevD","rotaD1",]
        defaults = [0.1, 7, 1, 20, 210]

        for i, label in enumerate(labels):
            tk.Label(inner_disk_frame, text=label).grid(row=i, column=0, sticky="w")
            entry = tk.Entry(inner_disk_frame)
            entry.insert(0, str(defaults[i]))
            entry.grid(row=i, column=1, padx=5)
            setattr(self, variable_names[i], entry)  # Store the Entry widget

    # Create a frame for outer disk parameters        
    def create_outer_disk_frame(self,parent_frame):
        outer_disk_frame = tk.LabelFrame(parent_frame, text="Outer Disk Parameters")
        outer_disk_frame.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")

        # Input labels and entries (example parameters)
        labels = ["Peak Surface Density [g/cm^2]","Ring distance [au]","Radial width [au]","Rotation [deg]","Asymmetry","Azimuthal width","Ellipse","Semi-major axis","Eccentricity"]
        variable_names = ["sigmadc2","ringcau2","ringwau2","rotaD2","OUTERASSYMETRIC","azimwid2","ELLIPSE","a_semi","e"]
        defaults = [10, 88, 9, 95, True, 0.75, True, 88, 0.2]

        for i, label in enumerate(labels):
            tk.Label(outer_disk_frame, text=label).grid(row=i, column=0, sticky="w")

            if label in ["Asymmetry", "Ellipse"]:
                # Create a checkbox for boolean parameter
                var = tk.BooleanVar(value=defaults[i])
                checkbutton = tk.Checkbutton(outer_disk_frame, variable=var)
                checkbutton.grid(row=i, column=1, padx=5)
                setattr(self, variable_names[i], var)  # Store the BooleanVar 
            else:
                # Create an entry for other parameters
                entry = tk.Entry(outer_disk_frame)
                entry.insert(0, str(defaults[i]))
                entry.grid(row=i, column=1, padx=5)
                setattr(self, variable_names[i], entry)  # Store the Entry widget
    
    # Function to create the plot parameter frame
    def create_plot_para_frame(self, parent_frame):
        plot_para = tk.LabelFrame(parent_frame, text="Image Parameter")
        plot_para.grid(row=4, column=0, padx=5, pady=5, sticky="nsew")
        labels = ['Wavelength [um]', 'Inclination [deg]', 'Position Angle [deg]', 'Index [0-8]']
        variable_names = ['wavelength', 'incl_value', 'posang_value','index']
        defaults = [2.2, 50, 7.4, 0]
        for i, label in enumerate(labels):
            tk.Label(plot_para, text=label).grid(row=i, column=0, sticky="w")
            entry = tk.Entry(plot_para)
            entry.insert(0, str(defaults[i]))
            entry.grid(row=i, column=1, padx=5)
            setattr(self, variable_names[i], entry)  # Store the Entry widget


    '''Running the simulation'''

    # Function deciding for sacttering or thermal re-emission of starlight by dust grains
    def update_mode(self):
        if self.scattering_var.get() == 1:
            self.run_simulation_scattering()
        else:
            self.run_simulation_thermal()

    # Funtction to run the scattering simulation
    def run_simulation_scattering(self):
        # Get values from input fields
        nr        = int(self.nr.get())
        nphi      = int(round(float(self.nphi.get())))
        ntheta    = int(round(float(self.ntheta.get())))
        rin       = int(round(float(self.rin.get()))) * au
        rout      = int(round(float(self.rout.get()))) * au 
        #INNER DISK
        sigmad0   = float(self.sigmad0.get())
        ringcau1  = float(self.ringcau1.get())
        elevD     = float(self.elevD.get())
        rotaD1    = float(self.rotaD1.get())

        #TESTING
        INNERASSYMETRIC = False 
        azipos1   = pi
        azimwid1  = 0.45

        #OUTER DISK
        sigmadc2  = float(self.sigmadc2.get())
        ringcau2  = float(self.ringcau2.get())
        ringwau2  = float(self.ringwau2.get())
        rotaD2    = float(self.rotaD2.get())

        OUTERASSYMETRIC = bool(self.OUTERASSYMETRIC.get())
        azimwid2  = float(self.azimwid2.get())


        ELLIPSE = bool(self.ELLIPSE.get())
        a_semi  = float(self.a_semi.get())  * au            # [cm] semi-major axis
        e       = float(self.e.get())                     # eccentricity
        
        nphot    = float(self.nphot.get())                  # Photones for Monte Carlo Simulation  

        wavelength  = float(self.wavelength.get())
        incl_value = float(self.incl_value.get())
        posang_value = float(self.posang_value.get())
        plot_index = int(self.index.get())

        ''' Make the r coordinates'''


        nlev_rin = 20                                       # Grid refinement at the inner edge: nr of cycles
        nspan_rin= 5                                        # Grid refinement at the inner edge: nr of cells each cycle
        ri       = np.logspace(np.log10(rin),np.log10(rout),nr+1)
        ri       = tools.grid_refine_inner_edge(ri,nlev_rin,nspan_rin)   # Refinement at inner edge
        rc       = 0.5 * ( ri[:-1] + ri[1:] )
        nr       = len(rc)           # Recompute nr, because of refinement at inner edge


        ''' Make the theta and phi coordinates '''


        thetaup  = 0.1
        thetai   = np.linspace(thetaup,np.pi-thetaup,ntheta+1)
        phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
        # thetac   = 0.5 * ( thetai[0:ntheta] + thetai[1:ntheta+1] )
        # phic     = 0.5 * ( phii[0:nphi] + phii[1:nphi+1] )
        thetac   = 0.5 * ( thetai[:-1] + thetai[1:] )
        phic     = 0.5 * ( phii[:-1] + phii[1:] )


        ''' Make the 3-D-grid '''


        qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
        r        = qq[0]   # The radial grid of the analytic disk model
        tt       = qq[1]
        pp       = qq[2]

        # ttdeg    = np.rad2deg(tt)
        # ppdeg    = np.rad2deg(pp)


        "Calculating each component of the coordinates (r, phi, theta) of the rotated frame that belongs to the inner disk"


        elevR    = np.deg2rad(elevD)        # Elevation of the disk in RAD
        newr     = r
        coselev  = (-np.sin(tt)*np.sin(pp)*np.sin(elevR))+(np.cos(tt)*np.cos(elevR))

        newtt       = np.arccos(coselev)
        sinnewpp    = np.sin(tt)*np.sin(pp)*np.cos(elevR)+np.cos(tt)*np.sin(elevR)/(np.sqrt(1-coselev**2))
        sinnewtt    = np.sin(newtt)
        cosnewpp    = (np.sin(tt)*np.cos(pp))/sinnewtt
        newpp       = np.arctan2(sinnewpp,cosnewpp)
        newpp       = np.where(newpp < 0, newpp + 2 * pi, newpp) #shifting atan2 domain from [-pi,pi] towards [0,2pi]


        '''Rearranging the assymmetrical shape'''


        dphi  =  newpp - azipos1

        # dphi = np.where(dphi >= pi, dphi + 2*pi, dphi)
        # dphi = np.where(dphi < 0, dphi + 2*pi, dphi)
        # dphi[dphi < 0] += 2 * pi 
        # if ((dphi>pi).all) or ((dphi<-pi).all):
        #      dphi += 2*pi

        ''' Make the dust density model '''

        #INNER DISK
        q_inner  = 0.25
        flang    = 0.1                         # The assumed constant radiative incidence angle
        lstar    = 4*pi*rstar**2*ss*tstar**4   # Stellar luminosity L = sigma*A*T^4
        firr     = flang*lstar/(4*pi*r**2)     # Irradiative flux F = L/A with A=4pi*r^2
        tmid     = (firr/ss)**q_inner          # Estimate of midplane temperature (^0.25)
        # ymid     = tmid_0 * (r/)
        cs1      = np.sqrt(kb*tmid/(2.3*mp))   # Isothermal sound speed at midplane
        omk      = np.sqrt(GG*mstar/r**3)      # The Kepler angular frequency
        hp1      = cs1/omk                      # The pressure scale height
        hpr1     = hp1/r                        # The dimensionless hp

        #OUTER DISK
        q_outer  = 0.25
        flang    = 0.1                         # The assumed constant radiative incidence angle
        lstar    = 4*pi*rstar**2*ss*tstar**4   # Stellar luminosity L = sigma*A*T^4
        firr     = flang*lstar/(4*pi*r**2)     # Irradiative flux F = L/A with A=4pi*r^2
        tmid2    = (firr/ss)**q_outer          # Estimate of midplane temperature (^0.25)

        cs2      = np.sqrt(kb*tmid2/(2.3*mp))   # Isothermal sound speed at midplane
        omk      = np.sqrt(GG*mstar/r**3)      # The Kepler angular frequency
        hp2      = cs2/omk                      # The pressure scale height
        hpr2     = hp2/r                        # The dimensionless hp

        ''' INNER DISK '''

        ringc1    = ringcau1*au                # Convert [au] to [cm]
        ringw1    = ringwau1*au                # Convert [au] to [cm]


        if INNERASSYMETRIC==True:
            sigmad    = sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0) * np.exp(-(((dphi)**2)/azimwid1**2)/2.0) # +  sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0) * np.exp(-(((newpp - 1.8*pi)**2)/azimwid1**2)/2.0)     # radial & azimuthal surface density profile for ring [g/cm^2]
        elif INNERASSYMETRIC==False:
            sigmad    = sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0)                                                  # radial surface density profile for ring [g/cm^2]
        
        newzr     = pi/2.e0 - newtt 

        rhod     = ( sigmad / (np.sqrt(2.e0*np.pi)*hp1) ) * np.exp(-(newzr**2/hpr1**2)/2.e0)                                   # volume density

        rotation1 = int(rotaD1/(360/nphi))          
        rhod     = np.roll(rhod, rotation1, axis=2)


        ''' OUTER DISK '''


        ringc2 = ringcau2*au                      # Radius of the ring, Convert [au] to [cm] 

        if ELLIPSE == True:
            ringc2_mod = a_semi*(1-e**2)/(1-e*np.cos(pp))
        else:
            ringc2_mod = ringc2

        ringw2  = ringwau2*au                   # Radial width of the ring, Convert [au] to [cm]

        if OUTERASSYMETRIC==True:
            sigmad2   = sigmadc2  * np.exp(-((newr-ringc2_mod)**2/ringw2**2)/2.0) * np.exp(-((pp-1.05*pi)**2/azimwid2**2)/2.0) # radial & azimuthal surface density profile for ring [g/cm^2]
        else:
            sigmad2   = sigmadc2  * np.exp(-((newr-ringc2_mod)**2/ringw2**2)/2.0)                                     # radial surface density profile for ring [g/cm^2]

        zr        = np.pi/2.e0 - tt                                                      


        rhod2     = ( sigmad2 / (np.sqrt(2.e0*np.pi)*hp2) ) * np.exp(-(zr**2/hpr2**2)/2.e0)                                   # volume density

        rotation2 = int(rotaD2/(360/nphi))                                                                                  # only if assymetrical shape 
        rhod2     = np.roll(rhod2, rotation2, axis=2)

    
        ''' ADDING OUTER & INNER DISK '''

        rhod     += rhod2

        #
        # Write the grid file
        #
        with open('amr_grid.inp','w+') as f:
            f.write('1\n')                       # iformat
            f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
            f.write('100\n')                     # Coordinate system: spherical
            f.write('0\n')                       # gridinfo
            f.write('1 1 1\n')                   # Include r,theta coordinates
            f.write('%d %d %d\n'%(nr,ntheta,nphi))  # Size of grid
            for value in ri:
                f.write('%21.14e\n'%(value))      # X coordinates (cell walls)
            for value in thetai:
                f.write('%21.14e\n'%(value))      # Y coordinates (cell walls)
            for value in phii:
                f.write('%21.14e\n'%(value))      # Z coordinates (cell walls)
                

        #
        # Write the density file
        #
        with open('dust_density.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
            f.write('1\n')                       # Nr of dust species
            data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%21.14e")
            f.write('\n')
    
        # ----------------------------------
        #
        # Wavelength [um] data and create wavelength_micron.inp file
        #
        # ----------------------------------

        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size

        #
        # Write the wavelength file
        #
        with open('wavelength_micron.inp','w+') as f:
            f.write('%d\n'%(nlam))
            for value in lam:
                f.write('%21.14e\n'%(value))


        #
        # Write the stars.inp file
        #
        with open('stars.inp','w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n'%(nlam))
            f.write('%21.14e %21.14e %21.14e %21.14e %21.14e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
            for value in lam:
                f.write('%21.14e\n'%(value))
            f.write('\n%21.14e\n'%(-tstar))
            

        #
        # Dust opacity control file
        #
        with open('dustopac.inp','w+') as f:
            f.write('2               Format number of this file\n')
            f.write('1               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('10              Way in which this dust species is read\n')#[MODIFIED] opacity with scattering information
            f.write('0               0=Thermal grain\n')
            f.write('pyrmg70         Extension of name of dustkappa_***.inp file\n')#[MODIFIED] opacity with scattering information
            f.write('----------------------------------------------------------------------------\n')
        #
        # Write the radmc3d.inp control file
        #
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 5\n') #[MODIFIED] radiative transfer with scattring
            f.write('iranfreqmode = 1\n')

        self.loading_label.config(text="Simulation running... 💫, don't interrupt")
        self.master.update()

        self.running = True
        thread = threading.Thread(target=self.run_radmc3d)
        thread.start()

        # Getting input parameters for imaging data 

        cmd = f'radmc3d image lambda {wavelength} npix 200 sizeau 258 incl {incl_value} posang {posang_value} nphot_scat 1000000 stokes setthreads 8'
        os.system(cmd)

        # incl_value = cmd.split("incl ")[1].split()[0]
        # posang_value = cmd.split("posang ")[1].split()[0]

        fnameread      = 'image.out'

        with open(fnameread,'r') as f:
            s        = f.readline()
            im_n      = f.readline()
            im_nx     = int(im_n.split()[0])
            im_ny     = int(im_n.split()[1])
            nlam      = int(f.readline())
            pixsize   = f.readline()
            pixsize_x = float(pixsize.split()[0])
            pixsize_y = float(pixsize.split()[1])
            wavelength_mic = np.zeros(nlam)
            im    = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
            rr_au = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
            phi   = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
            for ilam in range(nlam): wavelength_mic[ilam] = f.readline()
            for ilam in range(nlam):
                dummy = f.readline() #blank line
                for iy in range(im_ny):
                    for ix in range(im_nx):
                        image_temp = f.readline()
                        im[ix,iy,:,ilam] = [image_temp.split()[istokes] for istokes in [0,1,2,3]] #[MODIFIED] array dimension modified for Stokes I,Q,U,V
                        rr_au[ix,iy,ilam] = np.sqrt(((ix-im_nx/2)*pixsize_x)**2+((iy-im_ny/2)*pixsize_y)**2)/au
                        phi[ix,iy,ilam]   = np.arctan2((iy-im_ny/2)*pixsize_y,(ix-im_nx/2)*pixsize_x) 
        sizex_au = im_nx*pixsize_x/au
        sizey_au = im_ny*pixsize_y/au

        #polarized intensity is stored in the array im_PI
        im_PI = np.zeros((im_nx,im_ny,nlam,4))
        im_PI = ((im[:,:,1,0])**2.0 + (im[:,:,2,0])**2.0)**0.5
        d_pc = 136   # For IRS48 first measured 121 pc, which got updated to 136pc
        Bmaj_as        = 0.07 # Major axis of the beam in arcsec
        Bmin_as        = 0.07 # Minor axis of the beam in arcsec
        # wavelength2     = 0.00022  # wavelength in unit of cm

        Bmajinpix       = Bmaj_as*pixsize_x/au*d_pc
        Bmininpix       = Bmin_as*pixsize_x/au*d_pc
        pixelinbeam     = (np.pi/4.0/np.log(2.0))*Bmajinpix*Bmininpix
        # print('plotting with the unit of [Jy/beam]')
        # print('Bmaj, Bmin = ',  Bmaj_as, ',', Bmin_as, 'in arcseconds')
        conv           = pixsize_x * pixsize_y / (d_pc*pc)**2. * 1e23 * pixelinbeam
        # print('ref: how many pixels in one beam?:', pixelinbeam)
        # im_Jypbeam     = np.zeros((im_nx,im_ny,nlam))
        # im_Jypbeam     = im.copy()*conv
        im_PI_forplot  = im_PI.copy()*conv
        # im_forplot     = im_Jypbeam
        # pltunit        = '[Jy/beam]'

        # Determine the plot index based on simulation count or other criteria

        if self.single_plot_var.get():
             self.master.after(100, lambda: self.create_normplot(data=im_PI_forplot, 
                                                   wavelength=wavelength_mic[0], 
                                                   size=sizex_au, 
                                                   inc=incl_value, 
                                                   pos=posang_value,
                                                   ))
        else:
             self.master.after(100, lambda: self.create_normplot3x3(data=im_PI_forplot, 
                                                   wavelength=wavelength_mic[0], 
                                                   size=sizex_au, 
                                                   inc=incl_value, 
                                                   pos=posang_value,
                                                   plot_index=plot_index))

    # Function to run the thermal simulation  
    def run_simulation_thermal(self):
                # Get values from input fields
        nr        = int(self.nr.get())
        nphi      = int(round(float(self.nphi.get())))
        ntheta    = int(round(float(self.ntheta.get())))
        rin       = int(round(float(self.rin.get()))) * au
        rout      = int(round(float(self.rout.get()))) * au 
        #INNER DISK
        sigmad0   = float(self.sigmad0.get())
        ringcau1  = float(self.ringcau1.get())
        elevD     = float(self.elevD.get())
        rotaD1    = float(self.rotaD1.get())

        #TESTING
        INNERASSYMETRIC = False 
        azipos1   = pi
        azimwid1  = 0.45

        #OUTER DISK
        sigmadc2  = float(self.sigmadc2.get())
        ringcau2  = float(self.ringcau2.get())
        ringwau2  = float(self.ringwau2.get())
        rotaD2    = float(self.rotaD2.get())

        OUTERASSYMETRIC = bool(self.OUTERASSYMETRIC.get())
        azimwid2  = float(self.azimwid2.get())


        ELLIPSE = bool(self.ELLIPSE.get())
        a_semi  = float(self.a_semi.get())  * au            # [cm] semi-major axis
        e       = float(self.e.get())                     # eccentricity
        
        nphot    = float(self.nphot.get())                  # Photones for Monte Carlo Simulation  

        wavelength  = float(self.wavelength.get())
        incl_value = float(self.incl_value.get())
        posang_value = float(self.posang_value.get())
        plot_index = int(self.index.get())

        ''' Make the r coordinates'''


        nlev_rin = 20                                       # Grid refinement at the inner edge: nr of cycles
        nspan_rin= 5                                        # Grid refinement at the inner edge: nr of cells each cycle
        ri       = np.logspace(np.log10(rin),np.log10(rout),nr+1)
        ri       = tools.grid_refine_inner_edge(ri,nlev_rin,nspan_rin)   # Refinement at inner edge
        rc       = 0.5 * ( ri[:-1] + ri[1:] )
        nr       = len(rc)           # Recompute nr, because of refinement at inner edge


        ''' Make the theta and phi coordinates '''


        thetaup  = 0.1
        thetai   = np.linspace(thetaup,np.pi-thetaup,ntheta+1)
        phii     = np.linspace(0.e0,np.pi*2.e0,nphi+1)
        # thetac   = 0.5 * ( thetai[0:ntheta] + thetai[1:ntheta+1] )
        # phic     = 0.5 * ( phii[0:nphi] + phii[1:nphi+1] )
        thetac   = 0.5 * ( thetai[:-1] + thetai[1:] )
        phic     = 0.5 * ( phii[:-1] + phii[1:] )


        ''' Make the 3-D-grid '''


        qq       = np.meshgrid(rc,thetac,phic,indexing='ij')
        r        = qq[0]   # The radial grid of the analytic disk model
        tt       = qq[1]
        pp       = qq[2]

        # ttdeg    = np.rad2deg(tt)
        # ppdeg    = np.rad2deg(pp)


        "Calculating each component of the coordinates (r, phi, theta) of the rotated frame that belongs to the inner disk"


        elevR    = np.deg2rad(elevD)        # Elevation of the disk in RAD
        newr     = r
        coselev  = (-np.sin(tt)*np.sin(pp)*np.sin(elevR))+(np.cos(tt)*np.cos(elevR))

        newtt       = np.arccos(coselev)
        sinnewpp    = np.sin(tt)*np.sin(pp)*np.cos(elevR)+np.cos(tt)*np.sin(elevR)/(np.sqrt(1-coselev**2))
        sinnewtt    = np.sin(newtt)
        cosnewpp    = (np.sin(tt)*np.cos(pp))/sinnewtt
        newpp       = np.arctan2(sinnewpp,cosnewpp)
        newpp       = np.where(newpp < 0, newpp + 2 * pi, newpp) #shifting atan2 domain from [-pi,pi] towards [0,2pi]


        '''Rearranging the assymmetrical shape'''


        dphi  =  newpp - azipos1

        # dphi = np.where(dphi >= pi, dphi + 2*pi, dphi)
        # dphi = np.where(dphi < 0, dphi + 2*pi, dphi)
        # dphi[dphi < 0] += 2 * pi 
        # if ((dphi>pi).all) or ((dphi<-pi).all):
        #      dphi += 2*pi

        ''' Make the dust density model '''

        #INNER DISK
        q_inner  = 0.25
        flang    = 0.1                         # The assumed constant radiative incidence angle
        lstar    = 4*pi*rstar**2*ss*tstar**4   # Stellar luminosity L = sigma*A*T^4
        firr     = flang*lstar/(4*pi*r**2)     # Irradiative flux F = L/A with A=4pi*r^2
        tmid     = (firr/ss)**q_inner          # Estimate of midplane temperature (^0.25)
        # ymid     = tmid_0 * (r/)
        cs1      = np.sqrt(kb*tmid/(2.3*mp))   # Isothermal sound speed at midplane
        omk      = np.sqrt(GG*mstar/r**3)      # The Kepler angular frequency
        hp1      = cs1/omk                      # The pressure scale height
        hpr1     = hp1/r                        # The dimensionless hp

        #OUTER DISK
        q_outer  = 0.25
        flang    = 0.1                         # The assumed constant radiative incidence angle
        lstar    = 4*pi*rstar**2*ss*tstar**4   # Stellar luminosity L = sigma*A*T^4
        firr     = flang*lstar/(4*pi*r**2)     # Irradiative flux F = L/A with A=4pi*r^2
        tmid2    = (firr/ss)**q_outer          # Estimate of midplane temperature (^0.25)

        cs2      = np.sqrt(kb*tmid2/(2.3*mp))   # Isothermal sound speed at midplane
        omk      = np.sqrt(GG*mstar/r**3)      # The Kepler angular frequency
        hp2      = cs2/omk                      # The pressure scale height
        hpr2     = hp2/r                        # The dimensionless hp

        ''' INNER DISK '''

        ringc1    = ringcau1*au                # Convert [au] to [cm]
        ringw1    = ringwau1*au                # Convert [au] to [cm]


        if INNERASSYMETRIC==True:
            sigmad    = sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0) * np.exp(-(((dphi)**2)/azimwid1**2)/2.0) # +  sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0) * np.exp(-(((newpp - 1.8*pi)**2)/azimwid1**2)/2.0)     # radial & azimuthal surface density profile for ring [g/cm^2]
        elif INNERASSYMETRIC==False:
            sigmad    = sigmad0 * np.exp(-((newr-ringc1)**2/ringw1**2)/2.0)                                                  # radial surface density profile for ring [g/cm^2]
        
        newzr     = pi/2.e0 - newtt 

        rhod     = ( sigmad / (np.sqrt(2.e0*np.pi)*hp1) ) * np.exp(-(newzr**2/hpr1**2)/2.e0)                                   # volume density

        rotation1 = int(rotaD1/(360/nphi))          
        rhod     = np.roll(rhod, rotation1, axis=2)


        ''' OUTER DISK '''


        ringc2 = ringcau2*au                      # Radius of the ring, Convert [au] to [cm] 

        if ELLIPSE == True:
            ringc2_mod = a_semi*(1-e**2)/(1-e*np.cos(pp))
        else:
            ringc2_mod = ringc2

        ringw2  = ringwau2*au                   # Radial width of the ring, Convert [au] to [cm]

        if OUTERASSYMETRIC==True:
            sigmad2   = sigmadc2  * np.exp(-((newr-ringc2_mod)**2/ringw2**2)/2.0) * np.exp(-((pp-1.05*pi)**2/azimwid2**2)/2.0) # radial & azimuthal surface density profile for ring [g/cm^2]
        else:
            sigmad2   = sigmadc2  * np.exp(-((newr-ringc2_mod)**2/ringw2**2)/2.0)                                     # radial surface density profile for ring [g/cm^2]

        zr        = np.pi/2.e0 - tt                                                      


        rhod2     = ( sigmad2 / (np.sqrt(2.e0*np.pi)*hp2) ) * np.exp(-(zr**2/hpr2**2)/2.e0)                                   # volume density

        rotation2 = int(rotaD2/(360/nphi))                                                                                  # only if assymetrical shape 
        rhod2     = np.roll(rhod2, rotation2, axis=2)

    
        ''' ADDING OUTER & INNER DISK '''

        rhod     += rhod2

        #
        # Write the grid file
        #
        with open('amr_grid.inp','w+') as f:
            f.write('1\n')                       # iformat
            f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
            f.write('100\n')                     # Coordinate system: spherical
            f.write('0\n')                       # gridinfo
            f.write('1 1 1\n')                   # Include r,theta coordinates
            f.write('%d %d %d\n'%(nr,ntheta,nphi))  # Size of grid
            for value in ri:
                f.write('%21.14e\n'%(value))      # X coordinates (cell walls)
            for value in thetai:
                f.write('%21.14e\n'%(value))      # Y coordinates (cell walls)
            for value in phii:
                f.write('%21.14e\n'%(value))      # Z coordinates (cell walls)
                

        #
        # Write the density file
        #
        with open('dust_density.inp','w+') as f:
            f.write('1\n')                       # Format number
            f.write('%d\n'%(nr*ntheta*nphi))     # Nr of cells
            f.write('1\n')                       # Nr of dust species
            data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%21.14e")
            f.write('\n')
    
        # ----------------------------------
        #
        # Wavelength [um] data and create wavelength_micron.inp file
        #
        # ----------------------------------

        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size

        #
        # Write the wavelength file
        #
        with open('wavelength_micron.inp','w+') as f:
            f.write('%d\n'%(nlam))
            for value in lam:
                f.write('%21.14e\n'%(value))


        #
        # Write the stars.inp file
        #
        with open('stars.inp','w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n'%(nlam))
            f.write('%21.14e %21.14e %21.14e %21.14e %21.14e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
            for value in lam:
                f.write('%21.14e\n'%(value))
            f.write('\n%21.14e\n'%(-tstar))
            

        #
        # Dust opacity control file
        #
        with open('dustopac.inp','w+') as f:
            f.write('2               Format number of this file\n')
            f.write('1               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('1              Way in which this dust species is read\n')
            f.write('0               0=Thermal grain\n')
            f.write('silicate         Extension of name of dustkappa_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')
        #
        # Write the radmc3d.inp control file
        #
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 1\n') 
            f.write('iranfreqmode = 1\n')

        self.loading_label.config(text="Simulation running... 💫, don't interrupt")
        self.master.update()

        self.running = True
        thread = threading.Thread(target=self.run_radmc3d)
        thread.start()

        # Getting input parameters for imaging data 

        cmd = f'radmc3d image lambda {wavelength} npix 200 sizeau 356 incl {incl_value} posang {posang_value}'
        os.system(cmd)

        fnameread      = 'image.out'

        with open(fnameread,'r') as f:
            s        = f.readline()
            im_n      = f.readline()
            im_nx     = int(im_n.split()[0])
            im_ny     = int(im_n.split()[1])
            nlam      = int(f.readline())
            pixsize   = f.readline()
            pixsize_x = float(pixsize.split()[0])
            pixsize_y = float(pixsize.split()[1])
            wavelength_mic = np.zeros(nlam)
            im    = np.zeros((im_nx,im_ny,nlam))
            rr_au = np.zeros((im_nx,im_ny,nlam))
            phi   = np.zeros((im_nx,im_ny,nlam))
            for ilam in range(nlam): wavelength_mic[ilam] = f.readline()
            for ilam in range(nlam):
                dummy = f.readline() #blank line
                for iy in range(im_ny):
                    for ix in range(im_nx):
                        image_temp = f.readline()
                        im[ix,iy,ilam] = image_temp
                        rr_au[ix,iy,ilam] = np.sqrt(((ix-im_nx/2)*pixsize_x)**2+((iy-im_ny/2)*pixsize_y)**2)/au
                        phi[ix,iy,ilam]   = np.arctan2((iy-im_ny/2)*pixsize_y,(ix-im_nx/2)*pixsize_x) 
        sizex_au = im_nx*pixsize_x/au
        # sizey_au = im_ny*pixsize_y/au

        # im_Jypbeam : intensity with a unit of [Jy/beam]
        # beam size: Full Width of Half Maximum of the beam
        # beam size = pi/(4 log(2)) * Bmaj * Bmin
        d_pc = 136   # For IRS48 first measured 121 pc, which got updated to 136pc
        Bmaj_as        = 0.11 # Major axis of the beam in arcsec
        Bmin_as        = 0.072 # Minor axis of the beam in arcsec
        Bmajinpix       = Bmaj_as*pixsize_x/au*d_pc
        Bmininpix       = Bmin_as*pixsize_x/au*d_pc
        pixelinbeam    = (np.pi/4.0/np.log(2.0))*Bmajinpix*Bmininpix
        conv           = pixsize_x * pixsize_y / (d_pc*pc)**2. * 1e23 * pixelinbeam
        im_Jypbeam     = np.zeros((im_nx,im_ny,nlam))
        im_Jypbeam     = im.copy()*conv
        im_forplot     = im_Jypbeam

        # Determine the plot index based on simulation count or other criteria

        if self.single_plot_var.get():
             self.master.after(100, lambda: self.create_normplot(data=im_forplot, 
                                                   wavelength=wavelength_mic[0], 
                                                   size=sizex_au, 
                                                   inc=incl_value, 
                                                   pos=posang_value,
                                                   ))
        else:
             self.master.after(100, lambda: self.create_normplot3x3(data=im_forplot, 
                                                   wavelength=wavelength_mic[0], 
                                                   size=sizex_au, 
                                                   inc=incl_value, 
                                                   pos=posang_value,
                                                   plot_index=plot_index))

    # Function to run the thermal RADMC3D simulation
    def run_radmc3d(self):
        try:
            os.system('radmc3d mctherm')
            self.loading_label.config(text="Simulation completed successfully! 🎉")
        except Exception as e:
            self.loading_label.config(text=f"Simulation failed! 😔 Error: {e}")


    '''Creating the 3x3/single plots of the simulation data'''


    # Function to create the 3x3 plot of the simulation data
    def create_normplot3x3(self,data,wavelength,size,inc,pos,plot_index):
        
        if  self.scattering_var.get() == 1:                  # Scattering
            interpolation  = 'bicubic'
            cmap           = 'nipy_spectral'
            title          = f'Polarized Intensity of {float(wavelength)}$\mu m$' 

            simu_data = data[:,:].T
            simu_data = np.where(simu_data<10**(-15),10**-15 , simu_data.copy())   # cutting-off noise 

            if self.mask_var.get() == 1:

                simu_data = tools.circular_mask(simu_data, radius=23)  # masking image only if checkbox is checked
                nor_simu_data = simu_data / np.max(simu_data, axis=1)[100]     # normalization of image data 

            else:
                test = simu_data.copy()
                test[100, 80:120] = np.nan
                nor_simu_data = simu_data / np.nanmax(test, axis=1)[100]     # normalization of image data
            
            fig = plt.figure(figsize=(3.2, 2.5))
            ax = fig.add_axes([0.2, 0.15, 0.75, 0.75])  # Adjust the position of the plot within the figure
            im = ax.imshow(nor_simu_data,interpolation=interpolation,cmap=cmap,norm=None,origin='lower',extent= [-size/2.0,size/2.0,-size/2.0,size/2.0])
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('X [AU]', fontsize=6)
            ax.set_ylabel('Y [AU]', fontsize=6)
            ax.text(x=0, y=-110, s=f"Inclination: {inc}°, Position Angle: {pos}°", fontsize=5, bbox=dict(facecolor='white', alpha=1), ha='center')
            cbar=plt.colorbar(im, location='right')
            cbar.ax.tick_params(labelsize=7)
            tools.trace_outer_ring(nor_simu_data,ax=None, grid_x= 2, grid_y=0,int_threshold=0.2)  

        else:                                         # Thermal emission
            from matplotlib.colors import LogNorm
            intmax         =  np.max(data)            # maximum value
            intmin         =  np.max(data)*0.01       # 1% of the maximum value       
            interpolation  = 'bicubic'
            cmap           = 'magma'
            title          = f'Spectral irradiance of {float(wavelength)}$\mu m$'
            norm=LogNorm(vmin=intmin, vmax=intmax)
            #norm=None
            simu_data_org = data[:,:,0].T 

            # Plotting the image
            fig = plt.figure(figsize=(3.2, 2.6))
            ax = fig.add_axes([0.18, 0.15, 0.75, 0.75])  # Adjust the position of the plot within the figure
            im = ax.imshow(simu_data_org,interpolation=interpolation,cmap=cmap,norm=norm,origin='lower',extent= [-size/2.0,size/2.0,-size/2.0,size/2.0])
            ax.tick_params(axis='both', which='major', labelsize=7)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('X [AU]', fontsize=6)
            ax.set_ylabel('Y [AU]', fontsize=6)
            ax.text(x=0, y=-155, s=f"Inclination: {inc}°, Position Angle: {pos}°", fontsize=5, bbox=dict(facecolor='white', alpha=1), ha='center')
            cbar=plt.colorbar(im, location='right')
            cbar.ax.tick_params(labelsize=7)
            cbar.ax.set_title('Jy/beam', fontsize=7, y=-0.15)  # Label below the colorbar
            # tools.trace_outer_ring(simu_data_org,ax=None, grid_x= 2, grid_y=0,int_threshold=0.2)  

        # Place the plot in the GUI window depending on the plot index
        if self.plot_canvases[plot_index]:
            self.plot_canvases[plot_index].get_tk_widget().destroy()  # Remove old canvas

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame.winfo_children()[plot_index])  
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.plot_canvases[plot_index] = canvas # Store the new canvas in the list

        # Remove the loading label
        self.loading_label.pack_forget()

    # Function to create the single plot of the simulation data
    def create_normplot(self,data,wavelength,size,inc,pos):
        
        if self.scattering_var.get() == 1:                  # Scattering
            interpolation  = 'bicubic'
            cmap           = 'nipy_spectral'
            title          = f'Polarized Intensity of {float(wavelength)}$\mu m$' 

            simu_data = data[:,:].T
            simu_data = np.where(simu_data<10**(-15),10**-15 , simu_data.copy())   # cutting-off noise

            if self.mask_var.get() == 1:
                simu_data = tools.circular_mask(simu_data, radius=23)  # masking image only if checkbox is checked
                nor_simu_data = simu_data / np.max(simu_data, axis=1)[100]     # normalization of image data 

            else:
                test = simu_data.copy()
                test[100, 80:120] = np.nan
                nor_simu_data = simu_data / np.nanmax(test, axis=1)[100]     # normalization of image data  

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # Adjust the position of the plot within the figure
            im = ax.imshow(nor_simu_data,interpolation=interpolation,cmap=cmap,norm=None,origin='lower',extent= [-size/2.0,size/2.0,-size/2.0,size/2.0])
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('X [AU]', fontsize=8)
            ax.set_ylabel('Y [AU]', fontsize=8)
            ax.text(x=0, y=-120, s=f"Inclination: {inc}°, Position Angle: {pos}°", fontsize=7, bbox=dict(facecolor='white', alpha=1), ha='center')

            plt.subplots_adjust(right=0.9)

            cbar=plt.colorbar(im, location='right')
            cbar.ax.tick_params(labelsize=8)
            tools.trace_outer_ring(nor_simu_data,ax=None, grid_x= 2, grid_y=0,int_threshold=0.2)
        
        else:                                   # Thermal emission
            from matplotlib.colors import LogNorm
            intmax         =  np.max(data)            # maximum value
            intmin         =  np.max(data)*0.01       # 1% of the maximum value
            interpolation  = 'bicubic'
            cmap           = 'magma'
            title          = f'Spectral irradiance of {float(wavelength)}$\mu m$'
            norm=LogNorm(vmin=intmin, vmax=intmax)
            #norm=None
            simu_data_org = data[:,:,0].T     

            # Plotting the image
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # Adjust the position of the plot within the figure
            im = ax.imshow(simu_data_org,interpolation=interpolation,cmap=cmap,norm=norm,origin='lower',extent= [-size/2.0,size/2.0,-size/2.0,size/2.0])
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('X [AU]', fontsize=8)
            ax.set_ylabel('Y [AU]', fontsize=8)
            ax.text(x=0, y=-120, s=f"Inclination: {inc}°, Position Angle: {pos}°", fontsize=7, bbox=dict(facecolor='white', alpha=1), ha='center')
            plt.subplots_adjust(right=0.9)
            cbar=plt.colorbar(im, location='right',label='Jy/beam')
            cbar.ax.tick_params(labelsize=8)
            cbar.ax.set_title('Jy/beam', fontsize=7, y=-0.15)  # Label below the colorbar
            # tools.trace_outer_ring(simu_data_org,ax=None, grid_x= 2, grid_y=0,int_threshold=0.2)  

        # Place plot in the GUI
        for canvas in self.plot_canvases:
            canvas.get_tk_widget().destroy()  # Remove old canvas
        self.plot_canvases = []  # Clear the list of canvases

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.plot_canvases.append(canvas)

    #Function to save the plot
    def save_plot(self):

        if self.single_plot_var.get() == 0:     # 3x3 grid mode
            # Create a new figure to hold all the subplots
            grid_fig = plt.figure(figsize=(12, 10))  # Adjust size as needed
            # plot_index = int(self.index.get())

            # Iterate through the plot canvases and add them to the grid figure
            for i, canvas in enumerate(self.plot_canvases):
                fig = canvas.figure  # Get the Matplotlib figure from the canvas

                ax = grid_fig.add_subplot(3, 3, i+1)
                if fig.axes and fig.axes[0].images:  # Check if axes and images exist
                    im = ax.imshow(fig.axes[0].images[0].get_array(), origin='lower',
                            cmap=fig.axes[0].images[0].get_cmap(),interpolation='bicubic',
                            extent=fig.axes[0].images[0].get_extent())  # Include extent for correct aspect ratio
                    ax.tick_params(axis='both', which='major', labelsize=5)

                    # Copy title, labels, and colorbar from the original figure
                    ax.set_title(fig.axes[0].get_title(), fontsize=8)
                    ax.set_xlabel(fig.axes[0].get_xlabel(), fontsize=6)
                    ax.set_ylabel(fig.axes[0].get_ylabel(), fontsize=6)

                    # Get the colorbar from the original figure and add it to the subplot
                    cbar = plt.colorbar(im, location='right').ax.tick_params(labelsize=5)
                    if fig.axes[0].images[0].get_cmap() == 'magma':
                        cbar.ax.set_title('Jy/beam', fontsize=7, y=-0.15)
                else:
                    # Plot a white canvas if the subplot is empty
                    ax.axis('off')
                    pass
                    

                    
            # Ask the user for the location and name of the file to save
            file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                    filetypes=[("PNG Image", "*.png"),
                                                                ("JPEG Image", "*.jpg"),
                                                                ("All Files", "*.*")])
            if file_path:
                grid_fig.savefig(file_path)  # Save the combined figure

        else:
            # Get the figure that contains the single plot                
            file_path = filedialog.asksaveasfilename(defaultextension=".png", 
                                                    filetypes=[("PNG Image", "*.png"), 
                                                            ("JPEG Image", "*.jpg"),
                                                            ("All Files", "*.*")])
            if file_path:
                try:
                    plt.savefig(file_path)  # Save the figure using Matplotlib
                except Exception as e:
                    # Handle exception (e.g., show an error message)
                    print(f"Error saving plot: {e}")


    # Function to close the GUI window
    def abort_and_close(self):
        plt.close('all')  # Close all open Matplotlib figures
        self.master.destroy()  # Close the GUI window


# Assuming root is your Tkinter root window
root = tk.Tk()
# Make the window open in fullscreen
root.attributes('-fullscreen', True)

# Make the window resizable
root.resizable(True, True)

gui = DiskModelGUI(root)
root.mainloop()