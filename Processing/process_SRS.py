import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker, cm
import scipy.stats
import glob
from scipy import interpolate
import pandas as pd
import matplotlib.patches as mpatches
from mms_nirs.UCLN import DefaultValues
from mms_nirs.utils import ExtinctionCoefficients, calc_dpf,calc_attenuation_slope
from mms_nirs.SRS import srs_values

main_path = "/home/caredda/Videos/PROSPEKT/"


class display_config:
    def __init__(self, id_skin=0, id_adipose=0, id_muscle=0, id_sat_m=0, id_sat_p=0, id_HbT_m=0, id_HbT_p=0, id_mel=0):
        self.id_skin = id_skin
        self.id_adipose = id_adipose
        self.id_muscle = id_muscle
        self.id_sat_p = id_sat_p
        self.id_sat_m = id_sat_m
        self.id_HbT_m = id_HbT_m
        self.id_HbT_p = id_HbT_p
        self.id_mel = id_mel


def calc_conc(att, optode_dist, ext_coeffs_inv, wl_dep, dpf):
    return (np.matmul(ext_coeffs_inv, (att / wl_dep).T)) / (dpf * optode_dist)




def get_scanning_proba(id_subject, HbT_muscle, HbT_placenta, device_name, int_time_s, binning, detectors):

    #Thickness array
    skin_thickness_array_mm = np.array([1.9,2.6,3, 1])
    adipose_thickness_array_mm = np.array([3.2,4.1,5.6, 2.3])
    muscle_thickness_array_mm_=  np.array([7.6,10,14, 6.7])

    #HbT
    C_HbT_muscle_array = np.array([15,25,35,50])
    C_HbT_placenta_array = np.array([15,25,35,50])

    #find id HbT
    id_HbT_m = np.where(C_HbT_muscle_array == HbT_muscle)[0].item()
    id_HbT_p = np.where(C_HbT_placenta_array == HbT_placenta)[0].item()

    #Load data
    filename = main_path+"simulations/scanning_proba/"+device_name+"_ti_"+str(int_time_s)+"_binning_"+str(binning)+".npz"

    data = np.load(filename)
    S_placenta = data['S_placenta']
    detection_proba = data['Detection_proba']
    detection_proba[detection_proba>1] = 1
    detector_signal = data['detector_signal']
    noise_level = data['noise_level']
    SD_separation_cm = data['SD_separation_cm']
    t_skin = np.round(data['t_skin']*10)
    t_adipose = np.round(data['t_adipose']*10)
    t_muscle = np.round(data['t_muscle']*10)

    filename = main_path+"simulations/detection_proba/detection_probability_"+device_name+"_TI_"+str(int_time_s)+"s_binning_"+str(binning)+"_camera_model_0.npz"
    data = np.load(filename)
    min_signal_det = data['mu_p_min']

    #Get id patient
    for i in range(t_skin.shape[0]):
        if( (t_skin[i] == np.round(skin_thickness_array_mm[id_subject-1])) and
            (t_adipose[i] == np.round(adipose_thickness_array_mm[id_subject-1])) and
            (t_muscle[i] == np.round(muscle_thickness_array_mm_[id_subject-1]))):
            _id_subject = i
            break

    #get detectors id
    id_comb = np.array([],dtype=int)
    for det in detectors:
        id_comb = np.append(id_comb,np.where(SD_separation_cm==det)[0].item())

    #get value for the corresponding patient
    S_placenta = S_placenta[id_comb,_id_subject,:,:,id_HbT_m,id_HbT_p]
    detection_proba = detection_proba[id_comb,_id_subject,:,:,id_HbT_m,id_HbT_p]
    detector_signal = detector_signal[id_comb,_id_subject,:,:,id_HbT_m,id_HbT_p]
    noise_level=noise_level[id_comb]
    min_signal_det = min_signal_det[id_comb]

    return S_placenta, detection_proba, detector_signal, noise_level, min_signal_det



def get_SRS_data_Homogeneous_vol(detectors):

    #SatO2
    SatO2 = np.array([0.4,0.6,0.8])

    #Source detector separation
    SD_separations_in_cm = np.array([1,2,3,4])

    #get detectors id
    id_comb = np.array([],dtype=int)
    for det in detectors:
        id_comb = np.append(id_comb,np.where(SD_separations_in_cm==det)[0].item())


    path = main_path + "simulations/output_SatO2_lookup_table/"
    data = scipy.io.loadmat(path +"cst.mat")


    #Wavelength
    wavelength = np.squeeze(data['Lambdas'])
    WAVELENGTHS = np.arange(780,900+1)


    # Get absorption and extinction coefficient
    extinction_coefficients = (
        ExtinctionCoefficients.set_index("wavelength")
        .loc[WAVELENGTHS]
        .reset_index()
        .values)

    # Compute SRS with these chromophores
    chrom=np.array(["HbO2","HHb"])

    #Get chromophore names
    chrom_name = np.asarray(ExtinctionCoefficients.columns, dtype='<U18')

    #Only keep selected chromophores
    keep_id = np.array([],dtype=int)
    for i in range(chrom.shape[0]):
        keep_id = np.append(keep_id,np.where(chrom[i]==chrom_name)[0][0])
    extinction_coefficients = extinction_coefficients[:,keep_id]

    # Get MBLL Matrix
    ext_coeffs_inv = np.linalg.pinv(extinction_coefficients)

    #Init output
    SatO2_array = np.zeros((SatO2.shape[0]))


    for S in range(SatO2.shape[0]):

        #Load mat file
        data = scipy.io.loadmat(main_path + "simulations/output_homo/out_homo_Sat_"+str(SatO2[S])+".mat")

        #Load intensity
        _I = data['Diffuse_reflectance']


        #interpolate I
        I = np.zeros((_I.shape[0],WAVELENGTHS.shape[0]))
        for i in range(I.shape[0]):
            I[i,:] = interpolate.interp1d(wavelength,_I[i,:], kind='linear')(WAVELENGTHS)

        #Compute Attenuation
        A = np.zeros(I.shape)
        for i in range(I.shape[0]):
            A[i,:] = np.log10(1/I[i,:])


        #Compute StO2 with SRS algorith
        _StO2 = process_SRS(A,SD_separations_in_cm,id_comb,WAVELENGTHS,ext_coeffs_inv)
        SatO2_array[S] = _StO2


    SatO2_array[SatO2_array>100] = np.nan
    SatO2_array[SatO2_array<0] = np.nan

    return SatO2_array




def get_SRS_data(id_subject, HbT_muscle, HbT_placenta, f_melanosome, detectors):


    #Source detector separation
    SD_separations_in_cm = np.array([1,2,3,4])

    #get detectors id
    id_comb = np.array([],dtype=int)
    for det in detectors:
        id_comb = np.append(id_comb,np.where(SD_separations_in_cm==det)[0].item())


    path = main_path + "simulations/output_SatO2_lookup_table/"
    data = scipy.io.loadmat(path +"cst.mat")

    #SatO2
    SatO2_muscle = np.array([0.4,0.6,0.8])
    SatO2_placenta = np.array([0.4,0.6,0.8])

    #Wavelength
    wavelength = np.squeeze(data['Lambdas'])
    WAVELENGTHS = np.arange(780,900+1)


    # Get absorption and extinction coefficient
    extinction_coefficients = (
        ExtinctionCoefficients.set_index("wavelength")
        .loc[WAVELENGTHS]
        .reset_index()
        .values
    )

    # Compute SRS with these chromophores
    chrom=np.array(["HbO2","HHb"])

    #Get chromophore names
    chrom_name = np.asarray(ExtinctionCoefficients.columns, dtype='<U18')

    #Only keep selected chromophores
    keep_id = np.array([],dtype=int)
    for i in range(chrom.shape[0]):
        keep_id = np.append(keep_id,np.where(chrom[i]==chrom_name)[0][0])
    extinction_coefficients = extinction_coefficients[:,keep_id]

    # Get MBLL Matrix
    ext_coeffs_inv = np.linalg.pinv(extinction_coefficients)


    #Init output
    SatO2_array = np.zeros((SatO2_muscle.shape[0],SatO2_placenta.shape[0]))


    A_array = np.zeros((SatO2_muscle.shape[0],SatO2_placenta.shape[0],id_comb.shape[0],WAVELENGTHS.shape[0]))




    for S_m in range(SatO2_muscle.shape[0]):
        for S_p in range(SatO2_placenta.shape[0]):

            #Load mat file
            data = scipy.io.loadmat(path + "out_Subject_"+str(id_subject)+ "_Sat_muscle_"+str(SatO2_muscle[S_m])+"_St_placenta_"+str(SatO2_placenta[S_p])+ "_HbT_muscle_umol_"+str(HbT_muscle)+ "_HbT_placenta_umol_"+str(HbT_placenta)+ "_fmel_"+str(f_melanosome)+".mat")

            #Load intensity
            _I = data['Diffuse_reflectance']


            #interpolate I
            I = np.zeros((_I.shape[0],WAVELENGTHS.shape[0]))
            for i in range(I.shape[0]):
                I[i,:] = interpolate.interp1d(wavelength,_I[i,:], kind='linear')(WAVELENGTHS)

            #Compute Attenuation
            A = np.zeros(I.shape)
            for i in range(I.shape[0]):
                A[i,:] = np.log10(1/I[i,:])


            A_array[S_m,S_p,:,:] = A[id_comb,:]

            #Compute StO2 with SRS algorith
            _StO2 = process_SRS(A,SD_separations_in_cm,id_comb,WAVELENGTHS,ext_coeffs_inv)
            SatO2_array[S_m,S_p] = _StO2


    SatO2_array[SatO2_array>100] = np.nan
    SatO2_array[SatO2_array<0] = np.nan

    return SatO2_array, A_array


def process_SRS(A,SD_separations,comb,WAVELENGTHS,ext_coeffs_inv):

    #Calculate attenuation slop
    Attenuation_slope = calc_attenuation_slope(np.expand_dims(A[comb,:],axis=1),SD_separations[comb])


    #Calculate Concentration and StO2 with SRS
    _C, _StO2, _k_mua = srs_values(np.squeeze(Attenuation_slope),
                                    WAVELENGTHS,
                                    ext_coeffs_inv,
                                    SD_separations[comb].min(),
                                    SD_separations[comb].max())


    # SD_txt = ""
    # for j in range(len(SD_separations[comb])):
    #     SD_txt=SD_txt+str(int(SD_separations[comb][j]))+" "
    # SD_txt = SD_txt[0:-1] #remove last space

    return _StO2

def plot_SRS(SRS_array,detectors,conf = display_config(), mode='HbT',title="$SatO_2$ (SRS)"):


    #Init title
    title += " - Source/Detector separations "
    for i in range(detectors.shape[0]):
        title += str(detectors[i])+" mm "



    #Colormaps
    # cmap = cm.plasma
    cmap = cm.jet
    alpha = 1
    nb_row = 1

    ft_label = 12
    ft_title = 14

    plt.rcParams.update({'font.size': ft_label})


    #Mode HbT
    if(mode == "HbT"):

        title += "\nThickness skin "+str(skin_thickness[conf.id_skin]) + "mm - Thickness adipose "+ str(adipose_thickness[conf.id_adipose])+ " mm - Thickness muscle "+str(muscle_thickness[conf.id_muscle])+" mm - $SatO_2$ muscle "+str(int(SatO2_muscle[conf.id_sat_m]*100))+ "% - $SatO_2$ placenta "+str(int(SatO2_placenta[conf.id_sat_p]*100))

        #meshgrid(y,x)
        x = C_HbT_muscle
        y = C_HbT_placenta

        #labels
        x_label = "$HbT$\nmuscle (µM)"
        y_label = "$HbT$\nplacenta (µM)"

        #Select proba
        StO2 = SRS_array[conf.id_skin,conf.id_sat_m, conf.id_sat_p, :,:,:]

    #Mode SatO2
    if(mode == "SatO2"):

        title += "\n Thickness skin "+str(skin_thickness[conf.id_skin]) + "mm - Thickness adipose "+ str(adipose_thickness[conf.id_adipose])+ " mm - Thickness muscle "+str(muscle_thickness[conf.id_muscle])+" mm - $HbT$ muscle "+str(C_HbT_muscle[conf.id_HbT_m])+ "$\mu Mol$ - $HbT$ placenta "+str(C_HbT_placenta[conf.id_HbT_p])+ "$\mu Mol$"

        #meshgrid(y,x)
        x = SatO2_muscle*100
        y = SatO2_placenta*100

        #labels
        x_label = "$SatO_2$\nmuscle (%)"
        y_label = "$SatO_2$\nplacenta (%)"


        #Select proba
        StO2 = SRS_array[conf.id_skin,:,:, conf.id_HbT_m, conf.id_HbT_p,:]






    #Interpolation
    x_interp = np.linspace(x[0],x[-1],100)
    y_interp = np.linspace(y[0],y[-1],100)
    T_y,T_x = np.meshgrid(y_interp,x_interp)

    StO2_array = []
    for i in range(StO2.shape[-1]):
        interp = interpolate.RegularGridInterpolator((x, y),StO2[:,:,i],bounds_error=False, fill_value=None)
        StO2_array.append(interp((T_x,T_y)))

    StO2_array = np.asarray(StO2_array)





    #Calculate SatO2 errors
    E_muscle = []
    E_placenta = []
    if mode == "HbT":
        for i in range(len(StO2_array)):
            E_muscle.append(np.abs(100*(StO2_array[i]-SatO2_muscle[conf.id_sat_m]*100)/(SatO2_muscle[conf.id_sat_m]*100)))
            E_placenta.append(np.abs(100*(StO2_array[i]-SatO2_placenta[conf.id_sat_p]*100)/(SatO2_placenta[conf.id_sat_p]*100)))
    if mode == "SatO2":
        for i in range(len(StO2_array)):
            E_muscle.append(np.abs(100*(StO2_array[i]-T_x)/T_x))
            E_placenta.append(np.abs(100*(StO2_array[i]-T_y)/T_y))


    E_muscle = np.asarray(E_muscle)
    E_placenta = np.asarray(E_placenta)


    # Plot
    fig = plt.figure()
    plt.suptitle(title,fontsize=ft_title)


    #Plot Measured StO2
    vmin = np.nanmin(StO2_array)
    vmax = np.nanmax(StO2_array)

    for i in range(StO2_array.shape[0]):
        plt.subplot(3,StO2_array.shape[0],i+1)
        plt.title("Melanosome volume fraction "+str(f_melanosome[i])+"\n$SatO_2$")

        im=plt.pcolor(T_x,T_y,StO2_array[i,:,:],cmap=cmap,vmin = vmin,vmax=vmax)

        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        plt.colorbar(im)

    #Plot Error for muscle
    vmin = np.nanmin(E_muscle)
    vmax = np.nanmax(E_muscle)

    for i in range(StO2_array.shape[0]):
        plt.subplot(3,StO2_array.shape[0],StO2_array.shape[0]+i+1)
        plt.title("Melanosome volume fraction "+str(f_melanosome[i])+"\n$E_{SatO_2}$ for muscle")

        im=plt.pcolor(T_x,T_y,E_muscle[i,:,:],cmap=cmap,vmin = vmin,vmax=vmax)

        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        plt.colorbar(im)

    #Plot Error for placenta
    vmin = np.nanmin(E_placenta)
    vmax = np.nanmax(E_placenta)

    for i in range(StO2_array.shape[0]):
        plt.subplot(3,StO2_array.shape[0],2*StO2_array.shape[0]+i+1)
        plt.title("Melanosome volume fraction "+str(f_melanosome[i])+"\n$E_{SatO_2}$ for placenta")

        im=plt.pcolor(T_x,T_y,E_placenta[i,:,:],cmap=cmap,vmin = vmin,vmax=vmax)

        plt.xlabel(x_label,fontsize=ft_label)
        plt.ylabel(y_label,fontsize=ft_label)
        plt.colorbar(im)
    plt.show()




## Process SRS arrays for all cases

#Source detector separation
SD_separations_in_cm = np.array([1,2,3,4])

# Get all possible combination of 2 and more detectors
# id_comb = []
# for i in range(2,SD_separations.shape[0]+1):
#     id_comb += list(itertools.combinations(np.arange(0,SD_separations.shape[0]), i))
id_comb = [np.array([1,2]),np.array([2,3])] #2-3 cm, 3-4 cm


path = main_path + "simulations/output_SatO2_lookup_table/"
data = scipy.io.loadmat(path +"cst.mat")

#Wavelength
wavelength = np.squeeze(data['Lambdas'])
WAVELENGTHS = np.arange(wavelength[1],wavelength[-1]+1)

# Get absorption and extinction coefficient
extinction_coefficients = (
    ExtinctionCoefficients.set_index("wavelength")
    .loc[WAVELENGTHS]
    .reset_index()
    .values
)

#SatO2
SatO2_muscle = np.array([0.4,0.6,0.8])
SatO2_placenta = np.array([0.4,0.6,0.8])

#HbT
C_HbT_muscle = np.array([15,25,35,50])
C_HbT_placenta = np.array([15,25,35,50])


#melanosome volume fraction
f_melanosome = np.array([0.0255,0.155,0.305])

#thicknesses (25%, 50%, 75% percentile)
skin_thickness = np.round(np.array([1.9,2.6,3]))
adipose_thickness = np.round(np.array([3.2,4.1,5.6]))
muscle_thickness = np.round(np.array([7.6,10,14]))

dist_to_placenta = skin_thickness + adipose_thickness + muscle_thickness



# Compute SRS with these chromophores
chrom=np.array(["HHb", "HbO2"])

#Get chromophore names
chrom_name = np.asarray(ExtinctionCoefficients.columns, dtype='<U18')

#Only keep selected chromophores
keep_id = np.array([],dtype=int)
for i in range(chrom.shape[0]):
    keep_id = np.append(keep_id,np.where(chrom[i]==chrom_name)[0][0])
extinction_coefficients = extinction_coefficients[:,keep_id]
_extinction_coefficients = extinction_coefficients.copy()
_extinction_coefficients[:,0] = extinction_coefficients[:,1]
_extinction_coefficients[:,1] = extinction_coefficients[:,0]




# Get MBLL Matrix
ext_coeffs_inv = np.linalg.pinv(_extinction_coefficients)



#Init output
SatO2_array = np.zeros((adipose_thickness.shape[0], len(id_comb),
                        SatO2_muscle.shape[0],SatO2_placenta.shape[0],
                        C_HbT_muscle.shape[0],C_HbT_placenta.shape[0],f_melanosome.shape[0]))







for subject_id in range(SatO2_array.shape[0]):
    for S_m in range(SatO2_muscle.shape[0]):
        for S_p in range(SatO2_placenta.shape[0]):
            for HbT_m in range(C_HbT_muscle.shape[0]):
                for HbT_p in range(C_HbT_placenta.shape[0]):
                    for fmel in range(f_melanosome.shape[0]):
                        #Load mat file
                        data = scipy.io.loadmat(path + "out_Subject_"+str(subject_id+1)+ "_Sat_muscle_"+str(SatO2_muscle[S_m])+"_St_placenta_"+str(SatO2_placenta[S_p])+ "_HbT_muscle_umol_"+str(C_HbT_muscle[HbT_m])+ "_HbT_placenta_umol_"+str(C_HbT_placenta[HbT_p])+ "_fmel_"+str(f_melanosome[fmel])+".mat")

                         #Load intensity
                        _I = data['Diffuse_reflectance']


                        #interpolate I
                        I = np.zeros((_I.shape[0],WAVELENGTHS.shape[0]))
                        for i in range(I.shape[0]):
                            I[i,:] = interpolate.interp1d(wavelength,_I[i,:], kind='linear')(WAVELENGTHS)

                        #Compute Attenuation
                        A = np.zeros(I.shape)
                        for i in range(I.shape[0]):
                            A[i,:] = np.log10(1/I[i,:])


                        #Compute StO2 with SRS algorithm
                        for i in range(len(id_comb)):
                            _StO2 = process_SRS(A,SD_separations_in_cm,id_comb[i],WAVELENGTHS,ext_coeffs_inv)
                            SatO2_array[subject_id,i,S_m,S_p,HbT_m,HbT_p,fmel] = _StO2


SatO2_array[SatO2_array>100] = np.nan
SatO2_array[SatO2_array<0] = np.nan


## Plot SRS


#Source detector separation
SD_separations_in_cm = np.array([1,2,3,4])
id_comb = [np.array([1,2]),np.array([2,3])] #2-3 cm, 3-4 cm

plt.close('all')
id_SD = 0
detectors = SD_separations_in_cm[id_comb[id_SD]]*10

config = display_config()
# Thickness 25% (0),  50 % (1), 75 % percentile (2)
config.id_adipose = 0
config.id_skin = 0
config.id_muscle = 0

config.id_HbT_m = 1
config.id_HbT_p = 1
config.id_sat_p = 2
config.id_sat_m = 1


plot_SRS(SatO2_array[:,id_SD,:,:,:,:,:],detectors,config, mode='HbT')
# plot_SRS(SatO2_array[:,id_SD,:,:,:,:,:],detectors,config, mode='SatO2')


## Check Detla C values homogeneous phantom

detectors = np.array([3,4])
# detectors = np.array([2,3])


#SatO2
SatO2 = np.array([0.4,0.6,0.8])

#Source detector separation
SD_separations_in_cm = np.array([1,2,3,4])

#get detectors id
id_comb = np.array([],dtype=int)
for det in detectors:
    id_comb = np.append(id_comb,np.where(SD_separations_in_cm==det)[0].item())


path = main_path + "simulations/output_SatO2_lookup_table/"
data = scipy.io.loadmat(path +"cst.mat")


#Wavelength
wavelength = np.squeeze(data['Lambdas'])
WAVELENGTHS = np.arange(780,900+1)


# Get absorption and extinction coefficient
extinction_coefficients = (
    ExtinctionCoefficients.set_index("wavelength")
    .loc[WAVELENGTHS]
    .reset_index()
    .values)

# Compute SRS with these chromophores
chrom=np.array(["HbO2","HHb"])

#Get chromophore names
chrom_name = np.asarray(ExtinctionCoefficients.columns, dtype='<U18')

#Only keep selected chromophores
keep_id = np.array([],dtype=int)
for i in range(chrom.shape[0]):
    keep_id = np.append(keep_id,np.where(chrom[i]==chrom_name)[0][0])
extinction_coefficients = extinction_coefficients[:,keep_id]


# Get MBLL Matrix
ext_coeffs_inv = np.linalg.pinv(extinction_coefficients)

#Init output
SatO2_array = np.zeros((SatO2.shape[0]))
A_array = np.zeros((SatO2.shape[0],id_comb.shape[0],WAVELENGTHS.shape[0]))

for S in range(SatO2.shape[0]):

    #Load mat file
    data = scipy.io.loadmat(main_path + "simulations/output_homo/out_homo_Sat_"+str(SatO2[S])+".mat")

    #Load intensity
    _I = data['Diffuse_reflectance']


    #interpolate I
    I = np.zeros((_I.shape[0],WAVELENGTHS.shape[0]))
    for i in range(I.shape[0]):
        I[i,:] = interpolate.interp1d(wavelength,_I[i,:], kind='linear')(WAVELENGTHS)

    #Compute Attenuation
    A = np.zeros(I.shape)
    for i in range(I.shape[0]):
        A[i,:] = np.log10(1/I[i,:])

    A_array[S,:,:] = A[id_comb,:]

    #Compute StO2 with SRS algorith
    _StO2 = process_SRS(A,SD_separations_in_cm,id_comb,WAVELENGTHS,ext_coeffs_inv)
    SatO2_array[S] = _StO2


SatO2_array[SatO2_array>100] = np.nan
SatO2_array[SatO2_array<0] = np.nan

print(SatO2_array)



# Delta C

I = np.zeros((SatO2.shape[0],detectors.shape[0],WAVELENGTHS.shape[0]))

for S in range(SatO2.shape[0]):

    #Load mat file
    data = scipy.io.loadmat(main_path + "simulations/output_homo/out_homo_Sat_"+str(SatO2[S])+".mat")

    #Load intensity
    _I = data['Diffuse_reflectance']



    #interpolate I
    I_t = np.zeros((_I.shape[0],WAVELENGTHS.shape[0]))

    for i in range(I_t.shape[0]):
        I_t[i,:] = interpolate.interp1d(wavelength,_I[i,:], kind='linear')(WAVELENGTHS)
    I[S,:,:] = I_t[id_comb,:]


#Mean path
Mean_path = np.zeros((detectors.shape[0],WAVELENGTHS.shape[0]))
data = scipy.io.loadmat(main_path + "simulations/output_homo/out_homo_Sat_"+str(SatO2[0])+".mat")
mp = data['Mean_path']
#interpolate
mp_t = np.zeros((mp.shape[0],WAVELENGTHS.shape[0]))

for i in range(mp_t.shape[0]):
    mp_t[i,:] = interpolate.interp1d(wavelength,mp[i,:], kind='linear')(WAVELENGTHS)
Mean_path[:,:] = mp_t[id_comb,:]

#Delta A
dA = np.zeros(I.shape)
for S in range(SatO2.shape[0]):
    dA[S,:,:] = np.log10(np.divide(I[0,:,:],I[S,:,:]))



# Extinction coefficients
extinction_coefficients = DefaultValues().extinction_coefficients[:,0:2]

ext_coeffs_inv = np.linalg.pinv(extinction_coefficients)
dpf = 4.16
wl_dep = DefaultValues().wavelength_dependency

C_UCLN = []
for i in range(detectors.shape[0]):
    conc = calc_conc(dA[:,i,:], detectors[i], ext_coeffs_inv, wl_dep, dpf).T
    C_UCLN.append(conc*1000)

C_UCLN = np.asarray(C_UCLN)

# E_inv = get_MBLL_Matrix(extinction_coefficients,np.ones((extinction_coefficients.shape[0],)))
# dC = dA @ E_inv.T


#expected Delta C
C_HbO2 = 25 *SatO2
C_Hb = 25 *(1-SatO2)

expected_dC_HbO2 = C_HbO2 - C_HbO2[0]
expected_dC_Hb = C_Hb - C_Hb[0]



plt.close('all')
# plt.figure()
# plt.title("Ectinction coefficients")
# plt.plot(WAVELENGTHS,_extinction_coefficients[:,0],'r')
# plt.plot(WAVELENGTHS,_extinction_coefficients[:,1],'b')
# plt.show()



plt.figure()
plt.subplot(121)
plt.plot(A_array[:,0,:].T)
plt.subplot(122)
plt.plot(A_array[:,1,:].T)
plt.show()

plt.figure()
plt.title("$SatO_2$ measurements with SRS (homogeneous media)")
plt.plot(SatO2_array,label="mesured")
plt.plot(100*SatO2,label="expected")
plt.legend(loc="best")
plt.show()


plt.figure()
plt.suptitle("Concentration changes (homogeneous media)")
plt.subplot(121)
plt.title('$\Delta C_{HbO_2}$')
id_chrom = 0
id=0
plt.plot(C_UCLN[id,:,id_chrom],label="Detector "+str(detectors[id])+" cm")
id=1
plt.plot(C_UCLN[id,:,id_chrom],label="Detector "+str(detectors[id])+" cm")
plt.plot(expected_dC_HbO2,'k')

plt.subplot(122)
plt.title('$\Delta C_{Hb}$')
id_chrom = 1
id=0
plt.plot(C_UCLN[id,:,id_chrom],label="Detector "+str(detectors[id])+" cm")
id=1
plt.plot(C_UCLN[id,:,id_chrom],label="Detector "+str(detectors[id])+" cm")
plt.plot(expected_dC_Hb,'k')

plt.show()




## Get SRS values for fixed values (HbT, fmel and thickness)


#id_subject Thickness 25% (1),  50 % (2), 75 % percentile (3), 10 % percentile (4)
id_subject = 2
HbT_muscle = 25
HbT_placenta = 25
f_melanosome = 0.155

device_name="CYRIL"
# device_name="Mini CYRIL"

if device_name =="CYRIL":
    detectors = [np.array([2,3]), np.array([3,4])]
else:
    detectors = [np.array([3,4])]


#SatO2 array
SatO2_muscle_array = np.array([0.4,0.6,0.8])
SatO2_placenta_array = np.array([0.4,0.6,0.8])

#SatO2 array shape
#shape 0: SatO2 muscle
#shape 1: SatO2 placenta
SatO2_array = []
SatO2_array_homo = []
_SatO2_array = []
_SatO2_array_homo = []
A_array = []


for d in detectors:
    s, a = get_SRS_data(id_subject=id_subject,
                            HbT_muscle=HbT_muscle,
                            HbT_placenta=HbT_placenta,
                            f_melanosome=f_melanosome,
                            detectors=d)
    SatO2_array.append(s)
    A_array.append(a)

    SatO2_array_homo.append(get_SRS_data_Homogeneous_vol(detectors=d))

#Calculate Delta StO2
Delta_SatO2_muscle = np.array([])
Delta_SatO2_placenta = np.array([])
SatO2_muscle = np.array([])
SatO2_placenta = np.array([])

Delta_SatO2_SRS = []
for sat_m in range(SatO2_muscle_array.shape[0]):
    for sat_p in range(SatO2_placenta_array.shape[0]):
        Delta_SatO2_muscle = np.append(Delta_SatO2_muscle,SatO2_muscle_array[sat_m] - SatO2_muscle_array[0])
        Delta_SatO2_placenta = np.append(Delta_SatO2_placenta,SatO2_placenta_array[sat_p] - SatO2_placenta_array[0])

        SatO2_muscle = np.append(SatO2_muscle,SatO2_muscle_array[sat_m])
        SatO2_placenta = np.append(SatO2_placenta,SatO2_placenta_array[sat_p])

for i in range(len(SatO2_array)):
    temp = np.array([])
    temp2 = np.array([])

    for sat_m in range(SatO2_muscle_array.shape[0]):
        for sat_p in range(SatO2_placenta_array.shape[0]):
            temp = np.append(temp,SatO2_array[i][sat_m,sat_p] - SatO2_array[i][0,0])
            temp2 = np.append(temp2,SatO2_array[i][sat_m,sat_p])

    Delta_SatO2_SRS.append(temp)
    _SatO2_array.append(temp2)



#get scanning proba, detection proba and detector signal
det = np.array([],dtype=int)
for d in detectors:
   det = np.append(det,d)
det = np.unique(det)

Scanning_placenta, detection_proba, detector_signal, noise_level, min_signal_det = get_scanning_proba(id_subject=id_subject,
                                                                                      HbT_muscle= HbT_muscle,
                                                                                      HbT_placenta=HbT_placenta,
                                                                                      device_name=device_name,
                                                                                      int_time_s=1,
                                                                                      binning=10,
                                                                                      detectors=det)

#Calculate the contribution of placenta in the collected signal
prop_placenta_detector_signal = np.multiply(detector_signal,Scanning_placenta)

#Check if the contribution of placenta is higher than noise level
proba_contrib_placenta = np.zeros(prop_placenta_detector_signal.shape) #det, SatO2 muscle, SatO2 placenta

for d in range(proba_contrib_placenta.shape[0]):
    for id_sat_m in range(proba_contrib_placenta.shape[1]):
        for id_sat_p in range(proba_contrib_placenta.shape[2]):
            sample = np.random.normal(loc = prop_placenta_detector_signal[d,id_sat_m,id_sat_p],
                                    scale = noise_level[d],
                                    size = (100))

            # Process T-tests
            t,proba_contrib_placenta[d,id_sat_m,id_sat_p] = scipy.stats.ttest_1samp(sample,min_signal_det[d],
                                        alternative='greater')


